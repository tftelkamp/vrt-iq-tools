#!/usr/bin/env python3

# Copyright 2026 by Thomas Telkamp
#
# SPDX-License-Identifier: MIT

import sys
import os
import time
import numpy as np

from astropy.time import Time
from astropy.table import Table
from astropy.coordinates import AltAz, ICRS, ITRS, GCRS
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
from astropy.units import Quantity
import astropy.constants

import zmq
from argparse import ArgumentParser

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from antennas import ANTENNA_TABLE

parser = ArgumentParser(description="Fringe stopper for vrt_correlate")
parser.add_argument('-l', '--list', action='store_true', help='List all available antennas')

args = parser.parse_args()

if args.list:
	print("Available antennas:")
	for ant_name in sorted(ANTENNA_TABLE.keys()):
	    print(f"  - {ant_name}")
	exit(0)

port = 70001
context = zmq.Context()
socket = context.socket(zmq.ROUTER)
socket.bind("tcp://*:{0}".format(port))

initialized = False;

while True:
	parts = socket.recv_multipart()
	identity = parts[0]
	message = parts[1]

	received_fields = message.rstrip().split()

	print(received_fields)

	type = int(received_fields[0])

	if (type == 0):

		ant1_name = (received_fields[1]).decode("utf-8").rstrip()
		ant2_name = (received_fields[2]).decode("utf-8").rstrip()
		object_name = b' '.join(received_fields[3:]).decode("utf-8").rstrip()

		if ant1_name in ANTENNA_TABLE:
		    earth_loc1 = ANTENNA_TABLE[ant1_name]
		    ant1_xyz = earth_loc1.get_itrs()
		else:
			print(f"Error: Location '{ant1_name}' not found in the table.")
			print(f"\nAvailable locations:")
			for ant_name in sorted(ANTENNA_TABLE.keys()):
			    print(f"  - {ant_name}")
			to_send = "{0:d}".format(1);
			socket.send_multipart([identity, to_send.encode()])
			continue

		if ant2_name in ANTENNA_TABLE:
		    earth_loc2 = ANTENNA_TABLE[ant2_name]
		    ant2_xyz = earth_loc2.get_itrs()
		else:
			print(f"Error: Location '{ant2_name}' not found in the table.")
			print(f"\nAvailable locations:")
			for ant_name in sorted(ANTENNA_TABLE.keys()):
			    print(f"  - {ant_name}")
			to_send = "{0:d}".format(1);
			socket.send_multipart([identity, to_send.encode()])
			continue

		try:
			source = SkyCoord.from_name(object_name)
			print(source)
		except:
			print("Loading " + object_name + " failed")
			print("Not initialized!");
			to_send = "{0:d}".format(1);
			socket.send_multipart([identity, to_send.encode()])
			continue

		baseline = EarthLocation(ant1_xyz.x - ant2_xyz.x, ant1_xyz.y - ant2_xyz.y, ant1_xyz.z - ant2_xyz.z)

		# this north/east logic is copied from destevez
		north_radec = [source.ra.deg, source.dec.deg + 90]
		if north_radec[1] > 90:
		    north_radec[1] = 180 - north_radec[1]
		    north_radec[0] = 180 + north_radec[0]
		north = SkyCoord(ra = north_radec[0]*u.deg, dec = north_radec[1]*u.deg)

		baseline_itrs = baseline.get_itrs().cartesian.xyz

		initialized = True;

	elif (type == 1):

		if (not initialized):
			print("Not initialized!");
			to_send = "{0:d}".format(1);
			socket.send_multipart([identity, to_send.encode()])
			continue

		int_seconds = received_fields[1];
		frac_seconds = received_fields[2];

		unix_time = float(int_seconds)+float(frac_seconds)/1e12
		fringe_time = Time(unix_time, format='unix');

		print("Received request for:", fringe_time.isot)

		source_itrs = source.transform_to( ITRS(obstime = Time(fringe_time))).cartesian
		north_itrs = north.transform_to(ITRS(obstime = Time(fringe_time))).cartesian
		east_itrs = north_itrs.cross(source_itrs)

		baseline_projected_w = source_itrs.xyz.T.dot(baseline_itrs).value
		baseline_projected_w_min1sec = source.transform_to(ITRS(obstime = Time(fringe_time)-1*u.s)).cartesian.xyz.T.dot(baseline_itrs).value
		baseline_projected_w_plus1sec = source.transform_to(ITRS(obstime = Time(fringe_time)+1*u.s)).cartesian.xyz.T.dot(baseline_itrs).value
		baseline_projected_w_dot = (baseline_projected_w_plus1sec - baseline_projected_w_min1sec)/2.0

		baseline_projected_v = north_itrs.xyz.T.dot(baseline_itrs).value
		baseline_projected_u = east_itrs.xyz.T.dot(baseline_itrs).value

		to_send = "{0:d},{1:d},{2:d},{3:.12e},{4:.12e},{5:.12e},{6:.12e}".format(
			0,
			int(int_seconds),
			int(frac_seconds),
			baseline_projected_w,
			baseline_projected_w_dot,
			baseline_projected_u,
			baseline_projected_v
		)

		print("Responding:",to_send)

		socket.send_multipart([identity, to_send.encode()])


