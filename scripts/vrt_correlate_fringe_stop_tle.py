#!/usr/bin/env python3

# Copyright 2026 by Thomas Telkamp
#
# SPDX-License-Identifier: MIT

import sys
import os
import time

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import numpy as np

from astropy.time import Time
from astropy.table import Table
from astropy.coordinates import AltAz, ICRS, ITRS, GCRS
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
from astropy.units import Quantity
import astropy.constants

# Skyfield
from skyfield.api import EarthSatellite, load, wgs84, utc
from scipy.constants import c
from skyfield.api import utc

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
socket.bind ("tcp://*:{0}".format(port))

initialized = False;

while True:
	#  Wait for next request from client
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

		# Locations
		site1 = wgs84.latlon(earth_loc1.lat.deg, earth_loc1.lon.deg, earth_loc1.height.value)
		site2 = wgs84.latlon(earth_loc2.lat.deg, earth_loc2.lon.deg, earth_loc2.height.value)

		ts = load.timescale()

		print(object_name)

		url = 'https://celestrak.org/NORAD/elements/gp.php?FORMAT=TLE&CATNR={}'.format(object_name)
		filename = '/tmp/tle-CATNR-{}.txt'.format(object_name)
		try:
			satellites = load.tle_file(url, filename=filename)
		except:
			print("Loading " + object_name + " failed")
			print("Not initialized!");
			to_send = "{0:d}".format(1);
			socket.send_multipart([identity, to_send.encode()])
			continue

		print(satellites)

		source_sat = satellites[0]

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
		print(unix_time)
		fringe_time = Time(unix_time, format='unix');

		timestamp = Time(Time(fringe_time, precision=6).isot, precision=6)
		print("Received request for:", timestamp.isot)
		
		skyfield_ts = ts.from_astropy(fringe_time)

		state1 = (source_sat - site1).at(skyfield_ts).frame_latlon_and_rates(site1)
		state2 = (source_sat - site2).at(skyfield_ts).frame_latlon_and_rates(site2)

		range_site1 = state1[-4].m
		range_site2 = state2[-4].m

		v_site1 = state1[-1].m_per_s
		v_site2 = state2[-1].m_per_s

		baseline_projected_w = range_site2 - range_site1
		baseline_projected_w_dot = v_site2 - v_site1

		baseline_projected_v = 0; #north_itrs.xyz.T.dot(baseline_itrs).value
		baseline_projected_u = 0; #east_itrs.xyz.T.dot(baseline_itrs).value

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


