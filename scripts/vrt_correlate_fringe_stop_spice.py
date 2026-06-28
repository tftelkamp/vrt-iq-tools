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

import spiceypy as spice
import spiceypy

import zmq
from argparse import ArgumentParser

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from antennas import ANTENNA_TABLE

parser = ArgumentParser(description="Fringe stopper for vrt_correlate")
parser.add_argument('-l', '--list', action='store_true', help='List all available antennas')
parser.add_argument("-k", "--kernelfile", help="SPICE kernel file to use")
parser.add_argument("-d", "--kerneldir", help="SPICE default kernels directory", default="~/kernels")
args = parser.parse_args()

if args.list:
	print("Available antennas:")
	for ant_name in sorted(ANTENNA_TABLE.keys()):
	    print(f"  - {ant_name}")
	exit(0)

kernel_dir = args.kerneldir
kernel_file = args.kernelfile

# SPICE kernels
spiceypy.furnsh(kernel_dir+"/naif0012.tls")
spiceypy.furnsh(kernel_dir+"/de440s.bsp")
spiceypy.furnsh(kernel_dir+"/pck00011.tpc")
spiceypy.furnsh(kernel_dir+"/earth_latest_high_prec.bpc")
# spiceypy.furnsh(kernel_dir+"/mar097.bsp")
# spiceypy.furnsh(kernel_dir+"/jup346.bsp")
# spiceypy.furnsh(kernel_dir+"/dwingeloo.fk")
# spiceypy.furnsh(kernel_dir+"/dwingeloo.bsp")

if (kernel_file):
    spiceypy.furnsh(kernel_file)

port = 70001
context = zmq.Context()
socket = context.socket(zmq.ROUTER)
socket.bind ("tcp://*:{0}".format(port))

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

		[dim, radii] = spiceypy.bodvrd("EARTH", "RADII", 3)
		flattening = (radii[0]-radii[2])/radii[0];

		obs1 = spiceypy.georec(earth_loc1.lon.rad, earth_loc1.lat.rad, earth_loc1.height.value/1000, radii[0], flattening)
		obs2 = spiceypy.georec(earth_loc2.lon.rad, earth_loc2.lat.rad, earth_loc2.height.value/1000, radii[0], flattening)

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

		et = spiceypy.str2et( str(timestamp.isot) )
		try:
			[azlsta,lt] = spiceypy.azlcpo("ELLIPSOID", object_name, et, "CN+S", False, True, obs1, "EARTH", "ITRF93")
		except:
			print("Loading " + object_name + " failed")
			print("Not initialized!");
			to_send = "{0:d}".format(1);
			socket.send_multipart([identity, to_send.encode()])
			continue
		# [state, lt] = spiceypy.spkcpo(object_name, et, "J2000", "OBSERVER", "CN+S", obs1, "EARTH", "ITRF93")
		v1 = 1000*azlsta[3]
		d1 = 1000*azlsta[0]

		[azlsta,lt] = spiceypy.azlcpo("ELLIPSOID", object_name, et, "CN+S", False, True, obs2, "EARTH", "ITRF93")
		# [state, lt] = spiceypy.spkcpo(object_name, et, "J2000", "OBSERVER", "CN+S", obs2, "EARTH", "ITRF93")
		v2 = 1000*azlsta[3]
		d2 = 1000*azlsta[0]

		baseline_projected_w = d2 - d1
		baseline_projected_w_dot = v2 - v1

		baseline_projected_v = 0
		baseline_projected_u = 0

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


