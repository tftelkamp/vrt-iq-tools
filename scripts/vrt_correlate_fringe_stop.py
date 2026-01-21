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


ANTENNA_TABLE = {
    'dwingeloo': EarthLocation(lat=52.81213723180477, lon=6.396346463227839, height=70.26, ellipsoid="WGS84"),
    'stockert': EarthLocation(lat=50.56944039751571, lon=6.721943350231514, height=434, ellipsoid="WGS84"),
    'stockert_2m': EarthLocation(lat=50.569511621782546, lon=6.723179139703107, height=434, ellipsoid="WGS84"),
    'stockert_1m': EarthLocation(lat=50.56933956889883, lon=6.723980412727782, height=434, ellipsoid="WGS84"),
    'bochum': EarthLocation(lat=51.426990, lon=7.192566, height=205, ellipsoid="WGS84"),
    'debilt': EarthLocation(lat=52.1018688270963, lon=5.17973411317635, height=51, ellipsoid="WGS84"),
    'delft': EarthLocation(lat=51.999153, lon=4.373352, height=92, ellipsoid="WGS84"), # estimation
    'ata_1a': EarthLocation.from_geocentric(-2524036.0307912203*u.m, -4123528.101172219*u.m, 4147706.408318585*u.m),
    'ata_4j': EarthLocation.from_geocentric(-2523898.1150373477*u.m, -4123456.314794732*u.m, 4147860.3045849088*u.m)
}


cable_delay = 0
clock_offset = 0

parser = ArgumentParser(description="Fringe stopper for vrt_correlate")
parser.add_argument("-d", "--delay", help="fixed delay (ant1-ant2)", default=0, type=float)
parser.add_argument("-c", "--clock", help="clock offset (ant1-ant2)", default=0, type=float)
parser.add_argument('-l', '--list', action='store_true', help='List all available antennas')
parser.add_argument("-a1", "--ant1", help="antenna 1", type=str, required = True)
parser.add_argument("-a2", "--ant2", help="antenna 2", type=str, required = True)
parser.add_argument("-n", "--name", help="object name", type=str, required = True)

args = parser.parse_args()

if args.list:
	print("Available antennas:")
	for ant_name in sorted(ANTENNA_TABLE.keys()):
	    print(f"  - {ant_name}")
	exit(0)

ant1_name = args.ant1 #.lower()
ant2_name = args.ant2 #.lower()

if ant1_name in ANTENNA_TABLE:
    earth_loc1 = ANTENNA_TABLE[ant1_name]
    ant1_xyz = earth_loc1.get_itrs()
else:
    print(f"Error: Location '{args.ant1}' not found in the table.")
    print(f"\nAvailable locations:")
    for ant_name in sorted(ANTENNA_TABLE.keys()):
        print(f"  - {ant_name}")
    exit(1)

if ant2_name in ANTENNA_TABLE:
    earth_loc2 = ANTENNA_TABLE[ant2_name]
    ant2_xyz = earth_loc2.get_itrs()
else:
    print(f"Error: Location '{args.ant2}' not found in the table.")
    print(f"\nAvailable locations:")
    for ant_name in sorted(ANTENNA_TABLE.keys()):
        print(f"  - {ant_name}")
    exit(1)

cable_delay = args.delay
clock_offset = args.clock

source = SkyCoord.from_name(args.name)

baseline = EarthLocation(ant1_xyz.x - ant2_xyz.x, ant1_xyz.y - ant2_xyz.y, ant1_xyz.z - ant2_xyz.z)

print(source)

# this north/east logic is copied from destevez
north_radec = [source.ra.deg, source.dec.deg + 90]
if north_radec[1] > 90:
    north_radec[1] = 180 - north_radec[1]
    north_radec[0] = 180 + north_radec[0]
north = SkyCoord(ra = north_radec[0]*u.deg, dec = north_radec[1]*u.deg)

baseline_itrs = baseline.get_itrs().cartesian.xyz

port = 70001
context = zmq.Context()
socket = context.socket(zmq.REP)
socket.bind ("tcp://*:{0}".format(port))

while True:
	#  Wait for next request from client
	message = socket.recv()
	[int_seconds, frac_seconds] = message.rstrip().split()

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

	to_send = "{0:d},{1:d},{2:.12e},{3:.12e},{4:e},{5:e},{6:.12e},{7:.12e}".format(
		int(int_seconds),
		int(frac_seconds),
		baseline_projected_w,
		baseline_projected_w_dot,
		float(cable_delay),
		float(clock_offset),
		baseline_projected_u,
		baseline_projected_v
	)

	print("Responding:",to_send)

	socket.send_string(to_send)


