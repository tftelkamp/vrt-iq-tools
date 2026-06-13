#!/usr/bin/env python3

# Copyright 2026 by Thomas Telkamp
#
# SPDX-License-Identifier: MIT

from astropy.coordinates import EarthLocation
import astropy.units as u

ANTENNA_TABLE = {
    'dwingeloo': EarthLocation(lat=52.81213723180477, lon=6.396346463227839, height=70.26, ellipsoid="WGS84"),
    'stockert': EarthLocation(lat=50.56944039751571, lon=6.721943350231514, height=434, ellipsoid="WGS84"),
    'stockert_2m': EarthLocation(lat=50.569511621782546, lon=6.723179139703107, height=434, ellipsoid="WGS84"),
    'stockert_1m': EarthLocation(lat=50.56933956889883, lon=6.723980412727782, height=434, ellipsoid="WGS84"),
    'bochum': EarthLocation(lat=51.426990, lon=7.192566, height=205, ellipsoid="WGS84"),
    'debilt': EarthLocation(lat=52.1018688270963, lon=5.17973411317635, height=51, ellipsoid="WGS84"),
    'debilt_gps': EarthLocation(lat=52.102008951201654, lon=5.1796902174021175, height=51, ellipsoid="WGS84"),
    'delft': EarthLocation(lat=51.999153, lon=4.373352, height=92, ellipsoid="WGS84"),  # estimation
    'ata_1a': EarthLocation.from_geocentric(-2524036.0307912203*u.m, -4123528.101172219*u.m, 4147706.408318585*u.m),
    'ata_4j': EarthLocation.from_geocentric(-2523898.1150373477*u.m, -4123456.314794732*u.m, 4147860.3045849088*u.m)
}
