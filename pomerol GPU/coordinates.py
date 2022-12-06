#!/usr/bin/python3
# -*-coding:utf-8 -*

import sys as sys
import numpy as np
import astropy as apy
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from astropy import units as u

import scipy.integrate as int
import matplotlib
import matplotlib.pyplot as plt

import datetime as dt


DtoR = np.pi / 180.
RtoD = 1. / DtoR


grenoble = EarthLocation(lon=5.936935, lat=45.212343, height=600*u.m)

# time = Time("2019-01-16 4:00") #UT
time = dt.datetime.strptime("2019-01-16 4:00", "%Y-%m-%d %H:%M") #UT

ptcu = SkyCoord(obstime="2019-01-16 4:00", location=grenoble, alt=45*u.deg, az=0*u.deg, frame='altaz')

ptcu.galactic
