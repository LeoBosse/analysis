#!/usr/bin/python3
# -*-coding:utf-8 -*

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
import os
from subprocess import call
import datetime as dt

# import radis
import astropy.io.fits as fits
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_sun, get_body
from astropy.time import Time

from pybaselines import Baseline, utils

from spec_utils import *

import argparse


### Read the command line arguments. And store them in a dictionnary called args.
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--datafile', help='The data file to analyse. Relative path to /home/leob/iasb/analysis/data/skibotn_spectro/' )
args = parser.parse_args()


# file_name = 'data_100um_30sec_Sun Oct 8 2023_04.47.44_73.fits'
# file_name = 'data_30sec_100um_Mon Oct 9 2023_07.31.06_75.fits'

file_name = args.datafile  #'data_G1_30sec_100um_Tue Oct 10 2023_04.48.32_76.fits'

# file_name = 'calibration/G1_FVB_last_calib_Sat Oct 7 2023_17.11.32_52.fits'


dark_file_30 = 'darks/240x30sec_G1_FVB_Sat Oct 7 2023_15.56.51_37.fits'
dark_file_60 = 'darks/90x1min_G1_FVB_Sat Oct 7 2023_13.32.24_36.fits'

spec = Spectrum.FromFitsFile(file_name)


spec.ComputeDark(dark_file_30, dark_file_60)

spec.Subtract(spec.dark)

# serie.SetBaseLines()
# serie.data = serie.data - serie.baselines


### Plot Spectrum
spec.MakeFigure()


plt.legend()
plt.show()
