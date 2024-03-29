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
# file_name = args.datafile #'data_G1_30sec_100um_Tue Oct 10 2023_04.48.32_76.fits'



# file_names = [  'calibration/G1_FVB_last_calib_Sat Oct 7 2023_17.11.32_52.fits',
#                 'calibration/G1_FVB_last_calib_Sat Oct 7 2023_17.13.09_56.fits',
#                 'calibration/G1_FVB_last_calib_Sat Oct 7 2023_17.16.33_67.fits']

# file_names = [  'calibration/G1_FVB_last_calib_Sat Oct 7 2023_17.09.14_42.fits',
#                 'calibration/G1_FVB_last_calib_Sat Oct 7 2023_17.13.37_57.fits',
#                 'calibration/G1_FVB_last_calib_Sat Oct 7 2023_17.14.37_59.fits']

# file_names = ['/calibration/' + f for f in os.listdir(Serie.path + '/calibration/') if "last_calib" in f and ".fits" in f and int(f.split('.')[-2][3:]) < 53]

file_names = [  'calibration/G1_FVB_last_calib_Sat Oct 7 2023_17.07.51_38_0.1s_10um_553.64nm.fits',
                'calibration/G1_FVB_last_calib_Sat Oct 7 2023_17.08.30_40_0.01s_10um_553.64nm.fits',
                'calibration/G1_FVB_last_calib_Sat Oct 7 2023_17.09.14_42_0.001s_10um_553.64nm.fits',
                'calibration/G1_FVB_last_calib_Sat Oct 7 2023_17.10.21_46_0.0001s_10um_553.64nm.fits',
                'calibration/G1_FVB_last_calib_Sat Oct 7 2023_17.11.02_49_0.00001s_10um_553.64nm.fits']

dark_file_30 = Spectrum.path + 'darks/240x30sec_G1_FVB_Sat Oct 7 2023_15.56.51_37.fits'
dark_file_60 = Spectrum.path + 'darks/90x1min_G1_FVB_Sat Oct 7 2023_13.32.24_36.fits'


spectra = []
for n in file_names:
    spectra.append(Spectrum.FromFitsFile(n))


spectra = [s for s in spectra if np.max(s.data) < 65000]


# wavelength_shift = [0, -100, +100]
# color = ['g', 'b', 'r']
for i, s in enumerate(spectra):
    # s.color = color[i]
    # s.wavelengths += wavelength_shift[i]
    
    s.ComputeDark(dark_file_30, dark_file_60)
    s.Subtract(s.dark)

    # s.data /= np.average(s.data)


### Plot Spectrum
fig = plt.figure()
ax = fig.add_subplot(111)
for i, s in enumerate(spectra):
    # s.Plot(ax)
    plt.plot(s.wavelengths, s.data)


# data = np.array([s.data for s in spectra])
# plt.plot(spectra[0].wavelengths, np.average(data, axis=0), 'k')


fig = plt.figure()
for i in range(len(spectra) - 1):
    for j in range(i+1, len(spectra)):
        ratio = spectra[i].data / spectra[j].data
        ratio_theo = spectra[i].exposure / spectra[j].exposure
        
        plt.plot(spectra[i].wavelengths, ratio, label = f"{spectra[i].exposure} / {spectra[j].exposure}")


fig = plt.figure()
for i in range(len(spectra) - 1):
    for j in range(-10, 10):
        ratio = spectra[i].data[500+j] / spectra[i].data[800+j]
        
        plt.plot(spectra[i].exposure, ratio, '*')



plt.legend()
plt.show()
