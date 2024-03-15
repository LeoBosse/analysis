#!/usr/bin/python3
# -*-coding:utf-8 -*



#######################################################################
# Contains the PLIP data object.
# This is not useable as is, and should be develooped further once PLIP data is standard. Right now (2024), the format still moves regularly...


# Author: Léo Bosse
# License: the freeest one, do whatever with my bad code if it can be helpfull, nobody really cares!
#######################################################################




import numpy as np
import matplotlib.pyplot as plt
from utils import *
from scipy import signal
import sys as sys
import os
from subprocess import call

import h5py as h5

import datetime as dt
import scipy.io

# import observation as observation
# import observation as geometry
from vendange_configuration import *


class PLIP:
    def __init__(self, data_file_name):

        self.path = global_configuration.plip_data_path
        self.file_name = data_file_name

        self.timeshift = dt.timedelta(seconds=80)

        self.LoadData()
       

    def LoadData(self):
        with h5.File(self.path / self.file_name, 'r') as h5f:
            self.times = np.array([dt.datetime.fromtimestamp(t) - self.timeshift for t in h5f['times']]) 
            self.times += np.array([dt.timedelta(seconds = a) for a in range(len(self.times))])

            self.I = np.array(h5f['Intensity'])
            self.D = np.array(h5f['DoLP']) * 100 # / np.average(np.array(h5f['DoLP']) * 100)
            self.A = np.array(h5f['AoLP'])

            self.I_fine = np.array(h5f['Intensity_fine'])
            self.D_fine = np.array(h5f['DoLP_fine']) * 100
            self.A_fine = np.array(h5f['AoLP_fine'])

            self.I_laser = np.array(h5f['Intensity_Laser'])
            self.D_laser = np.array(h5f['DoLP_Laser']) * 100
            self.A_laser = np.array(h5f['AoLP_Laser'])