#!/usr/bin/python3
# -*-coding:utf-8 -*

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

import observation as observation
import observation as geometry
from vendange_configuration import *


class PLIP:
    def __init__(self, bottle, data_file_name):
        self.bottle = bottle
        self.path = global_configuration.plip_data_path
        self.file_name = data_file_name

        self.LoadData()

    def LoadData(self):
        with h5.File(self.path + self.file_name, 'r') as h5f:
            self.times = np.array([dt.datetime.fromtimestamp(t) for t in h5f['times']])
            self.I = np.array(h5f['Intensity_Laser'])
            self.D = np.array(h5f['DoLP_Laser']) * 100
            self.A = np.array(h5f['AoLP_Laser'])

        
