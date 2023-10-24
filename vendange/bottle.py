#!/usr/bin/python3
# -*-coding:utf-8 -*

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal
import scipy as sci
from scipy.interpolate import CubicSpline
import sys as sys
import os
from copy import deepcopy
from subprocess import call
import datetime as dt
import h5py as h5

import observation
import geometry
from utils import *
from vendange_configuration import *

import chaosmagpy as chaos


# Bottle classes.


class Bottle:
    def __init__(self, folder_name, auto=True, line=1, from_txt=False):
        self.folder = folder_name
        self.rotations = []

        self.AoBapp = False
        self.AoBlos = False
        self.AoRD = False
        self.time_stamp = ""

        self.from_txt = from_txt

        self.line = line

        self.valid = True

        try:  # If you specify the exact input file name (with .in)
            self.input_parameters = self.ReadInputFile(pwd_data + self.folder)
        except:
            try:  # If you only specify the folder, will look for a default input.in file
                self.input_parameters = self.ReadInputFile(pwd_data + self.folder + "/input.in")
            except:  # Don't use it, or verify it works first.
                self.input_parameters = self.ReadInputFile(pwd_data + self.folder + ".in")

        self.data_file_name = pwd_data + \
            self.input_parameters["data_files"].split(",")[0]
        self.instrument_name = self.input_parameters["instrument_name"]
        if self.instrument_name == "ptcu":
            self.instrument_name = "carmen"
        elif self.instrument_name == "ptcu_v2":
            self.instrument_name = "corbel"
        elif self.instrument_name == "corbel2":
            self.instrument_name = "corbel"

        try:
            self.observation_type = self.input_parameters["observation_type"]
        except:
            self.observation_type = "fixed"
        print("INFO: Observation type is", self.observation_type)

        self.NoVref = True
        try:
            if self.input_parameters["I_ref"] == "off":
                self.NoVref = True
        except:
            self.NoVref = False

        try:
            self.saving_name = self.input_parameters["saving_name"]
            self.saving_name += "_voie" + str(self.line)
        except:
            print("WARNING: No saving name specified. Will be saved in "
                  + self.data_file_name)
            self.saving_name = ""

        self.date, self.location, self.filters, self.azimut, self.elevation, self.com = "", "", [
            "0"], False, False, ""

        try:
            self.azimut = float(self.input_parameters["azimut"]) * DtoR
            self.elevation = float(self.input_parameters["elevation"]) * DtoR
        except:
            print("WARNING: No observation angles info in the input file.")

        try:
            self.time_zone = int(self.input_parameters["time_zone"])
            self.config_time_format = self.input_parameters["config_time_format"]
        except:
            self.time_zone = 0
            self.config_time_format = "UT"
        self.time_zone = dt.timedelta(hours=self.time_zone)

        # See PTCUBottle or SPPBottle for specific instructions

    def SetInfoFromDataFileName(self, f, data_f):
        self.folders, self.data_folders = f.split("/"), data_f.split("/")
        # i = len(data_folders) - 1
        self.folders_index = -3

        # self.instrument_name = self.folders[self.folders_index]

        ### Follows in SPPBottle or PTCUBottle

    def MiseEnBouteille(self):
        # Loading the data from the files, then creating a list of Rotation objects depending on the instrument name specified. Can be 'spp', 'spp2014', 'fake_spp', 'ptcu', 'fake_ptcu'
        if not self.from_txt:
            self.LoadData()
            
            self.SetJumps()
            
            if self.valid:
                self.CleanRotations()
            

            # Now that we have all our rotations, creating lists of usefull data for easier smoothing
                self.CreateLists()
                self.GetSmoothLists()
                self.TestResample()
                print(self.data['DoLP'])
        else:
            self.LoadFromTxt()

        if self.valid:
            self.Geometry()

    def ReadInputFile(self, filename):
        """Read a given input file and returns a dictionnary. First word = key, second = value. Remove empty lines, as many arguments as you want, in any order that you want."""
        with open(filename, "r") as f:
            input = f.readlines()
        input = [l.split(maxsplit=2) for l in input if l != "\n"]
        dict = {}
        for i in input:
            dict[i[0]] = i[1]
        # print(dict)
        return dict

    def SetJump(self, tail=False):
        if not tail:
            self.jump_unit = self.input_parameters["jump_unit"].lower()

        try:
            self.jump_mode = self.input_parameters["jump_mode"]
        except:
            self.jump_mode = "lenght"


        if not tail:
            param = 'head_jump'
        else:
            param = 'tail_jump'


        jump = self.input_parameters[param]

        if ';' in jump:
            jump = jump.split(';')[self.line-1]


        if ':' in jump:
            jump = dt.datetime.strptime(jump, '%Y%m%d-%H:%M:%S')
            jump -= self.datetime
        else:
            jump = float(jump)
            if "sec" in self.jump_unit:
                jump = dt.timedelta(seconds=jump)
            elif "min" in self.jump_unit:
                jump = dt.timedelta(minutes=jump)
            elif "h" in self.jump_unit:
                jump = dt.timedelta(hours=jump)


        if self.jump_mode == "length" or (tail and jump.total_seconds() == 0.):
            jump = self.data['Times'].iloc[-1] - jump


        return jump


    def SetJumps(self):
        self.head_jump = self.SetJump()
        self.tail_jump = self.SetJump(tail = True)

        # self.jump_unit = self.input_parameters["jump_unit"].lower()
        #
        # # try:
        # #     self.head_jump = float(self.input_parameters["head_jump"].split(';')[self.line-1])
        # # except:
        # #     self.head_jump = float(self.input_parameters["head_jump"])
        # # try:
        # #     self.tail_jump = float(self.input_parameters["tail_jump"].split(';')[self.line-1])
        # # except:
        # #     self.tail_jump = float(self.input_parameters["tail_jump"])
        #
        # self.head_jump = self.input_parameters["head_jump"]
        # if ';' in self.input_parameters["head_jump"]:
        #     self.head_jump = self.head_jump.split(';')[self.line-1]
        # if ':' in self.head_jump:
        #     self.head_jump = dt.strptime(self.head_jump, '%y%m%d-%H:%M:%S')
        #     self.head_jump -= self.rotations[0]
        # else:
        #     self.head_jump = float(self.input_parameters["head_jump"])
        #
        # self.tail_jump = self.input_parameters["tail_jump"]
        # if ';' in self.input_parameters["tail_jump"]:
        #     self.tail_jump = self.tail_jump.split(';')[self.line-1]
        #
        # if ':' in self.tail_jump:
        #     self.tail_jump = dt.strptime(self.tail_jump, '%y%m%d-%H:%M:%S')
        #     self.tail_jump -= self.rotations[0]
        # else:
        #     self.tail_jump = float(self.input_parameters["tail_jump"])
        #
        #
        # try:
        #     self.jump_mode = self.input_parameters["jump_mode"]
        # except:
        #     self.jump_mode = "lenght"
        #
        # # Until the bug at the begining of the observation is not fixed, delete head and tail after treatment. When fixed, we'll be able to do it before
        # # if self.jump_unit in ("seconds", "minutes", "hours"):
        # # if "sec" in self.jump_unit or "min" in self.jump_unit or "h" in self.jump_unit:
        # #     self.head_jump = 0
        # #     self.tail_jump = 0
        #
        # # and self.jump_unit in ["seconds", "minutes", "hours"]:
        # if self.jump_mode == "time":
        #     if "sec" in self.jump_unit:
        #         self.head_jump = dt.timedelta(seconds=self.head_jump)
        #         self.tail_jump = dt.timedelta(seconds=self.tail_jump)
        #     elif "min" in self.jump_unit:
        #         self.head_jump = dt.timedelta(minutes=self.head_jump)
        #         self.tail_jump = dt.timedelta(minutes=self.tail_jump)
        #     elif "h" in self.jump_unit:
        #         # print("DEBUG jumps", [r.time.total_seconds() for r in self.rotations[:10]])
        #         # print(time.timedelta(hours=float(self.input_parameters["head_jump"])), time.timedelta(hours=float(self.input_parameters["tail_jump"])))
        #         self.head_jump = dt.timedelta(hours=self.head_jump)
        #         self.tail_jump = dt.timedelta(hours=self.tail_jump)
        #
        # if self.jump_mode == "length" or self.tail_jump.total_seconds() == 0.:
        #     try:
        #         self.tail_jump = self.data['Times'][-1] - self.tail_jump
        #     except:
        #         self.tail_jump = self.rotations[-1].time - self.tail_jump
        #
        print('DEBUG JUMPS', self.head_jump, self.tail_jump)

    def CorrectDensity(self):
        try:
            use_density = int(self.input_parameters["density"].split(";")[self.line-1])
        except KeyError:
            print("Warning, no density paramters in the input files (exist only since 20220523). Not a probelm if you did not use a density filter during your mesurments.")
            return False

        print(f"use_density {use_density}")
        if not use_density:
            return False

        self.time_unit = self.input_parameters["density_time_unit"].lower()


        self.density_factor = int(self.input_parameters["density_factor"].split(";")[self.line-1])
        time_period = self.input_parameters["density_period"].split("_")[self.line-1].split(":")
        start = float(time_period[0])
        end = float(time_period[1])

        if "sec" in self.time_unit:
            start = dt.timedelta(seconds=start)
            end = dt.timedelta(seconds=end)
        elif "min" in self.time_unit:
            start = dt.timedelta(minutes=start)
            end = dt.timedelta(minutes=end)
        elif "h" in self.time_unit:
            start = dt.timedelta(hours=start)
            end = dt.timedelta(hours=end)

        self.density_start = start
        self.density_end = end

        return True

    def LoadData(self):
        # Loading the data files

        # See SPPBottle/PTCUBottle for first part

        # self.nb_rot = len(self.rotations)
        
        self.nb_rot = len(self.data)
        self.SetConfig()
        self.SetTimes()

    def SetConfig(self):
        try:
            self.Imin = self.input_parameters["Imin"].split(";")
            self.Imin = float(self.Imin[self.line - 1])
        except:
            self.Imin = 0.1
        try:
            self.Imax = self.input_parameters["Imax"].split(";")
            self.Imax = float(self.Imax[self.line - 1])
        except:
            self.Imax = 10000
        try:
            self.DoLPmax = self.input_parameters["DoLPmax"].split(";")
            self.DoLPmax = float(self.DoLPmax[self.line - 1])
        except:
            self.DoLPmax = 100
        try:
            self.DoLPmin = self.input_parameters["DoLPmin"].split(";")
            self.DoLPmin = float(self.DoLPmin[self.line - 1])
        except:
            self.DoLPmin = 0

    def CreateLists(self):
        # At this point, all times are expressed in seconds since the first rotation before deleting the head_jump
        # If SPP before 2019 03 10, all times are expressed in seconds since midnight of the first observation day.

        # print("DEBUG TIME LISTS: ", self.time, self.head_jump, self.rotations[0].time)
        # self.data['Times'] = np.array([r.time - self.head_jump for r in self.rotations])
        
        norm = self.data['Times'].iloc[0]
        self.data['Times'] -= norm
        
        # self.all_times_since_start = np.array([r.time  for r in self.rotations])
        # self.all_times_since_start = np.array([r.time - self.rotations[0].time for r in self.rotations])

        # self.all_datetimes = [time.gmtime(r.time + self.time + self.head_jump) for r in self.rotations]
        # self.all_datetimes = [time.gmtime(time.mktime(self.DateTime()) + r.time + self.time + self.head_jump) for r in self.rotations]

        # self.all_datetimes = [self.DateTime() + self.head_jump + t for t in self.data['Times']]
        # self.all_datetimes = [self.DateTime() + t for t in self.all_times_since_start]

        # print(self.head_jump, self.data['Times'][0], self.data['Times'][-1], self.all_datetimes[0], self.all_datetimes[-1])

        # for i in range(1, self.nb_rot):
        #     if self.data['Times'][i] < self.data['Times'][i-1]:
        #         self.data['Times'][i] += 24 * 3600

        # self.AoLP_correction = float(self.input_parameters["AoLP_correction"])*DtoR

        if not self.from_txt and self.instrument_name in ["corbel", "gdcu", "ptcu_v2"]:
            # if self.DateTime() < dt.datetime(year=2020, month=10, day=1):
            # self.AoLP_correction = 0  # USE only for data downloaded from the database BEFORE October 2020!!!! For everything observed after that or downloaded via KVA20, you have to apply the first correction below

            self.AoLP_correction = -(self.config['IS_PolarizerOffset' + str(self.line)]) * DtoR #This is the norm now (for all data downloaded after oct 2020)

            # self.AoLP_correction = (self.config['IS_PolarizerOffset' + str(self.line)] + 45) * DtoR ## This was the formula apllied to convert all data from before oct 2020 in the database.

            ##################################################
            ##################################################
            # ATTENTION, A UTILISER UNIQUEMENT POURT DES CAS PARTICULIER. DANS LE DOUTE, COMMENTEZ OU VERIFIER LE PARAMETRE 'AoLP_correction' DANS LE FICHIER D'INPUT.
            # print("self.AoLP_correction", self.AoLP_correction*RtoD)
            # self.AoLP_correction -= float(self.input_parameters["AoLP_correction"]) * DtoR
            # self.AoLP_correction -= float(self.input_parameters["AoLP_correction"].split(";")[self.line - 1]) * DtoR
            # print("self.AoLP_correction", self.AoLP_correction*RtoD)
            ##################################################
            ##################################################

        elif not self.from_txt and self.instrument_name == "spp":
            self.AoLP_correction = float(self.input_parameters["AoLP_correction"]) * DtoR
        else:
            self.AoLP_correction = 0
        print(self.config['IS_PolarizerOffset1'])
        print("AoLP correction:", self.line, self.AoLP_correction * RtoD)

        # self.data['V'] = np.array([r.V for r in self.rotations])
        # self.data['Vcos'] = np.array([r.Vcos for r in self.rotations])
        # self.data['Vsin'] = np.array([r.Vsin for r in self.rotations])
        # try:
        #     self.all_Iref = np.array([r.Iref for r in self.rotations])
        # except:
        #     self.NoVref = True
        #     print("WARNING: No reference canal")

        # self.all_TempPM = np.array([r.TempPM for r in self.rotations])
        # self.all_TempOptical = np.array(
        #     [r.TempOptical for r in self.rotations])
        # self.all_TempAmbiant = np.array(
        #     [r.TempAmbiant for r in self.rotations])

        # self.data['I0'] = np.array([r.I0 for r in self.rotations])
        # self.data['DoLP'] = np.array([r.DoLP for r in self.rotations])

        # AoLP_correction -133 en janvier 2019 -> angle de calib pas sur avant
        self.data['AoLP'] += self.AoLP_correction
        self.data['AoLP'] = SetAngleBounds(self.data['AoLP'], -np.pi / 2, np.pi / 2)

    # def CreateListsFromGivenData(self):
    #     ###At this point, all times are expressed in seconds since the first rotation before deleting the head_jump
    #     ###If SPP before 2019 03 10, all times are expressed in seconds since midnight of the first observation day.
    #     print("DEBUG TIME LISTS: ", self.time, self.head_jump, self.rotations[0].time)
    #     # self.data['Times'] = np.array([r.time - self.head_jump for r in self.rotations])
    #     self.data['Times'] = np.array([r.time for r in self.rotations])
    #     # self.all_times_since_start = np.array([r.time  for r in self.rotations])
    #     self.all_times_since_start = np.array([r.time - self.rotations[0].time for r in self.rotations])
    #
    #     # self.all_datetimes = [time.gmtime(r.time + self.time + self.head_jump) for r in self.rotations]
    #     # self.all_datetimes = [time.gmtime(time.mktime(self.DateTime()) + r.time + self.time + self.head_jump) for r in self.rotations]
    #
    #     self.all_datetimes = [self.DateTime() + t for t in self.data['Times']]
    #     # self.all_datetimes = [self.DateTime() + t for t in self.all_times_since_start]
    #
    #     print(self.head_jump, self.data['Times'][0], self.data['Times'][-1], self.all_datetimes[0], self.all_datetimes[-1])
    #
    #     # for i in range(1, self.nb_rot):
    #     #     if self.data['Times'][i] < self.data['Times'][i-1]:
    #     #         self.data['Times'][i] += 24 * 3600
    #
    #     self.AoLP_correction = float(self.input_parameters["AoLP_correction"])*DtoR
    #     print(self.AoLP_correction, self.AoLP_correction*RtoD)
    #
    #     self.data['V']      = np.array([r.V    for r in self.rotations])
    #     self.data['Vcos'] = np.array([r.Vcos for r in self.rotations])
    #     self.data['Vsin'] = np.array([r.Vsin for r in self.rotations])
    #     try:
    #         self.all_Iref = np.array([r.Iref for r in self.rotations])
    #     except:
    #         self.NoVref = True
    #         print("WARNING: No reference canal")
    #
    #     self.all_TempPM = np.array([r.TempPM for r in self.rotations])
    #     self.all_TempOptical = np.array([r.TempOptical for r in self.rotations])
    #     self.all_TempAmbiant = np.array([r.TempAmbiant for r in self.rotations])
    #
    #     self.data['I0']      = np.array([r.I0   for r in self.rotations])
    #     self.data['DoLP'] = np.array([r.DoLP for r in self.rotations])
    #
    #     ### AoLP_correction -133 en janvier 2019 -> angle de calib pas sur avant
    #     self.data['AoLP'] = np.array([r.AoLP for r in self.rotations]) + self.AoLP_correction
    #     self.data['AoLP'] = SetAngleBounds(self.data['AoLP'], -np.pi/2, np.pi/2)

    def SetSmoothLists(self):
        print("Set Smooth Lists")
        # Smoothing procedure. Can smooth for a given number of rotations (smoothing_unit==rotations) or a  given time period (smoothing_unit==seconds).
        # Average the data over smoothing_factor rotations
        self.smoothing_factor = float(self.input_parameters["smoothing_factor"])
        self.smoothing_unit = self.input_parameters["smoothing_unit"].lower()
        self.smoothing_method = self.input_parameters["smoothing_method"].lower(
        )

        self.nb_smooth_rot = self.nb_rot

        self.data["time_deltas"] = self.data['Times'].diff() #[self.data['Times'].iloc[i].total_seconds() - self.data['Times'].iloc[i - 1].total_seconds() for i in range(1, len(self.data['Times']))] # in seconds
        self.avg_dt = self.data["time_deltas"].mean().total_seconds() * 1000 #1000 * np.average(self.data["time_deltas"]).total_seconds()  # in millisec
        # self.time_deltas.append(self.avg_dt/1000)

        print("AVG DT (millisec)", self.avg_dt)

    def GetSmoothLists(self):

        self.SetSmoothLists()

        print("Smooting data over {} {}".format(
            self.smoothing_factor, self.smoothing_unit))
        self.data['smooth_V'] = self.GetSliddingAverage(
            self.data['V'],    self.data['Times'], self.smoothing_factor, self.smoothing_unit)
        self.data['smooth_Vcos'] = self.GetSliddingAverage(
            self.data['Vcos'], self.data['Times'], self.smoothing_factor, self.smoothing_unit)
        self.data['smooth_Vsin'] = self.GetSliddingAverage(
            self.data['Vsin'], self.data['Times'], self.smoothing_factor, self.smoothing_unit)
        if not self.NoVref:
            self.data['smooth_Iref'] = self.GetSliddingAverage(
                self.data['Iref'], self.data['Times'], self.smoothing_factor, self.smoothing_unit)
            self.Iref_average = self.data['Iref'].mean()

        self.nb_smooth_rot = len(self.data['smooth_V'])

        # Calculate the smooth I0, DoLP and AoLP
        self.data['smooth_I0'], self.data['smooth_DoLP'], self.data['smooth_AoLP'] = np.zeros(
            self.nb_smooth_rot), np.zeros(self.nb_smooth_rot), np.zeros(self.nb_smooth_rot)
        
        # for i in range(self.nb_smooth_rot):
        #     self.data['smooth_I0'].iloc[i], self.data['smooth_DoLP'].iloc[i], self.data['smooth_AoLP'].iloc[i] = Rotation.GetLightParameters(
        #         self.data['smooth_V'].iloc[i], self.data['smooth_Vcos'].iloc[i], self.data['smooth_Vsin'].iloc[i])
        self.data['smooth_I0'], self.data['smooth_DoLP'], self.data['smooth_AoLP'] = self.GetLightParameters(
                self.data['smooth_V'], self.data['smooth_Vcos'], self.data['smooth_Vsin'])

        self.data['smooth_AoLP'] = self.data['smooth_AoLP'] + self.AoLP_correction

        self.data['smooth_AoLP'] = SetAngleBounds(
            self.data['smooth_AoLP'], -np.pi / 2, np.pi / 2)

        self.V_average = self.data['V'].mean() #np.average(self.data['V'])
        self.Vcos_average = self.data['Vcos'].mean()
        self.Vsin_average = self.data['Vsin'].mean()
        self.I0_average, self.DoLP_average, self.AoLP_average = self.GetLightParameters(
            self.V_average, self.Vcos_average, self.Vsin_average)

        self.AoLP_average += self.AoLP_correction
        self.AoLP_average = SetAngleBounds(
            self.AoLP_average, -np.pi / 2, np.pi / 2)

        print("Computing Error bars...")
        # Get Variance
        if self.smoothing_unit == "seconds":
            smoothing_factor = self.smoothing_factor * 1000
        elif self.smoothing_unit == "minutes":
            smoothing_factor = self.smoothing_factor * 60 * 1000
        elif self.smoothing_unit == "hours":
            smoothing_factor = self.smoothing_factor * 3600 * 1000

        # self.std_I0 = 0         #np.sqrt(4 * self.data['V']     / 500)
        # self.std_smooth_I0 = 0# np.sqrt(4 * self.smooth_V  / smoothing_factor)
        #
        # self.std_DoLP =     0#    np.sqrt(4 * (1 + ((self.data['DoLP']/100)      ** 2 / 2)) / (self.data['I0']       * 500)) * 100
        # self.std_smooth_DoLP =0#     np.sqrt(4 * (1 + ((self.data['smooth_DoLP']/100) ** 2 / 2)) / (self.data['smooth_I0'] * smoothing_factor)) * 100
        #
        # self.std_AoLP =     0#    np.sqrt(1 / ((self.data['DoLP']/100)     ** 2 * self.data['I0']       * 500))
        # self.std_smooth_AoLP = 0#    np.sqrt(1 / ((self.data['smooth_DoLP']/100) ** 2 * self.data['smooth_I0'] * smoothing_factor))
        # self.smooth_AoLP_upper = 0#self.data['smooth_AoLP'] + self.std_smooth_AoLP
        # self.smooth_AoLP_lower = 0#self.data['smooth_AoLP'] - self.std_smooth_AoLP
        self.data['std_I0'] = np.sqrt(2 * self.data['I0'] / self.avg_dt)
        self.data['std_smooth_I0'] = np.sqrt(2 * self.data['smooth_I0'] / smoothing_factor)

        # print(self.data['V'])
        # print(self.std_I0)
        # print(self.std_smooth_I0)

        self.data['std_DoLP'] = np.sqrt(4 * (1 + ((self.data['DoLP'] / 100) ** 2 / 2)) / (self.data['I0'] * self.avg_dt)) * 100
        self.data['std_smooth_DoLP'] = np.sqrt(4 * (1 + ((self.data['smooth_DoLP'] / 100) ** 2 / 2)) / (self.data['smooth_I0'] * smoothing_factor)) * 100

        self.data['std_AoLP'] = np.sqrt(1 / ((self.data['DoLP'] / 100) ** 2 * self.data['I0'] * self.avg_dt))
        self.data['std_smooth_AoLP'] = np.sqrt(1 / ((self.data['smooth_DoLP'] / 100) ** 2 * self.data['smooth_I0'] * smoothing_factor))

        # for i in range(len(self.std_AoLP)):
        #     self.std_AoLP[i] = min(np.pi, self.std_AoLP[i])
        #     self.std_smooth_AoLP[i] = min(np.pi, self.std_smooth_AoLP[i])
        self.data['std_AoLP'] = np.clip(self.data['std_AoLP'], None, np.pi)
        self.data['std_smooth_AoLP'] = np.clip(self.data['std_smooth_AoLP'], None, np.pi)

        self.data['smooth_AoLP_upper'] = self.data['smooth_AoLP'] + self.data['std_smooth_AoLP']
        self.data['smooth_AoLP_lower'] = self.data['smooth_AoLP'] - self.data['std_smooth_AoLP']
        # print(self.data['smooth_DoLP'].shape, self.sm ooth_AoLP.shape, self.var_DoLP.shape, self.var_AoLP.shape)

        def SN(I, D, T): return D * np.sqrt(I * T) / 2

        self.data['all_SN'] = SN(self.data['I0'], self.data['DoLP'] / 100, self.avg_dt)
        self.data['smooth_SN'] = SN(self.data['smooth_I0'], self.data['smooth_DoLP'] / 100, smoothing_factor)


        self.data['I0_diff'] = np.gradient(self.data['smooth_I0'], 1.) #self.avg_dt)
        self.data['I0_diff_std'] = self.GetDiffErrors(smooth=True)
        self.data['all_I0_diff'] = np.gradient(self.data['I0'], self.avg_dt)
        self.data['all_I0_diff_std'] = self.GetDiffErrors(smooth=False)

        # self.smooth_I0_mean = np.mean(self.data['smooth_I0'])
        # self.smooth_I0_var = np.var(self.data['smooth_I0'])
        #
        # self.data['smooth_I0'] = (self.data['smooth_I0'] - self.smooth_I0_mean) / self.smooth_I0_var

        self.GetSliddingCorrCoef()

        # self.SetUnifyAngles()
        self.graph_angle_shift = 0

        print("Get Smooth Lists: DONE")

    def GetDiffErrors(self, smooth=True):
        if smooth:
            I = self.data['smooth_I0']
            I_err = self.data['std_smooth_I0']
            d = self.data['I0_diff']
        else:
            I = self.data['I0']
            I_err = self.data['std_I0']
            d = self.data['all_I0_diff']

        # print(I_err)

        d_after = np.append(I_err.iloc[1:], I_err.iloc[-1])
        d_before = np.append(I_err.iloc[0], I_err.iloc[:-1])
        I_after = np.append(I.iloc[1:], I.iloc[-1])
        I_before = np.append(I.iloc[0], I.iloc[:-1])

        diff_errors = abs(d) * .5 * np.sqrt(d_after / I_after**2 + d_before / I_before**2)

        # print(diff_errors)

        return diff_errors


    def SetUnifyAngles(self):
        print("Unifying AoLPs for least deviation")

        self.smooth_AoLP_upper = self.data['smooth_AoLP'] + self.std_smooth_AoLP
        self.smooth_AoLP_lower = self.data['smooth_AoLP'] - self.std_smooth_AoLP

        self.data['smooth_AoLP'], smooth_graph_angle_shift = UnifyAngles(
            self.data['smooth_AoLP'])
        self.data['AoLP'], all_graph_angle_shift = UnifyAngles(self.data['AoLP'])
        self.smooth_AoLP_upper, tmp = UnifyAngles(self.smooth_AoLP_upper)
        self.smooth_AoLP_lower, tmp = UnifyAngles(self.smooth_AoLP_lower)
        self.graph_angle_shift = max(
            all_graph_angle_shift, smooth_graph_angle_shift)

        if all_graph_angle_shift != self.graph_angle_shift:
            self.data['AoLP'], all_graph_angle_shift = UnifyAngles(
                self.data['AoLP'], manual_shift=self.graph_angle_shift)
            self.smooth_AoLP_upper, tmp = UnifyAngles(
                self.smooth_AoLP_upper, manual_shift=self.graph_angle_shift)
            self.smooth_AoLP_lower, tmp = UnifyAngles(
                self.smooth_AoLP_lower, manual_shift=self.graph_angle_shift)
        if smooth_graph_angle_shift != self.graph_angle_shift:
            self.data['smooth_AoLP'], smooth_graph_angle_shift = UnifyAngles(
                self.data['smooth_AoLP'], manual_shift=self.graph_angle_shift)
            self.smooth_AoLP_upper, tmp = UnifyAngles(
                self.smooth_AoLP_upper, manual_shift=self.graph_angle_shift)
            self.smooth_AoLP_lower, tmp = UnifyAngles(
                self.smooth_AoLP_lower, manual_shift=self.graph_angle_shift)

        if self.graph_angle_shift == 1:
            self.AoLP_average = SetAngleBounds(self.AoLP_average, 0, np.pi)
            self.smooth_AoLP_upper = SetAngleBounds(
                self.smooth_AoLP_upper, 0, np.pi)
            self.smooth_AoLP_lower = SetAngleBounds(
                self.smooth_AoLP_lower, 0, np.pi)
            # if self.AoLP_average < 0:
            #     self.AoLP_average += np.pi
        elif self.graph_angle_shift == 0:
            self.AoLP_average = SetAngleBounds(
                self.AoLP_average, -np.pi / 2, np.pi / 2)
            self.smooth_AoLP_upper = SetAngleBounds(
                self.smooth_AoLP_upper, -np.pi / 2, np.pi / 2)
            self.smooth_AoLP_lower = SetAngleBounds(
                self.smooth_AoLP_lower, -np.pi / 2, np.pi / 2)
            # if self.AoLP_average > np.pi/2:
            #     self.AoLP_average -= np.pi

        # for i, up in enumerate(self.smooth_AoLP_upper):
        #     if up < self.data['smooth_AoLP'][i]:


    def GetSliddingCorrCoef(self, window_size = 10, smooth = False):

        # Get the window size in number of points from the length in seconds and the average delta t.
        window_size = int(window_size * 1000 / self.avg_dt)
        if window_size % 2 == 1:
            window_size += 1

        # normalized numpy array of the window.
        window = np.ones(window_size) / window_size

        X_arr = self.data['all_I0_diff']
        Y_arr = self.data['smooth_DoLP']

        if not smooth:
            X_arr = self.data['all_I0_diff']
            Y_arr = self.data['DoLP']

        data = pd.DataFrame({'X':X_arr,
                             'Y':Y_arr})

        corr = data['X'].rolling(window_size, center=True).corr(data['Y'],method='pearson')

        # print('DEBUG ROLLING CORR', window_size)
        # print(corr[window_size:window_size+100])
        #
        # corr = np.zeros_like(self.data['Times'])
        # for it in range(len(self.data['Times'])):
        #     start = max(0, it - int(window_size / 2))
        #     end = min(len(self.data['Times']), it + int(window_size / 2))
        #
        #     corr[it] = np.corrcoef([X_arr[start:end], Y_arr[start:end]])[0, 1]
        #
        #     # print(corr[it] == data['X'][start:end].corr(data['Y'][start:end]))
        #
        # print(corr[window_size:window_size+100])
        return corr


    def GetFourierTransform(self, var, threshold_freq = 0.016):
        if var.lower() == "i0":
            signal = self.data['I0']
        elif var.lower() == "dolp":
            signal = self.data['DoLP']
        elif var.lower() == "aolp":
            signal = self.data['AoLP']

        ### Code copied from https://www.earthinversion.com/techniques/signal-denoising-using-fast-fourier-transform/
        delta = self.avg_dt/1000 #in sec
        n = len(signal)
        fhat = np.fft.fft(signal, n) #computes the fft
        # psd = fhat * np.conj(fhat) / n
        freq = (1/(delta*n)) * np.arange(n) #frequency array
        # idxs_half = np.arange(1, np.floor(n/2), dtype=np.int32) #first half index
        # psd_real = np.abs(psd[idxs_half]) #amplitude for first half


        ## Filter out noise
        # sort_psd = np.sort(psd_real)[::-1]
        # print(len(sort_psd))
        # threshold = sort_psd[threshold_freq]
        # psd_idxs = psd > threshold #array of 0 and 1
        # psd_clean = psd * psd_idxs #zero out all the unnecessary powers

        if threshold_freq > 0:
            psd_idxs = freq < threshold_freq #array of 0 and 1
        else:
            psd_idxs = np.ones_like(fhat)

        print(threshold_freq)
        print(freq)

        # f = plt.figure()
        # plt.plot(freq, psd)
        # plt.plot(freq, psd_clean)
        # plt.show()

        # print(threshold)
        # print(psd)
        # print(psd_clean)

        fhat_clean = psd_idxs * fhat #used to retrieve the signal

        signal_filtered = np.fft.ifft(fhat_clean) #inverse fourier transform

        return signal_filtered


    def GetGeometryAngles(self, obs, B_model=None):
        # print("DEBUG obs", obs)

        obs.SinglePointGeometry(B_model=B_model)
        AoBapp = obs.eta_chaos
        AoBapp_ortho = AoBapp + np.pi / 2
        AoBlos = obs.Blos

        AoRD = obs.GetRayleighAngle(obs.RD_src_azimut, obs.RD_src_elevation)
        AoRD_ortho = AoRD + np.pi / 2

        if self.graph_angle_shift == 1:
            if AoBapp < 0:
                AoBapp += np.pi
            if AoRD < 0:
                AoRD += np.pi
        if self.graph_angle_shift == 0 and AoBlos > np.pi / 2:
            if AoRD_ortho > np.pi / 2:
                AoRD_ortho -= np.pi
            if AoBlos > np.pi / 2:
                AoBlos -= np.pi
            if AoBapp_ortho > np.pi / 2:
                AoBapp_ortho -= np.pi

        return obs, AoBapp, AoBlos, AoRD, AoBapp_ortho, AoRD_ortho

    def Geometry(self, to_initiate=True):

        print("Geometry: Start")

        wd = os.getcwd()
        # os.chdir(pwd_src + "Geometry/Leo/src/")

        A_lon, A_lat = False, False
        A_lon, A_lat = geometry.GetLonLatFromName(self.location.lower())

        h = self.GetAltitude()

        try:
            self.source_azimut = float(
                self.input_parameters["pollution_source_azimut"]) * DtoR
            try:
                self.source_elevation = float(
                    self.input_parameters["pollution_source_elevation"]) * DtoR
            except:
                self.source_elevation = 0
        except:
            self.source_azimut = self.source_elevation = 0

        B_model = chaos.load_CHAOS_matfile(
            global_configuration.chaos_model_file)

        if self.observation_type == "fixed" and self.azimut is not None and self.elevation is not None:
            # print("DEBUG SOURCE:", self.source_azimut*RtoD, self.source_elevation*RtoD)
            try:
                # print("DEBUG: AZ/EL", self.azimut*RtoD, self.elevation*RtoD)
                self.observation, self.AoBapp, self.AoBlos, self.AoRD, self.AoBapp_ortho, self.AoRD_ortho = self.GetGeometryAngles(observation.ObservationPoint(
                    A_lon, A_lat, h, self.azimut, self.elevation, self.source_azimut, self.source_elevation, init_full=False), B_model=B_model)
            except:
                print("WARNING: No geometric angles where retrieved.")

        elif self.observation_type == "fixed_elevation_discrete_rotation":
            if to_initiate:
                self.discrete_rotation_elevation = float(self.input_parameters["discrete_rotation_elevation"]) * DtoR
                self.discrete_rotation_times = np.array([float(t) for t in self.input_parameters["discrete_rotation_times"].split("_")])
                self.discrete_rotation_azimuts = np.array([float(a) * DtoR for a in self.input_parameters["discrete_rotation_azimuts"].split("_")])

            ang_list = []
            for a in self.discrete_rotation_azimuts:
                ang_list.append(self.GetGeometryAngles(observation.ObservationPoint(A_lon, A_lat, h, a,
                                self.discrete_rotation_elevation, self.source_azimut, self.source_elevation, init_full=False), B_model=B_model))

            self.observation, self.AoBapp, self.AoBlos, self.AoRD, self.AoBapp_ortho, self.AoRD_ortho = zip(
                *ang_list)
            self.observation, self.AoBapp, self.AoBlos, self.AoRD, self.AoBapp_ortho, self.AoRD_ortho = np.array(self.observation), np.array(
                self.AoBapp), np.array(self.AoBlos), np.array(self.AoRD), np.array(self.AoBapp_ortho), np.array(self.AoRD_ortho)

            self.x_axis = self.data['Times']

        elif self.observation_type == "fixed_azimut_discrete_rotation":
            if to_initiate:
                self.discrete_rotation_azimut = float(
                    self.input_parameters["discrete_rotation_azimut"]) * DtoR
                self.discrete_rotation_times = [
                    float(t) for t in self.input_parameters["discrete_rotation_times"].split("_")]
                self.discrete_rotation_times_unit = self.input_parameters[
                    "discrete_rotation_times_unit"]
                self.discrete_rotation_elevations = [float(
                    a) * DtoR for a in self.input_parameters["discrete_rotation_elevations"].split("_")]

            ang_list = []
            for e in self.discrete_rotation_elevations:
                ang_list.append(self.GetGeometryAngles(observation.ObservationPoint(
                    A_lon, A_lat, h, self.discrete_rotation_azimut, e, self.source_azimut, self.source_elevation, init_full=False), B_model=B_model))

            self.observation, self.AoBapp, self.AoBlos, self.AoRD, self.AoBapp_ortho, self.AoRD_ortho = zip(
                *ang_list)
            self.observation, self.AoBapp, self.AoBlos, self.AoRD, self.AoBapp_ortho, self.AoRD_ortho = np.array(self.observation), np.array(
                self.AoBapp), np.array(self.AoBlos), np.array(self.AoRD), np.array(self.AoBapp_ortho), np.array(self.AoRD_ortho)

        elif self.observation_type == "fixed_elevation_continue_rotation":
            if to_initiate:
                self.continue_rotation_elevation = float(
                    self.input_parameters["continue_rotation_elevation"]) * DtoR
                self.continue_rotation_times = [
                    float(t) for t in self.input_parameters["continue_rotation_times"].split("_")]
                if self.continue_rotation_times[0] != 0:
                    self.continue_rotation_times.insert(0, 0)
                self.continue_rotation_times.append(
                    self.data['Times'].iloc[-1].total_seconds() / 60.)
                self.continue_rotation_times = np.array(
                    self.continue_rotation_times)
                self.nb_continue_rotation = len(
                    self.continue_rotation_times) - 1
                # self.rotation_direction = self.input_parameters["rotation_direction"]
            print("self.continue_rotation_times", self.continue_rotation_times)

            geo = geometry.Geometry("dummy", str(self.location), str(h), "e", str(
                self.continue_rotation_elevation * RtoD), str(self.source_azimut * RtoD), str(self.source_elevation * RtoD))
            try:
                self.azimut, self.observation = geo.FixedElevation(
                    direction=self.input_parameters["rotation_direction"], B_model=B_model)
            except:
                self.azimut, self.observation = geo.FixedElevation(
                    B_model=B_model)

            ang_list = [self.GetGeometryAngles(o, B_model=B_model)
                        for o in self.observation]
            self.observation, self.AoBapp, self.AoBlos, self.AoRD, self.AoBapp_ortho, self.AoRD_ortho = zip(
                *ang_list)
            self.observation, self.AoBapp, self.AoBlos, self.AoRD, self.AoBapp_ortho, self.AoRD_ortho = np.array(self.observation), np.array(
                self.AoBapp), np.array(self.AoBlos), np.array(self.AoRD), np.array(self.AoBapp_ortho), np.array(self.AoRD_ortho)

        os.chdir(wd)
        # print("Geometry:", self.AoBapp, self.AoBlos, self.AoRD)
        print("Geometry: DONE")

    def GetOneTimePeriod(self, start_time, end_time, smooth_data=True, direct=False):
        """Return the data points between the start and end times given."""
        start_time = dt.timedelta(minutes=start_time)
        end_time = dt.timedelta(minutes=end_time)

        if direct == True:
            I0 = self.data['I0'][(start_time < self.data['Times'])
                             & (self.data['Times'] < end_time)]
            DoLP = self.data['DoLP'][(start_time < self.data['Times'])
                                 & (self.data['Times'] < end_time)]
            AoLP = self.data['AoLP'][(start_time < self.data['Times'])
                                 & (self.data['Times'] < end_time)]
            nb_times = len(I0)
            return I0, DoLP, AoLP, nb_times
        else:
            V = self.data['V'].iloc[(start_time < self.data['Times'])
                           & (self.data['Times'] < end_time)]
            Vcos = self.data['Vcos'].iloc[(start_time < self.data['Times'])
                                 & (self.data['Times'] < end_time)]
            Vsin = self.data['Vsin'].iloc[(start_time < self.data['Times'])
                                 & (self.data['Times'] < end_time)]
            nb_times = len(V)
            return V, Vcos, Vsin, nb_times

    def GetAverageContinueRotations(self):
        """Returns a new bottle where all V, Vcos, Vsin are the average of the different (azimutal) rotations of this one. Is used for the mixer if the observation_type == fixed_elevation_continue_rotation."""

        # Get the smallest rotation (least number of points)
        min_nb_times = 10**100
        i = 0
        for i in range(0, self.nb_continue_rotation):
            start, end = self.continue_rotation_times[i], self.continue_rotation_times[i + 1]
            tmp_min_nb_times = self.GetOneTimePeriod(
                start, end, smooth_data=False)[3]
            if tmp_min_nb_times < min_nb_times:
                min_nb_times = tmp_min_nb_times
                i_min = i

        avg_V, avg_Vcos, avg_Vsin, nb_times = self.GetOneTimePeriod(
            self.continue_rotation_times[i_min], self.continue_rotation_times[i_min + 1], smooth_data=False)

        for i in range(1, self.nb_continue_rotation):
            start, end = self.continue_rotation_times[i], self.continue_rotation_times[i + 1]

            tmp_V, tmp_Vcos, tmp_Vsin, tmp_nb_times = self.GetOneTimePeriod(
                start, end, smooth_data=False)
            if tmp_nb_times > nb_times:
                tmp_V = ReduceList(tmp_V, nb_times)
                tmp_Vcos = ReduceList(tmp_Vcos, nb_times)
                tmp_Vsin = ReduceList(tmp_Vsin, nb_times)
            elif tmp_nb_times < nb_times:
                avg_V = ReduceList(avg_V, tmp_nb_times)
                avg_Vcos = ReduceList(avg_Vcos, tmp_nb_times)
                avg_Vsin = ReduceList(avg_Vsin, tmp_nb_times)
            avg_V += tmp_V
            avg_Vcos += avg_Vcos
            avg_Vsin += avg_Vsin
        avg_V /= self.nb_continue_rotation
        avg_Vcos /= self.nb_continue_rotation
        avg_Vsin /= self.nb_continue_rotation

        I0, DoLP, AoLP = Rotation.GetLightParameters(avg_V, avg_Vcos, avg_Vsin)

        new_bottle = deepcopy(self)

        new_bottle.continue_rotation_times = np.array(
            [self.continue_rotation_times[0], self.continue_rotation_times[1]])
        new_bottle.nb_continue_rotation = 1
        new_bottle.data['V'] = avg_V
        new_bottle.data['Vcos'] = avg_Vcos
        new_bottle.data['Vsin'] = avg_Vsin
        new_bottle.data['I0'] = I0
        new_bottle.data['DoLP'] = DoLP
        new_bottle.data['AoLP'] = AoLP
        new_bottle.all_times = np.arange(0)
        new_bottle.saving_name += "_rotations_average"

        for t in np.linspace(self.continue_rotation_times[0], self.continue_rotation_times[1], len(new_bottle.data['V'])):
            new_bottle.all_times = np.append(
                new_bottle.all_times, dt.timedelta(minutes=t))

        new_bottle.NoVref = True
        new_bottle.GetSmoothLists()
        new_bottle.Geometry(to_initiate=False)

        return new_bottle

    def GetWaveletTransform(self, values, smoothing_factor, unit = "seconds"):
        nb_values = len(values)
        if unit != "rotations":
            smoothing_factor *= 1000
            if unit == "minutes":
                smoothing_factor *= 60 #dt.timedelta(minutes = smoothing_factor)
            elif unit == "hours":
                smoothing_factor *= 3600 #dt.timedelta(hours = smoothing_factor)

            smoothing_factor = int(smoothing_factor / self.avg_dt)

        a = 0.75 # For a wavelet of around 12 seconds
        t0 = 0
        wave_fct = lambda t: -(t - t0) * np.exp(-(t-t0)**2 / (2*a)**2 / np.sqrt(2*np.pi) / a**3)

        window = np.linspace(-smoothing_factor/2, smoothing_factor/2, smoothing_factor)
        window = wave_fct(window)

        # plt.plot(np.linspace(-smoothing_factor/2, smoothing_factor/2, smoothing_factor), window)
        # plt.show()

        smooth_values = np.convolve(values, window, 'same')
        ### Correct the edge cases where the window does not overlapp completely.
        for i in range(int(smoothing_factor/2)):
            smooth_values[i]      *=  smoothing_factor / (smoothing_factor/2 + i)
            smooth_values[-1 - i] *=  smoothing_factor / (smoothing_factor/2 + i + 1) #Not sure why it needs a +1, but it works.

        return smooth_values


    def GetSliddingAverage(self, values, times, smoothing_factor, unit):
        """Given a list of values, a list of their associated times, the smoothing factor and its unit, return an array of smooth values. If unit==rotations each data is averaged with a window of smoothing_factor rotations centered on itself.
        If unit==seconds, minutes or hours, the smoothing factor is first converted to a number of rotations using the average delta t between all rotations points. Then each data is averaged with a window centered on itself as above. Both ways have the same results if each rotation has the same time step, but if the instrument has a probelm and each rotations are not equally spaced, this might cause unwanted effects.
        """
        nb_values = len(values)
        if unit != "rotations":
            smoothing_factor *= 1000
            if unit == "minutes":
                smoothing_factor *= 60 #dt.timedelta(minutes = smoothing_factor)
            elif unit == "hours":
                smoothing_factor *= 3600 #dt.timedelta(hours = smoothing_factor)

            smoothing_factor = int(smoothing_factor / self.avg_dt)
        if smoothing_factor % 2 == 1:
            smoothing_factor += 1

        window = np.ones(smoothing_factor) / smoothing_factor
        # smooth_values = signal.fftconvolve(values, window, 'same')
        smooth_values = np.convolve(values, window, 'same')
        ### Correct the edge cases where the window does not overlap completely.
        for i in range(int(smoothing_factor / 2)):
            smooth_values[i]      *=  smoothing_factor / (np.ceil(smoothing_factor / 2) + i)
            smooth_values[-1 - i] *=  smoothing_factor / (np.ceil(smoothing_factor / 2) + i + 1) #Not sure why it needs a +1, but it works.
        return smooth_values


    def PrintInfo(self):
        print("*********************************************")
        print("INFO:", self.data_file_name.split("/")[-3:])
        print("DoLP min:", np.min(self.data['smooth_DoLP']))
        print("DoLP max:", np.max(self.data['smooth_DoLP']))
        print("DoLP stddev:", np.std(self.data['smooth_DoLP']))
        print("AoLP stddev:", np.std(self.data['smooth_AoLP']) * RtoD)
        print("Averages: I0, DoLP, AoLP:", self.I0_average,
              self.DoLP_average, self.AoLP_average * RtoD)

        if self.observation_type == "fixed":
            print("AoBlos:", self.AoBlos * RtoD)
        print("*********************************************")

    def SaveTXT(self):
        print("Saving as .txt in", self.data_file_name
              + "/" + self.saving_name + '_results.txt')

    def ConvertRotTimeTogmTime(self, rot_time):
        t = dt.mktime(self.DateTime())
        tl = [dt.gmtime(t + rt) for rt in rot_time]
        return tl

    def GetRotationFrequency(self):
        all_freq = [1. / (self.data['Times'][i].total_seconds()
                          - self.data['Times'][i - 1].total_seconds()) for i in range(1, self.nb_rot)]
        return np.average(all_freq)

    def RenormTimes(self, shift):
        """Shift all times by a certain amount"""
        self.data['Times'] = np.array([t + shift for t in self.data['Times']])

    def DateTime(self, moment="start", delta=dt.timedelta(seconds=0), format="LT"):
        norm = dt.timedelta(seconds=0)
        if format == "LT" and self.config_time_format == "UT":
            norm = self.time_zone
        elif format == "UT" and self.config_time_format == "LT":
            norm = - self.time_zone

        if moment == "start":  # Datetime of first good rotation
            return self.data['Times'].iloc[0] + self.datetime + self.head_jump + delta + norm
            # return self.rotations[0].time + self.datetime + self.head_jump + delta
        elif moment == "end":  # Datetime of last good rotation
            return self.data['Times'].iloc[-1] + self.datetime + self.head_jump + delta + norm
            # return self.rotations[-1].time + self.datetime + self.head_jump + delta
        # Datetime written in the config file (== first bad rotation, before head_jump)
        elif moment == "config":
            return self.datetime + delta + norm
        # Datetime written in the config file (== first bad rotation, before head_jump)
        elif moment == "midnight":
            return dt.datetime(year=self.datetime.year, month=self.datetime.month, day=self.datetime.day) + delta + norm

    def GetInterpolation(self, time):
        """Return the linear interpolation of smooth data for the given times in seconds. Time must be in seconds since the start of this bottle."""
        all_times = list(map(lambda x: x.total_seconds(), self.data['Times']))

        # , left=0, right=0, period=None)
        interp_I0 = np.interp(time, all_times, self.data['smooth_I0'])
        # , left=0, right=0, period=None)
        interp_DoLP = np.interp(time, all_times, self.data['smooth_DoLP'])
        # , period=np.pi)
        interp_AoLP = np.interp(time, all_times, self.data['smooth_AoLP'])

        return interp_I0, interp_DoLP, interp_AoLP

    def __add__(self, bottle_to_add):

        print(f"Adding: {self} + {bottle_to_add}")

        if len(bottle_to_add.rotations) > 0:
            bottle_to_add.CleanRotations()

        norm = bottle_to_add.DateTime("start") - self.DateTime("start")
        for ir, r in enumerate(bottle_to_add.rotations):
            bottle_to_add.rotations[ir].time += norm
            # bottle_to_add.rotations[ir].time += self.rotations[-1].time

        if len(bottle_to_add.rotations) > 0:  # If not read from a file
            # print("DEBUG adding  with rotations")
            # print(len(self.rotations), len(bottle_to_add.rotations))
            self.rotations = np.append(self.rotations, bottle_to_add.rotations)

            # print(len(self.rotations), len(bottle_to_add.rotations))

            self.tail_jump = self.rotations[-1].time
            if not self.from_txt:
                # self.CleanRotations()

            # Now that we have all our rotations, creating lists of usefull data for easier smoothing
                self.CreateLists()
                self.GetSmoothLists()

        else:  # If read from a file (-f)
            # print("DEBUG adding  without rotations (from file)")
            bottle_to_add.all_times += bottle_to_add.DateTime("start") - self.DateTime("start")
            # bottle_to_add.all_times += self.data['Times'][-1]

            self.data['Times'] = np.append(self.data['Times'], bottle_to_add.all_times)
            self.nb_rot = len(self.data['Times'])

            self.data['smooth_V'] = np.append(self.data['smooth_V'], bottle_to_add.data['smooth_V'])
            self.data['smooth_Vcos'] = np.append(
                self.data['smooth_Vcos'], bottle_to_add.data['smooth_Vcos'])
            self.data['smooth_Vsin'] = np.append(
                self.data['smooth_Vsin'], bottle_to_add.data['smooth_Vsin'])

            self.nb_smooth_rot = len(self.data['smooth_V'])

            # Calculate the smooth I0, DoLP and AoLP
            self.data['smooth_I0'] = np.append(self.data['smooth_I0'], bottle_to_add.data['smooth_I0'])
            self.data['smooth_DoLP'] = np.append(
                self.data['smooth_DoLP'], bottle_to_add.data['smooth_DoLP'])
            self.data['smooth_AoLP'] = np.append(
                self.data['smooth_AoLP'], bottle_to_add.data['smooth_AoLP'])

            self.data['V'] = np.append(self.data['V'], bottle_to_add.data['V'])
            self.data['Vcos'] = np.append(self.data['Vcos'], bottle_to_add.data['Vcos'])
            self.data['Vsin'] = np.append(self.data['Vsin'], bottle_to_add.data['Vsin'])
            self.data['I0'] = np.append(self.data['I0'], bottle_to_add.data['I0'])
            self.data['DoLP'] = np.append(self.data['DoLP'], bottle_to_add.data['DoLP'])
            self.data['AoLP'] = np.append(self.data['AoLP'], bottle_to_add.data['AoLP'])

            self.V_average = np.average(self.data['V'])
            self.Vcos_average = np.average(self.data['Vcos'])
            self.Vsin_average = np.average(self.data['Vsin'])
            self.I0_average, self.DoLP_average, self.AoLP_average = Rotation.GetLightParameters(
                self.V_average, self.Vcos_average, self.Vsin_average)

            self.std_I0 = np.append(self.std_I0, bottle_to_add.std_I0)
            self.std_smooth_I0 = np.append(
                self.std_smooth_I0, bottle_to_add.std_smooth_I0)
            self.std_DoLP = np.append(self.std_DoLP, bottle_to_add.std_DoLP)
            self.std_smooth_DoLP = np.append(
                self.std_smooth_DoLP, bottle_to_add.std_smooth_DoLP)
            self.std_AoLP = np.append(self.std_AoLP, bottle_to_add.std_AoLP)
            self.std_smooth_AoLP = np.append(
                self.std_smooth_AoLP, bottle_to_add.std_smooth_AoLP)

            self.SetUnifyAngles()

            self.tail_jump = self.data['Times'][-1]

        self.graph_angle_shift = max(
            self.graph_angle_shift, bottle_to_add.graph_angle_shift)
        self.Geometry()

        return self


    def GetAltitude(self, filter = None):
        if filter is None:
            # print(self.filters)
            filter = self.filters

        altitude = 110

        if filter == "r": altitude = 220
        elif filter == "v": altitude = 110
        elif filter == "b" or "m": altitude = 100

        return altitude


    def TestResample(self):
        # f, axs = plt.subplots(3, sharex=True)
        
        N_pts = int(self.avg_dt)
        N_rot = len(self.data['I0'])
        filter_angles = np.linspace(0, N_rot * 2*np.pi, N_rot * N_pts) #List of angles (rad) for the polarising filter between 0 and 2π.

        print("Resampling signal")

        resamp_signal = CubicSpline(np.array(self.data['Times'].copy().dt.total_seconds()), np.array(self.data['I0'].copy()))
        resampled_t = np.linspace(self.data['Times'].iloc[0].total_seconds(), self.data['Times'].iloc[-1].total_seconds(), N_pts*N_rot, endpoint=True)
        resamp_signal = resamp_signal(resampled_t) / 2.
        

        def GetPola(V, Vcos, Vsin):
            """Given V, Vcos, Vsin, returns the initial intensity, DoLP and AoLP. This method is shared for spp and ptcu. It is also a static method, that can be called outside of the object. This way it can be used everywhere, each time you need I0, DoLP, AoLP to decrease the chances of making a mistake."""
            I0 = 2 * V
            DoLP = 2 * np.sqrt(Vcos**2 + Vsin**2) / V * 100
            AoLP = np.arctan2(Vsin, Vcos) / 2
            return abs(I0), abs(DoLP), AoLP

        def GetV(fake_signal):
            """return the average of a signal over 1 rotation"""
            integral = sum(fake_signal)
            return integral / N_pts

        def GetVcos(fake_signal):
            """Return the average of V*cos(2*theta) over 1 rotation"""
            x = fake_signal * np.cos(2 * filter_angles)
            return sum(x) / N_pts

        def GetVsin(fake_signal):
            """Return the average value of -V*sin(2*theta) over 1 rotation"""
            y = - fake_signal * np.sin(2 * filter_angles)
            return sum(y) / N_pts


        def GetStokesTime(fake_signal):
            """return the stokes parameters of a signal over N_rot rotation"""
            V = np.zeros(N_rot)
            Vcos = np.zeros(N_rot)
            Vsin = np.zeros(N_rot)
            for ir in range(N_rot):
                start_rot  = ir * N_pts
                end_rot    = start_rot + N_pts
                tmp_signal = fake_signal[start_rot:end_rot]

                V[ir] = sum(tmp_signal) / N_pts
                Vcos[ir] = sum(tmp_signal * np.cos(2 * filter_angles[start_rot:end_rot])) / N_pts
                Vsin[ir] = sum(tmp_signal * np.sin(2 * filter_angles[start_rot:end_rot])) / N_pts
            return V, Vcos, Vsin


        V, Vcos, Vsin = GetStokesTime(resamp_signal)
        I_list, DoLP_list, AoLP_list = GetPola(V, Vcos, Vsin)

        self.data['Vcos'] -= Vcos
        self.data['Vsin'] -= Vsin
        self.data['I0'], self.data['DoLP'], self.data['AoLP'] = GetPola(self.data['V'], self.data['Vcos'], self.data['Vsin'])
        self.data['AoLP'] += self.AoLP_correction
        self.data['AoLP']  = SetAngleBounds(self.data['AoLP'], -np.pi / 2, np.pi / 2)

        self.data['smooth_V']    = self.GetSliddingAverage(self.data['V'],    self.data['Times'], self.smoothing_factor, self.smoothing_unit)
        self.data['smooth_Vcos'] = self.GetSliddingAverage(self.data['Vcos'], self.data['Times'], self.smoothing_factor, self.smoothing_unit)
        self.data['smooth_Vsin'] = self.GetSliddingAverage(self.data['Vsin'], self.data['Times'], self.smoothing_factor, self.smoothing_unit)
        
        # Calculate the smooth I0, DoLP and AoLP
        # for i in range(self.nb_smooth_rot):
        #     self.data['smooth_I0'].iloc[i], self.data['smooth_DoLP'].iloc[i], self.data['smooth_AoLP'].iloc[i] = Rotation.GetLightParameters(
        #         self.data['smooth_V'].iloc[i], self.data['smooth_Vcos'].iloc[i], self.data['smooth_Vsin'].iloc[i])
        self.data['smooth_I0'], self.data['smooth_DoLP'], self.data['smooth_AoLP'] = Rotation.GetLightParameters(
                self.data['smooth_V'], self.data['smooth_Vcos'], self.data['smooth_Vsin'])
        self.data['smooth_AoLP'] = self.data['smooth_AoLP'] + self.AoLP_correction

        self.data['smooth_AoLP'] = SetAngleBounds(
            self.data['smooth_AoLP'], -np.pi / 2, np.pi / 2)

        # fig, axs = plt.subplots(3, sharex = True)
        # axs[0].plot(resampled_t, resamp_signal*2, ".")
        # axs[0].plot(self.data['Times'].dt.total_seconds(), self.data['I0'])
        # axs[0].plot(self.data['I0'], ".")
        # axs[0].plot(I_list, ".")
        # axs[0].plot(dI)
        # axs[1].plot(self.data['DoLP'], ".")
        # axs[1].plot(DoLP_list, ".")
        # axs[1].plot(dD, ".")
        # axs[2].plot(self.data['AoLP']*RtoD, ".")
        # axs[2].plot(AoLP_list*RtoD, ".")
        # axs[2].plot(dA*RtoD)



    def GetLightParameters(V, Vcos, Vsin):
            """Given V, Vcos, Vsin, returns the initial intensity, DoLP and AoLP. This method is shared for spp and ptcu. It is also a static method, that can be called outside of the object. This way it can be used everywhere, each time you need I0, DoLP, AoLP to decrease the chances of making a mistake."""
            I0 = 2 * V
            DoLP = 2 * np.sqrt(Vcos**2 + Vsin**2) / V * 100
            AoLP = np.arctan2(Vsin, Vcos) / 2
            return abs(I0), abs(DoLP), AoLP
    GetLightParameters = staticmethod(GetLightParameters)


    def IsGood(self):
        """Check if the rotation is a "good" rotation. the Vcos Vsin shouldn't be too high (<10), The DoLP < 30%, and all important parameters are a number (are defined)."""
        # good = True

        indexes = self.data[  (self.data['I0']   > self.Imax) 
                            | (self.data['I0']   < self.Imin)
                            | (self.data['DoLP'] > self.DoLPmax)
                            | (self.data['DoLP'] < self.DoLPmin)].index

        self.nb_bad_rot = len(indexes)
        print(indexes)
        self.data.drop(indexes, inplace=True)

        if len(self.data) == 0:
            self.valid = False
            
        # if True in np.isnan([self.V, self.Vcos, self.Vsin, self.I0, self.DoLP, self.AoLP]) :
        #     # print(np.isnan([self.V, self.Vcos, self.Vsin, self.I0, self.DoLP, self.AoLP]))
        #     good = False

#####################################################################################
###                                    Petit Cru                                      ###
#####################################################################################


class PTCUBottle(Bottle):
    def __init__(self, folder_name, auto=True, line=1, from_txt=False):
        Bottle.__init__(self, folder_name, auto, line=line, from_txt=from_txt)

        self.SetInfoFromDataFileName(self.data_file_name, pwd_data)

        if auto:
            self.MiseEnBouteille()

    def SetInfoFromDataFileName(self, f, data_f):
        Bottle.SetInfoFromDataFileName(self, f, data_f)

        i = self.folders_index
        date_location = self.folders[i + 1].split("_")
        print(self.folders[i + 1], date_location)
        self.date, self.location = date_location[0], date_location[1]
        rest = self.folders[i + 2].split("_")
        print(i, date_location, self.date, self.location, rest)
        # print(rest)
        # Filters: r,v,b,m pour rouge, vert, bleu, mauve
        # 0: no filters_list
        # o: orange 620
        # X, Y: 620nm et 413nm
        filters_list = ["r", "v", "b", "m", "0", "o", "X", "Y", "t"]
        # nb_filters = 2
        # if self.instrument_name == "gdcu":
        #     nb_filters = 4

        for r in rest:
            if np.all([a in filters_list for a in r]) and r != 'rot':
                if self.instrument_name == "carmen" or self.instrument_name == "fake_ptcu":
                    self.filters = r
                    if self.filters[1] in [0, "0"]:
                        self.NoVref = True
                elif self.instrument_name in ["corbel", "gdcu"]:
                    print(rest, r, self.line)
                    self.filters = r[self.line - 1]
                # print("FILTERS:", r, self.filters)

            elif r[0] == "a":
                self.azimut = float(r[1:]) * DtoR
            elif r[0] == "e":
                self.elevation = float(r[1:]) * DtoR
            else:
                self.com += "_" + r

        print("SetInfoFromDataFileName", self.instrument_name, self.date,
              self.location, self.filters, self.azimut * RtoD, self.elevation * RtoD, self.com)

    def LoadData(self):


        if self.instrument_name == "carmen" or self.instrument_name == "corbel" or self.instrument_name == "gdcu":
            self.raw_data, self.config = self.LoadPTCUData()
            self.nb_rot = len(self.raw_data)
            if self.nb_rot == 0:
                self.valid = False
            # self.nb_data_per_rot = len(self.raw_data[0])
            # if self.jump_mode == "length" or self.tail_jump == 0:
            #     self.tail_jump = len(self.raw_data) - self.tail_jump
            
            self.data = self.raw_data.copy()
            self.data = self.data.set_index('IDProcessedData')
            
            self.data.rename(columns={'Time': 'Times', 'PMAvg': 'V', 'PMCosAvg': 'Vcos', 'PMSinAvg': 'Vsin', 'REFAvg': 'Vref'}, inplace=True)
            
            self.data['Times'] = pd.to_timedelta(self.data['Times'], unit='ms')

            self.data['V'] = np.absolute(self.data['V'])
            self.data['I0'], self.data['DoLP'], self.data['AoLP'] = self.GetLightParameters(self.data['V'], self.data['Vcos'], self.data['Vsin'])
            if 'Vref' in self.data.columns:
                self.data['Iref'] = 2 * self.data['Vref']


            # print(self.data)
            # print(self.data.columns)
            
            # for r in self.raw_data[:]:
            #     self.rotations.append(PTCURotation(r))


        elif self.instrument_name == "fake_ptcu":
            self.raw_data, self.config = LoadPTCUTestData(
                nb_rot=1000, I0=100, D=0.10, phi=45 * RtoD)
            # if self.jump_mode == "length" or self.tail_jump == 0:
            #     self.tail_jump = len(self.raw_data) - self.tail_jump
            for r in self.raw_data[:]:
                self.rotations.append(PTCURotation(r))

        
        self.SetConfig()
        self.SetTimes()

        
        

    def GetTimeFromDateTime(self, date):
        """Return list self.data['Times'] shifted so that 0 is date."""
        shift = self.DateTime(moment="start") - date
        return self.data['Times'] + shift

    def SetTimeFromDateTime(self, date):
        self.date_time = date
        self.data['Times'] = self.GetTimeFromDateTime(date)

    def SetTimes(self):
        # try:
        #     self.time_stamp = str(self.config["Timestamp"]).split("\"")[1] #When importing config file from mysql workbench, timestamp is between "", but when using python script, there is no more ""
        # except:
        #     self.time_stamp = str(self.config["Timestamp"]).split("\'")[1]
        self.time_stamp = self.config["Timestamp"]
        print(self.time_stamp)
        # Date and time of first rotation (before deleting head_jump)

        self.datetime = dt.datetime.strptime(
            self.time_stamp, "%Y-%m-%d %H:%M:%S")
        print(self.datetime, self.datetime.strftime("%H:%M:%S"))
        # Time in sec since EPOCH
        # self.time = self.datetime.timestamp()

    # def DateTime(self, moment="start", delta=None, format="LT"):
    #     if not delta:
    #         delta = dt.timedelta(seconds=0)
    #
    #     norm = dt.timedelta(seconds=0)
    #     if format == "LT" and self.config_time_format == "UT":
    #         norm = self.time_zone
    #     elif format == "UT" and self.config_time_format == "LT":
    #         norm = - self.time_zone
    #
    #     if moment == "start": #Datetime of first good rotation
    #         return self.data['Times'][0] + self.datetime + self.head_jump + delta + norm
    #         # return self.rotations[0].time + self.datetime + self.head_jump + delta
    #     elif moment == "end": #Datetime of last good rotation
    #         return self.data['Times'][-1] + self.datetime + self.head_jump + delta + norm
    #         # return self.rotations[-1].time + self.datetime + self.head_jump + delta
    #     elif moment == "config": #Datetime written in the config file (== first bad rotation, before head_jump)
    #         return self.datetime + delta + norm

    def CleanRotations(self):
        
        self.nb_rot = len(self.data)

        print("Cleaning bad points...")

        # Put every rotation time in sec since first rotation
        norm = self.data['Times'].iloc[0]
        self.data['Times'] -= norm

        # Add 24h when next rotation is earlier than the previous one. Happens sometimes when we change the date
        one_day = dt.timedelta(days=1)
        for i in range(1, self.nb_rot):
            while self.data['Times'].iloc[i] < self.data['Times'].iloc[i - 1]:
                self.data['Times'].iloc[i] += one_day

        # Get and set the head and tail jumps
        # if self.jump_mode == "time": # and self.jump_unit in ["seconds", "minutes", "hours"]:
        #     self.head_jump = time.timedelta(seconds=float(self.input_parameters["head_jump"]))
        #     self.tail_jump = time.timedelta(seconds=float(self.input_parameters["tail_jump"]))
        #
        #     if self.jump_unit == "minutes":
        #         self.head_jump = time.timedelta(minutes=float(self.input_parameters["head_jump"]))
        #         self.tail_jump = time.timedelta(minutes=float(self.input_parameters["tail_jump"]))
        #     elif self.jump_unit == "hours":
        #         print("DEBUG jumps", [r.time.total_seconds() for r in self.rotations[:10]])
        #         print(time.timedelta(hours=float(self.input_parameters["head_jump"])), time.timedelta(hours=float(self.input_parameters["tail_jump"])))
        #         self.head_jump = time.timedelta(hours=float(self.input_parameters["head_jump"]))
        #         self.tail_jump = time.timedelta(hours=float(self.input_parameters["tail_jump"]))
        #
        #     if self.jump_mode == "length" or self.tail_jump.seconds == 0.:
        #         self.tail_jump = self.rotations[-1].time - self.tail_jump

        # Delete all rotations before the head and after the tail jump
        # self.rotations = [r for r in self.rotations if self.head_jump <= r.time <= self.tail_jump]
        # print(len(self.rotations), "good rotations in", self.nb_rot, ";", self.nb_rot - len(self.rotations), "deleted because of time jumps.")
        # self.nb_rot = len(self.rotations)

        indexes = self.data[  (self.data['Times']   < self.head_jump) 
                            | (self.data['Times']   > self.tail_jump)].index
        self.data.drop(indexes, inplace=True)
        self.nb_rot = len(self.data)

        # print(self.head_jump, self.tail_jump, self.rotations[0].time, self.rotations[-1].time)

        # Put every rotation time in sec since first rotation (after deleting head_jump)
        norm = self.data['Times'].iloc[0]
        self.data['Times'] -= norm
        # for i in range(len(self.rotations)):
        #     self.rotations[i].time -= norm

        # time = np.array([r.time.total_seconds()    for r in self.rotations])
        # all_V       = np.array([r.V    for r in self.rotations])
        # dV_dt = np.gradient(all_V, time )
        # avg_V = np.average(all_V)
        # print(avg_V)
        # self.rotations = [r for ir, r in enumerate(self.rotations) if abs(dV_dt[ir]) < 0.01 * avg_V]
        #
        # plt.plot(time, dV_dt)
        # plt.show()

        if self.CorrectDensity():
            print("Correcting density")
            # print(self.density_start, self.density_end)
            # print(self.rotations[0].time, self.rotations[-1].time)
            for ir, r in enumerate(self.rotations):
                if self.density_start <= r.time < self.density_end:
                    self.data['V'].iloc[ir]    *= self.density_factor
                    self.data['Vcos'].iloc[ir] *= self.density_factor
                    self.data['Vsin'].iloc[ir] *= self.density_factor
                    self.data['I0'].iloc[ir]   *= self.density_factor


        # self.rotations = [r for r in self.rotations if r.IsGood(
            # Imin=self.Imin, Imax=self.Imax, DoLPmin=self.DoLPmin, DoLPmax=self.DoLPmax)]
        self.IsGood()
        print(len(self.data), "good rotations in", self.nb_rot, ";",
              self.nb_rot - len(self.rotations), "deleted because of invalid data.")
        self.nb_rot = len(self.data)


    def SaveTXT(self):
        print("Saving as .txt in", self.data_file_name
              + "/" + self.saving_name + '_results.txt')

        # Default Format
        # times = [t.total_seconds(
        # ) * 1000 for t in self.GetTimeFromDateTime(self.DateTime(moment="config"))]
        # data = pd.DataFrame(np.array([times, self.data['V'], self.data['Vcos'], self.data['Vsin'], self.data['I0'], self.data['DoLP'], self.data['AoLP'], self.data['data['smooth_V']'], self.data['smooth_Vcos'], self.data['smooth_Vsin'], self.data['smooth_I0'],
        #                     self.data['smooth_DoLP'], self.data['smooth_AoLP'], self.std_I0, self.std_DoLP, self.std_AoLP, self.std_smooth_I0, self.std_smooth_DoLP, self.std_smooth_AoLP, self.all_SN, self.smooth_SN]).transpose())
        # data.columns = ["time", "V", "Vcos", "Vsin", "I0", "DoLP", "AoLP", "SV", "SVcos", "SVsin", "SI0",
        #                 "SDoLP", "SAoLP", "errI0", "errDoLP", "errAoLP", "errSI0", "errSDoLP", "errSAoLP", "SN", "SSN"]

        # Simple format for Jean
        # times = [t.total_seconds() / 3600. for t in self.GetTimeFromDateTime(self.DateTime(moment="midnight", format="LT"))]
        # data = pd.DataFrame(np.array([times, self.data['smooth_I0'], self.data['smooth_DoLP'], self.data['smooth_AoLP'], self.std_smooth_I0, self.std_smooth_DoLP, self.std_smooth_AoLP]).transpose())
        # data.columns = ["time","SI0","SDoLP","SAoLP","errSI0","errSDoLP","errSAoLP"]

        # data.to_csv(self.data_file_name + "/" + self.saving_name + '_results.txt', sep="\t", index=False)
        self.data.to_csv(self.data_file_name + "/" + self.saving_name + '_results.txt', sep="\t", index=False)

    def SaveHDF5(self):
        print("Saving as .hdf5 in", self.data_file_name
              + "/" + self.saving_name + '_results.hdf5')

        # times = [t.total_seconds(
        # ) * 1000 for t in self.GetTimeFromDateTime(self.DateTime(moment="config"))]

        # data = [times, self.data['V'], self.data['Vcos'], self.data['Vsin'], self.data['I0'], self.data['DoLP'], self.data['AoLP'], self.data['data['smooth_V']'], self.data['smooth_Vcos'], self.data['smooth_Vsin'], self.data['smooth_I0'],
        #         self.data['smooth_DoLP'], self.data['smooth_AoLP'], self.std_I0, self.std_DoLP, self.std_AoLP, self.std_smooth_I0, self.std_smooth_DoLP, self.std_smooth_AoLP, self.all_SN, self.smooth_SN]
        # columns = ["time", "V", "Vcos", "Vsin", "I0", "DoLP", "AoLP", "SV", "SVcos", "SVsin", "SI0",
        #            "SDoLP", "SAoLP", "errI0", "errDoLP", "errAoLP", "errSI0", "errSDoLP", "errSAoLP", "SN", "SSN"]

        # print(self.data)
        with h5.File(self.data_file_name + "/" + self.saving_name + '_results.hdf5', 'w') as f:

            f.create_group("/data")
            # print(self.data.columns)
            for ic, c in enumerate(self.data.columns):
                # print(ic, c)
                if c in ['time_deltas', 'Times']:
                    f.create_dataset("/data/" + c,  data=self.data[c].dt.total_seconds())
                elif c == 'Comment':
                    pass
                else:
                    # print(c, type(self.data[c].iloc[0]))
                    f.create_dataset("/data/" + c,  data=self.data[c])

            # print(self.config)
            # f.create_group("config")
            for k, v in self.config.items():
                f["data"].attrs[k] = v
                # f.create_dataset("/config/" + str(k), data = v)

    def LoadPTCUData(self):
        """Return an array of data and a dictionnary of config from a list of data files. The first is the data, the second is the config. Each line of the data is a rotation with 6 entries, the config dict is all the parameters of the observation, common to all rotations."""
        # DATA FILE
        # 0:     IDProcessedData,
        # 1:     IDConfiguration,
        # 2:     Time since begining of obseravtion (timestamp) in milliseconds,
        # 3:    PMAvg,
        # 4:    PMCosAvg,
        # 5:    PMSinAvg
        # 6:    IDTemperatures
        # 7:    TempPM
        # 8:    TempOptical
        # 9:    TempAmbiant
        # 10:    IDLiveComment
        # 11:    Comment
        # 12:    live_commentscol
        # CONFIG FILE
        # 0:        IDConfiguration,
        # 1:        Timestamp,
        # 2-3:        CM_Latitude, CM_Longitude,
        # 4-5:        CM_Elevation, CM_Azimuth,
        # 6:        CM_PolarizerSpeed,
        # 7:        CM_Tilt,
        # 8:        CM_Usage,
        # 9:        CM_Comments,
        # 10:        IS_PolarizerOffset1,
        # 11:        IS_PolarizerOffset2,
        # 12:        IS_PolarizerOffset3,
        # 13:        IS_PolarizerOffset4,
        # 14:        IS_ConverterOffset4,
        # 15:        IS_ConverterOffset3,
        # 16:        IS_ConverterOffset2,
        # 17:        IS_ConverterOffset1,
        # 18:        IS_EngineTicksPerTour,
        # 19:        IS_MotorGearedRatio,
        # 20:        IS_QuantumEfficiency,
        # 21:        IS_IpAddress

        # For old (< V3) csv files
        file_names = self.data_file_name
        line = self.line

        print("LOADING PTCU data: line", line)
        data_path = file_names
        # Loading the data from a binary file.
        # Output is one big 1-D array of length nb_data_tot = nb_toy * nb_data_per_rot
        data_file = data_path + "/data" + str(line) + ".csv"
        config_file = data_path + "/config.csv"

        # print(data_file, config_file)

        try:
            with open(data_file, "r") as f:
                first_line = f.readlines()[1]
            d = GetDelimiter(first_line)
        except:
            d = ""

        # raw_data = np.genfromtxt(data_file, delimiter=d, skip_header=1)
        
        raw_data = pd.read_csv(data_file, sep=",", na_values=['None'])
        # print(raw_data)

        # if self.instrument_name in ["corbel", "ptcu_v2", "gdcu"]:
        #     array_type = [('IDConfiguration',float),('Timestamp','S100'), ('CM_ID', 'S100'),('CM_Latitude',float),('CM_Longitude',float),('CM_Elevation',float),('CM_Azimuth',float),('CM_PolarizerSpeed',float),('CM_Tilt',float),('CM_Usage','S100'),('CM_Comments','S100'),('IS_PolarizerOffset1',float),('IS_PolarizerOffset2',float),('IS_PolarizerOffset3',float),('IS_PolarizerOffset4',float),('IS_ConverterGain1',float),('IS_ConverterGain2',float),('IS_ConverterGain3',float),('IS_ConverterGain4',float),("IS_ConverterOffset1", float),("IS_ConverterOffset2", float),("IS_ConverterOffset3", float),("IS_ConverterOffset4", float),('IS_MotorGearedRatio',float),('IS_QuantumEfficiency',float), ('IS_IpAddress',float)]
        # else:
        #     array_type = [('IDConfiguration',float),('Timestamp','S100'),('CM_Latitude',float),('CM_Longitude',float),('CM_Elevation',float),('CM_Azimuth',float),('CM_PolarizerSpeed',float),('CM_Tilt',float),('CM_Usage','S100'),('CM_Comments','S100'),('IS_PolarizerOffset1',float),('IS_PolarizerOffset2',float),('IS_PolarizerOffset3',float),('IS_PolarizerOffset4',float),("IS_ConverterOffset4", float),("IS_ConverterOffset3", float),("IS_ConverterOffset2", float),("IS_ConverterOffset1", float),('IS_EngineTicksPerTour',float),('IS_MotorGearedRatio',float),('IS_QuantumEfficiency',float), ('IS_IpAddress',float)]
        #
        # print("DEBUG", config_file, len(array_type))
        # configuration = np.genfromtxt(config_file, dtype=array_type, delimiter=",", skip_header=1)

        conf_df = pd.read_csv(config_file, sep=",")
        configuration = dict()
        # print(config_file)
        print('configuration', configuration)
        for n in conf_df.columns:
            configuration[n] = conf_df[n][0]

        print('configuration', configuration)

        # For new hdf5 files
        # file_names = self.data_file_name
        # line = self.line
        # print("LOADING PTCU data: line", line)
        # data_path = file_names
        #
        # data_file = data_path
        # config_file = data_path
        #
        # file = h5py.File(data_file, "r")
        # data = file[f"Data/Channel_{line-1}"]
        # raw_data = np.array([data["Time"], data["PMAvg"], data["PMCosAvg"], data["PMSinAvg"]])
        # configuration  = {'Timestamp': file.attrs["Timestamp"],
        # 'CM_Latitude': file.attrs["CM_Latitude"],
        # 'CM_Longitude': file.attrs["CM_Longitude"],
        # 'CM_Elevation': data.attrs["CM_Elevation"],
        # 'CM_Azimuth': data.attrs["CM_Azimuth"],
        # 'CM_PolarizerSpeed': data.attrs["CH_PolarizerSpeed"],
        # 'FL_Wavelength': data.attrs["FL_Wavelength"],
        # 'FL_Width': data.attrs["FL_Width"],
        # # 'CM_Tilt': ,
        # 'CM_Usage': file.attrs["CM_Usage"],
        # 'CM_Comments': file.attrs["CM_Comments"],
        # f'IS_ConverterGain{line}': data.attrs["CH_ConverterGain"],
        # f'IS_PolarizerOffset{line}': data.attrs["CH_PolarizationOffset"],
        # f"IS_ConverterOffset{line}": data.attrs["CH_ConverterOffset"],
        # # 'IS_EngineTicksPerTour':,
        # # 'IS_MotorGearedRatio':,
        # 'IS_QuantumEfficiency': data.attrs["PM_QuantumEfficiency"],
        # # 'IS_IpAddress':,
        # }

        print("PTCU DATA loaded")

        return raw_data, configuration

    def LoadFromHDF5(self):
        file_names = self.data_file_name
        line = self.line
        print("LOADING PTCU data: line", line)
        data_path = file_names

        data_file = data_path
        config_file = data_path

        file = h5py.File(data_file, "r")
        data = file[f"Data/Channel_{line-1}"]

        self.data['Times'] = np.array([dt.timedelta(milliseconds=t)
                                  for t in data["Time"]])
        self.nb_rot = len(self.data['Times'])

        self.SetJumps()

        self.data['smooth_V'] = np.array([d for d, t in zip(
            data["smooth_V"], self.data['Times']) if self.head_jump <= t < self.tail_jump])
        self.data['smooth_Vcos'] = np.array([d for d, t in zip(
            data["smooth_Vcos"], self.data['Times']) if self.head_jump <= t < self.tail_jump])
        self.data['smooth_Vsin'] = np.array([d for d, t in zip(
            data["smooth_Vsin"], self.data['Times']) if self.head_jump <= t < self.tail_jump])

        self.nb_smooth_rot = len(self.data['smooth_V'])

        # Calculate the smooth I0, DoLP and AoLP
        self.data['smooth_I0'] = np.array([d for d, t in zip(
            data["Smooth_I0"], self.data['Times']) if self.head_jump <= t < self.tail_jump])
        self.data['smooth_DoLP'] = np.array([d for d, t in zip(
            data["Smooth_DoLP"], self.data['Times']) if self.head_jump <= t < self.tail_jump])
        self.data['smooth_AoLP'] = np.array([d for d, t in zip(
            data["data['smooth_AoLP']"], self.data['Times']) if self.head_jump <= t < self.tail_jump])

        self.data['V'] = np.array([d for d, t in zip(
            data["PMAvg"], self.data['Times']) if self.head_jump <= t < self.tail_jump])
        self.data['Vcos'] = np.array([d for d, t in zip(
            data["PMCosAvg"], self.data['Times']) if self.head_jump <= t < self.tail_jump])
        self.data['Vsin'] = np.array([d for d, t in zip(
            data["PMSinAvg"], self.data['Times']) if self.head_jump <= t < self.tail_jump])
        self.data['I0'] = np.array([d for d, t in zip(
            data["I0"], self.data['Times']) if self.head_jump <= t < self.tail_jump])
        self.data['DoLP'] = np.array([d for d, t in zip(
            data["DoLP"], self.data['Times']) if self.head_jump <= t < self.tail_jump])
        self.data['AoLP'] = np.array([d for d, t in zip(
            data["AoLP"], self.data['Times']) if self.head_jump <= t < self.tail_jump])

        self.std_I0 = np.array([d for d, t in zip(
            data["errI0"], self.data['Times']) if self.head_jump <= t < self.tail_jump])
        self.std_smooth_I0 = np.array([d for d, t in zip(
            data["errSI0"], self.data['Times']) if self.head_jump <= t < self.tail_jump])

        self.std_DoLP = np.array([d for d, t in zip(
            data["errDoLP"], self.data['Times']) if self.head_jump <= t < self.tail_jump])
        self.std_smooth_DoLP = np.array([d for d, t in zip(
            data["errSDoLP"], self.data['Times']) if self.head_jump <= t < self.tail_jump])

        self.std_AoLP = np.array([d for d, t in zip(
            data["errAoLP"], self.data['Times']) if self.head_jump <= t < self.tail_jump])
        self.std_smooth_AoLP = np.array([d for d, t in zip(
            data["errSAoLP"], self.data['Times']) if self.head_jump <= t < self.tail_jump])

        norm = self.head_jump
        self.data['Times'] = np.array(
            [t - norm for t in self.data['Times'] if self.head_jump <= t < self.tail_jump])

        self.valid = True
        _, self.config = self.LoadPTCUData()

        self.SetTimes()
        self.SetSmoothLists()

        self.SetTimeFromDateTime(self.DateTime(moment="start"))

        self.graph_angle_shift = 0
        for a in self.data['smooth_AoLP']:
            if a > np.pi / 2.:
                self.graph_angle_shift = 1
                break
            if a < 0:
                break

    def LoadFromTxt(self):
        file_name = self.data_file_name + "/" + self.saving_name + '_results.txt'

        data = pd.read_csv(file_name, delimiter="\t")
        data.columns = [c.replace("#", "").replace(
            "(", "").replace(")", "").strip() for c in data.columns]
        print(data.columns)
        time_name = data.columns[0]

        self.data['Times'] = np.array([dt.timedelta(milliseconds=t)
                                  for t in data[time_name]])
        self.nb_rot = len(self.data['Times'])

        self.SetJumps()

        self.data['smooth_V'] = np.array([d for d, t in zip(
            data["SV"], self.data['Times']) if self.head_jump <= t < self.tail_jump])
        self.data['smooth_Vcos'] = np.array([d for d, t in zip(
            data["SVcos"], self.data['Times']) if self.head_jump <= t < self.tail_jump])
        self.data['smooth_Vsin'] = np.array([d for d, t in zip(
            data["SVsin"], self.data['Times']) if self.head_jump <= t < self.tail_jump])

        self.nb_smooth_rot = len(self.data['smooth_V'])

        # Calculate the smooth I0, DoLP and AoLP
        self.data['smooth_I0'] = np.array([d for d, t in zip(
            data["SI0"], self.data['Times']) if self.head_jump <= t < self.tail_jump])
        self.data['smooth_DoLP'] = np.array([d for d, t in zip(
            data["SDoLP"], self.data['Times']) if self.head_jump <= t < self.tail_jump])
        self.data['smooth_AoLP'] = np.array([d for d, t in zip(
            data["SAoLP"], self.data['Times']) if self.head_jump <= t < self.tail_jump])

        self.data['V'] = np.array([d for d, t in zip(
            data["V"], self.data['Times']) if self.head_jump <= t < self.tail_jump])
        self.data['Vcos'] = np.array([d for d, t in zip(
            data["Vcos"], self.data['Times']) if self.head_jump <= t < self.tail_jump])
        self.data['Vsin'] = np.array([d for d, t in zip(
            data["Vsin"], self.data['Times']) if self.head_jump <= t < self.tail_jump])
        self.data['I0'] = np.array([d for d, t in zip(
            data["I0"], self.data['Times']) if self.head_jump <= t < self.tail_jump])
        self.data['DoLP'] = np.array([d for d, t in zip(
            data["DoLP"], self.data['Times']) if self.head_jump <= t < self.tail_jump])
        self.data['AoLP'] = np.array([d for d, t in zip(
            data["AoLP"], self.data['Times']) if self.head_jump <= t < self.tail_jump])

        self.std_I0 = np.array([d for d, t in zip(
            data["errI0"], self.data['Times']) if self.head_jump <= t < self.tail_jump])
        self.std_smooth_I0 = np.array([d for d, t in zip(
            data["errSI0"], self.data['Times']) if self.head_jump <= t < self.tail_jump])

        self.std_DoLP = np.array([d for d, t in zip(
            data["errDoLP"], self.data['Times']) if self.head_jump <= t < self.tail_jump])
        self.std_smooth_DoLP = np.array([d for d, t in zip(
            data["errSDoLP"], self.data['Times']) if self.head_jump <= t < self.tail_jump])

        self.std_AoLP = np.array([d for d, t in zip(
            data["errAoLP"], self.data['Times']) if self.head_jump <= t < self.tail_jump])
        self.std_smooth_AoLP = np.array([d for d, t in zip(
            data["errSAoLP"], self.data['Times']) if self.head_jump <= t < self.tail_jump])

        self.all_SN = np.array([d for d, t in zip(
            data["SN"], self.data['Times']) if self.head_jump <= t < self.tail_jump])
        self.smooth_SN = np.array([d for d, t in zip(
            data["SSN"], self.data['Times']) if self.head_jump <= t < self.tail_jump])

        norm = self.head_jump
        self.data['Times'] = np.array(
            [t - norm for t in self.data['Times'] if self.head_jump <= t < self.tail_jump])

        self.valid = True
        _, self.config = self.LoadPTCUData()

        self.SetTimes()
        self.SetSmoothLists()

        self.SetTimeFromDateTime(self.DateTime(moment="start"))

        self.graph_angle_shift = 0
        for a in self.data['smooth_AoLP']:
            if a > np.pi / 2.:
                self.graph_angle_shift = 1
                break
            if a < 0:
                break

        # if self.instrument_name in ["corbel", "gdcu"]:
        #     self.data['AoLP']         -= 40 * DtoR
        #     self.data['smooth_AoLP']     -= 40 * DtoR


#####################################################################################
###                                    SPPBottle                                      ###
#####################################################################################


class SPPBottle(Bottle):
    def __init__(self, folder_name, auto=True, from_txt=False):
        Bottle.__init__(self, folder_name, auto, from_txt=from_txt)
        try:
            self.SetInfoFromDataFileName(self.data_file_name, pwd_data)
        except:
            self.location = ""
            print(
                "WARNING: Folders not standard. Could not retrieve info from folders name.")

        self.filters = ["r"]

        if auto:
            self.MiseEnBouteille()

    def SetInfoFromDataFileName(self, f, data_f):
        Bottle.SetInfoFromDataFileName(self, f, data_f)
        i = self.folders_index
        if len(self.folders[i + 1].split("_")) == 2:
            self.date, self.location = self.folders[i + 1].split("_")[0], self.folders[i + 1].split("_")[1]
        else:
            self.date = self.folders[i + 1]
        # print(self.instrument_name, self.date, self.location, self.filters, self.azimut*RtoD, self.elevation*RtoD, self.com)

        rest = self.folders[i + 2].split("_")
        print(rest)
        for r in rest:
            if len(r) == 2 and ((r[0] and r[1]) in ["r", "v", "b", "m", "0", "o"]):
                self.filters = r
            elif r[0] == "a":
                self.azimut = float(r[1:]) * DtoR
            elif r[0] == "e":
                self.elevation = float(r[1:]) * DtoR
            else:
                self.com += "_" + r

        # self.filters = ["r"]

        # if len(rest) <= 2:
        #     azimut, elevation = float(rest[0][1:]) * DtoR, float(rest[1][1:]) * DtoR
        # elif len(rest) > 2:
        #     filters, azimut, elevation = rest[0], float(rest[1][1:]) * DtoR, float(rest[2][1:]) * DtoR
        # else:
        #     com = folders[i+2].split("_")

        print(self.instrument_name, self.date, self.location, self.filters,
              self.azimut * RtoD, self.elevation * RtoD, self.com)

    def LoadData(self):
        # Loading the data files
        if self.instrument_name[:3] == "spp":
            self.raw_data = self.LoadSPPData(self.data_file_name)
            # if self.jump_mode == "length" or self.tail_jump.seconds == 0:
            #     self.tail_jump = len(self.raw_data) - self.tail_jump.seconds
            for r in self.raw_data[:]:
                self.rotations.append(SPPRotation(
                    r, instrument=self.instrument_name))

        elif self.instrument_name == "fake_spp":
            self.raw_data = self.LoadSPPTestData(1000, 100, 0.01, 45 * RtoD)
            # if self.jump_mode == "length" or self.tail_jump.seconds == 0:
            #     self.tail_jump = len(self.raw_data) - self.tail_jump.seconds
            for r in self.raw_data[:]:
                self.rotations.append(SPPRotation(r, fake="True"))

        Bottle.LoadData(self)

    # def DateTime(self, moment="start"):
    #     if moment == "start":
    #         return self.rotations[0].datetime
    #     if moment == "end":
    #         return self.rotations[-1].datetime

    def SetTimes(self):
        self.time_stamp = self.rotations[0].time_stamp
        # Time of obs in sec since EPOCH (except if before April 2019 -> in sec since midnight)
        self.time = self.rotations[0].time
        # Date and time of first rotation (before deleting head_jump)
        self.datetime = self.rotations[0].datetime
        print("DEBUG SetTimes: ", self.time_stamp, self.time, self.datetime)
        # self.datetime = time.strptime(self.time_stamp, "%Y%m%d")
        # self.time = time.mktime(self.datetime)

    def CleanRotations(self):
        self.nb_rot = len(self.rotations)

        for i in range(1, len(self.rotations)):
            while self.rotations[i].time < self.rotations[i - 1].time:
                self.rotations[i].time += dt.timedelta(seconds=24 * 3600)

        # and self.jump_unit in ["seconds", "minutes", "hours"]:
        if self.jump_mode == "time":
            self.head_jump = dt.timedelta(
                seconds=float(self.input_parameters["head_jump"]))
            self.tail_jump = dt.timedelta(
                seconds=float(self.input_parameters["tail_jump"]))

            if self.jump_unit == "minutes":
                self.head_jump = dt.timedelta(
                    minutes=float(self.input_parameters["head_jump"]))
                self.tail_jump = dt.timedelta(
                    minutes=float(self.input_parameters["tail_jump"]))
            elif self.jump_unit == "hours":
                self.head_jump = dt.timedelta(
                    hours=float(self.input_parameters["head_jump"]))
                self.tail_jump = dt.timedelta(
                    hours=float(self.input_parameters["tail_jump"]))

            if self.jump_mode == "length" or self.tail_jump.seconds == 0.:
                self.tail_jump = self.rotations[-1].time - self.tail_jump

            # print(self.head_jump, self.tail_jump, self.rotations[0].time, self.rotations[-1].time)
            print("Nb rotations before cleaning:", self.nb_rot)

            self.rotations = [
                r for r in self.rotations if self.head_jump <= r.time < self.tail_jump]

            print(self.nb_rot - len(self.rotations),
                  "bad rotations in", self.nb_rot)
            self.nb_rot = len(self.rotations)
            print(self.nb_rot, "good rotations left after jump cuts.")

            print(self.head_jump, self.tail_jump,
                  self.rotations[0].time, self.rotations[-1].time)

        self.rotations = [r for r in self.rotations if r.IsGood(
            Imin=self.Imin, Imax=self.Imax, DoLPmin=self.DoLPmin, DoLPmax=self.DoLPmax)]
        # print(self.Imin, self.Imax, self.DoLPmin, self.DoLPmax)

        print(self.nb_rot - len(self.rotations),
              "bad rotations in", self.nb_rot)
        self.nb_rot = len(self.rotations)
        print(self.nb_rot, "good rotations left after good rot cuts.")

    def SaveTXT(self):
        print("Saving as .txt in", self.data_file_name
              + "/" + self.saving_name + '_results.txt')
        times = [t.total_seconds() for t in self.data['Times']]
        np.savetxt("/".join(self.data_file_name.split("/")[:-1]) + "/" + self.saving_name + '_results.txt',                    np.array([times, self.data['smooth_V'], self.data['smooth_Vcos'], self.data['smooth_Vsin'],
                   self.data['smooth_I0'], self.data['smooth_DoLP'], self.data['smooth_AoLP']]).transpose(), delimiter="\t", header="time (ms)\tV\tVcos\tVsin\tI0 (mV)\tDoLP (percent)\tAoLP (deg)")

    def LoadSPPData(self, file_names):
        """Return a data array from a given data file (1 or more). Each line is a rotation, each being a list of data as written below."""
        # 0-5:         year, month, day, hour, min, sec
        # 6-9:         voltage, positive voltage, negative voltage, motor voltage
        # 10-11:     Temperature of base, of filter
        # 12-13:    elevation, azimuth (geog coord)
        # 14-93:    angle of polariser (degrees)
        # 94-173:    canal pola
        # 174-253:    canal total
        # 254-259:    V, Vcos, Vsin, I0, DoLP, AoLP
        nb_data_per_rot = 254
        # raw_data = np.array([])
        # data_path = "/home/bossel/These/Analysis/data/spp/"
        # Loading the data from a binary file.
        # Output is one big 1-D array of length nb_data_tot = nb_rot * nb_data_per_rot
        # for f in file_names:
        #     f = data_path + f + ".spp"
        #     raw_data = np.fromfile(f, dtype = np.int16)
        #     # raw_data = np.concatenate((raw_data, np.fromfile(f, dtype=np.int16)), axis=0)
        f = file_names + ".spp"
        raw_data = np.fromfile(f, dtype=np.int16)
        # raw_data = np.concatenate((raw_data, np.fromfile(f, dtype=np.int16)), axis=0)
        nb_data_tot = len(raw_data)
        print("nb_data_tot", nb_data_tot)

        # print(raw_data[1001 * 254 : 1002 * 254])
        # reshape the array to be (nb_pts, nb_data per pts)
        return raw_data.reshape(int(nb_data_tot / nb_data_per_rot), nb_data_per_rot)

    def LoadSPPTestData(self, nb_rot=100, I0=100, D=0.1, phi=0, file_name=''):
        """Return a test data array with the time, angle and fake voltage mesured every point as SPP would. The fake signal is a function of the form A*cos(theta-phi)+B where theta is the angle of the polariser and phi the angle of the polarisation"""

        # I0:  Intensity of the incoming light of the aurora (constant through time) arbitrary units
        # D:  Degree of polarized light (%)
        # phi:  Angle of the polarisation with respect to SPP (degrees)

        # Do not take into account the fact taht theta-phi cannot be > 90, because it would only change a sign, but it is squared in the end. So no difference. But I keep the ITheo function in case.
        # I = lambda theta, I0, D, phi: D * I0 * np.cos((theta - phi)*DtoR)**2 + (1-D) * I0 / 2 + D/5*np.random.normal(0, 1, theta.shape)
        def I(theta, I0, D, phi): return I0 + \
            np.random.normal(0, I0 * D, theta.shape)

        data = np.zeros((nb_rot, 254))
        da = 360. / 80
        angles = np.linspace(da, 360., 80) * DtoR
        dolps = np.linspace(0, 1., nb_rot)
        phis = np.linspace(-90, 90, nb_rot)
        for i in range(nb_rot):
            data[i][5] = i
            data[i][14:94] = angles
            data[i][94:174] = I(angles, I0, D, phi)
    #        data[i][94:174] = I(angles, I0, D, phis[i])
            #data[i][94:174] = I(angles, I0, D, i*360/nb_rot)
            #data[i][94:174] = I(angles, I0, dolps[i], phi)
            data[i][174:254] = [I0] * 80

            data[i][92:94] = [0, 0]
            data[i][172:174] = [0, 0]
            data[i][252:254] = [0, 0]
        return data
