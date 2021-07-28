#!/usr/bin/python3
# -*-coding:utf-8 -*

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from utils import *
from scipy import signal
import sys as sys
import os
from copy import deepcopy
from subprocess import call
import observation
import geometry
import datetime as dt
import chaosmagpy as chaos

### Bottle classes.

class Bottle:
	def __init__(self, folder_name, auto=True, line=1, from_txt = False):
		self.folder = folder_name
		self.rotations = []

		self.AoBapp = False
		self.AoBlos = False
		self.AoRD = False
		self.time_stamp = ""

		self.from_txt = from_txt

		self.line = line

		self.valid = True

		try: #If you specify the exact input file name (with .in)
			self.input_parameters = self.ReadInputFile(pwd_data + self.folder)
		except:
			try: #If you only specify the folder, will look for a default input.in file
				self.input_parameters = self.ReadInputFile(pwd_data + self.folder + "/input.in")
			except: #Don't use it, or verify it works first.
				self.input_parameters = self.ReadInputFile(pwd_src + "input_files/" + self.folder + ".in")

		self.data_file_name = pwd_data + self.input_parameters["data_files"].split(",")[0]
		self.instrument_name = self.input_parameters["instrument_name"]
		if self.instrument_name == "ptcu": self.instrument_name = "carmen"
		elif self.instrument_name == "ptcu_v2" : self.instrument_name = "corbel"


		try:
			self.observation_type = self.input_parameters["observation_type"]
		except:
			self.observation_type = "fixed"
		print("INFO: Observation type is", self.observation_type)


		try:
			if self.input_parameters["I_ref"] == "off":
				self.NoVref = True
		except:
			self.NoVref = False
		try:
			self.saving_name = self.input_parameters["saving_name"]
			self.saving_name += "_voie" + str(self.line)
		except:
			print("WARNING: No saving name specified. Will be saved in " + self.data_file_name)
			self.saving_name = ""

		self.date, self.location, self.filters, self.azimut, self.elevation, self.com = "", "", ["0"], False, False, ""

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
		self.time_zone = dt.timedelta(hours = self.time_zone)


		### See PTCUBottle or SPPBottle for specific instructions

	def SetInfoFromDataFileName(self, f, data_f):
		self.folders, self.data_folders = f.split("/"), data_f.split("/")
		# i = len(data_folders) - 1
		self.folders_index = -3

		# self.instrument_name = self.folders[self.folders_index]

		### Follows in SPPBottle or PTCUBottle

	def MiseEnBouteille(self):
		### Loading the data from the files, then creating a list of Rotation objects depending on the instrument name specified. Can be 'spp', 'spp2014', 'fake_spp', 'ptcu', 'fake_ptcu'
		if not self.from_txt:
			self.LoadData()
			self.SetJumps()
			if self.valid:
				self.CleanRotations()

			### Now that we have all our rotations, creating lists of usefull data for easier smoothing
				self.CreateLists()
				self.GetSmoothLists()
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


	def SetJumps(self):
		self.jump_unit = self.input_parameters["jump_unit"]
		self.head_jump = float(self.input_parameters["head_jump"])
		self.tail_jump = float(self.input_parameters["tail_jump"])
		try:
			self.jump_mode = self.input_parameters["jump_mode"]
		except:
			self.jump_mode = "lenght"
		#Until the bug at the begining of the observation is not fixed, delete head and tail after treatment. When fixed, we'll be able to do it before
		if self.jump_unit in ("seconds", "minutes", "hours"):
			self.head_jump = dt.timedelta(0)
			self.tail_jump = dt.timedelta(0)


		if self.jump_mode == "time": # and self.jump_unit in ["seconds", "minutes", "hours"]:
			self.head_jump = dt.timedelta(seconds=float(self.input_parameters["head_jump"]))
			self.tail_jump = dt.timedelta(seconds=float(self.input_parameters["tail_jump"]))
			if self.jump_unit == "minutes":
				self.head_jump = dt.timedelta(minutes=float(self.input_parameters["head_jump"]))
				self.tail_jump = dt.timedelta(minutes=float(self.input_parameters["tail_jump"]))
			elif self.jump_unit == "hours":
				# print("DEBUG jumps", [r.time.total_seconds() for r in self.rotations[:10]])
				# print(time.timedelta(hours=float(self.input_parameters["head_jump"])), time.timedelta(hours=float(self.input_parameters["tail_jump"])))
				self.head_jump = dt.timedelta(hours=float(self.input_parameters["head_jump"]))
				self.tail_jump = dt.timedelta(hours=float(self.input_parameters["tail_jump"]))

			if self.jump_mode == "length" or self.tail_jump.seconds == 0.:
				try:
					self.tail_jump = self.all_times[-1] - self.tail_jump
				except:
					self.tail_jump = self.rotations[-1].time - self.tail_jump


	def LoadData(self):
		# Loading the data files

		### See SPPBottle/PTCUBottle for first part

		self.nb_rot = len(self.rotations)
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
		###At this point, all times are expressed in seconds since the first rotation before deleting the head_jump
		###If SPP before 2019 03 10, all times are expressed in seconds since midnight of the first observation day.

		# print("DEBUG TIME LISTS: ", self.time, self.head_jump, self.rotations[0].time)
		# self.all_times = np.array([r.time - self.head_jump for r in self.rotations])
		norm = self.rotations[0].time
		self.all_times = np.array([r.time - norm for r in self.rotations])
		# self.all_times_since_start = np.array([r.time  for r in self.rotations])
		# self.all_times_since_start = np.array([r.time - self.rotations[0].time for r in self.rotations])

		# self.all_datetimes = [time.gmtime(r.time + self.time + self.head_jump) for r in self.rotations]
		# self.all_datetimes = [time.gmtime(time.mktime(self.DateTime()) + r.time + self.time + self.head_jump) for r in self.rotations]

		# self.all_datetimes = [self.DateTime() + self.head_jump + t for t in self.all_times]
		# self.all_datetimes = [self.DateTime() + t for t in self.all_times_since_start]

		# print(self.head_jump, self.all_times[0], self.all_times[-1], self.all_datetimes[0], self.all_datetimes[-1])

		# for i in range(1, self.nb_rot):
		# 	if self.all_times[i] < self.all_times[i-1]:
		# 		self.all_times[i] += 24 * 3600

		# self.AoLP_correction = float(self.input_parameters["AoLP_correction"])*DtoR

		if not self.from_txt and self.instrument_name in ["corbel", "gdcu", "ptcu_v2"]:
			# if self.DateTime() < dt.datetime(year=2020, month=10, day=1):
			self.AoLP_correction = 0 # USE only for data observed and downloaded BEFORE October 2020!!!! For everything observed after that or doawnloaded via KVA20, you have to apply the correction below
			# else:
				# self.AoLP_correction = -(self.config['IS_PolarizerOffset' + str(self.line)]) * DtoR #This is the nor now (after oct 2020)

			# self.AoLP_correction = (self.config['IS_PolarizerOffset' + str(self.line)] + 45) * DtoR ## This was the formula apllied to convert all data from before oct 2020 in the database.


		elif not self.from_txt and self.instrument_name == "spp":
			self.AoLP_correction = float(self.input_parameters["AoLP_correction"]) * DtoR
		else:
			self.AoLP_correction = 0
		print("AoLP correction:", self.line, self.AoLP_correction * RtoD)


		self.all_V 	  = np.array([r.V    for r in self.rotations])
		self.all_Vcos = np.array([r.Vcos for r in self.rotations])
		self.all_Vsin = np.array([r.Vsin for r in self.rotations])
		try:
			self.all_Iref = np.array([r.Iref for r in self.rotations])
		except:
			self.NoVref = True
			print("WARNING: No reference canal")

		self.all_TempPM = np.array([r.TempPM for r in self.rotations])
		self.all_TempOptical = np.array([r.TempOptical for r in self.rotations])
		self.all_TempAmbiant = np.array([r.TempAmbiant for r in self.rotations])

		self.all_I0   = np.array([r.I0   for r in self.rotations])
		self.all_DoLP = np.array([r.DoLP for r in self.rotations])

		### AoLP_correction -133 en janvier 2019 -> angle de calib pas sur avant
		self.all_AoLP = np.array([r.AoLP for r in self.rotations]) + self.AoLP_correction
		self.all_AoLP = SetAngleBounds(self.all_AoLP, -np.pi/2, np.pi/2)

	# def CreateListsFromGivenData(self):
	# 	###At this point, all times are expressed in seconds since the first rotation before deleting the head_jump
	# 	###If SPP before 2019 03 10, all times are expressed in seconds since midnight of the first observation day.
	# 	print("DEBUG TIME LISTS: ", self.time, self.head_jump, self.rotations[0].time)
	# 	# self.all_times = np.array([r.time - self.head_jump for r in self.rotations])
	# 	self.all_times = np.array([r.time for r in self.rotations])
	# 	# self.all_times_since_start = np.array([r.time  for r in self.rotations])
	# 	self.all_times_since_start = np.array([r.time - self.rotations[0].time for r in self.rotations])
	#
	# 	# self.all_datetimes = [time.gmtime(r.time + self.time + self.head_jump) for r in self.rotations]
	# 	# self.all_datetimes = [time.gmtime(time.mktime(self.DateTime()) + r.time + self.time + self.head_jump) for r in self.rotations]
	#
	# 	self.all_datetimes = [self.DateTime() + t for t in self.all_times]
	# 	# self.all_datetimes = [self.DateTime() + t for t in self.all_times_since_start]
	#
	# 	print(self.head_jump, self.all_times[0], self.all_times[-1], self.all_datetimes[0], self.all_datetimes[-1])
	#
	# 	# for i in range(1, self.nb_rot):
	# 	# 	if self.all_times[i] < self.all_times[i-1]:
	# 	# 		self.all_times[i] += 24 * 3600
	#
	# 	self.AoLP_correction = float(self.input_parameters["AoLP_correction"])*DtoR
	# 	print(self.AoLP_correction, self.AoLP_correction*RtoD)
	#
	# 	self.all_V 	 = np.array([r.V    for r in self.rotations])
	# 	self.all_Vcos = np.array([r.Vcos for r in self.rotations])
	# 	self.all_Vsin = np.array([r.Vsin for r in self.rotations])
	# 	try:
	# 		self.all_Iref = np.array([r.Iref for r in self.rotations])
	# 	except:
	# 		self.NoVref = True
	# 		print("WARNING: No reference canal")
	#
	# 	self.all_TempPM = np.array([r.TempPM for r in self.rotations])
	# 	self.all_TempOptical = np.array([r.TempOptical for r in self.rotations])
	# 	self.all_TempAmbiant = np.array([r.TempAmbiant for r in self.rotations])
	#
	# 	self.all_I0 	 = np.array([r.I0   for r in self.rotations])
	# 	self.all_DoLP = np.array([r.DoLP for r in self.rotations])
	#
	# 	### AoLP_correction -133 en janvier 2019 -> angle de calib pas sur avant
	# 	self.all_AoLP = np.array([r.AoLP for r in self.rotations]) + self.AoLP_correction
	# 	self.all_AoLP = SetAngleBounds(self.all_AoLP, -np.pi/2, np.pi/2)

	def SetSmoothLists(self):
		print("Set Smooth Lists")
		### Smoothing procedure. Can smooth for a given number of rotations (smoothing_unit==rotations) or a  given time period (smoothing_unit==seconds).
		self.smoothing_factor = int(self.input_parameters["smoothing_factor"]) #Average the data over smoothing_factor rotations
		self.smoothing_unit = self.input_parameters["smoothing_unit"].lower()
		self.smoothing_method = self.input_parameters["smoothing_method"].lower()

		self.nb_smooth_rot = self.nb_rot

		self.avg_dt = 1000 * np.average([self.all_times[i].total_seconds() - self.all_times[i-1].total_seconds() for i in range(1, len(self.all_times))]) #in millisec

		print("AVG DT (millisec)", self.avg_dt)

	def GetSmoothLists(self):

		self.SetSmoothLists()

		print("Smooting data over {} {}".format(self.smoothing_factor, self.smoothing_unit))
		self.smooth_V    = GetSliddingAverage(self.all_V,    self.all_times, self.smoothing_factor, self.smoothing_unit)
		self.smooth_Vcos = GetSliddingAverage(self.all_Vcos, self.all_times, self.smoothing_factor, self.smoothing_unit)
		self.smooth_Vsin = GetSliddingAverage(self.all_Vsin, self.all_times, self.smoothing_factor, self.smoothing_unit)
		if not self.NoVref:
			self.smooth_Iref = GetSliddingAverage(self.all_Iref, self.all_times, self.smoothing_factor, self.smoothing_unit)
			self.Iref_average = np.average(self.all_Iref)

		self.nb_smooth_rot = len(self.smooth_V)


		### Calculate the smooth I0, DoLP and AoLP
		self.smooth_I0, self.smooth_DoLP, self.smooth_AoLP = np.zeros(self.nb_smooth_rot), np.zeros(self.nb_smooth_rot), np.zeros(self.nb_smooth_rot)
		for i in range(self.nb_smooth_rot):
			self.smooth_I0[i], self.smooth_DoLP[i], self.smooth_AoLP[i] = Rotation.GetLightParameters(self.smooth_V[i], self.smooth_Vcos[i], self.smooth_Vsin[i])

		self.smooth_AoLP = self.smooth_AoLP + self.AoLP_correction

		self.smooth_AoLP = SetAngleBounds(self.smooth_AoLP, -np.pi/2, np.pi/2)

		self.V_average = np.average(self.all_V)
		self.Vcos_average = np.average(self.all_Vcos)
		self.Vsin_average = np.average(self.all_Vsin)
		self.I0_average, self.DoLP_average, self.AoLP_average = Rotation.GetLightParameters(self.V_average, self.Vcos_average, self.Vsin_average)

		self.AoLP_average += self.AoLP_correction
		self.AoLP_average = SetAngleBounds(self.AoLP_average, -np.pi/2, np.pi/2)


		print("Computing Error bars...")
		#Get Variance
		if self.smoothing_unit == "seconds": smoothing_factor = self.smoothing_factor * 1000
		elif self.smoothing_unit == "minutes": smoothing_factor = self.smoothing_factor * 60 * 1000
		elif self.smoothing_unit == "hours": smoothing_factor = self.smoothing_factor * 3600 * 1000



		# self.std_I0 = 0		 #np.sqrt(4 * self.all_V 	/ 500)
		# self.std_smooth_I0 = 0# np.sqrt(4 * self.smooth_V  / smoothing_factor)
		#
		# self.std_DoLP = 	0#	np.sqrt(4 * (1 + ((self.all_DoLP/100) 	 ** 2 / 2)) / (self.all_I0 	  * 500)) * 100
		# self.std_smooth_DoLP =0# 	np.sqrt(4 * (1 + ((self.smooth_DoLP/100) ** 2 / 2)) / (self.smooth_I0 * smoothing_factor)) * 100
		#
		# self.std_AoLP = 	0#	np.sqrt(1 / ((self.all_DoLP/100) 	** 2 * self.all_I0 	  * 500))
		# self.std_smooth_AoLP = 0#	np.sqrt(1 / ((self.smooth_DoLP/100) ** 2 * self.smooth_I0 * smoothing_factor))
		# self.smooth_AoLP_upper = 0#self.smooth_AoLP + self.std_smooth_AoLP
		# self.smooth_AoLP_lower = 0#self.smooth_AoLP - self.std_smooth_AoLP
		self.std_I0 = 		 np.sqrt(2 * self.all_I0 	/ self.avg_dt)
		self.std_smooth_I0 = np.sqrt(2 * self.smooth_I0  / smoothing_factor)

		# print(self.all_V)
		# print(self.std_I0)
		# print(self.std_smooth_I0)

		self.std_DoLP = 		np.sqrt(4 * (1 + ((self.all_DoLP/100) 	 ** 2 / 2)) / (self.all_I0 	  * self.avg_dt)) * 100
		self.std_smooth_DoLP = 	np.sqrt(4 * (1 + ((self.smooth_DoLP/100) ** 2 / 2)) / (self.smooth_I0 * smoothing_factor)) * 100

		self.std_AoLP = 		np.sqrt(1 / ((self.all_DoLP/100) 	** 2 * self.all_I0 	  * self.avg_dt))
		self.std_smooth_AoLP = 	np.sqrt(1 / ((self.smooth_DoLP/100) ** 2 * self.smooth_I0 * smoothing_factor))
		self.smooth_AoLP_upper = self.smooth_AoLP + self.std_smooth_AoLP
		self.smooth_AoLP_lower = self.smooth_AoLP - self.std_smooth_AoLP
		# print(self.smooth_DoLP.shape, self.sm ooth_AoLP.shape, self.var_DoLP.shape, self.var_AoLP.shape)

		SN = lambda I, D, T: D * np.sqrt(I * T) / 2

		self.all_SN = SN(self.all_I0, self.all_DoLP/100, self.avg_dt)
		self.smooth_SN = SN(self.smooth_I0, self.smooth_DoLP/100, smoothing_factor)


		self.SetUnifyAngles()
		# self.graph_angle_shift = 0

		print("Get Smooth Lists: DONE")

	def SetUnifyAngles(self):
		print("Unifying AoLPs for least deviation")

		self.smooth_AoLP_upper = self.smooth_AoLP + self.std_smooth_AoLP
		self.smooth_AoLP_lower = self.smooth_AoLP - self.std_smooth_AoLP

		self.smooth_AoLP, smooth_graph_angle_shift = UnifyAngles(self.smooth_AoLP)
		self.all_AoLP, all_graph_angle_shift = UnifyAngles(self.all_AoLP)
		self.smooth_AoLP_upper, tmp = UnifyAngles(self.smooth_AoLP_upper)
		self.smooth_AoLP_lower, tmp = UnifyAngles(self.smooth_AoLP_lower)
		self.graph_angle_shift = max(all_graph_angle_shift, smooth_graph_angle_shift)

		if all_graph_angle_shift != self.graph_angle_shift:
			self.all_AoLP, all_graph_angle_shift = UnifyAngles(self.all_AoLP, manual_shift = self.graph_angle_shift)
			self.smooth_AoLP_upper, tmp = UnifyAngles(self.smooth_AoLP_upper, manual_shift = self.graph_angle_shift)
			self.smooth_AoLP_lower, tmp = UnifyAngles(self.smooth_AoLP_lower, manual_shift = self.graph_angle_shift)
		if smooth_graph_angle_shift != self.graph_angle_shift:
			self.smooth_AoLP, smooth_graph_angle_shift = UnifyAngles(self.smooth_AoLP, manual_shift = self.graph_angle_shift)
			self.smooth_AoLP_upper, tmp = UnifyAngles(self.smooth_AoLP_upper, manual_shift = self.graph_angle_shift)
			self.smooth_AoLP_lower, tmp = UnifyAngles(self.smooth_AoLP_lower, manual_shift = self.graph_angle_shift)

		if self.graph_angle_shift == 1:
			self.AoLP_average = SetAngleBounds(self.AoLP_average, 0, np.pi)
			self.smooth_AoLP_upper = SetAngleBounds(self.smooth_AoLP_upper, 0, np.pi)
			self.smooth_AoLP_lower = SetAngleBounds(self.smooth_AoLP_lower, 0, np.pi)
			# if self.AoLP_average < 0:
			# 	self.AoLP_average += np.pi
		elif self.graph_angle_shift == 0:
			self.AoLP_average = SetAngleBounds(self.AoLP_average, -np.pi/2, np.pi/2)
			self.smooth_AoLP_upper = SetAngleBounds(self.smooth_AoLP_upper, -np.pi/2, np.pi/2)
			self.smooth_AoLP_lower = SetAngleBounds(self.smooth_AoLP_lower, -np.pi/2, np.pi/2)
			# if self.AoLP_average > np.pi/2:
			# 	self.AoLP_average -= np.pi

		# for i, up in enumerate(self.smooth_AoLP_upper):
		# 	if up < self.smooth_AoLP[i]:

	def GetGeometryAngles(self, obs, B_model = None):
		# print("DEBUG obs", obs)

		obs.SinglePointGeometry(B_model = B_model)
		AoBapp = obs.eta_chaos
		AoBapp_ortho = AoBapp + np.pi / 2
		AoBlos = obs.Blos

		AoRD = obs.GetRayleighAngle(obs.RD_src_azimut, obs.RD_src_elevation)
		AoRD_ortho = AoRD + np.pi/2

		if self.graph_angle_shift == 1:
			if AoBapp < 0:
				AoBapp += np.pi
			if AoRD < 0:
				AoRD += np.pi
		if self.graph_angle_shift == 0 and AoBlos > np.pi/2:
			if AoRD_ortho > np.pi/2:
				AoRD_ortho -= np.pi
			if AoBlos > np.pi/2:
				AoBlos -= np.pi
			if AoBapp_ortho > np.pi/2:
				AoBapp_ortho -= np.pi

		return obs, AoBapp, AoBlos, AoRD, AoBapp_ortho, AoRD_ortho

	def Geometry(self, to_initiate = True):

		print("Geometry: Start")

		wd = os.getcwd()
		# os.chdir(pwd_src + "Geometry/Leo/src/")

		A_lon, A_lat = False, False
		A_lon, A_lat = geometry.GetLonLatFromName(self.location.lower())

		if self.filters[0] == "r":
			h = 220
		elif self.filters[0] == "v":
			h = 110
		else:
			h = 90

		try:
			self.source_azimut = float(self.input_parameters["pollution_source_azimut"]) * DtoR
			try:
				self.source_elevation = float(self.input_parameters["pollution_source_elevation"]) * DtoR
			except:
				self.source_elevation = 0
		except:
			self.source_azimut = self.source_elevation = 0

		B_model = chaos.load_CHAOS_matfile('/home/bossel/These/Analysis/data/magn_field/CHAOS-7.4.mat')

		if self.observation_type == "fixed" and self.azimut is not None and self.elevation is not None:
			# print("DEBUG SOURCE:", self.source_azimut*RtoD, self.source_elevation*RtoD)
			try:
				# print("DEBUG: AZ/EL", self.azimut*RtoD, self.elevation*RtoD)
				self.observation, self.AoBapp, self.AoBlos, self.AoRD, self.AoBapp_ortho, self.AoRD_ortho = self.GetGeometryAngles(observation.ObservationPoint(A_lon, A_lat, h, self.azimut, self.elevation, self.source_azimut, self.source_elevation, init_full = False), B_model = B_model)
			except:
				print("WARNING: No geometric angles where retrieved.")

		elif self.observation_type == "fixed_elevation_discrete_rotation":
			if to_initiate:
				self.discrete_rotation_elevation 	= float(self.input_parameters["discrete_rotation_elevation"]) * DtoR
				self.discrete_rotation_times 		= np.array([float(t) for t in self.input_parameters["discrete_rotation_times"].split("_")])
				self.discrete_rotation_azimuts 		= np.array([float(a)*DtoR for a in self.input_parameters["discrete_rotation_azimuts"].split("_")])

			ang_list = []
			for a in self.discrete_rotation_azimuts:
				ang_list.append(self.GetGeometryAngles(observation.ObservationPoint(A_lon, A_lat, h, a, self.discrete_rotation_elevation, self.source_azimut, self.source_elevation, init_full = False), B_model = B_model))

			self.observation, self.AoBapp, self.AoBlos, self.AoRD, self.AoBapp_ortho, self.AoRD_ortho = zip(*ang_list)
			self.observation, self.AoBapp, self.AoBlos, self.AoRD, self.AoBapp_ortho, self.AoRD_ortho = np.array(self.observation), np.array(self.AoBapp), np.array(self.AoBlos), np.array(self.AoRD), np.array(self.AoBapp_ortho), np.array(self.AoRD_ortho)

			self.x_axis = self.all_times

		elif self.observation_type == "fixed_azimut_discrete_rotation":
			if to_initiate:
				self.discrete_rotation_azimut 	= float(self.input_parameters["discrete_rotation_azimut"]) * DtoR
				self.discrete_rotation_times 		= [float(t) for t in self.input_parameters["discrete_rotation_times"].split("_")]
				self.discrete_rotation_times_unit = self.input_parameters["discrete_rotation_times_unit"]
				self.discrete_rotation_elevations 		= [float(a)*DtoR for a in self.input_parameters["discrete_rotation_elevations"].split("_")]

			ang_list = []
			for e in self.discrete_rotation_elevations:
				ang_list.append(self.GetGeometryAngles(observation.ObservationPoint(A_lon, A_lat, h, self.discrete_rotation_azimut, e, self.source_azimut, self.source_elevation, init_full = False), B_model = B_model))

			self.observation, self.AoBapp, self.AoBlos, self.AoRD, self.AoBapp_ortho, self.AoRD_ortho = zip(*ang_list)
			self.observation, self.AoBapp, self.AoBlos, self.AoRD, self.AoBapp_ortho, self.AoRD_ortho = np.array(self.observation), np.array(self.AoBapp), np.array(self.AoBlos), np.array(self.AoRD), np.array(self.AoBapp_ortho), np.array(self.AoRD_ortho)


		elif self.observation_type == "fixed_elevation_continue_rotation":
			if to_initiate:
				self.continue_rotation_elevation 	= float(self.input_parameters["continue_rotation_elevation"]) * DtoR
				self.continue_rotation_times 		= [float(t) for t in self.input_parameters["continue_rotation_times"].split("_")]
				if self.continue_rotation_times[0] != 0: self.continue_rotation_times.insert(0, 0)
				self.continue_rotation_times.append(self.all_times[-1].total_seconds() / 60.)
				self.continue_rotation_times = np.array(self.continue_rotation_times)
				self.nb_continue_rotation			= len(self.continue_rotation_times) - 1
				# self.rotation_direction = self.input_parameters["rotation_direction"]
			print("self.continue_rotation_times", self.continue_rotation_times)

			geo = geometry.Geometry("dummy", str(self.location), str(h), "e", str(self.continue_rotation_elevation*RtoD), str(self.source_azimut*RtoD), str(self.source_elevation*RtoD))
			try:
				self.azimut, self.observation = geo.FixedElevation(direction = self.input_parameters["rotation_direction"], B_model = B_model)
			except:
				self.azimut, self.observation = geo.FixedElevation(B_model = B_model)


			ang_list = [self.GetGeometryAngles(o, B_model = B_model) for o in self.observation]
			self.observation, self.AoBapp, self.AoBlos, self.AoRD, self.AoBapp_ortho, self.AoRD_ortho = zip(*ang_list)
			self.observation, self.AoBapp, self.AoBlos, self.AoRD, self.AoBapp_ortho, self.AoRD_ortho = np.array(self.observation), np.array(self.AoBapp), np.array(self.AoBlos), np.array(self.AoRD), np.array(self.AoBapp_ortho), np.array(self.AoRD_ortho)

		os.chdir(wd)
		# print("Geometry:", self.AoBapp, self.AoBlos, self.AoRD)
		print("Geometry: DONE")


	def GetOneTimePeriod(self, start_time, end_time, smooth_data = True, direct = False):
		"""Return the data points between the start and end times given."""
		start_time = dt.timedelta(minutes = start_time)
		end_time = dt.timedelta(minutes = end_time)

		if direct == True:
			I0 = self.all_I0[(start_time < self.all_times) &  (self.all_times < end_time)]
			DoLP = self.all_DoLP[(start_time < self.all_times) &  (self.all_times < end_time)]
			AoLP = self.all_AoLP[(start_time < self.all_times) &  (self.all_times < end_time)]
			nb_times = len(I0)
			return I0, DoLP, AoLP, nb_times
		else:
			V = self.all_V[(start_time < self.all_times) &  (self.all_times < end_time)]
			Vcos = self.all_Vcos[(start_time < self.all_times) &  (self.all_times < end_time)]
			Vsin = self.all_Vsin[(start_time < self.all_times) &  (self.all_times < end_time)]
			nb_times = len(V)
			return V, Vcos, Vsin, nb_times



	def GetAverageContinueRotations(self):
		"""Returns a new bottle where all V, Vcos, Vsin are the average of the different (azimutal) rotations of this one. Is used for the mixer if the observation_type == fixed_elevation_continue_rotation."""

		#Get the smallest rotation (least number of points)
		min_nb_times = 10**100
		i = 0
		for i in range(0, self.nb_continue_rotation):
			start, end = self.continue_rotation_times[i], self.continue_rotation_times[i + 1]
			tmp_min_nb_times = self.GetOneTimePeriod(start, end, smooth_data = False)[3]
			if tmp_min_nb_times < min_nb_times:
				min_nb_times = tmp_min_nb_times
				i_min = i

		avg_V, avg_Vcos, avg_Vsin, nb_times = self.GetOneTimePeriod(self.continue_rotation_times[i_min], self.continue_rotation_times[i_min+1], smooth_data = False)

		for i in range(1, self.nb_continue_rotation):
			start, end = self.continue_rotation_times[i], self.continue_rotation_times[i + 1]

			tmp_V, tmp_Vcos, tmp_Vsin, tmp_nb_times  = self.GetOneTimePeriod(start, end, smooth_data = False)
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

		new_bottle  = deepcopy(self)

		new_bottle.continue_rotation_times = np.array([self.continue_rotation_times[0], self.continue_rotation_times[1]])
		new_bottle.nb_continue_rotation = 1
		new_bottle.all_V = avg_V
		new_bottle.all_Vcos = avg_Vcos
		new_bottle.all_Vsin = avg_Vsin
		new_bottle.all_I0 = I0
		new_bottle.all_DoLP = DoLP
		new_bottle.all_AoLP = AoLP
		new_bottle.all_times = np.arange(0)
		new_bottle.saving_name += "_rotations_average"

		for t in np.linspace(self.continue_rotation_times[0], self.continue_rotation_times[1], len(new_bottle.all_V)):
			new_bottle.all_times = np.append(new_bottle.all_times, dt.timedelta(minutes = t))

		new_bottle.NoVref = True
		new_bottle.GetSmoothLists()
		new_bottle.Geometry(to_initiate=False)

		return new_bottle

	def PrintInfo(self):
		print("*********************************************")
		print("INFO:", self.data_file_name.split("/")[-3:])
		print("DoLP min:", np.min(self.smooth_DoLP))
		print("DoLP max:", np.max(self.smooth_DoLP))
		print("DoLP stddev:", np.std(self.smooth_DoLP))
		print("AoLP stddev:", np.std(self.smooth_AoLP) * RtoD)
		if self.observation_type == "fixed":
			print("AoBlos:", self.AoBlos * RtoD)
		print("*********************************************")

	def SaveTXT(self):
		print("Saving as .txt in", self.data_file_name + "/" + self.saving_name + '_results.txt')

	def ConvertRotTimeTogmTime(self, rot_time):
		t = time.mktime(self.DateTime())
		tl = [time.gmtime(t + rt) for rt in rot_time]
		return tl

	def GetRotationFrequency(self):
		all_freq = [1. / (self.all_times[i].total_seconds() - self.all_times[i-1].total_seconds()) for i in range(1, self.nb_rot)]
		return np.average(all_freq)

	def RenormTimes(self, shift):
		"""Shift all times by a certain amount"""
		self.all_times = np.array([t + shift for t in self.all_times])


	def DateTime(self, moment="start", delta=dt.timedelta(seconds=0), format="LT"):
		norm = dt.timedelta(seconds=0)
		if format == "LT" and self.config_time_format == "UT":
			norm = self.time_zone
		elif format == "UT" and self.config_time_format == "LT":
			norm = - self.time_zone

		if moment == "start": #Datetime of first good rotation
			return self.all_times[0] + self.datetime + self.head_jump + delta + norm
			# return self.rotations[0].time + self.datetime + self.head_jump + delta
		elif moment == "end": #Datetime of last good rotation
			return self.all_times[-1] + self.datetime + self.head_jump + delta + norm
			# return self.rotations[-1].time + self.datetime + self.head_jump + delta
		elif moment == "config": #Datetime written in the config file (== first bad rotation, before head_jump)
			return self.datetime + delta + norm
		elif moment == "midnight": #Datetime written in the config file (== first bad rotation, before head_jump)
			return dt.datetime(year=self.datetime.year, month=self.datetime.month, day=self.datetime.day) + delta + norm

	def GetInterpolation(self, time):
		"""Return the linear interpolation of smooth data for the given times in seconds. Time must be in seconds since the start of this bottle."""
		all_times = list(map(lambda x: x.total_seconds(), self.all_times))

		interp_I0 = np.interp(time, all_times, self.smooth_I0)#, left=0, right=0, period=None)
		interp_DoLP = np.interp(time, all_times, self.smooth_DoLP)#, left=0, right=0, period=None)
		interp_AoLP = np.interp(time, all_times, self.smooth_AoLP)#, period=np.pi)

		return interp_I0, interp_DoLP, interp_AoLP


	def __add__(self, bottle_to_add):

		print("Adding:")


		norm = bottle_to_add.DateTime("config") - self.DateTime("config")
		for ir, r in enumerate(bottle_to_add.rotations):
			bottle_to_add.rotations[ir].time += norm
			# bottle_to_add.rotations[ir].time += self.rotations[-1].time

		if len(bottle_to_add.rotations) > 0: #If not read from a file

			# print(len(self.rotations), len(bottle_to_add.rotations))
			self.rotations = np.append(self.rotations, bottle_to_add.rotations)

			# print(len(self.rotations), len(bottle_to_add.rotations))

			self.tail_jump = self.rotations[-1].time
			if not self.from_txt:
				self.CleanRotations()

			### Now that we have all our rotations, creating lists of usefull data for easier smoothing
				self.CreateLists()
				self.GetSmoothLists()

		else: # If read from a file (-f)
			bottle_to_add.all_times += bottle_to_add.DateTime("start") - self.DateTime("start")
			# bottle_to_add.all_times += self.all_times[-1]

			self.all_times = np.append(self.all_times, bottle_to_add.all_times)
			self.nb_rot = len(self.all_times)

			self.smooth_V    = np.append(self.smooth_V, bottle_to_add.smooth_V)
			self.smooth_Vcos = np.append(self.smooth_Vcos, bottle_to_add.smooth_Vcos)
			self.smooth_Vsin = np.append(self.smooth_Vsin, bottle_to_add.smooth_Vsin)

			self.nb_smooth_rot = len(self.smooth_V)

			### Calculate the smooth I0, DoLP and AoLP
			self.smooth_I0 = np.append(self.smooth_I0, bottle_to_add.smooth_I0)
			self.smooth_DoLP = np.append(self.smooth_DoLP, bottle_to_add.smooth_DoLP)
			self.smooth_AoLP = np.append(self.smooth_AoLP, bottle_to_add.smooth_AoLP)

			self.all_V =  np.append(self.all_V, bottle_to_add.all_V)
			self.all_Vcos =  np.append(self.all_Vcos, bottle_to_add.all_Vcos)
			self.all_Vsin =  np.append(self.all_Vsin, bottle_to_add.all_Vsin)
			self.all_I0 =  np.append(self.all_I0, bottle_to_add.all_I0)
			self.all_DoLP =  np.append(self.all_DoLP, bottle_to_add.all_DoLP)
			self.all_AoLP =  np.append(self.all_AoLP, bottle_to_add.all_AoLP)

			self.V_average = np.average(self.all_V)
			self.Vcos_average = np.average(self.all_Vcos)
			self.Vsin_average = np.average(self.all_Vsin)
			self.I0_average, self.DoLP_average, self.AoLP_average = Rotation.GetLightParameters(self.V_average, self.Vcos_average, self.Vsin_average)

			self.std_I0 			= np.append(self.std_I0, bottle_to_add.std_I0)
			self.std_smooth_I0 		= np.append(self.std_smooth_I0, bottle_to_add.std_smooth_I0)
			self.std_DoLP 			= np.append(self.std_DoLP, bottle_to_add.std_DoLP)
			self.std_smooth_DoLP 	= np.append(self.std_smooth_DoLP, bottle_to_add.std_smooth_DoLP)
			self.std_AoLP 			= np.append(self.std_AoLP, bottle_to_add.std_AoLP)
			self.std_smooth_AoLP 	= np.append(self.std_smooth_AoLP, bottle_to_add.std_smooth_AoLP)

			self.SetUnifyAngles()


			self.tail_jump = self.all_times[-1]

		self.graph_angle_shift = max(self.graph_angle_shift, bottle_to_add.graph_angle_shift)
		self.Geometry()

		return self


#####################################################################################
###									Petit Cru									  ###
#####################################################################################


class PTCUBottle(Bottle):
	def __init__(self, folder_name, auto=True, line = 1, from_txt = False):
		Bottle.__init__(self, folder_name, auto, line = line, from_txt = from_txt)

		self.SetInfoFromDataFileName(self.data_file_name, pwd_data)

		if auto:
			self.MiseEnBouteille()


	def SetInfoFromDataFileName(self, f, data_f):
		Bottle.SetInfoFromDataFileName(self, f, data_f)

		i = self.folders_index
		date_location = self.folders[i+1].split("_")
		self.date, self.location = date_location[0], date_location[1]
		rest = self.folders[i+2].split("_")
		# print(rest)
		#Filters: r,v,b,m pour rouge, vert, bleu, mauve
		#0: no filters_list
		#o: orange 620
		#X, Y: 620nm et 413nm
		filters_list = ["r", "v", "b", "m", "0", "o", "X", "Y", "t"]
		nb_filters = 2
		if self.instrument_name == "gdcu":
			nb_filters = 4

		for r in rest:
			print("R", r)
			if len(r) == nb_filters and np.all([a in filters_list for a in r]):
				if self.instrument_name == "carmen" or self.instrument_name == "fake_ptcu":
					self.filters = r
					if self.filters[1] in [0, "0"]:
						self.NoVref = True
				elif self.instrument_name in ["corbel", "gdcu"]:
					self.filters = r[self.line - 1]
				# print("FILTERS:", r, self.filters)

			elif r[0] == "a":
				self.azimut = float(r[1:]) * DtoR
			elif r[0] == "e":
				self.elevation = float(r[1:]) * DtoR
			else:
				self.com += "_" + r

		print("SetInfoFromDataFileName", self.instrument_name, self.date, self.location, self.filters, self.azimut*RtoD, self.elevation*RtoD, self.com)

	def LoadData(self):
		if self.instrument_name == "carmen" or self.instrument_name == "corbel" or self.instrument_name == "gdcu":
			self.raw_data, self.config = self.LoadPTCUData()
			self.nb_rot = len(self.raw_data)
			if self.nb_rot == 0:
				self.valid = False
			# self.nb_data_per_rot = len(self.raw_data[0])
			# if self.jump_mode == "length" or self.tail_jump == 0:
			# 	self.tail_jump = len(self.raw_data) - self.tail_jump

			for r in self.raw_data[:]:
				self.rotations.append(PTCURotation(r))

		elif self.instrument_name == "fake_ptcu":
			self.raw_data, self.config = LoadPTCUTestData(nb_rot=1000, I0=100, D=0.10, phi=45*RtoD)
			# if self.jump_mode == "length" or self.tail_jump == 0:
			# 	self.tail_jump = len(self.raw_data) - self.tail_jump
			for r in self.raw_data[:]:
				self.rotations.append(PTCURotation(r))

		if self.valid:
			Bottle.LoadData(self)

	def GetTimeFromDateTime(self, date):
		"""Return list self.all_times shifted so that 0 is date."""
		shift = self.DateTime(moment="start") - date
		return self.all_times + shift

	def SetTimeFromDateTime(self, date):
		self.date_time = date
		self.all_times = self.GetTimeFromDateTime(date)

	def SetTimes(self):
		try:
			self.time_stamp = str(self.config["Timestamp"]).split("\"")[1] #When importing config file from mysql workbench, timestamp is between "", but when using python script, there is no more ""
		except:
			self.time_stamp = str(self.config["Timestamp"]).split("\'")[1]

		print(self.time_stamp)
		#Date and time of first rotation (before deleting head_jump)

		self.datetime = time.datetime.strptime(self.time_stamp, "%Y-%m-%d %H:%M:%S")
		print(self.datetime, self.datetime.strftime("%H:%M:%S"))
		#Time in sec since EPOCH
		# self.time = self.datetime.timestamp()


	# def DateTime(self, moment="start", delta=None, format="LT"):
	# 	if not delta:
	# 		delta = dt.timedelta(seconds=0)
	#
	# 	norm = dt.timedelta(seconds=0)
	# 	if format == "LT" and self.config_time_format == "UT":
	# 		norm = self.time_zone
	# 	elif format == "UT" and self.config_time_format == "LT":
	# 		norm = - self.time_zone
	#
	# 	if moment == "start": #Datetime of first good rotation
	# 		return self.all_times[0] + self.datetime + self.head_jump + delta + norm
	# 		# return self.rotations[0].time + self.datetime + self.head_jump + delta
	# 	elif moment == "end": #Datetime of last good rotation
	# 		return self.all_times[-1] + self.datetime + self.head_jump + delta + norm
	# 		# return self.rotations[-1].time + self.datetime + self.head_jump + delta
	# 	elif moment == "config": #Datetime written in the config file (== first bad rotation, before head_jump)
	# 		return self.datetime + delta + norm

	def CleanRotations(self):
		self.nb_rot = len(self.rotations)

		### Put every rotation time in sec since first rotation
		norm = self.rotations[0].time
		for i in range(len(self.rotations)):
			self.rotations[i].time -= norm

		### Add 24h when next rotation is earlier than the previous one. Happens sometimes when we change the date
		for i in range(1, len(self.rotations)):
			while self.rotations[i].time < self.rotations[i-1].time:
				self.rotations[i].time += dt.timedelta(day=1)

		### Get and set the head and tail jumps
		# if self.jump_mode == "time": # and self.jump_unit in ["seconds", "minutes", "hours"]:
		# 	self.head_jump = time.timedelta(seconds=float(self.input_parameters["head_jump"]))
		# 	self.tail_jump = time.timedelta(seconds=float(self.input_parameters["tail_jump"]))
		#
		# 	if self.jump_unit == "minutes":
		# 		self.head_jump = time.timedelta(minutes=float(self.input_parameters["head_jump"]))
		# 		self.tail_jump = time.timedelta(minutes=float(self.input_parameters["tail_jump"]))
		# 	elif self.jump_unit == "hours":
		# 		print("DEBUG jumps", [r.time.total_seconds() for r in self.rotations[:10]])
		# 		print(time.timedelta(hours=float(self.input_parameters["head_jump"])), time.timedelta(hours=float(self.input_parameters["tail_jump"])))
		# 		self.head_jump = time.timedelta(hours=float(self.input_parameters["head_jump"]))
		# 		self.tail_jump = time.timedelta(hours=float(self.input_parameters["tail_jump"]))
		#
		# 	if self.jump_mode == "length" or self.tail_jump.seconds == 0.:
		# 		self.tail_jump = self.rotations[-1].time - self.tail_jump

			### Delete all rotations before the head and after the tail jump
		self.rotations = [r for r in self.rotations if self.head_jump <= r.time < self.tail_jump]
		print(len(self.rotations), "good rotations in", self.nb_rot, ";", self.nb_rot - len(self.rotations), "deleted because of time jumps.")
		self.nb_rot = len(self.rotations)

		# print(self.head_jump, self.tail_jump, self.rotations[0].time, self.rotations[-1].time)


		### Put every rotation time in sec since first rotation (after deleting head_jump)
		norm = self.rotations[0].time
		for i in range(len(self.rotations)):
			self.rotations[i].time -= norm

		# time = np.array([r.time.total_seconds()    for r in self.rotations])
		# all_V 	  = np.array([r.V    for r in self.rotations])
		# dV_dt = np.gradient(all_V, time )
		# avg_V = np.average(all_V)
		# print(avg_V)
		# self.rotations = [r for ir, r in enumerate(self.rotations) if abs(dV_dt[ir]) < 0.01 * avg_V]
		#
		# plt.plot(time, dV_dt)
		# plt.show()

		self.rotations = [r for r in self.rotations if r.IsGood(Imin=self.Imin, Imax=self.Imax, DoLPmin=self.DoLPmin, DoLPmax=self.DoLPmax)]
		print(len(self.rotations), "good rotations in", self.nb_rot, ";", self.nb_rot - len(self.rotations), "deleted because of invalid data.")
		self.nb_rot = len(self.rotations)

	def SaveTXT(self):
		print("Saving as .txt in", self.data_file_name + "/" + self.saving_name + '_results.txt')

		### Default Format
		times = [t.total_seconds() * 1000 for t in self.GetTimeFromDateTime(self.DateTime(moment="config"))]
		data = pd.DataFrame(np.array([times, self.all_V, self.all_Vcos, self.all_Vsin, self.all_I0, self.all_DoLP, self.all_AoLP, self.smooth_V, self.smooth_Vcos, self.smooth_Vsin, self.smooth_I0, self.smooth_DoLP, self.smooth_AoLP, self.std_I0, self.std_DoLP, self.std_AoLP, self.std_smooth_I0, self.std_smooth_DoLP, self.std_smooth_AoLP]).transpose())
		data.columns = ["time","V","Vcos","Vsin","I0","DoLP","AoLP","SV","SVcos","SVsin","SI0","SDoLP","SAoLP","errI0","errDoLP","errAoLP","errSI0","errSDoLP","errSAoLP"]

		### Simple format for Jean
		# times = [t.total_seconds() / 3600. for t in self.GetTimeFromDateTime(self.DateTime(moment="midnight", format="LT"))]
		# data = pd.DataFrame(np.array([times, self.smooth_I0, self.smooth_DoLP, self.smooth_AoLP, self.std_smooth_I0, self.std_smooth_DoLP, self.std_smooth_AoLP]).transpose())
		# data.columns = ["time","SI0","SDoLP","SAoLP","errSI0","errSDoLP","errSAoLP"]

		data.to_csv(self.data_file_name + "/" + self.saving_name + '_results.txt', sep="\t", index=False)

	def LoadPTCUData(self):
		"""Return an array of data and a dictionnary of config from a list of data files. The first is the data, the second is the config. Each line of the data is a rotation with 6 entries, the config dict is all the parameters of the observation, common to all rotations."""
		###DATA FILE
		# 0: 	IDProcessedData,
		# 1: 	IDConfiguration,
		# 2: 	Time since begining of obseravtion (timestamp) in milliseconds,
		# 3:	PMAvg,
		# 4:	PMCosAvg,
		# 5:	PMSinAvg
		# 6:	IDTemperatures
		# 7:	TempPM
		# 8:	TempOptical
		# 9:	TempAmbiant
		# 10:	IDLiveComment
		# 11:	Comment
		# 12:	live_commentscol
		###CONFIG FILE
		# 0:		IDConfiguration,
		# 1:		Timestamp,
		# 2-3:		CM_Latitude, CM_Longitude,
		# 4-5:		CM_Elevation, CM_Azimuth,
		# 6:		CM_PolarizerSpeed,
		# 7:		CM_Tilt,
		# 8:		CM_Usage,
		# 9:		CM_Comments,
		# 10:		IS_PolarizerOffset1,
		# 11:		IS_PolarizerOffset2,
		# 12:		IS_PolarizerOffset3,
		# 13:		IS_PolarizerOffset4,
		# 14:		IS_ConverterOffset4,
		# 15:		IS_ConverterOffset3,
		# 16:		IS_ConverterOffset2,
		# 17:		IS_ConverterOffset1,
		# 18:		IS_EngineTicksPerTour,
		# 19:		IS_MotorGearedRatio,
		# 20:		IS_QuantumEfficiency,
		# 21:		IS_IpAddress

		###For old (< V3) csv files
		file_names = self.data_file_name
		line = self.line

		print("LOADING PTCU data: line", line)
		data_path = file_names
		#Loading the data from a binary file.
		#Output is one big 1-D array of length nb_data_tot = nb_toy * nb_data_per_rot
		data_file = data_path + "/data" + str(line) + ".csv"
		config_file = data_path + "/config.csv"


		print(data_file, config_file)

		try:
			with open(data_file, "r") as f:
				first_line = f.readlines()[1]
			d = GetDelimiter(first_line)
		except:
			d = ""

		raw_data = np.genfromtxt(data_file, delimiter = d, skip_header=1)


		if self.instrument_name in ["corbel", "ptcu_v2"]:
			array_type = [('IDConfiguration',float),('Timestamp','S100'), ('CM_ID', float),('CM_Latitude',float),('CM_Longitude',float),('CM_Elevation',float),('CM_Azimuth',float),('CM_PolarizerSpeed',float),('CM_Tilt',float),('CM_Usage','S100'),('CM_Comments','S100'),('IS_PolarizerOffset1',float),('IS_PolarizerOffset2',float),('IS_PolarizerOffset3',float),('IS_PolarizerOffset4',float),('IS_ConverterGain1',float),('IS_ConverterGain2',float),('IS_ConverterGain3',float),('IS_ConverterGain4',float),("IS_ConverterOffset1", float),("IS_ConverterOffset2", float),("IS_ConverterOffset3", float),("IS_ConverterOffset4", float),('IS_MotorGearedRatio',float),('IS_QuantumEfficiency',float), ('IS_IpAddress',float)]
		else:
			array_type = [('IDConfiguration',float),('Timestamp','S100'),('CM_Latitude',float),('CM_Longitude',float),('CM_Elevation',float),('CM_Azimuth',float),('CM_PolarizerSpeed',float),('CM_Tilt',float),('CM_Usage','S100'),('CM_Comments','S100'),('IS_PolarizerOffset1',float),('IS_PolarizerOffset2',float),('IS_PolarizerOffset3',float),('IS_PolarizerOffset4',float),("IS_ConverterOffset4", float),("IS_ConverterOffset3", float),("IS_ConverterOffset2", float),("IS_ConverterOffset1", float),('IS_EngineTicksPerTour',float),('IS_MotorGearedRatio',float),('IS_QuantumEfficiency',float), ('IS_IpAddress',float)]

		print("DEBUG", config_file, len(array_type))
		configuration = np.genfromtxt(config_file, dtype=array_type, delimiter=",", skip_header=1)


		### For new hdf5 files
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

		self.all_times = np.array([dt.timedelta(milliseconds=t) for t in data["Time"]])
		self.nb_rot = len(self.all_times)

		self.SetJumps()

		self.smooth_V    = np.array([d for d, t in zip(data["smooth_V"], self.all_times) if self.head_jump <= t < self.tail_jump])
		self.smooth_Vcos = np.array([d for d, t in zip(data["smooth_Vcos"], self.all_times) if self.head_jump <= t < self.tail_jump])
		self.smooth_Vsin = np.array([d for d, t in zip(data["smooth_Vsin"], self.all_times) if self.head_jump <= t < self.tail_jump])

		self.nb_smooth_rot = len(self.smooth_V)

		### Calculate the smooth I0, DoLP and AoLP
		self.smooth_I0 = np.array([d for d, t in zip(data["Smooth_I0"], self.all_times) if self.head_jump <= t < self.tail_jump])
		self.smooth_DoLP = np.array([d for d, t in zip(data["Smooth_DoLP"], self.all_times) if self.head_jump <= t < self.tail_jump])
		self.smooth_AoLP = np.array([d for d, t in zip(data["Smooth_AoLP"], self.all_times) if self.head_jump <= t < self.tail_jump])

		self.all_V =  np.array([d for d, t in zip(data["PMAvg"], self.all_times) if self.head_jump <= t < self.tail_jump])
		self.all_Vcos =  np.array([d for d, t in zip(data["PMCosAvg"], self.all_times) if self.head_jump <= t < self.tail_jump])
		self.all_Vsin =  np.array([d for d, t in zip(data["PMSinAvg"], self.all_times) if self.head_jump <= t < self.tail_jump])
		self.all_I0 =  np.array([d for d, t in zip(data["I0"], self.all_times) if self.head_jump <= t < self.tail_jump])
		self.all_DoLP =  np.array([d for d, t in zip(data["DoLP"], self.all_times) if self.head_jump <= t < self.tail_jump])
		self.all_AoLP =  np.array([d for d, t in zip(data["AoLP"], self.all_times) if self.head_jump <= t < self.tail_jump])

		self.std_I0 			= np.array([d for d, t in zip(data["errI0"], self.all_times) if self.head_jump <= t < self.tail_jump])
		self.std_smooth_I0 		= np.array([d for d, t in zip(data["errSI0"], self.all_times) if self.head_jump <= t < self.tail_jump])

		self.std_DoLP 			= np.array([d for d, t in zip(data["errDoLP"], self.all_times) if self.head_jump <= t < self.tail_jump])
		self.std_smooth_DoLP 	= np.array([d for d, t in zip(data["errSDoLP"], self.all_times) if self.head_jump <= t < self.tail_jump])

		self.std_AoLP 			= np.array([d for d, t in zip(data["errAoLP"], self.all_times) if self.head_jump <= t < self.tail_jump])
		self.std_smooth_AoLP 	= np.array([d for d, t in zip(data["errSAoLP"], self.all_times) if self.head_jump <= t < self.tail_jump])

		norm = self.head_jump
		self.all_times = np.array([t - norm for t in self.all_times if self.head_jump <= t < self.tail_jump])

		self.valid = True
		_, self.config = self.LoadPTCUData()

		self.SetTimes()
		self.SetSmoothLists()

		self.SetTimeFromDateTime(self.DateTime(moment="start"))

		self.graph_angle_shift = 0
		for a in self.smooth_AoLP:
			if a > np.pi/2.:
				self.graph_angle_shift = 1
				break
			if a < 0:
				break

	def LoadFromTxt(self):
		file_name = self.data_file_name + "/" + self.saving_name + '_results.txt'

		data = pd.read_csv(file_name, delimiter="\t")
		data.columns = [c.replace("#", "").replace("(", "").replace(")", "").strip() for c in data.columns]
		print(data.columns)
		time_name = data.columns[0]


		self.all_times = np.array([dt.timedelta(milliseconds=t) for t in data[time_name]])
		self.nb_rot = len(self.all_times)

		self.SetJumps()

		self.smooth_V    = np.array([d for d, t in zip(data["SV"], self.all_times) if self.head_jump <= t < self.tail_jump])
		self.smooth_Vcos = np.array([d for d, t in zip(data["SVcos"], self.all_times) if self.head_jump <= t < self.tail_jump])
		self.smooth_Vsin = np.array([d for d, t in zip(data["SVsin"], self.all_times) if self.head_jump <= t < self.tail_jump])

		self.nb_smooth_rot = len(self.smooth_V)

		### Calculate the smooth I0, DoLP and AoLP
		self.smooth_I0 = np.array([d for d, t in zip(data["SI0"], self.all_times) if self.head_jump <= t < self.tail_jump])
		self.smooth_DoLP = np.array([d for d, t in zip(data["SDoLP"], self.all_times) if self.head_jump <= t < self.tail_jump])
		self.smooth_AoLP = np.array([d for d, t in zip(data["SAoLP"], self.all_times) if self.head_jump <= t < self.tail_jump])

		self.all_V =  np.array([d for d, t in zip(data["V"], self.all_times) if self.head_jump <= t < self.tail_jump])
		self.all_Vcos =  np.array([d for d, t in zip(data["Vcos"], self.all_times) if self.head_jump <= t < self.tail_jump])
		self.all_Vsin =  np.array([d for d, t in zip(data["Vsin"], self.all_times) if self.head_jump <= t < self.tail_jump])
		self.all_I0 =  np.array([d for d, t in zip(data["I0"], self.all_times) if self.head_jump <= t < self.tail_jump])
		self.all_DoLP =  np.array([d for d, t in zip(data["DoLP"], self.all_times) if self.head_jump <= t < self.tail_jump])
		self.all_AoLP =  np.array([d for d, t in zip(data["AoLP"], self.all_times) if self.head_jump <= t < self.tail_jump])

		self.std_I0 			= np.array([d for d, t in zip(data["errI0"], self.all_times) if self.head_jump <= t < self.tail_jump])
		self.std_smooth_I0 		= np.array([d for d, t in zip(data["errSI0"], self.all_times) if self.head_jump <= t < self.tail_jump])

		self.std_DoLP 			= np.array([d for d, t in zip(data["errDoLP"], self.all_times) if self.head_jump <= t < self.tail_jump])
		self.std_smooth_DoLP 	= np.array([d for d, t in zip(data["errSDoLP"], self.all_times) if self.head_jump <= t < self.tail_jump])

		self.std_AoLP 			= np.array([d for d, t in zip(data["errAoLP"], self.all_times) if self.head_jump <= t < self.tail_jump])
		self.std_smooth_AoLP 	= np.array([d for d, t in zip(data["errSAoLP"], self.all_times) if self.head_jump <= t < self.tail_jump])

		norm = self.head_jump
		self.all_times = np.array([t - norm for t in self.all_times if self.head_jump <= t < self.tail_jump])


		self.valid = True
		_, self.config = self.LoadPTCUData()


		self.SetTimes()
		self.SetSmoothLists()

		self.SetTimeFromDateTime(self.DateTime(moment="start"))

		self.graph_angle_shift = 0
		for a in self.smooth_AoLP:
			if a > np.pi/2.:
				self.graph_angle_shift = 1
				break
			if a < 0:
				break

		# if self.instrument_name in ["corbel", "gdcu"]:
		# 	self.all_AoLP 		-= 40 * DtoR
		# 	self.smooth_AoLP 	-= 40 * DtoR






#####################################################################################
###									SPPBottle									  ###
#####################################################################################


class SPPBottle(Bottle):
	def __init__(self, folder_name, auto=True, from_txt = False):
		Bottle.__init__(self, folder_name, auto, from_txt=from_txt)
		try:
			self.SetInfoFromDataFileName(self.data_file_name, pwd_data)
		except:
			self.location = ""
			print("WARNING: Folders not standard. Could not retrieve info from folders name.")

		self.filters = ["r"]

		if auto:
			self.MiseEnBouteille()


	def SetInfoFromDataFileName(self, f, data_f):
		Bottle.SetInfoFromDataFileName(self, f, data_f)
		i = self.folders_index
		if len(self.folders[i+1].split("_")) == 2:
			self.date, self.location = self.folders[i+1].split("_")[0], self.folders[i+1].split("_")[1]
		else:
			self.date = self.folders[i+1]
		# print(self.instrument_name, self.date, self.location, self.filters, self.azimut*RtoD, self.elevation*RtoD, self.com)

		rest = self.folders[i+2].split("_")
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
		# 	azimut, elevation = float(rest[0][1:]) * DtoR, float(rest[1][1:]) * DtoR
		# elif len(rest) > 2:
		# 	filters, azimut, elevation = rest[0], float(rest[1][1:]) * DtoR, float(rest[2][1:]) * DtoR
		# else:
		# 	com = folders[i+2].split("_")

		print(self.instrument_name, self.date, self.location, self.filters, self.azimut*RtoD, self.elevation*RtoD, self.com)


	def LoadData(self):
		# Loading the data files
		if self.instrument_name[:3] == "spp":
			self.raw_data = self.LoadSPPData(self.data_file_name)
			# if self.jump_mode == "length" or self.tail_jump.seconds == 0:
			# 	self.tail_jump = len(self.raw_data) - self.tail_jump.seconds
			for r in self.raw_data[:]:
				self.rotations.append(SPPRotation(r, instrument = self.instrument_name))

		elif self.instrument_name == "fake_spp":
			self.raw_data = self.LoadSPPTestData(1000, 100, 0.01, 45*RtoD)
			# if self.jump_mode == "length" or self.tail_jump.seconds == 0:
			# 	self.tail_jump = len(self.raw_data) - self.tail_jump.seconds
			for r in self.raw_data[:]:
				self.rotations.append(SPPRotation(r, fake = "True"))

		Bottle.LoadData(self)

	# def DateTime(self, moment="start"):
	# 	if moment == "start":
	# 		return self.rotations[0].datetime
	# 	if moment == "end":
	# 		return self.rotations[-1].datetime


	def SetTimes(self):
		self.time_stamp  = self.rotations[0].time_stamp
		#Time of obs in sec since EPOCH (except if before April 2019 -> in sec since midnight)
		self.time 		 = self.rotations[0].time
		#Date and time of first rotation (before deleting head_jump)
		self.datetime = self.rotations[0].datetime
		print("DEBUG SetTimes: ", self.time_stamp, self.time, self.datetime)
		# self.datetime = time.strptime(self.time_stamp, "%Y%m%d")
		# self.time = time.mktime(self.datetime)

	def CleanRotations(self):
		self.nb_rot = len(self.rotations)

		for i in range(1, len(self.rotations)):
			while self.rotations[i].time < self.rotations[i-1].time:
				self.rotations[i].time += time.timedelta(seconds = 24 * 3600)

		if self.jump_mode == "time": # and self.jump_unit in ["seconds", "minutes", "hours"]:
			self.head_jump = time.timedelta(seconds=float(self.input_parameters["head_jump"]))
			self.tail_jump = time.timedelta(seconds=float(self.input_parameters["tail_jump"]))

			if self.jump_unit == "minutes":
				self.head_jump = time.timedelta(minutes=float(self.input_parameters["head_jump"]))
				self.tail_jump = time.timedelta(minutes=float(self.input_parameters["tail_jump"]))
			elif self.jump_unit == "hours":
				self.head_jump = time.timedelta(hours=float(self.input_parameters["head_jump"]))
				self.tail_jump = time.timedelta(hours=float(self.input_parameters["tail_jump"]))

			if self.jump_mode == "length" or self.tail_jump.seconds == 0.:
				self.tail_jump = self.rotations[-1].time - self.tail_jump

			# print(self.head_jump, self.tail_jump, self.rotations[0].time, self.rotations[-1].time)
			print("Nb rotations before cleaning:", self.nb_rot)

			self.rotations = [r for r in self.rotations if self.head_jump <= r.time < self.tail_jump]

			print(self.nb_rot - len(self.rotations), "bad rotations in", self.nb_rot)
			self.nb_rot = len(self.rotations)
			print(self.nb_rot, "good rotations left after jump cuts.")

			print(self.head_jump, self.tail_jump, self.rotations[0].time, self.rotations[-1].time)


		self.rotations = [r for r in self.rotations if r.IsGood(Imin=self.Imin, Imax=self.Imax, DoLPmin=self.DoLPmin, DoLPmax=self.DoLPmax)]
		# print(self.Imin, self.Imax, self.DoLPmin, self.DoLPmax)


		print(self.nb_rot - len(self.rotations), "bad rotations in", self.nb_rot)
		self.nb_rot = len(self.rotations)
		print(self.nb_rot, "good rotations left after good rot cuts.")


	def SaveTXT(self):
		print("Saving as .txt in", self.data_file_name + "/" + self.saving_name + '_results.txt')
		times = [t.total_seconds() for t in self.all_times]
		np.savetxt("/".join(self.data_file_name.split("/")[:-1]) + "/" + self.saving_name + '_results.txt', 				   np.array([times, self.smooth_V, self.smooth_Vcos, self.smooth_Vsin, self.smooth_I0, self.smooth_DoLP, self.smooth_AoLP]).transpose(), delimiter = "\t", header = "time (ms)\tV\tVcos\tVsin\tI0 (mV)\tDoLP (percent)\tAoLP (deg)")


	def LoadSPPData(self, file_names):
		"""Return a data array from a given data file (1 or more). Each line is a rotation, each being a list of data as written below."""
		# 0-5: 		year, month, day, hour, min, sec
		# 6-9: 		voltage, positive voltage, negative voltage, motor voltage
		# 10-11: 	Temperature of base, of filter
		# 12-13:	elevation, azimuth (geog coord)
		# 14-93:	angle of polariser (degrees)
		# 94-173:	canal pola
		# 174-253:	canal total
		# 254-259:	V, Vcos, Vsin, I0, DoLP, AoLP
		nb_data_per_rot  = 254
		# raw_data = np.array([])
		# data_path = "/home/bossel/These/Analysis/data/spp/"
		#Loading the data from a binary file.
		#Output is one big 1-D array of length nb_data_tot = nb_rot * nb_data_per_rot
		# for f in file_names:
		# 	f = data_path + f + ".spp"
		# 	raw_data = np.fromfile(f, dtype = np.int16)
		# 	# raw_data = np.concatenate((raw_data, np.fromfile(f, dtype=np.int16)), axis=0)
		f = file_names + ".spp"
		raw_data = np.fromfile(f, dtype = np.int16)
		# raw_data = np.concatenate((raw_data, np.fromfile(f, dtype=np.int16)), axis=0)
		nb_data_tot = len(raw_data)
		print("nb_data_tot", nb_data_tot)

		# print(raw_data[1001 * 254 : 1002 * 254])
		#reshape the array to be (nb_pts, nb_data per pts)
		return raw_data.reshape(int(nb_data_tot / nb_data_per_rot), nb_data_per_rot)


	def LoadSPPTestData(self, nb_rot = 100, I0=100, D=0.1, phi=0, file_name = ''):
		"""Return a test data array with the time, angle and fake voltage mesured every point as SPP would. The fake signal is a function of the form A*cos(theta-phi)+B where theta is the angle of the polariser and phi the angle of the polarisation"""

		#I0:  Intensity of the incoming light of the aurora (constant through time) arbitrary units
		#D:  Degree of polarized light (%)
		#phi:  Angle of the polarisation with respect to SPP (degrees)

		### Do not take into account the fact taht theta-phi cannot be > 90, because it would only change a sign, but it is squared in the end. So no difference. But I keep the ITheo function in case.
		# I = lambda theta, I0, D, phi: D * I0 * np.cos((theta - phi)*DtoR)**2 + (1-D) * I0 / 2 + D/5*np.random.normal(0, 1, theta.shape)
		I = lambda theta, I0, D, phi: I0 + np.random.normal(0, I0 * D, theta.shape)

		data = np.zeros((nb_rot, 254))
		da = 360./80
		angles = np.linspace(da, 360., 80) * DtoR
		dolps = np.linspace(0, 1., nb_rot)
		phis = np.linspace(-90, 90, nb_rot)
		for i in range(nb_rot):
			data[i][5] = i
			data[i][14:94] = angles
			data[i][94:174] = I(angles, I0, D, phi)
	#		data[i][94:174] = I(angles, I0, D, phis[i])
			#data[i][94:174] = I(angles, I0, D, i*360/nb_rot)
			#data[i][94:174] = I(angles, I0, dolps[i], phi)
			data[i][174:254] = [I0]*80

			data[i][92:94] = [0,0]
			data[i][172:174] = [0,0]
			data[i][252:254] = [0,0]
		return data
