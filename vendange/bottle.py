#!/usr/bin/python3
# -*-coding:utf-8 -*

import numpy as np
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

### Bottle classes.

class Bottle:
	def __init__(self, folder_name, auto=True, line=1):
		self.folder = folder_name
		self.rotations = []

		self.AoBapp = False
		self.AoBlos = False
		self.AoRD = False
		self.time_stamp = ""

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

		self.date, self.location, self.filters, self.azimut, self.elevation, self.com = "", "", "", False, False, ""

		try:
			self.azimut = float(self.input_parameters["azimut"]) * DtoR
			self.elevation = float(self.input_parameters["elevation"]) * DtoR
		except:
			print("WARNING: No observation angles info in the input file.")

		### See PTCUBottle or SPPBottle for specific instructions

	def SetInfoFromDataFileName(self, f, data_f):
		self.folders, self.data_folders = f.split("/"), data_f.split("/")
		# i = len(data_folders) - 1
		self.folders_index = -3

		# self.instrument_name = self.folders[self.folders_index]

		### Follows in SPPBottle or PTCUBottle

	def MiseEnBouteille(self):
		### Loading the data from the files, then creating a list of Rotation objects depending on the instrument name specified. Can be 'spp', 'spp2014', 'fake_spp', 'ptcu', 'fake_ptcu'
		self.SetJumps()
		self.LoadData()
		if self.valid:
			self.CleanRotations()

		### Now that we have all our rotations, creating lists of usefull data for easier smoothing
			self.CreateLists()
			self.GetSmoothLists()

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
			self.head_jump = time.timedelta(0)
			self.tail_jump = time.timedelta(0)
		# elif self.jump_unit == "minutes":
		# 	self.head_jump = 0
		# 	self.tail_jump = 0
		# elif self.jump_unit == "hours":
		# 	self.head_jump = 0
		# 	self.tail_jump = 0

	def LoadData(self):
		# Loading the data files

		### See SPPBottle/PTCUBottle for first part

		self.nb_rot = len(self.rotations)
		self.SetConfig()
		self.SetTimes()


	def SetConfig(self):
		try:
			self.Imin = float(self.input_parameters["Imin"])
		except:
			self.Imin = 0.1
		try:
			self.Imax = float(self.input_parameters["Imax"])
		except:
			self.Imax = 10000
		try:
			self.DoLPmax = float(self.input_parameters["DoLPmax"])
		except:
			self.DoLPmax = 100
		try:
			self.DoLPmin = float(self.input_parameters["DoLPmin"])
		except:
			self.DoLPmin = 0


	def CreateLists(self):
		###At this point, all times are expressed in seconds since the first rotation before deleting the head_jump
		###If SPP before 2019 03 10, all times are expressed in seconds since midnight of the first observation day.

		# print("DEBUG TIME LISTS: ", self.time, self.head_jump, self.rotations[0].time)
		# self.all_times = np.array([r.time - self.head_jump for r in self.rotations])
		self.all_times = np.array([r.time for r in self.rotations])
		# self.all_times_since_start = np.array([r.time  for r in self.rotations])
		self.all_times_since_start = np.array([r.time - self.rotations[0].time for r in self.rotations])

		# self.all_datetimes = [time.gmtime(r.time + self.time + self.head_jump) for r in self.rotations]
		# self.all_datetimes = [time.gmtime(time.mktime(self.DateTime()) + r.time + self.time + self.head_jump) for r in self.rotations]

		self.all_datetimes = [self.DateTime() + t for t in self.all_times]
		# self.all_datetimes = [self.DateTime() + t for t in self.all_times_since_start]

		print(self.head_jump, self.all_times[0], self.all_times[-1], self.all_datetimes[0], self.all_datetimes[-1])

		# for i in range(1, self.nb_rot):
		# 	if self.all_times[i] < self.all_times[i-1]:
		# 		self.all_times[i] += 24 * 3600

		# self.AoLP_correction = float(self.input_parameters["AoLP_correction"])*DtoR
		self.AoLP_correction = self.config['IS_PolarizerOffset' + str(self.line)] * DtoR
		print("AoLP correction:", self.AoLP_correction*RtoD)


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

	def GetSmoothLists(self):
		print("Get Smooth Lists")
		### Smoothing procedure. Can smooth for a given number of rotations (smoothing_unit==rotations) or a  given time period (smoothing_unit==seconds).
		self.smoothing_factor = int(self.input_parameters["smoothing_factor"]) #Average the data over smoothing_factor rotations
		self.smoothing_unit = self.input_parameters["smoothing_unit"].lower()
		self.smoothing_method = self.input_parameters["smoothing_method"].lower()

		self.nb_smooth_rot = int(self.nb_rot / self.smoothing_factor)

		self.avg_dt = 1000 * np.average([self.all_times[i].total_seconds() - self.all_times[i-1].total_seconds() for i in range(1, len(self.all_times))]) #in millisec

		print("AVG DT (millisec)", self.avg_dt)

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

		#Get Variance
		if self.smoothing_unit == "seconds": smoothing_factor = self.smoothing_factor*1000
		elif self.smoothing_unit == "minutes": smoothing_factor = self.smoothing_factor * 60*1000
		elif self.smoothing_unit == "hours": smoothing_factor = self.smoothing_factor * 3600*1000



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
		self.std_I0 = 		 np.sqrt(4 * self.all_V 	/ self.avg_dt)
		self.std_smooth_I0 = np.sqrt(4 * self.smooth_V  / smoothing_factor)

		self.std_DoLP = 		np.sqrt(4 * (1 + ((self.all_DoLP/100) 	 ** 2 / 2)) / (self.all_I0 	  * self.avg_dt)) * 100
		self.std_smooth_DoLP = 	np.sqrt(4 * (1 + ((self.smooth_DoLP/100) ** 2 / 2)) / (self.smooth_I0 * smoothing_factor)) * 100

		self.std_AoLP = 		np.sqrt(1 / ((self.all_DoLP/100) 	** 2 * self.all_I0 	  * self.avg_dt))
		self.std_smooth_AoLP = 	np.sqrt(1 / ((self.smooth_DoLP/100) ** 2 * self.smooth_I0 * smoothing_factor))
		self.smooth_AoLP_upper = self.smooth_AoLP + self.std_smooth_AoLP
		self.smooth_AoLP_lower = self.smooth_AoLP - self.std_smooth_AoLP
		# print(self.smooth_DoLP.shape, self.sm ooth_AoLP.shape, self.var_DoLP.shape, self.var_AoLP.shape)

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
			SetAngleBounds(self.AoLP_average, 0, np.pi)
			SetAngleBounds(self.smooth_AoLP_upper, 0, np.pi)
			SetAngleBounds(self.smooth_AoLP_lower, 0, np.pi)
			# if self.AoLP_average < 0:
			# 	self.AoLP_average += np.pi
		elif self.graph_angle_shift == 0:
			SetAngleBounds(self.AoLP_average, -np.pi/2, np.pi/2)
			SetAngleBounds(self.smooth_AoLP_upper, -np.pi/2, np.pi/2)
			SetAngleBounds(self.smooth_AoLP_lower, -np.pi/2, np.pi/2)
			# if self.AoLP_average > np.pi/2:
			# 	self.AoLP_average -= np.pi

		# for i, up in enumerate(self.smooth_AoLP_upper):
		# 	if up < self.smooth_AoLP[i]:




		print("Get Smooth Lists: DONE")


	def GetGeometryAngles(self, obs):
		obs.SinglePointGeometry()
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
		if self.location.lower() == "skibotn":
			A_lon = 20.24 * DtoR
			A_lat = 69.34 * DtoR
		elif self.location.lower() == "mens":
			A_lon = 5.76 * DtoR
			A_lat = 44.83 * DtoR
		elif self.location.lower() == "nyalesund":
			A_lon = 11.92288 * DtoR
			A_lat = 78.92320 * DtoR
		elif self.location.lower() == "vigan":
			A_lon = 3.504259 * DtoR
			A_lat = 44.039661 * DtoR
		elif self.location.lower() == "lagorge":
			A_lon = 5.936935 * DtoR
			A_lat = 45.212343 * DtoR
		elif self.location.lower() == "stveran":
			A_lon = 6.5430 * DtoR
			A_lat = 44.4156 * DtoR
		elif self.location.lower() == "sob":
			A_lon = -16.452403 * DtoR
			A_lat = 14.49496 * DtoR
		elif self.location.lower() == "skibotnsud":
			A_lon = 19.981116 * DtoR
			A_lat = 69.234956 * DtoR
		elif self.location.lower() == "kilpisjarvi":
			A_lon = 20.7846926 * DtoR
			A_lat = 69.0526675 * DtoR

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

		if self.observation_type == "fixed" and self.azimut and self.elevation:
			print("DEBUG SOURCE:", self.source_azimut*RtoD, self.source_elevation*RtoD)
			try:
				print("DEBUG: AZ/EL", self.azimut*RtoD, self.elevation*RtoD)
				self.observation, self.AoBapp, self.AoBlos, self.AoRD, self.AoBapp_ortho, self.AoRD_ortho = self.GetGeometryAngles(observation.ObservationPoint(A_lon, A_lat, h, self.azimut, self.elevation, self.source_azimut, self.source_elevation))
			except:
				print("WARNING: No geometric angles where retrieved.")

		elif self.observation_type == "fixed_elevation_discrete_rotation":
			if to_initiate:
				self.discrete_rotation_elevation 	= float(self.input_parameters["discrete_rotation_elevation"]) * DtoR
				self.discrete_rotation_times 		= np.array([float(t) for t in self.input_parameters["discrete_rotation_times"].split("_")])
				self.discrete_rotation_azimuts 		= np.array([float(a)*DtoR for a in self.input_parameters["discrete_rotation_azimuts"].split("_")])

			ang_list = []
			for a in self.discrete_rotation_azimuts:
				ang_list.append(self.GetGeometryAngles(observation.ObservationPoint(A_lon, A_lat, h, a, self.discrete_rotation_elevation, self.source_azimut, self.source_elevation)))

			self.observation, self.AoBapp, self.AoBlos, self.AoRD, self.AoBapp_ortho, self.AoRD_ortho = zip(*ang_list)
			self.observation, self.AoBapp, self.AoBlos, self.AoRD, self.AoBapp_ortho, self.AoRD_ortho = np.array(self.observation), np.array(self.AoBapp), np.array(self.AoBlos), np.array(self.AoRD), np.array(self.AoBapp_ortho), np.array(self.AoRD_ortho)

		elif self.observation_type == "fixed_azimut_discrete_rotation":
			if to_initiate:
				self.discrete_rotation_azimut 	= float(self.input_parameters["discrete_rotation_azimut"]) * DtoR
				self.discrete_rotation_times 		= [float(t) for t in self.input_parameters["discrete_rotation_times"].split("_")]
				self.discrete_rotation_times_unit = self.input_parameters["discrete_rotation_times_unit"]
				self.discrete_rotation_elevations 		= [float(a)*DtoR for a in self.input_parameters["discrete_rotation_elevations"].split("_")]

			ang_list = []
			for e in self.discrete_rotation_elevations:
				ang_list.append(self.GetGeometryAngles(observation.ObservationPoint(A_lon, A_lat, h, self.discrete_rotation_azimut, e, self.source_azimut, self.source_elevation)))

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
				self.azimut, self.obs = geo.FixedElevation(direction = self.input_parameters["rotation_direction"])
			except:
				self.azimut, self.obs = geo.FixedElevation()


			ang_list = [self.GetGeometryAngles(o) for o in self.obs]
			self.obs, self.AoBapp, self.AoBlos, self.AoRD, self.AoBapp_ortho, self.AoRD_ortho = zip(*ang_list)
			self.obs, self.AoBapp, self.AoBlos, self.AoRD, self.AoBapp_ortho, self.AoRD_ortho = np.array(self.obs), np.array(self.AoBapp), np.array(self.AoBlos), np.array(self.AoRD), np.array(self.AoBapp_ortho), np.array(self.AoRD_ortho)



		# try:
		# 	print(self.azimut, self.elevation)
		# 	self.observation = observation.ObservationPoint(A_lon, A_lat, h, self.azimut, self.elevation)
		# 	self.observation.SinglePointGeometry()
		# 	self.AoBapp = self.observation.eta_chaos
		# 	self.AoBapp_ortho = self.AoBapp + np.pi / 2
		# 	self.AoBlos = self.observation.Blos
		#
		# 	if self.graph_angle_shift == 1 and self.AoBapp < 0:
		# 		self.AoBapp += np.pi
		# 	if self.graph_angle_shift == 0 and self.AoBlos > np.pi/2:
		# 		self.AoBlos -= np.pi
		# 	if self.graph_angle_shift == 0 and self.AoBapp_ortho > np.pi/2:
		# 		self.AoBapp_ortho -= np.pi
		#
		# 	print("DEBUG B", self.AoBapp*RtoD)
		# except:
		# 	print("WARNING: Could not get the apparent angle of the magnetic field.")
		#
		# try:
		# 	self.source_azimut = float(self.input_parameters["pollution_source_azimut"]) * DtoR
		# 	try:
		# 		self.source_elevation = float(self.input_parameters["pollution_source_elevation"]) * DtoR
		# 	except:
		# 		self.source_elevation = 0
		# 	self.AoRD = self.observation.GetRayleighAngle(self.source_azimut, self.source_elevation)
		# 	self.AoRD_ortho = self.AoRD + np.pi/2
		#
		# 	if self.graph_angle_shift == 1:
		# 		if self.AoRD < 0:
		# 			self.AoRD += np.pi
		# 	elif self.graph_angle_shift == 0:
		# 		if self.AoRD_ortho > np.pi/2:
		# 			self.AoRD_ortho -= np.pi
		# except:
		# 	print("WARNING: No source pollution azimut to get the Angle of Rayleigh Diffusion (AoRD)")

		os.chdir(wd)
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



#####################################################################################
###									Petit Cru									  ###
#####################################################################################


class PTCUBottle(Bottle):
	def __init__(self, folder_name, auto=True, line = 1):
		Bottle.__init__(self, folder_name, auto, line = line)

		self.SetInfoFromDataFileName(self.data_file_name, pwd_data)

		if auto:
			self.MiseEnBouteille()


	def SetInfoFromDataFileName(self, f, data_f):
		Bottle.SetInfoFromDataFileName(self, f, data_f)

		i = self.folders_index
		date_location = self.folders[i+1].split("_")
		self.date, self.location = date_location[0], date_location[1]
		rest = self.folders[i+2].split("_")
		print(rest)
		#Filters: r,v,b,m pour rouge, vert, bleu, mauve
		#0: no filters_list
		#o: orange 620
		#X, Y: 620nm et 413nm
		filters_list = ["r", "v", "b", "m", "0", "o", "X", "Y"]
		nb_filters = 2
		if self.instrument_name == "gdcu":
			nb_filters = 4

		for r in rest:
			if len(r) == nb_filters and np.all([a in filters_list for a in r]):
				if self.instrument_name == "carmen" or self.instrument_name == "fake_ptcu":
					self.filters = r
					if self.filters[1] in [0, "0"]:
						self.NoVref = True
				elif self.instrument_name in ["corbel", "gdcu"]:
					self.filters = r[self.line - 1]
				print("FILTERS:", r, self.filters)

			elif r[0] == "a":
				self.azimut = float(r[1:]) * DtoR
			elif r[0] == "e":
				self.elevation = float(r[1:]) * DtoR
			else:
				self.com += "_" + r

		print(self.instrument_name, self.date, self.location, self.filters, self.azimut*RtoD, self.elevation*RtoD, self.com)

	def LoadData(self):
		if self.instrument_name == "carmen" or self.instrument_name == "corbel" or self.instrument_name == "gdcu":
			self.raw_data, self.config = LoadPTCUData(self.data_file_name, line = self.line)
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
		self.time = self.datetime.timestamp()

	def DateTime(self, moment="start"):
		if moment == "start":
			return self.rotations[0].time + self.datetime + self.head_jump
		if moment == "end":
			return self.rotations[-1].time + self.datetime + self.head_jump


	def CleanRotations(self):
		self.nb_rot = len(self.rotations)

		### Put every rotation time in sec since first rotation
		norm = self.rotations[0].time
		for i in range(len(self.rotations)):
			self.rotations[i].time -= norm

		### Add 24h when next rotation is earlier than the previous one. Happens sometimes when we change the date
		for i in range(1, len(self.rotations)):
			while self.rotations[i].time < self.rotations[i-1].time:
				self.rotations[i].time += 24 * 3600

		### Get and set the head and tail jumps
		if self.jump_mode == "time": # and self.jump_unit in ["seconds", "minutes", "hours"]:
			self.head_jump = time.timedelta(seconds=float(self.input_parameters["head_jump"]))
			self.tail_jump = time.timedelta(seconds=float(self.input_parameters["tail_jump"]))

			if self.jump_unit == "minutes":
				self.head_jump = time.timedelta(minutes=float(self.input_parameters["head_jump"]))
				self.tail_jump = time.timedelta(minutes=float(self.input_parameters["tail_jump"]))
			elif self.jump_unit == "hours":
				print("DEBUG jumps", [r.time.total_seconds() for r in self.rotations[:10]])
				print(time.timedelta(hours=float(self.input_parameters["head_jump"])), time.timedelta(hours=float(self.input_parameters["tail_jump"])))
				self.head_jump = time.timedelta(hours=float(self.input_parameters["head_jump"]))
				self.tail_jump = time.timedelta(hours=float(self.input_parameters["tail_jump"]))

			if self.jump_mode == "length" or self.tail_jump.seconds == 0.:
				self.tail_jump = self.rotations[-1].time - self.tail_jump

			### Delete all rotations before the head and after the tail jump
			self.rotations = [r for r in self.rotations if self.head_jump <= r.time < self.tail_jump]
			print(len(self.rotations), "good rotations in", self.nb_rot, ";", self.nb_rot - len(self.rotations), "deleted because of time jumps.")
			self.nb_rot = len(self.rotations)

			print(self.head_jump, self.tail_jump, self.rotations[0].time, self.rotations[-1].time)


			### Put every rotation time in sec since first rotation (after deleting head_jump)
			norm = self.rotations[0].time
			for i in range(len(self.rotations)):
				self.rotations[i].time -= norm

		self.rotations = [r for r in self.rotations if r.IsGood(Imin=self.Imin, Imax=self.Imax, DoLPmin=self.DoLPmin, DoLPmax=self.DoLPmax)]
		print(len(self.rotations), "good rotations in", self.nb_rot, ";", self.nb_rot - len(self.rotations), "deleted because of invalid data.")
		self.nb_rot = len(self.rotations)

	def SaveTXT(self):
		times = [t.total_seconds() for t in self.all_times]
		print("Saving as .txt in", self.data_file_name + "/" + self.saving_name + '_results.txt')

		np.savetxt(self.data_file_name + "/" + self.saving_name + '_results.txt', np.array([times, self.smooth_V, self.smooth_Vcos, self.smooth_Vsin, self.smooth_I0, self.smooth_DoLP, self.smooth_AoLP]).transpose(), delimiter = "\t", header = "time(s)\tSV\tSVcos\tSVsin\tSI0(mV)\tSDoLP(percent)\tSAoLP(deg)")



#####################################################################################
###									SPPBottle									  ###
#####################################################################################


class SPPBottle(Bottle):
	def __init__(self, folder_name, auto=True):
		Bottle.__init__(self, folder_name, auto)
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
			self.raw_data = LoadSPPData(self.data_file_name)
			# if self.jump_mode == "length" or self.tail_jump.seconds == 0:
			# 	self.tail_jump = len(self.raw_data) - self.tail_jump.seconds
			for r in self.raw_data[:]:
				self.rotations.append(SPPRotation(r, instrument = self.instrument_name))

		elif self.instrument_name == "fake_spp":
			self.raw_data = LoadSPPTestData(1000, 100, 0.01, 45*RtoD)
			# if self.jump_mode == "length" or self.tail_jump.seconds == 0:
			# 	self.tail_jump = len(self.raw_data) - self.tail_jump.seconds
			for r in self.raw_data[:]:
				self.rotations.append(SPPRotation(r, fake = "True"))

		Bottle.LoadData(self)

	def DateTime(self, moment="start"):
		if moment == "start":
			return self.rotations[0].datetime
		if moment == "end":
			return self.rotations[-1].datetime


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

			print(self.head_jump, self.tail_jump, self.rotations[0].time, self.rotations[-1].time)
			print(Rotation.nb_bad_rot, "bad rotations in", self.nb_rot)

			self.rotations = [r for r in self.rotations if self.head_jump <= r.time < self.tail_jump]
			self.nb_rot = len(self.rotations)

			print(Rotation.nb_bad_rot, "bad rotations in", self.nb_rot)
			print(self.head_jump, self.tail_jump, self.rotations[0].time, self.rotations[-1].time)


		print(Rotation.nb_bad_rot, "bad rotations in", self.nb_rot)

		self.rotations = [r for r in self.rotations if r.IsGood(Imin=self.Imin, Imax=self.Imax, DoLPmin=self.DoLPmin, DoLPmax=self.DoLPmax)]
		self.nb_rot = len(self.rotations)

		print(Rotation.nb_bad_rot, "bad rotations in", self.nb_rot)


	def SaveTXT(self):
		print("Saving as .txt in", self.data_file_name + "/" + self.saving_name + '_results.txt')
		times = [t.total_seconds() for t in self.all_times]
		np.savetxt("/".join(self.data_file_name.split("/")[:-1]) + "/" + self.saving_name + '_results.txt', 				   np.array([times, self.smooth_V, self.smooth_Vcos, self.smooth_Vsin, self.smooth_I0, self.smooth_DoLP, self.smooth_AoLP]).transpose(), delimiter = "\t", header = "time (ms)\tV\tVcos\tVsin\tI0 (mV)\tDoLP (percent)\tAoLP (deg)")
