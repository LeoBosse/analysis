#!/usr/bin/python3
# -*-coding:utf-8 -*

import numpy as np
import matplotlib.pyplot as plt
from utils import *
from scipy import signal
import sys as sys
import datetime as time

DtoR = np.pi / 180.
RtoD = 180. / np.pi

### Rotaion classes. Rotation object is the parent of SPPRotaion and PTCURotation objects. Both classes exist because of the initialisations that are a little different and SPPRotations needs more methods. But they still share some common things.

class Rotation:
	nb_bad_rot = 0
	def __init__(self, rotation_data, instrument="spp"):
		"""Initialise the Rotation class."""
		self.raw_data = rotation_data
		self.instrument = instrument.lower()

		self.V 	  	= 0
		self.Vcos 	= 0
		self.Vsin 	= 0

		self.I0		= 0
		self.DoLP	= 0
		self.AoLP	= 0

		self.is_good = True

		self.TempPM = 0
		self.TempOptical = 0
		self.TempAmbiant = 0

	def GetLightParameters(V, Vcos, Vsin):
		"""Given V, Vcos, Vsin, returns the initial intensity, DoLP and AoLP. This method is shared for spp and ptcu. It is also a static method, that can be called outside of the object. This way it can be used everywhere, each time you need I0, DoLP, AoLP to decrease the chances of making a mistake."""
		I0 = 2 * V
		DoLP = 2 * np.sqrt(Vcos**2 + Vsin**2) / V * 100
		AoLP = np.arctan2(Vsin, Vcos) / 2
		return abs(I0), abs(DoLP), AoLP
	GetLightParameters = staticmethod(GetLightParameters)

	def IsGood(self, Imin=0.1, Imax = 10000, DoLPmin=0, DoLPmax=100):
		"""Check if the rotation is a "good" rotation. the Vcos Vsin shouldn't be too high (<10), The DoLP < 30%, and all important parameters are a number (are defined)."""
		good = True
		### Usefull if real data because Very unlikely == bad
		if  self.I0 < Imin or self.I0 > Imax or self.DoLP >= DoLPmax or self.DoLP <= DoLPmin: #or abs(self.Vcos) >= 10. or abs(self.Vsin) >= 10.:
			 good = False

		if True in np.isnan([self.V, self.Vcos, self.Vsin, self.I0, self.DoLP, self.AoLP]) :
			# print(np.isnan([self.V, self.Vcos, self.Vsin, self.I0, self.DoLP, self.AoLP]))
			good = False

		if not good:
			Rotation.nb_bad_rot += 1

		# good = True
		return good

	def GetObservationTime(self, unit = "seconds"):
		"""Return the time of the rotation since the observation began. Several units can be specified: seconds (by default), milliseconds, ms, minutes, min, m, hours, hrs, h."""
		unit = unit.lower()
		if unit == "seconds":
			return self.time
		elif unit == "milliseconds" or unit == "ms":
			return self.time * 1000
		elif unit == "minutes" or unit == "min" or unit == "m":
			return self.time / 60
		elif unit == "hours" or unit == "hrs" or unit == "h":
			return self.time / 3600.
		else:
			print("WARNING: Time unit specified was incorrect. Using seconds by default.")
			return self.time

	def WriteTimeStamp(self):
		string = str(self.year) + "-" + str(self.month) + "-" + str(self.day)
		string += " "
		string += str(self.obs_loc_hour) + ":" + str(self.obs_loc_min) + ":" + str(self.obs_loc_sec)
		return string

	def GetTimeStamp(self, observation_start = 0, format = "local"):
		if format == "local":
			h = int(self.time / 3600) + self.obs_loc_hour
			if h > 23:
				h -= 24
				self.day += 1
			m = int(self.time / 60) + self.obs_loc_min
			s = int(self.time) + self.obs_loc_sec

		return self.WriteTimeStamp()


class PTCURotation(Rotation):
	def __init__(self, data):
		"""Initialise the PTCURotation class. Load important data from the input and calculate I0, DoLP, AoLP, then check if 'good'."""
		Rotation.__init__(self, data, instrument="PTCU" )
		###DATA INDICES
		# 0: 	IDConfiguration,
		# 1: 	Time,
		# 2: 	IDProcessedData,
		# 3:	PMAvg,
		# 4:	PMCosAvg,
		# 5:	PMSinAvg

		# self.config = config
		# ###CONFIG INDICES
		# # 0:		IDConfiguration,
		# # 1:		Timestamp,
		# # 2-3:		CM_Latitude, CM_Longitude,
		# # 4-5:		CM_Elevation, CM_Azimuth,
		# # 6:		CM_PolarizerSpeed,
		# # 7:		CM_Tilt,
		# # 8:		CM_Usage,
		# # 9:		CM_Comments,
		# # 10:		IS_PolarizerOffset1,
		# # 11:		IS_PolarizerOffset2,
		# # 12:		IS_PolarizerOffset3,
		# # 13:		IS_PolarizerOffset4,
		# # 14:		IS_ConverterOffset4,
		# # 15:		IS_ConverterOffset3,
		# # 16:		IS_ConverterOffset2,
		# # 17:		IS_ConverterOffset1,
		# # 18:		IS_EngineTicksPerTour,
		# # 19:		IS_MotorGearedRatio,
		# # 20:		IS_QuantumEfficiency,
		# # 21:		IS_IpAddress

		### Time since the begining of the observation (in time stamp) given in milliseconds, but stored in seconds, as all other times

		self.nb_data_per_rot = len(self.raw_data)
		# self.time				= time.timedelta(milliseconds = float(self.raw_data[1]))
		# print(self.nb_data_per_rot)
		# print(self.raw_data)
		if self.nb_data_per_rot == 10:
			self.time				= time.timedelta(milliseconds = float(self.raw_data[1]))
			self.V					= self.raw_data[2]
			self.Vcos				= self.raw_data[3]
			self.Vsin				= self.raw_data[4]
			self.Vref				= self.raw_data[5]
			self.Comment 			= self.raw_data[6]
			self.TempPM 			= self.raw_data[7]
			self.TempOptical 		= self.raw_data[8]
			self.TempAmbiant 		= self.raw_data[9]
		elif self.nb_data_per_rot == 9:
			self.time				= time.timedelta(milliseconds = float(self.raw_data[1]))
			self.V					= self.raw_data[2]
			self.Vcos				= self.raw_data[3]
			self.Vsin				= self.raw_data[4]
			self.Vref				= False
			self.Comment 			= self.raw_data[5]
			self.TempPM 			= self.raw_data[6]
			self.TempOptical 		= self.raw_data[7]
			self.TempAmbiant 		= self.raw_data[8]
		elif self.nb_data_per_rot == 13:
			self.time				= time.timedelta(milliseconds = float(self.raw_data[1]))
			self.V					= self.raw_data[3]
			self.Vcos				= self.raw_data[4]
			self.Vsin				= self.raw_data[5]
			self.Vref				= False
			self.IDTemperature		= self.raw_data[6]
			self.TempPM 			= self.raw_data[7]
			self.TempOptical 		= self.raw_data[8]
			self.TempAmbiant 		= self.raw_data[9]
			self.IDLiveComment 		= self.raw_data[10]
			self.Comment 			= self.raw_data[11]
			self.live_Commentscol 	= self.raw_data[12]
		elif self.nb_data_per_rot == 6:
			self.time				= time.timedelta(milliseconds = float(self.raw_data[2]))
			self.V					= self.raw_data[3]
			self.Vcos				= self.raw_data[4]
			self.Vsin				= self.raw_data[5]
			self.Vref				= False
			self.IDTemperature		= False
			self.TempPM 			= False
			self.TempOptical 		= False
			self.TempAmbiant 		= False
			self.IDLiveComment 		= False
			self.Comment 			= False
			self.live_Commentscol 	= False
		elif self.nb_data_per_rot == 4:
			self.time				= time.timedelta(milliseconds = float(self.raw_data[0]))
			self.V					= self.raw_data[1]
			self.Vcos				= self.raw_data[2]
			self.Vsin				= self.raw_data[3]
			self.Vref				= False
			self.IDTemperature		= False
			self.TempPM 			= False
			self.TempOptical 		= False
			self.TempAmbiant 		= False
			self.IDLiveComment 		= False
			self.Comment 			= False
			self.live_Commentscol 	= False
		# self.Vcos	= self.raw_data[4]
		# self.Vsin	= self.raw_data[5]

		# self.tilt 	= self.config["CM_Tilt"]

		### Get I0, DoLP, AoLP
		self.V = np.absolute(self.V)
		self.I0, self.DoLP,	self.AoLP = self.GetLightParameters(self.V, self.Vcos, self.Vsin)
		self.Iref = 2 * self.Vref
		### Check if 'good'
		self.is_good = self.IsGood()

class SPPRotation(Rotation):
	def __init__(self, rotation_data, fake = False, instrument = "spp"):
		"""Initialise the SPPRotation class. Load important data from the input and calculate directly V, Vcos, Vsin, I0, DoLP, AoLP. Clean the rotation and check if good."""
		Rotation.__init__(self, rotation_data, instrument=instrument)

		self.nb_data_per_rot 	= 254
		self.nb_pts_per_rot 	= 80

		self.date_stamp = ""
		self.time_stamp = ""
		self.year		= int(self.raw_data[0])
		self.date_stamp += str(self.year) + "-"
		self.month		= int(self.raw_data[1])
		self.date_stamp += str(self.month) + "-"
		self.day		= int(self.raw_data[2])
		self.date_stamp += str(self.day)
		self.hour		= int(self.raw_data[3])
		self.time_stamp += str(self.hour) + ":"
		self.min		= int(self.raw_data[4])
		self.time_stamp += str(self.min) + ":"
		self.sec		= int(self.raw_data[5])
		self.time_stamp += str(self.sec)

		### All the times are stored in seconds
		#Date and time of rotation
		self.datetime				= time.datetime.strptime(self.date_stamp + " " + self.time_stamp, "%Y-%m-%d %H:%M:%S")
		#Time of obs in sec since EPOCH (except if before April 2019)
		self.time_since_EPOCH		= time.timedelta(seconds = self.datetime.timestamp())
		#Time in sec since midnight of this day
		self.time_since_midnight 	=  self.datetime - time.datetime(self.year, self.month, self.day)

		if  self.datetime < time.datetime.strptime("2019-03-10", "%Y-%m-%d"): #I did this differently before this date, but all the input files are already done and I don't want to change everything...
			self.time = self.time_since_midnight
		else:
			self.time = self.time_since_EPOCH

		self.elevation	= self.raw_data[12]
		self.azimuth	= self.raw_data[13]

		### self.raw_data[14:94] ARE NOT THE CORRECT ANGLES? WEIRD VALUES INSIDE -> SET ANGLES BY HAND...
		# self.angles		= self.raw_data[14:94]
		# self.angles		= np.linspace(0, 355.5, self.nb_pts_per_rot) * DtoR
		self.angles		= np.linspace(4.5, 360, self.nb_pts_per_rot) * DtoR

		### in 2014, the pola and total canal were swapped -> instrument_name for 2014 spp data should be spp2014
		if self.instrument == "spp2014":
			self.I_pola		= self.raw_data[174:254]
			self.I_total	= self.raw_data[94:174]
		else:
			self.I_pola		= self.raw_data[94:174]
			self.I_total	= self.raw_data[174:254]

		### Clean the rotation
		self.CleanRotation()

		### Calculate V, Vcos, Vsin
		self.V 	  	= self.GetV()
		self.Vref	= self.GetVref()
		self.Vcos 	= self.GetVcos()
		self.Vsin 	= self.GetVsin()

		self.Iref = 2 * self.Vref
		### calculate I0, DoLP, AoLP
		self.I0, self.DoLP,	self.AoLP = self.GetLightParameters(self.V, self.Vcos, self.Vsin)

		### Check if the rotation is 'good'
		self.is_good = self.IsGood()

	def CleanRotation(self):
		"""Clean the rotation. Some points are sometime 0 at the end of a rotation. When it is the case, put the angle by hand, replace the intensity for both canal by the intensity of the previous point. I was first deleting them, but it is more complicated and DOESN'T WORK properly. This seems funny that only two or three points can mess with the results so much when removed. """
		# ### This commented section does not work, don't delete entries, but replaces them with neighbouring values. Maybe the Vcos and Vsin work differently when we don't consider a full rotation
		# to_del = []
		# for i, a in enumerate(self.angles):
		# 	if a == 0.:
		# 		self.nb_pts_per_rot -= 1
		# 		to_del.append(i)
		# self.angles  = np.delete(self.angles,  to_del)
		# self.I_pola  = np.delete(self.I_pola,  to_del)
		# self.I_total = np.delete(self.I_total, to_del)
		# self.nb_data_per_rot = self.nb_data_per_rot - 3 * len(to_del)

		### If the canal pola is 0, replace it with the previous value.
		for i, a in enumerate(self.I_pola):
			if a == 0.:
				if i == 0:
					self.angles[i]  = 4.5 * DtoR
					self.I_pola[i]  = self.I_pola[1]
					self.I_total[i] = self.I_total[1]
				else:
					self.angles[i]  = self.angles[i-1] + 4.5 * DtoR
					self.I_pola[i]  = self.I_pola[i-1]
					self.I_total[i] = self.I_total[i-1]


	def IsGood(self, Imin=0.1, Imax=10000, DoLPmin=0, DoLPmax=100):
		"""Check if the rotat
		ion is "good" like the Rotation class, but also check if the pola intensity is not 0 everywhere."""
		good = Rotation.IsGood(self, Imin, Imax, DoLPmin, DoLPmax)
		if good:
			if self.V == 0.:
				good = False
				Rotation.nb_bad_rot += 1
				return good
		return good

	def GetV(self):
		"""return the average of the signal this rotation"""
		integral = sum(self.I_pola)
		return integral / self.nb_pts_per_rot

	def GetVref(self):
		"""return the average of the signal this rotation"""
		integral = sum(self.I_total)
		return integral / self.nb_pts_per_rot

	def GetVcos(self):
		"""Return the average of V*cos(2*theta) of this rotation"""
		x = self.I_pola * np.cos(2 * self.angles)
		return sum(x) / self.nb_pts_per_rot

	def GetVsin(self):
		"""Return the average value of -V*sin(2*theta) of this rotation"""
		y = - self.I_pola * np.sin(2 * self.angles)
		return sum(y) / self.nb_pts_per_rot

	def __repr__(self):
		s ="SPP rotation.\n"
		s += str(self.raw_data)
		return s
