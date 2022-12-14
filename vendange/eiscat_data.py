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

class Eiscat:
	def __init__(self, bottle, hdf5 = False):

		self.valid = True
		self.path = global_configuration.eiscat_data_path
		# self.path = data_path + "eiscat/"

		self.data = []

		self.InitTime(bottle)
		self.SetPaths(bottle)

		if self.valid:
			self.GetDataFromFiles()

		if self.valid:
			print("EISCAT data exist! Number of points:", len(self.data))
			print(bottle.DateTime(format="UT").strftime("%Y-%m-%d"))
			# print(self.data)
		else:
			self.valid = False
			print("No EISCAT data available.")

	def SetPaths(self, bottle):
		self.folder = self.path
		if self.start_datetime.year == 2020:
			self.folder += "SP_FR_anadata/"

		self.folder += bottle.DateTime(format="UT").strftime("%Y-%m-%d")

		if self.start_datetime.year == 2019:
			self.folder += "_beata_10@uhfa/"
		elif self.start_datetime.year == 2020:
			self.folder += "_beata_60@vhf/"

		self.folders = [self.folder]
		if self.start_datetime.day != self.end_datetime.day:
			one_day = dt.timedelta(days=1)
			start_date = bottle.DateTime("start", format="UT")
			end_date = bottle.DateTime("end", format="UT")
			while start_date.day != end_date.day:
				start_date += one_day
				folder = self.path
				folder += start_date.strftime("%Y-%m-%d")
				folder += "_beata_10@uhfa/"
				self.folders.append(folder)
		# print(self.folders)
		self.file_list = []
		for f in self.folders:
			try:
				self.file_list.extend([f + fs for fs in os.listdir(f)])
			except:
				print("WARNING: No eiscat data files for this date.")

		if len(self.file_list) == 0:
			print("WARNING: No eiscat data files for this date.")
			self.valid = False

		if self.valid:
			self.files_list = [f for f in self.file_list if f[-4:] == ".mat"]


	def InitTime(self, bottle):
		self.start_datetime = bottle.DateTime(format="UT")
		self.end_datetime = bottle.DateTime("end", format="UT")

	def GetDataFromFile(self, file_name):
		fn = file_name[-12:-4]

		nb_seconds = int(fn)
		nb_seconds = dt.timedelta(seconds = nb_seconds)
		year = self.start_datetime.year

		datetime = dt.datetime(year, 1, 1) + nb_seconds

		data = False
		if self.start_datetime <= datetime < self.end_datetime:
			data = scipy.io.loadmat(file_name)

		return data

	def GetDataFromFiles(self):
		for f in self.files_list:
			data = self.GetDataFromFile(f)
			if data:
				self.data.append(data)

		if len(self.data) == 0:
			self.valid = False

	def GetParameter(self, parameter = "T_e", altitude = 110, time_divisor = 1):
		datetimes = []
		values = []
		errors = []

		tmp = abs(self.data[0]["r_h"] - altitude)
		alt_index = np.where(tmp == min(tmp))[0][0]

		if parameter == "N_e": param_index = 0
		elif parameter == "T_i": param_index = 1
		elif parameter == "T_e": param_index = 2
		elif parameter == "nu_in": param_index = 3
		elif parameter == "v_i": param_index = 4

		for d in self.data:
			time = dt.datetime(int(d["r_time"][0][0]), int(d["r_time"][0][1]), int(d["r_time"][0][2]),int( d["r_time"][0][3]), int(d["r_time"][0][4]), int(d["r_time"][0][5]))
			time = dt.timedelta(seconds = time.timestamp())
			datetimes.append(time.total_seconds() / time_divisor)


			values.append(d["r_param"][alt_index][param_index])
			errors.append(d["r_error"][alt_index][param_index])

			# if parameter == "T_e": #The value in the data correspond to T_e / T_i. Correcting for that
			# 	values[-1] *= d["r_param"][alt_index][1]


		zipped = zip(datetimes, values, errors)
		zipped = sorted(zipped, key=lambda x: x[0])
		to_del = []
		for i, (t, v, e) in enumerate(zipped):
			if np.isnan(t) or np.isnan(v) or np.isnan(e):
				to_del.insert(0, i)
		for id in to_del:
			del zipped[id]

		datetimes, values, errors = zip(*zipped)

		norm = datetimes[0]
		datetimes = [t - norm for t in datetimes]

		return datetimes, values, errors



class EiscatHDF5(Eiscat):

	def __init__(self, bottle, antenna = "tristatic"):
		self.antenna = antenna.lower()
		super().__init__(bottle, hdf5 = True)

	def SetPaths(self, bottle):

		self.folder = self.path #Init folder containing all Eiscat files for any dates

		self.folder += bottle.DateTime(format="UT").strftime("%Y%m%d") + "/"

		self.folders = []
		self.file_list = []


		one_day = dt.timedelta(days=1)
		start_date = bottle.DateTime("start", format="UT")
		end_date = bottle.DateTime("end", format="UT")
		while start_date.day <= end_date.day:
			folder = self.path
			folder += start_date.strftime("%Y%m%d") + "/"
			self.folders.append(folder)

			file_date = start_date.strftime("%Y-%m-%d")
			if   self.antenna == "tristatic":
				self.file_list.append(folder + f"MAD6529_{file_date}_beata_V139@kst.hdf5")

			elif self.antenna == "tromso":
				self.file_list.append(folder + f"MAD6400_{file_date}_beata_60@vhf.hdf5")

			elif self.antenna == "kiruna":
				self.file_list.append(folder + f"MAD6400_{file_date}_beata_V120@krva.hdf5")

			elif self.antenna == "sodankyla":
				self.file_list.append(folder + f"MAD6400_{file_date}_beata_V120@sdva.hdf5")

			elif self.antenna == "uhf_v":
				self.file_list.append(folder + f"MAD6502_{file_date}_beata_V225@uhf.hdf5")

			elif self.antenna == "uhf":
				self.file_list.append(folder + f"MAD6400_{file_date}_beata_ant@uhfa.hdf5")

			else:
				raise Exception("EISCAT antenna is incorrect. Please change it to one of the following: tristatic, tromso, kiruna, sodankila, uhf_v, uhf")

			start_date += one_day


		# print(self.file_list)

		if len(self.file_list) == 0:
			print("WARNING: No eiscat data HDF5 files for this date.")
			self.valid = False

		if self.valid:
			self.files_list = [f for f in self.file_list if f[-5:] == ".hdf5"]


	def GetDataFromFile(self, file_name):

		try:
			f = h5.File(file_name)
		except FileNotFoundError:
			print(f"ERROR: EISCAT data file not found: {file_name}")
			folder = '/'.join(file_name.split('/')[:-1])
			try:
				print(f"Available files in this directory: {os.listdir(folder)}")
			except FileNotFoundError:
				print(f"ERROR: EISCAT folder not found: {folder}")
				print(f"Available folders: {'/'.join(folder.split('/')[:-1])}")
			# self.valid = False
			return None


		d = np.array(f["Data/Table Layout"])

		data = []
		for i in range(len(d)):
			t = dt.datetime(d[i]["year"], d[i]["month"], d[i]["day"], d[i]["hour"], d[i]["min"], d[i]["sec"])
			if self.start_datetime <= t <= self.end_datetime:
				data.append(d[i])

		# print('EISCAT DATA COLUMNS', d[0].columns)
		return np.array(data)

	def GetDataFromFiles(self):
		self.data = np.array([])
		data_initialized = False

		for f in self.files_list:
			# print(f)
			data = self.GetDataFromFile(f)
			if data is None or data.size == 0:
				continue
			# print(data)
			if not data_initialized:
				# print("Init eiscat data")
				# print(self.data.shape)
				self.data = data
				data_initialized = True
				# print(self.data.shape)


			elif len(data[0]) == 37:
				# print(self.data.shape)
				# print(data.shape)
				self.data = np.append(self.data, data, axis = 0)

		if self.data is None:
			self.valid = False
		elif len(self.data) == 0:
			self.valid = False
		# print(self.data[0]["hour"], self.data[0]["min"], self.data[0]["gdalt"])
		# print(self.data[-1]["hour"], self.data[-1]["min"], self.data[-1]["gdalt"])
		# for i in range(len(self.data)):
		# 	print(self.data[i]["hour"], self.data[i]["min"], self.data[i]["gdalt"])


	def GetParameter(self, parameter = "ne", altitude = None, time_divisor = 1, time_format = "delta"):
		"""time_format= [delta, datetime]"""

		if altitude is not None:
			if altitude < 150:
				max_diff_altitude = altitude * 0.05 #km
			else:
				max_diff_altitude = altitude * 0.2 #km

		no_data = True
		while no_data:
			datetimes = []
			values = []
			errors = []
			for d in self.data:
				# print(altitude, d["gdalt"])

				tmpdt = dt.datetime(d["year"], d["month"], d["day"],d["hour"], d["min"], d["sec"])
				if (altitude is not None and abs(d["gdalt"] - altitude) >= max_diff_altitude) or (len(datetimes)>0 and tmpdt  == datetimes[-1]):
					# print("NO DATA :(")
					tmp = d["gdalt"]
					# print(f"NOT Taking data at altitude {tmp} ({altitude}+/-{max_diff_altitude}) at {tmpdt}")
					continue

				# tmp = d["gdalt"]
				# print(f"Taking data at altitude {tmp} (<{max_diff_altitude} of {altitude}) at {tmpdt}")
				# print("DATA!!!")
				datetimes.append(dt.datetime(d["year"], d["month"], d["day"],d["hour"], d["min"], d["sec"]))
				if time_format == "delta":
					datetimes[-1] = dt.timedelta(seconds = datetimes[-1].timestamp())
					datetimes[-1] = datetimes[-1].total_seconds() / time_divisor


				values.append(d[parameter])
				errors.append(d[f"d{parameter}"])

				# if parameter == "T_e": #The value in the data correspond to T_e / T_i. Correcting for that
				# 	values[-1] *= d["r_param"][alt_index][1]

			# If no data are found, maybe its because the constraint on the altitude of the parameter is too tight. So increase it until its meaningless and then give up.
			if len(values) < 1 and max_diff_altitude < altitude:
				max_diff_altitude *= 1.05
				continue
			elif len(values) >= 1 or max_diff_altitude >= altitude:
				no_data = False

		if not len(values):
			return [],[],[]
		else:
			print(f"EISCAT data taken at altitude {altitude} +/- {max_diff_altitude}.")

		zipped = zip(datetimes, values, errors)
		zipped = sorted(zipped, key=lambda x: x[0])
		# print("zipped", zipped)
		to_del = []
		for i, (t, v, e) in enumerate(zipped):
			if np.isnan(v) or np.isnan(e):
				to_del.insert(0, i)
		for id in to_del:
			del zipped[id]

		# print("zipped", zipped)

		datetimes, values, errors = zip(*zipped)

		if time_format == "delta":
			norm = datetimes[0]
			datetimes = [t - norm for t in datetimes]

		# print(datetimes, values, errors)

		return np.array(datetimes), np.array(values), np.array(errors)


	def GetUHFApparentAngle(self, bottle, time_divisor = 1, time_format = "delta"):

		altitude = bottle.GetAltitude()

		t, ve, dve = self.GetParameter(parameter = "vi1", altitude = altitude, time_divisor = time_divisor, time_format = time_format)
		t, vn, dvn = self.GetParameter(parameter = "vi2", altitude = altitude, time_divisor = time_divisor, time_format = time_format)
		t, vu, dvu = self.GetParameter(parameter = "vi3", altitude = altitude, time_divisor = time_divisor, time_format = time_format)

		nb_points = len(ve)

		# print(nb_points, ve)

		AoVi = np.zeros(nb_points)

		for i in range(nb_points):
			AoVi[i] = bottle.observation.GetApparentAngle([vu[i], ve[i], vn[i]]) + (np.pi/2)

		# print(f"Apparent angle of UHF Vi: {AoVi * RtoD}")

		AoVi = SetAngleBounds(AoVi, -np.pi/2 + np.pi/2 * bottle.graph_angle_shift, np.pi/2 + np.pi/2 * bottle.graph_angle_shift)

		print(f"Apparent angle of UHF Vi: {AoVi * RtoD}")

		return t, AoVi, None
