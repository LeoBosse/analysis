#!/usr/bin/python3
# -*-coding:utf-8 -*

import numpy as np
import matplotlib.pyplot as plt
from utils import *
from scipy import signal
import sys as sys
import os
from subprocess import call
import observation as observation
import observation as geometry

import datetime as time
import scipy.io

### Bottle classes.

class Eiscat:
	def __init__(self, bottle):

		self.valid = True
		self.path = "/home/bossel/These/Analysis/data/eiscat/"

		self.data = []

		self.start_datetime = bottle.DateTime(format="UT")
		self.end_datetime = bottle.DateTime("end", format="UT")

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
			one_day = time.timedelta(days=1)
			start_date = bottle.DateTime("start", format="UT")
			end_date = bottle.DateTime("end", format="UT")
			while start_date.day != end_date.day:
				start_date += one_day
				folder = self.path
				folder += start_date.strftime("%Y-%m-%d")
				folder += "_beata_10@uhfa/"
				self.folders.append(folder)
		print(self.folders)
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
			self.GetDataFromFiles()
			if len(self.data) == 0:
				self.valid = False
				print("No EISCAT data available.")
			else:
				print("EISCAT data exist! Number of points:", len(self.data))
				print(bottle.DateTime(format="UT").strftime("%Y-%m-%d"))


	def GetDataFromFile(self, file_name):
		fn = file_name[-12:-4]

		nb_seconds = int(fn)
		nb_seconds = time.timedelta(seconds = nb_seconds)
		year = self.start_datetime.year

		datetime = time.datetime(year, 1, 1) + nb_seconds

		data = False
		if self.start_datetime <= datetime < self.end_datetime:
			data = scipy.io.loadmat(file_name)

		return data

	def GetDataFromFiles(self):
		for f in self.files_list:
			data = self.GetDataFromFile(f)
			if data:
				self.data.append(data)

	def GetParameter(self, parameter = "T_e", altitude = 110, time_divisor = 1):
		datetimes = []
		values = []
		errors = []

		tmp = abs(self.data[0]["r_h"] - 110)
		alt_index = np.where(tmp == min(tmp))[0][0]

		if parameter == "N_e": param_index = 0
		elif parameter == "T_i": param_index = 1
		elif parameter == "T_e": param_index = 2
		elif parameter == "nu_in": param_index = 3
		elif parameter == "v_i": param_index = 4

		for d in self.data:
			dt = time.datetime(int(d["r_time"][0][0]), int(d["r_time"][0][1]), int(d["r_time"][0][2]),int( d["r_time"][0][3]), int(d["r_time"][0][4]), int(d["r_time"][0][5]))
			dt = time.timedelta(seconds = dt.timestamp())
			datetimes.append(dt.total_seconds() / time_divisor)


			values.append(d["r_param"][alt_index][param_index])
			errors.append(d["r_error"][alt_index][param_index])

			# if parameter == "T_e": #The value in the data correspond to T_e / T_i. Correcting for that
			# 	values[-1] *= d["r_param"][alt_index][1]


		zipped = zip(datetimes, values, errors)
		zipped = sorted(zipped, key=lambda x: x[0])
		to_del = []
		for i, (dt, v, e) in enumerate(zipped):
			if np.isnan(dt) or np.isnan(v) or np.isnan(e):
				to_del.insert(0, i)
		for id in to_del:
			del zipped[id]

		datetimes, values, errors = zip(*zipped)

		norm = datetimes[0]
		datetimes = [dt - norm for dt in datetimes]

		return datetimes, values, errors
