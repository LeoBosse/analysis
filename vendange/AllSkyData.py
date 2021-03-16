#!/usr/bin/python3
# -*-coding:utf-8 -*

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage{lmodern}']
from utils import *
from rotation import *
from bottle import *
from scipy import signal, io, misc
import sys
import os
from subprocess import call
import datetime as dtm

import imageio
# import skimage as ski



class AllSkyData:
	def __init__(self, bottle, filename=None):

		print("Creating AllSky Imager Data Object")

		self.data_path = "/home/bossel/These/Analysis/data/allsky_camera/tid.uio.no/plasma/aurora/"

		no_data = False

		self.file = self.data_path
		if filename:
			self.file += filename
		else:
			if   bottle.location.lower() == "nyalesund": self.location = "nya6"
			elif bottle.location.lower() == "skibotn": self.location = "skn4"
			else: self.location = "nan"

			if   bottle.filters and bottle.filters[0] in ["v", "m", "b"]: self.color = "5577"
			elif bottle.filters and bottle.filters[0] == "r": self.color = "6300"
			elif bottle.instrument_name == "spp": self.color = "6300"
			else: self.color = "nan"

			self.file += self.location + "/"
			self.file += self.color + "/"

		self.azimut = bottle.azimut * RtoD #Calib file coord are given in degrees!
		self.elevation = bottle.elevation * RtoD #Calib file coord are given in degrees!
		self.deg_delta = 1

		self.start = bottle.DateTime("start", format="UT")
		self.end   = bottle.DateTime("end", format="UT")
		print("DEBUG ALL SKY:", self.start, self.end)
		self.main_datetimes = []
		self.all_datetimes = []

		self.GetMainDateTimes()
		self.CleanMainDateTimes()

		self.GetAllDateTimes()
		self.CleanAllDateTimes()
		# print(self.all_datetimes)

		if self.all_datetimes:
			self.brightness = self.GetBrightnessList()

	def GetCalibrationFile(self, datetime):
		file = self.file
		file += str(datetime.year) + "/"
		file += datetime.strftime("%Y%m%d") + "/"
		file += "_".join((self.location, datetime.strftime("%Y%m%d"), self.color, "cal.dat"))

		return io.readsav(file)

	def GetMainDateTimes(self):
		t = self.start
		one_hour = dtm.timedelta(hours = 1)
		while t <= self.end:
			self.main_datetimes.append(t)
			t += one_hour
		if self.main_datetimes[-1].hour != self.end.hour:
			self.main_datetimes.append(self.main_datetimes[-1] + one_hour)

	def GetFolderName(self, datetime):
		file = self.file
		file += str(datetime.year) + "/"
		file += datetime.strftime("%Y%m%d") + "/"
		file += "ut"
		h = str(datetime.hour)
		if len(h) == 1:
			h = "0" + h
		file += h + "/"

		return file


	def GetDateTimeFromFile(self, file):
		file = file.split("/")[-1]
		datetime = file.split("_")[1:3]
		datetime = dtm.datetime.strptime("_".join(datetime), "%Y%m%d_%H%M%S")

		return datetime

	def GetAllDateTimes(self):
		# print(len(self.main_datetimes))
		for i, dt in enumerate(self.main_datetimes):
			folder = self.GetFolderName(dt)
			files_list = os.listdir(folder)
			files_list = [f for f in files_list if f[-4:] == ".png"]
			files_list = [f for f in files_list if len(f) == 33]

			for f in files_list:
				self.all_datetimes.append(self.GetDateTimeFromFile(f))

		self.all_datetimes.sort()


	def GetFileName(self, datetime):
		file = self.GetFolderName(datetime)
		file += "_".join((self.location, datetime.strftime("%Y%m%d_%H%M%S"), self.color, "cal.png"))
		return file

	def CleanMainDateTimes(self):
		self.main_datetimes = [dt for dt in self.main_datetimes if os.path.exists(self.GetFolderName(dt))]

	def CleanAllDateTimes(self):
		self.all_datetimes = [dt for dt in self.all_datetimes if self.start <= dt < self.end]


	def SetPixelArea(self, image, calibration):
		"""List of all pixels in the line of sight of the instrument"""
		self.pixel_area = []

		for i, lign in enumerate(image):
			for j, col in enumerate(lign):
				pix_azimut    = calibration["mazms"][i][j] #IN DEGREES
				pix_elevation = calibration["elevs"][i][j] #IN DEGREES
				if not np.isnan(pix_azimut) and not np.isnan(pix_elevation):
					# print(pix_azimut, pix_elevation)
					# print(self.deg_delta, self.azimut, self.elevation)

					if pix_azimut - self.deg_delta < self.azimut and self.azimut < pix_azimut + self.deg_delta and pix_elevation - self.deg_delta < self.elevation and self.elevation < pix_elevation + self.deg_delta:
						self.pixel_area.append((i, j))

		self.nb_pixels = len(self.pixel_area)


	def GetDataFromFile(self, file):
		self.image = imageio.imread(file)

		brightness = 0

		for pix in self.pixel_area:
			i, j = pix
			brightness += self.image[i][j]  #/ self.nb_pixels

		return brightness

	def GetBrightnessList(self):
		brightness_list = []

		datetime = self.all_datetimes[0]
		image_name = self.GetFileName(datetime)
		image = imageio.imread(image_name)
		calibration = self.GetCalibrationFile(datetime)
		# print("calibration", calibration)
		self.SetPixelArea(image, calibration)

		for dt in self.all_datetimes:
			if datetime.day != dt.day:
				datetime = dt
				image_name = self.GetFileName(datetime)
				image = imageio.imread(image_name)
				calibration = self.GetCalibrationFile(datetime)
				self.SetPixelArea(image, calibration)

			brightness_list.append(self.GetDataFromFile(self.GetFileName(dt)))

		return brightness_list


	def GetNormTimes(self, zero_datetime, divisor = 1):
		norm = zero_datetime
		timedelta_list = [dt - norm for dt in self.all_datetimes]

		return [t.total_seconds() / divisor for t in timedelta_list]
