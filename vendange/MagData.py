#!/usr/bin/python3
# -*-coding:utf-8 -*

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage{lmodern}']
from utils import *
from rotation import *
from bottle import *
from scipy import signal
import sys
import os
from subprocess import call
import datetime as time



class MagData:
	def __init__(self, bottle):

		print("Creating Magnetometer Data Object")

		self.data_path = "/home/bossel/These/Analysis/data/magnetometer/"

		self.file = self.data_path
		if bottle.location.lower() in ["tromso", "skibotn", "skibotnsud", "skibotnnord", "kilpisjarvi"]:
			self.file += "Tromso"
		elif bottle.location.lower() in ["nyalesund", "corbel"]:
			self.file += "Nyalesund"
		self.file += "/"

		self.additional_files = []
		if bottle.DateTime("start").day != bottle.DateTime("end").day:
			one_day = time.timedelta(days=1)
			start_date = bottle.DateTime("start")
			end_date = bottle.DateTime("end")
			while start_date.day != end_date.day:
				start_date += one_day
				self.additional_files.append(self.file + start_date.strftime("%Y%m%d"))

		self.file += bottle.DateTime().strftime("%Y%m%d")

		self.exist = True
		try:
			self.GetDataFromFile()
			print("INFO: Magnetic data available in file:", self.file)
		except:
			self.exist = False
			print("WARNING: No magnetometer data found for this observation")

	def GetDataFromFile(self):
		self.array_type = [('date','<U10'),('time','<U8'),('Dec',float),('Horiz',float),('Vert',float),('Incl',float),('Total',float)]
		self.data = np.genfromtxt(self.file, dtype=self.array_type, delimiter = [12, 8, 11, 10, 10, 10, 10], skip_header=7, skip_footer=1, names=None)

		for f in self.additional_files:
			print(f, len(self.data))
			self.data = np.concatenate((self.data, np.genfromtxt(f, dtype=self.array_type, delimiter = [12, 8, 11, 10, 10, 10, 10], skip_header=7, skip_footer=1, names=None)))
			print(f, len(self.data))

		self.datetime = ([time.datetime.strptime(d + " " + t, "%d/%m/%Y %H:%M:%S") for d, t in zip(self.data["date"], self.data["time"])])

		self.times = np.array([t.timestamp() for t in self.datetime]) # number of seconds since 01/01/1970
		# self.times = np.array([time.timedelta(seconds = t.timestamp()) for t in self.datetime]) # number of seconds since 01/01/1970
		# self.times_sec = [t.total_seconds() for t in self.times]


	def GetComponent(self, comp, divisor):
		# data = self.data["Dec"]
		# deriv = [(data[i+1] - data[i])**2 for i in range(len(data)-1)]
		# data = self.data["Horiz"]
		# deriv = [deriv[i] + (data[i+1] - data[i])**2 for i in range(len(data)-1)]
		# data = self.data["Vert"]
		# deriv = [deriv[i] + (data[i+1] - data[i])**2 for i in range(len(data)-1)]
		# deriv = [np.sqrt(d) for d in deriv]

		data = self.data[comp]
		times = self.GetNormTimes(divisor)

		return times, data

	def GetDerivative(self, comp, divisor):
		# data = self.data["Dec"]
		# deriv = [(data[i+1] - data[i])**2 for i in range(len(data)-1)]
		# data = self.data["Horiz"]
		# deriv = [deriv[i] + (data[i+1] - data[i])**2 for i in range(len(data)-1)]
		# data = self.data["Vert"]
		# deriv = [deriv[i] + (data[i+1] - data[i])**2 for i in range(len(data)-1)]
		# deriv = [np.sqrt(d) for d in deriv]

		data = self.data[comp]

		deriv = []
		for i in range(len(data)-1):
			if data[i+1] < 99999.9 and data[i] < 99999.9:
				deriv.append((data[i+1] - data[i])/(self.times_sec[i+1] - self.times_sec[i]))
			else:
				deriv.append(0)

		print([(data[i+1], data[i], self.times_sec[i+1], self.times_sec[i]) for i in range(len(data)-1) if deriv[i] > 100])

		times = self.GetNormTimes(divisor)
		del times[-1]

		return times, deriv

	def StripTime(self, start, end):
		# start = usefull_times[0]
		# end = usefull_times[-1]
		print("DEBUG MAGDATA TIME", start, end)
		print("DEBUG MAGDATA TIME", self.datetime[0], self.datetime[-1])
		print(len(self.data))

		self.data = np.array([d for i, d in enumerate(self.data) if start <= self.datetime[i] <= end], dtype = self.array_type)

		self.datetime = ([time.datetime.strptime(d + " " + t, "%d/%m/%Y %H:%M:%S") for d, t in zip(self.data["date"], self.data["time"])])

		self.times = ([time.timedelta(seconds = t.timestamp()) for t in self.datetime])
		self.times_sec = [t.total_seconds() for t in self.times]

		print(len(self.data))


	def GetNormTimes(self, divisor = 1):
		try:
			norm = self.times_sec[0]
		except:
			print("WARNING: GeoMag Data has no time list.")
			norm = 0.

		return [(t - norm) / divisor for t in self.times_sec]



	def __getattr__(self, name):
		if name in ["Dec", "Horiz", "Vert", "Total", 'Incl']:
			return self.data[name]

	def __repr__(self):
		return "Magnetic data object."
