#!/usr/bin/python3
# -*-coding:utf-8 -*

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
#matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage{lmodern}']
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=False)
from utils import *
from rotation import *
from bottle import *
from scipy import signal
import sys
import os
from subprocess import call
import datetime as time


class EqCurrent:
	def __init__(self, bottle = None, file = None, file_type = None):
		print("Initialazing Equivalent current object")

		self.path = "/home/bossel/These/Analysis/data/equivalent_currents/"
		self.file = file
		self.files = []
		self.valid = True

		self.file_type = file_type ## Used when given a file and not a bottle. If file_type is None, then equivalent current file from magnar. If == "digisonde", then digisonde file

		if bottle is not None:
			self.files.append(bottle.DateTime("start", format="UT").strftime("%Y%m%d") + "_")
			if bottle.DateTime("start").day != bottle.DateTime("end").day:
				one_day = time.timedelta(days=1)
				start_date = bottle.DateTime("start", format="UT")
				end_date = bottle.DateTime("end", format="UT")
				while start_date.day <= end_date.day:
					self.files.append(start_date.strftime("%Y%m%d") + "_")
					start_date += one_day

			print("self.files", self.files)

			for i in range(len(self.files)):
				if self.file_type is None:
					# self.files[i] = bottle.DateTime("start", format="UT").strftime("%Y%m%d") + "_"

					if bottle.location.lower() in ["tromso", "skibotn", "skibotnsud", "skibotnnord", "kilpisjarvi"]:
						self.files[i] += "TRO"
					elif bottle.location.lower() in ["nyalesund", "corbel"]:
						self.files[i] += "NAL"
					else:
						self.valid = False

					self.files[i] += "_pola.txt"

				elif self.file_type == "digisonde":
					# self.files[i] = bottle.DateTime("start", format="UT").strftime("%Y%m%d") + "_"
					self.files[i] += "digisonde.txt"


			self.start_time = bottle.DateTime("start", format="UT")
			self.end_time = bottle.DateTime("end", format="UT")


			if self.files:
				self.file = self.files[0]
			else:
				self.valid = False
				self.file = [""]


			print("TIMES: ", self.start_time, self.end_time)

		elif self.file is not None:
			self.file = file
			self.files = [file]
			self.start_time = None
			self.end_time = None



		# print(self.file)

		# print(self.start_time, self.end_time)

		self.x_axis = []

		if self.valid:
			self.saving_name = self.file.split(".")[0]
			if self.file != "fake":
				try:
					self.LoadData()
				except:
					print("Warning: Current data not valid")
					self.valid = False
			else:
				self.LoadFakeData()


	def LoadFakeData(self):

		N = 1000

		seconds = np.linspace(0, N, N)
		now = time.datetime.now()
		datetime = np.array([now + time.timedelta(seconds=s) for s in seconds])


		Jn = np.cos(seconds * 2. * np.pi / N)
		Je = np.sin(seconds * 2. * np.pi / N)

		self.data = pd.DataFrame({	"seconds": seconds,
								"datetime": datetime,
								"Jn": Jn,
								"Je": Je})

		self.data["J_norm"] = np.sqrt(self.data["Jn"] ** 2 + self.data["Je"] ** 2)

	def LoadData(self):
		if not self.file_type:
			data = pd.DataFrame(columns=["date", "time", "Jn", "Je"])
			# print(data)
			# print(self.files)
			for f in self.files:
				try:
					data = data.append(pd.read_csv(self.path + f, delimiter=" ", names=["date", "time", "Jn", "Je"], skiprows=1))
				except:
					print("Equivalent current file not valid.")
					self.valid = False
					return 0
			# 	print(data)
			# print(data)
			data = data.reset_index(drop=True)
			# print(data)

			data["Ju"] = np.zeros(data.index.stop)
			data["datetime"] = [time.datetime.strptime(d + " " + t, "%Y-%m-%d %H:%M:%S") for d, t in zip(data["date"], data["time"])]
			# print(data)

		elif self.file_type == "digisonde":
			# print(self.path + self.file)
			data = pd.DataFrame(columns=["date", "time", "Jn", "Je", "Ju"])
			# print(data)
			for f in self.files:
				data = data.append(pd.read_fwf(self.path + f, header=0, usecols=(5, 7, 8, 10, 16), names=("date", "time", "Jn", "Je", "Ju")))
				# print(data)
			data = data.reset_index(drop=True)
			data["datetime"] = [time.datetime.strptime(d, "%Y.%m.%d%H:%M:%S") for d in data["date"] + data["time"]]
			# print(data)

		if self.start_time is None:
			print(data.index, data.index.start, data.index.stop-1)
			self.start_time = data["datetime"][data.index.start]
			self.end_time = data["datetime"][data.index.stop-1]

		# print(data)
		data["seconds"] = [(dt - self.start_time).total_seconds() for dt in data["datetime"]]
		# print(data)

		data["J_norm"] = np.sqrt(data["Jn"] ** 2 + data["Je"] ** 2 + data["Ju"] ** 2)
		# print(data)

		self.data = data

		# self.MakePlot(show=True, app_angle=False, div=3600)

		del data["date"]
		del data["time"]

		data = data[0 <= data["seconds"]]
		data = data[data["seconds"] <= (self.end_time - self.start_time).total_seconds()]

		data = data.reset_index(drop=True)

		self.data = data

		# self.MakePlot(show=True, app_angle=False, xaxis = True, div=3600)

	def GetApparentAngle(self, obs, shift=0):
		app_angle = []
		for i in self.data.index:
			app_angle.append(obs.GetApparentAngle([self.data["Ju"][i], self.data["Je"][i], self.data["Jn"][i]], Blos=True))

		app_angle, J_los_angle = zip(*app_angle)
		app_angle = np.array(app_angle)
		# print(J_los_angle)
		if not shift:
			self.data["AoJapp"] = app_angle
			self.data["AoJapp"] = SetAngleBounds(np.array(app_angle), -np.pi/2, np.pi/2)
		else:
			self.data["AoJapp"] = app_angle
			self.data["AoJapp"] = SetAngleBounds(np.array(app_angle), 0, np.pi)

		# self.data["AoJapp"]	+= np.pi/2
		self.data["AoJapp"] = SetAngleBounds(np.array(app_angle), -np.pi/2, np.pi/2)

		self.data["AoJlos"] = J_los_angle

	def FindJup(self, obs, AoLP_array):
		Ce, Se, Ca, Sa = obs.GetTrigo()

		tanA = np.tan(AoLP_array)

		Ju = lambda Je, Jn, tanA: (Je * (Se*Sa - Ca / tanA) + Jn * (Se*Ca + Sa / tanA)) / Ce
		Ju_list = Ju(self.data["Je"], self.data["Jn"], tanA)

		self.data["Ju"] = np.where(abs(Ju_list) < 12000, Ju_list, np.zeros_like(self.data["Je"]))
		# self.data["Ju"] = Ju_list

	def GetNormTimes(self, divisor = 1):
		norm = self.data["seconds"][0]
		self.x_axis = [(t - norm) / divisor for t in self.data["seconds"]]

		return self.x_axis

	def MakePlot(self, show=False, app_angle = True, xaxis=False, div=1):
		f1, axs1 =  plt.subplots(4, sharex=True, figsize=(16, 8))
		f2, axs2 =  plt.subplots(1, sharex=True, figsize=(16, 8))

		if not self.x_axis or xaxis:
			self.x_axis = self.GetNormTimes(divisor=div)

		axs1[0].plot(self.x_axis, self.data["Jn"])
		axs1[1].plot(self.x_axis, self.data["Je"])
		axs1[2].plot(self.x_axis, self.data["Ju"])
		axs1[3].plot(self.x_axis, self.data["J_norm"])
		if app_angle:
			# comp = (self.data["Jn"] / self.data["Je"]) * RtoD
			# axs2.plot(self.data["seconds"], comp, color = "blue", label = "AoJapp - AoJlos")

			axs2.plot(self.x_axis, self.data["AoJapp"] * RtoD, color = "black", label = "AoJapp")
			axs2.set_ylabel("AoJapp")
			axs22 = plt.twinx(axs2)
			axs22.plot(self.x_axis, self.data["AoJlos"] * RtoD, color = "red", label = "AoJlos")
			axs22.set_ylabel("AoJlos")


			f2.legend()

		axs1[0].set_ylabel("Jn")
		axs1[1].set_ylabel("Je")
		axs1[2].set_ylabel("Ju")
		axs1[3].set_ylabel("J norm")

		if show:
			plt.show()


		print("Saving currents in", self.path + "/" + self.saving_name + ".png")
		f1.savefig(self.path + "/" + self.saving_name + '_currents.png', bbox_inches='tight')
		f2.savefig(self.path + "/" + self.saving_name + '_angles.png', bbox_inches='tight')


	def SaveTXT(self):
		### Default Format
		print("Saving currents in", self.path + "/" + self.saving_name + ".txt")
		self.data.to_csv(self.path + "/" + self.saving_name + '_DONE.txt', sep=",", index=False)

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
			start_date = bottle.DateTime("start", format="UT")
			end_date = bottle.DateTime("end", format="UT")
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
			# print(f, len(self.data))
			self.data = np.concatenate((self.data, np.genfromtxt(f, dtype=self.array_type, delimiter = [12, 8, 11, 10, 10, 10, 10], skip_header=7, skip_footer=1, names=None)))
			# print(f, len(self.data))

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

		# print([(data[i+1], data[i], self.times_sec[i+1], self.times_sec[i]) for i in range(len(data)-1) if deriv[i] > 100])

		times = self.GetNormTimes(divisor)
		del times[-1]

		return times, deriv

	def StripTime(self, start, end):
		# start = usefull_times[0]
		# end = usefull_times[-1]
		# print("DEBUG MAGDATA TIME", start, end)
		# print("DEBUG MAGDATA TIME", self.datetime[0], self.datetime[-1])
		# print(len(self.data))

		self.data = np.array([d for i, d in enumerate(self.data) if start <= self.datetime[i] <= end], dtype = self.array_type)

		self.datetime = ([time.datetime.strptime(d + " " + t, "%d/%m/%Y %H:%M:%S") for d, t in zip(self.data["date"], self.data["time"])])

		self.times = ([time.timedelta(seconds = t.timestamp()) for t in self.datetime])
		self.times_sec = [t.total_seconds() for t in self.times]

		# print(len(self.data))


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
