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
from scipy import signal
import sys
import os
import wget
from subprocess import call
import datetime as dt

from utils import *
from rotation import *
from bottle import *
from vendange_configuration import *


class EqCurrent:
	def __init__(self, bottle = None, file = None, file_type = None):
		print("Initialazing Equivalent current object")

		self.path = global_configuration.eq_current_data_path
		# self.path = "/home/bossel/These/Analysis/data/equivalent_currents/"
		self.file = file
		self.files = []
		self.valid = True

		self.file_type = file_type ## Used when given a file and not a bottle. If file_type is None, then equivalent current file from magnar. If == "digisonde", then digisonde file

		self.SetFileNames(bottle=bottle)

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

	def SetFileNames(self, bottle=None):

		date_format = "%Y%m%d_"
		if bottle.DateTime("start", format="UT").year >= 2022:
			date_format = "%Y-%m-%d_"

		if bottle is not None:
			self.files.append(bottle.DateTime("start", format="UT").strftime(date_format))
			if bottle.DateTime("start").day != bottle.DateTime("end").day:
				one_day = dt.timedelta(days=1)
				start_date = bottle.DateTime("start", format="UT")
				end_date = bottle.DateTime("end", format="UT")
				start_date += one_day
				while start_date.day <= end_date.day:
					self.files.append(start_date.strftime(date_format))
					start_date += one_day

			# print("self.files", self.files)

			for i in range(len(self.files)):
				if self.file_type is None:
					# self.files[i] = bottle.DateTime("start", format="UT").strftime("%Y%m%d") + "_"
					if bottle.location.lower() in ["tromso", "skibotn", "skibotnsud", "skibotnnord", "kilpisjarvi"]:
						self.files[i] += "TRO_pola"
					elif bottle.location.lower() in ["nyalesund", "corbel"]:
						self.files[i] += "NAL_pola"
					elif "finland" in bottle.location.lower():
						self.files[i] += "Location_"

						if int(bottle.azimut*RtoD)%360 == -42%360:
							self.files[i] += "0"
						elif int(bottle.azimut*RtoD)%360 == -162%360:
								self.files[i] += "1"
						elif int(bottle.azimut*RtoD)%360 == 48%360:
								self.files[i] += "4"
						elif int(bottle.azimut*RtoD)%360 == 210%360:
								self.files[i] += "5"
						elif int(bottle.azimut*RtoD)%360 == 260%360:
								self.files[i] += "6"
						else:
							print(f"Warning, azimut {bottle.azimut*RtoD} not recognized for equivalent current in finland")
							self.valid = False
					else:
						self.valid = False

					self.files[i] += ".txt"

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


			# print("TIMES: ", self.start_time, self.end_time)

		elif self.file is not None:
			self.file = file
			self.files = [file]
			self.start_time = None
			self.end_time = None

	def LoadFakeData(self):

		N = 1000

		seconds = np.linspace(0, N, N)
		now = dt.datetime.now()
		datetime = np.array([now + dt.timedelta(seconds=s) for s in seconds])


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
					print(f"Equivalent current file valid. {self.path + f}")
				except:
					print(f"Equivalent current file not valid. {self.path + f}")
					# return 0
			if data.empty:
				print("No equivalent current data retrieved from files")
				self.valid = False
			# 	print(data)
			# print(data)
			data = data.reset_index(drop=True)
			# print(data)

			data["Ju"] = np.zeros(data.index.stop)
			data["datetime"] = [dt.datetime.strptime(d + " " + t, "%Y-%m-%d %H:%M:%S") for d, t in zip(data["date"], data["time"])]
			# print(data)

		elif self.file_type == "digisonde":
			# print(self.path + self.file)
			data = pd.DataFrame(columns=["date", "time", "Jn", "Je", "Ju"])
			# print(data)
			for f in self.files:
				data = data.append(pd.read_fwf(self.path + f, header=0, usecols=(5, 7, 8, 10, 16), names=("date", "time", "Jn", "Je", "Ju")))
				# print(data)
			data = data.reset_index(drop=True)
			data["datetime"] = [dt.datetime.strptime(d, "%Y.%m.%d%H:%M:%S") for d in data["date"] + data["time"]]
			# print(data)

		if self.start_time is None:
			# print(data.index, data.index.start, data.index.stop-1)
			self.start_time = data["datetime"][data.index.start]
			self.end_time = data["datetime"][data.index.stop-1]

		# print(data)
		data["seconds"] = [(dt - self.start_time).total_seconds() for dt in data["datetime"]]
		# print(data)

		data["J_norm"] = np.sqrt(data["Jn"] ** 2 + data["Je"] ** 2 + data["Ju"] ** 2)
		data["Jaz"] = np.arctan2(data["Je"], data["Jn"])
		data["Jel"] = np.arcsin(data["Ju"] / data["J_norm"])
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


	def Polarisation(self, obs, DoLP_max = 100, DoLP_mode="rayleigh", AoLP_mode="para"):
		DoLPs, AoLPs = [], []
		for ju, je, jn in zip(self.data["Ju"], self.data["Je"], self.data["Jn"]):
			current = [ju, je, jn]
			DoLP, AoLP = EqCurrent.GetPolarisation(current, obs, DoLP_max=DoLP_max, DoLP_mode=DoLP_mode, AoLP_mode=AoLP_mode)

			DoLPs.append(DoLP)
			AoLPs.append(AoLP)

		self.data["DoLP"] = DoLPs
		self.data["AoLP"] = AoLPs



	@staticmethod
	def GetPolarisation(current, obs, DoLP_max = 100, DoLP_mode="rayleigh", AoLP_mode="para"):
		"""From a current (up-east-north array), and an observation object returns the DoLP (0-1) and AoLP (radians) of the polarization.
		Polarization modes define DoLP and AoLP.
		DoLP_mode: 	rayleigh: DoLP = sin**2 / (1+cos**2) of the "scattering angle", i.e. the angle between the current and the line of sight
					sin, sin2, cos, cos2: use a simple trigo function (squared if 2) on the "scattering angle"
		AoLP_mode: 	para= The direction of pola is parallel to the current
					perp= The direction of pola is perpendicular to the current and the line of sight (approx. like rayleigh scattering)
		DoLP_max: in % (1-100) The maximum DoLP of the induced polarisation. The DoLP_mode function is multiplied by a factor DoLP_max/100.
		"""

		AoLP_mode = AoLP_mode.lower()
		DoLP_mode = DoLP_mode.lower()

		AoLP = None
		DoLP = None

		if AoLP_mode in ["para", "parallel"]:
			AoLP, AoBlos = obs.GetApparentAngle(current, Blos=True)
		elif AoLP_mode in ["perp", "perpendicular", "rayleigh"]:
			los = [obs.los_uen[0,0], obs.los_uen[1,0], obs.los_uen[2,0]]
			pola_direction = np.cross(los, current)

			AoLP, AoBlos = obs.GetApparentAngle(pola_direction, Blos=True)

		if DoLP_mode in ["rayleigh"]:
			DoLP = np.sin(AoBlos)**2 / (1 + np.cos(AoBlos)**2)
		elif DoLP_mode in ["sin"]:
			DoLP = np.sin(AoBlos)
		elif DoLP_mode in ["sin2"]:
			DoLP = np.sin(AoBlos)**2
		elif DoLP_mode in ["cos"]:
			DoLP = np.cos(AoBlos)
		elif DoLP_mode in ["cos2"]:
			DoLP = np.cos(AoBlos)**2

		DoLP *= DoLP_max / 100.

		return DoLP, AoLP #[0-1] and radians


	def GetApparentAngle(self, obs, shift=0, Jup_mode=False):
		self.data["AoJapp"], self.data["AoJlos"] = EqCurrent.ApparentAngle(self.data["Ju"], self.data["Je"], self.data["Jn"],  obs, shift)

		if Jup_mode == "perp":
			self.data["AoJapp"] += np.pi/2
			if not shift:
				# current.data["AoJapp"] = app_angle
				app_angle = SetAngleBounds(self.data["AoJapp"], -np.pi/2, np.pi/2)
			else:
				# current.data["AoJapp"] = app_angle
				app_angle = SetAngleBounds(self.data["AoJapp"], 0, np.pi)


	@staticmethod
	def ApparentAngle(cu, ce, cn, obs, shift=0):
		app_angle = []
		for i in range(len(cu)):
			app_angle.append(obs.GetApparentAngle([cu[i], ce[i], cn[i]], Blos=True))



		# print(app_angle)
		app_angle, J_los_angle = zip(*app_angle)
		app_angle = np.array(app_angle) #+ (np.pi/2)
		# print(J_los_angle)
		if not shift:
			# current.data["AoJapp"] = app_angle
			app_angle = SetAngleBounds(np.array(app_angle), -np.pi/2, np.pi/2)
		else:
			# current.data["AoJapp"] = app_angle
			app_angle = SetAngleBounds(np.array(app_angle), 0, np.pi)

		# current.data["AoJapp"]	+= np.pi/2
		app_angle = SetAngleBounds(np.array(app_angle), -np.pi/2, np.pi/2)

		# current.data["AoJlos"] = J_los_angle

		return app_angle, J_los_angle


	def FindJup(self, obs, AoLP_array, mode = "para"):
		"""
		From an observation object (pointing direction), a horizontal current and an AoLP, computes the vertical current needed so that the apparent angle of the current is aligned on the AoLP.
		The vertical component is directly added to the EqCurrent object data array. Update the norm and the elevation of the current.
		"""

		print("Computing vertical current...")

		Ce, Se, Ca, Sa = obs.GetTrigo()

		tanA = np.tan(AoLP_array)

		### Transform the current vector from the up-east-north coordinates at the emission (H) to the uen coordinates at the instrument in A, then in xyz of the instrument.
		R_HA = obs.GetTransormationMatrixHA()
		R_AI = obs.GetRotMatrixAI()
		R_HI = R_AI @ R_HA

		### Transforms the current vector in H to the instrument coordiates
		J_H = np.array((np.zeros(len(self.data["Je"])), self.data["Je"], self.data["Jn"]))
		J_I = R_HI @ J_H
		Jx, Jy, Jz = J_I[0], J_I[1], J_I[2]

		### Transforms the up normal vector (1,0,0) in H to the instrument coordiates
		Up_H = np.array([[1], [0], [0]])
		Up_I = R_HI @ Up_H
		Ux, Uy, Uz = Up_I[0], Up_I[1], Up_I[2]


		if mode.lower() == "para":
		### We want to find x in J + x * Up. x is the same in H or I coordinat system.
		### We want that tan(AoJ) = tan(AoLP) = tanA = (Jy + x * Upy) / ((Jz + x * Upz)) (in the I coordinates)
		### So that x = (tanA * Jz - Jy) / (Uy - tanA * Uz)
			Findx = lambda Jy, Jz, Uy, Uz, tanA: (tanA * Jz - Jy) / (Uy - tanA * Uz)
		elif mode.lower() == "perp":
			Findx = lambda Jy, Jz, Uy, Uz, tanA: - (tanA * Jy + Jz) / (Uz + tanA * Uy)
		else:
			raise Exception("Error: Incorrect mode to compute Jup for equivalent currents. Correct strings are: perp or para")

		### Now in the Reference of H, the the new J is just (x, Je, Jn)
		Ju_list = Findx(Jy, Jz, Uy, Uz, tanA)


		### Old formula when considering that the A and H coordinates are equal (instrument points to the zenith for example)
		# Ju = lambda Je, Jn, tanA: (Je * (Se*Sa - Ca / tanA) + Jn * (Se*Ca + Sa / tanA)) / Ce
		# Ju_list = Ju(J_Ae, J_An, tanA)

		max_Ju = 2 * max(max(abs(self.data["Je"])), max(abs(self.data["Jn"])))

		self.data["Ju"] = np.where(abs(Ju_list) < max_Ju, Ju_list, np.zeros_like(self.data["Je"]))
		# self.data["Ju"] = Ju_list

		self.data["J_norm"] = np.sqrt(self.data["Jn"] ** 2 + self.data["Je"] ** 2 + self.data["Ju"] ** 2)
		self.data["Jel"] = np.arcsin(self.data["Ju"] / self.data["J_norm"])

	def GetNormTimes(self, divisor = 1, format = "delta"):
		if format == "delta":
			norm = self.data["seconds"][0]
			self.x_axis = np.array([(t - norm) / divisor for t in self.data["seconds"]])
		else:
			self.x_axis = self.data["datetime"]

		return self.x_axis

	def GetInterpolation(self, new_time, obs, divisor=1, shift=0):

		new_data = pd.DataFrame()

		new_data["seconds"] = new_time

		old_time = self.GetNormTimes(divisor = divisor)
		# print("new_time", new_time)
		# print("old_time", old_time)
		new_data["Ju"] = np.interp(new_time, old_time, self.data["Ju"])
		new_data["Je"] = np.interp(new_time, old_time, self.data["Je"])
		new_data["Jn"] = np.interp(new_time, old_time, self.data["Jn"])

		new_data["J_norm"] = np.sqrt(new_data["Jn"] ** 2 + new_data["Je"] ** 2 + new_data["Ju"] ** 2)
		new_data["Jaz"] = np.arctan2(new_data["Je"], new_data["Jn"])
		new_data["Jel"] = np.arcsin(new_data["Ju"] / new_data["J_norm"])

		new_data["AoJapp"], new_data["AoJlos"] = EqCurrent.ApparentAngle(new_data["Ju"], new_data["Je"], new_data["Jn"], obs, shift)

		return new_data

	def MakePlot(self, show = False, app_angle = True, xaxis = False, div = 1, coords = "uen"):
		f1, axs1 =  plt.subplots(3, sharex=True, figsize=(16, 8))
		f2, axs2 =  plt.subplots(1, sharex=True, figsize=(16, 8))

		# if not self.x_axis or xaxis:
		# 	self.x_axis = self.GetNormTimes(divisor=div)

		if coords == "uen":
			axs1[0].plot(self.x_axis, self.data["Jn"], "*")
			axs1[1].plot(self.x_axis, self.data["Je"], "*")
			axs1[2].plot(self.x_axis, self.data["Ju"], "*")
		elif coords == "azel":
			axs1[0].plot(self.x_axis, self.data["Jaz"] * RtoD, "*")
			axs1[1].plot(self.x_axis, self.data["Jel"] * RtoD, "*")
			axs1[2].plot(self.x_axis, self.data["J_norm"], "*")

		if app_angle:
			# comp = (self.data["Jn"] / self.data["Je"]) * RtoD
			# axs2.plot(self.data["seconds"], comp, color = "blue", label = "AoJapp - AoJlos")

			axs2.plot(self.x_axis, self.data["AoJapp"] * RtoD, color = "black", label = "AoJapp")
			axs2.set_ylabel("AoJapp")
			axs22 = plt.twinx(axs2)
			axs22.plot(self.x_axis, self.data["AoJlos"] * RtoD, color = "red", label = "AoJlos")
			axs22.set_ylabel("AoJlos")


			f2.legend()

		if coords == "uen":
			axs1[0].set_ylabel("Jn")
			axs1[1].set_ylabel("Je")
			axs1[2].set_ylabel("Ju")
		elif coords == "azel":
			axs1[0].set_ylabel("Jaz")
			axs1[1].set_ylabel("Jel")
			axs1[2].set_ylabel("J norm")

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
	def __init__(self, location, start, end, force_download = True):

		print("Creating Magnetometer Data Object")
		self.force_download = force_download
		self.location, self.start, self.end = location, start, end

		self.data_path = global_configuration.magnetometer_data_path
		# self.data_path = "/home/bossel/These/Analysis/data/magnetometer/"

		self.file = self.data_path
		if location.lower() in ["tromso", 'lacsud', "skibotn", "skibotnsud", "skibotnnord", "kilpisjarvi"]:
			self.file += "Tromso"
		elif location.lower() in ["nyalesund", "corbel"]:
			self.file += "Nyalesund"
		self.file += "/"


		self.additional_files = []
		if start.day != end.day:
			one_day = dt.timedelta(days=1)
			start_date = start
			end_date = end
			while start_date.day != end_date.day:
				start_date += one_day
				self.additional_files.append(self.file + start_date.strftime("%Y%m%d"))
				


		self.file += start.strftime("%Y%m%d")
		self.exist = True

		if self.force_download:
			self.DownloadFiles()

		try:
			self.GetDataFromFile()
			self.StripTime()
			print("INFO: Magnetic data available in file:", self.file)
		except:
			self.exist = False
			print(f"WARNING: No magnetometer data found for this observation ({self.file})")


	@classmethod
	def FromBottle(cls, bottle, force_download = True):
		return MagData(bottle.location, bottle.DateTime("start", format="UT"), bottle.DateTime("end", format="UT"))

	@classmethod
	def FromSpectroSerie(cls, serie, force_download = True):
		return MagData('skibotn', serie.times[0], serie.times[-1], force_download = force_download)


	def DownloadFiles(self):
		for i, f in enumerate([self.file] + self.additional_files):
			if not self.CheckDataFileExists(f):
				self.DownloadDataFile(f)

	def CheckDataFileExists(self, filename):
		if os.path.isfile(filename):
			return True
		else:
			return False
	
	def DownloadDataFile(self, save_path):
		if 'Tromso' in save_path:
			site = 'tro2a'
		elif 'Nyalesund' in save_path:
			site = 'nal1a'
		date = dt.datetime.strptime(save_path.split('/')[-1], "%Y%m%d")
		year  = date.year
		month = date.month
		day   = date.day
		password = "ResUseNoCom"
		res = "10sec"

		url = f"http://flux.phys.uit.no/cgi-bin/mkascii.cgi?site={site}&year={year}&month={month}&day={day}&res={res}&pwd={password}&format=html&comps=DHZ&getdata=+Get+Data+"
		
		print(f"Downloading magnetometer data files from {url} in {self.data_path}.")
		wget.download(url, save_path)


	def GetDataFromFile(self):
		self.array_type = [('date','<U10'),('time','<U8'),('Dec',float),('Horiz',float),('Vert',float),('Incl',float),('Total',float)]

		self.data = np.genfromtxt(self.file, dtype=self.array_type, delimiter = [12, 8, 11, 10, 10, 10, 10], skip_header=7, skip_footer=1, names=None)

		for f in self.additional_files:
			print(f, len(self.data))
			self.data = np.concatenate((self.data, np.genfromtxt(f, dtype=self.array_type, delimiter = [12, 8, 11, 10, 10, 10, 10], skip_header=7, skip_footer=1, names=None)))
			print(f, len(self.data))

		self.datetime = ([dt.datetime.strptime(d + " " + t, "%d/%m/%Y %H:%M:%S") for d, t in zip(self.data["date"], self.data["time"])])

		self.times = np.array([t.timestamp() for t in self.datetime]) # number of seconds since 01/01/1970
		# self.times = np.array([dt.timedelta(seconds = t.timestamp()) for t in self.datetime]) # number of seconds since 01/01/1970
		# self.times_sec = [t.total_seconds() for t in self.times]


	def GetComponent(self, comp, divisor, use_datetime=True):
		# data = self.data["Dec"]
		# deriv = [(data[i+1] - data[i])**2 for i in range(len(data)-1)]
		# data = self.data["Horiz"]
		# deriv = [deriv[i] + (data[i+1] - data[i])**2 for i in range(len(data)-1)]
		# data = self.data["Vert"]
		# deriv = [deriv[i] + (data[i+1] - data[i])**2 for i in range(len(data)-1)]
		# deriv = [np.sqrt(d) for d in deriv]

		data = self.data[comp]
		times = self.GetNormTimes(divisor, use_datetime = use_datetime)

		return times, data

	def GetDerivative(self, comp, divisor, use_datetime=True):
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

		times = self.GetNormTimes(divisor, use_datetime=use_datetime)
		del times[-1]

		return times, deriv

	def StripTime(self):
		# start = usefull_times[0]
		# end = usefull_times[-1]
		# print("DEBUG MAGDATA TIME", start, end)
		# print("DEBUG MAGDATA TIME", self.datetime[0], self.datetime[-1])
		# print(len(self.data))

		self.data = np.array([d for i, d in enumerate(self.data) if self.start <= self.datetime[i] <= self.end], dtype = self.array_type)

		self.datetime = ([dt.datetime.strptime(d + " " + t, "%d/%m/%Y %H:%M:%S") for d, t in zip(self.data["date"], self.data["time"])])

		self.times = ([dt.timedelta(seconds = t.timestamp()) for t in self.datetime])
		self.times_sec = [t.total_seconds() for t in self.times]

		# print(len(self.data))


	def GetNormTimes(self, divisor = 1, use_datetime = True):
		if use_datetime:
			return self.datetime

		try:
			norm = self.times_sec[0]
		except:
			print("WARNING: GeoMag Data has no time list.")
			norm = 0.

		return [(t - norm) / divisor for t in self.times_sec]


	def MakeFigure(self):
		fig, axs = plt.subplots(3, 1, sharex=True)
		axs = np.append(axs, [axs[0].twinx(), axs[1].twinx()])
		axs[0].plot(self.datetime, self.Dec, 'r', label='Declination')
		axs[3].plot(self.datetime, self.Incl, 'b', label='Inclination')
		axs[1].plot(self.datetime, self.Horiz, 'r', label='Horizontal')
		axs[4].plot(self.datetime, self.Vert, 'b', label='Vertical')
		axs[2].plot(self.datetime, self.Total, 'k', label = 'Total')
		
		lines, labels = np.empty_like(axs), np.empty_like(axs)
		for i, a in enumerate(axs):
			lines[i], labels[i] = a.get_legend_handles_labels()
		
		axs[0].legend(lines[0] + lines[3], labels[0] + labels[3], loc=0)
		axs[1].legend(lines[1] + lines[4], labels[1] + labels[4], loc=0)
		axs[2].legend(loc=0)

	def __getattr__(self, name):
		if name in ["Dec", "Horiz", "Vert", "Total", 'Incl']:
			return self.data[name]

	def __repr__(self):
		return "Magnetic data object."
