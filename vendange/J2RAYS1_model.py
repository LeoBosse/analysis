#!/usr/bin/python3
# -*-coding:utf-8 -*

import numpy as np
import pandas as pd
import matplotlib as mpl
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
import scipy.optimize as opt #import curve_fit


class Model:
	def __init__(self, data, obs_type = "fixed"):
		self.data = data

		self.obs_type = obs_type

		# print("DEBUG MODEL COLUMN NAME", self.data.columns)
		if self.data.shape[0] == 1:
			self.x_type = "index"
			self.x_axis = np.arange(1)

		elif "datetime" in self.data.columns and self.obs_type == "fixed":
			self.x_type = "datetime"
			try:
				self.data["datetime"] = [dt.datetime.strptime(t, "%Y%m%d-%H%M%S") for t in self.data["datetime"]]
				self.data["datetime"] = pd.to_datetime(self.data["datetime"])
			except:
				pass


		elif "azimuts" in self.data.columns:
			self.x_type = "azimuts"
			self.x_axis = self.data["azimuts"]

		# print("DEBUG INIT 1", self.data["RSV"], self.data["RSVcos"], self.data["RSDoLP"])
		self.SetMandatoryColumns()
		# print("DEBUG INIT 2", self.data["RSV"], self.data["RSVcos"], self.data["RSDoLP"], self.data["DVcos"], self.data["DDoLP"])
		self.UpdateFromFDA()
		# print("DEBUG INIT 3", self.data["RSV"], self.data["RSVcos"], self.data["RSDoLP"], self.data["DVcos"], self.data["DDoLP"])
		self.SetRSParam()
		# print("DEBUG INIT 4", self.data["RSV"], self.data["RSVcos"], self.data["RSDoLP"], self.data["DVcos"], self.data["DDoLP"])

		self.flux_average = np.average(self["I0"])

		self.bottle_function = np.average
		self.bottle_average = 1
		self.bottle_factor = 1

		# self.MakePlot()
		# plt.show()

	@staticmethod
	def InitFromBottle(bottle, RS_file, interp = None, time_divisor = 1, x_axis = None, shift = 0):
		filename = RS_file

		data = pd.read_csv(filename)

		if bottle.observation_type == "fixed_elevation_continue_rotation":
			old = data
			for i in range(1, bottle.nb_continue_rotation):
				new = old.copy()
				new["azimuts"] += i * 360
				data = pd.concat([data, new], ignore_index=True)

			data["azimuts"] += shift

		model_object = Model(data, obs_type = bottle.observation_type)

		# model_object.bottle_function = np.average
		model_object.bottle_factor = model_object.bottle_function(bottle.smooth_I0) / model_object.bottle_function(model_object["I0"])
		model_object.SetArbitraryUnits(model_object.bottle_factor)

		if model_object.x_type == "datetime":
			model_object.data = model_object.data.loc[(model_object.data["datetime"] > bottle.DateTime("start", format = "UT")) & (model_object.data["datetime"] < bottle.DateTime("end", format = "UT"))]
			model_object.x_axis = np.array([nt.total_seconds() / time_divisor for nt in model_object.data["datetime"] - bottle.DateTime("start", format = "UT")])
			# model_object.x_axis += shift


		elif model_object.x_type == "index" and x_axis is not None:
			model_object.single_data = model_object.data
			model_object.data = pd.concat([model_object.data] * int(len(x_axis)/100. + 1), ignore_index=True)
			model_object.x_axis = x_axis[::100]

		if model_object.obs_type == "fixed_elevation_discrete_rotation":
			# print("DEBUG MODEL fixed_elevation_discrete_rotation")
			model_object.x_type = "datetime"
			model_object.data["datetime"] = np.interp(model_object.data["azimuts"] * DtoR, bottle.discrete_rotation_azimuts, bottle.discrete_rotation_times)
			model_object.x_axis = model_object.data["datetime"]
			# print(model_object.data["azimuts"], bottle.discrete_rotation_azimuts, bottle.discrete_rotation_times, model_object.x_axis)

		return model_object


	def SetMandatoryColumns(self):
		length = len(self.data["DoLP"])

		for name in ["V", "Vcos", "Vsin", "RSI0", "RSDoLP", "RSAoRD", "RSV", "RSVcos", "RSVsin", "DI0", "DDoLP", "DAoLP", "DV", "DVcos", "DVsin"]:
			if name not in self.data.columns:
				self.data[name] = np.zeros(length)

	def GetVParam(self, F, D, A):
		V = F / 2.
		Vc = F * D/100. * np.cos(2 * A*DtoR) / 4.
		Vs = F * D/100. * np.sin(2 * A*DtoR) / 4.
		return V, Vc, Vs

	def GetFDAParam(self, V, Vc, Vs):
		F = V * 2.
		D = 2 * np.sqrt(Vc**2 + Vs**2) / V * 100
		A = np.arctan2(Vs, Vc) / 2 * RtoD
		return F, D, A

	def SetRSParam(self):
		self.data["RSV"] = self.data["V"] - self.data["DV"]
		self.data["RSVcos"] = self.data["Vcos"] - self.data["DVcos"]
		self.data["RSVsin"] = self.data["Vsin"] - self.data["DVsin"]

		self.data["RSI0"] = self.data["RSV"] * 2.
		self.data["RSDoLP"] = 2. * 100. * np.sqrt(self.data["RSVcos"]**2 + self.data["RSVsin"]**2) / self.data["RSV"]
		self.data["RSAoRD"] = np.arctan2(self.data["RSVsin"], self.data["RSVcos"]) / 2. * RtoD

	@staticmethod
	def GetVParamDF(data):
		data["V"] = data["I0"] / 2.
		data["Vcos"] = data["I0"] * (data["DoLP"]/100.) * np.cos(2 * data["AoRD"]*DtoR) / 4.
		data["Vsin"] = data["I0"] * (data["DoLP"]/100.) * np.sin(2 * data["AoRD"]*DtoR) / 4.

		data["DV"] = data["DI0"] / 2.
		data["DVcos"] = data["DI0"] * (data["DDoLP"]/100.) * np.cos(2 * data["DAoLP"]*DtoR) / 4.
		data["DVsin"] = data["DI0"] * (data["DDoLP"]/100.) * np.sin(2 * data["DAoLP"]*DtoR) / 4.

		data["RSV"] = data["RSI0"] / 2. #data["V"] - data["DV"]
		data["RSVcos"] = data["RSI0"] * (data["RSDoLP"]/100.) * np.cos(2 * data["RSAoRD"]*DtoR) / 4. #data["Vcos"] - data["DVcos"]
		data["RSVsin"] = data["RSI0"] * (data["RSDoLP"]/100.) * np.sin(2 * data["RSAoRD"]*DtoR) / 4. #data["Vsin"] - data["DVsin"]

		return data

	def SetVParam(self):
		self.data = Model.GetVParamDF(self.data)

	@staticmethod
	def GetFDAParamDF(data):
		data["I0"] = data["V"] * 2.
		data["DoLP"] = np.divide(2 * 100 * np.sqrt(data["Vcos"]**2 + data["Vsin"]**2), data["V"], where=data["V"] != 0, out = np.zeros(len(data["I0"])))
		data["AoRD"] = np.arctan2(data["Vsin"], data["Vcos"]) / 2. * RtoD

		data["DI0"] = data["DV"] * 2.

		data["DDoLP"] = np.divide(2. * 100. * np.sqrt(data["DVcos"]**2 + data["DVsin"]**2), data["DV"], where=data["DV"] != 0, out = np.zeros(len(data["I0"])))
		data["DAoLP"] = np.arctan2(data["DVsin"], data["DVcos"]) / 2. * RtoD

		data["RSI0"] = data["RSV"] * 2.
		data["RSDoLP"] = np.divide(2. * 100. * np.sqrt(data["RSVcos"]**2 + data["RSVsin"]**2), data["RSV"], where=data["RSV"] != 0, out = np.zeros(len(data["I0"])))
		data["RSAoRD"] = np.arctan2(data["RSVsin"], data["RSVcos"]) / 2. * RtoD

		return data

	def SetFDAParam(self):
		self.data = Model.GetFDAParamDF(self.data)

	def UpdateFromFDA(self):
		self.SetVParam()
		self.SetFDAParam()

	def UpdateFromVParam(self):
		self.SetFDAParam()
		self.SetVParam()


	def AddFDA(self, addF, addD, addA):
		avg_I0 = self.bottle_function(self["I0"])


		addV, addVcos, addVsin = self.GetVParam(addF, addD, addA)
		self["V"] 		+= addV
		self["Vcos"]	+= addVcos
		self["Vsin"]	+= addVsin

		self.UpdateFromVParam()

		self.bottle_factor = self.bottle_factor * avg_I0 / self.bottle_function(self["I0"])
		self.SetArbitraryUnits(avg_I0 / self.bottle_function(self["I0"]))

	def SetArbitraryUnits(self, flux_factor = 1, bottle = None):
		if bottle is not None:
			flux_factor = self.bottle_function(bottle.smooth_I0) / self.bottle_function(self["I0"])

		self["I0"] *= flux_factor
		self["DI0"] *= flux_factor
		self["RSI0"] *= flux_factor
		self["V"] *= flux_factor
		self["Vcos"] *= flux_factor
		self["Vsin"] *= flux_factor
		self["DV"] *= flux_factor
		self["DVcos"] *= flux_factor
		self["DVsin"] *= flux_factor
		self["RSV"] *= flux_factor
		self["RSVcos"] *= flux_factor
		self["RSVsin"] *= flux_factor


	def __add__(self, to_add):
		return Model(self.Add(to_add.data), obs_type = self.obs_type)

	def __radd__(self, to_add):
		return self + to_add

	def __mul__(self, f):
		prod = Model(self.data.copy(), obs_type = self.obs_type)
		# print("DEBUG PROD 1", prod["RSV"], prod["RSVcos"], prod["RSDoLP"])
		prod.data["I0"] 	*= f
		prod.data["DI0"] 	*= f
		prod.data["RSI0"] 	*= f
		prod.UpdateFromFDA()
		# print("DEBUG PROD 2", prod["RSV"], prod["RSVcos"], prod["RSDoLP"])
		return prod

	def __rmul__(self, f):
		return self * f

	def Add(self, add_dataframe, f1=1, f2=1):
		sum = self.data.copy()
		for n in ["V", "DV", "RSV"]:
			sum[n] = self.data[n] * f1 + add_dataframe[n] * f2

		for n in ["Vcos", "Vsin",  "DVcos", "DVsin", "RSVcos", "RSVsin"]:
			sum[n] += add_dataframe[n]

		sum = Model.GetFDAParamDF(sum)

		return sum

	def __getitem__(self, index):
		return self.data[index]

	def __setitem__(self, index, val):
		self.data[index] = val

	# @staticmethod
	# def InterpolateData(new_x_axis, x_axis, data):
	# 	new_data = pd.DataFrame()
	# 	print(len(new_x_axis), len(x_axis), len(data["V"]))
	# 	new_data["V"] = np.interp(new_x_axis, x_axis, data["V"])
	# 	new_data["Vcos"] = np.interp(new_x_axis, x_axis, data["Vcos"])
	# 	new_data["Vsin"] = np.interp(new_x_axis, x_axis, data["Vsin"])
	# 	new_data["DV"] = np.interp(new_x_axis, x_axis, data["DV"])
	# 	new_data["DVcos"] = np.interp(new_x_axis, x_axis, data["DVcos"])
	# 	new_data["DVsin"] = np.interp(new_x_axis, x_axis, data["DVsin"])
	#
	# 	return new_data
	#
	# def Interpolate(self, new_x_axis):
	# 	self.data = Model.InterpolateData(new_x_axis, self.x_axis, self.data)
	# 	self.x_axis = new_x_axis
	# 	self.UpdateFromVParam()

	def Copy(self):
		return Model(self.data.copy(), obs_type = self.obs_type)


	def FindDirectPola(self, bottle):
		x_axis = self.x_axis

		b_V = np.interp(x_axis, self.x_axis_list, bottle.smooth_V)
		b_Vcos = np.interp(x_axis, self.x_axis_list, bottle.smooth_Vcos)
		b_Vsin = np.interp(x_axis, self.x_axis_list, bottle.smooth_Vsin)

		diff_V = b_V - model["RSV"]
		diff_Vcos = b_Vcos - model["RSVcos"]
		diff_Vsin = b_Vsin - model["RSVsin"]

		diff_I0, diff_DoLP, diff_AoLP = Rotation.GetLightParameters(diff_V, diff_Vcos, diff_Vsin)



	def AddDirectPola(self, DoLP, AoLP):

		# print(self.data["RSI0"], self.data["RSDoLP"], self.data["RSAoRD"])

		new_model = self.Copy()

		DI0 = new_model.data["DI0"]
		DDoLP = DoLP #in % = [0, 100]
		DAoLP = AoLP #in degrees
		# print("DI0, DDoLP, DVsin", DI0, DDoLP, DAoLP)

		DV, DVcos, DVsin = FDAtoV(DI0, DDoLP / 100., DAoLP * DtoR)

		# print("DV, DVcos, DVsin", DV, DVcos, DVsin)

		new_model.data["DV"] = DV
		new_model.data["DVcos"] = DVcos
		new_model.data["DVsin"] = DVsin

		new_model.data["V"] = DV + new_model.data["RSV"]
		new_model.data["Vcos"] = DVcos + new_model.data["RSVcos"]
		new_model.data["Vsin"] = DVsin + new_model.data["RSVsin"]

		# print(new_model["RSV"], new_model["RSVcos"], new_model["RSVsin"])
		# print(new_model["V"], new_model["Vcos"], new_model["Vsin"])

		new_model.SetFDAParam()

		# print(new_model["I0"], new_model["DoLP"], new_model["AoRD"])

		return new_model



	@staticmethod
	def FitModelToFlux(x_axis, ptcu_flux, m1, m2, f1=1, f2=1):

		# print("DEBUG", m2.data.to_numpy(), m2.data.to_numpy()[0])
		# print("DEBUG", m1.data.to_numpy(), m1.data.to_numpy()[0])
		# print("DEBUG SHAPE", m1.data.shape, m2.data.shape)

		if m2.x_type == "index":
			# print("DEBUG", m2.data.values, m2.data.values[0])
			m1.data.reset_index(inplace=True, drop=True)
			m2.data = pd.concat([m2.single_data] * (len(m1.data)), ignore_index=True)
			m2.data.reset_index(inplace=True, drop=True)
			m2.x_axis = m1.x_axis
		elif m1.x_type == "index" and  m2.x_type != "index":
			m2.data.reset_index(inplace=True, drop=True)
			m1.data = pd.concat([m1.single_data] * (len(m2.data)), ignore_index=True)
			m1.data.reset_index(inplace=True, drop=True)
			m1.x_axis = m2.x_axis


		# print("DEBUG SHAPE", m1.data.shape, m2.data.shape)


		# print(m1.data, m2.data)
		ptcu_flux = np.interp(m1.x_axis, x_axis, ptcu_flux)

		avg_data = np.average(ptcu_flux)
		avg_m1 = np.average(m1.data["I0"])
		avg_m2 = np.average(m2.data["I0"])

		print(m1.data)
		print(m2.data)

		def least_square_flux(p, fact_1, fact_2):
			comblin = p[0] * fact_1 * m1.data["I0"] + p[1] * fact_2 * m2.data["I0"]
			# print(m1.data["I0"], m2.data["I0"])
			# print(comblin)
			least_square = sum((comblin - ptcu_flux) ** 2)
			return least_square

		p0 = [avg_data / avg_m1, avg_data / avg_m2]
		optimize = opt.minimize(least_square_flux, p0, args=(f1, f2))

		print("DEBUG FIT:", optimize, optimize.x)
		print("Best fit contributions:", np.average(m1["I0"] * optimize.x[0]), np.average(m2["I0"] * optimize.x[1]))
		print("Models bottle factor", m1.bottle_factor, m2.bottle_factor)


		# print("DEBUG bestfit 0", m1.data["RSV"], m1.data["RSVcos"], m1.data["RSDoLP"], m1.data["DVcos"])
		# print("DEBUG bestfit 0", m2.data["RSV"], m2.data["RSVcos"], m2.data["RSDoLP"], m2.data["DVcos"], m2.data["DV"])
		best = m1 * optimize.x[0]
		# print("DEBUG bestfit 1", best.data["RSV"], best.data["RSVcos"], best.data["RSDoLP"], best.data["DVcos"])
		best = best + m2 * optimize.x[1]
		# print("DEBUG bestfit 2", best.data["RSV"], best.data["RSVcos"], best.data["RSDoLP"], best.data["DVcos"])
		best.x_axis = m1.x_axis
		return best


	@staticmethod
	def FitModelToFluxPlusIso(x_axis, ptcu_flux, m1, m2, f1=1, f2=1):


		if m2.x_type == "index":
			# print("DEBUG", m2.data.values, m2.data.values[0])
			m1.data.reset_index(inplace=True, drop=True)
			m2.data = pd.concat([m2.single_data] * (len(m1.data)), ignore_index=True)
			m2.data.reset_index(inplace=True, drop=True)
			m2.x_axis = m1.x_axis
		elif m1.x_type == "index" and  m2.x_type != "index":
			m2.data.reset_index(inplace=True, drop=True)
			m1.data = pd.concat([m1.single_data] * (len(m2.data)), ignore_index=True)
			m1.data.reset_index(inplace=True, drop=True)
			m1.x_axis = m2.x_axis
		# print(m1.data, m2.data)
		ptcu_flux = np.interp(m1.x_axis, x_axis, ptcu_flux)


		avg_data = np.average(ptcu_flux)
		avg_m1 = np.average(m1.data["I0"])
		avg_m2 = np.average(m2.data["I0"])

		print(m1.data)
		print(m2.data)

		def least_square_flux(p, fact_1, fact_2):
			comblin = p[0] * fact_1 * m1.data["I0"] + p[1] * fact_2 * m2.data["I0"] + p[2]
			# print(m1.data["I0"], m2.data["I0"])
			# print(comblin)
			least_square = sum((comblin - ptcu_flux) ** 2)
			return least_square

		p0 = [avg_data / avg_m1, avg_data / avg_m2, 0]
		optimize = opt.minimize(least_square_flux, p0, args=(f1, f2))

		print("DEBUG FIT:", optimize, optimize.x)
		print("Best fit contributions:", np.average(m1["I0"] * optimize.x[0]), np.average(m2["I0"] * optimize.x[1]), optimize.x[2])

		best = m1 * optimize.x[0] + m2 * optimize.x[1]
		best.AddFDA(optimize.x[2], 0, 0)
		best.x_axis = m1.x_axis
		return best


	@staticmethod
	def FitModelToFluxAndDoLP(x_axis, ptcu_flux, ptcu_DoLP, m1, m2, f1=1, f2=1):

		ptcu_flux = np.interp(m1.x_axis, x_axis, ptcu_flux)
		ptcu_DoLP = np.interp(m1.x_axis, x_axis, ptcu_DoLP)

		avg_data = np.average(ptcu_flux)
		avg_m1 = np.average(m1.data["I0"])
		avg_m2 = np.average(m2.data["I0"])

		# add_model_flux = lambda i, p1, p2: p1 * f1 * m1.data["I0"][i] + p2 * f2 * m2.data["I0"][i]
		# param, param_cov = opt.curve_fit(add_model_flux, m1.x_axis, ptcu_flux)

		def least_square_flux(p, fact_1, fact_2):
			least_square = 0
			# comblin = (p[0] * fact_1 * m1.data["I0"] + p[1] * fact_2 * m2.data["I0"])
			comblin = (p[0] * fact_1 * m1 + fact_2 * m2)

			# for i in range(0, len(m1.x_axis)):
			flux_weight, DoLP_weight = 10, 1
			least_square = sum((comblin["I0"] * avg_data / np.average(comblin["I0"]) - ptcu_flux) ** 2) / np.average(ptcu_flux) * flux_weight + sum((comblin["DoLP"] - ptcu_DoLP) ** 2) / np.average(ptcu_DoLP) * DoLP_weight
			return least_square

		p0 = [avg_m2 / avg_m1]
		optimize = opt.minimize(least_square_flux, p0, args=(f1, f2))

		print("DEBUG FIT:", optimize, optimize.x)

		return m1 * optimize.x[0] + m2

	@staticmethod
	def FitModelToFluxRatio(x_axis, ptcu_flux, m1, m2, f1=1, f2=1):

		ptcu_flux = np.interp(m1.x_axis, x_axis, ptcu_flux)

		ratio = lambda d: (max(d) - min(d)) / min(d)
		avg_data = np.average(ptcu_flux)
		data_ratio = ratio(ptcu_flux)
		avg_m1 = np.average(m1.data["I0"])
		avg_m2 = np.average(m2.data["I0"])

		# add_model_flux = lambda i, p1, p2: p1 * f1 * m1.data["I0"][i] + p2 * f2 * m2.data["I0"][i]
		# param, param_cov = opt.curve_fit(add_model_flux, m1.x_axis, ptcu_flux)

		def least_square_flux_ratio(p, fact_1, fact_2):
			comblin = p[0] * fact_1 * m1.data["I0"] + (1 - p[0]) * fact_2 * m2.data["I0"]
			return abs(ratio(comblin) - data_ratio)

		p0 = [avg_m2 / avg_m1, 0]
		optimize = opt.minimize(least_square_flux_ratio, p0, args=(f1, f2), bounds = [(0, 1), (None, None)])

		param = optimize.x

		print("DEBUG FIT:", optimize, optimize.x)

		best = param[0] * f1 * m1 + (1 - param[0]) * f2 * m2
		# best.data["I0"] += param[1]
		best.UpdateFromFDA()

		return best

	@staticmethod
	def FitModelToDoLP(x_axis, ptcu_DoLP, m1, m2, f1=1, f2=1):

		ptcu_DoLP = np.interp(m1.x_axis, x_axis, ptcu_DoLP)

		avg_data = np.average(ptcu_DoLP)
		avg_m1 = np.average(m1.data["I0"])
		avg_m2 = np.average(m2.data["I0"])


		print("AVG", avg_data, avg_m1, avg_m2)
		m1.MakePlot()
		m2.MakePlot()

		# add_model_flux = lambda i, p1, p2: p1 * f1 * m1.data["I0"][i] + p2 * f2 * m2.data["I0"][i]
		# param, param_cov = opt.curve_fit(add_model_flux, m1.x_axis, ptcu_flux)

		def least_square_dolp(p, fact_1, fact_2):
			# least_square = 0
			comblin = m1 * p[0] * fact_1 + m2 * p[1] * fact_2
			comblin.UpdateFromFDA()
			least_square = sum((comblin.data["DoLP"] - ptcu_DoLP) ** 2)
			return least_square

		def match_variations_dolp(p, fact_1, fact_2):
			comblin = m1 * p[0] * fact_1 + m2 * fact_2
			return sum((np.gradient(ptcu_DoLP) - np.gradient(comblin.data["DoLP"])) ** 2)

		p0 = [avg_m2 / avg_m1]
		optimize = opt.minimize(match_variations_dolp, p0, args=(f1, f2))
		# optimize = opt.minimize(match_variations_dolp, p0, args=(f1, f2), bounds=[(avg_m2 / avg_m1 * 1e-5, avg_m2 / avg_m1 * 1e5), (-5, 5), (-5, 5)])

		print("DEBUG FIT:", optimize)
		p = optimize.x
		best = m1 * p[0] * f1 + m2 * f2
		# comblin.UpdateFromFDA()

		# best = m1 * 1 + m2 * 8e10 2.11825726e-10
		# best = m1 * 1 + m2 * 1e8

		best.MakePlot()
		plt.show()

		return best

	@staticmethod
	def FitModelToAoLP(x_axis, ptcu_AoLP, m1, m2, f1=1, f2=1):

		ptcu_AoLP = np.interp(m1.x_axis, x_axis, ptcu_AoLP)

		avg_data = np.average(ptcu_AoLP)
		avg_m1 = np.average(m1.data["I0"])
		avg_m2 = np.average(m2.data["I0"])

		# add_model_flux = lambda i, p1, p2: p1 * f1 * m1.data["I0"][i] + p2 * f2 * m2.data["I0"][i]
		# param, param_cov = opt.curve_fit(add_model_flux, m1.x_axis, ptcu_flux)

		def least_square_AoRD(p, fact_1, fact_2):
			model_sum = p[0] * fact_1 * m1 + fact_2 * m2
			return sum((model_sum.data["AoRD"] - ptcu_AoLP)**2)

		p0 = [avg_m2 / avg_m1]
		optimize = opt.minimize(least_square_AoRD, p0, args=(f1, f2), bounds = [(p0[0] * 1e-1, p0[0] * 1e1)])

		print("DEBUG FIT:", optimize)
		p = optimize.x
		best = m1 * p[0] * f1 + m2 * f2
		# comblin.UpdateFromFDA()

		# best = m1 * 1 + m2 * 8e10 2.11825726e-10
		# best = m1 * 1 + m2 * 1e8

		best.MakePlot()
		plt.show()

		return best






	# def FitDDoLPToDoLP(self, bottle, sky, grd, f1=1, f2=1):
	#
	# 	def add_model_DDoLP(i, p1, p2):
	# 		sky_DDoLP = p1 * sky["DDoLP"][i] + p2
	# 		sky_DV, sky_DVc, sky_DVs = GetVParam(sky["DI0"][i], p1 * sky["DDoLP"][i] + p2, sky["DAoLP"][i])
	#
	# 		sky_V, sky_Vc, sky_Vs = sky_DV + sky["RSV"], sky_DVc + sky["RSVc"], sky_DVs + sky["RSVs"]
	# 		sum_V, sum_Vc, sumVs = f1 * sky_V + f2 * grd["V"], sky_Vc + grd["Vc"], sky_Vs + grd["Vs"]
	#
	# 		sum_F, sum_D, sum_A = GetFDAParam(sum_V, sum_Vc, sumVs)
	# 	return sum_D
	#
	# 	param, param_cov = opt.curve_fit(add_model_DDoLP, np.arange(0, len(sky["DoLP"]), 1), data["SDoLP"])
	#
	# 	return param, param_cov


	def MakePlot(self):

		f, axs = plt.subplots(3, sharex=True)

		axs[0].plot(self.x_axis, self.data["I0"])
		axs[1].plot(self.x_axis, self.data["DoLP"])
		axs[2].plot(self.x_axis, self.data["AoRD"] * RtoD)
