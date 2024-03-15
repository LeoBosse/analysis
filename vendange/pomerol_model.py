#!/usr/bin/python3
# -*-coding:utf-8 -*


#######################################################################
# Contains the POMEROL model object.
# Methods are available to fit the model to the data, add and subtract moels from each other.
# This object reads the polarisation computed by a run of the POMEROL model (light pollutio calculator). 

# Author: Léo Bosse
# License: the freeest one, do whatever with my bad code if it can be helpfull, nobody really cares!
#######################################################################





import numpy as np
import pandas as pd
# import pymc3 as pm
import arviz as az
import matplotlib as mpl
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
from subprocess import call
import datetime as time
import scipy.optimize as opt #import curve_fit

from utils import *
from rotation import *
from bottle import *
from observation import *
from MagData import *

from vendange_configuration import *


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
			model_object.data.reset_index(drop=True, inplace=True)
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

	def SetXAxisAndLength(self, new_x_axis):
		new_len = len(new_x_axis)
		old_len = len(self.x_axis)
		if new_len <= old_len:
			self.data = self.data.iloc[:new_len]
		else:
			self.data = pd.DataFrame([list(self.data.iloc[0])] * new_len, columns = self.data.columns)
		self.x_axis = new_x_axis



	def SetMandatoryColumns(self):
		"""Call this after initialasing a new model. Set the missing columns as numpy zeros array."""
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


	def AddFDA(self, addF, addD, addA, direct=False):
		avg_I0 = self.bottle_function(self["I0"])

		addV, addVcos, addVsin = self.GetVParam(addF, addD, addA)

		self["V"] 		+= addV
		self["Vcos"]	+= addVcos
		self["Vsin"]	+= addVsin
		if direct:
			self["DV"] 		+= addV
			self["DVcos"]	+= addVcos
			self["DVsin"]	+= addVsin

		self.UpdateFromVParam()

		self.bottle_factor = self.bottle_factor * avg_I0 / self.bottle_function(self["I0"])
		self.SetArbitraryUnits(avg_I0 / self.bottle_function(self["I0"]))

	def SetArbitraryUnits(self, flux_factor = 1, bottle = None):
		"""Given a bottle, multiplies the fluxes of the model so that the flux units match the bottle units.
		Uses the function defined in __init__ to match the average/min/max flux of the bottle data with the average/min/max flux of the model."""

		if bottle is not None:
			flux_factor = self.bottle_function(bottle.smooth_I0) / self.bottle_function(self["I0"])

		for name in ["I0", "DI0", "RSI0", "V", "Vcos", "Vsin", "RSV", "RSVcos", "RSVsin", "DV", "DVcos", "DVsin"]:
			self[name] *= flux_factor


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

	# def FactorDoLP(self, DoLP_factor):
	#
	# 	new = Model(self.data.copy(), obs_type = self.obs_type)
	#
	# 	new.data["DoLP"] 	*= DoLP_factor
	# 	new.data["RSDoLP"] 	*= DoLP_factor
	# 	new.data["DDoLP"] 	*= DoLP_factor
	# 	new.UpdateFromFDA()
	#
	# 	return new



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
		new =  Model(self.data.copy(), obs_type = self.obs_type)
		new.flux_average = self.flux_average
		new.bottle_function = self.bottle_function
		new.bottle_average = self.bottle_average
		new.bottle_factor = self.bottle_factor
		# new.x_axis = self.x_axis
		return new


	# def FindDirectPola(self, bottle):
	# 	x_axis = self.x_axis
	#
	# 	b_V = np.interp(x_axis, self.x_axis_list, bottle.smooth_V)
	# 	b_Vcos = np.interp(x_axis, self.x_axis_list, bottle.smooth_Vcos)
	# 	b_Vsin = np.interp(x_axis, self.x_axis_list, bottle.smooth_Vsin)
	#
	# 	diff_V = b_V - model["RSV"]
	# 	diff_Vcos = b_Vcos - model["RSVcos"]
	# 	diff_Vsin = b_Vsin - model["RSVsin"]
	#
	# 	diff_I0, diff_DoLP, diff_AoLP = Rotation.GetLightParameters(diff_V, diff_Vcos, diff_Vsin)

	@staticmethod
	def BuildCompleteModel(p, bottle, grd_model, sky_model, DoLP_mode, AoLP_mode, eq_current = None):
		"""Give the parameters p: list of [current's azimut, current's elevation, current's max DoLP, ground proportion] in radians and % (0-100). Give a bottle, a ground and a sky model. The DoLP and AoLP mode for polarisation (see EqCurrent.GetPolarisation() for details).
		Returns a new model as a combination of the ground and sky model given in input. And a direct polarisation as defined by the current and modes."""
		if eq_current is None:
			if len(p) == 4:
				c_az, c_el, c_d, grd_c = p
				iso, iso_dir = 0, 0
			elif len(p) == 6:
				c_az, c_el, c_d, grd_c, iso, iso_dir = p

			current = np.array([np.sin(c_el),
								np.cos(c_el) * np.sin(c_az),
								np.cos(c_el) * np.cos(c_az)])

			current_x_axis = None
			# print("DEBUG simple current")

			# new_grd_model = grd_model.FactorDoLP(grd_D)
		else:
			if len(p) == 2:
				c_d, grd_c = p
				iso, iso_dir = 0, 0
			elif len(p) == 4:
				c_d, grd_c, iso, iso_dir = p

			current = eq_current
			current_x_axis = current["seconds"]
			# print("DEBUG complexe current")

		if iso != 0:
			sky_model = sky_model.Copy()
			sky_model.AddFDA(iso * iso_dir, 0, 0, direct=True)
			sky_model.AddFDA(iso * (1. - iso_dir), 0, 0, direct=False)
			sky_model.SetArbitraryUnits(bottle = bottle)

		model = grd_c * grd_model + (100 - grd_c) * sky_model
		# print(model["I0"], grd_model["I0"], sky_model["I0"])
		model.x_axis = grd_model.x_axis
		model.SetMandatoryColumns()
		model.SetArbitraryUnits(bottle = bottle)

		# print(model["I0"], model["DoLP"], model["AoRD"])

		DoLPs, AoLPs = model.GetDirectPola(current, bottle, DoLP_max = c_d, DoLP_mode=DoLP_mode, AoLP_mode=AoLP_mode, current_x_axis = current_x_axis)
		model = model.AddDirectPola(DoLPs * 100, AoLPs * RtoD)

		model.x_axis = grd_model.x_axis

		# if iso != 0:
		# 	model.AddFDA(iso, 0, 0)
		# 	model.SetArbitraryUnits(bottle = bottle)

		# print(model["I0"], model["DoLP"], model["AoRD"])

		return model

	def GetBottleLeastSquare(self, bottle, data_x, flux_factor=1, dolp_factor=1, aolp_factor=1):
		diff_flux, diff_DoLP, diff_AoLP = 0, 0, 0

		if flux_factor != 0:
			diff_flux = flux_factor * self.LeastSquare(data_x, bottle.smooth_I0, self.data["I0"], x2=self.x_axis, sigma1=bottle.std_smooth_I0, sigma2=0.1 * self.data["I0"])

		if dolp_factor != 0:
			diff_DoLP = dolp_factor * self.LeastSquare(data_x, bottle.smooth_DoLP, self.data["DoLP"], x2=self.x_axis, sigma1=bottle.std_smooth_DoLP)

		if aolp_factor != 0:
			diff_AoLP = aolp_factor * self.LeastSquare(data_x, bottle.smooth_AoLP, self.data["AoRD"]*DtoR, x2=self.x_axis, sigma1=bottle.std_smooth_AoLP, angle = True)

		return [diff_flux, diff_DoLP, diff_AoLP]


	@staticmethod
	def MCMC(bottle, data_x, grd_model, sky_model, DoLP_mode="rayleigh", AoLP_mode="para", ISO=False, least_square_factors = [1, 1, 1]):

		with pm.Model() as mcmc_model:
			# Priors for unknown model parameters
			current_az = pm.Uniform("current_az", lower = 0, upper = 2*np.pi)
			current_el = pm.Uniform("current_el", lower = 0, upper = np.pi/2)
			current_Dmax = pm.Uniform("current_Dmax", lower = 0, upper = 100)
			grd_prop = pm.Uniform("grd_prop", lower = 0, upper = 100)

			p = [current_az, current_el, current_Dmax, grd_prop]
			# Expected values of outcome
			POMEROL_model = Model.BuildCompleteModel(p, bottle, grd_model, sky_model, DoLP_mode, AoLP_mode)

			pred_flux = POMEROL_model["I0"]
			# pred_DoLP = POMEROL_model["DoLP"]
			# pred_AoLP = POMEROL_model["AoRD"]


			sigma = pm.HalfNormal("sigma", sigma=1)
			# Likelihood (sampling distribution) of observations
			flux_obs = pm.Normal("flux_obs", mu=pred_flux, sigma=sigma, observed=bottle.smooth_I0)

			trace = pm.sample(500, return_inferencedata=False)
			az.plot_trace(trace);

	@staticmethod
	def ManuallyFindBestDirectPola(bottle, data_x, grd_model, sky_model, DoLP_mode="rayleigh", AoLP_mode="para", ISO=False, least_square_factors = [1, 1, 1]):
		print("DEBUG MANUAL BEST MODEL:")
		Naz = 10
		Nel = 10
		Nd  = 10
		Ngrd = 10

		current_az 		= np.linspace(90*DtoR, 260*DtoR, Naz)
		current_el 		= np.linspace(0*DtoR, 30*DtoR, Nel)
		current_Dmax 	= np.linspace(0, 20, Nd)
		grd_prop		= np.linspace(0, 50, Ngrd)
		# current_az 		= np.linspace(0, 2 * np.pi, Naz)
		# current_el 		= np.linspace(0, np.pi/2., Nel)
		# current_Dmax 	= np.linspace(0, 30, Nd)
		# grd_prop		= np.linspace(0, 50, Ngrd)

		total_diff 	= np.zeros((Naz, Nel, Nd, Ngrd))
		min_total_diff 	= 1e60
		best_total_param 	= None
		best_total_index 	= None
		best_total_model 	= None

		fF, fD, fA = least_square_factors

		for ia, c_az in enumerate(current_az):
			for ie, c_el in enumerate(current_el):
				for id, c_d in enumerate(current_Dmax):
					for ig, grd_c in enumerate(grd_prop):
						p = [c_az, c_el, c_d, grd_c]
						model = Model.BuildCompleteModel(p, bottle, grd_model, sky_model, DoLP_mode, AoLP_mode)

						total_diff[ia, ie, id, ig] = sum(model.GetBottleLeastSquare(bottle, data_x, flux_factor = fF, dolp_factor = fD, aolp_factor = fA))

						if total_diff[ia, ie, id, ig] < min_total_diff:
							min_total_diff = total_diff[ia, ie, id, ig]
							best_total_param = p
							best_total_index = [ia, ie, id, ig]
							best_total_model = model
							print("New Best param az, el, dmax, grd_prop:", c_az*RtoD, c_el*RtoD, c_d, grd_c)
							print("new min_total_diff:", min_total_diff)

		print("FINAL DEBUG MANUAL BEST MODEL:")
		c_az, c_el, c_d, grd_c = best_total_param
		print("Best param az, el, dmax, grd_prop:", c_az*RtoD, c_el*RtoD, c_d, grd_c)
		print("min_total_diff:", min_total_diff)
		return best_total_model, total_diff, min_total_diff, best_total_param, best_total_index


	@staticmethod
	def FindBestDirectPola(bottle, data_x, grd_model, sky_model, DoLP_mode="rayleigh", AoLP_mode="para", ISO=False, least_square_factors = [1, 1, 1], eq_current = None):

		print("Fitting POMEROL model to data...")

		fF, fD, fA = least_square_factors

		def GetDiff(p, bottle, data_x, grd_model, sky_model, DoLP_mode, AoLP_mode):
			model = Model.BuildCompleteModel(p, bottle, grd_model, sky_model, DoLP_mode, AoLP_mode, eq_current = eq_current)
			least_square = model.GetBottleLeastSquare(bottle, data_x, flux_factor = fF, dolp_factor = fD, aolp_factor = fA)
			return np.sum(least_square)


		az_0, 	az_bnds 	= 240 * DtoR, 	(0*DtoR, 360*DtoR) #2 * np.pi)
		el_0, 	el_bnds 	= 20 * DtoR, 	(0, 90*DtoR)#(0, np.pi / 2.)
		d_0,  	d_bnds  	= 8, 	(0, 50)
		grd_0, 	grd_bnds 	= 20,  	(0, 50)
		iso_0, iso_bnds		= 0, 	(0, None)
		iso_dir_0, iso_dir_bnds		= 1., 	(0., 1.)
		# grd_D_0,grd_D_bnds	= 1, 	(0, 1)

		if eq_current is None:
			p0 = [az_0, el_0, d_0, grd_0]#, grd_D_0]
			bounds = [az_bnds, el_bnds, d_bnds, grd_bnds]#, grd_D_bnds]
		else:
			p0 = [d_0, grd_0]#, grd_D_0]
			bounds = [d_bnds, grd_bnds]#, grd_D_bnds]

		if ISO:
			p0.append(iso_0)
			p0.append(iso_dir_0)
			bounds.append(iso_bnds)
			bounds.append(iso_dir_bnds)

		methods = ["TNC"]#, "L-BFGS-B"] ### See https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html

		i = 0
		optimize = opt.minimize(GetDiff, p0, method = methods[i], args=(bottle, data_x, grd_model, sky_model, DoLP_mode, AoLP_mode), bounds=bounds, options={"maxiter": 500})
		print(f"DEBUG FIT: ({methods[i]})\n", optimize)
		while (not optimize.success) and (i+1 < len(methods)):
			i += 1
			optimize = opt.minimize(GetDiff, p0, method = methods[i], args=(bottle, data_x, grd_model, sky_model, DoLP_mode, AoLP_mode), bounds=bounds)
			print(f"DEBUG FIT: ({methods[i]})\n", optimize)

		if eq_current is None:
			if not ISO:
				az_best, el_best, d_best, grd_best = optimize.x
			else:
				az_best, el_best, d_best, grd_best, iso_best, iso_dir_best = optimize.x
		else:
			if not ISO:
				d_best, grd_best = optimize.x
			else:
				d_best, grd_best, iso_best, iso_dir_best = optimize.x

		best = Model.BuildCompleteModel(optimize.x, bottle, grd_model, sky_model, DoLP_mode, AoLP_mode, eq_current = eq_current)


		diff_flux, diff_DoLP, diff_AoLP = best.GetBottleLeastSquare(bottle, data_x, flux_factor = fF, dolp_factor = fD, aolp_factor = fA)

		if eq_current is None:
			current = np.array([m.sin(el_best),
			m.cos(el_best) * m.sin(az_best),
			m.cos(el_best) * m.cos(az_best)])
			print(f"BEST PARAM: Current: (az, el) = ({az_best * RtoD}, {el_best * RtoD})")
			print(f"BEST PARAM: Current: (UEN) = ({current[0]}, {current[1]}, {current[2]})")
		print(f"BEST PARAM: Max direct pola (%) = {d_best}")
		print(f"BEST PARAM: Ground proportion (%) = {grd_best}")
		if ISO:
			print(f"BEST PARAM: Total bkgrd (AU) = {iso_best}")
			print(f"BEST PARAM: Polar bkgrd (% of total bkgrd) = {iso_dir_best*100}")

		print(f"BEST PARAM: least square F, D, A = {diff_flux}, {diff_DoLP}, {diff_AoLP}")

		return best, optimize.x

	@staticmethod
	def MakeParamSpacePlot(best_param, free_param_index, bottle, data_x, grd_model, sky_model, DoLP_mode="rayleigh", AoLP_mode="para", ISO=False, least_square_factors = [1, 1, 1], eq_current = None):
		"""
		Compute, plot and show a parameter space 2D cut around a set of best fit parameters.
		Inputs:
		best_param: the list of best fit parameters around which computing the 2D map
		free_param_index: the indexes of the "best_param" list to use as free parameters. Should be listed in increasing order.
							0: azimuth, 1: elevation, 2: Max DolP, 3: ground proportion
		"""
		print("Plotting parameter map...")

		if eq_current is None:
			if len(best_param) == 6:
				b_az, b_el, b_d, b_g, b_iso, b_iso_dir = best_param
			elif len(best_param) == 4:
				b_az, b_el, b_d, b_g = best_param
				b_iso, b_iso_dir = 0, 0
		else:
			if len(best_param) == 4:
				b_d, b_g, b_iso, b_iso_dir = best_param
				b_az, b_el = 0, 0
			elif len(best_param) == 2:
				b_d, b_g = best_param
				b_az, b_el, b_iso, b_iso_dir = 0, 0, 0, 0


		current_az 		= np.array([b_az])
		current_el 		= np.array([b_el])
		current_Dmax 	= np.array([b_d])
		grd_prop		= np.array([b_g])
		iso_list		= np.array([b_iso])
		iso_dir_list	= np.array([b_iso_dir])

		N = [100, 10, 20, 20, 10, 10]
		def GetAxis(index):
			if 0 == index:
				axis_list	= np.linspace(0*DtoR, 360*DtoR, N[0])
				# axis_title 	= "Azimut (°)"
				axis_title 	= "Azimuth (°)"
			elif 1 == index:
				axis_list 	= np.linspace(0*DtoR, 90*DtoR, N[1])
				# axis_title 	= "Élévation (°)"
				axis_title 	= "Elevation (°)"
			elif 2 == index:
				axis_list 	= np.linspace(0, 40, N[2])
				# axis_title 	= "Maximum DoLP (%)"
				axis_title 	= "Max DoLP (%)"
			elif 3 == index:
				axis_list	= np.linspace(0, 50, N[3])
				# axis_title 	= "Proportion Sol (%)"
				axis_title 	= "Ground proportion (%)"
			elif 4 == index:
				axis_list	= np.linspace(0, 2*np.max(bottle.smooth_I0), N[4])
				# axis_title 	= "Background isotropique (AU)"
				axis_title 	= "Isotropic background (AU)"
			elif 5 == index:
				axis_list	= np.linspace(0, 1, N[5])
				# axis_title 	= "Background isotropique polarisé (%)"
				axis_title 	= "Pola Isotropic background (%)"
			return axis_list, axis_title

		if 0 in free_param_index: current_az = GetAxis(0)[0]
		if 1 in free_param_index: current_el = GetAxis(1)[0]
		if 2 in free_param_index: current_Dmax = GetAxis(2)[0]
		if 3 in free_param_index: grd_prop = GetAxis(3)[0]
		if 4 in free_param_index: iso_list = GetAxis(4)[0]
		if 5 in free_param_index: iso_dir_list = GetAxis(5)[0]

		y_axis, y_title = GetAxis(free_param_index[0])
		x_axis, x_title = GetAxis(free_param_index[1])

		Naz		= len(current_az)
		Nel		= len(current_el)
		Nd		= len(current_Dmax)
		Ngrd	= len(grd_prop)
		Niso	= len(iso_list)
		Niso_dir	= len(iso_dir_list)

		total_diff 	= np.zeros((Naz, Nel, Nd, Ngrd, Niso, Niso_dir))

		fF, fD, fA = least_square_factors

		# print("DEBUG", current_Dmax, Nd)

		for ia, c_az in enumerate(current_az):
			for ie, c_el in enumerate(current_el):
				for id, c_d in enumerate(current_Dmax):
					for ig, grd_c in enumerate(grd_prop):
						for iiso, iso in enumerate(iso_list):
							for iisodir, isodir in enumerate(iso_dir_list):
								if eq_current is None:
									p = [c_az, c_el, c_d, grd_c, iso, isodir]
								else:
									p = [c_d, grd_c, iso, isodir]

								# print(c_az, c_el, c_d, grd_c)
								model = Model.BuildCompleteModel(p, bottle, grd_model, sky_model, DoLP_mode = DoLP_mode, AoLP_mode = AoLP_mode, eq_current = eq_current)

								total_diff[ia, ie, id, ig, iiso, iisodir] = sum(model.GetBottleLeastSquare(bottle, data_x, flux_factor = fF, dolp_factor = fD, aolp_factor = fA))


		map = total_diff[:, :, :, :, :, :].reshape(N[free_param_index[0]], N[free_param_index[1]])
		if free_param_index == (0, 1):
			free_param_index = (1, 0)
			# print("map shape", map.shape)
			map = map.transpose()
			# print("map shape", map.shape)
			y_axis, y_title = GetAxis(free_param_index[0])
			x_axis, x_title = GetAxis(free_param_index[1])

		extent = [x_axis[0], x_axis[-1], y_axis[0], y_axis[-1]]
		if free_param_index[0] in [0, 1]:
			extent[2] = y_axis[0]*RtoD
			extent[3] = y_axis[-1]*RtoD
		if free_param_index[1] in [0, 1]:
			extent[0] = x_axis[0]*RtoD
			extent[1] = x_axis[-1]*RtoD

		fig, axs = plt.subplots()
		pos = axs.imshow(map, origin='lower', extent = extent)
		axs.contour(x_axis*RtoD, y_axis*RtoD, map, 1e5 * np.array([1, 2, 3]), cmap='Wistia') #iso is a list of the altitudes isobar to plot (in km)
		# axs.contour(x_axis, y_axis, map, 1e6 * np.array([0.05, 0.1, 0.5, 1, 3]), cmap='Wistia') #iso is a list of the altitudes isobar to plot (in km)
		plt.xlabel(x_title)
		plt.ylabel(y_title)
		fig.colorbar(pos, ax = axs)
		plt.show()

	@staticmethod
	def LeastSquare(x1, f1, f2, x2=None, sigma1=None, sigma2=None, angle=False):
		"""Given a function 1 (absciss=x1, ordinates=f1) and a function 2 (ordinates f2, absciss optional x2). And error bars on f1 and f2, compute the square difference between both functions.
		If the function unit is in radians, specify angle=True to correctly compute the angle difference in [-pi/2: pi/2].
		"""

		if x2 is None: #Set x2 if not specified
			x2 = x1

		### Linear interpolate the smallest function (least defined points) to match the longest function x axis.
		if len(x1) > len(x2):
			f2 = np.interp(x1, x2, f2)
			if sigma2 is not None:
				sigma2 = np.interp(x1, x2, sigma2)
			x2 = x1
		elif len(x1) < len(x2):
			f1 = np.interp(x2, x1, f1)
			if sigma1 is not None:
				sigma1 = np.interp(x2, x1, sigma1)
			x1 = x2

		### Simple differences
		if angle: # If functions units are angles, the difference is set between [-pi/2: pi/2].
			diff_square = GetAnglesDifference(f1, f2)
		else:
			diff_square = f1 - f2

		### Differences square
		diff_square *= diff_square

		### Divide by the error squared if given.

		error_divisor = np.zeros_like(diff_square)
		if sigma1 is not None:
			error_divisor += sigma1**2
		if sigma2 is not None:
			error_divisor += sigma2**2

		if sum(error_divisor) > 0:
			diff_square /= error_divisor

		return np.sum(diff_square)


	def GetDirectPola(self, current, bottle, DoLP_mode="rayleigh", AoLP_mode="para", DoLP_max = 100, current_x_axis = None):
		"""Given a current (in up-east-north coord), a bottle, and the polarisation mode, returns a list of DoLPs and AoLP for each observations (each row of the model). See EqCurrent.GetPolarisation() for details on polarisation modes.
		(5/04/2021): Might not work for discrete rotations...
		current_x_axis: None if given only one vector current. If current is a list of vectors, current_x_axis is the time/azimuth of each current vector.
		"""

		pola_list = []

		if self.x_type == "azimuts": # Continue rotation observation
			obs_list = [ObservationPoint(A_lon=0, A_lat=0, observed_altitude=110, azimuth=az*DtoR, elevation=bottle.elevation, init_full = False) for az in self.data["azimuts"]]

			for obs in obs_list:
				pola_list.append(EqCurrent.GetPolarisation(current, obs, DoLP_mode=DoLP_mode, AoLP_mode=AoLP_mode, DoLP_max = DoLP_max))

		elif self.x_type == "datetime": # Fixed observation
			bottle.observation
			if current_x_axis is None:
				pola_list.append(EqCurrent.GetPolarisation(current, bottle.observation, DoLP_mode=DoLP_mode, AoLP_mode=AoLP_mode, DoLP_max = DoLP_max))

				pola_list = pola_list * len(self.x_axis) # If fixed observation, do the computation once, and multiply to match the length of the model x_axis.

			else:
				# cu, ce, cn = zip(*current)
				# ncu = np.interp(self.x_axis, current_x_axis, cu)
				# nce = np.interp(self.x_axis, current_x_axis, ce)
				# ncn = np.interp(self.x_axis, current_x_axis, cn)
				#
				# current = np.array(zip(ncu, nce, ncn))

				for cur in zip(current["Ju"] / current["J_norm"], current["Je"] / current["J_norm"], current["Jn"] / current["J_norm"]):
					# print(cur)
					pola_list.append(EqCurrent.GetPolarisation(cur, bottle.observation, DoLP_mode=DoLP_mode, AoLP_mode=AoLP_mode, DoLP_max = DoLP_max))



		#Unzip the DoLP and AoLP results
		DoLP_list, AoLP_list = zip(*pola_list)
		DoLP_list, AoLP_list = np.array(DoLP_list), np.array(AoLP_list)

		# print(DoLP_list*100, AoLP_list*RtoD)

		return DoLP_list, AoLP_list #[0-1] and radians



	def AddDirectPola(self, DoLP, AoLP):
		"""Add direct polarisation using the flux of the self model. DoLP given in % (0-100). AoLP given in degrees.
		DoLP and AoLP are lists as long as the model x_axis
		"""

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

		# print(m1.data)
		# print(m2.data)

		def least_square_flux(p, fact_1, fact_2):
			comblin = p[0] * fact_1 * m1.data["I0"] + p[1] * fact_2 * m2.data["I0"]
			# print(m1.data["I0"], m2.data["I0"])
			# print(comblin)
			least_square = sum((comblin - ptcu_flux) ** 2)
			return least_square

		p0 = [avg_data / avg_m1, avg_data / avg_m2]
		optimize = opt.minimize(least_square_flux, p0, args=(f1, f2))

		# print("DEBUG bestfit 0", m1.data["RSV"], m1.data["RSVcos"], m1.data["RSDoLP"], m1.data["DVcos"])
		# print("DEBUG bestfit 0", m2.data["RSV"], m2.data["RSVcos"], m2.data["RSDoLP"], m2.data["DVcos"], m2.data["DV"])
		best = m1 * optimize.x[0] + m2 * optimize.x[1]
		# print("DEBUG bestfit 1", best.data["RSV"], best.data["RSVcos"], best.data["RSDoLP"], best.data["DVcos"])
		# print("DEBUG bestfit 2", best.data["RSV"], best.data["RSVcos"], best.data["RSDoLP"], best.data["DVcos"])
		best.x_axis = m1.x_axis

		print("DEBUG FIT:", optimize, optimize.x)
		c1, c2, s = np.average(m1["I0"] * optimize.x[0]), np.average(m2["I0"] * optimize.x[1]), np.average(best["I0"])
		print(f"Best fit contributions: , {c1} ({c1/s}), {c2} ({c2/s}), tot: {s}")
		# print("Models bottle factor", m1.bottle_factor, m2.bottle_factor)
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

		# print(m1.data)
		# print(m2.data)

		def least_square_flux(p, fact_1, fact_2):
			comblin = p[0] * fact_1 * m1.data["I0"] + p[1] * fact_2 * m2.data["I0"] + p[2]
			# print(m1.data["I0"], m2.data["I0"])
			# print(comblin)
			least_square = sum((comblin - ptcu_flux) ** 2)
			return least_square

		p0 = [avg_data / avg_m1, avg_data / avg_m2, 0]
		optimize = opt.minimize(least_square_flux, p0, args=(f1, f2))


		best = m1 * optimize.x[0] + m2 * optimize.x[1]
		tmp_best_avg = np.average(best["I0"]) + optimize.x[2]
		best.AddFDA(optimize.x[2], 0, 0)
		best.x_axis = m1.x_axis

		print("DEBUG FIT:", optimize, optimize.x)
		c1, c2, i, s = np.average(m1["I0"] * optimize.x[0]), np.average(m2["I0"] * optimize.x[1]), optimize.x[2], tmp_best_avg
		print(f"Best fit contributions: , {c1} ({c1/s}), {c2} ({c2/s}), {i} ({i/s}), tot: {s}")

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
