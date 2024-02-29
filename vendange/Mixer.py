#!/usr/bin/python3
# -*-coding:utf-8 -*

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
# plt.ion()
from matplotlib import rc
#matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage{lmodern}']
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=False)
mpl.rcParams['font.size'] = 16
mpl.rcParams['svg.fonttype'] = 'none'
import datetime as dt

import chaosmagpy as chaos
import scipy as sci
from scipy import signal
import sys
import os
from subprocess import call

import observation as obs
from geometry import GetLonLatFromName, UENToAzEl
from itertools import combinations

from utils import *
from rotation import *
from bottle import *
from AllSkyData import *
from eiscat_data import *
from MagData import *
from plip import *
from pomerol_model import *
from taster import *

from vendange_configuration import *

class Mixer:
	def __init__(self, bottles, comp_bottles=[], mode = "default"):

		self.bottles = bottles
		self.comp_bottles = comp_bottles
		self.mag_data = None
		self.pola_color = None
		
		self.mode = mode

		for ib, bottle in enumerate(self.bottles):

			self.SetGraphParameter(bottle)
			# self.MakeXAxis(bottle)
			
			
			self.LoadExternData(bottle)


			if self.comp_bottles:
				self.TriangulateAoLP(bottle, self.comp_bottles)
			# self.MakeFigure()
			# self.MakePlots(bottle)
			# self.MakeCleanCorrelationPlots(bottle, None, smooth=True)

			# self.MakeFigure()
			# self.MakePlots(bottle)
			# # self.MakeCleanCorrelationPlots(bottle, None, smooth=True)

			# # self.MakeFFT(bottle)

			# if self.show_correlations:
			# 	# self.MakeCleanCorrelationPlots(bottle, None, smooth=True, COMP = "DoLP")
			# 	# self.MakeCleanCorrelationPlots(bottle, None, smooth=True, COMP = "AoLP")
			# 	self.MakeCleanCorrelationPlots(bottle, None, smooth=True)
			# 	# self.MakeCleanCorrelationPlots(bottle, None, smooth=False)

			# # for m in self.pomerol_models:
			# # 	self.SubtractModel(bottle, m)
			# # if self.show_RS_model:
			# # 	self.CompareRayleigh(bottle, self.pomerol_model)
			# if self.show_SN:
			# 	self.MakeI0SNratio(bottle)
			# 	self.MakeSNratio(bottle)
			# if self.make_optional_plots:
			# 	if self.show_currents and self.eq_currents.valid and ib == 0:
			# 		self.eq_currents.MakePlot(coords = "uen")
			# 		self.eq_currents.MakePlot(coords = "azel")

			# if self.comp_bottles:
			# 	print('DEBUG COMP BOTTLES', self.comp_bottles)

			# 	# for b in self.comp_bottles:
			# 	# 	b.TestResample()
				
			# 	### If the same number of channels in bottles and comp_bottles, plot them one for one
			# 	### If not, then plot all comp_bottles on top of of every bottle.
			# 	comp_bottles_to_show = [ib]
			# 	if len(self.comp_bottles) != len(self.bottles):
			# 		comp_bottles_to_show = range(len(self.comp_bottles))
				
			# 	for ibc in comp_bottles_to_show:

			# 		self.SetGraphParameter(self.comp_bottles[ibc], comp=True)
			# 		self.MakeXAxis(self.comp_bottles[ibc])
			# 		self.MakePlots(self.comp_bottles[ibc], comp=True)
					
			# 		if self.make_optional_plots:
			# 			if self.show_SN:
			# 				self.MakeSNratio(self.comp_bottles[ibc])

			# 			self.SetGraphParameter(bottle)
			# 			self.MakeXAxis(bottle)

			# 			self.CompareBottles(self.comp_bottles[ibc], bottle)

			# # if self.make_optional_plots:
			# # 	self.MakeFFT(bottle)

			# if self.make_optional_plots and self.mag_data is not False:
			# 	if ib == 0:
			# 		self.MakeMagDataPlots(bottle)
			# 	self.MakeCorrelationPlots(bottle)

			# if self.show_eiscat and self.make_optional_plots and self.eiscat_data.valid:
			# 	self.MakeEiscatDataPlot(bottle)

			# if self.make_optional_plots and bottle.observation_type == "fixed_elevation_continue_rotation":
			# 	self.CompareAnglesPlots(bottle)

			# 	if len(bottle.continue_rotation_times) > 2:
			# 		print("AVERAGE ROTATIONS")
			# 		bottle = bottle.GetAverageContinueRotations()
			# 		self.SetGraphParameter(bottle)
			# 		self.MakeXAxis(bottle)
			# 		self.MakeFigure()
			# 		self.MakePlots(bottle)
			# 		self.MakeSNratio(bottle)
			# 		self.CompareAnglesPlots(bottle)



			# self.MakeCoherencePLots()
	def LoadExternData(self, bottle):

		if self.show_mag_data:
			# self.LoadMagData(bottles[0])
			self.LoadMagData(bottle)

		if self.show_eiscat:
			self.LoadEiscatData(bottle)
		
		if self.show_allsky:
			self.LoadASCData(bottle)

		if self.show_currents:
			self.LoadCurrentData(bottle)
		
		if self.show_RS_model and self.RS_model_files:
			self.LoadPOMEROLModel(bottle)

		self.LoadPLIPData(bottle)


	def LoadPLIPData(self, bottle):
		# self.plip_data = PLIP(bottle, "data_NOT_ROT.h5")
		self.plip_data = PLIP("data_vh-+.h5")
		

	def LoadMagData(self, bottle):
		self.mag_data = MagData.FromBottle(bottle)
		if not self.mag_data.exist:
			self.mag_data = None
		else:
			self.mag_data.StripTime()
			self.mag_data
		# if self.mag_data is not False:
		# 	self.mag_data.StripTime(bottle.DateTime("start"), bottle.DateTime("end"))

	def LoadASCData(self, bottle):
		self.allsky_data_available = False
		if not bottle.observation_type == "fixed":
			return 
		self.allsky_data = AllSkyData(bottle)
		self.allsky_data_available = bool(self.allsky_data.all_datetimes)

	def LoadEiscatData(self, bottle):
		self.eiscat_data = Eiscat(bottle)
		if not self.eiscat_data.valid:
			# print("EISCAT HDF5")
			self.eiscat_data = EiscatHDF5(bottle, antenna = self.eiscat_type)
			print("EISCAT data valid? ", self.eiscat_data.valid)


	def LoadCurrentData(self, bottle):
		self.eq_currents = EqCurrent(bottle, file_type=None)
		print("Equivalent current is valid? ", bool(self.eq_currents.valid))
		
		if not self.eq_currents.valid:
			return

		self.eq_currents.GetApparentAngle(bottle.observation, Jup_mode = self.compute_Jup)#, shift=bottle.graph_angle_shift)
		self.eq_currents.Polarisation(bottle.observation, DoLP_max = 1, DoLP_mode="rayleigh", AoLP_mode="perp")
		
		if not self.eq_currents.valid:
			return

		if self.compute_Jup and bottle.observation_type == "fixed":
			AoLPs = [bottle.GetInterpolation(t)[2] for t in self.eq_currents.GetNormTimes()]
			self.eq_currents.FindJup(bottle.observation, AoLPs, mode = self.compute_Jup)
	
			self.eq_currents.GetApparentAngle(bottle.observation, Jup_mode = self.compute_Jup)#, shift=bottle.graph_angle_shift)
	
			self.eq_currents.Polarisation(bottle.observation, DoLP_max = 1, DoLP_mode="rayleigh", AoLP_mode="perp")
	
			self.eq_currents.SaveTXT()


	def LoadPOMEROLModel(self, bottle):
		
		self.pomerol_models = [Model.InitFromBottle(bottle, f, time_divisor = self.divisor, x_axis = self.x_axis_list, shift = self.shift_model) for i, f in enumerate(self.RS_model_files[bottle.line - 1])]
		# print(self.pomerol_models[1].data.index)
		if bottle.observation_type == "fixed": #fixed obs -> set the xaxis of the ground (should be at index 0) model to the sky model (should be at index 1)
			self.pomerol_models[0].SetXAxisAndLength(self.pomerol_models[1].x_axis)


		# self.pomerol_models[0]["DoLP"] = self.pomerol_models[0]["DoLP"]  / 10.
		# self.pomerol_models[0].SetVParam()
		# self.pomerol_models[1]["DoLP"] = self.pomerol_models[1]["DoLP"]  * 2
		# self.pomerol_models[1].SetVParam()


		####Find the best model with equivalent current object
		# best_model, best_param = Model.FindBestDirectPola(bottle, self.x_axis_list, self.pomerol_models[0], self.pomerol_models[1], DoLP_mode = "rayleigh", AoLP_mode = "perp", ISO=True, least_square_factors = [1, 1, 1], eq_current = self.eq_currents.GetInterpolation(self.pomerol_models[0].x_axis, bottle.observation, divisor=self.divisor))
		# self.pomerol_models.append(best_model)
		####Find the best model with a constant and uniform current
		# best_model, best_param = Model.FindBestDirectPola(bottle, self.x_axis_list, self.pomerol_models[0], self.pomerol_models[1], DoLP_mode = "rayleigh", AoLP_mode = "perp", ISO=False, least_square_factors = [1, 1, 1], eq_current = None)
		# self.pomerol_models.append(best_model)

		### Buikd a model ans compute its square difference
		# model = Model.BuildCompleteModel((0, 0, 0, 30), bottle, self.pomerol_models[0], self.pomerol_models[1], DoLP_mode = "rayleigh", AoLP_mode = "perp")
		# least_square = model.GetBottleLeastSquare(bottle, self.x_axis_list)
		# self.pomerol_models.append(model)
		# print("least square 20%", least_square, np.sum(least_square))


		### For 20200224 ski sud green best model
		# best_param = (120*DtoR, 0*DtoR, 1, 0)

		### For 20200227 green best model
		# best_param = (230*DtoR, 25*DtoR, 5, 20)
		# best_param = (50*DtoR, -20*DtoR, 5, 20)

		### For 20200227 purple best model
		# best_param = (230*DtoR, 20*DtoR, 11.1, 10.7, 6.24, 0.35)

		### For 20190307 green best model
		# best_param = (244*DtoR, 5*DtoR, 1.3, 7.1)#, 6.24, 0.35) ### rotation e45
		# best_param = (4.49678563, 0.37226045, 0.        , 9.13981123) # a164e45 with free current
		# best_param = (0.3388598, 8.3487529) # with Equivalent current
		### For 20190307 purple best model
		# best_param = (249*DtoR, 2.3*DtoR, 3, 10, 0.00093955034896, 1.) ### rotation e45
		# best_param = (4.29988166e+00, 3.44626175e-01, 6.36127072e+00, 1.72867582e+01, 8.01634464e-04, 9.99956375e-01)  # with free current direction
		# best_param = (4.89073016, 21.40499469, 74.189954  ,  1.) # with Equivalent current
		### For 20190307 blue best model
		# best_param = (232*DtoR, 0*DtoR, 5.9, 26, 4.27494363, 0.98) ### rotation e45
		# best_param = (50*DtoR, 10*DtoR, 6, 25, 20, 0.98)
		# best_param = (3.89422949e+00, 3.81413283e-01, 4.09738726e+00, 1.29324498e+01, 8.13699802e-09, 9.90971225e-01)  # with free current
		# best_param = (2.90791549e+00, 1.97557372e+01, 9.26151992e-15, 1.00000000e+00)  # with Equivalent current


		# Model.MakeParamSpacePlot(best_param, (0,1), bottle, self.x_axis_list, self.pomerol_models[0], self.pomerol_models[1], DoLP_mode = "rayleigh", AoLP_mode = "perp", ISO=False, least_square_factors = [1, 1, 1])#, eq_current = self.eq_currents.GetInterpolation(self.pomerol_models[0].x_axis, bottle.observation, divisor=self.divisor))


		#Build a model with equivalent current object
		# self.pomerol_models.append(Model.BuildCompleteModel(best_param, bottle, self.pomerol_models[0], self.pomerol_models[1], DoLP_mode = "rayleigh", AoLP_mode = "perp", eq_current = self.eq_currents.GetInterpolation(self.pomerol_models[0].x_axis, bottle.observation, divisor=self.divisor)))
		# #Build a model with a constant and uniform current
		# best_param = (3.89422949e+00, 3.81413283e-01, 4.09738726e+00, 1.29324498e+01, 8.13699802e-09, 9.90971225e-01)  # with free current
		# self.pomerol_models.append(Model.BuildCompleteModel(best_param, bottle, self.pomerol_models[0], self.pomerol_models[1], DoLP_mode = "rayleigh", AoLP_mode = "perp", eq_current = None))

		# print(self.pomerol_models[-1].data)

		# for m in self.pomerol_models:
		# 	print(m.x_axis)



		# Model.MCMC(bottle, self.x_axis_list, self.pomerol_models[0], self.pomerol_models[1], DoLP_mode="rayleigh", AoLP_mode="para", ISO=False, least_square_factors = [1, 1, 1])


		# best_total_model, total_diff, min_total_diff, best_total_param, best_total_index = Model.ManuallyFindBestDirectPola(bottle, self.x_axis_list, self.pomerol_models[0], self.pomerol_models[1], DoLP_mode="rayleigh", AoLP_mode="para", ISO=False, least_square_factors = [1, 1, 1])


		# best_total_model, total_diff, min_total_diff, best_total_param = Model.ManuallyFindBestDirectPola(bottle, self.x_axis_list, self.pomerol_models[0], self.pomerol_models[1], DoLP_mode = "rayleigh", AoLP_mode = "perp", ISO=False)
		# self.pomerol_models.append(best_total_model)


		# if self.AddDirectPola and len(self.pomerol_models) >= 2:
		# 	print("Adding direct polarisation to model...")
		# 	if self.fit_RS_model_to_flux: # Model 0 and 1 are grd and sky. 2 is the direct light only and 3 is the Linear Comb of 0 and 1
		# 		direct_model = self.pomerol_models[2]
		# 		base_model = self.pomerol_models[3]
		# 	else: # Model 0 is the diffusion model. 1 is the direct only model
		# 		base_model = self.pomerol_models[1]
		# 		direct_model = self.pomerol_models[2]
		#
		# 	for dolp_factor in self.AddDirectPola:
		# 		self.pomerol_models.append(base_model.AddDirectPola(direct_model["DDoLP"] * dolp_factor, direct_model["DAoLP"]))
		# 		self.pomerol_models[-1].SetMandatoryColumns()
		# 		self.pomerol_models[-1].SetArbitraryUnits(bottle = bottle)
		
		
		if self.add_model and len(self.pomerol_models) >= 2:
			m1, m2 = self.pomerol_models[0], self.pomerol_models[1]
			print("Adding model together:", np.average(m1["I0"]), np.average(m2["I0"]))
			for a in self.add_model:
			# for a, b in self.add_model:

				sum = a * m1 + (1-a) * m2
				sum.SetMandatoryColumns()
				sum.SetArbitraryUnits(bottle = bottle)
				sum.x_axis = m1.x_axis
				self.pomerol_models.append(sum)

		if self.fit_RS_model_to_flux and len(self.pomerol_models) >= 2:
			print("Fitting J2RAYS-1 models to data...")
			if self.fit_func == "GS":
				best_model = Model.FitModelToFlux(self.x_axis_list, bottle.data['smooth_I0'], self.pomerol_models[0], self.pomerol_models[1])
			elif self.fit_func == "GSK":
				best_model = Model.FitModelToFluxPlusIso(self.x_axis_list, bottle.data['smooth_I0'], self.pomerol_models[0], self.pomerol_models[1])

			best_model.SetMandatoryColumns()
			best_model.SetArbitraryUnits(bottle = bottle)
			self.pomerol_models.append(best_model)

		if self.addFDA_to_model != []:
			max_data, min_data = np.max(bottle.data['smooth_I0']), np.min(bottle.data['smooth_I0'])
			ratio_data = max_data / min_data
			for i, m in enumerate(self.pomerol_models):
				try:
					addF, addD, addA = self.addFDA_to_model[bottle.line - 1][i]

					max_model, min_model = np.max(m["I0"]) , np.min(m["I0"])
					ratio_model = max_model / min_model
					shift = (ratio_data * min_model - max_model) / (1 - ratio_data)

					print(f"Adding constant background to model...: {addF} (old data unit), {addF / m.bottle_factor}, {addD}, {addA}")
					print(f"OLD: Min model, Max model, Ratio model, Shift to match data (model unit): {min_model / m.bottle_factor}, {max_model / m.bottle_factor}, {ratio_model}, {shift / m.bottle_factor}")

					m.AddFDA(addF, addD, addA)

					max_model, min_model = np.max(m["I0"]) , np.min(m["I0"])
					ratio_model = max_model / min_model
					shift = (ratio_data * min_model - max_model) / (1 - ratio_data)
					print(f"NEW: Min model, Max model, Ratio model, Shift to match data (model unit): {min_model / m.bottle_factor}, {max_model / m.bottle_factor}, {ratio_model}, {shift / m.bottle_factor}")


				except:
					print("WARNING: DID NOT add FDA to model number:", bottle.line - 1, i)
					pass
					# addF, addD, addA = 0, 0, 0

	def MakeXAxis(self, bottle):

		self.divisor = 1.
		delta = bottle.data['Times'].iloc[-1] - bottle.data['Times'].iloc[0]
		if delta > dt.timedelta(hours=2):
			self.divisor = 3600.
			if   self.langue == "en" : self.xlabel = "Time (hours)"
			elif self.langue == "fr" : self.xlabel = "Durée (heure)"
		elif delta > dt.timedelta(minutes=2):
			self.divisor = 60.
			if   self.langue == "en" : self.xlabel = "Time (minutes)"
			elif self.langue == "fr" : self.xlabel = "Durée (minutes)"
		else:
			if   self.langue == "en" : self.xlabel = "Time (seconds)"
			elif self.langue == "fr" : self.xlabel = "Durée (secondes)"

		# self.x_axis_list = np.array([t.total_seconds() for t in bottle.all_times_since_start]) / self.divisor
		norm = bottle.data['Times'].iloc[0].total_seconds()
		self.x_axis_list = np.array([t.total_seconds() - norm for t in bottle.data['Times']]) / self.divisor

		if self.use_24h_time_format and bottle.observation_type == "fixed":
			self.x_axis_list = bottle.DateTime("start", format=self.time_format) + bottle.data['Times']
			# day = self.x_axis_list[0].day
			# GetHour = lambda t: t.hour + t.minute / 60. + t.second / 3600.
			# GetDay = lambda d: d.day
			# self.x_axis_list = np.where(np.vectorize(GetDay)(self.x_axis_list) > day, np.vectorize(GetHour)(self.x_axis_list)+24, np.vectorize(GetHour)(self.x_axis_list))
			# self.x_axis_list = np.array([GetHour(t) for t in self.x_axis_list])
			self.xlabel = self.time_label

			if self.divisor == 3600.:
				self.xaxis_time_format = '%H:%M'
			else:
				self.xaxis_time_format = '%H:%M:%S'

			# plt.rcParams['axes.xmargin'] = 0


		self.x_time_list = self.x_axis_list

		if self.xaxis_azimut and bottle.observation_type == "fixed_elevation_continue_rotation":
			### Treat continue rotations differently to have all rotations appear the same length. x axis is coded by azimuth, not time. So if the rotations are done by hand, they will appear to have all the same length!
			self.x_axis_list = np.array(())
			self.x_time_ticks_label = [bottle.DateTime("start", format=self.time_format)]
			self.x_time_ticks_pos = [0]
			for ir in range(bottle.nb_continue_rotation):
				start 	= dt.timedelta(minutes=bottle.continue_rotation_times[ir])
				end 	= dt.timedelta(minutes=bottle.continue_rotation_times[ir+1])
				nb_points = len([t for t in bottle.data['Times'] if start <= t <= end ])

				self.x_axis_list = np.append(self.x_axis_list, np.linspace(360 * ir, 360 * (ir+1), nb_points))

				self.x_time_ticks_pos.append(360 * (ir+1))
				# self.x_time_ticks_pos.append(self.x_time_ticks_pos[-1] + 360)
				self.x_time_ticks_label.append(bottle.DateTime("start", format=self.time_format, delta=end))

			self.x_time_ticks_label = [x.strftime("%H:%M:%S") for x in self.x_time_ticks_label]

			# print(len(self.x_axis_list), len(bottle.data['Times']))
			self.xlabel = "Azimuth"

			if bottle.nb_continue_rotation <= 2:
				self.x_axis_ticks_pos = np.arange(0, 360 * bottle.nb_continue_rotation, 45)
				self.x_axis_ticks_label = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"] *  bottle.nb_continue_rotation
			else:
				self.x_axis_ticks_pos = np.arange(0, 360 * bottle.nb_continue_rotation, 90)
				self.x_axis_ticks_label = ["N", "E", "S", "W"] *  bottle.nb_continue_rotation


		else:
			# delta = bottle.data['Times'].iloc[-1] - bottle.data['Times'].iloc[0]
			# if delta > dt.timedelta(hours=2):
			# 	self.divisor = 3600.
			# 	self.xlabel = "Time (hours)"
			# elif delta > dt.timedelta(minutes=2):
			# 	self.divisor = 60.
			# 	self.xlabel = "Time (minutes)"
			# else:
			# 	self.xlabel = "Time (seconds)"
			#
			# # self.x_axis_list 	   = np.array([t.total_seconds() for t in bottle.all_times_since_start]) / self.divisor
			# norm = bottle.data['Times'].iloc[0].total_seconds()
			# # self.x_axis_list = np.array([t.total_seconds() - norm for t in bottle.data['Times']]) / self.divisor
			# self.x_axis_list = np.array([t.total_seconds() - norm for t in bottle.data['Times']]) / self.divisor

			if self.xaxis_azimut and bottle.observation_type == "fixed_elevation_discrete_rotation":
				self.xlabel = "Azimuth (°)"
				self.x_axis_ticks_pos = bottle.discrete_rotation_times[::9]
				self.x_axis_ticks_label = np.round(bottle.discrete_rotation_azimuts * RtoD % 360)[::9]
				self.x_axis_ticks_label = [int(x) for x in self.x_axis_ticks_label]

				self.x_time_ticks_pos = []
				self.x_time_ticks_label = []
				for it, t in enumerate(self.x_axis_ticks_pos):
					self.x_time_ticks_pos.append(t)
					self.x_time_ticks_label.append(bottle.DateTime("start", format=self.time_format) + dt.timedelta(minutes=t))

				self.x_time_ticks_label = [x.strftime("%H:%M:%S") for x in self.x_time_ticks_label]

		if self.max_error_bars:
			self.error_step = len(self.x_axis_list) // self.max_error_bars + 1
		else:
			self.error_step = 1

	def MakeFigure(self):
		if self.mode == "default":
			self.f1, (self.ax1, self.ax2, self.ax3) = plt.subplots(3, sharex=True, figsize=(16, 8))
			# f1, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True, figsize=(16, 8))

			self.ax1_lines = []
			self.ax2_lines = []
			self.ax3_lines = []

			self.ax1.set_ylabel("Flux (AU)")
			self.ax2.set_ylabel("DoLP (%)")
			self.ax3.set_ylabel("AoLP (°)")
			# ax4.set_ylabel("Temperature (Celsius)")

			# ax1.set_xlabel(self.xlabel)
			# ax2.set_xlabel(self.xlabel)
			self.ax3.set_xlabel(self.xlabel)

			if self.show_grid_lines:
				for ax in [self.ax1, self.ax2, self.ax3]:
					ax.set_axisbelow(True)
					ax.grid(True, zorder=-10)

			if self.use_24h_time_format:
				# self.ax3.xaxis.set_major_locator(mpl.dates.DayLocator())
				# self.ax3.xaxis.set_minor_locator(mpl.dates.HourLocator(range(0, 25, 1)))
				# self.ax3.xaxis.set_major_formatter(mpl.dates.DateFormatter('%Y-%m-%d'))
				self.f1.autofmt_xdate()

	def SetGraphParameter(self, bottle, comp = False):
		"""
		Use the attributes defined here to combine your plots, hide them, add error bars or other instrument data.
		For a color control, go to the SetColors() function.
		"""
		self.langue = "en" #Language used for the plots. "fr" or "en"

		self.marker_size = 1 # Size of the points used for all plots
		self.single_star_size = self.marker_size*5 # When plotting the apparent angle of B AoBapp or the light pollution Rayleighj angle AoRD, control the size of the ztar markers.
		self.marker = "." #Linestyle of the polarisation parameters. "." or "none" for point cloud, "-" or "solid" for solid lines.  Refer to https://matplotlib.org/stable/gallery/lines_bars_and_markers/linestyles.html or https://www.geeksforgeeks.org/linestyles-in-matplotlib-python/ for more.
		self.linestyle = "none" #Linestyle of the polarisation parameters. "none" for point cloud, "-" or "solid" for solid lines.  Refer to https://matplotlib.org/stable/gallery/lines_bars_and_markers/linestyles.html or https://www.geeksforgeeks.org/linestyles-in-matplotlib-python/ for more.
		self.show_Ipola = False # Show the plot for I0 * DoLP, i.e. the flux of polarized light
		self.show_Iref = False # For CarmenCru only. Show the reference channel with no polarizing lens

		self.show_time = not comp #Don't touch!
		self.time_format = "LT" #LT or UT. Self explainatory. Control the time format of the x-axis
		self.time_label = "LT" # The title of the time x-axis
		self.use_24h_time_format = 1

		self.show_raw_data = 1 and len(self.comp_bottles) == 0 # Show the data with no slidding average. All rotations of the polarizing filter. In black
		self.show_smooth_data = 1 # Show smoothed data (averaged over the time window defined in the input file)

		self.show_error_bars 		= 1 # Show error bars for the raw cru data. in grey
		self.show_smooth_error_bars = 1 # Show error bars for the smooth cru data. in grey
		self.max_error_bars = 10000 #If there are too many data, takes very long to plot error bars. If != 0, will plot max_error_bars error bars once every x points.

		self.show_avg_I0   = False #Show the flux average over the whole observation
		self.show_avg_Iref = False #Only for Carmen Cru. Show the reference flux average over the whole observation
		self.show_avg_DoLP = False #Show the DoLP for a smoothed data point obtained by averaging all data over the whle observation. It is different from the DoLP average because it is done on the V, Vcos and Vsin values.
		self.show_avg_AoLP = False #Show the AoLP for a smoothed data point obtained by averaging all data over the whle observation. It is different from the AoLP average because it is done on the V, Vcos and Vsin values.

		self.xaxis_azimut = True #only for almucantars. If True will show the azimut on the xaxis instead of time .

		self.show_legend = False #Show the legend bow over the plots

		self.show_correlations = 0

		self.show_allsky = False # If available, plot the allsky camera flux over the cru flux.

		self.show_plip = False

		self.show_eiscat = 0 # If available, plot the eiscat Ne over the Cru flux. Over things are possible if you want, just search for "PlotEiscat" function (called in the Mixer.MakePlot()) and have fun :)
		self.eiscat_type = "uhf" #Initally for March 2022 data. Choose the type of hdf5 files containing eiscat data (Possibilities for VHF mode: tromso, sodankyla, kiruna. For UHF mode: uhf, uhf_v)

		self.show_mag_data 	= 1 # If available, plot the magnetometer data. (field strength or its derivative, or orientation)
		self.B_component	= 'Horiz' # B field component to plot (Dec, Horiz, Vert, 'Incl' or 'Total')
		self.show_AoBapp 	= 0 # If True, show the apparent angle of the magnetic field computed in bottle.py from CHAOS model
		self.show_AoRD	 	= 0 # If True, show the AoLP produced by a point source defined in the input.in file.

		self.show_currents	= 1 # If available, show the apparent angle of the equivalent currents from Magnar.
		self.compute_Jup	= "" #False, "para" or "perp". Will add a vertical component to the equivalent current so that it match the AoLP. If para: the AoLP is parallel to the current. If perp, the AoLP is perpendicular to the current.

		self.show_grid_lines = True # Just to have a nicer grpah. self explainatory

		self.make_optional_plots = 0 # If True, will plot a lots of optional plots showing all kind of things. See the end of the __init__() function where it is used.
		self.show_SN = 0 # If make_optional_plots is True, plot the graph of the signal to noise equivalent defined in appendix of (Bosse et al. 2020) or in bottle.GetSmoothLists() as SN(I, DoLP, Period): DoLP * np.sqrt(I * Period) / 2. where I is the flux, DoLP the DoLP and Period the time of the averaging window.

		### The following paramters are used when comparing the data with the POMEROL model.
		self.show_RS_model		= False #Show or not the POMEROL model graphs.
		self.fit_RS_model_to_flux = False
		self.fit_func = "GS" #"GS" for Ground+Sky. "GSK" for Ground+Sky+background

		self.add_model = [] #Will add the 1st two models m1 and m2 together with a linear combination a*m1 + (1-a)m2
		# self.add_model = [(0.001, 0.999), (0.075, 0.925), (0.05, 0.95), (0.1, 0.90)]

		self.AddDirectPola = [] #[0.15]

		self.show_model_list = []

		self.addFDA_to_model = []

		self.adapt_flux_scale 	= True
		self.shift_model = 0
		self.model_colors = ["black", "red", "blue", "magenta", "purple", "xkcd:mustard", "green", "cyan"] * 10
		self.model_symbols = ["*", "+", "x", "1", "2", "3", "4", ".", "s", "<", "^", ">", "X", "D"] * 10



		# tmp_path = "/home/leob/These/Documentation/Mes_articles/POMEROL/new_pictures/FIG8/"
		# self.RS_model_files	= [	tmp_path + "v_grd_only.txt",
		# 						tmp_path + "v_sky_alb.txt",
		# 						tmp_path + "v_sky_albx1+v_grd_onlyx2e-09",
		# 						#tmp_path + "v_sky_albx1+v_direct_only_Bx1"],
		# 						#tmp_path + "m_sky_albx1+m_grd_onlyx2e-09"],
		# 						tmp_path + "b_sky_albx1+b_grd_onlyx5e-10",
		# 						tmp_path + "b_sky_albx1+b_grd_onlyx1e-10"]

		# tmp_path = "/home/leob/These/Documentation/Mes_articles/POMEROL/new_pictures/FIG7/"
		# self.RS_model_files	= [	tmp_path + "m_rot_e45_grd_only.txt",
		# 						tmp_path + "m_rot_e45_sky_albedo.txt",
		# 						tmp_path + "m_rot_e45_grd_onlyx1+m_rot_e45_sky_albedox12000000000.0",
		# 						tmp_path + "m_rot_e45_grd_onlyx1+m_rot_e45_sky_albedox12000000000.0"]

		# tmp_path = "/home/leob/These/Documentation/Mes_articles/POMEROL/new_pictures/FIG6/"
		# self.RS_model_files	= [	tmp_path + "o_grd_only_norm0.txt",
		# 						tmp_path + "o_grd_only_norm0.txt"]#,
								# tmp_path + "t_sky_alb.txt"]

		tmp_path = "/home/leob/These/Documentation/Mes_articles/POMEROL/new_pictures/FIG5/"
		# self.RS_model_files	= [	[tmp_path + "2021_lagorge_v_grd_only.txt", tmp_path + "2021_lagorge_v_grd_only_aero_mar.txt", tmp_path + "2021_lagorge_v_grd_only_FULL_aero_mar.txt", tmp_path + "2021_lagorge_v_grd_only_FULL_aero_mar_4000.txt", tmp_path + "2021_lagorge_v_grd_only_FULL_aero_mar_10000.txt", tmp_path + "2021_lagorge_v_grd_only_FULL_aero_mar_n4000_250nm.txt", tmp_path + "2021_lagorge_v_grd_only_FULL_aero_mar_n4000_200nm.txt", tmp_path + "2021_lagorge_v_grd_only_FULL_aero_mar_n4000_150nm.txt"]]#,#,

		### Bosse et al (POMEROL) Fig 10+
		# self.RS_model_files	= [	[tmp_path + "2021_lagorge_b_grd_only_FULL_aero_mar_n4000_150nm.txt"],
		# 						[tmp_path + "2021_lagorge_o_grd_only_FULL_aero_mar_n4000_150nm.txt"],
		# 						[tmp_path + "2021_lagorge_m_grd_only_FULL_aero_mar_n4000_150nm.txt"],
		# 						[tmp_path + "2021_lagorge_t_grd_only_FULL_aero_mar_n4000_150nm.txt"]]
		# self.addFDA_to_model = [[(10, 0, 0)], [(21, 0, 0)], [(2, 0, 0)], [(20, 0, 0)]] # 20210120_Lagorge_bomt_rot_lente NEW AEROSOLS

		### Bosse et al 2021 Fig 9
		# self.RS_model_files	= [	[tmp_path + "2021_lagorge_v_grd_only_FULL_aero_mar_n4000_150nm.txt"]]
		# self.addFDA_to_model = [[(190, 0, 0)]]
		### Bosse et al 2021 Fig 8
		self.RS_model_files	= [	[tmp_path + "2021_lagorge_v_grd_only.txt", tmp_path + "2021_lagorge_v_grd_only.txt", tmp_path + "2021_lagorge_v_grd_only_FULL_aero_mar_n4000_150nm.txt", tmp_path + "2021_lagorge_v_grd_only_FULL_aero_rur_1000.txt", tmp_path + "2021_lagorge_v_grd_only_FULL_aero_rur_500.txt"]]
		self.addFDA_to_model = [[(0, 0, 0), (67, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0)]]

		# tmp_path = "/home/leob/These/Documentation/Mes_articles/courants/POMEROL/"
		# # self.RS_model_files	= [[tmp_path + "20200224_skibotnSud_v_e66_MAR1.txt",
		# # 							tmp_path + "uni_sky_noDoLP.txt"]]
		# self.RS_model_files	= [[tmp_path + "20200224_skibotnSud_v_e30_MAR1.txt",
		# 						tmp_path + "uni_sky_noDoLP.txt"]]
		# # self.show_model_list = [0, 2]


		# tmp_path = "/home/leob/These/Analysis/results/rayleigh/Sob/"
		# self.RS_model_files = [ [tmp_path + "20191021_m_rot_e45_grd_Noaer.txt"]]
		# self.RS_model_files = [ [tmp_path + "20191022_r_rot_e45_grd_Noaer.txt"]]
		# self.RS_model_files = [ [tmp_path + "20191023_b_rot_e45_grd_Noaer.txt"]]
		# self.RS_model_files = [ [tmp_path + "20191023_v_rot_e45_grd_Noaer.txt"]]
		# self.RS_model_files = [ [tmp_path + "20191022_o_rot_e45_grd_Noaer.txt"]]

		# self.RS_model_files = [ [tmp_path + "20191022_o_rot_e45_grd_Noaer.txt",
		# 						 tmp_path + "20191022_o_rot_e45_grd_aerMAR.txt",
		# 						 # tmp_path + "20191022_o_rot_e45_grd_aerRUR.txt",
		# 						 # tmp_path + "20191022_o_rot_e45_grd_aerURB.txt",
		# 						 # tmp_path + "20191022_o_rot_e45_grd_aerDES2.txt",
		# 						 # tmp_path + "20191022_o_rot_e45_grd_aerDES3.txt",
		# 						 tmp_path + "20191022_o_rot_e45_grd_aerDES.txt",
		# 						 tmp_path + "20191022_o_rot_e45_grd_aerDES4.txt",
		# 						 tmp_path + "20191022_o_rot_e45_grd_aerDES5.txt"],
		# 						[]]
		# self.addFDA_to_model = [[(0, 0, 0),
		# 						 (9.3, 0, 0),
		# 						 # (38.5, 0, 0),
		# 						 # (76.3, 0, 0),
		# 						 # (497, 0, 0),
		# 						 # (429, 0, 0),
		# 						 (25.5, 0, 0),
		# 						 (332, 0, 0),
		# 						 (220, 0, 0)]]  * len(self.RS_model_files)



		# tmp_path = "/home/leob/These/Documentation/Mes_articles/auroral_pola/"
		### 20200227 XYvm e52 Direct Pola
		# self.RS_model_files	= [[tmp_path + "grd_only/o_e52_grd_only_aero_MAR1.txt",
		# 						tmp_path + "sky_only/uni_sky_o_e52_aero_MAR1_albedo.txt"],
		# 						[tmp_path + "grd_only/t_e52_grd_only_aero_MAR1.txt",
		# 						tmp_path + "sky_only/uni_sky_t_e52_aero_MAR1_albedo.txt"],
		# 						[tmp_path + "grd_only/v_e52_grd_only_aero_MAR1.txt",
		# 						tmp_path + "sky_only/20200227_v_e52_sky_only_aero_MAR1_albedo_skySHIFT20.txt"],
		# 						# tmp_path + "sky_only/20200227_v_e52_sky_only_aero_MAR1_albedo_Dpola-0.310.5_perp_rayleigh.txt"],
		# 						[tmp_path + "grd_only/m_e52_grd_only_aero_MAR1.txt",
		# 						tmp_path + "sky_only/20200227_m_e52_sky_only_aero_MAR1_albedo_skySHIFT20.txt"]]#,
		# 						# tmp_path + "sky_only/20200227_v_e52_sky_only_aero_MAR1_albedo_Dpola-0.310.5_perp_rayleigh.txt"]]

		# For turquoise bottle (2nd)
		# self.AddDirectPola = [0.4]
		# self.add_model = [(0.25, 0.75)]
		# self.show_model_list = [-1]
		# self.fit_RS_model_to_flux = True
		# self.fit_func = "GS" #"GSK" for Ground+Sky+background

		# For green bottle (3rd)
		# self.AddDirectPola = [0.4]
		# self.add_model = [(0.35, 0.65)]
		# # self.fit_RS_model_to_flux = True
		# # self.fit_func = "GS" #"GSK" for Ground+Sky+background
		# self.show_model_list = [0, 3, 4]

		# For mauve bottle (4rd)
		# self.AddDirectPola = [0.9]
		# self.add_model = [(0.25, 0.75)]
		# self.show_model_list = [-1]#[0, 2, 3, 4]
		# self.addFDA_to_model = [[(0, 0, 0), (6, 0, 0), (0, 0, 0), (6, 0, 0), (5.2, 0, 0)]]  * len(self.RS_model_files)

		### 20200227 XYvm e52 best fit
		# self.RS_model_files	= [[tmp_path + "grd_only/o_e52_grd_only_aero_MAR1.txt",
		# 						tmp_path + "sky_only/uni_sky_o_e52_aero_MAR1_albedo.txt"],
		# 						[tmp_path + "grd_only/t_e52_grd_only_aero_MAR1.txt",
		# 						tmp_path + "sky_only/uni_sky_t_e52_aero_MAR1_albedo.txt"],
		# 						# tmp_path + "t_e52_grd_only_aero_MAR1x1702.6351068301055+uni_sky_t_e52_aero_MAR1_albedox63627.19976961118"], #MANUAL fit because of the East flux jump. Coeff are : 4335.3242154348345 * 0.39273535777746454 and 25450.87990784447 * 2.5.
		# 						# [tmp_path + "grd_only/t_e52_grd_only_aero_MAR1.txt",
		# 						# tmp_path + "sky_only/uni_sky_t_e52_aero_MAR1_albedo.txt"],
		# 						[tmp_path + "grd_only/v_e52_grd_only_aero_MAR1.txt",
		# 						tmp_path + "sky_only/20200227_v_e52_sky_only_aero_MAR1_albedo.txt"],
		# 						[tmp_path + "grd_only/m_e52_grd_only_aero_MAR1.txt",
		# 						tmp_path + "sky_only/20200227_m_e52_sky_only_aero_MAR1_albedo.txt"]]

		### 20200227 XYvm e52 test aerosols
		# self.RS_model_files	= [[tmp_path + "grd_only/o_e52_grd_only_NO_aero.txt",
		# 						# tmp_path + "grd_only/o_e52_grd_only_NO_aero.txt",
		# 						# tmp_path + "grd_only/o_e52_grd_only_aero_arctic.txt",
		# 						# tmp_path + "grd_only/o_e52_grd_only_aero_MAR1.txt",
		# 						tmp_path + "grd_only/o_e52_grd_only_aero_MAR1.txt",
		# 						# tmp_path + "grd_only/o_e52_grd_only_aero_MAR1.txt"],
		# 						# tmp_path + "grd_only/o_e52_grd_only_aero_MAR2.txt",
		# 						tmp_path + "grd_only/o_e52_grd_only_aero_3mid.txt",
		# 						tmp_path + "grd_only/o_e52_grd_only_aero_2high.txt"],
		# 						[tmp_path + "grd_only/t_e52_grd_only_aero_MAR1.txt",
		# 						tmp_path + "grd_only/o_e52_grd_only_aero_MAR1.txt",
		# 						tmp_path + "grd_only/o_e52_grd_only_aero_MAR3.txt"],
		# 						[tmp_path + "grd_only/v_e52_grd_only_aero_MAR1.txt",
		# 						tmp_path + "grd_only/o_e52_grd_only_aero_MAR1.txt",
		# 						tmp_path + "grd_only/o_e52_grd_only_aero_MAR3.txt",
		# 						tmp_path + "grd_only/o_e52_grd_only_aero_3mid.txt",
		# 						tmp_path + "grd_only/o_e52_grd_only_aero_2high.txt"],
		# 						[tmp_path + "grd_only/m_e52_grd_only_aero_MAR1.txt",
		# 						tmp_path + "grd_only/o_e52_grd_only_aero_MAR1.txt",
		# 						tmp_path + "grd_only/o_e52_grd_only_aero_MAR3.txt",
		# 						tmp_path + "grd_only/o_e52_grd_only_aero_3mid.txt",
		# 						tmp_path + "grd_only/o_e52_grd_only_aero_2high.txt"]]
		# self.RS_model_files	= [[tmp_path + "CL/o_e52_grd_only_aero_MAR1x1+uni_sky_o_e52_aero_MAR1_albedox6.0"],
								# [],[],[]]

		# ### 20190307 br rot e45
		# self.RS_model_files	= [[tmp_path + "grd_only/b_grd_only_aero_MAR1.txt",
		# 						tmp_path + "sky_only/20190307_b_sky_only_aero_MAR1_albedo_skySHIFT20.txt"],
		# 						# tmp_path + "sky_only/20190307_b_sky_only_aero_MAR1_albedo_Dpola-0.210_para_Rayleigh.txt",],
		# 						[],[],[]]
		# self.AddDirectPola = [0.25]
		# self.add_model = [(0.15, 0.85)]
		# self.show_model_list = []# [1, 3, 4]
		# self.addFDA_to_model = [[(10, 0, 0), (25, 0, 0), (25, 0, 0), (25, 0, 0), (25, 0, 0), (25, 0, 0), (25, 0, 0)]]  * len(self.RS_model_files)

		### 20190307 mr rot e45
		# self.RS_model_files	= [[tmp_path + "grd_only/m_grd_only_aero_MAR1.txt",
		# 						tmp_path + "sky_only/20190307_m_sky_only_aero_MAR1_albedo_skySHIFT20.txt"],#,
		# 						# tmp_path + "sky_only/20190307_m_sky_only_aero_MAR1_albedo_gazm.txt"],
		# 						# tmp_path + "sky_only/20190307_m_sky_only_aero_MAR1_Dpola-0.110_Rayleigh.txt"],
		# 						[],[],[]]
		# # self.AddDirectPola = [0.175]
		# self.add_model = [(0.075, 0.925)]
		# self.show_model_list = [1, 3, 4]

		# ### 20190307 vr rot e45
		# self.RS_model_files	= [[tmp_path + "grd_only/v_grd_only_aero_MAR1.txt",
		# 						tmp_path + "sky_only/20190307_v_sky_only_aero_MAR1_albedo_skySHIFT20.txt"],
		# 						# tmp_path + "sky_only/20190307_v_sky_only_aero_MAR1_albedo_Dpola-0.210_Rayleigh.txt"],
		# 						[],[],[]]
		# self.AddDirectPola = [0.04]
		# self.add_model = [(0.05, 0.95)]
		# self.show_model_list = [1, 3, 4]

		## 20190307 vr a164 e45
		# self.RS_model_files	= [[tmp_path + "grd_only/v_a164_grd_only_aero_MAR1.txt",
		# 						tmp_path + "sky_only/v_a164_sky_only_aero_MAR1_step1.txt"],
		# 						[],[],[]]
		# self.add_model = [0.09]
		# self.show_model_list = [-3, -2, -1]
		## 20190307 br a164 e45
		# self.RS_model_files	= [[tmp_path + "grd_only/b_a164_grd_only_aero_MAR1.txt",
		# 						tmp_path + "sky_only/b_a164_sky_only_aero_MAR1_step1.txt"],
		# 						[],[],[]]
		# self.add_model = [0.15, 0.292]
		# self.show_model_list = [-4, -3, -2, -1]
		# ## 20190307 mr a164 e45
		# self.RS_model_files	= [[tmp_path + "grd_only/m_a164_grd_only_aero_MAR1.txt",
		# 						tmp_path + "sky_only/m_a164_sky_only_aero_MAR1_step1.txt"],
		# 						[],[],[]]
		# self.add_model = [0.2, 0.438]
		# self.show_model_list = [-4, -3, -2, -1]

		### 20200227 XYvm e52 Direct Pola test
		# self.RS_model_files	= [[],
		# 					[],
		# 					[tmp_path + "grd_only/v_e52_grd_only_aero_MAR1.txt",
		# 					 tmp_path + "sky_only/20200227_v_e52_sky_only_aero_MAR1_albedo.txt"],
		# 					 # tmp_path + "sky_only/20200227_v_e52_sky_only_aero_MAR1_albedo_Dpola201_sin2.txt"],
		# 					[tmp_path + "grd_only/m_e52_grd_only_aero_MAR1.txt",
		# 					 tmp_path + "sky_only/20200227_m_e52_sky_only_aero_MAR1_albedo.txt"]]
							 # tmp_path + "sky_only/20200227_v_e52_sky_only_aero_MAR1_albedo_Dpola201_sin2.txt"]]
		# self.RS_model_files	= [[],
		# 						[],
		# 						[tmp_path + "sky_only/20200227_v_e52_sky_only_aero_MAR1_albedo.txt",
		# 						 tmp_path + "sky_only/20200227_v_e52_sky_only_aero_MAR1_albedo_DpolaNW_cos.txt"],
		# 						[tmp_path + "sky_only/20200227_m_e52_sky_only_aero_MAR1_albedo.txt",
		# 						 tmp_path + "sky_only/20200227_v_e52_sky_only_aero_MAR1_albedo_DpolaNW_cos.txt"]]

		# self.addFDA_to_model = [[(0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0)]]  * len(self.RS_model_files)
		# self.addFDA_to_model = np.array([[(0.3, 0, 0), (2, 0, 0), (20, 0, 0), (18, 0, 0), (21, 0, 0)]])*1.5 # 20200227 rotation XYvm à skibotn pour test d'aerosols fit le median
		# self.addFDA_to_model = np.array([[(0.26, 0, 0), (2.14, 0, 0), (33.2, 0, 0), (30, 0, 0), (43, 0, 0)]]) # 20200227 rotation XYvm à skibotn pour test d'aerosols fit le min
		# self.addFDA_to_model = np.array([[(0, 0, 0),(0.26, 0, 0),(0, 0, 0),(43.3, 0, 0)]]) # MORE 20200227 rotation XYvm à skibotn pour test d'aerosols fit le min
		# self.addFDA_to_model = np.array([[(0, 0, 0), (0, 0, 0), (48, 0, 0)]]) # MORE 20200227 rotation XYvm à skibotn pour test d'aerosols fit le min
		# self.addFDA_to_model = np.array([[(43.3, 0, 0)]]) # 20200227 rotation XYvm à skibotn pour best fit (fit le min)

		# self.addFDA_to_model = [[], [], [], [(0, 0, 0), (0, 0, 0), (0, 0, 0)]] # 20200227 rotation XYvm à skibotn pour test de background dans le mauve (fct:average)

		# self.addFDA_to_model = [[(190, 0, 0)]] #Bosse et al 2021 Fig 9
		# self.addFDA_to_model = [[(0, 0, 0), (66.662415, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0)]] # 20210120 rotation verte à lagorge avec 3 profils d'aerosols et background. Bosse et al 2021 Fig 8

		# self.addFDA_to_model = [[(10, 0, 0)], [(21, 0, 0)], [(2, 0, 0)], [(20, 0, 0)]] # 20210120_Lagorge_bomt_rot_lente NEW AEROSOLS


		# self.SetColors(bottle, comp=comp)

		font = {'size'   : 24}
		matplotlib.rc('font', **font)

		# plt.tight_layout(pad=0, h_pad=None, w_pad=None, rect=None)


	def TriangulateAoLP(self, bottle, comp_bottles):
		
		all_bottles = comp_bottles.copy()
		all_bottles.append(bottle)
		
		for b in all_bottles:
			aolps = np.asarray(b.data['smooth_AoLP'])
			b.data['Pperp'] = [b.observation.GetPolaPlane(a, mode='perp', coord='A') for a in aolps] #Norm of the plane that is perp to the los and the AoLP (cross product)
			b.data['Ppara'] = [b.observation.GetPolaPlane(a, mode='para', coord='A') for a in aolps] #Norm of the plane that is parallel to the AoLP (contains AoLP and los -> norm is the cross product)
		
		# All possible pairings of the instruments. (1, 2)==(2, 1)
		pairings = combinations(all_bottles, 2)

		perp_intersections_uen  = []
		perp_intersections_azel = []
		para_intersections_uen  = []
		para_intersections_azel = []

		# Compute the intersection of the plane for each pair of observation
		for i, pair in enumerate(pairings):
			perp_uen = np.array([np.cross(p0, p1) for p0, p1 in zip(pair[0].data['Pperp'], pair[1].data['Pperp'])])
			perp_uen /= np.linalg.norm(perp_uen)

			para_uen = np.array([np.cross(p0, p1) for p0, p1 in zip(pair[0].data['Ppara'], pair[1].data['Ppara'])])
			para_uen /= np.linalg.norm(para_uen)

			perp_azel = np.array([UENToAzEl(uen) for uen in perp_uen])
			para_azel = np.array([UENToAzEl(uen) for uen in para_uen])

			perp_intersections_uen.append(perp_uen)
			para_intersections_uen.append(para_uen)
			perp_intersections_azel.append(perp_azel)
			para_intersections_azel.append(para_azel)


		if len(all_bottles) > 2:
			perp_diff = []
			para_diff = []
			perp_pairings = combinations(perp_intersections_uen, 2)
			para_pairings = combinations(para_intersections_uen, 2)

			print('Angles between every intersections of the perpendicular planes')
			for i, pair in enumerate(perp_pairings):
				angle_diff = np.arccos(np.dot(pair[0], pair[1]))
				perp_diff.append(angle_diff)
				print(angle_diff*RtoD)

			print('Angles between every intersections of the parallel planes')
			for i, pair in enumerate(para_pairings):
				angle_diff = np.arccos(np.dot(pair[0], pair[1]))
				para_diff.append(angle_diff)
				print(angle_diff*RtoD)

			perp_error = np.average(perp_diff)
			para_error = np.average(para_diff)
			if perp_error > para_error:
				print(f'PARALLEL currents are closer, with an average difference angle of {para_error*RtoD}')
			elif perp_error < para_error:
				print(f'PERPENDICULAR currents are closer, with an average difference angle of {perp_error*RtoD}')
			else:
				print('Parallel and perpendicular have exactly the same average difference! Hard to tell which is better...')

			print('All intersections for perp current')
			# print(perp_intersections_uen)
			for i in range(len(perp_intersections_uen)):
				print(perp_intersections_uen[i])

			print('All intersections for para current')
			# print(para_intersections_uen)
			for i in range(len(para_intersections_uen)):
				print(para_intersections_uen[i])



		# perp_current = np.average(perp_intersections_uen, axis=0)
		# perp_current /= np.linalg.norm(perp_current)
		# print("Average perpendicular current (uen)", perp_current)
		# print("Average perpendicular current (azel)", np.array(UENToAzEl(perp_current))*RtoD)
		
		

		# para_current = np.average(para_intersections_uen, axis=0)
		# para_current /= np.linalg.norm(para_current)
		# print("Average parallel current (uen)", para_current)
		# print("Average parallel current (azel)", np.array(UENToAzEl(para_current))*RtoD)


		fig, axs = plt.subplots(2,1, sharex=True)
		axs[0].plot(perp_intersections_azel[0][:, 0]*RtoD)
		axs[1].plot(perp_intersections_azel[0][:, 1]*RtoD)

		fig, axs = plt.subplots(2,1, sharex=True)
		axs[0].plot(para_intersections_azel[0][:, 0]*RtoD)
		axs[1].plot(para_intersections_azel[0][:, 1]*RtoD)
		
		return perp_intersections_azel
		

	# def SetColors(self, bottle, comp = False):

	# 	self.all_I0_color = "xkcd:black"
	# 	self.smooth_I0_color = "xkcd:red"

	# 	self.raw_error_bars_color = "grey"

	# 	self.smooth_ref_color = "xkcd:green"
	# 	self.AoBapp_color = "xkcd:turquoise"
	# 	self.AoRD_color = "xkcd:hot pink"
	# 	self.currents_color = "xkcd:mustard" #"xkcd:lilac"
	# 	if comp:
	# 		self.AoBapp_color = "xkcd:blue"
	# 		self.AoRD_color = "xkcd:salmon"

	# 	self.EISCAT_color = "xkcd:mustard"
	# 	self.mag_color = "xkcd:mustard"

	# 	self.all_sky_color = "xkcd:orange"

	# 	pola_to_line_color = {	"r": ["xkcd:red","xkcd:orange"],
	# 							"v": ["xkcd:green","xkcd:lime green"],
	# 							"b": ["xkcd:blue","xkcd:bright"],
	# 							"m": ["xkcd:purple","xkcd:lavender"],
	# 							"o": ["xkcd:orange","xkcd:orange"],
	# 							"t": ["xkcd:turquoise", "xkcd:turquoise"],
	# 							"X": ["xkcd:orange","xkcd:orange"],
	# 							"Y": ["xkcd:turquoise", "xkcd:turquoise"]}

	# 	### xkcd color guide: https://xkcd.com/color/rgb/
	# 	print("MIXER FILTER:", bottle.filters)
	# 	if bottle.filters:

	# 		if comp and self.pola_color == bottle.filters[0]: #If the last bottle color is the same as the comp bottle, then use a weird color for the comp
	# 			self.smooth_I0_color = "xkcd:olive"
	# 		else:
	# 			self.pola_color = bottle.filters[0]
	# 			if self.pola_color not in pola_to_line_color:
	# 				self.smooth_I0_color = "xkcd:black"
	# 			else:
	# 				self.smooth_I0_color = pola_to_line_color[self.pola_color][0]

	# 			# if 	 self.pola_color == "r": self.smooth_I0_color = "xkcd:red"
	# 			# elif self.pola_color == "v": self.smooth_I0_color = "xkcd:green"
	# 			# elif self.pola_color == "b": self.smooth_I0_color = "xkcd:blue"
	# 			# elif self.pola_color == "m": self.smooth_I0_color = "xkcd:purple"
	# 			# elif self.pola_color == "o": self.smooth_I0_color = "xkcd:orange"
	# 			# elif self.pola_color == "t": self.smooth_I0_color = "xkcd:turquoise"
	# 			# elif self.pola_color == "X": self.smooth_I0_color = "xkcd:orange"
	# 			# elif self.pola_color == "Y": self.smooth_I0_color = "xkcd:turquoise"
	# 			# else: self.smooth_I0_color = "red"

	# 		self.smooth_error_bars_color = self.smooth_I0_color

	# 		if bottle.instrument_name == "carmen" and bottle.filters[1] != 0 and not bottle.NoVref:
	# 			self.ref_color = bottle.filters[1]
	# 			self.smooth_ref_color = 'black'
	# 			if self.ref_color in pola_to_line_color:
	# 				if self.ref_color == self.pola_color:
	# 					self.smooth_ref_color = pola_to_line_color[self.ref_color][1]
	# 				else:
	# 					self.smooth_ref_color = pola_to_line_color[self.ref_color][0]

	# 			# if 	 self.ref_color == "r" and self.pola_color != "r": self.smooth_ref_color = "xkcd:red"
	# 			# elif self.ref_color == "r" and self.pola_color == "r": self.smooth_ref_color = "xkcd:orange"
	# 			# elif self.ref_color == "v" and self.pola_color != "v": self.smooth_ref_color = "xkcd:green"
	# 			# elif self.ref_color == "v" and self.pola_color == "v": self.smooth_ref_color = "xkcd:lime green"
	# 			# elif self.ref_color == "b" and self.pola_color != "b": self.smooth_ref_color = "xkcd:blue"
	# 			# elif self.ref_color == "b" and self.pola_color == "b": self.smooth_ref_color = "xkcd:bright blue"
	# 			# elif self.ref_color == "m" and self.pola_color != "m": self.smooth_ref_color = "xkcd:purple"
	# 			# elif self.ref_color == "m" and self.pola_color == "m": self.smooth_ref_color = "xkcd:lavender"
	# 			# elif self.ref_color == "o" and self.pola_color != "o": self.smooth_ref_color = "xkcd:orange"
	# 			# elif self.ref_color == "o" and self.pola_color == "o": self.smooth_ref_color = "xkcd:orange"
	# 			# else: self.smooth_ref_color = "green"


	# 	self.all_SN_color = "black"
	# 	self.smooth_SN_color = self.smooth_I0_color


	# def PlotFlux(self, ax, bottle, ax_lines=[]):

	# 	if self.show_raw_data:
	# 		if not self.show_error_bars:
	# 			l_all_I0, = ax.plot(self.x_axis_list, bottle.data['I0'],  color = self.all_I0_color, marker = self.marker, linestyle = self.linestyle, markersize=self.marker_size, label="Intensity", zorder=0)
	# 		else:
	# 			l_all_I0 = ax.errorbar(self.x_axis_list, bottle.data['I0'], yerr = bottle.data['std_I0'],  ecolor=self.raw_error_bars_color, color = self.all_I0_color, marker = self.marker, linestyle = self.linestyle, markersize=self.marker_size, label="Intensity", zorder=0, errorevery=self.error_step)


	# 	if self.show_smooth_data and not self.show_smooth_error_bars:
	# 		# l_smooth_I0, = ax.plot(self.x_axis_list, bottle.data['smooth_I0'], fmt=".", color = self.smooth_I0_color, linestyle = 'none', markersize=self.marker_size, label="Smooth Intensity (" + str(bottle.smoothing_factor) + " " + str(bottle.smoothing_unit)[:3] + ")", zorder=1)
	# 		y = bottle.data['smooth_I0']
	# 		if self.show_Ipola:
	# 			y *= bottle.data['smooth_DoLP'] / 100.
	# 		l_smooth_I0, = ax.plot(self.x_axis_list, y, color = self.smooth_I0_color, marker = self.marker, linestyle = self.linestyle, markersize=self.marker_size, label="Smooth Intensity (" + str(bottle.smoothing_factor) + " " + str(bottle.smoothing_unit)[:3] + ")", zorder=1)
	# 	elif self.show_smooth_data:
	# 		# l_smooth_I0 = ax.errorbar(self.x_axis_list, bottle.data['smooth_I0'], yerr = bottle.data['std_smooth_I0'], fmt=".", color = self.smooth_I0_color, ecolor=self.smooth_error_bars_color, linestyle = 'none', markersize=self.marker_size, label="Smooth Intensity (" + str(bottle.smoothing_factor) + " " + str(bottle.smoothing_unit)[:3] + ")", zorder=1, errorevery=self.error_step)
	# 		y = bottle.data['smooth_I0'] 
	# 		if self.show_Ipola:
	# 			y *= bottle.data['smooth_DoLP'] / 100.
	# 		l_smooth_I0 = ax.errorbar(self.x_axis_list, y, yerr = bottle.data['std_smooth_I0'],  color = self.smooth_I0_color, ecolor=self.smooth_error_bars_color, marker = self.marker, linestyle = self.linestyle, markersize=self.marker_size, label="Smooth Intensity (" + str(bottle.smoothing_factor) + " " + str(bottle.smoothing_unit)[:3] + ")", zorder=1, errorevery=self.error_step)


	# 	if self.show_avg_I0:
	# 		l_avg_I0, = ax.plot(self.x_axis_list,[bottle.I0_average] * len(self.x_axis_list), color = self.smooth_I0_color, label="Avg: " + str(bottle.I0_average)[:4])
	# 		ax_lines.append([l_avg_I0, l_avg_I0.get_label()])

	# 	if self.show_raw_data:
	# 		ax_lines.append([l_all_I0, l_all_I0.get_label()])
	# 	if  self.show_smooth_data:
	# 		ax_lines.append([l_smooth_I0, l_smooth_I0.get_label()])



	# 	# ### Graph the I0 * DoLP line on the first subplot
	# 	# ax4 = ax.twinx()
	# 	# # l_all_IDoLP, = ax4.plot(self.x_axis_list, bottle.data['I0'] * bottle.data['DoLP'], "k.", linestyle = 'none', markersize=self.marker_size, label="All Intensity * DoLP", zorder=2)
	# 	# smooth_IDoLP, = ax4.plot(self.x_axis_list, bottle.data['smooth_I0'] * bottle.data['smooth_DoLP'] / 100. * bottle.data['smooth_AoLP']
	# 	#
	# 	# cd, ".", color = "xkcd:hot pink", linestyle = 'none', markersize=self.marker_size, label="Smooth Intensity * DoLP", zorder=2)
	# 	#
	# 	# ax4.set_ylabel("Ref intensity")



	# 	if not bottle.NoVref and bottle.instrument_name == "carmen" and self.show_Iref:
	# 		ax3 = ax.twinx()
	# 		# ax13.set_xlim(ax1.get_xlim())
	# 		# xticks = plt.xticks()[0] * self.divisor
	# 		# xticks = [t + bottle.time + bottle.head_jump for t in xticks]
	# 		# ax12.set_xticklabels([dt.strftime("%H:%M:%S", dt.localtime(st)) for st in xticks])

	# 		# l_Iref, = ax3.plot(self.x_axis_list[1:], bottle.all_Iref[1:], ".", color = "black", linestyle = 'none', markersize=self.marker_size, label="Ref Intensity", zorder=2)
	# 		# ax_lines.append([l_Iref, l_Iref.get_label()])

	# 		l_smooth_Iref, = ax3.plot(self.x_axis_list, bottle.smooth_Iref, color = self.smooth_ref_color, marker = self.marker, linestyle = self.linestyle, markersize=self.marker_size, label="Smooth Ref Intensity (" + str(bottle.smoothing_factor) + " " + str(bottle.smoothing_unit)[:3] + ")", zorder=2)
	# 		ax_lines.append([l_smooth_Iref, l_smooth_Iref.get_label()])

	# 		if self.show_avg_Iref:
	# 			l_avg_Iref, = ax3.plot(self.x_axis_list, [bottle.Iref_average] * len(self.x_axis_list), color = self.smooth_ref_color, label="Avg Ref Intensity " + str(bottle.Iref_average)[:4], zorder=2)
	# 			ax_lines.append([l_avg_Iref, l_avg_Iref.get_label()])
	# 		ax3.set_ylabel("Ref intensity")


	# def PlotFluxDerivative(self, ax, bottle, ax_lines=[]):
	# 	l_smooth_I0diff, = ax.plot(self.x_axis_list, bottle.I0_diff, color = 'xkcd:red', marker = self.marker, linestyle = self.linestyle, markersize=self.marker_size, label="Flux Derivative (" + str(bottle.smoothing_factor) + " " + str(bottle.smoothing_unit)[:3] + ")", zorder=1)


	# def PlotDoLP(self, ax, bottle, ax_lines=[]):

	# 	# if self.show_error_bars:
	# 	# 	smooth_yerr = bottle.var_smooth_DoLP
	# 	# 	yerr = bottle.var_DoLP
	# 	# 	yerr = None
	# 	# else:
	# 	# 	smooth_yerr = None

	# 	if self.show_raw_data:
	# 		if not self.show_error_bars:
	# 			l_all_DoLP, = ax.plot(self.x_axis_list, bottle.data['DoLP'], color = self.all_I0_color, marker = self.marker, linestyle = self.linestyle, markersize=self.marker_size, label="DoLP", zorder=0)
	# 			ax_lines.append([l_all_DoLP, l_all_DoLP.get_label()])
	# 		else:
	# 			(l_all_DoLP, _, _) = ax.errorbar(self.x_axis_list, bottle.data['DoLP'], yerr = bottle.data['std_DoLP'],  ecolor=self.raw_error_bars_color, color = self.all_I0_color, marker = self.marker, linestyle = self.linestyle, markersize=self.marker_size, label="DoLP", zorder=0, errorevery=self.error_step)

	# 	if self.show_smooth_data and not self.show_smooth_error_bars:
	# 		l_smooth_DoLP, = ax.plot(self.x_axis_list, bottle.data['smooth_DoLP'], color = self.smooth_I0_color, marker = self.marker, linestyle = self.linestyle, markersize=self.marker_size, label="Smooth DoLP (" + str(bottle.smoothing_factor) + " " + str(bottle.smoothing_unit)[:3] + ")", zorder=1)
	# 	elif  self.show_smooth_data:
	# 		(l_smooth_DoLP, _, _) = ax.errorbar(self.x_axis_list, bottle.data['smooth_DoLP'], yerr = bottle.data['std_smooth_DoLP'],  color = self.smooth_I0_color, ecolor=self.smooth_error_bars_color, marker = self.marker, linestyle = self.linestyle, markersize=self.marker_size, label="Smooth DoLP (" + str(bottle.smoothing_factor) + " " + str(bottle.smoothing_unit)[:3] + ")", zorder=1, errorevery=self.error_step)

	# 	ax.set_ylim(bottom = 0)

	# 	if self.show_raw_data:
	# 		ax_lines.append([l_all_DoLP, l_all_DoLP.get_label()])
	# 	if  self.show_smooth_data:
	# 		ax_lines.append([l_smooth_DoLP, l_smooth_DoLP.get_label()])

	# 	# if bottle.location.lower() == "skibotn" and bottle.filters == "br" and bottle.DateTime().date() == dt.date(2019, 3, 7):
	# 	# 	ax.set_ylim((0, 5))

	# 	if self.show_avg_DoLP:
	# 		l_avg_DoLP, = ax.plot(self.x_axis_list,[bottle.DoLP_average] * len(self.x_axis_list), color = self.smooth_I0_color, label="Avg: " + str(bottle.DoLP_average)[:4])
	# 		ax_lines.append([l_avg_DoLP, l_avg_DoLP.get_label()])



	# def PlotB(self, ax, bottle, component, ax_lines=[], derivative=False):
	# 	if not self.mag_data:
	# 		return
	# 	# ax1 = ax.twinx()

	# 	# available compinents:
	# 	# "Dec"
	# 	# "Horiz"
	# 	# "Vert"

	# 	# print(self.mag_data.times[0], self.mag_data.times[-1])
	# 	# print(self.mag_data.GetNormTimes(self.divisor)[0], self.mag_data.GetNormTimes(self.divisor)[-1])
	# 	label = component
	# 	if component in ['Dec' or 'Incl']:
	# 		label += ' (deg)'
	# 	else:
	# 		label += ' (nT)'

	# 	if not derivative:
	# 		t, d = self.mag_data.GetComponent(component, self.divisor, use_datetime = self.use_24h_time_format)
	# 	else:
	# 		t, d = self.mag_data.GetDerivative(component, self.divisor, use_datetime = self.use_24h_time_format)
	# 	l_mag_data, = ax.plot(t, d, "orange", label=label)
	# 	# l_mag_data, = ax.plot(t, d, self.mag_color, label="dB/dt (nT/s)", zorder=2, linewidth=2)
	# 	# print(t, d)
	# 	# Md, md = max(d), min(d)
	# 	# ax.set_ylim(min(-20, md, -abs(Md)), max(20, Md, abs(md)))
	# 	ax_lines.append([l_mag_data, l_mag_data.get_label()])
	# 	ax.set_ylabel(label)
	# 	# ax.set_ylabel("dB/dt (nT/s)")



	# def PlotAllSky(self, ax, bottle, ax_lines=[]):

	# 	if self.allsky_data_available and self.show_allsky:
	# 		ax4 = ax.twinx()
	# 		offset = 10
	# 		# new_fixed_axis = ax4.get_grid_helper().new_fixed_axis
	# 		# ax4.axis["right"] = new_fixed_axis(loc="right", axes=par2, offset=(offset, 0))
	# 		#
	# 		# ax4.axis["right"].toggle(all=True)

	# 		l_ASI, = ax4.plot(self.allsky_data.GetNormTimes(bottle.DateTime("start", format=self.time_label), self.divisor), self.allsky_data.brightness, "*", color = "orange", linestyle = 'none', markersize=self.marker_size, label="AllSKy Imager", zorder=2)
	# 		ax_lines.append([l_ASI, l_ASI.get_label()])
	# 		# print(ax4.get_yticklabels())
	# 		# ax4.set_yticklabels([l.get_text() for l in ax4.get_yticklabels()], horizontalalignment = "left")
	# 		# ax4.tick_params(direction='in', labelright=True, pad=-5)
	# 		ax4.set_ylabel("ASI Brightness")


	# def PlotEiscat(self, bottle, parameter, ax, ax_lines):

	# 	# Parameters for VHF MAD6400_beata_60 hdf5 files
	# 	# dtype=[('year', '<i8'), ('month', '<i8'), ('day', '<i8'), ('hour', '<i8'), ('min', '<i8'), ('sec', '<i8'), ('recno', '<i8'), ('kindat', '<i8'), ('kinst', '<i8'), ('ut1_unix', '<f8'), ('ut2_unix', '<f8'), ('azm', '<f8'), ('elm', '<f8'), ('hsa', '<f8'), ('tfreq', '<f8'), ('power', '<f8'), ('systmp', '<f8'), ('range', '<f8'), ('gdalt', '<f8'), ('ne', '<f8'), ('ti', '<f8'), ('tr', '<f8'), ('co', '<f8'), ('vo', '<f8'), ('pm', '<f8'), ('po+', '<f8'), ('dne', '<f8'), ('dti', '<f8'), ('dtr', '<f8'), ('dco', '<f8'), ('dvo', '<f8'), ('dpm', '<f8'), ('dpo+', '<f8'), ('dwn', '<f8'), ('ddc', '<f8'), ('gfit', '<f8'), ('chisq', '<f8')])

	# 	# Parameters for UHF MAD6502_beata_ant hdf5 files
	# 	# dtype=[('year', '<i8'), ('month', '<i8'), ('day', '<i8'), ('hour', '<i8'), ('min', '<i8'), ('sec', '<i8'), ('recno', '<i8'), ('kindat', '<i8'), ('kinst', '<i8'), ('ut1_unix', '<f8'), ('ut2_unix', '<f8'), ('hsa', '<f8'), ('tfreq', '<f8'), ('azm', '<f8'), ('elm', '<f8'), ('power', '<f8'), ('systmp', '<f8'), ('range', '<f8'), ('gdalt', '<f8'), ('ne', '<f8'), ('ti', '<f8'), ('tr', '<f8'), ('co', '<f8'), ('vo', '<f8'), ('pm', '<f8'), ('po+', '<f8'), ('dne', '<f8'), ('dti', '<f8'), ('dtr', '<f8'), ('dco', '<f8'), ('dvo', '<f8'), ('dpm', '<f8'), ('dpo+', '<f8'), ('dwn', '<f8'), ('ddc', '<f8'), ('gfit', '<f8'), ('chisq', '<f8')])

	# 	# Parameters for UHF MAD6502_beata_V225  hdf5 files
	# 	# dtype=[('year', '<i8'), ('month', '<i8'), ('day', '<i8'), ('hour', '<i8'), ('min', '<i8'), ('sec', '<i8'), ('recno', '<i8'), ('kindat', '<i8'), ('kinst', '<i8'), ('ut1_unix', '<f8'), ('ut2_unix', '<f8'), ('gdalt', '<f8'), ('gdlat', '<f8'), ('glon', '<f8'), ('vi1', '<f8'), ('vi2', '<f8'), ('vi3', '<f8'), ('dvi1', '<f8'), ('dvi2', '<f8'), ('dvi3', '<f8')])

	# 	if not self.eiscat_data.valid:
	# 		return

	# 	time_format = "delta"
	# 	time_delta = dt.timedelta(hours=0)
	# 	if self.use_24h_time_format:
	# 		time_format = "datetime"
	# 		if self.time_format == "LT":
	# 			time_delta = dt.timedelta(hours=1)

	# 	# if self.eiscat_type == "uhf_v":
	# 	# 	parameter = "AoVi"

	# 	unit = r"m$^{-3}$"
	# 	if 'v' in parameter:
	# 		unit = r'm/s'
	# 	elif parameter == "AoVi":
	# 		unit = r'deg'

	# 	altitude = 110 # bottle.GetAltitude()
	# 	if self.langue == "en":
	# 		label = parameter.replace("_", "").capitalize() + " at " + str(altitude) + f" km ({unit})"
	# 	else:
	# 		label = parameter.replace("_", "").capitalize() + " à " + str(altitude) + f" km ({unit})"

	# 	if self.eiscat_type == "uhf_v" and parameter.lower() == "aovi":
	# 		t, d, e = self.eiscat_data.GetUHFApparentAngle(bottle, time_divisor = self.divisor, time_format = time_format)
	# 		d *= RtoD
	# 		# l_eiscat_data, capline, barlinecol = ax.errorbar(t+time_delta, d*RtoD, color = self.EISCAT_color, fmt = ".", label = label, zorder=2, markersize=self.marker_size*2)

	# 	else:
	# 		t, d, e = self.eiscat_data.GetParameter(parameter, altitude, time_divisor = self.divisor, time_format = time_format)
		
	# 	l_eiscat_data, capline, barlinecol = ax.errorbar(t + time_delta, d, yerr = e, color = self.EISCAT_color, fmt = ".", label = label, zorder=2, markersize=self.marker_size*2)

	# 	# Md, md = max(d), min(d)

	# 	if bottle.location.lower() == "skibotn" and bottle.filters == "vr" and bottle.DateTime().date() == dt.date(2019, 3, 7):
	# 		ax.set_ylim((10**11, 4.1*10**11))
	# 	elif bottle.location.lower() == "skibotn" and bottle.filters == "vr" and bottle.DateTime().date() == dt.date(2019, 3, 8):
	# 		ax.set_ylim((-3*10**10, 5*10**10))

	# 	ax_lines.append([l_eiscat_data, label])
	# 	ax.set_ylabel(label)


	# # ### Plot a second eiscat parameter onto the bottom graph
	# # if self.eiscat_type != "uhf_v":
	# # 	self.ax221 = self.ax3.twinx()
	# # 	parameter = "vo"
	# # else:
	# # 	parameter = "AoVi"
	# #
	# # altitude = bottle.GetAltitude()
	# # if self.langue == "en":
	# # 	label = parameter.replace("_", "") + " at " + str(altitude) + r" km (m/s)"
	# # else:
	# # 	label = parameter.replace("_", "") + " à " + str(altitude) + r" km (m/s)"
	# #
	# # if self.eiscat_type != "uhf_v":
	# # 	t, d, e = self.eiscat_data.GetParameter(parameter, altitude, time_divisor = self.divisor, time_format = time_format)
	# # 	l_eiscat_data, capline, barlinecol = self.ax221.errorbar(t+time_delta, d, yerr = e, color = self.EISCAT_color, fmt = ".", label = label, zorder=2)
	# # 	# Md, md = max(d), min(d)
	# # 	if bottle.location.lower() == "skibotn" and bottle.filters == "mr" and bottle.DateTime().date() == dt.date(2019, 3, 8):
	# # 		self.ax221.set_ylim((-180, 48))
	# # 	# elif bottle.location.lower() == "skibotn" and bottle.filters == "vr" and bottle.DateTime().date() == dt.date(2019, 3, 8):
	# # 	# 	self.ax221.set_ylim((-3*10**10, 5*10**10))
	# # 	self.ax221.set_ylabel(label)
	# # else:
	# # 	t, d, e = self.eiscat_data.GetUHFApparentAngle(bottle, time_divisor = self.divisor, time_format = time_format)
	# # 	l_eiscat_data, capline, barlinecol = self.ax3.errorbar(t+time_delta, d*RtoD, color = self.EISCAT_color, fmt = ".", label = label, zorder=2, markersize=self.marker_size*2)
	# #
	# # 	self.ax3_lines.append([l_eiscat_data, label])


	# def PlotAoLP(self, ax, bottle, ax_lines=[]):

	# 	if self.show_raw_data:
	# 		if not self.show_error_bars:
	# 			l_all_AoLP, = ax.plot(self.x_axis_list, bottle.data['AoLP'] * RtoD, color = self.all_I0_color, marker = self.marker, linestyle = self.linestyle, markersize=self.marker_size, label="AoLP", zorder=0)

	# 	if self.show_smooth_data and not self.show_smooth_error_bars:
	# 		l_smooth_AoLP, = ax.plot(self.x_axis_list, bottle.data['smooth_AoLP'] * RtoD, color = self.smooth_I0_color, marker = self.marker, linestyle = self.linestyle, markersize=self.marker_size, label="Smooth AoLP (" + str(bottle.smoothing_factor) + " " + str(bottle.smoothing_unit)[:3] + ")", zorder=1)

	# 	if (self.show_error_bars and self.show_raw_data) or self.show_smooth_error_bars:
	# 		# ax.fill_between(self.x_axis_list, bottle.data['AoLP'] * RtoD - bottle.data['std_AoLP'], bottle.data['AoLP'] * RtoD + bottle.data['std_AoLP'], color = "grey")
	# 		# ax.fill_between(self.x_axis_list, bottle.smooth_AoLP_lower * RtoD, bottle.smooth_AoLP_upper * RtoD, color = "yellow", alpha = 0.5)

	# 		ls_all = dict()
	# 		ls_smooth = dict()
	# 		if bottle.graph_angle_shift == 1:
	# 			min, mid, max = 90, 180, 270#0, 90, 180
	# 		elif bottle.graph_angle_shift == 0:
	# 			min, mid, max = -90, 0, 90

	# 		for x, s_error, s_angle, error, angle in zip(self.x_axis_list[::self.error_step], bottle.data['std_smooth_AoLP'].iloc[::self.error_step]*RtoD, bottle.data['smooth_AoLP'].iloc[::self.error_step]*RtoD, bottle.data['std_AoLP'].iloc[::self.error_step]*RtoD, bottle.data['AoLP'].iloc[::self.error_step]*RtoD):

	# 			if (self.show_error_bars and self.show_raw_data):
	# 				if angle > mid and angle + error > max:
	# 					temp = angle + error - 180
	# 					ls_all.update({x:[min, temp]})
	# 				elif angle < mid and angle - error < min:
	# 					temp = angle - error + 180
	# 					ls_all.update({x:[max, temp]})
	# 			if (self.show_smooth_data and self.show_smooth_error_bars):
	# 				if s_angle > mid and s_angle + s_error > max:
	# 					temp = s_angle + s_error - 180
	# 					ls_smooth.update({x:[min, temp]})
	# 				elif s_angle < mid and s_angle - s_error < min:
	# 					temp = s_angle - s_error + 180
	# 					ls_smooth.update({x:[max, temp]})

	# 		ax.set_ylim(min, max)
	# 		if (self.show_error_bars and self.show_raw_data):
	# 			l_all_AoLP = ax.errorbar(self.x_axis_list, bottle.data['AoLP'] * RtoD, yerr = bottle.data['std_AoLP'] * RtoD,  ecolor=self.raw_error_bars_color, color = self.all_I0_color, marker = self.marker, linestyle = self.linestyle, markersize=self.marker_size, label="AoLP", zorder=0, errorevery=self.error_step)
	# 			for i, a in ls_all.items():
	# 				ax.vlines(i, a[0], a[1], colors = self.raw_error_bars_color, zorder=0) #"green", linewidth=5)#

	# 		if self.show_smooth_data and self.show_smooth_error_bars:
	# 			l_smooth_AoLP = ax.errorbar(self.x_axis_list, bottle.data['smooth_AoLP'] * RtoD, yerr = bottle.data['std_smooth_AoLP'] * RtoD,  color = self.smooth_I0_color, ecolor=self.smooth_error_bars_color, marker = self.marker, linestyle = self.linestyle, markersize=self.marker_size, label="AoLP", zorder=1, errorevery=self.error_step)
	# 			# plt.errorbar(list(ls.keys()), [-90, 90], yerr=list(ls.values()), fmt='C0 ')
	# 			for i, a in ls_smooth.items():
	# 				ax.vlines(i, a[0], a[1], colors = self.smooth_I0_color, zorder=1) #"green", linewidth=5)#

	# 	if (self.show_error_bars and self.show_raw_data):
	# 		ax_lines.append([l_all_AoLP, l_all_AoLP.get_label()])

	# 	if  self.show_smooth_data and self.show_smooth_error_bars:
	# 		ax_lines.append([l_smooth_AoLP, l_smooth_AoLP.get_label()])

	# 	if self.show_avg_AoLP:
	# 		l_avg_AoLP, = ax.plot(self.x_axis_list,[bottle.AoLP_average * RtoD] * len(self.x_axis_list), color = self.smooth_I0_color, label="Avg: " + str(bottle.AoLP_average * RtoD)[:4])
	# 		ax_lines.append([l_avg_AoLP, l_avg_AoLP.get_label()])


	# def PlotCurrents(self, bottle):

	# 	if self.show_currents and self.eq_currents.valid:
	# 		time_format = "delta"
	# 		time_delta = 0
	# 		if self.use_24h_time_format:
	# 			time_format = "datetime"
	# 			if self.time_format == "LT":
	# 				time_delta = dt.timedelta(hours=1)

	# 		self.ax12 = self.ax1.twinx()
	# 		if self.show_eiscat and self.eiscat_data.valid: # Shift the y axis to not overlap with eiscat plot
	# 			self.ax12.spines["right"].set_position(("axes", 1.075))

	# 		# print(self.divisor, time_format,  time_delta)
	# 		# print(type(self.eq_currents.GetNormTimes(self.divisor, format=time_format)[0]))
	# 		# print(type(time_delta))
	# 		l_AoJnorm, = self.ax12.plot(self.eq_currents.GetNormTimes(self.divisor, format=time_format) + time_delta, self.eq_currents.data["J_norm"], "*", color = self.currents_color, label="J_norm")

	# 		self.ax12.set_ylabel(r"J (A/m)")

	# 		# self.ax23 = self.ax2.twinx()
	# 		# l_AoJlos, = self.ax23.plot(self.eq_currents.GetNormTimes(self.divisor), self.eq_currents.data["AoJlos"] * RtoD, "*", color = self.currents_color, label="AoJlos")
	# 		# # l_AoJlos, = self.ax23.plot(self.eq_currents.GetNormTimes(self.divisor), np.sin(self.eq_currents.data["AoJlos"]) * self.eq_currents.data["J_norm"], color = self.currents_color, label="AoJlos")
	# 		# self.ax2_lines.append([l_AoJlos, l_AoJlos.get_label()])

	# 		AoJapp = self.eq_currents.data["AoJapp"]
	# 		# if bottle.graph_angle_shift == 1:
	# 		# 	AoJapp = SetAngleBounds(self.eq_currents.data["AoJapp"], 0, np.pi)
	# 		AoJappperp = self.eq_currents.data["AoLP"]
	# 		# if bottle.graph_angle_shift == 1:
	# 		# 	AoJappperp = SetAngleBounds(self.eq_currents.data["AoLP"], 0, np.pi)

	# 		AoJapp = SetAngleBounds(AoJapp, -np.pi/2 + np.pi/2 * bottle.graph_angle_shift, np.pi/2 + np.pi/2 * bottle.graph_angle_shift)
	# 		AoJapp90 = SetAngleBounds(AoJapp+np.pi/2, -np.pi/2 + np.pi/2 * bottle.graph_angle_shift, np.pi/2 + np.pi/2 * bottle.graph_angle_shift)
	# 		AoJappperp = SetAngleBounds(AoJappperp, -np.pi/2 + np.pi/2 * bottle.graph_angle_shift, np.pi/2 + np.pi/2 * bottle.graph_angle_shift)

	# 		l_AoJapp, = self.ax3.plot(self.eq_currents.GetNormTimes(self.divisor, format=time_format) + time_delta, AoJapp * RtoD, "*", color = self.currents_color, label="AoJapp")
	# 		l_AoJapp, = self.ax3.plot(self.eq_currents.GetNormTimes(self.divisor, format=time_format) + time_delta, AoJapp90 * RtoD, "+", color = "xkcd:pink", label="AoJapp")

	# 		l_AoJapp, = self.ax3.plot(self.eq_currents.GetNormTimes(self.divisor), AoJappperp * RtoD, "*", color = "red", label="AoJapp_perp")

	# 		# self.ax3_lines.append([l_AoJapp, l_AoJapp.get_label()])

	# 	# self.ax3.set_ylim(45, 60)


	# def PlotAoRD(self, ax, bottle, ax_lines=[]):

	# 	if bottle.observation_type == "fixed":
	# 		if bottle.AoRD is not False and self.show_AoRD:
	# 			if bottle.graph_angle_shift == 1:
	# 				AoRD = SetAngleBounds(bottle.AoRD, 0, np.pi)
	# 			elif bottle.graph_angle_shift == 0:
	# 				AoRD = SetAngleBounds(bottle.AoRD, -np.pi/2, np.pi/2)
	# 			l_AoRD, = ax.plot(self.x_axis_list,[bottle.AoRD * RtoD] * len(self.x_axis_list), linewidth=3, color=self.AoRD_color, label="AoRD: " + str(bottle.AoRD*RtoD)[:5], zorder=2)
	# 			ax_lines.append([l_AoRD, l_AoRD.get_label()])
	# 			# l54, = ax.plot(self.x_axis_list,[bottle.AoRD_ortho * RtoD] * len(self.x_axis_list), ":g", linewidth=self.marker_size, label="AoRD ortho: " + str(bottle.AoRD_ortho*RtoD)[:5])
	# 		if (bottle.AoBapp and bottle.AoBlos) is not False and self.show_AoBapp:
	# 			if bottle.graph_angle_shift == 1:
	# 				AoBapp = SetAngleBounds(bottle.AoBapp, 0, np.pi)
	# 			elif bottle.graph_angle_shift == 0:
	# 				AoBapp = SetAngleBounds(bottle.AoBapp, -np.pi/2, np.pi/2)
	# 			l_AoBapp, = ax.plot(self.x_axis_list, [AoBapp * RtoD] * len(self.x_axis_list), linewidth=3, color = self.AoBapp_color, label="AoBapp: " + str(bottle.AoBapp*RtoD)[:5], zorder=2)
	# 			ax_lines.append([l_AoBapp, l_AoBapp.get_label()])
	# 			# l_AoBapp_ortho, = ax.plot(self.x_axis_list,[bottle.AoBapp_ortho * RtoD] * len(self.x_axis_list), ":b", linewidth=self.marker_size, label="AoBapp ortho: " + str(bottle.AoBapp_ortho*RtoD)[:5])
	# 			# ax_lines.append([l_AoBapp_ortho, l_AoBapp_ortho.get_label()])

	# 			# l_AoBlos, = ax.plot(self.x_axis_list,[bottle.AoBlos * RtoD] * len(self.x_axis_list), "orange", label="AoBlos: " + str(bottle.AoBlos*RtoD)[:5])
	# 			# ax_lines.append([l_AoBlos, l_AoBlos.get_label()])

	# 	elif self.xaxis_azimut and bottle.observation_type == "fixed_elevation_discrete_rotation":
	# 		# print("DEBUG plot discrete")
	# 		if self.show_AoRD:
	# 			l_AoRD, = ax.plot(bottle.discrete_rotation_times, bottle.AoRD * RtoD, "k*", markersize=self.single_star_size*1.5, label="AoRD", zorder=2)
	# 			l_AoRD, = ax.plot(bottle.discrete_rotation_times, bottle.AoRD * RtoD, "*", color = self.AoRD_color, markersize=self.single_star_size, label="AoRD", zorder=3)
	# 			ax_lines.append([l_AoRD, l_AoRD.get_label()])
	# 		if self.show_AoBapp:
	# 			l_AoBapp, = ax.plot(bottle.discrete_rotation_times, bottle.AoBapp * RtoD, "k*", markersize=self.single_star_size*1.5, label="AoBapp", zorder=2)
	# 			l_AoBapp, = ax.plot(bottle.discrete_rotation_times, bottle.AoBapp * RtoD, "*", color = self.AoBapp_color, markersize=self.single_star_size, label="AoBapp", zorder=3)
	# 			ax_lines.append([l_AoBapp, l_AoBapp.get_label()])

	# 		rot = 0
	# 		if len(self.x_axis_ticks_pos) > 15:
	# 			rot = 60

	# 		ax.set_xticks(self.x_axis_ticks_pos)
	# 		ax.set_xticklabels(self.x_axis_ticks_label, rotation = rot)
	# 		# ax.xticks(rotation = 30)

	# 	elif self.xaxis_azimut and bottle.observation_type == "fixed_azimut_discrete_rotation":
	# 		# print("DEBUG plot discrete")
	# 		if self.show_AoRD:
	# 			l_AoRD, = ax.plot(bottle.discrete_rotation_times, bottle.AoRD * RtoD, "*", color = self.AoRD_color, markersize=self.marker_size*4, label="AoRD", zorder=2)
	# 			ax_lines.append([l_AoRD, l_AoRD.get_label()])
	# 		if self.show_AoBapp:
	# 			l_AoBapp, = ax.plot(bottle.discrete_rotation_times, bottle.AoBapp * RtoD, "*", color = self.AoBapp_color, markersize=self.marker_size*4, label="AoBapp", zorder=2)
	# 			ax_lines.append([l_AoBapp, l_AoBapp.get_label()])

	# 	elif self.xaxis_azimut and bottle.observation_type == "fixed_elevation_continue_rotation":

	# 		ax1_lim = self.ax1.get_ylim()
	# 		self.ax1.set_ylim(ax1_lim)
	# 		ax2_lim = self.ax2.get_ylim()
	# 		ax3_lim = ax.get_ylim()
	# 		for i in range(bottle.nb_continue_rotation):
	# 			# print("DEBUG ROT TIMES", i * 360 + (bottle.source_azimut * RtoD)%360)

	# 			if not self.show_RS_model:
	# 				self.ax1.plot([i * 360 + (bottle.source_azimut * RtoD)%360, i * 360 + (bottle.source_azimut * RtoD)%360], ax1_lim, "--k")
	# 				self.ax2.plot([i * 360 + (bottle.source_azimut * RtoD)%360, i * 360 + (bottle.source_azimut * RtoD)%360], ax2_lim, "--k")
	# 				ax.plot([i * 360 + (bottle.source_azimut * RtoD)%360, i * 360 + (bottle.source_azimut * RtoD)%360], ax3_lim, "--k")
	# 				if i != 0: # Draw red lines to delimit the rotations
	# 					# print(i * 360)
	# 					self.ax1.plot([i * 360, i * 360], ax1_lim, "r")
	# 					self.ax2.plot([i * 360, i * 360], ax2_lim, "r")
	# 					ax.plot([i * 360, i * 360], ax3_lim, "r")

	# 			if self.show_AoRD:
	# 				l_AoRD, = ax.plot(np.linspace(i * 360, (i+1) * 360, len(bottle.AoRD)), bottle.AoRD * RtoD, "*", color = self.AoRD_color, markersize=self.marker_size*2, label="AoRD", zorder=2)
	# 				ax_lines.append([l_AoRD, l_AoRD.get_label()])
	# 			if self.show_AoBapp:
	# 				l_AoBapp, = ax.plot(np.linspace(i * 360, (i+1) * 360, len(bottle.AoBapp)), bottle.AoBapp * RtoD, "*", color = self.AoBapp_color, markersize=self.marker_size*2, label="AoBapp", zorder=2)
	# 				ax_lines.append([l_AoBapp, l_AoBapp.get_label()])


	# 		ax.set_xticks(self.x_axis_ticks_pos)
	# 		ax.set_xticklabels(self.x_axis_ticks_label)

	# 		# x_axis, RD_diff, Bapps_diff = self.CompareAngles()
	# 		# if ax3_lim[1] > 100:
	# 		# 	print("DEBUG!!!!!!")
	# 		# 	RD_diff, Bapps_diff  = UnifyAngles(RD_diff, 1)[0], UnifyAngles(Bapps_diff, 1)[0]
	# 		# ax.plot(x_axis, RD_diff * RtoD, "*", color = "xkcd:black", markersize=self.marker_size, label="AoRD")

	# 		# ax_lines.append([l_AoRD, l_AoRD.get_label()])


	# def PlotRSModel(self, bottle):

	# 	# self.ax11 = self.ax1.twinx()
	# 	# self.ax22 = self.ax2.twinx()

	# 	tmp_xshift = 0
	# 	maxDoLP = np.max(bottle.data['smooth_DoLP'])

	# 	if self.show_model_list:
	# 		self.pomerol_models = [self.pomerol_models[i] for i in self.show_model_list]

	# 	for im, models in enumerate(self.pomerol_models):
	# 		# print(self.pomerol_model.x_axis)
	# 		# self.ax11.yaxis.set_visible(False)
	# 		self.ax1.plot(models.x_axis + tmp_xshift, models["I0"], "*", color = self.model_colors[im], marker = self.model_symbols[im], markersize=self.marker_size)
	# 		self.ax2.plot(models.x_axis + tmp_xshift, models["DoLP"], "*", color = self.model_colors[im], marker = self.model_symbols[im], markersize=self.marker_size)
	# 		# self.ax2.set_ylim((0, max(np.max(bottle.data['smooth_DoLP']), np.max(models.data["DoLP"]))))
	# 		if np.max(models.data["DoLP"]) > maxDoLP:
	# 			maxDoLP = np.max(models.data["DoLP"])
	# 			self.ax2.set_ylim((0, 1.1 * maxDoLP))

	# 		self.ax3.plot(models.x_axis + tmp_xshift, RtoD * SetAngleBounds(DtoR * models.data["AoRD"], -np.pi/2 + np.pi/2 * bottle.graph_angle_shift, np.pi/2 + np.pi/2 * bottle.graph_angle_shift), "*", color = self.model_colors[im], marker = self.model_symbols[im], markersize=self.marker_size)

	# 	if self.adapt_flux_scale:
	# 		border = 1.1
	# 		max_data, min_data = np.max(bottle.data['smooth_I0']), np.min(bottle.data['smooth_I0'])
	# 		ratio_data = max_data / min_data
	# 		print(f"Min data, Max Data, Ratio data: {min_data}, {max_data}, {ratio_data}")
	# 		nt, nb = [max_data], [min_data]
	# 		for models in self.pomerol_models:
	# 			max_model, min_model = np.max(models.data["I0"]) , np.min(models.data["I0"])
	# 			ratio_model = max_model / min_model
	# 			shift = (ratio_data * min_model - max_model) / (1 - ratio_data)
	# 			# shift = (ratio_model * min_data - max_data) / (1 - ratio_model)
	# 			print(f"Min model, Max model, Ratio model, Shift to match data (data unit): {min_model}, {max_model}, {ratio_model}, {shift}")
	# 			print(f"Min model, Max model, Ratio model, Shift to match data (model unit): {min_model / models.bottle_factor}, {max_model / models.bottle_factor}, {ratio_model}, {shift / models.bottle_factor}")
	# 			nt.append(max_model)
	# 			nb.append(min_model)

	# 		self.ax1.set_ylim(np.min(nb) / border, np.max(nt) * border)

	# 			# nt, nb = [np.max(bottle.data['smooth_DoLP'])], [np.min(bottle.data['smooth_DoLP'])]
	# 			# for models in self.pomerol_models:
	# 			# 	nt.append(np.max(models.data["DoLP"]))
	# 			# 	nb.append(np.min(models.data["DoLP"]))
	# 			# self.ax2.set_ylim(np.min(nb) / border, np.max(nt) * border)

	# 		if bottle.observation_type == "fixed_elevation_continue_rotation":
	# 			ax1_lim = self.ax1.get_ylim()
	# 			self.ax1.set_ylim(ax1_lim)
	# 			ax2_lim = self.ax2.get_ylim()
	# 			ax3_lim = self.ax3.get_ylim()
	# 			for i in range(bottle.nb_continue_rotation):
	# 				# print("DEBUG ROT TIMES", i * 360 + (bottle.source_azimut * RtoD)%360)

	# 				self.ax1.plot([i * 360 + (bottle.source_azimut * RtoD)%360, i * 360 + (bottle.source_azimut * RtoD)%360], ax1_lim, "--k")
	# 				self.ax2.plot([i * 360 + (bottle.source_azimut * RtoD)%360, i * 360 + (bottle.source_azimut * RtoD)%360], ax2_lim, "--k")
	# 				self.ax3.plot([i * 360 + (bottle.source_azimut * RtoD)%360, i * 360 + (bottle.source_azimut * RtoD)%360], ax3_lim, "--k")
	# 				if i != 0: # Draw red lines to delimit the rotations
	# 					# print(i * 360)
	# 					self.ax1.plot([i * 360, i * 360], ax1_lim, "r")
	# 					self.ax2.plot([i * 360, i * 360], ax2_lim, "r")
	# 					self.ax3.plot([i * 360, i * 360], ax3_lim, "r")

	# 			# if r1 >= r11:
	# 			# 	# self.ax1.set_ylim(0, top*border)
	# 			# 	# self.ax11.set_ylim(0, nb * r1*border)
	# 			# 	self.ax1.set_ylim(bot/border, top*border)
	# 			# 	self.ax11.set_ylim(nt / r1 /border, nt*border)
	# 			# 	# self.ax11.set_ylim(nb/border, nb * r1*border)
	# 			# 	print(f"NEW RATIO 2 {r1}, from {nb} to {nb * r1}")
	# 			# else:
	# 			# 	# self.ax1.set_ylim(0, top * 1.1)
	# 			# 	# self.ax11.set_ylim(0, nt * 1.1)
	# 			# 	# self.ax1.set_ylim(0, bot * r11*border)
	# 			# 	# self.ax11.set_ylim(0, nt*border)
	# 			# 	self.ax1.set_ylim(nt/r11/border, nt *border)
	# 			# 	# self.ax1.set_ylim(bot/border, bot * r11*border)
	# 			# 	self.ax11.set_ylim(nb/border, nt*border)
	# 			# 	print(f"NEW RATIO 1 {r11}, from {bot} to {bot * r11}")

	# 			# self.axs[11].xlim((0, ))


	# def MakePlots(self, bottle, comp=False):
	# 	#Plotting the mean intensity I0, DoLP, AoLP for each rotation and more!

	# 	# if self.show_error_bars:
	# 	# 	smooth_yerr = bottle.var_smooth_I0
	# 	# 	yerr = bottle.var_I0
	# 	# 	yerr = None
	# 	# else:
	# 	# 	smooth_yerr = None

	# 	print("START PLOTTING")

	# 	###Flux
	# 	if not comp:
	# 		self.PlotFlux(self.ax1, bottle, ax_lines=self.ax1_lines)
	# 		# self.PlotFluxDerivative(self.ax2.twinx(), bottle, ax_lines=self.ax1_lines)
	# 	else:
	# 		self.PlotFlux(self.ax1.twinx(), bottle, ax_lines=self.ax1_lines)

	# 	self.PlotAllSky(self.ax1, bottle, ax_lines=self.ax1_lines)
	# 	###DoLP
	# 	# if not comp:
	# 	self.PlotDoLP(self.ax2, bottle, ax_lines=self.ax2_lines)
	# 	# else:
	# 	# 	self.PlotDoLP(self.ax2.twinx(), bottle, ax_lines=self.ax2_lines)

	# 	if self.show_mag_data:
	# 		self.PlotB(self.ax2.twinx(), bottle, component = self.B_component, ax_lines=self.ax2_lines)
	# 		# self.PlotB(self.ax2.twinx(), bottle, component = self.B_component, ax_lines=self.ax2_lines, derivative = True)
	# 		# self.PlotB(self.ax2.twinx(), bottle, component = 'Horiz', ax_lines=self.ax2_lines)
	# 		# self.PlotB(self.ax3.twinx(), bottle, component = 'Vert', ax_lines=self.ax3_lines)

	# 	if self.show_eiscat:
	# 		print('plotting eiscat')

	# 		if self.eiscat_type == "uhf_v":
	# 			self.PlotEiscat(bottle, 'AoVi', self.ax3.twinx(), self.ax3_lines)
	# 		else:
	# 			self.PlotEiscat(bottle, 'ne', self.ax1.twinx(), self.ax1_lines)
	# 			# self.PlotEiscat(bottle, 'co', self.ax2.twinx(), self.ax2_lines)
	# 			self.PlotEiscat(bottle, 'vo', self.ax3.twinx(), self.ax3_lines)

	# 		# self.PlotEiscat(bottle, 'AoVi', self.ax3, self.ax3_lines)
	# 		# self.PlotEiscat(bottle, 'vo', self.ax3.twinx(), self.ax3_lines)
	# 		print('DONE plotting eiscat')

	# 	self.PlotAoLP(self.ax3, bottle, ax_lines=self.ax3_lines)

	# 	self.PlotCurrents(bottle)

	# 	self.PlotAoRD(self.ax3, bottle, ax_lines=self.ax3_lines)

	# 	if self.show_RS_model:
	# 		self.PlotRSModel(bottle)

	# 	# self.ax4.plot(self.x_axis_list, bottle.all_TempPM, "k.", linestyle = 'none', markersize=self.marker_size, label="PM")
	# 	# self.ax4.plot(self.x_axis_list, bottle.all_TempOptical, "r.", linestyle = 'none', markersize=self.marker_size, label="Optical")
	# 	# self.ax4.plot(self.x_axis_list, bottle.all_TempAmbiant, "b.", linestyle = 'none', markersize=self.marker_size, label="Ambiant")
	# 	# self.ax2.plot(self.x_axis_list,[DoLP_average] * nb_rot, "b", label="Avg: " + str(DoLP_average))


	# 	self.f1.subplots_adjust(hspace=0)


	# 	###Set title
	# 	# self.f1.suptitle(bottle.saving_name.replace("_", " "))

	# 	# plt.setp([a.get_xticklabels() for a in f1.axes[:-1]], visible=False)
	# 	if bottle.observation_type == "fixed":
	# 		plt.minorticks_on()


	# 	if self.use_24h_time_format and bottle.observation_type == "fixed":
	# 		self.ax3.xaxis.set_major_formatter(mdates.DateFormatter(self.xaxis_time_format))

	# 	# if self.show_time:
	# 	#
	# 	# 	self.ax12 = self.ax1.twiny()
	# 	# 	self.ax12.set_xlim(self.ax1.get_xlim())
	# 	#
	# 	# 	# print("DEBUG*********************************************************************************************")
	# 	# 	# print("DEBUG*********************************************************************************************")
	# 	#
	# 	# 	if self.xaxis_azimut and bottle.observation_type in ["fixed_elevation_continue_rotation", "fixed_elevation_discrete_rotation"]:
	# 	# 		xticks = self.x_time_ticks_pos
	# 	# 		plt.xticks(ticks=xticks, labels=self.x_time_ticks_label)
	# 	#
	# 	# 	else:
	# 	# 		print(plt.xticks()[0], self.divisor)
	# 	# 		xticks = plt.xticks()[0] * self.divisor
	# 	# 		xticks = [bottle.DateTime("start", format=self.time_format, delta=dt.timedelta(seconds = t)) for t in xticks]
	# 	# 		# xticks = [bottle.DateTime("start", format=self.time_format) + dt.timedelta(seconds = t) for t in xticks]
	# 	# 		# xticks = [bottle.DateTime() + dt.timedelta(seconds = t) + bottle.head_jump for t in xticks]
	# 	#
	# 	# 	print(xticks)
	# 	#
	# 	# 	if bottle.observation_type == "fixed":
	# 	# 		self.ax12.set_xticklabels([st.strftime("%H:%M") for st in xticks])
	# 	# 	self.ax12.set_xlabel(self.time_label)

	# 	if self.show_legend:
	# 		self.ax1.legend(list(zip(*self.ax1_lines))[0], list(zip(*self.ax1_lines))[1])
	# 		self.ax2.legend(list(zip(*self.ax2_lines))[0], list(zip(*self.ax2_lines))[1])
	# 		self.ax3.legend(list(zip(*self.ax3_lines))[0], list(zip(*self.ax3_lines))[1], loc = "lower center")


	# 	print("Saving graphs in", bottle.data_file_name + "/" + bottle.saving_name + '_graphs.png')
	# 	if bottle.instrument_name in ["carmen", "corbel", "gdcu"]:
	# 		plt.savefig(bottle.data_file_name + "/" + bottle.saving_name + '_graphs.png', bbox_inches='tight')
	# 		# plt.savefig(bottle.data_file_name + "/" + bottle.saving_name + '_graphs.eps', bbox_inches='tight')
	# 	else:
	# 		plt.savefig("/".join(bottle.data_file_name.split("/")[:-1]) + "/" + bottle.saving_name + '_graphs.png', bbox_inches='tight')
	# 		# plt.savefig("/".join(bottle.data_file_name.split("/")[:-1]) + "/" + bottle.saving_name + '_graphs.eps', bbox_inches='tight')


	# def CompareRayleigh(self, bottle, model):

	# 	data_I = np.interp(model.x_axis, self.x_axis_list, bottle.data['smooth_I0'])
	# 	data_D = np.interp(model.x_axis, self.x_axis_list, bottle.data['smooth_DoLP'])
	# 	data_A = np.interp(model.x_axis, self.x_axis_list, bottle.data['smooth_AoLP'])

	# 	ratio = lambda x, y: x / y
	# 	diff = lambda x, y: x - y

	# 	fig, axs = plt.subplots(3, figsize=(50, 20))

	# 	axs[0].plot(model.x_axis, ratio(data_I, model.flux), "*")
	# 	axs[1].plot(model.x_axis, ratio(data_D, model.DoLP), "*")
	# 	axs[2].plot(model.x_axis, GetAnglesDifference(data_A, model.AoLP)*RtoD, "*")

	# 	axs[0].set_ylabel("Ratio")
	# 	axs[1].set_ylabel("Ratio")
	# 	axs[2].set_ylabel("Diff")
	# 	fig.suptitle("Data /- Model")

	# 	if bottle.instrument_name in ["carmen", "corbel", "gdcu"]:
	# 		plt.savefig(bottle.data_file_name + "/" + bottle.saving_name + '_RSmodel_comp.png', bbox_inches='tight')
	# 	else:
	# 		plt.savefig("/".join(bottle.data_file_name.split("/")[:-1]) + "/" + bottle.saving_name + '_RSmodel_comp.png', bbox_inches='tight')


	# def MakeSNratio(self, bottle):
	# 	if not self.show_SN:
	# 		return 0

	# 	f5, ax = plt.subplots(1, figsize=(10, 20))

	# 	ax.plot(self.x_axis_list, bottle.all_SN, ".", color = self.all_SN_color, label="all SN")
	# 	ax.set_ylabel("Raw Data")
	# 	ax1 = plt.twinx(ax)
	# 	ax1.plot(self.x_axis_list, bottle.smooth_SN, ".", color = self.smooth_SN_color, label="smooth SN")
	# 	ax1.set_ylabel("Smoothed Data")

	# 	ax.set_xlabel(self.xlabel)

	# 	plt.minorticks_on()

	# 	ax.set_title("SN ratio")
	# 	print("Saving graphs in", bottle.data_file_name + "/" + bottle.saving_name + '_SNratio.png')
	# 	if bottle.instrument_name in ["carmen", "corbel", "gdcu"]:
	# 		plt.savefig(bottle.data_file_name + "/" + bottle.saving_name + '_SNratio.png', bbox_inches='tight')
	# 		# plt.savefig(bottle.data_file_name + "/" + bottle.saving_name + '_graphs.eps', bbox_inches='tight')
	# 	else:
	# 		plt.savefig("/".join(bottle.data_file_name.split("/")[:-1]) + "/" + bottle.saving_name + '_SNratio.png', bbox_inches='tight')
	# 		# plt.savefig("/".join(bottle.data_file_name.split("/")[:-1]) + "/" + bottle.saving_name + '_graphs.eps', bbox_inches='tight')


	def TestResample(self, bottle):
		f, axs = plt.subplots(3, sharex=True)
		
		N_pts = int(bottle.avg_dt)
		N_rot = len(bottle.data['I0'])
		filter_angles = np.linspace(0, N_rot * 2*np.pi, N_rot * N_pts) #List of angles (rad) for the polarising filter between 0 and 2π.

		resamp_signal, resampled_t = sci.signal.resample(bottle.data['I0'], N_rot * N_pts, t = self.x_axis_list)
		resamp_signal /= 2.


		def GetPola(V, Vcos, Vsin):
			"""Given V, Vcos, Vsin, returns the initial intensity, DoLP and AoLP. This method is shared for spp and ptcu. It is also a static method, that can be called outside of the object. This way it can be used everywhere, each time you need I0, DoLP, AoLP to decrease the chances of making a mistake."""
			I0 = 2 * V
			DoLP = 2 * np.sqrt(Vcos**2 + Vsin**2) / V * 100
			AoLP = np.arctan2(Vsin, Vcos) / 2
			return abs(I0), abs(DoLP), AoLP

		def GetV(fake_signal):
			"""return the average of a signal over 1 rotation"""
			integral = sum(fake_signal)
			return integral / N_pts

		def GetVcos(fake_signal):
			"""Return the average of V*cos(2*theta) over 1 rotation"""
			x = fake_signal * np.cos(2 * filter_angles)
			return sum(x) / N_pts

		def GetVsin(fake_signal):
			"""Return the average value of -V*sin(2*theta) over 1 rotation"""
			y = - fake_signal * np.sin(2 * filter_angles)
			return sum(y) / N_pts


		def GetStokesTime(fake_signal):
			"""return the stokes parameters of a signal over N_rot rotation"""
			V = np.zeros(N_rot)
			Vcos = np.zeros(N_rot)
			Vsin = np.zeros(N_rot)
			for ir in range(N_rot):
				start_rot = ir * N_pts
				end_rot = start_rot + N_pts
				tmp_signal = fake_signal[start_rot:end_rot]

				V[ir] = sum(tmp_signal) / N_pts
				Vcos[ir] = sum(tmp_signal * np.cos(2 * filter_angles[start_rot:end_rot])) / N_pts
				Vsin[ir] = sum(tmp_signal * np.sin(2 * filter_angles[start_rot:end_rot])) / N_pts

			return V, Vcos, Vsin


		V, Vcos, Vsin = GetStokesTime(resamp_signal)
		I_list, DoLP_list, AoLP_list = GetPola(V, Vcos, Vsin)

		dV, dVcos, dVsin = bottle.all_V, bottle.all_Vcos - Vcos, bottle.all_Vsin - Vsin
		dI, dD, dA =  GetPola(dV, dVcos, dVsin)

		axs[0].plot(self.x_axis_list, bottle.data['I0'])
		axs[0].plot(self.x_axis_list, I_list)
		# axs[0].plot(self.x_axis_list, dI)
		axs[1].plot(self.x_axis_list, bottle.data['DoLP'])
		axs[1].plot(self.x_axis_list, DoLP_list)
		# axs[1].plot(self.x_axis_list, dD)
		axs[2].plot(self.x_axis_list, bottle.data['AoLP']*RtoD)
		axs[2].plot(self.x_axis_list, AoLP_list*RtoD)
		# axs[2].plot(self.x_axis_list, dA*RtoD)



	# def MakeI0SNratio(self, bottle):
	# 	f6, ax = plt.subplots(1, figsize=(10, 20))

	# 	# print(bottle.data['I0'],  bottle.data['std_I0'], bottle.data['I0'] / bottle.data['std_I0'])
	# 	# print(bottle.data['smooth_I0'], bottle.data['std_smooth_I0'], bottle.data['smooth_I0'] / bottle.data['std_smooth_I0'])

	# 	ax.plot(self.x_axis_list, bottle.data['I0'] / bottle.data['std_I0'], "--", color = self.all_SN_color, label="I0 SN")
	# 	ax.set_ylabel("Raw")
	# 	ax1 = plt.twinx(ax)
	# 	ax1.plot(self.x_axis_list, bottle.data['smooth_I0'] / bottle.data['std_smooth_I0'], "--", color = self.smooth_SN_color, label="smooth I0 SN")
	# 	ax1.set_ylabel("Smooth")

	# 	ax.set_title("Intensity SN ratio")
	# 	print("Saving graphs in", bottle.data_file_name + "/" + bottle.saving_name + '_I0SNratio.png')
	# 	if bottle.instrument_name in ["carmen", "corbel", "gdcu"]:
	# 		plt.savefig(bottle.data_file_name + "/" + bottle.saving_name + '_I0SNratio.png', bbox_inches='tight')
	# 		# plt.savefig(bottle.data_file_name + "/" + bottle.saving_name + '_graphs.eps', bbox_inches='tight')
	# 	else:
	# 		plt.savefig("/".join(bottle.data_file_name.split("/")[:-1]) + "/" + bottle.saving_name + '_I0SNratio.png', bbox_inches='tight')
	# 		# plt.savefig("/".join(bottle.data_file_name.split("/")[:-1]) + "/" + bottle.saving_name + '_graphs.eps', bbox_inches='tight')


	# def DrawEiscatParam(self, axe, param, alt):
	# 	label = param.replace("_", "\_")

	# 	t, d, e = self.eiscat_data.GetParameter(param, alt, time_divisor = self.divisor)
	# 	l_eiscat_data, capline, barlinecol = axe.errorbar(t, d, yerr = e, fmt = ".", label = label)
	# 	avg = np.average(d)
	# 	axe.plot(t, [avg] * len(t), label = "Average: " + str(avg)[:4])

	# def MakeEiscatDataPlot(self, bottle):
	# 	f4, ax = plt.subplots(4, figsize=(10, 20), sharex = True)

	# 	altitude = bottle.filters[0]
	# 	if altitude == "r": altitude = 220
	# 	elif altitude == "v": altitude = 110
	# 	elif altitude == "b" or "m": altitude = 85

	# 	self.DrawEiscatParam(ax[0], "ne", altitude)
	# 	self.DrawEiscatParam(ax[1], "ti", altitude)
	# 	self.DrawEiscatParam(ax[2], "tr", altitude)
	# 	# self.DrawEiscatParam(ax[3], "nu_in", altitude)
	# 	self.DrawEiscatParam(ax[3], "vo", altitude)
	# 	ax[-1].set_xlabel(self.xlabel)

	# 	for a in ax:
	# 		a.legend()

	# 	print("Saving graphs in", bottle.data_file_name + "/" + bottle.saving_name + '_eiscat.png')
	# 	if bottle.instrument_name in ["carmen", "corbel", "gdcu"]:
	# 		plt.savefig(bottle.data_file_name + "/" + bottle.saving_name + '_eiscat.png', bbox_inches='tight')
	# 	else:
	# 		plt.savefig("/".join(bottle.data_file_name.split("/")[:-1]) + "/" + bottle.saving_name + '_eiscat.png', bbox_inches='tight')


	# def MakeMagDataPlots(self, bottle):
	# 	f3, ax = plt.subplots(3, figsize=(10, 20), sharex = True)

	# 	ax[0].plot(self.mag_data.times_sec, self.mag_data.Total, "k", label="Total")
	# 	# ax[1].plot(self.mag_data.times_sec, self.mag_data.Dec, "r", label="Dec")
	# 	ax[1].plot(self.mag_data.times_sec, self.mag_data.Horiz, "b", label="Horiz")
	# 	# ax[3].plot(self.mag_data.times_sec, self.mag_data.Vert, "g", label="Vert")
	# 	ax[2].plot(self.mag_data.times_sec, self.mag_data.Incl, "y", label="Incl")
	# 	ax[-1].set_xlabel(self.xlabel)
	# 	for a in ax:
	# 		a.legend()
	# 	print("Saving graphs in", bottle.data_file_name + "/" + bottle.saving_name + '_magneto.png')
	# 	if bottle.instrument_name in ["carmen", "corbel", "gdcu"]:
	# 		plt.savefig(bottle.data_file_name + "/" + bottle.saving_name + '_magneto.png', bbox_inches='tight')
	# 	else:
	# 		plt.savefig("/".join(bottle.data_file_name.split("/")[:-1]) + "/" + bottle.saving_name + '_magneto.png', bbox_inches='tight')



	# def CompareAngles(self, bottle):
	# 	if bottle.observation_type == "fixed_elevation_continue_rotation":
	# 		AoRDs = bottle.AoRD
	# 		AoBapps = bottle.AoBapp
	# 		for i in range(bottle.nb_continue_rotation - 1):
	# 			AoRDs = np.append(AoRDs, bottle.AoRD)
	# 			AoBapps = np.append(AoBapps, bottle.AoBapp)

	# 		nb_AoLP = len(bottle.data['smooth_AoLP'])
	# 		nb_AoRD = len(AoRDs)
	# 		nb_AoBapp = len(AoBapps)
	# 		if nb_AoLP > nb_AoRD:
	# 			x_axis = ReduceList(self.x_axis_list, nb_AoRD)
	# 			AoLPs = ReduceList(bottle.data['smooth_AoLP'], nb_AoRD)
	# 		elif nb_AoLP < nb_AoRD:
	# 			x_axis = self.x_axis_list
	# 			AoLPs = bottle.data['smooth_AoLP']
	# 			AoRDs = ReduceList(AoRDs, nb_AoLP)
	# 			AoBapps = ReduceList(AoBapps, nb_AoLP)
	# 		else:
	# 			x_axis = self.x_axis_list
	# 			AoLPs = bottle.data['smooth_AoLP']

	# 		RD_diff = GetAnglesDifference(AoLPs, AoRDs)
	# 		Bapps_diff = GetAnglesDifference(AoLPs, AoBapps)

	# 		return x_axis, RD_diff, Bapps_diff


	# def CompareAnglesPlots(self, bottle):
	# 	f4, ax = plt.subplots(1, figsize=(10, 20), sharex = True)

	# 	x_axis, RD_diff, Bapps_diff = self.CompareAngles(bottle)

	# 	ax.plot(x_axis, RD_diff * RtoD, "*", color = self.AoRD_color, label="Total")
	# 	ax.plot(x_axis, Bapps_diff * RtoD, "*", color = self.AoBapp_color, label="Total")

	# 	print("Saving graphs in", bottle.data_file_name + "/" + bottle.saving_name + '_angles_comp.png')
	# 	if bottle.instrument_name in ["carmen", "corbel", "gdcu"]:
	# 		plt.savefig(bottle.data_file_name + "/" + bottle.saving_name + '_angles_comp.png', bbox_inches='tight')
	# 	else:
	# 		plt.savefig("/".join(bottle.data_file_name.split("/")[:-1]) + "/" + bottle.saving_name + '_angles_comp.png', bbox_inches='tight')


	# def MakeSmartCorrelationPlots(self, bottle, diff_thresholds=None, smooth=True, COMP = None):

	# 	nb_subplots = 3

	# 	gridspec = [nb_subplots - 1].extend([1]*(nb_subplots-1))

	# 	f2, axs = plt.subplots(ncols=1, nrows=nb_subplots, figsize=(20, 20), gridspec_kw = {'hspace':0, 'height_ratios':gridspec})

	# 	axs[2].get_shared_x_axes().join(axs[1], axs[2])
	# 	if nb_subplots > 3:
	# 		axs[3].get_shared_x_axes().join(axs[1], axs[3])

	# 	colors = ["b", 'k', "r", 'y', 'g']
	# 	masks = []

	# 	if not COMP:
	# 		COMP = "dolp"


	# 	X_array = bottle.I0_diff
	# 	Y_array = bottle.data['smooth_DoLP']
	# 	if COMP.lower() == "aolp":
	# 		Y_array = bottle.data['smooth_AoLP'] * RtoD
	# 	if not smooth:
	# 		# mask_array = bottle.GetWaveletTransform(bottle.all_I0_diff, 12, unit = "seconds") #bottle.all_I0_diff
	# 		X_array = bottle.all_I0_diff
	# 		Y_array = bottle.data['DoLP']
	# 		if COMP.lower() == "aolp":
	# 			Y_array = bottle.data['AoLP'] * RtoD

	# 	mask_array = bottle.GetSliddingCorrCoef(window_size = 2, smooth = smooth) #bottle.all_I0_diff
	# 	# mask_array = bottle.GetSliddingCorrCoef(window_size = 2, smooth = smooth) #bottle.all_I0_diff
	# 	if not diff_thresholds:
	# 		# diff_thresholds = np.percentile(mask_array, [25, 75])
	# 		diff_thresholds = [-0.5, 0.5]


	# 	diff_thresholds = np.sort(diff_thresholds)
	# 	# print(diff_thresholds)
	# 	masks.append(mask_array < diff_thresholds[0])
	# 	for id in range(1, len(diff_thresholds)):
	# 		masks.append((diff_thresholds[id-1] <= mask_array) * (mask_array < diff_thresholds[id]))
	# 	masks.append(mask_array >= diff_thresholds[-1])

	# 	subplots_counter = 0

	# 	for im, m in enumerate(masks):
	# 		# if im != 1:
	# 		X = np.ma.masked_equal(X_array * m, 0)
	# 		Y = np.ma.masked_equal(Y_array * m, 0)
	# 		axs[subplots_counter].plot(X, Y, ".", color = colors[im], alpha = 1)

	# 	# axs[subplots_counter].plot(mask_array, bottle.data['DoLP'], ".")

	# 	axs[subplots_counter].xaxis.tick_top()
	# 	axs[subplots_counter].set_xlabel("Dérivée Intensity")
	# 	axs[subplots_counter].xaxis.set_label_position('top')
	# 	if COMP.lower() == "aolp":
	# 		axs[subplots_counter].set_ylabel("AoLP (deg)")
	# 	else:
	# 		axs[subplots_counter].set_ylabel("DoLP (%)")



	# 	#Slidding correlation subplot
	# 	# subplots_counter += 1
	# 	#
	# 	# axs[subplots_counter].plot(self.x_axis_list, bottle.GetSliddingCorrCoef(window_size = 2, smooth = smooth), "-")
	# 	# axs[subplots_counter].plot(self.x_axis_list, mask_array, ".")
	# 	# #
	# 	# axs[subplots_counter].fill_between([self.x_axis_list[0], self.x_axis_list[-1]], [diff_thresholds[0], diff_thresholds[0]], np.min(mask_array), color = colors[0], zorder=-100, alpha = 0.3)
	# 	# for it in range(1, len(diff_thresholds)):
	# 	# 	axs[subplots_counter].fill_between([self.x_axis_list[0], self.x_axis_list[-1]], [diff_thresholds[it], diff_thresholds[it]], diff_thresholds[it-1], color = colors[it], zorder=-100, alpha = 0.3)
	# 	#
	# 	# axs[subplots_counter].fill_between([self.x_axis_list[0], self.x_axis_list[-1]], [diff_thresholds[-1], diff_thresholds[-1]], np.max(mask_array), color = colors[len(diff_thresholds)], zorder=-100, alpha = 0.3)
	# 	# # axs[subplots_counter].set_ylabel("dI")
	# 	# axs[subplots_counter].set_ylabel("Corr")


	# 	#Flux and DoLP subplot
	# 	subplots_counter += 1

	# 	ax33 = axs[subplots_counter].twinx()
	# 	axs[subplots_counter].plot(self.x_axis_list, Y_array, "r-")
	# 	if not smooth:
	# 		ax33.plot(self.x_axis_list, bottle.data['I0'], "b-")
	# 	else:
	# 		ax33.plot(self.x_axis_list, bottle.data['smooth_I0'], "b-")
	# 	if COMP.lower() == "aolp":
	# 		axs[subplots_counter].set_ylabel("I,  AoLP")
	# 	else:
	# 		axs[subplots_counter].set_ylabel("I,  DoLP")

	# 	#dFlux/dt and DoLP subplot
	# 	subplots_counter += 1

	# 	ax44 = axs[subplots_counter].twinx()
	# 	axs[subplots_counter].plot(self.x_axis_list, Y_array, "r-")
	# 	ax44.plot(self.x_axis_list, X_array, "b-")
	# 	if COMP.lower() == "aolp":
	# 		axs[subplots_counter].set_ylabel("dI,  AoLP")
	# 	else:
	# 		axs[subplots_counter].set_ylabel("dI,  DoLP")
	# 	axs[subplots_counter].set_xlabel(self.xlabel)


	# 	# f2.suptitle(f"smooth:{smooth}; {COMP.upper()}")

	# 	print("Saving correlation in", bottle.data_file_name + "/" + bottle.saving_name + '_smart_correlations.png')
	# 	if bottle.instrument_name in ["carmen", "corbel", "gdcu"]:
	# 		plt.savefig(bottle.data_file_name + "/" + bottle.saving_name + '_S' + str(smooth) + COMP.upper() + '_smart_correlation.png', bbox_inches='tight')
	# 	else:
	# 		plt.savefig("/".join(bottle.data_file_name.split("/")[:-1]) + '_S' + str(smooth) + "/" + bottle.saving_name + '_smart_correlation.png', bbox_inches='tight')

	# def MakeCleanCorrelationPlots(self, bottle, diff_thresholds=None, smooth=True):

	# 		nb_subplots = 2

	# 		gridspec = [1, 1]

	# 		f2, axs = plt.subplots(ncols=1, nrows=nb_subplots, figsize=(20, 20), gridspec_kw = {'hspace':0, 'height_ratios':gridspec})
	# 		subplots_counter = 0

	# 		for ax in axs:
	# 			ax.set_axisbelow(True)
	# 			ax.grid(True, zorder=-10)

	# 		colors = ["b", 'k', "r", 'y', 'g']

	# 		for icomp, COMP in enumerate(['dolp', 'aolp']):
	# 			masks = []

	# 			if smooth:
	# 				X_array = bottle.I0_diff
	# 				X_error = bottle.I0_diff_std
	# 				Y_array = bottle.data['smooth_DoLP']
	# 				Y_error = bottle.data['std_smooth_DoLP']
	# 				if COMP.lower() == "aolp":
	# 					Y_array = bottle.data['smooth_AoLP'] * RtoD
	# 					Y_error = bottle.data['std_smooth_AoLP'] * RtoD
	# 			elif not smooth:
	# 				# mask_array = bottle.GetWaveletTransform(bottle.all_I0_diff, 12, unit = "seconds") #bottle.all_I0_diff
	# 				X_array = bottle.all_I0_diff
	# 				X_error = bottle.all_I0_diff_std
	# 				Y_array = bottle.data['DoLP']
	# 				Y_error = bottle.data['std_DoLP']
	# 				if COMP.lower() == "aolp":
	# 					Y_array = bottle.data['AoLP'] * RtoD
	# 					Y_error = bottle.data['std_AoLP'] * RtoD

	# 			mask_array = bottle.GetSliddingCorrCoef(window_size = 2, smooth = smooth) #bottle.all_I0_diff
	# 			# mask_array = bottle.GetSliddingCorrCoef(window_size = 2, smooth = smooth) #bottle.all_I0_diff
	# 			if diff_thresholds is None:
	# 				# diff_thresholds = np.percentile(mask_array, [25, 75])
	# 				diff_thresholds = [-0.5, 0.5]


	# 			diff_thresholds = np.sort(diff_thresholds)
	# 			# print(diff_thresholds)
	# 			masks.append(mask_array < diff_thresholds[0])
	# 			for id in range(1, len(diff_thresholds)):
	# 				masks.append((diff_thresholds[id-1] <= mask_array) * (mask_array < diff_thresholds[id]))
	# 			masks.append(mask_array >= diff_thresholds[-1])

	# 			for im, m in enumerate(masks):
	# 				# if im != 1:
	# 				X = np.ma.masked_equal(X_array * m, 0)
	# 				Y = np.ma.masked_equal(Y_array * m, 0)
	# 				Xerr = np.ma.masked_equal(X_error * m, 0)
	# 				Yerr = np.ma.masked_equal(Y_error * m, 0)
	# 				axs[subplots_counter].errorbar(X, Y, yerr = Yerr, xerr = Xerr, fmt=".", color = colors[im], alpha = 1)

	# 			if COMP.lower() == "aolp":
	# 				axs[subplots_counter].set_ylabel("AoLP (deg)")
	# 			else:
	# 				axs[subplots_counter].set_ylabel("DoLP (%)")

	# 			subplots_counter += 1

	# 		axs[0].xaxis.tick_top()
	# 		axs[1].set_xlabel("Radiance time derivative (AU)")
	# 		axs[1].set_ylim(-90, 90)


	# 		print("Saving correlation in", bottle.data_file_name + "/" + bottle.saving_name + '_clean_correlations.png')
	# 		if bottle.instrument_name in ["carmen", "corbel", "gdcu"]:
	# 			plt.savefig(bottle.data_file_name + "/" + bottle.saving_name + '_S' + str(smooth) + '_clean_correlation.png', bbox_inches='tight')
	# 		else:
	# 			plt.savefig("/".join(bottle.data_file_name.split("/")[:-1]) + '_S' + str(smooth) + "/" + bottle.saving_name + '_clean_correlation.png', bbox_inches='tight')


	# def MakeCorrelationPlots(self, bottle):
	# 	f2, (ax1, ax2) = plt.subplots(2, figsize=(10, 20))

	# 	ax1.plot(bottle.data['I0'], bottle.data['DoLP'], "k+", label="Raw")
	# 	ax1.plot(bottle.data['smooth_I0'], bottle.data['smooth_DoLP'], "r+", label="Smoothed (" + str(bottle.smoothing_factor) + " " + str(bottle.smoothing_unit)[:3] + ")")
	# 	ax1.set_xlabel("Intensity (mV)")
	# 	ax1.set_ylabel("DoLP (%)")
	# 	ax1.legend()

	# 	if not bottle.NoVref and bottle.instrument_name=="carmen":
	# 		ax2.plot(bottle.data['I0'], bottle.all_Iref, "k+", label="Raw")
	# 		ax2.plot(bottle.data['smooth_I0'], bottle.smooth_Iref, "r+", label="Smoothed (" + str(bottle.smoothing_factor) + " " + str(bottle.smoothing_unit)[:3] + ")")
	# 		ax2.set_xlabel("Ref Intensity (mV)")
	# 		ax2.set_ylabel("Pola Intensity (mV)")

	# 		ax2.legend()

	# 	print("Saving correlation in", bottle.data_file_name + "/" + bottle.saving_name + '_correlations.png')
	# 	if bottle.instrument_name in ["carmen", "corbel", "gdcu"]:
	# 		plt.savefig(bottle.data_file_name + "/" + bottle.saving_name + '_correlation.png', bbox_inches='tight')
	# 	else:
	# 		plt.savefig("/".join(bottle.data_file_name.split("/")[:-1]) + "/" + bottle.saving_name + '_correlation.png', bbox_inches='tight')


	# def GetCoherence(self, data1, data2, freq):
	# 	NN = 2**9
	# 	f, P12 = signal.csd(data1, data2, freq, detrend='linear',nperseg=NN)
	# 	f, P11 = signal.csd(data1, data1, freq, detrend='linear',nperseg=NN)
	# 	f, P22 = signal.csd(data2, data2, freq, detrend='linear',nperseg=NN)
	# 	Nf = f.shape[0]

	# 	C = np.abs(P12)**2/(P11*P22)

	# 	ss = np.imag(P12[:])/np.abs(P12)
	# 	cc = np.real(P12[:])/np.abs(P12)
	# 	angle_12 = np.arctan2(ss,cc)

	# 	return f, C, angle_12

	# def MakeCoherencePLots(self, bottle):
	# 	f_coherence, axs = plt.subplots(2, figsize=(10, 20))

	# 	ax_coherence = axs[0]
	# 	ax_angle = axs[1]

	# 	f, C, A = self.GetCoherence(bottle.data['smooth_I0'], bottle.data['smooth_DoLP'], bottle.GetRotationFrequency())
	# 	ax_coherence.semilogx( (1./f), C, "k", label="Intensity-DoLP")
	# 	ax_angle.semilogx( (1./f), A * RtoD, "k", label="Intensity-DoLP")

	# 	f, C, A= self.GetCoherence(bottle.data['smooth_I0'], bottle.data['smooth_AoLP'], bottle.GetRotationFrequency())
	# 	ax_coherence.semilogx( (1./f), C, "r", label="Intensity-AoLP")
	# 	ax_angle.semilogx( (1./f), A * RtoD, "r", label="Intensity-AoLP")

	# 	f, C, A = self.GetCoherence(bottle.data['smooth_DoLP'], bottle.data['smooth_AoLP'], bottle.GetRotationFrequency())
	# 	ax_coherence.semilogx( (1./f), C, "b", label="DoLP-AoLP")
	# 	ax_angle.semilogx( (1./f), A * RtoD, "b", label="DoLP-AoLP")

	# 	ax_coherence.set_xlabel("Period (sec)")
	# 	ax_coherence.set_ylabel("Coherence")
	# 	# ax_coherence.set_xlim((-1,1))
	# 	ax_coherence.set_ylim((0,1))
	# 	ax_coherence.legend()

	# 	ax_angle.set_xlabel("Period (sec)")
	# 	ax_angle.set_ylabel("Shift angle")
	# 	# ax_angle.set_xlim((-1,1))
	# 	ax_angle.set_ylim((-180,180))
	# 	ax_angle.legend()


	# 	print("Saving coherence in", bottle.data_file_name + "/" + bottle.saving_name + '_coherence.png')
	# 	if bottle.instrument_name in ["carmen", "corbel", "gdcu"]:
	# 		plt.savefig(bottle.data_file_name + "/" + bottle.saving_name + '_coherence.png', bbox_inches='tight')
	# 	else:
	# 		plt.savefig("/".join(bottle.data_file_name.split("/")[:-1]) + "/" + bottle.saving_name + 'coherence.png', bbox_inches='tight')



	# def MakeFFT(self, bottle):
	# 	self.f1, (self.ax1, self.ax2, self.ax3) = plt.subplots(3, sharex=True, figsize=(16, 8))
	# 	self.ax1.plot(self.x_axis_list, bottle.data['I0'], "k.")
	# 	self.ax1.plot(self.x_axis_list, bottle.GetFourierTransform("I0"), "r.")
	# 	self.ax2.plot(self.x_axis_list, bottle.data['DoLP'], "k.")
	# 	self.ax2.plot(self.x_axis_list, bottle.GetFourierTransform("DoLP"), "r.")
	# 	self.ax3.plot(self.x_axis_list, bottle.data['AoLP'], "k.")
	# 	self.ax3.plot(self.x_axis_list, bottle.GetFourierTransform("AoLP"), "r.")
	# 	#
	# 	# self.ax1.set_ylabel("Intensity FFT")
	# 	# self.ax2.set_ylabel("DoLP FFT")
	# 	# self.ax3.set_ylabel("AoLP FFT")
	# 	#
	# 	# xaxis = np.fft.fftfreq(self.x_axis_list.shape[-1], d = 1000/20.)
	# 	#
	# 	# FFT_I0 = np.fft.fft(bottle.data['I0'])
	# 	# FFT_DoLP = np.fft.fft(bottle.data['DoLP'])
	# 	# FFT_AoLP = np.fft.fft(bottle.data['AoLP'])
	# 	#
	# 	# smooth_FFT_I0 = np.fft.fft(bottle.data['data['smooth_I0']'])
	# 	# smooth_FFT_DoLP = np.fft.fft(bottle.data['smooth_DoLP'])
	# 	# smooth_FFT_AoLP = np.fft.fft(bottle.data['smooth_AoLP'])
	# 	#
	# 	# self.ax1.plot(xaxis, FFT_I0, xaxis, smooth_FFT_I0)
	# 	# self.ax2.plot(xaxis, FFT_DoLP, xaxis, smooth_FFT_DoLP)
	# 	# self.ax3.plot(xaxis, FFT_AoLP, xaxis, smooth_FFT_AoLP)

	# def CompareBottles(self, b1, b2):
	# 	f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, figsize=(16, 8))

	# 	if b1.nb_rot > b2.nb_rot:
	# 		all_times = list(map(lambda x: x.total_seconds(), b2.data['Times']))
	# 		# print("DEBUG all_times", len(all_times), all_times[0])
	# 		I0_1, DoLP_1, AoLP_1 = b1.GetInterpolation(all_times)
	# 		I0_2, DoLP_2, AoLP_2 = b2.data['smooth_I0'], b2.data['smooth_DoLP'], b2.data['smooth_AoLP']
	# 	else:
	# 		all_times = list(map(lambda x: x.total_seconds(), b1.data['Times']))
	# 		# print("DEBUG all_times", len(all_times), all_times[0])
	# 		I0_2, DoLP_2, AoLP_2 = b2.GetInterpolation(all_times)
	# 		I0_1, DoLP_1, AoLP_1 = b1.data['smooth_I0'], b1.data['smooth_DoLP'], b1.data['smooth_AoLP']

	# 	delta = lambda a, b: (b-a) / a * 100
	# 	ax1.plot(self.x_axis_list, delta(I0_1, I0_2))
	# 	ax2.plot(self.x_axis_list, delta(DoLP_1, DoLP_2))
	# 	# ax3.plot(all_times, AoLP_1*RtoD)
	# 	# ax3.plot(all_times, AoLP_2*RtoD)
	# 	# ax3.plot(self.x_axis_list, delta(AoLP_1, AoLP_2))
	# 	ax3.plot(self.x_axis_list, (AoLP_1 - AoLP_2)*RtoD)

	# 	ax1.set_ylabel("I0")
	# 	ax1.grid(which="both")

	# 	ax2.set_ylabel("DoLP")
	# 	ax2.grid(which="both")

	# 	ax3.set_ylabel("AoLP")
	# 	ax3.set_xlabel(self.xlabel)
	# 	ax3.grid(which="both")

	# 	print("Saving bottle comparison in", b1.data_file_name + "/" + b1.saving_name + '_bottle_comparison.png')
	# 	if b1.instrument_name in ["carmen", "corbel", "gdcu"]:
	# 		plt.savefig(b1.data_file_name + "/" + b1.saving_name + '_bottle_comparison.png', bbox_inches='tight')
	# 	else:
	# 		plt.savefig("/".join(b1.data_file_name.split("/")[:-1]) + "/" + b1.saving_name + 'bottle_comparison.png', bbox_inches='tight')



	def GetDataFromTxt(self, filename, t_name=None, I_name=None, DoLP_name=None, AoLP_name=None, **kwargs):

		data = pd.read_csv(filename, **kwargs)
		# print("DEBUG NEW")
		# print(data)
		# print(data.columns)
		if t_name:
			t = data[t_name]
			if "time" in t_name:
				try:
					t = np.array([dt.datetime.strptime(t, "%Y%m%d-%H%M%S") for t in t])
				except:
					pass

		if I_name: 		I = data[I_name]
		if DoLP_name: 	DoLP = data[DoLP_name]
		if AoLP_name: 	AoLP = data[AoLP_name]

		# if len(t) == 1:
		# 	t = self.x_axis_list
		# 	I = np.array([I[0]] * len(t))
		# 	DoLP = np.array([DoLP[0]] * len(t))
		# 	AoLP = np.array([AoLP[0]] * len(t))

		return t, I, DoLP, AoLP


	def GetVParamFromLightParam(self, I0, DoLP, AoLP):
		"""Returns  the V parameters for a given radiant flux (any unit), a DoLP (between 0 and 1) and an AoLP in radians"""
		V = I0 / 2.
		Vcos = I0 * DoLP * np.cos(2 * AoLP) / 4.
		Vsin = I0 * DoLP * np.sin(2 * AoLP) / 4.
		return V, Vcos, Vsin

	def GetLightParameters(self, V, Vcos, Vsin):
		"""Given V, Vcos, Vsin, returns the initial intensity, DoLP and AoLP. This method is shared for spp and ptcu. It is also a static method, that can be called outside of the object. This way it can be used everywhere, each time you need I0, DoLP, AoLP to decrease the chances of making a mistake."""
		I0 = 2 * V
		DoLP = 2 * np.sqrt(Vcos**2 + Vsin**2) / V
		AoLP = np.arctan2(Vsin, Vcos) / 2
		return abs(I0), abs(DoLP), AoLP

	def AddLights(self, I1, D1, A1, I2, D2, A2, t1 = None, t2 = None):

		V1, Vc1, Vs1 = self.GetVParamFromLightParam(I1, D1, A1)
		V2, Vc2, Vs2 = self.GetVParamFromLightParam(I2, D2, A2)

		if t1 is None and t2 is None:
			return self.GetLightParameters(V1+V2, Vc1+Vc2, Vs1+Vs2)
		elif t1 is None:
			V1  = V1  + np.zeros_like(V2)
			Vc1 = Vc1 + np.zeros_like(Vc2)
			Vs1 = Vs1 + np.zeros_like(Vs2)
		elif t2 is None:
			V2  = V2  + np.zeros_like(V1)
			Vc2 = Vc2 + np.zeros_like(Vc1)
			Vs2 = Vs2 + np.zeros_like(Vs1)
		else:
			V2  = np.interp(t1, t2, V2)
			Vc2 = np.interp(t1, t2, Vc2)
			Vs2 = np.interp(t1, t2, Vs2)

		return self.GetLightParameters(V1+V2, Vc1+Vc2, Vs1+Vs2)

	def SubtractModel(self, bottle, model):
		# x_axis = self.x_axis_list
		#
		# m_V = np.interp(x_axis, model.x_axis, model["RSV"])
		# m_Vcos = np.interp(x_axis, model.x_axis, model["RSVcos"])
		# m_Vsin = np.interp(x_axis, model.x_axis, model["RSVsin"])
		#
		# diff_V = bottle.smooth_V - m_V
		# diff_Vcos = bottle.smooth_Vcos - m_Vcos
		# diff_Vsin = bottle.smooth_Vsin - m_Vsin


		x_axis = model.x_axis

		b_V = np.interp(x_axis, self.x_axis_list, bottle.smooth_V)
		b_Vcos = np.interp(x_axis, self.x_axis_list, bottle.smooth_Vcos)
		b_Vsin = np.interp(x_axis, self.x_axis_list, bottle.smooth_Vsin)

		diff_V = b_V - model["RSV"]
		diff_Vcos = b_Vcos - model["Vcos"]
		diff_Vsin = b_Vsin - model["Vsin"]

		diff_I0, diff_DoLP, diff_AoLP = Rotation.GetLightParameters(diff_V, diff_Vcos, diff_Vsin)

		f, axs = plt.subplots(3, sharex = True)
		axs[0].plot(x_axis, diff_I0, "k")
		axs[1].plot(x_axis, diff_DoLP, "b")
		axs[2].plot(x_axis, diff_AoLP * RtoD, "r")