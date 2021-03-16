#!/usr/bin/python3
# -*-coding:utf-8 -*

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
#matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage{lmodern}']
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=False)
matplotlib.rcParams['font.size'] = 16
matplotlib.rcParams['svg.fonttype'] = 'none'
import datetime as dt

import chaosmagpy as chaos

from utils import *
from rotation import *
from bottle import *
from AllSkyData import *
from eiscat_data import *
from MagData import *
from J2RAYS1_model import *
from scipy import signal
import sys
import os
from subprocess import call

class Mixer:
	def __init__(self, bottles, mag_data=False, comp_bottles=[], mode = "default"):

		self.bottles = bottles
		self.comp_bottles = comp_bottles
		self.mag_data = mag_data

		for ib,  bottle in enumerate(self.bottles):
			if self.mag_data is not False:
				self.mag_data.StripTime(bottle.DateTime("start"), bottle.DateTime("end"))

			self.allsky_data_available = False
			if bottle.observation_type == "fixed":
				self.allsky_data = AllSkyData(bottle)
				self.allsky_data_available = bool(self.allsky_data.all_datetimes)

			self.eiscat_data = Eiscat(bottle)
			# print("DEBUG:", self.eiscat_data.valid )

			self.eq_currents = EqCurrent(bottle, file_type=None)
			if self.eq_currents.valid and bottle.observation_type == "fixed":

				AoLPs = [bottle.GetInterpolation(t)[2] for t in self.eq_currents.GetNormTimes()]
				# self.eq_currents.FindJup(bottle.observation, AoLPs)

				self.eq_currents.GetApparentAngle(bottle.observation)#, shift=bottle.graph_angle_shift)

				self.eq_currents.SaveTXT()

			self.SetGraphParameter(bottle)
			self.MakeXAxis(bottle)

			if self.show_RS_model and self.RS_model_files:
				self.J2RAYS1_models = [Model.InitFromBottle(bottle, f, time_divisor = self.divisor, x_axis = self.x_axis_list, shift = self.shift_model) for i, f in enumerate(self.RS_model_files[bottle.line - 1])]

				# self.J2RAYS1_models[0]["DoLP"] = self.J2RAYS1_models[0]["DoLP"]  / 10.
				# self.J2RAYS1_models[0].SetVParam()
				# self.J2RAYS1_models[1]["DoLP"] = self.J2RAYS1_models[1]["DoLP"]  * 2
				# self.J2RAYS1_models[1].SetVParam()

				if self.fit_RS_model_to_flux and len(self.J2RAYS1_models) >= 2:
					print("Fitting J2RAYS-1 models to data...")
					# best_model = Model.FitModelToFlux(self.x_axis_list, bottle.smooth_I0, self.J2RAYS1_models[0], self.J2RAYS1_models[1])
					best_model = Model.FitModelToFluxPlusIso(self.x_axis_list, bottle.smooth_I0, self.J2RAYS1_models[0], self.J2RAYS1_models[1])

					# best_model = Model.FitModelToDoLP(self.x_axis_list, bottle.smooth_DoLP, self.J2RAYS1_models[0], self.J2RAYS1_models[1])
					# best_model = Model.FitModelToFluxAndDoLP(self.x_axis_list, bottle.smooth_I0, bottle.smooth_DoLP, self.J2RAYS1_models[0], self.J2RAYS1_models[1])

					# best_model = Model.FitModelToAoLP(self.x_axis_list, bottle.smooth_AoLP, self.J2RAYS1_models[0], self.J2RAYS1_models[1])
					# best_model = Model.FitModelToFluxRatio(self.x_axis_list, bottle.smooth_I0, self.J2RAYS1_models[0], self.J2RAYS1_models[1])
					# self.J2RAYS1_models = [best_model]
					best_model.SetMandatoryColumns()
					best_model.SetArbitraryUnits(bottle = bottle)
					self.J2RAYS1_models.append(best_model)
					# print("DEBUG best_model", self.J2RAYS1_models[0].data["DV"], self.J2RAYS1_models[0].data["DVcos"])#, "RSDoLP", "RSAoRD", "RSV", "RSVcos", "RSVsin"])
					# print("DEBUG best_model", self.J2RAYS1_models[1].data["DV"], self.J2RAYS1_models[1].data["DVcos"])#, "RSDoLP", "RSAoRD", "RSV", "RSVcos", "RSVsin"])
					# print("DEBUG best_model", best_model.data["RSV"], best_model.data["RSVcos"], best_model.data["RSDoLP"])
					# self.J2RAYS1_models.append((self.J2RAYS1_models[0] * param[0]) + (self.J2RAYS1_models[1] * param[1]))

				if self.AddDirectPola and len(self.J2RAYS1_models) >= 2:
					print("Adding direct polarisation to model...")
					if self.fit_RS_model_to_flux: # MOdel 0 and 1 are grd and sky. 2 is the direct light only and 3 is the Linear Comb of 0 and 1
						direct_model = self.J2RAYS1_models[2]
						base_model = self.J2RAYS1_models[3]
					else: # Model 0 is the diffusion model. 1 is the direct only model
						base_model = self.J2RAYS1_models[0]
						direct_model = self.J2RAYS1_models[1]

					for dolp_factor in self.AddDirectPola:
						self.J2RAYS1_models.append(base_model.AddDirectPola(direct_model["DDoLP"] * dolp_factor, direct_model["DAoLP"]))
						self.J2RAYS1_models[-1].SetMandatoryColumns()
						self.J2RAYS1_models[-1].SetArbitraryUnits(bottle = bottle)

				for i, m in enumerate(self.J2RAYS1_models):
					try:
						addF, addD, addA = self.addFDA_to_model[bottle.line - 1][i]
						m.AddFDA(addF, addD, addA)
						print("DEBUG: add FDA to model", addF, addD, addA)
					except:
						print("DEBUG: DID NOT add FDA to model", bottle.line - 1, i)
						pass
						# addF, addD, addA = 0, 0, 0
			self.mode = mode

			self.MakeFigure()
			self.MakePlots(bottle)
			for m in self.J2RAYS1_models:
				self.SubtractModel(bottle, m)
			# if self.show_RS_model:
			# 	self.CompareRayleigh(bottle, self.J2RAYS1_model)
			if self.make_optional_plots:
				self.MakeSNratio(bottle)
				self.MakeI0SNratio(bottle)
				if self.eq_currents.valid:
					self.eq_currents.MakePlot()

			if self.comp_bottles:
				self.SetGraphParameter(self.comp_bottles[ib], comp=True)
				self.MakeXAxis(self.comp_bottles[ib])
				self.MakePlots(self.comp_bottles[ib])
				if self.make_optional_plots:
					self.MakeSNratio(self.comp_bottles[ib])

					self.SetGraphParameter(bottle)
					self.MakeXAxis(bottle)

					self.CompareBottles(self.comp_bottles[ib], bottle)

			if self.make_optional_plots:
				self.MakeFFT(bottle)

			if self.make_optional_plots and self.mag_data is not False:
				self.MakeMagDataPlots(bottle)
				self.MakeCorrelationPlots(bottle)

			if self.make_optional_plots and self.eiscat_data.valid:
				self.MakeEiscatDataPlot(bottle)

			if self.make_optional_plots and bottle.observation_type == "fixed_elevation_continue_rotation":
				self.CompareAnglesPlots(bottle)

				if len(bottle.continue_rotation_times) > 2:
					print("AVERAGE ROTATIONS")
					bottle = bottle.GetAverageContinueRotations()
					self.SetGraphParameter(bottle)
					self.MakeXAxis(bottle)
					self.MakeFigure()
					self.MakePlots(bottle)
					self.MakeSNratio(bottle)
					self.CompareAnglesPlots(bottle)



			# self.MakeCoherencePLots()


	def MakeXAxis(self, bottle):

		self.divisor = 1.
		if self.xaxis_azimut and bottle.observation_type == "fixed_elevation_continue_rotation":
			self.x_axis_list = np.array(())
			for ir in range(bottle.nb_continue_rotation):
				start = dt.timedelta(minutes=bottle.continue_rotation_times[ir])
				end = dt.timedelta(minutes=bottle.continue_rotation_times[ir+1])
				nb_points = len([t for t in bottle.all_times if start <= t <= end ])

				self.x_axis_list = np.append(self.x_axis_list, np.linspace(360 * ir, 360 * (ir+1), nb_points))
			# print(len(self.x_axis_list), len(bottle.all_times))
			self.xlabel = "Azimuth"

			if bottle.nb_continue_rotation <= 2:
				self.x_axis_ticks_pos = np.arange(0, 360 * bottle.nb_continue_rotation, 45)
				self.x_axis_ticks_label = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"] *  bottle.nb_continue_rotation
			else:
				self.x_axis_ticks_pos = np.arange(0, 360 * bottle.nb_continue_rotation, 90)
				self.x_axis_ticks_label = ["N", "E", "S", "W"] *  bottle.nb_continue_rotation

		else:
			delta = bottle.all_times[-1] - bottle.all_times[0]
			if delta > time.timedelta(hours=2):
				self.divisor = 3600.
				self.xlabel = "Time (hours)"
			elif delta > time.timedelta(minutes=2):
				self.divisor = 60.
				self.xlabel = "Time (minutes)"
			else:
				self.xlabel = "Time (seconds)"

			# self.x_axis_list 	   = np.array([t.total_seconds() for t in bottle.all_times_since_start]) / self.divisor
			norm = bottle.all_times[0].total_seconds()
			# self.x_axis_list = np.array([t.total_seconds() - norm for t in bottle.all_times]) / self.divisor
			self.x_axis_list = np.array([t.total_seconds() - norm for t in bottle.all_times]) / self.divisor

			if self.xaxis_azimut and bottle.observation_type == "fixed_elevation_discrete_rotation":
				self.xlabel = "Azimuth (°)"
				self.x_axis_ticks_pos = bottle.discrete_rotation_times[::9]
				self.x_axis_ticks_label = np.round(bottle.discrete_rotation_azimuts * RtoD % 360)[::9]
				self.x_axis_ticks_label = [int(x) for x in self.x_axis_ticks_label]

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

	def SetGraphParameter(self, bottle, comp = False):
		self.marker_size = 3
		self.single_star_size = self.marker_size*5

		self.show_legend = False

		self.show_allsky = False
		self.show_eiscat = False

		self.show_mag_data 	= False
		self.show_AoBapp 	= False
		self.show_AoRD	 	= False
		self.show_currents	= True


		self.make_optional_plots = False

		# tmp_path = "/home/bossel/These/Documentation/Mes_articles/scattering/new_pictures/FIG8/"
		# self.RS_model_files	= [	tmp_path + "v_grd_only.txt",
		# 						tmp_path + "v_sky_alb.txt",
		# 						tmp_path + "v_sky_albx1+v_grd_onlyx2e-09",
		# 						tmp_path + "v_sky_albx1+v_direct_only_Bx1"]#,
								# tmp_path + "m_sky_albx1+m_grd_onlyx2e-09"]#,
								# tmp_path + "b_sky_albx1+b_grd_onlyx5e-10",
								# tmp_path + "b_sky_albx1+b_grd_onlyx1e-10"]

		# tmp_path = "/home/bossel/These/Documentation/Mes_articles/scattering/new_pictures/FIG7/"
		# self.RS_model_files	= [	tmp_path + "m_rot_e45_grd_only.txt",
		# 						tmp_path + "m_rot_e45_sky_albedo.txt",
		# 						tmp_path + "m_rot_e45_grd_onlyx1+m_rot_e45_sky_albedox12000000000.0",
		# 						tmp_path + "m_rot_e45_grd_onlyx1+m_rot_e45_sky_albedox12000000000.0"]

		# tmp_path = "/home/bossel/These/Documentation/Mes_articles/scattering/new_pictures/FIG6/"
		# self.RS_model_files	= [	tmp_path + "o_grd_only_norm0.txt",
		# 						tmp_path + "o_grd_only_norm0.txt"]#,
								# tmp_path + "t_sky_alb.txt"]

		# tmp_path = "/home/bossel/These/Documentation/Mes_articles/scattering/new_pictures/FIG5/"
		# self.RS_model_files	= [	[tmp_path + "2021_lagorge_v_grd_only.txt", tmp_path + "2021_lagorge_v_grd_only_aero_mar.txt", tmp_path + "2021_lagorge_v_grd_only_FULL_aero_mar.txt", tmp_path + "2021_lagorge_v_grd_only_FULL_aero_mar_4000.txt", tmp_path + "2021_lagorge_v_grd_only_FULL_aero_mar_10000.txt", tmp_path + "2021_lagorge_v_grd_only_FULL_aero_mar_n4000_250nm.txt", tmp_path + "2021_lagorge_v_grd_only_FULL_aero_mar_n4000_200nm.txt", tmp_path + "2021_lagorge_v_grd_only_FULL_aero_mar_n4000_150nm.txt"]]#,#,
		### Bosse et al 2021 Fig 10+
		# self.RS_model_files	= [	[tmp_path + "2021_lagorge_b_grd_only_FULL_aero_mar_n4000_150nm.txt"],
		# 						[tmp_path + "2021_lagorge_o_grd_only_FULL_aero_mar_n4000_150nm.txt"],
		# 						[tmp_path + "2021_lagorge_m_grd_only_FULL_aero_mar_n4000_150nm.txt"],
		# 						[tmp_path + "2021_lagorge_t_grd_only_FULL_aero_mar_n4000_150nm.txt"]]

		### Bosse et al 2021 Fig 9
		# self.RS_model_files	= [	[tmp_path + "2021_lagorge_v_grd_only_FULL_aero_mar_n4000_150nm.txt"]]
		### Bosse et al 2021 Fig 8
		# self.RS_model_files	= [	[tmp_path + "2021_lagorge_v_grd_only.txt", tmp_path + "2021_lagorge_v_grd_only.txt", tmp_path + "2021_lagorge_v_grd_only_FULL_aero_mar_n4000_150nm.txt", tmp_path + "2021_lagorge_v_grd_only_FULL_aero_rur_1000.txt", tmp_path + "2021_lagorge_v_grd_only_FULL_aero_rur_500.txt"]]


		tmp_path = "/home/bossel/These/Documentation/Mes_articles/auroral_pola/"
		### 20200227 XYvm e52 best fit
		# self.RS_model_files	= [[tmp_path + "grd_only/o_e52_grd_only_aero_1low.txt",
		# 						tmp_path + "sky_only/uni_sky_o_e52_aero_1low_albedo.txt"],
		# 						[tmp_path + "grd_only/t_e52_grd_only_aero_1low.txt",
		# 						tmp_path + "sky_only/uni_sky_t_e52_aero_1low_albedo.txt",
		# 						tmp_path + "t_e52_grd_only_aero_1lowx1702.6351068301055+uni_sky_t_e52_aero_1low_albedox63627.19976961118"], #MANUAL fit because of the East flux jump. Coeff are : 4335.3242154348345 * 0.39273535777746454 and 25450.87990784447 * 2.5.
		# 						# [tmp_path + "grd_only/t_e52_grd_only_aero_1low.txt",
		# 						# tmp_path + "sky_only/uni_sky_t_e52_aero_1low_albedo.txt"],
		# 						[tmp_path + "grd_only/v_e52_grd_only_aero_1low.txt",
		# 						tmp_path + "sky_only/20200227_v_e52_sky_only_aero_1low_albedo.txt"],
		# 						[tmp_path + "grd_only/m_e52_grd_only_aero_1low.txt",
		# 						tmp_path + "sky_only/20200227_m_e52_sky_only_aero_1low_albedo.txt"]]

		### 20200227 XYvm e52 test aerosols
		# self.RS_model_files	= [[tmp_path + "grd_only/o_e52_grd_only_NO_aero.txt",
		# 						tmp_path + "grd_only/o_e52_grd_only_aero_arctic.txt",
		# 						tmp_path + "grd_only/o_e52_grd_only_aero_1low.txt",
		# 						tmp_path + "grd_only/o_e52_grd_only_aero_3mid.txt",
		# 						tmp_path + "grd_only/o_e52_grd_only_aero_2high.txt"],
		# 						[tmp_path + "grd_only/t_e52_grd_only_aero_1low.txt"],
		# 						[tmp_path + "grd_only/v_e52_grd_only_aero_1low.txt"],
		# 						[tmp_path + "grd_only/m_e52_grd_only_aero_1low.txt"]]
		# self.RS_model_files	= [[tmp_path + "CL/o_e52_grd_only_aero_1lowx1+uni_sky_o_e52_aero_1low_albedox6.0"],
								# [],[],[]]

		# ### 20190307 br rot e45
		# self.shift_model = -30
		# self.RS_model_files	= [[tmp_path + "grd_only/b_grd_only_aero_1low.txt",
		# 						tmp_path + "sky_only/20190307_b_sky_only_aero_1low_albedo.txt"],
		# 						[],[],[]]
		### 20190307 mr rot e45
		# self.shift_model = -30
		# self.RS_model_files	= [[tmp_path + "grd_only/m_grd_only_aero_1low.txt",
		# 						tmp_path + "sky_only/20190307_m_sky_only_aero_1low_albedo.txt"],
		# 						[],[],[]]
		### 20190307 vr rot e45
		# self.RS_model_files	= [[tmp_path + "grd_only/v_grd_only_aero_1low.txt",
		# 						tmp_path + "sky_only/20190307_v_sky_only_aero_1low_albedo.txt"],
		# 						# tmp_path + "sky_only/20190307_v_sky_only_aero_1low_NOalbedo.txt"],
		# 						[],[],[]]
		## 20190307 vr a164 e45
		# self.RS_model_files	= [[tmp_path + "grd_only/v_a164_grd_only_aero_1low.txt",
		# 						tmp_path + "sky_only/v_a164_sky_only_aero_1low.txt"],
		# 						[],[],[]]
		## 20190307 br a164 e45
		# self.RS_model_files	= [[tmp_path + "grd_only/b_a164_grd_only_aero_1low.txt",
		# 						tmp_path + "sky_only/b_a164_sky_only_aero_1low.txt"],
		# 						[],[],[]]
		# ## 20190307 mr a164 e45
		# self.RS_model_files	= [[tmp_path + "grd_only/m_a164_grd_only_aero_1low.txt",
		# 						tmp_path + "sky_only/m_a164_sky_only_aero_1low.txt"],
		# 						[],[],[]]

		### 20200227 XYvm e52 Direct Pola test
		self.RS_model_files	= [[],
								[],
								[tmp_path + "grd_only/v_e52_grd_only_aero_1low.txt",
								 tmp_path + "sky_only/20200227_v_e52_sky_only_aero_1low_albedo.txt",
								 tmp_path + "sky_only/20200227_v_e52_sky_only_aero_1low_albedo_DpolaNW_cos.txt"],
								[tmp_path + "grd_only/m_e52_grd_only_aero_1low.txt",
								 tmp_path + "sky_only/20200227_m_e52_sky_only_aero_1low_albedo.txt",
								 tmp_path + "sky_only/20200227_v_e52_sky_only_aero_1low_albedo_DpolaNW_cos.txt"]]
		# self.RS_model_files	= [[],
		# 						[],
		# 						[tmp_path + "sky_only/20200227_v_e52_sky_only_aero_1low_albedo.txt",
		# 						 tmp_path + "sky_only/20200227_v_e52_sky_only_aero_1low_albedo_DpolaNW_cos.txt"],
		# 						[tmp_path + "sky_only/20200227_m_e52_sky_only_aero_1low_albedo.txt",
		# 						 tmp_path + "sky_only/20200227_v_e52_sky_only_aero_1low_albedo_DpolaNW_cos.txt"]]

		self.show_RS_model		= True
		self.fit_RS_model_to_flux = True

		self.AddDirectPola = [1, 2]

		self.adapt_flux_scale 	= True
		self.shift_model = 0
		self.model_colors = ["black", "red", "blue", "magenta", "purple",  "yellow", "pink", "orange", "magenta", "cyan"] * 10
		self.model_symbols = ["*", "+", "x", "1", "2", "3", "4", ".", "s", "<", "^", ">", "X", "D"] * 10
		self.addFDA_to_model = [[(0, 0, 0), (0, 0, 0)]]  * len(self.RS_model_files)
		# self.addFDA_to_model = np.array([[(0.3, 0, 0), (2, 0, 0), (20, 0, 0), (18, 0, 0), (21, 0, 0)]])*1.5 # 20200227 rotation XYvm à skibotn pour test d'aerosols fit le median
		# self.addFDA_to_model = np.array([[(0.26, 0, 0), (2.14, 0, 0), (33.2, 0, 0), (30, 0, 0), (43, 0, 0)]]) # 20200227 rotation XYvm à skibotn pour test d'aerosols fit le min
		# self.addFDA_to_model = [[], [], [], [(0, 0, 0), (0, 0, 0), (0, 0, 0)]] # 20200227 rotation XYvm à skibotn pour test de background dans le mauve (fct:average)

		# self.addFDA_to_model = [[(190, 0, 0)]] #Bosse et al 2021 Fig 9
		# self.addFDA_to_model = [[(0, 0, 0), (66.662415, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0)]] # 20210120 rotation verte à lagorge avec 3 profils d'aerosols et background. Bosse et al 2021 Fig 8

		# self.addFDA_to_model = [[(10, 0, 0)], [(21, 0, 0)], [(2, 0, 0)], [(20, 0, 0)]] # 20210120_Lagorge_bomt_rot_lente NEW AEROSOLS

		self.show_Ipola = False
		self.show_Iref = False

		self.show_time = not comp
		self.time_format = "UT" #LT or UT
		self.time_label = "UTC"


		self.show_raw_data = False
		self.show_smooth_data = True

		self.show_error_bars 		= False
		self.show_smooth_error_bars = True
		self.max_error_bars = 10000 #If there are too many data, takes very long to plot error bars. If != 0, will plot max_error_bars error bars once every x points.
		self.show_SN = False

		self.show_avg_I0 = False
		self.show_avg_Iref = False
		self.show_avg_DoLP = False
		self.show_avg_AoLP = False

		#If True Will show the azimut on the xaxis instead of time for rotations
		self.xaxis_azimut = True

		self.SetColors(bottle, comp=comp)

		font = {'weight' : 'bold',
        'size'   : 24}
		matplotlib.rc('font', **font)

		# plt.tight_layout(pad=0, h_pad=None, w_pad=None, rect=None)


	def SetColors(self, bottle, comp = False):

		self.all_I0_color = "xkcd:black"
		self.smooth_I0_color = "xkcd:red"

		self.raw_error_bars_color = "grey"

		self.smooth_ref_color = "xkcd:green"
		self.AoBapp_color = "xkcd:turquoise"
		self.AoRD_color = "xkcd:hot pink"
		self.currents_color = "xkcd:blue"
		if comp:
			self.AoBapp_color = "xkcd:blue"
			self.AoRD_color = "xkcd:salmon"

		self.EISCAT_color = "xkcd:mustard"
		self.mag_color = "xkcd:mustard"

		self.all_sky_color = "xkcd:orange"

		### xkcd color guide: https://xkcd.com/color/rgb/
		print("MIXER FILTER:", bottle.filters)
		if bottle.filters:
			self.pola_color = bottle.filters[0]
			if comp:
				self.smooth_I0_color = "xkcd:black"
			else:
				if 	 self.pola_color == "r": self.smooth_I0_color = "xkcd:red"
				elif self.pola_color == "v": self.smooth_I0_color = "xkcd:green"
				elif self.pola_color == "b": self.smooth_I0_color = "xkcd:blue"
				elif self.pola_color == "m": self.smooth_I0_color = "xkcd:purple"
				elif self.pola_color == "o": self.smooth_I0_color = "xkcd:orange"
				elif self.pola_color == "t": self.smooth_I0_color = "xkcd:turquoise"
				elif self.pola_color == "X": self.smooth_I0_color = "xkcd:orange"
				elif self.pola_color == "Y": self.smooth_I0_color = "xkcd:turquoise"
				else: self.smooth_I0_color = "red"

			self.smooth_error_bars_color = self.smooth_I0_color

			if bottle.instrument_name == "carmen" and bottle.filters[1] != 0:
				self.ref_color = bottle.filters[1]
				if 	 self.ref_color == "r" and self.pola_color != "r": self.smooth_ref_color = "xkcd:red"
				elif self.ref_color == "r" and self.pola_color == "r": self.smooth_ref_color = "xkcd:orange"
				elif self.ref_color == "v" and self.pola_color != "v": self.smooth_ref_color = "xkcd:green"
				elif self.ref_color == "v" and self.pola_color == "v": self.smooth_ref_color = "xkcd:lime green"
				elif self.ref_color == "b" and self.pola_color != "b": self.smooth_ref_color = "xkcd:blue"
				elif self.ref_color == "b" and self.pola_color == "b": self.smooth_ref_color = "xkcd:bright blue"
				elif self.ref_color == "m" and self.pola_color != "m": self.smooth_ref_color = "xkcd:purple"
				elif self.ref_color == "m" and self.pola_color == "m": self.smooth_ref_color = "xkcd:lavender"
				elif self.ref_color == "o" and self.pola_color != "o": self.smooth_ref_color = "xkcd:orange"
				elif self.ref_color == "o" and self.pola_color == "o": self.smooth_ref_color = "xkcd:orange"
				else: self.smooth_ref_color = "green"


		self.all_SN_color = "black"
		self.smooth_SN_color = self.smooth_I0_color

	def MakePlots(self, bottle):
		#Plotting the mean intensity I0, DoLP, AoLP for each rotation and more!

		# if self.show_error_bars:
		# 	smooth_yerr = bottle.var_smooth_I0
		# 	yerr = bottle.var_I0
		# 	yerr = None
		# else:
		# 	smooth_yerr = None

		print("START PLOTTING")

		if self.show_raw_data:
			if not self.show_error_bars:
				l_all_I0, = self.ax1.plot(self.x_axis_list, bottle.all_I0, ".", color = self.all_I0_color, linestyle = 'none', markersize=self.marker_size, label="Intensity", zorder=0)
			else:
				l_all_I0 = self.ax1.errorbar(self.x_axis_list, bottle.all_I0, yerr = bottle.std_I0, fmt=".", ecolor=self.raw_error_bars_color, color = self.all_I0_color , linestyle = 'none', markersize=self.marker_size, label="Intensity", zorder=0, errorevery=self.error_step)


		if self.show_smooth_data and not self.show_smooth_error_bars:
			# l_smooth_I0, = self.ax1.plot(self.x_axis_list, bottle.smooth_I0, ".", color = self.smooth_I0_color, linestyle = 'none', markersize=self.marker_size, label="Smooth Intensity (" + str(bottle.smoothing_factor) + " " + str(bottle.smoothing_unit)[:3] + ")", zorder=1)
			y = bottle.smooth_I0
			if self.show_Ipola:
				y *= bottle.smooth_DoLP / 100.
			l_smooth_I0, = self.ax1.plot(self.x_axis_list, y, ".", color = self.smooth_I0_color, linestyle = 'none', markersize=self.marker_size, label="Smooth Intensity (" + str(bottle.smoothing_factor) + " " + str(bottle.smoothing_unit)[:3] + ")", zorder=1)
		elif self.show_smooth_data:
			# l_smooth_I0 = self.ax1.errorbar(self.x_axis_list, bottle.smooth_I0, yerr = bottle.std_smooth_I0, fmt=".", color = self.smooth_I0_color, ecolor=self.smooth_error_bars_color, linestyle = 'none', markersize=self.marker_size, label="Smooth Intensity (" + str(bottle.smoothing_factor) + " " + str(bottle.smoothing_unit)[:3] + ")", zorder=1, errorevery=self.error_step)
			y = bottle.smooth_I0
			if self.show_Ipola:
				y *= bottle.smooth_DoLP / 100.
			l_smooth_I0 = self.ax1.errorbar(self.x_axis_list, y, yerr = bottle.std_smooth_I0, fmt=".", color = self.smooth_I0_color, ecolor=self.smooth_error_bars_color, linestyle = 'none', markersize=self.marker_size, label="Smooth Intensity (" + str(bottle.smoothing_factor) + " " + str(bottle.smoothing_unit)[:3] + ")", zorder=1, errorevery=self.error_step)


		if self.show_avg_I0:
			l_avg_I0, = self.ax1.plot(self.x_axis_list,[bottle.I0_average] * len(self.x_axis_list), color = self.smooth_I0_color, label="Avg: " + str(bottle.I0_average)[:4])
			self.ax1_lines.append([l_avg_I0, l_avg_I0.get_label()])

		if self.show_raw_data:
			self.ax1_lines.append([l_all_I0, l_all_I0.get_label()])
		if  self.show_smooth_data:
			self.ax1_lines.append([l_smooth_I0, l_smooth_I0.get_label()])


		# ### Graph the I0 * DoLP line on the first subplot
		# self.ax14 = self.ax1.twinx()
		# # l_all_IDoLP, = self.ax14.plot(self.x_axis_list, bottle.all_I0 * bottle.all_DoLP, "k.", linestyle = 'none', markersize=self.marker_size, label="All Intensity * DoLP", zorder=2)
		# smooth_IDoLP, = self.ax14.plot(self.x_axis_list, bottle.smooth_I0 * bottle.smooth_DoLP / 100. * bottle.smooth_AoLP
		#
		# cd, ".", color = "xkcd:hot pink", linestyle = 'none', markersize=self.marker_size, label="Smooth Intensity * DoLP", zorder=2)
		#
		# self.ax14.set_ylabel("Ref intensity")



		if not bottle.NoVref and bottle.instrument_name == "carmen" and self.show_Iref:
			self.ax13 = self.ax1.twinx()
			# ax13.set_xlim(ax1.get_xlim())
			# xticks = plt.xticks()[0] * self.divisor
			# xticks = [t + bottle.time + bottle.head_jump for t in xticks]
			# ax12.set_xticklabels([time.strftime("%H:%M:%S", time.localtime(st)) for st in xticks])

			# l_Iref, = self.ax13.plot(self.x_axis_list[1:], bottle.all_Iref[1:], ".", color = "black", linestyle = 'none', markersize=self.marker_size, label="Ref Intensity", zorder=2)
			# self.ax1_lines.append([l_Iref, l_Iref.get_label()])

			l_smooth_Iref, = self.ax13.plot(self.x_axis_list, bottle.smooth_Iref, ".", color = self.smooth_ref_color, linestyle = 'none', markersize=self.marker_size, label="Smooth Ref Intensity (" + str(bottle.smoothing_factor) + " " + str(bottle.smoothing_unit)[:3] + ")", zorder=2)
			self.ax1_lines.append([l_smooth_Iref, l_smooth_Iref.get_label()])

			if self.show_avg_Iref:
				l_avg_Iref, = self.ax13.plot(self.x_axis_list, [bottle.Iref_average] * len(self.x_axis_list), color = self.smooth_ref_color, label="Avg Ref Intensity " + str(bottle.Iref_average)[:4], zorder=2)
				self.ax1_lines.append([l_avg_Iref, l_avg_Iref.get_label()])
			self.ax13.set_ylabel("Ref intensity")



		if self.allsky_data_available and self.show_allsky:
			self.ax14 = self.ax1.twinx()
			offset = 10
			# new_fixed_axis = self.ax14.get_grid_helper().new_fixed_axis
			# self.ax14.axis["right"] = new_fixed_axis(loc="right", axes=par2, offset=(offset, 0))
			#
			# self.ax14.axis["right"].toggle(all=True)

			l_ASI, = self.ax14.plot(self.allsky_data.GetNormTimes(bottle.DateTime("start", format=self.time_label), self.divisor), self.allsky_data.brightness, "*", color = "orange", linestyle = 'none', markersize=self.marker_size, label="AllSKy Imager", zorder=2)
			self.ax1_lines.append([l_ASI, l_ASI.get_label()])
			# print(self.ax14.get_yticklabels())
			# self.ax14.set_yticklabels([l.get_text() for l in self.ax14.get_yticklabels()], horizontalalignment = "left")
			# self.ax14.tick_params(direction='in', labelright=True, pad=-5)
			self.ax14.set_ylabel("ASI Brightness")



		# if self.show_error_bars:
		# 	smooth_yerr = bottle.var_smooth_DoLP
		# 	yerr = bottle.var_DoLP
		# 	yerr = None
		# else:
		# 	smooth_yerr = None

		if self.show_raw_data:
			if not self.show_error_bars:
				l_all_DoLP, = self.ax2.plot(self.x_axis_list, bottle.all_DoLP, ".", color = self.all_I0_color, linestyle = 'none', markersize=self.marker_size, label="DoLP", zorder=0)
				self.ax2_lines.append([l_all_DoLP, l_all_DoLP.get_label()])
			else:
				(l_all_DoLP, _, _) = self.ax2.errorbar(self.x_axis_list, bottle.all_DoLP, yerr = bottle.std_DoLP, fmt=".", ecolor=self.raw_error_bars_color, color = self.all_I0_color, linestyle = 'none', markersize=self.marker_size, label="DoLP", zorder=0, errorevery=self.error_step)

		if self.show_smooth_data and not self.show_smooth_error_bars:
			l_smooth_DoLP, = self.ax2.plot(self.x_axis_list, bottle.smooth_DoLP, ".", color = self.smooth_I0_color, linestyle = 'none', markersize=self.marker_size, label="Smooth DoLP (" + str(bottle.smoothing_factor) + " " + str(bottle.smoothing_unit)[:3] + ")", zorder=1)
		elif  self.show_smooth_data:
			(l_smooth_DoLP, _, _) = self.ax2.errorbar(self.x_axis_list, bottle.smooth_DoLP, yerr = bottle.std_smooth_DoLP, fmt=".", color = self.smooth_I0_color, ecolor=self.smooth_error_bars_color, linestyle = 'none', markersize=self.marker_size, label="Smooth DoLP (" + str(bottle.smoothing_factor) + " " + str(bottle.smoothing_unit)[:3] + ")", zorder=1, errorevery=self.error_step)

		self.ax2.set_ylim(bottom = 0)

		if self.show_raw_data:
			self.ax2_lines.append([l_all_DoLP, l_all_DoLP.get_label()])
		if  self.show_smooth_data:
			self.ax2_lines.append([l_smooth_DoLP, l_smooth_DoLP.get_label()])

		if bottle.location.lower() == "skibotn" and bottle.filters == "br" and bottle.DateTime().date() == dt.date(2019, 3, 7):
			self.ax2.set_ylim((0, 5))

		if self.show_avg_DoLP:
			l_avg_DoLP, = self.ax2.plot(self.x_axis_list,[bottle.DoLP_average] * len(self.x_axis_list), color = self.smooth_I0_color, label="Avg: " + str(bottle.DoLP_average)[:4])
			self.ax2_lines.append([l_avg_DoLP, l_avg_DoLP.get_label()])

		if self.mag_data and self.show_mag_data:
			self.ax21 = self.ax2.twinx()

			# print(self.mag_data.times[0], self.mag_data.times[-1])
			# print(self.mag_data.GetNormTimes(self.divisor)[0], self.mag_data.GetNormTimes(self.divisor)[-1])
			component = "Horiz"
			t, d = self.mag_data.GetComponent(component, self.divisor)
			# t, d = self.mag_data.GetDerivative(component, self.divisor)
			l_mag_data, = self.ax21.plot(t, d, "orange", label="B (nT)")
			# l_mag_data, = self.ax21.plot(t, d, self.mag_color, label="dB/dt (nT/s)", zorder=2, linewidth=2)
			# print(t, d)
			# Md, md = max(d), min(d)
			# self.ax21.set_ylim(min(-20, md, -abs(Md)), max(20, Md, abs(md)))
			self.ax2_lines.append([l_mag_data, l_mag_data.get_label()])
			self.ax21.set_ylabel("B (nT)")
			# self.ax21.set_ylabel("dB/dt (nT/s)")


		if self.eiscat_data.valid and self.show_eiscat:
			self.ax22 = self.ax1.twinx()

			parameter = "N_e"
			altitude = bottle.filters[0]
			if altitude == "r": altitude = 220
			elif altitude == "v": altitude = 110
			elif altitude == "b" or "m": altitude = 85
			label = parameter.replace("_", "\_") + " at " + str(altitude) + " km (EISCAT)"
			t, d, e = self.eiscat_data.GetParameter(parameter, altitude, time_divisor = self.divisor)
			l_eiscat_data, capline, barlinecol = self.ax22.errorbar(t, d, yerr = e, color = self.EISCAT_color, fmt = ".", label = label, zorder=2)
			# Md, md = max(d), min(d)
			if bottle.location.lower() == "skibotn" and bottle.filters == "vr" and bottle.DateTime().date() == dt.date(2019, 3, 7):
				self.ax22.set_ylim((10**11, 4.1*10**11))
			elif bottle.location.lower() == "skibotn" and bottle.filters == "vr" and bottle.DateTime().date() == dt.date(2019, 3, 8):
				self.ax22.set_ylim((-3*10**10, 5*10**10))
			self.ax1_lines.append([l_eiscat_data, label])
			self.ax22.set_ylabel(label)


		if self.show_raw_data:
			if not self.show_error_bars:
				l_all_AoLP, = self.ax3.plot(self.x_axis_list, bottle.all_AoLP * RtoD, ".", color = self.all_I0_color, linestyle = 'none', markersize=self.marker_size, label="AoLP", zorder=0)

		if self.show_smooth_data and not self.show_smooth_error_bars:
			l_smooth_AoLP, = self.ax3.plot(self.x_axis_list, bottle.smooth_AoLP * RtoD, ".", color = self.smooth_I0_color, linestyle = 'none', markersize=self.marker_size, label="Smooth AoLP (" + str(bottle.smoothing_factor) + " " + str(bottle.smoothing_unit)[:3] + ")", zorder=1)

		if (self.show_error_bars and self.show_raw_data) or self.show_smooth_error_bars:
			# self.ax3.fill_between(self.x_axis_list, bottle.all_AoLP * RtoD - bottle.std_AoLP, bottle.all_AoLP * RtoD + bottle.std_AoLP, color = "grey")
			# self.ax3.fill_between(self.x_axis_list, bottle.smooth_AoLP_lower * RtoD, bottle.smooth_AoLP_upper * RtoD, color = "yellow", alpha = 0.5)

			ls_all = dict()
			ls_smooth = dict()
			if bottle.graph_angle_shift == 1:
				min, mid, max = 0, 90, 180
			elif bottle.graph_angle_shift == 0:
				min, mid, max = -90, 0, 90

			for x, s_error, s_angle, error, angle in zip(self.x_axis_list[::self.error_step], bottle.std_smooth_AoLP[::self.error_step]*RtoD, bottle.smooth_AoLP[::self.error_step]*RtoD, bottle.std_AoLP[::self.error_step]*RtoD, bottle.all_AoLP[::self.error_step]*RtoD):

				if (self.show_error_bars and self.show_raw_data):
					if angle > mid and angle + error > max:
						temp = angle + error - 180
						ls_all.update({x:[min, temp]})
					elif angle < mid and angle - error < min:
						temp = angle - error + 180
						ls_all.update({x:[max, temp]})
				if (self.show_smooth_data and self.show_smooth_error_bars):
					if s_angle > mid and s_angle + s_error > max:
						temp = s_angle + s_error - 180
						ls_smooth.update({x:[min, temp]})
					elif s_angle < mid and s_angle - s_error < min:
						temp = s_angle - s_error + 180
						ls_smooth.update({x:[max, temp]})

			self.ax3.set_ylim(min, max)
			if (self.show_error_bars and self.show_raw_data):
				l_all_AoLP = self.ax3.errorbar(self.x_axis_list, bottle.all_AoLP * RtoD, yerr = bottle.std_AoLP * RtoD, fmt=".", ecolor=self.raw_error_bars_color, color = self.all_I0_color, linestyle = 'none', markersize=self.marker_size, label="AoLP", zorder=0, errorevery=self.error_step)
				for i, a in ls_all.items():
					self.ax3.vlines(i, a[0], a[1], colors = self.raw_error_bars_color, zorder=0) #"green", linewidth=5)#

			if self.show_smooth_data and self.show_smooth_error_bars:
				l_smooth_AoLP = self.ax3.errorbar(self.x_axis_list, bottle.smooth_AoLP * RtoD, yerr = bottle.std_smooth_AoLP * RtoD, fmt=".", color = self.smooth_I0_color, ecolor=self.smooth_error_bars_color, linestyle = 'none', markersize=self.marker_size, label="AoLP", zorder=1, errorevery=self.error_step)
				# plt.errorbar(list(ls.keys()), [-90, 90], yerr=list(ls.values()), fmt='C0 ')
				for i, a in ls_smooth.items():
					self.ax3.vlines(i, a[0], a[1], colors = self.smooth_I0_color, zorder=1) #"green", linewidth=5)#

		if (self.show_error_bars and self.show_raw_data):
			self.ax3_lines.append([l_all_AoLP, l_all_AoLP.get_label()])

		if  self.show_smooth_data and self.show_smooth_error_bars:
			self.ax3_lines.append([l_smooth_AoLP, l_smooth_AoLP.get_label()])

		if self.show_avg_AoLP:
			l_avg_AoLP, = self.ax3.plot(self.x_axis_list,[bottle.AoLP_average * RtoD] * len(self.x_axis_list), color = self.smooth_I0_color, label="Avg: " + str(bottle.AoLP_average * RtoD)[:4])
			self.ax3_lines.append([l_avg_AoLP, l_avg_AoLP.get_label()])


		if self.show_currents and self.eq_currents.valid:
			self.ax23 = self.ax2.twinx()
			l_AoJlos, = self.ax23.plot(self.eq_currents.GetNormTimes(self.divisor), self.eq_currents.data["AoJlos"] * RtoD, "*", color = self.currents_color, label="AoJlos")
			# l_AoJlos, = self.ax23.plot(self.eq_currents.GetNormTimes(self.divisor), np.sin(self.eq_currents.data["AoJlos"]) * self.eq_currents.data["J_norm"], color = self.currents_color, label="AoJlos")
			self.ax2_lines.append([l_AoJlos, l_AoJlos.get_label()])

			l_AoJapp, = self.ax3.plot(self.eq_currents.GetNormTimes(self.divisor), self.eq_currents.data["AoJapp"] * RtoD, "*", color = self.currents_color, label="AoJapp")
			self.ax3_lines.append([l_AoJapp, l_AoJapp.get_label()])


		# self.ax3.set_ylim(45, 60)


		if bottle.observation_type == "fixed":
			if bottle.AoRD is not False and self.show_AoRD:
				if bottle.graph_angle_shift == 1:
					AoRD = SetAngleBounds(bottle.AoRD, 0, np.pi)
				elif bottle.graph_angle_shift == 0:
					AoRD = SetAngleBounds(bottle.AoRD, -np.pi/2, np.pi/2)
				l_AoRD, = self.ax3.plot(self.x_axis_list,[bottle.AoRD * RtoD] * len(self.x_axis_list), linewidth=3, color=self.AoRD_color, label="AoRD: " + str(bottle.AoRD*RtoD)[:5], zorder=2)
				self.ax3_lines.append([l_AoRD, l_AoRD.get_label()])
				# l54, = self.ax3.plot(self.x_axis_list,[bottle.AoRD_ortho * RtoD] * len(self.x_axis_list), ":g", linewidth=self.marker_size, label="AoRD ortho: " + str(bottle.AoRD_ortho*RtoD)[:5])
			if (bottle.AoBapp and bottle.AoBlos) is not False and self.show_AoBapp:
				if bottle.graph_angle_shift == 1:
					AoBapp = SetAngleBounds(bottle.AoBapp, 0, np.pi)
				elif bottle.graph_angle_shift == 0:
					AoBapp = SetAngleBounds(bottle.AoBapp, -np.pi/2, np.pi/2)
				l_AoBapp, = self.ax3.plot(self.x_axis_list, [AoBapp * RtoD] * len(self.x_axis_list), linewidth=3, color = self.AoBapp_color, label="AoBapp: " + str(bottle.AoBapp*RtoD)[:5], zorder=2)
				self.ax3_lines.append([l_AoBapp, l_AoBapp.get_label()])
				# l_AoBapp_ortho, = self.ax3.plot(self.x_axis_list,[bottle.AoBapp_ortho * RtoD] * len(self.x_axis_list), ":b", linewidth=self.marker_size, label="AoBapp ortho: " + str(bottle.AoBapp_ortho*RtoD)[:5])
				# self.ax3_lines.append([l_AoBapp_ortho, l_AoBapp_ortho.get_label()])

				# l_AoBlos, = self.ax3.plot(self.x_axis_list,[bottle.AoBlos * RtoD] * len(self.x_axis_list), "orange", label="AoBlos: " + str(bottle.AoBlos*RtoD)[:5])
				# self.ax3_lines.append([l_AoBlos, l_AoBlos.get_label()])

		elif self.xaxis_azimut and bottle.observation_type == "fixed_elevation_discrete_rotation":
			# print("DEBUG plot discrete")
			if self.show_AoRD:
				l_AoRD, = self.ax3.plot(bottle.discrete_rotation_times, bottle.AoRD * RtoD, "k*", markersize=self.single_star_size*1.5, label="AoRD", zorder=2)
				l_AoRD, = self.ax3.plot(bottle.discrete_rotation_times, bottle.AoRD * RtoD, "*", color = self.AoRD_color, markersize=self.single_star_size, label="AoRD", zorder=3)
				self.ax3_lines.append([l_AoRD, l_AoRD.get_label()])
			if self.show_AoBapp:
				l_AoBapp, = self.ax3.plot(bottle.discrete_rotation_times, bottle.AoBapp * RtoD, "k*", markersize=self.single_star_size*1.5, label="AoBapp", zorder=2)
				l_AoBapp, = self.ax3.plot(bottle.discrete_rotation_times, bottle.AoBapp * RtoD, "*", color = self.AoBapp_color, markersize=self.single_star_size, label="AoBapp", zorder=3)
				self.ax3_lines.append([l_AoBapp, l_AoBapp.get_label()])

			rot = 0
			if len(self.x_axis_ticks_pos) > 15:
				rot = 60

			self.ax3.set_xticks(self.x_axis_ticks_pos)
			self.ax3.set_xticklabels(self.x_axis_ticks_label, rotation = rot)
			# self.ax3.xticks(rotation = 30)

		elif self.xaxis_azimut and bottle.observation_type == "fixed_azimut_discrete_rotation":
			# print("DEBUG plot discrete")
			if self.show_AoRD:
				l_AoRD, = self.ax3.plot(bottle.discrete_rotation_times, bottle.AoRD * RtoD, "*", color = self.AoRD_color, markersize=self.marker_size*4, label="AoRD", zorder=2)
				self.ax3_lines.append([l_AoRD, l_AoRD.get_label()])
			if self.show_AoBapp:
				l_AoBapp, = self.ax3.plot(bottle.discrete_rotation_times, bottle.AoBapp * RtoD, "*", color = self.AoBapp_color, markersize=self.marker_size*4, label="AoBapp", zorder=2)
				self.ax3_lines.append([l_AoBapp, l_AoBapp.get_label()])

		elif self.xaxis_azimut and bottle.observation_type == "fixed_elevation_continue_rotation":

			ax1_lim = self.ax1.get_ylim()
			self.ax1.set_ylim(ax1_lim)
			ax2_lim = self.ax2.get_ylim()
			ax3_lim = self.ax3.get_ylim()
			for i in range(bottle.nb_continue_rotation):
				# print("DEBUG ROT TIMES", i * 360 + (bottle.source_azimut * RtoD)%360)

				if not self.show_RS_model:
					self.ax1.plot([i * 360 + (bottle.source_azimut * RtoD)%360, i * 360 + (bottle.source_azimut * RtoD)%360], ax1_lim, "--k")
					self.ax2.plot([i * 360 + (bottle.source_azimut * RtoD)%360, i * 360 + (bottle.source_azimut * RtoD)%360], ax2_lim, "--k")
					self.ax3.plot([i * 360 + (bottle.source_azimut * RtoD)%360, i * 360 + (bottle.source_azimut * RtoD)%360], ax3_lim, "--k")
					if i != 0: # Draw red lines to delimit the rotations
						# print(i * 360)
						self.ax1.plot([i * 360, i * 360], ax1_lim, "r")
						self.ax2.plot([i * 360, i * 360], ax2_lim, "r")
						self.ax3.plot([i * 360, i * 360], ax3_lim, "r")

				if self.show_AoRD:
					l_AoRD, = self.ax3.plot(np.linspace(i * 360, (i+1) * 360, len(bottle.AoRD)), bottle.AoRD * RtoD, "*", color = self.AoRD_color, markersize=self.marker_size, label="AoRD", zorder=2)
					self.ax3_lines.append([l_AoRD, l_AoRD.get_label()])
				if self.show_AoBapp:
					l_AoBapp, = self.ax3.plot(np.linspace(i * 360, (i+1) * 360, len(bottle.AoBapp)), bottle.AoBapp * RtoD, "*", color = self.AoBapp_color, markersize=self.marker_size, label="AoBapp", zorder=2)
					self.ax3_lines.append([l_AoBapp, l_AoBapp.get_label()])


			self.ax3.set_xticks(self.x_axis_ticks_pos)
			self.ax3.set_xticklabels(self.x_axis_ticks_label)

			# x_axis, RD_diff, Bapps_diff = self.CompareAngles()
			# if ax3_lim[1] > 100:
			# 	print("DEBUG!!!!!!")
			# 	RD_diff, Bapps_diff  = UnifyAngles(RD_diff, 1)[0], UnifyAngles(Bapps_diff, 1)[0]
			# self.ax3.plot(x_axis, RD_diff * RtoD, "*", color = "xkcd:black", markersize=self.marker_size, label="AoRD")

			# self.ax3_lines.append([l_AoRD, l_AoRD.get_label()])

		# self.ax4.plot(self.x_axis_list, bottle.all_TempPM, "k.", linestyle = 'none', markersize=self.marker_size, label="PM")
		# self.ax4.plot(self.x_axis_list, bottle.all_TempOptical, "r.", linestyle = 'none', markersize=self.marker_size, label="Optical")
		# self.ax4.plot(self.x_axis_list, bottle.all_TempAmbiant, "b.", linestyle = 'none', markersize=self.marker_size, label="Ambiant")
		# self.ax2.plot(self.x_axis_list,[DoLP_average] * nb_rot, "b", label="Avg: " + str(DoLP_average))


		if self.show_RS_model:
			# self.ax11 = self.ax1.twinx()
			# self.ax22 = self.ax2.twinx()

			tmp_xshift = 0

			maxDoLP = np.max(bottle.smooth_DoLP)
			for im, models in enumerate(self.J2RAYS1_models):
				# print(self.J2RAYS1_model.x_axis)
				# self.ax11.yaxis.set_visible(False)
				self.ax1.plot(models.x_axis + tmp_xshift, models["I0"], "*", color = self.model_colors[im], marker = self.model_symbols[im])
				self.ax2.plot(models.x_axis + tmp_xshift, models["DoLP"], "*", color = self.model_colors[im], marker = self.model_symbols[im])
				# self.ax2.set_ylim((0, max(np.max(bottle.smooth_DoLP), np.max(models.data["DoLP"]))))
				if np.max(models.data["DoLP"]) > maxDoLP:
					maxDoLP = np.max(models.data["DoLP"])
					self.ax2.set_ylim((0, 1.1 * maxDoLP))

				self.ax3.plot(models.x_axis + tmp_xshift, RtoD * SetAngleBounds(DtoR * models.data["AoRD"], -np.pi/2 + np.pi/2 * bottle.graph_angle_shift, np.pi/2 + np.pi/2 * bottle.graph_angle_shift), "*", color = self.model_colors[im], marker = self.model_symbols[im])

			if self.adapt_flux_scale:
				border = 1.1
				max_data, min_data = np.max(bottle.smooth_I0), np.min(bottle.smooth_I0)
				ratio_data = max_data / min_data
				print(f"Min data, Max Data, Ratio data: {min_data}, {max_data}, {ratio_data}")
				nt, nb = [max_data], [min_data]
				for models in self.J2RAYS1_models:
					max_model, min_model = np.max(models.data["I0"]) , np.min(models.data["I0"])
					ratio_model = max_model / min_model
					shift = (ratio_data * min_model - max_model) / (1 - ratio_data)
					# shift = (ratio_model * min_data - max_data) / (1 - ratio_model)
					print(f"Min model, Max model, Ratio model, Shift to match data: {min_model}, {max_model}, {ratio_model}, {shift}")
					nt.append(max_model)
					nb.append(min_model)

				self.ax1.set_ylim(np.min(nb) / border, np.max(nt) * border)

					# nt, nb = [np.max(bottle.smooth_DoLP)], [np.min(bottle.smooth_DoLP)]
					# for models in self.J2RAYS1_models:
					# 	nt.append(np.max(models.data["DoLP"]))
					# 	nb.append(np.min(models.data["DoLP"]))
					# self.ax2.set_ylim(np.min(nb) / border, np.max(nt) * border)

				if bottle.observation_type == "fixed_elevation_continue_rotation":
					ax1_lim = self.ax1.get_ylim()
					self.ax1.set_ylim(ax1_lim)
					ax2_lim = self.ax2.get_ylim()
					ax3_lim = self.ax3.get_ylim()
					for i in range(bottle.nb_continue_rotation):
						# print("DEBUG ROT TIMES", i * 360 + (bottle.source_azimut * RtoD)%360)

						self.ax1.plot([i * 360 + (bottle.source_azimut * RtoD)%360, i * 360 + (bottle.source_azimut * RtoD)%360], ax1_lim, "--k")
						self.ax2.plot([i * 360 + (bottle.source_azimut * RtoD)%360, i * 360 + (bottle.source_azimut * RtoD)%360], ax2_lim, "--k")
						self.ax3.plot([i * 360 + (bottle.source_azimut * RtoD)%360, i * 360 + (bottle.source_azimut * RtoD)%360], ax3_lim, "--k")
						if i != 0: # Draw red lines to delimit the rotations
							# print(i * 360)
							self.ax1.plot([i * 360, i * 360], ax1_lim, "r")
							self.ax2.plot([i * 360, i * 360], ax2_lim, "r")
							self.ax3.plot([i * 360, i * 360], ax3_lim, "r")

					# if r1 >= r11:
					# 	# self.ax1.set_ylim(0, top*border)
					# 	# self.ax11.set_ylim(0, nb * r1*border)
					# 	self.ax1.set_ylim(bot/border, top*border)
					# 	self.ax11.set_ylim(nt / r1 /border, nt*border)
					# 	# self.ax11.set_ylim(nb/border, nb * r1*border)
					# 	print(f"NEW RATIO 2 {r1}, from {nb} to {nb * r1}")
					# else:
					# 	# self.ax1.set_ylim(0, top * 1.1)
					# 	# self.ax11.set_ylim(0, nt * 1.1)
					# 	# self.ax1.set_ylim(0, bot * r11*border)
					# 	# self.ax11.set_ylim(0, nt*border)
					# 	self.ax1.set_ylim(nt/r11/border, nt *border)
					# 	# self.ax1.set_ylim(bot/border, bot * r11*border)
					# 	self.ax11.set_ylim(nb/border, nt*border)
					# 	print(f"NEW RATIO 1 {r11}, from {bot} to {bot * r11}")

					# self.axs[11].xlim((0, ))

		self.f1.subplots_adjust(hspace=0)

		###Set title
		# self.f1.suptitle(bottle.saving_name.replace("_", " "))

		# plt.setp([a.get_xticklabels() for a in f1.axes[:-1]], visible=False)
		if bottle.observation_type == "fixed":
			plt.minorticks_on()

		if self.show_time:
			self.ax12 = self.ax1.twiny()
			self.ax12.set_xlim(self.ax1.get_xlim())
			xticks = plt.xticks()[0] * self.divisor
			xticks = [bottle.DateTime("start", format=self.time_format) + time.timedelta(seconds = t) for t in xticks]
			# xticks = [bottle.DateTime() + time.timedelta(seconds = t) + bottle.head_jump for t in xticks]
			self.ax12.set_xticklabels([st.strftime("%H:%M") for st in xticks])
			self.ax12.set_xlabel(self.time_label)

		if self.show_legend:
			self.ax1.legend(list(zip(*self.ax1_lines))[0], list(zip(*self.ax1_lines))[1])
			self.ax2.legend(list(zip(*self.ax2_lines))[0], list(zip(*self.ax2_lines))[1])
			self.ax3.legend(list(zip(*self.ax3_lines))[0], list(zip(*self.ax3_lines))[1], loc = "lower center")


		print("Saving graphs in", bottle.data_file_name + "/" + bottle.saving_name + '_graphs.png')
		if bottle.instrument_name in ["carmen", "corbel", "gdcu"]:
			plt.savefig(bottle.data_file_name + "/" + bottle.saving_name + '_graphs.png', bbox_inches='tight')
			# plt.savefig(bottle.data_file_name + "/" + bottle.saving_name + '_graphs.eps', bbox_inches='tight')
		else:
			plt.savefig("/".join(bottle.data_file_name.split("/")[:-1]) + "/" + bottle.saving_name + '_graphs.png', bbox_inches='tight')
			# plt.savefig("/".join(bottle.data_file_name.split("/")[:-1]) + "/" + bottle.saving_name + '_graphs.eps', bbox_inches='tight')


	def CompareRayleigh(self, bottle, model):

		data_I = np.interp(model.x_axis, self.x_axis_list, bottle.smooth_I0)
		data_D = np.interp(model.x_axis, self.x_axis_list, bottle.smooth_DoLP)
		data_A = np.interp(model.x_axis, self.x_axis_list, bottle.smooth_AoLP)

		ratio = lambda x, y: x / y
		diff = lambda x, y: x - y

		fig, axs = plt.subplots(3, figsize=(50, 20))

		axs[0].plot(model.x_axis, ratio(data_I, model.flux), "*")
		axs[1].plot(model.x_axis, ratio(data_D, model.DoLP), "*")
		axs[2].plot(model.x_axis, GetAnglesDifference(data_A, model.AoLP)*RtoD, "*")

		axs[0].set_ylabel("Ratio")
		axs[1].set_ylabel("Ratio")
		axs[2].set_ylabel("Diff")
		fig.suptitle("Data /- Model")

		if bottle.instrument_name in ["carmen", "corbel", "gdcu"]:
			plt.savefig(bottle.data_file_name + "/" + bottle.saving_name + '_RSmodel_comp.png', bbox_inches='tight')
		else:
			plt.savefig("/".join(bottle.data_file_name.split("/")[:-1]) + "/" + bottle.saving_name + '_RSmodel_comp.png', bbox_inches='tight')


	def MakeSNratio(self, bottle):
		if not self.show_SN:
			return 0

		f5, ax = plt.subplots(1, figsize=(10, 20))

		ax.plot(self.x_axis_list, bottle.all_SN, "--", color = self.all_SN_color, label="all SN")
		ax.set_ylabel("Raw")
		ax1 = plt.twinx(ax)
		ax1.plot(self.x_axis_list, bottle.smooth_SN, "--", color = self.smooth_SN_color, label="smooth SN")
		ax1.set_ylabel("Smooth")

		ax.set_title("SN ratio")
		print("Saving graphs in", bottle.data_file_name + "/" + bottle.saving_name + '_SNratio.png')
		if bottle.instrument_name in ["carmen", "corbel", "gdcu"]:
			plt.savefig(bottle.data_file_name + "/" + bottle.saving_name + '_SNratio.png', bbox_inches='tight')
			# plt.savefig(bottle.data_file_name + "/" + bottle.saving_name + '_graphs.eps', bbox_inches='tight')
		else:
			plt.savefig("/".join(bottle.data_file_name.split("/")[:-1]) + "/" + bottle.saving_name + '_SNratio.png', bbox_inches='tight')
			# plt.savefig("/".join(bottle.data_file_name.split("/")[:-1]) + "/" + bottle.saving_name + '_graphs.eps', bbox_inches='tight')

	def MakeI0SNratio(self, bottle):
		f6, ax = plt.subplots(1, figsize=(10, 20))

		print(bottle.all_I0,  bottle.std_I0, bottle.all_I0 / bottle.std_I0)
		print(bottle.smooth_I0, bottle.std_smooth_I0, bottle.smooth_I0 / bottle.std_smooth_I0)

		ax.plot(self.x_axis_list, bottle.all_I0 / bottle.std_I0, "--", color = self.all_SN_color, label="I0 SN")
		ax.set_ylabel("Raw")
		ax1 = plt.twinx(ax)
		ax1.plot(self.x_axis_list, bottle.smooth_I0 / bottle.std_smooth_I0, "--", color = self.smooth_SN_color, label="smooth I0 SN")
		ax1.set_ylabel("Smooth")

		ax.set_title("Intensity SN ratio")
		print("Saving graphs in", bottle.data_file_name + "/" + bottle.saving_name + '_I0SNratio.png')
		if bottle.instrument_name in ["carmen", "corbel", "gdcu"]:
			plt.savefig(bottle.data_file_name + "/" + bottle.saving_name + '_I0SNratio.png', bbox_inches='tight')
			# plt.savefig(bottle.data_file_name + "/" + bottle.saving_name + '_graphs.eps', bbox_inches='tight')
		else:
			plt.savefig("/".join(bottle.data_file_name.split("/")[:-1]) + "/" + bottle.saving_name + '_I0SNratio.png', bbox_inches='tight')
			# plt.savefig("/".join(bottle.data_file_name.split("/")[:-1]) + "/" + bottle.saving_name + '_graphs.eps', bbox_inches='tight')


	def DrawEiscatParam(self, axe, param, alt):
		label = param.replace("_", "\_")

		t, d, e = self.eiscat_data.GetParameter(param, alt, time_divisor = self.divisor)
		l_eiscat_data, capline, barlinecol = axe.errorbar(t, d, yerr = e, fmt = ".", label = label)
		avg = np.average(d)
		axe.plot(t, [avg] * len(t), label = "Average: " + str(avg)[:4])

	def MakeEiscatDataPlot(self, bottle):
		f4, ax = plt.subplots(4, figsize=(10, 20), sharex = True)

		altitude = bottle.filters[0]
		if altitude == "r": altitude = 220
		elif altitude == "v": altitude = 110
		elif altitude == "b" or "m": altitude = 85

		self.DrawEiscatParam(ax[0], "N_e", altitude)
		self.DrawEiscatParam(ax[1], "T_i", altitude)
		self.DrawEiscatParam(ax[2], "T_e", altitude)
		# self.DrawEiscatParam(ax[3], "nu_in", altitude)
		self.DrawEiscatParam(ax[3], "v_i", altitude)
		ax[-1].set_xlabel(self.xlabel)

		for a in ax:
			a.legend()

		print("Saving graphs in", bottle.data_file_name + "/" + bottle.saving_name + '_eiscat.png')
		if bottle.instrument_name in ["carmen", "corbel", "gdcu"]:
			plt.savefig(bottle.data_file_name + "/" + bottle.saving_name + '_eiscat.png', bbox_inches='tight')
		else:
			plt.savefig("/".join(bottle.data_file_name.split("/")[:-1]) + "/" + bottle.saving_name + '_eiscat.png', bbox_inches='tight')


	def MakeMagDataPlots(self, bottle):
		f3, ax = plt.subplots(3, figsize=(10, 20), sharex = True)

		ax[0].plot(self.mag_data.times_sec, self.mag_data.Total, "k", label="Total")
		# ax[1].plot(self.mag_data.times_sec, self.mag_data.Dec, "r", label="Dec")
		ax[1].plot(self.mag_data.times_sec, self.mag_data.Horiz, "b", label="Horiz")
		# ax[3].plot(self.mag_data.times_sec, self.mag_data.Vert, "g", label="Vert")
		ax[2].plot(self.mag_data.times_sec, self.mag_data.Incl, "y", label="Incl")
		ax[-1].set_xlabel(self.xlabel)
		for a in ax:
			a.legend()
		print("Saving graphs in", bottle.data_file_name + "/" + bottle.saving_name + '_magneto.png')
		if bottle.instrument_name in ["carmen", "corbel", "gdcu"]:
			plt.savefig(bottle.data_file_name + "/" + bottle.saving_name + '_magneto.png', bbox_inches='tight')
		else:
			plt.savefig("/".join(bottle.data_file_name.split("/")[:-1]) + "/" + bottle.saving_name + '_magneto.png', bbox_inches='tight')



	def CompareAngles(self, bottle):
		if bottle.observation_type == "fixed_elevation_continue_rotation":
			AoRDs = bottle.AoRD
			AoBapps = bottle.AoBapp
			for i in range(bottle.nb_continue_rotation - 1):
				AoRDs = np.append(AoRDs, bottle.AoRD)
				AoBapps = np.append(AoBapps, bottle.AoBapp)

			nb_AoLP = len(bottle.smooth_AoLP)
			nb_AoRD = len(AoRDs)
			nb_AoBapp = len(AoBapps)
			if nb_AoLP > nb_AoRD:
				x_axis = ReduceList(self.x_axis_list, nb_AoRD)
				AoLPs = ReduceList(bottle.smooth_AoLP, nb_AoRD)
			elif nb_AoLP < nb_AoRD:
				x_axis = self.x_axis_list
				AoLPs = bottle.smooth_AoLP
				AoRDs = ReduceList(AoRDs, nb_AoLP)
				AoBapps = ReduceList(AoBapps, nb_AoLP)
			else:
				x_axis = self.x_axis_list
				AoLPs = bottle.smooth_AoLP

			RD_diff = GetAnglesDifference(AoLPs, AoRDs)
			Bapps_diff = GetAnglesDifference(AoLPs, AoBapps)

			return x_axis, RD_diff, Bapps_diff


	def CompareAnglesPlots(self, bottle):
		f4, ax = plt.subplots(1, figsize=(10, 20), sharex = True)

		x_axis, RD_diff, Bapps_diff = self.CompareAngles(bottle)

		ax.plot(x_axis, RD_diff * RtoD, "*", color = self.AoRD_color, label="Total")
		ax.plot(x_axis, Bapps_diff * RtoD, "*", color = self.AoBapp_color, label="Total")

		print("Saving graphs in", bottle.data_file_name + "/" + bottle.saving_name + '_angles_comp.png')
		if bottle.instrument_name in ["carmen", "corbel", "gdcu"]:
			plt.savefig(bottle.data_file_name + "/" + bottle.saving_name + '_angles_comp.png', bbox_inches='tight')
		else:
			plt.savefig("/".join(bottle.data_file_name.split("/")[:-1]) + "/" + bottle.saving_name + '_angles_comp.png', bbox_inches='tight')

	def MakeCorrelationPlots(self, bottle):
		f2, (ax1, ax2) = plt.subplots(2, figsize=(10, 20))

		ax1.plot(bottle.all_I0, bottle.all_DoLP, "k+", label="Raw")
		ax1.plot(bottle.smooth_I0, bottle.smooth_DoLP, "r+", label="Smoothed (" + str(bottle.smoothing_factor) + " " + str(bottle.smoothing_unit)[:3] + ")")
		ax1.set_xlabel("Intensity (mV)")
		ax1.set_ylabel("DoLP (%)")
		ax1.legend()

		if not bottle.NoVref and bottle.instrument_name=="carmen":
			ax2.plot(bottle.all_I0, bottle.all_Iref, "k+", label="Raw")
			ax2.plot(bottle.smooth_I0, bottle.smooth_Iref, "r+", label="Smoothed (" + str(bottle.smoothing_factor) + " " + str(bottle.smoothing_unit)[:3] + ")")
			ax2.set_xlabel("Ref Intensity (mV)")
			ax2.set_ylabel("Pola Intensity (mV)")

			ax2.legend()

		print("Saving correlation in", bottle.data_file_name + "/" + bottle.saving_name + '_correlations.png')
		if bottle.instrument_name in ["carmen", "corbel", "gdcu"]:
			plt.savefig(bottle.data_file_name + "/" + bottle.saving_name + '_correlation.png', bbox_inches='tight')
		else:
			plt.savefig("/".join(bottle.data_file_name.split("/")[:-1]) + "/" + bottle.saving_name + '_correlation.png', bbox_inches='tight')


	def GetCoherence(self, data1, data2, freq):
		NN = 2**9
		f, P12 = signal.csd(data1, data2, freq, detrend='linear',nperseg=NN)
		f, P11 = signal.csd(data1, data1, freq, detrend='linear',nperseg=NN)
		f, P22 = signal.csd(data2, data2, freq, detrend='linear',nperseg=NN)
		Nf = f.shape[0]

		C = np.abs(P12)**2/(P11*P22)

		ss = np.imag(P12[:])/np.abs(P12)
		cc = np.real(P12[:])/np.abs(P12)
		angle_12 = np.arctan2(ss,cc)

		return f, C, angle_12

	def MakeCoherencePLots(self, bottle):
		f_coherence, axs = plt.subplots(2, figsize=(10, 20))

		ax_coherence = axs[0]
		ax_angle = axs[1]

		f, C, A = self.GetCoherence(bottle.smooth_I0, bottle.smooth_DoLP, bottle.GetRotationFrequency())
		ax_coherence.semilogx( (1./f), C, "k", label="Intensity-DoLP")
		ax_angle.semilogx( (1./f), A * RtoD, "k", label="Intensity-DoLP")

		f, C, A= self.GetCoherence(bottle.smooth_I0, bottle.smooth_AoLP, bottle.GetRotationFrequency())
		ax_coherence.semilogx( (1./f), C, "r", label="Intensity-AoLP")
		ax_angle.semilogx( (1./f), A * RtoD, "r", label="Intensity-AoLP")

		f, C, A = self.GetCoherence(bottle.smooth_DoLP, bottle.smooth_AoLP, bottle.GetRotationFrequency())
		ax_coherence.semilogx( (1./f), C, "b", label="DoLP-AoLP")
		ax_angle.semilogx( (1./f), A * RtoD, "b", label="DoLP-AoLP")

		ax_coherence.set_xlabel("Period (sec)")
		ax_coherence.set_ylabel("Coherence")
		# ax_coherence.set_xlim((-1,1))
		ax_coherence.set_ylim((0,1))
		ax_coherence.legend()

		ax_angle.set_xlabel("Period (sec)")
		ax_angle.set_ylabel("Shift angle")
		# ax_angle.set_xlim((-1,1))
		ax_angle.set_ylim((-180,180))
		ax_angle.legend()


		print("Saving coherence in", bottle.data_file_name + "/" + bottle.saving_name + '_coherence.png')
		if bottle.instrument_name in ["carmen", "corbel", "gdcu"]:
			plt.savefig(bottle.data_file_name + "/" + bottle.saving_name + '_coherence.png', bbox_inches='tight')
		else:
			plt.savefig("/".join(bottle.data_file_name.split("/")[:-1]) + "/" + bottle.saving_name + 'coherence.png', bbox_inches='tight')


	def MakeFFT(self, bottle):
		self.f1, (self.ax1, self.ax2, self.ax3) = plt.subplots(3, sharex=True, figsize=(16, 8))

		self.ax1.set_ylabel("Intensity FFT")
		self.ax2.set_ylabel("DoLP FFT")
		self.ax3.set_ylabel("AoLP FFT")

		xaxis = np.fft.fftfreq(self.x_axis_list.shape[-1], d = 1000/20.)

		FFT_I0 = np.fft.fft(bottle.all_I0)
		FFT_DoLP = np.fft.fft(bottle.all_DoLP)
		FFT_AoLP = np.fft.fft(bottle.all_AoLP)

		smooth_FFT_I0 = np.fft.fft(bottle.smooth_I0)
		smooth_FFT_DoLP = np.fft.fft(bottle.smooth_DoLP)
		smooth_FFT_AoLP = np.fft.fft(bottle.smooth_AoLP)

		self.ax1.plot(xaxis, FFT_I0, xaxis, smooth_FFT_I0)
		self.ax2.plot(xaxis, FFT_DoLP, xaxis, smooth_FFT_DoLP)
		self.ax3.plot(xaxis, FFT_AoLP, xaxis, smooth_FFT_AoLP)

	def CompareBottles(self, b1, b2):
		f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, figsize=(16, 8))

		if b1.nb_rot > b2.nb_rot:
			all_times = list(map(lambda x: x.total_seconds(), b2.all_times))
			# print("DEBUG all_times", len(all_times), all_times[0])
			I0_1, DoLP_1, AoLP_1 = b1.GetInterpolation(all_times)
			I0_2, DoLP_2, AoLP_2 = b2.smooth_I0, b2.smooth_DoLP, b2.smooth_AoLP
		else:
			all_times = list(map(lambda x: x.total_seconds(), b1.all_times))
			# print("DEBUG all_times", len(all_times), all_times[0])
			I0_2, DoLP_2, AoLP_2 = b2.GetInterpolation(all_times)
			I0_1, DoLP_1, AoLP_1 = b1.smooth_I0, b1.smooth_DoLP, b1.smooth_AoLP

		delta = lambda a, b: (b-a) / a * 100
		ax1.plot(self.x_axis_list, delta(I0_1, I0_2))
		ax2.plot(self.x_axis_list, delta(DoLP_1, DoLP_2))
		# ax3.plot(all_times, AoLP_1*RtoD)
		# ax3.plot(all_times, AoLP_2*RtoD)
		# ax3.plot(self.x_axis_list, delta(AoLP_1, AoLP_2))
		ax3.plot(self.x_axis_list, (AoLP_1 - AoLP_2)*RtoD)

		ax1.set_ylabel("I0")
		ax1.grid(which="both")

		ax2.set_ylabel("DoLP")
		ax2.grid(which="both")

		ax3.set_ylabel("AoLP")
		ax3.set_xlabel(self.xlabel)
		ax3.grid(which="both")

		print("Saving bottle comparison in", b1.data_file_name + "/" + b1.saving_name + '_bottle_comparison.png')
		if b1.instrument_name in ["carmen", "corbel", "gdcu"]:
			plt.savefig(b1.data_file_name + "/" + b1.saving_name + '_bottle_comparison.png', bbox_inches='tight')
		else:
			plt.savefig("/".join(b1.data_file_name.split("/")[:-1]) + "/" + b1.saving_name + 'bottle_comparison.png', bbox_inches='tight')



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
