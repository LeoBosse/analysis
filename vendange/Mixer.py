#!/usr/bin/python3
# -*-coding:utf-8 -*

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage{lmodern}']
matplotlib.rcParams['font.size'] = 16
matplotlib.rcParams['svg.fonttype'] = 'none'
import datetime as dt


from utils import *
from rotation import *
from bottle import *
from AllSkyData import *
from eiscat_data import *
from scipy import signal
import sys
import os
from subprocess import call



class Mixer:
	def __init__(self, bottles, mag_data=False, comp_bottle=False, mode = "default"):

		self.bottles = bottles
		self.comp_bottle = comp_bottle
		self.mag_data = mag_data

		for bottle in self.bottles:
			if self.mag_data is not False:
				self.mag_data.StripTime(bottle.DateTime("start"), bottle.DateTime("end"))

			self.allsky_data = AllSkyData(bottle)
			self.allsky_data_available = bool(self.allsky_data.all_datetimes)

			self.eiscat_data = Eiscat(bottle)
			print("DEBUG:", self.eiscat_data.valid )


			self.mode = mode

			self.SetGraphParameter(bottle)
			self.MakeXAxis(bottle)
			self.MakeFigure()
			self.MakePlots(bottle)

			self.MakeFFT(bottle)

			if self.mag_data is not False:
				self.MakeMagDataPlots(bottle)
			self.MakeCorrelationPlots(bottle)

			if self.eiscat_data.valid:
				self.MakeEiscatDataPlot(bottle)

			if bottle.observation_type == "fixed_elevation_continue_rotation":
				self.CompareAnglesPlots(bottle)

				if len(bottle.continue_rotation_times) > 2:
					print("AVERAGE ROTATIONS")
					bottle = bottle.GetAverageContinueRotations()
					self.SetGraphParameter(bottle)
					self.MakeXAxis(bottle)
					self.MakeFigure()
					self.MakePlots(bottle)
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
			print(len(self.x_axis_list), len(bottle.all_times))
			self.xlabel = "Azimut"

			if bottle.nb_continue_rotation <= 2:
				self.x_axis_ticks_pos = np.arange(0, 360 * bottle.nb_continue_rotation, 45)
				self.x_axis_ticks_label = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"] *  bottle.nb_continue_rotation
			else:
				self.x_axis_ticks_pos = np.arange(0, 360 * bottle.nb_continue_rotation, 90)
				self.x_axis_ticks_label = ["N", "E", "S", "W"] *  bottle.nb_continue_rotation

		else:
			if bottle.all_times[-1] > time.timedelta(hours=2):
				self.divisor = 3600.
				self.xlabel = "Time (hours)"
			elif bottle.all_times[-1] > time.timedelta(minutes=2):
				self.divisor = 60.
				self.xlabel = "Time (minutes)"
			else:
				self.xlabel = "Time (seconds)"

			# self.x_axis_list 	   = np.array([t.total_seconds() for t in bottle.all_times_since_start]) / self.divisor
			norm =  bottle.all_times[0].total_seconds()
			self.x_axis_list = np.array([t.total_seconds() - norm for t in bottle.all_times]) / self.divisor

			if self.xaxis_azimut and bottle.observation_type == "fixed_elevation_discrete_rotation":
				self.xlabel = "Azimut (degrees)"
				self.x_axis_ticks_pos = bottle.discrete_rotation_times
				print(bottle.discrete_rotation_azimuts * RtoD)
				self.x_axis_ticks_label = np.round(bottle.discrete_rotation_azimuts * RtoD % 360)
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

			self.ax1.set_ylabel("Intensity (mV)")
			self.ax2.set_ylabel("DoLP (\%)")
			self.ax3.set_ylabel("AoLP (degrees)")
			# ax4.set_ylabel("Temperature (Celsius)")

			# ax1.set_xlabel(self.xlabel)
			# ax2.set_xlabel(self.xlabel)
			self.ax3.set_xlabel(self.xlabel)

	def SetGraphParameter(self, bottle):
		self.marker_size = 3
		self.single_star_size = self.marker_size*5

		self.show_error_bars 		= True
		self.show_smooth_error_bars = True
		self.max_error_bars = 10000 #If there are too many data, takes very long to plot error bars. If != 0, will plot max_error_bars error bars once every x points.

		#If True Will show the azimut on the xaxis instead of time for rotations
		self.xaxis_azimut = True

		self.all_I0_color = "xkcd:black"
		self.smooth_I0_color = "xkcd:red"

		self.smooth_ref_color = "xkcd:green"
		self.AoBapp_color = "xkcd:turquoise"
		self.AoRD_color = "xkcd:hot pink"
		self.EISCAT_color = "xkcd:mustard"
		self.mag_color = "xkcd:mustard"

		### xkcd color guide: https://xkcd.com/color/rgb/
		print("MIXER FILTER:", bottle.filters)
		if bottle.filters:
			self.pola_color = bottle.filters[0]
			if 	 self.pola_color == "r": self.smooth_I0_color = "xkcd:red"
			elif self.pola_color == "v": self.smooth_I0_color = "xkcd:green"
			elif self.pola_color == "b": self.smooth_I0_color = "xkcd:blue"
			elif self.pola_color == "m": self.smooth_I0_color = "xkcd:purple"
			elif self.pola_color == "o": self.smooth_I0_color = "xkcd:orange"
			elif self.pola_color == "X": self.smooth_I0_color = "xkcd:orange"
			elif self.pola_color == "Y": self.smooth_I0_color = "xkcd:turquoise"
			else: self.smooth_I0_color = "red"

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



	def MakePlots(self, bottle):
		#Plotting the mean intensity I0, DoLP, AoLP for each rotation and more!

		# if self.show_error_bars:
		# 	smooth_yerr = bottle.var_smooth_I0
		# 	yerr = bottle.var_I0
		# 	yerr = None
		# else:
		# 	smooth_yerr = None

		print("START PLOTTING")

		if not self.show_error_bars:
			l_all_I0, = self.ax1.plot(self.x_axis_list, bottle.all_I0, ".", color = self.all_I0_color, linestyle = 'none', markersize=self.marker_size, label="Intensity", zorder=0)
		else:
			l_all_I0 = self.ax1.errorbar(self.x_axis_list, bottle.all_I0, yerr = bottle.std_I0, fmt=".", ecolor="grey", color = self.all_I0_color , linestyle = 'none', markersize=self.marker_size, label="Intensity", zorder=0, errorevery=self.error_step)

		if not self.show_smooth_error_bars:
			l_smooth_I0, = self.ax1.plot(self.x_axis_list, bottle.smooth_I0, ".", color = self.smooth_I0_color, linestyle = 'none', markersize=self.marker_size, label="Smooth Intensity (" + str(bottle.smoothing_factor) + " " + str(bottle.smoothing_unit)[:3] + ")", zorder=1)
		else:
			l_smooth_I0 = self.ax1.errorbar(self.x_axis_list, bottle.smooth_I0, yerr = bottle.std_smooth_I0, fmt=".", color = self.smooth_I0_color, linestyle = 'none', markersize=self.marker_size, label="Smooth Intensity (" + str(bottle.smoothing_factor) + " " + str(bottle.smoothing_unit)[:3] + ")", zorder=1, errorevery=self.error_step)


		self.ax1_lines.append([l_all_I0, l_all_I0.get_label()])
		self.ax1_lines.append([l_smooth_I0, l_smooth_I0.get_label()])
		# l_avg_I0, = self.ax1.plot(self.x_axis_list,[bottle.I0_average] * len(self.x_axis_list), color = self.smooth_I0_color, label="Avg: " + str(bottle.I0_average)[:4])
		# self.ax1_lines.append([l_avg_I0, l_avg_I0.get_label()])


		# ### Graph the I0 * DoLP line on the first subplot
		# self.ax14 = self.ax1.twinx()
		# # l_all_IDoLP, = self.ax14.plot(self.x_axis_list, bottle.all_I0 * bottle.all_DoLP, "k.", linestyle = 'none', markersize=self.marker_size, label="All Intensity * DoLP", zorder=2)
		# smooth_IDoLP, = self.ax14.plot(self.x_axis_list, bottle.smooth_I0 * bottle.smooth_DoLP / 100. * bottle.smooth_AoLP
		#
		# cd, ".", color = "xkcd:hot pink", linestyle = 'none', markersize=self.marker_size, label="Smooth Intensity * DoLP", zorder=2)
		#
		# self.ax14.set_ylabel("Ref intensity")



		if not bottle.NoVref and bottle.instrument_name == "carmen":
			self.ax13 = self.ax1.twinx()
			# ax13.set_xlim(ax1.get_xlim())
			# xticks = plt.xticks()[0] * self.divisor
			# xticks = [t + bottle.time + bottle.head_jump for t in xticks]
			# ax12.set_xticklabels([time.strftime("%H:%M:%S", time.localtime(st)) for st in xticks])

			# l_Iref, = self.ax13.plot(self.x_axis_list[1:], bottle.all_Iref[1:], ".", color = "black", linestyle = 'none', markersize=self.marker_size, label="Ref Intensity", zorder=2)
			# self.ax1_lines.append([l_Iref, l_Iref.get_label()])

			l_smooth_Iref, = self.ax13.plot(self.x_axis_list, bottle.smooth_Iref, ".", color = self.smooth_ref_color, linestyle = 'none', markersize=self.marker_size, label="Smooth Ref Intensity (" + str(bottle.smoothing_factor) + " " + str(bottle.smoothing_unit)[:3] + ")", zorder=2)
			self.ax1_lines.append([l_smooth_Iref, l_smooth_Iref.get_label()])

			# l_avg_Iref, = self.ax13.plot(self.x_axis_list, [bottle.Iref_average] * len(self.x_axis_list), color = self.smooth_ref_color, label="Avg Ref Intensity " + str(bottle.Iref_average)[:4], zorder=2)
			# self.ax1_lines.append([l_avg_Iref, l_avg_Iref.get_label()])
			self.ax13.set_ylabel("Ref intensity")



		# if self.allsky_data_available:
		# 	self.ax14 = self.ax1.twinx()
		# 	offset = 10
		# 	# new_fixed_axis = self.ax14.get_grid_helper().new_fixed_axis
		# 	# self.ax14.axis["right"] = new_fixed_axis(loc="right", axes=par2, offset=(offset, 0))
		# 	#
		# 	# self.ax14.axis["right"].toggle(all=True)
		#
		# 	l_ASI, = self.ax14.plot(self.allsky_data.GetNormTimes(bottle.DateTime(), self.divisor), self.allsky_data.brightness, "*", color = "orange", linestyle = 'none', markersize=self.marker_size, label="AllSKy Imager", zorder=2)
		# 	self.ax1_lines.append([l_ASI, l_ASI.get_label()])
		# 	print(self.ax14.get_yticklabels())
		# 	# self.ax14.set_yticklabels([l.get_text() for l in self.ax14.get_yticklabels()], horizontalalignment = "left")
		# 	# self.ax14.tick_params(direction='in', labelright=True, pad=-5)
		# 	self.ax14.set_ylabel("ASI Brightness")



		# if self.show_error_bars:
		# 	smooth_yerr = bottle.var_smooth_DoLP
		# 	yerr = bottle.var_DoLP
		# 	yerr = None
		# else:
		# 	smooth_yerr = None
		if not self.show_error_bars:
			l_all_DoLP, = self.ax2.plot(self.x_axis_list, bottle.all_DoLP, ".", color = self.all_I0_color, linestyle = 'none', markersize=self.marker_size, label="DoLP", zorder=0)
			self.ax2_lines.append([l_all_DoLP, l_all_DoLP.get_label()])
		else:
			(l_all_DoLP, _, _) = self.ax2.errorbar(self.x_axis_list, bottle.all_DoLP, yerr = bottle.std_DoLP, fmt=".", ecolor="grey", color = self.all_I0_color, linestyle = 'none', markersize=self.marker_size, label="DoLP", zorder=0, errorevery=self.error_step)

		if not self.show_smooth_error_bars:
			l_smooth_DoLP, = self.ax2.plot(self.x_axis_list, bottle.smooth_DoLP, ".", color = self.smooth_I0_color, linestyle = 'none', markersize=self.marker_size, label="Smooth DoLP (" + str(bottle.smoothing_factor) + " " + str(bottle.smoothing_unit)[:3] + ")", zorder=1)
		else:
			self.ax2_lines.append([l_all_DoLP, l_all_DoLP.get_label()])
			(l_smooth_DoLP, _, _) = self.ax2.errorbar(self.x_axis_list, bottle.smooth_DoLP, yerr = bottle.std_smooth_DoLP, fmt=".", color = self.smooth_I0_color, linestyle = 'none', markersize=self.marker_size, label="Smooth DoLP (" + str(bottle.smoothing_factor) + " " + str(bottle.smoothing_unit)[:3] + ")", zorder=1, errorevery=self.error_step)

		self.ax2.set_ylim(bottom = 0)
		self.ax2_lines.append([l_all_DoLP, l_all_DoLP.get_label()])
		self.ax2_lines.append([l_smooth_DoLP, l_smooth_DoLP.get_label()])



		if bottle.location.lower() == "skibotn" and bottle.filters == "br" and bottle.DateTime().date() == dt.date(2019, 3, 7):
			self.ax2.set_ylim((0, 5))
		# l_avg_DoLP, = self.ax2.plot(self.x_axis_list,[bottle.DoLP_average] * len(self.x_axis_list), color = self.smooth_I0_color, label="Avg: " + str(bottle.DoLP_average)[:4])
		# self.ax2_lines.append([l_avg_DoLP, l_avg_DoLP.get_label()])

		if self.mag_data:
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


		if self.eiscat_data.valid:
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


		if not self.show_error_bars:
			l_all_AoLP, = self.ax3.plot(self.x_axis_list, bottle.all_AoLP * RtoD, ".", color = self.all_I0_color, linestyle = 'none', markersize=self.marker_size, label="AoLP", zorder=0)

		if not self.show_smooth_error_bars:
			l_smooth_AoLP, = self.ax3.plot(self.x_axis_list, bottle.smooth_AoLP * RtoD, ".", color = self.smooth_I0_color, linestyle = 'none', markersize=self.marker_size, label="Smooth AoLP (" + str(bottle.smoothing_factor) + " " + str(bottle.smoothing_unit)[:3] + ")", zorder=1)

		if self.show_error_bars or self.show_smooth_error_bars:
			# self.ax3.fill_between(self.x_axis_list, bottle.all_AoLP * RtoD - bottle.std_AoLP, bottle.all_AoLP * RtoD + bottle.std_AoLP, color = "grey")
			# self.ax3.fill_between(self.x_axis_list, bottle.smooth_AoLP_lower * RtoD, bottle.smooth_AoLP_upper * RtoD, color = "yellow", alpha = 0.5)

			ls_all = dict()
			ls_smooth = dict()
			if bottle.graph_angle_shift == 1:
				min, mid, max = 0, 90, 180
			elif bottle.graph_angle_shift == 0:
				min, mid, max = -90, 0, 90

			for x, s_error, s_angle, error, angle in zip(self.x_axis_list[::self.error_step], bottle.std_smooth_AoLP[::self.error_step]*RtoD, bottle.smooth_AoLP[::self.error_step]*RtoD, bottle.std_AoLP[::self.error_step]*RtoD, bottle.all_AoLP[::self.error_step]*RtoD):

				if self.show_error_bars:
					if angle > mid and angle + error > max:
						temp = angle + error - 180
						ls_all.update({x:[min, temp]})
					elif angle < mid and angle - error < min:
						temp = angle - error + 180
						ls_all.update({x:[max, temp]})
				if self.show_smooth_error_bars:
					if s_angle > mid and s_angle + s_error > max:
						temp = s_angle + s_error - 180
						ls_smooth.update({x:[min, temp]})
					elif s_angle < mid and s_angle - s_error < min:
						temp = s_angle - s_error + 180
						ls_smooth.update({x:[max, temp]})

			self.ax3.set_ylim(min, max)
			if self.show_error_bars:
				l_all_AoLP = self.ax3.errorbar(self.x_axis_list, bottle.all_AoLP * RtoD, yerr = bottle.std_AoLP * RtoD, fmt=".", ecolor="grey", color = self.all_I0_color, linestyle = 'none', markersize=self.marker_size, label="AoLP", zorder=0, errorevery=self.error_step)
				for i, a in ls_all.items():
					self.ax3.vlines(i, a[0], a[1], colors = "grey", zorder=0) #"green", linewidth=5)#
			if self.show_smooth_error_bars:
				l_smooth_AoLP = self.ax3.errorbar(self.x_axis_list, bottle.smooth_AoLP * RtoD, yerr = bottle.std_smooth_AoLP * RtoD, fmt=".", color = self.smooth_I0_color, linestyle = 'none', markersize=self.marker_size, label="AoLP", zorder=1, errorevery=self.error_step)
				# plt.errorbar(list(ls.keys()), [-90, 90], yerr=list(ls.values()), fmt='C0 ')
				for i, a in ls_smooth.items():
					self.ax3.vlines(i, a[0], a[1], colors = self.smooth_I0_color, zorder=1) #"green", linewidth=5)#

		self.ax3_lines.append([l_all_AoLP, l_all_AoLP.get_label()])
		self.ax3_lines.append([l_smooth_AoLP, l_smooth_AoLP.get_label()])
		# l_avg_AoLP, = self.ax3.plot(self.x_axis_list,[bottle.AoLP_average * RtoD] * len(self.x_axis_list), color = self.smooth_I0_color, label="Avg: " + str(bottle.AoLP_average * RtoD)[:4])
		# self.ax3_lines.append([l_avg_AoLP, l_avg_AoLP.get_label()])

		if bottle.observation_type == "fixed":
			if bottle.AoRD is not False:
				if bottle.graph_angle_shift == 1:
					AoRD = SetAngleBounds(bottle.AoRD, 0, np.pi)
				elif bottle.graph_angle_shift == 0:
					AoRD = SetAngleBounds(bottle.AoRD, -np.pi/2, np.pi/2)
				l_AoRD, = self.ax3.plot(self.x_axis_list,[bottle.AoRD * RtoD] * len(self.x_axis_list), linewidth=3, color=self.AoRD_color, label="AoRD: " + str(bottle.AoRD*RtoD)[:5], zorder=2)
				self.ax3_lines.append([l_AoRD, l_AoRD.get_label()])
				# l54, = self.ax3.plot(self.x_axis_list,[bottle.AoRD_ortho * RtoD] * len(self.x_axis_list), ":g", linewidth=self.marker_size, label="AoRD ortho: " + str(bottle.AoRD_ortho*RtoD)[:5])
			if (bottle.AoBapp and bottle.AoBlos) is not False:
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
			print("DEBUG plot discrete")
			# l_AoRD, = self.ax3.plot(bottle.discrete_rotation_times, bottle.AoRD * RtoD, "k*", markersize=self.marker_size*(5+2), label="AoRD", zorder=2)
			# l_AoRD, = self.ax3.plot(bottle.discrete_rotation_times, bottle.AoRD * RtoD, "*", color = self.AoRD_color, markersize=self.marker_size*(5-2), label="AoRD", zorder=3)
			# self.ax3_lines.append([l_AoRD, l_AoRD.get_label()])
			#
			# l_AoBapp, = self.ax3.plot(bottle.discrete_rotation_times, bottle.AoBapp * RtoD, "k*", markersize=self.marker_size*(5+2), label="AoBapp", zorder=2)
			# l_AoBapp, = self.ax3.plot(bottle.discrete_rotation_times, bottle.AoBapp * RtoD, "*", color = self.AoBapp_color, markersize=self.marker_size*(5-2), label="AoBapp", zorder=3)
			# self.ax3_lines.append([l_AoBapp, l_AoBapp.get_label()])

			rot = 0
			if len(self.x_axis_ticks_pos) > 15:
				rot = 60

			self.ax3.set_xticks(self.x_axis_ticks_pos)
			self.ax3.set_xticklabels(self.x_axis_ticks_label, rotation = rot)
			# self.ax3.xticks(rotation = 30)

		elif self.xaxis_azimut and bottle.observation_type == "fixed_azimut_discrete_rotation":
			print("DEBUG plot discrete")
			l_AoRD, = self.ax3.plot(bottle.discrete_rotation_times, bottle.AoRD * RtoD, "*", color = self.AoRD_color, markersize=self.marker_size*4, label="AoRD", zorder=2)
			self.ax3_lines.append([l_AoRD, l_AoRD.get_label()])
			l_AoBapp, = self.ax3.plot(bottle.discrete_rotation_times, bottle.AoBapp * RtoD, "*", color = self.AoBapp_color, markersize=self.marker_size*4, label="AoBapp", zorder=2)
			self.ax3_lines.append([l_AoBapp, l_AoBapp.get_label()])

		elif self.xaxis_azimut and bottle.observation_type == "fixed_elevation_continue_rotation":

			ax1_lim = self.ax1.get_ylim()
			self.ax1.set_ylim(ax1_lim)
			ax2_lim = self.ax2.get_ylim()
			ax3_lim = self.ax3.get_ylim()
			for i in range(bottle.nb_continue_rotation):
				print("DEBUG ROT TIMES", i * 360 + (bottle.source_azimut * RtoD)%360)

				self.ax1.plot([i * 360 + (bottle.source_azimut * RtoD)%360, i * 360 + (bottle.source_azimut * RtoD)%360], ax1_lim, "--k")
				self.ax2.plot([i * 360 + (bottle.source_azimut * RtoD)%360, i * 360 + (bottle.source_azimut * RtoD)%360], ax2_lim, "--k")
				self.ax3.plot([i * 360 + (bottle.source_azimut * RtoD)%360, i * 360 + (bottle.source_azimut * RtoD)%360], ax3_lim, "--k")
				if i != 0: # Draw red lines to delimit the rotations
					print(i * 360)
					self.ax1.plot([i * 360, i * 360], ax1_lim, "r")
					self.ax2.plot([i * 360, i * 360], ax2_lim, "r")
					self.ax3.plot([i * 360, i * 360], ax3_lim, "r")


				l_AoRD, = self.ax3.plot(np.linspace(i * 360, (i+1) * 360, len(bottle.AoRD)), bottle.AoRD * RtoD, "*", color = self.AoRD_color, markersize=self.marker_size, label="AoRD", zorder=2)
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

		self.f1.subplots_adjust(hspace=0)

		###Set title
		# self.f1.suptitle(bottle.saving_name.replace("_", " "))


		# plt.setp([a.get_xticklabels() for a in f1.axes[:-1]], visible=False)
		if bottle.observation_type == "fixed":
			plt.minorticks_on()

		self.ax12 = self.ax1.twiny()
		self.ax12.set_xlim(self.ax1.get_xlim())
		xticks = plt.xticks()[0] * self.divisor
		xticks = [bottle.DateTime() + time.timedelta(seconds = t) for t in xticks]
		# xticks = [bottle.DateTime() + time.timedelta(seconds = t) + bottle.head_jump for t in xticks]
		self.ax12.set_xticklabels([st.strftime("%H:%M") for st in xticks])
		self.ax12.set_xlabel("UT")

		# self.ax1.legend(list(zip(*self.ax1_lines))[0], list(zip(*self.ax1_lines))[1])
		# self.ax2.legend(list(zip(*self.ax2_lines))[0], list(zip(*self.ax2_lines))[1])
		# self.ax3.legend(list(zip(*self.ax3_lines))[0], list(zip(*self.ax3_lines))[1], loc = "lower center")

		print("Saving graphs in", bottle.data_file_name + "/" + bottle.saving_name + '_graphs.png')
		if bottle.instrument_name in ["carmen", "corbel", "gdcu"]:
			plt.savefig(bottle.data_file_name + "/" + bottle.saving_name + '_graphs.png', bbox_inches='tight')
			# plt.savefig(bottle.data_file_name + "/" + bottle.saving_name + '_graphs.eps', bbox_inches='tight')
		else:
			plt.savefig("/".join(bottle.data_file_name.split("/")[:-1]) + "/" + bottle.saving_name + '_graphs.png', bbox_inches='tight')
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
		ax1.set_ylabel("DoLP (\%)")
		ax1.legend()

		if  not bottle.NoVref:
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
