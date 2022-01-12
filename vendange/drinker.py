#!/usr/bin/python3
# -*-coding:utf-8 -*

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage{lmodern}']
from scipy import signal
import sys
import os
from subprocess import call

from utils import *
from rotation import *
from bottle import *


class Mixer:
	def __init__(self, bottle, comp_bottles = False, mode = "default"):
		self.bottle = bottle
		self.secondary_bottle = comp_bottles
		self.mode = mode

		self.MakeXAxis()
		self.MakeFigure()
		self.SetGraphParameter()

		self.MakePlots()
		self.MakeCorrelationPlots()


	def MakeXAxis(self):
		self.divisor = 1.
		if self.bottle.all_times[-1] > 3600 * 2:
			self.divisor = 3600.
			self.xlabel = "Time (hours)"
		elif self.bottle.all_times[-1] > 60:
			self.divisor = 60.
			self.xlabel = "Time (minutes)"
		else:
			self.xlabel = "Time (seconds)"

		self.x_axis_list 	   = np.array(self.bottle.all_times) / self.divisor
		self.x_axis_list_stamp = [t + self.bottle.time + self.bottle.head_jump for t in self.x_axis_list]
		self.x_axis_list_stamp = [time.strftime("%H:%M:%S", time.gmtime(st)) for st in self.x_axis_list_stamp]

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

	def SetGraphParameter(self):
		self.marker_size = 3

		self.smooth_I0_color = "red"
		self.smooth_ref_color = "green"
		### xkcd color guide: https://xkcd.com/color/rgb/
		if self.bottle.filters:
			self.pola_color = self.bottle.filters[0]
			self.ref_color = self.bottle.filters[1]
			if 	 self.pola_color == "r": self.smooth_I0_color = "xkcd:red"
			elif self.pola_color == "v": self.smooth_I0_color = "xkcd:green"
			elif self.pola_color == "b": self.smooth_I0_color = "xkcd:blue"
			elif self.pola_color == "m": self.smooth_I0_color = "xkcd:purple"
			else: self.smooth_I0_color = "red"
			if 	 self.ref_color == "r" and self.pola_color != "r": self.smooth_ref_color = "xkcd:red"
			elif self.ref_color == "r" and self.pola_color == "r": self.smooth_ref_color = "xkcd:orange"
			elif self.ref_color == "v" and self.pola_color != "v": self.smooth_ref_color = "xkcd:green"
			elif self.ref_color == "v" and self.pola_color == "v": self.smooth_ref_color = "xkcd:lime green"
			elif self.ref_color == "b" and self.pola_color != "b": self.smooth_ref_color = "xkcd:blue"
			elif self.ref_color == "b" and self.pola_color == "b": self.smooth_ref_color = "xkcd:bright blue"
			elif self.ref_color == "m" and self.pola_color != "m": self.smooth_ref_color = "xkcd:purple"
			elif self.ref_color == "m" and self.pola_color == "m": self.smooth_ref_color = "xkcd:lavender"
			else: self.smooth_ref_color = "green"

	def MakePlots(self):
		#Plotting the mean intensity I0, DoLP, AoLP for each rotation
		# if self.bottle.NoVref:
		l_all_I0, = self.ax1.plot(self.x_axis_list, self.bottle.all_I0, "k.", linestyle = 'none', markersize=self.marker_size, label="Intensity")
		self.ax1_lines.append([l_all_I0, l_all_I0.get_label()])

		if not self.bottle.NoVref:
			self.ax13 = self.ax1.twinx()
			# ax13.set_xlim(ax1.get_xlim())
			# xticks = plt.xticks()[0] * self.divisor
			# xticks = [t + self.bottle.time + self.bottle.head_jump for t in xticks]
			# ax12.set_xticklabels([time.strftime("%H:%M:%S", time.localtime(st)) for st in xticks])
			l_smooth_Iref, = self.ax13.plot(self.x_axis_list, self.bottle.smooth_Iref, ".", color = self.smooth_ref_color, linestyle = 'none', markersize=self.marker_size, label="Smooth Ref Intensity (" + str(self.bottle.smoothing_factor) + " " + str(self.bottle.smoothing_unit)[:3] + ")")
			self.ax1_lines.append([l_smooth_Iref, l_smooth_Iref.get_label()])
			l_avg_Iref, = self.ax13.plot(self.x_axis_list, [self.bottle.Iref_average] * len(self.x_axis_list), color = self.smooth_ref_color, label="Avg Ref Intensity " + str(self.bottle.Iref_average)[:4])
			self.ax1_lines.append([l_avg_Iref, l_avg_Iref.get_label()])
			self.ax13.set_ylabel("Ref intensity")

		l_smooth_I0, = self.ax1.plot(self.x_axis_list, self.bottle.smooth_I0, ".", color = self.smooth_I0_color, linestyle = 'none', markersize=self.marker_size, label="Smooth Intensity (" + str(self.bottle.smoothing_factor) + " " + str(self.bottle.smoothing_unit)[:3] + ")")
		self.ax1_lines.append([l_smooth_I0, l_smooth_I0.get_label()])
		l_avg_I0, = self.ax1.plot(self.x_axis_list,[self.bottle.I0_average] * len(self.x_axis_list), color = self.smooth_I0_color, label="Avg: " + str(self.bottle.I0_average)[:4])
		self.ax1_lines.append([l_avg_I0, l_avg_I0.get_label()])

		l_all_DoLP, = self.ax2.plot(self.x_axis_list, self.bottle.all_DoLP, "k.", linestyle = 'none', markersize=self.marker_size, label="DoLP")
		self.ax2_lines.append([l_all_DoLP, l_all_DoLP.get_label()])
		l_smooth_DoLP, = self.ax2.plot(self.x_axis_list, self.bottle.smooth_DoLP, ".", color = self.smooth_I0_color, linestyle = 'none', markersize=self.marker_size, label="Smooth DoLP (" + str(self.bottle.smoothing_factor) + " " + str(self.bottle.smoothing_unit)[:3] + ")")
		self.ax2_lines.append([l_smooth_DoLP, l_smooth_DoLP.get_label()])
		l_avg_DoLP, = self.ax2.plot(self.x_axis_list,[self.bottle.DoLP_average] * len(self.x_axis_list), color = self.smooth_I0_color, label="Avg: " + str(self.bottle.DoLP_average)[:4])
		self.ax2_lines.append([l_avg_DoLP, l_avg_DoLP.get_label()])

		l_all_AoLP, = self.ax3.plot(self.x_axis_list, self.bottle.all_AoLP * RtoD, "k.", linestyle = 'none', markersize=self.marker_size,  label="AoLP")
		self.ax3_lines.append([l_all_AoLP, l_all_AoLP.get_label()])
		l_smooth_AoLP, = self.ax3.plot(self.x_axis_list, self.bottle.smooth_AoLP * RtoD, ".", color = self.smooth_I0_color, linestyle = 'none', markersize=self.marker_size, label="Smooth AoLP (" + str(self.bottle.smoothing_factor) + " " + str(self.bottle.smoothing_unit)[:3] + ")")
		self.ax3_lines.append([l_smooth_AoLP, l_smooth_AoLP.get_label()])
		l_avg_AoLP, = self.ax3.plot(self.x_axis_list,[self.bottle.AoLP_average * RtoD] * len(self.x_axis_list), color = self.smooth_I0_color, label="Avg: " + str(self.bottle.AoLP_average * RtoD)[:4])
		self.ax3_lines.append([l_avg_AoLP, l_avg_AoLP.get_label()])
		if self.bottle.AoRD is not False:
			l_AoRD, = self.ax3.plot(self.x_axis_list,[self.bottle.AoRD * RtoD] * len(self.x_axis_list), "g", label="AoRD: " + str(self.bottle.AoRD*RtoD)[:5])
			self.ax3_lines.append([l_AoRD, l_AoRD.get_label()])
			# l54, = ax3.plot(self.x_axis_list,[self.bottle.AoRD_ortho * RtoD] * len(self.x_axis_list), ":g", linewidth=self.marker_size, label="AoRD ortho: " + str(self.bottle.AoRD_ortho*RtoD)[:5])
		if (self.bottle.AoBapp and self.bottle.AoBlos) is not False:
			l_AoBapp, = self.ax3.plot(self.x_axis_list,[self.bottle.AoBapp * RtoD] * len(self.x_axis_list), "b", label="AoBapp: " + str(self.bottle.AoBapp*RtoD)[:5])
			self.ax3_lines.append([l_AoBapp, l_AoBapp.get_label()])
			l_AoBapp_ortho, = self.ax3.plot(self.x_axis_list,[self.bottle.AoBapp_ortho * RtoD] * len(self.x_axis_list), ":b", linewidth=self.marker_size, label="AoBapp ortho: " + str(self.bottle.AoBapp_ortho*RtoD)[:5])
			self.ax3_lines.append([l_AoBapp_ortho, l_AoBapp_ortho.get_label()])

			l_AoBlos, = self.ax3.plot(self.x_axis_list,[self.bottle.AoBlos * RtoD] * len(self.x_axis_list), "orange", label="AoBlos: " + str(self.bottle.AoBlos*RtoD)[:5])
			self.ax3_lines.append([l_AoBlos, l_AoBlos.get_label()])

		# self.ax4.plot(self.x_axis_list, self.bottle.all_TempPM, "k.", linestyle = 'none', markersize=self.marker_size, label="PM")
		# self.ax4.plot(self.x_axis_list, self.bottle.all_TempOptical, "r.", linestyle = 'none', markersize=self.marker_size, label="Optical")
		# # self.ax4.plot(self.x_axis_list, self.bottle.all_TempAmbiant, "b.", linestyle = 'none', markersize=self.marker_size, label="Ambiant")
		# # self.ax2.plot(self.x_axis_list,[DoLP_average] * nb_rot, "b", label="Avg: " + str(DoLP_average))

		self.f1.subplots_adjust(hspace=0)
		self.f1.suptitle(self.bottle.saving_name.replace("_", " "))
		# plt.setp([a.get_xticklabels() for a in f1.axes[:-1]], visible=False)
		plt.minorticks_on()

		self.ax12 = self.ax1.twiny()
		self.ax12.set_xlim(self.ax1.get_xlim())
		xticks = plt.xticks()[0] * self.divisor
		xticks = [t + self.bottle.time + self.bottle.head_jump for t in xticks]
		self.ax12.set_xticklabels([time.strftime("%H:%M:%S", time.localtime(st)) for st in xticks])
		self.ax12.set_xlabel("UT")

		self.ax1.legend(list(zip(*self.ax1_lines))[0], list(zip(*self.ax1_lines))[1])
		self.ax2.legend(list(zip(*self.ax2_lines))[0], list(zip(*self.ax2_lines))[1])
		self.ax3.legend(list(zip(*self.ax3_lines))[0], list(zip(*self.ax3_lines))[1])

		print("Saving graphs in", self.bottle.data_file_name + "/" + self.bottle.saving_name + '_graphs.png')
		if self.bottle.instrument_name == "ptcu":
			plt.savefig(self.bottle.data_file_name + "/" + self.bottle.saving_name + '_graphs.png', bbox_inches='tight')
		else:
			plt.savefig("/".join(self.bottle.data_file_name.split("/")[:-1]) + "/" + self.bottle.saving_name + '_graphs.png', bbox_inches='tight')


	def MakeCorrelationPlots(self):
		f2, (ax1, ax2) = plt.subplots(2, figsize=(10, 20))

		ax1.plot(self.bottle.all_I0, self.bottle.all_DoLP, "k+", label="Raw")
		ax1.plot(self.bottle.smooth_I0, self.bottle.smooth_DoLP, "r+", label="Smoothed (" + str(self.bottle.smoothing_factor) + " " + str(self.bottle.smoothing_unit)[:3] + ")")
		ax1.set_xlabel("Intensity (mV)")
		ax1.set_ylabel("DoLP (\%)")
		ax1.legend()

		if  not self.bottle.NoVref:
			ax2.plot(self.bottle.all_I0, self.bottle.all_Iref, "k+", label="Raw")
			ax2.plot(self.bottle.smooth_I0, self.bottle.smooth_Iref, "r+", label="Smoothed (" + str(self.bottle.smoothing_factor) + " " + str(self.bottle.smoothing_unit)[:3] + ")")
			ax2.set_xlabel("Ref Intensity (mV)")
			ax2.set_ylabel("Pola Intensity (mV)")

			ax2.legend()

		print("Saving correlation in", self.bottle.data_file_name + "/" + self.bottle.saving_name + '_correlations.png')
		if self.bottle.instrument_name == "ptcu":
			plt.savefig(self.bottle.data_file_name + "/" + self.bottle.saving_name + '_correlation.png', bbox_inches='tight')
		else:
			plt.savefig("/".join(self.bottle.data_file_name.split("/")[:-1]) + "/" + self.bottle.saving_name + '_correlation.png', bbox_inches='tight')
