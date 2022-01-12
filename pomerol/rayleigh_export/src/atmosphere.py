#!/usr/bin/python3
# -*-coding:utf-8 -*

import sys as sys
import numpy as np
import time as tm
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Circle
from matplotlib.patches import Ellipse
from matplotlib.patches import Arrow
# from matplotlib.lines import mlines

import osgeo.gdal as gdal
gdal.UseExceptions()  # not required, but a good idea

import imageio

from observation import *
from rayleigh_utils import *

class Atmosphere:
	def __init__(self, in_dict):

		self.h_r_min			= float(in_dict["RS_min_altitude"])  			# Minimum altitude of the scatering layer
		self.h_r_max			= float(in_dict["RS_max_altitude"]) 			# Maximum altitude of the scatering layer
		# if self.h_r_max == self.h: 		# h_r_max must be != h or it crashes...
		# 	self.h_r_max += 0.01
		self.d_los = float(in_dict["resolution_along_los"]) 	# Lenth bin along line of sight


		self.d_0 	= 2.5468 * 10 ** 25 	#=1.2250 kg/m3	sea level density
		self.T 		= 288.15 				# K, sea level standard temperature
		self.g 		= 9.80665 				# m/s², earth-surface gravitational acceleration
		self.M 		= 0.0289654 			# kg/mol, molar mass of dry air
		self.R 		= 8.31447 				# J/(mol·K), ideal (universal) gas constant

	def GetVolume(self, AR, ouv_pc):
		"""Get the volume of a truncated cone of length d_los, and opening angle ouv_pc, centered at distance AR."""
		V = (np.pi * np.tan(ouv_pc) ** 2) / 3

		h1, h2 = AR - self.d_los, AR + self.d_los

		V *= (h2 - h1)
		V *= ((h2 + h1) ** 2 - h1 * h2)

		# print("V", V * (10 ** 9))
		return V * (10 ** 9) # in m3. (multiply by e9 because every distances is in km)


	def GetRSCrossSectionParticle(self, theta, wl = 557.7, Dp = 0.315, n = 1.000271375):
		"""Get Rayleigh scattering cross section for 1 angle. (See wikipedia or Seinfeld, John H. and Pandis, Spyros N. (2006) Atmospheric Chemistry and Physics, 2nd Edition, John Wiley and Sons, New Jersey, Chapter 15.1.1, ISBN 0471720186).
		Dp : particle size
		wl : wavelength of the light
		theta : scattering angle
		n : refractive index"""

		cs  = ((Dp * 10 ** (-9)) ** 6) / 8
		cs *= (np.pi / (wl * 10 ** (-9))) ** 4
		cs *= ((n ** 2 - 1) / (n ** 2 + 2)) ** 2
		cs *= (1 + np.cos(theta) ** 2)

		# print("cs", cs)

		return cs


	def GetRSCrossSectionMolecule(self, theta, wl = 557.7, alpha=2.3732 * 10**(-39)):
		"""Get Rayleigh scattering cross section for 1 angle. (See wikipedia).
		wl : wavelength of the light
		theta : scattering angle
		alpha : polarizability (in SI units !!! alpha(cm**3) = 8.988*10**15 alpha(SI) """

		cs  = (8 * np.pi ** 4 * alpha ** 2) / (wl ** 4)
		cs *= (1 + np.cos(theta) ** 2)

		return cs

	def GetParticuleDensity(self, h):
		"""Should give the particle density at altitude h in particules/m3."""

		h /= 10**3 #alt in meters

		A = (- self.g * self.M * h) / (self.R * self.T)

		d = self.d_0 * np.exp(A) #in particles / m3

		return d
