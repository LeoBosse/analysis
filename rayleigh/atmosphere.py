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

		Ang3toMeter3 = 10 ** (-30)
		self.names = ["N2", "O2", "H2O", "CO2", "O3"]
		self.polarizability = np.array([1.710, 1.562, 1.501, 2.507, 3.079]) * Ang3toMeter3

		self.profile_name = in_dict["Atmospheric_profile"] #name of the atmospheric profile file to use
		self.LoadAllProfiles(self.profile_name)

		# self.MakePlots()
		# plt.show()



	def LoadAllProfiles(self, filename):
		"""Load all atmospheric profiles from a file into a dictionnary. Download the file from: http://eodg.atm.ox.ac.uk/RFM/atm/ (MIPAS model 2001 works for sure)
		The pressure is converted in Pascals (1 mb = 100 Pa)
		Unit of density is part per million per volume, convert to fraction per volume"""
		self.profiles = dict()
		unit_factor = 1
		with open(filename, 'r') as f:
			inprofile = False #To know if we read a header or title line
			for line in f: #Read all lines one by one
				if line[0] == "*": #If title line. As the form: *NAME [UNIT]
					inprofile = True #We start to read a profile
					name = line.split()[0][1:] #Get the name of the profile
					if name == "END": break
					unit = line.split()[1][1:-1] #Get the unit
					if   unit == "mb": unit_factor = 100 #If unit pressure is millibars, convert to pascals
					elif unit == "ppmv": unit_factor = 10**(-6) #If unit of density is part per million per volume, convert to fraction per volume
					else: unit_factor = 1
					self.profiles[name] = [] #Create a new list in our dictionnary.
				elif inprofile: #If reading a profile line
					for word in line.split(): #Get all values one by one
						self.profiles[name].append(float(word) * unit_factor) #Save values in the dictionnary

		self.nb_profile_alt = len(self.profiles["HGT"]) #Number of points in the profiles

	def GetProfileValue(self, alt, name):
		"""From any altitude in km, gives the profile value with linear approx."""

		if alt in self.profiles["HGT"]:
			i = self.profiles["HGT"].index(alt)
			return self.profiles[name][i]
		else:
			i = 0
			while i < self.nb_profile_alt and alt > self.profiles["HGT"][i]:
				i += 1

			dalt = 0
			dprof = 0
			d = 0
			value = 0
			ip, im = i, i - 1
			if ip < self.nb_profile_alt:
				dalt = self.profiles["HGT"][ip] - self.profiles["HGT"][im]
				dprof = self.profiles[name][ip] - self.profiles[name][im]
				d = dprof / dalt
				value = self.profiles[name][im]

			return value + d * (alt - self.profiles["HGT"][im])


	def GetVolume(self, AR, ouv_pc, da = None, dr=None, unit="m", type="cone"):
		"""Get the volume of a truncated cone of length d_los, and half opening angle ouv_pc.
		Unit defines the unit you want it in. can be "cm", "m" or "km". (the cube is implicit)"""
		if dr is None:
			dr = self.d_los

		if type == "cone":
			V = (np.pi * np.tan(ouv_pc) ** 2) / 3

			h1, h2 = AR - dr/2, AR + dr/2

			V *= dr
			V *= (h1 ** 2 + h2 ** 2 + h1 * h2)
		elif type == "pole":

			R1, R2 = AR - dr/2., AR + dr/2.
			V = 2. * np.pi / 3. * (1 - np.cos(ouv_pc)) * (R2 ** 3 - R1 ** 3)

		elif type == "slice" and len(ouv_pc) == 2:

			R1, R2 = AR - dr/2., AR + dr/2.

			V = da / 3. * (R2 ** 3 - R1 ** 3) * (np.cos(ouv_pc[0]) - np.cos(ouv_pc[1]))

		else:
			raise ValueError("Incorrect arguments passed to atmosphere.GetVolume()")



		# Every length is given in km
		if unit == "m":
			V *= 10 ** 9
		elif unit == "cm":
			V *= 10 ** 15
		elif unit == "km":
			V *= 1
		else:
			V *= 1
			print("WARNING: Unit for diffusion volume is not correct. Set to km by default.")

		# print("V", V * (10 ** 9))

		return V


	def GetRSVolumeCS(self, wl, alt):
		"""Return the Rayleigh scattering volume cross section as calculated in Bucholtz 95 in km-1.
		"""

		beta_s = self.GetSquareLawFit(wl, "Volume CS") #in cm-1
		P = self.GetProfileValue(alt, "PRE") #in Pa
		T = self.GetProfileValue(alt, "TEM") #in K

		return beta_s * (P / 101325) * (288.15 / T)

	def GetRSOpticalDepth(self, wl, E_alt, R_alt):
		"""Return VERTICAL optical depth between two altitudes as calculated in Bucholtz 95."""
		tau_E = self.GetSquareLawFit(wl, "Optical Depth")
		P_E = self.GetProfileValue(E_alt, "PRE")
		tau_E *= (P_E / 101325)

		tau_R = self.GetSquareLawFit(wl, "Optical Depth")
		P_R = self.GetProfileValue(R_alt, "PRE")
		tau_R *= (P_R / 101325)

		tau_ER = abs(tau_E - tau_R)
		return tau_ER

	def GetSquareLawFit(self, wl, purpose):
		"""Return the result of the square law fit function presented in Bucholtz 1995. Purpose can be:
		"Single CS" = For single air particule rayleigh cross section in cm2.
		"Volume CS" = For volume cross rayleigh section in km-1.
		"Optical Depth" = For the rayleigh optical depth.
		Value of parameter A for optical depth is given for a Sub arctic winter atmosphere model."""

		A, B, C, D = self.GetSquareLawParam(wl, purpose)
		wl /= 1000  #wavelength given in nanometers, converted to micrometer
		E = - (B + C * wl + D / wl)
		return A * wl ** E


	def GetRSPhaseFunction(self, wl, theta):
		if   wl == 391.4:
			gamma = 1.499 * 0.01
		elif wl == 427.8:
			gamma = 1.483 * 0.01
		elif wl == 557.7:
			gamma = 1.442 * 0.01
		elif wl == 630:
			gamma = 1.413 * 0.01

		### Simple approx
		# A = 3. / 4
		# B = 1 + np.cos(theta) ** 2

		### Chandrasekhar formula
		A = 3 / (4 + 8 * gamma)
		B = (1 + 3 * gamma) + (1 - gamma) * np.cos(theta) ** 2

		return A * B

	def GetSquareLawParam(self, wl, purpose):
		"""Return the necessary parameters for calculating the result of the square law fit function presented in Bucholtz 1995. Purpose can be:
		"Single CS" = For single air particule rayleigh cross section.
		"Volume CS" = For volume rayleigh cross section.
		"Optical Depth" = For the rayleigh optical depth.
		Value of parameter A for optical depth is given for a Sub arctic winter atmosphere model."""
		if wl < 500:
			if purpose == "Single CS":
				A = 3.01577 * 10 ** (-28) #for cm2
			elif purpose == "Volume CS":
				A = 7.68246 * 10 ** (-4) #for km-1
			elif purpose == "Optical Depth":
				A = 6.49997 * 10 ** (-3)
			B = 3.55212
			C = 1.35579
			D = 0.11563
		else:
			if purpose == "Single CS":
				A = 4.01061 * 10 ** (-28) #for cm2
			elif purpose == "Volume CS":
				A = 10.21675 * 10 ** (-4) #for km-1
			elif purpose == "Optical Depth":
				A = 8.64145 * 10 ** (-3)
			B = 3.99668
			C = 0.00110298
			D = 0.0271393

		return A, B, C, D


	def MakePlots(self):
		f1, (ax1, ax2, ax3) = plt.subplots(3, sharey = True, figsize=(16, 8))
		ax1 = plt.subplot(131)
		ax2 = plt.subplot(132)
		ax3 = plt.subplot(133)

		ax1.plot(self.profiles["TEM"], self.profiles["HGT"])
		ax2.plot(self.profiles["PRE"], self.profiles["HGT"])
		ax3.plot(self.profiles["N2"], self.profiles["HGT"], label="N2")
		ax3.plot(self.profiles["O2"], self.profiles["HGT"], label="O2")
		ax3.plot(self.profiles["H2O"], self.profiles["HGT"], label="H2O")
		ax3.plot(self.profiles["CO2"], self.profiles["HGT"], label="CO2")
		ax3.plot(self.profiles["O3"], self.profiles["HGT"], label="O3")

		ax3.legend()
		ax3.set_xscale('log')

		ax1.set_title("Temperature")
		ax2.set_title("Pressure")
		ax3.set_title("Species")

	#
	# def GetRSCrossSection(self, alt, theta, wl = 557.7):
	# 	k = 2 * np.pi / (wl * 10 ** (-9))
	# 	alpha = self.GetPolarizability(alt)
	# 	A_ptcu = 5 * 10 ** (-4)
	#
	# 	cs = (1 + np.cos(theta) ** 2) * k ** 4 * alpha ** 2 / 2 / A_ptcu
	#
	# 	return cs


	# def GetPolarizability(self, alt):
	#
	# 	fraction 	= np.array([self.GetProfileValue(alt, n) for n in self.names])
	#
	# 	polarizability =  np.sum(self.polarizability * fraction) / np.sum(fraction)
	#
	# 	return polarizability
	# def GetEffCrossSection(self, alt, wl = 557.7):
	# 	"""Give the effective cross section of the atmosphere at a given altitude. From Born&Wolf (edition 7), formula 103 (p785) and 85 (p778). In the case a<<wl, only l=1 is necessary."""
	# 	a = 0.1 #Sphere radius
	# 	n = 1.00027716 #Complex refractive index
	#
	# 	eB_1 = 1j * (2 * np.pi * a / wl) ** 3 * (n**2 - 1) / (n**2 + 2)
	# 	mB_1 = 1j * (2 * np.pi * a / wl) ** 5 * (n**2 - 1) / 30
	# 	eB_2 = - (2 * np.pi * a / wl) ** 5 / 36 * (n**2 - 1) / (n**2 + 3/2)
	#
	# 	sum = - 2 * (eB_1 + mB_1) + 1j * 6 * eB_2
	# 	# print(sum, sum.real)
	#
	# 	Q = wl ** 2 / (2 * np.pi) * sum.real
	#
	# 	return Q
	#
	# def GetCExt(self, alt):
	# 	# print(self.GetEffCrossSection(alt),  self.GetParticuleDensity(alt), self.GetEffCrossSection(alt) * self.GetParticuleDensity(alt))
	# 	return self.GetEffCrossSection(alt) * self.GetParticuleDensity(alt)

	# def GetRSOpticalDepth(self, wl, alt):
	# 	tau_s = self.GetSquareLawFit(wl, "Optical Depth")
	# 	P = self.GetProfileValue(alt, "PRE")
	# 	return tau_s * (P / 101300)



	# def GetScattered(self, wl, alt, theta):
	# 	return self.GetRSVolumeCS(wl, alt) * self.GetRSPhaseFunction(wl, theta) / 4 / np.pi


	# def GetRSCrossSectionParticle(self, theta, wl = 557.7, Dp = 0.315, n = 1.000271375):
	# 	"""Get Rayleigh scattering cross section for 1 angle. (See wikipedia or Seinfeld, John H. and Pandis, Spyros N. (2006) Atmospheric Chemistry and Physics, 2nd Edition, John Wiley and Sons, New Jersey, Chapter 15.1.1, ISBN 0471720186).
	# 	Dp : particle size
	# 	wl : wavelength of the light
	# 	theta : scattering angle
	# 	n : refractive index"""
	#
	# 	cs  = ((Dp * 10 ** (-9)) ** 6) / 8
	# 	cs *= (np.pi / (wl * 10 ** (-9))) ** 4
	# 	cs *= ((n ** 2 - 1) / (n ** 2 + 2)) ** 2
	# 	cs *= (1 + np.cos(theta) ** 2)
	#
	# 	# print("cs", cs)
	#
	# 	return cs
	#
	#
	# def GetRSCrossSectionMolecule(self, theta, wl = 557.7, alpha=2.3732 * 10**(-39)):
	# 	"""Get Rayleigh scattering cross section for 1 angle. (See wikipedia).
	# 	wl : wavelength of the light
	# 	theta : scattering angle
	# 	alpha : polarizability (in SI units !!! alpha(cm**3) = 8.988*10**15 alpha(SI) """
	#
	# 	cs  = (8 * np.pi ** 4 * alpha ** 2) / (wl ** 4)
	# 	cs *= (1 + np.cos(theta) ** 2)
	#
	# 	return cs

	# def GetParticuleDensity(self, h):
	# 	"""Should give the particle density at altitude h in particules/m3."""
	#
	# 	T = self.GetProfileValue(h, "TEM")
	# 	P = self.GetProfileValue(h, "PRE")
	# 	k_B = 1.380649 * 10**(-23)
	#
	# 	d = P / (k_B * T) #in particles / m3
	#
	# 	return d
