#!/usr/bin/python3
# -*-coding:utf-8 -*

from mpi4py import MPI
mpi_comm = MPI.COMM_WORLD
mpi_rank = mpi_comm.Get_rank()
mpi_size = mpi_comm.Get_size()
mpi_name = mpi_comm.Get_name()

import sys as sys
import os as os
import subprocess as subp
import numpy as np
import pandas as pd
import time as tm
import scipy.constants as cst
import scipy.integrate as integrate
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Circle
from matplotlib.patches import Ellipse
from matplotlib.patches import Arrow
matplotlib.rcParams['text.usetex'] = True
# from matplotlib.lines import mlines

# import osgeo.gdal as gdal
# gdal.UseExceptions()  # not required, but a good idea

import imageio
# import numba as nba

from observation import *
from rayleigh_utils import *
# from pomerol.cythonized import *

class Atmosphere:
	def __init__(self, in_dict):

		self.wavelength = float(in_dict["wavelength"]) #in nanometers

		self.show_cross_section = False
		if "show_cross_section" in in_dict:
			self.show_cross_section = bool(int(in_dict["show_cross_section"]))

		self.show_atmosphere = False
		if "show_atmosphere" in in_dict:
			self.show_atmosphere = bool(int(in_dict["show_atmosphere"]))

		self.show_aer_size = False
		if "show_aer_size" in in_dict:
			self.show_aer_size = bool(int(in_dict["show_aer_size"]))


		self.h_r_min			= float(in_dict["RS_min_altitude"])  			# Minimum altitude of the scatering layer
		self.h_r_max			= float(in_dict["RS_max_altitude"]) 			# Maximum altitude of the scatering layer
		# if self.h_r_max == self.h: 		# h_r_max must be != h or it crashes...
		# 	self.h_r_max += 0.01
		self.Nlos = int(in_dict["Nb_points_along_los"]) 	# Lenth bin along line of sight
		# self.d_los = float(in_dict["resolution_along_los"]) 	# Lenth bin along line of sight

		# self.d_0 	= 2.5468 * 10 ** 25 	#=1.2250 kg/m3	sea level density
		# self.T 		= 288.15 				# K, sea level standard temperature
		# self.g 		= 9.80665 				# m/s², earth-surface gravitational acceleration
		# self.M 		= 0.0289654 			# kg/mol, molar mass of dry air
		# self.R 		= 8.31447 				# J/(mol·K), ideal (universal) gas constant

		self.use_ozone 		= bool(int(in_dict["use_ozone"]))
		self.use_aerosol 	= bool(int(in_dict["use_aerosol"]))

		Ang3toMeter3 = 10 ** (-30)
		# self.names = ["N2", "O2", "H2O", "CO2", "O3"]
		# self.polarizability = np.array([1.710, 1.562, 1.501, 2.507, 3.079]) * Ang3toMeter3

		self.profile_name = in_dict["Atmospheric_profile"] #name of the atmospheric profile file to use
		self.LoadAllProfiles(self.profile_name)


		self.P0 = 101325 # in Pa
		self.T0 = 288.15 # in Kelvin

		self.beta_0 = self.GetSquareLawFit(self.wavelength, "Volume CS") / 1 #in km-1
		self.tau_0 = self.GetSquareLawFit(self.wavelength, "Optical Depth")

		self.profiles["beta_ray"] = np.array([self.RSVolumeCS(a) for a in self.profiles["HGT"]])
		self.profiles["tau_ray"]  = np.array([self.GetRSAbsorptionCoeff(a) for a in self.profiles["HGT"]])

		self.profiles["sca_angle"] 		= np.linspace(0, np.pi, 100)
		self.profiles["ray_Phase_Fct"]  = np.array([self.RSPhaseFunction(self.wavelength, t) for t in self.profiles["sca_angle"]])
		self.profiles["ray_Phase_DoLP"]  = np.array([self.RSPhaseFunctionDoLP(t, f=1) for t in self.profiles["sca_angle"]])

		self.LoadO3AbsorptionCS()

		self.profiles["beta_aer"] = np.zeros_like(self.profiles['HGT'])
		self.profiles["aer_Phase_Fct"] = np.zeros_like(self.profiles['sca_angle'])
		self.profiles["aer_Phase_DoLP"] = np.zeros_like(self.profiles['sca_angle'])
		if self.use_aerosol:
			self.aer_complexity = int(in_dict['aer_complexity'])
			if self.aer_complexity > 0:
				self.aerosol = Aerosol(in_dict)

				self.profiles["aer_Phase_Fct"], self.profiles["aer_Phase_DoLP"]  = zip(*[self.AerosolPhaseFunction(t) for t in self.profiles["sca_angle"]])
				self.profiles["aer_Phase_Fct"], self.profiles["aer_Phase_DoLP"] = np.array(self.profiles["aer_Phase_Fct"]), np.array(self.profiles["aer_Phase_DoLP"])

				# self.aerosol.PlotPhaseFunction(["1low", "2high"], [391.4, 557.7, 630.0])
				# plt.show()

			self.GetAerosolProfil(in_dict)
			self.profiles["beta_aer"] = np.array([self.AerosolCS(a) for a in self.profiles["HGT"]])


		self.profiles['total_absorption'] = np.zeros_like(self.profiles['HGT'])
		self.profiles['total_absorption'] += self.profiles["tau_ray"]
		if self.use_ozone:
			self.profiles['total_absorption'] += self.tau_O3_list
		if self.use_aerosol:
			self.profiles['total_absorption'] += self.tau_aer_list



		if mpi_rank == 0:
			print("tau_0", self.tau_0)
			print("DEBUG ATMOSPHERE TOTAL OD (0-120km):", self.GetO3Absorbtion(0, 120), self.GetRSOpticalDepth(0, 120), self.GetAerosolsAbsorbtion(0, 120), self.GetO3Absorbtion(0, 120)+ self.GetRSOpticalDepth(0, 120)+ self.GetAerosolsAbsorbtion(0, 120))

		# self.depola = 0
		# if self.wavelength == 630:
		# 	self.depola = 2.787e-2
		# elif self.wavelength in [557.7, 550, 500]:
		# 	self.depola = 2.843e-2
		# elif self.wavelength == 427.8:
		# 	self.depola = 2.923e-2
		# elif self.wavelength == 391.4:
		# 	self.depola = 2.954e-2


		if mpi_rank == 0 and self.show_cross_section:
			self.MakeCrossSectionPlot()
			plt.show()
		if mpi_rank == 0 and self.show_atmosphere:
			self.MakePlots()
			plt.show()
		if mpi_rank == 0 and self.show_aer_size and self.aer_complexity > 0:
			self.MakeAerSizePlot()
			plt.show()



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
						self.profiles[name].append(float(word.strip(",; ")) * unit_factor) #Save values in the dictionnary

		if "win.atm" in filename:
			for io3, o3 in enumerate(self.profiles["O3"]):
				self.profiles["O3"][io3] = o3 * 2

		self.profiles["DEN"] = 1e9 * np.array(self.profiles["PRE"]) / np.array(self.profiles["TEM"]) / cst.value("Boltzmann constant") #in molecules/km-3

		self.nb_profile_alt = len(self.profiles["HGT"]) #Number of points in the profiles
		self.profiles["delta_z"] = [self.profiles["HGT"][i+1] - self.profiles["HGT"][i] for i in range(self.nb_profile_alt-1)]
		self.profiles["delta_z"].append(self.profiles["delta_z"][-1])

		for k, v in self.profiles.items():
			self.profiles[k] = np.array(v)

	@GlobalTimer
	def GetProfileValue(self, alt, name):
		"""From any altitude in km, gives the profile value with linear approx."""
		# index = np.searchsorted(self.profiles["HGT"], alt)
		# return self.profiles[name][index]
		return np.interp(alt, self.profiles["HGT"], self.profiles[name])
		# return FastInterp(alt, self.profiles["HGT"], self.profiles[name])

		#
		# if alt in self.profiles["HGT"]:
		# 	i = self.profiles["HGT"].index(alt)
		# 	return self.profiles[name][i]
		# else:
		# 	i = 0
		# 	while i < self.nb_profile_alt and alt > self.profiles["HGT"][i]:
		# 		i += 1
		#
		# 	dalt = 0
		# 	dprof = 0
		# 	d = 0
		# 	value = 0
		# 	ip, im = i, i - 1
		# 	if ip < self.nb_profile_alt:
		# 		dalt = self.profiles["HGT"][ip] - self.profiles["HGT"][im]
		# 		dprof = self.profiles[name][ip] - self.profiles[name][im]
		# 		d = dprof / dalt
		# 		value = self.profiles[name][im]
		#
		# 	return value + d * (alt - self.profiles["HGT"][im])

	def SetObservation(self, obs, ouv_pc):

		self.instrument_altitude = obs.A_alt

		if self.h_r_min <= self.instrument_altitude:
			self.range_min = 0.
		else:
			self.range_min = GetArbitraryDistance(self.instrument_altitude + RT, self.h_r_min + RT, obs.e)

		if self.h_r_max > self.instrument_altitude:
			self.range_max = GetArbitraryDistance(self.instrument_altitude + RT, self.h_r_max + RT, obs.e)
		else:
			raise Exception("Maximum atmospheric height < instrument altitude. The geometry can not handle it.")

		#numbers of scattering point along the line of sight
		self.los_length = self.range_max - self.range_min

		# if self.Nlos > 1:
		#array of range for each scattering points (distance along los from the instrument)

		### Linear dlos (all bins along line of sight of equal length)
		self.range_list = np.linspace(self.range_min, self.range_max, self.Nlos + 1) #in km

		### Quadratic dlos : Along the line of sight, bins close to the instrument smaller than bins away
		# self.range_list = np.array([x**2 for x in np.linspace(np.sqrt(self.range_min), np.sqrt(self.range_max), self.Nlos + 1)]) #in km


		self.mid_range_list = np.array([(self.range_list[i+1] + self.range_list[i]) / 2. for i in range(0, self.Nlos)]) #in km
		self.d_los_list = np.array([self.range_list[i+1] - self.range_list[i] for i in range(0, self.Nlos)]) #in km

		self.mid_altitudes_list = [GetArbitraryAltitude(RT + obs.A_alt, ab, obs.e) for ab in self.mid_range_list]

		self.volumes = np.array([self.GetVolume(r, ouv_pc, dr = dr, unit="km", type = "cone") for r, dr in zip(self.mid_range_list, self.d_los_list)])

		delta_z = np.abs(self.instrument_altitude - self.mid_altitudes_list)

		if self.use_ozone:
			los_O3_abs = np.array([self.GetO3Absorbtion(self.instrument_altitude, h)  for h in self.mid_altitudes_list])
			self.los_O3_transmittance = np.exp(- los_O3_abs * self.mid_range_list / delta_z)
		else:
			self.los_O3_transmittance = 1

		los_RS_abs = np.array([self.GetRSOpticalDepth(self.instrument_altitude, h)  for h in self.mid_altitudes_list])
		self.los_RS_transmittance = np.exp(- los_RS_abs * self.mid_range_list / delta_z)

		if self.use_aerosol:
			los_aer_abs = np.array([self.GetAerosolsAbsorbtion(self.instrument_altitude, h) for h in self.mid_altitudes_list])
			self.los_aer_transmittance = np.exp(- los_aer_abs * self.mid_range_list / delta_z)
		else:
			self.los_aer_transmittance = 1

		self.los_transmittance = self.los_O3_transmittance * self.los_RS_transmittance * self.los_aer_transmittance

		# for it, t in enumerate(self.los_transmittance):
		# 	print(los_RS_abs[it], los_O3_abs[it], los_aer_abs[it])

	# #@timer
	@GlobalTimer
	def GetVolume(self, AR, ouv_pc, da = None, dr=None, index = None , unit="km", type="cone"):
		"""Get the volume of a truncated cone of length d_los, and half opening angle ouv_pc.
		Unit defines the unit you want it in. can be "cm", "m" or "km". (the cube is implicit)"""

		if dr is None and index is None:
			raise ValueError("Incorrect arguments passed to atmosphere.GetVolume(). Must give at least one of dr or index parameter")

		if type == "index":
			V = self.volumes[index]

		elif type == "cone":
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

	def GetRSVolumeCS(self, alt):
		return self.GetProfileValue(alt, "beta_ray")

	#@timer
	def RSVolumeCS(self, alt):
		"""Return the Rayleigh scattering volume cross section as calculated in Bucholtz 95 in km-1.
		"""

		# beta_s = self.GetSquareLawFit(wl, "Volume CS") #in km-1
		P = self.GetProfileValue(alt, "PRE") #in Pa
		T = self.GetProfileValue(alt, "TEM") #in K

		# print("alt, P, T", alt, P, T)

		beta_z = self.beta_0 * (P / self.P0) * (self.T0 / T) # in km-1

		# if mpi_rank == 0:
		# 	print("beta_z", alt, beta_z)

		return beta_z  # in km-1

	def GetRSAbsorptionCoeff(self, alt):
		P = self.GetProfileValue(alt, "PRE")
		tau = self.tau_0 * (P / self.P0)
		return tau

	#@timer
	def GetRSOpticalDepth(self, E_alt, R_alt):
		"""Return VERTICAL optical depth between two altitudes as calculated in Bucholtz 95."""
		# tau_E = self.GetSquareLawFit(wl, "Optical Depth")
		# P_E = self.GetProfileValue(E_alt, "PRE")
		# tau_E = self.tau_0 * (P_E / self.P0)
		tau_E = self.GetRSAbsorptionCoeff(E_alt)

		# tau_R = self.GetSquareLawFit(wl, "Optical Depth")
		# P_R = self.GetProfileValue(R_alt, "PRE")
		# tau_R = self.tau_0 * (P_R / self.P0)
		tau_R = self.GetRSAbsorptionCoeff(R_alt)

		tau_ER = abs(tau_E - tau_R)
		return tau_ER

	def GetSquareLawFit(self, wl, purpose):
		"""Return the result of the square law fit function presented in Bucholtz 1995. Purpose can be:
		"Single CS" = For single air particule rayleigh cross section in cm2.
		"Volume CS" = For volume cross rayleigh section in km-1.
		"Optical Depth" = For the rayleigh optical depth.
		Value of parameter A for optical depth is given for a Sub arctic winter atmosphere model."""

		A, B, C, D = self.GetSquareLawParam(wl, purpose)
		wl /= 1000.  #wavelength given in nanometers, converted to micrometer
		E = - (B + C * wl + D / wl)
		return A * wl ** E

	def GetRSPhaseFunction(self, wl, theta):
		return np.interp(theta, self.profiles["sca_angle"], self.profiles["ray_Phase_Fct"])

	def RSPhaseFunctionDoLP(self, RD_angle, f = 1):
		return np.sin(RD_angle)**2 / (f + np.cos(RD_angle)**2)

	def GetRSPhaseFunctionDoLP(self, RD_angle):
		return np.interp(RD_angle, self.profiles["sca_angle"], self.profiles["ray_Phase_DoLP"])

	#@timer
	def RSPhaseFunction(self, wl, theta):
		### Simple approx
		A = 3. / 4
		B = 1 + np.cos(theta) ** 2


		### Chandrasekhar formula
		# if   wl == 391.4:
		# 	gamma = 1.499 * 0.01
		# elif wl == 427.8:
		# 	gamma = 1.483 * 0.01
		# elif wl in [557.7, 550, 500]:
		# 	gamma = 1.442 * 0.01
		# elif wl == 630:
		# 	gamma = 1.413 * 0.01
		# A = 3 / (4 + 8 * gamma)
		# B = (1 + 3 * gamma) + (1 - gamma) * np.cos(theta) ** 2

		# return 1
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

	def SetO3CS(self):
		self.O3CS = np.interp(self.wavelength, self.O3CS_list["wavelength"], self.O3CS_list["cross_section"]) #in cm2/molecules
		self.O3CS *= 1e-10 #in km2 / molecules

		# self.tauO3_0 = self.GetO3Absorbtion(self.h_r_min, self.h_r_max)

		# print("self.tauO3_0", self.tauO3_0, np.exp(-self.tauO3_0))

	def LoadO3AbsorptionCS(self):
		self.O3CS_list = pd.read_fwf("seO3_Burrows273-air-1.data", names=["wavelength", "cross_section"], skiprows=8) #in cm2/molecules

		self.SetO3CS() #in km2/molecules

		d = lambda i: self.profiles["O3"][i] * self.profiles["DEN"][i] * self.profiles["delta_z"][i] #[%] * [km-3] * [km]
		a = lambda d: np.sum(d) * self.O3CS #* (z_max - z_min)
		tau = lambda z_min, z_max: a([d(ih) for ih, h in enumerate(self.profiles["HGT"]) if z_min <= h <= z_max])

		self.tau_O3_list = np.array([tau(0, h) for h in self.profiles["HGT"]])

		self.tauO3_0 = self.GetO3Absorbtion(self.h_r_min, self.h_r_max)

		# # self.GetO3Absorbtion(0, 120)
		# # print("O3 column density (km-2):", sum(densities))
		#
		# x = np.array([h for h in self.profiles["HGT"] if h <= 60])
		# e = np.pi/4 #np.pi/4
		# yO3 = np.array([np.exp(-self.GetO3Absorbtion(0, h) / np.sin(e)) for h in x])
		# yRS = np.array([np.exp(-self.GetRSOpticalDepth(0, h) / np.sin(e)) for h in x])
		# ytot = yO3 * yRS
		# O3 = np.array(self.profiles["O3"][:len(x)]) * np.array(self.profiles["DEN"][:len(x)])
		#
		# if mpi_rank == 0:
		# 	fig, ax1 = plt.subplots()
		# 	ax1.plot(x, yO3, "b", label="Tr O3")
		# 	ax1.plot(x, yRS, "r", label="Tr Rayleigh")
		# 	ax1.plot(x, ytot, "k", label="Tr Total")
		# 	ax11 = ax1.twinx()
		# 	ax11.plot(x, O3 * 1e-15, label="O3 nb density (cm-2)")
		# 	ax1.set_xlabel("Altitude  (km)")
		# 	ax1.set_ylabel("Transmittance")
		# 	ax11.set_ylabel("Nb Density (cm-2)")
		# 	ax1.legend()
		# 	ax11.legend(loc="center right")
		# 	plt.title("Transmittance at 557.7nm, 45")
		# 	plt.show()
		# mpi_comm.Barrier()

	#@timer
	def GetO3Absorbtion(self, z_min, z_max):
		"""Returns the vertical absorption due to O3 between altitude zmin and zmax."""
		# if z_min > z_max:
		# 	z_min, z_max = z_max, z_min

		if self.use_ozone:
			absorbtion = np.interp(z_max, self.profiles["HGT"], self.tau_O3_list) - np.interp(z_min, self.profiles["HGT"], self.tau_O3_list)
			absorbtion = abs(absorbtion)
		else:
			absorbtion = 0
		# # d = lambda i: self.profiles["O3"][i] * self.profiles["PRE"][i] / cst.Boltzmann / self.profiles["TEM"][i]
		# d = lambda i: self.profiles["O3"][i] * self.profiles["DEN"][i] * self.profiles["delta_z"][i] #[%] * [km-3] * [km]
		#
		# densities = [d(ih) for ih, h in enumerate(self.profiles["HGT"]) if z_min <= h <= z_max]
		#
		# # print("O3 column density (km-2):", z_min, z_max, sum(densities))
		#
		# absorbtion = np.sum(densities) * self.O3CS #* (z_max - z_min)

		return absorbtion

	def MakeCrossSectionPlot(self):

		f, ax = plt.subplots(1)
		wavelengths = np.linspace(350, 700, 1000)
		beta_0 = [self.GetSquareLawFit(wl, "Volume CS") for wl in wavelengths] #in km-1

		C_O3 = np.interp(wavelengths, self.O3CS_list["wavelength"], self.O3CS_list["cross_section"]) #cm2
		# C_O3 =

		ax.plot(wavelengths, beta_0, "-k")
		ax1 = ax.twinx()
		ax1.plot(wavelengths, C_O3, "--k")

		ax.set_xlabel("Wavelength (nm)")
		ax.set_ylabel(r"Rayleigh volume scattering coefficient $\beta$ (km$^{-1}$)")
		ax1.set_ylabel(r"Ozone absorption cross section $\sigma$ (cm$^{2}$)")


	def MakeAerSizePlot(self):
		f1, axs = plt.subplots(1, sharey = True, figsize=(16, 8))

		r_list = np.logspace(np.log10(self.aerosol.r_min), np.log10(self.aerosol.r_max), num=100)
		# rn = rn0*exp(-(log(r/rm))**2/(2*(log(sig))**2))/r/log(sig)/(2*pi)**0.5 pour r<rmax
		distrib = self.aerosol.rn0 * np.exp(-(np.log(r_list / self.aerosol.rmg))**2 / (2 * self.aerosol.ln_sigma**2)) / r_list / self.aerosol.ln_sigma / np.sqrt(2 * np.pi)

		axs.plot(r_list, distrib)
		axs.loglog()

	def MakePlots(self):
		f1, axs = plt.subplots(3, sharey = True, figsize=(16, 8))
		axs[0] = plt.subplot(131)
		axs[1] = plt.subplot(132)
		axs[2] = plt.subplot(133)
		# axs[3] = plt.subplot(144)


		axs[0].plot(self.profiles["TEM"], self.profiles["HGT"])

		axs[1].plot(self.profiles["PRE"], self.profiles["HGT"])
		axs[1].set_xscale('log')
		# ax3.plot(self.profiles["N2"], self.profiles["HGT"], label="N2")
		# axs[2].plot(np.array(self.profiles["O3"]) * self.profiles["DEN"] * 1e-15, self.profiles["HGT"], label="O3")
		axs[2].plot(self.tau_aer_list, self.profiles["HGT"], label="O3")
		# axs[2].plot(np.array(self.profiles["AER"]), self.profiles["HGT"], label="O3")
		# axs[2].set_xscale('log')
		#
		# axs[0].set_yscale('log')
		# axs[1].set_yscale('log')
		# axs[2].set_yscale('log')


		# ax3.plot(self.profiles["H2O"], self.profiles["HGT"], label="H2O")
		# ax3.plot(self.profiles["CO2"], self.profiles["HGT"], label="CO2")
		# ax3.plot(self.profiles["O3"], self.profiles["HGT"], label="O3")
		#
		# ax3.legend()
		# axs[2].set_xscale('log')
		for profil_name in ["urb", "mar", "con"]:
			exctinction_profile = []
			if profil_name in ["urban", "urb"]:
				profil_lim_alt 	= [2, 12, 20, 30] # in km: Upper limits of layers of different exctinction coefficients.
				profil_ext		= np.array([0.5, 0.0025, 2.18e-4, 3.32e-5, 0]) # in km-1: Exctinction coefficient of the different altitude layers. index 0 correspond to the layer below altitude profil_lim_alt[0], index 1 correspond to the layer between altitude profil_lim_alt[0] and altitude profil_lim_alt[1],... index -1 correspond to the layer above altitude profil_lim_alt[-1]
			elif profil_name in ["continental", "con"]:
				profil_lim_alt 	= [2, 12, 20, 30] # in km: Upper limits of layers of different exctinction coefficients.
				profil_ext		= [0.1, 0.0025, 2.18e-4, 3.32e-5, 0] # in km-1: Exctinction coefficient of the different altitude layers. index 0 correspond to the layer below altitude profil_lim_alt[0], index 1 correspond to the layer between altitude profil_lim_alt[0] and altitude profil_lim_alt[1],... index -1 correspond to the layer above altitude profil_lim_alt[-1]
			elif profil_name in ["maritime", "mar"]:
				profil_lim_alt 	= [2, 12, 20, 30] # in km: Upper limits of layers of different exctinction coefficients.
				profil_ext		= [0.025, 0.0025, 2.18e-4, 3.32e-5, 0] # in km-1: Exctinction coefficient of the different altitude layers. index 0 correspond to the layer below altitude profil_lim_alt[0], index 1 correspond to the layer between altitude profil_lim_alt[0] and altitude profil_lim_alt[1],... index -1 correspond to the layer above altitude profil_lim_alt[-1]
			elif profil_name in ["maritime", "test1"]:
				profil_lim_alt 	= [2, 12, 20, 30] # in km: Upper limits of layers of different exctinction coefficients.
				profil_ext		= [0.05, 0.0025, 2.18e-4, 3.32e-5, 0] # in km-1: Exctinction coefficient of the different altitude layers. index 0 correspond to the layer below altitude profil_lim_alt[0], index 1 correspond to the layer between altitude profil_lim_alt[0] and altitude profil_lim_alt[1],... index -1 correspond to the layer above altitude profil_lim_alt[-1]
			elif profil_name in ["maritime", "test2"]:
				profil_lim_alt 	= [2, 12, 20, 30] # in km: Upper limits of layers of different exctinction coefficients.
				profil_ext		= [0.25, 0.0025, 2.18e-4, 3.32e-5, 0] # in km-1: Exctinction coefficient of the different altitude layers. index 0 correspond to the layer below altitude profil_lim_alt[0], index 1 correspond to the layer between altitude profil_lim_alt[0] and altitude profil_lim_alt[1],... index -1 correspond to the layer above altitude profil_lim_alt[-1]
			elif profil_name in ["maritime", "test3"]:
				profil_lim_alt 	= [2, 12, 20, 30] # in km: Upper limits of layers of different exctinction coefficients.
				profil_ext		= [0.001, 0.0025, 2.18e-4, 3.32e-5, 0] # in km-1: Exctinction coefficient of the different altitude layers. index 0 correspond to the layer below altitude profil_lim_alt[0], index 1 correspond to the layer between altitude profil_lim_alt[0] and altitude profil_lim_alt[1],... index -1 correspond to the layer above altitude profil_lim_alt[-1]
			elif profil_name in ["maritime", "test4"]:
				profil_lim_alt 	= [2, 12, 20, 30] # in km: Upper limits of layers of different exctinction coefficients.
				profil_ext		= [1, 0.0025, 2.18e-4, 3.32e-5, 0] # in km-1: Exctinction coefficient of the different altitude layers. index 0 correspond to the layer below altitude profil_lim_alt[0], index 1 correspond to the layer between altitude profil_lim_alt[0] and altitude profil_lim_alt[1],... index -1 correspond to the layer above altitude profil_lim_alt[-1]
			elif profil_name in ["maritime", "test5"]:
				profil_lim_alt 	= [2, 12, 20, 30] # in km: Upper limits of layers of different exctinction coefficients.
				profil_ext		= [0., 0.0025, 2.18e-4, 3.32e-5, 0] # in km-1: Exctinction coefficient of the different altitude layers. index 0 correspond to the layer below altitude profil_lim_alt[0], index 1 correspond to the layer between altitude profil_lim_alt[0] and altitude profil_lim_alt[1],... index -1 correspond to the layer above altitude profil_lim_alt[-1]

			for ih, h in enumerate(self.profiles["HGT"]):
				index = np.searchsorted(profil_lim_alt, h)
				exctinction_profile.append(profil_ext[index])

			exctinction_profile = np.array(exctinction_profile)
			# axs[3].plot(exctinction_profile, self.profiles["HGT"])

		f1.subplots_adjust(wspace=0)

		axs[1].set(yticklabels=[])
		axs[2].set(yticklabels=[])


		# axs[0].set_title("Temperature")
		# axs[1].set_title("Pressure")
		# axs[2].set_title("O3")

		axs[0].set_ylabel("Altitude (km)")
		axs[0].set_xlabel("Temperature (K)")
		axs[1].set_xlabel("Pressure (Pa)")
		axs[2].set_xlabel(r'O3 (cm$^{-3}$)')
		# axs[3].set_xlabel("Caer (km{-1})")
		# axs[3].semilogx()

	def GetAerosolProfil(self, in_dict):
		self.profiles["AER"] = np.zeros(len(self.profiles["HGT"]))
		self.aer_model = in_dict["aer_model"]

		if self.aer_complexity > 0:
			# Hn = [1, 2, 5, 10] #km
			# n0 = [3500, 3500, 3500, 3500] #cm-3
			# nB = [300, 300, 300, 300] #cm-3
			#
			# for Hn, n0, nB in zip(Hn, n0, nB):
			# 	p = n0 * (np.exp(- np.array(self.profiles["HGT"]) / Hn) + (nB/n0)) * (1e5)**3
			#
			# plt.loglog()
			# plt.legend()
			# plt.show()

			Hn = self.aerosol.Hn #float(in_dict['aer_Hn'])  #0.7 #km
			n0 = self.aerosol.n0 #float(in_dict['aer_n0']) #cm-3
			nB = self.aerosol.nB #float(in_dict['aer_nBK']) #cm-3
			v = Hn / abs(Hn)

			p = lambda z: n0 * (np.exp(-z / Hn) + (nB/n0)**v)**v * (1e5)**3

			self.aerosol_profil_param = [Hn, n0, nB]

			self.profiles["AER"] = np.zeros_like(self.profiles["HGT"])
			for i, alt in enumerate(self.profiles["HGT"]):
				if alt <= self.aerosol.max_alt:
					self.profiles["AER"][i] = p(alt) #in km-3
				else:
					break


			# f = plt.figure()
			# plt.loglog(self.profiles["AER"], self.profiles["HGT"], label = f"{Hn} {n0} {nB}")
			# plt.show()

			# print("AEROSOLS", self.profiles["AER"])

			# self.aerosol_mean_size = float(in_dict['aer_radius']) #in nanometers
			# self.aerosol_Qext = float(in_dict['aer_Qext']) #* (1 + (2 * np.pi * self.aerosol_mean_size / self.wavelength) ** (-2/3))
			# self.aerosol_Qabs = float(in_dict['aer_Qabs'])
			# self.aerosol_Csca = (self.aerosol_Qext - self.aerosol_Qabs) * np.pi * (self.aerosol_mean_size * 1e-9)**2 #in m2
			# self.aerosol_Csca *= 1e-6 #in km2
			# self.aerosol_Cext = self.aerosol_Qext * np.pi * (self.aerosol_mean_size * 1e-9)**2 #in m2
			# self.aerosol_Cext *= 1e-6 #in km2


			d = lambda i: self.profiles["AER"][i] * self.profiles["delta_z"][i] #[km-3] * [km]
			a = lambda d: np.sum(d) * self.aerosol.c_ext
			tau = lambda z_min, z_max: a([d(ih) for ih, h in enumerate(self.profiles["HGT"]) if z_min <= h <= z_max])

		elif self.aer_complexity == 0:
			self.aerosol_profile_name = in_dict['aer_ext_profile_name'].lower()

			print(f"Using simple WMO aerosol profile: {self.aerosol_profile_name}")
			if self.aerosol_profile_name in ["urban", "urb"]:
				profil_lim_alt 	= [0, 2, 12, 20, 30] # in km: lower limits of layers of different exctinction coefficients.
				profil_ext		= np.array([0.5, 0.0025, 2.18e-4, 3.32e-5, 0]) # in km-1: Exctinction coefficient of the different altitude layers. index 0 correspond to the layer below altitude profil_lim_alt[0], index 1 correspond to the layer between altitude profil_lim_alt[0] and altitude profil_lim_alt[1],... index -1 correspond to the layer above altitude profil_lim_alt[-1]
			elif self.aerosol_profile_name in ["continental", "con"]:
				profil_lim_alt 	= [0, 2, 12, 20, 30] # in km: lower limits of layers of different exctinction coefficients.
				profil_ext		= [0.1, 0.0025, 2.18e-4, 3.32e-5, 0] # in km-1: Exctinction coefficient of the different altitude layers. index 0 correspond to the layer below altitude profil_lim_alt[0], index 1 correspond to the layer between altitude profil_lim_alt[0] and altitude profil_lim_alt[1],... index -1 correspond to the layer above altitude profil_lim_alt[-1]
			elif self.aerosol_profile_name in ["maritime", "mar"]:
				profil_lim_alt 	= [0, 2, 12, 20, 30] # in km: lower limits of layers of different exctinction coefficients.
				profil_ext		= [0.025, 0.0025, 2.18e-4, 3.32e-5, 0] # in km-1: Exctinction coefficient of the different altitude layers. index 0 correspond to the layer below altitude profil_lim_alt[0], index 1 correspond to the layer between altitude profil_lim_alt[0] and altitude profil_lim_alt[1],... index -1 correspond to the layer above altitude profil_lim_alt[-1]
			elif self.aerosol_profile_name in ["pust"]:
				profil_lim_alt 	= [0, 5.5, 120] # in km: lower limits of layers of different exctinction coefficients.
				profil_ext		= [0.01, 0.005, 0] # in km-1: Exctinction coefficient of the different altitude layers. index 0 correspond to the layer below altitude profil_lim_alt[0], index 1 correspond to the layer between altitude profil_lim_alt[0] and altitude profil_lim_alt[1],... index -1 correspond to the layer above altitude profil_lim_alt[-1]
			elif self.aerosol_profile_name in ["pust2"]:
				profil_lim_alt 	= [0, 10] # in km: lower limits of layers of different exctinction coefficients.
				profil_ext		= [0.015, 0] # in km-1: Exctinction coefficient of the different altitude layers. index 0 correspond to the layer below altitude profil_lim_alt[0], index 1 correspond to the layer between altitude profil_lim_alt[0] and altitude profil_lim_alt[1],... index -1 correspond to the layer above altitude profil_lim_alt[-1]
			elif self.aerosol_profile_name in ["test1"]:
				profil_lim_alt 	= [0, 2, 12, 20, 30] # in km: lower limits of layers of different exctinction coefficients.
				profil_ext		= [0.05, 0.0025, 2.18e-4, 3.32e-5, 0] # in km-1: Exctinction coefficient of the different altitude layers. index 0 correspond to the layer below altitude profil_lim_alt[0], index 1 correspond to the layer between altitude profil_lim_alt[0] and altitude profil_lim_alt[1],... index -1 correspond to the layer above altitude profil_lim_alt[-1]
			elif self.aerosol_profile_name in ["test2"]:
				profil_lim_alt 	= [0, 2, 12, 20, 30] # in km: lower limits of layers of different exctinction coefficients.
				profil_ext		= [0.25, 0.0025, 2.18e-4, 3.32e-5, 0] # in km-1: Exctinction coefficient of the different altitude layers. index 0 correspond to the layer below altitude profil_lim_alt[0], index 1 correspond to the layer between altitude profil_lim_alt[0] and altitude profil_lim_alt[1],... index -1 correspond to the layer above altitude profil_lim_alt[-1]
			elif self.aerosol_profile_name in ["test3"]:
				profil_lim_alt 	= [0, 2, 12, 20, 30] # in km: lower limits of layers of different exctinction coefficients.
				profil_ext		= [0.001, 0.0025, 2.18e-4, 3.32e-5, 0] # in km-1: Exctinction coefficient of the different altitude layers. index 0 correspond to the layer below altitude profil_lim_alt[0], index 1 correspond to the layer between altitude profil_lim_alt[0] and altitude profil_lim_alt[1],... index -1 correspond to the layer above altitude profil_lim_alt[-1]
			elif self.aerosol_profile_name in ["test4"]:
				profil_lim_alt 	= [0, 2, 12, 20, 30] # in km: lower limits of layers of different exctinction coefficients.
				profil_ext		= [1, 0.0025, 2.18e-4, 3.32e-5, 0] # in km-1: Exctinction coefficient of the different altitude layers. index 0 correspond to the layer below altitude profil_lim_alt[0], index 1 correspond to the layer between altitude profil_lim_alt[0] and altitude profil_lim_alt[1],... index -1 correspond to the layer above altitude profil_lim_alt[-1]
			elif self.aerosol_profile_name in ["test5"]:
				profil_lim_alt 	= [0, 2, 12, 20, 30] # in km: lower limits of layers of different exctinction coefficients.
				profil_ext		= [0., 0.0025, 2.18e-4, 3.32e-5, 0] # in km-1: Exctinction coefficient of the different altitude layers. index 0 correspond to the layer below altitude profil_lim_alt[0], index 1 correspond to the layer between altitude profil_lim_alt[0] and altitude profil_lim_alt[1],... index -1 correspond to the layer above altitude profil_lim_alt[-1]

			# self.exctinction_profile = []
			# for ih, h in enumerate(self.profiles["HGT"]):
			# 	index = np.searchsorted(profil_lim_alt, h)
			# 	self.exctinction_profile.append(profil_ext[index-1])
			# self.exctinction_profile = np.array(self.exctinction_profile)
			self.exctinction_profile = np.interp(self.profiles["HGT"], profil_lim_alt, profil_ext)

			self.single_scatter_albedo = float(in_dict["aer_single_scatter_albedo"])
			self.aer_phase_fct_asym_g = float(in_dict["aer_phase_fct_asym_g"])

			d = lambda i: self.exctinction_profile[i] * self.profiles["delta_z"][i] #[km-1] * [km]
			tau = lambda z_min, z_max: np.sum([d(ih) for ih, h in enumerate(self.profiles["HGT"]) if z_min <= h <= z_max])

		self.tau_aer_list = np.array([tau(0, h) for h in self.profiles["HGT"]])

		if mpi_rank == 0: print("DEBUG AEROSOL Optical Depth (0-120km)", self.GetAerosolsAbsorbtion(0, 120))

		if "total_OD" in in_dict:
			self.total_OD_goal = float(in_dict["total_OD"])
			if self.total_OD_goal != 0:
				correction = self.total_OD_goal / self.GetAerosolsAbsorbtion(0, 120)
				self.tau_aer_list *= correction
				if self.aer_complexity > 0:
					self.profiles["AER"] *= correction

		if mpi_rank == 0: print("DEBUG AEROSOL Optical Depth (0-120km)", self.GetAerosolsAbsorbtion(0, 120))

		# f, ax = plt.subplots(1)
		# # ax.plot(self.exctinction_profile, self.profiles["HGT"])
		# # ax.plot(self.better_exctinction_profile, self.profiles["HGT"])
		# # ax1 = ax.twiny()
		# ax.plot(self.tau_aer_list, self.profiles["HGT"], "r")
		# plt.show()

	def GetAerosolCS(self, alt):
		return self.GetProfileValue(alt, "beta_aer")

	def AerosolCS(self, z):
		if self.use_aerosol:
			if self.aer_complexity > 0:
				aer_density = self.GetProfileValue(z, "AER") #km-3
				cross_section = self.aerosol.c_sca * aer_density #[km2 * km-3] = [km-1]
			else:
				cross_section = self.single_scatter_albedo * np.interp(z, self.profiles["HGT"], self.exctinction_profile)
		else:
			cross_section = 0

		return cross_section

	def GetAerosolsAbsorbtion(self, z_min, z_max):
		"""Returns the vertical absorption due to aerosols between altitude zmin and zmax. """

		if self.use_aerosol:
			absorbtion = np.interp(z_max, self.profiles["HGT"], self.tau_aer_list) - np.interp(z_min, self.profiles["HGT"], self.tau_aer_list)
			absorbtion = abs(absorbtion)
		else:
			absorbtion = 0
		return absorbtion


	def GetAerosolPhaseFunctionDoLP(self, a):
		dolp = np.interp(a, self.profiles["sca_angle"], self.profiles["aer_Phase_DoLP"])
		return dolp

	def GetAerosolPhaseFunction(self, a):
		pf   = np.interp(a, self.profiles["sca_angle"], self.profiles["aer_Phase_Fct"])
		dolp = np.interp(a, self.profiles["sca_angle"], self.profiles["aer_Phase_DoLP"])
		return pf, dolp


	def AerosolPhaseFunction(self, a):
		if self.aer_complexity == 0:
			g = self.aer_phase_fct_asym_g # g > O == front scattering predominant, g < 0 == backscattering predominant
			###################################################################
			### Here, a = 0 is FRONT-scattering and a=pi is BACK-scattering ###
			###################################################################
			Ptot = (1 - g**2) / (1 + g**2 - 2 * g * np.cos(a))**(3/2) #Heney-Greenstein Function
			DoLP = 0
			AoLP = 0

		else:
			Ptot = np.interp(a, self.aerosol.sca_data["theta"], self.aerosol.sca_data["i"])
			# For Mie scattering, a DoLP < 0 = AoLP is parallel to scatering plane ! So 90 degrees from Rayleigh scattering.
			DoLP = np.interp(a, self.aerosol.sca_data["theta"], self.aerosol.sca_data["d"])

		return Ptot, DoLP


class Aerosol:
	def __init__(self, in_dict, run_mie=True):

		self.in_dict = in_dict

		self.path = in_dict["src_path"] + "mie/"

		self.model = in_dict["aer_model"].lower()
		self.name = in_dict["aer_name"].lower()

		### Size distribution fct is:
		### rn = self.rn0 * np.exp(-(np.log(r_aer/self.rmg))**2 / (2 * self.ln_sigma**2)) / r_aer / self.ln_sigma / (2 * np.pi)**0.5
		###
		self.wavelength = float(in_dict["wavelength"])

		if self.model == "manual":
			self.index_re = float(in_dict["aer_nr"])
			self.index_im = float(in_dict["aer_ni"])
			self.rn0 = float(in_dict["aer_rn0"])
			self.rmg = float(in_dict["aer_rmg"])
			self.ln_sigma = float(in_dict["aer_ln_sigma"])
			self.r_min = float(in_dict["aer_r_min"])
			self.r_max = float(in_dict["aer_r_max"])

			self.Hn = float(in_dict['aer_Hn'])  #km
			self.n0 = float(in_dict['aer_n0']) #cm-3
			self.nB = float(in_dict['aer_nBK']) #cm-3

			self.max_alt = float(in_dict['aer_max_alt']) #km

		elif self.model != "default":
			self.SetModel()
			self.name = self.GetFileName()

		else:
			self.name = self.GetFileName()
			#Colette example values
			self.index_re = 1.43
			self.index_im = -2.0E-03 #<0
			self.rn0 = 1
			self.rmg = 0.1
			self.ln_sigma = 0.5
			self.r_min = 0.001
			self.r_max = 50

			self.Hn = 0#km
			self.n0 = 1 #cm-3
			self.nB = 0#cm-3

			self.max_alt = 12 #km


		self.Compute()


	def Compute(self):
		input_file_exists = os.path.isfile(self.path + self.name)
		result_file_exists = os.path.isfile(self.path + "res_" + self.name)
		if not result_file_exists or not input_file_exists:
			# if not self.standard:
			mpi_comm.Barrier()
			if mpi_rank == 0:
			# for wl in [620, 557.7, 427.8, 413, 391.4]:
				self.WriteInputFile()
				self.RunMieCode()
			mpi_comm.Barrier()

		self.LoadOpticalProperties()

	def GetFileName(self, model = None, wl = None):
		if model is None:
			model = self.model
		if wl is None:
			wl = self.wavelength

		return model + "_" + str(wl)


	def SetModel(self):
		self.standard = True

		if mpi_rank == 0:
			print("Aerosol Model:", self.model)

		if self.model in ["maritime", "1-low", "1low"]:
			self.model = "1low"

			self.index_re = 1.45
			self.index_im = -0.0035 #<0

			self.rn0 = 1
			self.rmg = 0.15 #0.133 #micrometers
			self.ln_sigma = 0.29
			self.r_min = 0.001
			self.r_max = 50

			self.Hn = 0.440  	#km
			self.n0 = 4000   	#cm-3
			self.nB = 10  		#cm-3
			self.max_alt = 12 #km

		elif self.model in ["mar1"]:
			self.index_re = 1.6
			self.index_im = -0.0035 #<0

			self.rn0 = 1
			self.rmg = 0.15 #0.133 #micrometers
			self.ln_sigma = 0.29
			self.r_min = 0.001
			self.r_max = 50

			self.Hn = 0.440  	#km
			self.n0 = 4000   	#cm-3
			self.nB = 10  		#cm-3
			self.max_alt = 12 #km
		elif self.model in ["mar2"]:
			self.index_re = 1.6
			self.index_im = -0.0035 #<0

			self.rn0 = 1
			self.rmg = 0.15 #0.133 #micrometers
			self.ln_sigma = 0.29
			self.r_min = 0.001
			self.r_max = 50

			self.Hn = 0.440  	#km
			self.n0 = 3000   	#cm-3
			self.nB = 10  		#cm-3
			self.max_alt = 12 #km

		elif self.model in ["mar3"]:
			self.index_re = 1.6
			self.index_im = -0.0035 #<0

			self.rn0 = 1
			self.rmg = 0.5 #0.133 #micrometers
			self.ln_sigma = 0.29
			self.r_min = 0.001
			self.r_max = 50

			self.Hn = 0.440  	#km
			self.n0 = 100   	#cm-3
			self.nB = 10  		#cm-3
			self.max_alt = 12 #km

		elif self.model in ["mar4"]:
			self.index_re = 1.6
			self.index_im = -0.005 #<0

			self.rn0 = 1
			self.rmg = 0.15 #0.133 #micrometers
			self.ln_sigma = 0.29
			self.r_min = 0.001
			self.r_max = 50

			self.Hn = 0.440  	#km
			self.n0 = 3000   	#cm-3
			self.nB = 10  		#cm-3
			self.max_alt = 12 #km
		elif self.model in ["mar5"]:
			self.index_re = 1.6
			self.index_im = -0.008 #<0

			self.rn0 = 1
			self.rmg = 0.15 #0.133 #micrometers
			self.ln_sigma = 0.29
			self.r_min = 0.001
			self.r_max = 50

			self.Hn = 0.440  	#km
			self.n0 = 3000   	#cm-3
			self.nB = 10  		#cm-3
			self.max_alt = 12 #km

		elif self.model in ["urban", "2-high", "2high"]:
			self.model = "2high"
			self.index_re = 1.61
			self.index_im = -0.03 #<0

			self.rn0 = 1
			self.rmg = 0.557 #micrometers
			self.ln_sigma = 0.266
			self.r_min = 0.001
			self.r_max = 50

			self.Hn = 0.50  	#km
			self.n0 = 1000   	#cm-3
			self.nB = 1  		#cm-3
			self.max_alt = 12 #km

		elif self.model in ["rural", "3-mid", "3mid"]:
			self.model = "3mid"
			self.index_re = 1.61
			self.index_im = -0.03 #<0

			self.rn0 = 1
			self.rmg = 0.557 #micrometers
			self.ln_sigma = 0.266
			self.r_min = 0.001
			self.r_max = 50

			self.Hn = 0.50  	#km
			self.n0 = 500   	#cm-3
			self.nB = 1  		#cm-3
			self.max_alt = 12 #km

		elif self.model in ["arctic"]:
			self.index_re = 1.45
			self.index_im = -0.03 #<0

			self.rn0 = 1
			self.rmg = 0.070 #micrometers
			self.ln_sigma = 0.245
			self.r_min = 0.001
			self.r_max = 50

			self.Hn = 100  	#km
			self.n0 = 200   	#cm-3
			self.nB = 200  		#cm-
			self.max_alt = 12 #km

		elif self.model in ["desert", "des"]:
			self.index_re = 1.53
			self.index_im = -0.008 #<0

			self.rn0 = 1
			self.rmg = 0.188 #micrometers
			self.ln_sigma = 0.77
			self.r_min = 0.001
			self.r_max = 50

			self.Hn = 2	 		#km
			self.n0 = 150   	#cm-3
			self.nB = 0  		#cm-3
			self.max_alt = 12 #km

		elif self.model in ["des2"]:
			self.index_re = 1.53
			self.index_im = -0.008 #<0

			self.rn0 = 1
			self.rmg = 0.188 #micrometers
			self.ln_sigma = 2
			self.r_min = 0.001
			self.r_max = 50

			self.Hn = 2	 		#km
			self.n0 = 150   	#cm-3
			self.nB = 0  		#cm-3
			self.max_alt = 12 #km
		elif self.model in ["des3"]:
			self.index_re = 1.53
			self.index_im = -0.1 #<0

			self.rn0 = 1
			self.rmg = 0.188 #micrometers
			self.ln_sigma = 2
			self.r_min = 0.001
			self.r_max = 50

			self.Hn = 2	 		#km
			self.n0 = 150   	#cm-3
			self.nB = 0  		#cm-3
			self.max_alt = 12 #km
		elif self.model in ["des4"]:
			self.index_re = 1.53
			self.index_im = -0.008 #<0

			self.rn0 = 1
			self.rmg = 0.188 #micrometers
			self.ln_sigma = 0.77
			self.r_min = 0.001
			self.r_max = 50

			self.Hn = 2	 		#km
			self.n0 = 1000   	#cm-3
			self.nB = 0  		#cm-3
			self.max_alt = 12 #km
		elif self.model in ["des5"]:
			self.index_re = 1.53
			self.index_im = -0.008 #<0

			self.rn0 = 1
			self.rmg = 0.188 #micrometers
			self.ln_sigma = 0.77
			self.r_min = 0.001
			self.r_max = 50

			self.Hn = 2	 		#km
			self.n0 = 500   	#cm-3
			self.nB = 0  		#cm-3
			self.max_alt = 12 #km
		else:
			self.standard = False

	def WriteInputFile(self):
		with open(self.path + self.name, "w") as f:
			f.write(f"{self.wavelength / 1000.} {self.index_re} {self.index_im} ! lambda (micrometre), indice reel nr, complexe ni (<0)\n")
			f.write(f"{self.rn0} {self.rmg} {self.ln_sigma} {self.r_min} {self.r_max} ! No,rmg,ln(sigmag),Rmin,Rmax: in micrometer\n")

	def RunMieCode(self):
		print(f"Running aerosol optical properties: name {self.name}. May take 10sec...")
		subp.run(f"{self.path}pmie < {self.path}{self.name} > {self.path}res_{self.name}", shell=True)

	def LoadOpticalProperties(self):
		with open("./mie/res_" + self.name, "r") as f:
			for il, line in enumerate(f, 1):
				if il == 16:
					# print(line.strip().split())
					self.c_ext = float(line.strip().split()[2]) * 1e-18 #microm**2 to km**2
					self.c_sca = float(line.strip().split()[6]) * 1e-18 #microm**2 to km**2
				elif il == 18:
					# print(line.strip().split())
					self.ssa = float(line.strip().split()[2])
					self.g = float(line.strip().split()[-1].split(":")[-1])
					break

		self.sca_data = pd.read_fwf(f"{self.path}res_{self.name}", skiprows = 22, names = ["j", "theta", "mu", "i", "q", "u", "d"], usecols=["theta", "mu", "i", "q", "u", "d"], infer_nrows=5000)
		self.sca_data = self.sca_data.dropna()

		# self.sca_data["DoLP"] = np.sqrt(self.sca_data["q"]**2 + self.sca_data["u"]**2) / self.sca_data["i"]
		# self.sca_data["AoLP"] = 0.5 * np.arctan2(self.sca_data["u"], self.sca_data["q"])
		# self.sca_data["DoLP"] = np.sqrt(self.sca_data["q"]**2 + self.sca_data["u"]**2) / self.sca_data["i"]
		# self.sca_data["AoLP"] = 0.5 * np.arctan2(self.sca_data["u"], self.sca_data["q"])
		self.sca_data["AoLP"] = np.where(np.sign(self.sca_data["d"]) >= 0, 0, np.pi/2) # Sign of DoLP indicates the AoLP: >0 AoLP perpendicular to duffusion plane. <0: AoLP parallel to diffusion plane.
		self.sca_data["d"] /= 100. # Normalize DoLP to 0-1 and get rid of the sign.
		# self.sca_data["d"] *= np.sign(self.sca_data["d"]) / 100. # Normalize DoLP to 0-1 and get rid of the sign.

		# self.sca_data["q"] = np.array([min(q, 0) for q in self.sca_data["q"]])
		# self.sca_data["d"] = - self.sca_data["q"] / self.sca_data["i"]

		self.sca_data["theta"] *= DtoR
		self.sca_data["delta_theta"] = np.gradient(self.sca_data["theta"])
		self.sca_data["delta_mu"] = np.gradient(self.sca_data["mu"])


		# if mpi_rank == 0:
		# 	print("AEROSOL phase function")
		#
		# 	print(self.sca_data["theta"])
		# 	print(self.sca_data["i"])
		# 	print(self.sca_data["d"])


		# f, axs = plt.subplots(3) #, sharex=True)
		#
		# axs[0].plot(self.sca_data["theta"], self.sca_data["i"], "k")
		# axs[0].plot(self.sca_data["theta"], self.sca_data["q"], "r")
		# axs[0].plot(self.sca_data["theta"], self.sca_data["u"], "b")
		# axs[0].plot(self.sca_data["theta"], np.sqrt(self.sca_data["q"]**2 + self.sca_data["u"]**2), "g")
		# axs[1].plot(self.sca_data["theta"], self.sca_data["d"]*100, "b")
		# rayleigh_dolp = lambda a: np.sin(a)**2 / (1 + np.cos(a)**2)
		# axs[1].plot(self.sca_data["theta"], rayleigh_dolp(self.sca_data["theta"])*100, "r--")
		# # axs[2].plot(self.sca_data["theta"], np.sqrt(self.sca_data["q"]**2 + self.sca_data["u"]**2) / self.sca_data["i"], "b")
		# # axs[2].plot(self.sca_data["theta"], self.sca_data["AoLP"], "b")
		# r_aer = np.logspace(np.log(self.r_min), np.log(self.r_max))
		# rn = self.rn0 * np.exp(-(np.log(r_aer/self.rmg))**2 / (2 * self.ln_sigma**2)) / r_aer / self.ln_sigma / (2 * np.pi)**0.5
		# axs[2].loglog(r_aer, rn, "b")

		# print(self.c_ext, self.c_sca, self.ssa, self.pi_0)
		# print(sum(self.sca_data["i"] * self.sca_data["delta_theta"]))
		# print(sum(-self.sca_data["i"] * self.sca_data["delta_mu"] / 2))
		# plt.show()





	def PlotPhaseFunction(self, aero_type_list, wavelengths):

		fig, axs = plt.subplots(nrows=len(aero_type_list), ncols=len(aero_type_list), sharex = True)

		rayleigh = lambda t: 3/4 * (1 + np.cos(t)**2)# / 4 / np.pi
		rayleigh_DoLP = lambda t: np.sin(t)**2 / (1 + np.cos(t)**2)

		for iax, aero_type in enumerate(aero_type_list):

			phase_fct = []
			for wl in wavelengths:
				aero = Aerosol(self.in_dict)
				aero.model = aero_type
				aero.wavelength = wl
				aero.name = aero.GetFileName()

				aero.SetModel()

				aero.Compute()

				phase_fct.append([aero.model, aero.wavelength, aero.sca_data["theta"], aero.sca_data["i"], aero.sca_data["d"]])

			if mpi_rank==0:

				for mod, wl, t, i, d in phase_fct:
					if wl == 391.4:
						c = "xkcd:purple"
					elif wl == 557.7:
						c = "xkcd:green"
					elif wl == 620:
						c = "xkcd:orange"
					elif wl == 630:
						c = "xkcd:red"

					axs[iax, 0].plot(np.array(t)*RtoD, i, "-", color=c)
					axs[iax, 0].plot(np.array(t)*RtoD, [rayleigh(theta) for theta in np.array(t)], "--k")

					axs[iax, 1].plot(np.array(t)*RtoD, d, "-", color=c)
					axs[iax, 1].plot(np.array(t)*RtoD, [rayleigh_DoLP(theta) for theta in np.array(t)], "--k")

					# print(t, i)
					print("INTEGRAL OF Phi(theta):", integrate.quad(lambda x: np.sin(x) * np.interp(x, t, i), 0, np.pi)[0] * 2 * np.pi, 4*np.pi)

					# ax.plot(np.array(t)*RtoD, i, color=c, linewidth=2, marker=m, markevery=10, markersize=7)
					# ax.plot(np.array(t)*RtoD, d, color=c, linewidth=2, marker=m, linestyle="dashed", markevery=10, markersize=7)


			axs[iax, 0].set_ylabel(aero_type)
		# axs[0, 0].set_title("Phase function")
		# axs[0, 1].set_title("DoLP")
		axs[0, 0].set_title("Fonction de phase")
		axs[0, 1].set_title("DoLP")

		axs[-1, 1].set_xlabel("Angle de diffusion (°)")
		axs[-1, 0].set_xlabel("Angle de diffusion (°)")
		# axs[-1, 0].set_xlabel("Scattering angle (°)")
		# axs[-1, 1].set_xlabel("Scattering angle (°)")

		axs[0, 0].semilogy()
		axs[1, 0].semilogy()

		# plt.show()





















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
