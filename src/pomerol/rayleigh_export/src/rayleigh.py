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
from sky_map import *
from ground_map import *
from atmosphere import *

class Rayleigh:
	def __init__(self, in_dict):
		"""Initialize the rayleigh object. It will be able to calculate a whole bunch of things!"""

		self.path = in_dict["path"]

		### Get all the azimuts and elevations of the instruments.
		self.a_pc_list		= np.array(in_dict["azimuts"].split(";"))  # List of azimuts for the observations
		if self.a_pc_list[0] == "all": 	self.a_pc_list = np.arange(0, 360, 30) * DtoR
		else: 						self.a_pc_list = np.array([float(a) for a in self.a_pc_list]) * DtoR
		self.e_pc_list		= np.array(in_dict["elevations"].split(";"))  # List of elevation for the observations
		if self.e_pc_list[0] == "all": 	self.e_pc_list = np.arange(10, 90, 30) * DtoR
		else: 						self.e_pc_list = np.array([float(e) for e in self.e_pc_list]) * DtoR
		self.a_pc, self.e_pc = self.a_pc_list[0], self.e_pc_list[0]
		self.Nb_a_pc, self.Nb_e_pc = len(self.a_pc_list), len(self.e_pc_list)

		self.is_single_observation = False
		if self.Nb_a_pc == self.Nb_e_pc == 1:
			self.is_single_observation = True

		### Init the atmosphere object
		self.atmosphere = Atmosphere(in_dict)

 		# Angle de demi-ouverture de l'instrument
		self.ouv_pc			= float(in_dict["instrument_opening_angle"]) * DtoR

		### Initialize the skymap. Even if we don't use it, some parameters must be initialized
		self.sky_map = SkyMap(in_dict)
		if self.sky_map.exist:
			self.sky_map.LoadSkyEmmisionsCube(self.Nb_a_pc, self.Nb_e_pc)
		self.has_sky_emission = self.sky_map.exist
		self.is_time_dependant = self.sky_map.is_time_dependant

		###Show the skymaps if they change in time or we move the instrument
		if self.has_sky_emission and (self.is_time_dependant or not self.is_single_observation):
			self.sky_map.MakeSkyCubePlot(self.a_pc_list, self.e_pc_list, self.ouv_pc)
			plt.show()
			plt.close("all")
		print("Has_sky_emission:", self.has_sky_emission)

		### Initialize the groundmap. Even if we don't use it, some parameters must be initialized
		self.ground_map = GroundMap(in_dict)
		if self.ground_map.exist:
			self.ground_map.LoadGroundEmmisionsMap(self.Nb_a_pc, self.Nb_e_pc)
		self.has_ground_emission = self.ground_map.exist
		print("Has_ground_emission:", self.has_ground_emission)


		if not self.is_single_observation and self.is_time_dependant:
			print("Cannot compute several observation directions with a time changing sky. Please, my code is already messy enough as it is...")
			raise SystemExit

		### Initialize lists to save the final observable of each observation
		self.I_direct_list 	= np.zeros((self.sky_map.Nt, len(self.e_pc_list), len(self.a_pc_list)))
		self.I_list 	= np.zeros((self.sky_map.Nt, len(self.e_pc_list), len(self.a_pc_list)))
		self.DoLP_list 	= np.zeros((self.sky_map.Nt, len(self.e_pc_list), len(self.a_pc_list)))
		self.AoRD_list 	= np.zeros((self.sky_map.Nt, len(self.e_pc_list), len(self.a_pc_list)))

		print("Initialization DONE")


	def GetGeometryFromAzEl(self, a_rd, e_rd, alt):
		"""From known geometry parameters : a, e of emission and alt of scatering , return missing parameters: Distance between emission and scattering and angle of scattering.
		A instrument positions pointing to (a_pc, e_pc)
		E emission points (a_rd, e_rd, h)
		R rayleigh diffusion point (a_pc, e_pc, atmosphere.h_r)"""

		e_rd = max(e_rd, 10**(-30)) # Avoid weird behaviour and division by 0

		v_rd_u = Getuen(a_rd, e_rd) # Vector from observer to emission point
		_, AE = self.obs.GetAH(elevation = e_rd, azimut = a_rd, altitude = self.sky_map.h)  # Distance between observer and emisison point

		RAE = GetAngle(self.v_pc_u, v_rd_u) # Angle between line of sight and emission

		_, AR = self.obs.GetAH(altitude = alt) # Distance between observer and scattering point

		RE = GenPythagore(AR, AE, RAE) # Distance between emission and scaterring point
		ARE = np.arcsin(AE * np.sin(RAE) / RE) # Scattering angle

		return AR, RE, ARE

	def GetGeometryFromAzDist(self, a_rd, d, h_r):
		"""From known geometry parameters : a, d, alt of emission and alt of scatering , return missing parameters: Distance between emission and scattering and angle of scattering.
		A instrument positions pointing to (a_pc, e_pc)
		E emission points (a_rd, e_rd, h)
		R rayleigh diffusion point (a_pc, e_pc, h_r)"""

		d = max(d, 10**(-30)) # Avoid weird behaviour and division by 0

		v_rd_u = Getuen(a_rd, 0.) # Vector from observer to emission point (flat ground, emission at elevation 0)
		AE = d * RT # Distance between observer and emisison point. (exact, following the curvature (~flat earth))

		RAE = GetAngle(self.v_pc_u, v_rd_u) # Angle between line of sight and emission

		_, AR = self.obs.GetAH(altitude = h_r) # Distance between observer and scattering point

		RE = GenPythagore(AR, AE, RAE) # Distance between emission and scaterring point
		ARE = np.arcsin(AE * np.sin(RAE) / RE) # Scattering angle

		return AR, RE, ARE

	def GetScattered(self, I0, AR, ER, RD_angle, alt, elevation = 0):
		"""Given an initial intensity of a source and some geometrical parameter, returns the intensity mesured at the instrument and its DoLP.
		Input parameters: elevation, altitude of source, scattering angle, distance between emission and scattering."""

		w_I = I0 * self.atmosphere.GetVolume(AR, self.ouv_pc) * self.atmosphere.GetParticuleDensity(alt) * self.atmosphere.GetRSCrossSectionParticle(RD_angle)

		if (ER and AR) != 0:
			w_I /= ER ** 2 	# Inverse square law
			w_I /= AR ** 2 	# Inverse square law
		else:
			w_I = 0

		w_DoLP = (1 - np.cos(RD_angle)**2) / (1 + np.cos(RD_angle)**2) # DoLP dependance on scattering angle

		return w_I, w_DoLP

	def GetDirectIntensity(self):
		"""BUGGED !!! Get the direct intensity going in the instrument whithout scattering. BUGGED"""
		direct_intensity = 0

		for ie, e in enumerate(self.sky_map.cube[self.time]):
			# print("DIRECT LIGHT EL", self.sky_map.elevations[ie], self.ouv_pc,  self.e_pc)
			if self.sky_map.elevations[ie] - self.ouv_pc < self.e_pc < self.sky_map.elevations[ie] + self.sky_map.de + self.ouv_pc:
				for ia, I in enumerate(e):
					# print("DIRECT LIGHT AZ", self.sky_map.azimuts[ia], self.a_pc)
					if self.sky_map.azimuts[ia] - self.ouv_pc < self.a_pc < self.sky_map.azimuts[ia] + self.sky_map.da + self.ouv_pc:
						# print("DIRECT LIGHT")
						direct_intensity += I / self.obs.GetAH(azimut=self.a_pc, elevation=self.e_pc, altitude=self.sky_map.h)[1] ** 2
		return direct_intensity

	def ComputeAllMaps(self):
		"""Will compute what the instrument sees for all the maps and all the observations directions"""

		if self.has_ground_emission: ### Compute how the ground map looks through the instrument for every observations directions
			self.is_ground_emission = True
			for ia_pc in range(self.Nb_a_pc):
				for ie_pc in range(self.Nb_e_pc):
					self.ia_pc, self.ie_pc = ia_pc, ie_pc
					self.a_pc, self.e_pc = self.a_pc_list[ia_pc], self.e_pc_list[ie_pc]
					self.time = 0
					self.SingleComputation()

		if self.has_sky_emission: ### Compute how the sky maps looks through the instrument for every time and observations directions
			self.is_ground_emission = False
			for t in range(self.sky_map.Nt):
				# self.SetIMap(time = t)
				for ia_pc in range(self.Nb_a_pc):
					for ie_pc in range(self.Nb_e_pc):
						self.ia_pc, self.ie_pc = ia_pc, ie_pc
						self.a_pc, self.e_pc = self.a_pc_list[ia_pc], self.e_pc_list[ie_pc]
						self.time = t
						self.SingleComputation()

		self.GetLightParametersList()

	def SingleComputation(self):
		"""Will compute and show what the instrument sees for one given emission map (I_map) and one observation"""
		self.ComputeMaps()
		self.MakePlots()

	def GetLightParametersList(self):
		"""When the contributions of all maps for all times and all observations are done, compute the intensity, DoLP, AoLP lists.
		the lists have the shape: (time, e_pc, a_pc)"""
		for t in range(self.sky_map.Nt):
			for ie_pc in range(self.Nb_e_pc):
				for ia_pc in range(self.Nb_a_pc):

					if self.has_sky_emission and self.has_ground_emission: #If sky and ground exist
						self.I_list[t, ie_pc, ia_pc] = np.sum(self.sky_map.scattering_map[t, ie_pc, ia_pc].flatten()) + np.sum(self.ground_map.scattering_map[ie_pc, ia_pc].flatten())
						self.I_direct_list[t, ie_pc, ia_pc] = self.GetDirectIntensity()
					elif self.has_ground_emission: #If sky doesn't exist
						self.I_list[t, ie_pc, ia_pc] = np.sum(self.ground_map.scattering_map[ie_pc, ia_pc].flatten())
					elif self.has_sky_emission:  #If ground doesn't exist
						self.I_list[t, ie_pc, ia_pc] = np.sum(self.sky_map.scattering_map[t, ie_pc, ia_pc].flatten())
						self.I_direct_list[t, ie_pc, ia_pc] = self.GetDirectIntensity()

					self.GetLightParameters(ground = self.has_ground_emission, sky = self.has_sky_emission, time = t, ie_pc = ie_pc, ia_pc = ia_pc)
					self.DoLP_list[t, ie_pc, ia_pc] = self.DoLP
					self.AoRD_list[t, ie_pc, ia_pc] = self.AoRD

					print("t, a_pc, e_pc, I0, DoLP, AoRD:")
					print(t, self.a_pc_list[ia_pc]*RtoD, self.e_pc_list[ie_pc]*RtoD, self.I_list[t, ie_pc, ia_pc], self.DoLP_list[t, ie_pc, ia_pc], self.AoRD_list[t, ie_pc, ia_pc] * RtoD)

	def ComputeMaps(self):
		"""Will initialize the last parameters that depend on the rest and will compute the contribution of each pixels from the map we set. We set the map by setting self.is_ground_emission to True or False. If False, then compute the sky map at the time set"""
		self.Nalt = int((self.atmosphere.h_r_max - self.atmosphere.h_r_min) / (self.atmosphere.d_los * np.sin(self.e_pc)))		# Number of bins along the line of sight between atmosphere.h_r_min and atmosphere.h_r_max of length atmosphere.d_los
		self.dh = int((self.atmosphere.h_r_max - self.atmosphere.h_r_min) / self.Nalt) #Height of each bin
		self.altitudes = np.linspace(self.atmosphere.h_r_min, self.atmosphere.h_r_max, self.Nalt) #List of all altitudes

		### Get the base name for the saving file -> unique for each parameter set
		self.save_name = "results/" + self.ground_map.location + "_" + self.sky_map.mode + "_a" + str(np.round(self.a_pc*RtoD, 0)) + "_e" + str(np.round(self.e_pc*RtoD, 0)) + "_h" + str(np.round(self.sky_map.h, 0)) + "_" + str(np.round(self.atmosphere.h_r_min, 0)) + "-" + str(np.round(self.atmosphere.h_r_max, 0)) + "_dlos" + str(np.round(self.atmosphere.d_los))

		if self.ground_map.radius > 0:
			self.save_name += "_D" + str(np.round(self.ground_map.radius * RT, 0))

		print("*******************************************************************************************************************************")
		self.already_done = False
		if self.already_done:
			print("Observation already done. Loading old results from file", self.save_name)
		else:
			self.obs = ObservationPoint(self.ground_map.A_lon, self.ground_map.A_lat, self.sky_map.h, self.a_pc, self.e_pc) # Usefull observation object
			print("STARTING calculation for obs:")
			print(self.obs)
			self.v_pc_u = Getuen(self.a_pc, self.e_pc) # Vector of observation (line of sight) in UEN system


			if self.is_ground_emission == True:
				self.ComputeGroundMaps()
			else:
				self.ComputeSkyMaps()


		print("Computing DONE")


	def ComputeSkyMaps(self):
		"""Compute the contribution of the sky map at the time set"""
		# try:
		# 	self.sky_scattering_map = np.loadtxt(self.path + self.save_name + "_scattering_map.cvs")
		# 	self.sky_DoLP_map = np.loadtxt(self.path + self.save_name + "_DoLP_map.cvs")
		# 	self.sky_total_scattering_map = np.loadtxt(self.path + self.save_name + "_total_scattering_map.cvs")
		# 	self.sky_AoRD_map = np.loadtxt(self.path + self.save_name + "_AoRD_map.cvs")
		# 	self.sky_already_done = True
		# except:
		# 	# Initialisation of all the maps (np.ndarray)
		# 	self.sky_scattering_map 		= np.zeros(self.sky_maps_shape) # Intensity from (e, a) reaching us
		# 	self.sky_DoLP_map				= np.zeros(self.sky_maps_shape) # Polarized intensity from (e, a) reaching us
		# 	self.sky_total_scattering_map 	= np.zeros(self.sky_maps_shape) # DoLP of scattered light from (e,a)
		# 	self.sky_AoRD_map 				= np.zeros(self.sky_maps_shape) # Angle of polaisation of light from (e,a)
		# 	self.sky_already_done = False


		count = 0
		start_time = tm.time()
		N = self.sky_map.N * self.Nalt
		print(N, "bins to compute")
		for ie, e in enumerate(self.sky_map.elevations):
			cos_e = np.cos(e)
			for ia, a in enumerate(self.sky_map.azimuts):
				if self.sky_map.cube[self.time, ie, ia] > 0:
					for ialt, alt in enumerate(self.altitudes): #for every altitude between minimum and maximum scattering altitude

						AR, RE, RD_angle = self.GetGeometryFromAzEl(a, e, alt)
						w_I, w_DoLP = self.GetScattered(self.sky_map.cube[self.time, ie, ia], AR, RE, RD_angle, alt, elevation = e)

						w_I *= cos_e # Bining effect, bin around e=90 are "smaller" if sky emission

						self.sky_map.scattering_map[self.time, self.ie_pc, self.ia_pc, ie, ia] += w_I # Intensity of light scattered from a given (e, a)
						self.sky_map.DoLP_map[self.time, self.ie_pc, self.ia_pc, ie, ia] += w_DoLP * w_I # Intensity of polarized light scattered from a given (e, a). Given a point source, the AoRD is the same for every altitude -> the addition makes sens

						count += 1
						if count == N // 10:
							print("\t", 9 * (tm.time() - start_time), "seconds left...")
						elif count == N // 2:
							print("\t", (tm.time() - start_time), "seconds left...")

					self.sky_map.total_scattering_map = self.sky_map.DoLP_map / self.sky_map.scattering_map # DoLP of a given (e, a)

				self.sky_map.AoRD_map[self.time, self.ie_pc, self.ia_pc, ie, ia] = self.obs.GetRayleighAngle(a, e)

	def ComputeGroundMaps(self):
		"""Compute the contribution of the ground map at the time set"""
		# try:
		# 	self.ground_map.total_scattering_map = np.loadtxt(self.path + self.save_name + "_scattering_map.cvs")
		# 	self.ground_map.DoLP_map = np.loadtxt(self.path + self.save_name + "_DoLP_map.cvs")
		# 	self.ground_map.total_scattering_map = np.loadtxt(self.path + self.save_name + "_total_scattering_map.cvs")
		# 	self.ground_map.AoRD_map = np.loadtxt(self.path + self.save_name + "_AoRD_map.cvs")
		# 	self.ground_already_done = True
		# except:
		# 	# Initialisation of all the maps (np.ndarray)
		# 	self.ground_map.total_scattering_map 		= np.zeros(self.ground_maps_shape) # Intensity from (e, a) reaching us
		# 	self.ground_map.DoLP_map				= np.zeros(self.ground_maps_shape) # Polarized intensity from (e, a) reaching us
		# 	self.ground_map.total_scattering_map 	= np.zeros(self.ground_maps_shape) # DoLP of scattered light from (e,a)
		# 	self.ground_map.AoRD_map 				= np.zeros(self.ground_maps_shape) # Angle of polaisation of light from (e,a)
		# 	self.ground_already_done = False

		count = 0
		start_time = tm.time()

		self.N = len(self.ground_map.longitudes) * len(self.ground_map.latitudes) * self.Nalt
		print(self.N, "bins to compute")
		for ilat, lat in enumerate(self.ground_map.latitudes):
			for ilon, lon in enumerate(self.ground_map.longitudes):
				if self.ground_map.I_map[ilat, ilon] > 0:
					a_rd, e_rd = LonLatToAzDist(lon, lat, self.ground_map.A_lon, self.ground_map.A_lat)
					for ialt, alt in enumerate(self.altitudes): #for every altitude between minimum and maximum scattering altitude
						AR, RE, RD_angle = self.GetGeometryFromAzDist(a_rd, e_rd, alt)
						w_I, w_DoLP = self.GetScattered(self.ground_map.I_map[ilat, ilon], AR, RE, RD_angle, alt)

						self.ground_map.scattering_map[self.ie_pc, self.ia_pc, ilat, ilon] += w_I # Intensity of light scattered from a given (e, a)
						self.ground_map.DoLP_map[self.ie_pc, self.ia_pc, ilat, ilon] += (w_DoLP * w_I) # Intensity of polarized light scattered from a given (e, a). Given a point source, the AoRD is the same for every altitude -> the addition makes sens

						count += 1
						if count == self.N // 10:
							print("\t", 9 * (tm.time() - start_time), "seconds left...")
						elif count == self.N // 2:
							print("\t", (tm.time() - start_time), "seconds left...")


					self.ground_map.total_scattering_map = self.ground_map.DoLP_map / self.ground_map.scattering_map # DoLP of a given (e, a)

				self.ground_map.AoRD_map[self.ie_pc, self.ia_pc, ilat, ilon] = self.obs.GetRayleighAngle(a_rd, 0) # e_rd is NOT the elevation if ground emission


	def GetLightParameters(self, ground = False, sky = False, time = None, ie_pc = None, ia_pc = None):
		"""Once the contributions of emmision maps are computed, compute the I, DOLP, AoLP and the AoLP histogram."""

		if time == None: time = self.time
		if ie_pc == None: ie_pc = self.ie_pc
		if ia_pc == None: ia_pc = self.ia_pc


		###Compute the AoLP contribution histogram.
		self.N_bins = 180
		self.bins = np.linspace(-np.pi/2, np.pi/2, self.N_bins + 1, endpoint=True)
		self.width = np.pi / self.N_bins

		self.hst = np.zeros(self.N_bins)

		if sky:
			sky_hst, b = np.histogram(self.sky_map.AoRD_map[time, ie_pc, ia_pc, :, :], bins=self.bins, weights=self.sky_map.DoLP_map[time, ie_pc, ia_pc, :, :], density = True)
			self.hst += sky_hst
		if ground:
			ground_hst, b = np.histogram(self.ground_map.AoRD_map[ie_pc, ia_pc, :, :], bins=self.bins, weights=self.ground_map.DoLP_map[ie_pc, ia_pc, :, :], density = True)
			self.hst += ground_hst

		self.hst /= sum(self.hst)

		###Simulate the instrument with a rotating polarizing filter to get V, Vcos, Vsin and then I, DoLP, AoLP
		Ns = 2 * len(self.bins)
		rs_signal = np.zeros(Ns)
		filter_orientation = np.linspace(0, 2 * np.pi, Ns)
		for i_f, f in enumerate(filter_orientation):
			for ihist, hist in enumerate(self.hst):
				rs_signal[i_f] += hist * np.cos(b[ihist] - f) ** 2

		self.V = np.average(rs_signal)
		self.Vcos = np.average(rs_signal * np.cos(2 * filter_orientation))
		self.Vsin = np.average(rs_signal * np.sin(2 * filter_orientation))

		# print("V, Vcos, Vsin", self.V, self.Vcos, self.Vsin)

		self.I0 = 2 * self.V
		self.DoLP = 100 * 2 * np.sqrt(self.Vcos ** 2 + self.Vsin ** 2) / self.V
		self.AoRD = np.arctan2(self.Vsin, self.Vcos) / 2

		# print(self.I0, self.DoLP, self.AoRD)

		return self.I0, self.DoLP, self.AoRD


	def MakePlots(self):
		################################################################################
		###	PLOTTING
		################################################################################

		################################################################################
		###	Polar plots of intensity

		if self.is_ground_emission:
			self.MakeGroundMapPlots()
		else:
			self.MakeSkyMapPlots()

		self.MakeAoLPHist(ground = self.is_ground_emission, sky = not self.is_ground_emission)

		self.MakeAoLPMap()

		if self.is_single_observation and not self.is_time_dependant:
			plt.show()

		plt.close("all")


	def MakeSkyMapPlots(self):
		"""Make plots of the sky emission map contributions. Plot the initial intensity, the scattered intensity, the DOLP..."""
		f1, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex = True, sharey = True, figsize=(16, 8))
		ax1 = plt.subplot(221, projection='polar')
		ax2 = plt.subplot(222, projection='polar')
		ax3 = plt.subplot(223, projection='polar')
		ax4 = plt.subplot(224, projection='polar')

		# So that e = 90 is in the center of the plot if the emission layer is NOT the ground
		el = 90 - (self.sky_map.elevations) * RtoD
		# i1 = ax1.pcolormesh(azimuts, el, scattering_map)
		i1 = ax1.pcolormesh(self.sky_map.azimuts, el, self.sky_map.cube[self.time, :, :])
		i2 = ax2.pcolormesh(self.sky_map.azimuts, el, self.sky_map.scattering_map[self.time, self.ie_pc, self.ia_pc, :, :])
		i3 = ax3.pcolormesh(self.sky_map.azimuts, el, self.sky_map.DoLP_map[self.time, self.ie_pc, self.ia_pc, :, :])
		# i4 = ax4.pcolormesh(azimuts, el, total_scattering_map)
		i4 = ax4.pcolormesh(self.sky_map.azimuts, el, self.sky_map.AoRD_map[self.time, self.ie_pc, self.ia_pc, :, :]*RtoD, cmap=cm.twilight)

		cbar1 = f1.colorbar(i1, extend='both', spacing='proportional', shrink=0.9, ax=ax1)
		cbar1.set_label('Initial emission map')
		cbar2 = f1.colorbar(i2, extend='both', spacing='proportional', shrink=0.9, ax=ax2)
		cbar2.set_label('Total intensity map')
		cbar3 = f1.colorbar(i3, extend='both', spacing='proportional', shrink=0.9, ax=ax3)
		cbar3.set_label('Polarised intensity map')
		cbar4 = f1.colorbar(i4, extend='both', spacing='proportional', shrink=0.9, ax=ax4)
		cbar4.set_label('DoLP')

		f1.suptitle("Relevant angles map at skibotn, with light pollution source at -45deg in azimut")

		for a in [ax1, ax2, ax3, ax4]:
			a.set_theta_zero_location("N")
			a.set_theta_direction(-1)
			a.set_thetamin(self.sky_map.I_zone_a_min * RtoD)
			a.set_thetamax(self.sky_map.I_zone_a_max * RtoD)

			a.set_rlim(0,90, 1)
			a.set_yticks(np.arange(0, 90, 20))
			a.set_yticklabels(a.get_yticks()[::-1])   # Change the labels

			a.add_artist(Ellipse((self.a_pc, 90 - self.e_pc * RtoD), width=self.ouv_pc, height=self.ouv_pc*RtoD, color="red"))

		plt.savefig(self.path + self.save_name + '_skymaps.png', bbox_inches='tight')

	def MakeGroundMapPlots(self):
		"""Make plots of the ground emission map contributions. Plot the initial intensity, the scattered intensity, the DOLP..."""
		f2, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex = True, sharey = True, figsize=(16, 8))
		ax1 = plt.subplot(221)
		ax2 = plt.subplot(222)
		ax3 = plt.subplot(223)
		ax4 = plt.subplot(224)

		i1 = ax1.pcolormesh((self.ground_map.longitudes-self.ground_map.A_lon) * RT, (self.ground_map.latitudes-self.ground_map.A_lat) * RT, self.ground_map.I_map)
		i2 = ax2.pcolormesh((self.ground_map.longitudes-self.ground_map.A_lon) * RT, (self.ground_map.latitudes-self.ground_map.A_lat) * RT, self.ground_map.total_scattering_map[self.ie_pc, self.ia_pc, :, :])		# Intensity from (e, a) reaching us
		i3 = ax3.pcolormesh((self.ground_map.longitudes-self.ground_map.A_lon) * RT, (self.ground_map.latitudes-self.ground_map.A_lat) * RT, self.ground_map.DoLP_map[self.ie_pc, self.ia_pc, :, :])				# Polarized intensity from (e, a) reaching us
		i4 = ax4.pcolormesh((self.ground_map.longitudes-self.ground_map.A_lon) * RT, (self.ground_map.latitudes-self.ground_map.A_lat) * RT, self.ground_map.total_scattering_map[self.ie_pc, self.ia_pc, :, :])	# DoLP of scattered light from (e,a)
		# i4 = ax4.pcolormesh((longitudes-A_lon) * RT, (latitudes-A_lat) * RT, AoRD_map * RtoD, cmap=cm.twilight)			# AoLP of scattered light from (e,a)


	# Angle of polaisation of light from (e,a)
		cbar1 = f2.colorbar(i1, extend='both', spacing='proportional', shrink=0.9, ax=ax1)
		cbar1.set_label('Initial emission map')
		cbar2 = f2.colorbar(i2, extend='both', spacing='proportional', shrink=0.9, ax=ax2)
		cbar2.set_label('Total intensity map')
		cbar3 = f2.colorbar(i3, extend='both', spacing='proportional', shrink=0.9, ax=ax3)
		cbar3.set_label('Polarised intensity map')
		cbar4 = f2.colorbar(i4, extend='both', spacing='proportional', shrink=0.9, ax=ax4)
		cbar4.set_label('DoLP')

		f2.suptitle("Relevant angles map at skibotn, with light pollution source at -45deg in azimut")

		dy = np.cos(self.a_pc) * self.atmosphere.h_r_max / np.tan(self.e_pc)
		dx = np.sin(self.a_pc) * self.atmosphere.h_r_max / np.tan(self.e_pc)
		for a in [ax1, ax2, ax3, ax4]:
			a.add_artist(Arrow(0, 0, dx, dy, color="red"))

		plt.savefig(self.path + self.save_name + '_groundmaps.png', bbox_inches='tight')


	def MakeAoLPHist(self, ground = False, sky = False):
		"""Make an pyplot histogram of all AoLP contributionsThe histogram is calculated in GetLightParameters()
		Need a call to plt.show() after calling this function."""
		f3, ax = plt.subplots(1, figsize=(16, 8))

		ax = plt.subplot(111, projection='polar')
		ax.set_theta_zero_location("N")
		ax.set_theta_direction(-1)
		ax.set_thetamin(-90)
		ax.set_thetamax(90)

		I0, DoLP, AoRD = self.GetLightParameters(ground=ground, sky=sky)

		bars = ax.bar(self.bins[:self.N_bins], self.hst, width=self.width)

		ax.set_title("Weighted AoRD: I0 = " + str(np.round(I0, 1)) + " DoLP = " + str(np.round(DoLP, 1)) + " AoRD = " + str(np.round(AoRD*RtoD, 1)))

		ax.plot([AoRD, AoRD], [0, max(self.hst)], "r")

		plt.savefig(self.path + self.save_name + '_AoRD_hist.png', bbox_inches='tight')


	def MakeAoLPMap(self):
		"""Make a pyplot quiver of the AoLP contribution of each emission point. Vertical means AoLP=0, horyzontal AoLP = 90.
		Need a call to plt.show() after calling this function."""
		f4, (ax1) = plt.subplots(1, sharex=True, figsize=(16, 8))

		ax1.set_xlabel("Azimuts")
		ax1.set_ylabel("Elevations")

		if not self.is_ground_emission:
			X, Y = self.sky_map.azimuts * RtoD, self.sky_map.elevations * RtoD
			U, V = np.sin(self.sky_map.AoRD_map[self.time, self.ie_pc, self.ia_pc, :, :]) * self.sky_map.DoLP_map[self.time, self.ie_pc, self.ia_pc, :, :], np.cos(self.sky_map.AoRD_map[self.time, self.ie_pc, self.ia_pc, :, :]) * self.sky_map.DoLP_map[self.time, self.ie_pc, self.ia_pc, :, :]
			M = self.sky_map.DoLP_map[self.time, self.ie_pc, self.ia_pc, :, :]
		else:
			X, Y = self.ground_map.longitudes * RtoD, self.ground_map.latitudes * RtoD
			U, V = np.sin(self.ground_map.AoRD_map[self.ie_pc, self.ia_pc, :, :]) * self.ground_map.DoLP_map[self.ie_pc, self.ia_pc, :, :], np.cos(self.ground_map.AoRD_map[self.ie_pc, self.ia_pc, :, :]) * self.ground_map.DoLP_map[self.ie_pc, self.ia_pc, :, :]
			M = self.ground_map.DoLP_map[self.ie_pc, self.ia_pc, :, :]


		q = ax1.quiver(X, Y, U, V, M, pivot="middle", headwidth=0, headlength=0, headaxislength=0)
		ax1.quiverkey(q, 0.5, 1.1, 1, "AoRD")
		if not self.is_ground_emission:
			ax1.add_artist(Circle((self.a_pc * RtoD, self.e_pc * RtoD), radius=self.ouv_pc * RtoD, color="red"))

		plt.savefig(self.path + self.save_name + '_AoRD_map.png', bbox_inches='tight')

		################################################################################
		###	AoRD map as stream plot (continuous lines)
	# def MakeAoLPStram(self):
		# f5, (ax1) = plt.subplots(1, sharex=True, figsize=(16, 8))
		#
		# ax1.set_xlabel("Azimuts")
		# ax1.set_ylabel("Elevations")
		#
		# if h != 0: 	X, Y = azimuts * RtoD, elevations * RtoD
		# else:		X, Y = longitudes * RtoD, latitudes * RtoD
		# U, V = np.sin(AoRD_map) * DoLP_map, np.cos(AoRD_map) * DoLP_map
		# M = DoLP_map
		#
		# ax1.streamplot(X, Y, U, V, color=M, density=2, arrowstyle="-")
		# if h != 0: 	ax1.add_artist(Circle((a_pc * RtoD, e_pc * RtoD), radius=ouv_pc * RtoD, color="red"))

		# np.savetxt(self.path + self.save_name + "_scattering_map.cvs", self.scattering_map)
		# np.savetxt(self.path + self.save_name + "_DoLP_map.cvs", self.DoLP_map)
		# np.savetxt(self.path + self.save_name + "_total_scattering_map.cvs", self.total_scattering_map)
		# np.savetxt(self.path + self.save_name + "_AoRD_map.cvs", self.AoRD_map)


	def MakeSummaryPlot(self):

		# print("Azimut, Elevation, I, DoLP, AoLP")
		# for t in range(self.sky_map.Nt):
		# 	for ia, a in enumerate(self.a_pc_list):
		# 		for ie, e in enumerate(self.e_pc_list):
		# 			print(a*RtoD, e*RtoD, self.I_list[t][ie][ia], self.DoLP_list[t][ie][ia], self.AoRD_list[t][ie][ia]*RtoD)

		if not self.is_time_dependant:
			if self.is_single_observation:
				self.MakeAoLPHist(ground = self.has_ground_emission, sky = self.has_sky_emission)

			else:
				f, axs = plt.subplots(3, sharex = True)
				axs[0] = plt.subplot(311)
				axs[1] = plt.subplot(312)
				axs[2] = plt.subplot(313)

				if self.Nb_e_pc == 1:
					# print("PLOT", self.a_pc_list*RtoD, self.I_list[0, 0, :].tolist(), self.DoLP_list[0, 0, :].tolist(), self.AoRD_list[0, 0, :].tolist())
					axs[0].plot(self.a_pc_list*RtoD, self.I_list[0, 0, :].tolist())
					axs[1].plot(self.a_pc_list*RtoD, self.DoLP_list[0, 0, :].tolist())
					axs[2].plot(self.a_pc_list*RtoD, (self.AoRD_list[0, 0, :]*RtoD).tolist())
					self.save_name = 'results/' + self.ground_map.location + "_rot_e" + str(self.e_pc_list[0]*RtoD) + ".png"
				elif self.Nb_a_pc == 1:
					axs[0].plot(self.e_pc_list*RtoD, self.I_list[0, :, 0].tolist())
					axs[1].plot(self.e_pc_list*RtoD, self.DoLP_list[0, :, 0].tolist())
					axs[2].plot(self.e_pc_list*RtoD, (self.AoRD_list[0, :, 0]*RtoD).tolist())
					self.save_name = 'results/' + self.ground_map.location + "_rot_a" + str(self.a_pc_list[0]*RtoD) + ".png"
				else:
					axs[0].imshow(self.I_list[0,:,:])
					axs[1].imshow(self.DoLP_list[0,:,:])
					axs[2].imshow(self.AoRD_list[0,:,:]*RtoD)
					self.save_name = 'results/' + self.ground_map.location + "_skymap" + ".png"

		else:
			f, axs = plt.subplots(3, sharex = True)
			axs[0] = plt.subplot(311)
			axs[1] = plt.subplot(312)
			axs[2] = plt.subplot(313)

			self.sky_map.MakeSkyCubePlot(self.a_pc_list, self.e_pc_list, self.ouv_pc)
			axs[0].set_yscale('log')
			# axs[0].plot(range(self.sky_map.Nt), (self.I_list[:, 0, 0] / (self.I_direct_list[:, 0, 0] + self.I_list[:, 0, 0])).tolist())
			# print(self.I_list[:, 0, 0].tolist())
			# print(self.I_direct_list[:, 0, 0].tolist())
			axs[0].plot(range(self.sky_map.Nt), self.I_list[:, 0, 0].tolist())
			# axs[0].plot(range(self.sky_map.Nt), self.I_direct_list[:, 0, 0].tolist())
			axs[1].plot(range(self.sky_map.Nt), self.DoLP_list[:, 0, 0].tolist())
			axs[2].plot(range(self.sky_map.Nt), (self.AoRD_list[:, 0, 0]*RtoD).tolist())
			self.save_name = 'results/' + self.ground_map.location + "_" + self.sky_map.mode + "_time_series_a" + str(self.a_pc_list[0]*RtoD) + "_e" + str(self.e_pc_list[0]*RtoD) + ".png"

		# plt.savefig(self.path + self.save_name, bbox_inches='tight')


if __name__ == "__main__":

	#Get the input file from the command arguments
	all_time_start = tm.time()
	arguments = sys.argv
	nb_args = len(arguments)

	try:
		in_dict = ReadInputFile("./input_files/" + arguments[1] + ".in")
		print("Correct input file in use:", arguments[1])
	except:
		in_dict = ReadInputFile("./input_files/RS_default.in")
		print("WARNING: Wrong or no input file specified, default in use.")

	#Init the rayleigh object
	rayleigh = Rayleigh(in_dict)
	rayleigh.ComputeAllMaps()
	rayleigh.MakeSummaryPlot()

	print("ALL TIME SINCE START:", tm.time() - all_time_start)

	plt.show()
