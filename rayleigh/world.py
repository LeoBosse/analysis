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

from tqdm import tqdm

from observation import *
from rayleigh_utils import *
from sky_map import *
from ground_map import *
from atmosphere import *


class World:
	def __init__(self, in_dict):
		self.wavelength = float(in_dict["wavelength"])

		### Get all the azimuts and elevations of the instruments.
		self.a_pc_list		= np.array(in_dict["azimuts"].split(";"))  # List of azimuts for the observations
		if self.a_pc_list[0] == "all":
			self.a_pc_list = np.arange(0, 360, 10) * DtoR
		else:
			self.a_pc_list = np.array([float(a) for a in self.a_pc_list]) * DtoR

		self.e_pc_list		= np.array(in_dict["elevations"].split(";"))  # List of elevation for the observations
		if self.e_pc_list[0] == "all":
			self.e_pc_list = np.linspace(10, 90, 9, endpoint=True) * DtoR
		else:
			self.e_pc_list = np.array([float(e) for e in self.e_pc_list]) * DtoR

		self.a_pc, self.e_pc = self.a_pc_list[0], self.e_pc_list[0]

		self.Nb_a_pc, self.Nb_e_pc = len(self.a_pc_list), len(self.e_pc_list)

		### Init the atmosphere object
		self.atmosphere = Atmosphere(in_dict)

		# Angle de demi-ouverture de l'instrument (1deg pour les crus)
		self.ouv_pc		= float(in_dict["instrument_opening_angle"]) * DtoR
		# Aire de l'instrument donnée en cm^2, transformée en m^2
		self.PTCU_area	= float(in_dict["instrument_area"]) * 10 ** (-4)

		### Initialize the skymap. Even if we don't use it, some parameters must be initialized
		self.sky_map = SkyMap(in_dict)
		if self.sky_map.exist:
			self.sky_map.LoadSkyEmmisionsCube(self.Nb_a_pc, self.Nb_e_pc)
		self.has_sky_emission = self.sky_map.exist
		self.is_time_dependant = self.sky_map.is_time_dependant

		self.is_single_observation = False
		if self.Nb_a_pc == self.Nb_e_pc == 1:
			self.is_single_observation = True

		###Show the skymaps if they change in time or we move the instrument
		if self.has_sky_emission and (self.is_time_dependant or not self.is_single_observation):
			self.sky_map.MakeSkyCubePlot(self.a_pc_list, self.e_pc_list, self.ouv_pc)
			plt.show()
			# plt.close("all")
		print("Has_sky_emission:", self.has_sky_emission)

		### Initialize the groundmap. Even if we don't use it, some parameters must be initialized
		self.ground_map = GroundMap(in_dict)
		if self.ground_map.exist:
			if float(in_dict["point_src_I0"]) > 0:
				I0 = float(in_dict["point_src_I0"])
				az = float(in_dict["point_src_az"]) * DtoR
				dist = float(in_dict["point_src_dist"])
				self.ground_map.LoadPointSourceMap(I0, az, dist, self.Nb_a_pc, self.Nb_e_pc)
			else:
				self.ground_map.LoadGroundEmmisionsMap(self.Nb_a_pc, self.Nb_e_pc)
		self.has_ground_emission = self.ground_map.exist
		print("Has_ground_emission:", self.has_ground_emission)


		### Initialize the atlitude maps
		self.alt_map = ElevationMap(in_dict)
		self.use_elevation_map = bool(int(in_dict["use_alt_map"]))
		# self.alt_map.PlotMap()

		self.instrument_altitude = self.alt_map.GetAltitudeFromLonLat(self.ground_map.A_lon, self.ground_map.A_lat) / 1000.#In km



	def SetObservation(self, a_pc, e_pc):

		#Create an observation object
		self.obs = ObservationPoint(self.ground_map.A_lon, self.ground_map.A_lat, self.sky_map.h, a_pc, e_pc, A_alt = self.instrument_altitude) # Usefull observation object

		range_min, range_max = self.atmosphere.h_r_min / np.sin(e_pc) , self.atmosphere.h_r_max / np.sin(e_pc)

		#numbers of scattering point along the line of sight
		los_length = range_max - range_min

		if los_length >= self.atmosphere.d_los:
			#array of range for each scattering points (distance along los from the instrument)
			dlos_list = np.arange(range_min + self.atmosphere.d_los / 2,  range_max - self.atmosphere.d_los / 2, self.atmosphere.d_los)
		else:
			dlos_list = [2.]
			# dlos_list = [(range_min + range_max) / 2.]

		print("DEBUG dlos list:", dlos_list)

		#array of lon, lat, alt (above sea level) for each scattering point
		self.dlos_list = [self.obs.GetPCoordinatesFromRange(d) for d in dlos_list]

		#list of absolute altitudes (above sea level) for each scattering point
		self.altitudes = [d[2] for d in self.dlos_list]
		# self.altitudes = [self.alt_map.GetAltitudeAboveGround(lon, lat, a) for lon, lat, a in dlos_list]
		self.Nalt = len(self.altitudes)


		# self.Nalt = int((self.atmosphere.h_r_max - self.atmosphere.h_r_min) / (self.atmosphere.d_los * np.sin(e_pc)))		# Number of bins along the line of sight between atmosphere.h_r_min and atmosphere.h_r_max of length atmosphere.d_los
		# if self.Nalt == 0:
		# 	sys.exit("ERROR: Nalt == 0")
		#
		# self.dh = int((self.atmosphere.h_r_max - self.atmosphere.h_r_min) / self.Nalt) #Height of each bin
		# self.altitudes = np.linspace(self.atmosphere.h_r_min, self.atmosphere.h_r_max, self.Nalt) #List of all altitudes


		self.v_pc_u = Getuen(a_pc, e_pc) # Vector of observation (line of sight) in UEN system


	# Print iterations progress
	def Progress(self, iteration, total, prefix = '', suffix = '', decimals = 1, length = 50, fill = '█', printEnd = "\r"):
		"""
		Call in a loop to create terminal progress bar
		@params:
		    iteration   - Required  : current iteration (Int)
		    total       - Required  : total iterations (Int)
		    prefix      - Optional  : prefix string (Str)
		    suffix      - Optional  : suffix string (Str)
		    decimals    - Optional  : positive number of decimals in percent complete (Int)
		    length      - Optional  : character length of bar (Int)
		    fill        - Optional  : bar fill character (Str)
		    printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
		"""
		percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
		filledLength = int(length * iteration // total)
		bar = fill * filledLength + '-' * (length - filledLength)

		sys.stdout.write('\r%s |%s| %s%% %s %s' % (prefix, bar, percent, suffix, printEnd))
		sys.stdout.flush()
		# Print New Line on Complete
		if iteration == total:
			print("\n")

	def ComputeSkyMaps(self, time, ia_pc, ie_pc):
		"""Compute the contribution of the sky map at the time set"""
		# try:
		# 	self.sky_scattering_map = np.loadtxt(self.path + self.save_name + "_scattering_map.cvs")
		# 	self.sky_DoLP_map = np.loadtxt(self.path + self.save_name + "_DoLP_map.cvs")
		# 	self.sky_total_scattering_map = np.loadtxt(self.path + self.save_name + "_total_scattering_map.cvs")
		# 	self.sky_AoRD_map = np.loadtxt(self.path + self.save_name + "_AoRD_map.cvs")
		# 	self.sky_already_done = True
		# except:
		# 	# Initialisation of all the maps (np.ndarray)
		# 	self.sky_scattering_map 		= np.zeros(self.world.sky_maps_shape) # Intensity from (e, a) reaching us
		# 	self.sky_DoLP_map				= np.zeros(self.world.sky_maps_shape) # Polarized intensity from (e, a) reaching us
		# 	self.sky_total_scattering_map 	= np.zeros(self.world.sky_maps_shape) # DoLP of scattered light from (e,a)
		# 	self.sky_AoRD_map 				= np.zeros(self.world.sky_maps_shape) # Angle of polaisation of light from (e,a)
		# 	self.sky_already_done = False


		count = 0
		start_time = tm.time()
		N = self.sky_map.N * self.Nalt
		print(N, "bins to compute")
		with tqdm(total=N, file=sys.stdout) as pbar:
			for ie, e in enumerate(self.sky_map.mid_elevations):
				for ia, a in enumerate(self.sky_map.mid_azimuts):
					if self.sky_map.cube[time, ie, ia] >= 0:
						for ialt, alt in enumerate(self.altitudes): #for every altitude between minimum and maximum scattering altitude

							sca_from_E_is_visible = True
							if self.use_elevation_map:
								sca_from_E_is_visible = self.alt_map.IsVisible(sca_lon, sca_lat, alt1=sca_alt, lon2=lon, lat2=lat, dlos = 0.5)

							if sca_from_E_is_visible:
								I0, w_DoLP = self.ComputeSingleRSSkyPointSource(time, ia, a, ie, e, alt)
								self.sky_map.total_scattering_map[time, ie_pc, ia_pc, ie, ia] += I0 # Intensity of light scattered from a given (e, a)
								self.sky_map.scattering_map[time, ie_pc, ia_pc, ie, ia] += w_DoLP * I0 # Intensity of polarized light scattered from a given (e, a). Given a point source, the AoRD is the same for every altitude -> the addition makes sens

							# count += 1
							# self.Progress(count, N, suffix="of sky point sources done")
							pbar.update(1)
							# if count == N // 10:
							# 	print("\t", 9 * (tm.time() - start_time), "seconds left...")
							# elif count == N // 2:
							# 	print("\t", (tm.time() - start_time), "seconds left...")

						self.sky_map.DoLP_map = self.sky_map.scattering_map / self.sky_map.total_scattering_map # DoLP of a given (e, a)

					self.sky_map.AoRD_map[time, ie_pc, ia_pc, ie, ia] = self.obs.GetRayleighAngle(a, e)

	def ComputeGroundMaps(self, time, ia_pc, ie_pc):
		"""Compute the contribution of the ground map at the time set"""
		# try:
		# 	self.world.ground_map.total_scattering_map = np.loadtxt(self.path + self.save_name + "_scattering_map.cvs")
		# 	self.world.ground_map.DoLP_map = np.loadtxt(self.path + self.save_name + "_DoLP_map.cvs")
		# 	self.world.ground_map.total_scattering_map = np.loadtxt(self.path + self.save_name + "_total_scattering_map.cvs")
		# 	self.world.ground_map.AoRD_map = np.loadtxt(self.path + self.save_name + "_AoRD_map.cvs")
		# 	self.ground_already_done = True
		# except:
		# 	# Initialisation of all the maps (np.ndarray)
		# 	self.world.ground_map.total_scattering_map 		= np.zeros(self.world.ground_maps_shape) # Intensity from (e, a) reaching us
		# 	self.world.ground_map.DoLP_map				= np.zeros(self.world.ground_maps_shape) # Polarized intensity from (e, a) reaching us
		# 	self.world.ground_map.total_scattering_map 	= np.zeros(self.world.ground_maps_shape) # DoLP of scattered light from (e,a)
		# 	self.world.ground_map.AoRD_map 				= np.zeros(self.world.ground_maps_shape) # Angle of polaisation of light from (e,a)
		# 	self.ground_already_done = False

		count = 0
		start_time = tm.time()

		N = len(self.ground_map.longitudes) * len(self.ground_map.latitudes) * self.Nalt
		print(N, "bins to compute")
		with tqdm(total=N, file=sys.stdout) as pbar:
			for ilat, lat in enumerate(self.ground_map.latitudes):
				for ilon, lon in enumerate(self.ground_map.longitudes):
					if self.ground_map.I_map[ilat, ilon] > 0:
						a_rd, e_rd = LonLatToAzDist(lon, lat, self.ground_map.A_lon, self.ground_map.A_lat)
						for sca_lon, sca_lat, sca_alt in self.dlos_list: #for every altitude between minimum and maximum scattering altitude

							print("############### DEBUG dlos:", )
							sca_from_E_is_visible = True
							if self.use_elevation_map:
								sca_from_E_is_visible = self.alt_map.IsVisible(sca_lon, sca_lat, alt1=sca_alt, lon2=lon, lat2=lat, dlos = 0.5)

							if sca_from_E_is_visible:

								I0, w_DoLP = self.ComputeSingleRSGroundPointSource(ilon, ilat, a_rd, e_rd, sca_alt)

								self.ground_map.total_scattering_map[ie_pc, ia_pc, ilat, ilon] += I0 # Intensity of light scattered from a given (e, a)
								self.ground_map.scattering_map[ie_pc, ia_pc, ilat, ilon] += (w_DoLP * I0) # Intensity of polarized light scattered from a given (e, a). Given a point source, the AoRD is the same for every altitude -> the addition makes sens

							# count += 1
							# self.Progress(count, N, suffix="of ground point sources done")
							pbar.update(1)


						# self.ground_map.DoLP_map = self.ground_map.scattering_map / ( self.ground_map.total_scattering_map + self.ground_map.scattering_map) # DoLP of a given (e, a)
						self.ground_map.DoLP_map = self.ground_map.scattering_map / self.ground_map.total_scattering_map # DoLP of a given (e, a)

						self.ground_map.AoRD_map[ie_pc, ia_pc, ilat, ilon] = self.obs.GetRayleighAngle(a_rd, 0) # e_rd is NOT the elevation if ground emission



	def ComputeSingleRSGroundPointSource(self, ilon, ilat, a_rd, e_rd, alt):

		I0 = self.ground_map.I_map[ilat, ilon]

		AR, RE, RD_angle = self.GetGeometryFromAzDist(a_rd, e_rd, alt)
		print("DEBUG: AR, RE, RD_angle, alt", AR, RE, RD_angle*RtoD, alt)

		print("DEBUG I0 1:", I0, I0*self.ground_map.GetArea(ilat) / RE ** 2)
		I0 *= self.ground_map.GetArea(ilat) / RE ** 2

		if alt != 0:
			opt_depth = self.atmosphere.GetRSOpticalDepth(self.wavelength, 0, alt) * RE / alt
			print("DEBUG opt_depth: ER", opt_depth, 1-opt_depth, np.exp(-opt_depth))
		else:
			opt_depth = 0
		I0 *= np.exp(-opt_depth)

		print("DEBUG I0 2:", I0)
		print("DEBUG scattered:", self.GetScattered(I0, AR, RE, RD_angle, alt)[0]/I0)

		I0, w_DoLP = self.GetScattered(I0, AR, RE, RD_angle, alt)

		print("DEBUG I0 3:", I0)

		if alt != 0:
			opt_depth = self.atmosphere.GetRSOpticalDepth(self.wavelength, 0, alt) * AR / alt
			print("DEBUG opt_depth: RA", opt_depth, 1-opt_depth, np.exp(-opt_depth))
		else:
			opt_depth = 0

		I0 *= np.exp(-opt_depth)
		print("DEBUG I0 4:", I0)

		return I0, w_DoLP


	def ComputeSingleRSSkyPointSource(self, time, ia, a, ie, e, alt):
		I0 = self.sky_map.cube[time, ie, ia]
		I0 *= np.cos(e) # Bining effect, bin around e=90 are "smaller" if sky emission

		AR, RE, RD_angle = self.GetGeometryFromAzEl(a, e, alt)

		if alt != 0:
			opt_depth = self.atmosphere.GetRSOpticalDepth(self.wavelength, self.sky_map.h, alt) * RE / alt
		else:
			opt_depth = 0

		I0 *= np.exp(-opt_depth)


		I0, w_DoLP = self.GetScattered(I0, AR, RE, RD_angle, alt, elevation = e)

		if alt != 0:
			opt_depth = self.atmosphere.GetRSOpticalDepth(self.wavelength, 0, alt) * AR / alt
		else:
			opt_depth = 0
		# print("opt_depth, p.exp(-opt_depth)", opt_depth, np.exp(-opt_depth))
		I0 *= np.exp(-opt_depth)

		return I0, w_DoLP

	def GetGeometryFromAzEl(self, a_rd, e_rd, alt):
		"""From known geometry parameters : a, e of emission and alt of scatering, return missing parameters: Distance between emission and scattering and angle of scattering.
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

		# e_E = np.arcsin(abs(alt - self.sky_map.h) / RE) #Angle of ER segment with horizontal plane (E is at a given altitude, AE is NOT horyzontal)

		return AR, RE, ARE #, e_E

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

		# AER = np.pi - ARE - self.e_pc #Angle of ER segment with horizontal plane (E is on the ground here so AE is horizontal)

		return AR, RE, ARE#, AER

	def GetScattered(self, I0, AR, ER, RD_angle, alt, elevation = 0):
		"""Given an initial intensity of a source and some geometrical parameter, returns the intensity mesured at the instrument and its DoLP.
		Input parameters: elevation, altitude of scattering, scattering angle, distance between emission and scattering."""

		V = self.atmosphere.GetVolume(AR, self.ouv_pc, unit="km") #in m3
		Crs = self.atmosphere.GetRSVolumeCS(self.wavelength, alt) #in cm-1 * 100 = m-1
		P = self.atmosphere.GetRSPhaseFunction(self.wavelength, RD_angle)

		# print("DEBUG V, Crs, P", V, Crs, P)

		if AR != 0:
			I0 *= Crs * P * V * self.PTCU_area / (AR * 1) ** 2 / 4 / np.pi
			# print("alt, V, Crs, P, omega", alt, V, Crs, P, self.PTCU_area / (AR*1000) ** 2)
		else:
			I0 = 0
			# print("WARNING!!!! I0==0")

		DoLP = np.sin(RD_angle)**2 / (1 + np.cos(RD_angle)**2) # DoLP dependance on scattering angle

		# print("DEBUG DOLP:", w_DoLP * 100)
		return I0, DoLP


	def GetDirect(self, t, ia, ie):
		"""Returns direct intensity (flux in nW) coming in the instrument"""
		a = self.a_pc_list[ia]
		e = self.e_pc_list[ie]

		return self.sky_map.GetFlux(a, e, self.ouv_pc, t = t, area = self.PTCU_area)
