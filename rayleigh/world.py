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
			# self.a_pc_list = np.arange(179.5, 180.5, 0.005) * DtoR
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

		self.max_angle_discretization = in_dict["max_angle_discretization"]
		if self.max_angle_discretization != "None":
			self.max_angle_discretization = float(self.max_angle_discretization) * DtoR
		else:
			self.max_angle_discretization = None

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

		#Show the skymaps if they change in time or we move the instrument
		# if self.has_sky_emission and (self.is_time_dependant or not self.is_single_observation):
		# self.sky_map.MakeSkyCubePlot(self.a_pc_list, self.e_pc_list, self.ouv_pc)
		# plt.show()
		# plt.close("all")

		print("Has_sky_emission:", self.has_sky_emission)

		### Initialize the groundmap. Even if we don't use it, some parameters must be initialized
		self.ground_map = GroundMap(in_dict, self.Nb_a_pc, self.Nb_e_pc)

		self.has_ground_emission = self.ground_map.exist
		print("Has_ground_emission:", self.has_ground_emission)


		### Initialize the atlitude maps
		self.alt_map = ElevationMap(in_dict)
		self.use_elevation_map = bool(int(in_dict["use_alt_map"]))
		# self.alt_map.PlotMap()

		self.instrument_altitude = self.alt_map.GetAltitudeFromLonLat(self.ground_map.A_lon, self.ground_map.A_lat) / 1000. #In km


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
			dlos_list = np.array([90. / np.sin(e_pc)])
			# dlos_list = [(range_min + range_max) / 2.]

		# print("DEBUG dlos list:", dlos_list)

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
	# def Progress(self, iteration, total, prefix = '', suffix = '', decimals = 1, length = 50, fill = '█', printEnd = "\r"):
	# 	"""
	# 	Call in a loop to create terminal progress bar
	# 	@params:
	# 	    iteration   - Required  : current iteration (Int)
	# 	    total       - Required  : total iterations (Int)
	# 	    prefix      - Optional  : prefix string (Str)
	# 	    suffix      - Optional  : suffix string (Str)
	# 	    decimals    - Optional  : positive number of decimals in percent complete (Int)
	# 	    length      - Optional  : character length of bar (Int)
	# 	    fill        - Optional  : bar fill character (Str)
	# 	    printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
	# 	"""
	# 	percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
	# 	filledLength = int(length * iteration // total)
	# 	bar = fill * filledLength + '-' * (length - filledLength)
	#
	# 	sys.stdout.write('\r%s |%s| %s%% %s %s' % (prefix, bar, percent, suffix, printEnd))
	# 	sys.stdout.flush()
	# 	# Print New Line on Complete
	# 	if iteration == total:
	# 		print("\n")

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
								V, Vcos, Vsin = self.ComputeSingleRSSkyPointSource(time, ia, a, ie, e, alt)
								self.sky_map.V_map[time, ie_pc, ia_pc, ie, ia] += V
								self.sky_map.Vcos_map[time, ie_pc, ia_pc, ie, ia] += Vcos
								self.sky_map.Vsin_map[time, ie_pc, ia_pc, ie, ia] += Vsin
								# self.sky_map.total_scattering_map[time, ie_pc, ia_pc, ie, ia] += I0 # Intensity of light scattered from a given (e, a)
								# self.sky_map.scattering_map[time, ie_pc, ia_pc, ie, ia] += DoLP * I0 # Intensity of polarized light scattered from a given (e, a). Given a point source, the AoRD is the same for every altitude -> the addition makes sens

							pbar.update(1)

					# 	self.sky_map.DoLP_map = self.sky_map.scattering_map / self.sky_map.total_scattering_map # DoLP of a given (e, a)
					#
					# self.sky_map.AoRD_map[time, ie_pc, ia_pc, ie, ia] = AoLP

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

		N = len(self.ground_map.mid_longitudes) * len(self.ground_map.mid_latitudes) * self.Nalt
		print(N, "bins to compute")
		with tqdm(total=N, file=sys.stdout) as pbar:
			for ilat, lat in enumerate(self.ground_map.mid_latitudes):
				for ilon, lon in enumerate(self.ground_map.mid_longitudes):
					if self.ground_map.I_map[ilat, ilon] > 0:
						a_rd, e_rd = LonLatToAzDist(lon, lat, self.ground_map.A_lon, self.ground_map.A_lat)
						for sca_lon, sca_lat, sca_alt in self.dlos_list: #for every altitude between minimum and maximum scattering altitude

							sca_from_E_is_visible = True
							if self.use_elevation_map:
								sca_from_E_is_visible = self.alt_map.IsVisible(sca_lon, sca_lat, alt1=sca_alt, lon2=lon, lat2=lat, dlos = 0.5)

							if sca_from_E_is_visible:

								V, Vcos, Vsin = self.ComputeSingleRSGroundPointSource(ilon, ilat, a_rd, e_rd, sca_alt)

								self.ground_map.V_map[ie_pc, ia_pc, ilat, ilon] += V
								self.ground_map.Vcos_map[ie_pc, ia_pc, ilat, ilon] += Vcos
								self.ground_map.Vsin_map[ie_pc, ia_pc, ilat, ilon] += Vsin

								# self.ground_map.total_scattering_map[ie_pc, ia_pc, ilat, ilon] += I0 # Intensity of light scattered from a given (e, a)
								# self.ground_map.scattering_map[ie_pc, ia_pc, ilat, ilon] += (w_DoLP * I0) # Intensity of polarized light scattered from a given (e, a). Given a point source, the AoRD is the same for every altitude -> the addition makes sens

							# count += 1
							# self.Progress(count, N, suffix="of ground point sources done")
							pbar.update(1)


						# self.ground_map.DoLP_map = self.ground_map.scattering_map / self.ground_map.total_scattering_map # DoLP of a given (e, a)
						#
						# self.ground_map.AoRD_map[ie_pc, ia_pc, ilat, ilon] = AoLP


	def CutVolume(self, ER, V, max_angle = None):

		if not max_angle or max_angle == 0:
			return 1

		condition = max_angle

		a = np.arctan2(np.cbrt(V), ER)

		N = (4 * int(a / condition)) ** 2 + 1

		# print("Nb Volumes:", N)
		return N


	def ComputeSingleRSGroundPointSource(self, ilon, ilat, a_E, e_E, alt):

		I0 = self.ground_map.I_map[ilat, ilon]

		AR, RE, RD_angle = self.GetGeometryFromAzDist(a_E, e_E, alt)
		vol = self.atmosphere.GetVolume(AR, self.ouv_pc, unit="km", type="pole")

		cut_V = self.CutVolume(RE, vol, max_angle = self.max_angle_discretization)

		# print("DEBUG: AR, RE, RD_angle, alt", AR, RE, RD_angle*RtoD, alt)

		geo_list = self.GetDiscreteRSPoint(self.a_pc, self.e_pc, alt, AR, RE, a_E, e_E, "ground", N=cut_V)

		V_T, Vcos_T, Vsin_T = 0, 0, 0

		vol_T = 0

		for a_A, e_A, r, AR, RE, RD_angle, alt, AoLP, dvol, dr in geo_list:

			# print(a_A*RtoD, e_A*RtoD, r, AR, RE, RD_angle*RtoD, alt, AoLP*RtoD, da*RtoD, dr)

			# print("DEBUG I0 1:", I0, I0*self.ground_map.GetArea(ilat) / RE ** 2)
			I0 = self.ground_map.I_map[ilat, ilon]
			I0 *= self.ground_map.GetArea(ilat) / RE ** 2

			if alt != 0:
				opt_depth = self.atmosphere.GetRSOpticalDepth(self.wavelength, 0, alt) * RE / alt
				# print("DEBUG opt_depth: ER", opt_depth, 1-opt_depth, np.exp(-opt_depth))
			else:
				opt_depth = 0

			I0 *= np.exp(-opt_depth)

			# print("DEBUG I0 2:", I0)
			# print("DEBUG scattered:", self.GetScattered(I0, AR, RE, RD_angle, alt)[0]/I0)

			# Vol = self.atmosphere.GetVolume(AR, da, dr=dr, unit="km") #in km3
			vol_T += dvol
			Crs = self.atmosphere.GetRSVolumeCS(self.wavelength, alt) #in km-1
			P = self.atmosphere.GetRSPhaseFunction(self.wavelength, RD_angle)

			# print("DEBUG V, Crs, P", V, Crs, P)

			if AR != 0:
				I0 *= Crs * P * dvol * self.PTCU_area / AR ** 2 / 4 / np.pi
				# print("alt, V, Crs, P, omega", alt, V, Crs, P, self.PTCU_area / (AR*1000) ** 2)
			else:
				I0 = 0
				# print("WARNING!!!! I0==0")

			DoLP = np.sin(RD_angle)**2 / (1 + np.cos(RD_angle)**2) # DoLP dependance on scattering angle


			# print("DEBUG I0 3:", I0)

			if alt != 0:
				opt_depth = self.atmosphere.GetRSOpticalDepth(self.wavelength, 0, alt) * AR / alt
				# print("DEBUG opt_depth: RA", opt_depth, 1-opt_depth, np.exp(-opt_depth))
			else:
				opt_depth = 0

			I0 *= np.exp(-opt_depth)
			# print("DEBUG I0 4:", I0)

			V, Vcos, Vsin = self.GetVParamFromLightParam(I0, DoLP, AoLP)

			V_T += V
			Vcos_T += Vcos
			Vsin_T += Vsin

		# print(V_T, Vcos_T, Vsin_T)
		if not abs(vol_T - vol) < 10**-10:
			print(len(geo_list), vol_T, vol, vol_T /vol, abs(vol_T - vol) < 10**-10)

		# V_T *= vol / vol_T
		# Vcos_T *= vol / vol_T
		# Vsin_T *= vol / vol_T

		return V_T, Vcos_T, Vsin_T

	def ComputeSingleRSSkyPointSource(self, time, ia_E, a_E, ie_E, e_E, alt):
		"""t: time of the map
		ia_E, a_E, ie_E, e_E: az and el and their indexes of the emission point E."""

		I0 = self.sky_map.cube[time, ie_E, ia_E]

		AR, RE, RD_angle = self.GetGeometryFromAzEl(a_E, e_E, alt)
		vol = self.atmosphere.GetVolume(AR, self.ouv_pc, unit="km", type="pole")

		cut_V = self.CutVolume(RE, vol, max_angle = self.max_angle_discretization)

		#Get a list of all small scattering volumes in the big troncated cone. 1 point if distance from emission to cone is 10 km
		geo_list = self.GetDiscreteRSPoint(self.a_pc, self.e_pc, alt, AR, RE, a_E, e_E, "sky", N = cut_V)

		V_T, Vcos_T, Vsin_T = 0, 0, 0

		vol_T = 0

		for a_A, e_A, r, AR, RE, RD_angle, alt, AoLP, dvol, dr in geo_list:

			# print(a_A*RtoD, e_A*RtoD, r, AR, RE, RD_angle*RtoD, alt, AoLP*RtoD, da*RtoD, dr)

			I0 = self.sky_map.cube[time, ie_E, ia_E]

			I0 *= self.sky_map.GetPixelArea(ie_E) / RE ** 2

			if alt != 0:
				opt_depth = self.atmosphere.GetRSOpticalDepth(self.wavelength, self.sky_map.h, alt) * RE / alt
			else:
				opt_depth = 0

			I0 *= np.exp(-opt_depth)

			# Vol = self.atmosphere.GetVolume(AR, da, dr=dr, unit="km") #in km3
			vol_T += dvol
			Crs = self.atmosphere.GetRSVolumeCS(self.wavelength, alt) #in km-1
			P = self.atmosphere.GetRSPhaseFunction(self.wavelength, RD_angle)

			# print("DEBUG V, Crs, P", V, Crs, P)

			if AR != 0:
				I0 *= Crs * P * dvol * self.PTCU_area / AR ** 2 / 4 / np.pi
				# print("alt, V, Crs, P, omega", alt, V, Crs, P, self.PTCU_area / (AR*1000) ** 2)
			else:
				I0 = 0
				# print("WARNING!!!! I0==0")

			DoLP = np.sin(RD_angle)**2 / (1 + np.cos(RD_angle)**2) # DoLP dependance on scattering angle


			if alt != 0:
				opt_depth = self.atmosphere.GetRSOpticalDepth(self.wavelength, 0, alt) * AR / alt
			else:
				opt_depth = 0
			# print("opt_depth, p.exp(-opt_depth)", opt_depth, np.exp(-opt_depth))
			I0 *= np.exp(-opt_depth)

			V, Vcos, Vsin = self.GetVParamFromLightParam(I0, DoLP, AoLP)

			V_T += V
			Vcos_T += Vcos
			Vsin_T += Vsin

		# print(V, Vcos, Vsin, I0, DoLP, AoLP)
		return V_T, Vcos_T, Vsin_T

	def GetDiscreteRSPoint(self, a, e, alt, AR, RE, a_E, e_E, map_type, N = 100):
		"""For a given observation direction and emission point, returns N points of the scattering volume (troncated cone).
		For each returned point:
		a, e: azimut elevation of the subvolume (close to a_pc, e_pc)
		r: range of the subvolume (distance to the instrument in km)
		AR, RE, RD_ange: distances (km) and angle
		alt: altitude of the subvolume in km
		AoLP: AoLP of the diffused light"""

		### Coordinates in the cylindrical reference frame.
		if N > 1:
			n = int(np.cbrt(N)) + 1
			ap_list, dap = np.linspace(0, 2 * np.pi, n, retstep=True, endpoint=True)
			ep_list, dep = np.linspace(0, self.ouv_pc, n, retstep=True, endpoint=True)
			r_list, dr = np.linspace(AR - self.atmosphere.d_los/2., AR + self.atmosphere.d_los/2., n, retstep=True)

			ap_list_mid = np.array([ia + dap / 2. for ia in ap_list[:-1]])
			ep_list_mid = np.array([ie + dep / 2. for ie in ep_list[:-1]])
			r_list_mid = np.array([ir + dr / 2. for ir in r_list[:-1]])
		else:
			n = 1
			ap_list, dap = np.array([0, 2 * np.pi]), 2 * np.pi
			ep_list, dep = np.array([0, self.ouv_pc]), self.ouv_pc
			r_list, dr = np.array([AR]), self.atmosphere.d_los

			ap_list_mid = np.array([0])
			ep_list_mid = np.array([0])
			r_list_mid = np.array([AR])


		### Returned list of interesting numbers for each subvolumes
		geo_list = np.zeros((len(ap_list_mid) * len(ep_list_mid) * len(r_list), 10))

		i = 0

		### Rotation matrix from instrument ref frame to UpEastNorth
		ca, sa = np.cos(a), np.sin(a)
		ce, se = np.cos(e), np.sin(e)
		Ria = np.array([	[se, 		0,		ce],
							[	 ce * sa,	-ca,	-se * sa],
							[	 ce * ca,	sa,		-se * ca]])


		# print("ap_list")
		# print(ap_list* RtoD)

		###Corrdinates in cylindrical of emission point
		E_u, E_e, E_n = np.sin(e_E), np.cos(e_E) * np.sin(a_E), np.cos(e_E) * np.cos(a_E)
		E_xI, E_yI, E_zI = np.dot(Ria.transpose(), (E_u, E_e, E_n))
		E_ap = np.arctan2(E_yI, E_zI)
		# print(E_ap * RtoD)
		ap_list = (ap_list + E_ap) % (2*np.pi)
		ap_list_mid = (ap_list_mid + E_ap) % (2*np.pi)
		# print(ap_list* RtoD)


		### For each subvolume
		for iep, ep in enumerate(ep_list_mid):
			cep, sep = np.cos(ep), np.sin(ep)
			for iap, ap in enumerate(ap_list_mid):

				cap, sap = np.cos(ap), np.sin(ap)

				xI, yI, zI = cep, sap * sep, cap * sep

				up, east, north = np.dot(Ria, ((xI), (yI), (zI)))

				az = np.arctan2(east, north)
				if az < 0:
					az += 2 *np.pi
				el = np.arcsin(up)

				for r in r_list_mid:
					dvol = self.atmosphere.GetVolume(r, [ep_list[iep], ep_list[iep + 1]], da = dap, dr=dr, unit="km", type="slice")

					obs = ObservationPoint(self.ground_map.A_lon, self.ground_map.A_lat, self.sky_map.h, az, el, A_alt = self.instrument_altitude)

					alt = obs.GetPCoordinatesFromRange(r)[2]
					AoLP = obs.GetRayleighAngle(a_E, e_E)

					if map_type == "ground":
						AR, RE, RD_angle = self.GetGeometryFromAzDist(a_E, e_E, alt, v_pc_u = Getuen(az, el), obs=obs)
					else:
						AR, RE, RD_angle = self.GetGeometryFromAzEl(a_E, e_E, alt, v_pc_u = Getuen(az, el), obs=obs)

					geo_list[i] = az, el, r, AR, RE, RD_angle, alt, AoLP, dvol, dr
					i += 1
		return geo_list

	def GetVParamFromLightParam(self, I0, DoLP, AoLP):
		V = I0 / 2.
		Vcos = I0 * DoLP * np.cos(2 * AoLP) / 4.
		Vsin = I0 * DoLP * np.sin(2 * AoLP) / 4.
		return V, Vcos, Vsin

	def MakeAoRDHist(self, A, I=None, D=None, Ipola=None, Inonpola=None):
		if Ipola is None and (I is not None and D is not None):
			Ipola = I * D
		if Inonpola is None and (I is not None and (D is not None or Ipola is not None)):
			Inonpola = I - Ipola

		###Compute the AoLP contribution histogram.
		N_bins = 180
		#bins is the bins limits, so 181 values
		bins, width = np.linspace(-np.pi/2, np.pi/2, N_bins + 1, endpoint=True, retstep=True)
		bins, width = np.linspace(-np.pi/2 - width/2, np.pi/2 + width/2, N_bins + 1, endpoint=True, retstep=True)
		mid_bins = [(bins[i+1] + bins[i])/2. for i in range(len(bins)-1)]

		#One value for each bin = 180 values
		hst = np.zeros(N_bins)
		# print("DEBUG ROT", self.N_bins, len(self.bins), len(self.mid_bins), len(self.hst))

		###Density=False: hst bins units are similar to intensities. Sum of hst does not depends on N_bins
		hst, bins = np.histogram(A, bins = bins, weights = Ipola, density = False)

		return hst, mid_bins

	def DrawAoLPHist(self, hst, bins, I0, DoLP, AoRD, double=True, save=False, save_name=None, show=True):
		"""Make an pyplot histogram of all AoLP contributionsThe histogram is calculated in GetLightParameters()
		Need a call to plt.show() after calling this function."""
		f3, ax = plt.subplots(1, figsize=(16, 8))

		ax = plt.subplot(111, projection='polar')
		ax.set_theta_zero_location("N")
		ax.set_theta_direction(-1)
		if double:
			ax.set_thetamin(-180)
			ax.set_thetamax(180)
		else:
			ax.set_thetamin(-90)
			ax.set_thetamax(90)

		# I0, DoLP, AoRD = self.I_list[0,0,0], self.DoLP_list[0,0,0], self.AoRD_list[0,0,0]
		# I0, DoLP, AoRD = self.GetLightParameters(ground=ground, sky=sky)

		N_bins = len(hst)
		width = bins[1] - bins[0]

		if not double:
			# bins = bins[:N_bins]
			h = hst
		else:
			bins = np.append(bins[:N_bins], 180 * DtoR + np.array(bins))
			# bins = np.append(bins[:N_bins], 180 * DtoR + np.array(bins)[:N_bins])
			h = np.append(hst, hst)

		bars = ax.bar(bins, h, width=width)

		ax.set_title("Weighted AoRD: I0 = " + str(np.format_float_scientific(I0, precision=3)) + " DoLP = " + str(np.round(DoLP, 1)) + " AoRD = " + str(np.round(AoRD*RtoD, 1)))

		if not double:
			ax.plot([AoRD, AoRD], [0, max(hst)], "r")
		else:
			ax.plot([AoRD, AoRD+np.pi], [max(hst), max(hst)], "r")

		# if show:
		# 	plt.show()

		if save and save_name is not None:
			plt.savefig(save_name, bbox_inches='tight')


	def LockInAmplifier(self, hst, mid_bins, InonPola):
		###Simulate the instrument with a rotating polarizing filter to get V, Vcos, Vsin and then I, DoLP, AoLP
		Ns = 1000 #5 * len(self.bins)
		filter_orientation, df = np.linspace(0, 2 * np.pi, Ns, endpoint=False, retstep=True)

		rs_signal = np.zeros(Ns) #+ 0.5 * self.InonPola * len(self.hst) * self.width * df

		N_bins = len(hst)

		for i_f, f in enumerate(filter_orientation):
			for ihist, hist in enumerate(hst):
				theta = mid_bins[ihist]
				# theta = b[(ihist+1)%len(self.hst)] - b[ihist]
				rs_signal[i_f] += (hist * np.cos(theta - f) ** 2 + 0.5 * InonPola / N_bins)
				# rs_signal[i_f] += hist * np.cos(b[ihist] - f) ** 2 # + 0.5 * self.InonPola

		# plt.plot(filter_orientation, rs_signal)
		# plt.plot(filter_orientation, [0.5 * self.InonPola * len(self.hst)]*Ns)
		# plt.show()

		V = np.average(rs_signal)
		Vcos = np.average(rs_signal * np.cos(2 * filter_orientation))
		Vsin = np.average(rs_signal * np.sin(2 * filter_orientation))

		# print("DEBUG I0:", 2*self.V, self.I0, abs(2*self.V-self.I0)<10**-15)

		# print("V, Vcos, Vsin", self.V, self.Vcos, self.Vsin)

		I0 = 2 * V
		DoLP = 100 * 2 * np.sqrt(Vcos ** 2 + Vsin ** 2) / V #in %
		AoRD = np.arctan2(Vsin, Vcos) / 2

		# print("DEBUG LIGHT PARAM:", 2 * self.V, self.I0, self.DoLP, self.AoRD)
	# print("TIME:", tm.time() - tmp_start)

		return V, Vcos, Vsin, I0, DoLP, AoRD




	def GetGeometryFromAzEl(self, a_E, e_E, alt, v_pc_u = None, obs = None):
		"""From known geometry parameters : a_E, e_E of emission point E and alt of scattering, return missing parameters: Distance between emission and scattering and angle of scattering.
		A instrument positions pointing to (a_pc, e_pc)
		E emission points (a_E, e_E, h)
		R rayleigh diffusion point (a_pc, e_pc, atmosphere.h_r)"""

		if obs is None:
			obs = self.obs

		e_E = max(e_E, 10**(-30)) # Avoid weird behaviour and division by 0

		v_rd_u = Getuen(a_E, e_E) # Vector from observer to emission point
		_, AE = obs.GetAH(elevation = e_E, azimut = a_E, altitude = self.sky_map.h)  # Distance between observer and emisson point

		if v_pc_u is None:
			v_pc_u = self.v_pc_u

		RAE = GetAngle(v_pc_u, v_rd_u) # Angle between line of sight and emission

		_, AR = obs.GetAH(altitude = alt) # Distance between observer and scattering point

		RE = GenPythagore(AR, AE, RAE) # Distance between emission and scaterring point
		ARE = np.arcsin(AE * np.sin(RAE) / RE) # Scattering angle

		# e_E = np.arcsin(abs(alt - self.sky_map.h) / RE) #Angle of ER segment with horizontal plane (E is at a given altitude, AE is NOT horyzontal)

		return AR, RE, ARE #, e_E

	def GetGeometryFromAzDist(self, a_rd, d, h_r, v_pc_u=None, obs = None):
		"""From known geometry parameters : a, d, alt of emission and alt of scatering , return missing parameters: Distance between emission and scattering and angle of scattering.
		A instrument positions pointing to (a_pc, e_pc)
		E emission points (a_rd, e_rd, h)
		R rayleigh diffusion point (a_pc, e_pc, h_r)"""

		if obs is None:
			obs = self.obs

		d = max(d, 10**(-30)) # Avoid weird behaviour and division by 0

		v_rd_u = Getuen(a_rd, 0.) # Vector from observer to emission point (flat ground, emission at elevation 0)
		AE = d * RT # Distance between observer and emisison point. (exact, following the curvature (~flat earth))

		if v_pc_u is None:
			v_pc_u = self.v_pc_u
		RAE = GetAngle(v_pc_u, v_rd_u) # Angle between line of sight and emission

		_, AR = obs.GetAH(altitude = h_r) # Distance between observer and scattering point

		RE = GenPythagore(AR, AE, RAE) # Distance between emission and scaterring point
		ARE = np.arcsin(AE * np.sin(RAE) / RE) # Scattering angle

		# AER = np.pi - ARE - self.e_pc #Angle of ER segment with horizontal plane (E is on the ground here so AE is horizontal)

		return AR, RE, ARE#, AER

	def GetScattered(self, I0, AR, RD_angle, alt):
		"""Given an initial intensity of a source and some geometrical parameter, returns the intensity mesured at the instrument and its DoLP.
		Input parameters: elevation, altitude of scattering, scattering angle, distance between emission and scattering."""

		V = self.atmosphere.GetVolume(AR, self.ouv_pc, unit="km") #in km3
		Crs = self.atmosphere.GetRSVolumeCS(self.wavelength, alt) #in km-1
		P = self.atmosphere.GetRSPhaseFunction(self.wavelength, RD_angle)

		# print("DEBUG V, Crs, P", V, Crs, P)

		if AR != 0:
			I0 *= Crs * P * V * self.PTCU_area / AR ** 2 / 4 / np.pi
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

		if self.has_sky_emission:
			opt_depth = self.atmosphere.GetRSOpticalDepth(self.wavelength, self.sky_map.h, 0) / np.sin(e)
			direct = self.sky_map.GetFlux(a, e, self.ouv_pc, t = t, area = self.PTCU_area)
			print("opt_depth, direct:")
			print(opt_depth, direct)
			return direct * np.exp(-opt_depth)
		else:
			return 0
