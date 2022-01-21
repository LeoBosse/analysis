#!/usr/bin/python3
# -*-coding:utf-8 -*

from mpi4py import MPI
mpi_comm = MPI.COMM_WORLD
mpi_rank = mpi_comm.Get_rank()
mpi_size = mpi_comm.Get_size()
mpi_name = mpi_comm.Get_name()

import sys as sys
import numpy as np
import scipy.integrate as integrate
import time as tm
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from matplotlib.patches import Circle
from matplotlib.patches import Ellipse
from matplotlib.patches import Arrow
# from matplotlib.lines import mlines

mpl.rcParams.update({'font.size': 30})

import astropy as apy
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from astropy import units as u

import chaosmagpy as chaos

# import osgeo.gdal as gdal
# gdal.UseExceptions()  # not required, but a good idea

import imageio

from tqdm import tqdm

from pomerol_configuration import *
from observation import *
from rayleigh_utils import *
from sky_map import *
from ground_map import *
from atmosphere import *


class World:
	def __init__(self, in_dict):
		self.wavelength = float(in_dict["wavelength"])
		self.src_path = in_dict["src_path"]

		self.show_ground_albedo = False
		self.show_sky_cube = False
		self.show_altitude = False
		if "show_ground_albedo" in in_dict:
			self.show_ground_albedo = bool(int(in_dict["show_ground_albedo"]))
		if "show_sky_cube" in in_dict:
			self.show_sky_cube = bool(int(in_dict["show_sky_cube"]))
		if "show_altitude" in in_dict:
			self.show_altitude = bool(int(in_dict["show_altitude"]))

		# Angle de demi-ouverture de l'instrument (1 deg pour les crus)
		self.ouv_pc		= float(in_dict["instrument_opening_angle"]) * DtoR
		self.inst_solid_angle = 2 * np.pi * (1 - np.cos(self.ouv_pc))
		# Aire de l'instrument donnée en cm^2, transformée en km^2
		self.PTCU_area	= float(in_dict["instrument_area"]) * 10 ** (-10) #in km**2

		self.use_analytic = in_dict["use_analytic"].lower() == "true"

		self.SetFluxUnit 		= lambda f: f

		self.flux_unit 			= in_dict["flux_unit"]
		self.eye_efficiency_data = pd.read_csv(self.src_path + "linCIE2008v2e_fine.csv", names=["wl", "func"])
		self.eye_efficiency		= np.interp(self.wavelength, self.eye_efficiency_data["wl"], self.eye_efficiency_data["func"])
		self.WattToLumen 		= lambda f: f * 683.002 * self.eye_efficiency
		self.LumenToCandela 	= lambda f: f / self.inst_solid_angle
		self.CandelaToCdm2 		= lambda f: f / (self.PTCU_area * 1e6)
		self.LumenToCdm2 		= lambda f: f / self.inst_solid_angle / (self.PTCU_area * 1e6)
		self.WattToCdm2 		= lambda f: self.LumenToCdm2(self.WattToLumen(f))
		self.WattToLambert 		= lambda f: (np.pi / 1e4) * self.WattToCdm2(f)
		self.LambertToCdm2 		= lambda f: (1e4 / np.pi) * f
		self.Cdm2ToWatt 		= lambda f: f * self.inst_solid_angle * (self.PTCU_area * 1e6) / 683.002 / self.eye_efficiency
		self.LambertToWatt 		= lambda f: self.Cdm2ToWatt(self.LambertToCdm2(f))
		self.LambertToMag 		= lambda f: 26.331 - 1.0857 * np.log(f)
		self.WattToMag 			= lambda f: self.LambertToMag(self.WattToLambert(f))
		self.S10TonanoWatt		= lambda f: self.LambertToWatt(0.22 * f) #1.962e-6 * f #

		# print(self.eye_efficiency, self.S10TonanoWatt(1))


		if "W" not in self.flux_unit:
			if self.flux_unit == "ph":
				self.SetFluxUnit = lambda f: f * self.wavelength * 1e-18 / (6.62607015e-34 * 299792458) # f[nW] * lambda[m]/(hc) * 1e-9[W/nW]
			# elif self.flux_unit == "mag":
			# 	self.SetFluxUnit = lambda f: 41.438 - 1.0857 * np.log(f * self.wavelength * 1e-18 / (6.62607015e-34 * 299792458)) # f[nW] * lambda[m]/(hc) * 1e-9[W/nW]

			elif self.flux_unit in ["lm", "cd", "cd/m2", "mag"]: #unit == lm or cd*
				data = pd.read_csv(self.src_path + "/linCIE2008v2e_fine.csv", names=["wl", "func"])
				if self.flux_unit == "lm":
					self.SetFluxUnit = lambda f: f * 683.002 * np.interp(self.wavelength, data["wl"], data["func"])
					self.flux_unit = "nlm"

				elif self.flux_unit == "cd":
					self.SetFluxUnit = lambda f: f * 683.002 * np.interp(self.wavelength, data["wl"], data["func"]) / self.inst_solid_angle
					self.flux_unit = "ncd"

				elif self.flux_unit == "cd/m2":
					print("EYE DATA",np.interp(self.wavelength, data["wl"], data["func"]))
					self.SetFluxUnit = lambda f: f * 683.002 * np.interp(self.wavelength, data["wl"], data["func"]) / self.inst_solid_angle / (self.PTCU_area * 1e6)
					self.flux_unit = "ncd/m2"

				elif self.flux_unit == "mag":
					self.SetFluxUnit = lambda f: 26.331 - 1.0857 * np.log(f * 683.002 * np.interp(self.wavelength, data["wl"], data["func"]) / self.inst_solid_angle / (self.PTCU_area * 1e6) * np.pi * 1e-4) # f[nW] * lambda[m]/(hc) * 1e-9[W/nW]

		else:
			self.flux_unit = "nW"



		### Get all the azimuts and elevations of the instruments.
		self.a_pc_list		= np.array(in_dict["azimuts"].split(";"))  # List of azimuts for the observations
		if self.a_pc_list[0] == "all":
			# self.a_pc_list = np.arange(179.5, 180.5, 0.005) * DtoR
			# self.a_pc_list = np.arange(0, 360, 10) * DtoR
			self.a_pc_list = np.arange(0, 360, 5) * DtoR
			# self.a_pc_list = np.arange(-10, 10, 1) * DtoR
		elif len(self.a_pc_list) == 3:
			self.a_pc_list = np.linspace(float(self.a_pc_list[0]), float(self.a_pc_list[1]), int(self.a_pc_list[2])) * DtoR
		else:
			self.a_pc_list = np.array([float(a) for a in self.a_pc_list]) * DtoR

		self.e_pc_list		= np.array(in_dict["elevations"].split(";"))  # List of elevation for the observations
		if self.e_pc_list[0] == "all":
			self.e_pc_list = np.linspace(1, 90, 20, endpoint=True) * DtoR
			# self.e_pc_list = np.linspace(35, 55, 21, endpoint=True) * DtoR
		elif len(self.e_pc_list) == 3:
			self.e_pc_list = np.linspace(float(self.e_pc_list[0]), float(self.e_pc_list[1]), int(self.e_pc_list[2])) * DtoR
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

		### Initialize the skymap. Even if we don't use it, some parameters must be initialized
		self.sky_map = SkyMap(in_dict, self.Nb_a_pc, self.Nb_e_pc)

		self.has_sky_emission = self.sky_map.exist
		self.is_time_dependant = self.sky_map.exist and self.sky_map.is_time_dependant

		self.start_datetime = dt.datetime.strptime(in_dict["start_datetime"], "%Y%m%d-%H:%M:%S")
		self.end_datetime = in_dict["end_datetime"].lower()
		if self.end_datetime in ["none", "0", "false"]:
			self.end_datetime = self.start_datetime
			self.times = [self.start_datetime]
		elif self.sky_map.mode == "movie":
			self.start_datetime = self.sky_map.times [0]
			self.end_datetime = self.sky_map.times[-1]
			self.times = self.sky_map.times
		else:
			self.end_datetime = dt.datetime.strptime(self.end_datetime, "%Y%m%d-%H:%M:%S")
			N_steps = max(self.Nb_a_pc, self.Nb_e_pc, self.sky_map.Nt)
			delta_t = (self.end_datetime - self.start_datetime) / N_steps
			self.times = np.array([self.start_datetime + i * delta_t for i in range(N_steps)])


		self.direct_light_mode = in_dict["direct_mode"].lower()
		if self.direct_light_mode not in ["none", "only", "add"]:
			self.direct_light_mode = "add"

		print(pomerol_configuration.chaos_model_file)
		self.B_model = chaos.load_CHAOS_matfile(pomerol_configuration.chaos_model_file)


		# print(self.times)

		self.is_single_observation = False
		if self.Nb_a_pc == self.Nb_e_pc == 1:
			self.is_single_observation = True

		if mpi_rank==0: print("Has_sky_emission:", self.has_sky_emission)
		#Show the skymaps if they change in time or we move the instrument
		# if self.has_sky_emission and (self.is_time_dependant or not self.is_single_observation):
		if self.has_sky_emission and self.show_sky_cube and mpi_rank==0:
			self.sky_map.MakeSkyCubePlot(self.a_pc_list, self.e_pc_list, self.ouv_pc)



		### Initialize the groundmap. Even if we don't use it, some parameters must be initialized
		self.ground_map = GroundMap(in_dict, self.Nb_a_pc, self.Nb_e_pc)

		self.has_ground_emission = self.ground_map.exist
		if mpi_rank == 0: print("Has_ground_emission:", self.has_ground_emission)

		### Initialize the atlitude maps
		self.alt_map = ElevationMap(in_dict)
		self.use_elevation_map = bool(int(in_dict["use_alt_map"]))
		if self.show_altitude:
			self.alt_map.PlotMap()

		self.instrument_altitude = self.alt_map.GetAltitudeFromLonLat(self.ground_map.A_lon, self.ground_map.A_lat) / 1000. #In km

		self.albedo_mode = in_dict["ground_albedo"].lower()
		if mpi_rank == 0 and self.show_ground_albedo:
			if "altitude" in self.albedo_mode:
				layers = np.array([float(a)*1000 for a in self.albedo_mode.split("_")[2:-1:2]])
				self.MakeGroundMapPlot(iso=layers)
			else:
				self.MakeGroundMapPlot()
			plt.show()

		if self.has_ground_emission and self.has_sky_emission and self.albedo_mode != "none":
			self.albedo_map = np.zeros_like(self.ground_map.I_map)
			if "constant" in self.albedo_mode:
				a = float(self.albedo_mode.split("_")[1])
				self.albedo_map += a

			elif "altitude" in self.albedo_mode:
				layers = np.array([float(a) for a in self.albedo_mode.split("_")[2:-1:2]])
				albedo = np.array([float(a) for a in self.albedo_mode.split("_")[1::2]])
				# print("layers, albedo", layers, albedo)
				for id, d in enumerate(self.ground_map.mid_distances):
					for ia, a in enumerate(self.ground_map.mid_azimuts):
						lon, lat = AzDistToLonLat(a, d, origin_lon = self.ground_map.A_lon, origin_lat = self.ground_map.A_lat)
						alt = self.alt_map.GetAltitudeFromLonLat(lon, lat)

						layer = int(np.searchsorted(layers, alt))
						self.albedo_map[id, ia] = albedo[layer]

			if mpi_rank == 0:
				print("Calculating sky on ground reflection with albedo", self.albedo_mode)
			self.ground_map.ProlongateTime(self.sky_map.Nt)

			reciever_cube = np.empty_like(self.ground_map.cube)

			for t in range(self.sky_map.Nt)[mpi_rank::mpi_size]:
				# t = i * mpi_size + mpi_rank
				self.ground_map.cube[t, :, :] += self.ground_map.I_map
				flat_map = self.sky_map.GetFlatMap(self.ground_map.A_lon, self.ground_map.A_lat, self.ground_map.mid_azimuts, self.ground_map.mid_distances, time = t, Nmax = 1000, tau0 = self.atmosphere.tau_0 + self.atmosphere.GetO3Absorbtion(0, self.sky_map.h))
				# flat_map, lon, lat = self.sky_map.GetFlatMap(self.ground_map.A_lon, self.ground_map.A_lat, self.ground_map.radius, time = t, Nmax = 1000, tau0 = self.atmosphere.tau_0)
				# flat_map, lon, lat = self.sky_map.GetFlatMap(self.ground_map.A_lon, self.ground_map.A_lat, mid_lon = self.ground_map.mid_longitudes, mid_lat = self.ground_map.mid_latitudes, time = t, Nmax = 10000, tau0 = self.atmosphere.tau_0)

				self.ground_map.cube[t, :, :] += flat_map * self.albedo_map

			mpi_comm.Allreduce(self.ground_map.cube, reciever_cube, op = MPI.SUM)
			self.ground_map.cube = reciever_cube


		self.add_B_pola = bool(int(in_dict["add_B_pola"]))


		# if mpi_rank == 0:
		# 	self.ground_map.MakePlots(0, 0)
		# 	self.sky_map.MakeSkyCubePlot(self.a_pc_list, self.e_pc_list, self.ouv_pc)
			# plt.show()


	# def GetGalacticCoord(self, az, el, time):
	# 	return GetGalacticCoord(self.ground_map.A_lon, self.ground_map.A_lat, self.instrument_altitude*1000)
		# location = EarthLocation(lon=self.ground_map.A_lon, lat=self.ground_map.A_lat, height=self.instrument_altitude*1000*u.m)
		# ptcu = SkyCoord(obstime = time, location=location, alt=el*u.rad, az=az*u.rad, frame='altaz')
		# Glon, Glat = ptcu.galactic.l.value, ptcu.galactic.b.value
		#
		# print("GLon, Glat", Glon, Glat)
		#
		# return Glon, Glat #in degrees

	def GetPolaFromB(self, DoLP_max = 100):
		"""Computes the apparent angle of the mag field AoBapp for self.observation and its angle with the line of sight AoBlos.
		Then computes a DoLP between 0 (AoBlos=0) and DoLP_max (AoBlos=90) in [0, 1] and AoBapp in radians."""

		### Chaos 6 B field
		# B_chaos = self.obs.GetBatP(B_model = self.B_model)[1]

		### MANUAL B/current
		B_up 	= 0
		B_east 	= 1
		B_north = 0
		B_chaos = np.array([0, 0, B_up, B_east, B_north])

		los = [self.obs.los_uen[0,0], self.obs.los_uen[1,0], self.obs.los_uen[2,0]]
		B = B_chaos[2:]

		# pola_direction = np.cross(self.obs.los_uen, [[B_up], [B_east], [B_north]])
		pola_direction = np.cross(los, B)

		AoBapp, AoBlos = self.obs.GetApparentAngle(B_chaos[2:], Blos=True)
		pola_app_angle, pola_los_angle = self.obs.GetApparentAngle(pola_direction, Blos=True)

		DoLP = np.sin(AoBlos)**2 / (1 + np.cos(AoBlos)**2)

		# return DoLP_max / 100 , AoBapp
		return DoLP_max * DoLP / 100., pola_app_angle



	def GetStarlight(self, time, ie_pc, ia_pc):
		time = self.times[time]
		el = self.e_pc_list[ie_pc]
		az = self.a_pc_list[ia_pc]

		return GetStarlight(self.ground_map.A_lon, self.ground_map.A_lat, self.instrument_altitude*1000, time, el, az)
		# Glon, Glat = self.GetGalacticCoord(az, el, time)
		#
		# # print(Glon, Glat)
		#
		# b = np.array([0, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80])
		# Jv = np.array([371.5, 246.5, 175.3, 136.5, 102.9, 68.7, 51.5, 41.3, 35.5, 32.5, 31.1]) #in S_10 units
		# Jv = self.S10TonanoWatt(Jv)
		#
		# return np.interp(abs(Glat), b, Jv)


	def SetObservation(self, a_pc, e_pc):
		#Create an observation object
		self.obs = ObservationPoint(self.ground_map.A_lon, self.ground_map.A_lat, self.sky_map.h, a_pc, e_pc, A_alt = self.instrument_altitude, init_full=False) # Usefull observation object

		self.atmosphere.SetObservation(self.obs, self.ouv_pc)



		# range_min, range_max = self.atmosphere.h_r_min / np.sin(e_pc) , self.atmosphere.h_r_max / np.sin(e_pc)
		#
		# #numbers of scattering point along the line of sight
		# los_length = range_max - range_min
		#
		# if self.atmosphere.Nlos > 1:
		# 	#array of range for each scattering points (distance along los from the instrument)
		# 	# dlos_list = np.arange(range_min + self.atmosphere.d_los / 2,  range_max - self.atmosphere.d_los / 2, self.atmosphere.d_los)
		#
		# 	self.range_list = np.array([x**2 for x in np.linspace(np.sqrt(range_min), np.sqrt(range_max), self.atmosphere.Nlos + 1)]) #in km
		# 	self.mid_range_list = np.array([(range_list[i+1] + range_list[i]) / 2. for i in range(0, self.atmosphere.Nlos)]) #in km
		# 	self.d_los_list = np.array([range_list[i+1] - range_list[i] for i in range(0, self.atmosphere.Nlos)]) #in km
		#
		# else:
		# 	dlos_list = np.array([90. / np.sin(e_pc)])
		# 	# dlos_list = [(range_min + range_max) / 2.]

		# print("DEBUG dlos list:", dlos_list)

		#array of lon, lat, alt (above sea level) for each scattering point
		self.dlos_list = [self.obs.GetPCoordinatesFromRange(d) for d in self.atmosphere.mid_range_list]

		#list of absolute altitudes (above sea level) for each scattering point
		self.altitudes = np.array([d[2] for d in self.dlos_list])
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

		# mpi_size = MPI.COMM_WORLD.Get_size()
		# mpi_rank = MPI.COMM_WORLD.Get_rank()
		# mpi_name = MPI.Get_processor_name()

		count = 0
		start_time = tm.time()
		N = self.sky_map.N * self.Nalt / mpi_size
		if mpi_rank==0: print(N, "bins to compute per processor")
		# with tqdm(total=N, file=sys.stdout) as pbar:
		for ie, e in enumerate(self.sky_map.mid_elevations):
			for ia, a in enumerate(self.sky_map.mid_azimuts):
				if self.sky_map.cube[time, ie, ia] > 0:

					if self.use_analytic and mpi_rank == 0:
						V, Vcos, Vsin = self.ComputeSingleRSSkyPointSourceAnal(time, ia, a, ie, e, self.sky_map.h)

						self.sky_map.V_map[time, ie_pc, ia_pc, ie, ia] += V
						self.sky_map.Vcos_map[time, ie_pc, ia_pc, ie, ia] += Vcos
						self.sky_map.Vsin_map[time, ie_pc, ia_pc, ie, ia] += Vsin

					else:
						for ialt, alt in enumerate(self.altitudes[mpi_rank::mpi_size]): #for every altitude between minimum and maximum scattering altitude
							ialt = ialt * mpi_size + mpi_rank
							# print("DEBUG MPI:", mpi_size, mpi_rank, ialt)

							sca_from_E_is_visible = True
							# if self.use_elevation_map:
							# 	sca_from_E_is_visible = self.alt_map.IsVisible(sca_lon, sca_lat, alt1=sca_alt, lon2=lon, lat2=lat, dlos = 0.5)

							if sca_from_E_is_visible:
								V, Vcos, Vsin = self.ComputeSingleRSSkyPointSource(time, ia, a, ie, e, ialt)
								# print("DEBUG Vs:", V, Vcos, Vsin)
								self.sky_map.V_map[time, ie_pc, ia_pc, ie, ia]    += V
								self.sky_map.Vcos_map[time, ie_pc, ia_pc, ie, ia] += Vcos
								self.sky_map.Vsin_map[time, ie_pc, ia_pc, ie, ia] += Vsin
								# self.sky_map.total_scattering_map[time, ie_pc, ia_pc, ie, ia] += I0 # Intensity of light scattered from a given (e, a)
								# self.sky_map.scattering_map[time, ie_pc, ia_pc, ie, ia] += DoLP * I0 # Intensity of polarized light scattered from a given (e, a). Given a point source, the AoRD is the same for every altitude -> the addition makes sens

								# pbar.update(1)

						# 	self.sky_map.DoLP_map = self.sky_map.scattering_map / self.sky_map.total_scattering_map # DoLP of a given (e, a)
						#
						# self.sky_map.AoRD_map[time, ie_pc, ia_pc, ie, ia] = AoLP

	def ComputeGroundMaps(self, time, ia_pc, ie_pc):
		"""Compute the contribution of the ground map at a given time"""
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

		start_time = tm.time()



		# mpi_size = MPI.COMM_WORLD.Get_size()
		# mpi_rank = MPI.COMM_WORLD.Get_rank()
		# mpi_name = MPI.Get_processor_name()

		# file_object = open('/home/bossel/These/Analysis/results/rayleigh/Aerosols/ptsrc_L100_a180_d5_ptcu_a180_e45.txt', mode='w')
		# print("source luminance (nW/m2/sr),source surface (m2),AR (km),Altitude (km),ER (km),Theta (degree),180 - theta,beta(z) (km-1),P_rs (theta),C_aer x SSA (km-1),P_aer (theta),Vol (km3),instrument surface (km2),opt_depth RS along ER,opt depth O3 along ER,opt depth aer along ER,transmittance along ER,transmittance along line of sight,F_RS (nW),F_aer (nW),DoLP RS,DoLP aer,F total (nW),DoLP total", file = file_object)
		# file_object.close()

		N = len(self.ground_map.mid_azimuts) * len(self.ground_map.mid_distances) * self.Nalt / mpi_size
		if mpi_rank==0: print(N, "bins to compute per processor")
		# with tqdm(total=N, file=sys.stdout) as pbar:
		for iaz, az in enumerate(self.ground_map.mid_azimuts):
			for idist, dist in enumerate(self.ground_map.mid_distances):
				if self.ground_map.cube[time, idist, iaz] > 0:
					a_rd, e_rd = az, DistToAngle(dist, lat = self.ground_map.A_lat, az = az)
					E_lon, E_lat = AzDistToLonLat(az, dist,  origin_lon = self.ground_map.A_lon, origin_lat = self.ground_map.A_lat)
					# print(a_rd * RtoD, e_rd * RT)
					# print(self.dlos_list)

					if self.use_analytic and mpi_rank == 0:
						V, Vcos, Vsin = self.ComputeSingleRSGroundPointSourceAnal(time, iaz, idist, a_rd, e_rd, 0)

						self.ground_map.V_map[time, ie_pc, ia_pc, idist, iaz] += V
						self.ground_map.Vcos_map[time, ie_pc, ia_pc, idist, iaz] += Vcos
						self.ground_map.Vsin_map[time, ie_pc, ia_pc, idist, iaz] += Vsin
					else:
						count_visible = 0
						for icount, sca_pos in enumerate(self.dlos_list[mpi_rank::mpi_size]): #for every altitude between minimum and maximum scattering altitude
							icount = icount * mpi_size + mpi_rank
							sca_lon, sca_lat, sca_alt = sca_pos
							# print("DEBUG MPI:", mpi_size, mpi_rank, icount*mpi_size+mpi_rank)
							# print(sca_lon, sca_lat, sca_alt)

							sca_from_E_is_visible = True
							if self.use_elevation_map and count_visible == 0: #If we take the altitude into account and if we are still in a mountains shadow. As soon as the first (lowest) scattering point is visible, we don't check the higher ones.
								sca_from_E_is_visible = self.alt_map.IsVisible(sca_lon, sca_lat, alt1=sca_alt, lon2=E_lon, lat2=E_lat, dlos = None)

							if sca_from_E_is_visible:

								if count_visible == 0:
									# print(f"First visible point nb {icount} at alt {sca_alt}")
									count_visible += 1
									if self.ground_map.mountain_shadow_heights[time, ie_pc, ia_pc, idist, iaz] == 0 or self.ground_map.mountain_shadow_heights[time, ie_pc, ia_pc, idist, iaz] < sca_alt:
										self.ground_map.mountain_shadow_heights[time, ie_pc, ia_pc, idist, iaz] = sca_alt

								V, Vcos, Vsin = self.ComputeSingleRSGroundPointSourceNormal(time, iaz, idist, a_rd, e_rd, icount)
								# if sca_alt == 0:
								# 	print(iaz, idist, a_rd, e_rd, sca_alt)

								self.ground_map.V_map[time, ie_pc, ia_pc, idist, iaz] += V
								self.ground_map.Vcos_map[time, ie_pc, ia_pc, idist, iaz] += Vcos
								self.ground_map.Vsin_map[time, ie_pc, ia_pc, idist, iaz] += Vsin

								# self.ground_map.total_scattering_map[ie_pc, ia_pc, ilat, ilon] += I0 # Intensity of light scattered from a given (e, a)
								# self.ground_map.scattering_map[ie_pc, ia_pc, ilat, ilon] += (w_DoLP * I0) # Intensity of polarized light scattered from a given (e, a). Given a point source, the AoRD is the same for every altitude -> the addition makes sens

								# count += 1
								# self.Progress(count, N, suffix="of ground point sources done")
								# pbar.update(1)


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

	def ComputeSingleRSGroundPointSourceAnal(self, time, iaz, idist, a_E, e_E, alt):
		I0 = self.ground_map.cube[time, idist, iaz] #[nW/m2/sr]

		AR, RE, RD_angle, RAE, alt_E = self.GetGeometryFromAzDist(a_E, e_E, alt, ret_RAE=True)
		AE = AngleToDist(e_E, lat = self.ground_map.A_lat, az = a_E)

		# print("AR, RE, AE", AR, RE, AE)

		V_T, Vcos_T, Vsin_T = 0, 0, 0

		ER = lambda r: np.sqrt(r**2 + AE**2 - 2 * r * AE * np.cos(RAE))
		theta = lambda r: np.arcsin(AE * np.sin(RAE) / ER(r))
		sinREA = lambda r: r * np.sin(RAE) / ER(r)
		V = lambda r, alpha: 2 * np.pi * (1 - np.cos(alpha)) * r ** 2

		beta0 = self.atmosphere.beta_0
		tau0 = self.atmosphere.tau_0

		# beta = lambda beta0, r: beta0 * np.exp(H * r)
		# tau = lambda r: tau0 * np.exp(H * r) * (ER(r) / r + 1) / np.sin(self.e_pc)
		altitude = lambda r: self.instrument_altitude + r * np.sin(self.e_pc)
		pressure = lambda a: self.atmosphere.GetProfileValue(a, "PRE")
		temperature = lambda a: self.atmosphere.GetProfileValue(a, "TEM")
		beta = lambda beta0, r: beta0 * pressure(altitude(r)) * self.atmosphere.T0 / self.atmosphere.P0 / temperature(altitude(r))
		tau = lambda r: self.atmosphere.GetRSOpticalDepth(0, altitude(r)) * ER(r) / altitude(r) + self.atmosphere.GetRSOpticalDepth(self.instrument_altitude, altitude(r)) / np.sin(self.e_pc)
		# (tau0 / self.atmosphere.P0) * ((pressure(altitude(r)) - pressure(0))/sinREA(r) + (pressure(altitude(r)) - pressure(self.instrument_altitude))/np.sin(self.e_pc))

		K0 = (3./8) * self.ground_map.GetArea(idist) * self.PTCU_area * self.ground_map.cube[time, idist, iaz] * (1 - np.cos(self.ouv_pc)) * beta0 * self.atmosphere.T0 / self.atmosphere.P0
		K2 = (AE * np.sin(RAE)) ** 2


		dF = lambda r: K0 * (2 - K2 / ER(r)**2) * np.exp(-tau(r)) * pressure(altitude(r)) / ER(r)**2 / temperature(altitude(r))
		# dF = lambda r: K0 * np.exp(H * r) * (2 - K2 / ER(r)**2) * np.exp(-tau(r)) / ER(r)**2
		DoLP = lambda r: np.sin(theta(r))**2 / (1 + np.cos(theta(r))**2)
		dFpola = lambda r: dF(r) * DoLP(r)
		dFNpola = lambda r: dF(r) * (1 - DoLP(r))

		min_range = max(0, (self.atmosphere.h_r_min - self.instrument_altitude) / np.sin(self.e_pc))
		max_range = (self.atmosphere.h_r_max - self.instrument_altitude) / np.sin(self.e_pc)
		F = integrate.quad(dF, min_range, max_range)#, points=[0])
		Fpola = integrate.quad(dFpola, min_range, max_range)#, points=[0])


		AoLP = self.obs.GetRayleighAngle(a_E, e_E)


		V_T, Vcos_T, Vsin_T = self.GetVParamFromLightParam(F[0], Fpola[0] / F[0], AoLP)
		return V_T, Vcos_T, Vsin_T


	def ComputeSingleRSGroundPointSourceNormal(self, time, iaz, idist, a_E, e_E, ialt):
		"""a_E, e_E is the azimut and distance (in radians) of the source seen from the instrument."""

		I0 = self.ground_map.cube[time, idist, iaz] #[nW/m2/sr]

		alt = self.altitudes[ialt]

		AR, RE, RD_angle, RAE, alt_E = self.GetGeometryFromAzDist(a_E, e_E, alt, ret_RAE=True)
		# print(AR, RE, RD_angle  *  RtoD, RAE * RtoD)
		# AE = self.ground_map.AngleToDist(e_E, lat = self.ground_map.A_lat, az = a_E)


		vol = self.atmosphere.GetVolume(AR, self.ouv_pc, index = ialt, unit="km", type="index") #km3
		# if mpi_rank == 0:
		# print(AR, vol)

		cut_V = self.CutVolume(RE, vol, max_angle = self.max_angle_discretization)

		# print("DEBUG: AR, RE, RD_angle, alt", AR, RE, RD_angle*RtoD, alt)

		if cut_V > 1:
			geo_list = self.GetDiscreteRSPoint(self.a_pc, self.e_pc, alt, AR, RE, a_E, e_E, "ground", N = cut_V)
		else:
			obs = ObservationPoint(self.ground_map.A_lon, self.ground_map.A_lat, self.sky_map.h, self.a_pc, self.e_pc, A_alt = self.instrument_altitude, init_full=False)
			alt = obs.GetPCoordinatesFromRange(AR)[2]
			AoLP = obs.GetRayleighAngle(a_E, e_E)

			geo_list = [[self.a_pc, self.e_pc, AR, AR, RE, RD_angle, alt, AoLP, vol, self.atmosphere.d_los_list[ialt]]]

		V_T, Vcos_T, Vsin_T = 0, 0, 0

		vol_T = 0

		# vol_T = np.sum([g[8] for g in geo_list])
		# if not abs(vol_T - vol) < 10**-10:
		# 	print(vol_T, vol, vol_T / vol, abs(vol_T - vol) < 10**-10)
		# file_object = open('/home/bossel/These/Analysis/results/rayleigh/Aerosols/ptsrc_L100_a180_d5_ptcu_a180_e45.txt', mode='a')
		# print("source luminance (nW/m2/sr),source surface (m²),AR (km),Altitude (km),ER (km),Theta (°),180 – theta,beta(z) (km-1),P_rs (theta),C_aer x SSA (km-1),P_aer (theta),Vol (km³),instrument surface (km²),opt_depth RS along ER,opt depth O3 along ER,opt depth aer along ER,transmittance along ER,transmittance along line of sight,F_RS (nW),F_aer (nW),DoLP RS,DoLP aer,F total (nW),DoLP total", file = file_object)
		# file_object.close()
		# file_object = open('/home/bossel/These/Analysis/results/rayleigh/Aerosols/ptsrc_L100_a180_d5_ptcu_a180_e45.txt', mode='a')

		for a_A, e_A, r, AR, RE, RD_angle, alt, AoLP, dvol, dr in geo_list:

			# print(AR, RD_angle*RtoD)

			if RE == 0:
				print("WARNING: ER = 0 !!!")
				continue
			# print(a_A*RtoD, e_A*RtoD, r, AR, RE, RD_angle*RtoD, alt, AoLP*RtoD, da*RtoD, dr)

			# print("DEBUG I0 1:", I0, I0*self.ground_map.GetArea(ilat) / RE ** 2)
			I0 = self.ground_map.cube[time, idist, iaz] #[nW/m2/sr]
			# print(f"I0 {I0}")
			I0 *= self.ground_map.GetArea(idist) # [nW / m2 / sr * m2] = [nW / sr]
			# print(f"I0*A_E {I0}")


			if RE > 0:
				I0 /= RE ** 2 # [nW / m2 / sr * m2 / km2] = [nW / km2]
			else:
				I0 = 0

			opt_depth = 0
			O3_abs = 0
			aer_abs = 0
			if alt != 0:
				delta_z = abs(alt_E - alt)
				opt_depth = self.atmosphere.GetRSOpticalDepth(alt_E, alt) * RE / delta_z
				if self.atmosphere.use_ozone:
					O3_abs = self.atmosphere.GetO3Absorbtion(alt_E, alt) * RE / delta_z
				if self.atmosphere.use_aerosol:
					aer_abs = self.atmosphere.GetAerosolsAbsorbtion(alt_E, alt) * RE / delta_z
				# print("DEBUG opt_depth: ER", opt_depth, 1-opt_depth, np.exp(-opt_depth))

			# if mpi_rank == 0: print(O3_abs, np.exp(- O3_abs))
			I0 *= np.exp(- opt_depth - O3_abs - aer_abs) # [nW / km2]

			# print("DEBUG I0 2:", I0)
			# print("DEBUG scattered:", self.GetScattered(I0, AR, RE, RD_angle, alt)[0]/I0)

			# Vol = self.atmosphere.GetVolume(AR, da, dr=dr, unit="km") #in km3
			vol_T += dvol
			Crs = self.atmosphere.GetRSVolumeCS(alt) #in km-1
			P = self.atmosphere.GetRSPhaseFunction(self.wavelength, RD_angle)# in sr

			Caer, Paer, DoLP_aer = 0, 0, 0
			if self.atmosphere.use_aerosol:
				Caer = self.atmosphere.GetAerosolCS(alt) #in km-1
				Paer, DoLP_aer = self.atmosphere.GetAerosolPhaseFunction(RD_angle) # Paer in sr. pi=FRONTscattering, 0=backscattring
			# print("DEBUG V, Crs, P", V, Crs, P)

			I0_rs = 0
			I0_aer = 0
			if AR != 0:
				# if mpi_rank == 0:
				# 	print(AR, dvol)
				I0_aer = I0 * Caer * Paer * dvol * self.PTCU_area / AR ** 2 / 4 / np.pi # [nW / km2] * [km-1 * km3 * km2 * km-2] = [nW]
				I0_rs = I0 * Crs * P * dvol * self.PTCU_area / AR ** 2 / 4 / np.pi # [nW / km2] * [km-1 * km3 * km2 * km-2] = [nW]
				# print("alt, V, Crs, P, omega", alt, V, Crs, P, self.PTCU_area / (AR*1000) ** 2)

			f = 1 #(1 + self.atmosphere.depola) / (1 - self.atmosphere.depola)
			DoLP_rs = np.sin(RD_angle)**2 / (f + np.cos(RD_angle)**2) # DoLP dependance on scattering angle for Rayleigh Scattering


			I0_aer *= self.atmosphere.los_transmittance[ialt]
			I0_rs *= self.atmosphere.los_transmittance[ialt]
			# I0 *= self.atmosphere.los_transmittance[ialt]


			# Option 1: Mix of flux and DoLP when using aerosols
			I0 = I0_rs + I0_aer
			if I0 != 0:
				DoLP = (I0_rs * DoLP_rs + I0_aer * DoLP_aer) / I0
			else:
				DoLP = 0

			if DoLP < 0: # For Mie scattering, a DoLP < 0 = AoLP is parallel to scatering plane ! So 90 degrees from Rayleigh scattering.
				AoLP += np.pi/2
				DoLP *= -1


			# Option 2: Mix of flux and DoLP when using aerosols
			# I0 = (Crs * I0_rs + Caer * I0_aer) / (Caer + Crs)
			# DoLP = (Crs * I0_rs * DoLP_rs + Caer * I0_aer * DoLP_aer) / (Caer + Crs) / I0


			# with open('/home/bossel/These/Analysis/results/rayleigh/Aerosols/ptsrc_L100_a180_d5_ptcu_a180_e45.txt', mode='w') as file_object:
			# print(self.ground_map.cube[time, idist, iaz], self.ground_map.GetArea(idist), AR, alt, RE, RD_angle*RtoD, 180-RD_angle*RtoD, Crs, P, Caer, Paer, dvol, self.PTCU_area, opt_depth, O3_abs, aer_abs, np.exp(- opt_depth - O3_abs - aer_abs), self.atmosphere.los_transmittance[ialt], I0_rs, I0_aer, DoLP_rs, DoLP_aer, I0, DoLP, sep=",", file=file_object)

			# V, Vcos, Vsin = self.GetVParamFromLightParam(I0, DoLP, AoLP) #Not using aerosols
			V, Vcos, Vsin = self.GetVParamFromLightParam(I0, DoLP, AoLP) #Using aerosols

			V_T += V
			Vcos_T += Vcos
			Vsin_T += Vsin

		# print(V_T, Vcos_T, Vsin_T)

		# if not abs(vol_T - vol) < 10**-10:
		# 	print(vol_T, vol, vol_T / vol, abs(vol_T - vol) < 10**-10)

		# V_T *= vol / vol_T
		# Vcos_T *= vol / vol_T
		# Vsin_T *= vol / vol_T

		# file_object.close()
		return V_T, Vcos_T, Vsin_T


	def ComputeSingleRSSkyPointSourceAnal(self, time, ia_E, a_E, ie_E, e_E, alt):
		I0 = self.sky_map.cube[time, ie_E, ia_E] # [nW / m2/ sr]

		AR, RE, RD_angle, RAE  = self.GetGeometryFromAzEl(a_E, e_E, alt, ret_RAE=True)
		AE = GenPythagore(AR, RE, RD_angle)

		# vol = self.atmosphere.GetVolume(AR, self.ouv_pc, unit="km", type="pole")

		# cut_V = self.CutVolume(RE, vol, max_angle = self.max_angle_discretization)

		# print("DEBUG: AR, RE, RD_angle, alt", AR, RE, RD_angle*RtoD, alt)

		# geo_list = self.GetDiscreteRSPoint(self.a_pc, self.e_pc, alt, AR, RE, a_E, e_E, "ground", N=cut_V)

		# V_T, Vcos_T, Vsin_T = 0, 0, 0
		#
		# vol_T = 0

		ER = lambda r: np.sqrt(r**2 + AE**2 - 2 * r * AE * np.cos(RAE))
		theta = lambda r: np.arcsin(AE * np.sin(RAE) / ER(r))
		V = lambda r, alpha: (2./3) * np.pi * (1 - np.cos(alpha)) * r ** 2

		wl  = self.atmosphere.wavelength / 1000.
		A = 7.68246 * 10 ** (-4) #for km-1
		B = 3.55212
		C = 1.35579
		D = 0.11563
		beta0 = A * wl ** (- (B + C * wl + D / wl))
		# print("beta0", beta0)

		A = 6.49997 * 10 ** (-3)
		B = 3.55212
		C = 1.35579
		D = 0.11563
		tau0 = A * wl ** (- (B + C * wl + D / wl))

		g = 9.80665
		M = 0.02896968
		T0 = 288.16
		R0 = 8.314462618
		H = - g * M / T0 / R0 * np.sin(self.e_pc)
		# print("H", H)

		beta = lambda beta0, r: beta0 * np.exp(H * r)
		tau = lambda r: tau0 * np.exp(H * r) * (ER(r) / r + 1) / np.sin(self.e_pc)

		K0 = beta0 * self.sky_map.GetPixelArea(ie_E) * self.PTCU_area * self.sky_map.cube[time, ie_E, ia_E] * (1 - np.cos(self.ouv_pc)) / 8
		K2 = (AE * np.sin(RAE)) ** 2


		dF = lambda r: K0 * np.exp(H * r) * (2 - K2 / ER(r)**2) * np.exp(-tau(r)) / ER(r)**2
		DoLP = lambda r: np.sin(theta(r))**2 / (1 + np.cos(theta(r))**2)
		dFpola = lambda r: dF(r) * DoLP(r)
		dFNpola = lambda r: dF(r) * (1 - DoLP(r))

		min_range = self.atmosphere.h_r_min / np.sin(self.e_pc)
		max_range = self.atmosphere.h_r_max / np.sin(self.e_pc)
		F = integrate.quad(dF, min_range, max_range, points=[0])
		Fpola = integrate.quad(dFpola, min_range, max_range, points=[0])

		obs = ObservationPoint(self.ground_map.A_lon, self.ground_map.A_lat, self.sky_map.h, self.a_pc, self.e_pc, A_alt = self.instrument_altitude, init_full=False)

		AoLP = obs.GetRayleighAngle(a_E, e_E)


		V_T, Vcos_T, Vsin_T = self.GetVParamFromLightParam(F[0], Fpola[0] / F[0], AoLP)
		return V_T, Vcos_T, Vsin_T


	def ComputeSingleRSSkyPointSource(self, time, ia_E, a_E, ie_E, e_E, ialt):
		"""t: time of the map
		ia_E, a_E, ie_E, e_E: az and el and their indexes of the emission point E."""

		I0 = self.sky_map.cube[time, ie_E, ia_E] # [nW / m2/ sr]

		# print("DEBUG I0", I0)

		alt = self.altitudes[ialt]

		AR, RE, RD_angle = self.GetGeometryFromAzEl(a_E, e_E, alt)
		vol = self.atmosphere.GetVolume(AR, self.ouv_pc, index = ialt, unit="km", type="index")
		# vol = self.atmosphere.volumes[ialt]

		cut_V = self.CutVolume(RE, vol, max_angle = self.max_angle_discretization)

		#Get a list of all small scattering volumes in the big troncated cone. 1 point if distance from emission to cone is 10 km
		if cut_V > 1:
			geo_list = self.GetDiscreteRSPoint(self.a_pc, self.e_pc, alt, AR, RE, a_E, e_E, "sky", N = cut_V)
		else:
			obs = ObservationPoint(self.ground_map.A_lon, self.ground_map.A_lat, self.sky_map.h, self.a_pc, self.e_pc, A_alt = self.instrument_altitude, init_full = False)
			alt = obs.GetPCoordinatesFromRange(AR)[2]
			AoLP = obs.GetRayleighAngle(a_E, e_E)

			geo_list = [[self.a_pc, self.e_pc, AR, AR, RE, RD_angle, alt, AoLP, vol, self.atmosphere.d_los_list[ialt]]]


		V_T, Vcos_T, Vsin_T = 0, 0, 0

		vol_T = 0

		for a_A, e_A, r, AR, RE, RD_angle, alt, AoLP, dvol, dr in geo_list:

			# print(a_A*RtoD, e_A*RtoD, r, AR, RE, RD_angle*RtoD, alt, AoLP*RtoD, da*RtoD, dr)

			I0 = self.sky_map.cube[time, ie_E, ia_E] # [nW / m2/ sr]

			if RE > 0:
				I0 *= self.sky_map.GetPixelArea(ie_E) / RE ** 2 # [nW / m2]
			else:
				I0 = 0

			opt_depth = 0
			O3_abs = 0
			aer_abs = 0
			if alt != 0:
				ER_eff = RE
				delta_z = abs(self.sky_map.h - alt)
				if self.atmosphere.h_r_max < self.sky_map.h:
					ER_eff = GetERForVeryBigE_alt(alt, self.sky_map.h, RE, self.atmosphere.h_r_max)
					delta_z = abs(self.atmosphere.h_r_max - alt)
				opt_depth = self.atmosphere.GetRSOpticalDepth(self.sky_map.h, alt)# * ER_eff / delta_z
				if self.atmosphere.use_ozone:
					O3_abs = self.atmosphere.GetO3Absorbtion(self.sky_map.h, alt)# * ER_eff / delta_z
				if self.atmosphere.use_aerosol:
					aer_abs = self.atmosphere.GetAerosolsAbsorbtion(self.sky_map.h, alt)# * ER_eff / delta_z

				# if mpi_rank	==0 : print("alt, opt_depth", alt, np.exp(-opt_depth), np.exp(-O3_abs), np.exp(-aer_abs))

			I0 *= np.exp(- (opt_depth - O3_abs - aer_abs) * ER_eff / delta_z) # [nW / m2]

			# Vol = self.atmosphere.GetVolume(AR, da, dr=dr, unit="km") #in km3
			vol_T += dvol
			Crs = self.atmosphere.GetRSVolumeCS(alt) #in km-1
			P = self.atmosphere.GetRSPhaseFunction(self.wavelength, RD_angle)

			# print("DEBUG V, Crs, P", V, Crs, P)
			Caer, Paer, DoLP_aer = 0, 0, 0
			if self.atmosphere.use_aerosol:
				Caer = self.atmosphere.GetAerosolCS(alt) #in km-1
				Paer, DoLP_aer = self.atmosphere.GetAerosolPhaseFunction(RD_angle) # in sr. pi=FRONTscattering, 0=backscattring
			# print("DEBUG V, Crs, P", V, Crs, P)

			I0_rs = 0
			I0_aer = 0
			if AR != 0:
				# if self.atmosphere.use_aerosol:
				I0_aer = I0 * Caer * Paer * dvol * self.PTCU_area / AR ** 2 / 4 / np.pi # [nW / km2] * [km-1 * km3 * km2 * km-2] = [nW]
				I0_rs = I0 * Crs * P * dvol * self.PTCU_area / AR ** 2 / 4 / np.pi # [nW / m2] * [km-1 * km3 * m2 * km-2] = [nW]
				# print("alt, V, Crs, P, omega", alt, V, Crs, P, self.PTCU_area / (AR*1000) ** 2)

			f = 1 #(1 + self.atmosphere.depola) / (1 - self.atmosphere.depola)
			DoLP_rs = np.sin(RD_angle)**2 / (f + np.cos(RD_angle)**2) # DoLP dependance on scattering angle
			# DoLP_aer = 0

			I0_aer *= self.atmosphere.los_transmittance[ialt]
			I0_rs *= self.atmosphere.los_transmittance[ialt]
			# I0 *= self.atmosphere.los_transmittance[ialt]

			# Option 1: Mix of flux and DoLP when using aerosols
			I0 = I0_rs + I0_aer
			if I0 != 0:
				DoLP = (I0_rs * DoLP_rs + I0_aer * DoLP_aer) / I0
			else:
				DoLP = 0

			if DoLP < 0: # For Mie scattering, a DoLP < 0 = AoLP is parallel to scatering plane ! So 90 degrees from Rayleigh scattering.
				AoLP += np.pi/2
				DoLP *= -1
			# Option 2: Mix of flux and DoLP when using aerosols
			# I0 = (Crs * I0_rs + Caer * I0_aer) / (Caer + Crs)
			# DoLP = (Crs * I0_rs * DoLP_rs + Caer * I0_aer * DoLP_aer) / (Caer + Crs) / I0


			# if alt != 0:
			# 	opt_depth = self.atmosphere.GetRSOpticalDepth(self.instrument_altitude, alt) * AR / alt
			# 	O3_abs = self.atmosphere.GetO3Absorbtion(self.instrument_altitude, alt) * AR / alt
			# else:
			# 	opt_depth = 0
			# 	O3_abs = 0
			# # print("opt_depth, p.exp(-opt_depth)", opt_depth, np.exp(-opt_depth))
			# I0 *= np.exp(- opt_depth - O3_abs)

			V, Vcos, Vsin = self.GetVParamFromLightParam(I0, DoLP, AoLP)

			V_T += V
			Vcos_T += Vcos
			Vsin_T += Vsin

		# print(V, Vcos, Vsin, I0, DoLP, AoLP)
		return V_T, Vcos_T, Vsin_T

	# @timer
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
			AR_index = np.searchsorted(self.atmosphere.range_list, AR, side = "right")
			r_min, r_max = self.atmosphere.range_list[AR_index-1], self.atmosphere.range_list[AR_index]
			r_list, dr = np.linspace(r_min, r_max, n, retstep=True)
			# if mpi_rank == 0: print(r_list)
			ap_list_mid = np.array([ia + dap / 2. for ia in ap_list[:-1]])
			ep_list_mid = np.array([ie + dep / 2. for ie in ep_list[:-1]])
			r_list_mid = np.array([ir + dr / 2. for ir in r_list[:-1]])
			# if mpi_rank == 0: print(r_list_mid)



			### Returned list of interesting numbers for each subvolumes
			geo_list = np.zeros((len(ap_list_mid) * len(ep_list_mid) * len(r_list_mid), 10))

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

						obs = ObservationPoint(self.ground_map.A_lon, self.ground_map.A_lat, self.sky_map.h, az, el, A_alt = self.instrument_altitude, init_full=False)

						alt = obs.GetPCoordinatesFromRange(r)[2]
						AoLP = obs.GetRayleighAngle(a_E, e_E)

						if map_type == "ground":
							AR, RE, RD_angle, alt_E = self.GetGeometryFromAzDist(a_E, e_E, alt, v_pc_u = Getuen(az, el), obs=obs)
						else:
							AR, RE, RD_angle = self.GetGeometryFromAzEl(a_E, e_E, alt, v_pc_u = Getuen(az, el), obs=obs)

						geo_list[i] = az, el, r, AR, RE, RD_angle, alt, AoLP, dvol, dr
						i += 1
		return geo_list


	# @timer
	def GetVParamFromLightParam(self, I0, DoLP, AoLP):
		"""Returns  the V parameters for a given radiant flux (any unit), a DoLP (between 0 and 1) and an AoLP in radians"""
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
		# N_bins = 180
		# #bins is the bins limits, so 181 values
		# bins, width = np.linspace(-np.pi/2, np.pi/2, N_bins + 1, endpoint=True, retstep=True)
		# bins, width = np.linspace(-np.pi/2 - width/2, np.pi/2 + width/2, N_bins + 1, endpoint=True, retstep=True)
		# mid_bins = [(bins[i+1] + bins[i])/2. for i in range(len(bins)-1)]

		if self.has_ground_emission and not self.has_sky_emission:
			# bins, width = self.ground_map.azimuts, self.ground_map.daz
			N_bins = int(self.ground_map.Naz / 2) + 1
		elif self.has_sky_emission and not self.has_ground_emission:
			# bins, width = self.sky_map.azimuts, self.sky_map.da
			N_bins = int(self.sky_map.Na / 2) + 1
		elif self.sky_map.da >= self.ground_map.daz:
			# bins, width = self.sky_map.azimuts, self.sky_map.da
			N_bins = int(self.sky_map.Na / 2) + 1
		elif self.ground_map.daz >= self.sky_map.da:
			# bins, width = self.ground_map.azimuts, self.ground_map.daz
			N_bins = int(self.ground_map.Naz / 2) + 1
		else:
			print("Error: No ground nor sky emmision. Should terminate before...")


		bins, width = np.linspace(-np.pi/2, np.pi/2, N_bins + 1, endpoint=True, retstep=True)
		bins, width = np.linspace(-np.pi/2 - width/2, np.pi/2 + width/2, N_bins + 1, endpoint=True, retstep=True)


		# print("DEBUG", N_bins, bins)

		#One value for each bin = 180 values
		hst = np.zeros(N_bins)
		# print("DEBUG ROT", self.N_bins, len(self.bins), len(self.mid_bins), len(self.hst))

		###Density=False: hst bins units are similar to intensities. Sum of hst does not depends on N_bins
		###Density=True: hst bins are NOT intensities! But a probability density function (sum height*width == 1)
		hst, bins = np.histogram(A, bins = bins, weights = Ipola, density = True)

		return hst, bins

	def DrawAoLPHist(self, hst, bins, I0, DoLP, AoRD, double=True, save="", show=True, metadata=None):
		"""Make an pyplot histogram of all AoLP contributionsThe histogram is calculated in GetLightParameters()
		Need a call to plt.show() after calling this function."""
		f3, ax = plt.subplots(1, figsize=(16, 8))
		# font = {'size'   : 24}
		# matplotlib.rc('font', **font)

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

		# print("DEBUG ", hst, bins)
		N_bins = len(hst)
		width = bins[1] - bins[0]
		bins = [(bins[i+1] + bins[i])/2. for i in range(len(bins)-1)]

		if not double:
			# bins = bins[:N_bins]
			h = hst
		else:
			bins = np.append(bins[:N_bins], 180 * DtoR + np.array(bins))
			# bins = np.append(bins[:N_bins], 180 * DtoR + np.array(bins)[:N_bins])
			h = np.append(hst, hst)

		# print(h)
		bars = ax.bar(bins, h, width=width, bottom = 0)

		ax.set_title("Weighted AoRD: I0 = " + str(np.format_float_scientific(I0, precision=3)) + " DoLP = " + str(np.round(DoLP, 2)) + " AoRD = " + str(np.round(AoRD*RtoD, 1)))

		if not double:
			ax.plot([AoRD, AoRD], [0, max(hst)], "r")
		else:
			ax.plot([AoRD, AoRD+np.pi], [max(hst), max(hst)], "r")

		# if show:
		# 	plt.show()

		if save:
			plt.savefig(save + '_hist.png', bbox_inches='tight', metadata=metadata)


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




	def GetGeometryFromAzEl(self, a_E, e_E, alt, v_pc_u = None, obs = None, ret_RAE=False):
		"""From known geometry parameters : a_E, e_E of emission point E and alt of scattering, return missing parameters: Distance between emission and scattering and angle of scattering.
		A instrument positions pointing to (a_pc, e_pc)
		E emission points (a_E, e_E, h)
		R rayleigh diffusion point (a_pc, e_pc, atmosphere.h_r)"""

		if obs is None:
			obs = self.obs

		# if e_E == 0:
		# 	e_E = 10**(-30) # Avoid weird behaviour and division by 0

		v_rd_u = Getuen(a_E, e_E) # Vector from observer to emission point
		_, AE = obs.GetAH(elevation = e_E, azimut = a_E, altitude = self.sky_map.h)  # Distance between observer and emisson point

		# if mpi_rank==0: print("AE", AE)

		if v_pc_u is None:
			v_pc_u = self.v_pc_u

		RAE = GetAngle(v_pc_u, v_rd_u) # Angle between line of sight and emission

		# if mpi_rank==0: print("RAE", RAE*RtoD)

		_, AR = obs.GetAH(altitude = alt) # Distance between observer and scattering point

		ER_vect = [AR * pc - AE * rd for pc, rd in zip(v_pc_u, v_rd_u)]
		ER_vect /= np.sqrt(sum([x**2 for x in ER_vect]))
		AER = GetAngle(-np.array(v_rd_u), ER_vect) # Angle between line of sight and emission

		# if mpi_rank==0: print('AER', AER*RtoD)

		RE = GenPythagore(AR, AE, RAE) # Distance between emission and scaterring point
		ARE = np.arcsin(AE * np.sin(RAE) / RE) # Scattering angle


		if RAE + AER < np.pi/2: #Check that the sum of the angles is 180 degrees! Mandatory if phase function is not symmetrical for back- and front- scattering. ARE = 0 -> BACK scattering. ARE=pi -> FRONT-scattring
			ARE = np.pi - ARE
		ARE = np.pi - ARE

		# if mpi_rank==0: print("RE, ARE", RE, ARE*RtoD)

		if ret_RAE:
			return AR, RE, ARE, RAE
		else:
			return AR, RE, ARE #, e_E
	# @timer
	def GetGeometryFromAzDist(self, a_rd, d, h_r, v_pc_u=None, obs = None, ret_RAE=False):
		"""From known geometry parameters : a, d, alt of emission and alt of scatering , return missing parameters: Distance between emission and scattering and angle of scattering.
		A instrument positions pointing to (a_pc, e_pc)
		E emission points (a_rd, e_rd, h)
		R rayleigh diffusion point (a_pc, e_pc, h_r)"""

		if obs is None:
			obs = self.obs

		d = max(d, 10**(-30)) # Avoid weird behaviour and division by 0

		alt_A = self.instrument_altitude
		lon_E, lat_E = AzDistToLonLat(a_rd, d)
		alt_E = self.alt_map.GetAltitudeFromLonLat(lon_E, lat_E)
		delta_alt = alt_E - alt_A

		AE_horiz = AngleToDist(d, lat = self.ground_map.A_lat, az = a_rd) # Distance between observer and emission point. Both points are at see level
		AE = np.sqrt(AE_horiz**2 + delta_alt**2) #REal AE distance when AE are not at the same altitude

		E_el = np.arctan2(delta_alt, AE_horiz) #From instrument in A, elevation at which we see the emission in E. (<0 if alt_E < alt_A)
		v_rd_u = Getuen(a_rd, E_el) # Vector from observer to emission point
		# AE = d * RT # Distance between observer and emisison point. (exact, following the curvature (~flat earth))

		if v_pc_u is None:
			v_pc_u = self.v_pc_u

		RAE = GetAngle(v_pc_u, v_rd_u) # Angle between line of sight and emission
		_, AR = obs.GetAH(altitude = h_r) # Distance between observer and scattering point

		ER_vect = [AR * pc - AE * rd for pc, rd in zip(v_pc_u, v_rd_u)]
		ER_vect /= np.sqrt(sum([x**2 for x in ER_vect]))
		AER = GetAngle(-np.array(v_rd_u), ER_vect) # Angle between line of sight and emission


		RE = GenPythagore(AR, AE, RAE) # Distance between emission and scaterring point
		ARE = np.arcsin(AE * np.sin(RAE) / RE) # Scattering angle.

		if RAE + AER < np.pi/2: #Check that the sum of the angles is 180 degrees! Mandatory if phase function is not symmetrical for back- and front- scattering. ARE = 0 -> BACK scattering. ARE=pi -> FRONT-scattring
			ARE = np.pi - ARE
		ARE = np.pi - ARE

		# AER = np.pi - ARE - self.e_pc #Angle of ER segment with horizontal plane (E is on the ground here so AE is horizontal)

		if ret_RAE:
			return AR, RE, ARE, RAE, alt_E
		else:
			return AR, RE, ARE, alt_E

	# def GetScattered(self, I0, AR, RD_angle, alt):
	# 	"""Given an initial intensity of a source and some geometrical parameter, returns the intensity mesured at the instrument and its DoLP.
	# 	Input parameters: elevation, altitude of scattering, scattering angle, distance between emission and scattering."""
	#
	# 	V = self.atmosphere.GetVolume(AR, self.ouv_pc, unit="km") #in km3
	# 	Crs = self.atmosphere.GetRSVolumeCS(self.wavelength, alt) #in km-1
	# 	P = self.atmosphere.GetRSPhaseFunction(self.wavelength, RD_angle)
	#
	# 	# print("DEBUG V, Crs, P", V, Crs, P)
	#
	# 	if AR != 0:
	# 		I0 *= Crs * P * V * self.PTCU_area / AR ** 2 / 4 / np.pi #Unit
	#
	# 		# print("alt, V, Crs, P, omega", alt, V, Crs, P, self.PTCU_area / (AR*1000) ** 2)
	# 	else:
	# 		I0 = 0
	# 		# print("WARNING!!!! I0==0")
	#
	# 	DoLP = np.sin(RD_angle)**2 / (1 + np.cos(RD_angle)**2) # DoLP dependance on scattering angle
	#
	# 	# print("DEBUG DOLP:", w_DoLP * 100)
	# 	return I0, DoLP


	def GetDirect(self, t, ia, ie):
		"""Returns direct intensity (flux in nW) coming in the instrument"""
		a = self.a_pc_list[ia]
		e = self.e_pc_list[ie]

		if self.has_sky_emission and self.direct_light_mode != "none":
			opt_depth  = self.atmosphere.GetRSOpticalDepth(self.sky_map.h, 0)
			opt_depth += self.atmosphere.GetO3Absorbtion(0, self.sky_map.h)
			opt_depth /= np.sin(e)
			direct = self.sky_map.GetFlux(a, e, self.ouv_pc, t = t) * self.PTCU_area * 1e6 * self.inst_solid_angle
			print("opt_depth, direct:")
			print(opt_depth, direct)
			return direct * np.exp(-opt_depth)
		else:
			return 0



	def MakeGroundMapPlot(self, iso = None):
		f = plt.figure(figsize=(16, 8))
		ax = plt.subplot(111, projection='polar')

		i1 = ax.pcolormesh(self.ground_map.azimuts, self.ground_map.distances, (self.ground_map.cube[0]))
		# cbar1 = f.colorbar(i1, extend='both', spacing='proportional', shrink=0.9, ax=ax)
		# cbar1.set_label('Initial emission map')

		print("DEBUG: CONVERT ALT MAP POLAR")
		polar_alt_map = self.ground_map.CarthToPolar(self.alt_map.map, self.alt_map.mid_longitudes, self.alt_map.mid_latitudes)

		adapt_map = np.append(polar_alt_map, np.vstack(polar_alt_map[:, 0]), axis=1)
		adapt_dist = self.ground_map.mid_distances #np.append(self.ground_map.mid_distances, self.ground_map.mid_distances[0])
		adapt_az = np.append(self.ground_map.mid_azimuts, self.ground_map.mid_azimuts[0]+2*np.pi)

		# contourset = ax.contourf(self.ground_map.mid_azimuts, self.ground_map.mid_distances, polar_alt_map, iso, cmap='Wistia') #iso is a list of the altitudes isobar to plot (in km)
		contourset = ax.contour(adapt_az, adapt_dist, adapt_map, iso, cmap='Wistia') #iso is a list of the altitudes isobar to plot (in m)

		# ax.add_artist(Circle((0, 0), 3, color="red", zorder = 1000))
		ax.plot([0, np.pi], [3, 3], "r")
		ax.plot([-np.pi/2, np.pi/2], [3, 3], "r")

		# ax.tick_params(axis='y', colors='green')
		ax.get_yaxis().set_visible(False)

		ax.set_theta_zero_location("N")
		ax.set_theta_direction(-1)
		ax.set_ylim(0, None)
