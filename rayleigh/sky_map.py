 #!/usr/bin/python3
# -*-coding:utf-8 -*

from mpi4py import MPI
mpi_comm = MPI.COMM_WORLD
mpi_rank = mpi_comm.Get_rank()
mpi_size = mpi_comm.Get_size()
mpi_name = mpi_comm.Get_name()

import sys as sys
import os as os
import numpy as np
import time as tm
import datetime as dt
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Circle
from matplotlib.patches import Ellipse
from matplotlib.patches import Arrow
# from matplotlib.lines import mlines

# import osgeo.gdal as gdal
# gdal.UseExceptions()  # not required, but a good idea

import imageio

from observation import *
from rayleigh_utils import *
from AllSkyData import *

class SkyMap:
	def __init__(self, in_dict, Nb_a_pc, Nb_e_pc):

		self.in_dict = in_dict

		self.path = in_dict["sky_path"]
		self.file = in_dict["sky_file"]

		self.h = float(in_dict["emission_altitude"])

		self.exist = True
		self.mode = in_dict["sky_mode"].lower()

		self.is_point_src = False

		if "band_e" in self.mode:
			self.band_elevation = float(self.mode[-2:]) * DtoR
			# print(self.mode[-2:], self.band_elevation)
			self.mode = self.mode[:-2]
		elif self.mode == "uniform_0" or self.mode == "none":
			self.mode = "none"
			self.exist = False

		# Used if mode == image or movie
		self.image_file = in_dict["sky_path"] + in_dict["sky_file"]
		self.image_wavelength = float(in_dict["sky_wavelength"]) * 10**-9
		self.RtonW = 6.62607015*10**(-34) * 299792458 / self.image_wavelength * 10**9

		self.cube_is_done = False
		self.Nt = 1
		self.times = [dt.datetime(year=1, month=1, day=1, second = 1)]
		if self.exist:
			self.LoadSkyEmmisionsCube(Nb_a_pc, Nb_e_pc)
			self.is_time_dependant = self.Nt > 1
			# print("DEBUG SKYMAP", self.Na, self.Ne, len(self.azimuts), len(self.elevations), len(self.mid_azimuts), len(self.mid_elevations))
			# # print("DEBUG SKYMAP", self.azimuts*RtoD, self.elevations*RtoD)
			# print("DEBUG SKYMAP: elevations", self.elevations*RtoD)
			# print("DEBUG SKYMAP: azimuts", self.azimuts*RtoD)
			# print("DEBUG SKYMAP: mid_elevations", self.mid_elevations*RtoD)
			# print("DEBUG SKYMAP: mid_azimuts", self.mid_azimuts*RtoD)


	def LoadSkyEmmisionsCube(self, Nb_a_pc, Nb_e_pc):
		"""In the case where h!=0, load a cube of data corresponding to the state of the sky at different times, depending on the type of sky we want.
		Unit in [nW / m2/ sr]"""

		if mpi_rank == 0: print("Loading the sky emission cube...")

		#Zone of emission
		self.I_zone_a_min = 0 	* DtoR	#float(in_dict["I_zone_a_min"])	* DtoR
		self.I_zone_a_max = 360 * DtoR	#float(in_dict["I_zone_a_max"]) 	* DtoR
		self.I_zone_e_min = 0 	* DtoR	#float(in_dict["I_zone_e_min"]) 	* DtoR
		self.I_zone_e_max = 90 	* DtoR	#float(in_dict["I_zone_e_max"]) 	* DtoR

		### Initialize some parameters for the case h!=0
		self.N = float(self.in_dict["Nb_emission_points"]) # Number of bins to compute (if h!=0)
		if self.N > 1:
			self.Na, self.Ne = int(2 * np.sqrt(self.N)), int(np.sqrt(self.N) / 2) # Numbers of bins in azimut/elevation. Na = 4*Ne and Na*Ne = N
		else:
			self.Na, self.Ne = 1, 1

		self.da = (self.I_zone_a_max - self.I_zone_a_min) / self.Na # Length of an azimut bin
		self.de = (self.I_zone_e_max - self.I_zone_e_min) / self.Ne # Length of an elevation bin

		self.mid_azimuts = np.linspace(self.I_zone_a_min, self.I_zone_a_max, self.Na)
		self.mid_elevations = np.linspace(self.I_zone_e_min, self.I_zone_e_max, self.Ne)

		self.azimuts = np.append(self.mid_azimuts - self.da/2., self.mid_azimuts[-1] + self.da/2.)
		self.elevations = np.append(self.mid_elevations - self.de/2., self.mid_elevations[-1] + self.de/2.)

		self.Nt = 1
		self.times = [dt.datetime(year=1, month=1, day=1, second = 1)]
		self.maps_shape = (self.Ne, self.Na)
		self.cube = np.zeros((self.Nt, self.Ne, self.Na))

		if "band" in self.mode:
			nb_sub_bands = int(self.mode[:1])

		if "moving_ew_band" in self.mode:
			self.Nt = int(self.in_dict["Nb_emission_maps"])
			self.times = [dt.datetime(year=1, month=1, day=1, second = s+1) for s in range(self.Nt)]

			self.cube = np.zeros((self.Nt, self.Ne, self.Na))
			for ie, e in enumerate(np.linspace(0, np.pi, self.Nt + 2)[1:-1]): #elevation of the auroral band center
				self.cube[ie, :, :] = self.GetBandSkyMap(a_band = 0, e_band = e, length = 1000, thickness = 50, height = 40, band_I = 100, nb_sub_bands = nb_sub_bands)

		elif "band_e" in self.mode:
			self.cube[0, :, :] = self.GetBandSkyMap(a_band = 0, e_band = self.band_elevation, length = 1000, thickness = 10, height = 50, band_I = 100, nb_sub_bands = nb_sub_bands)

		elif "uniform" in self.mode:
			self.cube[:, :, :] += float(self.mode.split("_")[1])

		elif "spot" in self.mode:

			m = self.mode.split("_")
			I = float(m[1][1:])
			a = float(m[2][1:]) * DtoR
			e = float(m[3][1:]) * DtoR

			self.Na, self.Ne = 1, 1
			self.da, self.de = DtoR, DtoR
			self.mid_azimuts = np.array([a])
			self.mid_elevations = np.array([e])
			self.azimuts = np.append(self.mid_azimuts - self.da/2., self.mid_azimuts[-1] + self.da/2.) #Not sure this is used anymore
			self.elevations = np.append(self.mid_elevations - self.de/2., self.mid_elevations[-1] + self.de/2.) #Not sure this is used anymore
			self.maps_shape = (self.Ne, self.Na)
			self.cube = np.zeros((self.Nt, self.Ne, self.Na))

			# self.cube[0, :, :] = self.GetBandSkyMap(a_band = a, e_band = e, length = 50, thickness = 50, height = 50, band_I = 100, nb_sub_bands = 1)
			self.cube[0, 0, 0] = I
			print("self.cube", self.cube)

			self.is_point_src = True

		elif self.mode == "image":
			self.cube[0, :, :] = self.LoadAllSkyImage()

		elif self.mode == "movie" or self.mode == "film":
			self.LoadAllSkyMovie()

		self.data_shape = (self.Nt, Nb_e_pc, Nb_a_pc, len(self.elevations), len(self.azimuts))

		self.scattering_map 		= np.zeros(self.data_shape) # Polarized intensity from (e, a) reaching us
		self.DoLP_map				= np.zeros(self.data_shape) # DoLP of scattered light from (e,a)
		self.total_scattering_map 	= np.zeros(self.data_shape) # Intensity from (e, a) reaching us
		self.AoRD_map 				= np.zeros(self.data_shape) # Angle of polaisation of light from (e,a)
		self.V_map 					= np.zeros(self.data_shape)
		self.Vcos_map 				= np.zeros(self.data_shape)
		self.Vsin_map 				= np.zeros(self.data_shape)

		self.V_total	= np.zeros((self.Nt, Nb_e_pc, Nb_a_pc))
		self.Vcos_total	= np.zeros((self.Nt, Nb_e_pc, Nb_a_pc))
		self.Vsin_total	= np.zeros((self.Nt, Nb_e_pc, Nb_a_pc))

		self.I0_total	= np.zeros((self.Nt, Nb_e_pc, Nb_a_pc))
		self.DoLP_total	= np.zeros((self.Nt, Nb_e_pc, Nb_a_pc))
		self.AoLP_total	= np.zeros((self.Nt, Nb_e_pc, Nb_a_pc))

		self.cube_is_done = True

		if mpi_rank == 0: print("Sky cube loaded!")

	def LoadAllSkyMovie(self):
		if mpi_rank==0: print("Importing all sky movie data...")

		ls = sorted([f for f in os.listdir(self.image_file) if ".png" in f and len(f) == 33])

		ls = ls[::2]
		self.times = ["_".join(n.split("_")[1:3]) for n in ls]
		self.times = [dt.datetime.strptime(n, "%Y%m%d_%H%M%S") for n in self.times]
		# print(ls)

		self.Nt = len(ls)
		self.cube = np.zeros((self.Nt, self.Ne, self.Na))

		for i, f in enumerate(ls):
			if mpi_rank==0: print(str(i+1) + "/" + str(len(ls)) + ": Importing all sky image data from " + self.image_file + f)
			self.cube[i, :, :] = self.LoadAllSkyImage(filename = self.image_file + f, verbose=False)


	def LoadAllSkyImage(self, filename = "", verbose=True):

		if not filename:
			filename = self.image_file

		if mpi_rank==0 and verbose: print("Importing all sky image data from: " + filename)

		image = imageio.imread(filename)
		file = filename.split("/")[-1]
		calib_file = "/".join(filename.split("/")[:-2]) + "/" + file[:14] + file[21:-3] + "dat"
		calib = io.readsav(calib_file)

		# print("self.maps_shape", self.maps_shape)
		# print(self.azimuts * RtoD)
		# print(self.elevations * RtoD)

		map = np.zeros(self.maps_shape)
		div = np.zeros(self.maps_shape)

		for i, line in enumerate(image):
			for j, pix in enumerate(line):
				if pix > 0:
					# print("PIX", pix)
					pix_azimut    = -(float(calib["gazms"][i][j]) * DtoR + np.pi) % (2*np.pi) # GIVEN IN DEGREES
					pix_elevation = float(calib["elevs"][i][j]) * DtoR # GIVEN IN DEGREES

					if pix_azimut and pix_elevation:
						ia, ie = self.GetPixFromAzEl(pix_azimut, pix_elevation)

						# print("AA", pix_azimut*RtoD, ia)
						# print("EE", pix_elevation*RtoD, ie)

						if ia is not None and ie is not None:
							map[ie, ia] += pix * self.RtonW
							div[ie, ia] += 1



		# return map / div #[nW / m2/ sr]
		return np.divide(map, div, out = np.zeros_like(map), where = div != 0) #[nW / m2/ sr]


	def MakeSkyCubePlot(self, a_pc_list, e_pc_list, ouv_pc):
		# Plot the sky_cube as fct of time
		nb_rows, nb_cols = int(np.ceil(np.sqrt(self.Nt))), int(np.sqrt(self.Nt))
		f, axs = plt.subplots(nb_rows, nb_cols)
		c = 0

		def MakeSubplot(ax, pos, nb_rows, nb_cols, ouv_pc):
			ax = plt.subplot(nb_rows, nb_cols, pos+1, projection='polar')
			el = 90 - (self.elevations) * RtoD
			i1 = ax.pcolormesh(self.azimuts, el, self.cube[pos, :, :])
			cbar1 = plt.colorbar(i1, extend='both', spacing='proportional', shrink=0.9, ax=ax)
			ax.set_theta_zero_location("N")
			ax.set_theta_direction(-1)

			Nb_a_pc, Nb_e_pc = len(a_pc_list), len(e_pc_list)

			if Nb_a_pc == 1 and Nb_e_pc == 1:
				ax.add_artist(Ellipse((a_pc_list[0], (90 - e_pc_list[0]) * RtoD), width=ouv_pc, height=ouv_pc*RtoD, color="red"))
			elif Nb_a_pc > 1 and Nb_e_pc > 1:
				ax.plot(a_pc_list, 90 - e_pc_list * RtoD, "r")
			elif Nb_a_pc > 1:
				ax.plot(a_pc_list, [90 - e_pc_list[0] * RtoD] * Nb_a_pc, "r")
			elif Nb_e_pc > 1:
				ax.plot([a_pc_list] * Nb_e_pc, 90 - e_pc_list * RtoD, "r")
			# col.set_rlim(0,90, 1)
			# col.set_yticks(np.arange(0, 90, 20))
			# col.set_yticklabels(a.get_yticks()[::-1])   # Change the labels

		if nb_rows > 1:
			for row in axs:
				if nb_cols > 1:
					for col in row:
						if c < self.Nt:
							MakeSubplot(col, c, nb_rows, nb_cols, ouv_pc)
						c+=1
		elif nb_cols > 1:
			for col in axs:
				if c < self.Nt:
					MakeSubplot(col, c, nb_rows, nb_cols, ouv_pc)
				c+=1

		else:
			MakeSubplot(axs, c, 1, 1, ouv_pc)


	def GetBandSkyMap(self, a_band = 0, e_band = np.pi/2, length = 1000, thickness = 10, height = 50, band_I = 100, nb_sub_bands = 1):
		map = np.zeros(self.maps_shape)
		for nsb in range(nb_sub_bands):
			band_thickness = thickness / (RT) #Thickness of the band in radians, put number in km
			sub_band_thickness = thickness / (nb_sub_bands * RT) #Thickness of the band in radians, put number in km
			band_length = length / RT #Length of the band in radians, put number in km

			o = ObservationPoint(0, 0, self.h, a_band, e_band)
			lon_band_center, lat_band_center = o.P_lon, o.P_lat - band_thickness / 2 + (nsb + 1) * band_thickness / (nb_sub_bands + 1)

			def InBand(az, el):
				in_band = 0
				for alt in range(height):
					o = ObservationPoint(0, 0, self.h + alt, az, el)
					observed_lon, observed_lat = o.P_lon, o.P_lat
					if lon_band_center - band_length / 2 < observed_lon < lon_band_center + band_length / 2 and lat_band_center - sub_band_thickness / 2 < observed_lat < lat_band_center + sub_band_thickness / 2:
						in_band = o.AH_norm
						break
				return in_band

			for iel, el in enumerate(map):
				for iaz, az in enumerate(el):
					AH = InBand(self.mid_azimuts[iaz], self.mid_elevations[iel])
					if AH:
						map[iel][iaz] += band_I # / AH ** 2
		return map

	def GetPixFromEl(self, pix_e):
		"""Return the elevation index of the pixel containing a given elevation in radians. For loop on self.elevations, but works also for """
		for ie, e in enumerate(self.elevations):
			if pix_e < e:
				pix_ie = (ie - 1) #% self.Ne
				return pix_ie
		# if not np.isnan(pix_e):
		# 	return int(((pix_e - self.elevations[0]) / self.de) % self.Ne)
		# else:
		return None

	def GetPixFromAz(self, pix_a):
		"""Return the azimut index of the pixel containing a given azimut in radians."""
		mod = 2 * np.pi
		pix_a %= mod
		for ia in range(self.Na):
			if self.azimuts[ia] % mod <= pix_a < self.azimuts[ia + 1] % mod:
				return ia #% (self.Na)
			elif self.azimuts[ia + 1] >= 2*np.pi and self.azimuts[ia] <= pix_a:
				return ia
			elif self.azimuts[ia] < 0 and self.azimuts[ia + 1] > pix_a:
				return ia
		# if not np.isnan(pix_a):
		# 	return int(((pix_a - self.azimuts[0]) / self.da) % self.Na)
		# else:
		return None

	def GetPixFromAzEl(self, a, e):
		return self.GetPixFromAz(a), self.GetPixFromEl(e)

	def GetPixelSolidAngleFromiEl(self, pix_ie):
		"""Returns the solid angle from an elevation index"""
		return self.da * self.de * np.cos(self.mid_elevations[pix_ie])

	def GetPixelSolidAngle(self, pix_e):
		"""Returns the solid angle from a elevation in radian"""
		ie = self.GetPixFromEl(pix_e)
		return GetPixelSolidAngleFromiEl(ie)

	@staticmethod
	def GetArea(em, eM, h, Na):
		"""Returns the area in m**2"""
		om, oM = ObservationPoint(0, 0, h, 0, em), ObservationPoint(0, 0, h, 0, eM)
		true_el = lambda e, o: np.arcsin(o.AH_norm * np.cos(e) / (RT + h))
		tem, teM = true_el(em, om), true_el(eM, oM)

		area = 2 * np.pi * (RT + h)**2 * (np.cos(teM) - np.cos(tem)) / Na

		return area * 1e6 #in m**2


	def GetPixelArea(self, ie):
		"""Return the flattened area in m**2 of a sky pixel from its elevation index."""
		if not self.is_point_src:
			e = self.mid_elevations[ie] #Get pixel mid elevation
			em, eM = self.elevations[ie], self.elevations[ie + 1] #Get pixel min em and Max eM elevations

			area = self.GetArea(em, eM, self.h, self.Na)
			# om, oM = ObservationPoint(0, 0, self.h, 0, em), ObservationPoint(0, 0, self.h, 0, eM)
			# true_el = lambda e, o: np.arcsin(o.AH_norm * np.cos(e) / (RT + self.h))
			# tem, teM = true_el(em, om), true_el(eM, oM)
			#
			# area = 2 * np.pi * (RT + self.h)**2 * (np.cos(teM) - np.cos(tem)) / self.Na * 1e6
		else:
			area = 1

		### Deprecated version with square areas
		# de = self.de / 2. # Get pixel half elevation difference
		# AEm, AE = self.h / np.sin(em), self.h / np.sin(e) #Get lenght of segment AE and AEm for an given emission altitude
		# Hm, HM = AEm * np.sin(de) / np.sin(e), AE * np.tan(de) / np.sin(e) #Get length of segment H
		#
		# area = (Hm + HM) * self.da * AE * np.cos(e) #km2

		return area # in m2


	def GetFlux(self, az, el, r, t=0, area=1):
		"""Return the flux in nW/sr/m2 for a given az, el and opening angle (all in radians).
		The pixel units are in nW/sr/m2.
		t is the time of the map
		area is the area of the detector. Set to 1 by default."""

		if "uniform" in self.mode: #If uniform sky, juste return the brigthness in nW/sr/m2
			print("DEBUG DIRECT", float(self.mode.split("_")[1]))
			return float(self.mode.split("_")[1])


		b = 0
		tot_area = 0

		ia_min, ie_min = self.GetPixFromAzEl(az - r, el - r)
		ia_max, ie_max = self.GetPixFromAzEl(az + r, el + r)
		# ia_min, ia_max = min(ia_min, ia_max), max(ia_min, ia_max)
		# ie_min, ie_max = min(ie_min, ie_max), max(ie_min, ie_max)

		# print(ia_min, ie_min, ia_max, ie_max)
		# print(az*RtoD, el*RtoD)
		# print(ia_min, ie_min)
		# print(ia_max, ie_max)

		if (ia_max is None or ie_max is None) and (ia_min is None or ie_min is None): # If the min and max pixel do not contain the instrument
			# print("No direct flux")
			return 0

		if ia_max is None or ie_max is None:
			ia_max, ie_max = ia_min, ie_min
		if ia_min is None or ie_min is None:
			ia_min, ie_min = ia_max, ie_max

		if ia_max >= ia_min:
			ia_list = np.arange(ia_min, ia_max + 1)
			# print("COOL")
		else:
			ia_list = np.arange(ia_min, ia_max + (self.Na) + 1) % (self.Na)
			# print("NOT COOL")

		ie_list = np.arange(ie_min, ie_max + 1)
		# print(ia_list, ie_list)

		mod = 2 * np.pi
		for ia in ia_list:
			shift = 0
			if self.azimuts[(ia + 1)]%mod <= self.azimuts[ia]%mod:
				shift = mod / 2.

			pix_am = max((self.azimuts[ia] + shift) % mod, (az - r + shift) % mod)
			pix_aM = min((self.azimuts[(ia + 1)] + shift) % mod, (az + r + shift) % mod)
			pix_da = pix_aM - pix_am
			pix_da %= mod

			# print("DEBUG DIRECT:", pix_da%mod*RtoD)
			# print((ia)%(self.Na-1), (ia + 1)%(self.Na-1))
			# print("AZ", (self.azimuts[(ia + 1)%(self.Na-1)] + shift) % mod*RtoD, (self.azimuts[ia] + shift) % mod*RtoD)
			# print("AZ", (az + r + shift) % mod*RtoD, (az - r + shift) % mod*RtoD)

			for ie in ie_list:
				# print("EL", self.elevations[ie + 1]*RtoD, (el + r)*RtoD, self.elevations[ie]*RtoD, (el - r)*RtoD)
				pix_em = max(self.elevations[ie], el - r)
				pix_eM = min(self.elevations[ie + 1], el + r)
				pix_de = pix_eM - pix_em
				# print(self.cube[t, ie, ia], self._GetPixelSolidAngleFromiEl(ie), self.cube[t, ie, ia] * self._GetPixelSolidAngleFromiEl(ie))

				pix_area = self.GetArea(pix_em, pix_eM, self.h, 2 * np.pi / pix_da)
				b += self.cube[t, ie, ia] * pix_area
				tot_area += pix_area
				# b += self.cube[t, ie, ia] * pix_da * pix_de * np.cos(self.mid_elevations[ie])


		# if mpi_rank==0: print("DEBUG DIRECT:", b / tot_area)

		return b / tot_area


	def SetLightParameters(self):

		self.total_scattering_map 	= 2 * self.V_map

		self.DoLP_map				= 2 * np.sqrt(self.Vcos_map**2 + self.Vsin_map**2) * 100 # DoLP of scattered light from (e,a)
		self.DoLP_map 				= np.divide(self.DoLP_map, self.V_map, out = np.zeros_like(self.DoLP_map), where = self.V_map != 0)

		self.AoRD_map 				= np.arctan2(self.Vsin_map, self.Vcos_map) / 2. # Angle of polaisation of light from (e,a)

		self.scattering_map 		= self.total_scattering_map * self.DoLP_map / 100. # Polarized intensity from (e, a) reaching us

		# print("self.scattering_map", self.scattering_map)
		# print(np.where(self.scattering_map < 0))

		for it, t in enumerate(self.V_map):
			for ie, e in enumerate(t):
				for ia, a in enumerate(e):
					self.V_total[it, ie, ia] = np.sum(self.V_map[it, ie, ia, :, :])
					self.Vcos_total[it, ie, ia] = np.sum(self.Vcos_map[it, ie, ia, :, :])
					self.Vsin_total[it, ie, ia] = np.sum(self.Vsin_map[it, ie, ia, :, :])

					# self.I0_total[it, ie, ia] = 2 * self.V_total[it, ie, ia]
					# self.DoLP_total[it, ie, ia] = 100 * 2 * np.sqrt(self.Vcos_total[it, ie, ia] ** 2 + self.Vsin_total[it, ie, ia] ** 2) / self.V_total[it, ie, ia] #in %
					# self.AoLP_total[it, ie, ia] = np.arctan2(self.Vsin_total[it, ie, ia], self.Vcos_total[it, ie, ia]) / 2

		self.I0_total = 2 * self.V_total
		self.DoLP_total = 100 * 2 * np.sqrt(self.Vcos_total ** 2 + self.Vsin_total ** 2) / self.V_total[it, ie, ia] #in %
		self.AoLP_total = np.arctan2(self.Vsin_total, self.Vcos_total) / 2


	def GetFlatMap(self, A_lon, A_lat, mid_lon = None, mid_lat = None, time = 0, Nmax = 1000):
		if mid_lon is None or mid_lat is None:
			lon_m = AzEltoLonLat(-90*DtoR, 0, self.h, A_lon, A_lat)[0]
			lon_M = AzEltoLonLat(90*DtoR, 0, self.h, A_lon, A_lat)[0]
			lat_m = AzEltoLonLat(180*DtoR, 0, self.h, A_lon, A_lat)[1]
			lat_M = AzEltoLonLat(0*DtoR, 0, self.h, A_lon, A_lat)[1]

			Nlon = Nlat = int(np.sqrt(Nmax))
			mid_longitudes, dlon = np.linspace(lon_m, lon_M, Nlon, retstep = True)
			mid_latitudes, dlat  = np.linspace(lat_m, lat_M, Nlat, retstep = True)

		else:
			Nlon = len(mid_lon)
			Nlat = len(mid_lat)
			mid_longitudes = mid_lon
			mid_latitudes  = mid_lat

		flat_map = np.zeros((Nlat, Nlon))

		for ie, e in enumerate(self.mid_elevations):
			for ia, a in enumerate(self.mid_azimuts):
				sky_lon, sky_lat = AzEltoLonLat(a, e, self.h, A_lon, A_lat)
				I0 = self.cube[time, ie, ia]
				sky_area = self.GetPixelArea(ie)
				for ilon, ground_lon in enumerate(mid_longitudes):
					for ilat, ground_lat in enumerate(mid_latitudes):
						L = np.sqrt((sky_lon - ground_lon)**2 + (sky_lat - ground_lat)**2) * RT * np.cos(sky_lat)
						theta = np.arctan(L / self.h)
						I = I0 * sky_area * np.cos(theta)**3 / (2 * np.pi * self.h**2 * 1e6) #nW/m2/sr

						flat_map[ilat, ilon] += I

		# if mpi_rank == 0: plt.imshow(flat_map)
		# if mpi_rank == 0: plt.show()

		return flat_map, mid_longitudes, mid_latitudes
