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

class SkyMap:
	def __init__(self, in_dict):

		self.path = in_dict["sky_path"]
		self.file = in_dict["sky_file"]

		self.h = float(in_dict["emission_altitude"])

		self.exist = True
		self.mode = in_dict["sky_mode"].lower()

		if "band_e" in self.mode:
			self.band_elevation = float(self.mode[-2:]) * DtoR
			# print(self.mode[-2:], self.band_elevation)
			self.mode = self.mode[:-2]
		elif self.mode == "uniform_0" or self.mode == "none":
			self.mode = "none"
			self.exist = False


		#Zone of emission
		self.I_zone_a_min = float(in_dict["I_zone_a_min"])	* DtoR
		self.I_zone_a_max = float(in_dict["I_zone_a_max"]) 	* DtoR
		self.I_zone_e_min = float(in_dict["I_zone_e_min"]) 	* DtoR
		self.I_zone_e_max = float(in_dict["I_zone_e_max"]) 	* DtoR


		### Initialize some parameters for the case h!=0
		self.N = float(in_dict["Nb_emission_points"]) # Number of bins to compute (if h!=0)
		self.Na, self.Ne = int(2 * np.sqrt(self.N)), int(np.sqrt(self.N) / 2) # Numbers of bins in azimut/elevation. Na = 4*Ne and Na*Ne = N
		self.da = (self.I_zone_a_max - self.I_zone_a_min) / self.Na # Length of an azimut bin
		self.de = (self.I_zone_e_max - self.I_zone_e_min) / self.Ne # Length of an elevation bin

		if self.exist and "moving" in self.mode:
			self.Nt = int(in_dict["Nb_emission_maps"])
		else:
			self.Nt = 1
		self.is_time_dependant = self.Nt > 1

		# self.azimuts = np.linspace(self.I_zone_a_min + self.da/2., self.I_zone_a_max + self.da/2, self.Na) # array of all azimuts
		# self.elevations = np.linspace(self.I_zone_e_min + self.de/2., self.I_zone_e_max - self.de/2., self.Ne-1) # array of all elevations
		self.azimuts = np.linspace(self.I_zone_a_min + self.da/2, self.I_zone_a_max + self.da/2, self.Na, endpoint=True) # array of all azimuts
		self.elevations = np.linspace(self.I_zone_e_min, self.I_zone_e_max, self.Ne, endpoint=True) # array of all elevations

		self.mid_azimuts = self.azimuts + self.da/2.
		self.mid_elevations = self.elevations + self.de/2.

		# print("DEBUG SKYMAP: elevations", self.elevations*RtoD)
		# print("DEBUG SKYMAP: azimuts", self.azimuts*RtoD)

		self.maps_shape = (len(self.elevations), len(self.azimuts))
		# self.sky_I_map = np.zeros(self.maps_shape) # Intensity of the emission at point (e,a)

		self.cube_is_done = False



	def LoadSkyEmmisionsCube(self, Nb_a_pc, Nb_e_pc):
		"""In the case where h!=0, load a cube of data corresponding to the state of the sky at different times, depending on the type of sky we want"""

		print("Loading the sky emission cube...")

		self.cube = np.zeros((self.Nt, self.maps_shape[0], self.maps_shape[1]))
		self.cube_is_done = True
		if "band" in self.mode:
			nb_sub_bands = int(self.mode[:1])

		if "moving_ew_band" in self.mode:
			for ie, e in enumerate(np.linspace(0, np.pi, self.Nt + 2)[1:-1]): #elevation of the auroral band center
				self.cube[ie, :, :] = self.GetBandSkyMap(a_band = 0, e_band = e, length = 1000, thickness = 50, height = 40, band_I = 100, nb_sub_bands = nb_sub_bands)

		elif "band_e" in self.mode:
			self.cube[0, :, :] = self.GetBandSkyMap(a_band = 0, e_band = self.band_elevation, length = 1000, thickness = 10, height = 50, band_I = 100, nb_sub_bands = nb_sub_bands)

		elif self.mode[:7] == "uniform":
			self.cube[:, :, :] += float(self.mode[8:])

		elif self.mode[:4] == "spot":
			a = int(self.mode[6:8]) * DtoR
			e = int(self.mode[-2:]) * DtoR
			self.cube[0, :, :] = self.GetBandSkyMap(a_band = a, e_band = e, length = 50, thickness = 50, height = 50, band_I = 100, nb_sub_bands = 1)

		self.data_shape = (self.Nt, Nb_e_pc, Nb_a_pc, len(self.elevations), len(self.azimuts))

		self.scattering_map 		= np.zeros(self.data_shape) # Polarized intensity from (e, a) reaching us
		self.DoLP_map				= np.zeros(self.data_shape) # DoLP of scattered light from (e,a)
		self.total_scattering_map 	= np.zeros(self.data_shape) # Intensity from (e, a) reaching us
		self.AoRD_map 				= np.zeros(self.data_shape) # Angle of polaisation of light from (e,a)

		self.cube_is_done = True


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
				ax.plot([a_pc_list] * Nb_e_pc, 90 - e_pc_list[0] * RtoD, "r")
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
						map[iel][iaz] += band_I / AH ** 2
		return map

	def GetPixFromEl(self, pix_e):
		"""Return the index of pixels containing a given elevation in radians."""
		for ie, e in enumerate(self.elevations):
			if pix_e < e:
				pix_ie = (ie - 1)%self.Ne
				return pix_ie

	def GetPixFromAz(self, pix_a):
		"""Return the index of pixels containing a given elevation in radians."""
		if pix_a < 0:
			pix_a += 360*DtoR
		for ia, a in enumerate(self.azimuts):
			if pix_a < a:
				pix_ia = (ia - 1)%self.Na
				return pix_ia

	def GetPixFromAzEl(self, a, e):
		return self.GetPixFromAz(a), self.GetPixFromEl(e)

	def GetPixelSolidAngleFromiEl(self, pix_ie):
		"""Returns the solid angle from an elevation index"""
		return self.da * self.de * np.cos(self.mid_elevations[pix_ie])

	def GetPixelSolidAngle(self, pix_e):
		"""Returns the solid angle from a elevation in radian"""
		ie = self.GetPixFromEl(pix_e)
		return GetPixelSolidAngleFromiEl(ie)


	def GetFlux(self, az, el, r, t=0, area=1):
		"""Return the flux in nW for a given az, el and opening angle (all in radians).
		The pixel units are in nW/sr/m2.
		t is the time of the map
		area is the area of the detector. Set to 1 by default."""
		b = 0

		ia_min, ie_min = self.GetPixFromAzEl(az - r, el - r)
		ia_max, ie_max = self.GetPixFromAzEl(az + r, el + r)
		# ia_min, ia_max = min(ia_min, ia_max), max(ia_min, ia_max)
		# ie_min, ie_max = min(ie_min, ie_max), max(ie_min, ie_max)

		# print(ia_min, ie_min, ia_max, ie_max)
		if ia_max >= ia_min:
			ia_list = np.arange(ia_min, ia_max + 1)
		else:
			ia_list = np.arange(ia_min, ia_max + self.Na + 1) % self.Na

		ie_list = np.arange(ie_min, ie_max + 1)
		print(ia_list, ie_list)

		mod = 360*DtoR
		for ia in ia_list:
			shift = 0
			if self.azimuts[(ia + 1)%(self.Na-1)]%mod < self.azimuts[ia]%mod:
				shift = mod / 2.
			pix_da = min((self.azimuts[(ia + 1)%(self.Na-1)] + shift) % mod, (az + r + shift) % mod) - max((self.azimuts[ia] + shift) % mod, (az - r + shift) % mod)
			pix_da %= mod
			print("DEBUG DIRECT:", pix_da%mod*RtoD)
			print((ia)%(self.Na-1), (ia + 1)%(self.Na-1))
			print("AZ", self.azimuts[(ia)%(self.Na-1)]%mod*RtoD, (self.azimuts[(ia + 1)%(self.Na-1)]%mod*RtoD))
			print("AZ", self.azimuts[(ia + 1)%(self.Na-1)]%mod*RtoD, (az + r)%mod*RtoD, self.azimuts[ia]%mod*RtoD, (az - r)%mod*RtoD)
			for ie in ie_list:
				print("EL", self.elevations[ie + 1]*RtoD, (el + r)*RtoD, self.elevations[ie]*RtoD, (el - r)*RtoD)
				pix_de = min(self.elevations[ie + 1], el + r) - max(self.elevations[ie], el - r)
				# print(self.cube[t, ie, ia], self._GetPixelSolidAngleFromiEl(ie), self.cube[t, ie, ia] * self._GetPixelSolidAngleFromiEl(ie))
				b += self.cube[t, ie, ia] * pix_da * pix_de * np.cos(self.mid_elevations[ie])


		# print("DEBUG DIRECT:", b, area, b * area)

		return b * area
