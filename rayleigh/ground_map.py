#!/usr/bin/python3
# -*-coding:utf-8 -*

import sys as sys
import numpy as np
import cv2
import time as tm
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Circle
from matplotlib.patches import Ellipse
from matplotlib.patches import Arrow
# from matplotlib.lines import mlines

import imageio

from observation import *
from rayleigh_utils import *



class ElevationMap():
	def __init__(self, in_dict):
		self.path = in_dict["alt_map_path"]

		self.location = in_dict["location"]
		self.A_lon, self.A_lat = GetLonLatFromName(self.location)
		self.A_lon, self.A_lat = self.A_lon, self.A_lat
		print(self.location, self.A_lon, self.A_lat)

		self.radius = float(in_dict["alt_map_radius"]) / RT
		self.N_bins_max = int(in_dict["ground_N_bins_max"])

		print(self.radius)

		self.LoadData()

		self.exist = (self.radius > 0)


	def GetNameFromLonLat(self, lon, lat):
		lon = int(lon)
		lat = int(lat)
		lon5 = lon - lon % 5
		lat5 = lat - lat % 5

		folder_name = self.FileNameCoordinates(lon5, lat5) + "_" + self.FileNameCoordinates(lon5+5, lat5+5)
		file_name = "/ALPSMLC30_" + self.FileNameCoordinates(lon, lat) + "_DSM.tif"

		return folder_name + file_name

	def FileNameCoordinates(self, lon, lat):
		name = ""
		if lat >= 0:
			name += "N"
		else:
			name += "S"
		name += str(int(lat)).zfill(3)

		if lon >= 0:
			name += "E"
		else:
			name += "W"
		name += str(int(lon)).zfill(3)

		return name

	def LoadData(self):
		min_lon, min_lat = self.A_lon - self.radius, self.A_lat - self.radius
		max_lon, max_lat = self.A_lon + self.radius, self.A_lat + self.radius
		print(min_lon, min_lat, max_lon, max_lat)

		lon_name_list = range(int(min_lon * RtoD), int(max_lon * RtoD) + 1,  1)
		lat_name_list = range(int(max_lat * RtoD), int(min_lat * RtoD) - 1, -1)

		self.map = None
		for n_lon in lon_name_list:
			lon_band = None
			old_tile_nb_lon = 0
			for n_lat in lat_name_list:

				file_name = self.GetNameFromLonLat(n_lon, n_lat)
				print(file_name)
				tile = gdal.Open(self.path + file_name, gdal.GA_ReadOnly) #get raster objsect from tif file

				nb_lon = tile.RasterXSize
				nb_lat = tile.RasterYSize

				tile_info = tile.GetGeoTransform()

				print("tile info", tile_info)
				origin_lon = tile_info[0]	 * DtoR
				origin_lat = tile_info[3]	 * DtoR
				pixel_width = tile_info[1]	 * DtoR
				pixel_height = tile_info[5]	 * DtoR

				LonToCol = lambda l: int((l - origin_lon) / pixel_width)
				LatToRow = lambda l: int((l - origin_lat) / pixel_height)

				lon_min = max(self.A_lon - self.radius, origin_lon)
				lon_max = min(self.A_lon + self.radius, origin_lon + pixel_width * nb_lon)

				lat_min = max(self.A_lat - self.radius, origin_lat + pixel_height * nb_lat)
				lat_max = min(self.A_lat + self.radius, origin_lat)

				nb_lon  = int((lon_max - lon_min) / pixel_width)
				nb_lat  = int((lat_min - lat_max) / pixel_height)

				print("GROUND MAP", nb_lon, nb_lat)

				self.longitudes, self.dlon = np.linspace(lon_min, lon_max, nb_lon, retstep=True) # list of pixel longitudes
				self.latitudes, self.dlat = np.linspace(lat_max, lat_min, nb_lat, retstep=True) # list of pixel latitudes

				tile = tile.GetRasterBand(1).ReadAsArray(LonToCol(lon_min), LatToRow(lat_max), nb_lon, nb_lat)

				print(tile)
				tile_nb_lat, new_tile_nb_lon = tile.shape

				if old_tile_nb_lon != 0:
					if old_tile_nb_lon != new_tile_nb_lon:
						tile = cv2.resize(tile, dsize=(old_tile_nb_lon, tile_nb_lat), interpolation = cv2.INTER_CUBIC)

				else:
					old_tile_nb_lon = new_tile_nb_lon

				if lon_band is not None:
					lon_band = np.append(lon_band, tile, axis=0)
				else:
					lon_band = tile


			if self.map is not None and lon_band.size > 0:
				self.map = np.append(self.map, lon_band, axis=1)
			else:
				self.map = lon_band

			print("MAP SHAPE", self.map.shape)


		self.map = np.flip(self.map, axis=0)
		nb_lat, nb_lon = self.map.shape

		self.longitudes, self.dlon = np.linspace(min_lon, max_lon, nb_lon, endpoint=True, retstep = True)
		self.latitudes, self.dlat = np.linspace(max_lat, min_lat, nb_lat, endpoint=True, retstep = True)

		print(self.longitudes, self.latitudes)



		if nb_lon * nb_lat > self.N_bins_max and self.N_bins_max > 0:
			print("GROUND MAP TOO BIG", nb_lon * nb_lat, self.N_bins_max)
			new_nb_lon, new_nb_lat = int(np.sqrt(self.N_bins_max)), int(np.sqrt(self.N_bins_max))
			print("GROUND MAP", new_nb_lon, new_nb_lat)

			self.map = cv2.resize(self.map, dsize=(new_nb_lon, new_nb_lat), interpolation = cv2.INTER_CUBIC)
			nb_lat, nb_lon = self.map.shape

			self.longitudes, self.dlon = np.linspace(min_lon, max_lon, nb_lon, retstep=True) # list of pixel longitudes
			self.latitudes, self.dlat = np.linspace(max_lat, min_lat, nb_lat, retstep=True) # list of pixel latitudes



		self.exist = True
		if not self.map.shape:
			self.exist = True


	def PlotMap(self):
		f = plt.figure()
		plt.pcolormesh(self.map)
		# plt.pcolormesh(np.flip(self.map, axis=0))
		# plt.pcolormesh(np.flip(np.rot90(self.map, k=-1)))
		plt.show()


class GroundMap:
	def __init__(self, in_dict):

		self.path = in_dict["ground_path"]
		self.file = in_dict["ground_file"]

		self.location = in_dict["location"]
		self.A_lon, self.A_lat 	= GetLonLatFromName(self.location)

		self.radius = float(in_dict["ground_emission_radius"]) / RT
		self.N_bins_max = int(in_dict["ground_N_bins_max"])

		self.exist = ((self.radius > 0 and self.file) or float(in_dict["point_src_I0"]) > 0)
		self.is_point_source = float(in_dict["point_src_I0"]) > 0


		self.src_I0, self.src_az, self.src_dist = None, None, None

	def LoadGroundEmmisionsMap(self, Nb_a_pc, Nb_e_pc):
		"""In the case where h==0, load the emission map from the geotiff files"""
		map = gdal.Open(self.path + self.file, gdal.GA_ReadOnly)
		nb_lon = map.RasterXSize
		nb_lat = map.RasterYSize

		nb_bands = map.RasterCount

		map_info = map.GetGeoTransform()
		origin_lon = map_info[0]	 * DtoR
		origin_lat = map_info[3]	 * DtoR
		pixel_width = map_info[1]	 * DtoR
		pixel_height = map_info[5]	 * DtoR

		LonToCol = lambda lon: int((lon - origin_lon) / pixel_width)
		LatToRow = lambda lat: int((lat - origin_lat) / pixel_height)

		lon_min = max(self.A_lon - self.radius, origin_lon)
		lon_max = min(self.A_lon + self.radius, origin_lon + pixel_width * nb_lon)

		lat_min = max(self.A_lat - self.radius, origin_lat + pixel_height * nb_lat)
		lat_max = min(self.A_lat + self.radius, origin_lat)

		nb_lon  = int( 2 * self.radius / pixel_width)
		nb_lat  = int(-2 * self.radius / pixel_height)

		print("GROUND MAP", nb_lon, nb_lat)

		self.longitudes, self.dlon = np.linspace(lon_min, lon_max, nb_lon, retstep=True) # list of pixel longitudes
		self.latitudes, self.dlat = np.linspace(lat_max, lat_min, nb_lat, retstep=True) # list of pixel latitudes

		map_band = map.GetRasterBand(1)
		self.I_map = map_band.ReadAsArray(LonToCol(lon_min), LatToRow(lat_max), nb_lon, nb_lat) # Emission map we will use

		if nb_lon * nb_lat > self.N_bins_max and self.N_bins_max > 0:
			print("GROUND MAP TOO BIG", nb_lon * nb_lat, self.N_bins_max)
			new_nb_lon, new_nb_lat = int(np.sqrt(self.N_bins_max)), int(np.sqrt(self.N_bins_max))
			print("GROUND MAP", new_nb_lon, new_nb_lat)

			self.longitudes, self.dlon = np.linspace(lon_min, lon_max, new_nb_lon, retstep=True) # list of pixel longitudes
			self.latitudes, self.dlat = np.linspace(lat_max, lat_min, new_nb_lat, retstep=True) # list of pixel latitudes

			self.I_map = cv2.resize(self.I_map, dsize=(new_nb_lon, new_nb_lat), interpolation = cv2.INTER_CUBIC)


		self.maps_shape = (Nb_e_pc, Nb_a_pc, len(self.latitudes), len(self.longitudes))

		self.scattering_map 		= np.zeros(self.maps_shape) # Intensity from (e, a) reaching us
		self.DoLP_map				= np.zeros(self.maps_shape) # Polarized intensity from (e, a) reaching us
		self.total_scattering_map 	= np.zeros(self.maps_shape) # DoLP of scattered light from (e,a)
		self.AoRD_map 				= np.zeros(self.maps_shape) # Angle of polaisation of light from (e,a)

	def GetArea(self, ilat):
		"""Return the area of a pixel on the map in km**2. If we use a point source, the area is set to one."""
		if not self.is_point_source:
			return RT ** 2 * abs(self.dlat) * abs(self.dlon) * np.cos(self.latitudes[ilat])
		else:
			return 1


	def LoadPointSourceMap(self, src_I0, src_az, src_dist, Nb_a_pc, Nb_e_pc):

		self.src_I0, self.src_az, self.src_dist = src_I0, src_az, src_dist

		src_lon, src_lat = AzDistToLonLat(src_az, src_dist, self.A_lon, self.A_lat)

		self.longitudes = np.linspace(src_lon, src_lon, 1) # list of pixel longitudes
		self.latitudes = np.linspace(src_lat, src_lat, 1) # list of pixel latitudes

		self.I_map = np.ones((1, 1)) * src_I0

		self.maps_shape = (Nb_e_pc, Nb_a_pc, 1, 1)

		self.scattering_map 		= np.zeros(self.maps_shape) # Polarized intensity from (e, a) reaching us
		self.DoLP_map				= np.zeros(self.maps_shape) # DoLP of scattered light from (e,a)
		self.total_scattering_map 	= np.zeros(self.maps_shape) # Intensity from (e, a) reaching us
		self.AoRD_map 				= np.zeros(self.maps_shape) # Angle of polaisation of light from (e,a)
