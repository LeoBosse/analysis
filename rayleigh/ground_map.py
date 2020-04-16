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

		# self.in_dict = in_dict

		self.path = in_dict["alt_map_path"]

		self.location = in_dict["location"]
		self.A_lon, self.A_lat = GetLonLatFromName(self.location)
		# self.A_lon, self.A_lat = self.A_lon, self.A_lat
		print(self.location, self.A_lon, self.A_lat)

		self.nb_pix_max = int(in_dict["alt_map_N_bins_max"])
		self.resolution = float(in_dict["alt_map_resolution"]) / RT #in radian
		if self.resolution == 0:
			self.resolution = (1/3600)*DtoR # 1 arcsec

		self.radius = float(in_dict["alt_map_radius"]) / RT
		self.N_bins_max = int(in_dict["ground_N_bins_max"])
		self.lon_min, self.lon_max = self.A_lon - self.radius, self.A_lon + self.radius
		self.lat_min, self.lat_max = self.A_lat - self.radius, self.A_lat + self.radius


		print(self.radius)

		self.SetMapProperties()

		self.LoadMultipleTiles()

		self.A_alt = self.GetAltitudeFromLonLat(self.A_lon, self.A_lat)

		self.exist = (self.radius > 0)

	def SetMapProperties(self):

		self.nb_pix_lon = int((self.lon_max - self.lon_min) / self.resolution) + 1
		self.nb_pix_lat = int((self.lat_max - self.lat_min) / self.resolution) + 1
		self.nb_pix 	= self.nb_pix_lon * self.nb_pix_lat

		if self.nb_pix > self.nb_pix_max and self.nb_pix_max > 0:
			self.nb_pix_lon = int(np.sqrt(self.nb_pix_max))
			self.nb_pix_lat = int(np.sqrt(self.nb_pix_max))
			self.nb_pix 	= self.nb_pix_lon * self.nb_pix_lat
			self.resolution = (self.lon_max - self.lon_min) / self.nb_pix_lon

		self.lon_min -= self.lon_min % self.resolution
		self.lon_max += self.resolution - self.lon_max % self.resolution

		self.lat_min -= self.lat_min % self.resolution
		self.lat_max += self.resolution - self.lat_max % self.resolution


		self.map = np.zeros((self.nb_pix_lat, self.nb_pix_lon))

		self.longitudes = np.linspace(self.lon_min, self.lon_max, self.nb_pix_lon)
		self.latitudes  = np.linspace(self.lat_max, self.lat_min, self.nb_pix_lon)

		self.lon_name_list = range(int(self.lon_min * RtoD), int(self.lon_max * RtoD) + 1,  1)
		self.lat_name_list = range(int(self.lat_max * RtoD), int(self.lat_min * RtoD) - 1, -1)

		self.nb_tile_lon = len(self.lon_name_list)
		self.nb_tile_lat = len(self.lat_name_list)

		self.nb_tiles = self.nb_tile_lon * self.nb_tile_lat

		self.tiles_coord = []

		for n_lon in self.lon_name_list:
			for n_lat in self.lat_name_list:
				self.tiles_coord.append((n_lon, n_lat))

	def IsInMap(self, lon, lat):
		print(self.lon_min ,lon ,self.lon_max ,self.lat_min ,lat ,self.lat_max)

		if self.lon_min < lon < self.lon_max and self.lat_min < lat < self.lat_max:
			return True
		else:
			return False


	def GetRelativeAltitude(self, lon, lat, unit=1000): #Elevation map is in meters
	"""Altitude of given coord above the instrument in A (center of the map).
	return result(in meter)/unit. (unit=1000:result in km)"""
		if self.IsInMap(lon, lat):
			l, c = self.GetIndexFromLonLat(lon, lat)
			return (self.map[l, c] - self.A_alt) / unit
		else:
			return - self.A_alt / unit

	def GetAltitudeFromLonLat(self, lon, lat, unit=1000):
		"""Altitude of given coord above sea level.
		return result(in meter)/unit. (unit=1000:result in km)"""
		if self.IsInMap(lon, lat):
			l, c = self.GetIndexFromLonLat(lon, lat)
			return self.map[l, c] / unit
		else:
			return 0

	def GetAltitudeAboveGround(self, lon, lat, alt, unit=1000):
		"""For a given altitude, returns the relative altitude to the given coord altitude.
		return result(in meter)/unit. (unit=1000:result in km)"""
		if self.IsInMap(lon, lat):
			l, c = self.GetIndexFromLonLat(lon, lat)
			return alt - self.map[l, c] / unit
		else:
			return alt


	def IsVisible(self, lon1, lat1, alt1=None, lon2=None, lat2=None, alt2=None, dlos=0.05):
		"""Return True if two points are visible. (if a straight line between them is always above ground.).Else returns False.
		dlos is the distance between each points along the line. (too big and it is long, too low and you will miss collisions.)
		If alt1==None, point 1 on ground.
		If point 2 not given, use instrument coord and alt.
		If point 2 coord are given but alt2==None, use point 2 on ground."""
		if not lon2:
			lon2 = self.A_lon
			lat2 = self.A_lat
			alt2 = self.A_alt
		if not alt2:
			alt2 = self.GetAltitudeFromLonLat(lon2, lat2)
		if not alt1:
			alt1 = self.GetAltitudeFromLonLat(lon1, lat1)

		







	def GetIndexFromLonLat(self, lon, lat):
		line = int((self.lat_max - lat) / self.resolution)
		col  = int((lon - self.lon_min) / self.resolution)
		return line, col


	def GetLonLatFromindex(self, l, c):
		lon = self.longitudes[c]
		lat = self.latitude[l]

		return lon, lat


	def LoadSingleTile(self, n_lon, n_lat):
		file_name = self.GetNameFromLonLat(n_lon, n_lat)
		print(file_name)

		try:
			tile = gdal.Open(self.path + file_name, gdal.GA_ReadOnly) #get raster objsect from tif file
		except:
			print("WARNING: Elevation tile not found:", self.path + file_name)
			return np.zeros((1,1)), self.A_lon, self.A_lat

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

		tile = tile.GetRasterBand(1).ReadAsArray(LonToCol(lon_min), LatToRow(lat_max), nb_lon, nb_lat)

		if pixel_width != self.resolution:
			nb_lon  = int((lon_max - lon_min) / self.resolution)
			nb_lat  = int((lat_max - lat_min) / self.resolution)

			tile = cv2.resize(tile, dsize=(nb_lon, nb_lat), interpolation = cv2.INTER_CUBIC)

		return tile, lon_min, lat_max

	def LoadMultipleTiles(self):

		for tile_lon, tile_lat in self.tiles_coord:
			tile, lon_min, lat_max = self.LoadSingleTile(tile_lon, tile_lat)

			tile_line, tile_col = self.GetIndexFromLonLat(lon_min, lat_max)
			nb_lat, nb_lon = tile.shape

			print("tile.shape", tile.shape)
			print("map chunk shape", self.map[tile_line : tile_line + nb_lat, tile_col : tile_col + nb_lon].shape)

			print(self.map.shape, tile_line, tile_line + nb_lat, tile_col, tile_col + nb_lon)

			self.map[tile_line : tile_line + nb_lat, tile_col : tile_col + nb_lon] = tile


	def PlotMap(self):
		f = plt.figure()
		plt.pcolormesh((self.longitudes - self.A_lon) * RT, (self.latitudes - self.A_lat) * RT, self.map)
		plt.colorbar()
		# plt.pcolormesh(np.flip(self.map, axis=0))
		# plt.pcolormesh(np.flip(np.rot90(self.map, k=-1)))
		plt.show()

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
