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

class GroundMap:
	def __init__(self, in_dict):

		self.path = in_dict["ground_path"]
		self.file = in_dict["ground_file"]

		self.location = in_dict["location"]
		self.A_lon, self.A_lat 	= GetLonLatFromName(self.location)

		self.radius = float(in_dict["ground_emission_radius"]) / RT

		self.exist = (self.radius > 0 and self.file)




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

		self.longitudes = np.linspace(lon_min, lon_max, nb_lon) # list of pixel longitudes
		self.latitudes = np.linspace(lat_max, lat_min, nb_lat) # list of pixel latitudes

		map_band = map.GetRasterBand(1)
		self.I_map = map_band.ReadAsArray(LonToCol(lon_min), LatToRow(lat_max), nb_lon, nb_lat) # Emission map we will use

		self.maps_shape = (Nb_e_pc, Nb_a_pc, len(self.latitudes), len(self.longitudes))

		self.scattering_map 		= np.zeros(self.maps_shape) # Intensity from (e, a) reaching us
		self.DoLP_map				= np.zeros(self.maps_shape) # Polarized intensity from (e, a) reaching us
		self.total_scattering_map 	= np.zeros(self.maps_shape) # DoLP of scattered light from (e,a)
		self.AoRD_map 				= np.zeros(self.maps_shape) # Angle of polaisation of light from (e,a)

	def LoadPointSourceMap(self, src_I0, src_lon, src_lat, Nb_a_pc, Nb_e_pc):

		self.longitudes = np.linspace(src_lon, src_lon, 1) # list of pixel longitudes
		self.latitudes = np.linspace(src_lat, src_lat, 1) # list of pixel latitudes

		self.I_map = np.array([src_I0])

		self.maps_shape = (Nb_e_pc, Nb_a_pc, 1, 1)

		self.scattering_map 		= np.zeros(self.maps_shape) # Intensity from (e, a) reaching us
		self.DoLP_map				= np.zeros(self.maps_shape) # Polarized intensity from (e, a) reaching us
		self.total_scattering_map 	= np.zeros(self.maps_shape) # DoLP of scattered light from (e,a)
		self.AoRD_map 				= np.zeros(self.maps_shape) # Angle of polaisation of light from (e,a)
