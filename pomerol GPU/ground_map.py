#!/usr/bin/python3
# -*-coding:utf-8 -*
from mpi4py import MPI
mpi_comm = MPI.COMM_WORLD
mpi_rank = mpi_comm.Get_rank()
mpi_size = mpi_comm.Get_size()
mpi_name = mpi_comm.Get_name()

import sys as sys
import numpy as np
import cv2
import time as tm
import scipy as sci
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Circle
from matplotlib.patches import Ellipse
from matplotlib.patches import Arrow
# from matplotlib.lines import mlines

import imageio

import osgeo.gdal as gdal
gdal.UseExceptions()  # not required, but a good idea

from pomerol_configuration import *
from observation import *
from rayleigh_utils import *



class ElevationMap():
	def __init__(self, in_dict):

		# self.in_dict = in_dict

		self.path = pomerol_configuration.alt_map_path

		self.location = in_dict["location"]
		self.A_lon, self.A_lat = GetLonLatFromName(self.location)
		# self.A_lon, self.A_lat = self.A_lon, self.A_lat
		# print(self.location, self.A_lon, self.A_lat)

		self.nb_pix_max = int(in_dict["alt_map_N_bins_max"])
		self.resolution = float(in_dict["alt_map_resolution"]) #in km
		if self.resolution == 0:
			self.resolution = 0.031 # 1 arcsec = 31m
		self.lon_resolution = DistToAngle(self.resolution, lat = self.A_lat, az = 90*DtoR)
		self.lat_resolution = DistToAngle(self.resolution, lat = self.A_lat, az =  0*DtoR)

		self.radius = float(in_dict["alt_map_radius"])
		self.N_bins_max = int(in_dict["ground_N_bins_max"])
		self.lon_min, self.lon_max = self.A_lon - DistToAngle(self.radius, lat = self.A_lat, az = 90*DtoR), self.A_lon + DistToAngle(self.radius, lat = self.A_lat, az = 90*DtoR)
		self.lat_min, self.lat_max = self.A_lat - DistToAngle(self.radius), self.A_lat + DistToAngle(self.radius)

		# print("ALT MAP BOUNDS", self.lon_min*RtoD, self.lon_max*RtoD, self.lat_min*RtoD, self.lat_max*RtoD)

		# print(self.radius)

		self.SetMapProperties()

		# print("ALT MAP BOUNDS", self.lon_min*RtoD, self.lon_max*RtoD, self.lat_min*RtoD, self.lat_max*RtoD)
		self.LoadMultipleTiles()
		# print("ALT MAP BOUNDS", self.lon_min*RtoD, self.lon_max*RtoD, self.lat_min*RtoD, self.lat_max*RtoD)

		self.A_alt = self.GetAltitudeFromLonLat(self.A_lon, self.A_lat)

		self.min_alt = np.min(self.map)
		self.max_alt = np.max(self.map)

		self.exist = (self.radius > 0)

	def SetMapProperties(self):

		### Compute numbers of lines, col and pixels.
		self.nb_pix_lon = int((self.lon_max - self.lon_min) / self.lon_resolution) + 1
		self.nb_pix_lat = int((self.lat_max - self.lat_min) / self.lat_resolution) + 1
		self.nb_pix 	= self.nb_pix_lon * self.nb_pix_lat
		# print("ALT MAP RESOLUTION lon, lat", self.lon_resolution, self.lat_resolution)

		### If a pix number maximum is specified (> 0), limit the map to that number
		if self.nb_pix > self.nb_pix_max and self.nb_pix_max > 0:
			self.nb_pix_lon = int(np.sqrt(self.nb_pix_max))
			self.nb_pix_lat = int(np.sqrt(self.nb_pix_max))
			self.nb_pix 	= self.nb_pix_lon * self.nb_pix_lat
			self.lon_resolution = (self.lon_max - self.lon_min) / self.nb_pix_lon
			self.lat_resolution = (self.lat_max - self.lat_min) / self.nb_pix_lat
			self.resolution = min(AngleToDist(self.lon_resolution), AngleToDist(self.lat_resolution))

		if mpi_rank == 0:
			print("ALT MAP RESOLUTION lon, lat", self.resolution, AngleToDist(self.lon_resolution), AngleToDist(self.lat_resolution))

		self.lon_min -= (self.lon_min % self.lon_resolution)
		self.lon_max -= (self.lon_max % self.lon_resolution)
		# self.lon_max += self.lon_resolution - self.lon_max % self.lon_resolution

		self.lat_min -= (self.lat_min % self.lat_resolution)
		self.lat_max -= (self.lat_max % self.lat_resolution)
		# self.lat_max += self.lat_resolution - self.lat_max % self.lat_resolution


		self.map = np.zeros((self.nb_pix_lat, self.nb_pix_lon))

		### Computes longitudes and latitudes of pixel borders and centers
		self.longitudes = np.linspace(self.lon_min, self.lon_max, self.nb_pix_lon + 1)
		self.latitudes  = np.linspace(self.lat_max, self.lat_min, self.nb_pix_lat + 1)

		self.mid_longitudes = np.array([(self.longitudes[i] + self.longitudes[i+1]) / 2 for i in range(self.nb_pix_lon)])
		self.mid_latitudes = np.array([(self.latitudes[i] + self.latitudes[i+1]) / 2 for i in range(self.nb_pix_lat)])

		### lon and lat values as int to get the map tiles. These are the coordinates of the upper left corner of the tile (north west corner)
		self.lon_name_list = range(int(np.floor(self.lon_min * RtoD)), int(np.floor(self.lon_max * RtoD)) + 1,  1)
		self.lat_name_list = range(int(np.floor(self.lat_max * RtoD)), int(np.floor(self.lat_min * RtoD)) - 1, -1)

		self.nb_tile_lon = len(self.lon_name_list)
		self.nb_tile_lat = len(self.lat_name_list)

		self.nb_tiles = self.nb_tile_lon * self.nb_tile_lat

		self.tiles_coord = []

		for n_lon in self.lon_name_list:
			for n_lat in self.lat_name_list:
				self.tiles_coord.append((n_lon, n_lat))

	def IsInMap(self, lon, lat):
		# print(self.lon_min ,lon ,self.lon_max ,self.lat_min ,lat ,self.lat_max)
		if self.lon_min < lon < self.lon_max and self.lat_min < lat < self.lat_max:
			return True
		else:
			return False

	def GetAltitudeFromLonLat(self, lon, lat, unit=1000):
		"""Altitude of given coord above sea level.
		return result(in meter)/unit. (unit=1000:result in km)"""
		if self.IsInMap(lon, lat):
			# print("IS in map")
			l, c = self.GetIndexFromLonLat(lon, lat)
			return self.map[l, c] / unit
		else:
			# print("NOT in map")
			return 0

	def GetRelativeAltitude(self, lon, lat, unit=1000): #Elevation map is in meters
		"""Altitude of given coord above the instrument in A (center of the map).
		return result(in meter)/unit. (unit=1000:result in km)."""
		return self.GetAltitudeFromLonLat(lon, lat, unit) - self.A_alt / unit

	def GetAltitudeAboveGround(self, lon, lat, alt, unit=1000):
		"""For a given altitude, returns the relative altitude to the given coord altitude.
		return result(in meter)/unit. (unit=1000:result in km)"""
		return alt - self.GetAltitudeFromLonLat(lon, lat, unit)

	def IsAboveGround(self, lon, lat, alt, unit=1000):
		# print("IsAboveGround", self.GetAltitudeAboveGround(lon, lat, alt, unit))
		return self.GetAltitudeAboveGround(lon, lat, alt, unit) >= 0

	def IsVisible(self, lon1, lat1, alt1=None, lon2=None, lat2=None, alt2=None, dlos=None, ret_coll=False):
		"""Return True if two points are visible  (ie if point 1 is visible from point 2, ie if a straight line between them is always above ground.).Else returns False.
		dlos is the distance between each points along the line. (too big and it is long, too low and you will miss collisions.)
		If alt1 is None, point 1 on ground.
		If point 2 not given, use instrument coord and alt.
		If point 2 coord are given but alt2 is None, use point 2 on ground."""

		if not lon2:
			lon2 = self.A_lon
			lat2 = self.A_lat
			alt2 = self.A_alt
		if not alt2:
			alt2 = self.GetAltitudeFromLonLat(lon2, lat2)
		if not alt1:
			alt1 = self.GetAltitudeFromLonLat(lon1, lat1)

		if dlos is None:
			dlos = self.resolution

		is_visible = True

		### Check if end points are above ground
		if not self.IsAboveGround(lon1, lat1, alt1):
			# print("Point 1 is under ground")
			if ret_coll:
				return False, lon1, lat1, alt1
			else:
				return False
		if not self.IsAboveGround(lon2, lat2, alt2):
			# print("Point 2 is under ground")
			if ret_coll:
				return False, lon2, lat2, alt2
			else:
				return False

		if alt1 < alt2:
			obs = ObservationToPoint(lon1, lat1, lon2, lat2, alt2, A_alt = alt1, init_full = False)
		else:
			obs = ObservationToPoint(lon2, lat2, lon1, lat1, alt1, A_alt = alt2, init_full = False)

		### Check every dlos km along path if pass below ground
		# old_alt = 0
		for r in np.arange(dlos/2, obs.AH_norm, dlos): #for every range
			lon, lat, alt = obs.GetPCoordinatesFromRange(r)
			# print(lon, lat, alt, self.GetAltitudeFromLonLat(lon, lat), alt - self.GetAltitudeFromLonLat(lon, lat))
			if alt > self.max_alt:
				break
			if alt < self.GetAltitudeFromLonLat(lon, lat):
				is_visible = False
				break

		### Take care of ret_coll argument
		if ret_coll:
			if not is_visible:
				return is_visible, lon, lat, alt
			else:
				return is_visible, None, None, None
		else:
			return is_visible


	def GetIndexFromLonLat(self, lon, lat):
		line = int((self.lat_max - lat) / self.lat_resolution)
		col  = int((lon - self.lon_min) / self.lon_resolution)
		return line, col


	def GetLonLatFromindex(self, l, c):
		lon = self.longitudes[c]
		lat = self.latitude[l]

		return lon, lat


	def LoadSingleTile(self, n_lon, n_lat):
		"""Load a single tile and returns it, as well as its coordinates. Input parameters are int for lon and lat of the tile."""

		file_name = self.GetNameFromLonLat(n_lon, n_lat)
		# if mpi_rank==0: print(file_name)

		try:
			tile = gdal.Open(self.path + file_name, gdal.GA_ReadOnly) #get raster objsect from tif file
		except:
			print("WARNING: Elevation tile not found:", self.path + file_name)
			return np.zeros((1,1)), self.A_lon, self.A_lat

		# Shape of the map, width and height
		nb_lon = tile.RasterXSize
		nb_lat = tile.RasterYSize

		### Get info on the tile:
		tile_info = tile.GetGeoTransform()
		# if mpi_rank==0: print("tile info", tile_info)
		origin_lon = tile_info[0]	 * DtoR #western border
		origin_lat = tile_info[3]	 * DtoR #northern border
		pixel_width = tile_info[1]	 * DtoR #pixel size in lon > 0
		pixel_height = tile_info[5]	 * DtoR #pixel size in lat < 0 because origin is lat_max !!!

		# if mpi_rank==0: print(origin_lon*RtoD, origin_lat*RtoD, pixel_width*RtoD, pixel_height*RtoD)

		### Function to pass from lon, lat to row and column index
		LonToCol = lambda l: int((l - origin_lon) / pixel_width)
		LatToRow = lambda l: int((l - origin_lat) / pixel_height) #both are < 0, so row is OK

		# Find min and max values for lon and lat depending on the instrument position and radius of the map.
		# if mpi_rank==0: print("test", self.A_lon*RtoD - DistToAngle(self.radius, lat=self.A_lat, az = 90*DtoR)*RtoD, origin_lon*RtoD, self.A_lon*RtoD + DistToAngle(self.radius, lat=self.A_lat, az = 90*DtoR)*RtoD, origin_lon*RtoD + pixel_width * nb_lon*RtoD)
		# if mpi_rank==0: print("Alon, Alat", self.A_lon*RtoD, self.A_lat*RtoD, DistToAngle(self.radius, lat=self.A_lat, az = 90*DtoR)*RtoD, DistToAngle(self.radius)*RtoD)
		lon_min = max(self.A_lon - DistToAngle(self.radius, lat=self.A_lat, az = 90*DtoR), origin_lon)
		lon_max = min(self.A_lon + DistToAngle(self.radius, lat=self.A_lat, az = 90*DtoR), origin_lon + pixel_width * nb_lon)

		lat_min = max(self.A_lat - DistToAngle(self.radius), origin_lat + pixel_height * nb_lat)
		lat_max = min(self.A_lat + DistToAngle(self.radius), origin_lat)

		nb_lon  = int((lon_max - lon_min) / pixel_width)
		nb_lat  = int((lat_min - lat_max) / pixel_height) #inversed because pixel_height < 0

		# if mpi_rank==0:
		# 	print("LON", lon_min*RtoD, lon_max*RtoD)
		# 	print("LAT", lat_min*RtoD, lat_max*RtoD)
		# 	print(LonToCol(lon_min), LatToRow(lat_max), nb_lon, nb_lat)

		tile = tile.GetRasterBand(1).ReadAsArray(LonToCol(lon_min), LatToRow(lat_max), nb_lon, nb_lat)

		if pixel_width != self.lon_resolution:
			nb_lon  = int((lon_max - lon_min) / self.lon_resolution)
			nb_lat  = int((lat_max - lat_min) / self.lat_resolution)

			tile = cv2.resize(tile, dsize=(nb_lon, nb_lat), interpolation = cv2.INTER_CUBIC)

		return tile, lon_min, lat_max

	def LoadMultipleTiles(self):
		"""For a large map made of several tiles, load all accessible tiles and merge them into one map"""

		### Loop through all tiles via there coord in lon, lat (as int)
		for tile_lon, tile_lat in self.tiles_coord:
			tile, lon_min, lat_max = self.LoadSingleTile(tile_lon, tile_lat) #Get tile and position of the upper left corner, lon_min (left of the map, west), lat_max (top of the map, north)

			tile_line, tile_col = self.GetIndexFromLonLat(lon_min, lat_max) # map index of the tile upper left corner
			nb_lat, nb_lon = tile.shape

			# print("tile.shape", tile.shape)
			# print("map chunk shape", self.map[tile_line : tile_line + nb_lat, tile_col : tile_col + nb_lon].shape)
			#
			# print(self.map.shape, tile_line, tile_line + nb_lat, tile_col, tile_col + nb_lon)

			### Merge tile into the final map
			self.map[tile_line : tile_line + nb_lat, tile_col : tile_col + nb_lon] = tile


	def PlotMap(self, show=False, iso = None):
		f = plt.figure()
		plt.pcolormesh((self.longitudes - self.A_lon) * RT, (self.latitudes - self.A_lat) * RT, self.map)
		plt.colorbar()

		if iso is not None:
			plt.contour((self.mid_longitudes - self.A_lon) * RT, (self.mid_latitudes - self.A_lat) * RT, self.map, iso)

		# plt.pcolormesh(np.flip(self.map, axis=0))
		# plt.pcolormesh(np.flip(np.rot90(self.map, k=-1)))
		if show:
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
		name += str(int(abs(lat))).zfill(3)

		if lon >= 0:
			name += "E"
		else:
			name += "W"
		name += str(int(abs(lon))).zfill(3)

		return name

class GroundMap:
	def __init__(self, in_dict, Nb_a_pc, Nb_e_pc):

		self.show_ground_albedo = False
		if "show_ground_albedo" in in_dict:
			self.show_ground_albedo = bool(int(in_dict["show_ground_albedo"]))

		self.path = pomerol_configuration.ground_map_path
		self.file = in_dict["ground_file"]

		self.location = in_dict["location"]
		self.A_lon, self.A_lat 	= GetLonLatFromName(self.location)

		self.mode = in_dict["ground_mode"].lower()
		self.shader = in_dict["shader_mode"].lower()

		self.radius = float(in_dict["ground_emission_radius"]) #in km
		self.N_bins_max = int(in_dict["ground_N_bins_max"])

		# self.exist = ((self.radius > 0 and self.file) or float(in_dict["point_src_I0"]) > 0)
		# self.is_point_source = float(in_dict["point_src_I0"]) > 0

		self.exist = self.mode != "none"
		self.is_point_source = "point" in self.mode

		self.src_I0, self.src_az, self.src_dist = None, None, None

		self.cube = None
		self.Nt = 1

		if self.exist:
			self.LoadMap(Nb_a_pc, Nb_e_pc)

			# if mpi_rank == 0:
				# self.MakeInitialMapPlot()
			# 	# self.MakePlots(0, 0)
				# plt.show()

	# def DistToAngle(self, dist, lat = 0, az = 0, R = RT):
	# 	return (dist  / R) * np.sqrt(np.cos(az)**2 + np.sin(az)**2 / np.cos(lat)**2)
	#
	# def AngleToDist(self, angle, lat = 0, az = 0, R = RT):
	# 	return (angle * R) / np.sqrt(np.cos(az)**2 + np.sin(az)**2 / np.cos(lat)**2)


	def Generator(self, it, Nlos):
		"""Generator function used for the compute shader."""
		for ie in range(self.Ndist):
			for ia in range(self.Naz):
				for ilos in range(self.Nlos):
					yield self.cube[it, ie, ia] # yield the emitted flux
					yield 0. # yield the emitted DoLP (always zero for now)
					yield 0. # yield the emitted DoLP (always zero for now)


	def LoadMap(self, Nb_a_pc, Nb_e_pc):

		self.Naz, self.Ndist = int(np.sqrt(self.N_bins_max)), int(np.sqrt(self.N_bins_max))

		self.azimuts,   self.daz   = np.linspace(0, 2 * np.pi, self.Naz + 1, retstep=True)
		### ADAPTATIVE DISTANCES (square)
		self.distances = np.array([x**2 for x in np.linspace(0, np.sqrt(self.radius), self.Ndist + 1)]) #in km
		### CONSTANT DISTANCES (linear)
		# self.distances, self.ddist = np.linspace(0, self.radius, self.Ndist + 1, retstep=True)  #in km

		self.mid_distances	= np.array([(self.distances[i+1] + self.distances[i]) / 2. for i in range(0, self.Ndist)]) #in km
		self.mid_azimuts	= self.azimuts[:-1] + self.daz / 2.

		# print(self.mid_distances*RT)

		self.maps_shape = (self.Nt, Nb_e_pc, Nb_a_pc, self.Ndist, self.Naz)

		if "image" in self.mode:
			self.LoadGroundImageMap()
		elif "point" in self.mode:
			s = self.mode.split("_")
			I0 = float(s[1][1:])
			az = float(s[2][1:]) * DtoR
			dist = float(s[3][1:])
			# print(f"{I0}, {az}, {dist}")
			self.LoadPointSourceMap(I0, az, dist, Nb_a_pc, Nb_e_pc)
			self.is_point_source = True
		elif "uniform" in self.mode:
			I0 = float(self.mode.split("_")[1][1:])
			self.LoadUniformMap(I0)
		elif "gaussian" in self.mode:
			s = self.mode.split("_")
			I0 = float(s[1][1:])
			width = float(s[2][1:])
			az = float(s[3][1:]) * DtoR
			dist = float(s[4][1:])
			self.LoadGaussianMap(I0, width, az, dist, Nb_a_pc, Nb_e_pc)
		elif "none" in self.mode:
			self.distances = self.azimuts = []
			self.ddist = self.daz = 0
			self.exist = False

		else:
			print("WARNING: Ground map mode is not correct. Ground map not used.")
			self.exist = False


		# self.I_map[self.GetPixFromAzDist(0*DtoR, 0.1)] +=

		if self.exist:
			self.Ndist 					= len(self.distances) - 1
			self.Naz 					= len(self.azimuts) - 1

			self.mid_distances			= np.array([(self.distances[i+1] + self.distances[i]) / 2. for i in range(0, self.Ndist)])
			self.mid_azimuts			= self.azimuts[:-1] + self.daz / 2.

			self.ShadeMap()

			self.cube 					= np.zeros((self.Nt, self.Ndist, self.Naz))
			self.cube[0] 				+= self.I_map

			self.scattering_map 		= np.zeros(self.maps_shape) # Intensity from (e, a) reaching us
			self.DoLP_map				= np.zeros(self.maps_shape) # Polarized intensity from (e, a) reaching us
			self.total_scattering_map 	= np.zeros(self.maps_shape) # DoLP of scattered light from (e,a)
			self.AoRD_map 				= np.zeros(self.maps_shape) # Angle of polarisation of light from (e,a)
			self.V_map 					= np.zeros(self.maps_shape)
			self.Vcos_map 				= np.zeros(self.maps_shape)
			self.Vsin_map 				= np.zeros(self.maps_shape)

			self.mountain_shadow_heights = np.zeros(self.maps_shape)

			self.V_total	= np.zeros((self.Nt, Nb_e_pc, Nb_a_pc))
			self.Vcos_total	= np.zeros((self.Nt, Nb_e_pc, Nb_a_pc))
			self.Vsin_total	= np.zeros((self.Nt, Nb_e_pc, Nb_a_pc))

			self.I0_total	= np.zeros((self.Nt, Nb_e_pc, Nb_a_pc))
			self.DoLP_total	= np.zeros((self.Nt, Nb_e_pc, Nb_a_pc))
			self.AoLP_total	= np.zeros((self.Nt, Nb_e_pc, Nb_a_pc))


	def ShadeMap(self, **kwargs):
		if "half" in self.shader:
			mid_a = (float(self.shader.split("_")[1][1:]) * DtoR)%(2*np.pi)
			a1, a2 = (mid_a - np.pi/2)%(np.pi*2), (mid_a + np.pi/2)%(np.pi*2)
			# print("DEBUG condition", mid_a*RtoD, a1*RtoD, a2*RtoD, a1 >= a2)
			def condition(d, a, **kwargs):
				ret = 1
				if a > a1 or a < a2 if a1 > a2 else a1 < a < a2:
					ret = 0
				return ret

			# condition = lambda d, a: min(a1, a2) <= a < max(a1, a2)

		elif "random" in self.shader:
			percentage = float(self.shader.split("_")[1][:]) / 100.
			# condition = lambda d, a, id=0: np.random.random() < percentage
			def condition(d, a, **kwargs):
				ret = 1
				if np.random.random() < percentage:
					ret = 0
				return ret

		elif "in_dist" in self.shader:
			min_dist = float(self.shader.split("_")[2][1:])
			max_dist = float(self.shader.split("_")[3][1:])
			def condition(d, a, **kwargs):
				ret = 1
				if min_dist < d < max_dist:
					ret = 0
				return ret

		elif "out_dist" in self.shader:
			min_dist = float(self.shader.split("_")[2][1:])
			max_dist = float(self.shader.split("_")[3][1:])
			def condition(d, a, **kwargs):
				ret = 1
				if d < min_dist  or d > max_dist:
					ret = 0
				return ret

		elif "ring" in self.shader:
			index = int(self.shader.split("_")[1])
			def condition(d, a, **kwargs):
				ret = 1
				if id != index:
					ret = 0
				return ret
		# elif "altitude" in self.shader and "alt_map" in kwargs:
		# 	alt_map = kwargs["alt_map"]
		# 	def condition(d, a, **kwargs):
		# 		ret = 1
		# 		lon, lat = AzDistToLonLat(a, d, origin_lon = self.A_lon, origin_lat = self.A_lat)
		# 		print("COORD", a*RtoD, d, lon*RtoD, lat*RtoD)
		# 		alt = alt_map.GetAltitudeFromLonLat(lon, lat)
		# 		print(alt)
		# 		if alt < 1:
		# 			ret = 0
		# 		else:
		# 			ret = 1
		# 		return ret
		else:
			def condition(d, a, **kwargs):
				return 1

		for id, d in enumerate(self.mid_distances):
			for ia, a in enumerate(self.mid_azimuts):
				self.I_map[id, ia] *= condition(d, a, id = id, ia = ia)



	def LoadUniformMap(self, I0):
		self.I_map = I0 * np.ones((self.Ndist, self.Naz)) # [nW / m2/ sr]


	def LoadGroundImageMap(self):
		"""In the case where h==0, load the emission map from the geotiff files"""

		if mpi_rank == 0:
			print("Loading ground emission map data...")

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

		# print(self.radius, self.A_lat,  np.pi/2)
		lat_radius = DistToAngle(self.radius, lat = self.A_lat, az = 0)
		lon_radius = DistToAngle(self.radius, lat = self.A_lat, az = np.pi/2)

		# print("DEBUG RADIUS SKI", AngleToDist(lat_radius, lat = self.A_lat, az = 0), AngleToDist(lon_radius, lat = self.A_lat, az = 90*DtoR))

		lon_min = max(self.A_lon - lon_radius, origin_lon)
		lon_max = min(self.A_lon + lon_radius, origin_lon + pixel_width * nb_lon)

		# print("DEBUG lon min", self.A_lon - lon_radius, origin_lon)

		lat_min = max(self.A_lat - lat_radius, origin_lat + pixel_height * nb_lat)
		lat_max = min(self.A_lat + lat_radius, origin_lat)

		nb_lon  = int( 2 * lon_radius / pixel_width)
		nb_lat  = int(-2 * lat_radius / pixel_height)

		# print("GROUND MAP", nb_lon, nb_lat)

		longitudes, dlon = np.linspace(lon_min, lon_max, nb_lon, retstep=True) # list of pixel longitudes
		latitudes,  dlat = np.linspace(lat_max, lat_min, nb_lat, retstep=True) # list of pixel latitudes

		map_band = map.GetRasterBand(1)
		I_map = map_band.ReadAsArray(LonToCol(lon_min), LatToRow(lat_max), nb_lon, nb_lat) # Emission map we will use in #[nW / cm2/ sr]
		I_map *= 1e4 #in [nW / m2 / sr]

		# plt.pcolormesh(AngleToDist(longitudes-self.A_lon, lat = self.A_lat, az = np.pi/2), AngleToDist(latitudes-self.A_lat, lat = self.A_lat, az = 0), I_map)
		# plt.show()

		self.I_map = self.CarthToPolar(I_map, longitudes, latitudes)


		### The following is to use when the polar map has better resolution than the carthsian one. Some pixel may then be equal to zero because no carthesian pixel center fall in them.
		### It does not work for now !!!!
		# for id, d in enumerate(self.mid_distances):
		# 	for ia, a in enumerate(self.mid_azimuts):
		# 		if self.I_map[id, ia] == 0:
		# 			pol_lon, pol_lat = AzDistToLonLat(a, d, origin_lon = self.A_lon, origin_lat = self.A_lat)
		# 			pol_lon_index = (pol_lon - lon_min) / pixel_width
		# 			pol_lat_index = (pol_lat - lat_max) / pixel_height
		# 			self.I_map[id, ia] = sci.ndimage.map_coordinates(I_map, [[pol_lon_index], [pol_lat_index]])

	def CarthToPolar(self, carth_map, longitudes, latitudes):
		polar_map = np.zeros((self.Ndist, self.Naz))
		divisor_map = np.zeros((self.Ndist, self.Naz))

		for ilon, lon in enumerate(longitudes):
			for ilat, lat in enumerate(latitudes):
				az, dist = LonLatToAzDist(lon, lat, self.A_lon, self.A_lat)
				dist *= RT
				iaz, idist = self.GetPixFromAzDist(az, dist)
				if idist is not None and iaz is not None:
					# print(dist, idist, az*RtoD, iaz)
					polar_map[idist, iaz] += carth_map[ilat, ilon]
					divisor_map[idist, iaz] += 1

		polar_map = np.divide(polar_map, divisor_map, out = np.zeros_like(polar_map), where = divisor_map != 0)

		for idist in range(self.Ndist):
			for iaz in range(self.Naz):
				if divisor_map[idist, iaz] == 0:
					lon, lat = AzDistToLonLat(self.mid_azimuts[iaz], self.mid_distances[idist], origin_lon = self.A_lon, origin_lat = self.A_lat)
					ilon, ilat = min(np.searchsorted(longitudes, lon), len(longitudes)-1), min(np.searchsorted(latitudes[-1::-1], lat), len(latitudes)-1)

					# print(idist, iaz, self.mid_azimuts[iaz], self.mid_distances[idist], ilon, ilat, len(longitudes), len(latitudes), carth_map[ilat, ilon])

					polar_map[idist, iaz] = carth_map[ilat, ilon]
					# pp = polar_map[min(self.Ndist, 	idist + 1), (iaz + 1) % self.Naz]
					# pm = polar_map[min(self.Ndist, 	idist + 1), (iaz - 1) % self.Naz]
					# mp = polar_map[max(0, 			idist - 1), (iaz + 1) % self.Naz]
					# mm = polar_map[max(0, 			idist - 1), (iaz - 1) % self.Naz]
					# polar_map[idist, iaz] = (pp + pm + mp + mm) / 4.

		return polar_map

	def LoadPointSourceMap(self, src_I0, src_az, src_dist, Nb_a_pc, Nb_e_pc):

		self.src_I0, self.src_az, self.src_dist = src_I0, src_az, src_dist

		src_lon, src_lat = AzDistToLonLat(src_az, src_dist, self.A_lon, self.A_lat)

		self.distances, self.ddist = np.linspace(src_dist - 0.5, src_dist + 0.5, 2, retstep=True) # list of pixel longitudes
		self.azimuts, self.daz 	   = np.linspace(src_az - 0.5, src_az + 0.5, 2, retstep=True) # list of pixel latitudes

		# print(self.longitudes)
		self.I_map = np.ones((1, 1)) * src_I0 #[nW / m2 / sr]

		self.Naz = self.Ndist = 1

		self.maps_shape = (self.Nt, Nb_e_pc, Nb_a_pc, 1, 1)
		self.is_point_source = True

	def GetPixFromAz(self, az):
		# If the azimuth is in the map range
		az = az%(2*np.pi)
		if az < self.azimuts[0] or az > self.azimuts[-1]:
			return None

		iaz = int((az % (2*np.pi)) / self.daz) % self.Naz
		# print(az*RtoD, iaz, self.azimuts*RtoD)
		return iaz

	def GetPixFromDist(self, dist):
		# print(self.Ndist, self.ddist * RT, dist * RT, int(dist / self.ddist))
		# idist = np.around(dist / self.ddist) #ONLY IF DISTANCE BINS ARE CONSTANT!!!

		# If the distance is in the map range
		if dist < self.distances[0] or dist > self.distances[-1]:
			return None

		idist = np.searchsorted(self.distances, dist, side="right") - 1
		# print(dist, idist, self.distances)

		return int(idist)

	def GetPixFromAzDist(self, az, dist):
		return self.GetPixFromAz(az), self.GetPixFromDist(dist)


	def GetRadiantFluxFromAzDist(self, az, dist, t=0):
		""" Return radiant flux (nW/m2/sr) for a given azimut (rad) and distance (km). Return 0 if outside of map.
		"""
		iaz, idist = self.GetPixFromAzDist(az, dist)
		if iaz is None or idist is None:
			return 0
		else:
			return self.cube[t, idist, iaz]

	def GetRadiantIntensity(self, az, dist, t=0):
		""" Return radiant intensity (nW/sr) for a given azimut (rad) and distance (km). Return 0 if outside of map.
		"""
		iaz, idist = self.GetPixFromAzDist(az, dist)
		if iaz is None or idist is None:
			return 0
		else:
			return self.cube[t, idist, iaz] * self.GetArea(idist)


	#@timer
	def GetArea(self, idist):
		"""Return the area of a pixel on the map in m**2. If we use a point source, the area is set to one."""
		if not self.is_point_source:
			area = 0.5 * abs(self.daz)
			area *= (self.distances[idist+1]**2 - self.distances[idist]**2)
			area *= 1e6
			return area
		else:
			return 1

	def SetLightParameters(self):
		self.total_scattering_map 	= 2 * self.V_map

		self.DoLP_map				= 2 * np.sqrt(self.Vcos_map**2 + self.Vsin_map**2) / 100 # DoLP of scattered light from (e,a)
		self.DoLP_map				= np.divide(self.DoLP_map, self.V_map, out = np.zeros_like(self.DoLP_map), where = self.V_map != 0)

		self.AoRD_map 				= np.arctan2(self.Vsin_map, self.Vcos_map) / 2. # Angle of polaisation of light from (e,a)
		self.scattering_map 		= self.total_scattering_map * self.DoLP_map / 100. # Polarized intensity from (e, a) reaching us

		for it, t in enumerate(range(self.Nt)):
			for ie, e in enumerate(self.V_map[0]):
				for ia, a in enumerate(e):
					self.V_total[it, ie, ia] = np.sum(self.V_map[it, ie, ia, :, :])
					self.Vcos_total[it, ie, ia] = np.sum(self.Vcos_map[it, ie, ia, :, :])
					self.Vsin_total[it, ie, ia] = np.sum(self.Vsin_map[it, ie, ia, :, :])

					self.I0_total[it, ie, ia] = 2 * self.V_total[it, ie, ia]
					self.DoLP_total[it, ie, ia] = 100 * 2 * np.sqrt(self.Vcos_total[it, ie, ia] ** 2 + self.Vsin_total[it, ie, ia] ** 2) / self.V_total[it, ie, ia] #in %
					self.AoLP_total[it, ie, ia] = np.arctan2(self.Vsin_total[it, ie, ia], self.Vcos_total[it, ie, ia]) / 2


	def ProlongateTime(self, new_Nt):
		self.Nt = new_Nt
		self.cube = np.zeros((self.Nt, self.Ndist, self.Naz))

		Nb_e_pc, Nb_a_pc = self.scattering_map.shape[1], self.V_total.shape[2]
		self.maps_shape = (self.Nt, Nb_e_pc, Nb_a_pc, self.Ndist, self.Naz)

		self.scattering_map 		= np.zeros(self.maps_shape) # Intensity from (e, a) reaching us
		self.DoLP_map				= np.zeros(self.maps_shape) # Polarized intensity from (e, a) reaching us
		self.total_scattering_map 	= np.zeros(self.maps_shape) # DoLP of scattered light from (e,a)
		self.AoRD_map 				= np.zeros(self.maps_shape) # Angle of polaisation of light from (e,a)
		self.V_map 					= np.zeros(self.maps_shape)
		self.Vcos_map 				= np.zeros(self.maps_shape)
		self.Vsin_map 				= np.zeros(self.maps_shape)

		self.mountain_shadow_heights = np.zeros(self.maps_shape)

		self.V_total	= np.zeros((self.Nt, Nb_e_pc, Nb_a_pc))
		self.Vcos_total	= np.zeros((self.Nt, Nb_e_pc, Nb_a_pc))
		self.Vsin_total	= np.zeros((self.Nt, Nb_e_pc, Nb_a_pc))

		self.I0_total	= np.zeros((self.Nt, Nb_e_pc, Nb_a_pc))
		self.DoLP_total	= np.zeros((self.Nt, Nb_e_pc, Nb_a_pc))
		self.AoLP_total	= np.zeros((self.Nt, Nb_e_pc, Nb_a_pc))


	def MakeInitialMapPlot(self, e_pc=None, a_pc=None, atmosphere_max = None, save = False, metadata=None):
		f = plt.figure(figsize=(16, 8))
		ax = plt.subplot(111, projection='polar')

		i1 = ax.pcolormesh(self.azimuts, self.distances, (self.cube[0]))

		cbar1 = f.colorbar(i1, extend='both', spacing='proportional', shrink=0.9, ax=ax)
		cbar1.set_label('Initial emission map')

		ax.set_theta_zero_location("N")
		ax.set_theta_direction(-1)

		if e_pc is not None and a_pc is not None:
			if atmosphere_max:
				length = min(atmosphere_max / np.tan(e_pc), self.radius * RT)
			else:
				length = self.radius * RT

			# a.add_artist(Arrow(0, 0, a_pc, length, color="red", width = 0.3))
			a.plot([0, a_pc], [0, length], "r")

		if save:
			plt.savefig(save + '_ground_emission.png', bbox_inches='tight', metadata=metadata)

		return f, ax

	def MakePlots(self, ie_pc, ia_pc, e_pc=None, a_pc=None, atmosphere_max = None, save = False, metadata=None):
		f, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex = True, sharey = True, figsize=(16, 8))
		ax1 = plt.subplot(221, projection='polar')
		ax2 = plt.subplot(222, projection='polar')
		ax3 = plt.subplot(223, projection='polar')
		ax4 = plt.subplot(224, projection='polar')

		i1 = ax1.pcolormesh(self.azimuts, self.distances, (self.cube[0]))
		i2 = ax2.pcolormesh(self.azimuts, self.distances, (self.total_scattering_map[0, ie_pc, ia_pc, :, :]))		# Intensity from (e, a) reaching us
		i3 = ax3.pcolormesh(self.azimuts, self.distances, (self.scattering_map[0, ie_pc, ia_pc, :, :]))				# Polarized intensity from (e, a) reaching us
		i4 = ax4.pcolormesh(self.azimuts, self.distances, self.mountain_shadow_heights[0, ie_pc, ia_pc, :, :])	# DoLP of scattered light from (e,a)
		# i4 = ax4.pcolormesh(self.azimuts, self.distances, self.DoLP_map[0, ie_pc, ia_pc, :, :])	# DoLP of scattered light from (e,a)
		# i4 = ax4.pcolormesh((azimuts-A_lon) * RT, (distances-A_lat) * RT, AoRD_map * RtoD, cmap=cm.twilight)			# AoLP of scattered light from (e,a)


		# Angle of polarisation of light from (e,a)
		cbar1 = f.colorbar(i1, extend='both', spacing='proportional', shrink=0.9, ax=ax1)
		cbar1.set_label('Initial emission map')
		cbar2 = f.colorbar(i2, extend='both', spacing='proportional', shrink=0.9, ax=ax2)
		cbar2.set_label('Total intensity map')
		cbar3 = f.colorbar(i3, extend='both', spacing='proportional', shrink=0.9, ax=ax3)
		cbar3.set_label('Polarised intensity map')
		cbar4 = f.colorbar(i4, extend='both', spacing='proportional', shrink=0.9, ax=ax4)
		cbar4.set_label('1st alt on LoS (km)')

		# f2.suptitle("Relevant angles map at skibotn, with light pollution source at -45deg in azimut")

		for a in [ax1, ax2, ax3, ax4]:
			a.set_theta_zero_location("N")
			a.set_theta_direction(-1)

			# if e_pc is not None and a_pc is not None:
			# 	if atmosphere_max:
			# 		length = min(atmosphere_max / np.tan(e_pc), self.radius * RT)
			# 	else:
			# 		length = self.radius * RT
			#
			# 	# a.add_artist(Arrow(0, 0, a_pc, length, color="red", width = 0.3))
			# 	a.plot([0, a_pc], [0, length], "r")

		if save:
			plt.savefig(save + '_groundmaps.png', bbox_inches='tight', metadata=metadata)


# class OldGroundMap:
# 	def __init__(self, in_dict, Nb_a_pc, Nb_e_pc):
#
# 		self.path = in_dict["ground_path"]
# 		self.file = in_dict["ground_file"]
#
# 		self.location = in_dict["location"]
# 		self.A_lon, self.A_lat 	= GetLonLatFromName(self.location)
#
# 		self.mode = in_dict["ground_mode"].lower()
#
# 		self.radius = float(in_dict["ground_emission_radius"]) / RT #in radians
# 		self.N_bins_max = int(in_dict["ground_N_bins_max"])
#
# 		# self.exist = ((self.radius > 0 and self.file) or float(in_dict["point_src_I0"]) > 0)
# 		# self.is_point_source = float(in_dict["point_src_I0"]) > 0
#
# 		self.exist = self.mode != "none"
# 		self.is_point_source = "point" in self.mode
#
# 		self.src_I0, self.src_az, self.src_dist = None, None, None
#
# 		self.cube = None
# 		self.Nt = 1
#
# 		if self.exist:
# 			self.LoadMap(Nb_a_pc, Nb_e_pc)
# 			# if self.is_point_source:
# 			# 	I0 = float(in_dict["point_src_I0"])
# 			# 	az = float(in_dict["point_src_az"]) * DtoR
# 			# 	dist = float(in_dict["point_src_dist"])
# 			# 	self.LoadPointSourceMap(I0, az, dist, Nb_a_pc, Nb_e_pc)
# 			# elif self.is_uniform:
# 			# 	self.LoadUniformMap(Nb_a_pc, Nb_e_pc)
# 			# else:
# 			# 	self.LoadGroundEmmisionsMap(Nb_a_pc, Nb_e_pc)
#
# 	def LoadMap(self, Nb_a_pc, Nb_e_pc):
# 		if "image" in self.mode:
# 			self.LoadGroundImageMap(Nb_a_pc, Nb_e_pc)
# 		elif "point" in self.mode:
# 			s = self.mode.split("_")
# 			I0 = float(s[1][1:])
# 			az = float(s[2][1:]) * DtoR
# 			dist = float(s[3][1:])
# 			self.LoadPointSourceMap(I0, az, dist, Nb_a_pc, Nb_e_pc)
# 		elif "uniform" in self.mode:
# 			I0 = float(self.mode.split("_")[1][1:])
# 			self.LoadUniformMap(I0, Nb_a_pc, Nb_e_pc)
# 		elif "gaussian" in self.mode:
# 			s = self.mode.split("_")
# 			I0 = float(s[1][1:])
# 			width = float(s[2][1:])
# 			az = float(s[3][1:]) * DtoR
# 			dist = float(s[4][1:])
# 			self.LoadGaussianMap(I0, width, az, dist, Nb_a_pc, Nb_e_pc)
# 		elif "none" in self.mode:
# 			self.longitudes = self.latitudes = []
# 			self.dlon = self.dlat = 0
# 			self.exist = False
#
# 		else:
# 			print("WARNING: Ground map mode is not correct. Ground map not used.")
# 			self.exist = False
#
#
# 		if self.exist:
# 			self.mid_longitudes			= self.longitudes[:-1] + self.dlon/2.
# 			self.mid_latitudes			= self.latitudes[:-1] + self.dlat/2.
#
# 			self.Nlat 					= len(self.mid_latitudes)
# 			self.Nlon 					= len(self.mid_longitudes)
#
# 			self.cube 					= np.zeros((self.Nt, self.Nlat, self.Nlon))
# 			self.cube[0] 				+= self.I_map
#
# 			self.scattering_map 		= np.zeros(self.maps_shape) # Intensity from (e, a) reaching us
# 			self.DoLP_map				= np.zeros(self.maps_shape) # Polarized intensity from (e, a) reaching us
# 			self.total_scattering_map 	= np.zeros(self.maps_shape) # DoLP of scattered light from (e,a)
# 			self.AoRD_map 				= np.zeros(self.maps_shape) # Angle of polaisation of light from (e,a)
# 			self.V_map 					= np.zeros(self.maps_shape)
# 			self.Vcos_map 				= np.zeros(self.maps_shape)
# 			self.Vsin_map 				= np.zeros(self.maps_shape)
#
# 			self.V_total	= np.zeros((self.Nt, Nb_e_pc, Nb_a_pc))
# 			self.Vcos_total	= np.zeros((self.Nt, Nb_e_pc, Nb_a_pc))
# 			self.Vsin_total	= np.zeros((self.Nt, Nb_e_pc, Nb_a_pc))
#
# 			self.I0_total	= np.zeros((self.Nt, Nb_e_pc, Nb_a_pc))
# 			self.DoLP_total	= np.zeros((self.Nt, Nb_e_pc, Nb_a_pc))
# 			self.AoLP_total	= np.zeros((self.Nt, Nb_e_pc, Nb_a_pc))
#
#
# 	def LoadGaussianMap(self, I0, width, a, d, Nb_a_pc, Nb_e_pc):
# 		nb_lon, nb_lat = int(np.sqrt(self.N_bins_max)/2+1)*2, int(np.sqrt(self.N_bins_max)/2+1)*2
#
# 		lon_min, lon_max = self.A_lon - self.radius / np.cos(self.A_lat), self.A_lon + self.radius / np.cos(self.A_lat)
#
# 		lat_min, lat_max = self.A_lat - self.radius, self.A_lat + self.radius
#
# 		self.longitudes, self.dlon = np.linspace(lon_min, lon_max, nb_lon + 1, retstep=True) # list of pixel longitudes
# 		self.latitudes, self.dlat = np.linspace(lat_max, lat_min, nb_lat + 1, retstep=True) # list of pixel latitudes
#
# 		self.mid_longitudes	= self.longitudes[:-1] + self.dlon/2.
# 		self.mid_latitudes	= self.latitudes[:-1] + self.dlat/2.
#
# 		x, y = np.meshgrid(self.mid_longitudes, self.mid_latitudes)
#
# 		x0 = self.A_lon + d / RT * np.sin(a)
# 		y0 = self.A_lat + d / RT * np.cos(a)
#
# 		d = np.sqrt((x-x0)**2 + (y-y0)**2)
#
# 		self.I_map = I0 * np.exp(-0.5 * (d / width) ** 2) #[nW / m2/ sr]
#
# 		self.maps_shape = (self.Nt, Nb_e_pc, Nb_a_pc, nb_lat, nb_lon)
#
#
# 	def LoadUniformMap(self, I0, Nb_a_pc, Nb_e_pc):
# 		nb_lon, nb_lat = int(np.sqrt(self.N_bins_max)/2+1)*2, int(np.sqrt(self.N_bins_max)/2+1)*2
#
# 		lon_min, lon_max = self.A_lon - self.radius / np.cos(self.A_lat), self.A_lon + self.radius / np.cos(self.A_lat)
# 		lat_min, lat_max = self.A_lat - self.radius, self.A_lat + self.radius
#
# 		self.longitudes, self.dlon = np.linspace(lon_min, lon_max, nb_lon + 1, retstep=True) # list of pixel longitudes
# 		self.latitudes, self.dlat = np.linspace(lat_max, lat_min, nb_lat + 1, retstep=True) # list of pixel latitudes
#
# 		# print(self.longitudes * RT, self.latitudes * RT)
#
# 		self.I_map = I0 * np.ones((nb_lat, nb_lon)) # [nW / m2/ sr]
#
# 		# print(self.I_map)
#
# 		self.maps_shape = (self.Nt, Nb_e_pc, Nb_a_pc, nb_lat, nb_lon)
# 		# print(self.maps_shape)
#
#
# 	def LoadGroundImageMap(self, Nb_a_pc, Nb_e_pc):
# 		"""In the case where h==0, load the emission map from the geotiff files"""
# 		map = gdal.Open(self.path + self.file, gdal.GA_ReadOnly)
# 		nb_lon = map.RasterXSize
# 		nb_lat = map.RasterYSize
#
# 		nb_bands = map.RasterCount
#
# 		map_info = map.GetGeoTransform()
# 		origin_lon = map_info[0]	 * DtoR
# 		origin_lat = map_info[3]	 * DtoR
# 		pixel_width = map_info[1]	 * DtoR
# 		pixel_height = map_info[5]	 * DtoR
#
# 		LonToCol = lambda lon: int((lon - origin_lon) / pixel_width)
# 		LatToRow = lambda lat: int((lat - origin_lat) / pixel_height)
#
# 		lon_min = max(self.A_lon - self.radius / np.cos(self.A_lat), origin_lon)
# 		lon_max = min(self.A_lon + self.radius / np.cos(self.A_lat), origin_lon + pixel_width * nb_lon)
#
# 		lat_min = max(self.A_lat - self.radius, origin_lat + pixel_height * nb_lat)
# 		lat_max = min(self.A_lat + self.radius, origin_lat)
#
# 		nb_lon  = int( 2 * self.radius / pixel_width)
# 		nb_lat  = int(-2 * self.radius / pixel_height)
#
# 		# print("GROUND MAP", nb_lon, nb_lat)
#
# 		self.longitudes, self.dlon = np.linspace(lon_min, lon_max, nb_lon, retstep=True) # list of pixel longitudes
# 		self.latitudes, self.dlat = np.linspace(lat_max, lat_min, nb_lat, retstep=True) # list of pixel latitudes
#
# 		map_band = map.GetRasterBand(1)
# 		self.I_map = map_band.ReadAsArray(LonToCol(lon_min), LatToRow(lat_max), nb_lon, nb_lat) # Emission map we will use in #[nW / m2/ sr]
# 		self.I_map *= 1e4
#
# 		if nb_lon * nb_lat > self.N_bins_max and self.N_bins_max > 0:
# 			# print("GROUND MAP TOO BIG", nb_lon * nb_lat, self.N_bins_max)
# 			new_nb_lon, new_nb_lat = int(np.sqrt(self.N_bins_max)), int(np.sqrt(self.N_bins_max))
# 			# print("GROUND MAP", new_nb_lon, new_nb_lat)
#
# 			self.longitudes, self.dlon = np.linspace(lon_min, lon_max, new_nb_lon, retstep=True) # list of pixel longitudes
# 			self.latitudes, self.dlat = np.linspace(lat_max, lat_min, new_nb_lat, retstep=True) # list of pixel latitudes
#
# 			self.I_map = cv2.resize(self.I_map, dsize=(new_nb_lon, new_nb_lat), interpolation = cv2.INTER_CUBIC)
#
#
# 		self.maps_shape = (self.Nt, Nb_e_pc, Nb_a_pc, len(self.latitudes), len(self.longitudes))
# 		# print("GROUND MAP SHAPE",self.maps_shape )
# 		#
# 		# self.scattering_map 		= np.zeros(self.maps_shape) # Intensity from (e, a) reaching us
# 		# self.DoLP_map				= np.zeros(self.maps_shape) # Polarized intensity from (e, a) reaching us
# 		# self.total_scattering_map 	= np.zeros(self.maps_shape) # DoLP of scattered light from (e,a)
# 		# self.AoRD_map 				= np.zeros(self.maps_shape) # Angle of polaisation of light from (e,a)
# 		# self.V_map 					= np.zeros(self.maps_shape)
# 		# self.Vcos_map 				= np.zeros(self.maps_shape)
# 		# self.Vsin_map 				= np.zeros(self.maps_shape)
# 		#
# 		# self.V_total	= np.zeros((Nb_e_pc, Nb_a_pc))
# 		# self.Vcos_total	= np.zeros((Nb_e_pc, Nb_a_pc))
# 		# self.Vsin_total	= np.zeros((Nb_e_pc, Nb_a_pc))
# 		#
# 		# self.I0_total	= np.zeros((Nb_e_pc, Nb_a_pc))
# 		# self.DoLP_total	= np.zeros((Nb_e_pc, Nb_a_pc))
# 		# self.AoLP_total	= np.zeros((Nb_e_pc, Nb_a_pc))
#
# 	def LoadPointSourceMap(self, src_I0, src_az, src_dist, Nb_a_pc, Nb_e_pc):
#
# 		self.src_I0, self.src_az, self.src_dist = src_I0, src_az, src_dist
#
# 		src_lon, src_lat = AzDistToLonLat(src_az, src_dist, self.A_lon, self.A_lat)
#
# 		self.longitudes, self.dlon = np.linspace(src_lon - 0.5/RT, src_lon + 0.5/RT, 2, retstep=True) # list of pixel longitudes
# 		self.latitudes, self.dlat = np.linspace(src_lat - 0.5/RT, src_lat + 0.5/RT, 2, retstep=True) # list of pixel latitudes
#
# 		print(self.longitudes)
# 		self.I_map = np.ones((1, 1)) * src_I0 #[nW / m2 / sr]
#
# 		self.maps_shape = (self.Nt, Nb_e_pc, Nb_a_pc, 1, 1)
#
# 		# self.scattering_map 		= np.zeros(self.maps_shape) # Polarized intensity from (e, a) reaching us
# 		# self.DoLP_map				= np.zeros(self.maps_shape) # DoLP of scattered light from (e,a)
# 		# self.total_scattering_map 	= np.zeros(self.maps_shape) # Intensity from (e, a) reaching us
# 		# self.AoRD_map 				= np.zeros(self.maps_shape) # Angle of polaisation of light from (e,a)
# 		# self.V_map 					= np.zeros(self.maps_shape)
# 		# self.Vcos_map 				= np.zeros(self.maps_shape)
# 		# self.Vsin_map 				= np.zeros(self.maps_shape)
# 		#
# 		# self.V_total	= np.zeros((Nb_e_pc, Nb_a_pc))
# 		# self.Vcos_total	= np.zeros((Nb_e_pc, Nb_a_pc))
# 		# self.Vsin_total	= np.zeros((Nb_e_pc, Nb_a_pc))
# 		#
# 		# self.I0_total	= np.zeros((Nb_e_pc, Nb_a_pc))
# 		# self.DoLP_total	= np.zeros((Nb_e_pc, Nb_a_pc))
# 		# self.AoLP_total	= np.zeros((Nb_e_pc, Nb_a_pc))
#
# 	def GetArea(self, ilat):
# 		"""Return the area of a pixel on the map in m**2. If we use a point source, the area is set to one."""
# 		if not self.is_point_source:
# 			return (RT ** 2) * abs(self.dlat) * abs(self.dlon) * np.cos(self.mid_latitudes[ilat])**1 * 1e6
# 		else:
# 			return 1
#
# 	def SetLightParameters(self):
# 		self.total_scattering_map 	= 2 * self.V_map
#
# 		self.DoLP_map				= 2 * np.sqrt(self.Vcos_map**2 + self.Vsin_map**2) / 100 # DoLP of scattered light from (e,a)
# 		self.DoLP_map				= np.divide(self.DoLP_map, self.V_map, out = np.zeros_like(self.DoLP_map), where = self.V_map != 0)
#
# 		self.AoRD_map 				= np.arctan2(self.Vsin_map, self.Vcos_map) / 2. # Angle of polaisation of light from (e,a)
# 		self.scattering_map 		= self.total_scattering_map * self.DoLP_map / 100. # Polarized intensity from (e, a) reaching us
#
# 		for it, t in enumerate(range(self.Nt)):
# 			for ie, e in enumerate(self.V_map[0]):
# 				for ia, a in enumerate(e):
# 					self.V_total[it, ie, ia] = np.sum(self.V_map[it, ie, ia, :, :])
# 					self.Vcos_total[it, ie, ia] = np.sum(self.Vcos_map[it, ie, ia, :, :])
# 					self.Vsin_total[it, ie, ia] = np.sum(self.Vsin_map[it, ie, ia, :, :])
#
# 					self.I0_total[it, ie, ia] = 2 * self.V_total[it, ie, ia]
# 					self.DoLP_total[it, ie, ia] = 100 * 2 * np.sqrt(self.Vcos_total[it, ie, ia] ** 2 + self.Vsin_total[it, ie, ia] ** 2) / self.V_total[it, ie, ia] #in %
# 					self.AoLP_total[it, ie, ia] = np.arctan2(self.Vsin_total[it, ie, ia], self.Vcos_total[it, ie, ia]) / 2
#
#
# 	def ProlongateTime(self, new_Nt):
# 		self.Nt = new_Nt
# 		self.cube = np.zeros((self.Nt, self.Nlat, self.Nlon))
#
# 		Nb_e_pc, Nb_a_pc = self.scattering_map.shape[1], self.V_total.shape[2]
# 		self.maps_shape = (self.Nt, Nb_e_pc, Nb_a_pc, self.Nlat, self.Nlon)
#
# 		self.scattering_map 		= np.zeros(self.maps_shape) # Intensity from (e, a) reaching us
# 		self.DoLP_map				= np.zeros(self.maps_shape) # Polarized intensity from (e, a) reaching us
# 		self.total_scattering_map 	= np.zeros(self.maps_shape) # DoLP of scattered light from (e,a)
# 		self.AoRD_map 				= np.zeros(self.maps_shape) # Angle of polaisation of light from (e,a)
# 		self.V_map 					= np.zeros(self.maps_shape)
# 		self.Vcos_map 				= np.zeros(self.maps_shape)
# 		self.Vsin_map 				= np.zeros(self.maps_shape)
#
# 		self.V_total	= np.zeros((self.Nt, Nb_e_pc, Nb_a_pc))
# 		self.Vcos_total	= np.zeros((self.Nt, Nb_e_pc, Nb_a_pc))
# 		self.Vsin_total	= np.zeros((self.Nt, Nb_e_pc, Nb_a_pc))
#
# 		self.I0_total	= np.zeros((self.Nt, Nb_e_pc, Nb_a_pc))
# 		self.DoLP_total	= np.zeros((self.Nt, Nb_e_pc, Nb_a_pc))
# 		self.AoLP_total	= np.zeros((self.Nt, Nb_e_pc, Nb_a_pc))
#
#
#


if __name__ == "__main__":
	in_dict = ReadInputFile("./input_files/RS_default.in")

	alt_map = ElevationMap(in_dict)

	# lon1, lat1, alt1 = alt_map.A_lon + 0/RT, alt_map.A_lat - 50/RT, alt_map.A_alt+6
	# lon2, lat2, alt2 = alt_map.A_lon, alt_map.A_lat, alt_map.A_alt
	I_map = np.zeros(alt_map.map.shape)

	for ilon, lon in enumerate(alt_map.longitudes):
		for ilat, lat in enumerate(alt_map.latitudes):
			if alt_map.IsVisible(lon, lat):
				I_map[ilon, ilat] = alt_map[ilon, ilat]

	# print("is_visible:", is_visible, clon, clat, calt)

	# x1, y1 = (lon1 - alt_map.A_lon) * RT, (lat1 - alt_map.A_lat) * RT
	# x2, y2 = (lon2 - alt_map.A_lon) * RT, (lat2 - alt_map.A_lat) * RT

	alt_map.PlotMap()
	f = plt.figure()
	plt.pcolormesh((alt_map.longitudes - alt_map.A_lon) * RT, (alt_map.latitudes - alt_map.A_lat) * RT, I_map)
	plt.colorbar()

	a_pc, e_pc = in_dict["azimuts"], in_dict["elevations"]
	dy = np.cos(a_pc) * 10 + alt_map.A_lat * RT
	dx = np.sin(a_pc) * 10 + alt_map.A_lon * RT
	plt.plot(np.linspace(alt_map.A_lon * RT, dx, 100), np.linspace(alt_map.A_lat * RT, dy, 100), color="red")

	# plt.plot([lon1*RT, lon2*RT], [lat1*RT, lat2*RT], "r")
	# plt.plot([x1, x2], [y1, y2], "r")

	# if not is_visible:
	# 	xcoll, ycoll = (clon - alt_map.A_lon) * RT, (clat - alt_map.A_lat) * RT
	# 	plt.plot(xcoll, ycoll, "*g")

	plt.show()
