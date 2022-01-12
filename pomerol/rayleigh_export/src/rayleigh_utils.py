#!/usr/bin/python3
# -*-coding:utf-8 -*

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.patches import Ellipse

import osgeo.gdal as gdal
gdal.UseExceptions()  # not required, but a good idea

import sys as sys
from observation import *

DtoR = np.pi / 180.
RtoD = 1. / DtoR
RT = 6371. # km

# Initialisation of a few usefull geometrical functions
Getuen 			= lambda a, e: (np.sin(e), np.cos(e) * np.sin(a), np.cos(e) * np.cos(a)) # From azimut, elevation, return a vector in the Up-East-North coordinate system
GetAngle 		= lambda v1, v2: np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0)) # Return angle in radian between two vectors
GenPythagore	= lambda a, b, theta: np.sqrt(a**2 + b**2 - 2*a*b*np.cos(theta)) # Generalised Pythagorean theorem. From 2 sides and the angle between them, return the lengthof third side.

RadToKm			= lambda rad: rad * RT

def AzDistToLonLat(a, d, origin_lon = 0, origin_lat = 0):
	dist_east, dist_north = d * np.cos(0) * np.sin(a), d * np.cos(0) * np.cos(a)

	delta_lon = dist_east / RT * np.cos(origin_lat)
	delta_lat = dist_north / RT

	# print(dist_east, dist_north, delta_lon, delta_lat)

	return origin_lon + delta_lon, origin_lat + delta_lat

def GetRotMatrixAO(lonA, latA):
	"""Return a matrix 3x3. It is the rotation matrix to pass a vector expressed in the reference frame of a point A on the Earth (up-east-noth) to the reference frame of the Earth's center O. lonA latA are the longitude and latitude of the point A."""
	Clon, Slon = np.cos(lonA), np.sin(lonA)
	Clat, Slat = np.cos(latA), np.sin(latA)

	return     np.array([[	Clon*Clat,	-Slon,	-Clon*Slat],
						 [	Slon*Clat, 	Clon, 	-Slon*Slat],
						 [	Slat, 		0, 		Clat]])


def LonLatToAzEl(E_lon, E_lat, h, A_lon, A_lat):
	OE = (RT + h) * np.array([		[np.cos(E_lon) * np.cos(E_lat)],
						[np.sin(E_lon) * np.cos(E_lat)],
						[np.sin(E_lat)]
	])
	OA = RT * np.array([		[np.cos(A_lon) * np.cos(A_lat)],
						[np.sin(A_lon) * np.cos(A_lat)],
						[np.sin(A_lat)]
	])

	AE = OE - OA
	R_OA = GetRotMatrixAO(A_lon, A_lat).T

	AE_uen = np.dot(R_OA, AE).T[0]
	AE_uen = AE_uen / np.sum(AE_uen)

	azimut 		= np.arctan2(AE_uen[1], AE_uen[2])
	elevation 	= GetAngle(AE_uen, (0, AE_uen[1], AE_uen[2]))

	return azimut, elevation


def LonLatToAzDist(E_lon, E_lat, A_lon, A_lat):
	OE = np.array([		[np.cos(E_lon) * np.cos(E_lat)],
						[np.sin(E_lon) * np.cos(E_lat)],
						[np.sin(E_lat)]
	])
	OA = np.array([		[np.cos(A_lon) * np.cos(A_lat)],
						[np.sin(A_lon) * np.cos(A_lat)],
						[np.sin(A_lat)]
	])

	AE = OE - OA
	R_OA = GetRotMatrixAO(A_lon, A_lat).T

	AE_uen = np.dot(R_OA, AE).T[0]

	# print(AE_uen)

	azimut = np.arctan2(AE_uen[1], AE_uen[2])
	# print(azimut*RtoD, "\t", AE_norm)
	# azimut = GetAngle(np.array([1, 0, 0]), AE_uen.reshape(1,3))

	dist_in_rad = GetAngle(OE.T[0], OA.T[0])

	return azimut, dist_in_rad

# def AzDistFromLonLat(E_lon, E_lat, h, A_lon, A_lat):
# 	"""From known geometry parameters : lon, lat, alt of emission and position of the instrument, return missing parameters: Distance between emission and scattering and angle of scattering."""
#
# 	# OE = np.array([	np.cos(E_lon) * np.cos(E_lat),
# 	# 						np.sin(E_lon) * np.cos(E_lat),
# 	# 						np.sin(E_lat)
# 	# ])
# 	#
# 	# e_rd = GetAngle(obs.OA.reshape(1, 3) / RT, OE)
# 	a_rd, d = LonLatToAzDist(E_lon, E_lat, obs.lon, obs.lat)
# 	e_rd = d / RT
#
# 	return a_rd, e_rd
# 	# GetGeometryFromAzEl(a_rd, e_rd, h, h_r, v_pc_u, obs)


def GetGeometryFromAzEl(a_rd, e_rd, h, h_r, v_pc_u, obs):
	"""From known geometry parameters : a, e, alt of emission and alt of scatering , return missing parameters: Distance between emission and scattering and angle of scattering.
	A instrument positions pointing to (a_pc, e_pc)
	E emission points (a_rd, e_rd, h)
	R rayleigh diffusion point (a_pc, e_pc, h_r)"""

	e_rd = max(e_rd, 10**(-30)) # Avoid weird behaviour and division by 0

	v_rd_u = Getuen(a_rd, e_rd) # Vector from observer to emission point
	_, AE = obs.GetAH(elevation = e_rd, azimut = a_rd, altitude = h)  # Distance between observer and emisison point

	RAE = GetAngle(v_pc_u, v_rd_u) # Angle between line of sight and emission

	_, AR = obs.GetAH(altitude = h_r) # Distance between observer and scattering point

	RE = GenPythagore(AR, AE, RAE) # Distance between emission and scaterring point
	ARE = np.arcsin(AE * np.sin(RAE) / RE) # Scattering angle

	return AR, RE, ARE

def GetGeometryFromAzDist(a_rd, d, h_r, v_pc_u, obs):
	"""From known geometry parameters : a, d, alt of emission and alt of scatering , return missing parameters: Distance between emission and scattering and angle of scattering.
	A instrument positions pointing to (a_pc, e_pc)
	E emission points (a_rd, e_rd, h)
	R rayleigh diffusion point (a_pc, e_pc, h_r)"""

	d = max(d, 10**(-30)) # Avoid weird behaviour and division by 0

	v_rd_u = Getuen(a_rd, 0.) # Vector from observer to emission point (flat ground, emission at elevation 0)
	AE = d * RT # Distance between observer and emisison point. (exact, following the curvature (~flat earth))

	RAE = GetAngle(v_pc_u, v_rd_u) # Angle between line of sight and emission

	_, AR = obs.GetAH(altitude=h_r) # Distance between observer and scattering point

	RE = GenPythagore(AR, AE, RAE) # Distance between emission and scaterring point
	ARE = np.arcsin(AE * np.sin(RAE) / RE) # Scattering angle

	return AR, RE, ARE

def GetVolume(AR, d_los, ouv_pc):
	"""Get the volume of a truncated cone of kength d_los, and opening angle ouv_pc."""
	V = (np.pi * np.tan(ouv_pc) ** 2) / 3

	h1, h2 = AR - d_los, AR + d_los

	V *= (h2 - h1)
	V *= ((h2 + h1) ** 2 - h1 * h2)

	# print("V", V * (10 ** 9))
	return V * (10 ** 9) # Every distances is in km

def GetRSCrossSectionParticle(theta, wl = 557.7, Dp = 0.315, n = 1.000271375):
	"""Get Rayleigh scattering cross section for 1 angle. (See wikipedia or Seinfeld, John H. and Pandis, Spyros N. (2006) Atmospheric Chemistry and Physics, 2nd Edition, John Wiley and Sons, New Jersey, Chapter 15.1.1, ISBN 0471720186).
	Dp : particle size
	wl : wavelength of the light
	theta : scattering angle
	n : refractive index"""

	cs  = ((Dp * 10 ** (-9)) ** 6) / 8
	cs *= (np.pi / (wl * 10 ** (-9))) ** 4
	cs *= ((n ** 2 - 1) / (n ** 2 + 2)) ** 2
	cs *= (1 + np.cos(theta) ** 2)

	# print("cs", cs)

	return cs


def GetRSCrossSectionMolecule(theta, wl = 557.7, alpha=2.3732 * 10**(-39)):
	"""Get Rayleigh scattering cross section for 1 angle. (See wikipedia).
	wl : wavelength of the light
	theta : scattering angle
	alpha : polarizability (in SI units !!! alpha(cm**3) = 8.988*10**15 alpha(SI) """

	cs  = (8 * np.pi ** 4 * alpha ** 2) / (wl ** 4)
	cs *= (1 + np.cos(theta) ** 2)

	return cs

def GetParticuleDensity(h = 0):
	"""Should give the particle density at altitude h in particules/m3. For now it's constant"""
	return 2.737 * (10 ** 25)


def GetScattered(I0, AR, ER, RD_angle, ouv_pc, alt, d_los, elevation = 0):
	"""Given an initial intensity of a source and some geometrical parameter, returns the intensity mesured at the instrument and its DoLP.
	Input parameters: elevation, altitude of source, scattering angle, distance between emission and scattering."""

	w_I  = 	I0 * GetVolume(AR, d_los, ouv_pc) * GetParticuleDensity(alt) * GetRSCrossSectionParticle(RD_angle)

	if (ER and AR) != 0:
		w_I /= ER ** 2 	# Inverse square law
		w_I /= AR ** 2 	# Inverse square law
	else:
		w_I = 0

	w_DoLP = (1 - np.cos(RD_angle)**2) / (1 + np.cos(RD_angle)**2) # DoLP dependance on scattering angle

	return w_I, w_DoLP


def GetLonLatFromName(name, unit="radians"):
	name = name.lower()
	if name == "mens":
		A_lon = 5.76
		A_lat = 44.83
	elif name == "skibotn":
		A_lon = 20.363109
		A_lat = 69.348027
	elif name == "nyalesund":
		A_lon = 11.92288
		A_lat = 78.92320
	elif name == "vigan":
		A_lon = 3.504259
		A_lat = 44.039661
	elif name == "lagorge":
		A_lon = 5.936935
		A_lat = 45.212343
	elif name == "stveran":
		A_lon = 6.906917
		A_lat = 44.696926
	else:
		A_lon = float(name.split(";")[0])
		A_lat = float(name.split(";")[1])

	if unit == "radians":
		A_lon *= DtoR
		A_lat *= DtoR

	return A_lon, A_lat

def ReadInputFile(filename):
		"""Read a given input file and return a dictionnary. First word = key, second = value. Remove empty lines, as many arguments as you want, in any order that you want."""
		with open(filename, "r") as f:
			input = f.readlines()
		input = [l.split() for l in input if l != "\n"]
		dict = {}
		for i in input:
			dict[i[0]] = i[1]

		return dict
