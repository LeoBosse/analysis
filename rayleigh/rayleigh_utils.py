#!/usr/bin/python3
# -*-coding:utf-8 -*

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.patches import Ellipse

from time import perf_counter #For the timer decorator

# import osgeo.gdal as gdal
# gdal.UseExceptions()  # not required, but a good idea

import sys as sys
from observation import *

import astropy as apy
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from astropy import units as u

import numba as nba

DtoR = np.pi / 180.
RtoD = 1. / DtoR
RT = 6371. # km

# Initialisation of a few usefull geometrical functions
Getuen 			= lambda a, e: (np.sin(e), np.cos(e) * np.sin(a), np.cos(e) * np.cos(a)) # From azimut, elevation, return a vector in the Up-East-North coordinate system
GetAngle 		= lambda v1, v2: np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0)) # Return angle in radian between two vectors
GenPythagore	= lambda a, b, theta: np.sqrt(a**2 + b**2 - 2*a*b*np.cos(theta)) # Generalised Pythagorean theorem. From 2 sides and the angle between them, return the length of third side.

RadToKm			= lambda rad: rad * RT
# Timer decorator

def timer(fn):

    def inner(*args, **kwargs):
        start_time = perf_counter()
        to_execute = fn(*args, **kwargs)
        end_time = perf_counter()
        execution_time = end_time - start_time
        print('{0} took \t \t {1:.8f}s to execute'.format(fn.__name__, execution_time))
        return to_execute

    return inner

def GetGalacticCoord(lon, lat, height, az, el, time):

	location = EarthLocation(lon=lon, lat=lat, height=height*u.m)
	ptcu = SkyCoord(obstime = time, location=location, alt=el*u.rad, az=az*u.rad, frame='altaz')
	Glon, Glat = ptcu.galactic.l.value, ptcu.galactic.b.value

	# print("GLon, Glat", Glon, Glat)

	return Glon, Glat #in degrees

def GetStarlight(lon, lat, height, time, el, az):

	Glon, Glat = GetGalacticCoord(lon, lat, height, az, el, time)

	b = np.array([0, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80])
	Jv = np.array([371.5, 246.5, 175.3, 136.5, 102.9, 68.7, 51.5, 41.3, 35.5, 32.5, 31.1]) #in S_10 units
	# Jv = self.S10TonanoWatt(Jv)

	return np.interp(abs(Glat), b, Jv)


def GetERForVeryBigE_alt(R_alt, E_alt, RE, alt_max):
	"""When the emission point in the sky is higher than the maximum altitude of the atmosphere, it creates troubles when transforming the vertical optical depth to the effective one. We compute the length of the segment of the line joining R and E, from R to the top of the atmosphere.
	R_alt: altitude of R (the scattering point) in km
	E_alt: altitude of E (the scattering point) in km
	RE: distance between points E and R
	alt_max: maximum altitude of the atmosphere in km. should be < E_alt
	"""
	if E_alt < alt_max:
		return False

	# radius = RT + alt_max
	# xR, yR = 0, RT + R_alt
	# R_elevation = GetArbitraryElevation(R_alt+RT, E_alt+RT, RE)
	# tan_e = np.tan(R_elevation)

	R_elevation = GetArbitraryElevation(RT + R_alt, RT + E_alt, RE)
	ER_eff = GetArbitraryDistance(RT + alt_max, RT + alt_max, R_elevation)

	# if tan_e > 0:
	# 	# y = 0.5 * (yR/tan_e + np.sqrt(1./tan_e**2 + 4*(radius - yR/tan_e)) )
	# 	x = (-tan_e * yR) + np.sqrt(radius**2 * (1 + tan_e**2) - yR**2)
	# else:
	# 	# y = 0.5 * (yR/tan_e - np.sqrt(1./tan_e**2 + 4*(radius - yR/tan_e)) )
	# 	x = (-tan_e * yR) - np.sqrt(radius**2 * (1 + tan_e**2) - yR**2)
	# x /= (1 + tan_e**2)
	#
	# # x = (y - yR) / tan_e
	# y = tan_e * x + yR
	#
	# ER_eff = np.sqrt((x-xR)**2 + (y-yR)**2)

	# print("ER_eff", R_alt, R_elevation*RtoD, x-xR, y-yR, ER_eff)
	return ER_eff

def GetArbitraryElevation(rA, rB, AB):
	"""Two concentric circles of radius rA and rB, and two points A and B on each circle (A on rA and B on rB). Given rA, rB and the distance AB between the two points, returns the elevation angle at which we see the outer point from the inner point. The elevation is the angle of the segment AB with the tangent to the inner circle passing by the inner point. we use the law of cosine, or generalized pythagorean theorem: c**2 = a**2 + b**2 - 2ab*cos(C) with abc the side length of a triangle and C the angle between a and b."""
	if rB < rA:
		rA, rB = rB, rA
	sin_e = (rB**2 - AB**2 - rA**2) / (2 * AB * rA)
	e = np.arcsin(sin_e)
	# print("Arb el", sin_e, e*RtoD)
	return e


def GetArbitraryDistance(rA, rB, elevation):
	"""Two concentric circles of radius rA and rB, and two points A and B on each circle (A on rA and B on rB). Given rA, rB and the elevation angle at which we see the outer point from the inner point, returns the distance AB between the two points. The elevation is the angle of the segment AB with the tangent to the inner circle passing by the inner point. We use the law of cosine, or generalized pythagorean theorem: c**2 = a**2 + b**2 - 2ab*cos(C) with abc the side length of a triangle and C the angle between a and b."""
	if rB < rA:
		rA, rB = rB, rA
	AB = - rA * np.sin(elevation) + np.sqrt(rB**2 - rA**2 * np.cos(elevation)**2)
	return AB

def GetArbitraryAltitude(rA, AB, e):
	"""Two concentric circles of radius rA and rB, and two points A and B on each circle (A on rA and B on rB). Given rA, AB and the elevation angle at which we see the outer point from the inner point, returns the distance rB the radius of the outer circle. The elevation is the angle of the segment AB with the tangent to the inner circle passing by the inner point. We use the law of cosine, or generalized pythagorean theorem: c**2 = a**2 + b**2 - 2ab*cos(C) with abc the side length of a triangle and C the angle between a and b."""
	rB = rA**2 + AB**2 + 2 * rA * AB * np.sin(e)
	return np.sqrt(rB) - RT #Sutract the earth radius to obtain the altitude of point B



# @timer
# @nba.njit
def DistToAngle(dist, lat = 0, az = 0, R = RT):
	return (dist  / R) * np.sqrt(np.cos(az)**2 + np.sin(az)**2 / np.cos(lat)**2)

# @timer
# @nba.njit
def AngleToDist(angle, lat = 0, az = 0, R = RT):
	return (angle * R) / np.sqrt(np.cos(az)**2 + np.sin(az)**2 / np.cos(lat)**2)


# @timer
# @nba.njit
def AzDistToLonLat(a, d, origin_lon = 0, origin_lat = 0):
	"""Azimut a in radians, distance d in km"""
	dist_east, dist_north = d * np.sin(a), d * np.cos(a)
	# dist_east, dist_north = DistToAngle(d, lat=origin_lat, az = a), d * np.cos(a)

	delta_lon = dist_east / np.cos(origin_lat) / RT
	delta_lat = dist_north / RT

	delta_lon = DistToAngle(dist_east, lat=origin_lat, az = 90)
	delta_lat = DistToAngle(dist_north, lat=origin_lat, az = 0)

	# print("AzDistToLonLat", dist_east, dist_north, delta_lon, delta_lat, delta_lon*RT, delta_lat*RT)

	return origin_lon + delta_lon, origin_lat + delta_lat

def GetRotMatrixAO(lonA, latA):
	"""Return a matrix 3x3. It is the rotation matrix to pass a vector expressed in the reference frame of a point A on the Earth (up-east-noth) to the reference frame of the Earth's center O. lonA latA are the longitude and latitude of the point A."""
	Clon, Slon = np.cos(lonA), np.sin(lonA)
	Clat, Slat = np.cos(latA), np.sin(latA)

	return     np.array([[	Clon*Clat,	-Slon,	-Clon*Slat],
						 [	Slon*Clat, 	Clon, 	-Slon*Slat],
						 [	Slat, 		0, 		Clat]])


def AzEltoLonLat(a, e, h, A_lon, A_lat):
	obs = ObservationPoint(A_lon, A_lat, h, a, e, init_full=False)
	return obs.P_lon, obs.P_lat

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

	# dist_in_rad = AE_uen[1], AE_uen[2]

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
	"""From known geometry parameters : a, e, alt of emission and alt of scatering, return missing parameters: Distance between emission and scattering and angle of scattering.
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

	AER = np.pi - ARE - e_rd

	return AR, RE, ARE, AER

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

	AER = np.pi - ARE - e_rd

	return AR, RE, ARE, AER


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
		A_lon = 20.363641
		A_lat = 69.348151
	elif name == "skibotnsud":
		A_lon = 19.981116
		A_lat = 69.234956
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
