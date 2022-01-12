#!/usr/bin/python3
# -*-coding:utf-8 -*
# First attempt in recreating the geometry calculated by Jean to get the magnetic field projection on the sensor

from subprocess import call
import numpy as np
import math as m
import matplotlib.pyplot as plt
import sys as sys

from time import perf_counter #For the timer decorator

import chaosmagpy as chaos

def timer(fn):

    def inner(*args, **kwargs):
        start_time = perf_counter()
        to_execute = fn(*args, **kwargs)
        end_time = perf_counter()
        execution_time = end_time - start_time
        print('{0} took \t \t {1:.8f}s to execute'.format(fn.__name__, execution_time))
        return to_execute

    return inner


path = "/home/bossel/These/Analysis/src/geomagic"
DtoR = np.pi/ 180.
RtoD = 180. / np.pi
RT = 6371 # #km


class ObservationPoint:
	"""This is the class of ONE observation point, about which we know the position of the observer (A_lon, A_lat), the direction of observation (elevation, azimuth), and the height of the observed point. It can then give us the position of the observed point H, the magnetic field at this position (from CHAOS-6) and the apparent angle of this magnetic field on our captor. """

	#@timer
	def __init__(self, A_lon, A_lat, observed_altitude, azimuth, elevation, RD_src_azimut=None, RD_src_elevation=None, init_full = True, A_alt=0):
		"""Initiating the class. The position of the observed point H is calculated automatically. Give everything in radians !"""
		self.lon, self.lat = A_lon, A_lat
		self.e, self.a = elevation, azimuth
		self.h = observed_altitude
		self.A_alt = A_alt

		#Calculating the observed point position
		self.OA = self.GetOA() #Earth center O to observer A
		self.AH, self.AH_norm = self.GetAH() #Observer A to observed point H
		self.P_lon, self.P_lat = self.GetPCoordinates() #Observed point coord (P is H projected on ground)
		self.P_colat = self.GetCoLatitude(self.P_lat)
		self.RD_src_azimut = RD_src_azimut
		self.RD_src_elevation = RD_src_elevation
		self.AoRD = False

		self.B_igrf=[0,0,0,0,0]
		self.B_chaos=[0,0,0,0,0]
		self.B_model = None

		if init_full:
			self.B_model = chaos.load_CHAOS_matfile('/home/bossel/These/Analysis/data/magn_field/CHAOS-7.4.mat')

	def SinglePointGeometry(self, GetBatP=True, B_model = None): #, A_lon, A_lat, h, a, e):
		if GetBatP:
			self.GetBatP(B_model = B_model)
		self.GetEta()
		if self.RD_src_azimut is not None and self.RD_src_elevation is not None:
			self.GetRayleighAngle(self.RD_src_azimut, self.RD_src_elevation)


	def GetTrigo(self):
		Ce, Se = np.cos(self.e), np.sin(self.e)
		Ca, Sa = np.cos(self.a), np.sin(self.a)

		return Ce, Se, Ca, Sa

	def GetPCoordinatesFromRange(self, range):
		"""Get lon, lat, alt from range. Range is the distance to the instrument along the line of sight in km.
		alt is the absolute altitude (above see level)"""

		AH_vect = (self.AH / self.AH_norm) * range
		# print(AH_vect)
		lon, lat = self.GetPCoordinates(AH = AH_vect)

		# alt = self.A_alt + np.sin(self.e) * range
		alt = (self.A_alt+RT)**2 + range**2 + 2 * (self.A_alt+RT) * range * np.sin(self.e)
		alt = np.sqrt(alt) - RT

		return lon, lat, alt


	def GetRayleighAngle(self, source_azimut, source_elevation, unit="radians"): #obs_azimut, source_azimut, elevation, unit="radians"):
		if unit == "radians":
			eff_azimut = self.a - source_azimut
		elif unit == "degrees":
			source_azimut *= DtoR
			source_elevation *= DtoR
			eff_azimut = self.a - source_azimut

		src_u = np.sin(source_elevation)
		src_e = np.cos(source_elevation) * np.sin(source_azimut)
		src_n = np.cos(source_elevation) * np.cos(source_azimut)
		# src_los_eun = np.array([src_u, src_e, src_n])
		#
		# los_uen = self.los_uen.flatten()
		# n_uen = np.cross(los_uen, src_los_eun)
		# # print(source_azimut, source_elevation, n_uen)
		#
		Ce, Se = m.cos(self.e), m.sin(self.e)
		Ca, Sa = m.cos(self.a), m.sin(self.a)
		# Raspp = np.array([	[Se, 	Ce * Sa,	Ce * Ca],
		# 				 [	 0,		-Ca,		Sa],
		# 				 [	 Ce,	-Se * Sa,	-Se * Ca]])
		# n_ptcu = np.dot(Raspp, n_uen)

		# sin_AoRD = n_ptcu[1]
		# cos_AoRD = n_ptcu[2]
		sin_AoRD = Se*(src_n*Ca + src_e*Sa) - src_u*Ce
		cos_AoRD = src_n*Sa - src_e*Ca

		self.AoRD = m.atan2(sin_AoRD, cos_AoRD)
		# self.AoRD = np.arctan(np.sin(self.e) * np.cos(eff_azimut) / np.sin(eff_azimut))

		if self.AoRD < -np.pi/2:
			self.AoRD += np.pi
		elif self.AoRD > np.pi/2:
			self.AoRD -= np.pi

		return self.AoRD

	#@timer
	def GetOA(self, **kwargs):
		"""Return the vector OA in the reference frame of O. O:Centre of the Earth to A:observer"""
		if not kwargs:
			lon, lat, alt = self.lon, self.lat, self.A_alt
		else:
			lon, lat, alt = kwargs["lon"], kwargs["lat"], kwargs["alt"]

		OA = (RT + alt) * np.array([[	 m.cos(lat) * m.cos(lon)],
									[	 m.cos(lat) * m.sin(lon)],
									[	 m.sin(lat)]])
		return OA

	@staticmethod
	def GetAHNorm(elevation, altitude, A_alt=0):
		RT = 6371 #km
		return - (RT + A_alt) * m.sin(elevation) + m.sqrt(- m.cos(elevation)**2 * (RT + A_alt)**2 + (RT + altitude) ** 2)
		# return - (RT + A_alt) * m.sin(elevation) + m.sqrt(abs(- m.cos(elevation)**2 * (RT + A_alt)**2 + (RT + altitude) ** 2))

	#@timer
	def GetAH(self, **kwargs):
		"""Return the vector AH and its norm in the reference frame of O. From A:observer to H: the observed aurora"""

		#By default, use own elevation, azimut, altitude. But can compute it for a given set of parameters
		if "elevation" in kwargs: e = kwargs["elevation"]
		else : e = self.e
		if "azimut" in kwargs: a = kwargs["azimut"]
		else : a = self.a
		if "altitude" in kwargs: altitude = kwargs["altitude"]
		else : altitude = self.h

		###Lign of sight in reference frame of A (up-east-north and NOT e-n-u)
		los = np.array([[m.sin(e)],
						[m.cos(e) * m.sin(a)],
						[m.cos(e) * m.cos(a)]])
		self.los_uen = los
		###Get the rotation matrix to pass a vector from the reference frame of A to the reference frame of O
		Rotation = self.GetRotMatrixAO(self.lon, self.lat)

		###Pass the LOS in reference to O
		losO = np.dot(Rotation, los)

		#Norm of AH. Distance between observer and observed aurora
		# print(np.cos(e), RT + self.A_alt, - m.cos(e)**2 * (RT + self.A_alt)**2 + (RT + altitude) ** 2)

		AH_norm =  self.GetAHNorm(e, altitude, A_alt=self.A_alt)
		# AH_norm =  - (RT + self.A_alt) * m.sin(e) + m.sqrt(abs(- m.cos(e)**2 * (RT + self.A_alt)**2 + (RT + altitude) ** 2))
		AH = losO * AH_norm

		return AH, AH_norm

	def SpericalToCarth(r, lat, lon, colat=False):
		if colat:
			lat = np.pi/2 - lat

		x = np.cos(lon) * np.cos(lat)
		y = np.sin(lon) * np.cos(lat)
		z = np.sin(lat)

		return x, y, z


	def GetRotMatrixAO(self, lonA, latA, transpose=False):
		"""Return a matrix 3x3. It is the rotation matrix to pass a vector expressed in the reference frame of a point A on the Earth (up-east-north) to the reference frame of the Earth's center O. lonA latA are the longitude and latitude of the point A."""
		Clon, Slon = m.cos(lonA), m.sin(lonA)
		Clat, Slat = m.cos(latA), m.sin(latA)

		R_ao =   np.array([  [	Clon*Clat,	-Slon,	-Clon*Slat],
						     [	Slon*Clat, 	Clon, 	-Slon*Slat],
						     [	Slat, 		0, 		Clat]])

		if transpose:
			R_ao = np.transpose(R_ao)

		return R_ao

	def GetTransormationMatrixHA(self, transpose = False):
		R_HO = self.GetRotMatrixAO(self.P_lon, self.P_lat, transpose = transpose)
		R_OA = self.GetRotMatrixAO(self.lon, self.lat, transpose = not transpose)

		if not transpose:
			return R_OA @ R_HO #matrix multiplication
		else:
			return R_HO @ R_OA #matrix multiplication


	#@timer
	def GetPCoordinates(self, AH=None):
		"""Return the coordinates of the observed point H (or P). The position of the observer and the vector OA and AH should be calculated before."""

		if AH is None:
			AH = self.AH

		###OH: Vector from Earth's center to observed point in reference frame of Earth's center O
		OH = self.OA + AH
		# OH_norm = np.sqrt(sum([x**2 for x in OH.flat]))

		###Calculating longitude and latitude of P
		P_lon = m.atan2(OH[1][0], OH[0][0])
		P_lat = np.pi/2. - m.atan2(np.sqrt(OH[0][0]**2 + OH[1][0]**2), OH[2][0])

		return P_lon, P_lat

	def GetCoLatitude(self, latitude, unit="radian"):
		"""Return the colatitude corresponding to a certain latitude in the unit given as parameter. colat = 90-lat"""
		if unit =="degree":
			return 90 - latitude
		else:
			return np.pi/2 - latitude


	def GetPolaPlane(self, AoLP):
		"""Return the polarization plane for the observation in Earth centered coordinates. Plane coordinates in an array (a, b, c, d) such that aX+bY+cZ+d=0"""

		# A = (0, np.sin(AoLP), np.cos(AoLP))
		# los = (1, 0, 0)

		plane = np.array((0, -np.cos(AoLP), np.sin(AoLP))).transpose() #in instrument coordinates (x*,y*,z*). THis is the vector perpendicular to the plane
		print("plane", plane)

		R_IA = self.GetRotMatrixAI(transpose = True)
		plane = np.dot(R_IA, plane) #plane expressed in up-east-north coordinate
		print("plane", plane)

		R_A0 = self.GetRotMatrixAO(self.lon, self.lat)
		plane = np.dot(R_A0, plane) #plane expressed in earth centered ref frame
		print("plane", plane)

		return plane

	def GetBatP(self, time=None, B_model = None):
		"""Return the mag field of igrf and chaos at the observation point"""
# 		#Writing the input file for igrf and chaos code
# 		np.savetxt(path + "/B_model/theta_phi_H.dat", [[self.P_colat, self.P_lon]], fmt="%1.8e")
#
# 		#defining the stdout files for igrf and chaos (avoid all the text writen in the terminal)
# 		# f1 = open(path+"/src/igrf_auto.out", "w")
# 		f2 = open(path + "/B_model/chaos_auto.out", "w")
#
# 		#Calling igrf and chaos
# #		call(["./igrf_auto"], stdout = f1)
# 		call([path + "/B_model/chaos_auto"], stdout = f2)
# #		f1.close()
# 		f2.close()
#
# 		#Loading the results of the codes
# #		B_igrf = np.loadtxt("Bxyz_H_igrf.dat")
# 		B_chaos = np.loadtxt(path + "/B_model/Bxyz_H_chaos6.dat")
#
# 		#Rearanging the vectors. The unit vector used in igrf and chaos are not the same here.
# 		#north-east-down to up-east-north
# 		#The first two values are the longitude and latitude of the observed point
# 		# self.B_igrf = [B_igrf[0], B_igrf[1], -B_igrf[4], B_igrf[3], B_igrf[2]]
# 		# self.B_igrf=[0,0,0,0,0]
# 		self.B_chaos = [B_chaos[0], B_chaos[1], -B_chaos[4], B_chaos[3], B_chaos[2]]
# 		# self.B_igrf = [B_igrf[0], B_igrf[1], 0, 0, 1]
# 		# self.B_chaos = [B_chaos[0], B_chaos[1], 0, 0, 1]

		if B_model is None:
			if self.B_model is None:
				B_model = chaos.load_CHAOS_matfile('/home/bossel/These/Analysis/data/magn_field/CHAOS-7.4.mat')
			else:
				B_model = self.B_model
		else:
			self.B_model = B_model

		if time is None:
			time = chaos.data_utils.mjd2000(2019, 1, 1)

		B_radius, B_theta, B_phi = B_model.synth_values_tdep(time, RT + self.h, self.P_colat*RtoD, self.P_lon*RtoD) #given in up, south, east

		# print((-B_chaos[4]-B_radius)/B_radius,( B_chaos[3]-B_phi)/B_phi,( B_chaos[2]+B_theta)/B_theta)

		self.B_chaos = [self.P_colat*RtoD, self.P_lon*RtoD, B_radius, B_phi, -B_theta] #transformed in colat, lon, up, east, north

		return self.B_igrf, self.B_chaos

	def GetEta(self, quadrant = True):
		"""Return the angle between the apparent magnetic field (given by igrf and chaos) and the zaxis of SPP. Si Quadrant==True permet de renvoyer un angle entre -pi/2 et pi/2, au lieu qu'il soit entre -pi et pi"""

		#Get the apparent angle for igrf and chaos
		self.eta_igrf = self.GetApparentAngle(self.B_igrf[2:])
		self.eta_chaos = self.GetApparentAngle(self.B_chaos[2:])

		#If quadrant==True, angle between -pi/2 and pi/2, instead of -pi et pi
		if quadrant:
			if self.eta_igrf <= -np.pi/2.:
				self.eta_igrf += np.pi
			elif self.eta_igrf > np.pi/2.:
				self.eta_igrf -= np.pi

			if self.eta_chaos <= -np.pi/2.:
				 self.eta_chaos += np.pi
			elif self.eta_chaos > np.pi/2.:
				 self.eta_chaos -= np.pi
		# self.eta = self.eta_chaos #because I only use chaos, both models are very close together
		return self.eta_igrf, self.eta_chaos

	def GetRotMatrixAI(self, transpose = False):
		""""Return a matrix 3x3. It is the rotation matrix to pass a vector expressed in the reference frame of a point A on the Earth (up-east-north) to the reference frame of the instrument (line of sight-west-up) if instrument points northern horyzon."""

		Ce, Se, Ca, Sa = self.GetTrigo()

		Raspp = np.array([	[Se, 	Ce * Sa,	Ce * Ca],
						 [	 0,		-Ca,		Sa],
						 [	 Ce,	-Se * Sa,	-Se * Ca]])

		if transpose:
			Raspp = np.transpose(Raspp)

		return Raspp


	def GetApparentAngle(self, Bp, Blos=False):
		"""Return the apparent angle of the input vector Bp (horizontal(1,3)) and the z axis of SPP. Angle returned in radians. 0 is vertical, pi/2 horizontal."""
		#Reshaping of Bp horizontal to vertical (1,3) to (3,1)
		Bpu, Bpe, Bpn = Bp[0], Bp[1], Bp[2]
		Bp = np.array([[Bpu], [Bpe], [Bpn]])

		#All the usefull rotation matrices.
		Rpo = self.GetRotMatrixAO(self.P_lon, self.P_lat)
		Roa = self.GetRotMatrixAO(self.lon, self.lat).transpose()
		Raspp = self.GetRotMatrixAI()

		#Calculation of B in frame of SPP
		Bspp = np.dot(Raspp, np.dot(Roa, np.dot(Rpo, Bp)))

		#Calculation of the apparent angle between the apparent B and the SPP z axis.
		#The angle is 0 if the magnetic field is on the vertical axis pointing up, then turning clockwise
		BsppX, BsppY, BsppZ = Bspp[0][0], Bspp[1][0], Bspp[2][0]
		eta = m.atan2(BsppY, BsppZ)


		#Calculating the projection of B on the line of sight
		if np.linalg.norm(Bspp) != 0:
			# print(BsppX, np.linalg.norm(Bspp), abs(np.arccos(BsppX / np.linalg.norm(Bspp))))
			self.Blos = abs(np.arccos(BsppX / np.linalg.norm(Bspp)))
			if self.Blos > np.pi/2:
				self.Blos = np.pi - self.Blos
		else:
			self.Blos = 0

		if Blos:
			return eta, self.Blos
		else:
			return eta

	def GetAllParameters(self):
		"""Return a list of all interesting parmaters of the observation to print in a terminal"""
		list = [ ("lon", self.lon*RtoD),
				("lat",self.lat*RtoD),
				("elevation", self.e*RtoD),
				("azimuth",self.a*RtoD),
				("altitude", self.h),
				("Hlon",self.P_lon*RtoD),
				("Hlat",self.P_lat*RtoD),
				("Norme AH", self.AH_norm),
				("B chaos (up, east, north)", self.B_chaos[2:]),
				("B chaos norme",np.linalg.norm(self.B_chaos[2:])),
				("Eta (angle apparent)", self.eta_chaos*RtoD),
				("Blos (angle B-lign of sight)", self.Blos*RtoD),
				# ("Pollution source (a, e)", (self.RD_src_azimut * RtoD, self.RD_src_elevation * RtoD)),
				("AoRD", self.AoRD * RtoD)]
		if self.AoRD is False:
			list[-1] = ("AoRD", self.AoRD)
		return list

	def PrintAllParameters(self):
		"""Print the list of all interesting parameters in a terminal"""
		for p in self.GetAllParameters():
			print(p[0], "\t", p[1])


	def __repr__(self):
		s = "Observation at point ({0}, {1}), azimuth {2}, elevation {3}. Observed point at altitude {4}, position ({5}, {6})".format(self.lon*RtoD, self.lat*RtoD, self.a*RtoD, self.e*RtoD, self.h, self.P_lon*RtoD, self.P_lat*RtoD)
		return s

class ObservationToPoint(ObservationPoint):
	def __init__(self, A_lon, A_lat, B_lon, B_lat, B_alt, RD_src_azimut=None, RD_src_elevation=None, init_full = True, A_alt=0):

		self.lon, self.lat = A_lon, A_lat
		self.h = B_alt
		self.A_alt = A_alt

		OA = self.GetOA() #Earth center O to observer A
		OB = self.GetOA(lon=B_lon, lat=B_lat, alt=B_alt) #Earth center O to observed point B

		AB = OB - OA #in XYZ earth center
		AB = np.dot(self.GetRotMatrixAO(A_lon, A_lat).transpose(), AB) #UEN
		AB_norm = np.sqrt(AB[0]**2 + AB[1]**2 + AB[2]**2)

		elevation 	= np.arcsin(AB[0] / AB_norm)
		azimuth		= np.arctan2(AB[1], AB[2])

		# print("ObservationToPoint: az, el:", azimuth*RtoD, elevation*RtoD)

		ObservationPoint.__init__(self, A_lon, A_lat, B_alt, azimuth, elevation, RD_src_azimut, RD_src_elevation, init_full, A_alt)
