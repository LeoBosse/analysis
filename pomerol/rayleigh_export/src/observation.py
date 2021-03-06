#!/usr/bin/python3
# -*-coding:utf-8 -*
# First attempt in recreating the geometry calculated by Jean to get the magnetic field projection on the sensor

from subprocess import call
import numpy as np
import math as m
import matplotlib.pyplot as plt
import sys as sys


path = "../" #"/home/bossel/These/Analysis/src/Geometry/Leo"
DtoR = np.pi/ 180.
RtoD = 180. / np.pi
RT = 6371 # #km

class ObservationPoint:
	"""This is the class of ONE observation point, about which we know the position of the observer (A_lon, A_lat), the direction of observation (elevation, azimuth), and the height of the observed point. It can then give us the position of the observed point H, the magnetic field at this position (from CHAOS-6) and the apparent angle of this magnetic field on our captor. """

	def __init__(self, A_lon, A_lat, observed_altitude, azimuth, elevation, RD_src_azimut=None, RD_src_elevation=None, init_full = True):
		"""Initiating the class. The position of the observed point H is calculated automatically. Give everything in radians !"""
		self.lon, self.lat = A_lon, A_lat
		self.e, self.a = elevation, azimuth
		self.h = observed_altitude

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

	def SinglePointGeometry(self, GetBatP=True): #, A_lon, A_lat, h, a, e):
		if GetBatP:
			self.GetBatP()
		self.GetEta()
		if self.RD_src_azimut is not None and self.RD_src_elevation is not None:
			self.GetRayleighAngle(self.RD_src_azimut, self.RD_src_elevation)

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

	def GetOA(self):
		"""Return the vector OA in the reference frame of O. O:Centre of the Earth to A:observer"""
		OA = RT * np.array([[	 m.cos(self.lat) * m.cos(self.lon)],
							[	 m.cos(self.lat) * m.sin(self.lon)],
							[	 m.sin(self.lat)]])
		return OA

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
		AH_norm =  -RT * m.sin(e) + m.sqrt(RT**2 * m.sin(e)**2 + altitude**2 + 2*altitude*RT)
		AH = losO * AH_norm

		return AH, AH_norm

	def GetRotMatrixAO(self, lonA, latA):
		"""Return a matrix 3x3. It is the rotation matrix to pass a vector expressed in the reference frame of a point A on the Earth (up-east-noth) to the reference frame of the Earth's center O. lonA latA are the longitude and latitude of the point A."""
		Clon, Slon = m.cos(lonA), m.sin(lonA)
		Clat, Slat = m.cos(latA), m.sin(latA)

		return     np.array([[	Clon*Clat,	-Slon,	-Clon*Slat],
						     [	Slon*Clat, 	Clon, 	-Slon*Slat],
						     [	Slat, 		0, 		Clat]])

	def GetPCoordinates(self):
		"""Return the coordinates of the observed point H (or P). The position of the observer and the vector OA and AH should be calculated before."""

		###OH: Vector from Earth's center to observed point in reference frame of Earth's center O
		OH = self.OA + self.AH
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

	def GetBatP(self):
		"""Return the mag field of igrf and chaos at the observation point"""
		#Writing the input file for igrf and chaos code
		np.savetxt("theta_phi_H.dat", [[self.P_colat, self.P_lon]], fmt="%1.8e")

		#defining the stdout files for igrf and chaos (avoid all the text writen in the terminal)
		# f1 = open(path+"/src/igrf_auto.out", "w")
		f2 = open(path + "/src/chaos_auto.out", "w")

		#Calling igrf and chaos
#		call(["./igrf_auto"], stdout = f1)
		call([path + "/src/chaos_auto"], stdout = f2)
#		f1.close()
		f2.close()

		#Loading the results of the codes
#		B_igrf = np.loadtxt("Bxyz_H_igrf.dat")
		B_chaos = np.loadtxt(path + "/src/Bxyz_H_chaos6.dat")

		#Rearanging the vectors. The unit vector used in igrf and chaos are not the same here.
		#north-east-down to up-east-north
		#The first two values are the longitude and latitude of the observed point
		# self.B_igrf = [B_igrf[0], B_igrf[1], -B_igrf[4], B_igrf[3], B_igrf[2]]
		# self.B_igrf=[0,0,0,0,0]
		self.B_chaos = [B_chaos[0], B_chaos[1], -B_chaos[4], B_chaos[3], B_chaos[2]]
		# self.B_igrf = [B_igrf[0], B_igrf[1], 0, 0, 1]
		# self.B_chaos = [B_chaos[0], B_chaos[1], 0, 0, 1]

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

	def GetApparentAngle(self, Bp):
		"""Return the apparent angle of the input vector Bp (horizontal(1,3)) and the z axis of SPP"""
		#Reshaping of Bp horizontal to vertical (1,3) to (3,1)
		Bpu, Bpe, Bpn = Bp[0], Bp[1], Bp[2]
		Bp = np.array([[Bpu], [Bpe], [Bpn]])

		#All the usefull rotation matrices.
		Rpo = self.GetRotMatrixAO(self.P_lon, self.P_lat)
		Roa = self.GetRotMatrixAO(self.lon, self.lat).transpose()

		Ce, Se = m.cos(self.e), m.sin(self.e)
		Ca, Sa = m.cos(self.a), m.sin(self.a)
		Raspp = np.array([	[Se, 	Ce * Sa,	Ce * Ca],
						 [	 0,		-Ca,		Sa],
						 [	 Ce,	-Se * Sa,	-Se * Ca]])

		#Calculation of B in frame of SPP
		Bspp = np.dot(Raspp, np.dot(Roa, np.dot(Rpo, Bp)))

		#Calculation of the apparent angle between the apparent B and the SPP z axis.
		#The angle is 0 if the magnetic field is on the vertical axis pointing up, then turning clockwise
		BsppX, BsppY, BsppZ = Bspp[0][0], Bspp[1][0], Bspp[2][0]
		eta = m.atan2(BsppY, BsppZ)

		#Calculating the projection of B on the line of sight
		if np.linalg.norm(Bspp) != 0:
			self.Blos = np.arccos(BsppX / np.linalg.norm(Bspp))
		else:
			self.Blos = 0
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
				("B chaos norme",np.linalg.norm(self.B_chaos)),
				("Eta (angle apparent)", self.eta_chaos*RtoD),
				("Blos (angle B-lign of sight)", self.Blos*RtoD),
				("POllution source (a, e)", (self.RD_src_azimut * RtoD, self.RD_src_elevation * RtoD)),
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
