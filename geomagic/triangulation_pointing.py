#!/usr/bin/python3
# -*-coding:utf-8 -*

import numpy as np
import observation as obs
from geometry import GetLonLatFromName

DtoR = np.pi / 180.
RtoD = 1. / DtoR
RT = 6371 # #km

##################################
# From the position and pointing of one instrument, find the direction of a second instrument in a known location such that they point toward the same point in the ionosphere
##################################

def GetRotMatrixAO(lonA, latA):
	"""Return a matrix 3x3. It is the rotation matrix to pass a vector expressed in the reference frame of a point A on the Earth (up-east-noth) to the reference frame of the Earth's center O. lonA latA are the longitude and latitude of the point A."""
	Clon, Slon = np.cos(lonA), np.sin(lonA)
	Clat, Slat = np.cos(latA), np.sin(latA)

	return     np.array([[	Clon*Clat,	-Slon,	-Clon*Slat],
						 [	Slon*Clat, 	Clon, 	-Slon*Slat],
						 [	Slat, 		0, 		Clat]])


### Point 0 is the known instrument, defining where we observe
name0 = 'eiscat_tromso' #"skibotnsud"
lon0, lat0 =GetLonLatFromName(name0) # 19.2123*DtoR, 69.34808*DtoR #
h = 110 #km
a0 = -180 * DtoR
e0 = 78 * DtoR
obs0 = obs.ObservationPoint(lon0, lat0, h, a0, e0)

### Longitude and latitude of the observed point H
lonH, latH = obs0.P_lon, obs0.P_lat
print("lonH, latH", lonH, latH)

### Lon, lat of the instrument 1 we want to know the direction
name1 =  "skibotn" #"kilpisjarvi"
# lat1, lon1 =  69.12602264116852*DtoR, 19.058594352049994*DtoR   # 68.924606*DtoR, 20.935396*DtoR
lon1, lat1 = GetLonLatFromName(name1) #19.21*DtoR, 70.5*DtoR

### vector from point 1 to point H in carthesian coord, then in up-east-north coordinate at point 1
v_OH = np.array((np.cos(lonH) * np.cos(latH), np.sin(lonH) * np.cos(latH), np.sin(latH))) #XYZ
v_O1 = np.array((np.cos(lon1) * np.cos(lat1), np.sin(lon1) * np.cos(lat1), np.sin(lat1))) #XYZ

v_1H = v_OH - v_O1 #XYZ
v_1H = np.dot(GetRotMatrixAO(lon1, lat1).transpose(), v_1H) #UEN
v_1H_norm = np.sqrt(v_1H[0]**2 + v_1H[1]**2 + v_1H[2]**2) * RT


# v_O1 = np.array((np.cos(lon1) * np.cos(lat1), np.sin(lon1) * np.cos(lat1), np.sin(lat1))) #XYZ
# v_O0 = np.array((np.cos(lon0) * np.cos(lat0), np.sin(lon0) * np.cos(lat0), np.sin(lat0))) #XYZ
# v_0H = np.array((np.sin(e0), np.cos(e0) * np.sin(a0), np.cos(e0) * np.cos(a0))) #uen
# v_0H = np.dot(GetRotMatrixAO(lon0, lat0), v_1H) #XYZ
#
# v_1H = - RT * v_O1 + RT * v_O0 + obs0.AH_norm * v_0H #XYZ
# print(v_O1)
# print(v_O0)
# print(v_0H)
# print(v_1H)
# v_1H = np.dot(GetRotMatrixAO(lon1, lat1).transpose(), v_1H) #UEN
# v_1H_norm = np.sqrt(sum([x**2 for x in v_1H]))

# dlon = (lonH - lon1) * RT * np.cos(lat1)
# dlat = (latH - lat1) * RT
# AH = RT * np.arccos(np.sin(latH) * np.sin(lat1) + np.cos(latH) * np.cos(lat1) * np.cos(lonH - lon1)) #See https://en.wikipedia.org/wiki/Great-circle_distance
# print("dlon, dlat", dlon*RT, dlat*RT)

###Azimut of instrument 1
a1 = np.arctan2(v_1H[1], v_1H[2])
# AH = np.sqrt(dlon ** 2 + dlat ** 2)
# print("AH", AH, v_1H_norm, AH - v_1H_norm)

###Elevation of instrument 1
e1 = np.arctan2(h, v_1H_norm)
# e1 = np.arcsin(v_1H[0] / v_1H_norm)
# print(v_1H)


obs1 = obs.ObservationPoint(lon1, lat1, h, a1, e1)

print("h:", h)
print(name0)
print("\ta:", a0 * RtoD)
print("\te:", e0 * RtoD)
print(name1)
print("\ta:", a1 * RtoD)
print("\te:", e1 * RtoD)

print((lonH - obs1.P_lon)*RT, (latH - obs1.P_lat)*RT)
