#!/usr/bin/python3
# -*-coding:utf-8 -*

import numpy as np
import observation as obs
from geometry import GetLonLatFromName

DtoR = np.pi / 180.
RtoD = 1. / DtoR
RT = 6371 # #km


I0 = 100
DoLP = 1 / 100
AoLP = 30 * DtoR

name0 = "skibotn"
lon0, lat0 = 0, 0 #GetLonLatFromName(name0)
h = 110 #km
a0 = -97 * DtoR
e0 = 45 * DtoR

obs0 = obs.ObservationPoint(lon0, lat0, h, a0, e0)
print(obs0)

plane0 = obs0.GetPolaPlane(AoLP)
print(plane0)

I0 = 100
DoLP = 1 / 100
AoLP = 40 * DtoR

name1 = "skibotn"
lon1, lat1 = 0, np.pi/2 #GetLonLatFromName(name0)
h = 110 #km
a1 = 24 * DtoR
e1 = 45 * DtoR

obs1 = obs.ObservationPoint(lon1, lat1, h, a1, e1)
print(obs1)

plane1 = obs1.GetPolaPlane(AoLP)
print(plane1)

current = np.cross(plane0, plane1) #in XYZ earth centered
print(current)

R_OH = obs0.GetRotMatrixAO(obs0.P_lon, obs0.P_lat)
R_HI = obs0.GetRotMatrixAI()

current = np.dot(R_OH, current)
print("Current in Up-East-North", current * RtoD)
current = np.dot(R_HI, current)
print("Current in instrument coord",current * RtoD)
