#!/usr/bin/python3
# -*-coding:utf-8 -*
from observation import *
from subprocess import call
import numpy as np
import math as m
import matplotlib.pyplot as plt
import sys as sys


DtoR = np.pi/ 180.
RtoD = 180. / np.pi


N = 50

###Altitude of the observed emission
h = 110 #km

#Earth's radius
RT = 6371 #km


#A_lon = np.linspace(0, 2*np.pi, N)
#A_lat = np.linspace(-np.pi / 2., np.pi / 2., N)

###MENS Coordinates
A_lon = 5.76 * DtoR
A_lat = 45 * DtoR#44.83 * DtoR

###Skibotn Coordinates
#A_lon = 20.24 * DtoR
#A_lat = 69.34 * DtoR

fig = plt.figure()

ax1 = fig.add_subplot(321)
ax1.set_xlabel("elevation")
ax1.set_ylabel("H longitude")
ax2 = fig.add_subplot(322)
ax2.set_xlabel("elevation")
ax2.set_ylabel("H latitude")
ax3 = fig.add_subplot(323)
ax3.set_xlabel("elevation")
ax3.set_ylabel("Eta, projected angle")
ax4 = fig.add_subplot(324)
ax4.set_xlabel("azimut")
ax4.set_ylabel("Blos")
ax5 = fig.add_subplot(325)
ax5.set_xlabel("azimuth")
ax5.set_ylabel("eta for e=90")
ax6 = fig.add_subplot(326)
ax6.set_xlabel("azimuth")
ax6.set_ylabel("eta")

#First test: For a few azimuth direction, varying elevation from 0 to 90
azimuth = np.array([0, 90, 180, 270]) * DtoR
# azimuth = np.linspace(0, np.pi/2 * 4, 10)
elevation = np.linspace(2.*DtoR, np.pi/2, N)

#Initiating lists
obs = [0]*N
# P_lon_list = np.zeros(N)
# P_lat_list = np.zeros(N)
# eta_igrf_list = np.zeros(N)
# eta_chaos_list = np.zeros(N)

azimuth_colors = ['red', 'blue', 'green', 'black']

#Looping over all positions, each get its own observation point object and the output values are plotted
for j, a in enumerate(azimuth):
	print("Fixed azimut: ", a*RtoD)
	for i, e in enumerate(elevation):
		obs[i] = ObservationPoint(A_lon, A_lat, e, a, h)
		obs[i].GetBatP()
		obs[i].GetEta(quadrant=False)
	# ax1.plot(elevation*RtoD, [o.P_lon*RtoD for o in obs], color = azimuth_colors[j], label = "azimuth: " + str(a*RtoD))
	# ax2.plot(elevation*RtoD, [o.P_lat*RtoD for o in obs], color = azimuth_colors[j], label = "azimuth: " + str(a*RtoD))

	# ax3.plot(elevation*RtoD, [o.eta_igrf*RtoD for o in obs], "*", color = azimuth_colors[j])
	# ax3.plot(elevation*RtoD, [o.eta_chaos*RtoD for o in obs], color = azimuth_colors[j],  label = "azimuth: " + str(a*RtoD))
	ax3.plot(elevation*RtoD, [o.Blos*RtoD for o in obs], color = (a/2./np.pi, 0, 1-a/2./np.pi))#,  label = "azimuth: " + str(a*RtoD))
# ax4.plot(elevation*RtoD, [o.AH_norm for o in obs], label = "Distance AH")

	# P_lon_list, P_lat_list, eta_igrf_list, eta_chaos_list, obs = np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N), [0] * N
	# obs = [0] * N

# Second test, this time the opposite: we fix a few values of elevation and we turn around at all azimuth
elevation = np.linspace(2.*DtoR, np.pi/2, 10)
# elevation  = np.linspace(0, 90, 5)*DtoR
# elevation  = np.array([2., 45., 90.])*DtoR
azimuth = np.linspace(0, 2*np.pi, N)
#same loop than before
for e in elevation:
	print("Fixed elevation:", e*RtoD)
	for i, a in enumerate(azimuth):
		obs[i] = ObservationPoint(A_lon, A_lat, e, a, h)
		obs[i].GetBatP()
		obs[i].GetEta(quadrant=False)
	# ax6.clear()
	# ax6.plot(azimuth*RtoD, [o.eta_igrf*RtoD for o in obs], "*", color = (2*e/np.pi, 0, 1 - 2*e/np.pi))
	ax6.plot(azimuth*RtoD, [o.eta_chaos*RtoD for o in obs], color = (2*e/np.pi, 0, 1-2*e/np.pi))
	ax4.plot(azimuth*RtoD, [o.Blos*RtoD for o in obs], color = (2*e/np.pi, 0, 1-2*e/np.pi), label="e:"+str(e*RtoD))

	if e <= 2.*DtoR:
		ax1.plot(azimuth*RtoD, [o.P_lon*RtoD for o in obs], label = "azimuth, elevation: " + str(a*RtoD) + str(e*RtoD))
		ax2.plot(azimuth*RtoD, [o.P_lat*RtoD for o in obs], label = "azimuth, elevation: " + str(a*RtoD) + str(e*RtoD))

		# ax5.plot(azimuth*RtoD, [o.eta_igrf*RtoD for o in obs], "*", color = "blue",  label = "igrf, e="+str(e*RtoD))
		ax5.plot(azimuth*RtoD, [o.eta_chaos*RtoD for o in obs], color = "blue", label = "Chaos, e="+str(e*RtoD))
		# ax6.plot(azimuth*RtoD, [o.P_lat*RtoD for o in obs], color = "red")
		# ax6.plot(azimuth*RtoD, [o.P_lon*RtoD for o in obs], color = "red")
	if e == np.pi/2:
		# ax5.plot(azimuth*RtoD, [o.eta_igrf*RtoD for o in obs], "*", color = "red",  label = "igrf, e="+str(e*RtoD))
		ax5.plot(azimuth*RtoD, [o.eta_chaos*RtoD for o in obs], color = "red", label = "chaos, e="+str(e*RtoD))

		# ax6.plot(azimuth*RtoD, [o.P_lat*RtoD for o in obs], color = "blue")
		# ax6.plot(azimuth*RtoD, [o.P_lon*RtoD for o in obs], color = "blue")
	# plt.draw()

for a in [ax1, ax2, ax3, ax4, ax5, ax6]:
	a.legend()


plt.show()
