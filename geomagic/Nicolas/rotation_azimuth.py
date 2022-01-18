#########################################################
# estimate apparent angle of the local magnetic field in the SPP frame
#	 during one rotation of the polarimeter
#########################################################

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import imp
import numpy as np
from numpy.linalg import inv
import scipy.io as sio
import routines

################################################
# points:
#	O: center of the Earth
#	A: measurement point, on ground: altitude "alt", latitude "lat", longitude "lon"
#	H: target point of altitude "h", elevation angle "eps", azimuth "alpha" wrt south
#	P: intersection of OH and the Earth's surface
#
################################################
# frames (see Laundal & Richmond, Space Science Review, 2016) :
# 	"enu": East, North, Up (depends on local position)
#	"ecef": east-centered-earth-fixed
#	"SPP": orthonormal frame (a,b,c)
#	 	a = unit vector along AH
#	 	b = unit vector perp. to a, in the horizontal plane
#	 	c = axb
################################################
# entries:

h = 110 			# altitude of point H (km)
lat = 44.83 * np.pi/180 	# lalitude of instrument (point A)
lon = 5.76 * np.pi/180		# longitude of instrument (point A)
alt = 0.9			# altitude of point A
R0 = 6371.2			# Earth's radius (km)
eps = 45 * np.pi/180		# elevation angle, 0 < eps < pi/2

################################################
# spherical coordinates for A:

R_A = R0 + alt
theta_A = np.pi/2-lat
phi_A = lon

R_H = R0 + h

# vector AH in enu(A)

# estimate of distance AH (independent of the beamline azimuth) :
#	Al-Kashi theorem in triangle OAH
#	... AH^2 + 2*R_A*sin(eps)*AH + (alt**2+2*R0*alt-h**2-2*R0*h) = 0

aa = 1.0
bb = 2*R_A*np.sin(eps)
cc = alt**2+2*R0*alt-h**2-2*R0*h
delta = bb**2-4*aa*cc
AH = (-bb+np.sqrt(delta))/(2*aa)

################################################
# loop over angle alpha (between South and BeamLine, in horizontal plane):

dangle = 10
alpha_all = np.linspace(dangle,360,360/dangle)
alpha_all = alpha_all * np.pi/180
Nangle = alpha_all.shape[0]

position_H_all = np.zeros((Nangle,2))
a_enu_P_all = np.zeros((Nangle,3))
b_enu_P_all = np.zeros((Nangle,3))
c_enu_P_all = np.zeros((Nangle,3))
a_enu_A_all = np.zeros((Nangle,3))

for j in range(Nangle):
	alpha = alpha_all[j]
	# unit vector a_enu_A along (AH) :
	# NB: alpha clockwise, 0 towards South !
	a_enu_A = np.zeros((3,1))
	a_enu_A[0] = -np.cos(eps)*np.sin(alpha)	# east
	a_enu_A[1] = -np.cos(eps)*np.cos(alpha)	# north
	a_enu_A[2] = np.sin(eps)		# up
	AH_enu_A = AH * a_enu_A
	# OA vector in ecef frame :
	OA_ecef = np.zeros((3,1))
	OA_ecef[0] = R_A*np.sin(theta_A)*np.cos(phi_A)
	OA_ecef[1] = R_A*np.sin(theta_A)*np.sin(phi_A)
	OA_ecef[2] = R_A*np.cos(theta_A)
	# rotation matrix from ecef frame to enu(A)
	M_ecef2enuA = routines.rotation_ecef2enu(theta_A,phi_A)
	M_enuA2ecef = M_ecef2enuA.T
	# OH vector in ecef frame :
	AH_ecef = np.dot(M_enuA2ecef,AH_enu_A)
	OH_ecef = OA_ecef + AH_ecef
	# estimate of (theta_H,phi_H)
	OH = np.sqrt(OH_ecef[0]**2+OH_ecef[1]**2+OH_ecef[2]**2)
	theta_H = np.arccos(OH_ecef[2]/OH)
	phi_H = np.arctan2(OH_ecef[1],OH_ecef[0])
	position_H_all[j,0] = theta_H
	position_H_all[j,1] = phi_H
	# rotation matrix from ecef frame to enu(P)
	theta_P = theta_H
	phi_P = phi_H
	M_ecef2enuP = routines.rotation_ecef2enu(theta_P,phi_P)
	M_enuP2ecef = M_ecef2enuP.T
	#### vector a in frame enu(H) :
	a_ecef = np.dot(M_enuA2ecef,a_enu_A)
	a_enu_P = np.dot(M_ecef2enuP,a_ecef)
	a_enu_P_all[j,:] = a_enu_P[:,0]
	a_enu_A_all[j,:] = a_enu_A[:,0]
	# in the frame of the SPP : (a,b,c)
	#	b perp to a and in horizontal plane
	#	c = axb
	b_enu_A = np.zeros((3,1))
	b_enu_A[2] = 0						# u
	b_enu_A[1] = np.sin(alpha)	# 1/np.sqrt(1+(a_enu_A[1]/a_enu_A[0])**2) 	# n
	b_enu_A[0] = -np.cos(alpha)	# -a_enu_A[1]/a_enu_A[0]*b_enu_A[1]		# e
	c_enu_A = routines.cross_product(a_enu_A,b_enu_A)
	# rotate (a,b,c) from ecef frame to enu(A)
	a_ecef = np.dot(M_enuA2ecef,a_enu_A)
	b_ecef = np.dot(M_enuA2ecef,b_enu_A)
	c_ecef = np.dot(M_enuA2ecef,c_enu_A)
	# rotate (a,b,c) from ecef frame to enu(P)
	a_enu_P = np.dot(M_ecef2enuP,a_ecef)
	b_enu_P = np.dot(M_ecef2enuP,b_ecef)
	c_enu_P = np.dot(M_ecef2enuP,c_ecef)
	b_enu_P_all[j,:] = b_enu_P[:,0]
	c_enu_P_all[j,:] = c_enu_P[:,0]


plt.figure()
plt.subplot(3,1,1)
plt.plot(alpha_all*180/np.pi,a_enu_P_all[:,0],'-r',label='a(e)')
plt.plot(alpha_all*180/np.pi,a_enu_P_all[:,1],'-g',label='a(n)')
plt.plot(alpha_all*180/np.pi,a_enu_P_all[:,2],'-b',label='a(u)')
plt.grid()
plt.ylabel('a')
plt.subplot(3,1,2)
plt.plot(alpha_all*180/np.pi,b_enu_P_all[:,0],'-r',label='b(e)')
plt.plot(alpha_all*180/np.pi,b_enu_P_all[:,1],'-g',label='b(n)')
plt.plot(alpha_all*180/np.pi,b_enu_P_all[:,2],'-b',label='b(u)')
plt.grid()
plt.ylabel('b')
plt.subplot(3,1,3)
plt.plot(alpha_all*180/np.pi,c_enu_P_all[:,0],'-r',label='c(e)')
plt.plot(alpha_all*180/np.pi,c_enu_P_all[:,1],'-g',label='c(n)')
plt.plot(alpha_all*180/np.pi,c_enu_P_all[:,2],'-b',label='c(u)')
plt.grid()
plt.legend()
# plt.xlabel('alpha',fontsize=myfontsize)
plt.ylabel('c')
plt.show()

'''
plt.plot(alpha_all*180/np.pi,a_enu_P_all)
plt.plot(alpha_all*180/np.pi,a_enu_A_all,'--')
plt.show()

plt.plot(position_H_all[:,0]*180/np.pi,position_H_all[:,1]*180/np.pi)
plt.plot(lat*180/np.pi,lon*180/np.pi,'+')
plt.show()
'''
np.savetxt('theta_phi_H.dat',position_H_all, fmt='%1.8e')

'''
run here fortran code to calculate field model predictions at point H:

for CHAOS-6:
compile with 'gfortran pred_chaos6.f -o chaos.out'
run with
./chaos.out <<+
2018.5
110
+

for IGRF-12:
compile with 'gfortran igrf12_loopH.f -o igrf.out'
run with
./igrf.out <<+

2
1
2
2018.5
6481.2
+
'''

Bxyz1 = np.loadtxt('Bxyz_H_igrf.dat')
Bxyz2 = np.loadtxt('Bxyz_H_chaos6.dat')

nB1_enu_P_all = np.zeros((Nangle,3))
nB2_enu_P_all = np.zeros((Nangle,3))
angleU2B1_all = np.zeros((Nangle))
angleU2B2_all = np.zeros((Nangle))
eta1_all = np.zeros((Nangle))
eta2_all = np.zeros((Nangle))

for j in range(Nangle):
	# unit vectors b and c in enu_P
	b_enu_P = np.zeros((3,1))
	b_enu_P[:,0] = b_enu_P_all[j,:]
	c_enu_P = np.zeros((3,1))
	c_enu_P[:,0] = c_enu_P_all[j,:]
	#--- IGRF
	B_enu_P = np.zeros((3,1))
	B_enu_P[0,0] = Bxyz1[j,3]
	B_enu_P[1,0] = Bxyz1[j,2]
	B_enu_P[2,0] = -Bxyz1[j,4]
	norm_B_enu_P = np.linalg.norm(B_enu_P)
	nB_enu_P = B_enu_P/norm_B_enu_P # normalized vector B
	nB1_enu_P_all[j,:] = nB_enu_P[:,0]
	tmp = routines.dot_product(nB_enu_P,a_enu_P_all[j,:])
	angleU2B = np.arccos(tmp) # careful between 0 and pi !!
	angleU2B1_all[j] = angleU2B
	# angle eta between B and b in plane (b,c)
	Bb = routines.dot_product(nB_enu_P,b_enu_P)
	Bc = routines.dot_product(nB_enu_P,c_enu_P)
	nBb = Bb/np.sqrt(Bb**2+Bc**2)
	nBc = Bc/np.sqrt(Bb**2+Bc**2)
	cos_eta = nBb
	sin_eta = nBc
	if (sin_eta>=0):
		eta1_all[j] = np.arccos(cos_eta)
	#--- CHAOS-6
	B_enu_P = np.zeros((3,1))
	B_enu_P[0,0] = Bxyz2[j,3]
	B_enu_P[1,0] = Bxyz2[j,2]
	B_enu_P[2,0] = -Bxyz2[j,4]
	norm_B_enu_P = np.linalg.norm(B_enu_P)
	nB_enu_P = B_enu_P/norm_B_enu_P # normalized vector B
	nB2_enu_P_all[j,:] = nB_enu_P[:,0]
	tmp = routines.dot_product(nB_enu_P,a_enu_P_all[j,:])
	angleU2B = np.arccos(tmp) # careful between 0 and pi !!
	angleU2B2_all[j] = angleU2B
	# angle eta between B and b in plane (b,c)
	Bb = routines.dot_product(nB_enu_P,b_enu_P)
	Bc = routines.dot_product(nB_enu_P,c_enu_P)
	nBb = Bb/np.sqrt(Bb**2+Bc**2)
	nBc = Bc/np.sqrt(Bb**2+Bc**2)
	cos_eta = nBb
	sin_eta = nBc
	if (sin_eta>=0):
		eta2_all[j] = np.arccos(cos_eta)

# careful: alpha defined wrt north in draft paper !!
azimuth = (alpha_all+np.pi)-2*np.pi

fig=plt.figure(figsize=(16,7))
myfontsize=16

plt.plot(azimuth*180/np.pi,eta1_all*180/np.pi,'o-',linewidth=2,label='eta IGRF')
plt.plot(azimuth*180/np.pi,eta2_all*180/np.pi,'+--',linewidth=2,label='eta CHAOS')
plt.ylabel('eta (degree)',fontsize=myfontsize)
plt.xlabel('azimuth (degree)',fontsize=myfontsize)
plt.legend(fontsize=16)
plt.xticks(fontsize=myfontsize)#, rotation=90)
plt.yticks(fontsize=myfontsize)#, rotation=90)
plt.grid()

plt.show()

#############################################

fig=plt.figure(figsize=(16,7))
myfontsize=16

plt.subplot(1,2,1)

plt.plot(azimuth*180/np.pi,nB1_enu_P_all[:,0],'b',label='B e IGRF',linewidth=2)
plt.plot(azimuth*180/np.pi,nB1_enu_P_all[:,1],'g',label='B n IGRF',linewidth=2)
plt.plot(azimuth*180/np.pi,nB1_enu_P_all[:,2],'r',label='B u IGRF',linewidth=2)
'''
plt.plot(azimuth*180/np.pi,nB2_enu_P_all[:,0],'b',label='B e CHAOS6',linewidth=1)
plt.plot(azimuth*180/np.pi,nB2_enu_P_all[:,1],'g',label='B n CHAOS6',linewidth=1)
plt.plot(azimuth*180/np.pi,nB2_enu_P_all[:,2],'r',label='B u CHAOS6',linewidth=1)
'''
plt.plot(azimuth*180/np.pi,a_enu_P_all[:,0],'b--',label='AH e',linewidth=2)
plt.plot(azimuth*180/np.pi,a_enu_P_all[:,1],'g--',label='AH n',linewidth=2)
plt.plot(azimuth*180/np.pi,a_enu_P_all[:,2],'r--',label='AH u',linewidth=2)
plt.grid()
plt.xlabel('azimuth',fontsize=myfontsize)
plt.legend(fontsize=16)
plt.xticks(fontsize=myfontsize)#, rotation=90)
plt.yticks(fontsize=myfontsize)#, rotation=90)
plt.ylabel('unit vector components', fontsize=myfontsize)
#plt.ylim((-9,1))
#plt.xlim((-3.,0.))

plt.subplot(1,2,2)

plt.plot(azimuth*180/np.pi,angleU2B1_all*180/np.pi,linewidth=2,label='B IGRF')
plt.plot(azimuth*180/np.pi,angleU2B2_all*180/np.pi,linewidth=1,label='B CHAOS6')
plt.xlabel('azimuth',fontsize=myfontsize)
plt.ylabel('angle (B,AH)',fontsize=myfontsize)
plt.legend(fontsize=16)
plt.xticks(fontsize=myfontsize)#, rotation=90)
plt.yticks(fontsize=myfontsize)#, rotation=90)
plt.grid()

plt.show()
