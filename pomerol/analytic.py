#!/usr/bin/python3
# -*-coding:utf-8 -*

import sys as sys
import numpy as np
import scipy.integrate as int
import matplotlib
import matplotlib.pyplot as plt

DtoR = np.pi / 180.
RtoD = 1. / DtoR


wl = 0.3914 #micro m

# ai = 0	* DtoR
ai_list = np.linspace(0, 360, 37, endpoint=True)	* DtoR
# ai_list = np.array([0])	* DtoR

ei = 45 	* DtoR
alpha = 1 * DtoR
A_A = 20e-10 # km**2

ae = 180 	* DtoR
ee = 0  	* DtoR
AE = 5 #km
LE = 100 #nW/m**2/sr
A_E = 1 #m**2

min_alt = 0
min_range = min_alt / np.sin(ei)
max_alt = 100
max_range = max_alt / np.sin(ei)


F_list = []
Fpola_list = []
FNpola_list = []
D_list = []
# RAE_list = np.ones(len(ai_list)) * 90 * DtoR
# RAE_list = np.append(np.linspace(135, 45, 19, endpoint=True), np.linspace(50, 135, 18, endpoint=True)) * DtoR
RAE_list = np.arccos(np.cos(ei) * np.sin(ai_list) * np.sin(ae) + np.cos(ei) * np.cos(ai_list) * np.cos(ae))
print(RAE_list*RtoD)


### Loop over PTCU orientations
for ai, RAE in zip(ai_list, RAE_list):
	F_list.append([])
	Fpola_list.append([])
	FNpola_list.append([])
	D_list.append([])
	# RAE = 45 	* DtoR

	ER = lambda r: np.sqrt(r**2 + AE**2 - 2 * r * AE * np.cos(RAE))

	theta = lambda r: np.arcsin(AE * np.sin(RAE) / ER(r))

	V = lambda r, alpha: (2./3.) * np.pi * (1 - np.cos(alpha)) * r ** 2

	A = 7.68246 * 10 ** (-4) #for km-1
	B = 3.55212
	C = 1.35579
	D = 0.11563
	beta0 = A * wl ** (- (B + C * wl + D / wl))
	# print("beta0", beta0)

	A = 6.49997 * 10 ** (-3)
	B = 3.55212
	C = 1.35579
	D = 0.11563
	tau0 = A * wl ** (- (B + C * wl + D / wl))

	g = 9.80665
	M = 0.02896968
	T0 = 288.16
	R0 = 8.314462618
	H = - g * M / T0 / R0 * np.sin(ei)
	# print("H", H)

	beta = lambda beta0, r: beta0 * np.exp(H * r)
	tau = lambda r: tau0 * np.exp(H * r) * (ER(r) / r + 1) / np.sin(ei)

	# P = lambda theta: 3./4. * (1 + np.cos(theta)**2)

	K0 = 3 * beta0 * A_E * A_A * LE * (1 - np.cos(alpha)) / 8
	# print(1 - np.cos(alpha))
	# print("K0", K0)
	K2 = (AE * np.sin(RAE)) ** 2
	# print("K2", K2)

	dF = lambda r: K0 * np.exp(H * r) * (2 - K2 / ER(r)**2) * np.exp(-tau(r)) / ER(r)**2
	DoLP = lambda r: np.sin(theta(r))**2 / (1 + np.cos(theta(r))**2)
	dFpola = lambda r: dF(r) * DoLP(r)
	dFNpola = lambda r: dF(r) * (1 - DoLP(r))

	### Loop over PTCU range
	max_range_list = np.linspace(min_range, max_range, 100)[1:]
	max_range_list = [max_range_list[-1]]
	for M in max_range_list:
		F = int.quad(dF, min_range, M, points=[0])
		Fpola = int.quad(dFpola, min_range, M, points=[0])
		FNpola = int.quad(dFNpola, min_range, M, points=[0])
		# print("F (nW):", F[0])
		# print("DoLP (%):", Fpola[0] / F[0] * 100)
		if abs(Fpola[0] + FNpola[0] - F[0]) > F[0] / 1000.:
			print(Fpola[0] + FNpola[0] - F[0], abs(Fpola[0] + FNpola[0] - F[0]) < F[0] / 10.)

		F_list[-1].append(F[0])
		Fpola_list[-1].append(Fpola[0])
		FNpola_list[-1].append(FNpola[0])

		D_list[-1].append(Fpola[0] / F[0] * 100)


print("F", F_list)
print("Fpola", Fpola_list)
print("Funpola", FNpola_list)
print("DoLP", D_list)


font = {'weight' : 'bold',
'size'   : 24}
matplotlib.rc('font', **font)
### Plot almucantar for every max range
f, axs = plt.subplots(2, sharex=True, figsize=(16, 8))
axs[0].plot(ai_list * RtoD, F_list, label = str(ai))
axs[0].plot(ai_list * RtoD, Fpola_list, label = str(ai), color="g")
axs[0].plot(ai_list * RtoD, FNpola_list, label = str(ai), color = "r")
axs[1].plot(ai_list * RtoD, D_list, label = str(ai))
axs[0].set_ylabel("Flux (nW)")
axs[1].set_xlabel("Azimut (°)")
axs[1].set_ylabel("DoLP (\%)")
f.subplots_adjust(hspace=0)
# plt.show()


### Plot max range variations for every PTCU orientations
f, axs = plt.subplots(2, sharex=True)
for i, ai in enumerate(ai_list):
	axs[0].plot(max_range_list, F_list[i], label = str(ai))
	axs[1].plot(max_range_list, D_list[i], label = str(ai))
f.subplots_adjust(hspace=0)
plt.show()
