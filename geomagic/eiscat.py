#!/usr/bin/python3
# -*-coding:utf-8 -*

import numpy as np

DtoR = np.pi / 180.
RtoD = 1. / DtoR

### For a given observation with petit cru in Skibotn, compute the direction of observation of EISCAT radar near tromso.

a_pc = 164 * DtoR
e_pc = 90 * DtoR

h = 110. #km
print("h:", h)
AE = 51.66 #km

gamma = 2*np.pi - 60*DtoR - a_pc
print("gamma:", gamma * RtoD)

AH = h / np.tan(e_pc)
print("AH: ", AH)

EH = np.sqrt(AH**2 + AE**2 - 2 * AH * AE * np.cos(gamma))
print("EH: ", EH)

gamma2 = np.arcsin(AH * np.sin(gamma) / EH) * RtoD
print("gamma2: ", gamma2)

a_eiscat = (90 + 30 + gamma2)
e_eiscat = np.arctan(h / EH)

print("a_ptcu:", a_pc * RtoD)
print("e_ptcu: ", e_pc * RtoD)
print("a_eiscat:", a_eiscat)
print("e_eiscat: ", e_eiscat * RtoD)
