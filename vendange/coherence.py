#########################################################
# calculate coherence spectrum of pola series
#########################################################

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import imp
import numpy as np
from numpy.linalg import inv
import scipy.io as sio
import routines
import webgeodyn
from webgeodyn.processing import psd
import scipy.signal as signal
import math

#########################################################

#tmp = np.loadtxt('./br_a164_e45/20190307_skibotn_br_a164_e45_results.txt')
#tmp = np.loadtxt('./mr_a164_e45/20190307_skibotn_mr_a164_e45_results.txt')
#tmp = np.loadtxt('./vr_a164_e45/20190307_skibotn_vr_a164_e45_results.txt')

#tmp = np.loadtxt('./20190307_skibotn_mr_a164_e45_results.txt')
#tmp = np.loadtxt('./20190307_skibotn_vr_a164_e45_results.txt')
#tmp = np.loadtxt('./20190307_skibotn_br_a164_e45_results.txt')
tmp = np.loadtxt('./20190307_skibotn_rr_a164_e90_results.txt')

T = tmp[:,0]
dt = T[1]-T[0]
fs = 1/dt
NT = T.shape[0]

I = tmp[:,4]
D = tmp[:,5]
A = tmp[:,6]


fig=plt.figure(figsize=(17,9))

plt.subplot(3,1,1)
plt.plot(T, I,'r-')
plt.ylabel(('I (mV)'))
plt.grid()

plt.subplot(3,1,2)
plt.plot(T, D,'g-')
plt.xlabel(('time (units = min)'))
plt.ylabel(('DoLP (%)'))
plt.grid()

plt.subplot(3,1,3)
plt.plot(T, A*180/np.pi,'b-')
plt.xlabel(('time (units = min)'))
plt.ylabel(('AoLP (degrees)'))
plt.grid()

plt.show()



NN = 2**9
f, Pid = signal.csd(I, D, fs, detrend='linear',nperseg=NN)
f, Pia = signal.csd(I, A, fs, detrend='linear', nperseg=NN)
f, Pad = signal.csd(A, D, fs, detrend='linear', nperseg=NN)
f, Pii = signal.csd(I, I, fs, detrend='linear', nperseg=NN)
f, Pdd = signal.csd(D, D, fs, detrend='linear', nperseg=NN)
f, Paa = signal.csd(A, A, fs, detrend='linear', nperseg=NN)
Nf = f.shape[0]

Cid = np.abs(Pid)**2/(Pii*Pdd)
Cad = np.abs(Pad)**2/(Paa*Pdd)
Cia = np.abs(Pia)**2/(Pii*Paa)

ss = np.imag(Pid[:])/np.abs(Pid)
cc = np.real(Pid[:])/np.abs(Pid)
angle_id = np.arctan2(ss,cc)

ss = np.imag(Pia[:])/np.abs(Pia)
cc = np.real(Pia[:])/np.abs(Pia)
angle_ia = np.arctan2(ss,cc)

ss = np.imag(Pad[:])/np.abs(Pad)
cc = np.real(Pad[:])/np.abs(Pad)
angle_ad = np.arctan2(ss,cc)

#f, Cid2 = signal.coherence(I, D, fs, nperseg=1024)

period = 1./f
period_min = 1./f/60
period_hour = period/3600

fig=plt.figure(figsize=(9,9))

plt.subplot(2,1,1)
plt.plot(np.log10(period_min), Cid,'ro-',label='Cid')
plt.plot(np.log10(period_min), Cad,'go-',label='Cad')
plt.plot(np.log10(period_min), Cia,'bo-',label='Cia')
plt.xlabel(('period (log10, units = min)'))
plt.ylabel(('coherence'))
plt.legend()
plt.ylim((0,1))
plt.xlim((-1,1))
plt.grid()

plt.subplot(2,1,2)
plt.plot(np.log10(period_min), angle_id*180/np.pi,'ro-',label='Cid')
plt.plot(np.log10(period_min), angle_ad*180/np.pi,'go-',label='Cad')
plt.plot(np.log10(period_min), angle_ia*180/np.pi,'bo-',label='Cia')
plt.xlabel(('period (log10, units = min)'))
plt.ylabel(('phase'))
plt.legend()
plt.ylim((-180,180))
plt.xlim((-1,1))
plt.grid()

plt.show()
