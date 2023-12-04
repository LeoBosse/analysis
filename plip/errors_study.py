#!/usr/bin/python3
# -*-coding:utf-8 -*

import numpy as np
from scipy.interpolate import interpn
import scipy.stats as stats
from scipy.signal import fftconvolve
import matplotlib as mpl
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys
import os
from subprocess import call
import datetime as dt

# import radis
import astropy.io.fits as fits

from plip import *

data = PLIP("data_SUB_0.h5")

RtoD = 180./np.pi

def CreateArcimage(shape, limits, fct, center = (0, 0)):
    Y = np.linspace(limits[0], limits[1], shape[0])
    X = np.linspace(limits[0], limits[1], shape[1])
    xx, yy = np.meshgrid(X, Y)
    return fct(xx - center[0], yy - center[1])

def _rebin(arr, binX:int, binY:int, func='mean'):
    shape = (arr.shape[0]//binX, binX,
             arr.shape[1]//binY, binY)
    return getattr(getattr(arr.reshape(shape), func)(-1), func)(1) 

def BinArray(arr, binX, binY, func = 'mean'):
    shape = arr.shape
    if (shape[0]%binX != 0):
        rest = shape[0]%binX
        arr1 = _rebin(arr[:shape[0]-rest], binX, 1, func=func)
        arr2 = _rebin(arr[shape[0]-rest:], rest, 1, func=func)
        arr = np.concatenate((arr1,arr2), axis=0)
        
    if (shape[1]%binY != 0):
        rest = shape[1]%binY
        arr1 = _rebin(arr[:,:shape[1]-rest], 1, binY, func=func)
        arr2 = _rebin(arr[:,shape[1]-rest:], 1, rest, func=func)
        arr = np.concatenate((arr1,arr2), axis=1)
    
    if (shape[0]%binX == 0) and (shape[1]%binY == 0):
        arr = _rebin(arr, binX, binY, func=func)
    
    return arr
  
def ShapeArray(arr, new_shape, func = 'mean'):
    return BinArray(arr, arr.shape[0]//new_shape[0], arr.shape[1]//new_shape[1], func = func)
          

shape = 1000, 4000
limits = -2000, 2000
fct = lambda x, y: stats.norm.pdf(x, 20, 600) * 650 / np.max(stats.norm.pdf(x, 20, 600)) + 30
# fct = lambda x, y: stats.maxwell.pdf(x, -1000, 60) * 650  / np.max(stats.maxwell.pdf(x, 0, 60)) + 30
centers = [(3, 0), (0, 0), (0, 0), (0, 0)]
imgs = []
for i in range(4):
    imgs.append(CreateArcimage(shape, limits, fct, center = centers[i]))
    imgs[-1] = BinArray(imgs[-1], 1, 1, 'mean')


# imgs[0] *= 1.02
# imgs[2] *= 1.01
# imgs[3] *= 1.03

# shape = data.I_fine[16].shape
# imgs = []
# shift = 1
# for i in range(4):
#     imgs.append(data.I_fine[14])
#     if i == 0 and shift != 0:
#         imgs[-1] = np.append(imgs[-1][:, shift:], np.zeros((imgs[-1].shape[0], shift)), axis = 1)
    
#     imgs[-1] = BinArray(imgs[-1], 20, 20, 'mean')
    

I = np.sum(imgs, axis = 0) / 2
Q = 2 * (imgs[0] - imgs[1]) / I
U = 2 * (imgs[2] - imgs[3]) / I

D = np.sqrt(Q**2 + U**2)
A = .5 * np.arctan2(U, Q)

print(f"Max DoLP: {np.max(D)*100}%")


fig, axs = plt.subplots(1, 3, sharex= True, sharey= True)
imI = axs[0].imshow(I)
cbar = fig.colorbar(imI, ax = axs[0])
imD = axs[1].imshow(D*100)
cbar = fig.colorbar(imD, ax = axs[1])
imA = axs[2].imshow(A*RtoD)
cbar = fig.colorbar(imA, ax = axs[2])


f, axs = plt.subplots(3, sharex=True)

# axs[0].twinx().plot(imgs[0, :] - imgs[1, :], 'r')

axs[0].plot(I[0, :])
axs[1].plot(D[0, :]*100)
axs[2].plot(A[0, :]*RtoD)

plt.show()
    
    
    


