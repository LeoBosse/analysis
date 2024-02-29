#!/usr/bin/python3
# -*-coding:utf-8 -*

import numpy as np
import scipy as sci
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
DtoR = 1./RtoD

def CreateArcimage(shape, limits, fct, center = (0, 0)):
    """
    Create a 2d array of the given shape, using the 2d function fct and the axis space defined by limits and center.
    """
    Y = np.linspace(limits[0], limits[1], shape[0])
    X = np.linspace(limits[0], limits[1], shape[1])
    xx, yy = np.meshgrid(X, Y)
    return fct(xx - center[0], yy - center[1])

def _rebin(arr, binX:int, binY:int, func='mean'):
    """
    Bins a 2d array using reshaping magic. Don't ask. Bin size defined by binX and binY. func in ['sum', 'mean', 'average'] to choose the binning mode.
    shape of the array must be a multiple of the bin size!
    """
    shape = (arr.shape[0]//binX, binX,
             arr.shape[1]//binY, binY)
    return getattr(getattr(arr.reshape(shape), func)(-1), func)(1) 

def BinArray(arr, binX, binY, func = 'mean'):
    """
    Bins a 2d array using reshaping magic. Don't ask. Bin size defined by binX and binY. func in ['sum', 'mean', 'average'] to choose the binning mode.
    """
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
    """
    Bins a 2d array using reshaping magic. Don't ask. Shape of the new array defined by new_shape. func in ['sum', 'mean', 'average'] to choose the binning mode.
    """
    return BinArray(arr, arr.shape[0]//new_shape[0], arr.shape[1]//new_shape[1], func = func)
          
### Choose the size of your images.
shape = 1000, 4000 #Nb pixels
limits = -2000, 2000 #axis limits. 


### Define fonction used to create your images. Can be anything but must only depend on the x and y coordinates of the image (x and y are in [limits], not pixel number!).
# fct = lambda x, y: stats.norm.pdf(x, 20, 600) * 650 / np.max(stats.norm.pdf(x, 20, 600)) + 30
fct  = lambda x, y: stats.maxwell.pdf(x, -1000, 600) * 650  / np.max(stats.maxwell.pdf(x, -1000, 600)) + 30
star = lambda x, y: stats.norm.pdf(x, 0, 10) * stats.norm.pdf(y, 0, 10) * 650000

### Define the shift in x and y for each image. (in limit space, not pixel number!)
centers = [(0, 0), (0, 0), (0, 0), (0, 0)]

### Define the spread of the gaussian filter used for defocus of each image. In pixels! Can be an int or a len 2 tuple for asymetric filters.
defocus_width = [10, 0, 5, 0] # in pixels


### Create the 4 camera images, defocus and bin them.
imgs = []
for i in range(4):
    imgs.append(CreateArcimage(shape, limits, fct, center = centers[i]))
    imgs[-1] += CreateArcimage(shape, limits, star, center = (500, 0))
    imgs[-1] += CreateArcimage(shape, limits, star, center = (1500, 00))
    
    imgs[-1] = sci.ndimage.gaussian_filter(imgs[-1], defocus_width[i])
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
    


### Compute polarisation
I = np.sum(imgs, axis = 0) / 2
Q = 2 * (imgs[0] - imgs[1]) / I
U = 2 * (imgs[2] - imgs[3]) / I

D = np.sqrt(Q**2 + U**2)
A = .5 * np.arctan2(U, Q)

print(f"Max DoLP: {np.max(D)*100}%")



### Plot graphs
fig, axs = plt.subplots(2, 2, sharex= True, sharey= True)
axs[0, 0].imshow(imgs[0])
axs[0, 1].imshow(imgs[1])
axs[1, 0].imshow(imgs[2])
axs[1, 1].imshow(imgs[3])


fig, axs = plt.subplots(3, 1, sharex= True, sharey= True)
imI = axs[0].imshow(I)
cbar = fig.colorbar(imI, ax = axs[0])
imD = axs[1].imshow(D*100)
cbar = fig.colorbar(imD, ax = axs[1])
imA = axs[2].imshow(A*RtoD)
cbar = fig.colorbar(imA, ax = axs[2])


f, axs = plt.subplots(4, sharex=True)

# axs[0].twinx().plot(imgs[0, :] - imgs[1, :], 'r')

for c in range(4):
    axs[0].plot(imgs[i][shape[0]//2, :])  
    
axs[1].plot(I[shape[0]//2, :])
axs[2].plot(D[shape[0]//2, :]*100)
axs[3].plot(A[shape[0]//2, :]*RtoD)

plt.show()
    
    
    


