#!/usr/bin/python3
# -*-coding:utf-8 -*

import numpy as np
import matplotlib.pyplot as plt
from utils import *
from scipy import signal
import sys as sys
import os
from subprocess import call

import h5py as h5

import datetime as dt
import scipy.io
# import skimage as ski
from PIL import Image, ImageEnhance

from plip import *



def IQUfromIDA(I, D, A):
    Q = D * np.cos(2 * A) / 100.
    U = D * np.sin(2 * A) / 100.
    return I, Q, U
    
def MZPfromIQU(I, Q, U):
    Z = (I + Q) / 2.
    tmp = I/2. - Q/4. 
    M = tmp - U * np.sqrt(3)/4
    P = tmp + U * np.sqrt(3)/4
    return M, Z, P


def BetterContrast(imgs, m, M):
    m, M = np.percentile(imgs, (m, M))
    
    a = 1./(M-m)
    b = -a * m
    
    return imgs * a + b 
    
    

# data = PLIP("data_vh-+.h5")
data = PLIP("data_SUB_0.h5")

# MZP pola is similar to RGB color
# IQU pola is similar to Y, Cb, Cr color (similar to HSL) Y~brightness, Cb and Cr are color in [-1, 1]

I, Q, U = IQUfromIDA(data.I_fine, data.D_fine, data.A_fine)
I, Q, U = I[1:], Q[1:], U[1:]

contrast = 4000
M, Z, P = MZPfromIQU(I / contrast, Q, U)
MZP  = np.stack((M, Z, P), axis=3)

MZP = BetterContrast(MZP, 20, 92)

nb_imgs = I.shape[0]
nrows, ncols = int(np.floor(np.sqrt(nb_imgs))), int(np.ceil(np.sqrt(nb_imgs)))
fig_size_coef = 400
fig, axs = plt.subplots(nrows, ncols, figsize=(2 * fig_size_coef, 2 * fig_size_coef), sharex=True, sharey=True)

# I_vmax = np.percentile(I, 98)

for irow in range(nrows):
    for icol in range(ncols):
        i = icol + irow * ncols
        if i >= nb_imgs:
            break
        # im_I = axs[irow, icol].imshow(MZP[i])
        im_I = axs[irow, icol].imshow(I[i])


# fig.tight_layout()
# fig.subplots_adjust(right=0.9)
# cbar_ax_I = fig_I.add_axes([0.90, 0.05, 0.02, 0.9])
# fig_I.colorbar(im_I, cax=cbar_ax_I, orientation='vertical')
plt.show()