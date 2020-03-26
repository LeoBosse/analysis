#!/usr/bin/python3
# -*-coding:utf-8 -*

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage{lmodern}']
from utils import *
from rotation import *
from bottle import *
from Mixer import *
from MagData import *
from scipy import signal
import sys
import os
import subprocess as subp

# import Geometry.Leo.src.geometry as geometry
import observation as observation
# import Geometry.Leo.src.rayleigh as rayleigh

### Manage the arguments. Should use only one argument: the name of the input file, without the path nor the extensions. e.g. 'ptcu20181115' will use the input file 'src/input_files/ptcu20181115.in'. Can add the data files to use as well, but I'm not sure this is working properly.
arguments = sys.argv
nb_args = len(arguments)
if nb_args < 1:
	arguments.extend(input('Do not forget the arguments').split(' '))
nb_args = len(arguments)

print(arguments)

show_plots, arguments 		= FindArgument("-s", arguments, default_value = True, getindex = 1)
show_plots = bool(show_plots)
mag_data_file, arguments 	= FindArgument("-m", arguments, default_value = True, getindex = 0)
mag_data_file = bool(mag_data_file)
another_bottle, arguments 	= FindArgument("-b", arguments, getindex = 1)


instrument_name = GetInstrumentName(arguments[1])
bottles = []

# ls = subp.check_output("ls /home/bossel/These/Analysis/data/" + arguments[1] + "/data*.csv")
if instrument_name != "spp":
	try:
		ls = [f for f in os.listdir("/home/bossel/These/Analysis/data/" + arguments[1]) if "data" in f]
	except:
		ls = [f for f in os.listdir("/home/bossel/These/Analysis/data/" + "/".join(arguments[1].split("/")[:-1])) if "data" in f]
	print("ls", ls)
	nb_lines = len(ls)
else:
	nb_lines = 1

if instrument_name in ["ptcu", "gdcu", "ptcu_v2", "carmen", "corbel"]:
	for l in range(1, nb_lines+1):
		# bottle = PTCUBottle(arguments[1], line = l)
		# try:
		# 	print("***************************************")
		# 	print("New Bottle")
		# 	bottle = PTCUBottle(arguments[1], line = l)
		# except:
		# 	print("New Bottle FAILED")
		# 	break
		bottle = PTCUBottle(arguments[1], line = l)
		if bottle.valid:
			bottles.append(bottle)


elif instrument_name == "spp":
	bottle = SPPBottle(arguments[1])
	bottles.append(bottle)

for b in bottles:
	b.PrintInfo()
	b.SaveTXT()

print("Optional arguments:", show_plots, mag_data_file)

if mag_data_file:
	mag_data = MagData(bottles[0])
	if not mag_data.exist:
		mag_data = False
else:
	mag_data = False
# mag_data = False

if another_bottle:
	comp_bottle = Bottle(another_bottle)
else:
	comp_bottle = False


Cru = Mixer(bottles, mag_data=mag_data, comp_bottle=comp_bottle)

if show_plots:
	plt.show()
