#!/usr/bin/python3
# -*-coding:utf-8 -*

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
#matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage{lmodern}']
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=False)
from scipy import signal
import sys
import os
import subprocess as subp

from utils import *
from rotation import *
from bottle import *
from Mixer import *
from MagData import *
# import Geometry.Leo.src.geometry as geometry
import observation as observation
# import Geometry.Leo.src.rayleigh as rayleigh

from vendange_configuration import *

### Full path to where the data are stored (everything that comes before what is written in the argument)
data_path = global_configuration.data_path
# data_path = "/home/bossel/These/Analysis/data/"

### Manage the arguments. Should use only one argument: the name of the input file, without the path nor the extensions. e.g. 'ptcu20181115' will use the input file 'src/input_files/ptcu20181115.in'. Can add the data files to use as well, but I'm not sure this is working properly.
arguments = sys.argv
nb_args = len(arguments)
if nb_args < 1:
	arguments.extend(input('Do not forget the arguments').split(' '))
nb_args = len(arguments)

print(arguments)

show_plots, arguments 		= FindArgument("-s", arguments, default_value = True, getindex = 1)
show_plots = bool(show_plots)
mag_data_file, arguments 	= FindArgument("-m", arguments, default_value = False, getindex = 1)
mag_data_file = bool(mag_data_file)
comp_bottle, arguments 	= FindArgument("-b", arguments, default_value = False, getindex = 1)

to_add_bottle, arguments 	= FindArgument("-a", arguments, default_value = None, getindex = 1)
from_txt, arguments 		= FindArgument("-f", arguments, default_value = False, getindex = 0)


bottle_name = arguments[1]
instrument_name = GetInstrumentName(bottle_name)
bottles = []
comp_bottles = []

def GetNbLines(bottle_name):
	"""Find the number of chanels of the instrument (the number of dataX.csv in the folder)"""
	instrument_name = GetInstrumentName(bottle_name)
	print(instrument_name, bottle_name)
	if instrument_name not in ["spp", 'carmen']:
		try:
			print("TRY", os.listdir(data_path + bottle_name))
			ls = [f for f in os.listdir(data_path + bottle_name) if "data" in f and ".csv" in f and len(f) in [9, 10]]
		except:
			print("EXCEPT", os.listdir(data_path + "/".join(bottle_name.split("/")[:-1])))
			ls = [f for f in os.listdir(data_path + "/".join(bottle_name.split("/")[:-1])) if "data" in f and ".csv" in f and len(f) in [9, 10]]
		print("ls", ls)
		nb_lines = len(ls)
	else:
		###SHOULD BE 2. USE ONLY IF COMPARING 2 INSTRUMENTS WITH DIFFERENT NB_LINES !!!
		# nb_lines = 1
		###BY DEFAULT SHOULD BE 2
		nb_lines = 1

	return nb_lines

if instrument_name in ["ptcu", "gdcu", "ptcu_v2", "carmen", "corbel"]:
	nb_lines = GetNbLines(bottle_name)
	print("nb_lines", nb_lines)
	# for l in [2]:
	for l in range(1, nb_lines+1):
		print("##################################################################")
		print("##################################################################")
		print(bottle_name)
		# bottle = PTCUBottle(arguments[1], line = l)
		# try:
		# 	print("***************************************")
		# 	print("New Bottle")
		# 	bottle = PTCUBottle(arguments[1], line = l)
		# except:
		# 	print("New Bottle FAILED")
		# 	break
		bottle = PTCUBottle(bottle_name, line = l, from_txt = from_txt)
		if to_add_bottle:
			bottle = bottle + PTCUBottle(to_add_bottle, line = l, from_txt = from_txt)

		if bottle.valid:
			bottles.append(bottle)


elif instrument_name == "spp":
	bottle = SPPBottle(bottle_name)
	bottles.append(bottle)

if comp_bottle:
	nb_lines = GetNbLines(comp_bottle)
	for l in [1]:
	# for l in range(1, len(bottles) + 1):
		print("##################################################################")
		print("##################################################################")
		print(comp_bottle)
		comp_bottles.append(PTCUBottle(comp_bottle, line = l, from_txt = from_txt))
		comp_bottles[-1].SetTimeFromDateTime(bottles[0].DateTime(moment="start"))

if not from_txt:
	for b in bottles:
		b.PrintInfo()
		b.SaveTXT()
		b.SaveHDF5()

if mag_data_file:
	mag_data = MagData(bottles[0])
	if not mag_data.exist:
		mag_data = False
else:
	mag_data = False
# mag_data = True

Cru = Mixer(bottles, mag_data=mag_data, comp_bottles=comp_bottles)

if show_plots:
	plt.show()
