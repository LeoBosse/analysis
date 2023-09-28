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

import argparse

from tkinter import Tk
from tkinter.filedialog import askdirectory
from tkinter import filedialog 

from utils import *
from rotation import *
from bottle import *
from Mixer import *
from MagData import *
# import Geometry.Leo.src.geometry as geometry
import observation as observation
# import Geometry.Leo.src.rayleigh as rayleigh

from vendange_configuration import *



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


### Full path to where the data are stored (everything that comes before what is written in the argument)
data_path = global_configuration.data_path
# data_path = "/home/bossel/These/Analysis/data/"

### Manage the arguments. Should use only one argument: the name of the input file, without the path nor the extensions. e.g. 'ptcu20181115' will use the input file 'src/input_files/ptcu20181115.in'. Can add the data files to use as well, but I'm not sure this is working properly.
# arguments = sys.argv
# nb_args = len(arguments)
# if nb_args < 1:
# 	arguments.extend(input('Do not forget the arguments').split(' '))
# nb_args = len(arguments)
# print(arguments)

# show_plots, arguments 		= FindArgument("-s", arguments, default_value = True, getindex = 1)
# show_plots = bool(show_plots)
# # mag_data_file, arguments 	= FindArgument("-m", arguments, default_value = False, getindex = 1)
# # mag_data_file = bool(mag_data_file)
# comp_bottle, arguments 	= FindArgument("-b", arguments, default_value = False, getindex = 1)

# to_add_bottle, arguments 	= FindArgument("-a", arguments, default_value = None, getindex = 1)
# from_txt, arguments 		= FindArgument("-f", arguments, default_value = False, getindex = 0)


arg_parser = argparse.ArgumentParser(
                    #prog='Vendange',
                    description='Analysis of CRU data.',
					)

arg_parser.add_argument('bottle_name', help='The path from the data/ folder where to find the input.in file to use for this run. If given a folder, will look for a "input.in" file in it. If a file, must be an input file.')        # positional argument
arg_parser.add_argument('-s', '--show_plots',  action = 'store_false', help="If used, will not open nor show the graphs at the end of execution. Usefull for launching multiple runs via a bash script for exemple. All plots are saved no matter what.")
arg_parser.add_argument('-b', '--comp_bottle', default = None, help='Indicate an other input file path. This data will be plotted against the main one for easy comparison.')
arg_parser.add_argument('-f', '--from_file',  action = 'store_true', help='If used, will load the data from an existing processed data file. Reuse the data and do not compute everything again. Faster for re-plotting the same data again.')
arg_parser.add_argument('-a', '--append', default = None, help='An input.in file path to append to the main data. for exemple if several observations were run one after the other and you want to plot them all on the same time serie.')
arg_parser.add_argument('-l', '--lines', action='extend', default=None, type=str, help='The instrument channels you want to analyse, counting from 1. The analysis order will match the order given here. If not given, will automatically search and run all found bottles. e.g: 1234 or 231')
arg_parser.add_argument('-cl', '--comp_lines', action='extend', default=None, type=str, help='Similar to --lines, but applied to the comparison bottles.')

args = arg_parser.parse_args()

print(args)

bottle_name = args.bottle_name #arguments[1]

show_plots = args.show_plots
comp_bottle = args.comp_bottle
to_add_bottle = args.append
from_txt = args.from_file

if args.lines is not None:
	lines = [int(l) for l in args.lines]
	nb_lines = len(lines)
else:
	nb_lines = GetNbLines(bottle_name)
	lines = range(1, nb_lines+1)

if args.comp_lines is not None:
	comp_lines = [int(l) for l in args.comp_lines]
	nb_comp_lines = len(comp_lines)
else:
	nb_comp_lines = 0
	comp_lines = range(1, nb_comp_lines+1)

print(nb_lines, lines)

print(bottle_name)
instrument_name = GetInstrumentName(bottle_name)
print(instrument_name)
bottles = []
comp_bottles = []

if instrument_name in ["ptcu", "gdcu", "ptcu_v2", "carmen", "corbel"]:
	
	# print("nb_lines", nb_lines)
	# for l in [2]:
	# for l in [2, 1]:
	for l in lines:
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
		if to_add_bottle is not None:
			bottle = bottle + PTCUBottle(to_add_bottle, line = l, from_txt = from_txt)

		if bottle.valid:
			bottles.append(bottle)


elif instrument_name == "spp":
	bottle = SPPBottle(bottle_name)
	bottles.append(bottle)

if comp_bottle:
	# nb_lines = GetNbLines(comp_bottle)
	for l in comp_lines:
	# for l in range(1, len(bottles) + 1):
		print("##################################################################")
		print("##################################################################")
		print(comp_bottle)
		comp_bottles.append(PTCUBottle(comp_bottle, line = l, from_txt = from_txt))
		# comp_bottles[-1].SetTimeFromDateTime(bottles[0].DateTime(moment="start"))

if not from_txt:
	for b in bottles:
		b.PrintInfo()
		b.SaveTXT()
		b.SaveHDF5()

# if True :#mag_data_file:
# 	mag_data = MagData(bottles[0])
# 	if not mag_data.exist:
# 		mag_data = False
# else:
# 	mag_data = False

mixer = Mixer(bottles, comp_bottles=comp_bottles)

taster = Taster(mixer)


if show_plots:
	plt.show()
