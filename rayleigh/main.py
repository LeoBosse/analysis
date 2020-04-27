#!/usr/bin/python3
# -*-coding:utf-8 -*

import sys as sys
import numpy as np
import time as tm
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Circle
from matplotlib.patches import Ellipse
from matplotlib.patches import Arrow
# from matplotlib.lines import mlines

import osgeo.gdal as gdal
gdal.UseExceptions()  # not required, but a good idea

import imageio

from simulation import *
from observation import *
from rayleigh_utils import *
from sky_map import *
from ground_map import *
from atmosphere import *
from world import *
from input import *

def RunSimulation(file_name = "", show = True):

	all_time_start = tm.time()

	try:
		in_dict = ReadInputFile("./input_files/" + file_name + ".in")
		print("Correct input file in use:", file_name)
	except:
		in_dict = ReadInputFile("./input_files/RS_default.in")
		print("WARNING: Wrong or no input file specified, default in use.")

	#Init the rayleigh object
	simu = Simulation(in_dict)
	simu.ComputeAllMaps()
	simu.MakeSummaryPlot()


	# old = sys.stdout
	# f = open("log/systematic_results.csv", "a")
	# sys.stdout = f

	simu.PrintSystematicResults()

	# sys.stdout = old

	print("ALL TIME SINCE START:", tm.time() - all_time_start)

	if show:
		plt.show()

if __name__ == "__main__":

	#Get the input file from the command arguments
	arguments = sys.argv
	nb_args = len(arguments)

	if nb_args == 2:
		file_name = arguments[1]
	else:
		file_name = "./input_files/RS_default.in"
	RunSimulation(file_name)
