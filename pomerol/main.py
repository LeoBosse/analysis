#!/bira-iasb/softs/21g/py39/bin/python3
# -*-coding:utf-8 -*

from mpi4py import MPI
mpi_comm = MPI.COMM_WORLD
mpi_rank = mpi_comm.Get_rank()
mpi_size = mpi_comm.Get_size()
mpi_name = mpi_comm.Get_name()

import sys as sys
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Circle
from matplotlib.patches import Ellipse
from matplotlib.patches import Arrow
# from matplotlib.lines import mlines

# import osgeo.gdal as gdal
# gdal.UseExceptions()  # not required, but a good idea

import imageio

from pomerol_configuration import *
from simulation import *
from observation import *
from rayleigh_utils import *
from sky_map import *
from ground_map import *
from atmosphere import *
from world import *
from input import *

def RunSimulation(file_name = "", show = True, output_result_file=None, header = True):

	datetime_start = dt.datetime.now()

	try:
		if file_name[-3:] != ".in":
			file_name += ".in"
		in_dict = ReadInputFile("./input_files/" + file_name)
		if mpi_rank == 0: print("Correct input file in use:", file_name)
	except:
		in_dict = ReadInputFile("./input_files/RS_default.in")
		if mpi_rank == 0: print("WARNING: Wrong or no input file specified, default in use.")

	if mpi_rank == 0:
		print("All input parameter from the input file:")
		for k, v in in_dict.items():
			if k[0] != "#":
				print(f"{k}: {v}")

	#Init the rayleigh object
	simu = Simulation(in_dict)
	simu.ComputeAllMaps()

	if int(in_dict["use_MS"]):
		simu.ComputeMultipleScattering()

	if mpi_rank == 0:

		old = sys.stdout
		if output_result_file:
			print(f"Print results in {output_result_file}")
			f = open(output_result_file, "a")
			sys.stdout = f

		simu.PrintSystematicResults(header = header)
		sys.stdout = old

		print("Simulation ran in:", dt.datetime.now() - datetime_start)

		simu.MakeSummaryPlot()

		print("Global timers recap:")
		for timer in GlobalTimer.fn_list:
			print(timer)

		if show:
			plt.show()


if __name__ == "__main__":
	#Get the input file from the command arguments
	arguments = sys.argv
	nb_args = len(arguments)

	# mpi_comm = MPI.COMM_WORLD
	# mpi_rank = mpi_comm.Get_rank()
	# mpi_size = mpi_comm.Get_size()
	# mpi_name = mpi_comm.Get_name()
	print(mpi_size, mpi_rank, mpi_name)

	if nb_args == 2:
		file_name = arguments[1]
	else:
		file_name = "./input_files/RS_default.in"

	# if mpi_rank == 0:
	RunSimulation(file_name)
