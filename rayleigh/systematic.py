#!/usr/bin/python3
# -*-coding:utf-8 -*


from main import *
from rayleigh_input import *
import pandas as pd

import sys


def GetData():
	# data = np.genfromtxt("./log/systematic_results.csv", delimiter=",", comments = "t", names=['t', 'a_pc', 'e_pc', 'src_dist', 'dlos', 'min_alt', 'max_alt', 'I0', 'DoLP', 'AoRD'])
	data = pd.read_csv("./log/systematic_results.csv", sep=",", comment = "t", names=['t', 'a_pc', 'e_pc', 'src_dist', 'dlos', 'min_alt', 'max_alt', 'I0', 'DoLP', 'AoRD'])

	# nb_data = data.index.stop
	#
	# wl = [360] * , 557.7, 427.8, 391.4] * int(nb_data / 4)
	# if nb_data % 4 == 1: wl.append([360])
	# elif nb_data % 4 == 2: wl.extend([360, 557.7])
	# elif nb_data % 4 == 3: wl.extend([360, 557.7, 427.7])
	#
	# data['wl'] = pd.Series(wl, index=data.index)

	return data

def PlotData(data, x_axis, y_axis, **kwargs):
	if kwargs:
		x_ax, y_ax = [], []
		for i in range(len(data)):
			valid = True
			for k, v in kwargs.items():
				if data[k][i] != v:
					valid = False
					break
			if valid:
				print(i, data.loc[[i]])
				x_ax.append(data[x_axis][i])
				y_ax.append(data[y_axis][i])
	else:
		x_ax, y_ax = data[x_axis], data[y_axis]
	print(x_ax, y_ax)
	plt.plot(x_ax, y_ax)


def RunSystematicSimulation():
	azimuths = "all"
	elevations = "all"
	src_dist = [5]
	dlos = [1, 5]
	alt_min = [0, 0.5, 1, 5, 10, 50]
	alt_max = [0.5, 1, 5, 10, 50, 100]
	wavelengths = [630, 557.7, 427.8, 391.4]

	for d in src_dist:
		for dl in dlos:
			for am in alt_min:
				for aM in [a for a in alt_max if a > am + dl]:
					for wl in wavelengths:
						in_dict = {"azimuts": azimuths,
								"elevations": elevations,
								"point_src_dist": d,
								"resolution_along_los": dl,
								"RS_min_altitude": am,
								"RS_max_altitude": aM,
								"wavelength": wl
								}

						r_input = RayleighInput()

						r_input.Update(in_dict)

						file_name = r_input.WriteInputFile(file_name = "systemic/systemic", **in_dict)

						print(file_name)

						log_file_name = "./log/" + file_name + ".out"
						log_file = open(log_file_name, "w")

						terminal = sys.stdout
						sys.stdout = log_file
						sys.stderr = log_file

						RunSimulation(file_name, show = False)

						sys.stdout = terminal
						sys.stderr = terminal


if __name__ == "__main__":
	RunSystematicSimulation()
