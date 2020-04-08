#!/usr/bin/python3
# -*-coding:utf-8 -*


from main import *
from rayleigh_input import *
import pandas as pd

import sys


def GetData():
	# data = np.genfromtxt("./log/systematic_results.csv", delimiter=",", comments = "t", names=['t', 'a_pc', 'e_pc', 'src_dist', 'dlos', 'min_alt', 'max_alt', 'I0', 'DoLP', 'AoRD'])
	data = pd.read_csv("./log/systematic_results.csv", sep=",", comment = "t", names=['t', 'a_pc', 'e_pc', 'src_dist', 'dlos', 'min_alt', 'max_alt', "wl", 'I0', 'DoLP', 'AoRD'])

	data = data.drop('t', 1)

	data["Ipola"] = data["I0"] * data["DoLP"]
	data["Inonpola"] = data["I0"] * (1 - data["DoLP"])

	# nb_data = data.index.stop
	#
	# wl = [360] * , 557.7, 427.8, 391.4] * int(nb_data / 4)
	# if nb_data % 4 == 1: wl.append([360])
	# elif nb_data % 4 == 2: wl.extend([360, 557.7])
	# elif nb_data % 4 == 3: wl.extend([360, 557.7, 427.7])
	#
	# data['wl'] = pd.Series(wl, index=data.index)

	return data

# def PlotData(data, x_axis, y_axis, t=None, a_pc=None, e_pc=None, src_dist=None, dlos=None, min_alt=None, max_alt=None, wl=None):


def PlotData(data, x_axis_name, y_axis_name, **kwargs):
	missing = []
	for c in data.columns:
		if c not in kwargs.keys() and c not in [x_axis_name, y_axis_name, "I0", "DoLP", "AoRD", "Ipola", "Inonpola"]:
			missing.append(c)

	if kwargs:
		for k, v in kwargs.items():
			data = data[abs(data[k] - v) < 0.000001]
			data = data[data["I0"] < 100]
		data = data.reset_index()

	if missing:
		axs = PlotMissingData(data, missing, x_axis_name, y_axis_name, **kwargs)
	else:
		axs = SimplePlot(data, x_axis_name, y_axis_name, **kwargs)


	axs[-1].set_xlabel(x_axis_name.replace("_", "-"))
	axs[0].set_title(" ".join([k + "=" + str(v) for k, v in kwargs.items()]).replace("_", "-"))




def SimplePlot(data, x_axis_name, y_axis_name, **kwargs):

	if y_axis_name == "all":
		fig, axs = plt.subplots(3, sharex=True, figsize=(16, 8))
		fig.subplots_adjust(hspace=0)
		axs_names = ["I0", "DoLP", "Ipola"]
	else:
		fig, axs = plt.subplots(1, sharex=True, figsize=(16, 8))
		axs = [axs]
		axs_names = [y_axis_name]


	for iax, ax in enumerate(axs):
		x_axis, y_axis = [], []
		y_axis_name = axs_names[iax]
		for i in data.index:
			x_axis.append(data[x_axis_name][i])
			y_axis.append(data[y_axis_name][i])

		x_axis, y_axis = zip(*sorted(zip(x_axis, y_axis), key=lambda x: x[0]))
		ax.plot(x_axis, y_axis, "-*")
		x_axis, y_axis = [], []

		ax.set_ylabel(y_axis_name.replace("_", "-"))

	return axs

def PlotMissingData(data, missing, x_axis_name, y_axis_name, **kwargs):
	# missing = []
	# for c in data.columns:
	# 	if c not in kwargs.keys() and c not in [x_axis_name, y_axis_name, "I0", "DoLP", "AoRD", "Ipola", "Inonpola"]:
	# 		missing.append(c)

	# if kwargs:
	# 	for k, v in kwargs.items():
	# 		data = data[abs(data[k] - v) < 0.000001]
	# 		data = data[data["I0"] < 100]
	# 	data = data.reset_index()

	def GetMissingValues(i):
		missing_values = []
		for m in missing:
			missing_values.append(data[m][i])
		return missing_values

	if y_axis_name == "all":
		fig, axs = plt.subplots(3, sharex=True, figsize=(16, 8))
		fig.subplots_adjust(hspace=0)
		axs_names = ["I0", "DoLP", "Ipola"]
	else:
		fig, axs = plt.subplots(1, sharex=True, figsize=(16, 8))
		axs = [axs]
		axs_names = [y_axis_name]

	legend_title = "; ".join(missing).replace("_", "-")
	for iax, ax in enumerate(axs):
		old_missing_values = GetMissingValues(0)
		x_axis, y_axis = [], []
		y_axis_name = axs_names[iax]
		for i in data.index:
			missing_values = GetMissingValues(i)
			if old_missing_values == missing_values:
				x_axis.append(data[x_axis_name][i])
				y_axis.append(data[y_axis_name][i])
			else:
				ax.plot(x_axis, y_axis, "-*", label = "; ".join(str(v) for v in old_missing_values))
				x_axis, y_axis = [], []

			if i == data.index.stop - 1:
				ax.plot(x_axis, y_axis, "-*", label = "; ".join(str(v) for v in old_missing_values))

			old_missing_values = missing_values

		ax.set_ylabel(y_axis_name.replace("_", "-"))
		handles, labels = ax.get_legend_handles_labels()
		# sort both labels and handles by labels
		labels, handles = zip(*sorted(zip(labels, handles), key = lambda x: float(x[0])))

	axs[0].legend(handles, labels, title=legend_title, loc='upper left', bbox_to_anchor=(1, 1))

	return axs


def RunSystematicSimulation():
	azimuths = 0
	elevations = 90
	src_dist = [5]
	dlos = [0.1]
	alt_min = [0]
	alt_max = [6, 8, 9, 20, 30, 40]
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
	# RunSystematicSimulation()

	plt.close("all")
	data = GetData()
	PlotData(data, "max_alt", "all", a_pc = 0, e_pc = 90, src_dist=5, wl=557.7, min_alt=0, dlos=0.100)
	plt.show()
