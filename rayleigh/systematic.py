#!/usr/bin/python3
# -*-coding:utf-8 -*

from main import *
from rayleigh_input import *
import pandas as pd
import matplotlib as mpl

import sys


def GetData():
	# data = np.genfromtxt("./log/systematic_results.csv", delimiter=",", comments = "t", names=['t', 'a_pc', 'e_pc', 'src_dist', 'dlos', 'min_alt', 'max_alt', 'I0', 'DoLP', 'AoRD'])
	data = pd.read_csv("./log/systematic_results.csv", sep=",", comment = "t", names=['t', 'a_pc', 'e_pc', 'src_dist', 'dlos', 'min_alt', 'max_alt', "wl", "mad", 'I0', 'DoLP', 'AoRD', 'Idir'])

	data = data.drop('t', 1)

	data["Ipola"] = data["I0"] * data["DoLP"] / 100.
	data["Inonpola"] = data["I0"] * (1 - data["DoLP"] / 100.)

	data["mag"] = -2.5 * np.log10(data["I0"])

	try:
		data["mad"] = np.round((data["mad"] * RtoD) % 360)
	except:
		pass

	# nb_data = data.index.stop
	#
	# wl = [360] * , 557.7, 427.8, 391.4] * int(nb_data / 4)
	# if nb_data % 4 == 1: wl.append([360])
	# elif nb_data % 4 == 2: wl.extend([360, 557.7])
	# elif nb_data % 4 == 3: wl.extend([360, 557.7, 427.7])
	#
	# data['wl'] = pd.Series(wl, index=data.index)

	# print(data)
	return data

# def PlotData(data, x_axis, y_axis, t=None, a_pc=None, e_pc=None, src_dist=None, dlos=None, min_alt=None, max_alt=None, wl=None):


def PlotData(data, x_axis_name, y_axis_name, zaxis=None, **kwargs):
	font = {'size'   : 24}
	matplotlib.rc('font', **font)

	missing = []
	for c in data.columns:
		if ((c not in kwargs.keys() or (c in kwargs.keys() and kwargs[c] == "*")) and c not in [x_axis_name, y_axis_name, zaxis, "I0", "DoLP", "AoRD", "Ipola", "Inonpola", "mag", "Idir"]) :
			missing.append(c)
			# print(missing)

	if kwargs:
		for k, v in kwargs.items():
			if v == "*":
				continue
			try:
				data = data[abs(data[k] - v) < 0.000001]
			except:
				pass
				# print(data[k], v)
			# data = data[data["I0"] < 100]
		data = data.reset_index()
		if data.empty:
			raise ValueError('No simulation correspond to the following given parameters: {}'.format(kwargs))
		# print(data)

	if missing:
		axs = PlotMissingData(data, missing, x_axis_name, y_axis_name, **kwargs)
	elif not zaxis:
		axs = SimplePlot(data, x_axis_name, y_axis_name, **kwargs)
	else:
		axs = ParameterMap(data, x_axis_name, y_axis_name, zaxis, **kwargs)



	axs[-1].set_xlabel(x_axis_name.replace("_", "-"))
	# axs[0].set_title(" ".join([k + "=" + str(v) for k, v in kwargs.items()]).replace("_", "-"))


def ParameterMap(data, x_axis_name, y_axis_name, z_axis_name, **kwargs):
	fig, axs = plt.subplots(1, sharex=True, figsize=(16, 8))
	axs = [axs]
	axs_names = [y_axis_name]

	x_axis = list(set(data[x_axis_name]))
	y_axis = list(set(data[y_axis_name]))
	z_axis = data[z_axis_name]

	x_axis = sorted(x_axis)
	y_axis = sorted(y_axis)

	nb_x = len(x_axis)
	nb_y = len(y_axis)

	x_axis_bins = x_axis
	x_axis_bins.append(2 * x_axis[-1] - x_axis[-2])
	y_axis_bins = y_axis
	y_axis_bins.append(2 * y_axis[-1] - y_axis[-2])

	z_map = np.zeros((nb_y, nb_x))

	for i in data.index:
		x = data[x_axis_name][i]
		y = data[y_axis_name][i]

		x_index = x_axis.index(x)
		y_index = y_axis.index(y)

		z_map[y_index][x_index] = data[z_axis_name][i]

	print(z_map.shape[-1])


	m = axs[0].pcolormesh(x_axis_bins, y_axis_bins, z_map, norm=matplotlib.colors.LogNorm())

	axs[0].set_xlabel(x_axis_name.replace("_", "-"))
	axs[0].set_ylabel(y_axis_name.replace("_", "-"))

	cbar1 = fig.colorbar(m, extend='both', spacing='proportional', shrink=0.9, ax=axs[0])
	cbar1.set_label(z_axis_name)

	return axs


def SimplePlot(data, x_axis_name, y_axis_name, **kwargs):
	print("SimplePlot")

	### Experimentation to test effect of dlos parameter on flux
	# ranges = lambda dlos: np.arange(dlos/2, 100/np.sin(np.pi/4)-dlos/2, dlos)
	# def closest(R, ranges):
	# 	c = ranges[-1]
	# 	R = R / np.sin(np.pi/4.)
	# 	for r in ranges:
	# 		if abs(R-r) < c:
	# 			c = abs(R-r)
	# 	return c

	if y_axis_name == "all":
		fig, axs = plt.subplots(2, sharex=True, figsize=(16, 8))
		fig.subplots_adjust(hspace=0)
		axs_names = ["I0;Ipola;Inonpola", "DoLP"]
		# fig, axs = plt.subplots(3, sharex=True, figsize=(16, 8))
		# fig.subplots_adjust(hspace=0)
		# axs_names = ["I0", "DoLP", "Ipola"]
	else:
		fig, axs = plt.subplots(1, sharex=True, figsize=(16, 8))
		axs = [axs]
		axs_names = [y_axis_name]

	for iax, ax in enumerate(axs):
		for yax in axs_names[iax].split(";"):
			x_axis, y_axis = [], []
			y_axis_name = yax
			# y_axis_name = axs_names[iax]
			for i in data.index:
				x_axis.append(data[x_axis_name][i])
				y_axis.append(data[y_axis_name][i])

			x_axis, y_axis = zip(*sorted(zip(x_axis, y_axis), key=lambda x: x[0]))
			ax.plot(x_axis, y_axis, "-*", label=y_axis_name)

			### Experimentation to test effect of dlos parameter on flux
			# dist = [closest(90, ranges(dlos)) for dlos in x_axis]
			# a = axs[0].twinx()
			# print(x_axis, dist)
			# a.plot(x_axis, dist, "r-*")

			x_axis, y_axis = [], []


		if len(axs_names[iax].split(";")) > 1:
			ax.legend()
		ax.set_ylabel(axs_names[iax].replace("_", "-"))




	return axs

def PlotMissingData(data, missing, x_axis_name, y_axis_name, **kwargs):
	print("PlotMissingData")

	x_axis_shift = 0
	x_axis_mod = 0

	def GetXAxis(x, x_axis_shift=None, x_axis_mod=None):
		if x_axis_shift:
			x += x_axis_shift
		if x_axis_mod:
			x %= x_axis_mod
		return x

	GetMissingValues = lambda i, missing: [data[m][i] for m in missing]


	if y_axis_name == "all":
		fig, axs = plt.subplots(2, sharex=True, figsize=(16, 8))
		fig.subplots_adjust(hspace=0)
		axs_names = ["I0", "DoLP"]#, "Ipola"]
	else:
		fig, axs = plt.subplots(1, sharex=True, figsize=(16, 8))
		axs = [axs]
		axs_names = [y_axis_name]

	legend_title = "; ".join(missing).replace("_", "-")
	nb_plots = len(set(data[missing[0]]))

	missing_value_list = sorted(list(set(data[missing[0]])))

	for iax, ax in enumerate(axs):
		#### Code from here to next "####" comment works best for only 1 missing argument
		# Set the default color cycle
		colors = [(i, 0, 1-i) for i in np.linspace(0, 1, nb_plots)]
		ax.set_prop_cycle("c", colors)

		x_axis, y_axis = [], []
		y_axis_name = axs_names[iax]

		for m in missing_value_list:
			tmp = data[data[missing[0]] == m]
			x_axis = GetXAxis(tmp[x_axis_name], x_axis_shift, x_axis_mod)
			y_axis = tmp[y_axis_name]

			x_axis, y_axis = zip(*sorted(zip(x_axis, y_axis), key = lambda x: x[0]))

			ax.plot(x_axis, y_axis, "-*", label = str(m))

		#### Following code should work if several arguments are missing
		# old_missing_values = GetMissingValues(0, missing)
		# x_axis, y_axis = [], []
		# y_axis_name = axs_names[iax]
		# for i in data.index: #Every line of the table once
		# 	missing_values = GetMissingValues(i, missing) #Get value for changing paramter. missing should be of length 1
		# 	if old_missing_values == missing_values: #If value for changing parameter do not change, append the plot
		# 		x_axis.append(GetXAxis(data[x_axis_name][i], x_axis_shift, x_axis_mod))
		# 		y_axis.append(data[y_axis_name][i])
		# 	else:
		#
		# 		x_axis, y_axis = zip(*sorted(zip(x_axis, y_axis), key = lambda x: x[0]))
		# 		ax.plot(x_axis, y_axis, "-*", label = "; ".join(str(v) for v in old_missing_values))
		# 		x_axis, y_axis = [GetXAxis(data[x_axis_name][i], x_axis_shift, x_axis_mod)], [data[y_axis_name][i]]
		#
		# 	if i == data.index.stop - 1:
		# 		x_axis, y_axis = zip(*sorted(zip(x_axis, y_axis), key = lambda x: x[0]))
		# 		ax.plot(x_axis, y_axis, "-*", label = "; ".join(str(v) for v in old_missing_values))
		#
		# 	old_missing_values = missing_values

		ax.set_ylabel(y_axis_name.replace("_", "-"))
		handles, labels = ax.get_legend_handles_labels()
		# sort both labels and handles by labels
		try:
			labels, handles = zip(*sorted(zip(labels, handles), key = lambda x: float(x[0])))
		except:
			labels, handles = zip(*sorted(zip(labels, handles), key = lambda x: x[0]))

	axs[0].legend(handles, labels, title=legend_title)#, loc='upper left', bbox_to_anchor=(0.98, 1))



	return axs


def RunSystematicSimulation():
	"""Automatically runs simulations for the given paramter space.
	1: Define parameter space and check for loops include all wanted paramters.
	2: For every other parameters of the simulation, the default value is taken. default value is the value written in the input file: "RS_default.in"
	3: From the terminal run ./systematic.py
	4: To plot parameter comparisons run ./systematic.py with the right keyword arguments.
	"""

	### Define paramter space. All combinations of these parameters will be simulated.
	azimuths = 0
	elevations = 45
	dlos = [0.1]
	alt_min = [0] #, 1, 5, 10, 25, 50, 100]
	alt_max = [1, 100] #np.arange(5, 100, 5) #[91, 92, 93, 94]
	wavelengths =  [557.7]
	max_angle_discretization = [1] 	#[1, 2, 5, 10, 15, 20, 40, 50, 90, 180, 360]
	ground_emission_radius = [50]  # in km.
	ground_N_bins_max = [1000]

	for ger in ground_emission_radius:
		for gNbm in ground_N_bins_max:
			for mad in max_angle_discretization:
				for sa in src_az:
					for d in src_dist:
						for dl in dlos:
							for am in alt_min:
								for aM in [a for a in alt_max if a > am]:
									for wl in wavelengths:
										h = 90
										if wl == 630:
											h = 220
										elif wl == 557.7:
											h = 100

										in_dict = {	"azimuts": azimuths,
													"elevations": elevations,
													"resolution_along_los": dl,
													"RS_min_altitude": am,
													"RS_max_altitude": aM,
													"wavelength": wl,
													"emission_altitude": h,
													"max_angle_discretization": mad,
													"ground_emission_radius": ger,
													"ground_N_bins_max": gNbm

										r_input = RayleighInput()

										r_input.Update(in_dict)

										file_name = r_input.WriteInputFile(file_name = "systemic/systemic", **in_dict)

										print(file_name)

										log_file_name = "./log/" + file_name + ".out"
										log_file = open(log_file_name, "w")

										terminal = sys.stdout
										sys.stdout = log_file
										sys.stderr = log_file

										RunSimulation(file_name, show = False, output_result_file = "log/systematic_results.csv")

										sys.stdout = terminal
										sys.stderr = terminal


if __name__ == "__main__":
	#Get the input file from the command arguments
	arguments = sys.argv
	nb_args = len(arguments)

	if nb_args == 1:
		RunSystematicSimulation()
	else:
		plt.close("all")
		data = GetData()
		# PlotData(data, "src_dist", "max_alt", zaxis="I0", e_pc = 90, a_pc=0, src_dist="*", wl=391.4, min_alt=0, max_alt="*", dlos=0.1, mad=1)
		PlotData(data, "max_alt", "all", e_pc = 90, a_pc=180, src_dist="None", wl=391.4, min_alt=0, max_alt="*", dlos=0.1, mad="*")
		# PlotData(data, "max_alt", "all", e_pc = 90, a_pc=0, src_dist=5, wl=391.4, min_alt=0, max_alt="*", dlos=0.1, mad=1)
		plt.show()
