#!/usr/bin/python3
# -*-coding:utf-8 -*

from mpi4py import MPI
# from mpi4py.futures import MPIPoolExecutor
mpi_comm = MPI.COMM_WORLD
mpi_rank = mpi_comm.Get_rank()
mpi_size = mpi_comm.Get_size()
mpi_name = mpi_comm.Get_name()

from main import *
from rayleigh_input import *
import pandas as pd
import matplotlib as mpl

mpl.rc('text', usetex=False)
mpl.rcParams.update({'font.size': 24})

import sys

import itertools

def GetData(file = "./log/systematic_results.csv"):
	# data = np.genfromtxt("./log/systematic_results.csv", delimiter=",", comments = "t", names=['t', 'a_pc', 'e_pc', 'src_dist', 'dlos', 'min_alt', 'max_alt', 'I0', 'DoLP', 'AoRD'])

	data = pd.read_csv(file, sep = ",", index_col = False)

	try:
		wrong_items = list(data[data['src_path'] == "src_path"].index)
		data = pd.read_csv(file, sep = ",", index_col = False, skiprows = lambda x: x-1 in wrong_items)
	except:
		pass

	# data = pd.read_csv("./log/systematic_results.csv", sep=",", comment = "t", names=['t', 'a_pc', 'e_pc', 'src_dist', 'dlos', 'min_alt', 'max_alt', "wl", "mad", 'I0', 'DoLP', 'AoRD', 'Idir'])

	# data = data.drop('t', 1)

	# data = data[data["path"] != "path"]

	data.reindex()
	print(data.columns)

	data["Ipola"] = data["I0"] * data["DoLP"] / 100.
	data["Inonpola"] = data["I0"] * (1 - data["DoLP"] / 100.)

	data["mag"] = -2.5 * np.log10(data["I0"])

	try:
		data["max_angle_discretization"] = np.round((data["max_angle_discretization"] * RtoD) % 360)
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
	# font = {'size'   : 24}
	# matplotlib.rc('font', **font)

	# print(data[x_axis_name], data["wavelength"])

	free_list = [] # list of parrameter names (column name) of free parameters
	for c in data.columns:
		if ((c in kwargs.keys() and kwargs[c] == "*") and c not in [x_axis_name, y_axis_name, zaxis, "I0", "DoLP", "AoRD", "Ipola", "Inonpola", "mag", "Idir"]) :
			free_list.append(c)
			# print(free_list)

	if kwargs:
		for k, v in kwargs.items():
			if v == "*":
				continue
			if type(v) != type(str()) and type(v) != type(bool()):
				data = data[abs(data[k] - v) < 0.000001]
			else:
				# print(f"V is :{v}: string for key :{k}:")
				# print(data[k])
				data = data[data[k] == v]
				# print(data[k], v)
			# data = data[data["I0"] < 100]
		data = data.reset_index()
		if data.empty:
			raise ValueError('No simulation correspond to the following given parameters: {}'.format(kwargs))
		# print(data)


	if free_list:
		axs = PlotFreeData(data, free_list, x_axis_name, y_axis_name, **kwargs)
	elif not zaxis:
		axs = SimplePlot(data, x_axis_name, y_axis_name, **kwargs)
	else:
		axs = ParameterMap(data, x_axis_name, y_axis_name, zaxis, **kwargs)

	### Plot Rayleigh scattering DoLP function for an infinitely far away point source
	### Start
	# N=100
	# a = np.linspace(0, 360*DtoR, N)
	# theta = np.arccos(np.cos(45*DtoR) * np.cos(a))
	# y = np.sin(theta)**2 / (1 + np.cos(theta)**2)
	#
	# axs[1].plot(a*RtoD, y*100, "k--")
	### End


	for ax in axs:
		ax.xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator())
		ax.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator())
		ax.grid(True)


	# if "azimut" in x_axis_name.lower():
	# 	x_axis_name = "azimuth"
	axs[-1].set_xlabel(x_axis_name.capitalize().replace("_", "-") + " (°)")

	axs[0].set_title(" ".join([k + "=" + str(v) for k, v in kwargs.items() if v != "*"]).replace("_", "-"))



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
		axs_names = ["I0;Ipola;Inonpola", "DoLP", "AoRD"]
		fig, axs = plt.subplots(len(axs_names), sharex=True, figsize=(16, 8))
		fig.subplots_adjust(hspace=0)
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

			# x_axis, y_axis = zip(*sorted(zip(x_axis, y_axis), key=lambda x: x[0]))
			x_axis, y_axis = SortAxis(x_axis, y_axis, x_axis_name, y_axis_name)
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

def PlotFreeData(data, free_list, x_axis_name, y_axis_name, **kwargs):
	# free_list == list of parrameter names (column name) of free parameters
	print("PlotFreeData")

	x_axis_shift = 0
	x_axis_mod = 0

	def GetXAxis(x, x_axis_shift=None, x_axis_mod=None):
		if x_axis_shift:
			x += x_axis_shift
		if x_axis_mod:
			x %= x_axis_mod
		return x

	# GetFreeValues = lambda i, free: [data[m][i] for m in free_list]

	#Initialize axis and figure
	if y_axis_name == "all":
		axs_names = ["I0", "DoLP", "AoRD"]
		axs_y_names = ["Flux (nW)", "DoLP (%)", "AoLP (°)"]
		fig, axs = plt.subplots(len(axs_names), sharex=True, figsize=(16, 8))
		fig.subplots_adjust(hspace=0)
	else:
		axs_names = y_axis_name.split(";")
		axs_y_names = axs_names
		fig, axs = plt.subplots(len(axs_names), sharex=True, figsize=(16, 8))
		if len(axs_names) == 1: axs = [axs]

	# legend_title = "Distance (km)"
	legend_title = "; ".join(free_list).replace("_", "-")

	#Single free parameter case
	# free_value_list = list(set(data[free_list[0]]))
	#Multiple free parameter case
	free_value_list = [list(set(data[f])) for f in free_list] #list of list of values for each free parameter in free_list

	if type(free_value_list[0][0]) == type(str()):
		if "point" in free_value_list[0][0]:
			print("DEBUG")
			free_value_list[0] = sorted(free_value_list[0], key = lambda x: int(x.split("_")[3][1:]))
	else:
		free_value_list = [sorted(fvl) for fvl in free_value_list]

	############################################################################
	### ONLY FOR POMEROL ARTICLE FIG 4 !!!!!!!! TO DECOMMENT !!!!!!
	### WILL NOT PLOT THE LAST TWO FREE PARAMETER VALUES !!!!!
	free_value_list[0] = free_value_list[0][:-2]
	############################################################################

	# Number of possible combinaison between all free parameters (number of lines to plots)
	nb_free_param = [len(fvl) for fvl in free_value_list]
	nb_plots = np.product(nb_free_param)

	#Creates a list of all possible free parameters dictionnaries with all different permutations possible.
	permutations_dicts = [dict(zip(free_list, v)) for v in itertools.product(*free_value_list)]
	# print(permutations_dicts)

	# Loop over the subplots if we want I0, DoLP and AoLP for example.
	for iax, ax in enumerate(axs):
		#### Code from here to next "####" comment works best for only 1 or 2 free_list argument

		# Set the default color cycle for the free paramter with the highest number of values
		if len(free_list) == 2:
			colors = ["black", "red", "blue", "yellow", "pink", "purple", "orange", "turquoise"]
			colors = [(i, 0, 1 - i) for i in np.linspace(0, 1, max(nb_free_param))]
		elif len(free_list) == 1:
			colors = [(i, 0, 1 - i) for i in np.linspace(0, 1, max(nb_free_param))]
		# Set the default style cycle for the free paramter with the lowest number of values
		styles = ["-","--", "-."]
		markers = [None, "x"]
		# ax.set_prop_cycle("c", colors)
		# ax.set_prop_cycle("linestyle", styles)

		x_axis, y_axis = [], []
		y_axis_name = axs_names[iax]

		# for ip, p in enumerate(free_list):
		# 	line_style = styles[ip]
		# 	for im, m in enumerate(free_value_list[ip]):

		for perm in permutations_dicts:
			keys, values = zip(*perm.items())
			indexes = [free_value_list[i].index(value) for i, value in enumerate(values)]

			tmp = data.copy()
			for k, v in perm.items():
				tmp = tmp[tmp[k] == v]

			print(perm)

			x_axis = GetXAxis(tmp[x_axis_name], x_axis_shift, x_axis_mod)
			y_axis = tmp[y_axis_name]
			x_axis = x_axis.reset_index(drop = True)


			x_axis, y_axis = SortAxis(x_axis, y_axis, x_axis_name, y_axis_name)
			# print(x_axis, y_axis)
			# if x_axis_name == "shader_mode":
			# 	if "dist" in x_axis[0]:
			# 		x_axis, y_axis = zip(*sorted(zip(x_axis, y_axis), key = lambda x: float(x[0].split("_")[2][1:])))
			# 	elif "ring" in x_axis[0]:
			# 		x_axis, y_axis = zip(*sorted(zip(x_axis, y_axis), key = lambda x: int(x[0].split("_")[1])))
			# else:
			# 	x_axis, y_axis = zip(*sorted(zip(x_axis, y_axis), key = lambda x: x[0]))

			# print(x_axis, y_axis)
			if len(free_list) == 2:
				linestyle = styles[indexes[0]]
				markerstyle = markers[indexes[0]]
				color = colors[indexes[1]]
			elif len(free_list) == 1:
				linestyle = styles[0]
				markerstyle = markers[0]
				color = colors[indexes[0]]
			# print(linestyle, color)

			ax.plot(x_axis, y_axis, linestyle = linestyle, marker = markerstyle, markersize = 10, color = color, label = ";".join([f"{v}" for v in values]))

		ax.set_ylabel(axs_y_names[iax].replace("_", "-"))
		handles, labels = ax.get_legend_handles_labels()
		# sort both labels and handles by labels
		try:
			labels, handles = SortAxis(labels, handles, free_list[0], y_axis_name, label=True)
			# labels, handles = zip(*sorted(zip(labels, handles), key = lambda x: float(x[0])))
		except:
			labels, handles = SortAxis(labels, handles, free_list[0], y_axis_name, label=True)
			# labels, handles = zip(*sorted(zip(labels, handles), key = lambda x: x[0]))

	axs[0].legend(handles, labels, title=legend_title, loc='upper left', bbox_to_anchor=(1, 1), fontsize="small", title_fontsize="xx-small")
	# axs[0].set_yscale("log")
	axs[0].tick_params(axis="y", which='minor')
	# axs[0].grid(axis="y", b=True, which='both', color='0.65', linestyle='-')
	# if len(axs) > 1:
	# 	axs[1].tick_params(axis="y", which='minor')
	# 	axs[1].minorticks_on()
	# 	axs[1].grid(axis="y", b=True, which='both', color='0.65', linestyle='-')
	return axs

def SortAxis(x_axis, y_axis, xname, yname, label=False):

	# print("SortAxis...")
	if xname == "shader_mode":
		mode = x_axis[0]
		if "dist" in mode:
			x_axis = np.array([float(x.split("_")[2][1:]) for x in x_axis])
			# x_axis, y_axis = zip(*sorted(zip(x_axis, y_axis), key = lambda x: float(x[0].split("_")[2][1:])))
		elif "ring" in mode:
			x_axis = np.array([int(x.split("_")[1]) for x in x_axis])
	elif xname == "ground_mode":
		if "point" in x_axis[0]:
			x_axis = np.array([int(x.split("_")[3][1:]) for x in x_axis])
			x_axis, y_axis = zip(*sorted(zip(x_axis, y_axis), key = lambda x: x[0]))
	else:
		# print(x_axis, y_axis)
		if label:
			# x_axis, y_axis = zip(*sorted(zip(x_axis, y_axis), key = lambda x: float(x[0].split(";")[1])))
			x_axis, y_axis = zip(*sorted(zip(x_axis, y_axis), key = lambda x: x[0]))
		else:
			print(x_axis, y_axis, xname, yname)
			x_axis, y_axis = zip(*sorted(zip(x_axis, y_axis), key = lambda x: float(x[0])))


	return x_axis, y_axis

def RunSystematicSimulation():
	"""Automatically runs simulations for the given paramter space.
	1: Define parameter space in in_dict dictionnary.
	2: For every other parameters of the simulation, the default value is taken. default value is the value written in the input file: "RS_default.in"
	3: Define the name of the file in wich to save the results. (just before the dictionnary output_result_file = <your_file_name>)
	3: From the terminal run ./systematic.py
	4: To plot parameter comparisons run ./systematic.py -p <your_file_name>. Hardcoded in the PlotData() function at the end of this script are several other paramters to allow you to choose which paramter to let free in the plots. Some exemples are commented.
	"""

	output_result_file = "log/systematic_zenith_wavelength"
	### Define parameter space. All combinations of these parameters will be simulated.
	in_dict = {	"azimuts": [0],
				"elevations": [90],
				# "use_analytic": ["false", "true"],
				# "Nb_points_along_los": [2**i for i in range(1, 9)]
				# "use_ozone": [0, 1],
				# "use_aerosol": [1],
				# "aer_complexity": [1],
				# "aer_model": ["rural", "urban", "maritime"],#, "test1", "test2", "test3", "test4", "test5"],
				# "aer_single_scatter_albedo": [0.9],
				# "aer_phase_fct_asym_g": [0.7] #[-1; 1], g>0 -> front scatter dominant. g<0 -> back scatter dominant

				# "aer_Qext": [2],
				# "aer_Qabs": [0.8],
				# "aer_radius": [0, 100, 200, 300, 500],
				# "aer_Hn": 	[1],
				# "aer_n0": 	[20000],
				# "aer_nBK": [300],
				# "Nb_points_along_los": [100],
				# "RS_min_altitude": [0],# np.arange(0, 100, 1), #, 1, 5, 10, 25, 50, 100]
				# "RS_max_altitude": np.arange(1, 10, 1),
				# "use_analytic": ["true", "false"],
				# "Atmospheric_profile": [f"atm_profiles/{p}.atm" for p in ["day", "equ", "ngt", "sum", "win"]],
				"wavelength": np.append(np.linspace(500, 600, 11), np.linspace(600, 900, 7)) # [540] 
				# "emission_altitude": [110],
				# "Nb_emission_points": [100, 500, 1000, 2000, 5000],
				# "max_angle_discretization": [0], 	#[1, 2, 5, 10, 15, 20, 40, 50, 90, 180, 360],
				# "ground_mode": ["point_I100_a180_d5"],
				# "ground_mode": [f"point_I100_a180_d{i}" for i in np.arange(1, 100)],
				# "shader_mode": [f"ring_{i}" for i in np.arange(0, 100)],
				# "ground_emission_radius": [50],  # in km.,
				# "ground_N_bins_max": [20, 50]

			}

	#Creates a list of all possible input dictionnaries with all different permutations possible.
	keys, values = zip(*in_dict.items())
	permutations_dicts = [dict(zip(keys, v)) for v in itertools.product(*values)]

	count = 0
	for perm in permutations_dicts:

		# if perm["RS_max_altitude"] != perm["RS_min_altitude"] + 1:
		# 	continue
		# if perm["aer_Qabs"] > perm["aer_Qext"]:
		# 	continue
		# if perm["use_aerosol"] == 1 and perm["aer_complexity"] == 1:
		# 	perm["aer_name"] = perm["aer_model"] + "_test"

		#Initialize the input file used for this simulation
		r_input = RayleighInput()

		#Update input file with the actual permutation of input parameters defined in in_dict
		r_input.Update(perm)

		#Write the input file and get its unique name (unique for this set of permutations, might be the same than older runs.)
		file_name = r_input.WriteInputFile(file_name = "systematic/systematic", **perm)
		# print("file_name", file_name)

		if mpi_rank == 0: print(file_name)

		#Create the log file for this run from the input file name
		log_file_name = "./log/" + file_name + ".out"
		log_file = open(log_file_name, "w")

		#Shift output to the log file
		terminal = sys.stdout
		sys.stdout = log_file
		sys.stderr = log_file

		#Run simulation and print the results in log/systematic_results.csv
		# with MPIPoolExecutor() as executor:
		# 	executor.map(RunSimulation, file_name, show = False, output_result_file = "log/systematic_results.csv", header = not bool(count))
		RunSimulation(file_name, show = False, output_result_file = output_result_file, header = not bool(count))

		#Shift output back to the terminal
		sys.stdout = terminal
		sys.stderr = terminal

		count += 1


if __name__ == "__main__":
	#Get the input file from the command arguments

	### To launch with:  mpiexec -n 4 python3 -m mpi4py.futures ./systematic.py

	arguments = sys.argv
	nb_args = len(arguments)

	if nb_args == 1:
		RunSystematicSimulation()
	else:
		plt.close("all")
		if nb_args == 3 and arguments[1] == "-p":
			data = GetData(file=arguments[2])
		else:
			data = GetData()

		# PlotData(data, "azimuts", "all",  aer_model="*", use_aerosol = "*")
		# PlotData(data, "azimuts", "all", use_analytic="*", Nb_points_along_los="*")
		# PlotData(data, "azimuts", "all", ground_N_bins_max="*")
		PlotData(data, "wavelength", "all",  wavelength='*')
		# PlotData(data, "azimuts", "all", elevations="*")
		plt.show()
