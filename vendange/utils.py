#!/usr/bin/python3
# -*-coding:utf-8 -*

import numpy as np
import time as time
import datetime as dt


from rotation import *

from vendange_configuration import *


DtoR = np.pi / 180.
RtoD = 180. / np.pi
pwd = global_configuration.path
# pwd = "/home/bossel/These/Analysis/"
# pwd_data = pwd + "data/"
# pwd_src = pwd + "src/"
pwd_data = global_configuration.data_path
pwd_src = global_configuration.src_path


def VtoFDA(V, Vcos, Vsin):
	"""Return the flux F, the DoLP D [0-1] and AoLP A in radians from V, Vcos, Vsin"""
	F = 2*V
	DoLP = 2 * np.sqrt(Vcos**2 + Vsin**2) / V
	AoLP = np.arctan2(Vsin, Vcos) / 2 * RtoD

	return F, DoLP, AoLP

def FDAtoV(F, D, A):
	"""Return V, Vcos, Vsin from the flux F, the DoLP D [0-1] and AoLP A in radians"""

	V = F / 2.
	Vcos = F * D * np.cos(2 * A) / 4.
	Vsin = F * D * np.sin(2 * A) / 4.

	return V, Vcos, Vsin



def ReadInputFile(filename):
	"""Read a given input file and return a dictionnary. First word = key, second = value. Remove empty lines, as many arguments as you want, in any order that you want."""
	with open(filename, "r") as f:
		input = f.readlines()
	input = [l.split()[:2] for l in input if l != "\n"]
	dict = {}
	for i in input:
		dict[i[0]] = i[1]
	# print(dict)
	return dict

def LoadSPPData(file_names):
	"""Return a data array from a given data file (1 or more). Each line is a rotation, each being a list of data as written below."""
	# 0-5: 		year, month, day, hour, min, sec
	# 6-9: 		voltage, positive voltage, negative voltage, motor voltage
	# 10-11: 	Temperature of base, of filter
	# 12-13:	elevation, azimuth (geog coord)
	# 14-93:	angle of polariser (degrees)
	# 94-173:	canal pola
	# 174-253:	canal total
	# 254-259:	V, Vcos, Vsin, I0, DoLP, AoLP
	nb_data_per_rot  = 254
	# raw_data = np.array([])
	# data_path = "/home/bossel/These/Analysis/data/spp/"
	#Loading the data from a binary file.
	#Output is one big 1-D array of length nb_data_tot = nb_rot * nb_data_per_rot
	# for f in file_names:
	# 	f = data_path + f + ".spp"
	# 	raw_data = np.fromfile(f, dtype = np.int16)
	# 	# raw_data = np.concatenate((raw_data, np.fromfile(f, dtype=np.int16)), axis=0)
	f = file_names + ".spp"
	raw_data = np.fromfile(f, dtype = np.int16)
	# raw_data = np.concatenate((raw_data, np.fromfile(f, dtype=np.int16)), axis=0)
	nb_data_tot = len(raw_data)
	print("nb_data_tot", nb_data_tot)

	# print(raw_data[1001 * 254 : 1002 * 254])
	#reshape the array to be (nb_pts, nb_data per pts)
	return raw_data.reshape(int(nb_data_tot / nb_data_per_rot), nb_data_per_rot)

def LoadPTCUData(file_names, line = 1):
	"""Return an array of data and a dictionnary of config from a list of data files. The first is the data, the second is the config. Each line of the data is a rotation with 6 entries, the config dict is all the parameters of the observation, common to all rotations."""
	###DATA FILE
	# 0: 	IDProcessedData,
	# 1: 	IDConfiguration,
	# 2: 	Time since begining of obseravtion (timestamp) in milliseconds,
	# 3:	PMAvg,
	# 4:	PMCosAvg,
	# 5:	PMSinAvg
	# 6:	IDTemperatures
	# 7:	TempPM
	# 8:	TempOptical
	# 9:	TempAmbiant
	# 10:	IDLiveComment
	# 11:	Comment
	# 12:	live_commentscol
	###CONFIG FILE
	# 0:		IDConfiguration,
	# 1:		Timestamp,
	# 2-3:		CM_Latitude, CM_Longitude,
	# 4-5:		CM_Elevation, CM_Azimuth,
	# 6:		CM_PolarizerSpeed,
	# 7:		CM_Tilt,
	# 8:		CM_Usage,
	# 9:		CM_Comments,
	# 10:		IS_PolarizerOffset1,
	# 11:		IS_PolarizerOffset2,
	# 12:		IS_PolarizerOffset3,
	# 13:		IS_PolarizerOffset4,
	# 14:		IS_ConverterOffset4,
	# 15:		IS_ConverterOffset3,
	# 16:		IS_ConverterOffset2,
	# 17:		IS_ConverterOffset1,
	# 18:		IS_EngineTicksPerTour,
	# 19:		IS_MotorGearedRatio,
	# 20:		IS_QuantumEfficiency,
	# 21:		IS_IpAddress

	print("LOADING PTCU data: line", line)
	data_path = file_names
	#Loading the data from a binary file.
	#Output is one big 1-D array of length nb_data_tot = nb_toy * nb_data_per_rot
	data_file = data_path + "/data" + str(line) + ".csv"
	config_file = data_path + "/config.csv"

	print(data_file, config_file)
	try:
		with open(data_file, "r") as f:
			first_line = f.readlines()[1]
		d = GetDelimiter(first_line)
	except:
		d = ""
	# if d == ",":
	# 	print("Delimiter: ,")
	# elif d == " ":
	# 	print("Delimiter: space")
	# elif d == "\t":
	# 	print("Delimiter: tab")
	# print("genfromtxt")
	raw_data = np.genfromtxt(data_file, delimiter = d, skip_header=1)
	# print(raw_data[:50])
	# print("genfromtxt DONE")
	# print("genfromtxt")

	array_type = [('IDConfiguration',float),('Timestamp','S100'), ("CM_ID", float),('CM_Latitude',float),('CM_Longitude',float),('CM_Elevation',float),('CM_Azimuth',float),('CM_PolarizerSpeed',float),('CM_Tilt',float),('CM_Usage','S100'),('CM_Comments','S100'),('IS_PolarizerOffset1',float),('IS_PolarizerOffset2',float),('IS_PolarizerOffset3',float),('IS_PolarizerOffset4',float),('IS_ConverterOffset4',float),('IS_ConverterOffset3',float),('IS_ConverterOffset2',float),('IS_ConverterOffset1',float),('IS_EngineTicksPerTour',float),('IS_MotorGearedRatio',float),('IS_QuantumEfficiency',float),('IS_IpAddress',float)]
	# configuration = np.genfromtxt(config_file, dtype=array_type, delimiter=",", skip_header=1)
	configuration = pd.read_csv(config_file, sep=",")

	print(configuration)

	# print("genfromtxt DONE")

	print("PTCU DATA loaded")
	return raw_data, configuration

def GetDelimiter(line):
	if "," in line:
		return ","
	elif " " in line:
		return " "
	else:
		return "\t"

# def CleanRotation(rot, i_dict):
# 	nb_good_pts = i_dict["nb_pts_per_rot"]
# 	to_del = []
# 	for i in range(i_dict["angles"], i_dict["angles"] + i_dict["nb_pts_per_rot"]):
# 		if rot[i] == 0.:
# 			nb_good_pts -= 1
# 			to_del.extend([i, i + i_dict["nb_pts_per_rot"], i + 2*i_dict["nb_pts_per_rot"]])
# 	rot = np.delete(rot, to_del)
# 	return rot, nb_good_pts


def LoadSPPTestData(nb_rot = 100, I0=100, D=0.1, phi=0, file_name = ''):
	"""Return a test data array with the time, angle and fake voltage mesured every point as SPP would. The fake signal is a function of the form A*cos(theta-phi)+B where theta is the angle of the polariser and phi the angle of the polarisation"""

	#I0:  Intensity of the incoming light of the aurora (constant through time) arbitrary units
	#D:  Degree of polarized light (%)
	#phi:  Angle of the polarisation with respect to SPP (degrees)

	### Do not take into account the fact taht theta-phi cannot be > 90, because it would only change a sign, but it is squared in the end. So no difference. But I keep the ITheo function in case.
	# I = lambda theta, I0, D, phi: D * I0 * np.cos((theta - phi)*DtoR)**2 + (1-D) * I0 / 2 + D/5*np.random.normal(0, 1, theta.shape)
	I = lambda theta, I0, D, phi: I0 + np.random.normal(0, I0 * D, theta.shape)

	data = np.zeros((nb_rot, 254))
	da = 360./80
	angles = np.linspace(da, 360., 80) * DtoR
	dolps = np.linspace(0, 1., nb_rot)
	phis = np.linspace(-90, 90, nb_rot)
	for i in range(nb_rot):
		data[i][5] = i
		data[i][14:94] = angles
		data[i][94:174] = I(angles, I0, D, phi)
#		data[i][94:174] = I(angles, I0, D, phis[i])
		#data[i][94:174] = I(angles, I0, D, i*360/nb_rot)
		#data[i][94:174] = I(angles, I0, dolps[i], phi)
		data[i][174:254] = [I0]*80

		data[i][92:94] = [0,0]
		data[i][172:174] = [0,0]
		data[i][252:254] = [0,0]
	return data


def LoadPTCUTestData(nb_rot = 100, I0=100, D=0.1, phi=0, file_name = ''):
	"""Return a test data array with the time, V, Vcos, Vsin mesured every rotation as Petit Cru would. The fake signal is a function of the form A*cos(theta-phi)+B where theta is the angle of the polariser and phi the angle of the polarisation"""

	fake_config = {	"IDConfiguration":0,
					"Timestamp":"fake",
					"CM_Latitude":0,
					"CM_Longitude":0,
					"CM_Elevation":0,
					"CM_Azimuth":0,
					"CM_PolarizerSpeed":2,
					"CM_Tilt":0,
					"CM_Usage":"fake",
					"CM_Comments":"fake",
					"IS_PolarizerOffset1":0,
					"IS_PolarizerOffset2":0,
					"IS_PolarizerOffset3":0,
					"IS_PolarizerOffset4":0,
					"IS_ConverterOffset4":0,
					"IS_ConverterOffset3":0,
					"IS_ConverterOffset2":0,
					"IS_ConverterOffset1":0,
					"IS_EngineTicksPerTour":0,
					"IS_MotorGearedRatio":0,
					"IS_QuantumEfficiency":0,
					"IS_IpAddress":0}

	ptcu_data = np.zeros((nb_rot, 6))
	for i in range(nb_rot):
		I = I0  + np.random.normal(0, 1)
		A = phi + np.random.normal(0, 1) * DtoR
		D1= D + 1 / 100 * np.random.normal(0, 1)
		ptcu_data[i][0] = 0
		ptcu_data[i][1] = fake_config["IDConfiguration"]
		ptcu_data[i][2] = i
		ptcu_data[i][3] = I / 2.
		ptcu_data[i][4] =  (D1 * I / 4.) * np.cos(2 * A)
		ptcu_data[i][5] = -(D1 * I / 4.) * np.sin(2 * A)
	print(ptcu_data[0])

	return ptcu_data, fake_config

#
# def GetSmoothAverage(values, times, smoothing_factor, unit):
# 	nb_values = len(values)
# 	if unit == "rotations":
# 		nb_smooth_values = int(nb_values / smoothing_factor) + 1
#
# 		smooth_values = np.zeros(nb_smooth_values)
#
# 		for i in range(nb_values):
# 			smooth_values[int(i/smoothing_factor)] += float(values[i]) / smoothing_factor
#
# 	elif unit == "seconds":
# 		smooth_values = []
# 		i, j = 0, 0
# 		while i < nb_values:
# 			smooth_values.append(0)
# 			j = 0
# 			while 0 <= i + j < nb_values and abs(times[i + j] - times[i]) < smoothing_factor:
# 				smooth_values[-1] += values[i + j]
# 				j += 1
# 			smooth_values[-1] /= (j + 1)
# 			i += j + 1
#
# 	return smooth_values

def ReduceList(values, new_length):

	old_length = len(values)
	new_values = np.zeros(new_length)
	len_factor = old_length / new_length

	old_indexes = np.arange(0, old_length - 1, len_factor)

	for i in range(new_length):
		old_index = old_indexes[i]
		low_index, low_factor = int(old_index), 1 - old_index + int(old_index)
		high_index, high_factor = low_index + 1, 1 - low_factor

		new_values[i] = values[low_index] * low_factor + values[high_index] * high_factor

	# print("REDUCING:")
	# print(old_length, new_length)
	# print(values, new_values)

	return new_values

def GetAnglesDifference(A1, A2):

	diff = np.array(A2) - np.array(A1)

	# diff = np.where(diff <   0, - diff, diff)
	diff = np.where(diff >=   np.pi / 2, diff - np.pi, diff)
	diff = np.where(diff <=  -np.pi / 2, diff + np.pi, diff)
	# diff = np.where(diff <   0, - diff, diff)
	# diff = np.where(diff <= -np.pi, - diff - np.pi, diff)

	return diff


def GetSliddingAverage(values, times, smoothing_factor, unit):
	"""DEPRECATED, USE bottle.GetSliddingAverage() NOW! USING np.convolve() IT GOES MUCH FASTER.
	Given a list of values, a list of their associated times, the smoothing factor and its unit, return an array of smooth values. If unit==rotations each data is averaged with a window of smoothing_factor rotations centered on itself.
	If unit==seconds each data is averaged with a window of smoothing_factor seconds centered on itself. Both ways have the same results if each rotation has the same time step, but this is not True for the 20181115 PTCU data!!! This has still to be explained, but this is a temporary solution.
	"""
	nb_values = len(values)
	if unit == "rotations":
		window = np.ones(smoothing_factor) / smoothing_factor
		smooth_values = np.convolve(values, window, 'valid')
		# smooth_values = np.zeros(nb_values)
		# for i in range(nb_values):
		# 	n = 0.
		# 	for j in range(- int(smoothing_factor / 2.), int(smoothing_factor / 2.) + 1):
		# 		if 0 <= i + j < nb_values:
		# 			smooth_values[i] += float(values[i+j])
		# 			n += 1.
		# 	smooth_values[i] /= n
	else:
		if unit == "seconds":
			smoothing_factor = dt.timedelta(seconds = smoothing_factor)
		elif unit == "minutes":
			smoothing_factor = dt.timedelta(minutes = smoothing_factor)
		elif unit == "hours":
			smoothing_factor = dt.timedelta(hours = smoothing_factor)

		# smooth_values = np.zeros(nb_values)
		# for i in range(nb_values):
		# 	j, n, stop_up, stop_down = 1, 1., False, False
		# 	smooth_values[i] += values[i]
		# 	while not stop_up or not stop_down:
		# 		if i + j < nb_values and abs(times[i + j] - times[i]) < smoothing_factor / 2.:
		# 			smooth_values[i] += values[i + j]
		# 			n += 1.
		# 		else:
		# 			stop_up = True
		# 		if 0 <= i - j and abs(times[i - j] - times[i]) < smoothing_factor / 2.:
		# 			smooth_values[i] += values[i - j]
		# 			n += 1.
		# 		else:
		# 			stop_down = True
		# 		j += 1
		# 	smooth_values[i] = smooth_values[i] / n
		
	return smooth_values

def UnifyAngles(angles, manual_shift = -1):
	"""Given an array of angle between -90, 90, tries to unify the angle such that we don't have a jump from 90 to -90 when the average angle is around 90. The output can have angles > 90, but this is better to graph and read."""

	# for i in range(1, len(angles)):
	# 	if angles[i] <= angles[i-1] - 120 * DtoR:
	# 		new_angles[i] += np.pi
	# 		angles[i] += np.pi
	# 	elif angles[i] > angles[i-1] + 120 * DtoR:
	# 		new_angles[i] -= np.pi
	# 		angles[i] -= np.pi

	shift = 0
	if manual_shift == -1:
		new_angles = np.array(angles)
		angles_std = np.std(angles)
		for i in range(len(angles)):
			if angles[i] <= 0:
				new_angles[i] += np.pi
		new_angles_std = np.std(new_angles)

		if new_angles_std < angles_std:
			angles = new_angles
			shift = 1

	elif manual_shift == 1:
		for i in range(len(angles)):
			if angles[i] <= 0:
				angles[i] += np.pi
		shift = 1

	return angles, shift

def SetAngleBounds(angles, min, max, mod = np.pi, unit="radians"):
	bounds_size = max - min
	try:
		angles = angles.tolist()
		for i, a in enumerate(angles):
			while angles[i] <= min:
				angles[i] += mod
			while angles[i] > max:
				angles[i] -= mod

		angles = np.array(angles)
	except:
		while angles <= min:
			angles += mod
		while angles > max:
			angles -= mod

	return angles

def GetInstrumentName(file_name):
	if "ptcu" in str(file_name):
		return "ptcu"
	elif "spp" in str(file_name):
		return "spp"



def FindArgument(arg_string, arguments, default_value = False, getindex = 0):
	"""Search for the arg_string argument in the list of arguments passed on the command line (arguments). If found, returns the value of the parameter (or the argument following at getindex index).
	If not found, returns the default_value.
	Only works for bool arguments (true, false, 0, 1, existing or not)"""
	try:
		index = arguments.index(arg_string) + getindex
		result = arguments[index]
		# print(arg_string, index, result)
		# del arguments[index]
		# if getindex == 1:
		# 	del arguments[index - getindex]
	except:
		result = default_value
		return result, arguments

	if   result.lower() in ["false", "f"]:
		result = False
	elif result.lower() in ["true", "t", arg_string]:
		result = True
	elif arg_string not in  ["-a", "-b"]:
		result = bool(int(result))

	print(arg_string, result)
	return result, arguments


		# def GetInfoFromDataFileName(f, data_f):
		# 	folders, data_folders = f.split("/"), data_f.split("/")
		# 	# i = len(data_folders) - 1
		# 	i = -3
		# 	instrument_name = folders[i]
		# 	date, location, filters, azimut, elevation, com = "", "", "", 0, 0, ""
		# 	if instrument_name == "ptcu" or instrument_name == "ptcu_v2" :
		# 		date_location = folders[i+1].split("_")
		# 		date, location = date_location[0], date_location[1]
		# 		rest = folders[i+2].split("_")
		# 		print(rest)
		# 		for r in rest:
		# 			if len(r) == 2 and (r[0] and r[1]) in ["r", "v", "b", "m", "0", "o"]:
		# 				filters = r
		# 			elif r[0] == "a":
		# 				azimut = float(r[1:]) * DtoR
		# 			elif r[0] == "e":
		# 				elevation = float(r[1:]) * DtoR
		# 			else:
		# 				com += "_" + r
		# 	elif instrument_name == "gdcu":
		# 		date_location = folders[i+1].split("_")
		# 		date, location = date_location[0], date_location[1]
		# 		rest = folders[i+2].split("_")
		# 		print(rest)
		# 		for r in rest:
		# 			if len(r) == 4 and (r[0] and r[1] and r[2] and r[3]) in ["r", "v", "b", "m", "0", "o"]:
		# 				filters = r
		# 			elif r[0] == "a":
		# 				azimut = float(r[1:]) * DtoR
		# 			elif r[0] == "e":
		# 				elevation = float(r[1:]) * DtoR
		# 			else:
		# 				com += "_" + r
		#
		# 	elif instrument_name == "spp":
		# 		print(instrument_name, date, location, filters, azimut*RtoD, elevation*RtoD, com)
		# 		if len(folders[i+1].split("_")) == 2:
		# 			date, location = folders[i+1].split("_")[0], folders[i+1].split("_")[1]
		# 		else:
		# 			date = folders[i+1]
		# 		print(instrument_name, date, location, filters, azimut*RtoD, elevation*RtoD, com)
		#
		# 		rest = folders[i+2].split("_")
		# 		print(rest)
		# 		for r in rest:
		# 			if len(r) == 2 and ((r[0] and r[1]) in ["r", "v", "b", "m", "0"]):
		# 				filters = r
		# 			elif r[0] == "a":
		# 				azimut = float(r[1:]) * DtoR
		# 			elif r[0] == "e":
		# 				elevation = float(r[1:]) * DtoR
		# 			else:
		# 				com += "_" + r
		#
		# 		# if len(rest) <= 2:
		# 		# 	azimut, elevation = float(rest[0][1:]) * DtoR, float(rest[1][1:]) * DtoR
		# 		# elif len(rest) > 2:
		# 		# 	filters, azimut, elevation = rest[0], float(rest[1][1:]) * DtoR, float(rest[2][1:]) * DtoR
		# 		# else:
		# 		# 	com = folders[i+2].split("_")
		#
		# 	print(instrument_name, date, location, filters, azimut*RtoD, elevation*RtoD, com)
		# 	return instrument_name, date, location, filters, azimut, elevation, com
