#!/usr/bin/python3
# -*-coding:utf-8 -*

import sys as sys
import numpy as np
import pandas as pd
import scipy.integrate as int
import matplotlib as mpl
import matplotlib.pyplot as plt

import scipy.optimize as opt #import curve_fit

from bottle import *

DtoR = np.pi / 180.
RtoD = 1. / DtoR

def GetVParam(F, D, A):
	V = F / 2.
	Vc = F * D/100. * np.cos(2 * A*DtoR) / 4.
	Vs = F * D/100. * np.sin(2 * A*DtoR) / 4.
	return V, Vc, Vs

def GetFDAParam(V, Vc, Vs):
	F = V *2.
	D = 2 * np.sqrt(Vc**2 + Vs**2) / V * 100
	A = np.arctan2(Vs, Vc) / 2 * RtoD
	return F, D, A

def GetVParamDF(data):
	data = data.copy()
	data["V"] = data["I0"] / 2.
	data["Vcos"] = data["I0"] * data["DoLP"]/100. * np.cos(2 * data["AoRD"]*DtoR) / 4.
	data["Vsin"] = data["I0"] * data["DoLP"]/100. * np.sin(2 * data["AoRD"]*DtoR) / 4.

	data["DV"] = data["DI0"] / 2.
	data["DVcos"] = data["DI0"] * data["DDoLP"]/100. * np.cos(2 * data["DAoLP"]*DtoR) / 4.
	data["DVsin"] = data["DI0"] * data["DDoLP"]/100. * np.sin(2 * data["DAoLP"]*DtoR) / 4.

	data["RSV"] = data["V"] - data["DV"]
	data["RSVcos"] = data["Vcos"] - data["DVcos"]
	data["RSVsin"] = data["Vsin"] - data["DVsin"]

	return data

def GetFDAParamDF(data):
	data = data.copy()
	data["I0"] = data["V"] * 2.
	data["DoLP"] = 2 * np.sqrt(data["Vcos"]**2 + data["Vsin"]**2) / data["V"] * 100
	data["AoRD"] = np.arctan2(data["Vsin"], data["Vcos"]) / 2 * RtoD

	data["DI0"] = data["DV"] * 2.
	data["DDoLP"] = 2 * np.sqrt(data["DVcos"]**2 + data["DVsin"]**2) / data["DV"] * 100
	data["DAoLP"] = np.arctan2(data["DVsin"], data["DVcos"]) / 2 * RtoD

	data["RSI0"] = data["RSV"] * 2.
	data["RSDoLP"] = 2 * np.sqrt(data["RSVcos"]**2 + data["RSVsin"]**2) / data["RSV"] * 100
	data["RSAoRD"] = np.arctan2(data["RSVsin"], data["RSVcos"]) / 2 * RtoD

	return data

def UpdateDF(df):
	df = GetVParamDF(df)
	df = GetFDAParamDF(df)
	return df

def AddDF(df1, df2, f1=1, f2=1):
	sum = df1.copy()
	for n in ["V", "DV", "RSV"]:
		sum[n] = df1[n] * f1 + df2[n] * f2

	for n in ["Vcos", "Vsin",  "DVcos", "DVsin", "RSVcos", "RSVsin"]:
		sum[n] += df2[n]

	sum = GetFDAParamDF(sum)

	return sum

def SubtractDirect(df):
	df = GetVParamDF(df)
	df["V"] -= df["DV"]
	df["Vcos"] -= df["DVcos"]
	df["Vsin"] -= df["DVsin"]

	df["DV"] -= df["DV"]
	df["DVcos"] -= df["DVcos"]
	df["DVsin"] -= df["DVsin"]

	return GetFDAParamDF(df)



def Shift(data, shift, x_axis="azimuts"):
	data = UpdateDF(data)

	for n in ["V", "DV", "RSV", "Vcos", "Vsin",  "DVcos", "DVsin", "RSVcos", "RSVsin"]:
		data[n] = np.interp(data[x_axis] - shift, data[x_axis], data[n])

	# data[x_axis] += shift

	data = GetFDAParamDF(data)

	return data

def FindRatio(goal_r, x, y, delta = 0.1):

	avg_x, avg_y = np.average(x), np.average(y)

	def FindMinInRange(avg, scale):
		fy = avg * 10.**np.arange(-scale, scale, scale/1000)
		ratio = []
		diff = []
		min_diff = False
		min_factor = False
		for f in fy:
			CL = x + f * y
			r = max(CL) / min(CL)
			d = abs(goal_r - r)
			ratio.append(r)
			diff.append(d)
			if not min_diff or min_diff > d:
				min_diff = d
				min_factor = f
		return min_factor

	min_f = avg_x / avg_y
	for s in [100, 10, 1, 0.1, 0.01, 0.001]:
		min_f = FindMinInRange(min_f, s)

	CL = x + min_f * y
	print(min_f, max(CL) / min(CL))
	return min_f

#
# data_path = "/home/bossel/These/Analysis/data/ptcu/"
# data_file = "20200307_Skibotn/vr_rot_e45/"
#
# ptcu_data = pd.read_csv(data_path + data_file, sep="\t")
# bottle = Bottle(data_path + data_file, from_txt = True)
# #"fixed_elevation_discrete_rotation"
# #"fixed_elevation_continue_rotation"
# bottle_type = "fixed"



### File 1 should be the sky, file 2 the ground
path = "/home/bossel/These/Documentation/Mes_articles/auroral_pola/"
# path = "/home/bossel/These/Analysis/results/rayleigh/"


file_1 = path + "grd_only/" + "t_e52_grd_only_aero_1low.txt"
file_2 = path + "sky_only/" + "uni_sky_t_e52_aero_1low_albedo.txt"


m1 = pd.read_csv(file_1)
m2 = pd.read_csv(file_2)

### If a model is a single observation, repeat it to match the length of the other model
if m1.shape[0] == 1 and m2.shape[0] != 1:
	m1 = pd.concat([m1]*m2.shape[0], ignore_index=True)
elif m2.shape[0] == 1 and m1.shape[0] != 1:
	m2 = pd.concat([m2]*m1.shape[0], ignore_index=True)


x_type = "datetime"
# x_type = "azimuts"

# param, cov = FitModelToFlux(ptcu_data, m1, m2)
# f1, f2 = param[0], param[1]
# best_flux_model = AddDF(m1, m2, f1=f1, f2=f2)
#
# param, cov = FitDDoLPToDoLP(ptcu_data, m1, m2, f1=f1, f2=f2)
# m1["DDoLP"] = m1["DDoLP"] * param[0] + param[1]
# m1 = UpdateDF(m1)
# best_DoLP_model = AddDF(m1, m2, f1=f1, f2=f2)

model_list = []

### OLD manual search for best fit

factor_1 = np.array([4335.3242154348345 * 0.39273535777746454])
factor_2 = np.array([25450.87990784447 * 2.5])
# factor_2 = np.array([0.01, 0.1, 1, 10, 100]) * 1e9
# factor_2 = np.array([1. / np.average(m2["I0"])]) * 1.765e-3 ### pour rotation a la gorge

# m1 = SubtractDirect(m1)

# f2_shift = 0
# m2 = Shift(m2, f2_shift, x_axis = "azimuts")

DoLP_factor_2 = 1 #0.3

# factor_2 = np.array([1]) #* 1e10

# ratio_list = [1.25] #[1.75, 1.25, 3.60, 1.83]
# factor_2 = np.array([FindRatio(r, m1["I0"], m2["I0"]) for r in ratio_list])

for if1, f1 in enumerate(factor_1):
	for if2, f2 in enumerate(factor_2):
		model_1 = m1.copy()
		model_2 = m2.copy()

		model_1["I0"] *= f1
		model_1["DI0"] *= f1
		model_2["I0"] *= f2
		model_2["DI0"] *= f2

		# DoLP_I0_link = abs(np.gradient(model_2["I0"]))
		model_2["DoLP"] *= DoLP_factor_2 #* np.interp(DoLP_I0_link, [np.min(DoLP_I0_link), np.max(DoLP_I0_link)], [1, 0])

		model_1 = GetVParamDF(model_1)
		model_2 = GetVParamDF(model_2)

		model_sum = AddDF(model_1, model_2)

		model_sum = GetFDAParamDF(model_sum)

		model_list.append((f1, f2, model_sum))

		print(model_sum["V"], model_1["V"], model_2["V"])


f, axs = plt.subplots(3, sharex=True, figsize=(16, 8))

try:
	x_axis = [mpl.dates.datestr2num(t) for t in model_sum["datetime"]]
except:
	x_axis = model_sum["azimuts"]


if model_1.shape[0] == model_2.shape[0]:
	axs[0].plot(x_axis, model_1["I0"], "+r")
	axs[0].plot(x_axis, model_2["I0"], "+g")
	axs[1].plot(x_axis, model_1["DoLP"], "+r")
	axs[1].plot(x_axis, model_2["DoLP"], "+g")
	axs[2].plot(x_axis, model_1["AoRD"], "+r")
	axs[2].plot(x_axis, model_2["AoRD"], "+g")

for f1, f2, model in model_list:
	axs[0].plot(x_axis, model["I0"], "*", label = f"{f1}, {f2}")
	axs[1].plot(x_axis, model["DoLP"], "*", label = f"{f1}, {f2}")
	axs[2].plot(x_axis, model["AoRD"], "*", label = f"{f1}, {f2}")
	# axs[0].plot(x_axis, model["V"], "*")
	# axs[0].plot(x_axis, model_1["V"], "*")
	# # axs[0].plot(x_axis, model_2["V"], "*")
	# axs[1].plot(x_axis, model["Vcos"], "*")
	# axs[1].plot(x_axis, model_1["Vcos"], "*")
	# # axs[1].plot(x_axis, model_2["Vcos"], "*")
	# axs[2].plot(x_axis, model["Vsin"], "*")
	# axs[2].plot(x_axis, model_1["Vsin"], "*")
	# axs[2].plot(x_axis, model_2["Vsin"], "*")

	save1 = ".".join(file_1.split("/")[-1].split(".")[:-1])
	save2 = ".".join(file_2.split("/")[-1].split(".")[:-1])
	save_name = f"{path}{save1}x{f1}+{save2}x{f2}"
	print(save_name)
	model.to_csv(save_name, index=False)


axs[0].set_ylabel("Flux (nW)")
axs[1].set_ylabel("DoLP (\%)")
axs[2].set_ylabel("AoLP (deg)")
axs[2].set_xlabel("Azimut (°)")

plt.legend()
plt.show()
