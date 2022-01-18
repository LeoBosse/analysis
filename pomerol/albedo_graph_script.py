#!/usr/bin/python3
# -*-coding:utf-8 -*

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rc('text', usetex=False)
mpl.rcParams.update({'font.size': 30})


def Albedo(altitudes, fct_str):

	lim_layers = [float(a) for a in fct_str.split("_")[1:-1:2]]
	lim_albedos = [float(a) for a in fct_str.split("_")[0::2]]

	albedos = [lim_albedos[int(np.searchsorted(lim_layers, a))] for a in altitudes]

	return albedos


altitudes = np.linspace(0, 1, 1000)
str_list = ["0.2_0.2_0.5_0.5_0.8", "0.1_0.2_0.3_0.5_0.5", "0.5_0.2_0.7_0.5_0.9"]#, "1_100", "0.5_100"]
file_legends = ["Varying albedo", "Low albedo", "High albedo", "Constant albedo=1", "Constant albedo=0.5"]
file_colors = ["red", "green", "blue", "magenta", "purple"]
line_type = ["-", "-", "-", "--", "--"]
# zorder = [10, 10, 10, 2, 2]

for i, alb_str in enumerate(str_list):
	plt.plot(altitudes, Albedo(altitudes, alb_str), line_type[i], label = file_legends[i], color = file_colors[i], zorder = 1)



plt.xlabel("Altitudes (km)")
plt.ylabel("Albedo")

# plt.legend(prop={'size': 22})
plt.show()
