#!/usr/bin/python3
# -*-coding:utf-8 -*

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rc('text', usetex=False)
mpl.rcParams.update({'font.size': 30})

DtoR = np.pi / 180.
RtoD = 1. / DtoR



### FIGURE 6
# path = "/home/bossel/These/Documentation/Mes_articles/scattering/new_pictures/FIG6/"
# file_names = [	path + "o_grd_onlyx1+o_sky_albx1000000000.0",
# 				path + "o_grd_onlyx1+o_sky_albx10000000000.0",
# 				path + "o_grd_onlyx1+o_sky_albx100000000000.0",
# 				path + "o_grd_onlyx1+o_sky_albx20000000000.0",
# 				path + "o_grd_onlyx1+o_sky_albx30000000000.0",
# 				path + "o_grd_onlyx1+o_sky_albx40000000000.0"]
#
# file_legends = ["Non pola", "No albedo", "Up", "East", "North", "test"]
# direct_file = path + "o_grd_onlyx1+o_sky_albx1000000000.0"
# direct_data = pd.read_csv(direct_file)


### FIGURE 4
# path = "/home/bossel/These/Documentation/Mes_articles/scattering/new_pictures/"
# file_names = [	path + "FIG6_new/v_unisky_noB.txt",
# 				path + "FIG6_new/v_unisky_EW.txt",
#  				path + "FIG6_new/v_unisky_NS.txt",
#  				path + "FIG6_new/v_unisky_B.txt"]#,
# 				# path + "FIG6_new/v_unisky_albedo_.txt"]
# file_legends = ["Non pola", "East", "North", "B", "albedo"]
# # file_names = [	path + "uniform_sky_albedo_No_Directex1+uniform_sky_Directe_only_no_polax1",
# # 				path + "sky_only_no_albedo.txt",
# #  				path + "uniform_sky_albedo_No_Directex1+direct_only_Btruex1",
# #  				# path + "uniform_sky_albedo_No_Directex1+uniform_sky_Directe_only_Buen_100_max1x1",
# # 				path + "uniform_sky_albedo_No_Directex1+uniform_sky_Directe_only_Buen_010_max1x1",
# # 				path + "uniform_sky_albedo_No_Directex1+uniform_sky_Directe_only_Buen_001_max1x1"]
# # file_legends = ["Non pola", "No albedo", "Up", "East", "North"]
# direct_file = path + "FIG6_new/uniform_sky_Directe_only_no_pola.txt"
# direct_data = pd.read_csv(direct_file)

##FIGURE 3
path = "/home/bossel/These/Documentation/Mes_articles/scattering/new_pictures/FIG3/"
file_names = [path + "FIG3_synth_aurora_NOdirect.txt", path + "FIG3_synth_aurora_alb_NOdirect_0.2_0.5_0.5_1_0.8.txt", path + "FIG3_synth_aurora_alb_NOdirect_0.1_0.5_0.3_1_0.5.txt", path + "FIG3_synth_aurora_alb_NOdirect_0.5_0.5_0.7_1_0.9.txt"]#, path + "FIG3_synth_aurora_alb_NOdirect_1.txt", path + "FIG3_synth_aurora_alb_NOdirect_0.5.txt"]
file_legends = ["Without snow", "Varying albedo", "Low albedo", "High albedo", "Constant albedo=1", "Constant albedo=0.5"]
direct_file = path + "FIG3_synth_aurora_100_bkgrd_10_DIRECTonly.txt"
direct_data = pd.read_csv(direct_file)

file_colors = ["black", "red", 	"green", 	"blue", "black", 	"purple", 	"grey"]
file_styles = ["-", 	"-",	"-", 		"-", 	"--", 		"-", 		"-"]


fig, axs = plt.subplots(3, sharex=True)

average_I0 = np.average(direct_data["I0"])
# average_I0 = np.average(pd.read_csv(file_names[0])["I0"])


for i, files in enumerate(file_names):

	data = pd.read_csv(files)

	if "azimuts" in data.columns:
		x_type = "azimuts"
		x_label = "Azimuth (°)"
		x_unit = "deg"
	elif "datetime" in data.columns:
		x_type = "datetime"
		x_label = "Time"
		x_unit = ""


	axs[0].plot(data[x_type], data["I0"] / 1, linestyle = file_styles[i], color = file_colors[i], label = file_legends[i])
	# axs[0].plot(data[x_type], data["I0"] / average_I0, linestyle = file_styles[i], color = file_colors[i], label = file_legends[i])
	axs[1].plot(data[x_type], data["DoLP"], linestyle = file_styles[i], color = file_colors[i], label = file_legends[i])
	# axs[0].plot(data[x_type], (data["I0"]+direct_data["I0"]) / average_I0, "-", color = file_colors[i], label = file_legends[i])
	# axs[1].plot(data[x_type], data["DoLP"] * data["I0"]/(data["I0"]+direct_data["I0"]), "-", color = file_colors[i], label = file_legends[i])

	axs[2].plot(data[x_type], data["AoRD"], linestyle = file_styles[i], color = file_colors[i], label = file_legends[i])


axs[0].set_ylabel("Flux (nW)")
axs[0].ticklabel_format(axis = "y", style = "sci")
# axs[0].set_ylim(1.0105, 1.042)
axs[0].set_ylim(1.52e-4, 1.74e-4)
axs[1].set_ylabel("DoLP (%)")
axs[2].set_ylabel("AoLP (°)")
axs[2].set_xlabel(x_label)

# axs[0].legend(prop={'size': 16})
# axs[1].legend(prop={'size': 16})
# axs[2].legend(prop={'size': 16})
fig.tight_layout()

fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0)

plt.show()
