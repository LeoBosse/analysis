#!/usr/bin/python3
# -*-coding:utf-8 -*

try:
	from Geometry.Leo.src.observation import *
except:
	from observation import *

from subprocess import call
import numpy as np
import math as m
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib import cm
import sys as sys

import datetime as dt

import astropy.coordinates as coord
from astropy.time import Time
import astropy.units as u

import argparse


DtoR = np.pi/ 180	#Convert Degree to Radian
RtoD = 180. / np.pi	#Convert Radians to Degree

RT = 6367 #km #Earth's radius

# Default parameters
path = "/home/bossel/These/Analysis/src/geomagic"
data_out_path = path + "/results/data"
images_out_path = path + "/results/images"
outfile_name = "/eta_chaos"
mode = "fixed_e"

azimuth = np.arange(0, 360*DtoR, 1*DtoR)
elevation = 45 * DtoR

N = 10000
RD_src_azimut = None
RD_src_elevation = None


def UENToAzEl(uen):
	'''Converts up-east-north coordinates to azimuth and elevation angles.'''
	u, e, n = np.array(uen) / np.sqrt(np.sum([d**2 for d in uen]))
	
	az = np.arctan2(e, n)
	el = np.arcsin(u)
	
	return az, el


def GetSunAltAz(lon, lat, datetime):
	
	loc = coord.EarthLocation(lon=lon * u.rad,
							  lat=lat * u.rad)
	time = Time(datetime, format='datetime', scale='utc')

	altaz = coord.AltAz(location=loc, obstime=time)
	sun = coord.get_sun(time).transform_to(altaz)
	# print(sun.transform_to(altaz).alt)

	return sun.alt.value, sun.az.value


def GetLonLatFromName(name, unit="radians"):
	name = name.lower()
	# print("DEBUG LOCATION", name)
	if name == "mens":
		A_lon = 5.76
		A_lat = 44.83
	elif name == "skibotn":
		A_lon = 20.363109
		A_lat = 69.348027
	elif name == "nyalesund":
		A_lon = 11.92288
		A_lat = 78.92320
	elif name == "corbel":
		A_lon = 12.1197763
		A_lat = 78.9035015
	elif name == "vigan":
		A_lon = 3.504259
		A_lat = 44.039661
	elif name == "lagorge":
		A_lon = 5.936935
		A_lat = 45.212343
	elif name == "stveran":
		A_lon = 6.5430
		A_lat = 44.4156
	elif name == "sob":
		A_lon = -16.452403
		A_lat = 14.49496
	elif name == "eiscat_tromso":
		A_lon = 19.2147
		A_lat = 69.5839
	elif name == "eiscat_sodankyla":
		A_lon = 26.630207
		A_lat = 67.364918
	elif name == "eiscat_kiruna":
		A_lon = 20.4341
		A_lat = 67.8607
	elif name == "kilpisjarvi":
		A_lon = 20.7846926
		A_lat = 69.0526675
	elif name == "skibotnnord":
		A_lon = 20.385487
		A_lat = 69.54756
	elif name == "skibotnsud":
		A_lon = 19.981116
		A_lat = 69.234956
	elif name == "finland":
		A_lon = 20.935396
		A_lat = 68.924606
	elif name == "lacsud":
		A_lon = 19.05496
		A_lat = 69.13801
	# elif name == "helligskogen":
	# 	A_lon = 20.699313
	# 	A_lat = 69.206369
	# elif name == "kilpisjarvi":
	# 	A_lon = 20.781710
	# 	A_lat = 69.043909
	# elif name == "kafjorddalen":
	# 	A_lon = 20.935186
	# 	A_lat = 69.461868
	# elif name == "oteren":
	# 	A_lon = 19.880466
	# 	A_lat = 69.256451
	# elif name == "lokvollen":
	# 	A_lon = 20.548937
	# 	A_lat = 69.541251

	else:
		try:
			A_lon = float(name.split(";")[0])
			A_lat = float(name.split(";")[1])
		except:
			A_lon, A_lat = 0, 0

	if unit == "radians":
		A_lon *= DtoR
		A_lat *= DtoR

	return A_lon, A_lat
	

if __name__ == "__main__":
	nb_args = len(sys.argv)
	arguments = sys.argv


	arg_parser = argparse.ArgumentParser(
                    #prog='Vendange',
                    description='Geometry computations.',
					)

	arg_parser.add_argument('-l', '--location', default='0;0', help='location of observation. Name or lon;lat pair')        # positional argument
	arg_parser.add_argument('-alt', '--altitude', default = 110,  help="Altitude of the observed emission")
	arg_parser.add_argument('-az', '--azimut', default = 0,  help='Azimut of the observation.')
	arg_parser.add_argument('-el', '--elevation', default = 0,  help='Elevation of the observation.')
	arg_parser.add_argument('-src_el', '--source_elevation', default = None,  help='Elevation of the pollution source.')
	arg_parser.add_argument('-src_az', '--source_azimut', default = None,  help='Azimut of the pollution source.')
	arg_parser.add_argument('-dt', '--datetime', default = None,  help='Date and time of observation in UTC (YYYYMMDD-HH:MM:SS)')
	
	args = arg_parser.parse_args()

	A_lon, A_lat = GetLonLatFromName(args.location)
	outfile_name += "_" + args.location.lower()

	if args.datetime:
		obs_time = dt.datetime.strptime(args.datetime, "%Y%m%d-%H:%M:%S")
	else:
		obs_time = None

	h = float(args.altitude)
	outfile_name += "_h" + str(h)

	if args.elevation.lower() in ["rot", ''] and args.azimut.lower() in ["rot", '']:
		mode = "map"
		outfile_name += "_map"
	elif args.azimut.lower() in ["rot", '']:
		mode = "fixed_e"
		elevation = float(args.elevation) * DtoR
		outfile_name += "_e" + args.elevation
	elif args.elevation.lower() in ["rot", '']:
		mode = "fixed_a"
		azimuth = float(args.azimut) * DtoR
		outfile_name += "_a" + args.azimut
	else:
		azimuth = float(args.azimut) * DtoR
		elevation = float(args.elevation) * DtoR
		mode = "fixed_obs"
		outfile_name += "_a" + args.azimut
		outfile_name += "_e" + args.elevation

	if args.source_azimut and args.source_elevation:
		RD_src_azimut = float(args.source_azimut) * DtoR
		RD_src_elevation = float(args.source_azimut) * DtoR


	# entries = []

	# if nb_args > 1:
	# 	entries.append(arguments[1])
	# else:
	# 	entries.append(input("Enter the location of observation. Name or lon;lat pair: "))

	# A_lon, A_lat = GetLonLatFromName(entries[0])
	# outfile_name += "_" + entries[0].lower()
	
	

	# if nb_args > 2:
	# 	entries.append(arguments[2])
	# else:
	# 	entries.append(input("Enter the altitude of observation (km): "))
	# h = float(entries[1])
	# outfile_name += "_h" + entries[1]

	# if nb_args > 3:
	# 	entries.append(arguments[3])
	# else:
	# 	entries.append(input("Enter the azimut of observation, or the fixed parameter name (a, e) for rotations: "))
		
	# if entries[2] == "a":
	# 	mode = "fixed_a"
	# 	if nb_args > 4:
	# 		entries.append(arguments[4])
	# 	else:
	# 		entries.append(input("Enter the fixed azimut for the rotation: "))
	# 	azimuth = float(entries[3]) * DtoR
	# 	outfile_name += "_a" + entries[3]
	# elif entries[2] == "e":
	# 	mode = "fixed_e"
	# 	if nb_args > 4:
	# 		entries.append(arguments[4])
	# 	else:
	# 		entries.append(input("Enter the fixed elevation for the rotation: "))
	# 	elevation = float(entries[3]) * DtoR
	# 	outfile_name += "_e" + entries[3]
	# elif entries[2] == "map":
	# 	mode = "map"
	# 	if nb_args > 4:
	# 		entries.append(arguments[4])
	# 	else:
	# 		entries.append(input("Optional: Enter the Rayleigh Diff source azimut: "))
	# 	if entries[3]:
	# 		RD_src_azimut = float(entries[3]) * DtoR
	# 	if nb_args > 5:
	# 		entries.append(arguments[5])
	# 	else:
	# 		entries.append(input("Optional: Enter the Rayleigh Diff source elevation: "))
	# 	if entries[4]:
	# 		RD_src_elevation = float(entries[4]) * DtoR
	# else:
	# 	azimuth = float(entries[2]) * DtoR
	# 	mode = "fixed_obs"
	# 	if nb_args > 4:
	# 		entries.append(arguments[4])
	# 	else:
	# 		entries.append(input("Enter the elevation: "))
	# 	elevation = float(entries[3]) * DtoR
	# 	outfile_name += "_e" + entries[3]


	# if RD_src_azimut is None:
	# 	if nb_args > 5:
	# 		entries.append(arguments[5])
	# 	else:
	# 		entries.append(input("Optional: Enter the Rayleigh Diff source azimut: "))
	# 	if entries[4]:
	# 		RD_src_azimut = float(entries[4]) * DtoR

	# if nb_args > 6:
	# 	entries.append(arguments[6])
	# 	print("test")
	# else:
	# 	entries.append(input("Optional: Enter the Rayleigh Diff source elevation: "))
	# if entries[5]:
		# RD_src_elevation = float(entries[5]) * DtoR

	

	if mode == "fixed_obs":
		obs = ObservationPoint(A_lon, A_lat, h, azimuth, elevation, RD_src_azimut, RD_src_elevation)
		obs.SinglePointGeometry()
		obs.PrintAllParameters()
		
		if obs_time:
			sun_alt_I, sun_az_I = GetSunAltAz(obs.lon, obs.lat, obs_time)
			sun_alt_O, sun_az_O = GetSunAltAz(obs.P_lon, obs.P_lat, obs_time)
			print('Sun alt az (instrument): ', sun_alt_I, sun_az_I)
			print('Sun alt az (observation): ', sun_alt_O, sun_az_O)
			
			sol_zenith_angle = lambda sun_alt: (90 - sun_alt)*DtoR
			sol_zenith_angle_I = sol_zenith_angle(sun_alt_I)
			sol_zenith_angle_O = sol_zenith_angle(sun_alt_O)
			
			# print('Solar Zenith angle (instrument): ', sol_zenith_angle_I)
			# print('Solar Zenith angle (observation): ', sol_zenith_angle_O)
			
			el_margin = 0 * DtoR   # Count as illuminated if the sun is below the horizon by a given margin.
			
			illumination_height = lambda zen_angle, margin=0: RT * (np.cos(margin) / np.sin(zen_angle - margin) - 1)
			illumination_height_I = illumination_height(sol_zenith_angle_I, DtoR)
			illumination_height_O = illumination_height(sol_zenith_angle_O, DtoR)
			
			print(f'Illumination height above instrument (km) with a {el_margin*RtoD} degree margin: ', illumination_height_I)
			print(f'Illumination height above observation (km) with a {el_margin*RtoD} degree margin: ', illumination_height_O)
			
			
			ranges = range(0, 5000, 50)
			longs, lats, alts = zip(*[obs.GetPCoordinatesFromRange(r) for r in ranges])
			illumination_heights = np.array([illumination_height(sol_zenith_angle(GetSunAltAz(lon, lat, obs_time)[0]), el_margin) for lon, lat in zip(longs, lats)])
			illumination_index  = np.where(alts > illumination_heights)[0][0]
			illumination_alt = alts[illumination_index]
			print(f"Line of sight is illuminated by the sun at altitude {illumination_alt}, longitude {longs[illumination_index]*RtoD}, latitude {lats[illumination_index]*RtoD}, (with a {el_margin*RtoD} degree margin)")
			
			el_margin = 10 * DtoR 
			illumination_heights_sca = np.array([illumination_height(sol_zenith_angle(GetSunAltAz(lon, lat, obs_time)[0]), el_margin) for lon, lat in zip(longs, lats)])
			illumination_index_sca  = np.where(alts > illumination_heights)[0][0]
			illumination_alt_sca = alts[illumination_index]
			print(f"Line of sight is illuminated by the sun at altitude {illumination_alt_sca}, longitude {longs[illumination_index]*RtoD}, latitude {lats[illumination_index]*RtoD}, (with a {el_margin*RtoD} degree margin)")
			
				
			f, axs = plt.subplots(3, sharex= True)
			axs[0].set_ylabel('Longitude')
			axs[0].plot(ranges, np.array(longs)*RtoD)
			axs[1].set_ylabel('Latitude')
			axs[1].plot(ranges, np.array(lats)*RtoD)
			axs[2].set_ylabel('Altitude (km)')
			axs[2].plot(ranges, np.array(alts))
			axs[2].plot(ranges, illumination_heights, 'r')
			axs[2].plot(ranges, illumination_heights_sca, 'y')
			
			# axs[2].add_patch(patches.Rectangle((illumination_range, alts[0]), ranges[-1] -illumination_range, alts[-1], color="yellow"))
			
			axs[-1].set_xlabel('Range (km)')
			plt.show()

	#If we do several observations and want to plot the results.
	else:
		#We fix the azimuth, and look from the horyzon to the zenith
		if mode == "fixed_a":
			obs=[]
			elevations = np.linspace(0, 180*DtoR, N)
			chaos_input_list=[]
			#Looping over all positions, each get its own observation point object and the output values are plotted
			for i, e in enumerate(elevations):
				obs.append(ObservationPoint(A_lon, A_lat, h, azimuth, e, RD_src_azimut, RD_src_elevation))
				chaos_input_list.append([obs[-1].P_colat, obs[-1].P_lon])
			print("choas input list:", len(chaos_input_list))
			np.savetxt("theta_phi_H.dat", chaos_input_list, fmt="%1.8e")
			f2 = open(path + "/src/chaos_auto.out", "w")
			call([path + "/src/chaos_auto10000"], stdout = f2)
			f2.close()
			B_chaos = np.loadtxt(path + "/src/Bxyz_H_chaos6.dat")
			print("B_chaos:", len(B_chaos))
			for i in range(len(obs)):
				obs[i].B_chaos = [B_chaos[i][0], B_chaos[i][1], -B_chaos[i][4], B_chaos[i][3], B_chaos[i][2]]
				obs[i].SinglePointGeometry(GetBatP=False)

			ax1 = plt.subplot(111) #Initiating figure
			ax1.set_xlabel("elevation")
			ax1.set_ylabel("Eta, a=" + str(azimuth*RtoD) + " deg")
			ax1.plot(elevations*RtoD, [o.eta_chaos*RtoD for o in obs], color = "blue", label = "AoBapp using chaos")
			ax1.plot(elevations*RtoD, [o.Blos*RtoD for o in obs], color = "red", label = "AoBlos using chaos")
			if RD_src_azimut:
				ax1.plot(elevations*RtoD, [o.AoRD*RtoD for o in obs], color = "green", label = "AoRD")
			np.savetxt(data_out_path + outfile_name + ".dat", [[o.e*RtoD, o.eta_chaos*RtoD] for o in obs], header = "elevation (deg) eta (deg)")
		#	np.savetxt(outfile_name + "Bchaos.dat", [o.B_chaos[2:] for o in obs], header="Bup Beast Bnorth")

		#We fix the elevation, and look all around, from azimuth 0 to 360 North-East-South-West-North
		elif mode == "fixed_e":
			obs=[]
			azimuth = np.linspace(0, 360*DtoR, N)
			chaos_input_list=[]
			#Looping over all positions, each get its own observation point object and the output values are plotted
			for i, a in enumerate(azimuth):
				obs.append(ObservationPoint(A_lon, A_lat, h, a, elevation, RD_src_azimut, RD_src_elevation))
				chaos_input_list.append([obs[-1].P_colat, obs[-1].P_lon])
			print("choas input list:", len(chaos_input_list))
			np.savetxt("theta_phi_H.dat", chaos_input_list, fmt="%1.8e")
			f2 = open(path + "/src/chaos_auto.out", "w")
			call([path + "/src/chaos_auto10000"], stdout = f2)
			f2.close()
			B_chaos = np.loadtxt(path + "/src/Bxyz_H_chaos6.dat")
			print("B_chaos:", len(B_chaos))
			for i in range(len(obs)):
				obs[i].B_chaos = [B_chaos[i][0], B_chaos[i][1], -B_chaos[i][4], B_chaos[i][3], B_chaos[i][2]]
				obs[i].SinglePointGeometry(GetBatP=False)

			ax1 = plt.subplot(111) #Initiating figure
			ax1.set_xlabel("azimuth")
			ax1.set_ylabel("Eta, e=" + str(elevation*RtoD) + " deg")

			ax1.plot(azimuth*RtoD, [(o.eta_chaos*RtoD) for o in obs], color = "blue", label = "AoBapp")
			# ax1.plot(azimuth*RtoD, [(o.eta_chaos*RtoD)%180 for o in obs], color = "blue", label = "AoBapp")

			# ax1.plot(azimuth*RtoD, [o.Blos*RtoD for o in obs], color = "red", label = "AoBlos")
			# ax1.plot(azimuth*RtoD, [90]*len(obs), color="black")
			if RD_src_azimut:
				ax1.plot(azimuth*RtoD, [(o.AoRD*RtoD) for o in obs], color = "green", label = "AoRD")
				# ax1.plot(azimuth*RtoD, [(o.AoRD*RtoD)%180 for o in obs], color = "green", label = "AoRD")

				# ax1.plot(azimuth*RtoD, [(o.AoRD-o.eta_chaos)*RtoD for o in obs], color = "green", label = "using chaos")
			np.savetxt(data_out_path + outfile_name + ".dat", [[o.a*RtoD, o.eta_chaos*RtoD] for o in obs], header = "azimuth (deg) eta (deg)")
		#	np.savetxt(outfile_name + "Bchaos.dat", [o.B_chaos[2:] for o in obs], header="Bup Beast Bnorth")
		elif mode == "map":
			elevations = np.linspace(0, 90*DtoR, int(np.sqrt(N)/2))
			azimuth = np.linspace(0, 360*DtoR, int(2*np.sqrt(N)))
			obs = []
			AoBapp = np.zeros((len(elevations), len(azimuth)))
			AoRD = np.zeros(AoBapp.shape)
			AoBlos = np.zeros(AoBapp.shape)
			diff = np.zeros(AoBapp.shape)
			chaos_input_list = []
			#Looping over all positions, each get its own observation point object and the output values are plotted
			for i, e in enumerate(elevations):
				obs.append([])
				for j, a in enumerate(azimuth):
					obs[i].append(ObservationPoint(A_lon, A_lat, h, a, e, RD_src_azimut, RD_src_elevation))
					chaos_input_list.append([obs[-1][-1].P_colat, obs[-1][-1].P_lon])
			print("choas input list:", len(chaos_input_list))
			np.savetxt("theta_phi_H.dat", chaos_input_list, fmt="%1.8e")
			f2 = open(path + "/src/chaos_auto.out", "w")
			call([path + "/src/chaos_auto10000"], stdout = f2)
			f2.close()
			B_chaos = np.loadtxt(path + "/src/Bxyz_H_chaos6.dat")
			print("B_chaos:", len(B_chaos))
			for i, e in enumerate(elevations):
				for j, a in enumerate(azimuth):
					B_index = i*(len(azimuth)-1) + j
					obs[i][j].B_chaos = [B_chaos[B_index][0], B_chaos[B_index][1], -B_chaos[B_index][4], B_chaos[B_index][3], B_chaos[i][2]]
					obs[i][j].SinglePointGeometry(GetBatP=False)

					AoBapp[i][j] = obs[i][j].eta_chaos
					AoBlos[i][j] = int(89*DtoR <= obs[i][j].Blos <= 91*DtoR)
					AoRD[i][j] = obs[i][j].AoRD
					diff[i][j] = abs(AoBapp[i][j] - AoRD[i][j])
					if diff[i][j] > 90 * DtoR:
						diff[i][j] = 180 * DtoR - diff[i][j]


			f1, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, figsize=(16, 8))
			ax3.set_xlabel("azimuth")
			ax1.set_ylabel("elevation")
			ax2.set_ylabel("elevation")
			ax3.set_ylabel("elevation")

			i1 = ax1.imshow(AoBapp * RtoD, extent=[0, 360, 0, 90], origin = "lower", cmap=cm.twilight)
			i2 = ax2.imshow(AoBlos, extent=[0, 360, 0, 90], origin = "lower")
			i3 = ax3.imshow(diff * RtoD, extent=[0, 360, 0, 90], origin = "lower", cmap=cm.twilight)
			cbar1 = f1.colorbar(i1, extend='both', spacing='proportional', shrink=0.9, ax=ax1)
			cbar1.set_label('AoBapp')
			cbar2 = f1.colorbar(i2, extend='both', spacing='proportional', shrink=0.9, ax=ax2)
			cbar2.set_label('AoBlos')
			cbar3 = f1.colorbar(i3, extend='both', spacing='proportional', shrink=0.9, ax=ax3)
			cbar3.set_label('abs(AoBapp-AoRD)')

			f1.suptitle("Relevant angles map at skibotn, with light pollution source at -45deg in azimut")
			f1.subplots_adjust(hspace=0)

			# test = np.zeros((90, 360))
			# for i in range(90):
			# 	for j in range(360):
			# 		test[i][j] = j
			# i1 = plt.imshow(test, extent=[0, 360, 0, 90], origin = "lower")
			# cbar1 = plt.colorbar(i1, extend='both', spacing='proportional', shrink=0.9)


		plt.legend()

		plt.savefig(images_out_path + outfile_name + '.png')

		plt.show()




class Geometry:
	def __init__(self, *arguments):

		self.N = 50
		self.RD_src_azimut = False
		self.RD_src_elevation = 0

		nb_args = len(arguments)

		entries = []

		if nb_args > 1:
			entries.append(arguments[1])
		else:
			entries.append(input("Enter the location of observation. Name or lon;lat pair: "))


		self.outfile_name = "/eta_chaos"
		self.A_lon, self.A_lat = GetLonLatFromName(entries[0].lower())
		self.outfile_name += "_" + entries[0].lower()
		if nb_args > 2:
			entries.append(arguments[2])
		else:
			entries.append(input("Enter the altitude of observation (km): "))
		self.h = float(entries[1])
		self.outfile_name += "_h" + entries[1]

		if nb_args > 3:
			entries.append(arguments[3])
		else:
			entries.append(input("Enter the azimut of observation, or the fixed parameter name (a, e) for rotation: "))
		if entries[2] == "a":
			self.mode = "fixed_a"
			if nb_args > 4:
				entries.append(arguments[4])
			else:
				entries.append(input("Enter the fixed azimut for the rotation: "))
			self.azimuth = float(entries[3]) * DtoR
			self.outfile_name += "_a" + entries[3]
		elif entries[2] == "e":
			self.mode = "fixed_e"
			if nb_args > 4:
				entries.append(arguments[4])
			else:
				entries.append(input("Enter the fixed elevation for the rotation: "))
			self.elevation = float(entries[3]) * DtoR
			self.outfile_name += "_e" + entries[3]
		elif entries[2] == "map":
			self.mode = "map"
			if nb_args > 3:
				entries.append(arguments[4])
			else:
				entries.append(input("Optional: Enter the Rayleigh Diff source azimut: "))
			if entries[3]:
				self.RD_src_azimut = float(entries[3]) * DtoR
			if nb_args > 4:
				entries.append(arguments[5])
			else:
				entries.append(input("Optional: Enter the Rayleigh Diff source elevation: "))
			if entries[4]:
				self.RD_src_elevation = float(entries[4]) * DtoR
		else:
			self.azimuth = float(entries[2]) * DtoR
			self.mode = "fixed_obs"
			if nb_args > 4:
				entries.append(arguments[4])
			else:
				entries.append(input("Enter the elevation: "))
			self.elevation = float(entries[3]) * DtoR
			self.outfile_name += "_e" + entries[3]

		if self.RD_src_azimut is False:
			if nb_args > 5:
				entries.append(arguments[5])
			else:
				entries.append(input("Optional: Enter the Rayleigh Diff source azimut: "))
			if entries[4]:
				self.RD_src_azimut = float(entries[4]) * DtoR

		if nb_args > 6:
			entries.append(arguments[6])
			print("test")
		else:
			entries.append(input("Optional: Enter the Rayleigh Diff source elevation: "))
		if entries[5]:
			self.RD_src_elevation = float(entries[5]) * DtoR


	def FixedObservation(self, B_model = None):
		obs = ObservationPoint(self.A_lon, self.A_lat, self.h, self.azimuth, self.elevation, self.RD_src_azimut, self.RD_src_elevation, init_full=False)
		obs.SinglePointGeometry(B_model = B_model)
		obs.PrintAllParameters()


	# #If we do several observations and want to plot the results.
	# def FixedAzimut(self):
	# 	#We fix the azimuth, and look from the horyzon to the zenith
	# 	self.obs=[]
	# 	self.elevations = np.linspace(0, 90*DtoR, self.N)
	# 	self.chaos_input_list=[]
	# 	#Looping over all positions, each get its own observation point object and the output values are plotted
	# 	for i, e in enumerate(self.elevations):
	# 		obs.append(ObservationPoint(self.A_lon, self.A_lat, self.h, self.azimuth, e, self.RD_src_azimut, self.RD_src_elevation))
	# 		chaos_input_list.append([self.obs[-1].P_colat, self.obs[-1].P_lon])
	# 	print("choas input list:", len(self.chaos_input_list))
	# 	np.savetxt("theta_phi_H.dat", self.chaos_input_list, fmt="%1.8e")
	# 	f2 = open(path + "/src/chaos_auto.out", "w")
	# 	call([path + "/src/chaos_auto10000"], stdout = f2)
	# 	f2.close()
	# 	self.B_chaos = np.loadtxt(path + "/src/Bxyz_H_chaos6.dat")
	# 	print("B_chaos:", len(self.B_chaos))
	# 	for i in range(len(self.obs)):
	# 		self.obs[i].B_chaos = [self.B_chaos[i][0], self.B_chaos[i][1], -self.B_chaos[i][4], self.B_chaos[i][3], self.B_chaos[i][2]]
	# 		self.obs[i].SinglePointGeometry(GetBatP=False)
	#
	# 	ax1 = plt.subplot(111) #Initiating figure
	# 	ax1.set_xlabel("elevation")
	# 	ax1.set_ylabel("Eta, a=" + str(azimuth*RtoD) + " deg")
	# 	ax1.plot(elevations*RtoD, [o.eta_chaos*RtoD for o in obs], color = "blue", label = "using chaos")
	# 	# ax1.plot(elevations*RtoD, [o.Blos*RtoD for o in obs], color = "blue", label = "using chaos")
	# 	if RD_src_azimut:
	# 		ax1.plot(elevations*RtoD, [o.AoRD*RtoD for o in obs], color = "green", label = "AoRD")
	# 	np.savetxt(data_out_path + outfile_name + ".dat", [[o.e*RtoD, o.eta_chaos*RtoD] for o in obs], header = "elevation (deg) eta (deg)")
	# #	np.savetxt(outfile_name + "Bchaos.dat", [o.B_chaos[2:] for o in obs], header="Bup Beast Bnorth")


	def FixedElevation(self, direction="NESW", B_model = None):
		#We fix the elevation, and look all around, from azimuth 0 to 360 North-East-South-West-North
		self.obs=[]

		if direction == "NESW": self.azimuth = np.linspace(0, 360*DtoR, self.N)
		elif direction == "NWSE": self.azimuth = np.linspace(360*DtoR, 0, self.N)
		self.chaos_input_list=[]
		#Looping over all positions, each get its own observation point object and the output values are plotted
		for i, a in enumerate(self.azimuth):
			self.obs.append(ObservationPoint(self.A_lon, self.A_lat, self.h, a, self.elevation, self.RD_src_azimut, self.RD_src_elevation, init_full=False))
			self.obs[-1].SinglePointGeometry(B_model = B_model)
			# self.chaos_input_list.append([self.obs[-1].P_colat, self.obs[-1].P_lon])
		# print("choas input list:", len(self.chaos_input_list))
		# np.savetxt(path + "/B_model/theta_phi_H.dat", self.chaos_input_list, fmt="%1.8e")
		# f2 = open(path + "/B_model/chaos_auto.out", "w")
		# call([path + "/B_model/chaos_auto" + str(self.N)], stdout = f2)
		# f2.close()
		# print("/B_model/chaos_auto" + str(self.N))
		# self.B_chaos = np.loadtxt(path + "/B_model/Bxyz_H_chaos6.dat")
		# print("B_chaos:", len(self.obs), len(self.B_chaos))

		# time = chaos.data_utils.mjd2000(2019, 1, 1)
		# self.B_chaos = B_model.synth_values_tdep(time, RT + self.h, self.P_colat*RtoD, self.P_lon*RtoD) #given in up, south, east

		# for i in range(len(self.obs)):
		# 	self.obs[i].B_chaos = [self.B_chaos[i][0], self.B_chaos[i][1], -self.B_chaos[i][4], self.B_chaos[i][3], self.B_chaos[i][2]]
			# self.obs[i].SinglePointGeometry(GetBatP=False)


		self.x_axis = self.azimuth*RtoD
		self.x_label = "azimuth"
		self.y_label = "Eta, e=" + str(self.elevation*RtoD) + " deg"

		return self.x_axis, self.obs


	def Graph(self):
		ax1 = plt.subplot(111) #Initiating figure
		ax1.set_xlabel(self.x_label)
		ax1.set_ylabel(self.y_label)

		ax1.plot(self.x_axis, [(o.eta_chaos*RtoD) for o in self.obs], color = "blue", label = "AoBapp")
		# ax1.plot(self.x_axis, [(o.eta_chaos*RtoD)%180 for o in obs], color = "blue", label = "AoBapp")
		# ax1.plot(self.x_axis, [o.Blos*RtoD for o in self.obs], color = "red", label = "AoBlos")

		if self.RD_src_azimut:
			ax1.plot(self.x_axis, [(o.AoRD*RtoD) for o in self.obs], color = "green", label = "AoRD")
			# ax1.plot(self.x_axis, [(o.AoRD*RtoD)%180 for o in obs], color = "green", label = "AoRD")

			# ax1.plot(self.x_axis, [(o.AoRD-o.eta_chaos)*RtoD for o in obs], color = "green", label = "using chaos")
		np.savetxt(data_out_path + self.outfile_name + ".dat", [[o.a*RtoD, o.eta_chaos*RtoD] for o in self.obs], header = "azimuth (deg) eta (deg)")
	#	np.savetxt(outfile_name + "Bchaos.dat", [o.B_chaos[2:] for o in obs], header="Bup Beast Bnorth")

		plt.legend()
		plt.savefig(images_out_path + self.outfile_name + '.png')
		plt.show()

	# def Map(self):
	# 	elevations = np.linspace(0, 90*DtoR, int(np.sqrt(N)/2))
	# 	azimuth = np.linspace(0, 360*DtoR, int(2*np.sqrt(N)))
	# 	obs = []
	# 	AoBapp = np.zeros((len(elevations), len(azimuth)))
	# 	AoRD = np.zeros(AoBapp.shape)
	# 	diff = np.zeros(AoBapp.shape)
	# 	chaos_input_list = []
	# 	#Looping over all positions, each get its own observation point object and the output values are plotted
	# 	for i, e in enumerate(elevations):
	# 		obs.append([])
	# 		for j, a in enumerate(azimuth):
	# 			obs[i].append(ObservationPoint(A_lon, A_lat, h, a, e, RD_src_azimut, RD_src_elevation))
	# 			chaos_input_list.append([obs[-1][-1].P_colat, obs[-1][-1].P_lon])
	# 	print("choas input list:", len(chaos_input_list))
	# 	np.savetxt("theta_phi_H.dat", chaos_input_list, fmt="%1.8e")
	# 	f2 = open(path + "/src/chaos_auto.out", "w")
	# 	call([path + "/src/chaos_auto10000"], stdout = f2)
	# 	f2.close()
	# 	B_chaos = np.loadtxt(path + "/src/Bxyz_H_chaos6.dat")
	# 	print("B_chaos:", len(B_chaos))
	# 	for i, e in enumerate(elevations):
	# 		for j, a in enumerate(azimuth):
	# 			B_index = i*(len(azimuth)-1) + j
	# 			obs[i][j].B_chaos = [B_chaos[B_index][0], B_chaos[B_index][1], -B_chaos[B_index][4], B_chaos[B_index][3], B_chaos[i][2]]
	# 			obs[i][j].SinglePointGeometry(GetBatP=False)
	#
	# 			AoBapp[i][j] = obs[i][j].eta_chaos
	# 			AoRD[i][j] = obs[i][j].AoRD
	# 			diff[i][j] = abs(AoBapp[i][j] - AoRD[i][j])
	# 			if diff[i][j] > 90 * DtoR:
	# 				diff[i][j] = 180 * DtoR - diff[i][j]
	#
	#
	# 	f1, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, figsize=(16, 8))
	# 	ax3.set_xlabel("azimuth")
	# 	ax1.set_ylabel("elevation")
	# 	ax2.set_ylabel("elevation")
	# 	ax3.set_ylabel("elevation")
	#
	# 	i1 = ax1.imshow(AoBapp * RtoD, extent=[0, 360, 0, 90], origin = "lower")
	# 	i2 = ax2.imshow(AoRD * RtoD, extent=[0, 360, 0, 90], origin = "lower")
	# 	i3 = ax3.imshow(diff * RtoD, extent=[0, 360, 0, 90], origin = "lower")
	# 	cbar1 = f1.colorbar(i1, extend='both', spacing='proportional', shrink=0.9, ax=ax1)
	# 	cbar1.set_label('AoBapp')
	# 	cbar2 = f1.colorbar(i2, extend='both', spacing='proportional', shrink=0.9, ax=ax2)
	# 	cbar2.set_label('AoRD')
	# 	cbar3 = f1.colorbar(i3, extend='both', spacing='proportional', shrink=0.9, ax=ax3)
	# 	cbar3.set_label('abs(AoBapp-AoRD)')
	#
	# 	f1.suptitle("Relevant angles map at skibotn, with light pollution source at -45deg in azimut")
	# 	f1.subplots_adjust(hspace=0)
	#
	# 	# test = np.zeros((90, 360))
	# 	# for i in range(90):
	# 	# 	for j in range(360):
	# 	# 		test[i][j] = j
	# 	# i1 = plt.imshow(test, extent=[0, 360, 0, 90], origin = "lower")
	# 	# cbar1 = plt.colorbar(i1, extend='both', spacing='proportional', shrink=0.9)
	#
	#
	# 	plt.legend()
	#
	# 	plt.savefig(images_out_path + outfile_name + '.png')
	#
	# 	plt.show()
