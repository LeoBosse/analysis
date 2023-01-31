#!/usr/bin/python3
# -*-coding:utf-8 -*

from mpi4py import MPI
mpi_comm = MPI.COMM_WORLD
mpi_rank = mpi_comm.Get_rank()
mpi_size = mpi_comm.Get_size()
mpi_name = mpi_comm.Get_name()

import sys as sys
import os as os
import numpy as np
import time as tm
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Circle
from matplotlib.patches import Ellipse
from matplotlib.patches import Arrow
# from matplotlib.lines import Line2D

# import osgeo.gdal as gdal
# gdal.UseExceptions()  # not required, but a good idea

import imageio

from observation import *
from rayleigh_utils import *
from sky_map import *
from ground_map import *
from atmosphere import *
from world import *
from input import *
from multiple_scattering import *
# from multiple_scattering_python import *


class Simulation:
	def __init__(self, in_dict):
		"""Initialize the rayleigh object. It will be able to calculate a whole bunch of things!"""

		self.in_dict = in_dict

		self.save_individual_plots = bool(int(in_dict["saving_graphs"]))
		self.save_path = in_dict["save_path"]
		self.save_name = self.save_path + in_dict["save_name"]

		self.use_GPU = in_dict["use_GPU"].lower() in ['true', '1', 'yes']

		mpi_comm.Barrier()
		# print(mpi_rank)
		if mpi_rank == 0 and self.save_individual_plots:
			# print(mpi_rank)
			confirm = input(f"Saving graph under the name: {self.save_name}. Press ENTER to confirm or any caracter to infirm, then change the input file save configuration.")
			if confirm != "":
				raise Exception("Go change the input file save configuration and rerun! I'm waiting for you.")

			if not os.path.exists(self.save_path):
				print(f"Creating saving folder {self.save_path}")
				os.makedirs(self.save_path)
				if not self.save_individual_plots:
					self.save_path = self.save_name = ""

		mpi_comm.Barrier()

		self.world = World(in_dict)

		self.cst_flux = float(in_dict["cst_Flux"]) #in nW
		self.cst_DoLP = float(in_dict["cst_DoLP"]) / 100 #between 0 and 1
		self.cst_AoLP = float(in_dict["cst_AoLP"]) * DtoR # in radians

		self.cst_V, self.cst_Vcos, self.cst_Vsin = self.world.GetVParamFromLightParam(self.cst_flux, self.cst_DoLP, self.cst_AoLP)

		self.add_starlight = in_dict["add_starlight"].lower() in ["1", "true", "yes", "ok"]

		self.direct_light_mode = in_dict["direct_mode"].lower()
		if self.direct_light_mode not in ["none", "only", "add"]:
			self.direct_light_mode = "add"

		if not self.world.is_single_observation and self.world.is_time_dependant:
			if mpi_rank==0: print("Cannot compute several observation directions with a time changing sky. Please, my code is already messy enough as it is...")
			raise SystemExit

		list_shapes = (self.world.sky_map.Nt, len(self.world.e_pc_list), len(self.world.a_pc_list))
		### Initialize lists to save the final observable of each observation
		self.I_direct_list 	= np.zeros(list_shapes)
		self.I_list 		= np.zeros(list_shapes)
		self.IPola_list 	= np.zeros(list_shapes)
		self.InonPola_list 	= np.zeros(list_shapes)
		self.DoLP_list 		= np.zeros(list_shapes)
		self.AoRD_list 		= np.zeros(list_shapes)

		self.I_list_err 	= np.zeros(list_shapes)
		self.InonPola_list_err = np.zeros(list_shapes)
		self.IPola_list_err = np.zeros(list_shapes)
		self.DoLP_list_err 	= np.zeros(list_shapes)
		self.AoRD_list_err 	= np.zeros(list_shapes)

		self.MS_I0			= np.zeros(list_shapes)
		self.MS_DoLP		= np.zeros(list_shapes)
		self.MS_AoLP		= np.zeros(list_shapes)



		self.add_B_pola = bool(float(in_dict["add_B_pola"]))
		self.add_B_Flux = np.zeros(list_shapes)
		self.add_B_DoLP = np.zeros(list_shapes)
		self.add_B_AoLP = np.zeros(list_shapes)

		if mpi_rank==0: print("Initialization DONE")

	# def GatherResults(root = 0):
	# 	self.world.GatherResults()


	def ComputeAllMaps(self):
		"""Will compute what the instrument sees for all the maps and all the observations directions"""


		if self.world.has_sky_emission: ### Compute how the sky maps looks through the instrument for every time and observations directions
			self.is_ground_emission = False
			observations_number = self.world.sky_map.Nt * self.world.Nb_a_pc * self.world.Nb_e_pc
			count = 0

			if self.use_GPU:
				self.world.CreateShaderWrap('sky', self.world.sky_map.cube[0, :, :])

			for t in range(self.world.sky_map.Nt):
				# self.SetIMap(time = t)
				for ia_pc in range(self.world.Nb_a_pc):
					for ie_pc in range(self.world.Nb_e_pc):
						count += 1
						if mpi_rank == 0:
							print("*******************************************************************************************************************************")
							print(f"Starting calculation for SKY map {count}/{observations_number}:")

						self.ia_pc, self.ie_pc = ia_pc, ie_pc
						self.world.a_pc, self.world.e_pc = self.world.a_pc_list[ia_pc], self.world.e_pc_list[ie_pc]
						self.a_pc, self.e_pc = self.world.a_pc, self.world.e_pc
						self.time = t
						self.SingleComputation()
			
			if self.use_GPU:
				self.world.shader.Release()

			reciever_V = None
			reciever_Vcos = None
			reciever_Vsin = None
			if mpi_rank == 0:
				print("Stacking processor results...")
				reciever_V = np.empty_like(self.world.sky_map.V_map)
				reciever_Vcos = np.empty_like(self.world.sky_map.Vcos_map)
				reciever_Vsin = np.empty_like(self.world.sky_map.Vsin_map)
			mpi_comm.Reduce(self.world.sky_map.V_map, reciever_V, op = MPI.SUM, root=0)
			self.world.sky_map.V_map = reciever_V
			mpi_comm.Reduce(self.world.sky_map.Vcos_map, reciever_Vcos, op = MPI.SUM, root=0)
			self.world.sky_map.Vcos_map = reciever_Vcos
			mpi_comm.Reduce(self.world.sky_map.Vsin_map, reciever_Vsin, op = MPI.SUM, root=0)
			self.world.sky_map.Vsin_map = reciever_Vsin

			if mpi_rank == 0:
				self.world.sky_map.SetLightParameters()

		if self.world.has_ground_emission: ### Compute how the ground map looks through the instrument for every time and observations directions
			observations_number = self.world.sky_map.Nt * self.world.Nb_a_pc * self.world.Nb_e_pc
			count = 0
			self.is_ground_emission = True

			if self.use_GPU:
				self.world.CreateShaderWrap('ground', self.world.ground_map.cube[0, :, :])

			for t in range(self.world.sky_map.Nt):
				for ia_pc in range(self.world.Nb_a_pc):
					for ie_pc in range(self.world.Nb_e_pc):
						count += 1
						if mpi_rank==0:
							print("*******************************************************************************************************************************")
							print(f"STARTING calculation for GROUND map {count}/{observations_number}:")

						self.ia_pc, self.ie_pc = ia_pc, ie_pc
						self.world.a_pc, self.world.e_pc = self.world.a_pc_list[ia_pc], self.world.e_pc_list[ie_pc]
						self.a_pc, self.e_pc = self.world.a_pc, self.world.e_pc
						if t == 0 or (t > 0 and self.world.ground_map.Nt > 0):
							self.time = t
							self.SingleComputation()
			if self.use_GPU:
				self.world.shader.Release()

			reciever_V = None
			reciever_Vcos = None
			reciever_Vsin = None
			reciever_shade = None
			if mpi_rank == 0:
				print("Stacking processor results...")
				reciever_V = np.empty_like(self.world.ground_map.V_map)
				reciever_Vcos =  np.empty_like(self.world.ground_map.Vcos_map)
				reciever_Vsin = np.empty_like(self.world.ground_map.Vsin_map)
				reciever_shade = np.empty_like(self.world.ground_map.mountain_shadow_heights)
			mpi_comm.Reduce(self.world.ground_map.V_map, reciever_V, op = MPI.SUM, root=0)
			self.world.ground_map.V_map = reciever_V
			mpi_comm.Reduce(self.world.ground_map.Vcos_map, reciever_Vcos, op = MPI.SUM, root=0)
			self.world.ground_map.Vcos_map = reciever_Vcos
			mpi_comm.Reduce(self.world.ground_map.Vsin_map, reciever_Vsin, op = MPI.SUM, root=0)
			self.world.ground_map.Vsin_map = reciever_Vsin
			mpi_comm.Reduce(self.world.ground_map.mountain_shadow_heights, reciever_shade, op = MPI.MIN, root=0)
			self.world.ground_map.mountain_shadow_heights = reciever_shade

			if mpi_rank == 0:
				self.world.ground_map.SetLightParameters()

		if mpi_rank == 0:
			self.GetLightParametersList()

	def SingleComputation(self):
		"""Will compute and show what the instrument sees for one given I_map and one observation"""
		self.ComputeMaps()
		# self.MakePlots()

	def ComputeMaps(self):
		"""Will initialize the last parameters that depend on the rest and will compute the contribution of each pixels from the map we set. We set the map by setting self.is_ground_emission to True or False. If False, then compute the sky map at the time set"""

		self.world.SetObservation(self.a_pc, self.e_pc)

		# self.SetSaveName()

		start_time = dt.datetime.now()
		if mpi_rank==0:
			print(self.world.obs)
			print("Starting time: ", start_time.time().strftime("%H:%M:%S"))

		if self.is_ground_emission == True:

			if not self.use_GPU:
				print('Ground map commputation with CPU...')
				self.world.ComputeGroundMaps(self.time, self.ia_pc, self.ie_pc)
				# print(self.world.ground_map.V_map[self.time, self.ie_pc, self.ia_pc, 0, :])
				# self.world.ground_map.Vcos_map[self.time, self.ie_pc, self.ia_pc, 0, :],
				# self.world.ground_map.Vsin_map[self.time, self.ie_pc, self.ia_pc, 0, :])
				# test_V = self.world.ground_map.V_map[self.time, self.ie_pc, self.ia_pc].copy()
				# test_Vcos = self.world.ground_map.Vcos_map[self.time, self.ie_pc, self.ia_pc].copy()
				# test_Vsin = self.world.ground_map.Vsin_map[self.time, self.ie_pc, self.ia_pc].copy()
				# test_I, test_D, test_A = np.vectorize(self.GetIDAFromV)(test_V, test_Vcos, test_Vsin)
				# self.world.ground_map.V_map[self.time, self.ie_pc, self.ia_pc]    *= 0
				# self.world.ground_map.Vsin_map[self.time, self.ie_pc, self.ia_pc] *= 0
				# self.world.ground_map.Vcos_map[self.time, self.ie_pc, self.ia_pc] *= 0

			else:
				print('Ground map computation with GPU...')
				self.world.ComputeGroundMapsGPU(self.time, self.ia_pc, self.ie_pc)
				# print('GPU / CPU')
				# print(self.world.ground_map.V_map[self.time, self.ie_pc, self.ia_pc], test_V)
				# print(np.average(test_V),
				# 	  np.average(test_Vcos),
				# 	  np.average(test_Vsin))
				# print(np.average(self.world.ground_map.V_map[   self.time, self.ie_pc, self.ia_pc]),
					#   np.average(self.world.ground_map.Vcos_map[self.time, self.ie_pc, self.ia_pc]),
					#   np.average(self.world.ground_map.Vsin_map[self.time, self.ie_pc, self.ia_pc]))
				# print(np.average(test_V),
					#   np.average(test_Vcos),
					#   np.average(test_Vsin))
				# print(np.average(test_I),
					#   np.average(test_D),
					#   np.average(test_A)*RtoD)
				# print(100*np.max(abs((self.world.ground_map.V_map[   self.time, self.ie_pc, self.ia_pc]) - test_V)),
					#   100*np.max(abs((self.world.ground_map.Vcos_map[self.time, self.ie_pc, self.ia_pc]) - test_Vcos)),
				#   100*np.max(abs((self.world.ground_map.Vsin_map[self.time, self.ie_pc, self.ia_pc]) - test_Vsin)))
				# I, D, A = np.vectorize(self.GetIDAFromV)(self.world.ground_map.V_map[   self.time, self.ie_pc, self.ia_pc], self.world.ground_map.Vcos_map[   self.time, self.ie_pc, self.ia_pc], self.world.ground_map.Vsin_map[   self.time, self.ie_pc, self.ia_pc])
				# print(np.average(I),
					#   np.average(D),
					#   np.average(A)*RtoD)
				# print(100*np.max(abs((I) - test_I)),
				# 	  100*np.max(abs((D) - test_D)),
				# 	  100*np.max(abs((A) - test_A)))
		else:
			if self.add_B_pola:
				# print("DEBUG DIRECT ONLY")
				B_DoLP, B_AoLP = self.world.GetPolaFromB(DoLP_max = 20)
				# print("Simu", B_DoLP, B_AoLP)
				# print(self.add_B_DoLP.shape, self.add_B_AoLP.shape)
				self.add_B_DoLP[self.time, self.ie_pc, self.ia_pc] = B_DoLP
				self.add_B_AoLP[self.time, self.ie_pc, self.ia_pc] = B_AoLP

			if self.direct_light_mode not in ["only"]:
				
				# if not self.use_GPU:
				print("Sky map CPU")
				self.world.ComputeSkyMaps(self.time, self.ia_pc, self.ie_pc)
				print(self.world.sky_map.V_map[self.time, self.ie_pc, self.ia_pc, 0, :],
					  self.world.sky_map.Vcos_map[self.time, self.ie_pc, self.ia_pc, 0, :],
					  self.world.sky_map.Vsin_map[self.time, self.ie_pc, self.ia_pc, 0, :])
				test_V = self.world.sky_map.V_map[self.time, self.ie_pc, self.ia_pc].copy()
				test_Vcos = self.world.sky_map.Vcos_map[self.time, self.ie_pc, self.ia_pc].copy()
				test_Vsin = self.world.sky_map.Vsin_map[self.time, self.ie_pc, self.ia_pc].copy()
				test_I, test_D, test_A = np.vectorize(self.GetIDAFromV)(test_V, test_Vcos, test_Vsin)
				self.world.sky_map.V_map[self.time, self.ie_pc, self.ia_pc]    *= 0
				self.world.sky_map.Vsin_map[self.time, self.ie_pc, self.ia_pc] *= 0
				self.world.sky_map.Vcos_map[self.time, self.ie_pc, self.ia_pc] *= 0
				# else:


				print("Sky map GPU")
				self.world.ComputeSkyMapsGPU(self.time, self.ia_pc, self.ie_pc)
				print('GPU / CPU')
				print(self.world.sky_map.V_map[self.time, self.ie_pc, self.ia_pc], test_V)
				print(np.average(test_V),
					  np.average(test_Vcos),
					  np.average(test_Vsin))
				print(np.average(self.world.sky_map.V_map[   self.time, self.ie_pc, self.ia_pc]),
					  np.average(self.world.sky_map.Vcos_map[self.time, self.ie_pc, self.ia_pc]),
					  np.average(self.world.sky_map.Vsin_map[self.time, self.ie_pc, self.ia_pc]))
				print(np.average(test_V),
					  np.average(test_Vcos),
					  np.average(test_Vsin))
				print(np.average(test_I),
					  np.average(test_D),
					  np.average(test_A)*RtoD)
				print(100*np.max(abs((self.world.sky_map.V_map[   self.time, self.ie_pc, self.ia_pc]) - test_V)),
					  100*np.max(abs((self.world.sky_map.Vcos_map[self.time, self.ie_pc, self.ia_pc]) - test_Vcos)),
				  100*np.max(abs((self.world.sky_map.Vsin_map[self.time, self.ie_pc, self.ia_pc]) - test_Vsin)))
				I, D, A = np.vectorize(self.GetIDAFromV)(self.world.sky_map.V_map[   self.time, self.ie_pc, self.ia_pc], self.world.sky_map.Vcos_map[   self.time, self.ie_pc, self.ia_pc], self.world.sky_map.Vsin_map[   self.time, self.ie_pc, self.ia_pc])
				print(np.average(I),
					  np.average(D),
					  np.average(A)*RtoD)
				print(100*np.max(abs((I) - test_I)),
					  100*np.max(abs((D) - test_D)),
					  100*np.max(abs((A) - test_A)))

		if mpi_rank==0:
			print("Computing DONE in: ", dt.datetime.now() - start_time)


	def ComputeMultipleScattering(self):
		if mpi_rank == 0: print("Computing multiple scattering")
		instrument = (self.a_pc, self.e_pc, 0)
		mul_sca = MultipleScattering(self.in_dict, instrument, self.world.atmosphere)

		# mul_sca.MakeScatteringHistograms(1, 0)
		# mul_sca.MakeScatteringHistograms(0, 1)
		# mul_sca.MakeScatteringHistograms(1, 1)

		# plt.show()

		mul_sca.PropagateAll()


		if mpi_rank == 0:
			mul_sca.SetRaysFluxList(self.world.ground_map)
			mul_sca.GetTotalUnpolarisedFlux(self.world.ground_map)

			mul_sca.SetStokesParameters()

			V, Vc, Vs, DoLP, AoLP = mul_sca.GetTotalContribution()
			self.MS_I0[self.time, self.ie_pc, self.ia_pc] = V
			self.MS_DoLP[self.time, self.ie_pc, self.ia_pc] = DoLP
			self.MS_AoLP[self.time, self.ie_pc, self.ia_pc] = AoLP

			print("MS TotalContribution", V, Vc, Vs, DoLP, AoLP)


			mul_sca.MakeOriginPlot(self.world.ground_map)
			mul_sca.Make3DPlot()
			mul_sca.MakeAltitudeHistogram()
			mul_sca.MakeScatteringHistograms()

		# plt.show()


	def PrintSystematicResults(self, header = True):
		# mag = -2.5 * np.log10(self.I_list[0, 0, 0] * 10**-9 / (3.0128 * 10**28))
		# print("Magnitude", mag)

		if self.save_name:
			save_file = open(self.save_name + ".txt", "w")

		if header:
			str_header = ""
			for key in self.in_dict.keys():
				if key[0] != "#":
					str_header += key + ","

			str_header += "datetime,I0,DoLP,AoRD,DI0,DDoLP,DAoLP,MSI0,MSDoLP,MSAoLP"
			print(str_header)
			if self.save_name:
				print(str_header, file=save_file)


		for t in range(self.world.sky_map.Nt):
			for ie_pc in range(self.world.Nb_e_pc):
				for ia_pc in range(self.world.Nb_a_pc):

					values = ""
					for k, v in self.in_dict.items():
						if k[0] != "#":
							if k == "azimuts":
								values += f"{self.world.a_pc_list[ia_pc] * RtoD}"
							elif k == "elevations":
								values += f"{self.world.e_pc_list[ie_pc] * RtoD}"
							else:
								values += v
							values += ","

					values += self.world.sky_map.times[t].strftime("%Y%m%d-%H%M%S") + ","
					values += f"{self.I_list[t, ie_pc, ia_pc]},"
					values += f"{self.DoLP_list[t, ie_pc, ia_pc]},"
					values += f"{self.AoRD_list[t, ie_pc, ia_pc] * RtoD},"
					values += f"{self.I_direct_list[t, ie_pc, ia_pc]},"
					values += f"{self.add_B_DoLP[t, ie_pc, ia_pc] * 100},"
					values += f"{self.add_B_AoLP[t, ie_pc, ia_pc] * RtoD}"
					values += f"{self.MS_I0[t, ie_pc, ia_pc]}"
					values += f"{self.MS_DoLP[t, ie_pc, ia_pc]}"
					values += f"{self.MS_AoLP[t, ie_pc, ia_pc]}"

					print(values[:])
					if self.save_name:
						print(values[:], file = save_file)

		# print("t,datetime,a_pc,e_pc,src_dist,dlos,min_alt,max_alt,wavelength,max_angle_discr,I0,DoLP,AoRD,Direct_I0:")
		#
		# for t in range(self.world.sky_map.Nt):
		# 	for ie_pc in range(self.world.Nb_e_pc):
		# 		for ia_pc in range(self.world.Nb_a_pc):
		# 			print(t, self.world.sky_map.times[t].strftime("%Y%m%d-%H%M%S"), self.world.a_pc_list[ia_pc]*RtoD, self.world.e_pc_list[ie_pc]*RtoD, self.world.ground_map.src_dist, self.world.atmosphere.d_los, self.world.atmosphere.h_r_min, self.world.atmosphere.h_r_max, self.world.wavelength, self.world.max_angle_discretization, self.I_list[t, ie_pc, ia_pc], self.DoLP_list[t, ie_pc, ia_pc], self.AoRD_list[t, ie_pc, ia_pc] * RtoD, self.I_direct_list[t, ie_pc, ia_pc], sep = "," )


	def GetLightParametersList(self):
		"""When the contributions of all maps for all times and all observations are done, compute the intensity, DoLP, AoLP lists.
		the lists have the shape: (time, e_pc, a_pc)"""
		for t in range(self.world.sky_map.Nt):
			for ie_pc in range(self.world.Nb_e_pc):
				for ia_pc in range(self.world.Nb_a_pc):
					self.GetLightParameters(ground = self.world.has_ground_emission, sky = self.world.has_sky_emission, time = t, ie_pc = ie_pc, ia_pc = ia_pc)

					self.I_list[t, ie_pc, ia_pc] = self.world.SetFluxUnit(self.I0)
					self.InonPola_list[t, ie_pc, ia_pc] = self.world.SetFluxUnit(self.InonPola)
					self.IPola_list[t, ie_pc, ia_pc] = self.world.SetFluxUnit(self.I0 - self.InonPola)
					# self.I_list[t, ie_pc, ia_pc] = self.I0 * self.world.luminosity_efficiency
					# self.InonPola_list[t, ie_pc, ia_pc] = self.InonPola * self.world.luminosity_efficiency
					# self.IPola_list[t, ie_pc, ia_pc] = (self.I0 - self.InonPola) * self.world.luminosity_efficiency
					self.DoLP_list[t, ie_pc, ia_pc] = self.DoLP
					self.AoRD_list[t, ie_pc, ia_pc] = self.AoRD

					# shutter_time = 60000 * 1000 #in milliseconds
					#
					# self.I_list_err[t, ie_pc, ia_pc] = np.sqrt(2 * self.I0 / shutter_time)
					# self.InonPola_list_err[t, ie_pc, ia_pc] = np.sqrt(2 * self.InonPola	/ shutter_time)
					# self.IPola_list_err[t, ie_pc, ia_pc] = np.sqrt(2 * self.IPola_list[t, ie_pc, ia_pc]	/ shutter_time)
					# self.DoLP_list_err[t, ie_pc, ia_pc] = np.sqrt(4 * (1 + ((self.DoLP/100.) ** 2 / 2)) / (self.I0 * shutter_time)) * 100
					# self.AoRD_list_err[t, ie_pc, ia_pc] = np.sqrt(1 / ((self.DoLP/100.)** 2 * self.I0 * shutter_time))


		# plt.plot(range(self.world.Nb_a_pc), self.InonPola_list[t, ie_pc, :])
		# plt.show()

	# def SetSaveName(self):
	# 	### Get the base name for the saving file -> unique for each parameter set
	# 	self.save_name = "results/" + self.world.ground_map.location + "_" + self.world.sky_map.mode + "_a" + str(np.round(self.a_pc*RtoD, 0)) + "_e" + str(np.round(self.e_pc*RtoD, 0)) + "_h" + str(np.round(self.world.sky_map.h, 0)) + "_" + str(np.round(self.world.atmosphere.h_r_min, 0)) + "-" + str(np.round(self.world.atmosphere.h_r_max, 0)) + "_Nlos" + str(np.round(self.world.atmosphere.Nlos))
	#
	# 	if self.world.ground_map.radius > 0:
	# 		self.save_name += "_D" + str(np.round(self.world.ground_map.radius * RT, 0))

	def GetIDAFromV(self, V, Vc, Vs):
		I0 = 2 * V
		if V != 0:
			DoLP = 100 * 2 * np.sqrt(Vc**2 + Vs**2) / V
		else:
			DoLP = 0
		AoLP = np.arctan2(Vs, Vc) / 2.

		return I0, DoLP, AoLP


	def GetLightParameters(self, ground = False, sky = False, time = None, ie_pc = None, ia_pc = None):
		"""Once the contributions of emission maps are computed, compute the I, DOLP, AoLP and the AoLP histogram."""

		if time == None: time = self.time
		if ie_pc == None: ie_pc = self.ie_pc
		if ia_pc == None: ia_pc = self.ia_pc

		if sky and ground: #If sky and ground exist
			self.I_direct_list[time, ie_pc, ia_pc] = self.world.GetDirect(time, ia_pc, ie_pc)
			self.I0  = np.sum(self.world.sky_map.total_scattering_map[time, ie_pc, ia_pc].flatten())
			self.I0 += np.sum(self.world.ground_map.total_scattering_map[time, ie_pc, ia_pc].flatten())
			self.I0 += self.I_direct_list[time, ie_pc, ia_pc]
			self.InonPola = self.I0 - np.sum(self.world.ground_map.scattering_map[time, ie_pc, ia_pc].flatten()) - np.sum(self.world.sky_map.scattering_map[time, ie_pc, ia_pc].flatten())
		elif ground: #If sky doesn't exist
			self.I0 = np.sum(self.world.ground_map.total_scattering_map[time, ie_pc, ia_pc].flatten())
			self.InonPola = self.I0 - np.sum(self.world.ground_map.scattering_map[time, ie_pc, ia_pc].flatten())
			# print(self.I0, self.InonPola)
		elif sky:  #If ground doesn't exist
			self.I_direct_list[time, ie_pc, ia_pc] = self.world.GetDirect(time, ia_pc, ie_pc)
			self.I0  = np.sum(self.world.sky_map.total_scattering_map[time, ie_pc, ia_pc].flatten())
			self.I0 += self.I_direct_list[time, ie_pc, ia_pc]
			self.InonPola = self.I0 - np.sum(self.world.sky_map.scattering_map[time, ie_pc, ia_pc].flatten())


		direct_V, direct_Vcos, direct_Vsin = self.world.GetVParamFromLightParam(self.I_direct_list[time, ie_pc, ia_pc], 0, 0)
		if self.add_B_pola and sky:
			direct_V, direct_Vcos, direct_Vsin = self.world.GetVParamFromLightParam(self.I_direct_list[time, ie_pc, ia_pc], self.add_B_DoLP[time, ie_pc, ia_pc], self.add_B_AoLP[time, ie_pc, ia_pc])

		self.V, self.Vcos, self.Vsin = self.cst_V, self.cst_Vcos, self.cst_Vsin

		if self.add_starlight:
			starlight = self.world.GetStarlight(time, ie_pc, ia_pc) / 2
			print("starlight", starlight)
			self.V += starlight

		if ground: #If ground and ground exist
			self.V 		+= self.world.ground_map.V_total[time, ie_pc, ia_pc]
			self.Vcos 	+= self.world.ground_map.Vcos_total[time, ie_pc, ia_pc]
			self.Vsin 	+= self.world.ground_map.Vsin_total[time, ie_pc, ia_pc]
		if sky: #If sky and ground exist
			self.V 		+= self.world.sky_map.V_total[time, ie_pc, ia_pc] + direct_V
			self.Vcos 	+= self.world.sky_map.Vcos_total[time, ie_pc, ia_pc] + direct_Vcos
			self.Vsin 	+= self.world.sky_map.Vsin_total[time, ie_pc, ia_pc] + direct_Vsin

		self.I0, self.DoLP, self.AoRD = self.GetIDAFromV(self.V, self.Vcos, self.Vsin)
		self.IPola = self.I0 * self.DoLP / 100.
		self.InonPola = self.I0 - self.IPola

		# hst, bins = np.zeros(self.N_bins), np.zeros(self.N_bins)
		if sky:
			A = self.world.sky_map.AoRD_map[time, ie_pc, ia_pc, :, :]
			Ip = self.world.sky_map.scattering_map[time, ie_pc, ia_pc, :, :]

			sky_hst, self.bins = self.world.MakeAoRDHist(A, Ipola = Ip)

		if ground:
			A = self.world.ground_map.AoRD_map[time, ie_pc, ia_pc, :, :]
			Ip = self.world.ground_map.scattering_map[time, ie_pc, ia_pc, :, :]
			ground_hst, self.bins = self.world.MakeAoRDHist(A, Ipola = Ip)

		if sky:
			self.hst = sky_hst
			if ground:
				self.hst += ground_hst
		elif ground:
			self.hst = ground_hst

		return self.I0, self.DoLP, self.AoRD

	def MakePlots(self):
		################################################################################
		###	PLOTTING
		################################################################################

		################################################################################
		###	Polar plots of intensity

		font = {'weight' : 'bold', 'size'   : 24}
		matplotlib.rc('font', **font)

		if mpi_rank == 0:
			if self.save_name:
				print(f"Making and saving plots in {self.save_name}")
			else:
				print("Making and NOT saving plots...")

		if self.world.has_ground_emission and not self.world.ground_map.is_point_source:
			self.world.ground_map.MakePlots(self.ie_pc, self.ia_pc, self.e_pc, self.a_pc, 10, save = self.save_name, metadata=self.in_dict)
			# self.MakeGroundMapPlots()
		if self.world.has_sky_emission:
			self.MakeSkyMapPlots(save = self.save_name, metadata=self.in_dict)

		# self.world.MakeAoRDHist(ground = self.is_ground_emission, sky = not self.is_ground_emission, double=True)
		# if self.world.is_single_observation:
		self.world.DrawAoLPHist(self.hst, self.bins, self.I0, self.DoLP, self.AoRD, double=True, save = self.save_name, show=False, metadata=self.in_dict)

		# self.MakeAoLPMap()

		# if self.world.is_single_observation and not self.world.is_time_dependant:
		# 	plt.show()

		# plt.close("all")


	def MakeSkyMapPlots(self, save = "", metadata = None):
		"""Make plots of the sky emission map contributions. Plot the initial intensity, the scattered intensity, the DOLP..."""
		f1, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex = True, sharey = True, figsize=(16, 8))
		ax1 = plt.subplot(221, projection='polar')
		ax2 = plt.subplot(222, projection='polar')
		ax3 = plt.subplot(223, projection='polar')
		ax4 = plt.subplot(224, projection='polar')

		# So that e = 90 is in the center of the plot if the emission layer is NOT the ground
		el = 90 - (self.world.sky_map.elevations[:]) * RtoD
		# i1 = ax1.pcolormesh(azimuts, el, scattering_map)
		i1 = ax1.pcolormesh(self.world.sky_map.azimuts[:], el, self.world.sky_map.cube[self.time, :, :])
		i2 = ax2.pcolormesh(self.world.sky_map.azimuts[:], el, self.world.sky_map.total_scattering_map[self.time, self.ie_pc, self.ia_pc, :, :])
		i3 = ax3.pcolormesh(self.world.sky_map.azimuts[:], el, self.world.sky_map.scattering_map[self.time, self.ie_pc, self.ia_pc, :, :])
		# i4 = ax4.pcolormesh(azimuts, el, total_scattering_map)
		i4 = ax4.pcolormesh(self.world.sky_map.azimuts[:], el, self.world.sky_map.AoRD_map[self.time, self.ie_pc, self.ia_pc, :, :]*RtoD, cmap=cm.twilight)

		cbar1 = f1.colorbar(i1, extend='both', spacing='proportional', shrink=0.9, ax=ax1)
		cbar1.set_label('Initial emission map')
		cbar2 = f1.colorbar(i2, extend='both', spacing='proportional', shrink=0.9, ax=ax2)
		cbar2.set_label('Total intensity map')
		cbar3 = f1.colorbar(i3, extend='both', spacing='proportional', shrink=0.9, ax=ax3)
		cbar3.set_label('Polarised intensity map')
		cbar4 = f1.colorbar(i4, extend='both', spacing='proportional', shrink=0.9, ax=ax4)
		cbar4.set_label('AoRD')

		# f1.suptitle("Relevant angles map at skibotn, with light pollution source at -45deg in azimut")

		for a in [ax1, ax2, ax3, ax4]:
			a.set_theta_zero_location("N")
			a.set_theta_direction(-1)
			a.set_thetamin(self.world.sky_map.I_zone_a_min * RtoD)
			a.set_thetamax(self.world.sky_map.I_zone_a_max * RtoD)

			a.set_rlim(0,90, 1)
			a.set_yticks(np.arange(0, 90, 20))
			a.set_yticklabels(a.get_yticks()[::-1])   # Change the labels

			a.add_artist(Ellipse((self.a_pc, 90 - self.e_pc * RtoD), width=self.world.ouv_pc, height=self.world.ouv_pc*RtoD, color="red"))

		if save:
			plt.savefig(save + '_skymaps.png', bbox_inches='tight', metadata = metadata)

	# def MakeGroundMapPlots(self):
	# 	"""Make plots of the ground emission map contributions. Plot the initial intensity, the scattered intensity, the DOLP..."""
	# 	f2, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex = True, sharey = True, figsize=(16, 8))
	# 	ax1 = plt.subplot(221, projection='polar')
	# 	ax2 = plt.subplot(222, projection='polar')
	# 	ax3 = plt.subplot(223, projection='polar')
	# 	ax4 = plt.subplot(224, projection='polar')
	#
	# 	i1 = ax1.pcolormesh(self.world.ground_map.azimuts, self.world.ground_map.distances * RT, self.world.ground_map.I_map)
	# 	i2 = ax2.pcolormesh(self.world.ground_map.azimuts, self.world.ground_map.distances * RT, self.world.ground_map.total_scattering_map[0, self.ie_pc, self.ia_pc, :, :])		# Intensity from (e, a) reaching us
	# 	i3 = ax3.pcolormesh(self.world.ground_map.azimuts, self.world.ground_map.distances * RT, self.world.ground_map.scattering_map[0, self.ie_pc, self.ia_pc, :, :])				# Polarized intensity from (e, a) reaching us
	# 	i4 = ax4.pcolormesh(self.world.ground_map.azimuts, self.world.ground_map.distances * RT, self.world.ground_map.DoLP_map[0, self.ie_pc, self.ia_pc, :, :])	# DoLP of scattered light from (e,a)
	# 	# i4 = ax4.pcolormesh((azimuts-A_lon) * RT, (distances-A_lat) * RT, AoRD_map * RtoD, cmap=cm.twilight)			# AoLP of scattered light from (e,a)
	#
	#
	# # Angle of polaisation of light from (e,a)
	# 	cbar1 = f2.colorbar(i1, extend='both', spacing='proportional', shrink=0.9, ax=ax1)
	# 	cbar1.set_label('Initial emission map')
	# 	cbar2 = f2.colorbar(i2, extend='both', spacing='proportional', shrink=0.9, ax=ax2)
	# 	cbar2.set_label('Total intensity map')
	# 	cbar3 = f2.colorbar(i3, extend='both', spacing='proportional', shrink=0.9, ax=ax3)
	# 	cbar3.set_label('Polarised intensity map')
	# 	cbar4 = f2.colorbar(i4, extend='both', spacing='proportional', shrink=0.9, ax=ax4)
	# 	cbar4.set_label('DoLP')
	#
	# 	# f2.suptitle("Relevant angles map at skibotn, with light pollution source at -45deg in azimut")
	#
	# 	for a in [ax1, ax2, ax3, ax4]:
	# 		a.set_theta_zero_location("N")
	# 		a.set_theta_direction(-1)
	#
	# 		a.add_artist(Arrow(0, 0, self.a_pc, self.world.atmosphere.h_r_max / np.tan(self.e_pc), color="red", width = 0.3))
	#
	# 	if self.save_individual_plots:
	# 		plt.savefig(self.path + self.save_name + '_groundmaps.png', bbox_inches='tight')

	# def MakeGroundMapPlots(self):
	# 	"""Make plots of the ground emission map contributions. Plot the initial intensity, the scattered intensity, the DOLP..."""
	# 	f2, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex = True, sharey = True, figsize=(16, 8))
	# 	ax1 = plt.subplot(221)
	# 	ax2 = plt.subplot(222)
	# 	ax3 = plt.subplot(223)
	# 	ax4 = plt.subplot(224)
	#
	# 	i1 = ax1.pcolormesh((self.world.ground_map.longitudes-self.world.ground_map.A_lon) * RT, (self.world.ground_map.latitudes-self.world.ground_map.A_lat) * RT, self.world.ground_map.I_map)
	# 	i2 = ax2.pcolormesh((self.world.ground_map.longitudes-self.world.ground_map.A_lon) * RT, (self.world.ground_map.latitudes-self.world.ground_map.A_lat) * RT, self.world.ground_map.total_scattering_map[0, self.ie_pc, self.ia_pc, :, :])		# Intensity from (e, a) reaching us
	# 	i3 = ax3.pcolormesh((self.world.ground_map.longitudes-self.world.ground_map.A_lon) * RT, (self.world.ground_map.latitudes-self.world.ground_map.A_lat) * RT, self.world.ground_map.scattering_map[0, self.ie_pc, self.ia_pc, :, :])				# Polarized intensity from (e, a) reaching us
	# 	i4 = ax4.pcolormesh((self.world.ground_map.longitudes-self.world.ground_map.A_lon) * RT, (self.world.ground_map.latitudes-self.world.ground_map.A_lat) * RT, self.world.ground_map.DoLP_map[0, self.ie_pc, self.ia_pc, :, :])	# DoLP of scattered light from (e,a)
	# 	# i4 = ax4.pcolormesh((longitudes-A_lon) * RT, (latitudes-A_lat) * RT, AoRD_map * RtoD, cmap=cm.twilight)			# AoLP of scattered light from (e,a)
	#
	#
	# # Angle of polaisation of light from (e,a)
	# 	cbar1 = f2.colorbar(i1, extend='both', spacing='proportional', shrink=0.9, ax=ax1)
	# 	cbar1.set_label('Initial emission map')
	# 	cbar2 = f2.colorbar(i2, extend='both', spacing='proportional', shrink=0.9, ax=ax2)
	# 	cbar2.set_label('Total intensity map')
	# 	cbar3 = f2.colorbar(i3, extend='both', spacing='proportional', shrink=0.9, ax=ax3)
	# 	cbar3.set_label('Polarised intensity map')
	# 	cbar4 = f2.colorbar(i4, extend='both', spacing='proportional', shrink=0.9, ax=ax4)
	# 	cbar4.set_label('DoLP')
	#
	# 	# f2.suptitle("Relevant angles map at skibotn, with light pollution source at -45deg in azimut")
	#
	# 	dy = np.cos(self.a_pc) * self.world.atmosphere.h_r_max / np.tan(self.e_pc)
	# 	dx = np.sin(self.a_pc) * self.world.atmosphere.h_r_max / np.tan(self.e_pc)
	# 	# for a in [ax1, ax2, ax3, ax4]:
	# 	# 	a.add_artist(Arrow(0, 0, dx, dy, color="red"))
	#
	# 	if self.save_individual_plots:
	# 		plt.savefig(self.path + self.save_name + '_groundmaps.png', bbox_inches='tight')
	#


	def MakeAoLPMap(self):
		"""Make a pyplot quiver of the AoLP contribution of each emission point. Vertical means AoLP=0, horyzontal AoLP = 90.
		Need a call to plt.show() after calling this function."""
		f4, (ax1) = plt.subplots(1, sharex=True, figsize=(16, 8))

		if not self.is_ground_emission:
			ax1.set_xlabel("Azimuth")
			ax1.set_ylabel("Elevations")

			X, Y = self.world.sky_map.mid_azimuts * RtoD, self.world.sky_map.mid_elevations * RtoD
			U, V = np.sin(self.world.sky_map.AoRD_map[self.time, self.ie_pc, self.ia_pc, :, :]) * self.world.sky_map.DoLP_map[self.time, self.ie_pc, self.ia_pc, :, :], np.cos(self.world.sky_map.AoRD_map[self.time, self.ie_pc, self.ia_pc, :, :]) * self.world.sky_map.DoLP_map[self.time, self.ie_pc, self.ia_pc, :, :]
			M = self.world.sky_map.DoLP_map[self.time, self.ie_pc, self.ia_pc, :, :]
		else:
			ax1.set_xlabel("Longitudes")
			ax1.set_ylabel("Latitudes")
			X, Y = self.world.ground_map.mid_longitudes * RtoD, self.world.ground_map.mid_latitudes * RtoD
			U, V = np.sin(self.world.ground_map.AoRD_map[self.ie_pc, self.ia_pc, :, :]) * self.world.ground_map.DoLP_map[self.ie_pc, self.ia_pc, :, :], np.cos(self.world.ground_map.AoRD_map[self.ie_pc, self.ia_pc, :, :]) * self.world.ground_map.DoLP_map[self.ie_pc, self.ia_pc, :, :]
			M = self.world.ground_map.DoLP_map[self.ie_pc, self.ia_pc, :, :]


		q = ax1.quiver(X, Y, U, V, M, pivot="middle", headwidth=0, headlength=0, headaxislength=0)
		ax1.quiverkey(q, 0.5, 1.1, 1, "AoRD")
		if not self.is_ground_emission:
			ax1.add_artist(Circle((self.a_pc * RtoD, self.e_pc * RtoD), radius=self.world.ouv_pc * RtoD, color="red"))

		if self.save_individual_plots:
			plt.savefig(self.save_name + '_AoRD_map.png', bbox_inches='tight', metadata = self.in_dict)

		################################################################################
		###	AoRD map as stream plot (continuous lines)
	# def MakeAoLPStram(self):
		# f5, (ax1) = plt.subplots(1, sharex=True, figsize=(16, 8))
		#
		# ax1.set_xlabel("Azimuts")
		# ax1.set_ylabel("Elevations")
		#
		# if h != 0: 	X, Y = azimuts * RtoD, elevations * RtoD
		# else:		X, Y = longitudes * RtoD, latitudes * RtoD
		# U, V = np.sin(AoRD_map) * DoLP_map, np.cos(AoRD_map) * DoLP_map
		# M = DoLP_map
		#
		# ax1.streamplot(X, Y, U, V, color=M, density=2, arrowstyle="-")
		# if h != 0: 	ax1.add_artist(Circle((a_pc * RtoD, e_pc * RtoD), radius=ouv_pc * RtoD, color="red"))

		# np.savetxt(self.path + self.save_name + "_scattering_map.cvs", self.scattering_map)
		# np.savetxt(self.path + self.save_name + "_DoLP_map.cvs", self.DoLP_map)
		# np.savetxt(self.path + self.save_name + "_total_scattering_map.cvs", self.total_scattering_map)
		# np.savetxt(self.path + self.save_name + "_AoRD_map.cvs", self.AoRD_map)

	def MakeMapPlots(self, f, axs):
		# extent = RtoD * np.array([self.world.a_pc_list[0], self.world.a_pc_list[-1], self.world.e_pc_list[0], self.world.e_pc_list[-1]])
		# a1 = axs[0].imshow(self.I_list[0,:,:], origin="lower", extent=extent)
		# axs[0].get_xaxis().set_visible(False)
		# cbar1 = f.colorbar(a1, extend='both', spacing='proportional', shrink=0.9, ax=axs[0])
		# a2 = axs[2].imshow(self.DoLP_list[0,:,:], origin="lower", extent=extent, cmap=plt.get_cmap("YlOrRd"))
		# axs[2].get_xaxis().set_visible(False)
		# cbar1 = f.colorbar(a2, extend='both', spacing='proportional', shrink=0.9, ax=axs[2])
		# a3 = axs[4].imshow(self.AoRD_list[0,:,:]*RtoD, origin="lower", extent=extent, cmap=plt.get_cmap("twilight"))
		# #axs[4].get_xaxis().set_visible(False)
		# cbar1 = f.colorbar(a3, extend='both', spacing='proportional', shrink=0.9, ax=axs[4])
		# a4 = axs[1].imshow(self.DoLP_list[0,:,:] * np.cos(2 * self.AoRD_list[0,:,:])/4, origin="lower", extent=extent, cmap=plt.get_cmap("bwr"))
		# axs[1].get_xaxis().set_visible(False)
		# cbar1 = f.colorbar(a4, extend='both', spacing='proportional', shrink=0.9, ax=axs[1])
		# a5 = axs[3].imshow(self.DoLP_list[0,:,:] * np.sin(2 * self.AoRD_list[0,:,:])/4, origin="lower", extent=extent, cmap=plt.get_cmap("bwr"))

		# cbar1 = f.colorbar(a5, extend='both', spacing='proportional', shrink=0.9, ax=axs[3])

		a1 = axs[0].pcolormesh(self.world.a_pc_list, 90 - (RtoD*self.world.e_pc_list), (self.I_list[0,:,:]))
		cbar1 = f.colorbar(a1, extend='both', spacing='proportional', shrink=0.9, ax=axs[0])

		a2 = axs[2].pcolormesh(self.world.a_pc_list, 90 - (RtoD*self.world.e_pc_list), (self.DoLP_list[0,:,:]))
		cbar1 = f.colorbar(a2, extend='both', spacing='proportional', shrink=0.9, ax=axs[2])

		a3 = axs[4].pcolormesh(self.world.a_pc_list, 90 - (RtoD*self.world.e_pc_list), (self.AoRD_list[0,:,:]))
		cbar1 = f.colorbar(a3, extend='both', spacing='proportional', shrink=0.9, ax=axs[4])

		a4 = axs[1].pcolormesh(self.world.a_pc_list, 90 - (RtoD*self.world.e_pc_list), (self.DoLP_list[0,:,:] * np.cos(2 * self.AoRD_list[0,:,:])/4))
		cbar1 = f.colorbar(a4, extend='both', spacing='proportional', shrink=0.9, ax=axs[1])

		a5 = axs[3].pcolormesh(self.world.a_pc_list, 90 - (RtoD*self.world.e_pc_list), (self.DoLP_list[0,:,:] * np.sin(2 * self.AoRD_list[0,:,:])/4))
		cbar1 = f.colorbar(a4, extend='both', spacing='proportional', shrink=0.9, ax=axs[3])

		for i in range(len(axs)):
			axs[i].get_yaxis().set_visible(False)
			axs[i].set_theta_zero_location("N")
			axs[i].set_theta_direction(-1)
			axs[i].set_ylim(0, None)


	def MakeSingleParamterPlots(self, f, axs):
		if self.world.Nb_e_pc == 1:
			xaxis = self.world.a_pc_list*RtoD
			I = self.I_list[0, 0, :]
			Ierr = self.I_list_err[0, 0, :] / 2.
			InonPola = self.InonPola_list[0, 0, :]
			InonPolaerr = self.InonPola_list_err[0, 0, :] / 2.
			IPola = self.IPola_list[0, 0, :]
			IPolaerr = self.IPola_list_err[0, 0, :] / 2.
			D = self.DoLP_list[0, 0, :]
			Derr = self.DoLP_list_err[0, 0, :] / 2.
			A = self.AoRD_list[0, 0, :] * RtoD
			Aerr = self.AoRD_list_err[0, 0, :] / 2. * RtoD

		elif self.world.Nb_a_pc == 1:
			xaxis = self.world.e_pc_list*RtoD
			I = self.I_list[0, :, 0]
			Ierr = self.I_list_err[0, :, 0] / 2.
			InonPola = self.InonPola_list[0, :, 0]
			InonPolaerr = self.InonPola_list_err[0, :, 0] / 2.
			IPola = self.IPola_list[0, :, 0]
			IPolaerr = self.IPola_list_err[0, :, 0] / 2.
			D = self.DoLP_list[0, :, 0]
			Derr = self.DoLP_list_err[0, :, 0] / 2.
			A = self.I_list[0, :, 0]
			Aerr = self.I_list_err[0, :, 0] / 2.

		axs[0].plot(xaxis, I, label = "Total")
		axs[0].fill_between(xaxis, I - Ierr, I + Ierr, alpha = 0.2)

		axs[0].plot(xaxis, InonPola, "r", label = "Non polarized")
		axs[0].fill_between(xaxis, InonPola - InonPolaerr, InonPola + InonPolaerr, alpha = 0.2)

		axs[0].plot(xaxis, IPola, "g", label = "Polarized")
		axs[0].fill_between(xaxis, IPola - IPolaerr, IPola + IPolaerr, alpha = 0.2)

		axs[2].plot(xaxis, D)
		axs[2].fill_between(xaxis, D - Derr, D + Derr, alpha = 0.2)

		axs[4].plot(xaxis, A)
		axs[4].fill_between(xaxis, A - Aerr, A + Aerr, alpha = 0.2)


	def MakeSummaryFigure(self, projection=None):
		f, axs = plt.subplots(nrows=3, ncols=2, sharex = True, figsize=(16, 8), subplot_kw={'projection': projection})
		axs = axs.flatten()
		axs[0] = plt.subplot(321)
		axs[1] = plt.subplot(322)
		axs[2] = plt.subplot(323)
		axs[3] = plt.subplot(324)
		axs[4] = plt.subplot(325)

		return f, axs

	def MakeSummaryPlot(self):

		# print("Azimut, Elevation, I, DoLP, AoLP")
		# for t in range(self.world.sky_map.Nt):
		# 	for ia, a in enumerate(self.a_pc_list):
		# 		for ie, e in enumerate(self.e_pc_list):
		# 			print(a*RtoD, e*RtoD, self.I_list[t][ie][ia], self.DoLP_list[t][ie][ia], self.AoRD_list[t][ie][ia]*RtoD)

		font = {'weight' : 'bold', 'size'   : 24}
		matplotlib.rc('font', **font)


		print("RECAP:")
		# print(f"Minimum Flux: {np.min(self.I_list[:,:,:])} nW")
		max_index = np.unravel_index(np.argmax(self.DoLP_list[:,:,:]), self.DoLP_list.shape)
		max_DoLP = self.DoLP_list[max_index]
		print(f"Maximum DoLP: {max_DoLP}% at time {max_index[0]} azimuth {self.world.a_pc_list[max_index[2]]*RtoD} and elevation {self.world.e_pc_list[max_index[1]]*RtoD}")

		if not self.world.is_time_dependant:
			if self.world.is_single_observation:
				# self.MakeAoLPHist(ground = self.world.has_ground_emission, sky = self.world.has_sky_emission)
				self.MakePlots()
				# pass
			else:
				if self.world.Nb_e_pc > 1 and self.world.Nb_a_pc > 1:
					f, axs = self.MakeSummaryFigure(projection='polar')
					self.MakeMapPlots(f, axs)
				else:
					f, axs = self.MakeSummaryFigure()
					self.MakeSingleParamterPlots(f, axs)

				axs[0].set_ylabel("Intensity ({})".format(self.world.flux_unit))
				# axs[0].legend()
				axs[2].set_ylabel("DoLP (\%)")
				axs[4].set_ylabel("AoLP (°)")
				axs[1].set_ylabel("Q (\%)")
				axs[3].set_ylabel("U (\%)")

				# f.subplots_adjust(hspace=0)

				if self.save_name:
					plt.savefig(self.save_name + "_summary.png", metadata=self.in_dict)
		else:
			f, axs = plt.subplots(3, sharex = True, figsize=(16, 8))
			axs[0] = plt.subplot(311)
			axs[1] = plt.subplot(312)
			axs[2] = plt.subplot(313)
			axs[0].set_ylabel("Intensity (nW)")
			axs[1].set_ylabel("DoLP (\%)")
			axs[2].set_ylabel("AoLP (°)")

			# axs[0].set_yscale('log')
			# axs[0].plot(range(self.world.sky_map.Nt), (self.I_list[:, 0, 0] / (self.I_direct_list[:, 0, 0] + self.I_list[:, 0, 0])).tolist())
			# print(self.I_list[:, 0, 0].tolist())
			# print(self.I_direct_list[:, 0, 0].tolist())
			axs[0].plot(range(self.world.sky_map.Nt), self.I_list[:, 0, 0].tolist())
			# axs[0].plot(range(self.world.sky_map.Nt), self.I_direct_list[:, 0, 0].tolist())
			axs[1].plot(range(self.world.sky_map.Nt), self.DoLP_list[:, 0, 0].tolist())
			axs[2].plot(range(self.world.sky_map.Nt), (self.AoRD_list[:, 0, 0]*RtoD).tolist())
			if self.save_name:
				plt.savefig(self.save_name + "_summary.png", metadata = self.in_dict)

			self.world.sky_map.MakeSkyCubePlot(self.world.a_pc_list, self.world.e_pc_list, self.world.ouv_pc, save = self.save_name, metadata=self.in_dict)

		# plt.savefig(self.path + self.save_name, bbox_inches='tight')
