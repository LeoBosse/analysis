#!/usr/bin/python3
# -*-coding:utf-8 -*

import mysql.connector as mysql
import sys
import os.path
import numpy as np

class Input:
	def __init__(self, template, input_file="", path = "./"):
		self.template = template
		self.path = path

		if not input_file:
			input_file = self.template

		self.file_name = input_file

		# print(self.path, self.file_name)

		self.dictionnary, self.comments = self.ReadInputFile(self.path + self.file_name)

	def Update(self, dict):
		for key, value in dict.items():
			self.dictionnary[key] = value

	def ReadInputFile(self, file_name):
		"""Read a given input file and return a dictionnary. First word = key, second = value. Remove empty lines, as many arguments as you want, in any order that you want."""
		with open(file_name, "r") as f: #open the file
			input = f.readlines() #store all lines in input

		input = [l.split()[:] for l in input if l != "\n"] #get first two words (separated by a space) for all non empty line

		dict = {i[0]: i[1] for i in input if len(i)>1} #Create dictionnary with all parameters names and values.
		comments = {i[0]: " ".join(i[2:]) for i in input if len(i) > 2} #Create a dictionnary of all comments of the input file.

		return dict, comments

	def WriteInputFile(self, folder="", file_name = "input", **kwargs):
		if folder == "":
			folder = self.path

		if kwargs:
			file_name = file_name.split(".")[0]
			# print("file_name", file_name)
			for k, v in kwargs.items():
				file_name += "_" + k + str(v)
			# print("file_name", file_name)
			# file_name += ".in"

		with open(folder + file_name + ".in", "w") as f_in:
			for key, value in self.dictionnary.items():
				f_in.write(" ".join((key, str(value))))
				if self.comments.__contains__(key):
					f_in.write("\t" + self.comments[key])
				f_in.write("\n")

		return file_name


class VendangeInput(Input):
	def __init__(self, input_file = ""):
		template = "template.in"
		path = "/home/bossel/These/Analysis/src/vendange/input_files/"
		Input.__init__(self, template, input_file, path = path)


	def Update(self, dict):
		for key, value in dict.items():
			if key == "observation_type":
				self.SetObservationType(value)
			else:
				self.dictionnary[key] = value

	def SetPollutionSource(self, lieu = "", az = 0, el = 0):
		if lieu.lower() == "skibotn":
			az, el = -45, 0
		elif lieu.lower() == "nyalesund":
			az, el = 1, 0
		elif lieu.lower() == "corbel":
			az, el = -60, 0

		self.dictionnary["pollution_source_azimut"] = az
		self.dictionnary["pollution_source_elevation"] = el


	def SetObservationType(self, new_type, el = 0):
		old_type = self.dictionnary["observation_type"]
		if new_type == old_type:
			return True

		self.dictionnary["observation_type"] = new_type

		if new_type == "fixed_elevation_continue_rotation":
			self.dictionnary["continue_rotation_elevation"] = self.dictionnary.pop("#continue_rotation_elevation")
			self.dictionnary["continue_rotation_times"] = self.dictionnary.pop("#continue_rotation_times")
			if el != 0:
				self.dictionnary["continue_rotation_elevation"] = el
			if old_type == "fixed_elevation_discrete_rotation":
				self.dictionnary["#discrete_rotation_azimut"] = self.dictionnary.pop("discrete_rotation_azimut")
				self.dictionnary["#discrete_rotation_times"] = self.dictionnary.pop("discrete_rotation_times")
				self.dictionnary["#discrete_rotation_elevations"] = self.dictionnary.pop("discrete_rotation_elevations")

		elif new_type == "fixed_elevation_discrete_rotation":
			self.dictionnary["discrete_rotation_azimut"] = self.dictionnary.pop("#discrete_rotation_azimut")
			self.dictionnary["discrete_rotation_times"] = self.dictionnary.pop("#discrete_rotation_times")
			self.dictionnary["discrete_rotation_elevations"] = self.dictionnary.pop("#discrete_rotation_elevations")

			if old_type == "fixed_elevation_continue_rotation":
				self.dictionnary["#continue_rotation_elevation"] = self.dictionnary.pop("continue_rotation_elevation")
				self.dictionnary["#continue_rotation_times"] = self.dictionnary.pop("continue_rotation_times")

		elif new_type == "fixed":
			if old_type == "fixed_elevation_continue_rotation":
				self.dictionnary["#continue_rotation_elevation"] = self.dictionnary.pop("continue_rotation_elevation")
				self.dictionnary["#continue_rotation_times"] = self.dictionnary.pop("continue_rotation_times")
			if old_type == "fixed_elevation_discrete_rotation":
				self.dictionnary["#discrete_rotation_azimut"] = self.dictionnary.pop("discrete_rotation_azimut")
				self.dictionnary["#discrete_rotation_times"] = self.dictionnary.pop("discrete_rotation_times")
				self.dictionnary["#discrete_rotation_elevations"] = self.dictionnary.pop("discrete_rotation_elevations")
