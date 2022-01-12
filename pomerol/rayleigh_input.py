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


class RayleighInput(Input):
	def __init__(self, input_file = ""):
		template = "RS_default.in"
		path = "/home/bossel/These/Analysis/src/rayleigh/input_files/"

		Input.__init__(self, template, input_file, path = path)


	def WriteInputFile(self, file_name, **kwargs):
		return Input.WriteInputFile(self, folder = self.path, file_name = file_name, **kwargs)
