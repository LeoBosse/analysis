#!/usr/bin/python3
# -*-coding:utf-8 -*

import mysql.connector as mysql
import sys
import os.path
import numpy as np
from input import Input

class RayleighInput(Input):
	def __init__(self, input_file = ""):
		template = "RS_default.in"
		path = "/home/bossel/These/Analysis/src/rayleigh/input_files/"

		Input.__init__(self, template, input_file, path = path)


	def WriteInputFile(self, file_name, **kwargs):
		return Input.WriteInputFile(self, folder = self.path, file_name = file_name, **kwargs)
