#!/usr/bin/python3
# -*-coding:utf-8 -*



#######################################################################
# Contains a script to download automatically the CRU data from the instrument!
# Super helpful if the data was acquired following the corresponding protocol and naming conventions.
# Connects to the instrument, looks for the given observation ID, downloads it, make a sanity check on the data, and store it in the correct folder with a useable default input file.
# You can gain soooooo much time using this if you respect the correct naming convention when naming the observation campaign before hitting START on the instrument dashboard.
# Correct name:
# [YYYMMDD]_[Location]_[FF(FF)]_a[A]_e[E]_[comments]
# YYYMMDD: Date with year-month-day like 20241503 for March 15th, 2024
# Location: Location name. For example Skibotn. You can check the geomagic/geometry.py file for a list of registered locations. Unknown locations are also accepted, but some analysis optional calculations might not work.
# F: filter color by channel (for exemple: rovm for grand cru if filters are red-orange-green-purple)
# A: azimuth in degrees. can be positive or negative, any number of digits
# E: elevation in degrees. between 0 and 90, any number of digits
# comments: anything you want. often the polarizer rotation speed (e.g: 5hz)

# Author: Léo Bosse
# License: the freeest one, do whatever with my bad code if it can be helpfull, nobody really cares!
#######################################################################





import mysql.connector as mysql
import sys
import os.path
from pathlib import Path 
import numpy as np

import input

from vendange_configuration import *

path = global_configuration.data_path
# path = "/home/bossel/These/Analysis/data/"


arguments = sys.argv
nb_args = len(arguments)
print(arguments, nb_args)

instrument = arguments[1].lower()
IDConfiguration = [int(arguments[2])]
if nb_args == 4:
	IDConfiguration = range(int(arguments[2]), int(arguments[3]) + 1)
	print(f"4 arguments given. Will try to download all observations from ID {arguments[2]} to {arguments[3]} included.")
elif nb_args > 4:
	IDConfiguration = [int(arguments[i]) for i in range(2, nb_args)]
	print(f"4 arguments given. Will try to download all observations ID: {IDConfiguration}")




if instrument == "carmen":
	db = mysql.connect(
	  host = "152.77.134.83",
	  user = "root",
	  passwd = "Skib0tn2019"
	)
	db_name = "petitcru"
	folder = "ptcu"
	nb_voies = 1
	data_requete = lambda v, id: "SELECT IDProcessedData, Time, PMAvg, PMCosAvg, PMSinAvg, REFAvg, Comment, TempPM, TempOptical, TempAmbiant FROM processeddata{0} LEFT OUTER JOIN temperatures USING (Time, IDConfiguration) LEFT OUTER JOIN live_comments USING (Time, IDConfiguration) WHERE IDConfiguration = {1};".format(v, id)
elif instrument == "corbel":
	db = mysql.connect(
	  host = "152.77.134.82",
	  user = "root",
	  passwd = "Skib0tn2019"
	)
	db_name = "K2PI"
	nb_voies = 2
	data_requete = lambda v, id: "SELECT IDProcessedData, Time, PMAvg, PMCosAvg, PMSinAvg, Comment, Temp_PM, Temp_Optic, Temp_Ambiant FROM processeddata{0} LEFT OUTER JOIN temperatures USING (Time, IDConfiguration) LEFT OUTER JOIN live_comments USING (Time, IDConfiguration) WHERE IDConfiguration = {1};".format(v, id)
elif instrument == "corbel2":
	instrument = "corbel"
	db = mysql.connect(
	  host = "193.156.98.136",
	  user = "root",
	  passwd = "Skib0tn2019"
	)
	db_name = "K2PI"
	nb_voies = 2
	data_requete = lambda v, id: "SELECT IDProcessedData, Time, PMAvg, PMCosAvg, PMSinAvg, Comment, Temp_PM, Temp_Optic, Temp_Ambiant FROM processeddata{0} LEFT OUTER JOIN temperatures USING (Time, IDConfiguration) LEFT OUTER JOIN live_comments USING (Time, IDConfiguration) WHERE IDConfiguration = {1};".format(v, id)
elif instrument == "gdcu":
	db = mysql.connect(
	  host = "152.77.134.81",
	  user = "root",
	  passwd = "Skib0tn2019"
	)
	db_name = "K2PI"
	nb_voies = 4
	data_requete = lambda v, id: "SELECT IDProcessedData, Time, PMAvg, PMCosAvg, PMSinAvg, Comment, Temp_PM, Temp_Optic, Temp_Ambiant FROM processeddata{0} LEFT OUTER JOIN temperatures USING (Time, IDConfiguration) LEFT OUTER JOIN live_comments USING (Time, IDConfiguration) WHERE IDConfiguration = {1};".format(v, id)


cursor = db.cursor()

cursor.execute("USE {}".format(db_name))

# num_fields = len(cursor.description)
# field_names = [i[0] for i in cursor.description]
# requete = "select * from configuration;"
# cursor.execute(requete)
# config = cursor.fetchall()
# print(config)

for id in IDConfiguration:
	###Get Configuration table
	requete = "SELECT * FROM configuration WHERE IDConfiguration = {};".format(id)
	cursor.execute(requete)
	config = cursor.fetchone()
	columns = [i[0] for i in cursor.description]

	print("IDConfiguration {0}: {1}".format(id, config[columns.index("CM_Comments")]))


	### Test if at least one PM was recording
	valid_data = False
	all_data = []
	for v in range(1, nb_voies + 1):
		requete = data_requete(v, id)

		cursor.execute(requete)
		data_columns = [i[0] for i in cursor.description]
		# print(data_columns)
		comment_index = data_columns.index("Comment")
		# print(comment_index)

		all_data.append(cursor.fetchall())


		#Replace all "," in comments to avoid messing up the .csv file.
		for l in range(len(all_data[-1])):
			comment = all_data[-1][l][comment_index]
			if comment:
				print(all_data[-1][l][comment_index])
				all_data[-1][l][comment_index] = str(comment).replace(",", ";")

		###Check if the data is empty
		if len(all_data[-1]) > 0:
			valid_data = True

	# If all canals are empty, the data is invalid and not saved
	if not valid_data:
		print("IDConfiguration {0}: No valid data on any channels.".format(id))
		continue

	### Get saving folder names
	title = config[columns.index("CM_Comments")].split("_")
	# print(title)
	date = config[columns.index("Timestamp")]
	lieu = title[1].capitalize()
	i = 0
	try:
		if title[2].lower() in ["sud", "north"]:
			i = 1
			lieu += title[2].capitalize()

		speed = config[columns.index("CM_PolarizerSpeed")]
		try:
			filters = title[2+i]
			az = title[3+i]
			observation_type = "fixed"
			if az == "rot":
				observation_type = "fixed_elevation_continue_rotation"

			el = title[4+i]

			if len(title) > 5+i:
				com = title[5+i:]
		except:
			if len(title) > 2+i:
				observation_type = "fixed"
				el = 45
				com = ""
				print("IDConfiguration {0}: Wrong Title. But I know where to save it!".format(id))
				print("Wrong title is:" + config[columns.index("CM_Comments")])
			else:
				raise("ValueError")
	except:
		if config[columns.index("CM_Comments")] == "Scheduled Campaign":
			# For automatic observation when Corbel Cru was in Corbel (2021-22)
			observation_type = "fixed"
			el = 45
			az = 35
			filters = "vm"
			lieu = "Corbel"
			date = config[columns.index("Timestamp")]
			i = 0
			speed = 5
			com = ""
		else:
			print("IDConfiguration {0}: Wrong Title. I can't save it...".format(id))
			print("Wrong title is:" + config[columns.index("CM_Comments")])
			# continue


	# ### For saving nice titles.
	folder = "ptcu" + date.strftime("%Y%m%d") + "_" + lieu
	# subfolder = "_".join(title[2+i:])
	subfolder = f"{filters}_{az}_{el}_{int(speed)}hz"

	### For saving weird titles. Just save it as is, with the ID configuration
	# folder = "ptcu/"
	# subfolder = str(config[columns.index("IDConfiguration")]) + "_" + str(config[columns.index("CM_Comments")])

	###Create saving data folder
	os.makedirs(path / folder, exist_ok=True)
	try:
		os.makedirs(path / folder / subfolder)
	except:
		nb_doublons = len([d for d in os.listdir(path / folder) if d[:len(subfolder)] == subfolder])
		subfolder += "_" + str(nb_doublons)
		os.makedirs(path / folder / subfolder)


	##Make input file:
	input_file = input.VendangeInput()
	# input_file.SetPollutionSource(lieu = lieu)
	input_file.SetObservationType(observation_type, el = el)
	input_file.Update({	"instrument_name": instrument,
					"data_files": folder + subfolder,
					"saving_name": "_".join((folder.replace("/", "_"), subfolder)),
					"observation_type": observation_type,
					"rotation_per_sec": speed
				})

	input_file.WriteInputFile(folder = path / folder / subfolder)


	###Write config.csv
	np.savetxt(path / folder / subfolder / "config.csv", [columns, config], delimiter = ",", fmt='%s')

	# print(f"saving config in {path + folder + subfolder}...")

	###Get Data
	for v in range(1, nb_voies + 1):
		# requete = data_requete(v, id)
		#
		# cursor.execute(requete)
		# data = cursor.fetchall()
		# columns = [i[0] for i in cursor.description]


		###Write dataV.csv
		np.savetxt(path / folder / subfolder / "data{}.csv".format(v), all_data[v - 1], delimiter = ",", fmt='%s', header=",".join(data_columns), comments="")
