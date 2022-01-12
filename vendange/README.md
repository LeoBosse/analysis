This is a README file for the Vendange data analysis code.

################################################################################
INSTALLATION
################################################################################

Download the sources files from https://gricad-gitlab.univ-grenoble-alpes.fr/bossel/analysis.

In a terminal, go to this folder and run the following command:
$pip install -e .
This will install vendange on your computer, and take into account all changes you make on the code as you play with it.

Do the same with the Geomagic package. It is a requirement for vendange.
Geomagic takes care of all complex geometry calculation (angles, distances, heights for a spherical Earth).


Open vendange_configuration.py and change all paths definition with your own computer setup.


################################################################################
Cru data management
################################################################################

Vendange was developped over three years and had to be adapted to the changes made to the instruments.
In particular, it does not read the hdf5 format yet.
Here is how you MUST store your Cru data in order to make it work.

This is absolutely not optimized, nor robust. It should all be rewriten and simplified to allow to read and use the hdf5 file format used in the database V3. A lot of code is made to accomodate SPP data, which are not in use anymore and should be completely forgotten!

Create a folder data/ in which you will put all the data you want.
Create a subfolder ptcu/, or any name you want where you will store all data of the Cru instruments.
In this folder, create as many folders as needed of the form: YYYYMMDD_Location (date_location, for example 20190307_Skibotn). Each of these folders will contain the data of this date and location.
For each observation, create subfolders of the form: cc_aAAA_eEE were "cc" is the initial of the colour of the filters on each channel of the instrument (for example botm if you use Grand Cru with blue-orange-mauve-turquoise). Coded are red 630, green 557.7, blue 427.8, mauve 391.4, turquoise 413, orange 620
aAA is the azimuth of the instrument. 0=North, 90=East, -90, 270=West, 180=South.
eEE is the elevation of the instrument. 0=horizontal, 90=vertical, zenith.
So if the instrument was pointing to at elevation 45deg facing west, you write a-90_e45 or a270_e45
You can also replace the aAA by "rot" if you performed almucantars. So it becomes rot_e45 for example.
You can put anything after that do differentiate two otherwise identical observations for example:
botm_a-90_e30_cloudy or botm_a-90_e30_5hz
Inside this last subfolder, you put the data files. They consist of N files of the form dataX.csv for X from 1 to N where N is the number of channels of the instrument (4 for Grand Cru, 2 for Corbel and Carmen).
There is also a config.csv file which recap the configuration of the instrument (IDConfiguration, calibration angles, dark, title, etc...)
Lastly, you create an input.in file from a template. The generic template can be found in input_files/template.in along with two other templates for (continuous or discrete) almucantars.
This file will tell the analysis software almost all it needs to know to treat the data (time jumps at end and beginig of observation, averaging window, value limit for bad data, etc...).

This process is tedious, but can be automatized using the recolte.py script.
This script goes to the SQL database, download the data you want and save them in the appropriate format based on the observation name you gave when lauching the measurments.
So if the observation with an IDConfiguration of 9999 is called 20190307_Skibotn_botm_a180_e90_blabla, you can execute in a terminal the line:
$python recolte 9999
And it will create the folders for you.
It should look like this:

data/
|		ptcu/
|		|		20190307_Skibotn/
|		|		|		botm_a180_e90_blabla/
|		|		|		|		input.in
|		|		|		|		config.csv
|		|		|		|		data1.csv
|		|		|		|		data2.csv
|		|		|		|		data3.csv
|		|		|		|		data4.csv

################################################################################
Other data management
################################################################################

If you want to compare Cru data with magnetometer data for example, you can put these data in your data folder, next to your folder ptcu (or watever you called it in the previous section Cru data management).

Create your magnetometer_data/ folder, and use the get_magnetometer_data.sh script to download the data in this folder. Then, vendange should be able to read them and plot them as you want.

You can do the same for all-sky camera images using the get_data camera.sh script. This will create a long and complicated arborescence of folders, but don't worry, vendange is able to read them as is!

Your data folder should then look something like this:

data/
|		ptcu/
|		magnetometer/
|		allsky/

Finally, go to the vendange_configuration.py script and change the paths in the __init__() function.

For EISCAT data, it is also possible to plot them against Cru data, but no automatised script is available. Just put your data somewhere (in the data folder maybe?) and change the path in the eiscat_data.py class. For this, you will have to play with the strings to comprehend your data format.



################################################################################
Running Vendange
################################################################################

If everything is set up right, running vendange is easy.
Just run the following from a terminal from the vendange folder:

$./main <path of the observation>

<path of the observation> is the path where the data you are interested in are stored.
For example, following the previous example: ptcu/20190307_Skibotn/botm_a180_e90_blabla
$./main ptcu/20190307_Skibotn/botm_a180_e90_blabla

This will treat the data contained in this folder using the inut.in file. You could give the input file explicitly if you want like that:
$./main ptcu/20190307_Skibotn/botm_a180_e90_blabla/your_input_file.in

Now, you have to setup a few things for it to work. Make sure your data are stored with a readable path as explained before.
Then, make sure the input file is correct (paths, and sensible argument values).
All arguments are available in the input_files_template.in file with extensive comments.

Then, if you want to play with the graphs, go to the mixer.py module, and play with the variables defined in the Mixer.SetGraphParameter() function.
This should be sufficient to obtain all kind of graphs.



################################################################################
Special arguments
################################################################################
You can run vendange with some arguments in the command line.
They are a little bit buggy, and should be used with care.

-f
Read the data from a .csv file that was created in a previous analysis of the data. It contains smoothed data and everything. This gains time for very long observations with a lot of data because you don't compute the windowed average.

-s [True, true, t, 1, False, false, f, 0]
Show the plots by default (if omitted). If put to False, don't show the graphs but produce and save them anyway.

-b [comparison bottle path]
Add a second observation to compare with the first one. Indicate the path of the data for this 2nd observation as you would do when running vendange normally. This will plot both bottle on the same graph.

-a  [append bottle path]
Add a second observation to appende to the first one. Indicate the path of the data for this 2nd observation as you would do when running vendange normally. This will just append the 2nd to the 1st and treat them as a single joined observation.


################################################################################
Trouble shooting
################################################################################

There is a lot of bugs in this programm! I don't use chemicals product to kill them in my vines, so that they reproduce very fast...
Seriously, as previously said, this code was started without a global development in mind, and had to accomodate for updates during 3 years. I learned a lot during this time and would do a lot of things differently today.
All this to say it is not a good code and that I apologize in advance for all the errors you'll enjoy while using it.
  
