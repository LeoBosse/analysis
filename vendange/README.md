This is a README file for the Vendange data analysis code.

***
# INSTALLATION
***

Download the sources files from https://gricad-gitlab.univ-grenoble-alpes.fr/bossel/analysis.

In a terminal, go to this folder containing this README file and run the following command:   

        $pip install -e .

This will install Vendange on your computer, and take into account all changes you make on the code as you play with it.
If you play with different python environment, make sure to install it in the right one.

Do the same with the Geomagic package. It is a requirement for vendange.
Run:  

        $cd ../geomagic
        $pip install -e .

Geomagic takes care of all complex geometry calculation (angles, distances, heights for a spherical Earth).


Open `vendange_configuration.py` and change all paths definition with your own computer setup. 


***
# Cru data management
***

Vendange was developped over three years and had to be adapted to the changes made to the instruments.
In particular, it does not read the hdf5 format yet.
Here is how you MUST store your Cru data in order to make it work.

This is absolutely not optimized, nor robust. It should all be rewriten and simplified to allow to read and use the hdf5 file format used in the database V3. A lot of code is made to accomodate SPP data, which are not in use anymore and should be completely forgotten!

1. Create a folder data/ in which you will put all the data you want.

2. Create a subfolder ptcu/, or any name you want where you will store all data of the Cru instruments.

3. In this folder, create as many folders as needed of the form: `YYYYMMDD_Location` (date_location, for example `20190307_Skibotn`). Each of these folders will contain the data acquired on this date and location.

4. For each observation, create subfolders of the form: `cc_aAAA_eEE` were `cc` is the initials of the colour of the filters on each channel of the instrument (for example `botm` if you use Grand Cru with blue-orange-mauve-turquoise). hard coded in vendange are red 630, green 557.7, blue 427.8, mauve 391.4, turquoise 413, orange 620  
`aAA` is the azimuth of the instrument. 0=North, 90=East, -90=270=West, 180=South.  
`eEE` is the elevation of the instrument. 0=horizontal, 90=vertical, zenith.  
So if the instrument was pointing at elevation 45deg facing west, you write `a-90_e45` or `a270_e45`
You can also replace the `aAA` by `rot` if you performed almucantars. So it becomes `rot_e45` for example.
You can put anything after that do differentiate two otherwise identical observations for example: `botm_a-90_e30_cloudy` or `botm_a-90_e30_5hz`
5. Inside this last subfolder, you put the data files. They consist of N files of the form `dataX.csv` for X from 1 to N where N is the number of channels of the instrument (4 for Grand Cru, 2 for Corbel and Carmen).  
There is also a `config.csv` file which recap the configuration of the instrument (IDConfiguration, calibration angles, dark, title, etc...)
1. Lastly, you create an `input.in` file from a template. The generic template can be found in `input_files/template.in` along with two other templates for (continuous or discrete) almucantars. The default name of the input file is `input.in`, however, you can rename it or create a copy to test different analysis parameters. See the Running Vendange for more info on how to run none default input files.  
This file will tell the vendange analysis software almost all it needs to know to treat the data (time jumps at end and beginig of observation, averaging window, limit values for bad data, etc...).

This process is tedious, but can be automatized using the `recolte.py` script.
This script goes to the SQL database, download the data you want and save them in the appropriate format based on the observation name you gave when lauching the measurments.
So if the observation with an IDConfiguration of 9999 is called `20190307_Skibotn_botm_a180_e90_blabla`, you can execute in a terminal the line:  

        $python recolte.py 9999

And it will create the folders for you.
It should look like this:  

+ data/  
  + ptcu/  
    + 20190307_Skibotn/  
      + botm_a180_e90_blabla/  
        + input.in  
        + config.csv  
        + data1.csv  
        + data2.csv  
        + data3.csv  
        + data4.csv  



***
# Other data management
***

If you want to compare Cru data with magnetometer data for example, you can put these data in your data folder, next to your ptcu folder (or watever you called it in the previous section Cru data management).

Create your `magnetometer_data/` folder, and use the `get_magnetometer_data.sh` script to download the data in this folder. Then, vendange should be able to read them and plot them as you want.

You can do the same for all-sky camera images using the `get_data_camera.sh` script. This will create a long and complicated arborescence of folders, but don't worry, vendange is able to read them as is!

Your data folder should then look something like this:  

+ data/  
    + ptcu/  
    + magnetometer/  
    + allsky/  

Finally, go to the `vendange_configuration.py` script and change the paths in the `__init__()` function to tell vendange where these new data sources are located.

For EISCAT data, it is also possible to plot them against Cru data, but no automatised script is available. Just put your data somewhere (in the data folder maybe?) and change the path in the `eiscat_data.py` class. For this, you will have to play with the strings to comprehend your data format. From my experience, I had to rewrite the `eiscat_data.py` script each year to accomodate new or changing eiscat configurations and data format. So be prepared to read and write some python code!



***
# Running Vendange
***

If everything is set up right, running vendange is easy.
Just run the following from a terminal from the vendange folder:  

        $./main <path of the observation>

`<path of the observation>` is the path where the data you are interested in are stored.
For example, following the previous example: `ptcu/20190307_Skibotn/botm_a180_e90_blabla` 

        $./main ptcu/20190307_Skibotn/botm_a180_e90_blabla

This will treat the data contained in this folder using the `input.in` file. You could give the input file explicitly like that:  

        $./main ptcu/20190307_Skibotn/botm_a180_e90_blabla/your_input_file.in

This is usefull to test different analysis procedure, zoom in on interesting events and save the results to a different save path.
Now, you have to setup a few things for it to work. Make sure your data are stored with a readable path as explained before.
Then, make sure the input file is correct (paths, and sensible argument values).
All arguments are available in the `input_files_template.in` file with extensive comments.

Then, if you want to play with the graphs, go to the `mixer.py` module, and play with the variables defined in the `Mixer.SetGraphParameter()` function. I'm sorry, the plotting parameters are not defined outside of the python script, but hard coded in vendange... I tried to put all of them on the same place though.
This should be sufficient to obtain all kind of graphs.



***
# Special arguments
***

You can run vendange with some arguments in the command line.
They are a little bit buggy, and should be used with care.

`-f`  
Read the data from a .csv file that was created in a previous analysis of the data. It contains smoothed data and everything. This gains time for very long observations with a lot of data because you don't compute the windowed average.

`-s [True, true, t, 1, False, false, f, 0]`  
Show the plots by default (if omitted). If put to False, don't show the graphs but produce and save them anyway.

`-b [comparison bottle path]`  
Add a second observation to compare with the first one. Indicate the path of the data for this 2nd observation as you would do when running vendange normally. This will plot both bottle on the same graph.

`-a  [append bottle path]`  
Add a second observation to appende to the first one. Indicate the path of the data for this 2nd observation as you would do when running vendange normally. This will just append the 2nd to the 1st and treat them as a single joined observation.


***
# Trouble shooting
***

There is a lot of bugs in this programm! I don't use chemicals product to kill them in my vines, so that they reproduce very fast...
As previously said, this code was started without a global development in mind, and had to accomodate for instrumental updates, new models and needs during -3- 5 years. I learned a lot during this time and would do a lot of things differently today.
All this to say it is not a good code and that I apologize in advance for all the errors you'll enjoy while using it.
  
To help you on your quest, I will try to present shortly what is happening in vendange at the crucial steps. It may help you find understand the origin of a error you encounter.

1. `main.py`: Everything starts here. This scripts read the command line arguments, and go to the folder of the specified input file to find some infos. It counts the number of dataX.csv to infer the number of channels in the instrument. From this, it create a `Bottle` object (see `bottle.py`) for each one of them. If you want only a given channel, you can find the corresponding for loop and change the list of channels. Warning: the channels are numbered starting at 1, there is no channel 0. It will append additional bottles to the first one, or create 'comparison bottles' (comp_bottle), to plot in parallel for better comparisons. Then, creates and run the Mixer object, where all data types are mixed and plotted.

2. `bottle.py`: This is the main object containing the data and sorting, cleaning and computing all interesting values. A `Bottle` object is the parent of two `PTCUBottle` and `SPPBottle` objects corresponding respectively to the CRU instrument serie and the SPP instrument. SPP is long gone and in a museum. The `PTCUBottle` is what you will use. Here is the creation process of a bottle:
   1. Read and load the `input.in` file specified on the command line. Store all parameters in a dictionnary.
   2. Read and load the data from the `dataX.csv` file corresponding to the bottle channel. Also loads the `config.csv` file containing useful info on the instrument parameters of the observation.
   3. All data is now stored in a panda DataFrame (since sept 2023) for faster and easier code. Some computations might still be using old list methods which are not optimized. Feel free to hunt them down and correct.
   4. The times of each data point is converted to a `datetime.TimeDelta` object, and the start of the observation as a `datetime.DateTime` object.
   5. Clean the data by removing the time jumps at the beggining and the end of the data. Removes all points outside the intensity and DoLP bounds set in the `input.in` file.
   6. Creates columns in the `bottle.data` DataFrame for the Intensity, DoLP, AoLP (IDA) from the V, Vcos, Vsin measurements of the instruments. Apply a rolling window to compute the smoothed V and IDA parameters, uncertainties, derivatives...
   7. From the derivative of the intensity measurments, compute the polarisation due only to the change of intensity and removes it from the final polarisation. This prevents phantom pola due to very rapid changes in polarisation (pulsating aurora, car headlamps...)
   8. When needed, also computes the geometry of the observation with the magnetic field orientation, simple light pollution model of the AoLP, orientation of the instrument as a function of time for almucantars.
   9. It also saves the data uner a `.txt` and a `.hdf5` file format. It should also be able to load the data from any save file.

3.  `mixer.py`: Contains the `Mixer` object where other datasets are loaded, and the plots are drawn and saved.
    1. The `Mixer.__init__` method is a mess right now. Everything is happening inside it. Let's go through it anyway. The idea is to loop through each `Bottle` object to adapt the multiinstrumental comparisons and plotting. So it simply loops over every bottle created earlier, and do the following steps.
    2. Set the graph parameters in `Mixer.SetGraphParameter(bottle)`. Go here to change and customize by hand a lot of your graphs.
    3. Tries to load data from other instruments (all sky cameras, magnetometer, EISCAT, equivalent currents, Rayleigh scattering (RS) models from POMEROL...). Each instruments has its own object to read, load, clean and validate its data. issing data can be due to:
       1. none existence, have you downloaded it properly, and put in the correct folder specified in `vendange_configuration.py`?
       2. The time range of available data does not match the time range of the bottle. Check if you have data for 2 days in case it stops at midnight, check the UT or LT time format and time zones.
       3. Did you ask to show the data? go to `Mixer.SetGraphParameter(bottle)` and look for a parameter called `show_eiscat`, `show_magnetometer` or `show_currents`...
    4. Then plots the graphs with all the data asked for and available. `Mixer.MakePlots()` is the method plotting the main graph with intensity, DoLP and AoLP triple time series from top to bottom. Each value is plotted on a subplot in its own function. However, to create more figures with other data, set `Mixer.make_optional_plots`, `Mixer.show_correlations` or `Mixer.show_SN` to `True` in `Mixer.SetGraphParameter()`. All graphs are automatically saved in the data folder using the appendix defined in the `input.in` folder. Feel free to add your own plotting method for your needs. You can simply call it in the `Mixer.__init__()` method, next to already existing calls to plotting function.
