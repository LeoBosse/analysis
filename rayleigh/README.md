This is a README file for the POMEROL simulation code.

################################################################################
INSTALLATION
################################################################################

Download the sources files from https://gricad-gitlab.univ-grenoble-alpes.fr/bossel/analysis.
Place them all in a directory, and change the "src_path" in the RS_default.in input file (and all input files you'll use).

In a terminal, go to this folder and run the following command:
$pip install -e .
This will install POMEROL on your computer, and take into account all changes you make on the code as you play with it.

Do the same with the Geomagic package. It is a requirement for POMEROL.
It takes care of all complex geometry calculation (angles, distances, heights for a spherical Earth).

Then downlaod the data you need. (ground emission maps, all sky camera pictures, elevation map).
They are not included in the package for size and legal reasons.



################################################################################
RUNING YOUR FIRST MODEL
################################################################################

We will present a few steps to set up your first POMEROL model.
First, you need an input file. Tis is the file where you will define all the parameters of your model.
You can find a few example files in the folder input_files.

You can copy one of the templates, change its name and change the parameters.
Always put the extension ".in".
You can also organize your input files in subfolders inside the input_files folder. In this case, you must specify the your subfolders when calling POMEROL.

A description of the input parameters found in the input files can be found in the default and example input files. Each parameter is commented.


Once you have your input file ready, just run POMEROL in a terminal like this:
./main <input file name>
or
mpirun -n <Number of cores> <input file name>


In all cases, only specify the name, without the extension or the path. For example, if your input_files folder look like this:
input_file:
|		example_file_1.in
|		example_file_2.in
|		example_folder:
|		|		test1.in
|		|		test2.in

you can call:
./main example_file_2
or:
./main example_folder/test1


To run on multiple cores, POMEROL uses the mpi4py module (https://mpi4py.readthedocs.io/en/stable/).
To use this functionnality, launch the code with the following command:

mpirun -n <Number of cores> <input file name>


################################################################################
RUNING SYSTEMATIC PARAMETER STUDY
################################################################################

A tool is created to help you run automatically POMEROL. This is helpfull if you want to test the influence of one or two parameters one your results.
The code is contained on the systematic.py script
It is, for now, highly not user friendly... Sorry!

It works like this.

In the funciton RunSystematicSimulation() the dictionnary in_dict lets you define all the parameters you want to cycle through.
All possible combination of parameters will be launched automatically one after the other.
All the parameters from the RS_default.in input file can be specified, and the one you do not explicitely define here are taken from this default file.

The output_result_file will contain the results of all the runs you performed with this script. Each line contains as many paramters value as there is in RS_default.in file, plus the flux, DoLP, AoLP for total and direct light.

Once this is set, you can launch your automatic runs with:

./systematic

Once all your runs are done, a few basic tools are coded to plot the results.
run the following in a terminal:

./systematic -p <result_log_file>

the <result_log_file> argument is the file created when you launched your runs earlier and defined with the output_result_file variable.
To choose which parameter you wwant to plot against which, go near the end of systematic.py, in the if __name__ == "__main__": condition.
There, edit the PlotData() function as follows.
The second argument is the x-axis parameter, usually "azimuts", but can be anything existing in the RS_default.in file.
The third argument is the parameter(s) on the y-axis.  Again, can be anything from RS_default.in, but can also be "all". The "all" option will produce 3 plots with the Flux, DoLP and AoLP. Multiple plots can be achieved by specifying a list of parameters as the third argument. Finally, you can also specify multiple parameters separated by ";" so that they appear on the same graph.
All arguments after these are optional but highly recommanded if you want the script to work! You should specify all argument values to select only data where these arguments were set to a specific value. And you should set the other variable arguments to "*". Never put more than 2 free parameters, the code cannot handle it!

As an example, to study the effect of the first 10 kilometer of the atmosphere, you can set it like this:

output_result_file = "log/systematic_atm_max_height"
in_dict = {	"azimuts": ["all"],
						"elevations": [45],
						wavelength: [557.7, 630.0],
						"RS_max_altitude": np.arange(1, 10, 1),
						"ground_mode": ["point_I100_a180_d5"]
					}

run ./systematic. This will create a file with all model results for all azimuts and all atospheric max altitudes.
It also prints a long string for each POMEROL run it launches. It allows to follow the progress and estimates the time it needs.
It stores all the logs of each runs in a file build like this: log/<long printed string>.out so that you can check if there is an error somewhere.

Then, to plot the results, set PlotData() like this:

PlotData(data, "azimuts", "all", RS_max_altitude="*", wavelength=557.7)

This will produce a plot of almucantars at 45° of elevation, for a point source on ground 5km south of the instrument.

!!!! Be carefull with the name of the reult output file !!! It is not erased when launching a new batch of runs, so you might mix results that are very different !!!!

For more control on the plots or the way this script launches the model, you will have to go down the rabbit hole and code it yourself.
