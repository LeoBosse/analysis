README file to read the data in the .txt save files.

Files are classic csv files with a space " " as column delimiter. One line gives the data for one rotation of the polarisating filter.
Four files are contained in the zip archive, corresponding to the four observations described in the article. Their titles are self explanatory.

Data files titles:
blue_rotations_without_magnet.txt
blue_rotations_with_magnet.txt
purple_rotations_without_magnet.txt
purple_rotations_with_magnet.txt

Here is a short description of each columns:

time  : time since the begining of the observation in milliseconds. Might not start at zero if the observation was cleaned.

V     : raw data of the V polarisation parameter (see https://doi.org/10.1051/swsc/2020036 for details)
Vcos	: raw data of the Vcos polarisation parameter (see https://doi.org/10.1051/swsc/2020036 for details)
Vsin	: raw data of the Vsin polarisation parameter (see https://doi.org/10.1051/swsc/2020036 for details)

I0	  : raw flux data, computed from V, Vcos, Vsin as I0 = 2V
DoLP	: raw DoLP data, computed from V, Vcos, Vsin as DoLP = 200 * sqrt(Vcos^2 + Vsin^2) / V
AoLP	: raw AoLP data, computed from V, Vcos, Vsin as AoLP = arctan(Vsin, Vcos) / 2

SV	  : smoothed data. Slidding averaging window of 10 seconds over the V raw data
SVcos	: smoothed data. Slidding averaging window of 10 seconds over the Vcos raw data
SVsin	: smoothed data. Slidding averaging window of 10 seconds over the Vsin raw data

SI0		: smoothed flux data computed from SV, SVcos and SVsin. Same formula than the raw data.
SDoLP	: smoothed DoLP data computed from SV, SVcos and SVsin. Same formula than the raw data.
SAoLP	: smoothed AoLP data computed from SV, SVcos and SVsin. Same formula than the raw data.

errI0	  : error associated with the raw data parameter I0
errDoLP : error associated with the raw data parameter DoLP
errAoLP	: error associated with the raw data parameter AoLP

errSI0	: error associated with the smoothed data parameter SI0
errSDoLP: error associated with the smoothed data parameter SDoLP
errSAoLP: error associated with the smoothed data parameter SAoLP
