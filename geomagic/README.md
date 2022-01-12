This is a README file for the geomagic module.

################################################################################
INSTALLATION
################################################################################

Download the sources files from https://gricad-gitlab.univ-grenoble-alpes.fr/bossel/analysis.

In a terminal, go to this folder and run the following command:
$pip install -e .
This will install geomagic on your computer, and take into account all changes you make on the code as you play with it.
You should then be able to just import a module like this in your python code:
import observation
or like any other module like numpy or matplotlib.

Geomagic is a requirement for POMEROL and Vendange.
It takes care of all complex geometry calculation (angles, distances, heights for a spherical Earth) in the observation.py file, where the observation class is defined.
