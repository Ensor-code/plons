
<h1 align="center">
<img src="https://raw.githubusercontent.com/Ensor-code/plons/main/plons.png" width="300">
</h1>

![Build status](https://github.com/Ensor-code/plons/actions/workflows/build-and-test.yaml/badge.svg)
![Build status](https://github.com/Ensor-code/plons/actions/workflows/upload-to-pypi.yaml/badge.svg)
![Build status](https://github.com/Ensor-code/plons/actions/workflows/upload-to-anaconda.yaml/badge.svg)
![Build status](https://readthedocs.org/projects/plons/badge/?version=latest)


This is the README for the PLONS PLOtting tool for Nice Simulations. 
PLONS can be used to read in and analyse data of hydrodynamical simulations with the SPH code
[PHANTOM](https://phantomsph.bitbucket.io/) ([Price et al. 2018](https://ui.adsabs.harvard.edu/abs/2018PASA...35...31P/abstract)). 
This python package is currently tailored to single, binary and triple AGB wind models.
For a more detailed explanation of how to use PLONS, see ([readthedocs](https://plons.readthedocs.io/en/latest/)).

This tool is complementary to the SPH visualisation tool [SPLASH](https://users.monash.edu.au/~dprice/splash/) ([Price 2007](https://adsabs.harvard.edu/abs/2007PASA...24..159P)), but can be analysed more toroughly with PLONS.


PLONS pipeline
------------------

PLONS contains a pipeline that produces useful output files with information about your simulation, and currently contains the following functionalities:

1. 2D slice plots showing the density, temperature and velocity distribution of the last full dump of the model.​

2. 1D line plots (radial structure) of the global structure of the last dump of the model along the x-, y- and z-axes.​

3. Information and plots of the orbital evolution​

4. 1D spherical profiles, which are useful for single star models

To run this pipeline, use the script 'main.py' in plons/scripts/ .
To select a certain model and option, you can execute this script as python main.py -m ... -o ... (and filling in the model and option number instead of the ...).

Source file functions and example notebooks
----------
Next to the functionalities included in the pipeline, you can use Plons for more detailed analysises using the functions defined in the source files (in plons/src/plons/). 
Examples of how to use these functions can be found in some example notebooks (in docs/src/1_examples/), describing a more thorough analysis that makes use of line slices and sliceplots (/0_creating_figures_package), and an analysis of accretion disks around the companion (/2_accretion_disk).


Reading in PHANTOM data
-------

PHANTOM returns different types of useful data files.
Every certain timestep during the evolution of the model, dump files are output, typically named e.g. as 'wind_xxxxx'.
These files contain relevant data (such as position, velocity,
mass, density, energy, temperature) of the SPH particles in the model at a certain timestep.
The last 'x' rows in these datafiles correspond to the 'x' sink particles in the model.
Further, for every run, evolutionary '*.ev' are produced that contain information about the time evolution of usefull quantities such as the mass of a sink particle for example.


