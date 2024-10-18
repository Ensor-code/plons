Pipeline
########

PLONS contains a pipeline that produces useful output files with information about your simulation

To use the pipeline, the pipeline needs to be :ref:`installed <link-installation-pipeline>` first.

When running the pipeline, the user settings will be integrated, and models will be searched in the ``data_location``. All models found will be displayed, and the model for which the pipeline needs to run can be selected by typing the correct number. Multiple models can be selected by adding them with a space bar inbetween. When typing ``a`` or ``all``, all models will be used.

Later a selection of the pipeline can be chosen to be executed. It currently contains the following options:

1. 2D slice plots showing the density, temperature and velocity distribution of the last full dump of the model

2. 1D line plots (radial structure) of the global structure of the last dump of the model along the x-, y- and z-axes

3. Information and plots of the orbital evolution​

4. 1D spherical profiles, which are useful for single star models

To run this pipeline, use the script 'main.py' in plons/scripts/ .
To select a certain model and option, you can execute this script as python main.py -m ... -o ... (and filling in the model and option number instead of the ...).


Pipeline - More detailed information
------------------------------------

This pipeline loads the wind.in and wind.setup files first. This is the general setup information
of the model (configuration of the system, information of the thermodynamics, evolution time, wind
outflow setup,...). Next the last full dump is loaded, and using this data, some other useful
quantities are calculated (gas pressure/temperature, speed, sperical coordinates, sound speed of
the gas,...).  Lastly, the data of the sink particle(s) is loaded. From this file, also the period
and orbital velocity of the sink particles is calculated, if the system is a multiple, as a function of
time. The last entry in the sink data corresponds to the full dump loaded. Therefore the
position/velocity/mass/... of this last entry are also saved in the full dump data for easy usage.


(1) 2D slice plots 
___________________

This option creates 2D sliceplots in the orbital plane (face-on view) and a perpendicular meridional plane (edge-on view) 
at different scales. It shows the density, temperature and velocity distribution of your model in these 2 planes. 
To get the data in these infinitely thin slices, the smoothing kernel from PHANTOM is used, and the calculations are done on a grid of n*n
pixels, using the data of the m nearest neighbours (SPH particles) of those grid points. 
These plots can be useful to get a first visual representation of the morphology in the AGB wind of the simulation, 
and can be optimised according to your personal needs. 


(2) 1D radial structure (line slice) plots 
___________________________________________

This option creates 1D radial line profiles of the density, velocity and temperature distribution along the x-, y- and z-axes.
The data along these lines is achieved in the same way as explained in (1), but for a 1D grid. 
These plots can be usefull as a more detailed and quantitive representation of the morphology, for example to identify global assymetries and non-spherical density distributions.
For an example use case, see Appendix A of [Malfait et al. 2021](https://ui.adsabs.harvard.edu/abs/2021A%26A...652A..51M/abstract)


(3) Information of the orbital evolution.
_________________________________________

This option produces several plots with usefull information about the orbital properties and orbital evolution of the simulation,
such as the mass and angular momentum accretion by the companion.
In order to get this information, the evolutionary sink files are used (.ev), which give relevant quantities of the sink particles as a function of time.


(4) 1D spherical profiles for single star models
________________________________________________

This option should be used for single star simulations, and compares your 3D data to a 1D solution. It reads the wind_1D.dat data to plot the 1D solution.
