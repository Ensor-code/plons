


-------------------------------
    README PLONS
-------------------------------
![Build status](https://github.com/Ensor-code/plons/actions/workflows/build-and-test.yaml/badge.svg)
![Build status](https://github.com/Ensor-code/plons/actions/workflows/upload-to-pypi.yaml/badge.svg)
![Build status](https://github.com/Ensor-code/plons/actions/workflows/upload-to-anaconda.yaml/badge.svg)


This is the README for the PLONS PLOtting tool for Nice Simulations. It can be used to
read in and analyse data of hydrodynamical simulations with the SPH code PHANTOM 
(Price, D. J., Wurster, J., Tricco, T. S., et al. 2018, PASA, 35, e031). PHANTOM returns 
(different types of useful data files: 'wind_xxxxx'-files and '.ev'-files. The former can
also be visualised by the SPH visualisation tool SPLASH (Price, D. J. 2007, PASA, 24, 159; 
https://users.monash.edu.au/~dprice/splash/), and can be analysed more toroughly with PLONS.


PHANTOM
-------

PHANTOM returns (full) dump files every certain timestep during the evolution of the 
model. These are the wind_xxxxx-files and contain the relevant data (position, velocity, 
mass, density, energy, temperature) of the SPH particles in the model. The last x rows of 
the output are the x sink particles in the model. More detailed information of the sink 
particles (position, velocity, mass, accreted mass, spin, angular momentum) in function of 
evolution time of the model, can be retrieved in the .ev files.

PLONS - General
------------------

This pipeline is suited for single, binary and triple AGB wind models.

The following info can be attained:

(1) 2D slice plots of the global structure of the last dump full of the model.
(2) 1D line plots (radial structure) of the global structure of the last dump of the model along the x-, y- and z-axes.
(3) Information about the velocity related quantities of the model + Quantitative measurement of the degree of aspherical morphology: morphological parameters eta, Qp and epsilon.
(4) Cummulative mass fraction in function of the polar coordinate theta.
(5) Information of the orbital evolution.
(6) 1D spherical profiles for single star models



Pipeline - More detailed information
------------------------------------

This pipeline loads the wind.in and wind.setup files first. This is the general setup information 
of the model (configuration of the system, information of the thermodynamics, evolution time, wind 
outflow setup,...). Next the last full dump is loaded, and using this data, some other useful 
quantities are calculated (gas pressure/temperature, speed, sperical coordinates, sound speed of 
the gas,...).  Lastly, the data of the sink particle(s) is loaded. From this file, also the period 
and orbital velocity of the sink particles is calculated, if the system is a multiple, in function of 
time. The last entry in the sink data, corresponds to the full dump loaded. Therefore the 
position/velocity/mass/... of this last entry are also saved in the full dump data for easy usage.


(1)

The 2D slices in the plots are infinitely thin and calculated by using the smoothing kernel from
PHANTOM. The density, speed and temperature are calculated by this method on a grid of n*n 
pixels, using the data of the m nearest neighbours (SPH particles) of those grid points. These 
plots give a visual representation of the morphology present in the AGB wind of the simulation.


(2)

In 1D, the data along the x-, y- and z-axes is constructed in the same way as explained in (1), 
but for a 1D grid. This results in plots of the radial structure, which is an important addition 
to the 2D slice plots.


(3)

The terminal velocity in PHANTOM SPH code is not an input parameter, so it needs to be calculated. 
Due to the morphology in the outflow of the AGB star, at a certain radius a wide range of speeds 
will be obtained. Therefore calculating only one value for the terminal velocity is not possible. 
Instead, three values are obtained, a minimum, mean and maximum terminal velocity. The method used 
here is binning the speed in function of the radius of the SHP particle to the AGB star. Bins of 
1AU are constructed and per bin, the minimum, mean and maximum speed is calculated. Using speed of 
the outer 20% of the data, three values for the terminal velocity are obtained. 

NOTE: For some models (most often with slow initial wind velocity, low mass loss rate and/or low-mass 
companion), the interaction of the companion has not reached the outer boundary of the model. Therefore, 
this method will not lead to a correct terminal velocity value and thus the unphysical part of the data 
will be cut off. Visually, the speed profile in function of radius shows a strong, unphysical linear 
increase, starting with a clear twist.


In order to get some quantitative indication about the morphology of the models, several morphology 
parameters are currently in use: 
    - eta = v/v_orb (see Saladino, M. I., Pols, O. R., van der Helm, E., Pelupessy, I., & 
      Portegies Zwart, S. 2018, A&A, 618, A50; El Mellah, I., Bolte, J., Decin, L., Homan, 
      W., & Keppens, R. 2020, arXiv e-prints, arXiv:2001.04482)
    - Qp = p_comp/p_wind 
         = (M_comp v_orb) / (M_wind v_wind)
        (see Decin, L., Montarges, M., Richards, A. M. S., et al. 2020, Science, 369, 1497)
    - epsilon = e_kin/e_grav
              = (v_wind^2 a)/(24 G^3 M_comp^2 M_AGB)^(1/3)
    - Rcapt/a = (2 G M_comp)/(v_wind^2) * 1/a
Depending on the value of these parameters, the model is expected to show radial/EDE/complex morphology.

Rhill = a(1-e)(M_comp/ 3M_AGB)^(1/3)
Therefore, Rhill/Rcapt = epsilon, if for both radii the same v_wind is used.

VERY IMPORTANT NOTE: 
The parameter 'v_wind' is not unambiguously defined. Even more, if only a binary models is used, 
this parameter is most likely not to be constrained properly. For the different morphological parameters, 
multiple values for v_wind are used.
For eta, different velocities/speed v are used:
    - terminal velocity
    - speed of the wind at the location of the companion 
For Qp and epsilon, two different values are used as v_wind:
    - the mean velocity of the wind at the location is used as calculated from the binning
    - the average of the min and max velocity at the location of the companion is used
For Rcapt/a the initial wind velocity is used to get a rough indication.

Better options for v_wind would be to use the velocity of the corresponding single model as follows:
    - speed of the wind at the location of the companion
    - speed resulting from the vector sum of the speed of the wind at the location of the companion 
      and the orbital velocity of the AGB star.
These can be calculated using the output of this pipeline.


(4)

The cummulative mass fraction (CMF) is calculated in function of the polar angle theta (theta = 0.5pi 
is the orbital plane, theta = 0 the north pole), again by using the method of binning. Only the northern 
hemisphere is used, because of symmetry. Because of asymmetry along the x-axis, the positive (right) 
and negative (left) part according to x are also calculated seperately. Certainly for eccentric models 
this can differ. To calculate the CMF, only the selfsimilar part of the outflow it used, therefore a 
factor (called factor) is given. The inner factor * semi-major axis AU is left out in this calculation.

From this CMF, the angle is calculated where 25, 50 and 75% of the mass of the wind is present, called 
theta25/50/75 respectively. From this, the parameter delta is calculated, defined as 
(theta25-theta50)/(theta50-theta75). When this parameter is normalised to the value for the corresponding 
single model, it gives a measure for the degree of EDE.

Also the mean density is calculated using a similar approach as the CMF, using the method of binning. 
This mean density can easily be normalised to compare the morphology in different models.


(5)

In order to get information about the orbital evolution, the sink files are used (.ev), which gives
relevant quantities of the sink particles in function of time. This file returns plots of the orbit, 
evolution of orbital separation, orbital velocity and accreted mass of the companion.





