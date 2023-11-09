Pipeline
########

The following examples show how to create figures using the plons pipeline.

To use the pipeline, the pipeline needs to be :ref:`installed <link-installation-pipeline>` first.

When running the pipeline, the user settings will be integrated, and models will be searched in the ``data_location``. All models found will be displayed, and the model for which the pipeline needs to run can be selected by typing the correct number. Multiple models can be selected by adding them with a space bar inbetween. When typing ``a`` or ``all``, all models will be used.

Later a selection of the pipeline can be chosen. The following info can be attained:

1. 2D slice plots of the global structure of the last dump full of the model.
2. 1D line plots (radial structure) of the global structure of the last dump of the model along the x-, y- and z-axes.
3. Information about the velocity related quantities of the model + Quantitative measurement of the degree of aspherical morphology: morphological parameters eta, Qp and epsilon.
4. Cummulative mass fraction as a function of the polar coordinate theta.
5. Information of the orbital evolution.
6. 1D spherical profiles for single star models


Pipeline - More detailed information
------------------------------------

This pipeline loads the wind.in and wind.setup files first. This is the general setup information
of the model (configuration of the system, information of the thermodynamics, evolution time, wind
outflow setup,...). Next the last full dump is loaded, and using this data, some other useful
quantities are calculated (gas pressure/temperature, speed, sperical coordinates, sound speed of
the gas,...).  Lastly, the data of the sink particle(s) is loaded. From this file, also the period
and orbital velocity of the sink particles is calculated, if the system is a multiple, as a function of
time. The last entry in the sink data, corresponds to the full dump loaded. Therefore the
position/velocity/mass/... of this last entry are also saved in the full dump data for easy usage.


(1) 2D slice plots of the global structure of the last dump full of the model.
______________________________________________________________________________

The 2D slices in the plots are infinitely thin and calculated by using the smoothing kernel from
PHANTOM. The density, speed and temperature are calculated by this method on a grid of n*n
pixels, using the data of the m nearest neighbours (SPH particles) of those grid points. These
plots give a visual representation of the morphology present in the AGB wind of the simulation.


(2) 1D line plots (radial structure) of the global structure of the last dump of the model along the x-, y- and z-axes.
_______________________________________________________________________________________________________________________

In 1D, the data along the x-, y- and z-axes is constructed in the same way as explained in (1),
but for a 1D grid. This results in plots of the radial structure, which is an important addition
to the 2D slice plots.


(3) Information about the velocity related quantities of the model + Quantitative measurement of the degree of aspherical morphology: morphological parameters eta, Qp and epsilon.
___________________________________________________________________________________________________________________________________________________________________________________

The terminal velocity in PHANTOM SPH code is not an input parameter, so it needs to be calculated.
Due to the morphology in the outflow of the AGB star, at a certain radius a wide range of speeds
will be obtained. Therefore calculating only one value for the terminal velocity is not possible.
Instead, three values are obtained, a minimum, mean and maximum terminal velocity. The method used
here is binning the speed as a function of the radius of the SHP particle to the AGB star. Bins of
1AU are constructed and per bin, the minimum, mean and maximum speed is calculated. Using speed of
the outer 20% of the data, three values for the terminal velocity are obtained.

NOTE: For some models (most often with slow initial wind velocity, low mass loss rate and/or low-mass
companion), the interaction of the companion has not reached the outer boundary of the model. Therefore,
this method will not lead to a correct terminal velocity value and thus the unphysical part of the data
will be cut off. Visually, the speed profile as a function of radius shows a strong, unphysical linear
increase, starting with a clear twist.


In order to get some quantitative indication about the morphology of the models, several morphology
parameters are currently in use:

- :math:`\eta = v/v_{orb}` (see `Saladino et al. 2018 <https://ui.adsabs.harvard.edu/abs/2018A%26A...618A..50S/abstract>`_; `El Mellah et al. 2020 <https://ui.adsabs.harvard.edu/abs/2020A%26A...637A..91E/abstract>`_)
- :math:`Q_p = p_{comp}/p_{wind} = (M_{comp} v_{orb}) / (M_{wind} v_{wind})` (see `Decin, L. et al. 2020 <https://ui.adsabs.harvard.edu/abs/2020Sci...369.1497D/abstract>`_)
- :math:`\epsilon = e_{kin}/e_{grav} = (v_{wind}^2 a)/(24 G^3 M_{comp}^2 M_{AGB})^{1/3}` (see `Maes et al. 2021 <https://ui.adsabs.harvard.edu/abs/2021A%26A...652A..51M/abstract>`_; `Malfait et al. 2021 <https://ui.adsabs.harvard.edu/abs/2021A%26A...653A..25M/abstract>`_)
- :math:`R_{capt}/a = (2 G M_{comp})/(v_{wind}^2)/a`

Depending on the value of these parameters, the model is expected to show radial/EDE/complex morphology.

.. math:: 
    R_{hill} = a(1-e)(M_{comp}/ 3M_{AGB})^{1/3}

Therefore, Rhill/Rcapt = epsilon, if for both radii the same v_wind is used.

VERY IMPORTANT NOTE:
The parameter ':math:`v_{wind}`' is not unambiguously defined. Even more, if only a binary models is used,
this parameter is most likely not to be constrained properly. For the different morphological parameters,
multiple values for :math:`v_{wind}` are used.
For eta, different velocities/speed :math:`v` are used:

- terminal velocity
- speed of the wind at the location of the companion

For Qp and epsilon, two different values are used as v_wind:

- the mean velocity of the wind at the location is used as calculated from the binning
- the average of the min and max velocity at the location of the companion is used

For :math:`R_{capt}/a` the initial wind velocity is used to get a rough indication.

Better options for :math:`v_{wind}` would be to use the velocity of the corresponding single model as follows:
- speed of the wind at the location of the companion
- speed resulting from the vector sum of the speed of the wind at the location of the companion and the orbital velocity of the AGB star.

These can be calculated using the output of this pipeline.


(4) Cummulative mass fraction as a function of the polar coordinate theta.
__________________________________________________________________________

The cummulative mass fraction (CMF) is calculated as a function of the polar angle theta (:math:`\theta =  \pi/2`)
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

(5) Information of the orbital evolution.
_________________________________________

In order to get information about the orbital evolution, the sink files are used (.ev), which gives
relevant quantities of the sink particles as a function of time. This file returns plots of the orbit,
evolution of orbital separation, orbital velocity and accreted mass of the companion.


(6) 1D spherical profiles for single star models
________________________________________________

In order to check your simulations and compare the 3D simulations to a 1D solution, this routine reads wind_1D.dat to plot the 1D solution, as well as the 3D sph values.

