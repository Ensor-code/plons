.. _link-installation-pipeline:

Installation Pipeline
#####################

Next to the Plons package, an automated pipeline is provided. To work with this pipeline, the git repository needs to be cloned.

Within the git repository, there is a directory ``scripts``. In this directory the pipeline is situated, called ``main.py``.

To work with this pipeline, the script needs to be executed within this directory, or linking to this directory via

::

    (cd /home/matse/codes/plons/scripts; python main.py $@)

The first time to run this pipeline, a few questions will be asked:

::

    Please enter the prefix that has been used to name all the PHANTOM data files

representing the prefix given for the phantom simulations (prefix.in, prefix.setup, prefix_00000 ...)

::

    Enter the path where the PHANTOM outputs of the models are saved

This is the repository where the phantom models are saved. Where the pipeline will search for models

::

    Enter the path where the pipeline output should be saved

In this repository the figures constructed with the pipeline will be saved.

::

    Enter the path where PHANTOM is installed

This repository is the hard path to where phantom is installed.