.. _link-installation:

Installation
############

The easiest way to install plons is via a package installer like pip or conda. Use one of the following commands:

.. code-block:: shell
    
    pip install plons

.. code-block:: shell
    
    conda install -c ensor plons

In this way you will have the most recent stable version.

An other way to install plons is using the `GitHub <https://github.com/ensor-code/plons>`_ repository. Clone the repository, enter the directory, and install the package using:

.. code-block:: shell

    git clone git@github.com:Ensor-code/plons.git
    cd plons
    pip install .

If you have the github repository on your computer, you can also access our custom pipeline. Add to your .bashrc the following function:

.. code-block:: shell

    function plons() {
    (cd "path to plons"; python main.py $@) && notify-send "Plotting Tool" "Rendering figures done" -i "path to plons"/plons.png
    }

now you can use the plons pipeline as an executable