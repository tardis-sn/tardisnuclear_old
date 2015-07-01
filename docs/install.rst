************
Installation
************

``tardisnuclear`` has several dependencies that are easiest installed via
 `Anaconda <http://continuum.io/downloads>`_ and all further instructions will
 assume that you use that for your installation. 

We will first clone the github repository of tardisnuclear::

    git clone https://github.com/tardis-sn/tardisnuclear.git
    cd tardisnuclear

Next we will create an environment that has exactly the libraries needed for
tardisnuclear to run::

    conda create -n tardisnuclear python=2 --file conda-requirements

Then activate it::

    source activate tardisnuclear


A further important packate that needs to be installed is pyne::

    conda install -c cyclus -c pyne pyne=0.5.0

and then we are ready to install tardisnuclear with::

    python setup.py install
