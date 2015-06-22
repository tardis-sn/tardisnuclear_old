***************
Getting started
***************

The first useful thing is to create a lightcurve model for a Ni56 powered
explosion. Here the values for the different are given in units of solar
masses::

    >>> from tardisnuclear import models
    >>> sn_model = models.make_bolometric_model(ni56=1.0)
    >>> sn_model
    <BolometricLightCurve(ni56=1.0)>

It is however very likely that this will fail due to the inavailabilty of the
necessary nuclear data. The nuclear data.
