Install
=======

Via pip
-------

Run

    pip install oospectro


Note about dependencies
-----------------------

* numpy
* scipy
* scikit-image
* matplotlib

Installing scipy >= 1.1.0 will provide better performances.
If scipy < 1.1.0 is install, then a non-optimized copy of scipy algorithms is used internally.

The most experienced users can use a virtualenv to get a recent scipy version.
To install scipy from master, use

    pip install git+https://github.com/scipy/scipy.git

