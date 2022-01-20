Full class list and documentation
=================================

Primary Beams
-------------

Simulations may be run using :class:`pyuvdata.UVBeam` objects or :class:`pyuvsim.AnalyticBeam` objects.


.. autoclass:: pyuvsim.AnalyticBeam
    :members:

Antenna objects
---------------

.. autoclass:: pyuvsim.Antenna
    :members:

Baseline objects
----------------
.. autoclass:: pyuvsim.Baseline
    :members:


Telescopes
----------

Shared properties of all antennas, such as array center location and list of primary beam objects.

.. autoclass:: pyuvsim.Telescope
    :members:


.. autoclass:: pyuvsim.BeamList
    :members:

Sources
-------

See documentation for :class:`pyradiosky.SkyModel`.


Simulation setup
----------------

Simulations are run from .yaml files with specified keywords. See :ref:`usage`.

.. automodule:: pyuvsim.simsetup
    :members:

UV Simulation functions
-----------------------

Methods for running simulations.

.. automodule:: pyuvsim.uvsim
    :members:


MPI Tools
---------

Tools to initialize MPI and share data among processors.

.. automodule:: pyuvsim.mpi
    :members:

Profiling
---------

Developers will want to use line profiling tools to track down bottlenecks and
find room for performance improvements. These tools better allow profilers to
operate within the MPI environment.

.. automodule:: pyuvsim.profiling
    :members:
