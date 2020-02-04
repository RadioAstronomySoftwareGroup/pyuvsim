Full class list and documentation
=================================

Primary Beams
-------------

Simulations may be run using pyuvdata UVBeam objects or AnalyticBeam objects.


.. autoclass:: pyuvsim.AnalyticBeam
    :members:

Antenna objects
---------------

.. autoclass:: pyuvsim.Antenna
    :members:

Profiling
---------

Developers will want to use line profiling tools to track down bottlenecks and
find room for performance improvements. These tools better allow profilers to
operate within the MPI environment.

.. automodule:: pyuvsim.profiling
    :members:

Simulation setup
----------------

Simulations are run from .yaml files with specified keywords. See :ref:`usage`.

.. automodule:: pyuvsim.simsetup
    :members:

Sources
-------

See documentation for `pyradiosky`.

Telescopes
----------

Shared properties of all antennas, such as array center location and primary beam objects.

.. autoclass:: pyuvsim.Telescope
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
