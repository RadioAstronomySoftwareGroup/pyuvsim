Full class list and documentation
=================================

Primary Beams
-------------

Simulations may be run using pyuvdata UVBeam objects or AnalyticBeam objects.


.. automodule:: pyuvsim.analyticbeam
    :members:

Antenna objects
--------------

.. automodule:: pyuvsim.antenna
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

.. automodule:: pyuvsim.source
    :members:

Telescopes
----------

Shared properties of all antennas, such as array center location and primary beam objects.

.. automodule:: pyuvsim.telescope
    :members:

UV Simulation functions
-----------------------

Methods for setting up the simulation, handling concurrent MPI processes, and saving computed visibilities.

.. automodule:: pyuvsim.uvsim
    :members:
