.. _usage:

Running a simulation
====================

``pyuvsim`` takes a UVData object, instrument configuration settings, and a source catalog and effectively "fills" the UVData object with simulated data. The function ``uvsim.run_pyuvsim`` in uvsim.py accepts as input a path to a yaml ``obsparam`` file and writes out the data to ``uvfits``, ``miriad``, or ``uvh5`` format. Optionally, it can also skip writing the data out and just return a ``UVData`` object filled with simulated data. The default behavior is to write to ``uvfits``.

This is demonstrated in more detail in ``run_param_pyuvsim.py`` in the scripts directory. See :doc:`parameter_files` for information on the parameter files.

Under the hood, pyuvsim generates a ``pyuvdata.UVData`` object without data and then fills it with simulated data. The function ``pyuvsim.run_uvdata_uvsim`` provides this lower-level functionality if needed.

Using MPI
^^^^^^^^^

``pyuvsim`` is parallelized using the Message Passing Interface (MPI). To take full advantage of this, any wrapper must be run with ``mpirun``:

    .. code-block:: python

        # Running with 50 MPI processing units
        > mpirun -n 50 python run_param_pyuvsim obsparam_filename.yaml   # This will run a parameter file job with 10 processing units.


Further speedup is achieved through ``numpy``/``scipy`` internal threading. How effective this is depends on the linear algebra library that ``numpy`` is compiled on, which can be checked with ``numpy.show_config()``.

Enabling Profiling
^^^^^^^^^^^^^^^^^^

The ``line_profiler`` module provides runtime estimates on a line by line basis. It is built into ``pyuvsim`` to work within the MPI framework using the functions in ``pyuvsim.profiling``. To run a simulation with profiling enabled, run the command ``profiling.set_profiler()`` before starting ``run_uvsim()``. This function can be passed a list of functions you wish to profile (by name), as well as which rank to return data for (it will only profile one MPI rank at a time!).

Generating Config Files from Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The scripts ``uvfits_to_telescope_config.py`` and ``uvfits_to_config.py`` are provided for convenience. These will accept the path to any valid uvfits file as input, along with output file name options, and will generate telescope layout, telescope configuration, and an obsparam file for a simulation with the same time, frequency, and baseline structure as the data file.

Note that the generated configuration files will still need to be given paths to beam models and source catalogs before they can be used.


Parallelization
^^^^^^^^^^^^^^^

Given ``Npus`` MPI processes, and ``Nsrcs``, ``Nbls``, ``Ntimes``, ``Nfreqs`` sources, baselines, times, and frequencies (respectively), the choice of splitting the various axes goes like this:

    1. If ``Npus`` < ``Nsrcs`` and ``Npus > Nbls * Ntimes * Nfreqs``:
           a. The source axis will be split across MPI processes.
           b. Each process will receive all baselines, times, and frequencies.
    2. Otherwise:
           a. Each process will receive all sources, and times/frequencies/baselines will be split among MPI processes.
           b. If it is estimated that loading all sources on each process simultaneously will overrun the available memory, then the sources will be split into chunks for processing.

In each case, the source axis is handled through ``numpy``'s threading. It's recommended that jobs with especially large source axes be given more cores per MPI process (in SLURM, for instance, this is set by the ``--cpus-per-task`` option in ``sbatch`` or ``srun``). Usually around 2 to four cpus per process is sufficient. For large numbers of times/baselines/frequencies, however, running with more MPI processes offers a better speedup.

The source array is initially shared in immutable shared memory, and parts of it are copied for usage within each MPI process. Likewise, UVBeam-class beam objects are loaded on the root process only, and broadcast using shared memory. These measures prevent large data arrays from being copied over ``Npus`` times, which can cause an unacceptable memory bloat.
