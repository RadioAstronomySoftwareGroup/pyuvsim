For Developers
==============

The Radio Astronomy Software Group welcomes contributions from anyone. Full contribution guidelines and code of conduct may be found on the `RASG Github Pages <https://radioastronomysoftwaregroup.github.io>`_.

This page discusses more specific points to help developers contribute to ``pyuvsim``.

Computing Resources
-------------------

Instrument simulation is a computationally-demanding task, and ``pyuvsim`` is designed to run a brute-force calculation that requires much computing power. There is a continual effort to improve in performance, but for the foreseeable future the best results can be obtained using large, distributed computing clusters.

Much of the development and testing of ``pyuvsim`` was done on the Oscar computing cluster at Brown University's Center for Computation and Visualization (CCV_). This cluster is managed with the Simple Linux Utility Resource Manager (SLURM). Though ``pyuvsim`` can be run on most Linux or Mac-OS systems, its operation with large numbers of times/baselines/frequencies/sources or especially complicated beam objects has not been tested well. Testing within other cluster management environments is encouraged.

.. _CCV: https://docs.ccv.brown.edu/oscar/


Running and Comparing Reference Simulations
-------------------------------------------

The *reference simulations* comprise a set of standard simulations with well-defined, simple parameters that most instrument simulation tools should be capable of running. As such, they represent a standard point of comparison between ``pyuvsim`` and other simulators. They are also used to track that ``pyuvsim`` produces consistent results over time.

The first generation reference simulations run automatically upon pull request to ensure consistent runtime and output. For some types of pull requests, we ask that the 2nd generation reference simulations be re-run with the new code to confirm that any changes to results are detected. The older first generation reference simulation files are stored on the `Brown Digital Repository <https://repository.library.brown.edu/studio/collections/bdr:wte2qah8/>`_. Some old second generation reference simulations are available from this `Google Drive <https://drive.google.com/drive/folders/14hH-zBhHGddVacc0ncqRWq7ofhGLWfND?usp=drive_link>`_. If a change to the reference simulations is expected, it should be documented in the CHANGELOG and the new simulation files given to the project managers to store on the Brown Digital Repository. This process is documented in `Benchmarking <https://pyuvsim.readthedocs.io/en/latest/developers.html#benchmarking>`_. Reference simulation files are stored on the BDR with name corresponding to the name of the simulation. All versions of the output from the same reference simulation will have the same name, and their order is distinguishable by upload date, with the latest uploaded simulation of the same name taken as the standard for comparison. Specific information about an individual run can be found by loading the saved uvh5 file and viewing the attributes of interest (for example the history attribute of the object to see the versions of pyuvsim and pyuvdata used to run the simulation).

More information on checking the reference simulations can be found in the README in the ``reference_simulations`` directory in the repository.

We will periodically issue a new set of reference simulations as ``pyuvsim`` gains more capabilities.

For more details, see `reference_simulations/README.md <https://github.com/RadioAstronomySoftwareGroup/pyuvsim/tree/main/reference_simulations>`_.

Benchmarking
------------

Benchmarking Simulations
~~~~~~~~~~~~~~~~~~~~~~~~

The ``benchmarking`` directory contains tools to test the runtime and memory usage of large simulations. There is no requirement to check benchmarks for pull requests, but it's a good idea to make sure changes don't drastically alter the runtime. The file BENCHMARKS.log keeps a record of performance over time.

The README file in the ``benchmarking`` directory gives more details on how to do benchmarking.

Note that the benchmarking scripts are designed only for SLURM systems.

For more details, see `benchmarking/README.md <https://github.com/RadioAstronomySoftwareGroup/pyuvsim/tree/main/benchmarking>`_.

Running a Reference Simulation with pytest-benchmark
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To run a single core regression test of the reference simulations, you need to specify a reference simulation with the ``refsim`` flag and use ``benchmark-only``. Additionally, you need to use mpiexec to run pytest as follows:

    .. code-block:: python

        # use mpiexec to run pytest specifying one core
        > mpiexec -n 1 -np 1 pytest --refsim=1.1_baseline_number --benchmark-only

Here "1.1_baseline_number" would be the specific reference simulation being tested. You can use the ``refsim`` flag multiple times to parametrize multiple reference simulations: ``--refsim=refsim1 --refsim=refsim2``. Multiple simulations can also be specified with ``--refsim={refsim1,refsim2,...}``.

We run single core regression tests of the available reference simulations with pytest and pytest-benchmark via our github ci workflow on every pull request. We do so to ensure output and runtime consistency. As we only run the simulations with a single core, the benchmarking aspect of these tests is only relevant for linear operations and not a test of any parallelism. Historical simulation output is currently stored on the Brown Digital Repository.

The available ``refsim`` values are:

* 1.1_baseline_number
* 1.2_time_axis
* 1.3_frequency_axis
* 1.4_source_axis
* 1.5_uvbeam
* 1.6_healpix
* 1.7_multi_beam
* 1.8_lunar

These values, and additional custom options defined for pytest such as ``--savesim``, can be seen under the custom options section of the output of ``pytest --help``.

To merge a pull request that changes the output of the reference simulations, you should generate a new set of reference simulation outputs to upload to the Brown Digital Repository. The line to do this currently is:

    .. code-block:: python

        # use mpiexec to run pytest specifying one core
        > mpiexec -n 1 -np 1 pytest --refsim={1.1_baseline_number,1.2_time_axis,1.3_frequency_axis,1.4_source_axis,1.5_uvbeam,1.6_healpix,1.7_multi_beam,1.8_lunar} --benchmark-only --savesim

This will run all the active reference simulations, and write the new reference simulation output out in a ``new_data`` directory located in the current working directory. The files in the ``new_data`` directory should then be uploaded to the Brown Digital Repository to serve as the updated reference simulations. If you just want to run the reference simulations locally but not save anything, the ``--savesim`` flag can be dropped -- in this case it can be nice to use the ``-s`` flag to see all printed output.
