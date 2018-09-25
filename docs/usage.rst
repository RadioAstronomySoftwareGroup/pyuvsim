.. _usage:
Running a simulation
====================

``pyuvsim`` takes a UVData object, instrument configuration settings, and a source catalog and effectively "fills" the UVData object with simulated data. The function ``run_pyuvsim`` in uvsim.py


From a uvfits file
^^^^^^^^^^^^^^^^^^

This is demonstrated in ``run_pyuvsim.py`` in the scripts directory.

To simulate from a uvfits file:
    1. Read the file into a UVData object (``input_uv``)
    2. Define a ``beam_list`` and a ``beam_dict``.
    3. Provide the path to a catalog file, or specify a mock catalog.
    4. Run ``uvsim.run_uvsim()``:
        ::

             output_uv = run_pyuvsim(input_uv, beam_list, beam_dict=beam_dict)


From parameter files
^^^^^^^^^^^^^^^^^^^^

This is demonstrated in more detail in ``run_param_pyuvsim.py`` in the scripts directory. See :doc:`parameter_files` for information on the parameter files.

To simulate from parameter files:
    1.  Use ``initialize_uvdata_from_params`` to obtain a uvdata object, beam_list, beam_dict, and beam_ids:
        ::

            input_uv, beam_list, beam_dict, beam_ids = simsetup.initialize_uvdata_from_params(params)

        This function can accept ``params`` as either a parameter dictionary or a filename for a yaml file.
    2. Provide a catalog file path.
    3. Run ``uvsim.run_uvsim()`` as before.


Using MPI
^^^^^^^^^

``pyuvsim`` is parallelized using the Message Passing Interface (MPI). To take full advantage of this, any wrapper must be run with ``mpirun``:
    ::

        > mpirun -n 50 python run_param_pyuvsim -p obsparam_filename.yaml   # This will run a parameter file job with 10 processing units.
    
