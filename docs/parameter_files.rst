Parameter and configuration Files
=================================

When running simulations from yaml and csv files, there are four configuration files
that must be used.

The outermost parameter file is the `obsparam_*.yaml`, which is parsed by
``initialize_uvdata_from_params`` into a UVdata object, a source list, a beam dictionary,
and a beam list.

The antenna layout and telescope config yaml files determine the full properties of the
array, including location, beam models, layout, and naming.

The catalog text files give point source lists.


These files contain overall simulation parameters.
Passed into ``run_param_pyuvsim.py``

.. code-block:: yaml
    :caption: Example obsparam yaml file

    filing:
      outdir: '.'   #Output file directory
      outfile_prefix: 'sim' # Prefix for the output file name, separated by underscores
      outfile_suffix: 'results' # Suffix for output file name
      outfile_name: 'sim_results' # Alternatively, give the full name
      output_format: 'uvfits'  # Format for output. Default is 'uvfits', but measurement sets ('ms'), 'miriad' and 'uvh5' are also supported.
      clobber: False        # overwrite existing files. (Default False)
    freq:
      Nfreqs: 10    # Number of frequencies
      channel_width: 80000.0    # Frequency channel width
      end_freq: 100800000.0     # Start and end frequencies (Hz)
      start_freq: 100000000.0
      freq_array : [1.0000e+08,   1.0008e+08,   1.0016e+08, 1.0024e+08,
          1.0032e+08,   1.0040e+08,   1.0048e+08, 1.0056e+08,
          1.0064e+08, 1.0072e+08]
      bandwidth: 800000.0
    sources:
      catalog: '../pyuvsim/data/gleam_50srcs.vot'   # Path to catalog file (txt, vot, hdf5, etc.) readable with pyradiosky.
      spectral_type: flat # If using the GLEAM catalog, specify the spectral type (flat, subband or spectral_index). Defaults to flat.
      table_name: single  # Required for non-GLEAM VO table files
      id_column: name  # Required for non-GLEAM VO table files
      flux_columns: Si  # Required for non-GLEAM VO table files
      ra_column: RAJ2000  # Recommended for non-GLEAM VO table files
      dec_column: DEJ2000  # Recommended for non-GLEAM VO table files
      catalog: 'mock'       # Alternatively, use 'mock' to use a builtin catalog).
      mock_arrangement: 'zenith'    # If using the mock catalog, specify which one. Additional mock keywords are specified here.
    telescope:
      array_layout: 'triangle_bl_layout.csv'    # Antenna layout csv file
      telescope_config_name: '28m_triangle_10time_10chan.yaml'  # Telescope metadata file.
    time:
      Ntimes: 10        # Number of times.
      integration_time: 11.0  # Time step size  (seconds)
      start_time: 2457458.1738949567    # Start and end times (Julian date)
      end_time: 2457458.175168105
      duration_hours: 0.0276
    select: # limit which baselines are simulated. Use any UVData.select keywords (except polarizations) and/or redundant_threshold
      bls: [(1, 2), (3, 4), (5, 6)]
      ant_str: 'cross'
      antenna_nums: [1, 7, 9, 15]
      redundant_threshold: 0.1 # redundancy threshold in meters. Only simulate one baseline per redundant group

**Note** The example above is shown with all allowed keywords, but many of these are
redundant. This will be further explained below. Only one source catalog will be used
at a time.

Filing
^^^^^^
    Specifies where the results file will be output, what name the file should have,
    and whether or not to overwrite existing files. None of these parameters are required.

Frequency
^^^^^^^^^

    As is the standard with ``pyuvdata``, the frequency channel numbers refer to the
    **frequency at the channel center**. The ``bandpass`` refers to the total band
    covered by all channels, and the channel width is the separation between channel
    centers. Therefore, the following relations hold::

		bandpass = Nfreqs * channel_width
		bandpass = (end_freq + channel_width/2.) - (start_freq - channel_width/2.) = ( end_freq - start_freq) + channel_width
		start_freq = end_freq - bandpass + channel_width
		end_freq = start_freq + bandpass - channel_width


    Time and frequency structure may be defined with different combinations of keywords
    to suit the user's purposes. The user must specify sufficient information for the
    frequency array to be defined.

    Minimum frequency requirements:

    Specify bandwidth via one of the following combinations:

        * (``start_freq``, ``end_freq``)
        * (``channel_width``, ``Nfreqs``)
        * (``bandwidth``)

    Specify channel width via:

        * (``bandwidth``, ``Nfreqs``)
        * (``channel width``)

    Specify a reference frequency via:

        * (``start_freq``)
        * (``end_freq``)

    As long as one of the sets from each category above is met by the supplied
    keywords, the frequency array will be successfully built.
    You can also just give an explicit ``freq_array``.

    The ``channel_width`` should be specified as a scalar unless ``freq_array`` is specified,
    in which case ``channel_width`` can either be a scalar or an array of the same
    length as ``freq_array``.

    If you specify an explicit ``freq_array`` that is not evenly spaced or is only
    length one, you must specify the ``channel_width``, either as a single value (in Hz)
    or as an array of the same length as ``freq_array``.

Time
^^^^

    The time array is specified similarly. The entries in the ``time_array`` indicate the
    **center of each time step in Julian date**. The ``integration_time`` is the time
    step size in seconds. The user may also specify ``duration_hours`` or ``duration_days``
    to specify the total time covered by all time steps. The following relations among
    parameters hold::

        duration_hours = Ntimes * integration_time / (3600.)
        duration_days = duration_hours / 24.
        duration_days = (end_time - start_time) + integration_time / 86400
        start_time = end_time - duration_days + integration_time / 86400
        end_time = start_time + duration_days - integration_time / 86400

    The numerical factors are to convert among seconds, days, and hours. The user must
    specify sufficient information for the time array to be defined:

    Minimum time requirements:

    Total time:

        * (``start_time``, ``end_time``)
        * (``integration_time``, ``Ntimes``)
        * (``duration_hours``) or (``duration_days``)

    Time step:

        * (``duration_hours`` or ``duration_days``, ``Ntimes``)
        * (``integration_time``)

    Reference time:

        * (``start_time``)
        * (``end_time``)

    As long as one of the sets from each category above is met by the supplied keywords,
    the time array will be successfully built.



Telescope Configuration
^^^^^^^^^^^^^^^^^^^^^^^

    Under the telescope section, the keywords ``array_layout`` and ``telescope_config_name``
        give paths to, respectively, the array layout text file and the telescope metadata
        configuration yaml file. These path may either be absolute or specified relative
        to the location of the obsparam yaml file.

    Example array layout with four antennas:

    .. literalinclude:: example_configs/baseline_lite.csv

    Columns here provide, in order from left to right, the antenna name, antenna number,
    a beam ID number, and the antenna positions relative to the array center in
    east, north, up (ENU) in meters. The layout file has a corresponding telescope
    metadata file, shown below:

    .. literalinclude:: example_configs/bl_lite_mixed.yaml

    This yaml file provides the telescope name, location in latitude/longitude/altitude
    in degrees/degrees/meters (respectively), and the `beam dictionary`.
    In this case, beam_id == 0 is the UVBeam file hera.uvbeam, beam_id == 1
    is an Airy disk with diameter 16 m, beam_id == 2 is a Gaussian beam with
    sigma 0.03, beam_id == 3 is another Airy beam with diameter 12 m, 4 is a Gaussian
    with diameter 14 m, and 5 is a Gaussian with with diameter 12 m. When specifying a shape
    parameter for a specific beam_id, the beam type needs to be specified using
    the type keyword (rather than on the same line as the beam_id) and then
    the shape keyword can be specified in the next line at the same indent level.
    When no shape parameter is added in the beam_dictionary (as with 3), pyuvsim
    will look for a default parameter below. So in this case, the beam_id == 3
    ends up with a diameter of 12 m. The dictionary only needs to be as long as
    the number of unique beams used in the array, while the layout file specifies
    which antennas will use which beam type. This allows for a mixture of beams
    to be used, as in this example. Unassigned beams will be ignored (the given
    layout file does not use beams 2 or 3).

    `freq_interp_kind` sets the type of frequency interpolation for all UVBeam objects
    defined in the beam list (see documentation on UVBeam for options).

    When used with version of pyuvdata >= 2.0.1, the `spline_interp_opts` keyword lets
    the user set the order on the angular interpolating polynomial spline function. By default,
    it is cubic.

    Analytic beams may require additional parameters.

    - uniform = The same response in all directions. No additional parameters.
    - gaussian = Gaussian function shaped beam. Requires either an antenna diameter
      (in meters) or a standard deviation sigma (in radians). This standard deviation sets
      the width of the beam in zenith angle. Note that defining gaussian beams via `sigma`
      will be deprecated in the future.
    - airy = Airy disk (ie, diffraction pattern of a circular aperture). Requires an
      antenna diameter.

    Note that beams defined with an antenna diameter will be chromatic (their widths on
    the sky will change with frequency).

    The figure below shows the array created by these configurations, with beam type
    indicated by color.

    .. image:: Images/baseline_lite.png
	    :width: 600
	    :alt: Graphical depiction of the example antenna layout.

Telescopes on the Moon
~~~~~~~~~~~~~~~~~~~~~~
   If the ``lunarsky`` module is installed, the ``telescope_location`` can be interpreted as the
   lon/lat/alt of an observatory on the Moon, defined in the "Mean Earth/ Mean Rotation"
   frame (see documentation on ``lunarsky``). Setting the keyword ``world: moon`` in the
   telescope_config file enables this:

   .. literalinclude:: example_configs/tranquility_config.yaml


Sources
^^^^^^^
    Specify the path to a text catalog file via ``catalog``. The path can be given as an
    absolute path or relative to the location of the obsparam. This catalog should be
    readable with `pyradiosky`.

    An example catalog file:

    .. literalinclude:: ../pyuvsim/data/mock_catalog_heratext_2458098.27471265.txt
        :lines: 1-5

    The columns are:

        * ``SOURCE_ID`` : Identifier for the source
        * ``RA_J2000`` : Right ascension of source at J2000 epoch, in decimal degrees.
        * ``DEC_J2000`` : Declination of source at J2000 epoch, in decimal degrees.
        * ``FLUX``: Source stokes I brightness in Janskies.  (Currently only point sources are supported).
        * ``Frequency``: A reference frequency for the given flux. This will be used for spectral modeling.

    If the catalog is a GLEAM VO table file, optionally specify the ``spectral_type``
    as one of: ``flat``, ``subband`` or ``spectral_index``. If not specified it defaults to ``flat``.

    If the catalog is a different VO table file, several other keywords are required or recommended:
      * ``table_name`` : The name of the table to use from the file (required).
      * ``id_column`` : The name of the column to use for the source IDs (required).
      * ``flux_columns`` : One or a list of columns to use for the source fluxes
        (a list for fluxes at multiple frequencies) (required).
      * ``ra_column`` : The name of the column to use for the source RAs (recommended, defaults to ``RAJ2000``).
      * ``dec_column`` : The name of the column to use for the source Decs (recommended, defaults to ``DEJ2000``).

    Alternatively, you can specify a ``mock`` and provide the ``mock_arrangement``
    keyword to specify which mock catalog to generate. Available options are shown
    in the :func:`pyuvsim.simsetup.create_mock_catalog` docstring.

    Flux limits can be made by providing the keywords ``min_flux`` and ``max_flux``.
    These specify the min/max stokes I flux to choose from the catalog.

    The option ``horizon_buffer`` can be set (in radians) to adjust the tolerance on the
    coarse horizon cut. After reading in the catalog, ``pyuvsim`` roughly calculates the
    rise and set times (in local sidereal time, in radians) for each source. If the
    source never rises, it is excluded from the simulation, and if the source never sets
    its rise/set times are set to None. This calculation is less accurate than the
    astropy alt/az calculation used in the main task loop, so a "buffer" angle is added
    to the set lst (and subtracted from the rise lst) to ensure sources aren't
    accidentally excluded. Tests indicate that a 10 minute buffer is sufficient.
    Pyuvsim also excludes sources below the horizon after calculating their AltAz
    coordinates, which is more accurate. The coarse cut is only to reduce computational load.

Select
^^^^^^
    Specify keywords to select which baselines to simulate. The selection is done by
    UVData.select, so it can accept any keyword that function accepts, except ones that
    affect polarization because pyuvsim computes all polarizations.

    In addition to the UVData.select keywords, a ``redundant_threshold`` parameter can
    be specified. If it is present, only one baseline from each set of redundant
    baselines is simulated. The ``redundant_threshold`` specifies how different two
    baseline vectors can be to still be called redundant -- the magnitude of the vector
    differences must be less than or equal to the threshold. The vector differences are
    calculated for a phase center of zenith (i.e. in drift mode).
