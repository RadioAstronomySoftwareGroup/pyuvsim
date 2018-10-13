
Parameter and configuration Files
===================================

When running simulations from yaml and csv files, there are four configuration files that must be used.

The outermost parameter file is the `obsparam_*.yaml`, which is parsed by ``initialize_uvdata_from_params`` into a UVdata object, a source list, a beam dictionary, and a beam list.

The antenna layout and telescope config yaml files determine the full properties of the array, including location, beam models, layout, and naming.

The catalog text files give point source lists. 


These files contain overall simulation parameters.
Passed into ````run_param_pyuvsim.py````

.. code-block:: yaml
    :caption: Example obsparam yaml file

    filing:
      outdir: '.'   #Output file directory
      outfile_prefix: 'sim_' # Prefix for the output uvfits file.
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
      catalog: 'mock'       #Choice of catalog file name (or specify mock to use builtin catalogs).
      mock_arrangement: 'zenith'    # Choosing a mock layout
    telescope:
      array_layout: 'triangle_bl_layout.csv'    # Antenna layout csv file
      telescope_config_name: '28m_triangle_10time_10chan.yaml'  # Telescope metadata file.
    time:
      Ntimes: 10        # Number of times.
      integration_time: 11.0  # Time step size  (seconds)
      start_time: 2457458.1738949567    # Start and end times (Julian date)
      end_time: 2457458.175168105  
      duration_hours: 0.0276

**Note** The example above is shown with all allowed keywords, but many of these are redundant. This will be further explained below.


Time and frequency
^^^^^^^^^^^^^^^^^^

    Time and frequency structure may be defined with different combinations of keywords to suit the user's purposes. The user must specify sufficient information for the frequency and time arrays to be defined.
    
    Minimum frequency requirements:
        specify bandwidth via one of the following combinations:
            * (``start_freq``, ``end_freq``)
            * (``channel_width``, ``Nfreqs``)
            * (``bandwidth``)
    
        specify channel width via:
            * (``bandwidth``, ``Nfreqs``)
            * (``channel width``)
    
        specify a reference frequency via:
            * (``start_freq``)
            * (``end_freq``)
    
    As long as one of the sets from each category above is met by the supplied keywords, the frequency array will be successfully built.
    Likewise for the time array:
    
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
    
    You can also just give an explicit ``time_array`` or ``freq_array``.

    
Telescope Configuration
^^^^^^^^^^^^^^^^^^^^^^^

    Under the telescope section, the keywords ``array_layout`` and ``telescope_config_name`` give paths to, respectively, the array layout text file and the telescope metadata configuration yaml file. If the full path is not specified, pyuvsim will look in the same directory as the obsparam file.

    Example array layout:

    .. literalinclude:: ../pyuvsim/data/test_config/triangle_bl_layout.csv
        :end-before: 3

    Columns here provide, in order from left to right, the antenna name, antenna number, a beam ID number, and the antenna positions relative to the array center in eastings/northings/up-ings in meters. The layout file has a corresponding telescope metadata file, shown below:

    .. literalinclude:: ../pyuvsim/data/test_config/28m_triangle_10time_10chan.yaml
        :end-before: 3

    This yaml file provides the telescope name, location in latitude/longitude/altitude in degrees/degrees/meters (respectively), and the `beam dictionary`.

    The beam dictionary, along with the beam IDs in the layout file, allow one to specify different beam models be used for different antennas. The telescope configuration yaml file lists beam models to be used along with a "Beam ID". The corresponding Beam IDs in the layout file then tells pyuvsim which beam model to assign to each antenna.

Sources
^^^^^^^
    Specify the path to a text catalog file via ``catalog``.
    
    An example catalog file:

    .. literalinclude:: ../pyuvsim/data/mock_catalog_heratext_2458098.27471265.txt
        :end-before: 3
   
    The columns are:
        * ``SOURCE_ID`` : Identifier for the source
        * ``RA_J2000`` : Right ascension of source at J2000 epoch, in decimal degrees.
        * ``DEC_J2000`` : Declination of source at J2000 epoch, in decimal degrees.
        * ``FLUX``: Source stokes I brightness in Janskies.  (Currently only point sources are supported).
        * ``Frequency``: A reference frequency for the given flux. This will be used for spectral modeling. 
    
    Alternatively, you can specify a ``mock`` and provide the ``mock_arrangement`` keyword to specify which mock catalog to generate. Available options are shown in the ``create_mock_catalog`` doctstring:

    .. module:: pyuvsim

    .. autofunction:: create_mock_catalog


