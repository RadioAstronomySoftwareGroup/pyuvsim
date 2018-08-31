# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 2-clause BSD License

import numpy as np
import yaml
import os
import sys
import warnings
import astropy.units as units
from astropy.time import Time
from astropy.io.votable import parse
from astropy.coordinates import Angle, SkyCoord, EarthLocation, AltAz

from pyuvdata import UVBeam, UVData
import pyuvdata.utils as uvutils

from .source import Source
from .analyticbeam import AnalyticBeam
from .mpi import get_rank
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH


def strip_extension(filepath):
    if '.' not in filepath:
        return filepath, ''
    file_list = filepath.split('.')
    return ".".join(file_list[:-1]), '.' + file_list[-1]


def check_file_exists_and_increment(filepath):
    """
        Given filepath (path + filename), check if it exists. If so, add a _1
        at the end, if that exists add a _2, and so on.
    """
    if os.path.exists(filepath):
        filepath, ext = strip_extension(filepath)
        if not filepath.endswith("_0"):
            filepath += "_0" + ext
        else:
            filepath += ext
    else:
        return filepath
    n = 1
    while os.path.exists(filepath):
        filepath, ext = strip_extension(filepath)
        filepath = filepath[:-2] + "_" + str(n) + ext
        n += 1
    return filepath


def parse_layout_csv(layout_csv):
    """ Interpret the layout csv file """

    header = open(layout_csv, 'r').readline()
    header = [h.strip() for h in header.split()]
    dt = np.format_parser(['a10', 'i4', 'i4', 'f8', 'f8', 'f8'],
                          ['name', 'number', 'beamid', 'e', 'n', 'u'], header)

    return np.genfromtxt(layout_csv, autostrip=True, skip_header=1,
                         dtype=dt.dtype)


def read_votable_catalog(gleam_votable):
    """
    Creates a list of pyuvsim source objects from the GLEAM votable catalog.
    Despite the semi-standard votable format, there are enough differences that every catalog probably
    needs its own function.
    List of tested catalogs: GLEAM EGC catalog, version 2
    """
    class Found(Exception):
        pass

    resources = parse(gleam_votable).resources
    try:
        for rs in resources:
            for tab in rs.tables:
                if 'GLEAM' in tab.array.dtype.names:
                    raise Found
    except Found:
        table = tab

    data = table.array

    sourcelist = []
    for entry in data:
        source = Source(entry['GLEAM'], Angle(entry['RAJ2000'], unit=units.deg),
                        Angle(entry['DEJ2000'], unit=units.deg),
                        freq=(200e6 * units.Hz),
                        stokes=np.array([entry['Fintwide'], 0., 0., 0.]))
        sourcelist.append(source)
    return sourcelist


def read_text_catalog(catalog_csv):
    """
        Read in a text file of sources.
        Columns:
            Source_ID = source id
            ra_j2000  = right ascension at J2000 epoch, in decimal degrees
            dec_j2000 = declination at J2000 epoch, in decimal degrees
            flux_density_I = Stokes I flux density in Janskies
            frequency = reference frequency (for future spectral indexing) [Hz]
        For now, flat spectrum sources.
    """
    header = open(catalog_csv, 'r').readline()
    header = [h.strip() for h in header.split() if not h[0] == '[']  # Ignore units in header
    dt = np.format_parser(['a10', 'f8', 'f8', 'f8', 'f8'],
                          ['source_id', 'ra_j2000', 'dec_j2000', 'flux_density_I', 'frequency'], header)

    catalog_table = np.genfromtxt(catalog_csv, autostrip=True, skip_header=1,
                                  dtype=dt.dtype)

    catalog = []

    for si in xrange(catalog_table.size):
        catalog.append(Source(catalog_table['source_id'][si],
                              Angle(catalog_table['ra_j2000'][si], unit=units.deg),
                              Angle(catalog_table['dec_j2000'][si], unit=units.deg),
                              catalog_table['frequency'][si] * units.Hz,
                              [catalog_table['flux_density_I'][si], 0, 0, 0]))

    return catalog


def create_mock_catalog(time, arrangement='zenith', array_location=None, Nsrcs=None,
                        zen_ang=None, save=False, max_za=-1.0):
    """
        Create a mock catalog with test sources at zenith.

        arrangement = Choose test point source pattern (default = 1 source at zenith)

        Keywords:
            * Nsrcs = Number of sources to put at zenith
            * array_location = EarthLocation object.
            * zen_ang = For off-zenith and triangle arrangements, how far from zenith to place sources. (deg)
            * save = Save mock catalog as npz file.

        Accepted arrangements:
            * 'triangle' = Three point sources forming a triangle around the zenith
            * 'cross'    = An asymmetric cross
            * 'horizon'  = A single source on the horizon   ## TODO
            * 'zenith'   = Some number of sources placed at the zenith.
            * 'off-zenith' = A single source off zenith
            * 'long-line' = Horizon to horizon line of point sources
            * 'hera_text' = Spell out HERA around the zenith

    """

    if not isinstance(time, Time):
        time = Time(time, scale='utc', format='jd')

    if array_location is None:
        array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s',
                                       height=1073.)
    freq = (150e6 * units.Hz)

    if arrangement not in ['off-zenith', 'zenith', 'cross', 'triangle', 'long-line', 'hera_text']:
        raise KeyError("Invalid mock catalog arrangement: " + str(arrangement))

    mock_keywords = {'time': time.jd, 'arrangement': arrangement,
                     'array_location': repr((array_location.lat.deg, array_location.lon.deg, array_location.height.value))}

    if arrangement == 'off-zenith':
        if zen_ang is None:
            zen_ang = 5.0  # Degrees
        mock_keywords['zen_ang'] = zen_ang
        Nsrcs = 1
        alts = [90. - zen_ang]
        azs = [90.]   # 0 = North pole, 90. = East pole
        fluxes = [1.0]

    if arrangement == 'triangle':
        Nsrcs = 3
        if zen_ang is None:
            zen_ang = 3.0
        mock_keywords['zen_ang'] = zen_ang
        alts = [90. - zen_ang, 90. - zen_ang, 90. - zen_ang]
        azs = [0., 120., 240.]
        fluxes = [1.0, 1.0, 1.0]

    if arrangement == 'cross':
        Nsrcs = 4
        alts = [88., 90., 86., 82.]
        azs = [270., 0., 90., 135.]
        fluxes = [5., 4., 1.0, 2.0]

    if arrangement == 'zenith':
        if Nsrcs is None:
            Nsrcs = 1
        mock_keywords['Nsrcs'] = Nsrcs
        alts = np.ones(Nsrcs) * 90.
        azs = np.zeros(Nsrcs, dtype=float)
        fluxes = np.ones(Nsrcs) * 1 / Nsrcs
        # Divide total Stokes I intensity among all sources
        # Test file has Stokes I = 1 Jy

    if arrangement == 'long-line':
        if Nsrcs is None:
            Nsrcs = 10
        if max_za < 0:
            max_za = 85
        mock_keywords['Nsrcs'] = Nsrcs
        mock_keywords['zen_ang'] = zen_ang
        fluxes = np.ones(Nsrcs, dtype=float)
        zas = np.linspace(-max_za, max_za, Nsrcs)
        alts = 90. - zas
        azs = np.zeros(Nsrcs, dtype=float)
        inds = np.where(alts > 90.0)
        azs[inds] = 180.
        alts[inds] = 90. + zas[inds]

    if arrangement == 'hera_text':

        azs = np.array([-254.055, -248.199, -236.310, -225.000, -206.565,
                        -153.435, -123.690, -111.801, -105.945, -261.870,
                        -258.690, -251.565, -135.000, -116.565, -101.310,
                        -98.130, 90.000, 90.000, 90.000, 90.000, 90.000,
                        -90.000, -90.000, -90.000, -90.000, -90.000,
                        -90.000, 81.870, 78.690, 71.565, -45.000, -71.565,
                        -78.690, -81.870, 74.055, 68.199, 56.310, 45.000,
                        26.565, -26.565, -45.000, -56.310, -71.565])

        zas = np.array([7.280, 5.385, 3.606, 2.828, 2.236, 2.236, 3.606,
                        5.385, 7.280, 7.071, 5.099, 3.162, 1.414, 2.236,
                        5.099, 7.071, 7.000, 6.000, 5.000, 3.000, 2.000,
                        1.000, 2.000, 3.000, 5.000, 6.000, 7.000, 7.071,
                        5.099, 3.162, 1.414, 3.162, 5.099, 7.071, 7.280,
                        5.385, 3.606, 2.828, 2.236, 2.236, 2.828, 3.606, 6.325])

        alts = 90. - zas
        Nsrcs = zas.size
        fluxes = np.ones_like(azs)

    catalog = []

    source_coord = SkyCoord(alt=Angle(alts, unit=units.deg), az=Angle(azs, unit=units.deg),
                            obstime=time, frame='altaz', location=array_location)
    icrs_coord = source_coord.transform_to('icrs')

    ra = icrs_coord.ra
    dec = icrs_coord.dec
    for si in range(Nsrcs):
        catalog.append(Source('src' + str(si), ra[si], dec[si], freq, [fluxes[si], 0, 0, 0]))
    if get_rank() == 0 and save:
        np.savez('mock_catalog_' + arrangement, ra=ra.rad, dec=dec.rad, alts=alts, azs=azs, fluxes=fluxes)

    catalog = np.array(catalog)
    return catalog, mock_keywords


@profile
def initialize_uvdata_from_params(obs_params):
    """
        Construct a uvdata object from parameters in a valid yaml file.

        Arguments:
            obs_params: Either an obs_param file name or a dictionary of parameters read in.
                Any uvdata parameters may be passed in through here.
        Returns:
            uv_obj, beam_list, beam_dict, beam_ids
    """

    if isinstance(obs_params, str):
        with open(obs_params, 'r') as pfile:
            param_dict = yaml.safe_load(pfile)

        param_dict['config_path'] = obs_params
    else:
        param_dict = obs_params
    # Parse telescope parameters
    tele_params = param_dict['telescope']

    telescope_config_name = tele_params['telescope_config_name']
    layout_csv = tele_params['array_layout']
    if not os.path.isdir(param_dict['config_path']):
        param_dict['config_path'] = os.path.dirname(param_dict['config_path'])
    if not os.path.exists(telescope_config_name):
        telescope_config_name = os.path.join(param_dict['config_path'], telescope_config_name)
    if not os.path.exists(layout_csv):
        layout_csv = os.path.join(param_dict['config_path'], layout_csv)

    extra_keywords = {'obs_param_file': os.path.basename(param_dict['config_path']),
                      'telescope_config_file': tele_params['telescope_config_name'],
                      'antenna_location_file': tele_params['array_layout']}

    param_dict['extra_keywords'] = extra_keywords

    ant_layout = parse_layout_csv(layout_csv)

    with open(telescope_config_name, 'r') as yf:
        telparam = yaml.safe_load(yf)
        tloc = telparam['telescope_location'][1:-1]  # drop parens
        tloc = map(float, tloc.split(","))
        tloc[0] *= np.pi / 180.
        tloc[1] *= np.pi / 180.   # Convert to radians
        telparam['telescope_location'] = uvutils.XYZ_from_LatLonAlt(*tloc)

    param_dict.update(telparam)
    E, N, U = ant_layout['e'], ant_layout['n'], ant_layout['u']
    antnames = ant_layout['name']
    beam_ids = ant_layout['beamid']
    beam_list = []
    beam_dict = {}
    for beamID in np.unique(beam_ids):
        beam_model = telparam['beam_paths'][beamID]
        which_ants = antnames[np.where(beam_ids == beamID)]
        for a in which_ants:
            beam_dict[a] = beamID
        uvb = UVBeam()
        if beam_model in ['gaussian', 'uniform', 'airy']:
            # Identify analytic beams
            if beam_model == 'gaussian':
                if 'sigma' not in telparam:
                    raise KeyError("Missing sigma for gaussian beam.")
                beam = AnalyticBeam('gaussian', sigma=telparam['sigma'])
            elif beam_model == 'airy':
                if 'diameter' not in telparam:
                    raise KeyError("Missing diameter for airy beam")
                beam = AnalyticBeam('airy', diameter=telparam['diameter'])
            else:
                beam = AnalyticBeam('uniform')
            beam_list.append(beam)
            continue
        if not os.path.exists(beam_model):
            filename = beam_model
            path = os.path.join(SIM_DATA_PATH, filename)
            if not os.path.exists(path):
                raise OSError("Could not find beam file " + filename)
        else:
            path = beam_model   # beam_model = path to beamfits
        uvb.read_beamfits(path)
        beam_list.append(uvb)

    param_dict['Nants_data'] = antnames.size
    param_dict['Nants_telescope'] = antnames.size
    param_dict['antenna_names'] = np.array(antnames.tolist())
    param_dict['antenna_numbers'] = np.array(ant_layout['number'])
    antpos_enu = np.vstack((E, N, U)).T

    param_dict['antenna_positions'] = uvutils.ECEF_from_ENU(antpos_enu, *tloc) - param_dict['telescope_location']

    # Parse frequency structure
    freq_params = param_dict['freq']

    freq_keywords = ['freq_array', 'start_freq', 'end_freq', 'Nfreqs',
                     'channel_width', 'bandwidth']
    fa, sf, ef, nf, cw, bw = [fk in freq_params for fk in freq_keywords]
    kws_used = ", ".join(freq_params.keys())

    if fa:
        freq_arr = np.array(freq_params['freq_array'])
        freq_params['Nfreqs'] = freq_arr.size
        if freq_params['Nfreqs'] > 1:
            freq_params['channel_width'] = np.diff(freq_arr)[0]
        elif 'channel_width' not in freq_params:
            raise ValueError("Channel width must be specified "
                             "if freq_arr has length 1")
    else:
        if not nf:
            if not cw:
                raise ValueError("Either channel_width or Nfreqs "
                                 " must be included in parameters:" + kws_used)
            if sf and ef:
                freq_params['bandwidth'] = freq_params['end_freq'] - freq_params['start_freq']
                bw = True
            if bw:
                freq_params['Nfreqs'] = int(np.floor(freq_params['bandwidth']
                                                     / freq_params['channel_width']))
            else:
                raise ValueError("Either bandwidth or band edges "
                                 "must be specified: " + kws_used)

        if not cw:
            if not bw:
                raise ValueError("Either bandwidth or channel width"
                                 " must be specified: " + kws_used)
            freq_params['channel_width'] = (freq_params['bandwidth']
                                            / float(freq_params['Nfreqs']))

        if not bw:
            freq_params['bandwidth'] = (freq_params['channel_width']
                                        * freq_params['Nfreqs'])
            bw = True
        if not sf:
            if ef and bw:
                freq_params['start_freq'] = freq_params['end_freq'] - freq_params['bandwidth']
        if not ef:
            if sf and bw:
                freq_params['end_freq'] = freq_params['start_freq'] + freq_params['bandwidth']

        freq_arr = np.linspace(freq_params['start_freq'],
                               freq_params['end_freq'],
                               freq_params['Nfreqs'], endpoint=False)
    if freq_params['Nfreqs'] != 1:
        assert np.allclose(np.diff(freq_arr), freq_params['channel_width'] * np.ones(freq_params["Nfreqs"] - 1), atol=1.0)  # 1 Hz

    Nspws = 1 if 'Nspws' not in freq_params else freq_params['Nspws']
    freq_arr = np.repeat(freq_arr, Nspws).reshape(Nspws, freq_params['Nfreqs'])

    param_dict['channel_width'] = freq_params['channel_width']
    param_dict['freq_array'] = freq_arr
    param_dict['Nfreqs'] = freq_params['Nfreqs']

    # Parse time structure
    time_params = param_dict['time']

    optional_time_params = ['time_format', 'snapshot_length_hours']

    for k in optional_time_params:
        if k in time_params:
            param_dict[k] = time_params[k]

    time_keywords = ['start_time', 'end_time', 'Ntimes', 'integration_time',
                     'duration_hours', 'duration_days']
    st, et, nt, it, dh, dd = [tk in time_params for tk in time_keywords]
    kws_used = ", ".join(time_params.keys())
    daysperhour = 1 / 24.
    hourspersec = 1 / 60.**2
    dayspersec = daysperhour * hourspersec

    if dh and not dd:
        time_params['duration'] = time_params['duration_hours'] * daysperhour
        dd = True
    elif dd:
        time_params['duration'] = time_params['duration_days']

    if not nt:
        if not it:
            raise ValueError("Either integration_time or Ntimes must be "
                             "included in parameters:" + kws_used)
        if st and et:
            time_params['duration'] = time_params['end_time'] - time_params['start_time']
            dd = True
        if dd:
            time_params['Ntimes'] = int(np.round(time_params['duration']
                                                 / (time_params['integration_time']
                                                    * dayspersec))) + 1
        else:
            raise ValueError("Either duration or time bounds must be specified: "
                             + kws_used)

    if not it:
        if not dd:
            raise ValueError("Either duration or integration time "
                             "must be specified: " + kws_used)
        time_params['integration_time'] = (time_params['duration']
                                           / float(time_params['Ntimes']))
        if 'snapshot_length_hours' in time_params:
            print('Warning: Setting integration time from Ntimes for snapshots. '
                  'This may not be well-defined.')

    inttime_days = time_params['integration_time'] * 1 / (24. * 3600.)
    inttime_days = np.trunc(inttime_days * 24 * 3600) / (24. * 3600.)
    if not dd:
        time_params['duration'] = inttime_days * (time_params['Ntimes'])
        dd = True
    if not st:
        if et and dd:
            time_params['start_time'] = time_params['end_time'] - time_params['duration']
    if not et:
        if st and dd:
            time_params['end_time'] = time_params['start_time'] + time_params['duration']
    if not (st or et):
        raise ValueError("Either a start or end time must be specified: " + kws_used)

    time_arr = np.linspace(time_params['start_time'],
                           time_params['end_time'],
                           time_params['Ntimes'], endpoint=False)

    if time_params['Ntimes'] != 1:
        assert np.allclose(np.diff(time_arr), inttime_days * np.ones(time_params["Ntimes"] - 1), atol=1e-4)   # To nearest second

    Nbl = (param_dict['Nants_data'] + 1) * param_dict['Nants_data'] / 2

    time_arr = np.sort(np.tile(time_arr, Nbl))
    param_dict['integration_time'] = (np.ones_like(time_arr, dtype=np.float64)
                                      * time_params['integration_time'])
    param_dict['time_array'] = time_arr
    param_dict['Ntimes'] = time_params['Ntimes']
    param_dict['Nspws'] = 1
    param_dict['Npols'] = 4
    # Now make a UVData object with these settings built in.
    # The syntax below allows for other valid uvdata keywords to be passed
    #  without explicitly setting them here.

    uv_obj = UVData()
    for k in param_dict:
        if hasattr(uv_obj, k):
            setattr(uv_obj, k, param_dict[k])

    bls = np.array([uv_obj.antnums_to_baseline(uv_obj.antenna_numbers[j], uv_obj.antenna_numbers[i])
                    for i in range(0, uv_obj.Nants_data)
                    for j in range(i, uv_obj.Nants_data)])

    uv_obj.baseline_array = np.tile(bls, uv_obj.Ntimes)
    uv_obj.Nbls = bls.size
    uv_obj.Nblts = uv_obj.Nbls * uv_obj.Ntimes

    uv_obj.ant_1_array, uv_obj.ant_2_array = \
        uv_obj.baseline_to_antnums(uv_obj.baseline_array)

    return uv_obj, beam_list, beam_dict, beam_ids


@profile
def uvdata_to_telescope_config(uvdata_in, beam_filepath, layout_csv_name=None,
                               telescope_config_name=None,
                               return_names=False, path_out='.'):
    """
        From a uvfits file, generate telescope parameters files.
        Output config files are written to the current directory, unless keep_path is set.
        Arguments:
            uvdata_in (UVData): object to process
            path_out (str): Target directory for the config file.
            beam_filepath (str): Path to a beamfits file.
        Keywords:
            layout_csv_name (str, optional): The name for the antenna positions
                csv file (Default <telescope_name>_layout.csv)
            telescope_config_name: The name for the telescope config file
                (Default teleconfig_#number.yaml)
            return_names: Return the file names for loopback tests.
        Returns:
            if return_names, returns tuple (path, telescope_config_name, layout_csv_name)

    """

    if telescope_config_name is None:
        telescope_config_path = \
            check_file_exists_and_increment(os.path.join(path_out, 'telescope_config_'
                                                         + uvdata_in.telescope_name
                                                         + '.yaml'))
        telescope_config_name = os.path.basename(telescope_config_path)

    if layout_csv_name is None:
        layout_csv_path = check_file_exists_and_increment(os.path.join(path_out, uvdata_in.telescope_name + "_layout.csv"))
        layout_csv_name = os.path.basename(layout_csv_path)

    antpos_enu, antenna_numbers = uvdata_in.get_ENU_antpos()

    e, n, u = antpos_enu.T
    beam_ids = np.zeros_like(e).astype(int)
    col_width = max([len(name) for name in uvdata_in.antenna_names])
    header = ("{:" + str(col_width) + "} {:8} {:8} {:10} {:10} {:10}\n").format("Name", "Number", "BeamID", "E", "N", "U")

    with open(os.path.join(path_out, layout_csv_name), 'w') as lfile:
        lfile.write(header + '\n')
        for i in range(beam_ids.size):
            e, n, u = antpos_enu[i]
            beam_id = beam_ids[i]
            name = uvdata_in.antenna_names[i]
            num = uvdata_in.antenna_numbers[i]
            line = ("{:" + str(col_width) + "} {:8d} {:8d} {:10.4f} {:10.4f} {:10.4f}\n").format(name, num, beam_id, e, n, u)
            lfile.write(line)

    # Write the rest to a yaml file.
    yaml_dict = dict(
        telescope_name=uvdata_in.telescope_name,
        telescope_location=repr(uvdata_in.telescope_location_lat_lon_alt_degrees),
        Nants=uvdata_in.Nants_telescope,
        beam_paths={
            0: beam_filepath
        }
    )

    with open(os.path.join(path_out, telescope_config_name), 'w+') as yfile:
        yaml.dump(yaml_dict, yfile, default_flow_style=False)

    print('Path: {}, telescope_config: {}, layout: {}'.format(path_out, telescope_config_name,
                                                              layout_csv_name))

    if return_names:
        return path_out, telescope_config_name, layout_csv_name


@profile
def uvdata_to_config_file(uvdata_in, param_filename=None, telescope_config_name='',
                          layout_csv_name='', catalog='mock', path_out='.'):
    """
    Extract simulation configuration settings from uvfits.
    Arguments:
        uvdata_in (UVData) : uvdata object.

    Keywords:
        param_filename (str, optional) : output param file name, defaults to obsparam_#.yaml.
        telescope_config_name (str, optional) : Name of yaml file file. Defaults to blank string.
        layout_csv_name (str, optional) : Name of layout csv file. Defaults to blank string.
        catalog (str, optional): Path to catalog file, defaults to 'mock'.
        path_out (str, optional) : Where to put config files.

    Returns:
        Nothing

    """

    if param_filename is None:
        param_filename = check_file_exists_and_increment(os.path.join(path_out, 'obsparam.yaml'))

    freq_array = uvdata_in.freq_array[0, :].tolist()
    time_array = uvdata_in.time_array.tolist()
    integration_time_array = np.array(uvdata_in.integration_time)
    if np.max(integration_time_array) != np.min(integration_time_array):
        warnings.warn('The integration time is not constant. Using the shortest integration time')
    integration_time = float(np.min(integration_time_array))
    param_dict = dict(
        time=dict(
            start_time=time_array[0],
            end_time=time_array[-1],
            Ntimes=uvdata_in.Ntimes,
            integration_time=integration_time,
        ),
        freq=dict(
            start_freq=freq_array[0],
            end_freq=freq_array[-1] + uvdata_in.channel_width,
            channel_width=uvdata_in.channel_width,
            Nfreqs=uvdata_in.Nfreqs,
        ),
        sources=dict(
            catalog=catalog
        ),

        telescope=dict(
            telescope_config_name=telescope_config_name,
            array_layout=layout_csv_name
        ),

        filing=dict(
            outdir='.',
            outfile_name='',
            outfile_prefix=''
        )
    )

    if catalog == 'mock':
        param_dict['sources']['mock_arrangement'] = 'zenith'

    with open(os.path.join(path_out, param_filename), 'w') as yfile:
        yaml.dump(param_dict, yfile, default_flow_style=False)


def write_uvfits(uv_obj, param_dict, return_filename=False, dryrun=False):
    """
        Parse output file information from parameters and write uvfits to file.
    """
    param_dict = param_dict['filing']
    if 'outdir' not in param_dict:
        param_dict['outdir'] = '.'
    if 'outfile_name' not in param_dict or param_dict['outfile_name'] == '':
        outfile_prefix = ""
        outfile_suffix = "results"
        if 'outfile_prefix' in param_dict:
            outfile_prefix = param_dict['outfile_prefix']
        if 'outfile_suffix' in param_dict:
            outfile_suffix = param_dict['outfile_suffix']
        outfile_name = "_".join([outfile_prefix, outfile_suffix])
        outfile_name = os.path.join(param_dict['outdir'], outfile_name)
    else:
        outfile_name = os.path.join(param_dict['outdir'], param_dict['outfile_name'])
    print('Outfile path: ', outfile_name)
    if not outfile_name.endswith(".uvfits"):
        outfile_name = outfile_name + ".uvfits"

    if 'clobber' not in param_dict:
        outfile_name = check_file_exists_and_increment(outfile_name)

    if not os.path.exists(param_dict['outdir']):
        os.makedirs(param_dict['outdir'])

    if not dryrun:
        uv_obj.write_uvfits(outfile_name, force_phase=True, spoof_nonessential=True)
    if return_filename:
        return outfile_name
