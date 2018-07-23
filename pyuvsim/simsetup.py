# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 2-clause BSD License

from pyuvdata import UVBeam, UVData
import numpy as np
import yaml
import os
import pyuvdata.utils as uvutils
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
import pyuvsim
import astropy.units as units
from astropy.coordinates import Angle
# Utilities for setting up simulations for parameter files,
# and for generating parameter files from uvfits files


def strip_extension(filepath):
    if '.' not in filepath:
        return filepath, ''
    file_list = filepath.split('.')
    return ".".join(file_list[:-1]), '.' + file_list[-1]


def check_file_exists_and_increment(filepath):
    """
        Given filepath (path + filename), check if it exists. If so, add a _1
        at the end. etc.
    """
    if os.path.exists(filepath):
        filepath, ext = strip_extension(filepath)
        filepath += "_0" + ext
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


def point_sources_from_params(catalog_csv):
    """
        Read in a text file of sources.
        Columns:
            Source_ID = source id
            ra_j2000  = right ascension at J2000 epoch, in decimal hours
            dec_j2000 = declination at J2000 epoch, in decimal degrees
            flux_density_I = Stokes I flux density in Janskies
            frequency = reference frequency (for future spectral indexing) [Hz]
        For now, flat spectrum sources.
    """
    header = open(catalog_csv, 'r').readline()
    header = [h.strip() for h in header.split()]
    dt = np.format_parser(['a10', 'f8', 'f8', 'f8', 'f8'],
                          ['source_id', 'ra_j2000', 'dec_j2000', 'flux_density_I', 'frequency'], header)

    catalog_table = np.genfromtxt(catalog_csv, autostrip=True, skip_header=1,
                                  dtype=dt.dtype)

    catalog = []

    for si in xrange(catalog_table.size):
        catalog.append(pyuvsim.Source(catalog_table['source_id'][si],
                                      Angle(catalog_table['ra_j2000'][si], unit=units.hour),
                                      Angle(catalog_table['dec_j2000'][si], unit=units.deg),
                                      catalog_table['frequency'][si] * units.Hz,
                                      [catalog_table['flux_density_I'][si], 0, 0, 0]))

    return catalog


def initialize_uvdata_from_params(obs_params):
    """
        Construct a uvdata object from parameters in a valid yaml file.

        Arguments:
            obs_params: Either an obs_param file name or a dictionary of parameters read in.
                Any uvdata parameters may be passed in through here.
        Returns:
            uv_obj, beam_list, beam_ids
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

    if not os.path.exists(telescope_config_name):
        path = os.path.dirname(param_dict['config_path'])
        telescope_config_name = os.path.join(path, telescope_config_name)

    if not os.path.exists(layout_csv):
        path = os.path.dirname(param_dict['config_path'])
        layout_csv = os.path.join(path, layout_csv)

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
        if beam_model in ['gaussian', 'tophat']:
            # Identify analytic beams
            if beam_model == 'gaussian':
                try:
                    beam = pyuvsim.AnalyticBeam('gaussian', sigma=telparam['sigma'])
                except KeyError as err:
                    print("Missing sigma for gaussian beam.")
                    raise err
            else:
                beam = pyuvsim.AnalyticBeam('tophat')
            beam_list.append(beam)
            continue
        if not os.path.exists(beam_model):
            filename = beam_model
            path = os.path.join(SIM_DATA_PATH, filename)
            if not os.path.exists(path):
                raise OSError("Could not find file " + filename)
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
        try:
            if freq_params['Nfreqs'] > 1:
                freq_params['channel_width'] = np.diff(freq_arr)[0]
            else:
                freq_params['channel_width'] = freq_params['channel_width']
        except KeyError:
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
        try:
            assert np.allclose(np.diff(freq_arr), freq_params['channel_width'] * np.ones(freq_params["Nfreqs"] - 1), atol=1.0)  # 1 Hz
        except AssertionError as err:
            print freq_params
            print freq_params['Nfreqs']
            print freq_params['bandwidth'] / 2.
            print np.diff(freq_arr)[0]
            raise err

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
        raise ValueError("Either a start or end time must be specified" + kws_used)

    time_arr = np.linspace(time_params['start_time'],
                           time_params['end_time'],
                           time_params['Ntimes'], endpoint=False)

    if time_params['Ntimes'] != 1:
        try:
            assert np.allclose(np.diff(time_arr), inttime_days * np.ones(time_params["Ntimes"] - 1), atol=1e-4)   # To nearest second
        except AssertionError as err:
            print time_params
            print np.diff(time_arr)[0]
            raise err

    Nbl = (param_dict['Nants_data'] + 1) * param_dict['Nants_data'] / 2

    time_arr = np.sort(np.tile(time_arr, Nbl))
    param_dict['integration_time'] = time_params['integration_time']
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
    inttime_days = time_array[1] - time_array[0]
    param_dict = dict(
        time=dict(
            start_time=time_array[0],
            end_time=time_array[-1],
            Ntimes=uvdata_in.Ntimes,
            integration_time=uvdata_in.integration_time,
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
