from pyuvdata import UVBeam, UVData
import numpy as np
import yaml
import os
import pyuvdata.utils as uvutils

# Utilities for setting up simulations for parameter files,
# and for generating parameter files from uvfits files


def check_file_exists_and_increment(filepath):
    """
        Given filepath (filename + path), check if it exists. If so, add a _1
        at the end. etc.
    """
    if os.path.exists(filepath):
        filepath += "_0"
    else:
        return filepath
    n = 1
    while os.path.exists(filepath):
        filepath = filepath[:-2] + "_" + str(n)
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


def initialize_uvdata_from_params(param_dict):
    """
        Construct a uvdata object from parameters in a valid yaml file.

        Arguments:
            param : dict
                If dict, this is treated as the yaml file that has been read
        Returns:
            uv_obj, beam_list, beam_ids

        TODO:
            Verify the parameter file before reading.
            Supply defaults when possible.
    """

    # Parse telescope parameters
    tele_params = param_dict['telescope']

    teleyaml = tele_params['teleyaml']
    layout_csv = tele_params['array_layout']

    if not os.path.exists(teleyaml):
        path = os.path.dirname(param_dict['config_path'])
        teleyaml = os.path.join(path, teleyaml)

    if not os.path.exists(layout_csv):
        path = os.path.dirname(param_dict['config_path'])
        layout_csv = os.path.join(path, layout_csv)

    ant_layout = parse_layout_csv(layout_csv)

    with open(teleyaml, 'r') as yf:
            telparam = yaml.safe_load(yf)
            tloc = telparam['telescope_location'][1:-1]  # drop parens
            tloc = map(float, tloc.split(","))
            tloc[0] *= np.pi/180.
            tloc[1] *= np.pi/180.   # Convert to radians
            telparam['telescope_location'] = uvutils.XYZ_from_LatLonAlt(*tloc)

    param_dict.update(telparam)
    E, N, U = ant_layout['e'], ant_layout['n'], ant_layout['u']
    antnames = ant_layout['name']
    beam_ids = ant_layout['beamid']
    beam_list = []
    for beamID in np.unique(beam_ids):
        path = telparam['beam_paths'][beamID]
        uvb = UVBeam()
        uvb.read_beamfits(path)
        beam_list.append(uvb)
    param_dict['Nants_data'] = antnames.size
    param_dict['Nants_telescope'] = antnames.size
    param_dict['antenna_names'] = np.array(antnames.tolist())
    param_dict['antenna_numbers'] = np.array(ant_layout['number'])
    antpos_enu = np.vstack((E, N, U)).T

    param_dict['antenna_positions'] = uvutils.ECEF_from_ENU(antpos_enu.T, *tloc).T - param_dict['telescope_location']
    del param_dict['telescope']

    # Parse frequency structure
    freq_params = param_dict['freq']

    freq_keywords = ['freq_array', 'start_freq', 'end_freq', 'Nfreqs',
                     'channel_width', 'bandwidth', 'center_freq']
    fa, sf, ef, nf, cw, bw, cf = [fk in freq_params for fk in freq_keywords]
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
            raise ValueError("Channel width must be specified " +
                             "if freq_arr has length 1")
    else:
        if not nf:
            if not cw:
                raise ValueError("Either channel_width or Nfreqs " +
                                 " must be included in parameters:" + kws_used)
            if sf and ef:
                freq_params['bandwidth'] = freq_params['start_freq'] - freq_params['end_freq']
                bw = True
            if bw:
                freq_params['Nfreqs'] = int(np.floor(freq_params['bandwidth'] /
                                            freq_params['channel_width'])) + 1
            else:
                raise ValueError("Either bandwidth or band edges " +
                                 "must be specified: " + kws_used)

        if not cw:
            if not bw:
                raise ValueError("Either bandwidth or channel width" +
                                 " must be specified: " + kws_used)
            freq_params['channel_width'] = (freq_params['bandwidth'] /
                                            float(freq_params['Nfreqs']))
    
        if not bw:
            freq_params['bandwidth'] = (freq_params['channel_width'] *
                                        freq_params['Nfreqs'])
        if not sf:
            if ef and bw:
                freq_params['start_freq'] = freq_params['end_freq'] - freq_params['bandwidth']
            if cf and bw:
                freq_params['start_freq'] = freq_params['center_freq'] - freq_params['bandwidth']/2.
        if not ef:
            if sf and bw:
                freq_params['end_freq'] = freq_params['start_freq'] + freq_params['bandwidth']
            if cf and bw:
                freq_params['end_freq'] = freq_params['center_freq'] + freq_params['bandwidth']/2.

        freq_arr = np.linspace(freq_params['start_freq'], freq_params['end_freq'], freq_params['Nfreqs'])

    Nspws = 1 if 'Nspws' not in freq_params else freq_params['Nspws']
    freq_arr = np.repeat(freq_arr, Nspws).reshape(Nspws, freq_params['Nfreqs'])

    param_dict['channel_width'] = freq_params['channel_width']
    param_dict['freq_array'] = freq_arr
    param_dict['Nfreqs'] = freq_params['Nfreqs']
    del param_dict['freq']

    # Parse time structure
    time_params = param_dict['time']

    time_keywords = ['start_time', 'end_time', 'Ntimes', 'integration_time',
                     'duration_hours', 'duration_days']
    st, et, nt, it, dh, dd = [tk in time_params for tk in time_keywords]
    kws_used = ", ".join(time_params.keys())
    daysperhour = 1/24.
    hourspersec = 1/60.**2
    dayspersec = daysperhour*hourspersec

    if dh and not dd:
        time_params['duration'] = time_params['duration_hours'] * daysperhour
        dd = True
    elif dd:
        time_params['duration'] = time_params['duration_days']

    if not nt:
        if not it:
            raise ValueError("Either integration_time or Ntimes " +
                             "must be included in parameters:" + kws_used)
        if st and et:
            time_params['duration'] = time_params['start_time'] - time_params['end_time']
            dd = True
        if dd:
            time_params['Ntimes'] = int(np.floor(time_params['duration'] /
                                        (time_params['integration_time']
                                        * dayspersec))) + 1
        else:
            raise ValueError("Either duration or time bounds " +
                             "must be specified: " + kws_used)

    if not it:
        if not dd:
            raise ValueError("Either duration or integration time " +
                             "must be specified: " + kws_used)
        time_params['integration_time'] = (time_params['duration'] /
                                           float(time_params['Ntimes']))

    if not dd:
        time_params['duration'] = time_params['integration_time'] * time_params['Ntimes']
    if not st:
        if et and dd:
            time_params['start_time'] = time_params['end_time'] - time_params['duration']
    if not et:
        if st and dd:
            time_params['end_time'] = time_params['start_time'] + time_params['duration']
    if not (st or et):
        raise ValueError("Either a start or end time must be specified" + kws_used)

    time_arr = np.linspace(time_params['start_time'],
                           time_params['end_time'] +
                           time_params['integration_time']*dayspersec,
                           time_params['Ntimes'])

    Nbl = (param_dict['Nants_data']+1)*param_dict['Nants_data'] / 2

    time_arr = np.sort(np.tile(time_arr, Nbl))
    param_dict['integration_time'] = time_params['integration_time']
    param_dict['time_array'] = time_arr
    param_dict['Ntimes'] = time_params['Ntimes']
    del param_dict['time']

    # Now make a UVData object with these settings built in.
    # The syntax below allows for other valid uvdata keywords to be passed
    #  without explicitly setting them here.

    uv_obj = UVData()
    print param_dict
    for k in param_dict:
        print k, param_dict[k]
        if hasattr(uv_obj, k):
            setattr(uv_obj, k, param_dict[k])

    bls = np.array([uv_obj.antnums_to_baseline(j, i)
                   for i in range(0, uv_obj.Nants_data)
                   for j in range(i, uv_obj.Nants_data)])

    uv_obj.baseline_array = np.tile(bls, uv_obj.Ntimes)
    uv_obj.Nbls = bls.size
    uv_obj.Nblts = uv_obj.Nbls * uv_obj.Ntimes

    uv_obj.ant_1_array, uv_obj.ant_2_array = \
        uv_obj.baseline_to_antnums(uv_obj.baseline_array)

    return uv_obj, beam_list, beam_ids


def uvfits_to_teleyaml(file_in, beam_filepath, layout_csv_fname=None,
                       teleyaml_fname=None, keep_path=False,
                       return_names=False, path='.'):
    """
        From a uvfits file, generate telescope parameters files.
        Output config files are written to the current directory,
            unless keep_path is set.
        Arguments:
            file_in = A uvfits filename (required).
            beam_filepath = Path to a beamfits file (required).
            layout_fname = The name for the antenna positions csv file (defaults to telescope_name_layout.csv)
            teleyaml_fname = The name for the telescope yaml file (defaults to file_in.yaml)
            keep_path = Put config files in the same directory as the input uvfits file.
            return_names = Return the file names for loopback tests.

    """

    uv_fname = os.path.basename(file_in)
    uv_path = os.path.dirname(file_in)
    if keep_path:
        path = uv_path

    if teleyaml_fname is None:
        teleyaml_fname = '.'.join(uv_fname.split('.')[:-1]) + '.yaml'

    input_uv = UVData()
    input_uv.read_uvfits(os.path.join(uv_path, uv_fname), read_data=False)

    if layout_csv_fname is None:
        layout_csv_fname = input_uv.telescope_name + "_layout.csv"

    antpos_enu = uvutils.ENU_from_ECEF((input_uv.antenna_positions +
                                        input_uv.telescope_location).T,
                                       * input_uv.telescope_location_lat_lon_alt).T

    e, n, u = antpos_enu.T

    beam_ids = np.zeros_like(e).astype(int)
    col_width = max([len(name) for name in input_uv.antenna_names])
    header = ("{:"+str(col_width)+"} {:8} {:8} {:10} {:10} {:10}\n").format("Name", "Number", "BeamID", "E", "N", "U")

    with open(os.path.join(path, layout_csv_fname), 'w') as lfile:
        lfile.write(header+'\n')
        for i in range(beam_ids.size):
            e, n, u = antpos_enu[i]
            beam_id = beam_ids[i]
            name = input_uv.antenna_names[i]
            num = input_uv.antenna_numbers[i]
            line = ("{:"+str(col_width)+"} {:8d} {:8d} {:10.4f} {:10.4f} {:10.4f}\n").format(name, num, beam_id, e, n, u)
            lfile.write(line)

    # Write the rest to a yaml file.
    yaml_dict = dict(
        telescope_name=input_uv.telescope_name,
        telescope_location=repr(input_uv.telescope_location_lat_lon_alt_degrees),
        Nants=input_uv.Nants_telescope,
        beam_paths={
                0: beam_filepath
                    }
    )

    with open(os.path.join(path, teleyaml_fname), 'w+') as yfile:
        yaml.dump(yaml_dict, yfile, default_flow_style=False)

    print('Path: {}, teleyaml: {}, layout: {}'.format(path, teleyaml_fname,
                                                      layout_csv_fname))

    if return_names:
        return path, teleyaml_fname, layout_csv_fname


def uvfits_to_yaml(file_in, yaml_fname=None, keep_path=False, teleyaml_path='',
                   layout_csv_path='', catalog='mock', path='.'):
    """
    Extract simulation configuration settings from uvfits.
    Args:
        file_in = (string) uvfits file name
        yaml_fname (optional) = output param file name, defaults to file_in.yaml.
        keep_path (optional) = If true, place config file in the same directory as the uvfits.
        teleyaml_path (optional) = Path to teleyaml file. Defaults to ''
        layout_csv_path (optional) = Path to layout csv file. Defaults to ''
        catalog = Path to catalog file, defaults to 'mock'.

    """

    uv_fname = os.path.basename(file_in)
    uv_path = os.path.dirname(file_in)
    if keep_path:
        path = uv_path

    if yaml_fname is None:
        yaml_fname = 'config_'+'.'.join(uv_fname.split('.')[:-1]) + '.yaml'

    input_uv = UVData()
    input_uv.read_uvfits(os.path.join(uv_path, uv_fname), read_data=False)

    freq_array = input_uv.freq_array[0, :].tolist()
    time_array = input_uv.time_array.tolist()

    param_dict = dict(
        time=dict(
            start_time=time_array[0],
            end_time=time_array[-1],
            Ntimes=input_uv.Ntimes,
            integration_time=input_uv.integration_time,
        ),
        freq=dict(
            start_freq=freq_array[0],
            end_freq=freq_array[-1],
            channel_width=input_uv.channel_width,
            Nfreqs=input_uv.Nfreqs,
        ),

        # Params will still require a catalog and telescope parameters:
        sources=dict(
            catalog=catalog
            ),

        telescope=dict(
            teleyaml=teleyaml_path,
            array_layout=layout_csv_path
        ),

        filing=dict(
            outdir='.',
            outfile_name='',
            outfile_prefix=''
        )
    )

    if catalog == 'mock':
        param_dict['sources']['mock_arrangement'] = 'zenith'

    with open(os.path.join(path, yaml_fname), 'w') as yfile:
        yaml.dump(param_dict, yfile, default_flow_style=False)
