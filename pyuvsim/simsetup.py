### Utilities for setting up simulations for parameter files, and for generating parameter files from uvfits files

from pyuvdata.utils import ECEF_from_ENU, XYZ_from_LatLonAlt
from pyuvdata import UVBeam, UVData
import numpy as np, yaml, os


def check_file_exists_and_increment(filepath):
    """
        Given filepath (filename + path), check if it exists. If so, add a _1 at the end. etc.
    """
    if os.path.exists(filepath):
        filepath += "_0"
    else:
        return filepath
    n=1
    while os.path.exists(filepath):
        filepath = filepath[:-2] + "_" + str(n)
        n += 1
    return filepath

def parse_layout_csv(layout_csv):
    header = open(layout_csv,'r').readline()
    header = [h.strip() for h in header.split()]
    dt = np.format_parser(['a10','i4','i4','f8','f8','f8'],['name','number','beamid','e','n','u'], header)
    return np.genfromtxt(layout_csv,autostrip=True, skip_header=1,dtype=dt.dtype)

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

    ## Parse telescope parameters
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
            telescope_location_lat_lon_alt_degrees = map(float,tloc.split(","))
            tloc = telescope_location_lat_lon_alt_degrees
            tloc[0] *= np.pi/180.
            tloc[1] *= np.pi/180.   #Convert to radians
            telparam['telescope_location'] = XYZ_from_LatLonAlt(*tloc)
   
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
    antpos_enu = np.vstack((E,N,U)).T

    param_dict['antenna_positions'] = ECEF_from_ENU((antpos_enu).T, *tloc).T - param_dict['telescope_location']
    del param_dict['telescope']


    ## Parse frequency structure
    freq_params = param_dict['freq']

    freq_keywords = ['freq_array', 'start_freq', 'end_freq', 'Nfreqs', 'channel_width', 'bandwidth', 'center_freq']
    fa, sf, ef, nf, cw, bw, cf = [ fk in freq_params for fk in freq_keywords]
    kws_used = ", ".join(freq_params.keys())

    if fa:
        freq_arr = np.array(freq_params['freq_array'])
        freq_params['Nfreqs'] = freq_arr.size
        try:
            freq_params['channel_width'] = np.diff(freq_arr)[0] if freq_params['Nfreqs'] > 1 else freq_params['channel_width']
        except KeyError:
            raise ValueError("Channel width must be specified if freq_arr has length 1")
    else:
        if not nf:
            if not cw: raise ValueError("Either channel_width or Nfreqs must be included in parameters:" + kws_used)
            if sf and ef:
                freq_params['bandwidth'] = freq['start_freq'] - freq['end_freq']
                bw=True
            if bw: freq_params['Nfreqs'] = int(np.floor(freq_params['bandwidth'] / freq_params['channel_width'])) + 1
            else: raise ValueError("Either bandwidth or band edges must be specified: " + kws_used)
    
        if not cw:
            if not bw: raise ValueError("Either bandwidth or channel width must be specified: "+ kws_used)
            freq_params['channel_width'] = freq_params['bandwidth']/float(freq_params['Nfreqs'])
    
        if not bw: freq_params['bandwidth'] = freq_params['channel_width'] * freq_params['Nfreqs']
        if not sf:
            if ef and bw: freq_params['start_freq'] = freq_params['end_freq'] - freq_params['bandwidth']
            if cf and bw: freq_params['start_freq'] = freq_params['center_freq'] - freq_params['bandwidth']/2.
        if not ef:
            if sf and bw: freq_params['end_freq'] = freq_params['start_freq'] + freq_params['bandwidth']
            if cf and bw: freq_params['end_freq'] = freq_params['center_freq'] + freq_params['bandwidth']/2.

        freq_arr = np.linspace(freq_params['start_freq'], freq_params['end_freq'], freq_params['Nfreqs'])

    Nspws = 1 if not 'Nspws' in freq_params else freq_params['Nspws']
    freq_arr = np.repeat(freq_arr,Nspws).reshape(Nspws,freq_params['Nfreqs'])

    param_dict['channel_width'] = freq_params['channel_width']
    param_dict['freq_array'] = freq_arr
    param_dict['Nfreqs'] = freq_params['Nfreqs']
    del param_dict['freq']


    ## Parse time structure
    time_params = param_dict['time']

    time_keywords = ['start_time', 'end_time', 'Ntimes', 'integration_time', 'duration_hours', 'duration_days']
    st, et, nt, it, dh, dd = [ tk in time_params for tk in time_keywords]
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
        if not it: raise ValueError("Either integration_time or Ntimes must be included in parameters:" + kws_used)
        if st and et:
            time_params['duration'] = time_params['start_time'] - time_params['end_time']
            dd=True
        if dd: time_params['Ntimes'] = int(np.floor(time_params['duration'] / (time_params['integration_time'] * dayspersec))) + 1
        else: raise ValueError("Either duration or time bounds must be specified: " + kws_used)
    
    if not it:
        if not dd: raise ValueError("Either duration or integration time must be specified: "+ kws_used)
        time_params['integration_time'] = time_params['duration'] / float(time_params['Ntimes']) 
    
    if not dd: time_params['duration'] = time_params['integration_time'] * time_params['Ntimes']
    if not st:
        if et and dd: time_params['start_time'] = time_params['end_time'] - time_params['duration']
    if not et:
        if st and dd: time_params['end_time'] = time_params['start_time'] + time_params['duration']
    if not (st or et): raise ValueError("Either a start or end time must be specified")

    time_arr = np.linspace(time_params['start_time'], time_params['end_time'] + time_params['integration_time']*dayspersec, time_params['Ntimes'])

    Nbl = (param_dict['Nants_data']+1)*param_dict['Nants_data'] / 2

    time_arr =  np.sort(np.tile(time_arr,Nbl))
    param_dict['integration_time'] = time_params['integration_time']
    param_dict['time_array'] = time_arr
    param_dict['Ntimes'] = time_params['Ntimes']
    del param_dict['time']


    ## Now make a UVData object with these settings built in.
    ## The syntax below allows for other valid uvdata keywords to be passed in without explicitly setting them here.

    uv_obj = UVData()
    print param_dict
    for k in param_dict:
        print k, param_dict[k]
        if hasattr(uv_obj,k):
            setattr(uv_obj,k,param_dict[k])

    bls = np.array([uv_obj.antnums_to_baseline(j,i)#,attempt256=True)
               for i in range(0,uv_obj.Nants_data) 
               for j in range(i,uv_obj.Nants_data) ])
    
    uv_obj.baseline_array = np.tile(bls, uv_obj.Ntimes)
    uv_obj.Nbls = bls.size
    uv_obj.Nblts = uv_obj.Nbls * uv_obj.Ntimes

    uv_obj.ant_1_array, uv_obj.ant_2_array = \
          uv_obj.baseline_to_antnums(uv_obj.baseline_array)

    return uv_obj, beam_list, beam_ids
