### Utilities for setting up simulations for parameter files, and for generating parameter files from uvfits files

from pyuvdata.utils import ECEF_from_ENU, XYZ_from_LatLonAlt
from pyuvdata import UVBeam, UVData
import numpy as np, yaml, os


def check_required_params(param_dict):
    """
        Assert that the minimum parameters are provided for simulation.
    """
    assert(1==1)

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
        Construct a uvdata object from the parameters in param_dict.
        TODO:
            Verify the parameter file before reading.
            Supply defaults when possible.
    """

    ## Parse telescope parameters
    tele_params = param_dict['telescope']

    teleyaml = tele_params['teleyaml']
    layout_csv = tele_params['array_layout']
    
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
    for beamID in beam_ids:
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
    if 'Nfreqs' in freq_params:
        freq_arr = np.linspace(freq_params['start_freq'],freq_params['end_freq'],freq_params['Nfreqs'])
    else:
        freq_arr = np.arange(freq_params['start_freq'],freq_params['end_freq'],freq_params['channel_width'])

    Nspws = 1
    freq_arr = np.repeat(freq_arr,Nspws).reshape(Nspws,freq_params['Nfreqs'])

    if freq_params['Nfreqs'] > 1:
        assert np.diff(freq_arr)[0] == freq_params['channel_width']

    param_dict['channel_width'] = freq_params['channel_width']
    param_dict['freq_array'] = freq_arr
    param_dict['Nfreqs'] = freq_params['Nfreqs']
    del param_dict['freq']


    ## Parse time structure
    time_params = param_dict['time']
    if 'Ntimes' in time_params:
        time_arr = np.linspace(time_params['start_time'],time_params['end_time'],time_params['Ntimes'])
    else:
        time_arr = np.arange(time_params['start_time'],time_params['end_time'],time_params['integration_time'])

    if time_params['Ntimes'] > 1:
        assert np.diff(time_arr)[0] == time_params['integration_time']

    param_dict['channel_width'] = time_params['integration_time']

    Nbl = (param_dict['Nants_data']+1)*param_dict['Nants_data'] / 2

    param_dict['time_array'] = np.sort(np.tile(time_arr,Nbl))
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

    bls = np.array([uv_obj.antnums_to_baseline(j,i,attempt256=True)
               for i in range(0,uv_obj.Nants_data) 
               for j in range(i,uv_obj.Nants_data) ])
    
    uv_obj.baseline_array = np.tile(bls, uv_obj.Ntimes)
    uv_obj.Nbls = bls.size
    uv_obj.Nblts = uv_obj.Nbls * uv_obj.Ntimes

    uv_obj.ant_1_array, uv_obj.ant_2_array = \
          uv_obj.baseline_to_antnums(uv_obj.baseline_array)

    return uv_obj, beam_list, beam_ids
