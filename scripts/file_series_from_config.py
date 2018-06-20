
## Make a series of simulation files from a given configuration

from pyuvsim import simsetup
import argparse
import yaml
import datetime
import os
from astropy.time import Time
from astropy.units import degree
from astropy.coordinates import EarthLocation, Angle
from pyuvdata import utils as uvutils
import numpy as np
import copy


def set_other_required_params(input_uv):
    input_uv.Nblts = input_uv.Nbls * input_uv.Ntimes
    enu = uvutils.ENU_from_ECEF((input_uv.antenna_positions + input_uv.telescope_location).T, *input_uv.telescope_location_lat_lon_alt).T
    uvws = []
    ant_1_arr = []
    ant_2_arr = []
    for i in range(input_uv.Nants_data):
        for j in range(i,input_uv.Nants_data):
            uvws.append(enu[j,:] - enu[i,:])
            ant_2_arr.append(i)
            ant_1_arr.append(j)
   
    uvws = np.array(uvws)
    input_uv.ant_1_array = np.tile(ant_1_arr, input_uv.Ntimes).T
    input_uv.ant_2_array = np.tile(ant_2_arr, input_uv.Ntimes).T
    input_uv.uvw_array = np.tile(uvws.T, input_uv.Ntimes).T
    input_uv.polarization_array  =  np.arange(-5,-9,-1)
    input_uv.Npols = 4
    input_uv.Nspws = 1
    input_uv.set_drift()
    input_uv.history = ''
    if input_uv.instrument is None: input_uv.instrument=input_uv.telescope_name
    input_uv.data_array = np.zeros((input_uv.Nblts, input_uv.Nspws, input_uv.Nfreqs, input_uv.Npols), dtype=np.complex)
    input_uv.nsample_array = np.zeros((input_uv.Nblts, input_uv.Nspws, input_uv.Nfreqs, input_uv.Npols)) + 1
    input_uv.flag_array = np.zeros((input_uv.Nblts, input_uv.Nspws, input_uv.Nfreqs, input_uv.Npols)).astype(bool)
    input_uv.antenna_positions = (uvutils.ECEF_from_ENU(enu.T,*input_uv.telescope_location_lat_lon_alt).T - input_uv.telescope_location)
    input_uv.set_lsts_from_time_array()
    input_uv.phase_to_time(input_uv.time_array[0])
    input_uv.spw_array = [0]
    input_uv.vis_units = 'Jy'
    input_uv.object_name='zenith'


def set_ofilename(parms, uvd):
    ofile = 'empty' if not 'testfile_prefix' in parms else parms['testfile_prefix']
    if not 'snapshot_number' in parms: en = ''
    else: en = parms['snapshot_number']

    if not 'Nfiles' in parms: nout=1
    else: nout = parms['Nfiles']
    odir = '.' if not 'outdir' in parms else parms['outdir']
    if not os.path.exists(odir): os.makedirs(odir)
    if parms['time_format'] == 'lst':
            t0str = "{:09.6f}".format(uvd.lst_array[0]*(12/np.pi))
            tfstr = "{:09.6f}".format(uvd.lst_array[-1]*(12/np.pi))
            label = "{:03d}".format(en)+"_"+t0str + "-" + tfstr
    else:
            label=str(en)
    if nout > 1:
         ofiu = ofile.split(".")[0] +"_"  +  label + ".uvfits"
    elif nout == 1:
         ofiu = ofile.split(".")[0] + ".uvfits"
    return os.path.join(odir,ofiu)

parser = argparse.ArgumentParser(description=("A command-line script "
                                              "to generate simulation files from parameters. "))

parser.add_argument('-p','--paramsfile',dest='paramsfile', type=str, help='Parameter yaml file.')
args = vars(parser.parse_args())

if 'paramsfile' not in args:
    raise KeyError("Parameter file required")

with open(args['paramsfile'], 'r') as pfile:
    params = yaml.safe_load(pfile)

time_params = params['time']

if not 'time_format' in time_params:
    time_params['time_format'] = 'julian'

if time_params['time_format'] == 'lst':
    
    tele_params = params['telescope']

    telescope_config_name = tele_params['telescope_config_name']

    if not os.path.exists(telescope_config_name):
        path = os.path.dirname(params['config_path'])
        telescope_config_name = os.path.join(path, telescope_config_name)

    with open(telescope_config_name, 'r') as yf:
            telparam = yaml.safe_load(yf)

    lat, lon, alt = eval(telparam['telescope_location'])   #Degrees

    loca = EarthLocation(lat=Angle(lat, degree), lon=Angle(lon, degree) )
    year   = 2017  if not 'year' in time_params else int(time_params['year'])
    month  = 10    if not 'month' in time_params else int(time_params['month'])
    day    = 28    if not 'day' in time_params else int(time_params['day'])
    hour   = 0     if not 'hour' in time_params else int(time_params['hour'])
    minute = 0     if not 'minute' in time_params else int(time_params['minute'])
    second = 0     if not 'second' in time_params else int(time_params['second'])
    basetime = Time(datetime.datetime(year,month,day,hour,minute,second),location=loca)


    dstart = (time_params['start_time'] - basetime.sidereal_time('apparent').hour)
    dend   = (time_params['end_time']   - basetime.sidereal_time('apparent').hour)
    tzero = basetime.jd + (1/24.)*dstart
    tfin =  basetime.jd + (1/24.)*dend

if time_params['time_format'] == 'julian':
    tzero =  Time(time_params['start_time'],format='jd').jd
    tfin  =  Time(time_params['end_time'],format='jd').jd

time_params['start_time'] = tzero
time_params['end_time'] = tfin

file_params = params['filing']
params.update(file_params)
del params['filing']

params['config_path'] = args['paramsfile']   #Ensure path is present

if 'snapshot_length_hours' in time_params:
    ### Splitting the output into a series of files.
    length_days = float(time_params['snapshot_length_hours'])*1/24.
    total_dur_hours = (tfin - tzero)*24.  #Hours
    Nfiles = int(np.floor(total_dur_hours/float(time_params['snapshot_length_hours'])) + 1)
    try:
        fi = int(os.environ['SLURM_ARRAY_TASK_ID'])
        params['snapshot_number'] = fi
        serial=False
        time_params['start_time'] = tzero + length_days * fi
        time_params['end_time'] = time_params['start_time'] + length_days
        if time_params['end_time'] > tfin:
            time_params['end_time'] = tfin
            time_params['Ntimes'] = int((tfin - tzero)/float(time_params['integration_time']))
    except KeyError:
        print "Warning: Running in serial"
        serial=True
else:
    serial = True
    Nfiles = 1
    length_days = tfin - tzero

params['Nfiles'] = Nfiles

if not serial:
    params['time'] = time_params
    input_uv, beam_list, beam_ids = simsetup.initialize_uvdata_from_params(params)
    set_other_required_params(input_uv)
    input_uv.history = 'Made on '+str(datetime.date.today())+ ' from config file '+args['paramsfile']+' by Adam Lanman'
    opath = set_ofilename(params, input_uv)
    input_uv.write_uvfits(opath,spoof_nonessential=True)
else:
    for fi in range(Nfiles):
        time_params['start_time'] = tzero + length_days * fi
        time_params['end_time'] = time_params['start_time'] + length_days
        if time_params['end_time'] > tfin:
            time_params['end_time'] = tfin
            time_params['Ntimes'] = int((tfin - time_params['start_time'])*3600.*24./float(time_params['integration_time']))
        params['time'] = time_params
        params['snapshot_number'] = fi
        input_uv, beam_list, beam_ids = simsetup.initialize_uvdata_from_params(params)
        set_other_required_params(input_uv)
        input_uv.history = 'Made on '+str(datetime.date.today())+ ' from config file '+args['paramsfile']+' by Adam Lanman'
        opath = set_ofilename(params, input_uv)
        print opath
#        input_uv.write_uvfits(opath,spoof_nonessential=True)


