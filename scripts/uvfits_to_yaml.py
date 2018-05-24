### Take a uvfits file, and save sim parameters as a yaml file.

from pyuvdata import UVData
import sys, os, yaml
from pyuvdata.data import DATA_PATH

def instr_to_beam(inst):
    """
        Choose a beam model file based on an instrument name.
    """
    if inst.lower() == 'hera':
        return os.path.join(DATA_PATH, 'HERA_NicCST_150MHz.txt')
    if inst.lower() == 'mwa':
        raise ValueError("Need to define a default MWA beam model")
    else:
        print("Warning: No beam model for "+str(inst)+", defaulting to HERA.")
        return instr_to_beam("hera")

uv_fname = sys.argv[1]
yaml_fname = '.'.join(uv_fname.split('.')[:-1]) + '.yaml'   #Replace extension

input_uv = UVData()
input_uv.read_uvfits(uv_fname, read_data=False)

freq_array = input_uv.freq_array[0,:].tolist()
time_array = input_uv.time_array.tolist()

param_dict = dict(

    start_time = time_array[0],
    end_time = time_array[-1],
    Ntimes = input_uv.Ntimes,

    integration_time = input_uv.integration_time,
    
    start_freq = freq_array[0],
    end_freq = freq_array[-1],
    channel_width = input_uv.channel_width,
    Nfreqs = input_uv.Nfreqs,
    
    instrument = input_uv.instrument,
    beam_model = instr_to_beam(input_uv.instrument)

)

with open(yaml_fname, 'w') as yfile:
    yaml.dump(param_dict, yfile, default_flow_style=False)
