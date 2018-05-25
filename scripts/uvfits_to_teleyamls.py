

## Read in a uvfits file and automatically generate yaml files. 
## Will assume the same beam_id for all antennas for now.

from pyuvdata import UVData
from pyuvdata.utils import ENU_from_ECEF
from astropy.io import fits
import yaml, argparse
import sys, os, numpy as np
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH

parser = argparse.ArgumentParser(description=("Extracts antenna position info from uvfits."))

parser.add_argument('file_in', metavar='<FILE>', type=str, nargs='+')
parser.add_argument('-l','--layout_fname', type=str, default='')
parser.add_argument('-b','--beam_file', type=str, default=os.path.join(SIM_DATA_PATH, 'HERA_NicCST.uvbeam'))
parser.add_argument('-y','--teleyaml_fname', type=str, default='')

args = parser.parse_args()

layout_fname = args.layout_fname
teleyaml_fname = args.teleyaml_fname
beam_file = args.beam_file
uv_fname = args.file_in[0]

if teleyaml_fname == '':
    teleyaml_fname = '.'.join(os.path.basename(uv_fname).split('.')[:-1]) + '.yaml'   #Replace extension

input_uv = UVData()
input_uv.read_uvfits(uv_fname, read_data=False)

if layout_fname == '':
    layout_fname = input_uv.telescope_name + "_layout.csv"

antpos_enu = ENU_from_ECEF((input_uv.antenna_positions + input_uv.telescope_location).T, *input_uv.telescope_location_lat_lon_alt).T

beam_filename = beam_file    ## beamfits file. For now, just the cst file.

e,n,u = antpos_enu.T

beam_ids = np.zeros_like(e).astype(int) #np.repeat(np.array(['beam0']), e.size)
antpos_full_array = np.vstack([input_uv.antenna_names, input_uv.antenna_numbers, beam_ids, e, n, u])
delim = ' '
col_width = max([len(name) for name in input_uv.antenna_names])
header = ("{:"+str(col_width)+"} {:8} {:8} {:10} {:10} {:10}\n").format("Name", "Number", "BeamID", "E", "N", "U")

with open(layout_fname, 'w') as lfile:
    lfile.write(header+'\n')
    for i in range(beam_ids.size):
        e,n,u = antpos_enu[i]
        beam_id = beam_ids[i]
        name = input_uv.antenna_names[i]
        num = input_uv.antenna_numbers[i]
        line = ("{:"+str(col_width)+"} {:8d} {:8d} {:10.4f} {:10.4f} {:10.4f}\n").format(name, num, beam_id, e, n, u)
        lfile.write(line)

### Write the rest to a yaml file.
yaml_dict = dict(
    telescope_name = input_uv.telescope_name,
    telescope_location = repr(input_uv.telescope_location_lat_lon_alt_degrees),
    Nants = input_uv.Nants_telescope,
    beam_paths = {
            0 : beam_file 
                }
)

with open(teleyaml_fname, 'w+') as yfile:
    yaml.dump(yaml_dict, yfile, default_flow_style=False)


print(('Filenames: ', teleyaml_fname, layout_fname))
