
import pyuvsim
from pyuvdata import UVBeam, UVData
from astropy.time import Time
import numpy as np, os, yaml
import nose.tools as nt
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
from test_uvsim import create_zenith_source, beam_files


EW_uvfits_file = os.path.join(SIM_DATA_PATH, '28mEWbl_10time_10chan.uvfits')
herabeam_default = os.path.join(SIM_DATA_PATH, 'HERA_NicCST.uvbeam')
param_filenames = [os.path.join(SIM_DATA_PATH,'test_config', 'param_10time_10chan_{}.yaml'.format(x)) for x in range(5)]   # Five different test configs


def test_param_reader():
    for n in range(4):
        yield (check_param_reader, n)

def check_param_reader(config_num):
    """
        tests initialize_uvdata_from_params
    """

    param_filename = param_filenames[config_num]
    print param_filename
    hera_uv = UVData()
    hera_uv.read_uvfits(EW_uvfits_file)

    time = Time(hera_uv.time_array[0], scale='utc', format='jd')
    sources = np.array([create_zenith_source(time, 'zensrc')])

    beam = UVBeam()
    beam.read_beamfits(herabeam_default)
    beam_list = [beam]

    expected_uvtask_list = pyuvsim.uvdata_to_task_list(hera_uv, sources, beam_list)

#    for pf in param_filenames:
    with open(param_filename, 'r') as pf:
        param_dict = yaml.safe_load(pf)

    param_dict['config_path'] = param_filename    #Ensure path is present

    
    ## Check default configuration
    uv_obj, new_beam_list, beam_ids = pyuvsim.initialize_uvdata_from_params(param_dict)
    uvtask_list = pyuvsim.uvdata_to_task_list(uv_obj, sources, new_beam_list)


    ## Tasks are not ordered in UVTask lists, so need to sort them.
    ## This is enabled by the comparison operator in UVTask

    uvtask_list = sorted(uvtask_list)
    expected_uvtask_list = sorted(expected_uvtask_list)

    nt.assert_true(uvtask_list == expected_uvtask_list)


if __name__ == '__main__':
    test_param_reader()
