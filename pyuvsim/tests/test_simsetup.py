
import pyuvsim
from pyuvdata import UVBeam, UVData
from astropy.time import Time
import numpy as np, os, yaml, shutil, copy
import nose.tools as nt
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
from test_uvsim import create_zenith_source, beam_files


EW_uvfits_file = os.path.join(SIM_DATA_PATH, '28mEWbl_10time_10chan.uvfits')
herabeam_default = os.path.join(SIM_DATA_PATH, 'HERA_NicCST.uvbeam')
param_filenames = [os.path.join(SIM_DATA_PATH,'test_config', 'param_10time_10chan_{}.yaml'.format(x)) for x in range(5)]   # Five different test configs


def compare_dictionaries(dic1, dic2):
    """
        Recursively compare two dictionaries.
    """
    compare=True
    for k in dic1.keys():
        if isinstance(dic1[k],dict):
            compare *= compare_dictionaries(dic1[k],dic2[k])
        else:
            compare *= (dic1[k] == dic2[k])
    return bool(compare)


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

def test_uvfits_to_yaml():
    """
        Loopback test of reading parameters from uvfits file, generating uvfits file, and reading in again.
    """
    opath = 'uvfits_yaml_temp'
    param_filename = 'test_config.yaml'
    second_param_filename = 'test2_config.yaml'
    test_uv_name = 'test_uv.uvfits'
    telescope_config = 'test_telescope_config.yaml'
    if not os.path.exists(opath):
        os.makedirs(opath)        #Directory will be deleted when test completed.
    
    #Read uvfits file to params.
    path, telescope_config, layout_fname = pyuvsim.simsetup.uvfits_to_telescope_config(EW_uvfits_file, herabeam_default, telescope_config_name=telescope_config, path=opath, return_names=True)
    pyuvsim.simsetup.uvfits_to_config_file(EW_uvfits_file, config_filename = param_filename, telescope_config_path = os.path.join(path, telescope_config), layout_csv_path = os.path.join(path, layout_fname), path=opath)
    
    # From parameters, generate a uvdata object.

    with open(os.path.join(opath,param_filename), 'r') as pf:
        param_dict = yaml.safe_load(pf)
    param_dict['config_path'] = param_filename    #Ensure path is present

    orig_param_dict = copy.deepcopy(param_dict)   # The parameter dictionary gets modified in the function below.
    uv_obj, new_beam_list, beam_ids = pyuvsim.initialize_uvdata_from_params(param_dict)

    # Write uvfits file to temp dir. Need to spoof some required parameters.
    uv_obj.Npols = 4
    uv_obj.Nspws = 1
    uv_obj.data_array = np.zeros((uv_obj.Nblts, uv_obj.Nspws, uv_obj.Nfreqs, uv_obj.Npols),dtype=np.complex)
    uv_obj.flag_array = np.ones_like(uv_obj.data_array).astype(bool)
    uv_obj.nsample_array = np.ones_like(uv_obj.data_array).astype(float)
    uv_obj.history = "Temporary file"
    uv_obj.instrument = 'HERA'
    uv_obj.set_lsts_from_time_array()
    uv_obj.object_name = 'temp'
    uv_obj.set_phased()
    uv_obj.polarization_array = [-5,-6,-7,-8]
    uv_obj.spw_array = [0]
    uv_obj.uvw_array = np.ones((uv_obj.Nblts, 3))
    uv_obj.vis_units = "Jy"
    uv_obj.phase_center_dec = 0.0
    uv_obj.phase_center_ra = 0.0
    uv_obj.phase_center_epoch = 2000.0

    uv_obj.write_uvfits(os.path.join(opath, test_uv_name), spoof_nonessential=True)

    # Generate parameters from new uvfits and compare with old.
    path, telescope_config, layout_fname = pyuvsim.simsetup.uvfits_to_telescope_config(os.path.join(opath,test_uv_name), herabeam_default, telescope_config_name=telescope_config, path=opath, return_names=True)
    pyuvsim.simsetup.uvfits_to_config_file(EW_uvfits_file, config_filename = second_param_filename, telescope_config_path = os.path.join(path, telescope_config), layout_csv_path = os.path.join(path, layout_fname), path=opath)

    del param_dict
    with open(os.path.join(path,second_param_filename), 'r') as pf:
        param_dict = yaml.safe_load(pf)

    param_dict['config_path'] = param_filename

    nt.assert_true(compare_dictionaries(param_dict,orig_param_dict))

    shutil.rmtree(opath)

if __name__ == '__main__':
    test_uvfits_to_yaml()
