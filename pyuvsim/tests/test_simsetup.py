# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from __future__ import absolute_import, division, print_function

import numpy as np
import os
import yaml
import shutil
import copy
from six.moves import map, range, zip
import nose.tools as nt
from astropy.time import Time
from astropy.coordinates import Angle, SkyCoord, EarthLocation
from astropy import units


from pyuvdata import UVBeam, UVData

import pyuvsim
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
import pyuvsim.tests as simtest


EW_uvfits_file = os.path.join(SIM_DATA_PATH, '28mEWbl_10time_10chan.uvfits')
herabeam_default = os.path.join(SIM_DATA_PATH, 'HERA_NicCST.uvbeam')
param_filenames = [os.path.join(SIM_DATA_PATH, 'test_config', 'param_10time_10chan_{}.yaml'.format(x)) for x in range(5)]   # Five different test configs
longbl_uvfits_file = os.path.join(SIM_DATA_PATH, '5km_triangle_1time_1chan.uvfits')
triangle_uvfits_file = os.path.join(SIM_DATA_PATH, '28m_triangle_10time_10chan.uvfits')
GLEAM_vot = os.path.join(SIM_DATA_PATH, 'gleam_50srcs.vot')


def test_setup_airy():
    # what is this line supposed to be checking??
    pyuvsim.initialize_uvdata_from_params(os.path.join(SIM_DATA_PATH,
                                                       'simple_equator_sim_airy.yaml'))
    nt.assert_raises(KeyError, pyuvsim.initialize_uvdata_from_params,
                     os.path.join(SIM_DATA_PATH, 'simple_equator_sim_airy_broken.yaml'))


def test_mock_catalog_zenith_source():

    time = Time(2457458.65410, scale='utc', format='jd')

    array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s',
                                   height=1073.)
    freq = (150e6 * units.Hz)

    source_coord = SkyCoord(alt=Angle(90 * units.deg), az=Angle(0 * units.deg),
                            obstime=time, frame='altaz', location=array_location)
    icrs_coord = source_coord.transform_to('icrs')

    ra = icrs_coord.ra
    dec = icrs_coord.dec

    test_source = pyuvsim.Source('src0', ra, dec, freq, [1, 0, 0, 0])

    cat, mock_keywords = pyuvsim.create_mock_catalog(time, arrangement='zenith')
    cat_source = cat[0]

    nt.assert_equal(cat_source, test_source)


def test_mock_catalog_off_zenith_source():

    src_az = Angle('90.0d')
    src_alt = Angle('85.0d')
    src_za = Angle('90.0d') - src_alt

    time = Time(2457458.65410, scale='utc', format='jd')

    array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s',
                                   height=1073.)
    freq = (150e6 * units.Hz)

    source_coord = SkyCoord(alt=src_alt, az=src_az,
                            obstime=time, frame='altaz', location=array_location)
    icrs_coord = source_coord.transform_to('icrs')

    ra = icrs_coord.ra
    dec = icrs_coord.dec
    test_source = pyuvsim.Source('src0', ra, dec, freq, [1.0, 0, 0, 0])

    cat, mock_keywords = pyuvsim.create_mock_catalog(time, arrangement='off-zenith', zen_ang=src_za.deg)
    cat_source = cat[0]

    nt.assert_equal(cat_source, test_source)


def check_param_reader(config_num):
    """
        tests initialize_uvdata_from_params
    """

    param_filename = param_filenames[config_num]
    hera_uv = UVData()
    hera_uv.read_uvfits(triangle_uvfits_file)
    hera_uv.telescope_name = 'HERA'

    time = Time(hera_uv.time_array[0], scale='utc', format='jd')
    sources, _ = pyuvsim.create_mock_catalog(time, arrangement='zenith')

    beam0 = UVBeam()
    beam0.read_beamfits(herabeam_default)
    beam1 = pyuvsim.AnalyticBeam('uniform')
    beam2 = pyuvsim.AnalyticBeam('gaussian', sigma=0.02)
    beam3 = pyuvsim.AnalyticBeam('airy', diameter=14.6)
    beam_list = [beam0, beam1, beam2, beam3]

    beam_dict = {'ANT1': 0, 'ANT2': 1, 'ANT3': 2, 'ANT4': 3}
    expected_uvtask_list = pyuvsim.uvdata_to_task_list(hera_uv, sources, beam_list, beam_dict=beam_dict)

    # Check error conditions:
    if config_num == 0:
        with open(param_filename, 'r') as pfile:
            params_bad = yaml.safe_load(pfile)
        params_bad['config_path'] = os.path.join(SIM_DATA_PATH, "test_config")
        params_bad['telescope']['telescope_config_name'] = os.path.join(SIM_DATA_PATH, 'test_config', '28m_triangle_10time_10chan_nosigma.yaml')
        nt.assert_raises(KeyError, pyuvsim.initialize_uvdata_from_params, params_bad)
        params_bad['telescope']['telescope_config_name'] = os.path.join(SIM_DATA_PATH, 'test_config', '28m_triangle_10time_10chan_nodiameter.yaml')
        nt.assert_raises(KeyError, pyuvsim.initialize_uvdata_from_params, params_bad)
        params_bad['telescope']['telescope_config_name'] = os.path.join(SIM_DATA_PATH, 'test_config', '28m_triangle_10time_10chan_nofile.yaml')
        nt.assert_raises(OSError, pyuvsim.initialize_uvdata_from_params, params_bad)

        # Errors on frequency configuration
        with open(param_filename, 'r') as pfile:
            params_bad = yaml.safe_load(pfile)
        # Define freq_arr but not channel_width
        params_bad['config_path'] = os.path.join(SIM_DATA_PATH, "test_config")
        params_bad['freq']['freq_array'] = np.array([1e8])
        del params_bad['freq']['channel_width']
        nt.assert_raises(ValueError, pyuvsim.initialize_uvdata_from_params, params_bad)
        del params_bad['freq']['freq_array']

        # Don't define Nfreqs or channel_width
        del params_bad['freq']['Nfreqs']
        nt.assert_raises(ValueError, pyuvsim.initialize_uvdata_from_params, params_bad)

        # Define Nfreqs but not bandwidth
        del params_bad['freq']['end_freq']  # Can't make bandwidth without start and end
        params_bad['freq']['Nfreqs'] = 10
        nt.assert_raises(ValueError, pyuvsim.initialize_uvdata_from_params, params_bad)

        # Now check time configuration:
        with open(param_filename, 'r') as pfile:
            params_bad = yaml.safe_load(pfile)
        params_bad['config_path'] = os.path.join(SIM_DATA_PATH, "test_config")

        # Don't define start or end time:
        del params_bad['time']['end_time']
        del params_bad['time']['start_time']
        nt.assert_raises(ValueError, pyuvsim.initialize_uvdata_from_params, params_bad)

        # Don't define Ntimes or integration_time
        del params_bad['time']['Ntimes']
        del params_bad['time']['integration_time']
        nt.assert_raises(ValueError, pyuvsim.initialize_uvdata_from_params, params_bad)

        params_bad['time']['Ntimes'] = 10
        nt.assert_raises(ValueError, pyuvsim.initialize_uvdata_from_params, params_bad)

    # Check default configuration
    uv_obj, new_beam_list, new_beam_dict, beam_ids = pyuvsim.initialize_uvdata_from_params(param_filename)

    # write_uvfits test:
    with open(param_filename, 'r') as fhandle:
        param_dict = yaml.safe_load(fhandle)
    expected_ofilepath = pyuvsim.simsetup.write_uvfits(uv_obj, param_dict, return_filename=True, dryrun=True)
    ofilename = 'sim_results.uvfits'
    if config_num == 1:
        if os.path.isdir('tempdir'):
            os.rmdir('tempdir')
        ofilename = os.path.join('.', 'tempdir', ofilename)
    else:
        ofilename = os.path.join('.', ofilename)
    nt.assert_equal(ofilename, expected_ofilepath)

    uvtask_list = pyuvsim.uvdata_to_task_list(uv_obj, sources, new_beam_list, beam_dict=new_beam_dict)
    # Tasks are not ordered in UVTask lists, so need to sort them.
    # This is enabled by the comparison operator in UVTask
    uvtask_list = sorted(uvtask_list)
    expected_uvtask_list = sorted(expected_uvtask_list)

    nt.assert_true(uvtask_list == expected_uvtask_list)


# This loops through different config files and tests all of them the same way
def test_param_reader():
    for n in range(5):
        yield (check_param_reader, n)


def test_uvfits_to_config():
    """
        Loopback test of reading parameters from uvfits file, generating uvfits file, and reading in again.
    """
    opath = 'uvfits_yaml_temp'
    param_filename = 'test_config.yaml'
    second_param_filename = 'test2_config.yaml'
    telescope_config = 'test_telescope_config.yaml'
    if not os.path.exists(opath):
        os.makedirs(opath)        # Directory will be deleted when test completed.

    # Read uvfits file to params.
    uv0 = UVData()
    # uv0.read_uvfits(EW_uvfits_file)
    uv0.read_uvfits(longbl_uvfits_file)
    path, telescope_config, layout_fname = \
        pyuvsim.simsetup.uvdata_to_telescope_config(uv0, herabeam_default,
                                                    telescope_config_name=telescope_config,
                                                    path_out=opath, return_names=True)
    pyuvsim.simsetup.uvdata_to_config_file(uv0, param_filename=param_filename,
                                           telescope_config_name=os.path.join(path, telescope_config),
                                           layout_csv_name=os.path.join(path, layout_fname),
                                           path_out=opath)
    # From parameters, generate a uvdata object.

    with open(os.path.join(opath, param_filename), 'r') as pf:
        param_dict = yaml.safe_load(pf)
    param_dict['config_path'] = opath    # Ensure path is present

    orig_param_dict = copy.deepcopy(param_dict)   # The parameter dictionary gets modified in the function below.
    uv1, new_beam_list, new_beam_dict, beam_ids = pyuvsim.initialize_uvdata_from_params(param_dict)

    # Generate parameters from new uvfits and compare with old.
    path, telescope_config, layout_fname = \
        pyuvsim.simsetup.uvdata_to_telescope_config(uv1, herabeam_default,
                                                    telescope_config_name=telescope_config,
                                                    layout_csv_name=layout_fname,
                                                    path_out=opath, return_names=True)
    pyuvsim.simsetup.uvdata_to_config_file(uv1, param_filename=second_param_filename,
                                           telescope_config_name=os.path.join(path, telescope_config),
                                           layout_csv_name=os.path.join(path, layout_fname),
                                           path_out=opath)

    del param_dict
    with open(os.path.join(path, second_param_filename), 'r') as pf:
        param_dict = yaml.safe_load(pf)

    nt.assert_true(simtest.compare_dictionaries(param_dict, orig_param_dict))

    shutil.rmtree(opath)


def test_point_catalog_reader():
    catfile = os.path.join(SIM_DATA_PATH, 'test_config', 'pointsource_catalog.txt')
    catalog = pyuvsim.simsetup.read_text_catalog(catfile)

    with open(catfile, 'r') as fhandle:
        header = fhandle.readline()
    header = [h.strip() for h in header.split()]
    dt = np.format_parser(['a10', 'f8', 'f8', 'f8', 'f8'],
                          ['source_id', 'ra_j2000', 'dec_j2000', 'flux_density_I', 'frequency'], header)

    catalog_table = np.genfromtxt(catfile, autostrip=True, skip_header=1,
                                  dtype=dt.dtype)

    for src in catalog:
        nt.assert_true(src.name in catalog_table['source_id'])
        nt.assert_true(src.ra.deg in catalog_table['ra_j2000'])
        nt.assert_true(src.dec.deg in catalog_table['dec_j2000'])
        nt.assert_true(src.stokes[0] in catalog_table['flux_density_I'])
        nt.assert_true(src.freq.to("Hz").value in catalog_table['frequency'])


def test_read_gleam():

    sourcelist = pyuvsim.simsetup.read_votable_catalog(GLEAM_vot)

    nt.assert_equal(len(sourcelist), 50)


def test_file_namer():
    """
    File name incrementer utility
    """
    existing_file = param_filenames[0]
    new_filepath = pyuvsim.simsetup.check_file_exists_and_increment(existing_file)
    print(new_filepath)
    print(existing_file)
    nt.assert_true(new_filepath.endswith("_5.yaml"))    # There are four other of these param test files


def test_mock_catalogs():
    time = Time(2457458.1739, scale='utc', format='jd')
    cat1, mock_kwds1 = pyuvsim.simsetup.create_mock_catalog(time, 'zenith')
    cat2, mock_kwds2 = pyuvsim.simsetup.create_mock_catalog(time, 'off-zenith')
    cat3, mock_kwds3 = pyuvsim.simsetup.create_mock_catalog(time, 'cross')
    cat4, mock_kwds4 = pyuvsim.simsetup.create_mock_catalog(time, 'triangle')
    cat5, mock_kwds5 = pyuvsim.simsetup.create_mock_catalog(time, 'long-line')
    cat6, mock_kwds6 = pyuvsim.simsetup.create_mock_catalog(time, 'hera_text')

    nt.assert_true(len(cat1) == 1)
    nt.assert_true(len(cat2) == 1)
    nt.assert_true(len(cat3) == 4)
    nt.assert_true(len(cat4) == 3)
    nt.assert_true(len(cat5) == 10)
    nt.assert_true(len(cat6) == 43)
