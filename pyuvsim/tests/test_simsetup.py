# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import copy
import os
import shutil

import numpy as np
import pytest
import pyuvdata.tests as uvtest
import yaml
from astropy import units
from astropy.coordinates import Angle, SkyCoord, EarthLocation
from astropy.time import Time
from pyuvdata import UVBeam, UVData
import pyradiosky

import pyuvsim
import pyuvsim.tests as simtest
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH

herabeam_default = os.path.join(SIM_DATA_PATH, 'HERA_NicCST.uvbeam')

# Five different test configs
param_filenames = [
    os.path.join(SIM_DATA_PATH, 'test_config', 'param_10time_10chan_{}.yaml'.format(x))
    for x in range(6)
]

longbl_uvfits_file = os.path.join(SIM_DATA_PATH, '5km_triangle_1time_1chan.uvfits')
triangle_uvfits_file = os.path.join(SIM_DATA_PATH, '28m_triangle_10time_10chan.uvfits')
manytimes_config = os.path.join(
    SIM_DATA_PATH, 'test_config', 'param_100times_1.5days_triangle.yaml'
)
gleam_param_file = os.path.join(SIM_DATA_PATH, 'test_config', 'param_1time_1src_testgleam.yaml')


def test_mock_catalog_zenith_source():
    time = Time(2457458.65410, scale='utc', format='jd')

    array_location = EarthLocation(
        lat='-30d43m17.5s', lon='21d25m41.9s', height=1073.
    )

    source_coord = SkyCoord(
        alt=Angle(90 * units.deg), az=Angle(0 * units.deg),
        obstime=time, frame='altaz', location=array_location
    )
    icrs_coord = source_coord.transform_to('icrs')

    ra = icrs_coord.ra
    dec = icrs_coord.dec

    test_source = pyradiosky.SkyModel('src0', ra, dec, [1, 0, 0, 0], [1e8], 'flat')

    cat, mock_keywords = pyuvsim.create_mock_catalog(time, arrangement='zenith')

    assert cat == test_source


def test_mock_catalog_off_zenith_source():
    src_az = Angle('90.0d')
    src_alt = Angle('85.0d')

    time = Time(2457458.65410, scale='utc', format='jd')

    array_location = EarthLocation(
        lat='-30d43m17.5s', lon='21d25m41.9s', height=1073.
    )

    source_coord = SkyCoord(
        alt=src_alt, az=src_az, obstime=time, frame='altaz', location=array_location
    )
    icrs_coord = source_coord.transform_to('icrs')

    ra = icrs_coord.ra
    dec = icrs_coord.dec
    test_source = pyradiosky.SkyModel('src0', ra, dec, [1.0, 0, 0, 0], [1e8], 'flat')

    cat, mock_keywords = pyuvsim.create_mock_catalog(time, arrangement='off-zenith',
                                                     alt=src_alt.deg)

    assert cat == test_source


def test_catalog_from_params():
    # Pass in parameter dictionary as dict
    hera_uv = UVData()
    with pytest.warns(UserWarning) as telwarn:
        hera_uv.read_uvfits(triangle_uvfits_file)
    assert str(telwarn.pop().message).startswith('Telescope 28m_triangle_10time_10chan.yaml is not in known_telescopes.')

    source_dict = {}
    simtest.assert_raises_message(
        KeyError, "No catalog defined.",
        pyuvsim.simsetup.initialize_catalog_from_params,
        {'sources': source_dict}
    )
    pytest.raises(
        KeyError, pyuvsim.simsetup.initialize_catalog_from_params, {'sources': source_dict}
    )
    arrloc = '{:.5f},{:.5f},{:.5f}'.format(*hera_uv.telescope_location_lat_lon_alt_degrees)
    source_dict = {
        'catalog': 'mock',
        'mock_arrangement': 'zenith',
        'Nsrcs': 5,
        'time': hera_uv.time_array[0]
    }
    with pytest.warns(UserWarning) as warn:
        pyuvsim.simsetup.initialize_catalog_from_params({'sources': source_dict})
    assert str(warn.pop().message).startswith("No array_location specified. Defaulting to the HERA site.")
    catalog_uv, srclistname = pyuvsim.simsetup.initialize_catalog_from_params(
        {'sources': source_dict}, hera_uv
    )
    catalog_uv = pyradiosky.array_to_skymodel(catalog_uv)
    source_dict['array_location'] = arrloc
    del source_dict['time']

    simtest.assert_raises_message(TypeError, 'input_uv must be UVData object',
                                  pyuvsim.simsetup.initialize_catalog_from_params,
                                  {'sources': source_dict},
                                  input_uv='not_uvdata')
    simtest.assert_raises_message(ValueError,
                                  'input_uv must be supplied if using mock catalog without '
                                  'specified julian date',
                                  pyuvsim.simsetup.initialize_catalog_from_params,
                                  {'sources': source_dict})
    with pytest.warns(UserWarning) as warn:
     catalog_str, srclistname2 = pyuvsim.simsetup.initialize_catalog_from_params({'sources': source_dict}, hera_uv)
    assert str(warn.pop().message).startswith("Warning: No julian date given for mock catalog. Defaulting to first time step.")
    
    catalog_str = pyradiosky.array_to_skymodel(catalog_str)
    assert np.all(catalog_str == catalog_uv)

# parametrize will loop over all the give values
@pytest.mark.parametrize("config_num", [0, 2])
def test_param_reader(config_num):

    # Reading in various configuration files

    param_filename = param_filenames[config_num]
    hera_uv = UVData()
    with pytest.warns(UserWarning) as warn:
        hera_uv.read_uvfits(triangle_uvfits_file)
    assert str(warn.pop().message).startswith('Telescope 28m_triangle_10time_10chan.yaml is not in known_telescopes.')

    hera_uv.telescope_name = 'HERA'
    if config_num == 5:
        hera_uv.select(bls=[(0, 1), (1, 2)])

    time = Time(hera_uv.time_array[0], scale='utc', format='jd')
    sources, _ = pyuvsim.create_mock_catalog(time, arrangement='zenith', return_table=True)

    beam0 = UVBeam()
    beam0.read_beamfits(herabeam_default)
    beam1 = pyuvsim.AnalyticBeam('uniform')
    beam2 = pyuvsim.AnalyticBeam('gaussian', sigma=0.02)
    beam3 = pyuvsim.AnalyticBeam('airy', diameter=14.6)
    beam_list = [beam0, beam1, beam2, beam3]

    beam_dict = {'ANT1': 0, 'ANT2': 1, 'ANT3': 2, 'ANT4': 3}
    Ntasks = hera_uv.Nblts * hera_uv.Nfreqs
    taskiter = pyuvsim.uvdata_to_task_iter(range(Ntasks), hera_uv, sources,
                                           beam_list, beam_dict=beam_dict)
    expected_uvtask_list = list(taskiter)

    # Check error conditions:
    if config_num == 0:
        params_bad = pyuvsim.simsetup._config_str_to_dict(param_filename)
        bak_params_bad = copy.deepcopy(params_bad)

        # Missing config file info
        params_bad['config_path'] = os.path.join(
            SIM_DATA_PATH, 'nonexistent_directory', 'nonexistent_file'
        )
        simtest.assert_raises_message(
            ValueError, 'nonexistent_directory is not a directory',
            pyuvsim.initialize_uvdata_from_params, params_bad
        )

        params_bad['config_path'] = os.path.join(SIM_DATA_PATH, "test_config")
        params_bad['telescope']['array_layout'] = 'nonexistent_file'
        simtest.assert_raises_message(
            ValueError, 'nonexistent_file from yaml does not exist',
            pyuvsim.initialize_uvdata_from_params, params_bad
        )

        params_bad['telescope']['telescope_config_name'] = 'nonexistent_file'
        simtest.assert_raises_message(
            ValueError, 'telescope_config_name file from yaml does not exist',
            pyuvsim.initialize_uvdata_from_params, params_bad
        )

        # Missing beam keywords
        params_bad = copy.deepcopy(bak_params_bad)

        params_bad['config_path'] = os.path.join(SIM_DATA_PATH, "test_config")

        params_bad = copy.deepcopy(bak_params_bad)
        params_bad['telescope']['telescope_config_name'] = os.path.join(
            SIM_DATA_PATH, 'test_config', '28m_triangle_10time_10chan_gaussnoshape.yaml'
        )
        simtest.assert_raises_message(
            KeyError, 'Missing shape parameter for gaussian beam (diameter or sigma).',
            pyuvsim.initialize_uvdata_from_params, params_bad
        )
        params_bad['telescope']['telescope_config_name'] = os.path.join(
            SIM_DATA_PATH, 'test_config', '28m_triangle_10time_10chan_nodiameter.yaml'
        )
        simtest.assert_raises_message(
            KeyError, 'Missing diameter for airy beam.',
            pyuvsim.initialize_uvdata_from_params, params_bad
        )
        params_bad['telescope']['telescope_config_name'] = os.path.join(
            SIM_DATA_PATH, 'test_config', '28m_triangle_10time_10chan_nofile.yaml'
        )
        simtest.assert_raises_message(ValueError, 'Undefined beam model',
                                      pyuvsim.initialize_uvdata_from_params, params_bad)

    # Check default configuration
    uv_obj, new_beam_list, new_beam_dict = pyuvsim.initialize_uvdata_from_params(param_filename)
    for i, bm in enumerate(new_beam_list):
        new_beam_list[i] = pyuvsim.simsetup.beam_string_to_object(bm)
    # write_uvdata tests with different configs:
    with open(param_filename, 'r') as fhandle:
        param_dict = yaml.safe_load(fhandle)
    expected_ofilepath = pyuvsim.utils.write_uvdata(
        uv_obj, param_dict, return_filename=True, dryrun=True
    )
    ofilename = 'sim_results.uvfits'
    if config_num == 1:
        if os.path.isdir('tempdir'):
            os.rmdir('tempdir')
        ofilename = os.path.join('.', 'tempdir', ofilename)
    else:
        ofilename = os.path.join('.', ofilename)
    assert ofilename == expected_ofilepath

    Ntasks = uv_obj.Nblts * uv_obj.Nfreqs
    taskiter = pyuvsim.uvdata_to_task_iter(
        range(Ntasks), hera_uv, sources, beam_list, beam_dict=beam_dict
    )
    uvtask_list = list(taskiter)

    # Tasks are not ordered in UVTask lists, so need to sort them.
    uvtask_list = sorted(uvtask_list)
    expected_uvtask_list = sorted(expected_uvtask_list)
    assert uvtask_list == expected_uvtask_list


def test_tele_parser():
    """
    Check a variety of cases not already tested by param reader
    """
    # check no tele config passed
    tdict = dict(array_layout=os.path.join(SIM_DATA_PATH, 'test_layout_6ant.csv'))
    tel_error = 'If telescope_config_name not provided in `telescope` obsparam section, ' \
                'you must provide telescope_location'
    simtest.assert_raises_message(
        KeyError, tel_error, pyuvsim.simsetup.parse_telescope_params, tdict
    )
    tdict['telescope_location'] = '(-30.72152777777791, 21.428305555555557, 1073.0000000093132)'
    tel_error = 'If telescope_config_name not provided in `telescope` obsparam section, ' \
                'you must provide telescope_name'
    simtest.assert_raises_message(
        KeyError, tel_error, pyuvsim.simsetup.parse_telescope_params, tdict
    )

    tdict['telescope_name'] = 'tele'
    tpars, blist, bdict = pyuvsim.simsetup.parse_telescope_params(tdict)
    assert tpars['Nants_data'] == 6
    assert blist == []
    assert bdict == {}

    tdict.pop('array_layout')
    simtest.assert_raises_message(
        KeyError, 'array_layout must be provided.',
        pyuvsim.simsetup.parse_telescope_params, tdict
    )


def test_freq_parser():
    """
    Check all valid input parameter cases for frequencies.
    """

    fdict_base = dict(
        Nfreqs=10,
        channel_width=0.5,
        start_freq=0.0,
        end_freq=4.5,
        bandwidth=5.0)

    freq_array = np.linspace(
        fdict_base['start_freq'],
        fdict_base['start_freq'] + fdict_base['bandwidth'] - fdict_base['channel_width'],
        fdict_base['Nfreqs'], endpoint=True
    )

    fdict_base['freq_array'] = freq_array

    # As long as one tuple from each set is represented,
    # the param parser should work.

    bpass_kwd_combos = [
        ('start_freq', 'end_freq'), ('channel_width', 'Nfreqs'), ('bandwidth',)
    ]
    chwid_kwd_combos = [('bandwidth', 'Nfreqs'), ('channel_width',)]
    ref_freq_combos = [('start_freq',), ('end_freq',)]

    for bpass in bpass_kwd_combos:
        for chwid in chwid_kwd_combos:
            for ref in ref_freq_combos:
                keys = tuple(set(bpass + chwid + (ref)))  # Get unique keys
                subdict = {key: fdict_base[key] for key in keys}
                test = pyuvsim.parse_frequency_params(subdict)
                assert np.allclose(test['freq_array'][0], freq_array)

    # Now check error cases
    err_cases = [
        ('bandwidth',),
        ('start_freq', 'Nfreqs'),
        ('start_freq', 'channel_width'),
        ('start_freq', 'end_freq')
    ]
    err_mess = [
        'Either start or end frequency must be specified: bandwidth',
        'Either bandwidth or channel width must be specified: Nfreqs, start_freq',
        'Either bandwidth or band edges must be specified: channel_width, start_freq',
        'Either channel_width or Nfreqs  must be included in parameters:end_freq, '
        'start_freq'
    ]
    for ei, er in enumerate(err_cases):
        subdict = {key: fdict_base[key] for key in er}
        simtest.assert_raises_message(
            ValueError, err_mess[ei], pyuvsim.parse_frequency_params, subdict
        )

    subdict = {'freq_array': freq_array[0]}
    simtest.assert_raises_message(
        ValueError, 'Channel width must be specified if freq_arr has length 1',
        pyuvsim.parse_frequency_params, subdict
    )

    subdict = {'freq_array': freq_array[[0, 1, 4, 8]]}
    simtest.assert_raises_message(ValueError, 'Spacing in frequency array is uneven.',
                                  pyuvsim.parse_frequency_params,
                                  subdict)

    subdict = {'channel_width': 3.14, 'start_freq': 1.0, 'end_freq': 8.3}
    simtest.assert_raises_message(
        ValueError, 'end_freq - start_freq must be evenly divisible by channel_width',
        pyuvsim.parse_frequency_params, subdict
    )

    subdict = fdict_base.copy()
    subdict['Nfreqs'] = 7
    del subdict['freq_array']
    simtest.assert_raises_message(
        ValueError, 'Frequency array spacings are not equal to channel width.',
        pyuvsim.parse_frequency_params, subdict
    )


def test_time_parser():
    """
    Check a variety of cases for the time parser.
    """

    daysperhour = 1 / 24.
    dayspersec = 1 / (24 * 3600.)

    tdict_base = {
        'Ntimes': 24,
        'duration_hours': 0.9999999962747097 / daysperhour,
        'end_time': 2457458.9583333298,
        'integration_time': 3599.999986588955,
        'start_time': 2457458.0
    }

    inttime_days = tdict_base['integration_time'] * dayspersec
    time_array = np.linspace(
        tdict_base['start_time'] + inttime_days / 2.,
        tdict_base['start_time'] + tdict_base['duration_hours'] * daysperhour - inttime_days / 2.,
        tdict_base['Ntimes'], endpoint=True
    )

    tdict_base['time_array'] = time_array

    # As long as one tuple from each set is represented,
    # the param parser should work.

    bpass_kwd_combos = [
        ('start_time', 'end_time'),
        ('integration_time', 'Ntimes'),
        ('duration_hours',)
    ]
    chwid_kwd_combos = [('duration_hours', 'Ntimes'), ('integration_time',)]
    ref_freq_combos = [('start_time',), ('end_time',)]

    for bpass in bpass_kwd_combos:
        for chwid in chwid_kwd_combos:
            for ref in ref_freq_combos:
                keys = tuple(set(bpass + chwid + (ref)))  # Get unique keys
                subdict = {key: tdict_base[key] for key in keys}
                test = pyuvsim.parse_time_params(subdict)
                assert np.allclose(test['time_array'], time_array, atol=dayspersec)

    subdict = {'time_array': time_array}
    test = pyuvsim.parse_time_params(subdict)
    assert np.allclose(test['time_array'], time_array, atol=dayspersec)

    # Now check error cases
    err_cases = [
        ('duration_hours',),
        ('start_time', 'Ntimes'),
        ('start_time', 'integration_time'),
        ('start_time', 'end_time')
    ]
    err_mess = [
        'Start or end time must be specified: duration_hours',
        'Either duration or integration time must be specified: Ntimes, start_time',
        'Either duration or time bounds must be specified: integration_time, start_time',
        'Either integration_time or Ntimes must be included in parameters: end_time, '
        'start_time'
    ]

    for ei, er in enumerate(err_cases):
        subdict = {key: tdict_base[key] for key in er}
        simtest.assert_raises_message(ValueError, err_mess[ei], pyuvsim.parse_time_params, subdict)

    subdict = {'integration_time': 3.14, 'start_time': 10000.0, 'end_time': 80000.3, 'Ntimes': 30}
    simtest.assert_raises_message(
        ValueError, 'Calculated time array is not consistent with set '
        'integration_time', pyuvsim.parse_time_params, subdict
    )

    subdict = tdict_base.copy()
    subdict['Ntimes'] = 7
    del subdict['time_array']
    simtest.assert_raises_message(
        ValueError, 'Calculated time array is not consistent with set '
        'integration_time.', pyuvsim.parse_time_params, subdict
    )


def test_single_input_time():
    time_dict = pyuvsim.simsetup.time_array_to_params([1.0])
    assert time_dict['integration_time'] == 1.0


@pytest.fixture(scope='module')
def times_and_freqs():
    freqs = np.linspace(100, 200, 1024)
    times = np.linspace(2458570, 2458570 + 0.5, 239)
    # yield the time and frequency arrays to the tests
    # then delete after
    yield times, freqs

    del (times, freqs)


def test_freq_time_params_match(times_and_freqs):
    times, freqs = times_and_freqs
    time_dict = pyuvsim.simsetup.time_array_to_params(times)
    freq_dict = pyuvsim.simsetup.freq_array_to_params(freqs)
    ftest = pyuvsim.simsetup.parse_frequency_params(freq_dict)
    ttest = pyuvsim.simsetup.parse_time_params(time_dict)
    assert np.allclose(ftest['freq_array'], freqs)
    assert np.allclose(ttest['time_array'], times)


def test_uneven_time_array_to_params(times_and_freqs):
    times, freqs = times_and_freqs
    # Check that this works for unevenly-spaced times
    times = np.random.choice(times, 150, replace=False)
    times.sort()
    time_dict = pyuvsim.simsetup.time_array_to_params(times)
    ttest = pyuvsim.simsetup.parse_time_params(time_dict)
    assert np.allclose(ttest['time_array'], times)


def test_single_time_array_to_params(times_and_freqs):
    times, freqs = times_and_freqs
    # check Ntimes = 1 and Nfreqs = 1 case
    times = np.linspace(2458570, 2458570.5, 1)
    tdict = pyuvsim.simsetup.time_array_to_params(times)
    assert tdict['Ntimes'] == 1
    assert tdict['start_time'] == times


def test_single_freq_array_to_params(times_and_freqs):
    freqs = np.linspace(100, 200, 1)
    fdict = pyuvsim.simsetup.freq_array_to_params(freqs)
    assert fdict['Nfreqs'] == 1
    assert fdict['start_freq'] == freqs


def test_param_select_cross():
    param_filename = os.path.join(SIM_DATA_PATH, 'test_config', 'obsparam_mwa_nocore.yaml')
    param_dict = pyuvsim.simsetup._config_str_to_dict(param_filename)
    uv_obj_full, new_beam_list, new_beam_dict = pyuvsim.initialize_uvdata_from_params(param_dict)

    # test only keeping cross pols
    param_dict['select'] = {'ant_str': 'cross'}
    uv_obj_cross, new_beam_list, new_beam_dict = pyuvsim.initialize_uvdata_from_params(param_dict)
    uv_obj_cross2 = uv_obj_full.select(ant_str='cross', inplace=False)

    assert uv_obj_cross == uv_obj_cross2


def test_param_select_bls():
    param_filename = os.path.join(SIM_DATA_PATH, 'test_config', 'obsparam_mwa_nocore.yaml')
    param_dict = pyuvsim.simsetup._config_str_to_dict(param_filename)
    uv_obj_full, new_beam_list, new_beam_dict = pyuvsim.initialize_uvdata_from_params(param_dict)

    # test only keeping certain baselines
    param_dict['select'] = {'bls': '[(40, 41), (42, 43), (44, 45)]'}  # Test as string
    uv_obj_bls, new_beam_list, new_beam_dict = pyuvsim.initialize_uvdata_from_params(param_dict)

    uv_obj_bls2 = uv_obj_full.select(
        bls=[(40, 41), (42, 43), (44, 45)], inplace=False
    )
    uv_obj_bls.history, uv_obj_bls2.history = '', ''
    assert uv_obj_bls == uv_obj_bls2

    param_dict['object_name'] = 'foo'
    uv_obj_full, new_beam_list, new_beam_dict = pyuvsim.initialize_uvdata_from_params(param_dict)
    assert uv_obj_full.object_name == 'foo'


def test_param_select_redundant():
    param_filename = os.path.join(SIM_DATA_PATH, 'test_config', 'obsparam_hex37_14.6m.yaml')
    param_dict = pyuvsim.simsetup._config_str_to_dict(param_filename)
    uv_obj_full, new_beam_list, new_beam_dict = pyuvsim.initialize_uvdata_from_params(param_dict)

    # test only keeping one baseline per redundant group
    param_dict['select'] = {'redundant_threshold': 0.1}
    uv_obj_red, new_beam_list, new_beam_dict = pyuvsim.initialize_uvdata_from_params(param_dict)
    uv_obj_red2 = uv_obj_full.compress_by_redundancy(tol=0.1, inplace=False)
    uv_obj_red.history, uv_obj_red2.history = '', ''

    assert uv_obj_red == uv_obj_red2
    assert uv_obj_red.Nbls < uv_obj_full.Nbls


@pytest.mark.parametrize('case', np.arange(6))
def check_uvdata_keyword_init(case):
    base_kwargs = dict(
        antenna_layout_filepath=os.path.join(SIM_DATA_PATH, "test_config/triangle_bl_layout.csv"),
        telescope_location=(-30.72152777777791, 21.428305555555557, 1073.0000000093132),
        telescope_name="HERA", Nfreqs=10, start_freq=1e8, bandwidth=1e8, Ntimes=60,
        integration_time=100.0, start_time=2458101.0, polarization_array=['xx'],
        no_autos=True, write_files=False, run_check=True
    )

    if case == 0:
        # check it runs through
        uvd = pyuvsim.simsetup.initialize_uvdata_from_keywords(**base_kwargs)
        assert np.allclose(base_kwargs['telescope_location'],
                           uvd.telescope_location_lat_lon_alt_degrees)
        assert np.allclose(base_kwargs['integration_time'], uvd.integration_time)
        assert base_kwargs['telescope_name'] == uvd.telescope_name
        assert base_kwargs['start_freq'] == uvd.freq_array[0, 0]
        assert base_kwargs['start_time'] == uvd.time_array[0]
        assert base_kwargs['Ntimes'] == uvd.Ntimes
        assert base_kwargs['Nfreqs'] == uvd.Nfreqs
        assert base_kwargs['polarizations'] == uvd.get_pols()
        assert not np.any(uvd.ant_1_array == uvd.ant_2_array)

    elif case == 1:
        # check bls and antenna_nums selections work
        bls = [(1, 0), (2, 0), (3, 0)]
        new_kwargs = copy.deepcopy(base_kwargs)
        new_kwargs['bls'] = bls
        uvd = pyuvsim.simsetup.initialize_uvdata_from_keywords(**new_kwargs)
        antpairs = uvd.get_antpairs()
        assert antpairs == bls

    elif case == 2:
        # also check that '1' gets converted to [1]
        new_kwargs = copy.deepcopy(base_kwargs)
        new_kwargs['polarization_array'] = ['xx', 'yy']
        new_kwargs['no_autos'] = False
        new_kwargs['antenna_nums'] = '1'
        uvd = pyuvsim.simsetup.initialize_uvdata_from_keywords(**new_kwargs)

        assert uvd.Nbls == 1
        assert uvd.get_pols() == new_kwargs['polarizations']
    elif case == 3:
        # check time and freq array definitions supersede other parameters
        fa = np.linspace(100, 200, 11) * 1e6
        ta = np.linspace(2458101, 2458102, 21)
        new_kwargs = copy.deepcopy(base_kwargs)
        new_kwargs['freq_array'] = fa
        new_kwargs['time_array'] = ta
        uvd = pyuvsim.simsetup.initialize_uvdata_from_keywords(**new_kwargs)

        assert np.allclose(uvd.time_array[::uvd.Nbls], ta)
        assert np.allclose(uvd.freq_array[0], fa)
    elif case == 4:
        # test feeding array layout as dictionary
        uvd = pyuvsim.simsetup.initialize_uvdata_from_keywords(**base_kwargs)
        antpos, ants = uvd.get_ENU_antpos()
        antpos_d = dict(zip(ants, antpos))
        layout_fname = 'temp_layout.csv'
        obsparam_fname = 'temp_obsparam.yaml'

        new_kwargs = copy.deepcopy(base_kwargs)
        new_kwargs['output_layout_filename'] = layout_fname
        new_kwargs['output_yaml_filename'] = obsparam_fname
        new_kwargs['array_layout'] = antpos_d
        new_kwargs['path_out'] = simtest.TESTDATA_PATH
        new_kwargs['write_files'] = True

        uvd = pyuvsim.simsetup.initialize_uvdata_from_keywords(**new_kwargs)
        layout_path = os.path.join(simtest.TESTDATA_PATH, layout_fname)
        obsparam_path = os.path.join(simtest.TESTDATA_PATH, obsparam_fname)
        assert os.path.exists(layout_path)
        assert os.path.exists(obsparam_path)

        assert uvd.Nbls == 6
        assert uvd.Nants_data == 4
        ap, a = uvd.get_ENU_antpos()
        apd = dict(zip(a, ap))
        assert np.all([np.isclose(antpos_d[ant], apd[ant]) for ant in ants])

    elif case == 5:
        # Check defaults when writing to file.
        uvd = pyuvsim.simsetup.initialize_uvdata_from_keywords(**base_kwargs)
        antpos, ants = uvd.get_ENU_antpos()
        antpos_d = dict(zip(ants, antpos))
        # Checking -- Default to a copy of the original layout, if layout is provided.
        layout_fname = 'triangle_bl_layout.csv'
        obsparam_fname = 'obsparam.yaml'

        new_kwargs = copy.deepcopy(base_kwargs)
        new_kwargs['write_files'] = True

        pyuvsim.simsetup.initialize_uvdata_from_keywords(**new_kwargs)
        assert os.path.exists(layout_fname)
        assert os.path.exists(obsparam_fname)

        os.remove(layout_fname)
        os.remove(obsparam_fname)

        # Default if no antenna_layout_filepath is provided.
        layout_fname = 'antenna_layout.csv'
        new_kwargs.pop('antenna_layout_filepath')
        new_kwargs['array_layout'] = antpos_d
        new_kwargs['complete'] = True

        pyuvsim.simsetup.initialize_uvdata_from_keywords(**new_kwargs)
        assert os.path.exists(layout_fname)
        assert os.path.exists(obsparam_fname)

        os.remove(layout_fname)
        os.remove(obsparam_fname)


def test_uvfits_to_config():
    """
        Loopback test of reading parameters from uvfits file, generating uvfits file, and reading
        in again.
    """
    opath = 'uvfits_yaml_temp'
    param_filename = 'obsparam.yaml'
    second_param_filename = 'test2_config.yaml'
    if not os.path.exists(opath):
        os.makedirs(opath)  # Directory will be deleted when test completed.

    # Read uvfits file to params.
    uv0 = UVData()
    uv0.read_uvfits(longbl_uvfits_file)

    warningmessages = [
        'The default for the `center` keyword has changed. Previously it defaulted to True, '
        'using the median antennna location; now it defaults to False, using the '
        'telescope_location.',
        'The xyz array in ENU_from_ECEF is being interpreted as (Npts, 3). Historically this '
        'function has supported (3, Npts) arrays, please verify that array ordering is as '
        'expected.']
    path, telescope_config, layout_fname = \
            pyuvsim.simsetup.uvdata_to_telescope_config(uv0, herabeam_default,
            path_out=opath, return_names=True)

    uv0.integration_time[-1] += 2  # Test case of non-uniform integration times
    with pytest.warns(UserWarning) as warn:
        pyuvsim.simsetup.uvdata_to_config_file(uv0,
            telescope_config_name=os.path.join(path, telescope_config),
            layout_csv_name=os.path.join(path, layout_fname),
            path_out=opath
        )
    assert str(warn[0].message).startswith('The integration time is not constant. Using the shortest integration time')

    # From parameters, generate a uvdata object.
    param_dict = pyuvsim.simsetup._config_str_to_dict(os.path.join(opath, param_filename))

    orig_param_dict = copy.deepcopy(
        param_dict)  # The parameter dictionary gets modified in the function below.
    uv1, new_beam_list, new_beam_dict = pyuvsim.initialize_uvdata_from_params(param_dict)
    # Generate parameters from new uvfits and compare with old.
    path, telescope_config, layout_fname = \
        pyuvsim.simsetup.uvdata_to_telescope_config(uv1, herabeam_default,
            telescope_config_name=telescope_config,
            layout_csv_name=layout_fname,
            path_out=opath, return_names=True
        )
    pyuvsim.simsetup.uvdata_to_config_file(
        uv1, param_filename=second_param_filename,
        telescope_config_name=os.path.join(path, telescope_config),
        layout_csv_name=os.path.join(path, layout_fname),
        path_out=opath
    )

    del param_dict

    param_dict = pyuvsim.simsetup._config_str_to_dict(os.path.join(opath, second_param_filename))
    shutil.rmtree(opath)
    assert param_dict['obs_param_file'] == second_param_filename
    assert orig_param_dict['obs_param_file'] == param_filename
    orig_param_dict['obs_param_file'] = second_param_filename
    assert simtest.compare_dictionaries(param_dict, orig_param_dict)


def test_mock_catalogs():
    time = Time(2458098.27471265, scale='utc', format='jd')

    arrangements = ['off-zenith', 'zenith', 'cross', 'triangle', 'long-line', 'random', 'hera_text']

    cats = {}
    for arr in arrangements:
        # rseed is only used by the "random" mock catalog
        cat, mock_kwds = pyuvsim.simsetup.create_mock_catalog(time, arr, rseed=2458098)
        cats[arr] = cat

    # For each mock catalog, verify the Ra/Dec source positions against a text catalog.

    text_catalogs = {
        'cross': 'mock_cross_2458098.27471.txt',
        'hera_text': 'mock_hera_text_2458098.27471.txt',
        'long-line': 'mock_long-line_2458098.27471.txt',
        'off-zenith': 'mock_off-zenith_2458098.27471.txt',
        'triangle': 'mock_triangle_2458098.27471.txt',
        'random': 'mock_random_2458098.27471.txt',
        'zenith': 'mock_zenith_2458098.27471.txt'
    }
    simtest.assert_raises_message(
        KeyError, "Invalid mock catalog arrangement: invalid_catalog_name",
        pyuvsim.simsetup.create_mock_catalog, time, 'invalid_catalog_name'
    )

    for arr in arrangements:
        radec_catalog = pyradiosky.read_text_catalog(
            os.path.join(SIM_DATA_PATH, 'test_catalogs', text_catalogs[arr])
        )
        assert np.all(radec_catalog == cats[arr])

    cat, mock_kwds = pyuvsim.simsetup.create_mock_catalog(time, 'random', save=True)
    loc = eval(mock_kwds['array_location'])
    loc = EarthLocation.from_geodetic(loc[1], loc[0], loc[2])  # Lon, Lat, alt
    fname = 'mock_catalog_random.npz'
    alts_reload = np.load(fname)['alts']
    cat.update_positions(time, loc)
    alt, az = cat.alt_az
    assert np.all(alt > np.radians(30.))
    assert np.allclose(alts_reload, np.degrees(alt))
    os.remove(fname)


def test_keyword_param_loop():
    # Check that yaml/csv files made by intialize_uvdata_from_keywords will work
    # on their own.
    layout_fname = 'temp_layout_kwdloop.csv'
    obsparam_fname = 'temp_obsparam_kwdloop.yaml'
    path_out = simtest.TESTDATA_PATH
    antpos_enu = np.ones(30).reshape((10, 3))
    antnums = np.arange(10)
    antpos_d = dict(zip(antnums, antpos_enu))
    uvd = pyuvsim.simsetup.initialize_uvdata_from_keywords(
        array_layout=antpos_d,
        telescope_location=(-30.72152777777791, 21.428305555555557, 1073.0000000093132),
        telescope_name="HERA", Nfreqs=10, start_freq=1e8, bandwidth=1e8, Ntimes=60,
        integration_time=100.0, start_time=2458101.0, no_autos=True,
        path_out=path_out, antenna_layout_filepath=layout_fname, output_yaml_filename=obsparam_fname
    )

    uv2, _, _ = pyuvsim.simsetup.initialize_uvdata_from_params(
        os.path.join(path_out, obsparam_fname))

    uv2.extra_keywords = {}
    uvd.extra_keywords = {}  # These will not match

    assert uv2 == uvd


def test_multi_analytic_beams():
    # Test inline definitions of beam attributes.
    # eg. (in beam configuration file):
    #
    # beam_paths:
    #   0 : airy, diameter=14
    #   1 : airy, diameter=20
    #   2 : gaussian, sigma=0.5

    par_fname = os.path.join(simtest.TESTDATA_PATH, 'test_teleconfig.yaml')
    layout_fname = os.path.join(simtest.TESTDATA_PATH, 'test_layout_5ant.csv')

    telescope_location = (-30.72152777777791, 21.428305555555557, 1073.0000000093132)
    telescope_name = 'SKA'
    beam_specs = {0: 'airy, diameter=14', 1: 'airy, diameter=20', 2: 'gaussian, sigma=0.5'}
    expected = ['analytic_airy_diam_14', 'analytic_airy_diam_20', 'analytic_gaussian_sig_0.5']

    Nants = 5
    antenna_numbers = np.arange(Nants)
    antpos = np.zeros((Nants, 3))
    antpos[:, 0] = np.arange(Nants)
    names = antenna_numbers.astype(str)
    beam_ids = [0, 1, 2, 2, 0]
    pyuvsim.simsetup._write_layout_csv(layout_fname, antpos, names, antenna_numbers, beam_ids)

    # Write tele config to file.
    pdict = dict(telescope_location=str(telescope_location),
                 telescope_name=telescope_name, beam_paths=beam_specs)
    with open(par_fname, 'w') as yfile:
        yaml.dump(pdict, yfile, default_flow_style=False)

    # Check error condition:
    ambig_beam_specs = copy.copy(beam_specs)
    ambig_beam_specs[0] = 'airy, gaussian, diameter=14'
    bad_pdict = copy.copy(pdict)
    bad_pdict['beam_paths'] = ambig_beam_specs

    simtest.assert_raises_message(ValueError, "Ambiguous beam specification",
                                  pyuvsim.simsetup._construct_beam_list, beam_ids, bad_pdict)

    param_dict = {'telescope_config_name': par_fname, 'array_layout': layout_fname}

    pdict, beam_list, beam_dict = pyuvsim.simsetup.parse_telescope_params(
        param_dict, simtest.TESTDATA_PATH)

    for i, nm in enumerate(names):
        bid = beam_ids[i]
        assert beam_dict[nm] == bid
        assert beam_list[bid] == expected[bid]
