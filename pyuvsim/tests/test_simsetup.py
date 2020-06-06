# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import copy
import os
import shutil

import numpy as np
import pytest
import yaml
from astropy import units
from astropy.coordinates import Angle, SkyCoord, EarthLocation, Latitude, Longitude
from pyuvdata import UVBeam, UVData
import pyradiosky
from pyradiosky.data import DATA_PATH as SKY_DATA_PATH

import pyuvsim
import pyuvsim.tests as simtest
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH

from pyuvsim.astropy_interface import Time

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

    test_source = pyradiosky.SkyModel('src0', ra, dec, [1, 0, 0, 0], 'flat')

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
    test_source = pyradiosky.SkyModel('src0', ra, dec, [1.0, 0, 0, 0], 'flat')

    cat, mock_keywords = pyuvsim.create_mock_catalog(time, arrangement='off-zenith',
                                                     alt=src_alt.deg)

    assert cat == test_source


def test_catalog_from_params():
    # Pass in parameter dictionary as dict
    hera_uv = UVData()
    with pytest.warns(
        UserWarning,
        match='Telescope 28m_triangle_10time_10chan.yaml is not in known_telescopes.'
    ):
        hera_uv.read_uvfits(triangle_uvfits_file)

    source_dict = {}
    with pytest.raises(KeyError, match='No catalog defined.'):
        pyuvsim.simsetup.initialize_catalog_from_params({'sources': source_dict})

    arrloc = '{:.7f},{:.7f},{:.7f}'.format(*hera_uv.telescope_location_lat_lon_alt_degrees)
    source_dict = {
        'catalog': 'mock',
        'mock_arrangement': 'zenith',
        'Nsrcs': 5,
        'time': hera_uv.time_array[0]
    }
    with pytest.warns(
        UserWarning,
        match="No array_location specified. Defaulting to the HERA site."
    ):
        pyuvsim.simsetup.initialize_catalog_from_params({'sources': source_dict})

    catalog_uv, srclistname = pyuvsim.simsetup.initialize_catalog_from_params(
        {'sources': source_dict}, hera_uv
    )
    catalog_uv = catalog_uv.get_skymodel()
    source_dict['array_location'] = arrloc
    del source_dict['time']

    with pytest.raises(TypeError, match="input_uv must be UVData object"):
        pyuvsim.simsetup.initialize_catalog_from_params({'sources': source_dict},
                                                        input_uv='not_uvdata')

    with pytest.raises(ValueError, match="input_uv must be supplied if using mock catalog"):
        pyuvsim.simsetup.initialize_catalog_from_params({'sources': source_dict})

    with pytest.warns(
        UserWarning,
        match="Warning: No julian date given for mock catalog. Defaulting to first time step."
    ):
        catalog_str, srclistname2 = pyuvsim.simsetup.initialize_catalog_from_params(
            {'sources': source_dict}, hera_uv
        )

    catalog_str = catalog_str.get_skymodel()
    assert np.all(catalog_str == catalog_uv)


def test_vot_catalog():
    vot_param_filename = os.path.join(SIM_DATA_PATH, 'test_config', 'param_1time_1src_testvot.yaml')
    vot_catalog = (
        pyuvsim.simsetup.initialize_catalog_from_params(vot_param_filename)[0]
    ).get_skymodel()

    txt_param_filename = os.path.join(SIM_DATA_PATH, 'test_config', 'param_1time_1src_testcat.yaml')
    txt_catalog = (
        pyuvsim.simsetup.initialize_catalog_from_params(txt_param_filename)[0]
    ).get_skymodel()

    assert vot_catalog == txt_catalog


def test_gleam_catalog():
    gleam_param_filename = os.path.join(
        SIM_DATA_PATH, 'test_config', 'param_1time_1src_testgleam.yaml'
    )
    with pytest.warns(UserWarning, match="No spectral_type specified for GLEAM, using 'flat'."):
        gleam_catalog = (
            pyuvsim.simsetup.initialize_catalog_from_params(gleam_param_filename)[0]
        ).get_skymodel()

    # flux cuts applied
    assert gleam_catalog.Ncomponents == 23

    # no cuts
    with open(gleam_param_filename, 'r') as pfile:
        param_dict = yaml.safe_load(pfile)
    param_dict['config_path'] = os.path.dirname(gleam_param_filename)
    param_dict["sources"].pop("min_flux")
    param_dict["sources"].pop("max_flux")

    gleam_catalog = (
        pyuvsim.simsetup.initialize_catalog_from_params(param_dict)[0]
    ).get_skymodel()
    assert gleam_catalog.Ncomponents == 50


@pytest.mark.parametrize(
    "spectral_type",
    ["flat", "subband", "spectral_index"])
def test_gleam_catalog_spectral_type(spectral_type):
    gleam_param_filename = os.path.join(
        SIM_DATA_PATH, 'test_config', 'param_1time_1src_testgleam.yaml'
    )
    with open(gleam_param_filename, 'r') as pfile:
        param_dict = yaml.safe_load(pfile)
    param_dict['config_path'] = os.path.dirname(gleam_param_filename)
    param_dict["sources"].pop("min_flux")
    param_dict["sources"].pop("max_flux")
    param_dict["sources"]["spectral_type"] = spectral_type

    gleam_catalog = pyuvsim.simsetup.initialize_catalog_from_params(param_dict)[0]
    assert gleam_catalog.spectral_type == spectral_type
    assert gleam_catalog.Ncomponents == 50


@pytest.mark.parametrize(
    ("key_pop", "message"),
    [("ra_column", "No RA column name specified"),
     ("dec_column", "No Dec column name specified")])
def test_vot_catalog_warns(key_pop, message):
    vot_param_filename = os.path.join(SIM_DATA_PATH, 'test_config', 'param_1time_1src_testvot.yaml')

    vot_catalog = (
        pyuvsim.simsetup.initialize_catalog_from_params(vot_param_filename)[0]
    ).get_skymodel()

    with open(vot_param_filename, 'r') as pfile:
        param_dict = yaml.safe_load(pfile)
    param_dict['config_path'] = os.path.dirname(vot_param_filename)
    param_dict["sources"].pop(key_pop)

    with pytest.warns(UserWarning, match=message):
        vot_catalog2 = (
            pyuvsim.simsetup.initialize_catalog_from_params(param_dict)[0]
        ).get_skymodel()

    assert vot_catalog == vot_catalog2


@pytest.mark.parametrize(
    ("key_pop", "message"),
    [("table_name", "No VO table name specified"),
     ("id_column", "No ID column name specified"),
     ("flux_columns", "No Flux column names specified")])
def test_vot_catalog_error(key_pop, message):
    vot_param_filename = os.path.join(SIM_DATA_PATH, 'test_config', 'param_1time_1src_testvot.yaml')

    with open(vot_param_filename, 'r') as pfile:
        param_dict = yaml.safe_load(pfile)
    param_dict['config_path'] = os.path.dirname(vot_param_filename)
    param_dict["sources"].pop(key_pop)

    with pytest.raises(ValueError, match=message):
        pyuvsim.simsetup.initialize_catalog_from_params(param_dict)[0]


# parametrize will loop over all the give values
@pytest.mark.parametrize("config_num", [0, 2])
def test_param_reader(config_num):
    pytest.importorskip('mpi4py')
    # Reading in various configuration files

    param_filename = param_filenames[config_num]
    hera_uv = UVData()
    with pytest.warns(UserWarning) as warn:
        hera_uv.read_uvfits(triangle_uvfits_file)
    assert str(warn.pop().message).startswith('Telescope 28m_triangle_10time_10chan.yaml '
                                              'is not in known_telescopes.')

    hera_uv.telescope_name = 'HERA'
    if config_num == 5:
        hera_uv.select(bls=[(0, 1), (1, 2)])

    time = Time(hera_uv.time_array[0], scale='utc', format='jd')
    sources, _ = pyuvsim.create_mock_catalog(time, arrangement='zenith', return_data=True)

    beam0 = UVBeam()
    beam0.read_beamfits(herabeam_default)
    beam1 = pyuvsim.AnalyticBeam('uniform')
    beam2 = pyuvsim.AnalyticBeam('gaussian', sigma=0.02)
    beam3 = pyuvsim.AnalyticBeam('airy', diameter=14.6)
    beam_list = pyuvsim.BeamList([beam0, beam1, beam2, beam3])

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
        with pytest.raises(ValueError, match="nonexistent_directory is not a directory"):
            pyuvsim.simsetup.initialize_uvdata_from_params(params_bad)

        params_bad['config_path'] = os.path.join(SIM_DATA_PATH, "test_config")
        params_bad['telescope']['array_layout'] = 'nonexistent_file'
        with pytest.raises(ValueError, match="nonexistent_file from yaml does not exist"):
            pyuvsim.simsetup.initialize_uvdata_from_params(params_bad)

        params_bad['telescope']['telescope_config_name'] = 'nonexistent_file'
        with pytest.raises(ValueError, match="telescope_config_name file from yaml does not exist"):
            pyuvsim.simsetup.initialize_uvdata_from_params(params_bad)

        # Missing beam keywords
        params_bad = copy.deepcopy(bak_params_bad)

        params_bad['config_path'] = os.path.join(SIM_DATA_PATH, "test_config")

        params_bad = copy.deepcopy(bak_params_bad)
        params_bad['telescope']['telescope_config_name'] = os.path.join(
            SIM_DATA_PATH, 'test_config', '28m_triangle_10time_10chan_gaussnoshape.yaml'
        )
        with pytest.raises(KeyError,
                           match="Missing shape parameter for gaussian beam"):
            pyuvsim.simsetup.initialize_uvdata_from_params(params_bad)

        params_bad['telescope']['telescope_config_name'] = os.path.join(
            SIM_DATA_PATH, 'test_config', '28m_triangle_10time_10chan_nodiameter.yaml'
        )
        with pytest.raises(KeyError, match="Missing diameter for airy beam."):
            pyuvsim.simsetup.initialize_uvdata_from_params(params_bad)

        params_bad['telescope']['telescope_config_name'] = os.path.join(
            SIM_DATA_PATH, 'test_config', '28m_triangle_10time_10chan_nofile.yaml'
        )
        with pytest.raises(ValueError, match="Undefined beam model"):
            pyuvsim.simsetup.initialize_uvdata_from_params(params_bad)

    # Check default configuration
    uv_obj, new_beam_list, new_beam_dict = pyuvsim.initialize_uvdata_from_params(param_filename)
    new_beam_list.set_obj_mode()

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
    tdict = {'array_layout': os.path.join(SIM_DATA_PATH, 'test_layout_6ant.csv')}
    tel_error = ('If telescope_config_name not provided in `telescope` obsparam section, '
                 'you must provide telescope_location')

    with pytest.raises(KeyError, match=tel_error):
        pyuvsim.simsetup.parse_telescope_params(tdict)

    tdict['telescope_location'] = '(-30.72152777777791, 21.428305555555557, 1073.0000000093132)'
    tel_error = 'If telescope_config_name not provided in `telescope` obsparam section, ' \
                'you must provide telescope_name'
    with pytest.raises(KeyError, match=tel_error):
        pyuvsim.simsetup.parse_telescope_params(tdict)

    tdict['telescope_name'] = 'tele'
    tpars, blist, bdict = pyuvsim.simsetup.parse_telescope_params(tdict)

    assert tpars['Nants_data'] == 6
    assert len(blist) == 0

    tdict.pop('array_layout')
    with pytest.raises(KeyError, match="array_layout must be provided."):
        pyuvsim.simsetup.parse_telescope_params(tdict)


def test_freq_parser():
    """
    Check all valid input parameter cases for frequencies.
    """

    fdict_base = {
        "Nfreqs": 10,
        "channel_width": 0.5,
        "start_freq": 0.0,
        "end_freq": 4.5,
        "bandwidth": 5.0
    }

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
        with pytest.raises(ValueError, match=err_mess[ei]):
            pyuvsim.parse_frequency_params(subdict)

    subdict = {'freq_array': freq_array[0]}
    with pytest.raises(ValueError,
                       match='Channel width must be specified if freq_arr has length 1'):
        pyuvsim.parse_frequency_params(subdict)

    subdict = {'freq_array': freq_array[[0, 1, 4, 8]]}
    with pytest.raises(ValueError, match='Spacing in frequency array is uneven.'):
        pyuvsim.parse_frequency_params(subdict)

    subdict = {'channel_width': 3.14, 'start_freq': 1.0, 'end_freq': 8.3}
    with pytest.raises(ValueError,
                       match='end_freq - start_freq must be evenly divisible by channel_width'):
        pyuvsim.parse_frequency_params(subdict)

    subdict = fdict_base.copy()
    subdict['Nfreqs'] = 7
    del subdict['freq_array']
    with pytest.raises(ValueError,
                       match='Frequency array spacings are not equal to channel width.'):
        pyuvsim.parse_frequency_params(subdict)


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
        with pytest.raises(ValueError, match=err_mess[ei]):
            pyuvsim.parse_time_params(subdict)

    subdict = {'integration_time': 3.14, 'start_time': 10000.0, 'end_time': 80000.3, 'Ntimes': 30}
    with pytest.raises(ValueError,
                       match='Calculated time array is not consistent with set integration_time'):
        pyuvsim.parse_time_params(subdict)

    subdict = tdict_base.copy()
    subdict['Ntimes'] = 7
    del subdict['time_array']
    with pytest.raises(ValueError,
                       match='Calculated time array is not consistent with set integration_time'):
        pyuvsim.parse_time_params(subdict)


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
    base_kwargs = {
        "antenna_layout_filepath": os.path.join(SIM_DATA_PATH,
                                                "test_config/triangle_bl_layout.csv"),
        "telescope_location": (-30.72152777777791, 21.428305555555557, 1073.0000000093132),
        "telescope_name": "HERA",
        "Nfreqs": 10,
        "start_freq": 1e8,
        "bandwidth": 1e8,
        "Ntimes": 60,
        "integration_time": 100.0,
        "start_time": 2458101.0,
        "polarization_array": ['xx'],
        "no_autos": True,
        "write_files": False,
        "run_check": True
    }

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

    path, telescope_config, layout_fname = \
        pyuvsim.simsetup.uvdata_to_telescope_config(uv0, herabeam_default,
                                                    path_out=opath, return_names=True)

    uv0.integration_time[-1] += 2  # Test case of non-uniform integration times
    with pytest.warns(UserWarning) as warn:
        pyuvsim.simsetup.uvdata_to_config_file(
            uv0,
            telescope_config_name=os.path.join(path, telescope_config),
            layout_csv_name=os.path.join(path, layout_fname),
            path_out=opath
        )
    assert str(warn[0].message).startswith('The integration time is not constant. '
                                           'Using the shortest integration time')

    # From parameters, generate a uvdata object.
    param_dict = pyuvsim.simsetup._config_str_to_dict(os.path.join(opath, param_filename))

    orig_param_dict = copy.deepcopy(
        param_dict)  # The parameter dictionary gets modified in the function below.
    uv1, new_beam_list, new_beam_dict = pyuvsim.initialize_uvdata_from_params(param_dict)
    # Generate parameters from new uvfits and compare with old.
    path, telescope_config, layout_fname = \
        pyuvsim.simsetup.uvdata_to_telescope_config(
            uv1, herabeam_default,
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
    with pytest.raises(KeyError, match="Invalid mock catalog arrangement: invalid_catalog_name"):
        pyuvsim.create_mock_catalog(time, 'invalid_catalog_name')

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
    beam_specs = {0: {'type': 'airy', 'diameter': 14},
                  1: {'type': 'airy', 'diameter': 20},
                  2: {'type': 'gaussian', 'sigma': 0.5}}
    expected = ['analytic_airy_diam=14', 'analytic_airy_diam=20', 'analytic_gaussian_sig=0.5']

    Nants = 5
    antenna_numbers = np.arange(Nants)
    antpos = np.zeros((Nants, 3))
    antpos[:, 0] = np.arange(Nants)
    names = antenna_numbers.astype(str)
    beam_ids = [0, 1, 2, 2, 0]
    pyuvsim.simsetup._write_layout_csv(layout_fname, antpos, names, antenna_numbers, beam_ids)

    # Write tele config to file.
    pdict = {
        "telescope_location": str(telescope_location),
        "telescope_name": telescope_name,
        "beam_paths": beam_specs
    }
    with open(par_fname, 'w') as yfile:
        yaml.dump(pdict, yfile, default_flow_style=False)

    param_dict = {'telescope_config_name': par_fname, 'array_layout': layout_fname}

    pdict, beam_list, beam_dict = pyuvsim.simsetup.parse_telescope_params(
        param_dict, simtest.TESTDATA_PATH)

    for i, nm in enumerate(names):
        bid = beam_ids[i]
        assert beam_dict[nm] == bid
        assert beam_list[bid] == expected[bid]


def test_direct_fname():
    shutil.copyfile(
        os.path.join(SIM_DATA_PATH, "test_config", "28m_triangle_10time_10chan.yaml"),
        "28m_triangle_10time_10chan.yaml"
    )
    shutil.copyfile(
        os.path.join(SIM_DATA_PATH, "test_config", "param_100times_1.5days_triangle.yaml"),
        "param_100times_1.5days_triangle.yaml"
    )
    shutil.copyfile(
        os.path.join(SIM_DATA_PATH, "test_config", "triangle_bl_layout.csv"),
        "triangle_bl_layout.csv"
    )

    # This should now run without errors
    pyuvsim.simsetup.initialize_uvdata_from_params("param_100times_1.5days_triangle.yaml")

    os.remove("28m_triangle_10time_10chan.yaml")
    os.remove("param_100times_1.5days_triangle.yaml")
    os.remove("triangle_bl_layout.csv")


def test_beamopts_init():
    # Check that spline_interp_opts is passed along correctly to BeamList
    telescope_config_name = os.path.join(SIM_DATA_PATH, 'mwa128_config.yaml')
    with open(telescope_config_name, 'r') as yf:
        telconfig = yaml.safe_load(yf)
    telconfig['spline_interp_opts'] = {'kx' : 2, 'ky' : 2}
    beam_list = pyuvsim.simsetup._construct_beam_list(np.arange(1), telconfig)
    assert beam_list.spline_interp_opts is not None


def test_moon_lsts():
    # Check that setting lsts for a Moon simulation works as expected.
    pytest.importorskip('lunarsky')

    param_filename = os.path.join(SIM_DATA_PATH, 'test_config', 'obsparam_tranquility_hex.yaml')
    param_dict = pyuvsim.simsetup._config_str_to_dict(param_filename)
    uv_obj, beam_list, beam_dict = pyuvsim.initialize_uvdata_from_params(param_dict)
    assert 'world' in uv_obj.extra_keywords.keys()
    assert uv_obj.extra_keywords['world'] == 'moon'

    # Check ordering on lsts -- unique lsts should correspond with unique times.
    lsts, lst_inds = np.unique(uv_obj.lst_array, return_inverse=True)
    times, time_inds = np.unique(uv_obj.time_array, return_inverse=True)
    assert np.all(lst_inds == time_inds)

    # Confirm that the change in LST makes sense for a lunar month.
    dlst = np.degrees(lsts[1] - lsts[0]) / 15 * 3600  # Hours per step.
    tstep_per_lunar_month = 28 * 24 * 3600 / uv_obj.integration_time[0]
    sec_per_tstep = 24 * 3600 / tstep_per_lunar_month

    assert np.isclose(dlst, sec_per_tstep, rtol=1e-2)

    # Unset the lst array and confirm that the call from _complete_uvdata returns the same.
    backup_lst_array = uv_obj.lst_array.copy()
    uv_obj.lst_array = None

    new_obj = pyuvsim.simsetup._complete_uvdata(uv_obj)

    assert np.allclose(new_obj.lst_array, backup_lst_array)
    assert new_obj.check()


@pytest.mark.filterwarnings("ignore:The _ra parameters are not")
@pytest.mark.filterwarnings("ignore:The _dec parameters are not")
@pytest.mark.filterwarnings("ignore:Future equality does not pass")
def test_mock_catalog_moon():
    # A mock catalog made with a MoonLocation.
    pytest.importorskip('lunarsky')
    import lunarsky
    from pyuvsim.astropy_interface import Time

    time = Time.now()
    loc = lunarsky.MoonLocation.from_selenodetic(24.433333333, 0.687500000)
    mmock, mkwds = pyuvsim.simsetup.create_mock_catalog(time, 'hera_text', array_location=loc)
    eloc = EarthLocation.from_geodetic(24.433, 0.6875)
    emock, ekwds = pyuvsim.simsetup.create_mock_catalog(time, 'hera_text', array_location=eloc)

    assert mkwds['world'] == 'moon'
    assert ekwds['world'] == 'earth'

    # Simple check that the given lat/lon were interpreted differently in each call.
    assert mmock != emock


@pytest.fixture
def cat_with_some_pols():
    # Mock catalog with a couple sources polarized.
    Nsrcs = 30
    Nfreqs = 10
    freqs = np.linspace(100, 130, Nfreqs) * 1e6 * units.Hz

    pol_inds = range(12, 15)
    stokes = np.zeros((4, Nfreqs, Nsrcs))
    stokes[0, :, :] = 1.0

    stokes[1, :, pol_inds] = 0.2
    stokes[2, :, pol_inds] = 1.2
    stokes[3, :, pol_inds] = 0.3

    ra = Longitude(np.linspace(0, 2 * np.pi, Nsrcs), 'rad')
    dec = Latitude(np.linspace(-np.pi / 2, np.pi / 3, Nsrcs), 'rad')

    sky = pyradiosky.SkyModel(
        name=np.arange(Nsrcs).astype(str),
        ra=ra,
        dec=dec,
        stokes=stokes,
        spectral_type='full',
        freq_array=freqs
    )

    return sky


@pytest.mark.skipif("units.Quantity not in pyradiosky.SkyModel()._stokes.expected_type")
def test_skymodeldata_with_quantity_stokes(cat_with_some_pols):
    # Support for upcoming pyradiosky change setting SkyModel.stokes
    # to an astropy Quantity.
    sky = cat_with_some_pols
    sky.stokes *= units.Jy

    smd = pyuvsim.simsetup.SkyModelData(sky)
    assert np.all(sky.stokes.to('Jy').value[0] == smd.stokes_I)

    sky2 = smd.get_skymodel()
    assert sky2 == sky


@pytest.mark.parametrize('component_type', ['point', 'healpix'])
def test_skymodeldata(component_type, cat_with_some_pols):
    # Test that SkyModelData class can properly recreate a SkyModel and subselect.
    if component_type == 'point':
        sky = cat_with_some_pols
    else:
        pytest.importorskip('astropy-healpix')
        path = os.path.join(SKY_DATA_PATH, 'healpix_disk.hdf5')
        sky = pyradiosky.SkyModel()
        sky.read_healpix_hdf5(path)
    smd = pyuvsim.simsetup.SkyModelData(sky)

    assert (smd.ra == sky.ra.deg).all()
    assert (smd.dec == sky.dec.deg).all()
    assert (smd.stokes_I == sky.stokes[0]).all()

    if smd.polarized is not None:
        assert (smd.stokes_Q == sky.stokes[..., smd.polarized][1]).all()
        assert (smd.stokes_U == sky.stokes[..., smd.polarized][2]).all()
        assert (smd.stokes_V == sky.stokes[..., smd.polarized][3]).all()

    # Make skymodel from SkyModelData.
    sky1 = smd.get_skymodel()

    assert sky1 == sky

    # Now try with subselection:
    sky1_sub = smd.get_skymodel(range(8, 13))

    assert sky1.check()
    assert sky1_sub.check()
    assert sky1_sub.Ncomponents == 5
    if smd.polarized is not None:
        assert sky1_sub._n_polarized == 1


@pytest.mark.parametrize('inds', [range(30), range(5), range(9, 14)])
def test_skymodeldata_pol_select(inds, cat_with_some_pols):
    # When running SkyModelData.subselect, confirm that the
    # polarization array and Q, U, V are properly selected.

    smd = pyuvsim.simsetup.SkyModelData(cat_with_some_pols)
    sub_smd = smd.subselect(inds)

    test_q = np.zeros((smd.Nfreqs, smd.Ncomponents))
    temp = np.zeros((sub_smd.Nfreqs, sub_smd.Ncomponents))
    temp[..., sub_smd.polarized] = sub_smd.stokes_Q
    test_q[..., inds] = temp[()]

    full_q = np.zeros_like(test_q)
    full_q[..., smd.polarized] = smd.stokes_Q

    assert np.all(full_q[..., inds] == test_q[..., inds])


@pytest.mark.parametrize('inds', [range(30), range(5)])
def test_skymodeldata_attr_bases(inds, cat_with_some_pols):
    # Check that downselecting doesn't copy length Ncompnent arrays.

    smd = pyuvsim.simsetup.SkyModelData(cat_with_some_pols)
    smd_copy = smd.subselect(inds)
    assert smd_copy.ra.base is smd.ra.base
    assert smd_copy.dec.base is smd.dec.base
    assert smd_copy.stokes_I.base is smd.stokes_I.base
