# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import os
import copy
import shutil

import numpy as np
import pytest
import yaml
from astropy import units
from astropy.coordinates import Angle, SkyCoord, EarthLocation, Latitude, Longitude
from pyuvdata import UVBeam, UVData
import pyuvdata.tests as uvtest
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


def test_mock_catalog_zenith_source(hera_loc):
    time = Time(2457458.65410, scale='utc', format='jd')

    array_location = hera_loc

    source_coord = SkyCoord(
        alt=Angle(90 * units.deg), az=Angle(0 * units.deg),
        obstime=time, frame='altaz', location=array_location
    )
    icrs_coord = source_coord.transform_to('icrs')

    ra = icrs_coord.ra
    dec = icrs_coord.dec

    test_source = pyradiosky.SkyModel('src0', ra, dec, units.Quantity([1, 0, 0, 0], 'Jy'), 'flat')

    cat, mock_keywords = pyuvsim.create_mock_catalog(time, arrangement='zenith')

    assert cat == test_source


def test_mock_catalog_off_zenith_source(hera_loc):
    src_az = Angle('90.0d')
    src_alt = Angle('85.0d')

    time = Time(2457458.65410, scale='utc', format='jd')

    array_location = hera_loc

    source_coord = SkyCoord(
        alt=src_alt, az=src_az, obstime=time, frame='altaz', location=array_location
    )
    icrs_coord = source_coord.transform_to('icrs')

    ra = icrs_coord.ra
    dec = icrs_coord.dec
    test_source = pyradiosky.SkyModel('src0', ra, dec, units.Quantity([1.0, 0, 0, 0], "Jy"), 'flat')

    cat, mock_keywords = pyuvsim.create_mock_catalog(time, arrangement='off-zenith',
                                                     alt=src_alt.deg)

    assert cat == test_source


def test_mock_diffuse_maps_errors():
    analytic_diffuse = pyuvsim.simsetup.analytic_diffuse
    astropy_healpix = pyuvsim.simsetup.astropy_healpix
    if (analytic_diffuse is not None) and (astropy_healpix is not None):

        # Error cases:
        with pytest.raises(ValueError, match="Diffuse arrangement selected"):
            pyuvsim.simsetup.create_mock_catalog(Time.now(), arrangement='diffuse')

        with uvtest.check_warnings(UserWarning, match="No nside chosen"):
            pyuvsim.simsetup.create_mock_catalog(
                Time.now(), arrangement='diffuse', diffuse_model='monopole'
            )

    else:
        with pytest.raises(ValueError, match="analytic_diffuse and astropy_healpix"):
            pyuvsim.simsetup.create_mock_catalog(Time.now(), arrangement='diffuse')


@pytest.mark.filterwarnings("ignore:Input ra and dec parameters are being used instead")
@pytest.mark.parametrize('model', [
                         ['monopole', {}],
                         ['gauss', {'a': 0.05}],
                         ['polydome', {'n': 4}]
                         ])
@pytest.mark.parametrize('location', ['earth', 'moon'])
def test_mock_diffuse_maps(model, hera_loc, apollo_loc, location):
    analytic_diffuse = pytest.importorskip('analytic_diffuse')
    pytest.importorskip('astropy_healpix')
    if location == 'earth':
        loc = hera_loc
    else:
        loc = apollo_loc
    modname, modkwargs = model
    map_nside = 128
    t0 = Time.now()
    cat, kwds = pyuvsim.simsetup.create_mock_catalog(
        t0, arrangement='diffuse', array_location=loc,
        diffuse_model=modname, map_nside=map_nside,
        diffuse_params=modkwargs
    )

    cat.update_positions(t0, loc)

    modfunc = analytic_diffuse.get_model(modname)
    alt, az = cat.alt_az
    za = np.pi / 2 - alt

    vals = modfunc(az, za, **modkwargs)

    assert cat.nside == map_nside
    assert np.allclose(cat.stokes[0].to_value("K"), vals)


@pytest.mark.filterwarnings("ignore:LST values stored in this file are not self-consistent")
def test_catalog_from_params():
    # Pass in parameter dictionary as dict
    hera_uv = UVData()
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
    with uvtest.check_warnings(
        [UserWarning, PendingDeprecationWarning, DeprecationWarning],
        match=[
            "No array_location specified. Defaulting to the HERA site.",
            "initialize_catalog_from_params will return a SkyModel instance",
            "recarray flux columns will no longer be labeled",
        ]
    ):
        pyuvsim.simsetup.initialize_catalog_from_params({'sources': source_dict})

    catalog_uv, srclistname = pyuvsim.simsetup.initialize_catalog_from_params(
        {'sources': source_dict}, hera_uv, return_recarray=False
    )
    source_dict['array_location'] = arrloc
    del source_dict['time']

    with pytest.raises(TypeError, match="input_uv must be UVData object"):
        pyuvsim.simsetup.initialize_catalog_from_params({'sources': source_dict},
                                                        input_uv='not_uvdata')

    with pytest.raises(ValueError, match="input_uv must be supplied if using mock catalog"):
        pyuvsim.simsetup.initialize_catalog_from_params({'sources': source_dict})

    with uvtest.check_warnings(
        UserWarning,
        match="No julian date given for mock catalog. Defaulting to first time step."
    ):
        catalog_str, srclistname2 = pyuvsim.simsetup.initialize_catalog_from_params(
            {'sources': source_dict}, hera_uv, return_recarray=False
        )
    assert np.all(catalog_str == catalog_uv)


def test_vot_catalog():
    vot_param_filename = os.path.join(SIM_DATA_PATH, 'test_config', 'param_1time_1src_testvot.yaml')
    vot_catalog = (
        pyuvsim.simsetup.initialize_catalog_from_params(
            vot_param_filename, return_recarray=False
        )[0]
    )

    txt_param_filename = os.path.join(SIM_DATA_PATH, 'test_config', 'param_1time_1src_testcat.yaml')
    txt_catalog = (
        pyuvsim.simsetup.initialize_catalog_from_params(
            txt_param_filename, return_recarray=False
        )[0]
    )

    assert vot_catalog == txt_catalog


def test_gleam_catalog():
    gleam_param_filename = os.path.join(
        SIM_DATA_PATH, 'test_config', 'param_1time_1src_testgleam.yaml'
    )
    with uvtest.check_warnings(
        UserWarning, match="No spectral_type specified for GLEAM, using 'flat'."
    ):
        gleam_catalog = (
            pyuvsim.simsetup.initialize_catalog_from_params(
                gleam_param_filename, return_recarray=False
            )[0]
        )

    # flux cuts applied
    assert gleam_catalog.Ncomponents == 23

    # no cuts
    with open(gleam_param_filename, 'r') as pfile:
        param_dict = yaml.safe_load(pfile)
    param_dict['config_path'] = os.path.dirname(gleam_param_filename)
    param_dict["sources"].pop("min_flux")
    param_dict["sources"].pop("max_flux")

    with uvtest.check_warnings(
        UserWarning, match="No spectral_type specified for GLEAM, using 'flat'."
    ):
        gleam_catalog = (
            pyuvsim.simsetup.initialize_catalog_from_params(param_dict, return_recarray=False)[0]
        )
    assert gleam_catalog.Ncomponents == 50


def test_skyh5_catalog(tmp_path):
    gleam_filename = os.path.join(
        SIM_DATA_PATH, 'gleam_50srcs.vot'
    )
    skyobj = pyradiosky.SkyModel.from_gleam_catalog(gleam_filename)
    assert skyobj.Ncomponents == 50

    skyh5_file = os.path.join(tmp_path, 'gleam.skyh5')
    skyobj.write_skyh5(skyh5_file, clobber=True)

    starting_param_filename = os.path.join(
        SIM_DATA_PATH, 'test_config', 'param_1time_1src_testcat.yaml'
    )

    param_filename = os.path.join(
        tmp_path, 'param_test_skyh5_gleam.yaml'
    )

    with open(starting_param_filename, 'r') as yf:
        param_dict = yaml.safe_load(yf)

    param_dict['sources']['catalog'] = skyh5_file
    with open(param_filename, 'w') as yfile:
        yaml.dump(param_dict, yfile, default_flow_style=False)

    skyh5_catalog = (
        pyuvsim.simsetup.initialize_catalog_from_params(
            param_filename, return_recarray=False
        )[0]
    )
    assert skyh5_catalog.Ncomponents == 50

    # test works with `hdf5` extension:
    hdf5_file = os.path.join(tmp_path, 'gleam.hdf5')
    skyobj.write_skyh5(hdf5_file, clobber=True)

    param_dict['sources']['catalog'] = hdf5_file
    with open(param_filename, 'w') as yfile:
        yaml.dump(param_dict, yfile, default_flow_style=False)

    hdf5_catalog = (
        pyuvsim.simsetup.initialize_catalog_from_params(
            param_filename, return_recarray=False
        )[0]
    )
    assert hdf5_catalog == skyh5_catalog

    # test error with unknown extension
    h5_file = os.path.join(tmp_path, 'gleam.h5')
    skyobj.write_skyh5(hdf5_file, clobber=True)

    param_dict['sources']['catalog'] = h5_file
    with open(param_filename, 'w') as yfile:
        yaml.dump(param_dict, yfile, default_flow_style=False)

    with pytest.raises(ValueError, match="The file extension is not recognized as a valid one"):
        pyuvsim.simsetup.initialize_catalog_from_params(
            param_filename, return_recarray=False
        )


@pytest.mark.filterwarnings("ignore:Input ra and dec parameters are being used instead")
def test_healpix_catalog():
    pytest.importorskip('astropy_healpix')
    path = os.path.join(SKY_DATA_PATH, 'healpix_disk.skyh5')
    sky = pyradiosky.SkyModel()
    sky.read_skyh5(path)

    params = {'sources': {'catalog': path}}
    hpx_sky = pyuvsim.simsetup.initialize_catalog_from_params(
        params, return_recarray=False
    )[0]
    assert hpx_sky == sky


@pytest.mark.filterwarnings("ignore:This method reads an old 'healvis' style healpix HDF5 file")
def test_healpix_hdf5_catalog():
    """Test that an old healvis style hdf5 file works.

    Remove this when we require a newer version of pyradiosky.
    """
    pytest.importorskip('astropy_healpix')
    path = os.path.join(SKY_DATA_PATH, 'healpix_disk.hdf5')
    sky = pyradiosky.SkyModel()
    sky.read_healpix_hdf5(path)

    params = {'sources': {'catalog': path}}
    hpx_sky = pyuvsim.simsetup.initialize_catalog_from_params(
        params, return_recarray=False
    )[0]
    assert hpx_sky == sky


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

    gleam_catalog = pyuvsim.simsetup.initialize_catalog_from_params(
        param_dict, return_recarray=False
    )[0]
    assert gleam_catalog.spectral_type == spectral_type
    assert gleam_catalog.Ncomponents == 50


@pytest.mark.parametrize(
    ("key_pop", "message"),
    [("ra_column", "No ra_column specified"),
     ("dec_column", "No dec_column specified")])
def test_vot_catalog_warns(key_pop, message):
    vot_param_filename = os.path.join(SIM_DATA_PATH, 'test_config', 'param_1time_1src_testvot.yaml')

    vot_catalog = (
        pyuvsim.simsetup.initialize_catalog_from_params(
            vot_param_filename, return_recarray=False
        )[0]
    )
    with open(vot_param_filename, 'r') as pfile:
        param_dict = yaml.safe_load(pfile)
    param_dict['config_path'] = os.path.dirname(vot_param_filename)
    param_dict["sources"].pop(key_pop)

    with uvtest.check_warnings(UserWarning, match=message):
        vot_catalog2 = (
            pyuvsim.simsetup.initialize_catalog_from_params(param_dict, return_recarray=False)[0]
        )

    assert vot_catalog == vot_catalog2


@pytest.mark.parametrize(
    ("key_pop", "message"),
    [("table_name", "No table_name specified"),
     ("id_column", "No id_column specified"),
     ("flux_columns", "No flux_columns specified")])
def test_vot_catalog_error(key_pop, message):
    vot_param_filename = os.path.join(SIM_DATA_PATH, 'test_config', 'param_1time_1src_testvot.yaml')

    with open(vot_param_filename, 'r') as pfile:
        param_dict = yaml.safe_load(pfile)
    param_dict['config_path'] = os.path.dirname(vot_param_filename)
    param_dict["sources"].pop(key_pop)

    with pytest.raises(ValueError, match=message):
        pyuvsim.simsetup.initialize_catalog_from_params(param_dict, return_recarray=False)[0]


@pytest.mark.filterwarnings("ignore:LST values stored in this file are not self-consistent")
def test_param_reader():
    pytest.importorskip('mpi4py')
    # Reading in various configuration files

    param_filename = os.path.join(SIM_DATA_PATH, "test_config", "param_10time_10chan_0.yaml")
    hera_uv = UVData()
    hera_uv.read_uvfits(triangle_uvfits_file)
    # set missing x_orientation
    hera_uv.x_orientation = "east"

    hera_uv.unphase_to_drift()
    hera_uv.telescope_name = 'HERA'

    time = Time(hera_uv.time_array[0], scale='utc', format='jd')
    sources, _ = pyuvsim.create_mock_catalog(time, arrangement='zenith', return_data=True)

    beam0 = UVBeam()
    beam0.read_beamfits(herabeam_default)
    beam0.extra_keywords['beam_path'] = herabeam_default
    beam1 = pyuvsim.AnalyticBeam('uniform')
    beam2 = pyuvsim.AnalyticBeam('gaussian', sigma=0.02)
    beam3 = pyuvsim.AnalyticBeam('airy', diameter=14.6)
    beam_list = pyuvsim.BeamList([beam0, beam1, beam2, beam3])

    # To fill out other parameters in the UVBeam.
    beam_list.set_str_mode()
    beam_list.set_obj_mode()

    beam_dict = {'ANT1': 0, 'ANT2': 1, 'ANT3': 2, 'ANT4': 3}

    # Error conditions:
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
    assert uv_obj.x_orientation == "east"

    pyuvsim.simsetup._complete_uvdata(uv_obj, inplace=True)

    with open(param_filename, 'r') as fhandle:
        param_dict = yaml.safe_load(fhandle)
    expected_ofilepath = pyuvsim.utils.write_uvdata(
        uv_obj, param_dict, return_filename=True, dryrun=True
    )
    assert './sim_results.uvfits' == expected_ofilepath

    # Spoof attributes that won't match.
    uv_obj.antenna_names = uv_obj.antenna_names.tolist()
    uv_obj.antenna_diameters = hera_uv.antenna_diameters
    uv_obj.history = hera_uv.history

    uvfits_required_extra = [
        "_antenna_positions",
        "_gst0",
        "_rdate",
        "_earth_omega",
        "_dut1",
        "_timesys",
    ]
    for attr in uvfits_required_extra:
        param = getattr(uv_obj, attr)
        if param.value is None:
            param.value = param.spoof_val
            setattr(uv_obj, attr, param)

    assert new_beam_dict == beam_dict
    assert new_beam_list == beam_list

    # remove filename attribute to ensure equality
    hera_uv.filename = None
    assert uv_obj == hera_uv


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
def test_uvdata_keyword_init(case, tmpdir):
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
        assert base_kwargs['polarization_array'] == uvd.get_pols()
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
        assert uvd.get_pols() == new_kwargs['polarization_array']
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
        new_kwargs['path_out'] = str(tmpdir)
        new_kwargs['write_files'] = True

        uvd = pyuvsim.simsetup.initialize_uvdata_from_keywords(**new_kwargs)
        layout_path = str(tmpdir.join(layout_fname))
        obsparam_path = str(tmpdir.join(obsparam_fname))
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
    with uvtest.check_warnings(
        UserWarning,
        match='The integration time is not constant. '
        'Using the shortest integration time'
    ):
        pyuvsim.simsetup.uvdata_to_config_file(
            uv0,
            telescope_config_name=os.path.join(path, telescope_config),
            layout_csv_name=os.path.join(path, layout_fname),
            path_out=opath
        )

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


@pytest.mark.parametrize(
    "arrangement,text_cat",
    [
        ('cross', 'mock_cross_2458098.27471.txt'),
        ('hera_text', 'mock_hera_text_2458098.27471.txt'),
        ('long-line', 'mock_long-line_2458098.27471.txt'),
        ('off-zenith', 'mock_off-zenith_2458098.27471.txt'),
        ('triangle', 'mock_triangle_2458098.27471.txt'),
        ('random', 'mock_random_2458098.27471.txt'),
        ('zenith', 'mock_zenith_2458098.27471.txt'),
    ]
)
def test_mock_catalogs(arrangement, text_cat):
    time = Time(2458098.27471265, scale='utc', format='jd')
    cat, mock_kwds = pyuvsim.simsetup.create_mock_catalog(time, arrangement, rseed=2458098)

    # For each mock catalog, verify the Ra/Dec source positions against a saved text catalog.
    radec_catalog = pyradiosky.SkyModel()
    radec_catalog.read_text_catalog(
        os.path.join(SIM_DATA_PATH, 'test_catalogs', text_cat)
    )
    assert np.all(radec_catalog == cat)


def test_saved_mock_catalog():
    time = Time(2458098.27471265, scale='utc', format='jd')
    cat, mock_kwds = pyuvsim.simsetup.create_mock_catalog(time, 'random', Nsrcs=100, save=True)
    loc = eval(mock_kwds['array_location'])
    loc = EarthLocation.from_geodetic(loc[1], loc[0], loc[2])  # Lon, Lat, alt
    fname = 'mock_catalog_random.npz'
    alts_reload = np.load(fname)['alts']
    cat.update_positions(time, loc)
    alt, az = cat.alt_az
    os.remove(fname)
    assert np.allclose(alts_reload, np.degrees(alt))


@pytest.mark.parametrize('min_alt', [-20, 0, None, 50])
def test_randsource_minalt(min_alt):
    time = Time(2458098.27471265, scale='utc', format='jd')
    cat, mock_kwds = pyuvsim.simsetup.create_mock_catalog(
        time, 'random', Nsrcs=100, min_alt=min_alt
    )
    loc = eval(mock_kwds['array_location'])
    loc = EarthLocation.from_geodetic(loc[1], loc[0], loc[2])  # Lon, Lat, alt
    cat.update_positions(time, loc)
    alt, az = cat.alt_az
    if min_alt is None:
        min_alt = 30    # Checking default
    assert np.all(alt >= np.radians(min_alt))


def test_randsource_distribution():
    # Check that random sources are uniformly distributed per solid angle

    astropy_healpix = pytest.importorskip("astropy_healpix")
    Nsrcs = 40000
    time = Time(2458098.27471265, scale='utc', format='jd')
    cat, mock_kwds = pyuvsim.simsetup.create_mock_catalog(
        time, 'random', Nsrcs=Nsrcs, min_alt=-90, rseed=2458098
    )
    loc = eval(mock_kwds['array_location'])
    loc = EarthLocation.from_geodetic(loc[1], loc[0], loc[2])  # Lon, Lat, alt
    cat.update_positions(time, loc)
    alt, az = cat.alt_az

    # Bin into low-res HEALPix pixels.
    nside = 32
    npix = 12 * nside**2
    hp = astropy_healpix.HEALPix(nside)

    inds = hp.lonlat_to_healpix(Longitude(az, unit='rad'), Latitude(alt, unit='rad'))
    un, counts = np.unique(inds, return_counts=True)

    # counts should be Poisson-distributed with rate lambda = nsrcs / npix
    # variance and mean should be close to lambda
    lam = Nsrcs / npix
    assert np.isclose(np.mean(counts), lam, atol=1.0)
    assert np.isclose(np.var(counts), lam, atol=1.0)


def test_mock_catalog_error():
    time = Time(2458098.27471265, scale='utc', format='jd')
    with pytest.raises(KeyError, match="Invalid mock catalog arrangement: invalid_catalog_name"):
        pyuvsim.create_mock_catalog(time, 'invalid_catalog_name')


def test_keyword_param_loop(tmpdir):
    # Check that yaml/csv files made by intialize_uvdata_from_keywords will work
    # on their own.
    layout_fname = 'temp_layout_kwdloop.csv'
    obsparam_fname = 'temp_obsparam_kwdloop.yaml'
    path_out = str(tmpdir)
    # add some jiggle so you get non-zero uvws
    antpos_enu = (np.ones(30) + np.random.uniform(-10, 10, 30)).reshape((10, 3))
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


def test_multi_analytic_beams(tmpdir):
    # Test inline definitions of beam attributes.
    # eg. (in beam configuration file):
    #
    # beam_paths:
    #   0 : airy, diameter=14
    #   1 : airy, diameter=20
    #   2 : gaussian, sigma=0.5
    par_fname = str(tmpdir.join('test_teleconfig.yaml'))
    layout_fname = str(tmpdir.join('test_layout_5ant.csv'))

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
        param_dict, str(tmpdir))

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


def test_beamlist_init_errors():
    telescope_config_name = os.path.join(SIM_DATA_PATH, 'bl_lite_mixed.yaml')
    with open(telescope_config_name, 'r') as yf:
        telconfig = yaml.safe_load(yf)

    # The path for beam 0 is invalid, and it's not needed for this test.
    del telconfig['beam_paths'][0]
    beam_ids = np.arange(1, 6)

    bad_conf_0 = copy.deepcopy(telconfig)
    bad_conf_0['beam_paths'][1] = 1.35
    with pytest.raises(ValueError, match="Beam model is not properly specified"):
        pyuvsim.simsetup._construct_beam_list(beam_ids, bad_conf_0)

    bad_conf_1 = copy.deepcopy(telconfig)
    del bad_conf_1['beam_paths'][1]['type']

    with pytest.raises(ValueError, match="Beam model must have a"):
        pyuvsim.simsetup._construct_beam_list(beam_ids, bad_conf_1)

    bad_conf_2 = copy.deepcopy(telconfig)
    bad_conf_2['beam_paths'][1]['type'] = 'unsupported_type'
    with pytest.raises(ValueError, match="Undefined beam model type"):
        pyuvsim.simsetup._construct_beam_list(beam_ids, bad_conf_2)

    # Check that spline_interp_opts is passed along correctly to BeamList
    telconfig['spline_interp_opts'] = {'kx' : 2, 'ky' : 2}
    beam_list = pyuvsim.simsetup._construct_beam_list(beam_ids, telconfig)
    assert beam_list.spline_interp_opts is not None


def test_beamlist_init():
    telescope_config_name = os.path.join(SIM_DATA_PATH, 'bl_lite_mixed.yaml')
    with open(telescope_config_name, 'r') as yf:
        telconfig = yaml.safe_load(yf)

    telconfig['beam_paths'][0] = os.path.join(SIM_DATA_PATH, 'HERA_NicCST.uvbeam')

    beam_list = pyuvsim.simsetup._construct_beam_list(np.arange(6), telconfig)
    beam_list.set_obj_mode()

    # How the beam attributes should turn out for this file:
    assert isinstance(beam_list[0], UVBeam)
    assert beam_list[1].type == 'airy'
    assert beam_list[1].diameter == 16
    assert beam_list[2].type == 'gaussian'
    assert beam_list[2].sigma == 0.03
    assert beam_list[3].type == 'airy'
    assert beam_list[3].diameter == 12
    assert beam_list[4].type == 'gaussian'
    assert beam_list[4].diameter == 14
    assert beam_list[5].type == 'gaussian'
    assert beam_list[5].diameter == 12


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


@pytest.fixture(scope='module')
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
        stokes=stokes * units.Jy,
        spectral_type='full',
        freq_array=freqs
    )

    return sky


@pytest.mark.filterwarnings("ignore:Input ra and dec parameters are being used instead")
@pytest.mark.parametrize('unit', ['Jy', 'K'])
def test_skymodeldata_with_quantity_stokes(unit, cat_with_some_pols):
    # Support for upcoming pyradiosky change setting SkyModel.stokes
    # to an astropy Quantity.
    if unit == 'Jy':
        sky = cat_with_some_pols
    else:
        pytest.importorskip('analytic_diffuse')
        sky, _ = pyuvsim.simsetup.create_mock_catalog(
            Time.now(), arrangement='diffuse', diffuse_model='monopole', map_nside=16
        )
    if not isinstance(sky.stokes, units.Quantity):
        sky.stokes *= units.Unit(unit)

    smd = pyuvsim.simsetup.SkyModelData(sky)
    assert np.all(sky.stokes.to_value(unit)[0] == smd.stokes_I)
    assert smd.flux_unit == unit

    sky2 = smd.get_skymodel()
    if units.Quantity != pyradiosky.SkyModel()._stokes.expected_type:
        sky.stokes = sky.stokes.to_value(unit)
    assert sky2 == sky


@pytest.mark.filterwarnings("ignore:Input ra and dec parameters are being used instead")
@pytest.mark.parametrize('component_type', ['point', 'healpix'])
def test_skymodeldata(component_type, cat_with_some_pols):
    # Test that SkyModelData class can properly recreate a SkyModel and subselect.
    if component_type == 'point':
        sky = cat_with_some_pols
    else:
        pytest.importorskip('astropy_healpix')
        path = os.path.join(SKY_DATA_PATH, 'healpix_disk.skyh5')
        sky = pyradiosky.SkyModel()
        sky.read_skyh5(path)

    smd = pyuvsim.simsetup.SkyModelData(sky)

    assert (smd.ra == sky.ra.deg).all()
    assert (smd.dec == sky.dec.deg).all()

    if isinstance(sky.stokes, units.Quantity):
        smd.stokes_I *= units.Unit(smd.flux_unit)
        if smd.polarized is not None:
            smd.stokes_Q *= units.Unit(smd.flux_unit)
            smd.stokes_U *= units.Unit(smd.flux_unit)
            smd.stokes_V *= units.Unit(smd.flux_unit)

    assert (smd.stokes_I == sky.stokes[0]).all()

    if smd.polarized is not None:
        assert (smd.stokes_Q == sky.stokes[..., smd.polarized][1]).all()
        assert (smd.stokes_U == sky.stokes[..., smd.polarized][2]).all()
        assert (smd.stokes_V == sky.stokes[..., smd.polarized][3]).all()

    # Make skymodel from SkyModelData.
    sky1 = smd.get_skymodel()
    # history is not copied into SkyModelData.
    sky1.history = sky.history

    assert sky1 == sky

    # Now try with subselection:
    sky1_sub = smd.get_skymodel(range(8, 13))

    assert sky1.check()
    assert sky1_sub.check()
    assert sky1_sub.Ncomponents == 5
    if smd.polarized is not None:
        assert sky1_sub._n_polarized == 1


@pytest.mark.parametrize('inds', [range(30), range(5), np.arange(9, 14)])
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
    # Check that downselecting doesn't copy length-Ncomponent arrays.

    smd = pyuvsim.simsetup.SkyModelData(cat_with_some_pols)
    smd_copy = smd.subselect(inds)
    assert smd_copy.ra.base is smd.ra.base
    assert smd_copy.dec.base is smd.dec.base
    assert smd_copy.stokes_I.base is smd.stokes_I.base


def test_set_lsts_errors():
    # Error cases on set_lsts function.
    uv0 = UVData()
    uv0.read_uvfits(longbl_uvfits_file)
    uv0.lst_array = None

    uv0.extra_keywords['world'] = 'moon'
    if not pyuvsim.astropy_interface.hasmoon:
        with pytest.raises(ValueError, match="Cannot construct lsts for MoonLocation"):
            pyuvsim.simsetup._set_lsts_on_uvdata(uv0)

    uv0.extra_keywords['world'] = 'tatooine'
    with pytest.raises(ValueError, match="Invalid world tatooine."):
        pyuvsim.simsetup._set_lsts_on_uvdata(uv0)
