
import os
import copy

import numpy as np
from astropy.coordinates import EarthLocation

import pytest
from pyuvdata import UVBeam
import pyuvdata.tests as uvtest

import pyuvsim
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
from pyuvsim.telescope import BeamConsistencyError

herabeam_default = os.path.join(SIM_DATA_PATH, 'HERA_NicCST.uvbeam')


# Ignore warnings of pending sigma deprecation
@pytest.mark.filterwarnings('ignore:Achromatic gaussian')
@pytest.fixture(scope='module')
def beam_objs_main():
    uvb = UVBeam()
    uvb.read_beamfits(herabeam_default)
    uvb.extra_keywords['beam_path'] = herabeam_default

    uvb2 = uvb.copy()

    beams = [uvb, uvb2]

    beams.append(pyuvsim.AnalyticBeam('uniform'))
    diameter_m = 14.
    beams.append(pyuvsim.AnalyticBeam('airy', diameter=diameter_m))
    sigma = 0.03
    beams.append(pyuvsim.AnalyticBeam('gaussian', sigma=sigma))
    ref_freq, alpha = 100e6, -0.5
    beams.append(pyuvsim.AnalyticBeam('gaussian', sigma=sigma,
                 ref_freq=ref_freq, spectral_index=alpha))
    return beams


@pytest.mark.filterwarnings('ignore:Achromatic gaussian')
@pytest.fixture(scope='function')
def beam_objs(beam_objs_main):
    beams_copy = copy.deepcopy(beam_objs_main)

    return beams_copy


@pytest.mark.filterwarnings('ignore:Achromatic gaussian')
def test_convert_loop(beam_objs):
    beams = beam_objs
    beams[0].freq_interp_kind = 'linear'
    beams[1].freq_interp_kind = 'cubic'

    # Should warn about inconsistent params on UVBeams.
    with uvtest.check_warnings(UserWarning, match="Conflicting settings for"):
        beamlist = pyuvsim.BeamList(beams)

    # Convert beams to strings:
    # Fail, because UVBeams are inconsistent.
    with pytest.raises(ValueError, match='Conflicting settings for'):
        beamlist.set_str_mode()

    beams[1].freq_interp_kind = 'linear'

    beamlist.set_str_mode()

    # check that _obj_to_str on a string beam works
    beamlist2 = copy.deepcopy(beamlist)
    beamlist2._obj_to_str(beamlist2[0])

    assert beamlist2 == beamlist

    assert beamlist.uvb_params['freq_interp_kind'] == 'linear'

    for bs in beamlist:
        assert isinstance(bs, str)
    assert beamlist._obj_beam_list == []

    # Convert strings to beams. Need to set additional parameters for comparison.
    beamlist._set_params_on_uvbeams(beams)
    beamlist.set_obj_mode()
    for bi, b in enumerate(beamlist):
        assert b == beams[bi]

    assert beamlist._str_beam_list == []

    # Reset UVBeams
    beams[0].freq_interp_kind = None
    beams[1].freq_interp_kind = None


@pytest.mark.filterwarnings('ignore:Achromatic gaussian')
def test_object_mode(beam_objs):
    beams = beam_objs
    beamlist = pyuvsim.BeamList(beams)

    beamlist[0].freq_interp_kind = 'cubic'

    uvb = copy.deepcopy(beams[0])
    uvb.freq_interp_kind = 'quartic'

    # Warn if inserted object mismatches.
    with uvtest.check_warnings(UserWarning, match="Conflicting settings for"):
        beamlist.append(uvb)
    assert len(beamlist) == 7

    # Error if converting to string mode with mismatched keywords:
    with pytest.raises(ValueError, match='Conflicting settings '):
        beamlist.set_str_mode()

    beamlist._set_params_on_uvbeams(beamlist._obj_beam_list)

    # Error if converting to string mode without beam_paths:
    beamlist[0].extra_keywords.pop('beam_path')
    with pytest.raises(ValueError, match='Need to set '):
        beamlist.set_str_mode()

    # Insert string -- Converts to object
    new_anabeam = 'analytic_gaussian_sig=3.0'
    beamlist[-1] = new_anabeam

    assert isinstance(beamlist[-1], pyuvsim.AnalyticBeam)
    assert beamlist[-1].sigma == 3.0


@pytest.mark.filterwarnings('ignore:Achromatic gaussian')
def test_string_mode(beam_objs):
    beams = beam_objs
    beamlist = pyuvsim.BeamList(beams)
    beamlist.set_str_mode()

    uvb = beams[0]
    uvb.freq_interp_kind = 'quartic'

    with pytest.raises(ValueError, match='UVBeam parameters do not'):
        beamlist.append(uvb)

    uvb.freq_interp_kind = beamlist.uvb_params['freq_interp_kind']
    beamlist.append(uvb)

    assert isinstance(beamlist[-1], str)
    beamlist.set_obj_mode()

    # Check that parameters are set properly.
    try:
        new_pars = beamlist._scrape_uvb_params(beamlist._obj_beam_list, strict=True)
        assert new_pars == beamlist.uvb_params
    except ValueError:
        assert False


@pytest.mark.filterwarnings('ignore:Achromatic gaussian')
@pytest.mark.filterwarnings("ignore:Cannot check consistency of a string-mode BeamList")
def test_comparison(beam_objs):
    beamlist = pyuvsim.BeamList(beam_objs)
    beamlist.set_str_mode()

    beamlist2 = pyuvsim.BeamList(beamlist._str_beam_list)
    assert beamlist == beamlist2

    beamlist.set_obj_mode()
    beamlist2.set_obj_mode()
    assert beamlist == beamlist2


def test_no_overwrite(beam_objs):
    # Ensure UVBeam keywords are not overwritten by BeamList.uvb_params
    # while in object mode.
    newbeams = beam_objs
    beamlist = pyuvsim.BeamList(newbeams)
    assert beamlist.uvb_params['freq_interp_kind'] == 'cubic'

    uvb = newbeams[0].copy()
    uvb.freq_interp_kind = 'quintic'

    beamlist.append(uvb)
    assert uvb.freq_interp_kind == 'quintic'
    assert beamlist.uvb_params['freq_interp_kind'] == 'cubic'


@pytest.mark.filterwarnings('ignore:Achromatic gaussian')
def test_beamlist_errors(beam_objs):
    # make a copy to enable Telescope equality checking
    beams = copy.deepcopy(beam_objs)
    beamlist = pyuvsim.BeamList(beams, check=False)

    # Try to make a BeamList with a mixture of strings and objects.
    newlist = copy.deepcopy(beamlist._obj_beam_list)
    newlist[2] = beamlist._obj_to_str(newlist[2])
    with pytest.raises(ValueError, match='Invalid beam list:'):
        pyuvsim.BeamList(newlist)

    # Try to append an invalid beam path while in object mode.
    beam_path = 'invalid_file.uvbeam'
    with pytest.raises(ValueError, match='Invalid file path'):
        beamlist.append(beam_path)

    # test error on beams with different x_orientation
    beams[0].x_orientation = None
    with pytest.raises(
        BeamConsistencyError,
        match='x_orientation of beam 2 is not consistent with beam 1'
    ):
        pyuvsim.BeamList(beams, check=True)

    # test warning on beams with different x_orientation
    beams[0].x_orientation = None
    beams[1].x_orientation = None
    with uvtest.check_warnings(
        UserWarning,
        match="All polarized beams have x_orientation set to None. This will make it "
        "hard to interpret the polarizations of the simulated visibilities.",
    ):
        pyuvsim.BeamList(beams).x_orientation

    # Compare Telescopes with beamlists of different lengths

    del beams[0]
    array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s',
                                   height=1073.)
    tel0 = pyuvsim.Telescope('tel0', array_location, beams)
    tel1 = pyuvsim.Telescope('tel1', array_location, beam_objs)
    assert tel0 != tel1


def test_beamlist_consistency(beam_objs):

    # Does check, but raises no error
    pyuvsim.BeamList(beam_objs[:2])

    pyuvsim.BeamList(beam_objs[:3])


def test_beamlist_consistency_properties(beam_objs):
    beams = beam_objs
    beamlist = pyuvsim.BeamList(beams[:2])
    assert beamlist.x_orientation == beams[0].x_orientation


def test_beamlist_consistency_stringmode(beam_objs):
    beams = beam_objs
    beamlist = pyuvsim.BeamList(beams[:2])
    beamlist.set_str_mode()
    beamlist.check_consistency(force=True)
    assert beamlist.string_mode
    with uvtest.check_warnings(
        UserWarning, match="Cannot check consistency of a string-mode BeamList!"
    ):
        beamlist.check_consistency(force=False)

    beamlist = pyuvsim.BeamList(beams[:3], check=False)
    beamlist.set_str_mode()

    pyuvsim.BeamList(beamlist._str_beam_list, check=True, force_check=True)


def test_beam_basis_type(beam_objs):
    beamlist = pyuvsim.BeamList(beam_objs)

    basis_types = beamlist._get_beam_basis_type()

    assert basis_types == {index: "az_za" for index in range(len(beam_objs))}


@pytest.mark.filterwarnings("ignore:key beam_path in extra_keywords is longer")
@pytest.mark.filterwarnings('ignore:Achromatic gaussian')
def test_beam_basis_type_errors(beam_objs):
    beam_objs[0].pixel_coordinate_system = "orthoslant_zenith"
    beam_objs[0].check()

    beamlist = pyuvsim.BeamList(beam_objs, check=False)
    with pytest.raises(
        ValueError,
        match="pyuvsim currently only supports UVBeams with 'az_za' or "
            "'healpix' pixel coordinate systems."
    ):
        beamlist.check_consistency()

    # test with non-orthogonal basis vectors
    # first construct a beam with non-orthogonal basis vectors
    new_basis_vecs = np.zeros_like(beam_objs[1].basis_vector_array)
    new_basis_vecs[0, 0, :, :] = np.sqrt(0.5)
    new_basis_vecs[0, 1, :, :] = np.sqrt(0.5)
    new_basis_vecs[1, :, :, :] = beam_objs[1].basis_vector_array[1, :, :, :]
    new_data = np.zeros_like(beam_objs[1].data_array)
    # drop all the trailing colons in the slicing below
    new_data[0] = np.sqrt(2) * beam_objs[1].data_array[0]
    new_data[1] = beam_objs[1].data_array[1] - beam_objs[1].data_array[0]
    beam_objs[1].basis_vector_array = new_basis_vecs
    beam_objs[1].data_array = new_data
    beam_objs[1].check()

    beamlist = pyuvsim.BeamList(beam_objs[1:], check=False)
    with pytest.raises(
        ValueError,
        match="pyuvsim currently only supports beams with basis vectors that"
            "are aligned with the azimuth and zenith angle in each pixel."
            "Work is in progress to add other basis vector systems."
    ):
        beamlist.check_consistency()

    beamlist.set_str_mode()
    with pytest.raises(
        ValueError,
        match="Cannot get beam basis type from a string-mode BeamList."
    ):
        beamlist._get_beam_basis_type()


def test_empty_beamlist():
    a = pyuvsim.BeamList(check=False)
    assert a.x_orientation is None
    assert a.beam_type is None


@pytest.mark.filterwarnings("ignore:key beam_path in extra_keywords is longer than 8")
def test_powerbeam_consistency(beam_objs):
    newbeams = beam_objs[:2]
    for beam in newbeams:
        beam.efield_to_power()

    beamlist = pyuvsim.BeamList(newbeams)
    beamlist.check_consistency()


@pytest.mark.filterwarnings('ignore:Achromatic gaussian')
@pytest.mark.filterwarnings("ignore:key beam_path in extra_keywords is longer than 8")
def test_check_azza_full_sky(beam_objs):
    beam = beam_objs

    beamlist = pyuvsim.BeamList(beam)
    assert beamlist.check_all_azza_beams_full_sky()

    # Downselect a beam to a small sky area
    za_max = np.deg2rad(10.0)
    za_inds_use = np.nonzero(beam[0].axis2_array <= za_max)[0]
    beam[0].select(axis2_inds=za_inds_use)

    beamlist = pyuvsim.BeamList(beam)
    assert not beamlist.check_all_azza_beams_full_sky()
