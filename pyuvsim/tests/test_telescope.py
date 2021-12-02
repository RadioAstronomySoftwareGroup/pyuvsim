
import os
import copy
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
def beam_objs():
    uvb = UVBeam()
    uvb.read_beamfits(herabeam_default)
    uvb.extra_keywords['beam_path'] = herabeam_default

    uvb2 = UVBeam()
    uvb2.read_beamfits(herabeam_default)
    uvb2.extra_keywords['beam_path'] = herabeam_default

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
    newbeams = copy.deepcopy(beams)
    beamlist = pyuvsim.BeamList(newbeams)

    beamlist[0].freq_interp_kind = 'cubic'

    uvb = copy.deepcopy(newbeams[0])
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
    newbeams = copy.deepcopy(beam_objs)
    beamlist = pyuvsim.BeamList(newbeams)
    beamlist.set_str_mode()

    uvb = newbeams[0]
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
    newbeams = copy.deepcopy(beam_objs)
    beamlist = pyuvsim.BeamList(newbeams)
    assert beamlist.uvb_params['freq_interp_kind'] == 'cubic'

    uvb = copy.deepcopy(newbeams[0])
    uvb.freq_interp_kind = 'quintic'

    beamlist.append(uvb)
    assert uvb.freq_interp_kind == 'quintic'
    assert beamlist.uvb_params['freq_interp_kind'] == 'cubic'


def test_beamlist_errors(beam_objs):
    newbeams = copy.deepcopy(beam_objs)
    beamlist = pyuvsim.BeamList(newbeams, check=False)

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
    newbeams[0].x_orientation = None
    with pytest.raises(BeamConsistencyError):
        pyuvsim.BeamList(newbeams, check=True)

    # test warning on beams with different x_orientation
    newbeams[0].x_orientation = None
    newbeams[1].x_orientation = None
    with uvtest.check_warnings(
        UserWarning,
        match="All polarized beams have x_orientation set to None. This will make it "
        "hard to interpret the polarizations of the simulated visibilities.",
    ):
        pyuvsim.BeamList(newbeams).x_orientation

    # Compare Telescopes with beamlists of different lengths
    del newbeams[0]
    array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s',
                                   height=1073.)
    tel0 = pyuvsim.Telescope('tel0', array_location, newbeams)
    tel1 = pyuvsim.Telescope('tel1', array_location, beam_objs)
    assert tel0 != tel1


def test_beamlist_consistency(beam_objs):
    newbeams = copy.deepcopy(beam_objs)

    # Does check, but raises no error
    pyuvsim.BeamList(newbeams[:2])

    pyuvsim.BeamList(newbeams[:3])


def test_beamlist_consistency_properties(beam_objs):
    newbeams = copy.deepcopy(beam_objs)
    beamlist = pyuvsim.BeamList(newbeams[:2])
    assert beamlist.x_orientation == newbeams[0].x_orientation


def test_beamlist_consistency_stringmode(beam_objs):
    newbeams = copy.deepcopy(beam_objs)
    beamlist = pyuvsim.BeamList(newbeams[:2])
    beamlist.set_str_mode()
    beamlist.check_consistency(force=True)
    assert beamlist.string_mode
    with uvtest.check_warnings(
        UserWarning, match="Cannot check consistency of a string-mode BeamList!"
    ):
        beamlist.check_consistency(force=False)

    beamlist = pyuvsim.BeamList(newbeams[:3], check=False)
    beamlist.set_str_mode()

    pyuvsim.BeamList(beamlist._str_beam_list, check=True, force_check=True)
