
import os
import copy

from pyuvdata import UVBeam
import pytest

import pyuvsim
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH


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


def test_convert_loop(beam_objs):
    beams = beam_objs
    beams[0].freq_interp_kind = 'linear'
    beams[1].freq_interp_kind = 'cubic'

    # Should warn about inconsistent params on UVBeams.
    with pytest.warns(UserWarning) as confwarn:
        beamlist = pyuvsim.BeamList(beams)
    assert str(confwarn.pop().message).startswith("Conflicting settings for")

    # Convert beams to strings:
    # Fail, because UVBeams are inconsistent.
    with pytest.raises(ValueError, match='Conflicting settings for'):
        beamlist.set_str_mode()

    beams[1].freq_interp_kind = 'linear'

    beamlist.set_str_mode()

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
    with pytest.warns(UserWarning) as confwarn:
        beamlist.append(uvb)
    assert str(confwarn.pop().message).startswith("Conflicting settings for")
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
