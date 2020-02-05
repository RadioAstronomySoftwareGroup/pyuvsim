
import os

from pyuvdata import UVBeam

import pyuvsim
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH


herabeam_default = os.path.join(SIM_DATA_PATH, 'HERA_NicCST.uvbeam')


def test_beam_list():
    uvb = UVBeam()
    uvb.read_beamfits(herabeam_default)
    uvb.extra_keywords['beam_path'] = herabeam_default
    uvb.freq_interp_kind = 'linear'

    beams = [uvb]

    beams.append(pyuvsim.AnalyticBeam('uniform'))
    diameter_m = 14.
    beams.append(pyuvsim.AnalyticBeam('airy', diameter=diameter_m))
    sigma = 0.03
    beams.append(pyuvsim.AnalyticBeam('gaussian', sigma=sigma))
    ref_freq, alpha = 100e6, -0.5
    beams.append(pyuvsim.AnalyticBeam('gaussian', sigma=sigma,
                 ref_freq=ref_freq, spectral_index=alpha))

    beamlist = pyuvsim.BeamList(beams)

    # Convert beams to strings:
    beamlist.set_str_mode()
    for bs in beamlist:
        assert isinstance(bs, str)
    assert beamlist._obj_beam_list == []

    beamlist.set_obj_mode()

    for bi, b in enumerate(beamlist):
        assert b == beams[bi]

    assert beamlist._str_beam_list == []
