
import os

from pyuvdata import UVBeam
import pytest

import pyuvsim
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH


herabeam_default = os.path.join(SIM_DATA_PATH, 'HERA_NicCST.uvbeam')

# Ignore warnings of pending sigma deprecation
@pytest.mark.filterwarnings('ignore:Achromatic gaussian')
class TestBeamList():
    def setup(self):
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

        self.beams = beams

    def test_convert_loop(self):
        # Convert beams to strings:
        beamlist = pyuvsim.BeamList(self.beams)

        beamlist.set_str_mode()
        for bs in beamlist:
            assert isinstance(bs, str)
        assert beamlist._obj_beam_list == []

        # Convert strings to beams.
        beamlist.set_obj_mode()

        for bi, b in enumerate(beamlist):
            assert b == self.beams[bi]

        assert beamlist._str_beam_list == []

    def test_insert_beam(self):
        # Insert a new beam model into the array
        beamlist = pyuvsim.BeamList(self.beams)
        assert beamlist.uvb_params['freq_interp_kind'] == 'linear'

        uvb2 = UVBeam()
        uvb2.read_beamfits(herabeam_default)
        uvb2.extra_keywords['beam_path'] = herabeam_default
        uvb2.freq_interp_kind = 'quartic'

        beamlist[0] = uvb2

        assert uvb2 is not self.beams[0]
        assert beamlist.uvb_params['freq_interp_kind'] == 'quartic'

        beamlist.set_str_mode()

        # If inserting an object while in string mode, it should
        # convert the object on the fly.
        beamlist[0] = self.beams[0]

        # Since the new object has freq_interp_kind set to linear
        assert beamlist.uvb_params['freq_interp_kind'] == 'linear'

    def test_error(self):

        uvb2 = UVBeam()
        uvb2.read_beamfits(herabeam_default)
        uvb2.freq_interp_kind = 'quartic'

        new_beam_list = [self.beams[0], uvb2]

        pytest.raises(ValueError, pyuvsim.BeamList, new_beam_list)
