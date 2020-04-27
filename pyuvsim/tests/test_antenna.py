
import numpy as np
import os
import yaml
from astropy.coordinates import EarthLocation

import pyuvsim
import pyuvsim.tests as simtest
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH


def test_get_beam_jones():
    # Run get_beam_jones with spline options.
    array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s', height=1073.)
    beam0 = simtest.make_cst_beams()
    beam0.freq_interp_kind = 'cubic'
    telescope_config_name = os.path.join(SIM_DATA_PATH, 'mwa128_config.yaml')
    with open(telescope_config_name, 'r') as yf:
        telconfig = yaml.safe_load(yf)
    telconfig['spline_interp_opts'] = {'kx' : 1, 'ky' : 1}
    beam_list = pyuvsim.simsetup._construct_beam_list(np.arange(1), telconfig)
    beam_list.set_obj_mode()
    beam_list.append(beam0)

    assert beam0 is beam_list[-1]

    # Make antenna that uses beam #1
    antenna = pyuvsim.Antenna('ant1', 1, np.array([0, 10, 0]), 1)
    array = pyuvsim.Telescope('telescope_name', array_location, beam_list)

    altaz = [[0.0134], [1.0]]

    alts = np.linspace(0, np.pi / 4, 50)
    azs = np.linspace(0, 2 * np.pi, 50)

    alts, azs = np.meshgrid(alts, azs)

    altaz = np.zeros((50**2, 3))
    altaz[:, 0] = alts.flatten()
    altaz[:, 1] = azs.flatten()

    antenna.get_beam_jones(array, altaz, 150e6)
