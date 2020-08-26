
import numpy as np
import os
import yaml
from astropy import units
import pytest

import pyuvsim
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH


def test_jones_set_spline(cst_beam, hera_loc):
    # Run get_beam_jones with spline options.
    array_location = hera_loc
    beam0 = cst_beam.copy()
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

    array.beam_list.spline_interp_opts = None
    antenna.get_beam_jones(array, altaz, 150e6, interpolation_function='az_za_simple')


def test_jones_set_interp(cst_beam, hera_loc):
    # check setting the interpolation method

    array_location = hera_loc

    beam = cst_beam.copy()
    beam.freq_interp_kind = None

    beam_list = pyuvsim.BeamList([beam])
    antenna1 = pyuvsim.Antenna('ant1', 1, np.array([0, 10, 0]), 0)
    array = pyuvsim.Telescope('telescope_name', array_location, beam_list)
    source_altaz = np.array([[0.0], [np.pi / 4.]])
    freq = 123e6 * units.Hz

    with pytest.raises(ValueError, match='freq_interp_kind must be set'):
        antenna1.get_beam_jones(array, source_altaz, freq)

    jones = antenna1.get_beam_jones(array, source_altaz, freq, freq_interp_kind='cubic')
    assert beam.freq_interp_kind == 'cubic'
    jones0 = antenna1.get_beam_jones(array, source_altaz, freq)
    jones1 = antenna1.get_beam_jones(array, source_altaz, freq, freq_interp_kind='linear')
    assert beam.freq_interp_kind == 'linear'
    jones2 = antenna1.get_beam_jones(array, source_altaz, freq)

    assert (np.all(jones2 == jones0)
            and np.all(jones1 == jones)
            and np.all(jones1 == jones0))


def test_set_interps(cst_beam, hera_loc):
    array_location = hera_loc

    beam = cst_beam.copy()
    beam.interpolation_function = None

    beam_list = pyuvsim.BeamList([beam])
    antenna1 = pyuvsim.Antenna('ant1', 1, np.array([0, 10, 0]), 0)
    array = pyuvsim.Telescope('telescope_name', array_location, beam_list)
    source_altaz = np.array([[0.0], [np.pi / 4.]])
    freq = 123e6 * units.Hz

    with pytest.warns(UserWarning, match="UVBeam interpolation_function is not set"):
        antenna1.get_beam_jones(array, source_altaz, freq)

    assert beam.interpolation_function == 'az_za_simple'
