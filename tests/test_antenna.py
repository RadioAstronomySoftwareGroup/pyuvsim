# Copyright (c) 2021 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License
import copy
import os

import numpy as np
import pytest
import yaml
from astropy import units

from pyuvsim import simsetup
from pyuvsim.antenna import Antenna
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
from pyuvsim.telescope import BeamList, Telescope


def test_jones_set_spline(cst_beam, hera_loc):
    # Run get_beam_jones with spline options.
    array_location = hera_loc
    beam0 = cst_beam.copy()
    telescope_config_name = os.path.join(SIM_DATA_PATH, "mwa128_config.yaml")
    with open(telescope_config_name) as yf:
        telconfig = yaml.safe_load(yf)
    telconfig["spline_interp_opts"] = {"kx": 1, "ky": 1}
    telconfig["beam_paths"][1] = beam0

    Nfreqs = 10
    freqs = np.linspace(100, 130, Nfreqs) * 1e6 * units.Hz

    beam_list = simsetup._construct_beam_list(np.arange(2), telconfig, freq_array=freqs)
    assert len(beam_list) == 2

    assert beam0 == beam_list[-1].beam

    # Make antenna that uses beam #1
    antenna = Antenna("ant1", 1, np.array([0, 10, 0]), 1)
    array = Telescope("telescope_name", array_location, beam_list)

    altaz = [[0.0134], [1.0]]

    alts = np.linspace(0, np.pi / 4, 50)
    azs = np.linspace(0, 2 * np.pi, 50)

    alts, azs = np.meshgrid(alts, azs)

    altaz = np.zeros((50**2, 3))
    altaz[:, 0] = alts.flatten()
    altaz[:, 1] = azs.flatten()

    jones_matrix = antenna.get_beam_jones(
        array, altaz, 150e6, interpolation_function="az_za_simple"
    )

    # These are just values from a run, so this just tests for unexpected changes.
    expected_jones = np.array(
        [
            [
                [
                    -4.57061296e-04 - 3.88626249e-04j,
                    -7.28285993e-05 - 1.45479743e-04j,
                    -7.28285993e-05 - 1.45479743e-04j,
                ],
                [
                    -3.35886569e-02 - 1.83636844e-02j,
                    -3.36205621e-02 - 1.84105336e-02j,
                    -3.36205621e-02 - 1.84105336e-02j,
                ],
            ],
            [
                [
                    -1.04000784e-02 - 1.11629186e-02j,
                    -1.03973090e-02 - 1.11516998e-02j,
                    -1.03973090e-02 - 1.11516998e-02j,
                ],
                [
                    5.32870283e-04 + 1.16831373e-04j,
                    1.26946128e-06 - 1.22843330e-06j,
                    1.26946128e-06 - 1.22843330e-06j,
                ],
            ],
        ]
    )

    np.testing.assert_allclose(expected_jones, jones_matrix)


def test_jones_set_interp(cst_beam, hera_loc):
    # check setting the interpolation method

    array_location = hera_loc

    beam = cst_beam.copy()
    beam_list = BeamList([beam])
    antenna1 = Antenna("ant1", 1, np.array([0, 10, 0]), 0)
    array = Telescope("telescope_name", array_location, beam_list)
    source_altaz = np.array([[0.0], [np.pi / 4.0]])
    freq = 123e6 * units.Hz

    jones = antenna1.get_beam_jones(array, source_altaz, freq, freq_interp_kind="cubic")
    jones0 = antenna1.get_beam_jones(array, source_altaz, freq)
    jones1 = antenna1.get_beam_jones(
        array, source_altaz, freq, freq_interp_kind="linear"
    )
    jones2 = antenna1.get_beam_jones(array, source_altaz, freq)

    assert np.all(jones2 == jones0)
    assert np.all(jones1 == jones)
    assert np.all(jones1 == jones0)


@pytest.mark.parametrize(
    ("problem", "err_msg"),
    [
        (
            "beam_id",
            "This antenna beam_id is 1, which is too large "
            "for the beam_list, which has length 1.",
        ),
        ("beam_type", "Beam type must be efield!"),
        ("normalization", "UVBeams must be peak normalized."),
    ],
)
def test_get_beam_jones_errors(cst_beam, hera_loc, problem, err_msg):
    array_location = hera_loc

    beam = cst_beam.copy()
    if problem == "normalization":
        beam.data_normalization = "physical"
        peak_normalize = False
    else:
        peak_normalize = True

    if problem == "beam_type":
        beam.efield_to_power()
        beam_type = "power"
    else:
        beam_type = "efield"
    beam_list = BeamList([beam], beam_type=beam_type, peak_normalize=peak_normalize)
    if problem == "beam_id":
        beam_id = 1
    else:
        beam_id = 0

    antenna1 = Antenna("ant1", 1, np.array([0, 10, 0]), beam_id)
    array = Telescope("telescope_name", array_location, beam_list)
    source_altaz = np.array([[0.0], [np.pi / 4.0]])
    freq = 123e6 * units.Hz

    with pytest.raises(ValueError, match=err_msg):
        antenna1.get_beam_jones(array, source_altaz, freq, freq_interp_kind="cubic")


def test_ant_comparison():
    antenna1 = Antenna("ant1", 1, np.array([0, 10, 0]), 1)
    antenna2 = Antenna("ant2", 2, np.array([0, 20, 0]), 1)

    ant1_copy = copy.deepcopy(antenna1)

    assert antenna1 < antenna2
    assert antenna1 <= antenna2
    assert antenna1 <= antenna1
    assert antenna1 == ant1_copy

    assert antenna2 > antenna1
    assert antenna2 >= antenna2
    assert antenna1 >= antenna1
