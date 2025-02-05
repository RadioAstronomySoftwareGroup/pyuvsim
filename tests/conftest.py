# Copyright (c) 2019 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

"""Testing environment setup and teardown for pytest."""

import contextlib
import os

import numpy as np
import pytest
from astropy.coordinates import EarthLocation
from astropy.time import Time
from astropy.utils import iers
from pyuvdata import UVBeam

try:
    from pyuvsim import mpi
except ImportError:
    mpi = None

# this is pretty hacky, but want just one convenient place
# to set if we can use the parallel tests
from importlib.util import find_spec

pytest.pyuvsim_can_parallel = find_spec("pytest_mpi") is not None

issubproc = os.environ.get("TEST_IN_PARALLEL", 0)
with contextlib.suppress(ValueError):
    issubproc = bool(int(issubproc))


def pytest_collection_modifyitems(session, config, items):
    # Enforce that the profiler test is run last.

    if len(items) <= 1:
        return
    profiler_index = None
    for ii, it in enumerate(items):
        if "profiler" in it.name:
            profiler_index = ii
            break
    if profiler_index is not None:
        items.append(items.pop(profiler_index))  # Move profiler tests to the end.


def pytest_addoption(parser):
    parser.addoption(
        "--nompi", action="store_true", help="skip mpi-parallelized tests."
    )
    parser.addoption(
        "--savesim",
        action="store_true",
        default=False,
        help="saves reference simulation output to the current working directory",
    )
    parser.addoption(
        "--refsim",
        action="append",
        default=[],
        help="specify an available reference simulation to pass to test functions",
        choices=[
            "1.1_baseline_number",
            "1.2_time_axis",
            "1.3_frequency_axis",
            "1.4_source_axis",
            "1.5_uvbeam",
            "1.6_healpix",
            "1.7_multi_beam",
            "1.8_lunar",
        ],
    )


def pytest_generate_tests(metafunc):
    if "refsim" in metafunc.fixturenames:
        metafunc.parametrize("refsim", metafunc.config.getoption("refsim"))


def pytest_runtest_setup(item):
    if "parallel" in item.keywords:
        if mpi is None:
            pytest.skip("Need mpi4py to run parallelized tests.")
        elif item.config.getoption("nompi", False):
            pytest.skip("Skipping parallelized tests with --nompi option.")


@pytest.fixture(autouse=True, scope="session")
def _setup_and_teardown_package():
    # Do a calculation that requires a current IERS table. This will trigger
    # automatic downloading of the IERS table if needed, including trying the
    # mirror site in python 3 (but won't redownload if a current one exists).
    # If there's not a current IERS table and it can't be downloaded, turn off
    # auto downloading for the tests and turn it back on once all tests are
    # completed (done by extending auto_max_age).
    try:
        t1 = Time.now()
        t1.ut1  # noqa
    except Exception:
        iers.conf.auto_max_age = None

    # yield to allow tests to run
    yield

    iers.conf.auto_max_age = 30


@pytest.fixture
def goto_tempdir(tmpdir):
    # Run test within temporary directory.
    newpath = str(tmpdir)
    cwd = os.getcwd()
    os.chdir(newpath)

    yield newpath

    os.chdir(cwd)


@pytest.fixture(scope="session")
def cst_beam():
    beam = UVBeam()

    freqs = [150e6, 123e6]

    try:
        from pyuvdata.datasets import fetch_data

        beam_files = fetch_data(["hera_fagnoni_dipole_150", "hera_fagnoni_dipole_123"])
    except ImportError:
        # This can be removed once we require pyuvdata > 3.2.3
        from pyuvdata.data import DATA_PATH as UV_DATA_PATH

        cst_files = ["HERA_NicCST_150MHz.txt", "HERA_NicCST_123MHz.txt"]
        beam_files = [os.path.join(UV_DATA_PATH, "NicCSTbeams", f) for f in cst_files]

    kwargs = {
        "beam_type": "efield",
        "frequency": freqs,
        "telescope_name": "HERA",
        "feed_name": "PAPER",
        "feed_version": "0.1",
        "feed_pol": ["x"],
        "model_name": "E-field pattern - Rigging height 4.9m",
        "model_version": "1.0",
        "feed_angle": [np.pi / 2],
        "mount_type": "fixed",
    }
    beam.read_cst_beam(beam_files, **kwargs)
    beam.feed_angle[1] = 0.0
    beam.peak_normalize()
    return beam


@pytest.fixture(scope="session")
def hera_loc():
    return EarthLocation(lat="-30d43m17.5s", lon="21d25m41.9s", height=1073.0)


@pytest.fixture(scope="session")
def apollo_loc():
    with contextlib.suppress(ImportError):
        from lunarsky import MoonLocation

        return MoonLocation(lat=0.6875, lon=24.433, height=0)
    return None


@pytest.fixture(scope="session")
def savesim(request):
    return request.config.getoption("--savesim")
