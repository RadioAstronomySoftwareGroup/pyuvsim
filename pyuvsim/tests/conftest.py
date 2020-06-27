# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2019 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

"""Testing environment setup and teardown for pytest."""
import os
import warnings
from subprocess import check_output, CalledProcessError

import pytest
from astropy.time import Time
from astropy.utils import iers
from astropy.coordinates import EarthLocation
from pyuvdata import UVBeam
from pyuvdata.data import DATA_PATH


def pytest_collection_modifyitems(session, config, items):
    # Enforce that the profiler test is run last.

    if len(items) <= 1:
        return
    for ii, it in enumerate(items):
        if 'profiler' in it.name:
            break
    items.append(items.pop(ii))     # Move to the end.


def pytest_configure(config):
    """Register an additional marker."""
    config.addinivalue_line(
        "markers",
        "parallel: mark test to run in parallel"
    )


def pytest_runtest_call(item):
    # If a test is to be run in parallel, spawn a subprocess that runs it in parallel.
    # Avoids the bug that comes up when a test fails in an inconsistent mpi state.

    ispartest = item.get_closest_marker("parallel") is not None
    if not ispartest:
        return

    issubproc = os.environ.get('TEST_IN_PARALLEL', 0)
    try:
        issubproc = int(issubproc)
    except ValueError:
        pass

    os.environ['TEST_IN_PARALLEL'] = '1'

    call = ['mpirun', '--oversubscribe', '-n', "2",
            "pytest", "{:s}::{:s}".format(str(item.fspath), str(item.name))]
    call.extend(["--tb=no", '-q'])
    if not issubproc:
        try:
            check_output(call).decode('utf-8')
        except CalledProcessError:
            pass
        finally:
            del os.environ['TEST_IN_PARALLEL']


@pytest.fixture(autouse=True, scope="session")
def setup_and_teardown_package():
    # Do a calculation that requires a current IERS table. This will trigger
    # automatic downloading of the IERS table if needed, including trying the
    # mirror site in python 3 (but won't redownload if a current one exists).
    # If there's not a current IERS table and it can't be downloaded, turn off
    # auto downloading for the tests and turn it back on once all tests are
    # completed (done by extending auto_max_age).
    try:
        t1 = Time.now()
        t1.ut1
    except (Exception):
        iers.conf.auto_max_age = None

    # yield to allow tests to run
    yield

    iers.conf.auto_max_age = 30


@pytest.fixture(autouse=True)
def ignore_deprecation():
    warnings.filterwarnings('ignore', message='Achromatic gaussian beams will not be supported',
                            category=PendingDeprecationWarning)
    warnings.filterwarnings('ignore', message='"initialize_catalog_from_params will not return'
                                              ' recarray by default in the future.',
                            category=PendingDeprecationWarning)


@pytest.fixture(scope='session')
def cst_beam():
    beam = UVBeam()
    beam.freq_interp_kind = 'linear'

    freqs = [150e6, 123e6]

    cst_files = ['HERA_NicCST_150MHz.txt', 'HERA_NicCST_123MHz.txt']
    beam_files = [os.path.join(DATA_PATH, 'NicCSTbeams', f) for f in cst_files]
    beam.read_cst_beam(
        beam_files, beam_type='efield', frequency=freqs,
        telescope_name='HERA', feed_name='PAPER', feed_version='0.1', feed_pol=['x'],
        model_name='E-field pattern - Rigging height 4.9m', model_version='1.0'
    )
    beam.peak_normalize()
    return beam


@pytest.fixture(scope='session')
def hera_loc():
    return EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s', height=1073.)
