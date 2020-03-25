# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2019 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

"""Testing environment setup and teardown for pytest."""
import os
import shutil
import warnings

import pytest
import urllib
from astropy.time import Time
from astropy.utils import iers

from pyuvsim.data import DATA_PATH


@pytest.fixture(autouse=True, scope="session")
def setup_and_teardown_package():
    """Make data/test directory to put test output files in."""
    testdir = os.path.join(DATA_PATH, 'temporary_test_data/')
    if not os.path.exists(testdir):
        os.mkdir(testdir)

    # try to download the iers table. If it fails, turn off auto downloading for the tests
    # and turn it back on in teardown_package (done by extending auto_max_age)
    try:
        iers.IERS_A.open(iers.IERS_A_URL)
        t1 = Time.now()
        t1.ut1
    except urllib.error.URLError:
        iers.conf.auto_max_age = None

    # yield to allow tests to run
    yield

    iers.conf.auto_max_age = 30

    # clean up the test directory after
    if os.path.exists(testdir):
        shutil.rmtree(testdir)


@pytest.fixture(autouse=True)
def ignore_deprecation():
    warnings.filterwarnings('ignore', message='Achromatic gaussian beams will not be supported',
                            category=PendingDeprecationWarning)
