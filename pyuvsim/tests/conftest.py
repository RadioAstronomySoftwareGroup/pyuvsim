# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2019 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

"""Testing environment setup and teardown for pytest."""
from __future__ import absolute_import, division, print_function

import os
import pytest
import shutil
from pyuvsim.data import DATA_PATH


@pytest.fixture(autouse=True, scope="session")
def setup_and_teardown_package():
    """Make data/test directory to put test output files in."""
    testdir = os.path.join(DATA_PATH, 'temporary_test_data/')
    if not os.path.exists(testdir):
        print('making test directory')
        os.mkdir(testdir)

    # yield to allow tests to run
    yield

    # clean up the test directory after
    shutil.rmtree(testdir)
