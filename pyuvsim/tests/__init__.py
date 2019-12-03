# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import os

import numpy as np
import pytest
from pyuvdata import UVBeam
from pyuvdata.data import DATA_PATH

from pyuvsim.data import DATA_PATH as SIM_DATA_PATH

TESTDATA_PATH = os.path.join(SIM_DATA_PATH, 'temporary_test_data/')


# functions used by many tests


def compare_dictionaries(d1, d2):
    """Recursively compare dictionaries.

    Keys of each dict must match.
    Walks through two input dicts and compares each key.
    Makes calls to assert statements and np.allclose to compare values.
    """
    assert set(d1.keys()) == set(d2.keys())
    for key in d1:
        if isinstance(d1[key], (list)):
            assert d1[key] == list(d2[key]), (
                "key: {key} has type {key1_type} in d1 and {key2_type} in d2\n"
                "d1:  data has type {data1_type} and value {data1_val}\n"
                "d2:  data has type {data2_type} and value {data2_val}\n".format(
                    key=key,
                    key1_type=type(d1[key]),
                    key2_type=type(d2[key]),
                    data1_type=type(d1[key][0]),
                    data1_val=d1[key],
                    data2_type=type(d2[key][0]),
                    data2_val=d2[key]
                )
            )

        elif isinstance(d1[key], (np.ndarray)):
            if np.issubdtype(d1[key].dtype, np.string_):
                assert np.array_equal(d1[key], np.asarray(d2[key]))
            else:
                assert np.allclose(d1[key], np.asarray(d2[key])), (
                    "key: {key} has type {key1_type} in d1 and {key2_type} in d2\n"
                    "d1:  data has type {data1_type} and value {data1_val}\n"
                    "d2:  data has type {data2_type} and value {data2_val}\n".format(
                        key=key,
                        key1_type=type(d1[key]),
                        key2_type=type(d2[key]),
                        data1_type=type(d1[key][0]),
                        data1_val=d1[key],
                        data2_type=type(d2[key][0]),
                        data2_val=d2[key]
                    )
                )
        elif isinstance(d1[key], dict):
            compare_dictionaries(d1[key], d2[key])
        elif isinstance(d1[key], (float, np.float, np.float32)):
            assert np.allclose(d1[key], d2[key]), (
                "key: {key} has type {key1_type} in d1 and {key2_type} in d2\n"
                "d1:  data has type {data1_type} and value {data1_val}\n"
                "d2:  data has type {data2_type} and value {data2_val}\n".format(
                    key=key,
                    key1_type=type(d1[key]),
                    key2_type=type(d2[key]),
                    data1_type=type(d1[key][0]),
                    data1_val=d1[key],
                    data2_type=type(d2[key][0]),
                    data2_val=d2[key]
                )
            )
        else:
            assert d1[key] == d2[key], (
                "key: {key} has type {key1_type} in d1 and {key2_type} in d2\n"
                "d1:  data has type {data1_type} and value {data1_val}\n"
                "d2:  data has type {data2_type} and value {data2_val}\n".format(
                    key=key,
                    key1_type=type(d1[key]),
                    key2_type=type(d2[key]),
                    data1_type=type(d1[key][0]),
                    data1_val=d1[key],
                    data2_type=type(d2[key][0]),
                    data2_val=d2[key]
                )
            )
    return True


def make_cst_beams(freqs=None):
    beam = UVBeam()
    beam.freq_interp_kind = 'linear'

    if freqs is None:
        freqs = [150e6, 123e6]

    cst_files = ['HERA_NicCST_150MHz.txt', 'HERA_NicCST_123MHz.txt']
    beam_files = [os.path.join(DATA_PATH, 'NicCSTbeams', f) for f in cst_files]
    beam.read_cst_beam(
        beam_files, beam_type='efield', frequency=freqs,
        telescope_name='HERA', feed_name='PAPER', feed_version='0.1', feed_pol=['x'],
        model_name='E-field pattern - Rigging height 4.9m', model_version='1.0'
    )

    return beam


def assert_raises_message(exception_type, message, func, *args, **kwargs):
    """
    Check that the correct error message is raised.
    """
    nocatch = kwargs.pop('nocatch', False)
    if nocatch:
        func(*args, **kwargs)

    with pytest.raises(exception_type) as err:
        func(*args, **kwargs)

    try:
        assert message in str(err.value)
    except AssertionError as excp:
        print("{} not in {}".format(message, str(err.value)))
        raise excp
