# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import numpy as np

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
                    data2_val=d2[key],
                )
            )

        elif isinstance(d1[key], (np.ndarray)):
            if np.issubdtype(d1[key].dtype, np.bytes_):
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
                        data2_val=d2[key],
                    )
                )
        elif isinstance(d1[key], dict):
            compare_dictionaries(d1[key], d2[key])
        elif isinstance(d1[key], (float, np.float32)):
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
                    data2_val=d2[key],
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
                    data2_val=d2[key],
                )
            )
    return True
