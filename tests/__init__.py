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
                f"key: {key} has type {type(d1[key])} in d1 and {type(d2[key])} in d2\n"
                f"d1:  data has type {type(d1[key][0])} and value {d1[key]}\n"
                f"d2:  data has type {type(d2[key][0])} and value {d2[key]}\n"
            )

        elif isinstance(d1[key], (np.ndarray)):
            if np.issubdtype(d1[key].dtype, np.bytes_):
                assert np.array_equal(d1[key], np.asarray(d2[key]))
            else:
                assert np.allclose(d1[key], np.asarray(d2[key])), (
                    f"key: {key} has type {type(d1[key])} in d1 and {type(d2[key])} in d2\n"
                    f"d1:  data has type {type(d1[key][0])} and value {d1[key]}\n"
                    f"d2:  data has type {type(d2[key][0])} and value {d2[key]}\n"
                )
        elif isinstance(d1[key], dict):
            compare_dictionaries(d1[key], d2[key])
        elif isinstance(d1[key], float | np.float32):
            assert np.allclose(d1[key], d2[key]), (
                f"key: {key} has type {type(d1[key])} in d1 and {type(d2[key])} in d2\n"
                f"d1:  data has type {type(d1[key][0])} and value {d1[key]}\n"
                f"d2:  data has type {type(d2[key][0])} and value {d2[key]}\n"
            )
        else:
            assert d1[key] == d2[key], (
                f"key: {key} has type {type(d1[key])} in d1 and {type(d2[key])} in d2\n"
                f"d1:  data has type {type(d1[key][0])} and value {d1[key]}\n"
                f"d2:  data has type {type(d2[key][0])} and value {d2[key]}\n"
            )
    return True


def get_vizier_mirror():
    # Check to see if Vizier is avaliable, if not check if CfA mirror is avaliable
    import requests  # a dependency of astroquery

    default_vizier_mirror = "vizier.u-strasbg.fr"
    vizier_mirror = None
    path = "/viz-bin/votable"
    try:
        response = requests.head(
            f"https://{default_vizier_mirror}{path}", timeout=1.0, allow_redirects=True
        )
        response.raise_for_status()
    except (
        requests.exceptions.RequestException,
        requests.exceptions.ConnectionError,
        requests.exceptions.HTTPError,
    ):
        vizier_mirror = "vizier.cfa.harvard.edu"
        try:
            response = requests.head(
                f"https://{vizier_mirror}{path}", timeout=1.0, allow_redirects=True
            )
            response.raise_for_status()
        except (
            requests.exceptions.RequestException,
            requests.exceptions.ConnectionError,
            requests.exceptions.HTTPError,
        ) as e:
            raise requests.exceptions.RequestException(
                "Connection error w/ Vizier"
            ) from e

    return vizier_mirror
