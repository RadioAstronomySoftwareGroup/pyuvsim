# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License
"""Define namespace."""
from setuptools_scm import get_version
from pathlib import Path
from pkg_resources import get_distribution, DistributionNotFound
import warnings

from .branch_scheme import branch_scheme


try:  # pragma: nocover
    # get accurate version for developer installs
    version_str = get_version(Path(__file__).parent.parent, local_scheme=branch_scheme)

    __version__ = version_str

except (LookupError, ImportError):
    try:
        # Set the version automatically from the package details.
        __version__ = get_distribution(__name__).version
    except DistributionNotFound:  # pragma: nocover
        # package is not installed
        pass

# Filter distutils Deprecation Warning from pyuvdata
# (can be removed when we require pyuvdata > 2.2.6)
# needs to be done before the imports to work properly
warnings.filterwarnings("ignore", message="distutils Version classes are deprecated")

from .profiling import *  # noqa
from .uvsim import *  # noqa
from .simsetup import *  # noqa
from .analyticbeam import *  # noqa
from .antenna import *  # noqa
from .baseline import *  # noqa
from .telescope import *  # noqa
from .utils import *  # noqa
