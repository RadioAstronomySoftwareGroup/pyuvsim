# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License
"""Define namespace."""
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from setuptools_scm import get_version

from .branch_scheme import branch_scheme

try:
    # get accurate version for developer installs
    version_str = get_version(Path(__file__).parent.parent, local_scheme=branch_scheme)

    __version__ = version_str

except (LookupError, ImportError):
    try:
        # Set the version automatically from the package details.
        __version__ = version("pyuvsim")
    except PackageNotFoundError:  # pragma: nocover
        # package is not installed
        pass

from .analyticbeam import *  # noqa
from .antenna import *  # noqa
from .baseline import *  # noqa
from .profiling import *  # noqa
from .simsetup import *  # noqa
from .telescope import *  # noqa
from .utils import *  # noqa
from .uvsim import *  # noqa
