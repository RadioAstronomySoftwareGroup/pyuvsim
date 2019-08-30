# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from __future__ import absolute_import, division, print_function

from . import version
from .simsetup import *  # noqa
from .analyticbeam import *  # noqa
from .antenna import *  # noqa
from .baseline import *  # noqa
from .source import *  # noqa
from .telescope import *  # noqa
from .utils import *  # noqa

try:
    from .uvsim import *  # noqa
except ImportError:
    # can only be used in setup mode
    pass

try:
    from .profiling import *  # noqa
except ImportError:
    # can only be used in non-profiling mode.
    pass

__version__ = version.version
