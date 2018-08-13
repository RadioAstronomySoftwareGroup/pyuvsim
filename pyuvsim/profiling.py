# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 2-clause BSD License

"""
Use the line profiler when requested.

If the set_profiler method is called, the profile decorators
throughout the code will do line profiling on all instance methods.
"""

from __future__ import absolute_import, division, print_function

from line_profiler import LineProfiler
import inspect
import sys
import atexit

PY3 = sys.version_info[0] == 3

if PY3:
    import builtins
else:
    import __builtin__ as builtins

prof = LineProfiler()


def profile(cls):
    """ cls can be a class or function."""
    if inspect.isfunction(cls):
        return prof(cls)

    class wrapper(object):
        def __init__(self, *args, **kwargs):
            self.oInstance = cls(*args, **kwargs)

        def __getattribute__(self, s):
            try:
                x = super(wrapper, self).__getattribute__(s)
            except AttributeError:
                pass
            else:
                return x
            x = self.oInstance.__getattribute__(s)
#            if type(x) == type(self.__init__):  # Instance method
#                return prof(x)
#            else:
#                return x
            return prof(x)
    return wrapper


def set_profiler():
    """ If profiling is requested, then assign it to the builtins """
    builtins.__dict__['profile'] = profile
    atexit.register(prof.print_stats)


# By default, the profile decorator has no effect.
if 'profile' not in builtins.__dict__:
    builtins.__dict__['profile'] = lambda f: f
else:
    set_profiler()  # Will activate if run with kernprof
