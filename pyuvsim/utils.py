# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 2-clause BSD License

import time as pytime
import sys


class progsteps:
    """
        Similar to progress bar. Prints a percentage of completion.
        For when running in batch and progress bar doesn't work well.
    """
    def __init__(self, maxval=None):
        self.t0 = pytime.time()
        if maxval is None:
            raise ValueError("Maximum value is needed.")
        self.maxval = float(maxval)
        step = self.maxval * 0.01
        if step < 1.0:
            step = 1
        self.step = step

    def update(self, count):
        if count % self.step == 0:
            print("{:0.2f}% completed. {:0.3f} minutes elapsed \n".format(
                  (count / self.maxval) * 100., (pytime.time() - self.t0) / 60.))
            sys.stdout.flush()

    def finish(self):
        self.update(self.maxval)
