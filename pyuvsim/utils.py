import time as pytime
import sys

class progsteps:
    """
        Similar to progress bar. Prints a percentage of completion.
        For when running in batch and progress bar doesn't work well.
    """
    def __init__(self, maxval=None):
        self.t0 = time.time()
        if maxval is None:
            raise ValueError("Maximum value is needed.")
        self.maxval = float(maxval)
        step = self.maxval*0.01
        if step < 1.0: step = 1
        self.step = step
        

    def update(self, count):
        if count % self.step == 0:
            print("".join(map(str,[(count/float(tot)) * 100., "% completed. Elapsed: ", (pytime.time() - time0)/60., 'minutes \n'])))
            sys.stdout.flush()
