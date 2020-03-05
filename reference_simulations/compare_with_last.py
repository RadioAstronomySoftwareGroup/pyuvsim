#!/bin/python

# Compare the visibilities in a given set of simulation outputs to the corresponding
# values in a reference set.

import numpy as np
import sys
import glob
import warnings
import os

import h5py

# Loading data using h5py directly instead of reading the full uvh5 files.
#   (The 1.2 reference sims recalculate all lsts when reloaded, which takes a long time.)

ref_set_path = 'latest_ref_data/v1'

new_file_paths = sys.argv[1:]
old_file_paths = glob.glob(os.path.join(ref_set_path, '*.uvh5'))

new_file_names = [os.path.basename(fn) for fn in new_file_paths]
old_file_names = [os.path.basename(fn) for fn in old_file_paths]

new_files = dict(zip(new_file_names, new_file_paths))
old_files = dict(zip(old_file_names, old_file_paths))

for k in new_files.keys():
    if k not in old_files.keys():
        warnings.warn("File {} not in reference data.".format(k))
        continue
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')     # Ignore known_telescopes warning
        uv_old = h5py.File(old_files[k], 'r')['Data']
        uv_new = h5py.File(new_files[k], 'r')['Data']

    print(k, np.allclose(uv_old['visdata'][()], uv_new['visdata'][()]))
