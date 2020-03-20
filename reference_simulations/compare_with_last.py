#!/bin/python

# Compare the visibilities in a given set of simulation outputs to the corresponding
# values in a reference set.

import sys
import glob
import warnings
import os

from pyuvdata import UVData

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
    uv_old = UVData()
    uv_new = UVData()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')     # Ignore known_telescopes warning
        uv_old.read_uvh5(old_files[k], run_check_acceptability=False)
        uv_new.read_uvh5(new_files[k], run_check_acceptability=False)

    print(k, uv_old == uv_new)
