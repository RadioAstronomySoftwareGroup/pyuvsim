#!/bin/python
# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2020 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License
"""
Compare a run of a reference sim with the reference run.

Compare the visibilities in a given set of simulation outputs to the corresponding
values in a reference set.
"""

import argparse
import glob
import os
import warnings

from pyuvdata import UVData

parser = argparse.ArgumentParser(
    description="A script to compare new data files (uvh5) to downloaded archived data."
)
parser.add_argument('paths', nargs='+', type=str, help="Paths to new files, "
                                                       "to compare with archived data.")
parser.add_argument('-g', '--generation', type=int, help="Generation of reference sims "
                                                         "to download (1 or 2).", default=1)

args = parser.parse_args()

ref_set_path = os.path.join('latest_ref_data', f"gen{args.generation}")

new_file_paths = args.paths
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
    uv_new.history = uv_old.history
    print(k, uv_old == uv_new)
