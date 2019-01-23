# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from __future__ import absolute_import, division, print_function

import argparse
import os

from pyuvdata import UVData

from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
import pyuvsim.simsetup

# Read in a uvfits file and automatically generate yaml files.
# Will assume the same beam_id for all antennas for now.

parser = argparse.ArgumentParser(description=("Extracts antenna position info from uvfits."))

parser.add_argument('file_in', metavar='<FILE>', type=str, nargs='+')
parser.add_argument('-l', '--layout_csv_name', default=None)
parser.add_argument('-t', '--telescope_config_name', default=None)
parser.add_argument('-b', '--beam_filepath', type=str, default=os.path.join(SIM_DATA_PATH, 'HERA_NicCST.uvbeam'))

args = parser.parse_args()

uvd = UVData()
uvd.read(args.file_in[0])

pyuvsim.simsetup.uvdata_to_telescope_config(uvd, args.beam_filepath, layout_csv_name=args.layout_csv_name,
                                            telescope_config_name=args.telescope_config_name,
                                            return_names=False, path_out='.')
