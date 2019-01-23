# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from __future__ import absolute_import, division, print_function

import argparse

from pyuvdata import UVData

import pyuvsim.simsetup

# Take a uvfits file, and save sim parameters as a yaml file.

parser = argparse.ArgumentParser(description=("Generate basic simulation parameters from uvfits."))

parser.add_argument('file_in', metavar='<FILE>', type=str, nargs='+')
parser.add_argument('-p', '--param_filename', default=None)
parser.add_argument('-t', '--telescope_config_path', default='')
parser.add_argument('-l', '--layout_csv_path', default='')

args = parser.parse_args()

uvd = UVData()
uvd.read(args.file_in[0])

pyuvsim.simsetup.uvdata_to_config_file(uvd, param_filename=args.param_filename,
                                       telescope_config_name=args.telescope_config_path,
                                       layout_csv_name=args.layout_csv_path,
                                       catalog='mock', path_out='.')
