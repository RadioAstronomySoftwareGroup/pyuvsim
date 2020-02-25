#!/usr/bin/env python
# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import argparse
import os

import pyuvsim

parser = argparse.ArgumentParser(
    description="A command-line script to execute a pyuvsim simulation from a parameter file."
)
parser.add_argument('paramsfile', type=str, help='Parameter yaml file.', default=None)
parser.add_argument('--profile', type=str, help='Time profiling output file name.')
parser.add_argument('--raw_profile', help='Also save pickled LineStats data for line profiling.', action='store_true')

args = parser.parse_args()

if args.paramsfile is None:
    raise ValueError("Parameter file required")

if args.profile is not None:
    pyuvsim.profiling.set_profiler(outfile_prefix=args.profile, dump_raw=args.raw_profile)

if not os.path.isdir(os.path.dirname(args.paramsfile)):
    args.paramsfile = os.path.join('.', args.paramsfile)
pyuvsim.uvsim.run_uvsim(args.paramsfile)
