#!/usr/bin/env python
# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License
"""Run a pyuvsim simulation from a parameter file."""

import argparse
import os
import time as pytime
from datetime import datetime, timedelta

import pyuvsim

parser = argparse.ArgumentParser(
    description="A command-line script to execute a pyuvsim simulation from a parameter file."
)
parser.add_argument("paramsfile", type=str, help="Parameter yaml file.", default=None)
parser.add_argument("--profile", type=str, help="Time profiling output file name.")
parser.add_argument("--quiet", action="store_true", help="Suppress stdout printing.")
parser.add_argument(
    "--keep_nonroot_stdout",
    action="store_true",
    help="Do not redirect stdout on nonzero ranks to /dev/null.",
)
parser.add_argument(
    "--raw_profile",
    help="Also save pickled LineStats data for line profiling.",
    action="store_true",
)
parser.add_argument(
    "--backend",
    type=str,
    help="Backend task collection for simulation",
    choices=["rma", "send_recv"],
    default="rma",
)

parser.add_argument(
    "--progbar",
    type=str,
    help="Monitor and reporting module for simulation progress",
    choices=["progsteps", "tqdm"],
    default="progsteps",
)

args = parser.parse_args()

if args.paramsfile is None:
    raise ValueError("Parameter file required")

if args.profile is not None:
    pyuvsim.profiling.set_profiler(
        outfile_prefix=args.profile, dump_raw=args.raw_profile
    )

if not os.path.isdir(os.path.dirname(args.paramsfile)):
    args.paramsfile = os.path.join(".", args.paramsfile)

t0 = pytime.time()

block_nonroot_stdout = not args.keep_nonroot_stdout

pyuvsim.uvsim.run_uvsim(
    args.paramsfile,
    quiet=args.quiet,
    block_nonroot_stdout=block_nonroot_stdout,
    backend=args.backend,
    progbar=args.progbar,
)

if args.profile:
    dt = pytime.time() - t0
    maxrss = pyuvsim.mpi.get_max_node_rss()
    rtime = str(timedelta(seconds=dt))
    if isinstance(maxrss, float):
        print("\tRuntime: {} \n\tMaxRSS: {:.3f} GiB".format(rtime, maxrss))
    if hasattr(pyuvsim.profiling.prof, "meta_file"):
        with open(pyuvsim.profiling.prof.meta_file, "a") as afile:
            afile.write("Runtime \t {}\nMaxRSS \t {:.3f}\n".format(rtime, maxrss))
            afile.write("Date/Time \t {}".format(str(datetime.now())))
