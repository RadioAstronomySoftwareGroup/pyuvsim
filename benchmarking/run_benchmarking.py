#!/bin/env python
# Script to setup/run benchmarking jobs.

import argparse
import os
from datetime import datetime
from subprocess import check_call

from benchmark import make_benchmark_configuration, make_jobscript, settings_setup, update_runlog

now = datetime.now()

parser = argparse.ArgumentParser(
    description="A command-line script to run a benchmarking pyuvsim job."
)
parser.add_argument('yaml_file', type=str, help='Settings file.', default='settings.yaml')
parser.add_argument('--init_configs', action='store_true',
                    help='Make simulation config files (obsparam, etc.).')
parser.add_argument('--submit', action='store_true', help='Make jobscript and submit.')
parser.add_argument('--jobscript', action='store_true', help='Make jobscript only.')
parser.add_argument('--update_log', action='store_true',
                    help='Update the log of benchmarking runs.')
parser.add_argument('-o', '--outdir', type=str, help='Path for all configuration/results to go.',
                    default='{:4d}-{:02d}-{:02d}'.format(now.year, now.month, now.day))


args = parser.parse_args()

if args.submit:
    args.jobscript = True

if args.submit and args.update_log:
    raise ValueError("Cannot submit benchmarking job and update the log at the same time."
                     " Run the job first, then run again with ``--update_log`` when it's finished.")

settings = settings_setup(args.yaml_file, outdir=args.outdir)

obspath = os.path.join(args.outdir, settings['config_dir'], settings['obsparam_name'])
if args.submit and (not os.path.exists(obspath)):
    print("Existing simulation config not found. Making them now.")
    args.init_configs = True

if args.init_configs:
    make_benchmark_configuration(settings)

if args.jobscript:
    # Create a SLURM jobscript and submit it.
    make_jobscript(settings)

if args.submit:
    output = check_call(['sbatch', '-o', os.path.join(args.outdir, 'slurm_%A.out'),
                        settings['jobscript']])

if args.update_log:
    # Can only be run after the job is complete.
    update_runlog(settings, args.yaml_file)
