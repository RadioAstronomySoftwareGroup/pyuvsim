
# Script to setup/run benchmarking jobs.

import argparse
import yaml
import os
import shutil
from datetime import datetime
from subprocess import check_call

from benchmark import make_benchmark_configuration, make_jobscript

now = datetime.now()

parser = argparse.ArgumentParser(
    description="A command-line script to run a benchmarking pyuvsim job."
)
parser.add_argument('yaml_file', type=str, help='Settings file.', default='settings.yaml')
parser.add_argument('--init_configs', action='store_true',
                    help='Make simulation config files (obsparam, etc.).')
parser.add_argument('--submit', action='store_true', help='Make jobscript and submit.')
parser.add_argument('--jobscript', action='store_true', help='Make jobscript only.')
parser.add_argument('--cleanup', action='store_true', help='Delete configuration and output files.')
parser.add_argument('-o', '--outdir', type=str, help='Path for all configuration/results to go.',
                    default='{:4d}-{:02d}-{:02d}'.format(now.year, now.month, now.day))


args = parser.parse_args()

args.init_configs = True

if args.submit:
    args.jobscript = True

with open(args.yaml_file, 'r') as yfile:
    settings = yaml.safe_load(yfile)

for key in ['config_dir', 'profiles', 'data_out']:
    settings[key] = os.path.join(args.outdir, settings[key])

if args.init_configs:
    obspath = make_benchmark_configuration(
        settings['config_dir'],
        settings['data_out'],
        settings['profiles'],
        settings['obsparam_name'],
        settings['Nbls'],
        settings['Nfreqs'],
        settings['Ntimes'],
        settings['Nside']
    )

else:
    obspath = os.path.join(settings['config_dir'], settings['obsparam_name'])

if args.jobscript:
    # Create a SLURM jobscript and submit it.
    profile_path = os.path.join(settings['profiles'], settings['profile_prefix'])
    make_jobscript(obspath, profile_path,
                   settings['Ntasks'],
                   settings['Nnodes'],
                   settings['Ncpus_per_task'],
                   settings['walltime'],
                   settings['memory'])

if args.submit:
    output = check_call(['sbatch', '-o', os.path.join(args.outdir, 'slurm_%A.out'), 'jobscript.sh'])


if args.cleanup:
    # Clear all outputs.
    shutil.rmtree(settings['config_dir'])
    shutil.rmtree(settings['data_out'])
    shutil.rmtree(settings['profiles'])
