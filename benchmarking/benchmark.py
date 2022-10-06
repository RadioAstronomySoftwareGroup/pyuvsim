# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2022 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License
"""Utility functions for benchmarking."""
import os
import re
import socket

import h5py
import numpy as np
import yaml
from pyradiosky import write_healpix_hdf5

from pyuvsim.simsetup import _write_layout_csv, freq_array_to_params


def settings_setup(settings_file, outdir=None):
    """
    Parse settings file and make other variables.

    Parameters
    ----------
    settings_file : str
        Path to settings yaml file
    outdir : str
        Output file directory (optional)

    Returns
    -------
    settings : dict
        Dictionary of configuration and output file parameters.

    """
    with open(settings_file, 'r') as yfile:
        settings = yaml.safe_load(yfile)

    odir_keys = ['config_dir', 'profiles', 'data_out']
    if outdir is not None:
        for key in odir_keys:
            settings[key] = os.path.join(outdir, settings[key])

    if 'Nside' in settings.keys():
        settings['Nsrcs'] = 12 * settings['Nside']**2

    settings['hostname'] = socket.getfqdn()
    settings['settings_file'] = settings_file
    settings['teleconfig_file'] = 'tele_config.yaml'
    settings['layout_fname'] = 'layout.csv'
    settings['hpx_fname'] = 'skymodel.hdf5'
    settings['beamtype'] = 'gaussian'
    settings['beamshape'] = {'sigma' : 0.08449}
    settings['profile_path'] = os.path.join(settings['profiles'], settings['profile_prefix'])
    settings['jobscript'] = os.path.join(outdir, 'jobscript.sh')

    return settings


def make_benchmark_configuration(settings_dict):
    """
    Make configuration files and directories for benchmarking simulation.

    Creates input/output directories and configuration files for a benchmarking simulation.

    Parameters
    ----------
    settings_dict : dict
        Dictionary of parameters from benchmark.settings_setup

    """
    confdir = settings_dict['config_dir']
    outdir = settings_dict['data_out']
    Nside = settings_dict['Nside']
    teleconfig_file = settings_dict['teleconfig_file']
    beamshape = settings_dict['beamshape']
    beamtype = settings_dict['beamtype']
    layout_fname = settings_dict['layout_fname']
    hpx_fname = settings_dict['hpx_fname']

    Nfreqs = settings_dict['Nfreqs']
    Nbls = settings_dict['Nbls']
    Ntimes = settings_dict['Ntimes']

    odir_keys = ['config_dir', 'profiles', 'data_out']
    for key in odir_keys:
        odir = settings_dict[key]
        if not os.path.exists(odir):
            os.makedirs(odir)

    Nsrcs = 12 * Nside**2

    # ----------------
    # Telescope config
    # ----------------
    teleconfig_path = os.path.join(confdir, teleconfig_file)
    shapekey, shapeval = beamshape.popitem()
    teleconfig = {
        'beam_paths': {
            0: {'type': beamtype, shapekey: shapeval}
        },
        'telescope_location': '(-30.72153, 21.42831, 1073.00000)',
        'telescope_name': 'test_array'
    }

    with open(teleconfig_path, 'w') as yfile:
        yaml.dump(teleconfig, yfile, default_flow_style=False)

    # ----------------
    # Frequency setup
    # ----------------

    f0 = 1e8
    df = 195312.5
    freqs = np.arange(Nfreqs) * df + f0

    # ----------------
    # Array layout setup
    # ----------------

    max_enu = 200
    nearest_nants = np.floor((np.sqrt(1 + 8 * Nbls) + 3) / 2.) + 1
    Nants = int(nearest_nants)
    antpos_enu = np.random.uniform(-max_enu, max_enu, (Nants, 3))
    antpos_enu[:, 2] = 0.0   # Set heights to zero.
    antnums = np.arange(Nants)

    _write_layout_csv(os.path.join(confdir, layout_fname),
                      antpos_enu, antnums.astype('str'), antnums)

    blsel = []
    bi = 0
    for a1 in range(Nants):
        for a2 in range(a1, Nants):
            if bi >= Nbls:
                break
            blsel.append('({},{})'.format(a1, a2))
            bi += 1

    # ----------------
    # Sky Model setup
    #   Full frequency HEALPix shell.
    # ----------------

    fmin, fmax = 0.0, 1.0   # K (fluxes)
    skydat = np.random.uniform(fmin, fmax, (Nfreqs, Nsrcs))

    write_healpix_hdf5(os.path.join(confdir, hpx_fname), skydat, range(Nsrcs), freqs)

    # ----------------
    # Make config dictionaries
    # ----------------

    filedict = {
        'outdir': outdir,
        'outfile_name': 'benchmark',
        'output_format': 'uvh5'
    }

    freqdict = freq_array_to_params(freqs)

    srcdict = {
        'catalog': hpx_fname
    }

    teledict = {
        'array_layout': layout_fname,
        'telescope_config_name': teleconfig_file
    }

    timedict = {
        'Ntimes': Ntimes,
        'integration_time': 11.0,
        'start_time': 2458116.24485307
    }

    seldict = {
        'bls': '[' + ", ".join(blsel) + ']'
    }

    param_dict = {
        'filing': filedict,
        'freq': freqdict,
        'sources': srcdict,
        'telescope': teledict,
        'time': timedict,
        'select': seldict
    }

    obsparam_path = os.path.join(settings_dict['config_dir'], settings_dict['obsparam_name'])
    with open(obsparam_path, 'w') as yfile:
        yaml.dump(param_dict, yfile, default_flow_style=False)


def make_jobscript(settings_dict):
    """
    Write out a SLURM submittable jobscript.

    Parameters
    ----------
    settings_dict : dict
        Dictionary of parameters from benchmark.settings_setup

    """
    mem = settings_dict['MemoryLimit']
    walltime = settings_dict['walltime']
    Ncpus_per_task = settings_dict['Ncpus_per_task']
    Nnodes = settings_dict['Nnodes']
    Ntasks = settings_dict['Ntasks']
    profile_path = settings_dict['profile_path']
    obspath = os.path.join(settings_dict['config_dir'], settings_dict['obsparam_name'])

    script = "#!/bin/bash\n\n"
    script += "#SBATCH -J pyuvsim_benchmark\n"
    script += "#SBATCH --mem={}\n".format(mem)
    script += "#SBATCH --time={}\n".format(walltime)
    script += "#SBATCH --cpus-per-task={:d}\n".format(Ncpus_per_task)
    script += "#SBATCH --nodes={}-{}\n".format(Nnodes, Nnodes)
    script += "#SBATCH --ntasks={:d}\n".format(Ntasks)
    script += "#SBATCH -m cyclic\n\n"

    script += "srun --mpi=pmi2 python ../scripts/run_param_pyuvsim.py "\
              + "{} --profile='{}' --raw_profile".format(obspath, profile_path)

    with open(settings_dict['jobscript'], 'w') as jfile:
        jfile.write(script)


def update_runlog(settings_dict, logfile='BENCHMARKS.log'):
    """
    Read results from a benchmarking simulation and update a log file.

    Parameters
    ----------
    settings_dict : dict
        Dictionary of parameters from benchmark.settings_setup
    logfile : str
        Name of log file (Default: BENCHMARKS.log)

    """
    meta_file = settings_dict['profile_path'] + "_meta.out"

    data_file = os.path.join(settings_dict['data_out'], 'benchmark.uvh5')

    # Look for "pyuvsim verison: ??" but disregard any number of spaces between
    # "pyuvsim" "version" the ":" and the version itself.
    # \d+\.\d+\.\d+ ==> version ?.?.?
    # \.dev(?:\d)*\+g\w{7} ==> optional subversion string too ".dev##+g" followed by the hash
    # (\.\w+)? possibly a branch name too
    # (?:.....)? means an match this group 0 or one times (essentially making it optional)
    pattern = r"pyuvsim\s*version\s*:\s*(?P<VERSION>\d+\.\d+\.\d+(\.dev(?:\d)*\+g\w{7}(?:\.\w+)?)?)"

    with h5py.File(data_file, 'r') as dfile:
        hist = dfile['Header/history'][()].decode('UTF-8')

        match = re.search(pattern, hist)

        if match is not None:
            pyuvsim_version = match.group("VERSION")
        else:
            # try to just split it out of the history brute force.
            # diregard the last character (a '.' in a sentence).
            pyuvsim_version = hist.split("pyuvsim version: ")[1].split(" ")[0][:-1]

    header_vals = [
        "Date/Time", 'uvsim_version', 'HostName', 'SettingsFile',
        'cpus-per-task', 'Ntasks', 'Nnodes', 'MemLimit',
        'Ntimes', 'Nbls', 'Nfreqs', 'Nsrcs', 'Nsrcs_part',
        'MaxRSS [GiB]', 'Runtime'
    ]

    widths = [len(s) for s in header_vals]

    with open(meta_file, 'r') as mfile:
        lines = mfile.readlines()
    lines = (line.split() for line in lines)
    meta = {line[0]: '-'.join(line[1:]) for line in lines}

    results = [
        meta['Date/Time'],
        pyuvsim_version,
        settings_dict['hostname'],
        settings_dict['settings_file'],
        settings_dict['Ncpus_per_task'],
        settings_dict['Ntasks'],
        settings_dict['Nnodes'],
        settings_dict['MemoryLimit'],
        settings_dict['Ntimes'],
        settings_dict['Nbls'],
        settings_dict['Nfreqs'],
        settings_dict['Nsrcs'],
        meta['Nsrcs_loc'],
        meta['MaxRSS'],
        meta['Runtime'],
    ]

    results = [str(res) for res in results]
    for ii, wid in enumerate(widths):
        widths[ii] = max(len(results[ii]), wid)

    formats = ['{' + ': <{}'.format(w) + '}' for w in widths]

    header = '\t'.join([formats[i].format(hkey) for i, hkey in enumerate(header_vals)])

    if not os.path.exists(logfile):
        log = open(logfile, 'w')
        log.write(header)
    else:
        log = open(logfile, 'a')

    results = [formats[i].format(str(r)) for i, r in enumerate(results)]

    log.write('\n' + '\t'.join(results))
    log.close()
