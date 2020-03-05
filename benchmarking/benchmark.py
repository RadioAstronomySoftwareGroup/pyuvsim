import numpy as np
import yaml
import os
from pyradiosky import write_healpix_hdf5

from pyuvsim.simsetup import _write_layout_csv, freq_array_to_params


def make_benchmark_configuration(config_dir=None, data_out=None, profiles=None, obsparam_name=None,
                                 Nbls=None, Nfreqs=None, Ntimes=None, Nside=None):
    """
    Setup configuration files and directories for benchmarking simulation.

    """

    confdir = config_dir
    outdir = data_out
    profdir = profiles

    for odir in [confdir, outdir, profdir]:
        if not os.path.exists(odir):
            os.makedirs(odir)

    Nsrcs = 12 * Nside**2

    # ----------------
    # Copy telescope config
    # ----------------
    teleconfig_file = 'benchmark_tele_config.yaml'
    teleconfig_path = os.path.join(confdir, teleconfig_file)

    beamtype = 'gaussian'
    beamshape = "sigma=0.08449"

    teleconfig = {
        'beam_paths': {
            0: "{}, {}".format(beamtype, beamshape)
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

    layout_fname = 'benchmark_layout.csv'

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

    hpx_fname = 'benchmark_skymodel.hdf5'
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

    obspath = os.path.join(confdir, obsparam_name)

    with open(obspath, 'w') as yfile:
        print(obspath)
        yaml.dump(param_dict, yfile, default_flow_style=False)

    # Append a string giving the axis sizes:
    with open(obspath, 'a') as ofile:
        ofile.write(
            "#Ntimes={:d}, Nfreqs={:d}, Nbls={:d}, Nsrcs={:d}\n".format(
                Ntimes, Nfreqs, Nbls, Nsrcs
            )
        )

    return obspath


def make_jobscript(obspath, profile_path, Ntasks, Nnodes, Ncpus_per_task, walltime, mem):
    """Write out a submittable jobscript."""

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

    with open("jobscript.sh", 'w') as jfile:
        jfile.write(script)
