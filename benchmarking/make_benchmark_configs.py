import numpy as np
import yaml
import os
from shutil import copyfile
from pyradiosky import write_healpix_hdf5

from pyuvsim.simsetup import _write_layout_csv, freq_array_to_params

# To be set by calling function:

confdir = 'configdir'
outdir = 'simfiles'

obsparam_name = 'obsparam_benchmark.yaml'

if not os.path.exists(confdir):
    os.mkdir(confdir)


# Nbls = 150
# Nfreqs = 256
# Ntimes = 10
# Nside = 128
Nbls = 100
Nfreqs = 100
Ntimes = 70
Nside = 16


Nsrcs = 12 * Nside**2

# ----------------
# Copy telescope config
# ----------------
teleconfig_file = 'benchmark_analytic_config.yaml'
copyfile(teleconfig_file, os.path.join(confdir, teleconfig_file))


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

_write_layout_csv(os.path.join(confdir, layout_fname), antpos_enu, antnums.astype('str'), antnums)


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
    yaml.dump(param_dict, yfile, default_flow_style=False)

# Append a string giving the axis sizes:
with open(obspath, 'a') as ofile:
    ofile.write(
        "#Ntimes={:d}, Nfreqs={:d}, Nbls={:d}, Nsrcs={:d}\n".format(
            Ntimes, Nfreqs, Nbls, Nsrcs
        )
    )
print(obspath)
