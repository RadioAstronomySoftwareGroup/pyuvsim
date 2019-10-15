# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import os
import sys
import time as pytime

import numpy as np
try:
    import psutil
    HAVE_PSUTIL = True
except ImportError:
    HAVE_PSUTIL = False

from . import version as simversion


def get_version_string():
    version_string = ('Simulated with pyuvsim version: ' + simversion.version + '.')
    if simversion.git_hash:
        version_string += ('  Git origin: ' + simversion.git_origin
                           + '.  Git hash: ' + simversion.git_hash
                           + '.  Git branch: ' + simversion.git_branch
                           + '.  Git description: ' + simversion.git_description + '.')
    return version_string


class progsteps:
    """
    Similar to a progress bar, this prints a percentage of task completion.

    Parameters
    ----------

    maxval : int
        Maximum value to count to.

    """

    def __init__(self, maxval=None):
        self.t0 = pytime.time()
        if maxval is None:
            raise ValueError("Maximum value is needed.")
        self.maxval = float(maxval)
        step = self.maxval * 0.01
        if step < 1.0:
            step = 1
        self.step = step
        self.curval = -1
        self.remain = None

    def update(self, count):
        """
        Update the progress bar.

        Parameters
        ----------
        count : int
            Current value of counter
        """
        if count >= self.curval + self.step:
            doprint = False
            if not self.curval == count:
                doprint = True
                self.curval = count
            if doprint:
                dt = pytime.time() - self.t0
                frac_done = count / self.maxval
                self.remain = dt * (1 / frac_done - 1)
                print(("{:0.2f}% completed. {:0.3f} minutes elapsed."
                       + "{:0.3f} minutes remaining. \n").format(
                    frac_done * 100., dt / 60., self.remain / 60.))
                sys.stdout.flush()

    def finish(self):
        self.update(self.maxval)


def altaz_to_zenithangle_azimuth(altitude, azimuth):
    """
    Convert from astropy altaz convention to UVBeam az/za convention.

    Args:
        altitude: in radians
        azimuth: in radians in astropy convention: East of North (N=0, E=pi/2)

    Returns:
        zenith_angle in radians
        azimuth in radians in uvbeam convention: North of East(East=0, North=pi/2)
    """
    input_alt = np.array(altitude)
    input_az = np.array(azimuth)
    if input_alt.size != input_az.size:
        raise ValueError('number of altitude and azimuth values must match.')

    zenith_angle = np.pi / 2 - input_alt
    new_azimuth = np.pi / 2 - input_az

    if new_azimuth.size > 1:
        wh_neg = np.where(new_azimuth < -1e-9)
        if wh_neg[0].size > 0:
            new_azimuth[wh_neg] = new_azimuth[wh_neg] + np.pi * 2
    else:
        if new_azimuth < -1e-9:
            new_azimuth = new_azimuth + np.pi * 2

    return zenith_angle, new_azimuth


def zenithangle_azimuth_to_altaz(zenith_angle, azimuth):
    """
    Convert from astropy altaz convention to UVBeam az/za convention.

    Args:
        zenith_angle: in radians
        azimuth: in radians in uvbeam convention: North of East(East=0, North=pi/2)

    Returns:
        altitude in radians
        azimuth in radians in astropy convention: East of North (N=0, E=pi/2)
    """
    input_za = np.array(zenith_angle)
    input_az = np.array(azimuth)
    if input_za.size != input_az.size:
        raise ValueError('number of zenith_angle and azimuth values must match.')

    altitude = np.pi / 2 - input_za
    new_azimuth = np.pi / 2 - input_az

    if new_azimuth.size > 1:
        wh_neg = np.where(new_azimuth < -1e-9)
        if wh_neg[0].size > -1e-9:
            new_azimuth[wh_neg] = new_azimuth[wh_neg] + np.pi * 2
    else:
        if new_azimuth < -1e-9:
            new_azimuth = new_azimuth + np.pi * 2

    return altitude, new_azimuth


def strip_extension(filepath, ext=None):
    """ Remove extension from file. """
    if '.' not in filepath:
        return filepath, ''
    file_list = filepath.split('.')
    if ext is not None:
        return filepath[:-len(ext) - 1], '.' + ext
    ext = file_list[-1]
    # miriad files might not have an extension
    # limited list of recognized extensions
    if ext not in ['uvfits', 'uvh5', 'yaml']:
        return filepath, ''
    return ".".join(file_list[:-1]), '.' + file_list[-1]


def check_file_exists_and_increment(filepath, extension=None):
    """
        Given filepath (path + filename), check if it exists. If so, add a _1
        at the end, if that exists add a _2, and so on.
    """
    base_filepath, ext = strip_extension(filepath, extension)
    bf_list = base_filepath.split('_')
    if bf_list[-1].isdigit():
        base_filepath = '_'.join(bf_list[:-1])
    n = 0
    while os.path.exists(filepath):
        filepath = "{}_{}".format(base_filepath, n) + ext
        n += 1
    return filepath


def write_uvdata(uv_obj, param_dict, return_filename=False, dryrun=False, out_format=None):
    """
    Parse output file information from parameters and write uvfits to file.

    Args:
        uv_obj: UVData object to write out.
        param_dict: parameter dictionary defining output path, filename, and
                    whether or not to clobber.
        return_filename: (Default false) Return the file path
        dryrun: (Default false) Don't write to file.
        out_format: (Default uvfits) Write as uvfits/miriad/uvh5

    Returns:
        File path, if return_filename is True
    """
    if 'filing' in param_dict.keys():
        param_dict = param_dict['filing']
    if 'outdir' not in param_dict:
        param_dict['outdir'] = '.'
    if 'output_format' in param_dict:
        out_format = param_dict['output_format']
    elif out_format is None:
        out_format = 'uvfits'

    if 'outfile_name' not in param_dict or param_dict['outfile_name'] == '':
        outfile_prefix = ""
        outfile_suffix = "results"
        if 'outfile_prefix' in param_dict:
            outfile_prefix = param_dict['outfile_prefix']
        if 'outfile_suffix' in param_dict:
            outfile_suffix = param_dict['outfile_suffix']
        outfile_name = "_".join([outfile_prefix, outfile_suffix])
        outfile_name = os.path.join(param_dict['outdir'], outfile_name)
    else:
        outfile_name = os.path.join(param_dict['outdir'], param_dict['outfile_name'])

    if not os.path.exists(param_dict['outdir']):
        os.makedirs(param_dict['outdir'])

    if out_format == 'uvfits':
        if not outfile_name.endswith(".uvfits"):
            outfile_name = outfile_name + ".uvfits"

    if out_format == 'uvh5':
        if not outfile_name.endswith(".uvh5"):
            outfile_name = outfile_name + ".uvh5"

    noclobber = ('clobber' not in param_dict) or not bool(param_dict['clobber'])
    if noclobber:
        outfile_name = check_file_exists_and_increment(outfile_name)

    print('Outfile path: ', outfile_name)
    if not dryrun:
        if out_format == 'uvfits':
            uv_obj.write_uvfits(outfile_name, force_phase=True, spoof_nonessential=True)
        elif out_format == 'miriad':
            uv_obj.write_miriad(outfile_name, clobber=not noclobber)
        elif out_format == 'uvh5':
            uv_obj.write_uvh5(outfile_name)
        else:
            raise ValueError(
                "Invalid output format. Options are \" uvfits\", \"uvh5\", or \"miriad\"")
    if return_filename:
        return outfile_name


def get_avail_memory():
    """
    Method for estimating the virtual memory available (in bytes)
    on the current node to a running process.

    Currently only supports the SLURM array scheduler.

    If this is not called from within a SLURM task, it will estimate
    using psutils methods.
    """
    if not HAVE_PSUTIL:
        raise ImportError("You need psutils to estimate available memory. "
                          "Install it by running pip install pyuvsim[sim] "
                          "or pip install pyuvsim[all] if you also want "
                          "h5py and line_profiler installed.")

    slurm_key = 'SLURM_MEM_PER_NODE'
    if slurm_key in os.environ:
        return float(os.environ[slurm_key]) * 1e6  # MB -> B

    return psutil.virtual_memory().available


def iter_array_split(part_index, N, M):
    """
    Returns an iterator giving the indices of `part` below:
        part = np.array_split(np.arange(N), M)[part_index]

    This mimics the behavior of array_split without having to make
    the whole array that will be split.
    """

    Neach_section, extras = divmod(N, M)
    if part_index < extras:
        length = Neach_section + 1
        start = part_index * (length)
        end = start + length
    else:
        length = Neach_section
        start = extras * (Neach_section + 1) + (part_index - extras) * length
        end = start + length

    return range(start, end), end - start


def estimate_skymodel_memory_usage(Ncomponents, Nfreqs):
    """
    Estimate the memory footprint of a SkyModel by summing the sizes
    of the data types that go into SkyModel.

    This aims to anticipate the full memory required to handle a SkyModel
    class in simulation, accounting for its attributes as well as
    data generated while used.

    Parameters
    ----------

    Ncomponents : int
        Number of source components.
    Nfreqs : int
        Number of frequencies per source component.
        (size of SkyModel.freq_array)

    Returns
    -------
    mem_est : float
        Estimate of memory usage in bytes
    """

    base_float = [1.5]    # A float
    base_bool = [True]
    base_str = ["source_name"]

    Ncomp_attrs = {'ra': base_float, 'dec': base_float,
                   'alt_az': 2 * base_float, 'rise_lst': base_float, 'set_lst': base_float,
                   'pos_lmn': 3 * base_float, 'name': base_str, 'horizon_mask': base_bool}
    Ncomp_Nfreq_attrs = {'stokes': 4 * base_float,
                         'coherency_radec': 4 * base_float,
                         'coherency_local': 4 * base_float}

    mem_est = np.sum([sys.getsizeof(v) * Ncomponents for k, v in Ncomp_attrs.items()])
    mem_est += np.sum([sys.getsizeof(v) * Ncomponents * Nfreqs
                       for k, v in Ncomp_Nfreq_attrs.items()])
    return mem_est
