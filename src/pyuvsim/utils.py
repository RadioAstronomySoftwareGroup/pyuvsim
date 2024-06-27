# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License
"""Define various utility functions."""
import itertools
import os
import sys
import time as pytime
import warnings
from datetime import timedelta

import numpy as np
import psutil

from . import __version__


def get_version_string():
    """Get the current version string for pyuvsim."""
    return "Simulated with pyuvsim version: " + __version__ + "."


class progsteps:  # noqa This should be named with CapWords convention
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
        self.curval = 0
        self.remain = None

    def __enter__(self):
        """Enter a context manager."""
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Exit a context manager."""
        self.finish()

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
                print(
                    ("{:0.2f}% completed. {}  elapsed. " + "{} remaining. \n").format(
                        frac_done * 100.0,
                        str(timedelta(seconds=dt)),
                        str(timedelta(seconds=self.remain)),
                    ),
                    flush=True,
                )

    def finish(self):
        """Finalize the progress steps info."""
        self.update(self.maxval)


def altaz_to_zenithangle_azimuth(altitude, azimuth):
    """
    Convert from astropy altaz convention to UVBeam az/za convention.

    Parameters
    ----------
    altitude, azimuth: float or array of float
        altitude above horizon
        azimuth in radians in astropy convention: East of North (N=0, E=pi/2)

    Returns
    -------
    zenith_angle: float or array of float
        In radians
    azimuth: float or array of float
        In radians in uvbeam convention: North of East(East=0, North=pi/2)

    """
    input_alt = np.asarray(altitude)
    input_az = np.asarray(azimuth)
    if input_alt.size != input_az.size:
        raise ValueError("number of altitude and azimuth values must match.")

    zenith_angle = np.pi / 2 - input_alt
    new_azimuth = np.pi / 2 - input_az

    if new_azimuth.size > 1:
        wh_neg = np.where(new_azimuth < -1e-9)
        if wh_neg[0].size > 0:
            new_azimuth[wh_neg] = new_azimuth[wh_neg] + np.pi * 2
    elif new_azimuth.size == 1:
        if new_azimuth < -1e-9:
            new_azimuth = new_azimuth + np.pi * 2

    return zenith_angle, new_azimuth


def zenithangle_azimuth_to_altaz(zenith_angle, azimuth):
    """
    Convert from astropy altaz convention to UVBeam az/za convention.

    Parameters
    ----------
    zenith_angle: float, array_like of float
        Zenith angle in radians
    azimuth: float, array_like of float
        Azimuth in radians in uvbeam convention: North of East(East=0, North=pi/2)

    Returns
    -------
    altitude: array of float
        Altitude in radians
    azimuth: array of float
        In radians in astropy convention: East of North (N=0, E=pi/2)

    """
    input_za = np.array(zenith_angle)
    input_az = np.array(azimuth)
    if input_za.size != input_az.size:
        raise ValueError("number of zenith_angle and azimuth values must match.")

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


def check_file_exists_and_increment(filepath):
    """
    Check for a file and increment the name if it does to ensure a unique name.

    Given filepath (path + filename), check if it exists. If so, add a _1
    at the end, if that exists add a _2, and so on.

    """
    base_filepath, ext = os.path.splitext(filepath)
    bf_list = base_filepath.split("_")
    if bf_list[-1].isdigit():
        base_filepath = "_".join(bf_list[:-1])
    n = 0
    while os.path.exists(filepath):
        filepath = "{}_{}".format(base_filepath, n) + ext
        n += 1
    return filepath


def write_uvdata(
    uv_obj,
    param_dict,
    return_filename=False,
    dryrun=False,
    out_format=None,
    fix_autos=True,
    quiet=False,
):
    """
    Parse output file information from parameters and write out to a file.

    Parameters
    ----------
    uv_obj : UVData Object
        The object to be written out.
    param_dict : Dict
        parameter dictionary defining output path, filename, and whether or not to clobber.
    return_filename : Bool
        (Default false) Return the file path
    dryrun : Bool
        (Default false) Don't write to file.
    out_format : Str
        (Default uvh5) Write as uvfits/miriad/uvh5/ms
    fix_autos : bool
        If auto-correlations with imaginary values are found, fix those values so
        that they are real-only in data_array.
    quiet : bool
        If True, do not print anything to stdout.

    Returns
    -------
        File path, if return_filename is True

    """
    if "filing" in param_dict.keys():
        param_dict = param_dict["filing"]
    if "outdir" not in param_dict:
        param_dict["outdir"] = "."

    if "outfile_name" not in param_dict or param_dict["outfile_name"] == "":
        outfile_prefix = ""
        outfile_suffix = "results"
        if "outfile_prefix" in param_dict:
            outfile_prefix = param_dict["outfile_prefix"]
        if "outfile_suffix" in param_dict:
            outfile_suffix = param_dict["outfile_suffix"]
        outfile_name = "_".join([outfile_prefix, outfile_suffix])
        outfile_name = os.path.join(param_dict["outdir"], outfile_name)
    else:
        outfile_name = os.path.join(param_dict["outdir"], param_dict["outfile_name"])

    _, file_extension = os.path.splitext(outfile_name)

    if "output_format" in param_dict:
        out_format = param_dict["output_format"]
    elif file_extension in [".uvfits", ".uvh5", ".ms"]:
        out_format = file_extension[1:]
    elif out_format is None:
        # should be removed eventually. Maybe in v1.4? (it was added in 1.2.6)
        warnings.warn(
            "No out format specified for uvdata file. Defaulting to uvh5 (note "
            "this is a defaulting change, it used to default to uvfits)."
        )
        out_format = "uvh5"

    if not os.path.exists(param_dict["outdir"]):
        os.makedirs(param_dict["outdir"])

    if out_format in ["uvfits", "uvh5", "ms"]:
        if not outfile_name.endswith(f".{out_format}"):
            outfile_name = outfile_name + f".{out_format}"

    noclobber = ("clobber" not in param_dict) or not bool(param_dict["clobber"])
    if noclobber:
        outfile_name = check_file_exists_and_increment(outfile_name)

    if not quiet:
        print("Outfile path: ", outfile_name, flush=True)
    if not dryrun:
        if out_format == "uvfits":
            uv_obj.write_uvfits(outfile_name, force_phase=True, fix_autos=fix_autos)
        elif out_format == "miriad":
            uv_obj.write_miriad(
                outfile_name, clobber=not noclobber, fix_autos=fix_autos
            )
        elif out_format == "uvh5":
            uv_obj.write_uvh5(outfile_name, clobber=not noclobber, fix_autos=fix_autos)
        elif out_format == "ms":
            try:
                import casacore.tables  # noqa
                import casacore.tables.tableutil  # noqa
            except ImportError as error:  # pragma: no cover
                raise ImportError(
                    "casacore is not installed but is required for measurement set "
                    "functionality"
                ) from error
            uv_obj.write_ms(outfile_name, clobber=not noclobber, fix_autos=fix_autos)
        else:
            raise ValueError(
                "Invalid output format. Options are 'uvfits', 'uvh5', 'miriad' or 'ms'."
            )
    if return_filename:
        return outfile_name


def get_avail_memory():
    """
    Estimate the virtual memory available (in bytes).

    This gives the available memory on the current node to a running process.

    If this is not called from within a SLURM task, it will estimate
    using psutils methods.

    """
    slurm_key = "SLURM_MEM_PER_NODE"
    if slurm_key in os.environ:
        return float(os.environ[slurm_key]) * 1e6  # MB -> B

    return psutil.virtual_memory().available


def iter_array_split(part_index, N, M):
    """
    Return an iterator giving the indices of part of an array.

    part = np.array_split(np.arange(N), M)[part_index]

    This mimics the behavior of numpy.array_split without having to make
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
    Estimate the memory footprint of a :class:`pyradiosky.SkyModel`.

    By summing the sizes of the data types that go into SkyModel.

    This aims to anticipate the full memory required to handle a SkyModel
    class in simulation, accounting for its attributes as well as
    intermediate data generated.

    Parameters
    ----------
    Ncomponents : int
        Number of source components.
    Nfreqs : int
        Number of frequencies per source component.
        (size of :attr:`pyradiosky.SkyModel.freq_array`)

    Returns
    -------
    mem_est : float
        Estimate of memory usage in bytes

    """
    base_float = [1.5]  # A float
    base_bool = [True]
    base_str = ["source_name"]

    Ncomp_attrs = {
        "ra": base_float,
        "dec": base_float,
        "alt_az": 2 * base_float,
        "rise_lst": base_float,
        "set_lst": base_float,
        "pos_lmn": 3 * base_float,
        "name": base_str,
        "horizon_mask": base_bool,
    }
    Ncomp_Nfreq_attrs = {
        "stokes": 4 * base_float,
        "coherency_radec": 4 * base_float,
        "coherency_local": 4 * base_float,
    }

    mem_est = np.sum([sys.getsizeof(v) * Ncomponents for k, v in Ncomp_attrs.items()])
    mem_est += np.sum(
        [sys.getsizeof(v) * Ncomponents * Nfreqs for k, v in Ncomp_Nfreq_attrs.items()]
    )
    return mem_est


def _grouper_it(iterable, chunksize=1):
    """Chunk an iterator and return an iterator.

    Parameters
    ----------
    iterable : Iterable
        The iterable object to chunk
    chunksize : int
        size of chunks desired

    Returns
    -------
    iterable chunked into sizes
    """
    it = iter(iterable)
    while True:
        chunk = list(itertools.islice(it, chunksize))
        if not chunk:
            return
        yield chunk


def _chunked_iterator_product(iter1, iter2, chunksize1, chunksize2):
    """Iterate over the product of two chunked iterators.

    Parameters
    ----------
    iter1 : Iterable
        One iterator to chunk through.
    iter2 : Iterable
        The second iterator to chunk through.
    chunksize1 : int
        Chunk size for iter1
    chunksize2 : int
        Chunk size for iter2

    Returns
    -------
    An iterator over the chunked product of all combinations of iter1 and iter2

    """
    for i1 in _grouper_it(iter1, chunksize1):
        for i2 in _grouper_it(iter2, chunksize2):
            yield i1, i2
