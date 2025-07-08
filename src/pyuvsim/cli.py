# Copyright (c) 2025 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License
"""Command line scripts."""

import argparse
import os
import resource
import subprocess
import sys
import time as pytime
from datetime import datetime, timedelta

import numpy as np
import psutil
import yaml
from astropy import units as u
from astropy.coordinates import Angle, EarthLocation, SkyCoord
from astropy.time import Time
from astropy.units import Quantity
from numpy.lib.recfunctions import append_fields
from pyradiosky import SkyModel
from pyuvdata import UVBeam, UVData

from pyuvsim import mpi, profiling, simsetup, telescope, uvsim
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
from pyuvsim.simsetup import _parse_layout_csv as parse_csv


def run_pyuvsim(argv=None):
    """Run a pyuvsim simulation from a parameter file or a UVData, SkyModel and BeamList."""
    parser = argparse.ArgumentParser(
        description="A command-line script to execute a pyuvsim simulation from "
        "either a single parameter file or all of: UVData readable file, "
        "SkyModel readable file, and UVBeam readable file."
    )
    parser.add_argument("--param", type=str, help="Parameter yaml file.", default=None)
    parser.add_argument("--uvdata", type=str, help="UVData readable file.")
    parser.add_argument("--uvbeam", type=str, help="UVBeam readable file.")
    parser.add_argument("--skymodel", type=str, help="SkyModel readable file.")

    parser.add_argument("--profile", type=str, help="Time profiling output file name.")
    parser.add_argument(
        "--quiet", action="store_true", help="Suppress stdout printing."
    )
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

    args = parser.parse_args(argv)

    uvdata_set = np.array(
        [args.uvdata is not None, args.uvbeam is not None, args.skymodel is not None]
    )

    if args.param is None:
        if np.all(~uvdata_set) or not np.all(uvdata_set):
            raise ValueError(
                "Either pass a parameter file or all of: uvdata, skymodel, uvbeam files."
            )
    elif np.any(uvdata_set):
        raise ValueError(
            "A parameter file and was passed along with one or more of uvdata, "
            "skymodel, uvbeam. Either pass a parameter file or all of: uvdata, "
            "skymodel, uvbeam files."
        )

    if args.profile is not None:
        profiling.set_profiler(outfile_prefix=args.profile, dump_raw=args.raw_profile)

    t0 = pytime.time()

    block_nonroot_stdout = not args.keep_nonroot_stdout

    if args.param is not None:
        if not os.path.isdir(os.path.dirname(args.param)):
            args.param = os.path.join(".", args.param)

        uvsim.run_uvsim(
            args.param, quiet=args.quiet, block_nonroot_stdout=block_nonroot_stdout
        )
    else:
        uvd = UVData.from_file(args.uvdata)
        skymodel = SkyModel.from_file(args.skymodel)
        uvb = UVBeam.from_file(args.uvbeam)
        beam_list = telescope.BeamList([uvb])

        uvsim.run_uvdata_uvsim(
            input_uv=uvd,
            beam_list=beam_list,
            beam_dict=None,
            catalog=skymodel,
            quiet=args.quiet,
            block_nonroot_stdout=block_nonroot_stdout,
        )

    if args.profile:
        dt = pytime.time() - t0
        maxrss = mpi.get_max_node_rss()
        rtime = str(timedelta(seconds=dt))
        if isinstance(maxrss, float):
            print(f"\tRuntime: {rtime} \n\tMaxRSS: {maxrss:.3f} GiB")
        if hasattr(profiling.prof, "meta_file"):
            with open(profiling.prof.meta_file, "a") as afile:
                afile.write(f"Runtime \t {rtime}\nMaxRSS \t {maxrss:.3f}\n")
                afile.write(f"Date/Time \t {str(datetime.now())}")


def uvdata_to_config(argv=None):
    """Generate basic simulation parameters from a file readable by :class:`pyuvdata.UVData`."""
    # Take a uvdata readable file, and save sim parameters as a yaml file.

    parser = argparse.ArgumentParser(
        description=(
            "Generate basic simulation parameters from a data file readable by "
            ":class:`pyuvdata.UVData`."
        )
    )

    parser.add_argument("file_in", metavar="<FILE>", type=str, nargs="+")
    parser.add_argument("-p", "--param_filename", default=None)
    parser.add_argument("-t", "--telescope_config_path", default="")
    parser.add_argument("-l", "--layout_csv_path", default="")
    parser.add_argument("--outpath", default=".")

    args = parser.parse_args(argv)

    uvd = UVData()
    uvd.read(args.file_in[0])

    simsetup.uvdata_to_config_file(
        uvd,
        param_filename=args.param_filename,
        telescope_config_name=args.telescope_config_path,
        layout_csv_name=args.layout_csv_path,
        path_out=args.outpath,
    )


def uvdata_to_telescope_config(argv=None):
    """Extract antenna position info from a data file readable by :class:`pyuvdata.UVData`."""
    # Read in a uvdata readable file and automatically generate yaml files.
    # Will assume the same beam_id for all antennas for now.

    parser = argparse.ArgumentParser(
        description=(
            "Extracts antenna position info from a data file readable by "
            ":class:`pyuvdata.UVData`."
        )
    )

    parser.add_argument("file_in", metavar="<FILE>", type=str, nargs="+")
    parser.add_argument("-l", "--layout_csv_name", default=None)
    parser.add_argument("-t", "--telescope_config_name", default=None)
    parser.add_argument(
        "-b",
        "--beam_filepath",
        type=str,
        default=os.path.join(SIM_DATA_PATH, "HERA_NicCST.uvbeam"),
    )

    args = parser.parse_args(argv)

    uvd = UVData()
    uvd.read(args.file_in[0])

    simsetup.uvdata_to_telescope_config(
        uvd,
        args.beam_filepath,
        layout_csv_name=args.layout_csv_name,
        telescope_config_name=args.telescope_config_name,
        return_names=False,
        path_out=".",
    )


def removeallspaces(s):
    """Remove all spaces."""
    return "".join([c for c in s if c != " "])


def create_text_catalog(
    text: str,
    char_pix_height: int = 10,
    lat: float = -30,
    lon: float = 116,
    jd: float = 2468000.0,
    plot: bool = False,
    thresh: int = 140,
    verbose: int = 0,
):
    """
    Create a test pattern catalog that spells out something.

    The test pattern is centered at zenith for a particular location and time.
    The test pattern is generated by creating a bitmap image of text and then
    turning each pixel with text in it into a source.

    Parameters
    ----------
    text :  str
        Text to make into sources.
    char_pix_height : int
        How many sources high should a letter be.
    lat : float
        Latitude in degrees on which to center the text at zenith.
    lon : float
        Longitude in degrees on which to center the text at zenith.
    jd : float
        Julian date on which to center the text at zenith.
    plot : bool
        Option to plot the sources.
    thresh : float
        Threshold for how bright a pixel needs to be to be convert to a source.
        A number between 1 and 255.
    verbose : int
        How verbose to be, default is 0, max is 2.
    """
    import matplotlib.image as mpimg
    import matplotlib.pyplot as plt

    catname = removeallspaces(text)
    imgfname = catname + ".bmp"

    # first lets construct our image
    fontsize = 10  # this ends up being arbitrary here
    emsize = 1 / 72.0 * fontsize  # height of a letter in inches
    dpi = np.int(np.ceil(char_pix_height / emsize))
    # dpi = height of a letter in pixels/hight of letter in inches
    im_cmd = [
        "convert",
        "-background",
        "black",
        "-fill",
        "white",
        "-font",
        "Keyboard",  # still looking for optimal font
        "-pointsize",
        "10",
        "-units",
        "PixelsPerInch",
        "-density",
        str(dpi),
        "-trim",
        "+repage",
        "-antialias",
        "label:" + text.upper(),
        f"{imgfname}",
    ]
    if verbose > 1:
        print(im_cmd)
    if verbose > 0:
        print("generating image file", imgfname)
    p = psutil.Popen(im_cmd)
    if verbose > 1:
        print(p.communicate())
    p.wait(timeout=2)

    if verbose > 0:
        print("processing image into a catalog")
    im = mpimg.imread(imgfname)
    arr = im[:, :, 0]  # image is rgb but actually white, grab the r channel
    arr = arr[::-1, :]  # flip updown
    if verbose > 1:
        print(arr)
    # get the coordinates for each pixel
    nx, ny = arr.shape
    y, x = np.where(arr > thresh)
    dx = x - int(np.floor(nx / 2.0))
    dy = y - int(np.floor(ny / 2.0))

    pixsize_deg = 1.0

    zas = np.sqrt(dx**2 + dy**2) * pixsize_deg
    aztmp = np.arctan2(dy, dx) * 180.0 / np.pi
    azs = np.arctan2(dy, dx) * 180.0 / np.pi - 90.0

    alts = 90.0 - zas

    # convert from alt-az to ra/dec overhead
    time = Time(jd, scale="utc", format="jd")
    location = EarthLocation(lat=lat * u.degree, lon=lon * u.degree)

    source_coord = SkyCoord(
        alt=Angle(alts, unit=u.deg),
        az=Angle(azs, unit=u.deg),
        obstime=time,
        frame="altaz",
        location=location,
    )

    icrs_coord = source_coord.transform_to("icrs")
    if verbose > 0:
        print("generated", len(azs), "sources")
    if plot:
        zas_rad = np.radians(zas)
        aztmp_rad = np.radians(aztmp)
        plt.scatter(
            zas_rad * np.cos(aztmp_rad), zas_rad * np.sin(aztmp_rad), label="Original"
        )
        plt.legend()
        plt.show()

    n_srcs = len(icrs_coord)
    sky = SkyModel(
        name=[f"TEST{i}" for i in range(n_srcs)],
        skycoord=icrs_coord,
        stokes=Quantity(np.ones((4, n_srcs)), "Jy"),
        spectral_type="flat",
        component_type="point",
    )

    # save catalog
    catfile = catname + ".skyh5"
    sky.write_skyh5(catfile)
    print(f"saved catalog file to {catfile}")


def text_to_catalog(argv=None):
    """
    Create a test pattern catalog that spells out something.

    The test pattern is centered at zenith for a particular location and time.
    The test pattern is generated by creating a bitmap image of text and then
    turning each pixel with text in it into a source.
    """
    parser = argparse.ArgumentParser(
        description=("make a test pattern catalog that spells out something")
    )

    parser.add_argument(
        "-t",
        "--text",
        dest="text",
        default="EH",
        type=str,
        help="Text to make into sources",
    )
    parser.add_argument(
        "-n",
        dest="char_pix_height",
        type=int,
        default=10,
        help="how many sources high should a letter be",
    )
    parser.add_argument(
        "--jd",
        metavar="jd",
        type=float,
        default=2468000,
        help="Julian date on which to center the text at zenith.",
    )
    parser.add_argument(
        "--thresh",
        metavar="thresh",
        type=float,
        default=140,
        help="""Threshold for converting font pixels to sources.
                        number between 1 and 255.""",
    )
    parser.add_argument(
        "--array_location",
        type=str,
        default="-30,116",
        help="LAT,LON in degrees on which to center the text at zenith.",
    )
    parser.add_argument("--plot", action="store_true")

    parser.add_argument("-v", "--verbose", action="count", default=0)

    args = parser.parse_args(argv)

    lat, lon = map(float, args.array_location.split(","))

    create_text_catalog(
        text=args.text,
        char_pix_height=args.char_pix_height,
        lat=lat,
        lon=lon,
        jd=args.jd,
        plot=args.plot,
        thresh=args.thresh,
        verbose=args.verbose,
    )


def plot_csv_antpos(argv=None):
    """Plot antenna positions from a layout csv."""
    import matplotlib.pyplot as plt

    parser = argparse.ArgumentParser(
        description=("Plot antenna positions from a layout csv.")
    )
    parser.add_argument("layout_file", help="File to plot antenna positions from.")

    args = parser.parse_args(argv)

    antpos = parse_csv(args.layout_file)

    enu = antpos[["e", "n", "u"]]
    antnames = antpos["number"]

    plt.scatter(enu["e"], enu["n"])
    for i, _ in enumerate(enu):
        plt.annotate(antnames[i], (enu["e"][i], enu["n"][i]))

    plt.xlabel("East [m]")
    plt.ylabel("North [m]")
    plt.show()


def run_profile_pyuvsim(argv=None):
    """Run a pyuvsim simulation for profiling purposes."""
    parser = argparse.ArgumentParser(
        description=(
            "A command-line script to execute a pyuvsim simulation for profiling purposes."
        )
    )

    parser.add_argument("--Nsrcs", dest="Nsrcs", type=int, default=1)
    parser.add_argument("--Ntimes", dest="Ntimes", type=int, default=1)
    parser.add_argument("--Nfreqs", dest="Nfreqs", type=int, default=1)
    parser.add_argument("--Nbls", dest="Nbls", type=int, default=1)
    parser.add_argument("--beam", dest="beam", type=str, default="uniform")
    parser.add_argument(
        "--prof_out", dest="prof_out", type=str, default="time_profile.out"
    )
    parser.add_argument(
        "--mem_out", dest="mem_out", type=str, default="memory_usage.out"
    )
    parser.add_argument(
        "--time_out", dest="time_out", type=str, default="time_usage.out"
    )

    args = parser.parse_args(argv)

    paramsfile = os.path.join(SIM_DATA_PATH, "profiling_params.yaml")

    with open(paramsfile) as pfile:
        params = yaml.safe_load(pfile)

    params["config_path"] = os.path.dirname(paramsfile)

    min_alt = 70  # Degrees

    mpi.start_mpi()
    rank = mpi.get_rank()

    if rank == 0:
        t0 = pytime.time()

    beam_list = None
    beam_dict = None
    input_uv = UVData()
    mock_keywords = None
    catalog = None
    profiling.set_profiler(outfile_prefix=args.prof_out)

    if rank == 0:
        print(
            f"{args.Nfreqs} freqs, {args.Ntimes} times, {args.Nbls} bls, "
            "{args.Nsrcs} srcs, {args.beam} beam"
        )
        params["freq"]["Nfreqs"] = args.Nfreqs
        params["time"]["Ntimes"] = args.Ntimes
        params["sources"] = {"catalog": "mock"}

        input_uv, beam_list, beam_dict = simsetup.initialize_uvdata_from_params(
            params, return_beams=True
        )

        if input_uv.Nbls < args.Nbls:
            raise ValueError(
                f"Cannot profile for more than {input_uv.Nbls} baselines, requeted {args.Nbls}"
            )

        # Baseline selection:
        input_uv.baseline_array = np.repeat(
            input_uv.baseline_array[: args.Nbls], args.Ntimes
        )
        input_uv.ant_1_array, input_uv.ant_2_array = input_uv.baseline_to_antnums(
            input_uv.baseline_array
        )
        ants_new = np.unique(
            input_uv.ant_1_array.tolist() + input_uv.ant_2_array.tolist()
        )
        input_uv.antenna_numbers = ants_new
        input_uv.antenna_names = ants_new.astype(str)
        Nants = ants_new.size
        # For now, all use the same beam model
        beam_dict = dict(
            zip(input_uv.antenna_names, np.zeros(Nants, dtype=int), strict=False)
        )
        input_uv.antenna_positions = input_uv.antenna_positions[:Nants, :]
        input_uv.Nants_data = Nants
        input_uv.Nants_telescope = Nants

        # Time selection:
        inds = np.array(
            [np.arange(args.Nbls) + i * input_uv.Nbls for i in range(args.Ntimes)]
        ).flatten()
        input_uv.time_array = input_uv.time_array[inds]
        input_uv.lst_array = input_uv.lst_array[inds]
        input_uv.Nbls = args.Nbls
        input_uv.Nblts = args.Nbls * args.Ntimes

        # Beam selection:
        # Default is uniform
        if args.beam == "hera":
            # TODO: fix this. How?
            beamfile = "/users/alanman/data/alanman/NickFagnoniBeams/HERA_NicCST_fullfreq.uvbeam"
            beam_list = [beamfile]

        mock_keywords = {
            "mock_arrangement": "random",
            "Nsrcs": args.Nsrcs,
            "min_alt": min_alt,
            "time": input_uv.time_array[0],
        }
        print(f"Beam: {beam_list[0]}")
        params["sources"].update(**mock_keywords)

        # Catalog setup
        catalog = simsetup.initialize_catalog_from_params(params, return_catname=False)

    comm = mpi.world_comm
    input_uv = comm.bcast(input_uv, root=0)
    beam_list = comm.bcast(beam_list, root=0)
    beam_dict = comm.bcast(beam_dict, root=0)
    catalog = mpi.shared_mem_bcast(catalog, root=0)
    if rank == 0:
        print("Starting simulation.")
        sys.stdout.flush()

    uvsim.run_uvdata_uvsim(input_uv, beam_list, beam_dict=beam_dict, catalog=catalog)

    memory_usage_gb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
    comm.Barrier()

    memory_usage_gb = comm.gather(memory_usage_gb, root=0)

    if rank == 0:
        elapsed_time_sec = pytime.time() - t0
        print("Elapsed: " + str(elapsed_time_sec))
        with open(args.time_out, "w") as timefile:
            timefile.write(str(elapsed_time_sec))

        memory_usage_gb = np.min(memory_usage_gb)
        print("Mem_usage: " + str(memory_usage_gb))
        with open(args.mem_out, "w") as memfile:
            memfile.write(str(memory_usage_gb))


def summarize_profiling(argv=None):
    """Summarize profiling results in a table."""
    parser = argparse.ArgumentParser(
        description=("Summarize profiling results in a table.")
    )
    parser.add_argument("filename", help="File to get profiling results from.")
    args = parser.parse_args(argv)

    slurmids = []

    fname = args.filename

    with open(fname, "a+") as fhandle:
        header = fhandle.readline()

    header = [h.strip().upper() for h in header.split(",")]
    dt = np.format_parser(
        ["i4", "i4", "i4", "i4", "U8", "i4"],
        ["Nsrcs", "Ntimes", "Nfreqs", "Nbls", "beam", "slurm_id"],
        header,
    )

    # TODO: this is reading in data from a file. How does the file get generated?
    filedat = np.genfromtxt(
        fname, autostrip=True, skip_header=1, delimiter=",", dtype=dt.dtype
    )

    slurmids = filedat["slurm_id"].astype(str)
    slurmids = [sid + "_0.0" for sid in slurmids]

    p = subprocess.Popen(
        'sacct --jobs="'
        + ",".join(slurmids)
        + '" --format="JobID, Start, Elapsed, MaxRSS, NNodes, NTasks, NCPUS"',
        shell=True,
        stdout=subprocess.PIPE,
    )

    table = np.genfromtxt(
        p.stdout, dtype=None, names=True, comments="--", encoding=None
    )
    table["MaxRSS"] = (float(x[:-1]) * 1e3 / 1e9 for x in table["MaxRSS"])

    dt = table.dtype.descr
    ind = dt.index(("MaxRSS", table.dtype["MaxRSS"]))
    dt[ind] = ("MaxRSS (GB)", "f")
    ind = dt.index(("NTasks", table.dtype["NTasks"]))
    dt[ind] = ("NProcs", "i")  # So there's no confusion with slurm tasks vs. UVTasks.
    dt = np.dtype(dt)
    table = table.astype(dt)

    njobs = len(slurmids)
    Nsrcs, Nchans, Ntimes, Nbls = np.zeros((4, njobs)).astype(int)
    beam_type = []
    # They got resorted in the sacct task.
    for ent in filedat:
        sid = str(ent["slurm_id"]) + "_0.0"
        ind = np.where(table["JobID"] == sid)
        Nsrcs[ind] = ent["Nsrcs"]
        Nchans[ind] = ent["Nfreqs"]
        Ntimes[ind] = ent["Ntimes"]
        Nbls[ind] = ent["Nbls"]
        beam_type.append(ent["beam"])

    def hms2sec(hms):
        """Convert hour, minute, seconds to seconds."""
        h, m, s = map(float, hms.split(":"))
        return h * 60.0**2 + m * 60.0 + s

    runtime_sec = np.array(map(hms2sec, table["Elapsed"]))
    cores_per_node = table["NProcs"] / table["NNodes"]
    ntasks = Nsrcs * Ntimes * Nbls * Nchans

    table = append_fields(
        table,
        [
            "CoresPerNode",
            "Nbls",
            "Ntimes",
            "Nchan",
            "Nsrc",
            "Beam",
            "Ntasks",
            "Runtime_Seconds",
        ],
        [cores_per_node, Nbls, Ntimes, Nchans, Nsrcs, beam_type, ntasks, runtime_sec],
        usemask=False,
    )

    # rec2csv(table, "profiling_results_table.csv")
    np.savetxt("profiling_results_table.csv", table, delimiter=",")


def profiling_plots(argv=None):
    """Make plots of profiling results under different constraints."""
    import matplotlib.pyplot as plt

    parser = argparse.ArgumentParser(
        description=("Summarize profiling results in a table.")
    )
    parser.add_argument("filename", help="File to get profiling results from.")
    args = parser.parse_args(argv)

    fname = args.filename

    dat = np.genfromtxt(fname, names=True, max_rows=2, delimiter=",")

    fields = [
        "JobID",
        "Start",
        "MaxRSS_GB",
        "NNodes",
        "NProcs",
        "Nbls",
        "Ntimes",
        "Nchan",
        "Nsrc",
        "Beam",
        "Ntasks",
        "Runtime_Second",
    ]
    titles = [f.lower() for f in fields]
    fmts = ["U10", "U10", "f8", "i4", "i4", "i4", "i4", "i4", "i4", "U10", "i4", "f8"]
    dt = np.format_parser(fmts, dat.dtype.names, titles)

    dat = np.genfromtxt(
        sys.argv[1], autostrip=True, dtype=dt.dtype, delimiter=",", skip_header=1
    )

    Ncpus = np.unique(dat["NProcs"])
    beams = np.unique(dat["Beam"])
    # Ntasks, inv, counts = np.unique(dat["Ntasks"], return_inverse=True, return_counts=True)
    NNodes = np.unique(dat["NNodes"])
    markers = [".", "o", "+", "<", ">", "*", "v", "^", "h", "p"]
    cmap = plt.get_cmap("tab10")
    colors = [cmap(i) for i in range(len(Ncpus))]

    # params = ["Nbls", "Ntimes", "Nchan", "Nsrc"]
    _, axes = plt.subplots(nrows=1, ncols=len(beams))
    # handles, labels = [], []
    ymin, ymax = np.min(dat["Runtime_Seconds"]), np.max(dat["Runtime_Seconds"])
    ymax *= 1.20
    for nni in range(len(NNodes)):
        cond_n = dat["NNodes"] == NNodes[nni]
        for bi in range(len(beams)):
            cond_b = dat["Beam"] == beams[bi]
            for nc in range(len(Ncpus)):
                inds = np.where(cond_n & cond_b & (dat["NProcs"] == Ncpus[nc]))
                if len(inds[0]) == 0:
                    continue
                axes[bi].scatter(
                    dat["Ntasks"][inds],
                    dat["Runtime_Seconds"][inds],
                    label=f"{Ncpus[nc]:d} cores, {NNodes[nni]:d} Nodes",
                    color=colors[nc],
                    marker=markers[nni],
                )

            axes[bi].set_title(f"{beams[bi]} beam")
            axes[bi].set_ylabel("Runtime (seconds)")
            axes[bi].set_xlabel("Ntasks")
            axes[bi].set_ylim([ymin, ymax])

    def plot_handles(m, c):
        """Get some plot handles."""
        return plt.plot([], [], marker=m, color=c, ls="none")[0]

    colhandles = [plot_handles("s", colors[i]) for i in range(len(Ncpus))]
    markhandles = [plot_handles(markers[i], "k") for i in range(len(NNodes))]

    collabels = map("{:d} cpus".format, Ncpus.tolist())
    marklabels = map("{:d} nodes".format, NNodes.tolist())
    plt.figlegend(labels=collabels, handles=colhandles, loc="center right")
    leg2 = plt.figlegend(labels=marklabels, handles=markhandles, loc="lower right")
    plt.gca().add_artist(leg2)

    _, axes = plt.subplots(nrows=1, ncols=len(beams))
    ymin, ymax = np.min(dat["MaxRSS_GB"]), np.max(dat["MaxRSS_GB"])
    ymax *= 1.20
    for nni in range(len(NNodes)):
        cond_n = dat["NNodes"] == NNodes[nni]
        for bi in range(len(beams)):
            cond_b = dat["Beam"] == beams[bi]
            for nc in range(len(Ncpus)):
                inds = np.where(cond_n & cond_b & (dat["NProcs"] == Ncpus[nc]))
                if len(inds[0]) == 0:
                    continue
                axes[bi].scatter(
                    dat["Ntasks"][inds],
                    dat["MaxRSS_GB"][inds],
                    label=f"{Ncpus[nc]:d} cores, {NNodes[nni]:d} Nodes",
                    color=colors[nc],
                    marker=markers[nni],
                )

            axes[bi].set_title(f"{beams[bi]} beam")
            axes[bi].set_ylabel("MaxRSS (GB)")
            axes[bi].set_xlabel("Ntasks")
            axes[bi].set_ylim([ymin, ymax])
    plt.figlegend(labels=collabels, handles=colhandles, loc="center right")
    leg2 = plt.figlegend(labels=marklabels, handles=markhandles, loc="lower right")
    plt.gca().add_artist(leg2)

    plt.show()
