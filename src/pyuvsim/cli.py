# Copyright (c) 2025 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License
"""Command line scripts."""

import argparse
import time as pytime
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import psutil
from astropy import units as u
from astropy.coordinates import Angle, EarthLocation, SkyCoord
from astropy.time import Time
from astropy.units import Quantity
from pyradiosky import SkyModel
from pyuvdata import UVBeam, UVData

from pyuvsim import profiling, simsetup, telescope, utils, uvsim
from pyuvsim.simsetup import _parse_layout_csv as parse_csv


def run_pyuvsim(argv=None):
    """Run a pyuvsim simulation from a parameter file or a UVData, SkyModel and BeamList."""
    from pyuvsim import mpi

    parser = argparse.ArgumentParser(
        description="A command-line script to execute a pyuvsim simulation from "
        "either a single parameter file or all of: UVData readable file, "
        "SkyModel readable file, UVBeam readable file and output file."
    )
    parser.add_argument("--param", type=str, help="Parameter yaml file.", default=None)
    parser.add_argument("--uvdata", type=str, help="UVData readable file.")
    parser.add_argument("--uvbeam", type=str, help="UVBeam readable file.")
    parser.add_argument("--skymodel", type=str, help="SkyModel readable file.")
    parser.add_argument(
        "--outfile", type=str, help="File to save output UVData object to."
    )

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
        [
            args.uvdata is not None,
            args.uvbeam is not None,
            args.skymodel is not None,
            args.outfile is not None,
        ]
    )

    if args.param is None:
        if np.all(~uvdata_set) or not np.all(uvdata_set):
            raise ValueError(
                "Either pass a parameter file or all of: uvdata, skymodel, "
                "uvbeam files and outfile."
            )
    elif np.any(uvdata_set):
        raise ValueError(
            "A parameter file and was passed along with one or more of uvdata, "
            "skymodel, uvbeam and outfile. Either pass a parameter file or all of: "
            "uvdata, skymodel, uvbeam files and outfile."
        )

    if args.profile is not None:
        profiling.set_profiler(outfile_prefix=args.profile, dump_raw=args.raw_profile)

    t0 = pytime.time()

    block_nonroot_stdout = not args.keep_nonroot_stdout

    if args.param is not None:
        uvsim.run_uvsim(
            args.param, quiet=args.quiet, block_nonroot_stdout=block_nonroot_stdout
        )
    else:
        uvd = UVData.from_file(args.uvdata)
        skymodel = SkyModel.from_file(args.skymodel)
        uvb = UVBeam.from_file(args.uvbeam)
        beam_list = telescope.BeamList([uvb])

        uvd_out = uvsim.run_uvdata_uvsim(
            input_uv=uvd,
            beam_list=beam_list,
            beam_dict=None,
            catalog=skymodel,
            quiet=args.quiet,
            block_nonroot_stdout=block_nonroot_stdout,
        )
        pobj = Path(args.outfile)
        utils.write_uvdata(
            uvd_out,
            param_dict={"outdir": str(pobj.parent), "outfile_name": str(pobj.name)},
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
    parser.add_argument("-b", "--beamfile", type=str, default=None)
    parser.add_argument("-l", "--layout_csv_name", default=None)
    parser.add_argument("-t", "--telescope_config_name", default=None)

    args = parser.parse_args(argv)

    if args.beamfile is None:
        raise ValueError("beamfile must be passed.")

    uvd = UVData()
    uvd.read(args.file_in[0])

    simsetup.uvdata_to_telescope_config(
        uvd,
        args.beamfile,
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
    dpi = int(np.ceil(char_pix_height / emsize))
    # dpi = height of a letter in pixels/hight of letter in inches
    im_cmd = [
        "convert",
        "-background",
        "black",
        "-fill",
        "white",
        "-font",
        "C059-Bold",  # still looking for optimal font
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
        plt.savefig(f"{catname}.png", bbox_inches="tight")
        plt.clf()

    n_srcs = len(icrs_coord)
    sky = SkyModel(
        name=[f"text_{text}_{i}" for i in range(n_srcs)],
        skycoord=icrs_coord,
        stokes=Quantity(np.ones((4, 1, n_srcs)), "Jy"),
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
        "--lat",
        type=float,
        default=-30,
        help="Latitude in degrees on which to center the text at zenith.",
    )
    parser.add_argument(
        "--lon",
        type=float,
        default=116,
        help="Longitude in degrees on which to center the text at zenith.",
    )
    parser.add_argument("--plot", action="store_true")

    parser.add_argument("-v", "--verbose", action="count", default=0)

    args = parser.parse_args(argv)

    create_text_catalog(
        text=args.text,
        char_pix_height=args.char_pix_height,
        lat=args.lat,
        lon=args.lon,
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
    layout_path = Path(args.layout_file)
    plot_path = layout_path.stem + ".png"

    enu = antpos[["e", "n", "u"]]
    antnames = antpos["number"]

    plt.scatter(enu["e"], enu["n"])
    for i, _ in enumerate(enu):
        plt.annotate(antnames[i], (enu["e"][i], enu["n"][i]))

    plt.xlabel("East [m]")
    plt.ylabel("North [m]")
    plt.savefig(str(plot_path), bbox_inches="tight")
    plt.clf()
