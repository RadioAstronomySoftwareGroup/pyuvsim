# Copyright (c) 2025 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License
"""Command line scripts."""

import argparse
import os
import subprocess  # nosec
import time as pytime
import warnings
from argparse import RawTextHelpFormatter
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import psutil
from astropy import units as u
from astropy.config import get_cache_dir
from astropy.coordinates import Angle, EarthLocation, SkyCoord
from astropy.time import Time
from astropy.units import Quantity
from astropy.utils.data import (
    clear_download_cache,
    download_file,
    get_cached_urls,
    import_file_to_cache,
)
from pyradiosky import SkyModel
from pyradiosky.utils import download_gleam
from pyuvdata import UVBeam, UVData

from pyuvsim import profiling, simsetup, telescope, utils, uvsim
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
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
    parser.add_argument(
        "--uvdata", type=str, help="UVData readable file.", default=None
    )
    parser.add_argument(
        "--uvbeam", type=str, help="UVBeam readable file.", default=None
    )
    parser.add_argument(
        "--skymodel", type=str, help="SkyModel readable file.", default=None
    )
    parser.add_argument(
        "--outfile",
        type=str,
        help="File to save output UVData object to.",
        default=None,
    )

    parser.add_argument(
        "--profile", type=str, help="Time profiling output file name.", default=None
    )
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
        if not np.all(uvdata_set):
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
        if hasattr(profiling.prof, "meta_file"):
            with open(profiling.prof.meta_file, "a") as afile:
                afile.write(f"Runtime \t {rtime}\nMaxRSS \t {maxrss:.3f}\n")
                afile.write(f"Date/Time \t {str(datetime.now())}")


def run_param_pyuvsim(argv=None):
    """
    Run a pyuvsim simulation from a parameter file. Deprecated.

    Deprecated in favor of run_pyuvsim.
    """
    parser = argparse.ArgumentParser(
        description="A command-line script to execute a pyuvsim simulation from a parameter file."
    )
    parser.add_argument(
        "paramsfile", type=str, help="Parameter yaml file.", default=None
    )
    parser.add_argument(
        "--profile", type=str, help="Time profiling output file name.", default=None
    )
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

    warnings.warn(
        "This script is deprecated in favor of run_pyuvsim. It will be removed in version 1.6",
        DeprecationWarning,
    )

    args = parser.parse_args(argv)

    arglist = ["--param", args.paramsfile]
    if args.profile is not None:
        arglist.extend(["--profile", args.profile])
    if args.quiet:
        arglist.append("--quiet")
    if args.keep_nonroot_stdout:
        arglist.append("--keep_nonroot_stdout")
    if args.raw_profile:
        arglist.append("--raw_profile")

    run_pyuvsim(arglist)


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


def create_text_catalog(
    text: str,
    char_pix_height: int = 10,
    lat: float = -30,
    lon: float = 116,
    jd: float = 2460000.0,
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
    try:
        import matplotlib.image as mpimg
        import matplotlib.pyplot as plt
    except ImportError as ie:
        raise ImportError(
            "matplotlib must be installed to use create text catalogs."
        ) from ie

    try:
        subprocess.check_output(["convert", "--version"])  # nosec
    except (FileNotFoundError, subprocess.CalledProcessError) as err:
        raise RuntimeError(
            "ImageMagick must installed to create text catalogs"
        ) from err

    catname = "".join(text.split())
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
    del icrs_coord.obstime
    del icrs_coord.location
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
        default=2460000,
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


def im_to_catalog(argv=None):
    """
    Create a test pattern catalog that spells out something. Deprecated.

    Deprecated in favor of text_to_catalog.

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
        "-s",
        dest="spacing",
        type=float,
        default=1.0,
        help="Not used, retained for backwards compatibility.",
    )
    parser.add_argument(
        "-n",
        dest="char_pix_height",
        type=int,
        default=10,
        help="how many sources high should a letter be",
    )
    parser.add_argument("--jd", metavar="jd", type=float, default=2468000)
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
    args = parser.parse_args(argv)

    warnings.warn(
        "This script is deprecated in favor of text_to_catalog. It will be removed in version 1.6",
        DeprecationWarning,
    )

    arglist = [
        "-t",
        args.text,
        "-n",
        str(args.char_pix_height),
        "--jd",
        str(args.jd),
        "--thresh",
        str(args.thresh),
        "--lat",
        str(args.lat),
        "--lon",
        str(args.lon),
    ]
    if args.plot:
        arglist.append("--plot")
    arglist.append("-vv")

    text_to_catalog(arglist)


def plot_csv_antpos(argv=None):
    """Plot antenna positions from a layout csv."""
    try:
        import matplotlib.pyplot as plt
    except ImportError as ie:
        raise ImportError(
            "matplotlib must be installed to use plot antenna positions."
        ) from ie

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
    plt.gca().set_aspect("equal")
    plt.savefig(str(plot_path), bbox_inches="tight")
    plt.clf()


def download_gleam_vot(url="file://gleam.vot", row_limit=None):
    """
    Download gleam.vot to pyuvsim cache specified by astropy.

    Runs import_file_to_cache to format the installed file to be recognizable through astropy.
    Includes adding a fake url with which to specify the file. The fake download url acts as
    the key to load the file in astropy -- allowing filename specification in yaml files
    without any path details as astropy handles searching the cache with get_cached_urls
    and download_file.

    Parameters
    ----------
    url: str
        A fake url for the gleam.vot file.
    """
    # check if url already cached
    if url in get_cached_urls("pyuvsim"):
        warnings.warn(
            f"astropy cached url for {url} already exists! If you want to redownload,"
            " please clear the cached url."
        )
        return

    # around 1/4 of total rows print warning
    if row_limit is None or row_limit > 75000:
        print("gleam.vot download can take a while")

    # downloads gleam.vot using astroquery
    cache_dir = get_cache_dir("pyuvsim")
    download_gleam(path=cache_dir, filename="gleam.vot", row_limit=row_limit)

    path_to_file = os.path.join(cache_dir, "gleam.vot")

    # converts downloaded gleam.vot to format recognizable for astropy with fake url key
    if os.path.exists(path_to_file):
        print("running import_file_to_cache to convert to astropy downloaded format")
        path_to_file = import_file_to_cache(
            url, path_to_file, remove_original=True, pkgname="pyuvsim"
        )
    else:  # pragma: nocover
        warnings.warn(f"gleam.vot not downloaded to expected location: {path_to_file}")

    return path_to_file


def download_file_using_astropy(url):
    """
    Download file at url to pyuvsim cache specified by astropy using download_file.

    Parameters
    ----------
    url: str
        The url of the file to be downloaded.

    The download url acts as the key to load the file in astropy -- allowing filename
    specification in yaml files without any path details as astropy handles searching
    the cache with get_cached_urls and download_file.
    """
    # check if url already cached
    if url in get_cached_urls("pyuvsim"):
        warnings.warn(
            f"astropy cached url for {url} already exists! If you want to redownload,"
            " please clear the cached url."
        )

    # astropy file download method which defaults to cached file and hashes the file
    # should have a default timeout of 10 seconds
    filepath = download_file(url, cache=True, pkgname="pyuvsim")

    return filepath


def download_data_files(argv=None):
    """
    Large file downloading using the astropy cache to store beams and catalogs.

    Takes a command line argument of the file to download.

    Current options:
        - "gleam": the gleam catalog as a VOTable
        - "mwa": mwa uvbeam file
        - "healpix": gsm 2016 nside 128 healpix map saved as skyh5

    All files require astropy. Gleam additionally requires pyradiosky and astroquery.

    Downloads files to astropy default cache location for pyuvsim e.g. "get_cache_dir(pyuvsim)"

    The intention is to have all files downloaded under the astropy paradigm, and thus findable
    in cache throug the same paradigm with methods get_cached_urls and download_file. This allows
    us to list only the file url corresponding to the astropy cached file in either yaml instead
    of having to worry about relative or absolute path variables.
    """
    parser = argparse.ArgumentParser(
        description="A command-line script to download large data files "
        "for the reference simulations.",
        formatter_class=RawTextHelpFormatter,
    )

    # dictionary of keywords for downloadable files with url values (gleam we don't have a url)
    file_download_dict = {
        # mwa uvbeam
        "mwa": "http://ws.mwatelescope.org/static/mwa_full_embedded_element_pattern.h5",
        # gsm 2016 nside 128 healpix map
        "healpix": "https://repository.library.brown.edu/storage/bdr:eafzyycj/content/",
        "gleam": "",
    }

    # create some strings for formatted help output
    formatted_key_list = "\n\t-" + "\n\t-".join(file_download_dict.keys()) + "\n\n"
    keyword_options = " ".join(list(file_download_dict.keys())[:2])

    parser.add_argument(
        "files",
        nargs="*",  # '*' consumes zero or more arguments, storing them in a list
        default=file_download_dict.keys(),  # If no arguments are provided, use the full list
        help="List of files to download via keyword. Downloads all available files as default if "
        "no file arguments are given. Currently available file keywords:"
        f'{formatted_key_list}Sample use: "download_data_files {keyword_options}"',
    )

    parser.add_argument(
        "--row_limit",
        type=int,
        default=None,
        help="Optionally specify a maximum number of rows to download for GLEAM. May not always "
        "download the same set of rows so this option should not be used for anything requiring "
        "repeatablility!",
    )

    parser.add_argument(
        "--clear",
        action="store_true",
        default=False,
        help="If clear is passed, the pyuvsim cache is first wiped by calling "
        'astropy\'s clear_download_cache(pkgname="pyuvsim"). All other input '
        "is treated subsequently.",
    )

    args = parser.parse_args(argv)

    if args.clear:
        print("Clearing entirety of pyuvsim cache.")
        clear_download_cache(pkgname="pyuvsim")

    return_paths = []

    # loop through every file keyword in the list of filenames passed in from the command line
    # for each file keyword, check that it exists in the dictionary, then try to download
    for file in args.files:
        if file in file_download_dict:
            # call download_gleam_vot if gleam, otherwise download from url
            if file == "gleam":
                filepath = download_gleam_vot(row_limit=args.row_limit)
            else:
                filepath = download_file_using_astropy(file_download_dict[file])

            return_paths.append(filepath)
        else:
            warnings.warn(
                f'file "{file}" not found! Check the available keys in the file_download_dict'
                "or the listed current options at the top of the file."
            )

    # NOTE: this causes a print (maybe to stderr?) at the end which I don't know how to avoid
    return return_paths


def robust_response(url, n_retry=10):
    """
    Attempt to GET the url n_retry times if return code in range.

    returns the GET request
    Stackoverflow / Docs link for approach / Retry:
    https://stackoverflow.com/questions/15431044/can-i-set-max-retries-for-requests-request
    https://urllib3.readthedocs.io/en/latest/reference/urllib3.util.html#module-urllib3.util.retry
    """
    import requests
    from requests.adapters import HTTPAdapter, Retry

    s = requests.Session()
    retries = Retry(
        total=n_retry,
        backoff_factor=10,
        raise_on_status=True,
        status_forcelist=range(500, 600),
    )
    s.mount("http://", HTTPAdapter(max_retries=retries))

    return s.get(url)


def get_latest_api_response_pid(sim_name):
    """
    Get the pid of the most recent upload matching the sim_name on the Brown Digital Repository.

    Sends an api call to the "pyuvsim historical reference simulations" collection,
    then identifies and returns the latest uploaded object in the response.

    Link to BDR API DOCS:
    https://repository.library.brown.edu/studio/api-docs/

    Parameters
    ----------
    sim_name: str
        The simulation name to use as a search pattern

    Returns
    -------
    str
        the pid as a string
    """
    # api url
    api_url = (
        "https://repository.library.brown.edu/api/collections/bdr:wte2qah8/?q="
        f"primary_title:{sim_name}&wt=json&sort=object_last_modified_dsi desc&rows=1"
    )

    print(
        f"\nquerying BDR collection for latest item matching {sim_name} via api: {api_url}"
    )
    response = robust_response(api_url)

    # attempt to parse json assuming we have the result
    try:
        # get pid from first response item as we query by time desc and only get 1 result
        pid = response.json()["items"]["docs"][0]["pid"]
    except Exception as e:
        warnings.warn(f"Failed to parse BDR response for file PID with error: {e}")

    return pid


def download_ref_sims(argv=None):
    """Download latest reference simulations from the Brown Digital Repository using astropy."""
    ref_sims = [
        "1.1_baseline_number",
        "1.2_time_axis",
        "1.3_frequency_axis",
        "1.4_source_axis",
        "1.5_uvbeam",
        "1.6_healpix",
        "1.7_multi_beam",
        "1.8_lunar",
    ]

    # create parser that grabs all input command line arguments as file keywords
    parser = argparse.ArgumentParser(
        description="A command line script to download the latest reference simulations.",
        formatter_class=RawTextHelpFormatter,
    )

    # create some strings for formatted help output
    formatted_key_list = "\n\t-" + "\n\t-".join(ref_sims) + "\n\n"
    keyword_options = " ".join(ref_sims[:2])

    parser.add_argument(
        "sims",
        nargs="*",  # '*' consumes zero or more arguments, storing them in a list
        default=ref_sims,  # If no arguments are provided, use FULL_LIST as the default
        help="List of reference simulations to download via keyword. Downloads all available "
        "files if no arguments given. Currently available file keywords:"
        f'{formatted_key_list}Sample use: "download_ref_sims {keyword_options}"',
    )

    args = parser.parse_args(argv)

    return_paths = []

    for sim in args.sims:
        if sim in ref_sims:
            # for each sim: get pid, construct download url, download file and return path
            pid = get_latest_api_response_pid(sim)
            download_url = (
                f"https://repository.library.brown.edu/storage/{pid}/content/"
            )
            filepath = download_file_using_astropy(download_url)
            return_paths.append(filepath)
        else:
            warnings.warn(
                f'sim "{sim}" not found! Check the available keys in the ref_sims list'
                " or the listed current options at the top of the file."
            )

    # NOTE: this causes a print (maybe to stderr?) at the end which I don't know how to avoid
    return return_paths


def run_ref_sim(argv=None):
    """Run the reference simulations using a wrapper of run_pyuvsim."""
    ref_sims = {
        "1": "1.1_baseline_number",
        "2": "1.2_time_axis",
        "3": "1.3_frequency_axis",
        "4": "1.4_source_axis",
        "5": "1.5_uvbeam",
        "6": "1.6_healpix",
        "7": "1.7_multi_beam",
        "8": "1.8_lunar",
    }

    # create parser that grabs all input command line arguments as file keywords
    parser = argparse.ArgumentParser(
        description="A command line script to run the reference simulations.",
        formatter_class=RawTextHelpFormatter,
    )

    parser.add_argument(
        "version",
        type=str,
        help="Specify which first generation reference simulation you wish to run. Writes the "
        "output to `results_data` in the current working directory. Can be done "
        "with numbers (1-8) or by explicit name:"
        "\n\t1: 1.1_baseline_number"
        "\n\t2: 1.2_time_axis"
        "\n\t3: 1.3_freq_axis"
        "\n\t4: 1.4_source_axis"
        "\n\t5: 1.5_uvbeam"
        "\n\t6: 1.6_healpix"
        "\n\t7: 1.7_multi_beam"
        "\n\t8: 1.8_lunar\n\n"
        'Sample use: "mpiexec run_ref_sim 1" or "mpiexec run_ref_sim 1.1_baseline_number"',
    )

    args = parser.parse_args(argv)

    ver = args.version

    # if passed as a number (1-8), convert to the actual name
    if ver in ref_sims:
        ver = ref_sims[ver]

    if ver not in ref_sims.values():
        raise ValueError("No match found for Input reference simulation.")

    # yaml filepath construction
    yaml_filename = "obsparam_ref_" + ver + ".yaml"
    yaml_filepath = os.path.join(SIM_DATA_PATH, "test_config", yaml_filename)

    # use run_pyuvsim to run the simulation
    run_pyuvsim(["--param", yaml_filepath])
