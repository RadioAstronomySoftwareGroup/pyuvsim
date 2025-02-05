# Copyright (c) 2025 Radio Astronomy Software Group
# Licensed under the 2-clause BSD License
"""Command line scripts."""

import argparse
import os
import warnings

try:
    from astropy.config import get_cache_dir
    from astropy.utils.data import download_file, get_cached_urls, import_file_to_cache
    from pyradiosky.utils import download_gleam
except ImportError as e:
    raise ImportError(
        "astropy, pyradiosky, astroquery are required for full script functionality"
    ) from e

# TODO: should the bulk of this code be here?
# TODO: add verbosity as a lot of these should maybe not be warnings


def download_gleam_vot(url="file://gleam.vot"):
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

    warnings.warn("gleam.vot download can take a while")
    # downloads gleam.vot using astroquery
    cache_dir = get_cache_dir("pyuvsim")
    download_gleam(path=cache_dir, filename="gleam.vot")

    path_to_file = os.path.join(cache_dir, "gleam.vot")

    # converts downloaded gleam.vot to format recognizable for astropy with fake url key
    if os.path.exists(path_to_file):
        warnings.warn(
            "running import_file_to_cache to convert to astropy downloaded format"
        )
        import_file_to_cache(url, path_to_file, remove_original=True, pkgname="pyuvsim")
    else:
        warnings.warn(f"gleam.vot not downloaded to expected location: {path_to_file}")


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
        description="A command-line script to download the GLEAM vot table from Vizier."
    )

    # dictionary of keywords for downloadable files with url values (gleam we don't have a url)
    file_download_dict = {
        # mwa uvbeam
        "mwa": "http://ws.mwatelescope.org/static/mwa_full_embedded_element_pattern.h5",
        # gsm 2016 nside 128 healpix map
        "healpix": "https://repository.library.brown.edu/storage/bdr:eafzyycj/content/",
        "gleam": "",
    }

    parser.add_argument(
        "files",
        nargs="*",  # '*' consumes zero or more arguments, storing them in a list
        default=file_download_dict.keys(),  # If no arguments are provided, use the full list
        help="List of files to download via keyword. Downloads all available files if no arguments"
        f" given. Currently available file keywords: {' '.join(file_download_dict.keys())}."
        ' Sample use: "python3 download_data_files.py'
        f' {" ".join(list(file_download_dict.keys())[:2])}"',
    )

    args = parser.parse_args(argv)

    # loop through every file keyword in the list of filenames passed in from the command line
    # for each file keyword, check that it exists in the dictionary, then try to download
    for file in args.files:
        if file in file_download_dict:
            # call download_gleam_vot if gleam, otherwise download from url
            if file == "gleam":
                download_gleam_vot()
            else:
                download_file_using_astropy(file_download_dict[file])
        else:
            warnings.warn(
                f'file "{file}" not found! Check the available keys in the file_download_dict'
                "or the listed current options at the top of the file."
            )


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
        description="a list of reference simulations to download"
    )

    # TODO: Implement support for local file downloading (currently only allowed to download files
    #       to the astropy cache). This method currently only nicely functions as a function call
    #       that returns the url key of the astropy cachedd file so isn't very command line
    #       functional.

    parser.add_argument(
        "sims",
        nargs="*",  # '*' consumes zero or more arguments, storing them in a list
        default=ref_sims,  # If no arguments are provided, use FULL_LIST as the default
        help="List of reference simulations to download via keyword. Downloads all available"
        "files if"
        " no arguments given. Currently available file keywords:"
        f"{' '.join(ref_sims)}."
        f' Sample use: "python3 download_latest_reference_sims.py {" ".join(ref_sims)}"',
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

    return return_paths
