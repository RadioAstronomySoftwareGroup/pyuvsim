#!/bin/env python
# Copyright (c) 2022 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License
"""
An approach to large file downloading using the astropy cache to store beams and catalogs.

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

import argparse
import os

try:
    from astropy.config import get_cache_dir
    from astropy.utils.data import download_file, get_cached_urls, import_file_to_cache
    from pyradiosky.utils import download_gleam
except ImportError as e:
    raise ImportError(
        "astropy, pyradiosky, astroquery are required for full script functionality"
    ) from e


# get the cache directory for pyuvsim as prescribed by astropy
cache_dir = get_cache_dir("pyuvsim")


# dictionary of keywords for downloadable files with url values (gleam we don't have a url)
file_download_dict = {
    # mwa uvbeam
    "mwa": "http://ws.mwatelescope.org/static/mwa_full_embedded_element_pattern.h5",
    # gsm 2016 nside 128 healpix map
    "healpix": "https://repository.library.brown.edu/storage/bdr:eafzyycj/content/",
    # gleam sky catalog as VOTable
    "gleam": "",
}


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
        print(
            f"astropy cached url for {url} already exists! If you want to redownload,"
            " please clear the cached url."
        )
        return

    print("gleam.vot download can take a while")
    # downloads gleam.vot using astroquery
    download_gleam(path=cache_dir, filename="gleam.vot")

    path_to_file = os.path.join(cache_dir, "gleam.vot")

    # converts downloaded gleam.vot to format recognizable for astropy with fake url key
    if os.path.exists(path_to_file):
        print("running import_file_to_cache to convert to astropy downloaded format")
        import_file_to_cache(url, path_to_file, remove_original=True, pkgname="pyuvsim")
    else:
        print(f"gleam.vot not downloaded to expected location: {path_to_file}")


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
        print(
            f"astropy cached url for {url} already exists! If you want to redownload,"
            " please clear the cached url."
        )
        return

    # astropy file download method which defaults to cached file and hashes the file
    # should have a default timeout of 10 seconds
    download_file(url, cache=True, pkgname="pyuvsim")


def main():
    """
    Run the main method as the entrypoint of download_data_files.py when ran directly.

    Parses the arguments provided using argparse. The arguments should be be space separated keys
    corresponding to the file_download_dict. Calls download_from_instructions for each valid key.
    """
    # create parser that grabs all input command line arguments as file keywords
    parser = argparse.ArgumentParser(
        description="parse a list of filenames by keyword to download"
    )

    parser.add_argument(
        "files",
        nargs="*",  # '*' consumes zero or more arguments, storing them in a list
        default=file_download_dict.keys(),  # If no arguments are provided, use the full list
        help="List of files to download via keyword. Downloads all available files if no arguments"
        f" given. Currently available file keywords: {' '.join(file_download_dict.keys())}."
        ' Sample use: "python3 download_data_files.py'
        f' {" ".join(list(file_download_dict.keys())[:2])}"',
    )

    args = parser.parse_args()

    files = args.files

    # loop through every file keyword in the list of filenames passed in from the command line
    # for each file keyword, check that it exists in the dictionary, then try to download
    for file in files:
        if file in file_download_dict:
            # call download_gleam_vot if gleam, otherwise download from url
            if file == "gleam":
                download_gleam_vot()
            else:
                download_file_using_astropy(file_download_dict[file])
        else:
            print(
                f'file "{file}" not found! Check the available keys in the file_download_dict'
                "or the listed current options at the top of the file."
            )


if __name__ == "__main__":
    main()
