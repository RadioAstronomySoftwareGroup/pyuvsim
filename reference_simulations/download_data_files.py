#!/bin/env python
# Copyright (c) 2022 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License
"""
An extensible approach to large file downloading with variable methods -- beams and catalogs.

Takes a command line argument of the file to download.

Current options:
    - "gleam": the gleam catalog as a VOTable
    - "mwa": mwa uvbeam file
    - "healpix": gsm 2016 nside 128 healpix map saved as skyh5

gleam requires pyradiosky and astroquery, mwa and healpix require pooch.

downloads files to pooch default file download location: "~/.cache/pooch/"

approach inspired by: https://stackoverflow.com/questions/67790907/how-to-call-a-function-stored-as-variable-with-arguments-python-3
"""

import argparse

try:
    from pyradiosky.utils import download_gleam
    from pooch import retrieve, os_cache
except ImportError as e:
    raise ImportError(
        "pooch, pyradiosky, astroquery are required for full script functionality"
    ) from e

# TODO: add historical output as well (have them download to file_cache/reference_data or something

# TODO: need to figure out if there is a more proper "pooch" approach to downloading and storing
#       the files even when pooch isn't used to download the file

# TODO: maybe change to pyuvsim_reference_simulations?
# TODO: need to figure out how to properly configure this target via yaml files
file_cache = os_cache("pyuvsim")

# dictionary of keywords for downloadable files
# each downloadable file has a dictionary containing:
# - the function used to download the file
# - the args to pass to the function
# - the kwargs to pass to the function
file_download_dict = {
    # mwa uvbeam
    "mwa" : {
             "function": retrieve,
             "args": [],
             "kwargs": {
                        "url": "http://ws.mwatelescope.org/static/mwa_full_embedded_element_pattern.h5",
                        "fname": "mwa_full_embedded_element_pattern.h5",
                        "known_hash": "a7649c6e03b8128a1de4614c2f363af5fa44f3890ae27bf893d56ca337bc48ee",
                        #"progressbar": True,
                        "path": file_cache
             }
    },
    # gsm 2016 nside 128 healpix map
    "healpix" : {
             "function": retrieve,
             "args": [],
             "kwargs": {
                        "url": "https://repository.library.brown.edu/storage/bdr:eafzyycj/content/",
                        "fname": "gsm16_nside128_100mhz.skyh5",
                        "known_hash": "cedbbe5180830e69cefe72c3a1642a8205dc2d595e6bc0a35c08e5948b12290d",
                        #"progressbar": True,
                        "path": file_cache
             }
    },
    # gleam sky catalog as VOTable
    "gleam" : {
             "function": download_gleam,
             "args": [],
             "kwargs": {
                        "path": file_cache
             }
    }
}

def main():
    # create parser that grabs all input command line arguments as file keywords
    parser = argparse.ArgumentParser(description="parse a list of filenames by keyword to download")

    parser.add_argument(
        "files",
        nargs="*",  # '*' consumes zero or more arguments, storing them in a list
        default=file_download_dict.keys(),  # If no arguments are provided, use the full list
        help="List of files to download via keyword. Downloads all available files if no arguments"
             f" given. Currently available file keywords: {" ".join(file_download_dict.keys())}."
             f" Sample use: \"python3 download_data_files.py {" ".join(list(file_download_dict.keys())[:2])}\""
    )

    args = parser.parse_args()

    files = args.files

    # loop through every file keyword in the list of filenames passed in from the command line
    # for each file keyword, check that it exists in the dictionary, then try to download
    for file in files:
        if file in file_download_dict:
            # get download instructions for file
            download_instructions = file_download_dict[file]
            # call method which construct download call from instruction
            download_from_instructions(file, download_instructions)
        else:
            print(f"file \"{file}\" not found! Check the available keys in the file_download_dict"
                   "or the listed current options at the top of the file."
            )

def download_from_instructions(filename, instruction):
    # takes as input the filename and a dictionary containing three keys:
    # - function: the method to call to download the file(s)
    # - args: a list of arguments to pass to the method in order
    # - kwargs: a list of keyword arguments to pass to the method
    # then calls the method using the args and kwargs

    # parse download instructions
    download_method = instruction.get("function")
    args = instruction.get("args", [])
    kwargs = instruction.get("kwargs", {})

    if "path" in kwargs:
        print(f"downloading {filename} to {kwargs["path"]}")
    else:
        print(f"download path for {filename} not specified in dictionary")

    if filename == "gleam":
        print("gleam filename download can take a while")

    # call the method with args and kwargs to download the file
    download_method(*args, **kwargs)

if __name__ == "__main__":
    main()
