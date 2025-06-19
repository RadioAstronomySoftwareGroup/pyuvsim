# Copyright (c) 2020 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

"""File to download ref sims."""

import argparse
import json
import time

try:
    from pooch import os_cache, retrieve
except ImportError as e:
    raise ImportError("pooche required for full script functionality") from e


# to test for now but swap to new reference simulations after
ref_sims = [
    "1.1_uniform",
    "1.1_gauss",
    "1.1_mwa",
    "1.2_uniform",
    "1.2_gauss",
    "1.3_uniform",
    "1.3_gauss",
]


json_cache = os_cache("pyuvsim/json")


def main():
    """
    Run the main method as the entrypoint of download_ref_sims.py when ran directly.

    Parses the arguments provided using argparse. The arguments should be be space separated keys
    corresponding to reference simulation data files stored on the Brown Digital Repository.

    """
    # create parser that grabs all input command line arguments as file keywords
    parser = argparse.ArgumentParser(
        description="a list of reference simulations to download"
    )

    # TODO: flip flop on if default is caching or not, and similarly on if a script available for
    #       simulation comparison should grab cache or local
    # TODO: the comparison script should honestly search both
    parser.add_argument(
        "--cache",
        action="store_true",
        help="Caches files to pooch pyuvsim cache if passed, otherwise downloads to local"
        " directory.",
    )

    # TODO: fix this up a bit
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

    args = parser.parse_args()

    to_cache = args.cache
    sims = args.sims

    file_path = "latest_ref_sims"
    if to_cache:
        file_path = os_cache("pyuvsim/latest_ref_data")

    # TOOD: make this and other file main method
    for sim in sims:
        if sim in ref_sims:
            download_sim(sim, fpath=file_path)
        else:
            print(
                f'sim "{sim}" not found! Check the available keys in the ref_sims list'
                " or the listed current options at the top of the file."
            )


def download_sim(sim_name, fpath="."):
    """
    Download the historical reference simulations from the Brown Digital Repository.

    Sends an api call to the "pyuvsim historical reference simulations" collection,
    then identifies the latest uploaded object in the response. Downloads that object to the
    target directory. API response and download uses pooch.

    Parameters
    ----------
    sim_name: str
        A a string key matching the reference simulation data stored on the Brown Digital
        Repository.
    fpath: str
        The file path to download the file to.

    """
    # Link to BDR API DOCS:
    # https://repository.library.brown.edu/studio/api-docs/

    # api url
    api_url = (
        "https://repository.library.brown.edu/api/collections/bdr:wte2qah8/?q="
        f"primary_title:{sim_name}&wt=json&sort=object_last_modified_dsi desc&rows=1"
    )

    print(
        f"\nquerying BDR collection for latest item matching {sim_name} via api: {api_url}"
    )
    api_response_path = retrieve(url=api_url, path=json_cache, known_hash=None)
    # api_response_path = retrieve(url=api_url, path=json_cache, known_hash=None, progressbar=True)
    time.sleep(1)

    # attempt to parse json assuming we have the result
    try:
        # read downloaded file as a json and parse for pid
        with open(api_response_path) as file:
            api_response_json = json.load(file)
        pid = api_response_json["items"]["docs"][0]["pid"]
    except Exception as e:
        print(f"Failed to parse BDR response for file PID with error: {e}")

    # download url
    download_url = f"https://repository.library.brown.edu/storage/{pid}/content/"

    print(f"attempting download of requested file by http: {download_url}\n")
    print(f"downloading file to {fpath}")
    file_download_path = retrieve(
        url=download_url, fname=f"ref_{sim_name}.uvh5", path=fpath, known_hash=None
    )
    print(f"downloaded file to {file_download_path}")
    # file_download_path = retrieve(url=download_url, fname=f"ref_{sim_name}.uvh5", path=fpath,
    # known_hash=None, progressbar=True)
    time.sleep(1)


if __name__ == "__main__":
    main()
