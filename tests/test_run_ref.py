# Copyright (c) 2020 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import importlib
import os

import pytest
from pyuvdata import UVData

from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
#from pyuvsim.mpi import rank
import pyuvsim
from pyuvsim.uvsim import run_uvsim

hasbench = importlib.util.find_spec("pytest_benchmark") is not None

pytest.importorskip("mpi4py")  # noqa

# gets latest pid from api call
def get_latest_pid(response):
    response_content_arr = response.json()['response']['docs']

    # construct an array of just the timestamps for object creation from all the search results
    time_arr = [item['object_created_dsi'] for item in response_content_arr]

    if len(time_arr) == 0:
        print("fail")
        return

    # get the max timestamp, then get the index at which the max time occurs
    # this corresponds to the latest created object
    latest_item_pos = time_arr.index(max(time_arr))
    # get the pid of the latest item
    latest_pid = response_content_arr[latest_item_pos]['pid']
    
    return latest_pid


# TODO: fix up order of operations and comment the method better
def download_sim(target_dir, sim_name):
    import requests
    import urllib.request

    api_url="https://repository.library.brown.edu/api/search/?q="
    response=requests.get(api_url+sim_name)

    target_dir = os.path.join(target_dir, "results_data")

    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    # set filename with full filepath
    fname = "ref_" + sim_name + ".uvh5"
    fname = os.path.join(target_dir, fname)

    pid = get_latest_pid(response)

    # download url
    download_url = f"https://repository.library.brown.edu/storage/{pid}/content/"

    # download the file to the location
    urllib.request.urlretrieve(download_url, fname)

    # additionally download mwa uvbeam to SIM_DATA_PATH if necessary
    if "mwa" in sim_name:
        download_url = "http://ws.mwatelescope.org/static/mwa_full_embedded_element_pattern.h5"
        fname=os.path.join(SIM_DATA_PATH,"mwa_full_embedded_element_pattern.h5")
        print(fname)

        # download the file to the location
        if os.path.isfile(fname):
            print("skipping download mwa uvbeam file already downloaded")
        else:
            urllib.request.urlretrieve(download_url, fname)


# TODO: make more complete
def compare_uvh5(uv_ref, uv_new):
    # fix part(s) that should deviate
    # TODO: maybe assert that deviation occurred
    uv_new.history = uv_ref.history

    # perform equality check between historical and current reference
    # simulation output
    # TODO: implement better asserts
    #       include lower tolerance for deviations
    assert uv_new == uv_ref

def construct_filepaths(target_dir, sim):
    # construct filename and filepath of downloaded file
    uvh5_filename = "ref_" + sim + ".uvh5"
    uvh5_filepath = os.path.join(target_dir, "results_data", uvh5_filename)

    # construct filename and filepath of yaml file used
    # to run simulation
    yaml_filename = "obsparam_ref_" + sim + ".yaml"
    yaml_filepath = os.path.join(SIM_DATA_PATH, "test_config", yaml_filename)
    assert os.path.exists(yaml_filepath)

    return uvh5_filepath, yaml_filepath


@pytest.fixture
def goto_tempdir(tmpdir):
    # Run test within temporary directory.
    newpath = str(tmpdir)
    cwd = os.getcwd()
    os.chdir(newpath)

    yield newpath

    os.chdir(cwd)


# TODO: FIXME: ASK WHY MWA 1.2 BEAM HARD TO RUN?
# TODO: FIXME: (TECHNICALLY GOES ABOVE BUT HAVE IT EXPLICITLY SEARCH COLLECTION AND NOT GENERAL SEARCH API FOR BETTER RESULTS)

# list of sims to benchmark pedantically as they are too long to run many times
# TODO: if need be swap to a dictionary with more specific custom pedantic arguments
long_ref_sims = ["1.2_uniform", "1.2_gauss"]

@pytest.mark.skipif(not hasbench, reason="benchmark utility not installed")
def test_run_sim(benchmark, goto_tempdir, refsim):
    # download reference sim output to compare
    download_sim(goto_tempdir, refsim)

    # construct paths to necessary file using paramfile
    uvh5_filepath, yaml_filepath = construct_filepaths(goto_tempdir, refsim)

    # instantiate UVData object and read in the downloaded uvh5
    # file as the point of historical comparison
    uv_ref = UVData.from_file(uvh5_filepath)

    # benchmark uvsim.run_uvsim method with param_filename argument
    # runs multiple times to get benchmark data
    # outdir is given by the yaml file but should not write anything
    # and simply return a UVData object
    # uses long_ref_sims global array to determine if pedantic benchmarking should be used
    if refsim in long_ref_sims:
        uv_new = benchmark.pedantic(
            run_uvsim,   
            args=[yaml_filepath],   
            kwargs={'return_uv':True},   
            setup=None,   
            rounds=1,   
            warmup_rounds=0,   
            iterations=1) 
    else:
        uv_new = benchmark(run_uvsim, yaml_filepath, return_uv=True)

    # loading the file and comparing is only done on rank 0.
    # probably not necessary
    if pyuvsim.mpi.rank != 0:
        return

    compare_uvh5(uv_ref, uv_new)
