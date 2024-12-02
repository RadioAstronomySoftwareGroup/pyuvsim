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

#TODO: make downloading simulation output file and companion files more robust (add retries, etc...)

# gets latest pid from api call
def get_latest_pid(response):
    import requests

    # list of item dicts
    collection_response_items  = response.json()['items']['docs']

    if len(collection_response_items) == 0:
        return "No items found in collection matching query."

    # using "object_created_dsi" key to sort the items, so need to request that
    # for each item via the "json_uri", and get the pid as well to return
    print(f"requesting json_uri for each item in reponse, and parsing 'object_created_dsi' and 'pid'")
    json_uris = [ item['json_uri'] for item in collection_response_items ]
    object_created_dsis = []
    pids = []
    for uri in json_uris:
        # get response for item as json
        response_json = requests.get(uri).json()

        # get values from json
        time_created = response_json['object_created_dsi']
        item_pid = response_json['pid']

        # append to lists
        object_created_dsis.append(time_created)
        pids.append(item_pid)

    # get the max timestamp, then get the index at which the max time occurs
    # this corresponds to the latest created object
    latest_item_pos = object_created_dsis.index(max(object_created_dsis))
    
    # get the pid of the latest item by shared pos
    latest_pid = pids[latest_item_pos]

    print("returning pid of most recent file uploaded to the collection matching the query")
    
    return latest_pid


# TODO: fix up order of operations and comment the method better
def download_sim(target_dir, sim_name):
    import requests
    import urllib.request

    api_url="https://repository.library.brown.edu/api/collections/bdr:wte2qah8/?q="
    
    print(f"\nquerying BDR collection for items matching {sim_name} via api: {api_url+sim_name}")
    response=requests.get(api_url+sim_name)

    # parse out the latest file in the collection from the search result and return its pid
    pid = get_latest_pid(response)

    # append results_data to the target directory download path
    target_dir = os.path.join(target_dir, "results_data")

    # check if we need to make target_dir
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    # set filename with full filepath
    fname = "ref_" + sim_name + ".uvh5"
    fname = os.path.join(target_dir, fname)

    # download url
    download_url = f"https://repository.library.brown.edu/storage/{pid}/content/"

    print(f"attempting download of requested file by http: {download_url}\n")
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
    import numpy as np

    # print histories
    print("History 1:")
    print(uv_ref.history)
    print("")

    print("History 2:")
    print(uv_new.history)
    print("")

    ref_arr, new_arr = uv_ref.data_array, uv_new.data_array

    # mean diff
    mean_diff_of_vis = np.mean(np.abs(new_arr - ref_arr))
    mean_diff_of_abs = np.mean(np.abs(ref_arr) - np.abs(new_arr))

    # compute max absolute difference and max relative difference
    # should be identical to the allclose output for maximum absolute
    # and relative difference
    max_absolute_diff=np.amax(np.abs(new_arr - ref_arr))
    frac_diff=frac_diff=np.abs(new_arr - ref_arr)/np.abs(ref_arr)
    frac_diff[np.isnan(frac_diff)]=0
    max_relative_diff=np.amax(frac_diff)

    # generate a true/false ndarray for passing and failing cases
    # should match output of np.testing.assert_allclose 
    cases = np.abs(new_arr - ref_arr) <= (1e-8 + 1e-5 * np.abs(ref_arr))
    outcome = cases.all()
    num_mismatched = str(len(cases[cases == False])) + "/" + str(cases.size)

    # print some things for reference
    print(f"mean of abs of diff of visibilities 'mean(abs(old_data_arr - new_data_arr)): {mean_diff_of_vis}")
    print(f"mean of diff of abs of visibilities 'mean(abs(old_data_arr) - abs(new_data_arr)): {mean_diff_of_abs}")
    print(f"max_absolute_diff: {max_absolute_diff}")
    print(f"max_relative_diff: {max_relative_diff}")
    print(f"outcome: {outcome}")
    print(f"num_mismatched: {num_mismatched}")

    # perform basic equality check between historical and current reference
    # simulation output
    assert uv_new == uv_ref

    # perform more intensive check
    # TODO

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


# list of sims to benchmark pedantically as they are too long to run many times
# if need be can swap to a dictionary with more specific custom pedantic arguments
long_ref_sims = ["1.2_uniform", "1.2_gauss"]

# TODO: generic comment on top about how this is used (conftest command line / workflow setup)
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
