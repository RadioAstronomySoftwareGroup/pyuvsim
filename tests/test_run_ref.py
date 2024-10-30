# Copyright (c) 2020 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import importlib
import os

import pytest
from pyuvdata import UVData

from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
from pyuvsim.mpi import rank
from pyuvsim.uvsim import run_uvsim

hasbench = importlib.util.find_spec("pytest_benchmark") is not None

pytest.importorskip("mpi4py")  # noqa

# TODO: get running for current 8 sims (add necessary data files)
# AND DO THE FIX FOR THE GOOGLE DRIVE (UPLOAD THE FILES AND
# GET THEM TO DOWNLOAD (I GUESS TEMP JUST MAKE ANOTHER
# DICT OR FILE OR DATA STRUCTURE THAT CAN GRAB THE FILES)

# TODO: check syntax preference for global lists!
ci_ref_sims = [
    "1.1_uniform",
    "1.1_gauss",
    #"1.1_hera",
    #"1.1_gauss",
    #"1.2_uniform",
    #"1.2_hera",
    #"1.3_uniform",
    #"1.3_gauss",
]


# TODO: hardset it using dictionary or something
# TODO: swap to something more permanent!
@pytest.fixture
def download_sims():
    import requests

    target_dir = "results_data"

    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    # temporary way to store all the fileids
    # control which simulations via global ci_ref_sims array
    fids = {
        "1.1_uniform":"1V4RmmUGrx5iH-Zyuj45p1ag3vmMuUI2A",
        "1.1_gauss":"1gTj9BSunEoyde5Dq8a5IhUieMo9PuMuC",
        #"1.1_hera":"13B0YOSO9v4KsEGX9GVHuBEQrt5Tnjsxr",
        #"1.2_gauss":a,
        #"1.2_uniform":a,
        #"1.2_hera":a,
        #"1.3_uniform":a,
        #"1.3_gauss":a
    }


    urlbase = "https://drive.google.com/uc?export=download"
    
    # for each sim name in ci_ref_sims, checks that needs hera uvbeam and
    # has downloadable data (is a key to fids)
    download_hera_uvbeam = any(["hera" in sim for sim in ci_ref_sims if sim in fids])
    print(download_hera_uvbeam)

    # download the hera uvbeam data file if we need it
    if download_hera_uvbeam:
        print("skipping as file cannot currently be downloaded")

        #fid = "1lqkLlnB3uE17FcPB2GcJJQoCcn7-wVP8"
        #fname = os.path.join(target_dir, "HERA_NicCST_fullfreq.uvbeam")
        #r = requests.get(urlbase, params={"id": fid})
        # write the file
        #with open(fname, "wb") as ofile:
        #    ofile.write(r.content)

    # TODO: upload all other files and hardcode
    #       them via a dictionary with ci_ref_sims for keys
    for sim in ci_ref_sims:
        if sim not in fids.keys():
            print(f"key {sim} not found!")
            continue            

        # get id from dict
        fid = fids[sim]
        
        # get content
        r = requests.get(urlbase, params={"id": fid})

        # set filename with full filepath
        fname = "ref_" + sim + ".uvh5"
        fname = os.path.join(target_dir, fname)

        # write the file
        with open(fname, "wb") as ofile:
            print(f"writing {fname}")
            ofile.write(r.content)


@pytest.fixture
def goto_tempdir(tmpdir):
    # Run test within temporary directory.
    newpath = str(tmpdir)
    cwd = os.getcwd()
    os.chdir(newpath)

    yield newpath

    os.chdir(cwd)


# TODO: add all 8 yaml files
@pytest.mark.parametrize("paramfile", ci_ref_sims)
@pytest.mark.skipif(not hasbench, reason="benchmark utility not installed")
def test_run_11_uniform(benchmark, goto_tempdir, download_sims, paramfile):
    # download reference sim
    # download_sims(goto_tempdir)

    # construct filename and filepath of downloaded file
    uvh5_filename = "ref_" + paramfile + ".uvh5"
    uvh5_filepath = os.path.join(goto_tempdir, "results_data", uvh5_filename)

    # instantiate UVData object and read in the downloaded uvh5
    # file as the point of historical comparison
    uv_ref = UVData.from_file(uvh5_filepath)

    # construct filename and filepath of yaml file used
    # to run simulation
    yaml_filename = "obsparam_ref_" + paramfile + ".yaml"
    yaml_filepath = os.path.join(SIM_DATA_PATH, "test_config", yaml_filename)

    # benchmark uvsim.run_uvsim method with param_filename argument
    # runs around 10 times to get benchmark data
    # outdir is given by the yaml file but should be current working directory
    # for all the reference simulations
    # TODO: think more about filepath specification, also there was a warning
    # TODO: check if somehow making like 10 copies of output accidentally (whoops)
    benchmark(run_uvsim, yaml_filepath)

    print(f"\nfor run with {paramfile}:\n")
    print(f"\n\ncurrent working directory is: {os.getcwd()}\n\n")
    print(f"\n\nfiles: {os.listdir()}\n\n")

    # loading the file and comparing is only done on rank 0.
    if rank != 0:
        return

    # instantiate new UVData object from the benchmark run_uvsim output
    # as current point of comparison
    new_uvh5_filepath = os.path.join(goto_tempdir, uvh5_filename)
    uv_new = UVData.from_file(new_uvh5_filepath)

    # fix part(s) that should deviate
    # TODO: maybe assert that deviation occurred
    uv_new.history = uv_ref.history

    # perform equality check between historical and current reference
    # simulation output
    # TODO: implement better asserts
    #       include lower tolerance for deviations
    assert uv_new == uv_ref
