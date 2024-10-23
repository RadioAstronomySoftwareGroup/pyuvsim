# Copyright (c) 2020 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import os
#import warnings


import pytest
from pyuvdata import UVData

import pyuvsim
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH

import importlib
hasbench = importlib.util.find_spec("pytest_benchmark") is not None

pytest.importorskip("mpi4py")  # noqa

# TODO: get running for current 8 sims (add necessary data files)
# AND DO THE FIX FOR THE GOOGLE DRIVE (UPLOAD THE FILES AND 
# GET THEM TO DOWNLOAD (I GUESS TEMP JUST MAKE ANOTHER 
# DICT OR FILE OR DATA STRUCTURE THAT CAN GRAB THE FILES)

# TODO: check syntax preference for global lists!
ci_ref_sims = [
               "1.1_uniform",
               "1.1_hera",
               "1.1_gauss",
               "1.2_uniform",
               "1.2_hera",
               "1.3_uniform",
               "1.3_gauss"
              ]

# TODO: swap to something more permanent!
@pytest.fixture
def download_sims():
    import requests
    target_dir = "results_data" 

    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    urlbase = "https://drive.google.com/uc?export=download"

    dat = []

    fid = "1V4RmmUGrx5iH-Zyuj45p1ag3vmMuUI2A"
    fname = "ref_1.1_uniform.uvh5"
    # fileid name type size size_unit date time

    r = requests.get(urlbase, params={"id": fid})

    fname = os.path.join(target_dir, fname)

    with open(fname, "wb") as ofile:
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
    #download_sims(goto_tempdir)

    # construct filename and filepath of downloaded file
    uvh5_filename = "ref_" + paramfile + ".uvh5"
    uvh5_filepath = os.path.join(goto_tempdir, "results_data", uvh5_filename)

    # instantiate UVData object and read in the downloaded uvh5 
    # file as the point of historical comparison
    uv_ref = UVData()
    uv_ref.read_uvh5(filepath)

    # construct filename and filepath of yaml file used 
    # to run simulation
    yaml_filename = "obsparam_ref_" + paramfile + ".yaml"
    yaml_filepath = os.path.join(SIM_DATA_PATH, "test_config", yaml_filename)

    # benchmark uvsim.run_uvsim method with param_filename argument
    # runs around 10 times to get benchmark data
    # outdir is given by the yaml file but should be current working directory
    # for all the reference simulations
    # TODO: think more about filepath specification, also there was a warning I can get rid of probably
    # TODO: check if somehow making like 10 copies of output accidentally (whoops)
    benchmark(
        pyuvsim.uvsim.run_uvsim,
        yaml_filepath
    )

    # loading the file and comparing is only done on rank 0.
    if pyuvsim.mpi.rank != 0:
        return

    # instantiate new UVData object from the benchmark run_uvsim output
    # as current point of comparison
    path = goto_tempdir
    ofilepath = os.path.join(path, uvh5_filename)
    uv_new = UVData.from_file(ofilepath)

    # reset part(s) that should deviate 
    # TODO: maybe assert that deviation occurred 
    uv_new.history = uv_ref.history

    # perform equality check between historical and current reference
    # simulation output
    # TODO: implement better asserts
    #       include lower tolerance for deviations
    assert uv_new == uv_ref
