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

#paramfile (fix name) = ["a", "b",]
# FIXME: example for that in lunar stuff (pyuvdata in pyuvdata/test_utils/test_coordinates.py (maybe inaccurate locationing))

# conftest implementation -- setup and teardown happens here
# fixture can run previous fixture
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


#@pytest.mark.filterwarnings("ignore:antenna_diameters are not set")
#@pytest.mark.filterwarnings("ignore:Telescope Triangle is not in known_telescopes.")
#@pytest.mark.filterwarnings("ignore:Fixing auto-correlations to be be real-only")
#@pytest.mark.parallel(1)
# to test: call pytest with mpiexec / similar and see if runs
@pytest.mark.parametrize("paramfile", ["obsparam_ref_1.1_uniform.yaml"]) 
@pytest.mark.skipif(not hasbench, reason="benchmark utility not installed")
def test_run_11_uniform(benchmark, goto_tempdir, download_sims, paramfile):
    # download reference sim
    #download_sims(goto_tempdir)

    filepath = os.path.join(goto_tempdir, "results_data", "ref_1.1_uniform.uvh5")
    # Test vot and txt catalogs for parameter simulation
    # Compare to reference files.
    uv_ref = UVData()
    # TODO: implement methods to make this download from remote to appropriate path 
    # for now from google drive 

    uv_ref.read_uvh5(filepath)

    param_filename = os.path.join(SIM_DATA_PATH, "test_config", paramfile)
    # This test obsparam file has "single_source.txt" as its catalog.
    #with warnings.catch_warnings():
    #    warnings.simplefilter("ignore")
    #	# not sure
    benchmark(
        pyuvsim.uvsim.run_uvsim,
        param_filename
    )
    #pyuvsim.uvsim.run_uvsim(param_filename)

    # Loading the file and comparing is only done on rank 0.
    # maybe comment out
    if pyuvsim.mpi.rank != 0:
        return

    path = goto_tempdir
    ofilepath = os.path.join(path, "ref_1.1_uniform.uvh5")
    uv_new = UVData.from_file(ofilepath)

    # Reset parts that will deviate
    uv_new.history = uv_ref.history
    #uv_ref.dut1 = uv_new.dut1
    #uv_ref.gst0 = uv_new.gst0
    #uv_ref.rdate = uv_new.rdate

    # TODO: implement better asserts

    assert uv_new == uv_ref
