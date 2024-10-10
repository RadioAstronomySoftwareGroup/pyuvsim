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
@pytest.mark.parametrize("paramfile", ["obsparam_ref_1.1_uniform.yaml"]) 
@pytest.mark.parallel(1)
@pytest.mark.skipif(not hasbench, reason="benchmark utility not installed")
def test_run_11_uniform(benchmark, goto_tempdir, paramfile):
    # Test vot and txt catalogs for parameter simulation
    # Compare to reference files.
    uv_ref = UVData()
    uv_ref.read_uvh5("/oscar/data/jpober/mburdorf/pyuvsim/results_data/ref_1.1_uniform.uvh5")

    param_filename = os.path.join(SIM_DATA_PATH, "test_config", paramfile)
    # This test obsparam file has "single_source.txt" as its catalog.
    #with warnings.catch_warnings():
    #    warnings.simplefilter("ignore")
    #	# not sure
    benchmark(
        pyuvsim.uvsim.run_uvsim,
        param_filename
    )
    #benchmark(print, "hi")
    pyuvsim.uvsim.run_uvsim(param_filename)

    # Loading the file and comparing is only done on rank 0.
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

    assert uv_new == uv_ref

@pytest.mark.skipif(not hasbench, reason="benchmark utility not installed")
def test_example(benchmark): #goto_tempdir, paramfile):
    #paramfile = "obsparam_ref_1.1_uniform.yaml" 
    # Test vot and txt catalogs for parameter simulation
    # Compare to reference files.
    uv_ref = UVData()
    benchmark(uv_ref.read_uvh5, os.path.join("/oscar/data/jpober/mburdorf/pyuvsim/results_data/ref_1.1_uniform.uvh5"))
