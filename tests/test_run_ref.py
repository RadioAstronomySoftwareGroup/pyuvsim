# Copyright (c) 2020 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import importlib
import os

import pytest
from pyuvdata import UVData

import pyuvsim
from pyuvsim.cli import (
    download_data_files,
    download_ref_sims,
    get_latest_api_response_pid,
    run_ref_sim,
)
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
from pyuvsim.uvsim import run_uvsim

hasbench = importlib.util.find_spec("pytest_benchmark") is not None

pytest.importorskip("mpi4py")  # noqa


# ONLY RUN ON download_test workflow accomplished by requiring pytest-benchmark
@pytest.mark.skipif(not hasbench, reason="benchmark utility not installed")
@pytest.mark.filterwarnings("ignore:gleam.vot download")
@pytest.mark.parametrize(
    "files",
    [
        ("--clear", "--row_limit", "75001", "gleam", 1),
        ("mwa", 1),
        ("healpix", 1),
        ("fail", 1),
        ("gleam", 2),
    ],
)
def test_download_data_files(files):
    # NOTE: Not covering the case where files are already downloaded and the tests run.
    #       In that case additional warnings will occur.

    # check for duplicate download of gleam and bad file input to see warnings
    if files[-1] == 2:
        with pytest.warns(UserWarning, match=r"astropy cached url for .*"):
            download_data_files([*files[:-1]])
    elif files[0] == "fail":
        with pytest.warns(UserWarning, match=r'file "fail" not found! .*'):
            download_data_files([*files[:-1]])
    else:
        download_data_files([*files[:-1]])


@pytest.mark.filterwarnings("ignore:Some Stokes I values are negative")
@pytest.mark.filterwarnings("ignore:astropy cached url")
@pytest.mark.filterwarnings("ignore:The check_broadcast function")
@pytest.mark.parametrize(
    "ref_sim",
    [
        "1.1_baseline_number",
        "2",
        "1.3_frequency_axis",
        "4",
        "1.5_uvbeam",
        "6",
        "1.7_multi_beam",
        "8",
        "fail",
    ],
)
# ONLY RUN ON download_test workflow accomplished by requiring pytest-benchmark
@pytest.mark.skipif(not hasbench, reason="benchmark utility not installed")
def test_run_ref_sim_cli(ref_sim):
    if ref_sim == "4":
        # decided to clear the cache here for more consistent performance
        download_data_files(["--clear", "gleam", "--row_limit", "1000"])
    if ref_sim == "1.5_uvbeam":
        download_data_files(["mwa"])
    if ref_sim == "6":
        download_data_files(["healpix"])

    if ref_sim == "fail":
        with pytest.raises(ValueError, match=r"No match found for .*"):
            run_ref_sim([ref_sim])
    else:
        run_ref_sim([ref_sim])


# ONLY RUN ON download_test workflow accomplished by requiring pytest-benchmark
@pytest.mark.skipif(not hasbench, reason="benchmark utility not installed")
@pytest.mark.filterwarnings("ignore:astropy cached url")
@pytest.mark.parametrize(
    "ref_sim",
    [
        "1.1_baseline_number",
        "1.2_time_axis",
        "1.3_frequency_axis",
        "1.4_source_axis",
        "1.5_uvbeam",
        "1.6_healpix",
        "1.7_multi_beam",
        "1.8_lunar",
        "fail",
        "1.1_baseline_number",  # to check downloading a file that already exists
    ],
)
def test_download_ref_sims(ref_sim):
    # catch a hard to hit piece of code that requires failed downloading from BDR
    if ref_sim == "fail":
        with pytest.warns(UserWarning, match=r'sim "fail" not found! .*'):
            download_ref_sims([ref_sim])
        with (
            pytest.raises(UnboundLocalError) as e_info,
            pytest.warns(UserWarning, match=r"Failed to parse BDR response .*"),
        ):
            get_latest_api_response_pid(ref_sim)
        # after second with statement resolves we check the returned error
        assert (
            e_info.value.args[0] == "cannot access local variable 'pid' "
            "where it is not associated with a value"
        )
    else:
        download_ref_sims([ref_sim])


def compare_uvh5(uv_ref, uv_new):
    # takes as input two UVData objects, and computes relevant quantities for determining how
    # similar the data are. Prints the histories before setting them equal.
    import numpy as np

    # print histories
    print("History 1:")
    print(uv_ref.history)
    print("")

    print("History 2:")
    print(uv_new.history)
    print("")

    # set history to match so equality check doesn't fail
    uv_new.history = uv_ref.history

    # print extra_keywords
    print("extra_keywords 1:")
    print(uv_ref.extra_keywords)
    print("")

    print("extra_keywords 2:")
    print(uv_new.extra_keywords)
    print("")

    # set extra_keywords to match so equality check doesn't fail
    uv_new.extra_keywords = uv_ref.extra_keywords

    print("phase_center_catalog 1:")
    uv_ref.print_phase_center_info()
    print("")

    print("phase_center_catalog 2:")
    uv_new.print_phase_center_info()
    print("")

    # set phase_center_catalog name to match so equality check doesn't fail
    uv_ref.rename_phase_center(0, "unprojected")

    print("mount type 1:")
    print(uv_ref.telescope.mount_type)
    print("")

    print("mount type 2:")
    print(uv_new.telescope.mount_type)
    print("")

    # set telescope mount_type to match so equality check doesn't fail
    uv_ref.telescope.mount_type = uv_new.telescope.mount_type

    ref_arr, new_arr = uv_ref.data_array, uv_new.data_array

    # mean diff
    mean_diff_of_vis = np.mean(np.abs(new_arr - ref_arr))
    mean_diff_of_abs = np.mean(np.abs(ref_arr) - np.abs(new_arr))

    # compute max absolute difference and max relative difference
    # should be identical to the allclose output for maximum absolute
    # and relative difference
    max_absolute_diff = np.amax(np.abs(new_arr - ref_arr))
    frac_diff = np.abs(new_arr - ref_arr) / np.abs(ref_arr)
    frac_diff[np.isnan(frac_diff)] = 0
    max_relative_diff = np.amax(frac_diff)

    # using np.testing.assert_allclose defaults for comparison
    rtol = 1e-7
    atol = 0

    # generate a true/false ndarray for passing and failing cases
    # should match output of np.testing.assert_allclose
    cases = np.abs(new_arr - ref_arr) <= (atol + rtol * np.abs(ref_arr))
    outcome = cases.all()

    # get unique outcomes (true / false) and corresponding counts
    # then convert to dict and get result
    unique, counts = np.unique(cases, return_counts=True)
    outcome_dict = dict(zip(unique, counts, strict=False))

    # need to check that key exists
    if False in outcome_dict:
        num_mismatched = str(outcome_dict[False]) + "/" + str(cases.size)
    else:
        num_mismatched = "0" + "/" + str(cases.size)

    # print some things for reference
    print(
        f"mean of abs of diff of visibilities "
        f"'mean(abs(old_data_arr - new_data_arr))': {mean_diff_of_vis}"
    )
    print(
        f"mean of diff of abs of visibilities "
        f"'mean(abs(old_data_arr) - abs(new_data_arr))': {mean_diff_of_abs}"
    )
    print(f"max_absolute_diff: {max_absolute_diff}")
    print(f"max_relative_diff: {max_relative_diff}")
    print(f"outcome: {outcome}")
    print(f"num_mismatched: {num_mismatched}")

    # set UVParameter tols to local value for overloaded equality check
    uv_new._data_array.tols = (rtol, atol)

    # perform equality check between historical and current reference
    # simulation output
    assert uv_new == uv_ref


def construct_yaml_filepath(sim):
    # takes as input the sim name (NEEDS TO BE AN EXISTING SIM IN THE DATA DIRECTORY), then
    # constructs the expected yaml_filepath to run the simulation

    # yaml filepath construction
    yaml_filename = "obsparam_ref_" + sim + ".yaml"
    yaml_filepath = os.path.join(SIM_DATA_PATH, "test_config", yaml_filename)

    # if yaml_filepath does not exist then comparison should fail
    assert os.path.exists(yaml_filepath)

    return yaml_filepath


@pytest.mark.skipif(not hasbench, reason="benchmark utility not installed")
@pytest.mark.filterwarnings("ignore:astropy cached url for")
def test_run_sim(benchmark, refsim, savesim):
    # pytest method to benchmark reference simulations. currently only called on one reference
    # simulation at a time. takes as input the benchmark fixture, a fixture to generate a temporary
    # directory for testing, and a fixture defined in conftest.py that is used to parametrize this
    # method with a specific reference simulation name via command line input.
    if savesim:
        print(f"\nSaving reference simulation for {refsim} to cwd")

    # TODO: decide if better to have keyword as refsim name to begin with to avoid having to match
    large_data_files = {
        "1.4_source_axis": "gleam",
        "1.5_uvbeam": "mwa",
        "1.6_healpix": "healpix",
    }

    if refsim in large_data_files:
        keyword = large_data_files[refsim]
        download_data_files([keyword])

    # download reference sim output to compare
    uvh5_filepath = download_ref_sims([refsim])

    # construct paths to necessary file using paramfile
    yaml_filepath = construct_yaml_filepath(refsim)

    # instantiate UVData object and read in the downloaded uvh5
    # file as the point of historical comparison
    uv_ref = UVData.from_file(uvh5_filepath, file_type="uvh5")

    # benchmark uvsim.run_uvsim method with param_filename argument
    # runs multiple times to get benchmark data
    # outdir is given by the yaml file but should not write anything
    # and simply return a UVData object
    # currently always uses pedantic with 2 rounds 1 iteration
    uv_new = benchmark.pedantic(
        run_uvsim,
        args=[yaml_filepath],
        kwargs={"return_uv": True},
        setup=None,
        rounds=2,
        warmup_rounds=0,
        iterations=1,
    )

    # loading the file and comparing is only done on rank 0.
    # probably not necessary
    if pyuvsim.mpi.rank != 0:
        return

    # save the output to the temporary directory if savesim
    if savesim:
        savepath = os.path.join(os.getcwd(), "new_data")
        if not os.path.exists(savepath):
            os.makedirs(savepath)
        uv_new.write_uvh5(os.path.join(savepath, f"ref_{refsim}.uvh5"), clobber=True)

    # performs any assertions to confirm that the reference simulation output hasn't diverged
    compare_uvh5(uv_ref, uv_new)
