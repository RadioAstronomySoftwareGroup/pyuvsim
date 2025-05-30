# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License
import atexit
import os
import re
import shutil

import pytest
from numpy import unique
from pyuvdata.testing import check_warnings

import pyuvsim
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH


def profdata_dir_setup(tmpdir):
    # The pytest temporary test directory doesn't work with the profiling outputs
    # because they are only created after the exit functions are run, and pytest
    # deletes its temporary directory before then.
    # As a workaround, make a temporary directory for profiling data and register
    # its deletion in the exit functions.
    outpath = tmpdir.mkdir("temporary_profiling_data")
    atexit.register(shutil.rmtree, outpath)
    return outpath


@pytest.mark.filterwarnings("ignore:The mount_type parameter must be set for UVBeam")
def test_profiler(tmpdir):
    line_profiler = pytest.importorskip("line_profiler")
    outpath = profdata_dir_setup(tmpdir)
    testprof_fname = str(outpath.join("time_profile.out"))
    pyuvsim.profiling.set_profiler(outfile_prefix=testprof_fname, dump_raw=True)
    with check_warnings(UserWarning, match="Profiler already set"):
        pyuvsim.profiling.set_profiler(
            outfile_prefix=testprof_fname[:-4], dump_raw=True
        )
    param_filename = os.path.join(
        SIM_DATA_PATH, "test_config", "param_1time_1src_testcat.yaml"
    )
    pyuvsim.uvsim.run_uvsim(param_filename, return_uv=True)
    time_profiler = pyuvsim.profiling.get_profiler()
    time_profiler.disable_by_count()
    assert isinstance(time_profiler, line_profiler.LineProfiler)
    assert hasattr(time_profiler, "rank")
    assert hasattr(time_profiler, "meta_file")
    lstats = time_profiler.get_stats()
    assert len(lstats.timings) != 0
    func_names = [k[2] for k in lstats.timings]
    assert unique(func_names).tolist() == sorted(
        pyuvsim.profiling.default_profile_funcs
    )


def test_profiler_mock_import(tmpdir):
    try:
        from line_profiler import LineProfiler  # noqa
    except ImportError:
        outpath = profdata_dir_setup(tmpdir)
        testprof_fname = str(outpath.join("time_profile.out"))

        with pytest.raises(
            ImportError,
            match=re.escape(
                "You need mpi4py and line_profiler to use the "
                "profiling module. Install them both by running pip "
                "install pyuvsim[all]."
            ),
        ):
            pyuvsim.profiling.set_profiler(outfile_prefix=testprof_fname, dump_raw=True)

        lp = pyuvsim.profiling.LineProfiler()
        assert lp is None
