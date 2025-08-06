# Copyright (c) 2025 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import pytest

from pyuvsim import cli


def test_run_pyuvsim_errors():
    with pytest.raises(
        ValueError,
        match="Either pass a parameter file or all of: uvdata, skymodel, uvbeam "
        "files and outfile.",
    ):
        cli.run_pyuvsim(["--uvdata", "foo"])

    with pytest.raises(
        ValueError,
        match="A parameter file and was passed along with one or more of uvdata, "
        "skymodel, uvbeam and outfile. Either pass a parameter file or all of: "
        "uvdata, skymodel, uvbeam files and outfile.",
    ):
        cli.run_pyuvsim(["--uvdata", "foo", "--param", "bar"])


def test_uvdata_to_telescope_config_errors():
    with pytest.raises(ValueError, match="beamfile must be passed."):
        cli.uvdata_to_telescope_config(["foo"])
