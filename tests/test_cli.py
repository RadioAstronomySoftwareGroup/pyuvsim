# Copyright (c) 2025 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License
import subprocess  # nosec
from pathlib import Path

import pytest
from pyuvdata.telescopes import known_telescope_location

from pyuvsim import cli
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH

try:
    subprocess.check_output(["convert", "--version"])  # nosec
    img_mgk_installed = True
except FileNotFoundError:
    img_mgk_installed = False


def test_run_pyuvsim_errors():
    pytest.importorskip("mpi4py")
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


@pytest.mark.parametrize(["verbosity", "plot"], [(None, True), (1, False), (2, False)])
def test_text_to_catalog_basic(verbosity, plot, goto_tempdir):
    """This just tests that it runs, not that it's correct."""
    try:
        import matplotlib  # noqa
    except ImportError:
        with pytest.raises(
            ImportError,
            match="matplotlib must be installed to use create text catalogs.",
        ):
            cli.text_to_catalog(["-t", "R"])

    if not img_mgk_installed:
        with pytest.raises(
            RuntimeError, match="ImageMagick must installed to create text catalogs"
        ):
            cli.text_to_catalog(["-t", "R"])

    pytest.importorskip("matplotlib")
    mwa_location = known_telescope_location("mwa")
    command = [
        "text_to_catalog",
        "-t",
        "R",
        "-n",
        str(20),
        "--lat",
        str(mwa_location.lat.deg),
        "--lon",
        str(mwa_location.lon.deg),
    ]
    if verbosity is not None:
        command.append(f"-{'v' * verbosity}")
    if plot:
        command.append("--plot")

    output = subprocess.check_output(command)  # nosec
    if verbosity is None:
        assert output.decode("utf-8").startswith("saved catalog file to")
    elif verbosity == 1:
        assert output.decode("utf-8").startswith("generating image file")
        assert "saved catalog file to" in output.decode("utf-8")
    else:
        assert output.decode("utf-8").startswith(
            "['convert', '-background', 'black', '-fill', 'white'"
        )
        assert "generating image file" in output.decode("utf-8")
        assert "saved catalog file to" in output.decode("utf-8")

    path = Path(goto_tempdir)
    bmp_file = path / "R.bmp"
    skyh5_file = path / "R.skyh5"

    assert bmp_file.exists()
    assert skyh5_file.exists()
    if plot:
        plotfile = path / "R.png"
        assert plotfile.exists()


def test_plot_csv_antpos_basic(goto_tempdir):
    """This just tests that it runs, not that it's correct."""
    layout_file = Path(SIM_DATA_PATH) / "test_config" / "layout_hex37_14.6m.csv"
    try:
        import matplotlib  # noqa
    except ImportError:
        with pytest.raises(
            ImportError,
            match="matplotlib must be installed to use plot antenna positions.",
        ):
            cli.plot_csv_antpos([layout_file])

    pytest.importorskip("matplotlib")
    subprocess.check_output(["plot_csv_antpos", str(layout_file)])  # nosec
    plotfile = Path(goto_tempdir) / (layout_file.stem + ".png")
    assert plotfile.exists()
