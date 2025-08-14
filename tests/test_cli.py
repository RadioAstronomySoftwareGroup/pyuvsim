# Copyright (c) 2025 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License
import subprocess  # nosec
from pathlib import Path

import pytest
from pyuvdata.telescopes import known_telescope_location
from pyuvdata.testing import check_warnings

from pyuvsim import cli
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH

try:
    subprocess.check_output(["convert", "--version"])  # nosec
    img_mgk_installed = True
except (FileNotFoundError, subprocess.CalledProcessError):
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


@pytest.mark.parametrize(
    ["verbosity", "plot", "use_old"],
    [(None, True, False), (1, False, False), (2, True, True)],
)
def test_text_to_catalog_basic(verbosity, plot, use_old, goto_tempdir, capsys):
    """This just tests that it runs, not that it's correct."""
    try:
        import matplotlib  # noqa
    except ImportError:
        with pytest.raises(
            ImportError,
            match="matplotlib must be installed to use create text catalogs.",
        ):
            cli.text_to_catalog(["-t", "R"])

    pytest.importorskip("matplotlib")

    if not img_mgk_installed:
        with pytest.raises(
            RuntimeError, match="ImageMagick must installed to create text catalogs"
        ):
            cli.text_to_catalog(["-t", "R"])
    else:
        mwa_location = known_telescope_location("mwa")
        if use_old:
            command = []
        else:
            command = ["text_to_catalog"]
        command += [
            "-t",
            "R",
            "-n",
            str(20),
            "--jd",
            str(2460000),
            "--lat",
            str(mwa_location.lat.deg),
            "--lon",
            str(mwa_location.lon.deg),
        ]
        if verbosity is not None and not use_old:
            command.append(f"-{'v' * verbosity}")
        if plot:
            command.append("--plot")

        if use_old:
            with check_warnings(
                DeprecationWarning,
                match="This script is deprecated in favor of text_to_catalog. "
                "It will be removed in version 1.6",
            ):
                cli.im_to_catalog(command)
            output = capsys.readouterr().out
        else:
            output = subprocess.check_output(command).decode("utf-8")  # nosec
        if verbosity is None:
            assert output.startswith("saved catalog file to")
        elif verbosity == 1:
            assert output.startswith("generating image file")
            assert "saved catalog file to" in output
        else:
            assert output.startswith(
                "['convert', '-background', 'black', '-fill', 'white'"
            )
            assert "generating image file" in output
            assert "saved catalog file to" in output

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
