# Copyright (c) 2025 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License
import importlib
import os
import subprocess  # nosec
from pathlib import Path

import numpy as np
import pytest
from pyradiosky import SkyModel

# from pyuvdata.telescopes import known_telescope_location
from pyuvdata.testing import check_warnings

from pyuvsim import cli
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH

from . import get_vizier_mirror

try:
    subprocess.check_output(["convert", "--version"])  # nosec
    img_mgk_installed = True
except (FileNotFoundError, subprocess.CalledProcessError):
    img_mgk_installed = False

has_astroquery = importlib.util.find_spec("astroquery") is not None


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
        if use_old:
            command = []
        else:
            command = ["text_to_catalog"]
        # run command with --test as not testing image generation
        command += [
            "-t",
            "R",
            "-n",
            str(10),
            "--jd",
            str(2460000),
            "--lat",
            "-26.70331941",
            "--lon",
            "116.6708152",
            "--test",
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
        # just use the stored R.bmp file as we no longer test the image generation
        bmp_file = Path(os.path.join(SIM_DATA_PATH, "test_catalogs", "R.bmp"))
        skyh5_file = path / "R.skyh5"

        assert bmp_file.exists()
        assert skyh5_file.exists()

        if plot:
            plotfile = path / "R.png"
            assert plotfile.exists()

        # this next part could run only once but maybe better to check old and new?

        # load current catalog as SkyModel
        current_catalog = SkyModel()
        current_catalog.read_skyh5("R.skyh5")

        # load reference catalog as SkyModel
        reference_catalog = SkyModel()
        reference_catalog_path = os.path.join(SIM_DATA_PATH, "test_catalogs", "R.txt")
        reference_catalog.read(reference_catalog_path)

        # get lons and lats from current catalog
        # offset lons by 30 degrees as that is how reference catalog was made
        # round to 6 decimals as I believe other catalog was rounded
        lons = np.round(current_catalog.skycoord.data.lon.value + 30, 6)
        lats = np.round(current_catalog.skycoord.data.lat.value, 6)

        # get lons and lats from reference catalog
        old_lats = reference_catalog.skycoord.data.lat.value
        old_lons = reference_catalog.skycoord.data.lon.value

        # APPARENTLY ON SOME OS VERSIONS THE PRODUCED IMAGE HAS 19 POINTS!!!
        # SPECIAL CASING AGAINST THOSE VERSIONS FOR NOW
        # if current_catalog.Ncomponents == reference_catalog.Ncomponents:
        # check agreement of all of the source locations
        assert (lons == old_lons).all()
        assert (lats == old_lats).all()


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


# very limited set of small files, more comprehensive tests exist in test_run_ref.py
@pytest.mark.parametrize(
    ("args", "duplicate"), [(["--row_limit", "20", "gleam"], True), (["fail"], False)]
)
def test_download_data_files_small(args, duplicate):

    if "gleam" in args:
        pytest.importorskip("astroquery")
        vizier_mirror = get_vizier_mirror()
        if vizier_mirror is not None:
            args = args + ["--vizier_mirror", vizier_mirror]

    # check for duplicate download of gleam and bad file input to see warnings
    if duplicate:
        # run it twice to generate the warning. First time with clearing
        cli.download_data_files(args + ["--clear"])

        if vizier_mirror is None:
            # make sure we have coverage for supplying the vizier mirror
            args += ["--vizier_mirror", "vizier.u-strasbg.fr"]
        with pytest.warns(UserWarning, match=r"astropy cached url for .*"):
            cli.download_data_files(args)
    elif args[0] == "fail":
        with pytest.warns(UserWarning, match=r"file 'fail' not found! .*"):
            cli.download_data_files(args)
    else:
        cli.download_data_files(args)


@pytest.mark.skipif(has_astroquery, reason="Test only when astroquery not installed.")
def test_download_gleam_error_noastroquery():
    with pytest.raises(
        ImportError, match="astroquery must be installed to download GLEAM."
    ):
        cli.download_data_files(["--row_limit", "75001", "gleam"])
