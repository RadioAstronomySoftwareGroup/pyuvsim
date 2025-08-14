# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License
"""
Functions for setting up simulations based on input config files.

This provides a configuration interface of human-readable yaml and csv files to specify
all the required simulation parameters. This configuration interface is used by other
simulators as well.

This module contains methods to create configuration files from :class:`pyuvdata.UVData`
objects and empty :class:`pyuvdata.UVData` objects from configuration files.
"""

import ast
import contextlib
import copy
import logging
import os
import shutil
import warnings

import astropy.units as units
import numpy as np
import numpy.typing as npt
import yaml
from astropy.coordinates import (
    ICRS,
    AltAz,
    Angle,
    EarthLocation,
    Latitude,
    Longitude,
    SkyCoord,
)
from astropy.time import Time
from astropy.utils.data import cache_contents, is_url_in_cache
from pyradiosky import SkyModel
from pyuvdata import (
    AiryBeam,
    GaussianBeam,
    ShortDipoleBeam,
    Telescope,
    UniformBeam,
    UVBeam,
    UVData,
    utils as uvutils,
)
from pyuvdata.analytic_beam import AnalyticBeam

from pyuvsim.data import DATA_PATH as SIM_DATA_PATH

from .telescope import BeamList

try:
    import astropy_healpix
except ImportError:
    astropy_healpix = None
try:
    import analytic_diffuse
except ImportError:
    analytic_diffuse = None
try:
    from . import mpi
    from .mpi import get_rank
except ImportError:

    def get_rank():
        """Mock to prevent import errors."""
        return 0

    mpi = None


from .utils import check_file_exists_and_increment

__all__ = [
    "create_mock_catalog",
    "SkyModelData",
    "initialize_catalog_from_params",
    "initialize_uvdata_from_params",
    "initialize_uvdata_from_keywords",
    "uvdata_to_telescope_config",
    "uvdata_to_config_file",
]

logger = logging.getLogger(__name__)

# this dict can go away in version 1.5 when we require the !AnalyticBeam
# tag in yaml files to specify analytic beams
analytic_beams = {
    "airy": AiryBeam,
    "gaussian": GaussianBeam,
    "short_dipole": ShortDipoleBeam,
    "uniform": UniformBeam,
}


def _parse_layout_csv(layout_csv):
    """
    Interpret the layout csv file.

    Columns in the file provide, in order from left to right, the antenna name, antenna number,
    a beam ID number, and the antenna positions relative to the telescope location in
    east, north, up (ENU) in meters.
    See https://pyuvsim.readthedocs.io/en/latest/parameter_files.html for more details.

    Parameters
    ----------
    layout_csv : str
        Filename of a layout csv.

    """
    with open(layout_csv) as fhandle:
        header = fhandle.readline()

    header = [h.strip() for h in header.split()]
    str_format_code = "U"

    # get data types for each column
    dtypes = {
        "name": str_format_code + "10",
        "number": "i4",
        "beamid": "i4",
        "e": "f8",
        "n": "f8",
        "u": "f8",
    }

    # check columns in file
    lower_header = [col.lower() for col in header]
    columns = ["name", "number", "beamid", "e", "n", "u"]
    col_exist = [col for col in columns if col in lower_header]

    dt = np.rec.format_parser([dtypes[col] for col in col_exist], col_exist, header)

    return np.genfromtxt(layout_csv, autostrip=True, skip_header=1, dtype=dt.dtype)


def _write_layout_csv(
    filepath, antpos_enu, antenna_names, antenna_numbers, beam_ids=None
):
    """
    Write a layout csv file.

    Columns in the file provide, in order from left to right, the antenna name, antenna number,
    a beam ID number, and the antenna positions relative to the telescope location in
    east, north, up (ENU) in meters.
    See https://pyuvsim.readthedocs.io/en/latest/parameter_files.html for more details.

    Parameters
    ----------
    filepath : str
        Filename to save the file to.
    antpos_enu : array_like of float
        East, North, Up positions of antennas in meters relative to the telescope location.
    antenna_names : array_like of str
        Names of antennas.
    antenna_numbers : array_like of int
        Antenna numbers.
    beam_ids : array_like of int
        Beam ID numbers to associate each antenna with a beam from a BeamList object.
        Defaults to all zeros (all antennas have the same beam) if nothing is passed in.

    """
    col_width = max(len(name) for name in antenna_names)
    header = ("{:" + str(col_width) + "} {:8} {:8} {:10} {:10} {:10}\n").format(
        "Name", "Number", "BeamID", "E", "N", "U"
    )
    if beam_ids is None:
        beam_ids = np.zeros(len(antenna_names), dtype=int)
    with open(filepath, "w") as lfile:
        lfile.write(header + "\n")
        for i, (e, n, u) in enumerate(antpos_enu):
            beam_id = beam_ids[i]
            name = antenna_names[i]
            num = antenna_numbers[i]
            line = ("{:{}} {:8d} {:8d} {:10.4f} {:10.4f} {:10.4f}\n").format(
                name, col_width, num, beam_id, e, n, u
            )
            lfile.write(line)


def _config_str_to_dict(config_str):
    """
    Read a yaml file into a dictionary, also add the file path info to the dictionary.

    Parameters
    ----------
    config_str : str
        Filename of a configuration yaml file to read.

    """
    with open(config_str) as pfile:
        param_dict = yaml.safe_load(pfile)

    config_str = os.path.abspath(config_str)
    param_dict["config_path"] = os.path.dirname(config_str)
    param_dict["obs_param_file"] = os.path.basename(config_str)

    return param_dict


def _create_catalog_diffuse(
    map_nside, diffuse_model, diffuse_params, time, array_location
):
    """
    Make a :class:`pyradiosky.SkyModel` object from an analog diffuse map function.

    Only used for internal testing, should not be called by users.

    Parameters
    ----------
    map_nside : int
        HEALPix nside of desired map.
    diffuse_model : str
        Name of diffuse model.
    diffuse_params : dict
        Dictionary of parameters corresponding to the `diffuse_model`.
    time : :class:`astropy.time.Time` object
        Time to use in creating map (to assign RA/Decs based on Alt/Az).
    array_location : :class:`astropy.coordinates.EarthLocation` or :class:`lunarsky.MoonLocation`
        Location to use in creating map (to assign RA/Decs based on Alt/Az).

    Returns
    -------
    catalog : class:`pyradiosky.SkyModel`
        SkyModel object containing the diffuse map.
    icrs_coord : :class:`astropy.coordinates.SkyCoord`
        Astropy SkyCoord object containing the coordinates for the pixels in ICRS.
    alts :  array_like of float
        Altitudes of the pixels in radians.
    azs :  array_like of float
        Azimuths of the pixels in radians.
    fluxes : array_like of float
        Brightness temperature of the pixels in K.

    """
    # Make the map, calculate AltAz positions of pixels, evaluate the model.
    npix = 12 * map_nside**2
    pixinds = np.arange(npix)

    hpix = astropy_healpix.HEALPix(nside=map_nside, frame=ICRS())
    icrs_coord = hpix.healpix_to_skycoord(pixinds)

    use_altaz = True
    if not isinstance(array_location, EarthLocation):
        with contextlib.suppress(ImportError):
            from lunarsky import LunarTopo, MoonLocation, SkyCoord as LunarSkyCoord

            if isinstance(array_location, MoonLocation):
                localframe = LunarTopo(location=array_location, obstime=time)
                icrs_coord = LunarSkyCoord(icrs_coord)
                use_altaz = False

    if use_altaz:
        localframe = AltAz(location=array_location, obstime=time)

    source_coord = icrs_coord.transform_to(localframe)

    alts = source_coord.alt.rad
    azs = source_coord.az.rad

    modfunc = analytic_diffuse.get_model(diffuse_model)
    fluxes = modfunc(source_coord.az.rad, source_coord.zen.rad, **diffuse_params)

    stokes = np.zeros((4, 1, npix))
    stokes[0, :] = fluxes

    stokes *= units.K

    catalog = SkyModel(
        stokes=stokes,
        nside=map_nside,
        hpx_inds=pixinds,
        spectral_type="flat",
        frame="icrs",
    )

    return catalog, icrs_coord, alts, azs, fluxes


def _create_catalog_discrete(Nsrcs, alts, azs, fluxes, time, array_location):
    """
    Make a :class:`pyradiosky.SkyModel` object from a set of point sources.

    Only used for internal testing, should not be called by users.

    Parameters
    ----------
    Nsrcs : int
        Number of sources
    alts :  array_like of float
        Altitudes or Latitudes (depending on array_location) of the sources in radians.
    azs :  array_like of float
        Azimuths or Longitues (depending on array_location) of the sources in radians.
    fluxes : array_like of float
        Brightnesses of the sources in Jy.
    time : :class:`astropy.time.Time` object
        Time to use in defining the positions.
    array_location : :class:`astropy.coordinates.EarthLocation` or :class:`lunarsky.MoonLocation`
        Location to use in defining the positions.

    Returns
    -------
    catalog : class:`pyradiosky.SkyModel`
        SkyModel object containing the diffuse map.
    icrs_coord : :class:`astropy.coordinates.SkyCoord`
        Astropy SkyCoord object containing the coordinates for the sources in ICRS.

    """
    localframe = AltAz
    coord_class = SkyCoord
    if not isinstance(array_location, EarthLocation):
        with contextlib.suppress(ImportError):
            from lunarsky import LunarTopo, MoonLocation, SkyCoord as LunarSkyCoord

            if isinstance(array_location, MoonLocation):
                localframe = LunarTopo
                coord_class = LunarSkyCoord

    source_coord = coord_class(
        alt=Angle(alts, unit=units.deg),
        az=Angle(azs, unit=units.deg),
        obstime=time,
        frame=localframe,
        location=array_location,
    )
    icrs_coord = source_coord.transform_to("icrs")

    names = np.array(["src" + str(si) for si in range(Nsrcs)])
    stokes = np.zeros((4, 1, Nsrcs))
    stokes[0, :] = fluxes

    stokes *= units.Jy

    catalog = SkyModel(
        name=names, skycoord=icrs_coord, stokes=stokes, spectral_type="flat"
    )

    return catalog, icrs_coord


def create_mock_catalog(
    time,
    arrangement="zenith",
    array_location=None,
    Nsrcs=None,
    alt=None,
    save=False,
    min_alt=None,
    rseed=None,
    return_data=False,
    diffuse_model=None,
    diffuse_params=None,
    map_nside=None,
):
    """
    Create a mock catalog.

    SkyModels are defined in an AltAz frame at the given time, then returned in
    ICRS ra/dec coordinates. They always have a 'flat' spectral type.

    Parameters
    ----------
    time: float or :class:`astropy.time.Time` object
        Julian date, if a float it is interpreted as a JD.
    arrangement: str
        Point source pattern (default = 1 source at zenith).
        Accepted arrangements:

        * `triangle`:  Three point sources forming a triangle around the zenith
        * `cross`: An asymmetric cross
        * `zenith`: Some number of sources placed at the zenith.
        * `off-zenith`:  A single source off zenith
        * `long-line`:  Horizon to horizon line of point sources
        * `hera_text`:  Spell out HERA around the zenith
        * `random`:  Randomly distributed point sources near zenith
        * `diffuse`: Choice of solved models in the analytic_diffuse module.
    Nsrcs: int
        Number of sources to make (used for zenith, off-zenith, long-line, and random arrangements).
    array_location: :class:`astropy.coordinates.EarthLocation`
        Location of observer. Source positions are defined with respect to a particular zenith.
        Can also be a :class:`lunarsky.MoonLocation` object to make catalogs for
        observers on the Moon.
        Default = HERA site
    alt: float
        For off-zenith and triangle arrangements, altitude to place sources.
        In degrees.
    min_alt: float
        For random and long-line arrangements, minimum altitude at which to place sources.
        In degrees.
    save: bool
        Save mock catalog as npz file.
    rseed: int
        If using the random configuration, pass in a RandomState seed.
    return_data: bool
        If True, return a :class:`~.SkyModelData` object instead of
        :class:`pyradiosky.SkyModel`.
    diffuse_model: str
        If arrangement is 'diffuse', name of the diffuse model to generate.
        See documentation in `analytic_diffuse.models` for details.
    diffuse_params: dict
        If arrangement is 'diffuse', extra parameters accepted by the chosen model.
        See documentation in `analytic_diffuse.models` for details.
    map_nside: int
        If arrangement is 'diffuse', nside parameter for the healpix map to make.

    Returns
    -------
    :class:`pyradiosky.SkyModel` or :class:`~.SkyModelData`
        The catalog, as either a SkyModel or a SkyModelData (if `return_data` is True)
    dict
        A dictionary of keywords used to define the catalog.

    Notes
    -----
    Generating diffuse models requires the `analytic_diffuse` and `astropy_healpix` modules.

    """
    if isinstance(time, float):
        time = Time(time, scale="utc", format="jd")

    if array_location is None:
        array_location = EarthLocation(
            lat="-30d43m17.5s", lon="21d25m41.9s", height=1073.0
        )

    if arrangement not in [
        "off-zenith",
        "zenith",
        "cross",
        "triangle",
        "long-line",
        "hera_text",
        "random",
        "diffuse",
    ]:
        raise KeyError("Invalid mock catalog arrangement: " + str(arrangement))

    mock_keywords = {
        "time": time.jd,
        "arrangement": arrangement,
        "array_location": repr(
            (
                float(array_location.lat.deg),
                float(array_location.lon.deg),
                float(array_location.height.value),
            )
        ),
    }

    mock_keywords["world"] = "earth"
    if not isinstance(array_location, EarthLocation):
        with contextlib.suppress(ImportError):
            from lunarsky import MoonLocation

            if isinstance(array_location, MoonLocation):
                mock_keywords["world"] = "moon"

    if arrangement == "off-zenith":
        if alt is None:
            alt = 85.0  # Degrees
        mock_keywords["alt"] = alt
        Nsrcs = 1
        alts = [alt]
        azs = [90.0]  # 0 = North pole, 90. = East pole
        fluxes = [1.0]

    if arrangement == "triangle":
        Nsrcs = 3
        if alt is None:
            alt = 87.0  # Degrees
        mock_keywords["alt"] = alt
        alts = [alt, alt, alt]
        azs = [0.0, 120.0, 240.0]
        fluxes = [1.0, 1.0, 1.0]

    if arrangement == "cross":
        Nsrcs = 4
        alts = [88.0, 90.0, 86.0, 82.0]
        azs = [270.0, 0.0, 90.0, 135.0]
        fluxes = [5.0, 4.0, 1.0, 2.0]

    if arrangement == "zenith":
        if Nsrcs is None:
            Nsrcs = 1
        mock_keywords["Nsrcs"] = Nsrcs
        alts = np.ones(Nsrcs) * 90.0
        azs = np.zeros(Nsrcs, dtype=float)
        fluxes = np.ones(Nsrcs) * 1 / Nsrcs
        # Divide total Stokes I intensity among all sources
        # Test file has Stokes I = 1 Jy

    if arrangement == "random":
        if Nsrcs is None:
            Nsrcs = 1
        if min_alt is None:
            min_alt = 30  # Degrees
        mock_keywords["Nsrcs"] = Nsrcs
        np.random.seed(seed=rseed)
        mock_keywords["rseed"] = np.random.get_state()[1][0]

        # Necessary to get uniform distribution per solid angle.
        min_alt_rad = np.radians(min_alt)
        rv = np.random.uniform(-1, np.cos(min_alt_rad + np.pi / 2), Nsrcs)
        alts = np.degrees(np.arccos(rv) - np.pi / 2)
        azs = np.degrees(np.random.uniform(0, 2 * np.pi, Nsrcs))
        fluxes = np.ones(Nsrcs, dtype=float)

    if arrangement == "long-line":
        if Nsrcs is None:
            Nsrcs = 10
        if min_alt is None:
            min_alt = 5
        mock_keywords["Nsrcs"] = Nsrcs
        mock_keywords["alt"] = alt
        fluxes = np.ones(Nsrcs, dtype=float)
        if Nsrcs % 2 == 0:
            length = 180 - min_alt * 2
            spacing = length / (Nsrcs - 1)
            max_alt = 90.0 - spacing / 2
            alts = np.linspace(min_alt, max_alt, Nsrcs // 2)
            alts = np.append(alts, np.flip(alts, axis=0))
            azs = np.append(
                np.zeros(Nsrcs // 2, dtype=float) + 180.0,
                np.zeros(Nsrcs // 2, dtype=float),
            )
        else:
            alts = np.linspace(min_alt, 90, (Nsrcs + 1) // 2)
            alts = np.append(alts, np.flip(alts[1:], axis=0))
            azs = np.append(
                np.zeros((Nsrcs + 1) // 2, dtype=float) + 180.0,
                np.zeros((Nsrcs - 1) // 2, dtype=float),
            )

    if arrangement == "hera_text":
        # fmt: off
        azs = np.array([-254.055, -248.199, -236.310, -225.000, -206.565,
                        -153.435, -123.690, -111.801, -105.945, -261.870,
                        -258.690, -251.565, -135.000, -116.565, -101.310,
                        -98.130, 90.000, 90.000, 90.000, 90.000, 90.000,
                        -90.000, -90.000, -90.000, -90.000, -90.000,
                        -90.000, 81.870, 78.690, 71.565, -45.000, -71.565,
                        -78.690, -81.870, 74.055, 68.199, 56.310, 45.000,
                        26.565, -26.565, -45.000, -56.310, -71.565])

        zas = np.array([7.280, 5.385, 3.606, 2.828, 2.236, 2.236, 3.606,
                        5.385, 7.280, 7.071, 5.099, 3.162, 1.414, 2.236,
                        5.099, 7.071, 7.000, 6.000, 5.000, 3.000, 2.000,
                        1.000, 2.000, 3.000, 5.000, 6.000, 7.000, 7.071,
                        5.099, 3.162, 1.414, 3.162, 5.099, 7.071, 7.280,
                        5.385, 3.606, 2.828, 2.236, 2.236, 2.828, 3.606, 6.325])
        # fmt: on

        alts = 90.0 - zas
        Nsrcs = alts.size
        fluxes = np.ones_like(alts)

    if arrangement == "diffuse":
        if analytic_diffuse is None or astropy_healpix is None:
            raise ValueError(
                "analytic_diffuse and astropy_healpix "
                "modules are required to use diffuse mock catalog."
            )
        if diffuse_model is None:
            raise ValueError(
                "Diffuse arrangement selected, but diffuse model not chosen."
            )
        if diffuse_params is None:
            diffuse_params = {}
        if map_nside is None:
            warnings.warn("No nside chosen. Defaulting to 32")
            map_nside = 32

        mock_keywords["diffuse_model"] = diffuse_model
        mock_keywords["map_nside"] = map_nside
        mock_keywords["diffuse_params"] = repr(diffuse_params)

        catalog, icrs_coord, alts, azs, fluxes = _create_catalog_diffuse(
            map_nside, diffuse_model, diffuse_params, time, array_location
        )

    else:
        catalog, icrs_coord = _create_catalog_discrete(
            Nsrcs, alts, azs, fluxes, time, array_location
        )

    if return_data:
        catalog = SkyModelData(catalog)

    if get_rank() == 0 and save:
        np.savez(
            "mock_catalog_" + arrangement,
            ra=icrs_coord.ra.rad,
            dec=icrs_coord.dec.rad,
            alts=alts,
            azs=azs,
            fluxes=fluxes,
        )

    return catalog, mock_keywords


class SkyModelData:
    """
    Carries immutable SkyModel data in simple ndarrays.

    This is to facilitate sharing SkyModel objects in MPI, without
    excessive copying and memory bloat.

    When used with MPI, this can be initialized simultaneously on all
    processes such that sky_in is provided only on the root process.
    The data can then be shared across all processes by running the `share`
    method.

    Parameters
    ----------
    sky_in: :class:`pyradiosky.SkyModel`
        A valid SkyModel object.
    filename : str or list of str, optional
        The filename (or other string identifier) of the input catalog. This overrides
        the filename set on the sky_in object (if it has one). If not set, this defaults
        to the filename set on the sky_in object (if it has one)

    """

    filename = None
    pyradiosky_version_str = None
    Ncomponents = None
    component_type = None
    spectral_type = None
    ra = None
    dec = None
    name = None
    Nfreqs = None
    stokes_I = None  # noqa should be all lower case
    stokes_Q = None  # noqa should be all lower case
    stokes_U = None  # noqa should be all lower case
    stokes_V = None  # noqa should be all lower case
    freq_array = None
    freq_edge_array = None
    reference_frequency = None
    spectral_index = None
    polarized = None
    nside = None
    hpx_inds = None
    flux_unit = None

    put_in_shared = [
        "stokes_I",
        "stokes_Q",
        "stokes_U",
        "stokes_V",
        "polarized",
        "ra",
        "dec",
        "reference_frequency",
        "spectral_index",
        "hpx_inds",
    ]

    def __init__(self, sky_in=None, filename=None):
        # Collect relevant attributes.
        if sky_in is not None:
            if sky_in.name is not None:
                self.name = np.asarray(sky_in.name)
            else:
                self.name = None
            self.nside = sky_in.nside
            self.hpx_inds = sky_in.hpx_inds
            self.component_type = sky_in.component_type
            self.spectral_type = sky_in.spectral_type
            self.Ncomponents = sky_in.Ncomponents
            if sky_in.frame != "icrs":
                sky_in.transform_to(ICRS)
            sky_ra, sky_dec = sky_in.get_lon_lat()
            self.ra = sky_ra.deg
            self.dec = sky_dec.deg
            self.Nfreqs = sky_in.Nfreqs
            stokes_in = sky_in.stokes

            if isinstance(stokes_in, units.Quantity):
                self.flux_unit = stokes_in.unit.to_string()
                stokes_in = stokes_in.value

            self.stokes_I = stokes_in[0, ...]

            if sky_in._n_polarized > 0:
                self.polarized = sky_in._polarized
                Q, U, V = stokes_in[1:, :, sky_in._polarized]
                self.stokes_Q = Q
                self.stokes_U = U
                self.stokes_V = V

            if sky_in._freq_array.required:
                self.freq_array = sky_in.freq_array.to("Hz").value
            if sky_in.spectral_type == "subband" and sky_in.freq_edge_array is not None:
                self.freq_edge_array = sky_in.freq_edge_array.to("Hz").value

            if sky_in._reference_frequency.required:
                self.reference_frequency = sky_in.reference_frequency.to("Hz").value
            if sky_in._spectral_index.required:
                self.spectral_index = sky_in.spectral_index

            self.pyradiosky_version_str = sky_in.pyradiosky_version_str
        if filename is not None:
            if isinstance(filename, str):
                filename_use = [filename]
            else:
                filename_use = filename
            self.filename = filename_use
        elif sky_in is not None and sky_in.filename is not None:
            self.filename = sky_in.filename

    def subselect(self, inds):
        """
        Subselect, returning a new SkyModelData object.

        Parameters
        ----------
        inds: range or index array
            Indices to select along the component axis.

        Returns
        -------
        SkyModelData
            A new SkyModelData with Ncomp axes downselected.

        Notes
        -----
        If inds is a range object, this method will avoid copying data in numpy arrays,
        such that the returned SkyModelData object carries views into the current object's arrays.

        """
        new_sky = SkyModelData()

        new_sky.Ncomponents = len(inds)
        new_sky.nside = self.nside
        new_sky.component_type = self.component_type
        if self.name is not None:
            new_sky.name = self.name[inds]
        if isinstance(inds, range):
            new_sky.stokes_I = self.stokes_I[:, slice(inds.start, inds.stop, inds.step)]
        else:
            new_sky.stokes_I = self.stokes_I[:, inds]
        new_sky.ra = self.ra[inds]
        new_sky.dec = self.dec[inds]
        new_sky.Nfreqs = self.Nfreqs
        new_sky.spectral_type = self.spectral_type
        new_sky.flux_unit = self.flux_unit

        if self.reference_frequency is not None:
            new_sky.reference_frequency = self.reference_frequency[inds]
        if self.spectral_index is not None:
            new_sky.spectral_index = self.spectral_index[inds]
        if self.freq_array is not None:
            new_sky.freq_array = self.freq_array
        if self.freq_edge_array is not None:
            new_sky.freq_edge_array = self.freq_edge_array
        if self.hpx_inds is not None:
            new_sky.hpx_inds = self.hpx_inds[inds]

        if self.polarized is not None:
            sub_inds = np.isin(self.polarized, inds)
            new_sky.stokes_Q = self.stokes_Q[..., sub_inds]
            new_sky.stokes_U = self.stokes_U[..., sub_inds]
            new_sky.stokes_V = self.stokes_V[..., sub_inds]
            new_sky.polarized = np.where(np.isin(inds, self.polarized))[0]

        return new_sky

    def share(self, root=0):
        """
        Share across MPI processes (requires mpi4py to use).

        All attributes are put in shared memory.

        Parameters
        ----------
        root : int
            Root rank on COMM_WORLD, from which data will be broadcast.

        """
        if mpi is None:
            raise ImportError(
                "You need mpi4py to use this method. "
                "Install it by running pip install pyuvsim[sim] "
                "or pip install pyuvsim[all] if you also want the "
                "line_profiler installed."
            )
        mpi.start_mpi()
        self.Ncomponents = mpi.world_comm.bcast(self.Ncomponents, root=root)

        # Get list of attributes that are set.
        isset = None
        if mpi.rank == root:
            isset = [key for key, value in self.__dict__.items() if value is not None]
        isset = mpi.world_comm.bcast(isset, root=root)

        for key in isset:
            attr = getattr(self, key)
            if key in self.put_in_shared:
                val = mpi.shared_mem_bcast(attr, root=root)
            else:
                val = mpi.world_comm.bcast(attr, root=root)
            if val is not None:
                setattr(self, key, val)

        mpi.world_comm.Barrier()

    def get_skymodel(self, inds=None):
        """
        Initialize :class:`pyradiosky.SkyModel` from current settings.

        Parameters
        ----------
        inds : range or index array
            Indices to select along the component axis.

        """
        if inds is not None:
            obj = self.subselect(inds)
            return obj.get_skymodel()

        stokes_use = np.zeros((4, self.Nfreqs, self.Ncomponents), dtype=float)

        stokes_use[0, ...] = self.stokes_I

        if self.polarized is not None:
            stokes_use[1, :, self.polarized] = self.stokes_Q.T
            stokes_use[2, :, self.polarized] = self.stokes_U.T
            stokes_use[3, :, self.polarized] = self.stokes_V.T

        if self.flux_unit is not None:
            stokes_use *= units.Unit(self.flux_unit)

        other = {}
        if self.spectral_type in ["full", "subband"]:
            other["freq_array"] = self.freq_array * units.Hz
        if self.freq_edge_array is not None:
            other["freq_edge_array"] = self.freq_edge_array * units.Hz
        if self.spectral_type == "spectral_index":
            other["reference_frequency"] = self.reference_frequency * units.Hz
            other["spectral_index"] = self.spectral_index
        if self.component_type == "healpix":
            other["nside"] = self.nside
            other["hpx_inds"] = self.hpx_inds
            other["frame"] = "icrs"
        else:
            other["name"] = self.name
            other["skycoord"] = SkyCoord(
                ra=Longitude(self.ra, unit="deg"),
                dec=Latitude(self.dec, unit="deg"),
                frame="icrs",
            )
        other["filename"] = self.filename

        return SkyModel(stokes=stokes_use, spectral_type=self.spectral_type, **other)


def _sky_select_calc_rise_set(sky, source_params, telescope_lat_deg=None):
    """
    Apply flux and non-rising source cuts, calculate rise and set lsts.

    Parameters
    ----------
    sky : pyradiosky.Skymodel
        SkyModel object to apply cuts to.
    source_params : dict
        Dict specifying flux cut and horizon buffer parameters.
    telescope_lat_deg : float
        Latitude of telescope in degrees, used for horizon calculations.

    """
    # do other selections before non-rising to avoid warnings about nans and
    # negatives that are filtered in this step in the cut_nonrising step
    select_params = {}
    if "min_flux" in source_params:
        min_brightness = float(source_params["min_flux"]) * units.Jy
        select_params["min_brightness"] = min_brightness
    if "max_flux" in source_params:
        max_brightness = float(source_params["max_flux"]) * units.Jy
        select_params["max_brightness"] = max_brightness
    for par in ["non_nan", "non_negative"]:
        if par in source_params:
            select_params[par] = source_params[par]

    sky.select(**select_params)

    if telescope_lat_deg is not None:
        telescope_lat = Latitude(telescope_lat_deg * units.deg)
        sky.cut_nonrising(telescope_lat)
        calc_kwargs = {}
        if "horizon_buffer" in source_params:
            calc_kwargs["horizon_buffer"] = float(source_params["horizon_buffer"])
        sky.calculate_rise_set_lsts(telescope_lat, **calc_kwargs)

    return sky


def initialize_catalog_from_params(
    obs_params, input_uv=None, filetype=None, return_catname=False
):
    """
    Make catalog from parameter file specifications.

    Default behavior is to do coarse horizon cuts.

    Parameters
    ----------
    obs_params: str or dict
        Either an obsparam file name or a dictionary of parameters.
    input_uv: :class:`pyuvdata.UVData`
        Used to set location for mock catalogs and for horizon cuts.
    filetype : str
        One of ['skyh5', 'gleam', 'vot', 'text', 'fhd'] or None.
        If None, the code attempts to guess what the file type is.
    return_catname: bool
        Return the catalog name. Default is False.

    Returns
    -------
    sky: :class:`pyradiosky.SkyModel`
        Source catalog filled with data.
    source_list_name: str
        Catalog identifier for metadata. Only returned if return_catname is True.

    """
    if input_uv is not None and not isinstance(input_uv, UVData):
        raise TypeError("input_uv must be UVData object")

    if isinstance(obs_params, str):
        with open(obs_params) as pfile:
            param_dict = yaml.safe_load(pfile)

        param_dict["config_path"] = os.path.dirname(obs_params)
    else:
        param_dict = obs_params

    source_params = param_dict["sources"]
    if "catalog" in source_params:
        catalog = source_params["catalog"]
    else:
        raise KeyError("No catalog defined.")

    if catalog == "mock":
        mock_keywords = {"arrangement": source_params["mock_arrangement"]}
        extra_mock_kwds = [
            "time",
            "Nsrcs",
            "zen_ang",
            "save",
            "min_alt",
            "array_location",
            "diffuse_model",
            "diffuse_params",
            "map_nside",
        ]
        for k in extra_mock_kwds:
            if k in source_params:
                if k == "array_location":
                    # String -- lat, lon, alt in degrees
                    latlonalt = [float(s.strip()) for s in source_params[k].split(",")]
                    lat, lon, alt = latlonalt
                    mock_keywords[k] = EarthLocation.from_geodetic(lon, lat, alt)
                else:
                    mock_keywords[k] = source_params[k]

        # time, arrangement, array_location, save, Nsrcs, min_alt

        if "array_location" not in mock_keywords:
            if input_uv is not None:
                mock_keywords["array_location"] = input_uv.telescope.location
            else:
                warnings.warn(
                    "No array_location specified. Defaulting to the HERA site."
                )
        if "time" not in mock_keywords:
            if input_uv is not None:
                mock_keywords["time"] = input_uv.time_array[0]
                warnings.warn(
                    "No julian date given for mock catalog. Defaulting to first time step."
                )
            else:
                raise ValueError(
                    "input_uv must be supplied if using mock catalog without specified julian date"
                )

        time = mock_keywords.pop("time")

        sky, mock_keywords = create_mock_catalog(time, **mock_keywords)
        mock_keyvals = [str(key) + str(val) for key, val in mock_keywords.items()]
        source_list_name = "mock_" + "_".join(mock_keyvals)
    elif isinstance(catalog, str):
        # if catalog file is not found, first check relative to config_path,
        # then check astropy cache treating catalog input as url
        if not os.path.isfile(catalog):
            # create boolean to determine if cache should be checked
            check_cache = True
            if "config_path" in param_dict:
                relative_path = os.path.join(param_dict["config_path"], catalog)
                if os.path.isfile(relative_path):
                    catalog = relative_path
                    # if file is found no need to check cache
                    check_cache = False
            if check_cache and is_url_in_cache(catalog, pkgname="pyuvsim"):
                catalog = cache_contents("pyuvsim")[catalog]

        allowed_read_params = [
            "spectral_type",
            "table_name",
            "id_column",
            "ra_column",
            "dec_column",
            "lon_column",
            "lat_column",
            "frame",
            "flux_columns",
            "reference_frequency",
            "freq_array",
            "spectral_index_column",
        ]
        read_params = {}
        if filetype is not None:
            read_params["filetype"] = filetype
        elif "filetype" in source_params:
            read_params["filetype"] = source_params["filetype"]

        for param, value in source_params.items():
            if param in allowed_read_params:
                read_params[param] = value

        detect_gleam = filetype == "gleam" or (
            os.path.splitext(catalog)[1] == ".vot" and "gleam" in catalog.casefold()
        )

        if detect_gleam and "spectral_type" not in source_params:
            read_params["spectral_type"] = "subband"

        source_list_name = os.path.basename(catalog)
        # defer the acceptibility check until after selections are done.
        sky = SkyModel.from_file(catalog, **read_params, run_check_acceptability=False)

    if input_uv is not None:
        telescope_lat_deg = input_uv.telescope.location.lat.to("deg").value
    else:
        telescope_lat_deg = None
    # Do source selections, if any.
    sky = _sky_select_calc_rise_set(
        sky, source_params, telescope_lat_deg=telescope_lat_deg
    )
    # do the full check after selections to avoid warnings about nans and
    # negatives that are removed in the selection
    sky.check()

    # If the filename parameter is None (e.g. for mock skies)
    # add the source_list_name to the object so it can be put in the UVData history later
    if sky.filename is None:
        sky.filename = [source_list_name]

    if return_catname:
        return sky, source_list_name
    else:
        return sky


def _check_uvbeam_file(beam_file):
    if not os.path.exists(beam_file):
        testdata_beam_file = os.path.join(SIM_DATA_PATH, beam_file)
        if os.path.exists(testdata_beam_file):
            beam_file = testdata_beam_file
        else:
            raise ValueError(f"Unrecognized beam file or model: {beam_file}")
    return beam_file


def _construct_beam_list(
    beam_ids: npt.NDArray[int] | list[float] | tuple[float],
    telconfig: dict,
    freq_array: npt.NDArray[float],
    freq_range: (npt.NDArray[float] | list[float] | tuple[float] | None) = None,
):
    """
    Construct the BeamList from the telescope parameter dict.

    Parameters
    ----------
    beam_ids : list of int
        List of beam ID integers to include from the telconfig. All integers in
        this list must be included in the telconfig.
    telconfig : dict
        Telescope parameters dict, typically parsed from the telescope yaml file.
        See pyuvsim documentation for allowable keys.
        https://pyuvsim.readthedocs.io/en/latest/parameter_files.html#telescope-configuration
    freq_array : array-like of float
        Simulation frequency array.
    freq_range : array-like of float
        Range of frequencies to keep in the beam object, shape (2,). Should
        include all frequencies in the simulation with a buffer to support
        interpolation.

    Returns
    -------
    BeamList
        The constructed beam list object.

    """
    beam_list = []
    assert isinstance(beam_ids, np.ndarray | list | tuple), (
        "beam_ids must be a list, tuple or numpy array."
    )
    beam_ids = np.asarray(beam_ids, dtype=int)
    assert isinstance(telconfig, dict), "telconfig must be a dict."

    if freq_range is not None:
        assert isinstance(freq_range, np.ndarray | list | tuple), (
            "If passed, freq_range must be a list, tuple or numpy array"
        )
        freq_range = np.asarray(freq_range, dtype=float).squeeze()
        if freq_range.size != 2:
            raise ValueError("If passed, freq_range have 2 elements")

    # possible global shape options
    if "diameter" in telconfig or "sigma" in telconfig:
        raise ValueError(
            "Beam shape options diameter and sigma should be specified per beamID "
            "in the 'beam_paths' section not as globals. For examples see the "
            "parameter_files documentation."
        )

    select = telconfig.pop("select", {})
    freq_buffer = select.pop("freq_buffer", None)
    if freq_range is not None and "freq_range" not in select and freq_buffer is None:
        select["freq_range"] = freq_range
    if freq_buffer is not None and "freq_range" not in select:
        freq_arr_val = freq_array
        freq_range = (
            freq_arr_val.min() - freq_buffer,
            freq_arr_val.max() + freq_buffer,
        )
        select["freq_range"] = freq_range

    for beam_id in beam_ids:
        if beam_id not in telconfig["beam_paths"]:
            raise ValueError(f"beam_id {beam_id} is not listed in the telconfig.")
        beam_model = telconfig["beam_paths"][beam_id]

        if not isinstance(beam_model, dict | AnalyticBeam | UVBeam):
            raise ValueError(
                "Beam model is not properly specified in telescope config file. "
                f"It is specified as {beam_model}."
            )

        if not isinstance(beam_model, AnalyticBeam | UVBeam) and not (
            isinstance(beam_model, dict) and "filename" in beam_model
        ):
            warnings.warn(
                "Entries in 'beam_paths' should be specified using either the "
                "UVBeam or AnalyticBeam constructors or using a dict syntax for "
                "UVBeams. For examples see the parameter_files documentation. "
                "Specifying analytic beam without the AnalyticBeam constructors "
                "will cause an error in version 1.6",
                DeprecationWarning,
            )

        if isinstance(beam_model, AnalyticBeam | UVBeam):
            if len(select) > 0 and isinstance(beam_model, UVBeam):
                if "freq_range" in select:
                    freq_range = select.pop("freq_range")
                    freq_chans = np.nonzero(
                        (beam_model.freq_array >= freq_range[0])
                        & (beam_model.freq_array <= freq_range[1])
                    )[0]
                    select["freq_chans"] = freq_chans
                beam_model.select(**select)
            beam_list.append(beam_model)
        elif isinstance(beam_model, dict) and "filename" in beam_model:
            # this a UVBeam readable file with (possibly) some read kwargs
            beam_file = beam_model.pop("filename")
            read_kwargs = beam_model

            beam_file = _check_uvbeam_file(beam_file)
            uvb = UVBeam()

            if "freq_range" in select:
                read_kwargs["freq_range"] = select["freq_range"]

            uvb.read(beam_file, **read_kwargs)

            beam_list.append(uvb)
        else:
            if "type" in beam_model:
                beam_type = beam_model["type"]
            else:
                raise ValueError(
                    "Beam model must have either a 'filename' field for UVBeam "
                    "files or a 'type' field for analytic beams."
                )

            if beam_type not in analytic_beams:
                raise ValueError(f"Undefined beam model type: {beam_type}")

            this_beam_opts = {}
            if isinstance(beam_model, dict):
                for key in beam_model:
                    if key != "type":
                        this_beam_opts[key] = beam_model[key]

            # Gaussian beam requires either diameter or sigma
            # Airy beam requires diameter

            # Default to None for diameter and sigma.
            # Values in the "beam_paths" override globally-defined options.
            if beam_type == "gaussian":
                shape_opts = {"diameter": None, "sigma": None}
            elif beam_type == "airy":
                shape_opts = {"diameter": None}
            else:
                shape_opts = {}

            for opt in shape_opts:
                shape_opts[opt] = this_beam_opts.get(opt)

            ab = analytic_beams[beam_type](**shape_opts)

            beam_list.append(ab)

    bl_options = {}
    if "spline_interp_opts" in telconfig:
        bl_options["spline_interp_opts"] = telconfig["spline_interp_opts"]
    if "freq_interp_kind" in telconfig:
        bl_options["freq_interp_kind"] = telconfig["freq_interp_kind"]

    beam_list_obj = BeamList(beam_list, **bl_options)

    return beam_list_obj


def _get_tel_loc(tel_dict: dict):
    """Get the telescope location info out of a dict."""
    telescope_location_latlonalt = tel_dict["telescope_location"]
    if isinstance(telescope_location_latlonalt, str):
        telescope_location_latlonalt = ast.literal_eval(telescope_location_latlonalt)
    world = tel_dict.pop("world", None)
    # get the lunar ellipsoid. Default to None for earth or SPHERE for moon
    if world == "moon":
        ellipsoid = tel_dict.pop("ellipsoid", "SPHERE")
    else:
        ellipsoid = tel_dict.pop("ellipsoid", None)

    return telescope_location_latlonalt, world, ellipsoid


def _get_array_layout(tele_params, config_path):
    # get array layout
    if "array_layout" not in tele_params:
        raise KeyError("array_layout must be provided.")
    array_layout = tele_params.pop("array_layout")
    if not isinstance(array_layout, str | dict):
        raise ValueError(
            "array_layout must be a string or have options that parse as a dict."
        )
    if isinstance(array_layout, str):
        # Interpet as file path to layout csv file.
        layout_csv = array_layout
        # if array layout is a str, parse it as .csv filepath
        if isinstance(layout_csv, str):
            if not os.path.exists(layout_csv) and isinstance(config_path, str):
                layout_csv = os.path.join(config_path, layout_csv)
                if not os.path.exists(layout_csv):
                    raise ValueError(
                        f"layout_csv file {layout_csv} from yaml does not exist"
                    )
            ant_layout = _parse_layout_csv(layout_csv)
            east, north, up = ant_layout["e"], ant_layout["n"], ant_layout["u"]
            antnames = ant_layout["name"]
            antnums = np.array(ant_layout["number"])
            beam_ids = ant_layout["beamid"]
    elif isinstance(array_layout, dict):
        # Receiving antenna positions directly
        antnums = tele_params.pop("antenna_numbers")
        antnames = tele_params.pop("antenna_names")
        east, north, up = np.array([array_layout[an] for an in antnums]).T
        layout_csv = "user-fed dict"
        beam_ids = None

    return antnames, antnums, east, north, up, layout_csv, beam_ids


def parse_telescope_params(
    tele_params: dict,
    *,
    config_path: str = "",
    freq_array: npt.NDArray[float],
    freq_range: (npt.NDArray[float] | list[float] | tuple[float] | None) = None,
    return_beams: bool = True,
):
    """
    Parse the "telescope" section of obsparam.

    Parameters
    ----------
    tele_params : dict
        Telescope parameters
        See pyuvsim documentation for allowable keys.
        https://pyuvsim.readthedocs.io/en/latest/parameter_files.html#telescope-configuration
    config_path : str
        Path to directory holding configuration and layout files.
    freq_array : array-like of float
        Simulation frequency array.
    freq_range : float
        If given, select frequencies on reading the beam.
    return_beams : bool
        Option to return the beam_list and beam_dict.

    Returns
    -------
    param_dict : dict
        Parameters related to the telescope and antenna layout, to be included in the
        :class:`pyuvdata.UVData` object.

        * `Nants_data`: Number of antennas
        * `Nants_telescope`: Number of antennas
        * `antenna_names`: list of antenna names
        * `antenna_numbers`: corresponding list of antenna numbers
        * `antenna_positions`: Array of ECEF antenna positions
        * `telescope_location`: ECEF array center location
        * `telescope_location_lat_lon_alt`: Lat Lon Alt array center location
        * `telescope_config_file`: Path to configuration yaml file
        * `antenna_location_file`: Path to csv layout file
        * `telescope_name`: observatory name
        * `x_orientation`: old, prefer feed_array + feed_angle + mount_type.
            physical orientation of the dipole labelled as "x"
        * `feed_array`: array of feeds (e.g. ['x', 'y']). Must be
        * `feed_angle`: array of feed angles. Interpretation depends on mount_type
            see pyuvdata docs.
        * `mount_type`: antenna mount type, see pyuvdata docs. Defaults to "fixed"
            if feed_array and feed_angle are passed.
        * `world`: Either "earth" or "moon" (requires lunarsky).
    beam_list : :class:`pyuvsim.BeamList`
        This is a BeamList object, only returned if return_beams is True.
    beam_dict : dict
        Antenna numbers to beam indices, only returned if return_beams is True.

    """
    tele_params = copy.deepcopy(tele_params)
    # check for telescope config
    tele_config = "telescope_config_name" in tele_params
    if tele_config:
        # parse telescope config
        if not os.path.isdir(config_path):
            config_path = os.path.dirname(config_path)
            if not os.path.isdir(config_path):
                raise ValueError(f"config_path {config_path} is not a directory")
        telescope_config_name = tele_params["telescope_config_name"]
        if not os.path.exists(telescope_config_name):
            telescope_config_name = os.path.join(config_path, telescope_config_name)
            if not os.path.exists(telescope_config_name):
                raise ValueError(
                    f"telescope_config_name file from yaml does not exist: {telescope_config_name}"
                )
        with open(telescope_config_name) as yf:
            telconfig = yaml.safe_load(yf)

        telescope_location_latlonalt, world, ellipsoid = _get_tel_loc(telconfig)

        tele_params["telescope_name"] = telconfig["telescope_name"]
    else:
        # if not provided, get bare-minumum keys from tele_params
        if "telescope_location" not in tele_params:
            raise KeyError(
                "If telescope_config_name not provided in `telescope` obsparam section, "
                "you must provide telescope_location"
            )
        if "telescope_name" not in tele_params:
            raise KeyError(
                "If telescope_config_name not provided in `telescope` obsparam section, "
                "you must provide telescope_name"
            )

        telescope_location_latlonalt, world, ellipsoid = _get_tel_loc(tele_params)

    lat_rad = telescope_location_latlonalt[0] * np.pi / 180.0
    long_rad = telescope_location_latlonalt[1] * np.pi / 180.0
    alt = telescope_location_latlonalt[2]

    if world == "moon":
        frame = "mcmf"
    elif world == "earth" or world is None:
        frame = "itrs"
    else:
        raise ValueError(f"Invalid world {world}")
    tele_params["telescope_location"] = uvutils.XYZ_from_LatLonAlt(
        latitude=lat_rad,
        longitude=long_rad,
        altitude=alt,
        frame=frame,
        ellipsoid=ellipsoid,
    )

    telescope_name = tele_params["telescope_name"]

    antnames, antnums, east, north, up, layout_csv, beam_ids = _get_array_layout(
        tele_params=tele_params, config_path=config_path
    )

    # fill in outputs with just array info
    return_dict = {}
    return_dict["Nants_data"] = antnames.size
    return_dict["Nants_telescope"] = antnames.size
    return_dict["antenna_names"] = np.array(antnames.tolist())
    return_dict["antenna_numbers"] = np.array(antnums)
    antpos_enu = np.vstack((east, north, up)).T
    return_dict["antenna_positions"] = (
        uvutils.ECEF_from_ENU(
            antpos_enu,
            latitude=lat_rad,
            longitude=long_rad,
            altitude=alt,
            frame=frame,
            ellipsoid=ellipsoid,
        )
        - tele_params["telescope_location"]
    )
    if world is not None:
        return_dict["world"] = world
    return_dict["telescope_frame"] = frame
    return_dict["ellipsoid"] = ellipsoid

    return_dict["array_layout"] = layout_csv
    return_dict["telescope_location"] = np.asarray(tele_params["telescope_location"])
    return_dict["telescope_location_lat_lon_alt"] = np.asarray(
        telescope_location_latlonalt
    )
    return_dict["telescope_name"] = telescope_name

    if not tele_config:
        beam_info = False
        for key in ["feed_array", "feed_angle", "mount_type"]:
            if key in tele_params:
                return_dict[key] = tele_params[key]
                beam_info = True
        if not beam_info:
            # if no info on beams, just return what we have
            # if telescope has feed angle, warn
            msg = (
                "No beam information, so cannot determine telescope mount_type, "
                "feed_array or feed_angle. Specify a telescope config file in "
                "the obs param file to get beam information or, if calling "
                "from `initialize_uvdata_from_keywords`, specify mount_type, "
                "feed_array and feed_angle."
            )
            warnings.warn(msg)

        if not return_beams:
            return return_dict
        else:
            return return_dict, BeamList([]), {}

    # if provided, parse sections related to beam files and types
    return_dict["telescope_config_name"] = telescope_config_name
    beam_ids_inc = np.unique(beam_ids)
    beam_dict = {}

    beam_list = _construct_beam_list(
        beam_ids_inc, telconfig, freq_array=freq_array, freq_range=freq_range
    )

    # construct feed_angles if appropriate
    # use dtype=object so strings don't get truncated.
    mount_type = np.full((antnames.size,), None, dtype=object)
    feed_angle = np.zeros((antnames.size, beam_list[0].beam.Nfeeds), dtype=float)
    # feed_arrays have to be the same -- checked in consistency checker
    feed_array = np.repeat(
        beam_list[0].beam.feed_array[np.newaxis, :], antnames.size, axis=0
    )

    for beam_ind, beam_id in enumerate(beam_ids_inc):
        wh_this_beam = np.nonzero(beam_ids == beam_id)
        which_ants = antnames[wh_this_beam]
        for ant in which_ants:
            beam_dict[str(ant)] = beam_ind
        mount_type[wh_this_beam] = beam_list[beam_ind].beam.mount_type
        feed_angle[wh_this_beam, :] = beam_list[beam_ind].beam.feed_angle

    if np.any(mount_type):
        if not np.all(mount_type):
            mount_type[np.nonzero(mount_type == np.array(None))] = "other"
        return_dict["mount_type"] = mount_type
    return_dict["feed_angle"] = feed_angle
    return_dict["feed_array"] = feed_array
    _ = return_dict.pop("x_orientation", None)

    if not return_beams:
        return return_dict
    else:
        return return_dict, beam_list, beam_dict


def _setup_coord_arrays(coord_params, param_map, tols):
    kws_given = ", ".join(sorted([param_map[coord_key] for coord_key in coord_params]))
    if "coord_array" in coord_params:
        # If the coordinate array is provided just use it.
        kwds_used = ["coord_array"]
        coord_params["coord_array"] = np.asarray(coord_params["coord_array"])
        coord_params["coord_n"] = coord_params["coord_array"].size
        coord_params["coord_start"] = coord_params["coord_array"][0]
        coord_params["coord_end"] = coord_params["coord_array"][-1]
        if "coord_delta" in coord_params:
            # If the deltas are provided use them
            coord_params["coord_delta"] = np.atleast_1d(coord_params["coord_delta"])
            if coord_params["coord_delta"].size == 1:
                coord_params["coord_delta"] = np.full(
                    (coord_params["coord_n"],),
                    coord_params["coord_delta"][0],
                    dtype=float,
                )
            elif coord_params["coord_delta"].size != coord_params["coord_n"]:
                raise ValueError(
                    f"If {param_map['coord_delta']} has multiple elements, the "
                    f"{param_map['coord_delta']} must be the same length as "
                    f"{param_map['coord_array']}."
                )
            kwds_used.append("coord_delta")
        else:
            # Can only calculate the delta if coords are evenly spaced
            if coord_params["coord_n"] > 1:
                if not uvutils.tools._test_array_constant_spacing(
                    coord_params["coord_array"], tols=tols
                ):
                    raise ValueError(
                        f"{param_map['coord_delta']} must be specified if spacing "
                        f"in {param_map['coord_array']} is uneven."
                    )
                coord_params["coord_delta"] = np.full(
                    (coord_params["coord_n"],),
                    np.mean(np.diff(coord_params["coord_array"])),
                    dtype=float,
                )
            else:
                raise ValueError(
                    f"{param_map['coord_delta']} must be specified if "
                    f"{param_map['coord_array']} has length 1"
                )
        # length is the *full length including widths* so delta is needed to
        # calculate the lengths
        coord_params["coord_length"] = (
            coord_params["coord_array"][-1] + coord_params["coord_delta"][-1] * 0.5
        ) - (coord_params["coord_array"][0] - coord_params["coord_delta"][0] * 0.5)
    else:
        # The final linspace call needs coord_start, coord_end and coord_n.
        # Calculate any of these we need from the other available parameters
        # Also calculate the other parameters so we can test for consistency
        # on any extra parameters that are provided
        if "coord_start" not in coord_params and "coord_end" not in coord_params:
            raise ValueError(
                f"Either {param_map['coord_start']} or {param_map['coord_end']} "
                "must be specified. The parameters that were specified were: "
                + kws_given
            )
        if (
            "coord_delta" in coord_params
            and np.asarray(coord_params["coord_delta"]).size > 1
        ):
            raise ValueError(
                f"{param_map['coord_delta']} must be a scalar if "
                f"{param_map['coord_array']} is not specified"
            )
        if "coord_start" in coord_params and "coord_end" in coord_params:
            # Both start and end are specified
            kwds_used = ["coord_start", "coord_end"]
            if "coord_n" not in coord_params and "coord_delta" not in coord_params:
                raise ValueError(
                    f"If both {param_map['coord_start']} and {param_map['coord_end']} "
                    f"are specified, either {param_map['coord_n']} or "
                    f"{param_map['coord_delta']} must be specified as well. The "
                    "parameters that were specified were: " + kws_given
                )
            if "coord_n" in coord_params:
                kwds_used.append("coord_n")
                extent = coord_params["coord_end"] - coord_params["coord_start"]
                if coord_params["coord_n"] > 1:
                    coord_delta = extent / (coord_params["coord_n"] - 1)
                    coord_length = extent + coord_delta
                    coord_params["coord_delta"] = coord_delta
                elif (
                    "coord_delta" not in coord_params
                    and "coord_length" not in coord_params
                ):
                    raise ValueError(
                        f"If {param_map['coord_n']} is 1 then either "
                        f"{param_map['coord_delta']} or {param_map['coord_length']} "
                        "must be specified."
                    )
                elif "coord_delta" in coord_params:
                    kwds_used.remove("coord_end")
                    kwds_used.append("coord_delta")
                    coord_params["coord_delta"] = np.atleast_1d(
                        coord_params["coord_delta"]
                    )
                    coord_length = coord_params["coord_delta"][0]
                else:
                    kwds_used.remove("coord_end")
                    kwds_used.append("coord_length")
                    coord_length = coord_params["coord_length"]
                    coord_params["coord_delta"] = np.atleast_1d(
                        coord_params["coord_length"]
                    )
            else:
                kwds_used.append("coord_delta")
                coord_length = (
                    coord_params["coord_end"]
                    - coord_params["coord_start"]
                    + coord_params["coord_delta"]
                )
                coord_params["coord_n"] = float(
                    coord_length / coord_params["coord_delta"]
                )
            coord_params["coord_length"] = coord_length
        else:
            # Only one of start & end are specified
            other_specs = [
                "coord_n" in coord_params,
                "coord_delta" in coord_params,
                "coord_length" in coord_params,
            ]
            if np.sum(other_specs) < 2:
                raise ValueError(
                    f"If only one of {param_map['coord_start']} and "
                    f"{param_map['coord_end']} is specified, two more of "
                    f"{param_map['coord_n']}, {param_map['coord_delta']} and "
                    f"{param_map['coord_length']} must be specified as well. "
                    "The parameters that were specified were: " + kws_given
                )
            if "coord_n" in coord_params:
                kwds_used = ["coord_n"]
                if "coord_delta" in coord_params:
                    kwds_used.append("coord_delta")
                    coord_length = coord_params["coord_n"] * coord_params["coord_delta"]
                    coord_params["coord_length"] = coord_length
                else:
                    kwds_used.append("coord_length")
                    coord_params["coord_delta"] = (
                        coord_params["coord_length"] / coord_params["coord_n"]
                    )
            else:
                kwds_used = ["coord_delta", "coord_length"]
                coord_params["coord_n"] = (
                    coord_params["coord_length"] / coord_params["coord_delta"]
                )

            if "coord_start" in coord_params:
                kwds_used.append("coord_start")
                coord_params["coord_end"] = (
                    coord_params["coord_start"]
                    + coord_params["coord_length"]
                    - coord_params["coord_delta"]
                )
            else:
                kwds_used.append("coord_end")
                coord_params["coord_start"] = (
                    coord_params["coord_end"]
                    - coord_params["coord_length"]
                    + coord_params["coord_delta"]
                )

        if not np.isclose(coord_params["coord_n"], np.round(coord_params["coord_n"])):
            raise ValueError(
                f"{param_map['coord_end']} - {param_map['coord_start']} must be "
                f"evenly divisible by {param_map['coord_delta']}"
            )
        coord_params["coord_n"] = int(np.round(coord_params["coord_n"]))

        coord_params["coord_array"], coord_delta = np.linspace(
            coord_params["coord_start"],
            coord_params["coord_end"],
            coord_params["coord_n"],
            retstep=True,
            endpoint=True,
        )
        if coord_params["coord_n"] > 1:
            # Set the delta to be the linspace step size unless there's only
            # one value, in which case use the value that was passed in or
            # calculated earlier
            coord_params["coord_delta"] = np.full(
                (coord_params["coord_n"],), coord_delta, dtype=float
            )
        else:
            # ensure changes in start/end are detected
            coord_params["coord_start"] = coord_params["coord_array"][0]
            coord_params["coord_end"] = coord_params["coord_array"][-1]
            # ensure delta is an array
            coord_params["coord_delta"] = np.atleast_1d(coord_params["coord_delta"])
    coord_params["params_used"] = kwds_used
    return coord_params


def parse_frequency_params(freq_params):
    """
    Parse the "freq" section of obsparam.

    Parameters
    ----------
    freq_params : dict
        Dictionary of frequency parameters.
        See pyuvsim documentation for examples of allowable key combinations.
        https://pyuvsim.readthedocs.io/en/latest/parameter_files.html#frequency

    Returns
    -------
    dict
        Dictionary of :class:`pyuvdata.UVData` parameters related to frequency:

            * `channel_width`: (dtype float, ndarray, shape=(Nfreqs)) Frequency channel
              widths in Hz
            * `Nfreqs`: (int) Number of frequencies
            * `freq_array`: (dtype float, ndarray, shape=(Nfreqs)) Frequency channel
              centers in Hz

    """
    freq_params_orig = copy.deepcopy(freq_params)
    # remove Nspws because we don't use it we set it to 1 it at the end anyway
    if "Nspws" in freq_params_orig:
        del freq_params_orig["Nspws"]
    freq_params = copy.deepcopy(freq_params_orig)
    ftols = UVData()._freq_array.tols

    freq_param_map = {
        "coord_array": "freq_array",
        "coord_delta": "channel_width",
        "coord_start": "start_freq",
        "coord_end": "end_freq",
        "coord_n": "Nfreqs",
        "coord_length": "bandwidth",
    }

    freq_param_map_rev = {
        freq_key: coord_key for coord_key, freq_key in freq_param_map.items()
    }

    coord_params = {
        freq_param_map_rev[freq_key]: value for freq_key, value in freq_params.items()
    }
    coord_params = _setup_coord_arrays(
        coord_params, param_map=freq_param_map, tols=ftols
    )

    freq_params = {
        freq_key: coord_params[coord_key]
        for coord_key, freq_key in freq_param_map.items()
    }
    kwds_used = [freq_param_map[coord_key] for coord_key in coord_params["params_used"]]

    # Warn about inconsistencies in supplied parameters
    inconsistent_params = []
    for kwd in freq_params_orig:
        if kwd == "channel_width" or kwd == "freq_array":
            orig_arr = np.atleast_1d(freq_params_orig[kwd])
            if orig_arr.size != freq_params[kwd].size:
                orig_arr = np.full(
                    (freq_params["Nfreqs"],), freq_params_orig[kwd], dtype=float
                )
            if not np.allclose(
                orig_arr, freq_params[kwd], rtol=ftols[0], atol=ftols[1]
            ):
                inconsistent_params.append(kwd)
        else:
            if not np.isclose(
                freq_params_orig[kwd], freq_params[kwd], rtol=ftols[0], atol=ftols[1]
            ):
                inconsistent_params.append(kwd)

    if len(inconsistent_params) > 0:
        warnings.warn(
            f"The {', '.join(inconsistent_params)} is not consistent with "
            f"the {', '.join(kwds_used)} specified. Using the values calculated "
            f"from the {', '.join(kwds_used)} parameters. Input values were: "
            f"{[freq_params_orig[kwd] for kwd in inconsistent_params]}, calculated "
            f"values were: {[freq_params[kwd] for kwd in inconsistent_params]}."
        )

    # return the things needed for UVData objects
    return_dict = {}
    return_dict["Nfreqs"] = freq_params["Nfreqs"]
    return_dict["freq_array"] = freq_params["freq_array"]
    return_dict["channel_width"] = freq_params["channel_width"]
    return_dict["Nspws"] = 1

    return return_dict


def parse_time_params(time_params):
    """
    Parse the "time" section of obsparam.

    Parameters
    ----------
    time_params : dict
        Dictionary of time parameters
        See pyuvsim documentation for examples of allowable key combinations.
        https://pyuvsim.readthedocs.io/en/latest/parameter_files.html#time

    Returns
    -------
    dict
        Dictionary of :class:`pyuvdata.UVData` parameters related to time:

            * `integration_time`: (float) Time array spacing in seconds.
            * `Ntimes`: (int) Number of times
            * `start_time`: (float) Starting time in Julian Date
            * `time_array`: (dtype float, ndarray, shape=(Ntimes,)) Time step centers in JD.

    """
    init_time_params = copy.deepcopy(time_params)
    # Handle the time_offset: add it to all absolute time parameters unless
    # it looks like they already have it, in which case warn!
    if "time_offset" in init_time_params:
        for param in ["start_time", "end_time", "time_array"]:
            if param in init_time_params:
                if np.max(init_time_params[param]) > init_time_params["time_offset"]:
                    warnings.warn(
                        f"time_offset is present, but {param} is larger than "
                        f"time_offset, so not adding time_offset to {param}."
                    )
                else:
                    if isinstance(init_time_params[param], list):
                        init_time_params[param] = (
                            np.asarray(init_time_params[param])
                            + init_time_params["time_offset"]
                        ).tolist()
                    else:
                        init_time_params[param] += init_time_params["time_offset"]
        del init_time_params["time_offset"]

    time_params = copy.deepcopy(init_time_params)
    dtols = UVData()._time_array.tols
    stols = UVData()._integration_time.tols

    # deal with the different units for different time parameters
    daysperhour = 1 / 24.0
    hourspersec = 1 / 60.0**2
    dayspersec = daysperhour * hourspersec
    if "integration_time" in time_params:
        time_params["int_time_days"] = (
            np.asarray(time_params.pop("integration_time")) * dayspersec
        )

    if "duration_hours" in time_params:
        if "duration_days" in time_params:
            warnings.warn(
                "Both duration_hours and duration_days are specified, using duration_days."
            )
            time_params.pop("duration_hours")
        else:
            time_params["duration_days"] = (
                time_params.pop("duration_hours") * daysperhour
            )

    time_param_map = {
        "coord_array": "time_array",
        "coord_delta": "int_time_days",
        "coord_start": "start_time",
        "coord_end": "end_time",
        "coord_n": "Ntimes",
        "coord_length": "duration_days",
    }

    time_param_map_rev = {
        freq_key: coord_key for coord_key, freq_key in time_param_map.items()
    }

    coord_params = {
        time_param_map_rev[freq_key]: value for freq_key, value in time_params.items()
    }

    # fix the names in errors/warnings to be the original names
    time_param_map_for_func = copy.deepcopy(time_param_map)
    time_param_map_for_func["coord_delta"] = "integration_time"
    time_param_map_for_func["coord_length"] = "duration"

    coord_params = _setup_coord_arrays(
        coord_params, param_map=time_param_map_for_func, tols=dtols
    )

    time_params = {
        freq_key: coord_params[coord_key]
        for coord_key, freq_key in time_param_map.items()
    }
    kwds_used = [time_param_map[coord_key] for coord_key in coord_params["params_used"]]

    # Go back to the right units
    time_params["integration_time"] = time_params.pop("int_time_days") / dayspersec
    kwds_used = [
        "integration_time" if kwd == "int_time_days" else kwd for kwd in kwds_used
    ]
    if "duration_days" in kwds_used and "duration_days" not in init_time_params:
        kwds_used = [
            "duration_hours" if kwd == "duration_days" else kwd for kwd in kwds_used
        ]

    # Warn about inconsistencies in supplied parameters
    inconsistent_params = []
    for kwd in init_time_params:
        if kwd == "integration_time" or kwd == "time_array":
            orig_arr = np.atleast_1d(init_time_params[kwd])
            if orig_arr.size != time_params[kwd].size:
                orig_arr = np.full(
                    (time_params["Ntimes"],), init_time_params[kwd], dtype=float
                )
            if not np.allclose(
                orig_arr, time_params[kwd], rtol=stols[0], atol=stols[1]
            ):
                inconsistent_params.append(kwd)
        elif kwd == "duration_hours":
            orig_duration_days = init_time_params[kwd] * daysperhour
            if not np.isclose(
                orig_duration_days,
                time_params["duration_days"],
                rtol=dtols[0],
                atol=dtols[1],
            ):
                inconsistent_params.append(kwd)
                time_params[kwd] = time_params["duration_days"] / daysperhour
        else:
            if not np.isclose(
                init_time_params[kwd], time_params[kwd], rtol=dtols[0], atol=dtols[1]
            ):
                inconsistent_params.append(kwd)

    if len(inconsistent_params) > 0:
        warnings.warn(
            f"The {', '.join(inconsistent_params)} is not consistent with "
            f"the {', '.join(kwds_used)} specified. Using the values calculated "
            f"from the {', '.join(kwds_used)} parameters. Input values were: "
            f"{[init_time_params[kwd] for kwd in inconsistent_params]}, calculated "
            f"values were: {[time_params[kwd] for kwd in inconsistent_params]}."
        )

    # return the things needed for UVData objects
    return_dict = {}
    return_dict["integration_time"] = time_params["integration_time"]
    return_dict["time_array"] = time_params["time_array"]
    return_dict["Ntimes"] = time_params["Ntimes"]
    return_dict["start_time"] = time_params["start_time"]

    return return_dict


def _coord_arrays_to_params(coord_array, coord_delta, tols):
    if coord_array.size > 1 and (np.asarray(coord_delta)).size == 1:
        coord_delta = np.full_like(coord_array, coord_delta)

    if (
        not uvutils.tools._test_array_constant_spacing(coord_array, tols=tols)
        or not uvutils.tools._test_array_constant(coord_delta, tols=tols)
        or not uvutils.tools._test_array_consistent(coord_array, coord_delta, tols=tols)
        or coord_array.size == 1
    ):
        # arrays are not evenly spaced or deltas do not match diffs or only one
        # value in the array. Just use the arrays and deltas.
        return {
            "coord_array": coord_array.tolist(),
            "coord_delta": coord_delta.tolist(),
        }
    else:
        # Evenly spaced arrays and delta matches array diffs
        # Use start/end and number because these are easier for people to
        # read and understand
        return {
            "coord_start": float(coord_array[0]),
            "coord_end": float(coord_array[-1]),
            "coord_n": coord_array.size,
        }


def freq_array_to_params(freq_array, channel_width):
    """
    Get a set of parameters that can be used to generate a given frequency array.

    Parameters
    ----------
    freq_array : array of float
        Frequencies in Hz.
    channel_width : array of float
        Channel widths in Hz.

    Returns
    -------
    dict
        Dictionary of frequency parameters consistent with freq_array.
        (channel_width, Nfreqs, bandwidth, start_freq, end_freq)
        See pyuvsim documentation for details:
        https://pyuvsim.readthedocs.io/en/latest/parameter_files.html#frequency

    """
    freq_param_map = {
        "coord_array": "freq_array",
        "coord_delta": "channel_width",
        "coord_start": "start_freq",
        "coord_end": "end_freq",
        "coord_n": "Nfreqs",
        "coord_length": "bandwidth",
    }

    coord_params = _coord_arrays_to_params(
        coord_array=freq_array,
        coord_delta=channel_width,
        tols=UVData()._freq_array.tols,
    )

    return {
        freq_param_map[coord_key]: value for coord_key, value in coord_params.items()
    }


def time_array_to_params(time_array, integration_time):
    """
    Get a set of parameters that can be used to generate a given time array.

    Returns a dictionary of parameters that can be used by parse_time_params
    to obtain the time_array passed in.

    Parameters
    ----------
    time_array : array of float
        Julian dates.
    integration_time : array of float
        Integration times in seconds.

    Returns
    -------
    dict
        Dictionary of time parameters consistent with time_array.
        (integration_time, Ntimes, duration_days, start_time, end_times)
        See pyuvsim documentation for details:
        https://pyuvsim.readthedocs.io/en/latest/parameter_files.html#time

    """
    # times maybe supplied as length Nblts with lots of duplicates.
    # Find unique times.
    time_array, unique_inverse = np.unique(np.asarray(time_array), return_inverse=True)

    # If there are duplicate times, we need to downselect the integration time array
    # to match. Easy if integration times are all the same, hard otherwise.
    if time_array.size < unique_inverse.size:
        # there were some repeated times
        int_time_use = np.full_like(time_array, np.nan)
        if not uvutils.tools._test_array_constant(
            integration_time, tols=UVData()._integration_time.tols
        ):
            int_time_varies_per_time = False
            for tind in range(time_array.size):
                int_times_this_time = integration_time[
                    np.nonzero(unique_inverse == tind)
                ]
                if not uvutils.tools._test_array_constant(
                    int_times_this_time, tols=UVData()._integration_time.tols
                ):
                    int_time_varies_per_time = True
                    int_time_use[tind] = np.min(int_times_this_time)
                else:
                    int_time_use[tind] = np.mean(int_times_this_time)
            if int_time_varies_per_time:
                warnings.warn(
                    "integration time varies for unique times, using the shortest "
                    "integration time for each unique time."
                )
        else:
            int_time_use = np.full_like(time_array, np.mean(integration_time))
        integration_time = int_time_use

    # deal with the different units for different time parameters
    daysperhour = 1 / 24.0
    hourspersec = 1 / 60.0**2
    dayspersec = daysperhour * hourspersec
    int_time_days = integration_time * dayspersec

    time_param_map = {
        "coord_array": "time_array",
        "coord_delta": "int_time_days",
        "coord_start": "start_time",
        "coord_end": "end_time",
        "coord_n": "Ntimes",
        "coord_length": "duration_days",
    }

    coord_params = _coord_arrays_to_params(
        coord_array=time_array,
        coord_delta=int_time_days,
        tols=UVData()._time_array.tols,
    )

    time_params = {
        time_param_map[coord_key]: value for coord_key, value in coord_params.items()
    }

    # go back to the right units
    if "int_time_days" in time_params:
        time_params["integration_time"] = (
            np.asarray(time_params.pop("int_time_days")) / dayspersec
        ).tolist()

    # use time_offset because if these dicts are dumped to yamls we lose too
    # much precision on the absolute time parameters
    if "time_array" in time_params:
        # Set the offset as the JD at midnight of the first time
        time_params["time_offset"] = float(
            np.floor(time_params["time_array"][0] - 0.5) + 0.5
        )
        time_params["time_array"] = (
            np.asarray(time_params["time_array"]) - time_params["time_offset"]
        ).tolist()
    else:
        # Set the offset as the JD at midnight of the first time
        time_params["time_offset"] = float(
            np.floor(time_params["start_time"] - 0.5) + 0.5
        )
        time_params["start_time"] = (
            time_params["start_time"] - time_params["time_offset"]
        )
        time_params["end_time"] = time_params["end_time"] - time_params["time_offset"]

    return time_params


def subselect(uv_obj, param_dict):
    """
    Do selection on a UVData object.

    Parameters
    ----------
    uv_obj : UVData
        UVData object to do the selection on.
    param_dict : dict
        Dict of parameters that can include a `select` section.

    """
    # downselect baselines (or anything that can be passed to pyuvdata's select method)
    # Note: polarization selection is allowed here, but will cause an error if the incorrect pols
    # are passed to pyuvsim.
    if "select" not in param_dict:
        return

    valid_select_keys = [
        "antenna_nums",
        "antenna_names",
        "ant_str",
        "frequencies",
        "freq_chans",
        "times",
        "blt_inds",
    ]

    select_params = param_dict["select"]
    no_autos = bool(select_params.pop("no_autos", False))
    select_params = {k: v for k, v in select_params.items() if k in valid_select_keys}
    if "antenna_nums" in select_params:
        select_params["antenna_nums"] = list(map(int, select_params["antenna_nums"]))
    redundant_threshold = param_dict["select"].get("redundant_threshold", None)
    if len(select_params) > 0:
        uv_obj.select(**select_params)

    if no_autos:
        uv_obj.select(ant_str="cross")

    if redundant_threshold is not None:
        uv_obj.compress_by_redundancy(tol=redundant_threshold, use_grid_alg=True)


def set_ordering(uv_obj, param_dict, reorder_blt_kw):
    """
    Do conjugation/reordering on a UVData object.

    Parameters
    ----------
    uv_obj : UVData
        UVData object to reorder.
    param_dict : dict
        Dict of parameters that can include an `order` section.
    reorder_blt_kw : dict
        Dict giving order & minor order for blt reordering.

    """
    bl_conjugation_convention = None

    # this is the order required by the simulator. Might as well default to it.
    default_blt_order = ["time", "baseline"]
    if "ordering" in param_dict:
        ordering_dict = param_dict["ordering"]
        if "blt_order" not in ordering_dict:
            ordering_dict["blt_order"] = default_blt_order
    else:
        ordering_dict = {"blt_order": default_blt_order}

    bl_conjugation_convention = ordering_dict.pop("conjugation_convention", None)
    if reorder_blt_kw is None:
        blt_order = ordering_dict["blt_order"]

        if isinstance(blt_order, str):
            reorder_blt_kw = {"order": blt_order}
        if isinstance(blt_order, tuple | list | np.ndarray):
            reorder_blt_kw = {"order": blt_order[0]}
            if len(blt_order) > 1:
                reorder_blt_kw["minor_order"] = blt_order[1]
    if bl_conjugation_convention is None:
        warnings.warn(
            "The default baseline conjugation convention has changed. In the past "
            "it was 'ant2<ant1', it now defaults to 'ant1<ant2'. You can specify "
            "the baseline conjugation convention in `obs_param` by setting the "
            "obs_param['ordering']['conjugation_convention'] field. This warning will go "
            "away in version 1.5."
        )
    elif bl_conjugation_convention != "ant1<ant2":
        uv_obj.conjugate_bls(convention=bl_conjugation_convention)

    obj_blt_order = uv_obj.blt_order
    if obj_blt_order is not None:
        obj_blt_order = {"order": uv_obj.blt_order[0]}
        if len(uv_obj.blt_order) > 1:
            obj_blt_order["minor_order"] = uv_obj.blt_order[1]

    if reorder_blt_kw and reorder_blt_kw != obj_blt_order:
        uv_obj.reorder_blts(**reorder_blt_kw)


def _initialize_polarization_helper(param_dict, uvparam_dict, beam_list=None):
    # There does not seem to be any way to get polarization_array into uvparam_dict, so
    # let's add it explicitly.
    if param_dict.get("polarization_array", None) is not None:
        if beam_list is not None:
            polstr_list = uvutils.polnum2str(param_dict["polarization_array"])
            feed_set = set()
            for polstr in polstr_list:
                feed_set = feed_set.union(uvutils.pol.POL_TO_FEED_DICT[polstr])
            for feed in feed_set:
                beam_feeds = beam_list[0].feed_array
                if feed not in beam_feeds:
                    raise ValueError(
                        f"Specified polarization array {param_dict['polarization_array']} "
                        f"requires feeds {feed_set} but the beams have feeds {beam_feeds}."
                    )

        uvparam_dict["polarization_array"] = np.array(param_dict["polarization_array"])

    # Parse polarizations
    if uvparam_dict.get("polarization_array") is None:
        if beam_list is not None:
            # default polarization array based on the beam feeds
            uvparam_dict["polarization_array"] = uvutils.pol.convert_feeds_to_pols(
                feed_array=beam_list[0].feed_array, include_cross_pols=True
            )
        else:
            uvparam_dict["polarization_array"] = np.array([-5, -6, -7, -8])

    if "Npols" not in uvparam_dict:
        uvparam_dict["Npols"] = len(uvparam_dict["polarization_array"])

    return uvparam_dict


def initialize_uvdata_from_params(
    obs_params, return_beams=False, reorder_blt_kw=None, check_kw=None
):
    """
    Construct a :class:`pyuvdata.UVData` object from parameters in a valid yaml file.

    Sufficient information must be provided by the parameters to define time and frequency arrays
    and verify the channel widths and time steps. This will error if insufficient or incompatible
    parameters are defined.

    The parameter dictionary may contain any valid :class:`pyuvdata.UVData` attributes as well.

    If the polarization array is not specified, it defaults to (XX, XY, YX, YY).

    Parameters
    ----------
    obs_params : dict or str
        Either an obs_param file name or a dictionary of parameters read in.
        Additional :class:`pyuvdata.UVData` parameters may be passed in through here.
    return_beams : bool
        Option to return the beam_list and beam_dict. Default is False.
    reorder_blt_kw : dict (optional)
        Deprecated. Use obs_param['ordering']['blt_order'] instead.
        Keyword arguments to send to the ``uvdata.reorder_blts`` method at the end of
        the setup process. Typical parameters include "order" and "minor_order". Default
        values are ``order='time'`` and ``minor_order='baseline'``.
    check_kw : dict (optional)
        A dictionary of keyword arguments to pass to the :func:`uvdata.UVData.check`
        method after object creation. Caution: turning off critical checks can
        result in a UVData object that cannot be written to a file.


    Returns
    -------
    uv_obj : :class:`pyuvdata.UVData`
        Initialized UVData object.
    beam_list : :class:`pyuvsim.BeamList`
        List of beam specifiers as strings, only returned if return_beams is True.
    beam_dict : dict
        Map of antenna numbers to index in beam_list, only returned if return_beams is
        True.

    """
    logger.info("Starting SimSetup")

    if reorder_blt_kw is not None:
        warnings.warn(
            "The reorder_blt_kw parameter is deprecated in favor of setting "
            "obs_param['ordering']['blt_order']. This will become an error in "
            "version 1.5",
            DeprecationWarning,
        )

    uvparam_dict = {}  # Parameters that will go into UVData
    if isinstance(obs_params, str):
        param_dict = _config_str_to_dict(obs_params)  # Container for received settings.
    else:
        param_dict = copy.deepcopy(obs_params)

    logger.info("Finished reading obsparams")

    # Parse frequency structure
    freq_dict = param_dict["freq"]
    parsed_freq_dict = parse_frequency_params(freq_dict)
    uvparam_dict["freq_array"] = parsed_freq_dict["freq_array"]
    uvparam_dict["channel_width"] = parsed_freq_dict["channel_width"]

    freq_array = parsed_freq_dict["freq_array"]
    logger.info("Finished reading frequency info")

    # Parse telescope parameters
    tele_dict = param_dict["telescope"]
    beamselect = tele_dict.pop("select", {})
    if len(beamselect) > 0:
        warnings.warn(
            "Beam selections should be specified in the telescope "
            "configuration, not in the obsparam. This will become an error in "
            "version 1.5",
            DeprecationWarning,
        )

    freq_buffer = beamselect.pop("freq_buffer", None)
    if freq_buffer:
        freq_range = (freq_array.min() - freq_buffer, freq_array.max() + freq_buffer)
    else:
        freq_range = None

    ptp_ret = parse_telescope_params(
        tele_dict,
        config_path=param_dict["config_path"],
        freq_array=freq_array,
        freq_range=freq_range,
        return_beams=return_beams,
    )
    if return_beams:
        tele_params, beam_list, beam_dict = ptp_ret
    else:
        tele_params = ptp_ret

    logger.info("Finished Setup of BeamList")

    # Use extra_keywords to pass along required paths for file history.
    extra_keywords = {}
    if "obs_param_file" in param_dict:
        extra_keywords["obsparam"] = param_dict["obs_param_file"]
        if "telescope_config_name" in tele_dict and isinstance(
            tele_dict["telescope_config_name"], str
        ):
            extra_keywords["telecfg"] = tele_dict["telescope_config_name"]
        if "array_layout" in tele_dict and isinstance(tele_dict["array_layout"], str):
            extra_keywords["layout"] = tele_dict["array_layout"]

    if "world" in tele_params:
        extra_keywords["world"] = tele_params.get("world")
    uvparam_dict["extra_keywords"] = extra_keywords

    # Parse time structure
    time_dict = param_dict["time"]
    parsed_time_dict = parse_time_params(time_dict)
    uvparam_dict["times"] = parsed_time_dict["time_array"]
    uvparam_dict["integration_time"] = parsed_time_dict["integration_time"]
    logger.info("Finished Setup of Time Dict")

    # figure out polarization
    if return_beams:
        bl_pass = beam_list
    else:
        bl_pass = None
    uvparam_dict = _initialize_polarization_helper(
        param_dict, uvparam_dict, beam_list=bl_pass
    )

    # telescope frame is set from world in parse_telescope_params.
    # Can only be itrs or (if lunarsky is installed) mcmf
    if tele_params["telescope_frame"] == "itrs":
        telescope_location = EarthLocation.from_geocentric(
            *tele_params["telescope_location"], unit="m"
        )
    elif tele_params["telescope_frame"] == "mcmf":
        # to get here, lunarsky has to be installed, so don't need to test for it
        # It's checked in the earlier call to parse_telescope_params, which sets
        # tele_params["telescope_frame"] to either "itrs" or "mcmf"
        from lunarsky import MoonLocation

        telescope_location = MoonLocation.from_selenocentric(
            *tele_params["telescope_location"], unit="m"
        )
        telescope_location.ellipsoid = tele_params["ellipsoid"]

    if "cat_name" not in param_dict and "object_name" not in param_dict:
        cat_name = "unprojected"
    elif "object_name" in param_dict:
        cat_name = param_dict["object_name"]
    else:
        cat_name = param_dict["cat_name"]
    phase_center_catalog = {0: {"cat_name": cat_name, "cat_type": "unprojected"}}

    # Set the antpairs _before_ creating the uvdata object to conserve memory.
    antpairs = param_dict.get("select", {}).pop("bls", None)
    if isinstance(antpairs, str):
        antpairs = ast.literal_eval(antpairs)

    tel_init_params = {"location": telescope_location}

    telescope_param_map = {
        "telescope_name": "name",
        "antenna_names": "antenna_names",
        "antenna_numbers": "antenna_numbers",
        "antenna_positions": "antenna_positions",
        "antenna_diameters": "antenna_diameters",
        "x_orientation": "x_orientation",
        "feed_angle": "feed_angle",
        "feed_array": "feed_array",
        "mount_type": "mount_type",
        "instrument": "instrument",
    }
    for key, tel_key in telescope_param_map.items():
        if key in param_dict["telescope"]:
            tel_init_params[tel_key] = param_dict["telescope"][key]
        if key in tele_params:
            tel_init_params[tel_key] = tele_params[key]
    if "instrument" not in tel_init_params:
        tel_init_params["instrument"] = tel_init_params["name"]

    uv_obj = UVData.new(
        telescope=Telescope.new(**tel_init_params),
        phase_center_catalog=phase_center_catalog,
        vis_units="Jy",
        history="",
        do_blt_outer=True,
        time_axis_faster_than_bls=False,
        antpairs=antpairs,
        check_kw=check_kw,
        **uvparam_dict,
    )

    logger.info(f"  Baseline Array: {uv_obj.baseline_array.nbytes / 1024**3:.2f} GB")
    logger.info(f"      Time Array: {uv_obj.time_array.nbytes / 1024**3:.2f} GB")
    logger.info(f"Integration Time: {uv_obj.integration_time.nbytes / 1024**3:.2f} GB")
    logger.info(f"       Ant1Array: {uv_obj.ant_1_array.nbytes / 1024**3:.2f} GB")
    logger.info(f"BLT-ORDER: {uv_obj.blt_order}")
    logger.info("Initialized UVData object")

    subselect(uv_obj, param_dict)
    logger.info("After Select")

    set_ordering(uv_obj, param_dict, reorder_blt_kw)

    logger.info(f"BLT-ORDER: {uv_obj.blt_order}")
    logger.info("After Re-order BLTS")

    logger.info("After Check")
    logger.info(f"BLT-ORDER: {uv_obj.blt_order}")
    if return_beams:
        return uv_obj, beam_list, beam_dict
    else:
        return uv_obj


def _complete_uvdata(uv_in, inplace=False):
    """
    Initialize data-like arrays on a a :class:`pyuvdata.UVData` object.

    This will overwrite any existing data in `uv_in` unless inplace=True.

    Parameters
    ----------
    uv_in : :class:`pyuvdata.UVData` instance
        Usually an incomplete object, containing only metadata.
    inplace : bool, optional
        Whether to perform the filling on the passed object, or a copy.

    Returns
    -------
    :class:`pyuvdata.UVData` : object with initialized data-like arrays
        (if `inplace` is `True`, it is on a copy of the input). With zeroed
        data_array, no flags and nsample_array of all ones.

    """
    if not inplace:
        uv_obj = copy.deepcopy(uv_in)
    else:
        uv_obj = uv_in

    # moved to have meta-data compatible objects from the init
    # but retain this functionality for now.
    # in most cases should be a no-op anyway.
    if uv_obj.lst_array is None:
        uv_obj.set_lsts_from_time_array()

    # Clear existing data, if any.
    _shape = (uv_obj.Nblts, uv_obj.Nfreqs, uv_obj.Npols)
    uv_obj.data_array = np.zeros(_shape, dtype=complex)
    uv_obj.flag_array = np.zeros(_shape, dtype=bool)
    uv_obj.nsample_array = np.ones(_shape, dtype=float)

    return uv_obj


def initialize_uvdata_from_keywords(
    output_yaml_filename=None,
    antenna_layout_filepath=None,
    output_layout_filename=None,
    array_layout=None,
    telescope_location=None,
    telescope_name=None,
    feed_array=None,
    feed_angle=None,
    mount_type=None,
    Nfreqs=None,
    start_freq=None,
    bandwidth=None,
    freq_array=None,
    channel_width=None,
    Ntimes=None,
    integration_time=None,
    start_time=None,
    time_array=None,
    bls=None,
    antenna_nums=None,
    antenna_names=None,
    polarization_array=None,
    no_autos=False,
    redundant_threshold=None,
    conjugation_convention=None,
    blt_order=None,
    write_files=True,
    path_out=None,
    complete=False,
    check_kw: dict | None = None,
    **kwargs,
):
    """
    Initialize a :class:`pyuvdata.UVData` object from keyword arguments.

    Optionally, write out the configuration to YAML and CSV files such that
    :func:`~initialize_uvdata_from_params` will produce the same :class:`pyuvdata.UVData` object.

    Parameters
    ----------
    output_yaml_filename : str (optional)
        Specify filename for yaml file to write out.
        Defaults to obsparam.yaml
    antenna_layout_filepath : str (optional)
        Path to csv file of antenna positions (see documentation for details).
    output_layout_filename : str (optional)
        File name for antenna layout if writing files out.
        If unspecified, defaults to the basename of the antenna_layout_filepath.
        If that is not given, it defaults to "antenna_layout.csv"
        Will not overwrite existing files.
    array_layout : dictionary (required if antenna_layout_filepath not given)
        keys are integer antenna numbers, values are len-3 antenna positions
        in ENU coordinates [meters].
    antenna_names : list of str (optional)
        If unset, antenna names are assigned as "%s" % antnum.
    telescope_location : array of float, shape (3,)
        Telescope location on Earth in LatLonAlt coordinates [deg, deg, meters]
    telescope_name : str
        Name of telescope
    feed_array : array-like of str or None
        List of feeds for each antenna in the telescope, must be one of
        "x", "y", "l", "r". Shape (Nants, Nfeeds), dtype str. Can also be shape
        (Nfeeds,), in which case the same values for feed_array are used for
        all antennas in the object.
    feed_angle : array-like of float or None
        Orientation of the feed with respect to zenith (or with respect to north if
        pointed at zenith). Units is in rads, vertical polarization is nominally 0,
        and horizontal polarization is nominally pi / 2. Shape (Nants, Nfeeds),
        dtype float.  Can also be shape (Nfeeds,), in which case the same values for
        feed_angle are used for all antennas in the object.
    mount_type : str or array-like of str
        Antenna mount type, which describes the optics of the antenna in question.
        Supported options include: "alt-az" (primary rotates in azimuth and
        elevation), "equatorial" (primary rotates in hour angle and declination)
        "orbiting" (antenna is in motion, and its orientation depends on orbital
        parameters), "x-y" (primary rotates first in the plane connecting east,
        west, and zenith, and then perpendicular to that plane),
        "alt-az+nasmyth-r" ("alt-az" mount with a right-handed 90-degree tertiary
        mirror), "alt-az+nasmyth-l" ("alt-az" mount with a left-handed 90-degree
        tertiary mirror), "phased" (antenna is "electronically steered" by
        summing the voltages of multiple elements, e.g. MWA), "fixed" (antenna
        beam pattern is fixed in azimuth and elevation, e.g., HERA), and "other"
        (also referred to in some formats as "bizarre"). See the "Conventions"
        page of the documentation for further details. Shape (Nants,), dtype str.
        Can also provide a single string, in which case the same mount_type is
        used for all antennas in the object.
    Nfreqs : int
        Number of frequency channels
    start_freq : float
        Starting frequency [Hz]
    bandwidth : float
        Total frequency bandwidth of spectral window [Hz]
    channel_width: float
        Frequency channel spacing [Hz]
    freq_array : ndarray
        frequency array [Hz], if this is specified, it supersedes Nfreqs,
        start_freq, bandwidth
    Ntimes : int
        Number of integration bins
    integration_time : float
        Width of time bins [seconds]
    start_time : float
        Time of the first integration bin [Julian Date]
    time_array : ndarray
        time array [Julian Date]. If this is specified it supersedes values of Ntimes,
        start_time and integration_time
    bls : list
        List of antenna-pair tuples for baseline selection
    redundant_threshold: float
        Redundant baseline selection tolerance for selection [meters]
    antenna_nums : list
        List of antenna numbers to keep in array
    polarization_array : list
        List of polarization strings (or ints) to insert into object
    no_autos : bool
        If True, eliminate all auto correlations
    conjugation_convention : str, optional
        A convention for the directions of the baselines, options are:
        'ant1<ant2', 'ant2<ant1', 'u<0', 'u>0', 'v<0', 'v>0'
    blt_order : str or array-like of str, optional
        Specifies the ordering along the blt axis. Defaults to ['time', 'baseline'].
        If set to anything else, it is passed to :meth:`pyuvdata.UVData.reorder_blts`.
        A single string specifies the order, a list of two strings specifies the
        [order, minor_order]. See the docs on the UVData method for more information.
    write_files : bool
        If True, write out the parameter information to yaml files.
    path_out : str (optional)
        Path in which to place generated configuration files, if write_files is True.
        Defaults to current directory.
    complete : bool (optional)
        Whether to fill out the :class:`pyuvdata.UVData` object with its requisite
        data arrays, and check if it's all consistent.
    check_kw : dict (optional)
        A dictionary of keyword arguments to pass to the :func:`uvdata.UVData.check`
        method after object creation. Caution: turning off critical checks can
        result in a UVData object that cannot be written to a file.
    kwargs : dictionary
        Any additional valid :class:`pyuvdata.UVData` attribute to assign to object.

    Returns
    -------
    :class:`pyuvdata.UVData`
        Initialized based on keywords, with a zeroed data_array, no flags and
        nsample_array of all ones.

    """
    arrfile = antenna_layout_filepath is not None
    outfile = output_layout_filename is not None

    if array_layout is None and antenna_layout_filepath is None:
        raise ValueError(
            "Either array_layout or antenna_layout_filepath must be passed."
        )

    if (
        array_layout is None or isinstance(array_layout, str) or write_files
    ) and path_out is None:
        path_out = "."

    if write_files:
        if not outfile:
            if not arrfile:
                output_layout_filename = "antenna_layout.csv"
            else:
                output_layout_filename = os.path.basename(antenna_layout_filepath)
            outfile = True

        # Increment name appropriately:
        output_layout_filepath = os.path.join(path_out, output_layout_filename)
        output_layout_filename = os.path.basename(
            check_file_exists_and_increment(output_layout_filepath)
        )

        if output_yaml_filename is None:
            output_yaml_filename = "obsparam.yaml"
        output_yaml_filename = check_file_exists_and_increment(
            os.path.join(path_out, output_yaml_filename)
        )

        # Copying original file to new place, if it exists
        if antenna_layout_filepath is not None and os.path.exists(
            antenna_layout_filepath
        ):
            shutil.copyfile(
                antenna_layout_filepath, os.path.join(path_out, output_layout_filename)
            )

    antenna_numbers = None
    if isinstance(array_layout, dict):
        antenna_numbers = np.fromiter(array_layout.keys(), dtype=int)
        antpos_enu = array_layout.values()
        if antenna_names is None:
            antenna_names = antenna_numbers.astype("str")
        if write_files:
            _write_layout_csv(
                output_layout_filepath, antpos_enu, antenna_names, antenna_numbers
            )

    if array_layout is None:
        if write_files and output_layout_filename is not None:
            array_layout = output_layout_filename
        else:
            array_layout = antenna_layout_filepath

    freq_params = {
        "Nfreqs": Nfreqs,
        "start_freq": start_freq,
        "bandwidth": bandwidth,
        "freq_array": freq_array,
        "channel_width": channel_width,
    }
    time_params = {
        "Ntimes": Ntimes,
        "start_time": start_time,
        "integration_time": integration_time,
        "time_array": time_array,
    }
    selection_params = {
        "bls": bls,
        "redundant_threshold": redundant_threshold,
        "antenna_nums": antenna_nums,
        "no_autos": no_autos,
    }
    tele_params = {
        "telescope_location": repr(tuple(telescope_location)),
        "telescope_name": telescope_name,
        "feed_array": feed_array,
        "feed_angle": feed_angle,
        "mount_type": mount_type,
    }
    layout_params = {
        "antenna_names": antenna_names,
        "antenna_numbers": antenna_numbers,
        "array_layout": array_layout,
    }

    ordering_params = {
        "conjugation_convention": conjugation_convention,
        "blt_order": blt_order,
    }

    freq_params = {k: v for k, v in freq_params.items() if v is not None}
    time_params = {k: v for k, v in time_params.items() if v is not None}
    selection_params = {k: v for k, v in selection_params.items() if v is not None}
    tele_params = {k: v for k, v in tele_params.items() if v is not None}
    layout_params = {k: v for k, v in layout_params.items() if v is not None}
    ordering_params = {k: v for k, v in ordering_params.items() if v is not None}

    valid_param_names = [getattr(UVData(), param).name for param in UVData()]

    extra_kwds = {k: v for k, v in kwargs.items() if k in valid_param_names}

    # Convert str polarization array to int.
    if polarization_array is not None and type(polarization_array[0]) is not int:
        polarization_array = np.array(uvutils.polstr2num(polarization_array))

    if output_yaml_filename is None:
        output_yaml_filename = ""

    param_dict = {
        "time": time_params,
        "freq": freq_params,
        "select": selection_params,
        "ordering": ordering_params,
        "telescope": tele_params,
        "config_path": path_out,
        "polarization_array": polarization_array,
    }

    param_dict.update(**extra_kwds)

    if write_files:
        tele_params["array_layout"] = antenna_layout_filepath
        param_dict["telescope"] = tele_params
        writeable_param_dict = copy.deepcopy(param_dict)
        if polarization_array is not None:
            writeable_param_dict["polarization_array"] = polarization_array.tolist()
        with open(output_yaml_filename, "w") as yfile:
            yaml.dump(writeable_param_dict, yfile, default_flow_style=False)

    param_dict["obs_param_file"] = os.path.basename(output_yaml_filename)
    param_dict["telescope"].update(layout_params)
    uv_obj = initialize_uvdata_from_params(
        param_dict, return_beams=False, check_kw=check_kw
    )

    if complete:
        _complete_uvdata(uv_obj, inplace=True)

    return uv_obj


def uvdata_to_telescope_config(
    uvdata_in,
    beam_filepath,
    layout_csv_name=None,
    telescope_config_name=None,
    return_names=False,
    path_out=".",
):
    """
    Make telescope parameter files from a :class:`pyuvdata.UVData` object.

    Makes both a telescope_config file and a layout_csv file:

    * telescope_config: YAML file with telescope_location and telescope_name
        The beam list is spoofed, since that information cannot be found in a UVData object.
    * layout_csv: tab separated value file giving ENU antenna positions.
        Beam ID is spoofed as well.

    See https://pyuvsim.readthedocs.io/en/latest/parameter_files.html for details.

    Parameters
    ----------
    uvdata_in : pyuvdata.UVData
        UVData object for which to make config files.
    path_out : str
        Target directory for the config files.
    beam_filepath : str
        Path to a beamfits file.
    layout_csv_name : str
        The name for the antenna positions.
        Default <telescope_name>_layout.csv, where <telescope_name> is uvdata_in.telescope_name
    telescope_config_name : str
        The name for the telescope config file
        Default telescope_config_<telescope_name>.yaml
    return_names : bool
        Return the file names. Default is False. Used in tests.

    Returns
    -------
    tuple
        if return_names, returns (path, telescope_config_name, layout_csv_name)

    """
    tel_name = uvdata_in.telescope.name
    antpos_enu = uvdata_in.telescope.get_enu_antpos()
    ant_names = uvdata_in.telescope.antenna_names
    ant_numbers = uvdata_in.telescope.antenna_numbers
    tel_lla_deg = uvdata_in.telescope.location_lat_lon_alt_degrees

    # fix formating issue for numpy types in numpy>2.0
    tel_lla_deg = tuple(np.asarray(tel_lla_deg).tolist())

    if telescope_config_name is None:
        telescope_config_path = check_file_exists_and_increment(
            os.path.join(path_out, f"telescope_config_{tel_name}.yaml")
        )
        telescope_config_name = os.path.basename(telescope_config_path)

    if layout_csv_name is None:
        layout_csv_path = check_file_exists_and_increment(
            os.path.join(path_out, tel_name + "_layout.csv")
        )
        layout_csv_name = os.path.basename(layout_csv_path)

    _write_layout_csv(
        os.path.join(path_out, layout_csv_name), antpos_enu, ant_names, ant_numbers
    )

    # create a nearly empty beam object, only defining the filepath and
    # mount_type for the yaml dumper
    beam = UVBeam()
    beam.filename = [beam_filepath]
    beam.mount_type = uvdata_in.telescope.mount_type[0]

    # Write the rest to a yaml file.
    yaml_dict = {
        "telescope_name": tel_name,
        "telescope_location": repr(tel_lla_deg),
        "beam_paths": {0: beam},
    }

    with open(os.path.join(path_out, telescope_config_name), "w+") as yfile:
        yaml.safe_dump(yaml_dict, yfile, default_flow_style=False)

    logger.info(
        f"Path: {path_out}, telescope_config: {telescope_config_name}, layout: {layout_csv_name}"
    )

    if return_names:
        return path_out, telescope_config_name, layout_csv_name


def uvdata_to_config_file(
    uvdata_in,
    param_filename=None,
    telescope_config_name="",
    layout_csv_name="",
    catalog="mock",
    path_out=".",
):
    """
    Extract simulation configuration settings from a UVData object.

    When used with :func:`~uvdata_to_telescope_config`, this will produce all the necessary
    configuration yaml and csv file to make an "empty" :class:`pyuvdata.UVData` object comparable to
    the argument `uvdata_in`. The generated file will match `uvdata_in` in frequency, time,
    antenna positions, and uvw coordinates. Note that it may be different in the
    baseline conjugation convention and the ordering of the baseline-time axis,
    but those can be modified using the :meth:`pyuvdata.UVData.conjugate_bls` and
    :meth:`pyuvdata.UVData.reorder_blt` methods.

    Parameters
    ----------
    uvdata_in: :class:`pyuvdata.UVData`
        UVData object for which to make config files.
    param_filename: str
        output param file name, defaults to obsparam_#.yaml.
    telescope_config_name: str
        Name of telescope configuration yaml file. Defaults to blank string.
    layout_csv_name: str
        Name of antenna layout csv file. Defaults to blank string.
    catalog: str
        Path to catalog file. Defaults to 'mock'.
    path_out: str
        Where to put config files.

    """
    if param_filename is None:
        param_filename = check_file_exists_and_increment(
            os.path.join(path_out, "obsparam.yaml")
        )
        param_filename = os.path.basename(param_filename)

    tdict = time_array_to_params(
        time_array=uvdata_in.time_array, integration_time=uvdata_in.integration_time
    )
    fdict = freq_array_to_params(
        freq_array=uvdata_in.freq_array, channel_width=uvdata_in.channel_width
    )

    param_dict = {
        "time": tdict,
        "freq": fdict,
        "sources": {"catalog": catalog},
        "telescope": {
            "telescope_config_name": telescope_config_name,
            "array_layout": layout_csv_name,
        },
        "filing": {"outdir": ".", "outfile_name": "", "outfile_prefix": ""},
        # Add this in to prevent warnings and make it explicit, these are just
        # the default values.
        "ordering": {
            "conjugation_convention": "ant1<ant2",
            "blt_order": ["time", "baseline"],
        },
    }

    if catalog == "mock":
        param_dict["sources"]["mock_arrangement"] = "zenith"

    with open(os.path.join(path_out, param_filename), "w") as yfile:
        yaml.safe_dump(param_dict, yfile, default_flow_style=False)
