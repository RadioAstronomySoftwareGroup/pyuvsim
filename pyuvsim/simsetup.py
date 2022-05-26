# -*- mode: python; coding: utf-8 -*
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
import copy
import os
import shutil
import warnings

import numpy as np
import yaml
from packaging import version  # packaging is installed with setuptools
from astropy.coordinates import Angle, EarthLocation, Latitude, Longitude, AltAz, ICRS
import astropy.units as units
from pyuvdata import utils as uvutils, UVData
import pyradiosky

from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
from .analyticbeam import AnalyticBeam
from .telescope import BeamList
from .astropy_interface import SkyCoord, MoonLocation, Time, LunarTopo, hasmoon
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
    with open(layout_csv, 'r') as fhandle:
        header = fhandle.readline()

    header = [h.strip() for h in header.split()]
    str_format_code = 'U'

    # get data types for each column
    dtypes = {
        "name": str_format_code + '10',
        'number': 'i4',
        'beamid': 'i4',
        'e': 'f8',
        'n': 'f8',
        'u': 'f8'
    }

    # check columns in file
    lower_header = [col.lower() for col in header]
    columns = ['name', 'number', 'beamid', 'e', 'n', 'u']
    col_exist = [col for col in columns if col in lower_header]

    dt = np.format_parser([dtypes[col] for col in col_exist],
                          col_exist, header)

    return np.genfromtxt(layout_csv, autostrip=True, skip_header=1,
                         dtype=dt.dtype)


def _write_layout_csv(filepath, antpos_enu, antenna_names, antenna_numbers, beam_ids=None):
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
    with open(filepath, 'w') as lfile:
        lfile.write(header + '\n')
        for i, (e, n, u) in enumerate(antpos_enu):
            beam_id = beam_ids[i]
            name = antenna_names[i]
            num = antenna_numbers[i]
            line = (
                "{:{}} {:8d} {:8d} {:10.4f} {:10.4f} {:10.4f}\n").format(
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
    with open(config_str, 'r') as pfile:
        param_dict = yaml.safe_load(pfile)

    config_str = os.path.abspath(config_str)
    param_dict['config_path'] = os.path.dirname(config_str)
    param_dict['obs_param_file'] = os.path.basename(config_str)

    return param_dict


def _set_lsts_on_uvdata(uv_obj):
    """Set the LSTs on a :class:`pyuvdata.UVData` object, with handling for MoonLocations."""
    # If the telescope location is a MoonLocation,
    # then uv_obj.extra_keywords['world'] == 'moon'.
    world = uv_obj.extra_keywords.get('world', 'earth')

    if world == 'earth':
        uv_obj.set_lsts_from_time_array()
    elif world == 'moon':
        if not hasmoon:
            raise ValueError("Cannot construct lsts for MoonLocation without lunarsky module")
        un_jds, inv = np.unique(uv_obj.time_array, return_inverse=True)
        loc = MoonLocation(*uv_obj.telescope_location, unit='m')
        times = Time(un_jds, format='jd', scale='utc', location=loc)
        uv_obj.lst_array = times.sidereal_time('apparent').rad[inv]
    else:
        raise ValueError(f"Invalid world {world}.")
    return uv_obj


def _create_catalog_diffuse(
    map_nside, diffuse_model, diffuse_params, time, localframe, array_location
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

    if localframe == 'lunartopo':
        frame = LunarTopo(location=array_location, obstime=time)
    else:
        frame = AltAz(location=array_location, obstime=time)
    hpix = astropy_healpix.HEALPix(nside=map_nside, frame=ICRS())
    icrs_coord = hpix.healpix_to_skycoord(pixinds)
    source_coord = icrs_coord.transform_to(frame)

    alts = source_coord.alt.rad
    azs = source_coord.az.rad

    modfunc = analytic_diffuse.get_model(diffuse_model)
    fluxes = modfunc(source_coord.az.rad, source_coord.zen.rad, **diffuse_params)

    stokes = np.zeros((4, 1, npix))
    stokes[0, :] = fluxes

    stokes *= units.K

    if version.parse(pyradiosky.__version__) > version.parse("0.1.2"):
        catalog = pyradiosky.SkyModel(
            stokes=stokes,
            nside=map_nside,
            hpx_inds=pixinds,
            spectral_type='flat',
            frame="icrs"
        )
    else:
        catalog = pyradiosky.SkyModel(
            stokes=stokes,
            nside=map_nside,
            hpx_inds=pixinds,
            spectral_type='flat',
        )

    return catalog, icrs_coord, alts, azs, fluxes


def _create_catalog_discrete(Nsrcs, alts, azs, fluxes, time, localframe, array_location):
    """
    Make a :class:`pyradiosky.SkyModel` object from a set of point sources.

    Only used for internal testing, should not be called by users.

    Parameters
    ----------
    Nsrcs : int
        Number of sources
    alts :  array_like of float
        Altitudes or Latitudes (depending on localframe) of the sources in radians.
    azs :  array_like of float
        Azimuths or Longitues (depending on localframe) of the sources in radians.
    fluxes : array_like of float
        Brightnesses of the sources in Jy.
    time : :class:`astropy.time.Time` object
        Time to use in defining the positions.
    localframe : :class:`astropy.coordinates.BaseCoordinateFrame` class or str
        Frame that the source positions are defined in (e.g. altaz or ICRS)
    array_location : :class:`astropy.coordinates.EarthLocation` or :class:`lunarsky.MoonLocation`
        Location to use in defining the positions.

    Returns
    -------
    catalog : class:`pyradiosky.SkyModel`
        SkyModel object containing the diffuse map.
    icrs_coord : :class:`astropy.coordinates.SkyCoord`
        Astropy SkyCoord object containing the coordinates for the sources in ICRS.

    """
    source_coord = SkyCoord(alt=Angle(alts, unit=units.deg), az=Angle(azs, unit=units.deg),
                            obstime=time, frame=localframe, location=array_location)
    icrs_coord = source_coord.transform_to('icrs')

    names = np.array(['src' + str(si) for si in range(Nsrcs)])
    stokes = np.zeros((4, 1, Nsrcs))
    stokes[0, :] = fluxes

    stokes *= units.Jy

    if version.parse(pyradiosky.__version__) > version.parse("0.1.2"):
        catalog = pyradiosky.SkyModel(
            name=names,
            ra=icrs_coord.ra,
            dec=icrs_coord.dec,
            stokes=stokes,
            spectral_type='flat',
            frame="icrs"
        )
    else:
        catalog = pyradiosky.SkyModel(
            name=names,
            ra=icrs_coord.ra,
            dec=icrs_coord.dec,
            stokes=stokes,
            spectral_type='flat'
        )

    return catalog, icrs_coord


def create_mock_catalog(time, arrangement='zenith', array_location=None, Nsrcs=None,
                        alt=None, save=False, min_alt=None, rseed=None, return_data=False,
                        diffuse_model=None, diffuse_params=None, map_nside=None):
    """
    Create a mock catalog.

    SkyModels are defined in an AltAz frame at the given time, then returned in
    ICRS ra/dec coordinates.

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
        time = Time(time, scale='utc', format='jd')

    if array_location is None:
        array_location = EarthLocation(lat='-30d43m17.5s', lon='21d25m41.9s',
                                       height=1073.)

    if arrangement not in ['off-zenith', 'zenith', 'cross', 'triangle', 'long-line', 'hera_text',
                           'random', 'diffuse']:
        raise KeyError("Invalid mock catalog arrangement: " + str(arrangement))

    mock_keywords = {
        'time': time.jd, 'arrangement': arrangement,
        'array_location': repr(
            (array_location.lat.deg, array_location.lon.deg, array_location.height.value))
    }

    if isinstance(array_location, MoonLocation):
        localframe = 'lunartopo'
        mock_keywords['world'] = 'moon'
    else:
        localframe = 'altaz'
        mock_keywords['world'] = 'earth'

    if arrangement == 'off-zenith':
        if alt is None:
            alt = 85.0  # Degrees
        mock_keywords['alt'] = alt
        Nsrcs = 1
        alts = [alt]
        azs = [90.]  # 0 = North pole, 90. = East pole
        fluxes = [1.0]

    if arrangement == 'triangle':
        Nsrcs = 3
        if alt is None:
            alt = 87.0  # Degrees
        mock_keywords['alt'] = alt
        alts = [alt, alt, alt]
        azs = [0., 120., 240.]
        fluxes = [1.0, 1.0, 1.0]

    if arrangement == 'cross':
        Nsrcs = 4
        alts = [88., 90., 86., 82.]
        azs = [270., 0., 90., 135.]
        fluxes = [5., 4., 1.0, 2.0]

    if arrangement == 'zenith':
        if Nsrcs is None:
            Nsrcs = 1
        mock_keywords['Nsrcs'] = Nsrcs
        alts = np.ones(Nsrcs) * 90.
        azs = np.zeros(Nsrcs, dtype=float)
        fluxes = np.ones(Nsrcs) * 1 / Nsrcs
        # Divide total Stokes I intensity among all sources
        # Test file has Stokes I = 1 Jy

    if arrangement == 'random':
        if Nsrcs is None:
            Nsrcs = 1
        if min_alt is None:
            min_alt = 30  # Degrees
        mock_keywords['Nsrcs'] = Nsrcs
        np.random.seed(seed=rseed)
        mock_keywords['rseed'] = np.random.get_state()[1][0]

        # Necessary to get uniform distribution per solid angle.
        min_alt_rad = np.radians(min_alt)
        rv = np.random.uniform(-1, np.cos(min_alt_rad + np.pi / 2), Nsrcs)
        alts = np.degrees(np.arccos(rv) - np.pi / 2)
        azs = np.degrees(np.random.uniform(0, 2 * np.pi, Nsrcs))
        fluxes = np.ones(Nsrcs, dtype=float)

    if arrangement == 'long-line':
        if Nsrcs is None:
            Nsrcs = 10
        if min_alt is None:
            min_alt = 5
        mock_keywords['Nsrcs'] = Nsrcs
        mock_keywords['alt'] = alt
        fluxes = np.ones(Nsrcs, dtype=float)
        if Nsrcs % 2 == 0:
            length = 180 - min_alt * 2
            spacing = length / (Nsrcs - 1)
            max_alt = 90. - spacing / 2
            alts = np.linspace(min_alt, max_alt, Nsrcs // 2)
            alts = np.append(alts, np.flip(alts, axis=0))
            azs = np.append(np.zeros(Nsrcs // 2, dtype=float) + 180.,
                            np.zeros(Nsrcs // 2, dtype=float))
        else:
            alts = np.linspace(min_alt, 90, (Nsrcs + 1) // 2)
            alts = np.append(alts, np.flip(alts[1:], axis=0))
            azs = np.append(np.zeros((Nsrcs + 1) // 2, dtype=float) + 180.,
                            np.zeros((Nsrcs - 1) // 2, dtype=float))

    if arrangement == 'hera_text':
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

        alts = 90. - zas
        Nsrcs = alts.size
        fluxes = np.ones_like(alts)

    if arrangement == 'diffuse':
        if analytic_diffuse is None or astropy_healpix is None:
            raise ValueError("analytic_diffuse and astropy_healpix "
                             "modules are required to use diffuse mock catalog.")
        if (diffuse_model is None):
            raise ValueError("Diffuse arrangement selected, but diffuse model not chosen.")
        if diffuse_params is None:
            diffuse_params = {}
        if map_nside is None:
            warnings.warn("No nside chosen. Defaulting to 32")
            map_nside = 32

        mock_keywords['diffuse_model'] = diffuse_model
        mock_keywords['map_nside'] = map_nside
        mock_keywords['diffuse_params'] = repr(diffuse_params)

        catalog, icrs_coord, alts, azs, fluxes = _create_catalog_diffuse(
            map_nside, diffuse_model, diffuse_params, time, localframe, array_location
        )

    else:

        catalog, icrs_coord = _create_catalog_discrete(
            Nsrcs, alts, azs, fluxes, time, localframe, array_location
        )

    if return_data:
        catalog = SkyModelData(catalog)

    if get_rank() == 0 and save:
        np.savez('mock_catalog_' + arrangement,
                 ra=icrs_coord.ra.rad,
                 dec=icrs_coord.dec.rad,
                 alts=alts,
                 azs=azs,
                 fluxes=fluxes)

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
    reference_frequency = None
    spectral_index = None
    polarized = None
    nside = None
    hpx_inds = None
    flux_unit = None

    put_in_shared = ['stokes_I', 'stokes_Q', 'stokes_U', 'stokes_V', 'polarized',
                     'ra', 'dec', 'reference_frequency', 'spectral_index', 'hpx_inds']

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
            if hasattr(sky_in, "get_lon_lat"):
                if sky_in.frame != "icrs":
                    sky_in.transform_to(ICRS)
                sky_ra, sky_dec = sky_in.get_lon_lat()
                self.ra = sky_ra.deg
                self.dec = sky_dec.deg
            else:
                # backwards compatibility for pyradiosky < 0.1.3
                self.ra = sky_in.ra.deg
                self.dec = sky_in.dec.deg
            self.Nfreqs = sky_in.Nfreqs
            stokes_in = sky_in.stokes

            if isinstance(stokes_in, units.Quantity):
                if stokes_in.unit.is_equivalent("Jy"):
                    stokes_in = stokes_in.to_value("Jy")
                    self.flux_unit = "Jy"
                elif stokes_in.unit.is_equivalent("K"):
                    stokes_in = stokes_in.to_value("K")
                    self.flux_unit = "K"

            self.stokes_I = stokes_in[0, ...]

            if sky_in._n_polarized > 0:
                self.polarized = sky_in._polarized
                Q, U, V = stokes_in[1:, :, sky_in._polarized]
                self.stokes_Q = Q
                self.stokes_U = U
                self.stokes_V = V

            if sky_in._freq_array.required:
                self.freq_array = sky_in.freq_array.to("Hz").value
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
        elif hasattr(sky_in, "filename"):
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
        if self.hpx_inds is not None:
            new_sky.hpx_inds = self.hpx_inds[inds]

        if self.polarized is not None:
            sub_inds = np.in1d(self.polarized, inds)
            new_sky.stokes_Q = self.stokes_Q[..., sub_inds]
            new_sky.stokes_U = self.stokes_U[..., sub_inds]
            new_sky.stokes_V = self.stokes_V[..., sub_inds]
            new_sky.polarized = np.where(np.in1d(inds, self.polarized))[0]

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
            raise ImportError("You need mpi4py to use this method. "
                              "Install it by running pip install pyuvsim[sim] "
                              "or pip install pyuvsim[all] if you also want the "
                              "line_profiler installed.")
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

        ra_use = Longitude(self.ra, unit='deg')
        dec_use = Latitude(self.dec, unit='deg')
        stokes_use = np.zeros((4, self.Nfreqs, self.Ncomponents), dtype=float)

        stokes_use[0, ...] = self.stokes_I

        if self.polarized is not None:
            stokes_use[1, :, self.polarized] = self.stokes_Q.T
            stokes_use[2, :, self.polarized] = self.stokes_U.T
            stokes_use[3, :, self.polarized] = self.stokes_V.T

        if self.flux_unit is not None:
            stokes_use *= units.Unit(self.flux_unit)

        other = {}
        if self.spectral_type in ['full', 'subband']:
            other['freq_array'] = self.freq_array * units.Hz
        if self.spectral_type == 'spectral_index':
            other['reference_frequency'] = self.reference_frequency * units.Hz
            other['spectral_index'] = self.spectral_index
        if self.component_type == 'healpix':
            other['nside'] = self.nside
            other['hpx_inds'] = self.hpx_inds
        else:
            other['name'] = self.name

        empty_sky = pyradiosky.SkyModel()
        if hasattr(empty_sky, "filename"):
            # this only works for pyradiosky >= 0.1.3
            other["filename"] = self.filename

        return pyradiosky.SkyModel(
            ra=ra_use, dec=dec_use, stokes=stokes_use,
            spectral_type=self.spectral_type, **other
        )


def _sky_select_calc_rise_set(sky, source_params, telescope_lat_deg=None):
    """
    Apply flux and non-rising source cuts, calculate rise and set lsts.

    Parameters
    ----------
    source_params : dict
        Dict specifying flux cut and horizon buffer parameters.
    telescope_lat_deg : float
        Latitude of telescope in degrees, used for horizon calculations.

    """
    if hasattr(sky, "cut_nonrising"):

        if telescope_lat_deg is not None:
            telescope_lat = Latitude(telescope_lat_deg * units.deg)
            sky.cut_nonrising(telescope_lat)
            calc_kwargs = {}
            if "horizon_buffer" in source_params:
                calc_kwargs["horizon_buffer"] = float(source_params["horizon_buffer"])
            sky.calculate_rise_set_lsts(telescope_lat, **calc_kwargs)

        min_brightness = None
        if "min_flux" in source_params:
            min_brightness = float(source_params["min_flux"]) * units.Jy
        max_brightness = None
        if "max_flux" in source_params:
            max_brightness = float(source_params["max_flux"]) * units.Jy

        sky.select(min_brightness=min_brightness, max_brightness=max_brightness)

    else:
        # backwards compatibility for pyradiosky < 0.1.3

        # Parse source selection options
        select_options = ['min_flux', 'max_flux', 'horizon_buffer']

        # Put catalog selection cuts in source section.
        source_select_kwds = {}
        for key in select_options:
            if key in source_params:
                source_select_kwds[key] = float(source_params[key])
        if telescope_lat_deg is not None:
            source_select_kwds['latitude_deg'] = telescope_lat_deg

        sky.source_cuts(**source_select_kwds)

    return sky


def _read_pyradiosky_file(filename, source_params, filetype=None):
    """
    Read in a pyradiosky catalog file.

    Parameters
    ----------
    filename : str
        File to read in.
    source_params : dict
        Dictionary specifying options to pyradiosky read methods.
    filetype : str
        One of ['skyh5', 'gleam', 'vot', 'text', 'hdf5'] or None.
        If None, the code attempts to guess what the file type is.

    Returns:
    --------
    :class:`pyradiosky.SkyModel`
        SkyModel object based on filename

    """
    sky = pyradiosky.SkyModel()

    allowed_filetypes = ['skyh5', 'gleam', 'vot', 'text', 'hdf5']
    if filetype is not None:
        if filetype not in allowed_filetypes:
            raise ValueError(f"Invalid filetype. Filetype options are: {allowed_filetypes}")
    else:
        if filename.endswith("txt"):
            filetype = "text"
        elif filename.endswith('vot'):
            if "gleam" in filename.lower():
                filetype = "gleam"
            else:
                filetype = "vot"
        elif filename.endswith("skyh5"):
            filetype = "skyh5"
        elif filename.endswith("hdf5"):
            filetype = "hdf5"

    if filetype == "text":
        sky.read_text_catalog(filename)
    elif filetype == "gleam":
        if "spectral_type" in source_params:
            spectral_type = source_params["spectral_type"]
        else:
            warnings.warn("No spectral_type specified for GLEAM, using 'flat'.")
            spectral_type = "flat"
        sky.read_gleam_catalog(
            filename, spectral_type=spectral_type
        )
    elif filetype == "vot":
        vo_params = {}
        required_params = ["table_name", "id_column", "flux_columns"]
        warn_params = {"ra_column": "RAJ2000", "dec_column": "DEJ2000"}
        for param in required_params:
            if param not in source_params:
                raise ValueError(f"No {param} specified for {filename}.")
            vo_params[param] = source_params[param]
        for param, default in warn_params.items():
            if param not in source_params:
                warnings.warn(
                    f"No {param} specified for {filename}, using default: {default}."
                )
                vo_params[param] = default
            else:
                vo_params[param] = source_params[param]

        vo_params["reference_frequency"] = None
        sky.read_votable_catalog(filename, **vo_params)
    elif filetype == "skyh5":
        sky.read_skyh5(filename)
    elif filetype == "hdf5":
        # In principal, this could either be a skyh5 or an old-style healvis hdf5 file
        # Try skyh5 first, if that doesn't work, use the old method
        try:
            sky.read_skyh5(filename)
        except ValueError:
            sky.read_healpix_hdf5(filename)
    else:
        raise ValueError(
            "The file extension is not recognized as a valid one for SkyModel objects. "
            "Expected extensions include: '.txt', '.vot', '.skyh5' and '.hdf5'."
        )

    return sky


def initialize_catalog_from_params(
    obs_params, input_uv=None, return_recarray=True, filetype=None, return_catname=None
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
    return_recarray: bool
        Return a recarray instead of a :class:`pyradiosky.SkyModel` instance.
    filetype : str
        One of ['skyh5', 'gleam', 'vot', 'text', 'hdf5'] or None.
        If None, the code attempts to guess what the file type is.
    return_catname: bool
        Return the catalog name. Defaults to True currently but will default to False
        in version 1.4.

    Returns
    -------
    skydata: numpy.recarray or :class:`pyradiosky.SkyModel`
        Source catalog filled with data.
    source_list_name: str
        Catalog identifier for metadata. Only returned if return_catname is True.

    """
    if return_catname is None:
        warnings.warn(
            "The return_catname parameter currently defaults to True, but starting in"
            "version 1.4 it will default to False.",
            DeprecationWarning
        )
        return_catname = True

    if input_uv is not None and not isinstance(input_uv, UVData):
        raise TypeError("input_uv must be UVData object")

    if isinstance(obs_params, str):
        with open(obs_params, 'r') as pfile:
            param_dict = yaml.safe_load(pfile)

        param_dict['config_path'] = os.path.dirname(obs_params)
    else:
        param_dict = obs_params

    source_params = param_dict['sources']
    if 'catalog' in source_params:
        catalog = source_params['catalog']
    else:
        raise KeyError("No catalog defined.")

    if catalog == 'mock':
        mock_keywords = {'arrangement': source_params['mock_arrangement']}
        extra_mock_kwds = ['time', 'Nsrcs', 'zen_ang', 'save', 'min_alt',
                           'array_location', 'diffuse_model', 'diffuse_params', 'map_nside']
        for k in extra_mock_kwds:
            if k in source_params.keys():
                if k == 'array_location':
                    # String -- lat, lon, alt in degrees
                    latlonalt = [float(s.strip()) for s in source_params[k].split(',')]
                    lat, lon, alt = latlonalt
                    mock_keywords[k] = EarthLocation.from_geodetic(lon, lat, alt)
                else:
                    mock_keywords[k] = source_params[k]

        # time, arrangement, array_location, save, Nsrcs, min_alt

        if 'array_location' not in mock_keywords:
            if input_uv is not None:
                mock_keywords['array_location'] = EarthLocation.from_geocentric(
                    *input_uv.telescope_location, unit='m')
            else:
                warnings.warn("No array_location specified. Defaulting to the HERA site.")
        if 'time' not in mock_keywords:
            if input_uv is not None:
                mock_keywords['time'] = input_uv.time_array[0]
                warnings.warn(
                    "No julian date given for mock catalog. Defaulting to first time step."
                )
            else:
                raise ValueError(
                    "input_uv must be supplied if using mock catalog without specified julian date"
                )

        time = mock_keywords.pop('time')

        sky, mock_keywords = create_mock_catalog(time, **mock_keywords)
        mock_keyvals = [str(key) + str(val) for key, val in mock_keywords.items()]
        source_list_name = 'mock_' + "_".join(mock_keyvals)
    elif isinstance(catalog, str):
        if not os.path.isfile(catalog):
            catalog = os.path.join(param_dict['config_path'], catalog)

        if filetype is None:
            if "filetype" in source_params:
                filetype = source_params["filetype"]

        source_list_name = os.path.basename(catalog)
        sky = _read_pyradiosky_file(catalog, source_params, filetype=filetype)

    if input_uv is not None:
        telescope_lat_deg = input_uv.telescope_location_lat_lon_alt_degrees[0]
    else:
        telescope_lat_deg = None
    # Do source selections, if any.
    sky = _sky_select_calc_rise_set(
        sky, source_params, telescope_lat_deg=telescope_lat_deg
    )

    # If the filename parameter doesn't exist on the sky object
    # (because of an older version of pyradiosky) or if it is None (e.g. for mock skies)
    # add the source_list_name to the object so it can be put in the UVData history later
    if getattr(sky, "filename", None) is None:
        sky.filename = source_list_name

    if return_recarray:
        warnings.warn("initialize_catalog_from_params will return a SkyModel instance, not "
                      "a recarray, by default in version 1.3. Set keyword "
                      "`return_recarray=True` to ensure recarray is returned.",
                      DeprecationWarning)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sky = sky.to_recarray()

    if return_catname:
        return sky, source_list_name
    else:
        return sky


def _construct_beam_list(beam_ids, telconfig):
    """Construct BeamList."""
    beam_list = []
    for beamID in beam_ids:
        beam_model = telconfig['beam_paths'][beamID]

        if not isinstance(beam_model, (str, dict)):
            raise ValueError('Beam model is not properly specified in telescope config file.')

        # first check to see if the beam_model is a string giving a file location
        if (isinstance(beam_model, str)
                and (os.path.exists(beam_model)
                     or os.path.exists(os.path.join(SIM_DATA_PATH, beam_model)))):
            if os.path.exists(beam_model):
                beam_list.append(beam_model)
            else:
                beam_list.append(os.path.join(SIM_DATA_PATH, beam_model))
        # Failing that, try to parse the beam string as an analytic beam.
        else:
            if isinstance(beam_model, str):
                beam_type = beam_model
            elif 'type' in beam_model:
                beam_type = beam_model['type']
            else:
                raise ValueError("Beam model must have a 'type' field.")

            if beam_type not in AnalyticBeam.supported_types:
                raise ValueError("Undefined beam model type: {}".format(beam_type))

            this_beam_opts = {}
            if isinstance(beam_model, dict):
                for key in beam_model:
                    if key != 'type':
                        this_beam_opts[key] = beam_model[key]

            # Gaussian beam requires either diameter or sigma
            # Airy beam requires diameter

            # Default to None for diameter and sigma.
            # Values in the "beam_paths" override globally-defined options.
            shape_opts = {'diameter': None, 'sigma': None}

            for opt in shape_opts.keys():
                shape_opts[opt] = this_beam_opts.get(opt, None)

            if all(v is None for v in shape_opts.values()):
                for opt in shape_opts.keys():
                    shape_opts[opt] = telconfig.get(opt, None)

            diameter = shape_opts.pop('diameter')
            sigma = shape_opts.pop('sigma')

            if beam_type == 'uniform':
                beam_model = 'analytic_uniform'

            if beam_type == 'gaussian':
                if diameter is not None:
                    beam_model = '_'.join(['analytic_gaussian', 'diam=' + str(diameter)])
                elif sigma is not None:
                    beam_model = '_'.join(['analytic_gaussian', 'sig=' + str(sigma)])
                else:
                    raise KeyError("Missing shape parameter for gaussian beam (diameter or sigma).")
            if beam_type == 'airy':
                if diameter is not None:
                    beam_model = '_'.join(['analytic_airy', 'diam=' + str(diameter)])
                else:
                    raise KeyError("Missing diameter for airy beam.")

            beam_list.append(beam_model)

    beam_list_obj = BeamList(beam_list=beam_list)
    if 'spline_interp_opts' in telconfig.keys():
        beam_list_obj.spline_interp_opts = telconfig['spline_interp_opts']

    for key in beam_list_obj.uvb_params.keys():
        if key in telconfig:
            beam_list_obj.uvb_params[key] = telconfig[key]

    return beam_list_obj


def parse_telescope_params(tele_params, config_path=''):
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
        * `x_orientation`: physical orientation of the dipole labelled as "x"
        * `world`: Either "earth" or "moon".
    beam_list : :class:`pyuvsim.BeamList`
        This is a BeamList object in string mode.
        The strings provide either paths to beamfits files or the specifications to make
        a :class:`pyuvsim.AnalyticBeam`.
    beam_dict : dict
        Antenna numbers to beam indices

    """
    tele_params = copy.deepcopy(tele_params)
    # check for telescope config
    tele_config = 'telescope_config_name' in tele_params
    if tele_config:
        # parse telescope config
        if not os.path.isdir(config_path):
            config_path = os.path.dirname(config_path)
            if not os.path.isdir(config_path):
                raise ValueError('config_path {} is not a directory'.format(config_path))
        telescope_config_name = tele_params['telescope_config_name']
        if not os.path.exists(telescope_config_name):
            telescope_config_name = os.path.join(config_path, telescope_config_name)
            if not os.path.exists(telescope_config_name):
                raise ValueError('telescope_config_name file from yaml does not exist')
        with open(telescope_config_name, 'r') as yf:
            telconfig = yaml.safe_load(yf)
        telescope_location_latlonalt = ast.literal_eval(telconfig['telescope_location'])
        world = telconfig.pop('world', None)
        tele_params['telescope_name'] = telconfig['telescope_name']
    else:
        # if not provided, get bare-minumum keys from tele_params
        if 'telescope_location' not in tele_params:
            raise KeyError(
                "If telescope_config_name not provided in `telescope` obsparam section, "
                "you must provide telescope_location"
            )
        if 'telescope_name' not in tele_params:
            raise KeyError(
                "If telescope_config_name not provided in `telescope` obsparam section, "
                "you must provide telescope_name"
            )
        telescope_location_latlonalt = tele_params['telescope_location']
        if isinstance(telescope_location_latlonalt, str):
            telescope_location_latlonalt = ast.literal_eval(telescope_location_latlonalt)
        world = tele_params.pop('world', None)

    telescope_location = np.array(telescope_location_latlonalt)
    telescope_location[0] *= np.pi / 180.
    telescope_location[1] *= np.pi / 180.  # Convert to radians
    tele_params['telescope_location'] = uvutils.XYZ_from_LatLonAlt(*telescope_location)
    telescope_name = tele_params['telescope_name']

    # get array layout
    if 'array_layout' not in tele_params:
        raise KeyError('array_layout must be provided.')
    array_layout = tele_params.pop('array_layout')
    if isinstance(array_layout, str):
        # Interpet as file path to layout csv file.
        layout_csv = array_layout
        # if array layout is a str, parse it as .csv filepath
        if isinstance(layout_csv, str):
            if not os.path.exists(layout_csv):
                layout_csv = os.path.join(config_path, layout_csv)
                if not os.path.exists(layout_csv):
                    raise ValueError(
                        'layout_csv file {} from yaml does not exist'.format(layout_csv))
            ant_layout = _parse_layout_csv(layout_csv)
            E, N, U = ant_layout['e'], ant_layout['n'], ant_layout['u']
            antnames = ant_layout['name']
            antnums = np.array(ant_layout['number'])
    elif isinstance(array_layout, dict):
        # Receiving antenna positions directly
        antnums = tele_params.pop('antenna_numbers')
        antnames = tele_params.pop('antenna_names')
        E, N, U = np.array([array_layout[an] for an in antnums]).T
        layout_csv = 'user-fed dict'

    # fill in outputs with just array info
    return_dict = {}
    beam_list = BeamList([])
    beam_dict = {}
    return_dict['Nants_data'] = antnames.size
    return_dict['Nants_telescope'] = antnames.size
    return_dict['antenna_names'] = np.array(antnames.tolist())
    return_dict['antenna_numbers'] = np.array(antnums)
    antpos_enu = np.vstack((E, N, U)).T
    return_dict['antenna_positions'] = (
        uvutils.ECEF_from_ENU(antpos_enu, *telescope_location)
        - tele_params['telescope_location'])
    if world is not None:
        return_dict['world'] = world

    return_dict['array_layout'] = layout_csv
    return_dict['telescope_location'] = np.asarray(tele_params['telescope_location'])
    return_dict['telescope_location_lat_lon_alt'] = np.asarray(telescope_location_latlonalt)
    return_dict['telescope_name'] = telescope_name

    # if provided, parse sections related to beam files and types
    if not tele_config:
        beam_dict = {n: 0 for n in antnames}
        return return_dict, beam_list, beam_dict

    return_dict['telescope_config_name'] = telescope_config_name
    beam_ids = ant_layout['beamid']
    beam_dict = {}

    for beamID in np.unique(beam_ids):
        which_ants = antnames[np.where(beam_ids == beamID)]
        for a in which_ants:
            beam_dict[a] = beamID

    beam_list = _construct_beam_list(np.unique(beam_ids), telconfig)

    return return_dict, beam_list, beam_dict


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
    freq_keywords = ['freq_array', 'start_freq', 'end_freq', 'Nfreqs', 'channel_width', 'bandwidth']
    fa, sf, ef, nf, cw, bw = [fk in freq_params for fk in freq_keywords]
    kws_used = ", ".join(sorted(freq_params.keys()))
    freq_params = copy.deepcopy(freq_params)

    if fa:
        freq_arr = np.asarray(freq_params['freq_array'])
        freq_params['Nfreqs'] = freq_arr.size
        if 'channel_width' not in freq_params:
            if freq_params['Nfreqs'] > 1:
                freq_params['channel_width'] = np.diff(freq_arr)[0]
                if not np.allclose(np.diff(freq_arr), freq_params['channel_width']):
                    raise ValueError(
                        "Channel width must be specified if spacing in freq_array is uneven."
                    )
            else:
                raise ValueError("Channel width must be specified if freq_array has length 1")
    else:
        if not (sf or ef):
            raise ValueError('Either start or end frequency must be specified: ' + kws_used)
        if cw:
            if np.asarray(freq_params['channel_width']).size > 1:
                raise ValueError("channel_width must be a scalar if freq_array is not specified")

        if not nf:
            if not cw:
                raise ValueError("Either channel_width or Nfreqs "
                                 " must be included in parameters:" + kws_used)

            if sf and ef:
                freq_params['bandwidth'] = (
                    freq_params['end_freq']
                    - freq_params['start_freq']
                    + freq_params['channel_width'])

                bw = True
            if bw:
                Nfreqs = float(freq_params['bandwidth'] / freq_params['channel_width'])
                freq_params['Nfreqs'] = Nfreqs

            else:
                raise ValueError("Either bandwidth or band edges must be specified: " + kws_used)

        if not cw:
            if not bw:
                raise ValueError("Either bandwidth or channel width must be specified: " + kws_used)
            freq_params['channel_width'] = (freq_params['bandwidth'] / float(freq_params['Nfreqs']))

        if not bw:
            freq_params['bandwidth'] = freq_params['channel_width'] * freq_params['Nfreqs']
            bw = True

        if not sf:
            if ef and bw:
                freq_params['start_freq'] = (
                    freq_params['end_freq']
                    - freq_params['bandwidth']
                    + freq_params['channel_width'])

        if not ef:
            if sf and bw:
                freq_params['end_freq'] = (
                    freq_params['start_freq']
                    + freq_params['bandwidth']
                    - freq_params['channel_width'])

        if not np.isclose(freq_params['Nfreqs'] % 1, 0):
            raise ValueError("end_freq - start_freq must be evenly divisible by channel_width")
        freq_params['Nfreqs'] = int(freq_params['Nfreqs'])

        freq_arr = np.linspace(freq_params['start_freq'],
                               freq_params['end_freq'] + freq_params['channel_width'],
                               freq_params['Nfreqs'], endpoint=False)

    chan_width = np.atleast_1d(freq_params["channel_width"])
    if chan_width.size < freq_params['Nfreqs'] and chan_width.size == 1:
        chan_width = np.full((freq_params['Nfreqs'],), chan_width, dtype=float)

    return_dict = {}
    return_dict['Nfreqs'] = freq_params['Nfreqs']
    return_dict['freq_array'] = freq_arr
    return_dict['channel_width'] = chan_width
    return_dict['Nspws'] = 1

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
    return_dict = {}

    init_time_params = copy.deepcopy(time_params)
    time_params = copy.deepcopy(time_params)

    time_keywords = ['time_array', 'start_time', 'end_time', 'Ntimes', 'integration_time',
                     'duration_hours', 'duration_days']
    ta, st, et, nt, it, dh, dd = [tk in time_params for tk in time_keywords]
    kws_used = ", ".join(sorted(time_params.keys()))
    daysperhour = 1 / 24.
    hourspersec = 1 / 60. ** 2
    dayspersec = daysperhour * hourspersec

    if ta:
        # Time array is defined. Supersedes all other parameters:
        time_arr = time_params['time_array']
        time_params['Ntimes'] = len(time_arr)
        time_params['start_time'] = np.min(time_arr)
        time_params['integration_time'] = np.mean(np.diff(time_arr))

    else:
        if not (st or et):
            raise ValueError("Start or end time must be specified: " + kws_used)
        if dh and not dd:
            time_params['duration'] = time_params['duration_hours'] * daysperhour
            dd = True
        elif dd:
            time_params['duration'] = time_params['duration_days']

        if not nt:
            if not it:
                raise ValueError("Either integration_time or Ntimes must be "
                                 "included in parameters: " + kws_used)
            if st and et:
                time_params['duration'] = (
                    time_params['end_time']
                    - time_params['start_time']
                    + time_params['integration_time'] * dayspersec)

                dd = True
            if dd:
                time_params['Ntimes'] = int(
                    np.round(
                        time_params['duration'] / (time_params['integration_time'] * dayspersec))
                )
            else:
                raise ValueError("Either duration or time bounds must be specified: " + kws_used)

        if not it:
            if not dd:
                raise ValueError("Either duration or integration time "
                                 "must be specified: " + kws_used)
            # In seconds
            time_params['integration_time'] = \
                time_params['duration'] / dayspersec / float(time_params['Ntimes'])

        inttime_days = time_params['integration_time'] * dayspersec
        if not dd:
            time_params['duration'] = inttime_days * (time_params['Ntimes'])
            dd = True
        if not st:
            if et and dd:
                time_params['start_time'] = (
                    time_params['end_time']
                    - time_params['duration']
                    + inttime_days)

        if not et:
            if st and dd:
                time_params['end_time'] = (
                    time_params['start_time']
                    + time_params['duration']
                    - inttime_days)

        time_arr = np.linspace(time_params['start_time'],
                               time_params['end_time'] + inttime_days,
                               time_params['Ntimes'], endpoint=False)

        if time_params['Ntimes'] != 1:
            if not np.allclose(
                    np.diff(time_arr), inttime_days * np.ones(time_params["Ntimes"] - 1),
                    atol=dayspersec):  # To nearest second
                raise ValueError(
                    "Calculated time array is not consistent with set integration_time."
                    "\nInput parameters are: {}".format(str(init_time_params)))

    return_dict['integration_time'] = time_params['integration_time']
    return_dict['time_array'] = time_arr
    return_dict['Ntimes'] = time_params['Ntimes']
    return_dict['start_time'] = time_params['start_time']

    return return_dict


def freq_array_to_params(freq_array):
    """
    Get a set of parameters that can be used to generate a given frequency array.

    Parameters
    ----------
    freq_array : array of float
        Frequencies in Hz.

    Returns
    -------
    dict
        Dictionary of frequency parameters consistent with freq_array.
        (channel_width, Nfreqs, bandwidth, start_freq, end_freq)
        See pyuvsim documentation for details:
        https://pyuvsim.readthedocs.io/en/latest/parameter_files.html#frequency

    """
    freq_array = np.asarray(freq_array).ravel()

    fdict = {}
    if freq_array.size < 2:
        fdict['channel_width'] = 1.0
        fdict['Nfreqs'] = 1
        fdict['start_freq'] = freq_array.item(0)
        return fdict

    fdict['channel_width'] = np.diff(freq_array).item(0)
    fdict['Nfreqs'] = freq_array.size
    fdict['bandwidth'] = fdict['channel_width'] * fdict['Nfreqs']
    fdict['start_freq'] = freq_array.item(0)
    fdict['end_freq'] = freq_array.item(-1)

    return fdict


def time_array_to_params(time_array):
    """
    Get a set of parameters that can be used to generate a given time array.

    Returns a dictionary of parameters that can be used by parse_time_params
    to obtain the time_array passed in.

    Parameters
    ----------
    time_array : array of float
        Julian dates.

    Returns
    -------
    dict
        Dictionary of time parameters consistent with time_array.
        (integration_time, Ntimes, duration_days, start_time, end_times)
        See pyuvsim documentation for details:
        https://pyuvsim.readthedocs.io/en/latest/parameter_files.html#time

    """
    time_array = np.asarray(time_array)
    Ntimes_uniq = np.unique(time_array).size

    tdict = {}
    if Ntimes_uniq < 2:
        tdict['integration_time'] = 1.0
        tdict['Ntimes'] = 1
        tdict['start_time'] = time_array.item(0)
        return tdict

    dt = np.min(np.diff(time_array)).item()

    if not np.allclose(np.diff(time_array), np.ones(time_array.size - 1) * dt):
        tdict['time_array'] = time_array.tolist()

    tdict['integration_time'] = dt * (24. * 3600.)
    tdict['Ntimes'] = Ntimes_uniq
    tdict['duration_days'] = tdict['integration_time'] * tdict['Ntimes'] / (24. * 3600.)
    tdict['start_time'] = time_array.item(0)
    tdict['end_time'] = time_array.item(-1)

    return tdict


def initialize_uvdata_from_params(obs_params, return_beams=None):
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
        Option to return the beam_list and beam_dict. Currently defaults to True, but
        will default to False starting in version 1.4.

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
    if return_beams is None:
        warnings.warn(
            "The return_beams parameter currently defaults to True, but starting in"
            "version 1.4 it will default to False.",
            DeprecationWarning
        )
        return_beams = True

    uvparam_dict = {}  # Parameters that will go into UVData
    if isinstance(obs_params, str):
        param_dict = _config_str_to_dict(obs_params)  # Container for received settings.
    else:
        param_dict = copy.deepcopy(obs_params)

    # Parse telescope parameters
    tele_dict = param_dict['telescope']
    uvparam_dict.update(tele_dict)
    tele_params, beam_list, beam_dict = parse_telescope_params(tele_dict, config_path=param_dict[
        'config_path'])
    uvparam_dict.update(tele_params)
    uvparam_dict["x_orientation"] = beam_list.x_orientation

    # Use extra_keywords to pass along required paths for file history.
    extra_keywords = {}
    if 'obs_param_file' in param_dict:
        extra_keywords['obsparam'] = param_dict['obs_param_file']
        if (
            'telescope_config_name' in tele_dict
            and isinstance(tele_dict['telescope_config_name'], str)
        ):
            extra_keywords['telecfg'] = tele_dict['telescope_config_name']
        if (
            'array_layout' in tele_dict
            and isinstance(tele_dict['array_layout'], str)
        ):
            extra_keywords['layout'] = tele_dict['array_layout']

    if "world" in tele_params:
        extra_keywords['world'] = tele_params.get('world')
    uvparam_dict['extra_keywords'] = extra_keywords

    # Parse frequency structure
    freq_dict = param_dict['freq']
    uvparam_dict.update(parse_frequency_params(freq_dict))

    # Parse time structure
    time_dict = param_dict['time']
    uvparam_dict.update(parse_time_params(time_dict))

    # There does not seem to be any way to get polarization_array into uvparam_dict, so
    # let's add it explicitly.
    if "polarization_array" in param_dict:
        uvparam_dict['polarization_array'] = param_dict['polarization_array']

    # Parse polarizations
    if uvparam_dict.get('polarization_array', None) is None:
        uvparam_dict['polarization_array'] = np.array([-5, -6, -7, -8])
    if 'Npols' not in uvparam_dict:
        uvparam_dict['Npols'] = len(uvparam_dict['polarization_array'])

    if 'object_name' not in param_dict:
        tloc = EarthLocation.from_geocentric(*uvparam_dict['telescope_location'], unit='m')
        time = Time(uvparam_dict['time_array'][0], scale='utc', format='jd')
        src, _ = create_mock_catalog(time, arrangement='zenith', array_location=tloc)
        if 'sources' in param_dict:
            source_file_name = os.path.basename(param_dict['sources']['catalog'])
            uvparam_dict['object_name'] = '{}_ra{:.4f}_dec{:.4f}'.format(
                source_file_name, src.ra.deg[0], src.dec.deg[0]
            )
        else:
            uvparam_dict['object_name'] = 'Unspecified'
    else:
        uvparam_dict['object_name'] = param_dict['object_name']

    uv_obj = UVData()
    uv_obj._set_future_array_shapes()
    # use the __iter__ function on UVData to get list of UVParameters on UVData
    valid_param_names = [getattr(uv_obj, param).name for param in uv_obj]
    for k in valid_param_names:
        if k in param_dict:
            setattr(uv_obj, k, param_dict[k])
        if k in uvparam_dict:
            setattr(uv_obj, k, uvparam_dict[k])

    bls = np.array(
        [
            uv_obj.antnums_to_baseline(uv_obj.antenna_numbers[j], uv_obj.antenna_numbers[i])
            for i in range(0, uv_obj.Nants_data)
            for j in range(i, uv_obj.Nants_data)
        ]
    )

    uv_obj.baseline_array = np.tile(bls, uv_obj.Ntimes)
    uv_obj.Nbls = bls.size
    uv_obj.time_array = np.repeat(uv_obj.time_array, uv_obj.Nbls)
    uv_obj.integration_time = np.repeat(uv_obj.integration_time, uv_obj.Nbls * uv_obj.Ntimes)
    uv_obj.Nblts = uv_obj.Nbls * uv_obj.Ntimes

    uv_obj.ant_1_array, uv_obj.ant_2_array = uv_obj.baseline_to_antnums(uv_obj.baseline_array)

    # add other required metadata to allow select to work without errors
    # these will all be overwritten in uvsim._complete_uvdata, so it's ok to hardcode them here

    _set_lsts_on_uvdata(uv_obj)
    uv_obj.set_uvws_from_antenna_positions()
    uv_obj.history = ''

    uv_obj._set_drift()
    uv_obj.vis_units = 'Jy'

    uv_obj.instrument = uv_obj.telescope_name

    uv_obj.spw_array = np.array([0])
    if uv_obj.Ntimes == 1:
        uv_obj.integration_time = np.ones_like(uv_obj.time_array, dtype=np.float64)  # Second
    else:
        # Note: currently only support a constant spacing of times
        uv_obj.integration_time = (
            np.ones_like(uv_obj.time_array, dtype=np.float64)
            * np.diff(np.unique(uv_obj.time_array))[0]
            * (24. * 60 ** 2)  # Seconds
        )

    # select on object
    valid_select_keys = [
        'antenna_nums', 'antenna_names', 'ant_str', 'bls',
        'frequencies', 'freq_chans', 'times', 'blt_inds'
    ]

    # downselect baselines (or anything that can be passed to pyuvdata's select method)
    # Note: polarization selection is allowed here, but will cause an error if the incorrect pols
    # are passed to pyuvsim.
    if 'select' in param_dict:
        select_params = param_dict['select']
        no_autos = bool(select_params.pop('no_autos', False))
        select_params = {k: v for k, v in select_params.items() if k in valid_select_keys}
        if 'antenna_nums' in select_params:
            select_params['antenna_nums'] = list(map(int, select_params['antenna_nums']))
        redundant_threshold = param_dict['select'].get('redundant_threshold', None)
        if 'bls' in select_params:
            bls = select_params['bls']
            if isinstance(bls, str):
                # If read from file, this should be a string.
                bls = ast.literal_eval(bls)
                select_params['bls'] = bls
        if len(select_params) > 0:
            uv_obj.select(**select_params)

        if no_autos:
            uv_obj.select(ant_str='cross')

        if redundant_threshold is not None:
            uv_obj.compress_by_redundancy(tol=redundant_threshold)
    # we construct uvdata objects in (time, ant1) order
    # but the simulator will force (time, baseline) later
    # so order this now so we don't get any warnings.
    uv_obj.reorder_blts(order="time", minor_order="baseline")
    uv_obj.check()

    if return_beams:
        return uv_obj, beam_list, beam_dict
    else:
        return uv_obj


def _complete_uvdata(uv_in, inplace=False):
    """
    Fill out all required parameters of a :class:`pyuvdata.UVData` object.

    Ensure that it passes the :func:`pyuvdata.UVData.check()`.
    This will overwrite existing data in `uv_in`!

    Parameters
    ----------
    uv_in : :class:`pyuvdata.UVData` instance
        Usually an incomplete object, containing only metadata.
    inplace : bool, optional
        Whether to perform the filling on the passed object, or a copy.

    Returns
    -------
    :class:`pyuvdata.UVData` : filled/completed object (if `inplace` is `True`, it is
        the modified input). With zeroed data_array, no flags and nsample_array of all ones.

    """
    if not inplace:
        uv_obj = copy.deepcopy(uv_in)
    else:
        uv_obj = uv_in

    # moved to have meta-data compatible objects from the init
    # but retain this functionality for now.
    # in most cases should be a no-op anyway.
    if uv_obj.lst_array is None:
        _set_lsts_on_uvdata(uv_obj)

    # Clear existing data, if any.
    _shape = (uv_obj.Nblts, uv_obj.Nfreqs, uv_obj.Npols)
    uv_obj.data_array = np.zeros(_shape, dtype=complex)
    uv_obj.flag_array = np.zeros(_shape, dtype=bool)
    uv_obj.nsample_array = np.ones(_shape, dtype=float)

    uv_obj.extra_keywords = {}

    uv_obj.check()

    return uv_obj


def initialize_uvdata_from_keywords(
        output_yaml_filename=None, antenna_layout_filepath=None, output_layout_filename=None,
        array_layout=None, telescope_location=None, telescope_name=None, Nfreqs=None,
        start_freq=None, bandwidth=None, freq_array=None, channel_width=None, Ntimes=None,
        integration_time=None, start_time=None, time_array=None, bls=None,
        antenna_nums=None, antenna_names=None, polarization_array=None, no_autos=False,
        redundant_threshold=None, write_files=True, path_out=None, complete=False, **kwargs):
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
    write_files : bool
        If True, write out the parameter information to yaml files.
    path_out : str (optional)
        Path in which to place generated configuration files, if write_files is True.
        Defaults to current directory.
    complete : bool (optional)
        Whether to fill out the :class:`pyuvdata.UVData` object with its requisite
        data arrays, and check if it's all consistent.
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

    if write_files:
        if path_out is None:
            path_out = '.'
        if not outfile:
            if not arrfile:
                output_layout_filename = 'antenna_layout.csv'
            else:
                output_layout_filename = os.path.basename(antenna_layout_filepath)
            outfile = True

        # Increment name appropriately:
        output_layout_filepath = os.path.join(path_out, output_layout_filename)
        output_layout_filename = os.path.basename(
            check_file_exists_and_increment(output_layout_filepath, 'csv')
        )

        if output_yaml_filename is None:
            output_yaml_filename = 'obsparam.yaml'
        output_yaml_filename = check_file_exists_and_increment(
            os.path.join(path_out, output_yaml_filename), 'yaml'
        )

        if antenna_layout_filepath is not None:
            # Copying original file to new place, if it exists
            if os.path.exists(antenna_layout_filepath):
                shutil.copyfile(
                    antenna_layout_filepath,
                    os.path.join(path_out, output_layout_filename)
                )

    antenna_numbers = None
    if isinstance(array_layout, dict):
        antenna_numbers = np.fromiter(array_layout.keys(), dtype=int)
        antpos_enu = array_layout.values()
        if antenna_names is None:
            antenna_names = antenna_numbers.astype('str')
        if write_files:
            _write_layout_csv(
                output_layout_filepath, antpos_enu, antenna_names, antenna_numbers
            )

    if array_layout is None:
        if outfile:
            array_layout = output_layout_filename
        else:
            array_layout = antenna_layout_filepath

    freq_params = {
        'Nfreqs': Nfreqs, 'start_freq': start_freq, 'bandwidth': bandwidth,
        'freq_array': freq_array, 'channel_width': channel_width
    }
    time_params = {
        'Ntimes': Ntimes, 'start_time': start_time,
        'integration_time': integration_time, 'time_array': time_array
    }
    selection_params = {
        'bls': bls, 'redundant_threshold': redundant_threshold,
        'antenna_nums': antenna_nums, 'no_autos': no_autos
    }
    tele_params = {
        'telescope_location': repr(tuple(telescope_location)),
        'telescope_name': telescope_name
    }
    layout_params = {
        'antenna_names': antenna_names, 'antenna_numbers': antenna_numbers,
        'array_layout': array_layout
    }

    freq_params = {k: v for k, v in freq_params.items() if v is not None}
    time_params = {k: v for k, v in time_params.items() if v is not None}
    selection_params = {k: v for k, v in selection_params.items()
                        if v is not None}
    tele_params = {k: v for k, v in tele_params.items() if v is not None}
    layout_params = {k: v for k, v in layout_params.items() if v is not None}

    uv_obj = UVData()

    valid_param_names = [getattr(uv_obj, param).name for param in uv_obj]

    extra_kwds = {k: v for k, v in kwargs.items() if k in valid_param_names}

    # Convert str polarization array to int.
    if polarization_array is not None:
        if type(polarization_array[0]) is not int:
            polarization_array = np.array(uvutils.polstr2num(polarization_array))

    if output_yaml_filename is None:
        output_yaml_filename = ''

    param_dict = {
        'time': time_params,
        'freq': freq_params,
        'select': selection_params,
        'telescope': tele_params,
        'config_path': path_out,
        'polarization_array': polarization_array
    }

    param_dict.update(**extra_kwds)

    if write_files:
        tele_params['array_layout'] = antenna_layout_filepath
        param_dict['telescope'] = tele_params
        writeable_param_dict = copy.deepcopy(param_dict)
        if polarization_array is not None:
            writeable_param_dict['polarization_array'] = polarization_array.tolist()
        with open(output_yaml_filename, 'w') as yfile:
            yaml.dump(writeable_param_dict, yfile, default_flow_style=False)

    param_dict['obs_param_file'] = os.path.basename(output_yaml_filename)
    param_dict['telescope'].update(layout_params)
    uv_obj = initialize_uvdata_from_params(param_dict, return_beams=False)

    if complete:
        _complete_uvdata(uv_obj, inplace=True)

    return uv_obj


def uvdata_to_telescope_config(
        uvdata_in, beam_filepath, layout_csv_name=None, telescope_config_name=None,
        return_names=False, path_out='.'):
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
    if telescope_config_name is None:
        telescope_config_path = \
            check_file_exists_and_increment(
                os.path.join(
                    path_out, 'telescope_config_{}.yaml'.format(uvdata_in.telescope_name)
                )
            )
        telescope_config_name = os.path.basename(telescope_config_path)

    if layout_csv_name is None:
        layout_csv_path = check_file_exists_and_increment(
            os.path.join(path_out, uvdata_in.telescope_name + "_layout.csv"))
        layout_csv_name = os.path.basename(layout_csv_path)

    antpos_enu, antenna_numbers = uvdata_in.get_ENU_antpos()

    _write_layout_csv(
        os.path.join(path_out, layout_csv_name),
        antpos_enu, uvdata_in.antenna_names, uvdata_in.antenna_numbers
    )

    # Write the rest to a yaml file.
    yaml_dict = {
        "telescope_name": uvdata_in.telescope_name,
        "telescope_location": repr(tuple(uvdata_in.telescope_location_lat_lon_alt_degrees)),
        "Nants": uvdata_in.Nants_telescope,
        "beam_paths": {0: beam_filepath}
    }

    with open(os.path.join(path_out, telescope_config_name), 'w+') as yfile:
        yaml.dump(yaml_dict, yfile, default_flow_style=False)

    print('Path: {}, telescope_config: {}, layout: {}'.format(
        path_out, telescope_config_name, layout_csv_name)
    )

    if return_names:
        return path_out, telescope_config_name, layout_csv_name


def uvdata_to_config_file(uvdata_in, param_filename=None, telescope_config_name='',
                          layout_csv_name='', catalog='mock', path_out='.'):
    """
    Extract simulation configuration settings from uvfits.

    When used with :func:`~uvdata_to_telescope_config`, this will produce all the necessary
    configuration yaml and csv file to make an "empty" :class:`pyuvdata.UVData` object comparable to
    the argument `uvdata_in`. The generated file will match `uvdata_in` in frequency, time,
    antenna positions, and uvw coordinates.

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
        param_filename = check_file_exists_and_increment(os.path.join(path_out, 'obsparam.yaml'))
        param_filename = os.path.basename(param_filename)

    freq_array = uvdata_in.freq_array
    time_array = uvdata_in.time_array
    integration_time_array = np.array(uvdata_in.integration_time)
    if np.max(integration_time_array) != np.min(integration_time_array):
        warnings.warn('The integration time is not constant. Using the shortest integration time')

    tdict = time_array_to_params(time_array)
    fdict = freq_array_to_params(freq_array)

    tdict.pop('time_array', None)

    param_dict = {
        "time": tdict,
        "freq": fdict,
        "sources": {"catalog": catalog},
        "telescope": {
            "telescope_config_name": telescope_config_name,
            "array_layout": layout_csv_name
        },
        "filing": {
            "outdir": '.',
            "outfile_name": '',
            "outfile_prefix": ''
        }
    }

    if catalog == 'mock':
        param_dict['sources']['mock_arrangement'] = 'zenith'

    with open(os.path.join(path_out, param_filename), 'w') as yfile:
        yaml.dump(param_dict, yfile, default_flow_style=False)
