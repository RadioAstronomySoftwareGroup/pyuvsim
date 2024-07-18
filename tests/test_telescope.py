import copy
import os
import shutil
import warnings

import numpy as np
import pytest
from astropy.coordinates import EarthLocation
from pyuvdata import UVBeam

import pyuvsim
from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
from pyuvsim.telescope import BeamConsistencyError

try:
    from pyuvdata.testing import check_warnings
except ImportError:
    # this can be removed once we require pyuvdata >= v3.0
    from pyuvdata.tests import check_warnings

try:
    import mpi4py  # noqa

    has_mpi = True
except ImportError:
    has_mpi = False

herabeam_default = os.path.join(SIM_DATA_PATH, "HERA_NicCST.beamfits")


@pytest.fixture(scope="module")
def beam_objs_main():
    uvb = UVBeam()
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", "The shapes of several attributes will be changing"
        )
        uvb.read_beamfits(herabeam_default)
    if hasattr(uvb, "use_current_array_shapes"):
        uvb.use_future_array_shapes()

    uvb.extra_keywords["beam_path"] = herabeam_default

    # beams are always peak normalized inside BeamList
    uvb.peak_normalize()

    uvb2 = uvb.copy()

    beams = [uvb, uvb2]

    beams.append(pyuvsim.AnalyticBeam("uniform"))
    diameter_m = 14.0
    beams.append(pyuvsim.AnalyticBeam("airy", diameter=diameter_m))
    sigma = 0.03
    # Ignore warnings of pending sigma deprecation
    beams.append(pyuvsim.AnalyticBeam("gaussian", sigma=sigma))
    ref_freq, alpha = 100e6, -0.5
    beams.append(
        pyuvsim.AnalyticBeam(
            "gaussian", sigma=sigma, ref_freq=ref_freq, spectral_index=alpha
        )
    )
    return beams


@pytest.fixture
def beam_objs(beam_objs_main):
    beams_copy = copy.deepcopy(beam_objs_main)

    return beams_copy


def test_convert_loop(beam_objs):
    beams = beam_objs

    if hasattr(beams[0], "_freq_interp_kind"):
        # this can go away when we require pyuvdata version >= 2.4.2
        beams[0].freq_interp_kind = "linear"
        beams[1].freq_interp_kind = "cubic"

        # Should warn about inconsistent params on UVBeams.
        with check_warnings(UserWarning, match="Conflicting settings for"):
            beamlist = pyuvsim.BeamList(beams)

        # Convert beams to strings:
        # Fail, because UVBeams are inconsistent.
        with pytest.raises(ValueError, match="Conflicting settings for"):
            beamlist.set_str_mode()

        beams[1].freq_interp_kind = "linear"
    else:
        beamlist = pyuvsim.BeamList(beams)

    beamlist.set_str_mode()

    # check that _obj_to_str on a string beam works
    beamlist2 = copy.deepcopy(beamlist)
    beamlist2._obj_to_str(beamlist2[0])

    assert beamlist2 == beamlist

    for bs in beamlist:
        assert isinstance(bs, str)
    assert beamlist._obj_beam_list == []

    # Convert strings to beams.
    beamlist.set_obj_mode()
    for bi, b in enumerate(beamlist):
        assert b == beams[bi]

    assert beamlist._str_beam_list == []


@pytest.mark.filterwarnings("ignore:This method will be removed in version 3.0")
def test_force_future_shapes(beam_objs):
    if hasattr(beam_objs[0], "use_current_array_shapes"):
        beam_objs[0].use_current_array_shapes()

        beamlist = pyuvsim.BeamList(beam_objs)

        assert beamlist[0].future_array_shapes


def test_object_mode(beam_objs, tmp_path):
    beams = beam_objs
    beamlist = pyuvsim.BeamList(beams)

    beamfits_file = os.path.join(SIM_DATA_PATH, "HERA_NicCST.beamfits")

    new_beam_file = os.path.join(tmp_path, "HERA_NicCST.uvbeam")
    shutil.copyfile(beamfits_file, new_beam_file)

    uvb = copy.deepcopy(beams[0])
    uvb.extra_keywords["beam_path"] = new_beam_file
    if hasattr(beams[0], "_freq_interp_kind"):
        # this can go away when we require pyuvdata version >= 2.4.2
        beamlist[0].freq_interp_kind = "cubic"
        uvb.freq_interp_kind = "quartic"
        warn_type = UserWarning
        msg = "Conflicting settings for"
    else:
        warn_type = None
        msg = ""
    # Warn if inserted object mismatches.
    with check_warnings(warn_type, match=msg):
        beamlist.append(uvb, uvb_read_kwargs={"file_type": "beamfits"})
    assert len(beamlist) == 7

    if hasattr(beams[0], "_freq_interp_kind"):
        # Error if converting to string mode with mismatched keywords:
        with pytest.raises(ValueError, match="Conflicting settings "):
            beamlist.set_str_mode()
    else:
        # otherwise check that looping str/obj modes works
        beamlist.set_str_mode()
        with check_warnings(
            DeprecationWarning,
            match="This beamfits file does not have a '.fits' or '.beamfits' "
            "extension, so UVBeam does not recognize it as a beamfits file. "
            "Either change the file extension or specify the beam_type. This is "
            "currently handled with a try/except clause in pyuvsim, but this "
            "will become an error in version 1.4",
        ):
            beamlist.set_obj_mode()

    # Error if converting to string mode without beam_paths:
    beamlist[0].extra_keywords.pop("beam_path")
    with pytest.raises(ValueError, match="Need to set "):
        beamlist.set_str_mode()

    # Insert string -- Converts to object
    new_anabeam = "analytic_gaussian_sig=3.0"
    beamlist[-1] = new_anabeam

    assert isinstance(beamlist[-1], pyuvsim.AnalyticBeam)
    assert beamlist[-1].sigma == 3.0


def test_string_mode(beam_objs):
    beams = beam_objs
    beamlist = pyuvsim.BeamList(beams)
    beamlist.set_str_mode()

    uvb = beams[0]
    beamlist.append(uvb)

    assert isinstance(beamlist[-1], str)


def test_comparison(beam_objs):
    beamlist = pyuvsim.BeamList(beam_objs)
    beamlist.set_str_mode()

    beamlist2 = pyuvsim.BeamList(beamlist._str_beam_list)
    assert beamlist == beamlist2

    if not has_mpi:
        with pytest.raises(ImportError, match="You need mpi4py to use shared memory"):
            beamlist.set_obj_mode(use_shared_mem=True)

    beamlist.set_obj_mode()
    beamlist2.set_obj_mode()
    assert beamlist == beamlist2


def test_beamlist_errors(beam_objs):
    # make a copy to enable Telescope equality checking
    beams = copy.deepcopy(beam_objs)
    beamlist = pyuvsim.BeamList(beams, check=False)

    # Try to make a BeamList with a mixture of strings and objects.
    newlist = copy.deepcopy(beamlist._obj_beam_list)
    newlist[2] = beamlist._obj_to_str(newlist[2])
    with pytest.raises(ValueError, match="Invalid beam list:"):
        pyuvsim.BeamList(newlist)

    # Try to append an invalid beam path while in object mode.
    beam_path = "invalid_file.beamfits"
    with pytest.raises(ValueError, match="Invalid file path"):
        beamlist.append(beam_path)

    # test error on beams with different x_orientation
    beams[0].x_orientation = None
    with pytest.raises(
        BeamConsistencyError,
        match="x_orientation of beam 2 is not consistent with beam 1",
    ):
        pyuvsim.BeamList(beams, check=True)

    # test warning on beams with no x_orientation
    beams[0].x_orientation = None
    beams[1].x_orientation = None
    with check_warnings(
        UserWarning,
        match="All polarized beams have x_orientation set to None. This will make it "
        "hard to interpret the polarizations of the simulated visibilities.",
    ):
        assert pyuvsim.BeamList(beams).x_orientation is None

    # Compare Telescopes with beamlists of different lengths

    del beams[0]
    array_location = EarthLocation(lat="-30d43m17.5s", lon="21d25m41.9s", height=1073.0)
    tel0 = pyuvsim.Telescope("tel0", array_location, beams)
    tel1 = pyuvsim.Telescope("tel1", array_location, beam_objs)
    assert tel0 != tel1


def test_beamlist_consistency(beam_objs):
    # Does check, but raises no error
    pyuvsim.BeamList(beam_objs[:2])

    pyuvsim.BeamList(beam_objs[:3])


def test_beamlist_consistency_properties(beam_objs):
    beams = beam_objs
    beamlist = pyuvsim.BeamList(beams[:2])
    assert beamlist.x_orientation == beams[0].x_orientation


def test_beamlist_consistency_stringmode(beam_objs):
    beams = beam_objs
    beamlist = pyuvsim.BeamList(beams[:2])
    beamlist.set_str_mode()
    beamlist.check_consistency(force=True)
    assert beamlist.string_mode
    with check_warnings(
        UserWarning, match="Cannot check consistency of a string-mode BeamList!"
    ):
        beamlist.check_consistency(force=False)

    beamlist = pyuvsim.BeamList(beams[:3], check=False)
    beamlist.set_str_mode()

    beamlist2 = pyuvsim.BeamList(beamlist._str_beam_list, check=True, force_check=True)
    assert beamlist2 == beamlist


def test_beam_basis_type(beam_objs):
    beamlist = pyuvsim.BeamList(beam_objs)

    basis_types = beamlist._get_beam_basis_type()

    assert basis_types == dict.fromkeys(range(len(beam_objs)), "az_za")


@pytest.mark.filterwarnings("ignore:key beam_path in extra_keywords is longer")
def test_beam_basis_type_errors(beam_objs):
    beam_objs[0].pixel_coordinate_system = "orthoslant_zenith"
    beam_objs[0].check()

    beamlist = pyuvsim.BeamList(beam_objs, check=False)
    with pytest.raises(
        ValueError,
        match="pyuvsim currently only supports UVBeams with 'az_za' or "
        "'healpix' pixel coordinate systems.",
    ):
        beamlist.check_consistency()

    beamlist.set_str_mode()
    with pytest.raises(
        ValueError, match="Cannot get beam basis type from a string-mode BeamList."
    ):
        beamlist._get_beam_basis_type()


@pytest.mark.filterwarnings("ignore:key beam_path in extra_keywords is longer")
def test_beam_basis_non_orthogonal_error(beam_objs):
    # test with non-orthogonal basis vectors
    # first construct a beam with non-orthogonal basis vectors
    new_basis_vecs = np.zeros_like(beam_objs[1].basis_vector_array)
    new_basis_vecs[0, 0, :, :] = np.sqrt(0.5)
    new_basis_vecs[0, 1, :, :] = np.sqrt(0.5)
    new_basis_vecs[1, :, :, :] = beam_objs[1].basis_vector_array[1, :, :, :]
    new_data = np.zeros_like(beam_objs[1].data_array)
    # drop all the trailing colons in the slicing below
    new_data[0] = np.sqrt(2) * beam_objs[1].data_array[0]
    new_data[1] = beam_objs[1].data_array[1] - beam_objs[1].data_array[0]
    beam_objs[1].basis_vector_array = new_basis_vecs
    beam_objs[1].data_array = new_data
    beam_objs[1].check()

    beamlist = pyuvsim.BeamList(beam_objs, check=False)
    with pytest.raises(
        ValueError,
        match="pyuvsim currently only supports beams with basis vectors that"
        "are aligned with the azimuth and zenith angle in each pixel."
        "Work is in progress to add other basis vector systems.",
    ):
        beamlist.check_consistency()


def test_empty_beamlist():
    a = pyuvsim.BeamList(check=False)
    assert a.x_orientation is None
    assert a.beam_type is None


@pytest.mark.filterwarnings("ignore:key beam_path in extra_keywords is longer than 8")
@pytest.mark.filterwarnings("ignore:Fixing auto polarization power beams")
def test_powerbeam_consistency(beam_objs):
    newbeams = beam_objs[:2]
    for beam in newbeams:
        beam.efield_to_power()

    beamlist = pyuvsim.BeamList(newbeams)
    beamlist.check_consistency()


@pytest.mark.filterwarnings("ignore:key beam_path in extra_keywords is longer than 8")
def test_check_azza_full_sky(beam_objs):
    beam = beam_objs

    beamlist = pyuvsim.BeamList(beam)
    assert beamlist.check_all_azza_beams_full_sky()

    # Downselect a beam to a small sky area
    za_max = np.deg2rad(10.0)
    za_inds_use = np.nonzero(beam[0].axis2_array <= za_max)[0]
    beam[0].select(axis2_inds=za_inds_use)

    beamlist = pyuvsim.BeamList(beam)
    assert not beamlist.check_all_azza_beams_full_sky()
