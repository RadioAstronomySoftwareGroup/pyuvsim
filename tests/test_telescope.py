import copy
import os

import numpy as np
import pytest
from astropy.coordinates import EarthLocation
from pyuvdata import AiryBeam, GaussianBeam, ShortDipoleBeam, UniformBeam, UVBeam
from pyuvdata.testing import check_warnings

from pyuvsim.data import DATA_PATH as SIM_DATA_PATH
from pyuvsim.telescope import BeamConsistencyError, BeamList, Telescope

try:
    import mpi4py  # noqa

    has_mpi = True
except ImportError:
    has_mpi = False

herabeam_default = os.path.join(SIM_DATA_PATH, "HERA_NicCST.beamfits")


@pytest.fixture(scope="module")
def beam_objs_main():
    uvb = UVBeam()
    uvb.read_beamfits(herabeam_default)

    uvb2 = uvb.copy()

    beams = [uvb, uvb2]

    beams.append(UniformBeam())
    diameter_m = 14.0
    beams.append(AiryBeam(diameter=diameter_m))
    sigma = 0.03
    beams.append(GaussianBeam(sigma=sigma))
    reference_frequency = 100e6
    alpha = -0.5
    beams.append(
        GaussianBeam(
            sigma=sigma, reference_frequency=reference_frequency, spectral_index=alpha
        )
    )
    beams.append(ShortDipoleBeam())

    return beams


@pytest.fixture
def beam_objs(beam_objs_main):
    beams_copy = copy.deepcopy(beam_objs_main)

    return beams_copy


def test_comparison(beam_objs):
    beamlist = BeamList(beam_objs)

    beamlist2 = BeamList([bi.beam for bi in beamlist.beam_list])
    assert beamlist == beamlist2


def test_en_feed_beamlist_no_errors(beam_objs):
    beams = beam_objs

    beams[0].feed_array = ["e", "n"]
    beams[0].check()

    assert "n" in beams[0].feed_array
    BeamList([beams[0]])


def test_beamlist_errors(beam_objs):
    # make a copy to enable Telescope equality checking
    beams = copy.deepcopy(beam_objs)
    beam_list_init = BeamList(beam_objs)

    # test error on beams with different x_orientation
    beams[0].x_orientation = None
    with pytest.raises(
        BeamConsistencyError,
        match="x_orientation of beam 2 is not consistent with beam 1",
    ):
        BeamList(beams)

    # test warning on beams with no x_orientation
    beams[0].x_orientation = None
    beams[1].x_orientation = None
    with check_warnings(
        UserWarning,
        match="All polarized beams have x_orientation set to None. This will make it "
        "hard to interpret the polarizations of the simulated visibilities.",
    ):
        assert BeamList(beams[:2]).x_orientation is None

    # Compare Telescopes with beamlists of different lengths
    beam_list_init = BeamList(beam_objs)
    beams = copy.deepcopy(beam_objs)
    del beams[0]
    beam_list = BeamList(beams)
    array_location = EarthLocation(lat="-30d43m17.5s", lon="21d25m41.9s", height=1073.0)
    tel0 = Telescope("tel0", array_location, beam_list_init)
    tel1 = Telescope("tel1", array_location, beam_list)
    assert tel0 != tel1


def test_beamlist_consistency(beam_objs):
    # Does check, but raises no error
    BeamList(beam_objs[:2])

    BeamList(beam_objs[:3])


def test_beamlist_consistency_properties(beam_objs):
    beams = beam_objs
    beamlist = BeamList(beams[:2])
    assert beamlist.x_orientation == beams[0].x_orientation


def test_beam_basis_type(beam_objs):
    beamlist = BeamList(beam_objs)

    basis_types = beamlist._get_beam_basis_type()

    assert basis_types == dict.fromkeys(range(len(beam_objs)), "az_za")


@pytest.mark.filterwarnings("ignore:key beam_path in extra_keywords is longer")
def test_beam_basis_type_errors(beam_objs):
    beam_objs[0].pixel_coordinate_system = "orthoslant_zenith"
    beam_objs[0].check()

    with pytest.raises(
        ValueError,
        match="pyuvsim currently only supports UVBeams with 'az_za' or "
        "'healpix' pixel coordinate systems.",
    ):
        BeamList(beam_objs)


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

    with pytest.raises(
        ValueError,
        match="pyuvsim currently only supports beams with basis vectors that"
        "are aligned with the azimuth and zenith angle in each pixel."
        "Work is in progress to add other basis vector systems.",
    ):
        BeamList(beam_objs)


def test_empty_beamlist():
    blist = BeamList([])
    assert blist.x_orientation is None
    assert blist.beam_type == "efield"
    assert blist.data_normalization is None


@pytest.mark.filterwarnings("ignore:key beam_path in extra_keywords is longer than 8")
@pytest.mark.filterwarnings("ignore:Fixing auto polarization power beams")
def test_powerbeam_consistency(beam_objs):
    newbeams = copy.deepcopy(beam_objs)
    newbeams = newbeams[:2]
    for beam in newbeams:
        beam.efield_to_power()

    BeamList(newbeams, beam_type="power")


@pytest.mark.filterwarnings("ignore:key beam_path in extra_keywords is longer than 8")
def test_check_azza_full_sky(beam_objs):
    beams = beam_objs

    beamlist = BeamList(beams)
    assert beamlist.check_all_azza_beams_full_sky()

    # Downselect a beam to a small sky area
    za_max = np.deg2rad(10.0)
    za_inds_use = np.nonzero(beams[0].axis2_array <= za_max)[0]
    beams[0].select(axis2_inds=za_inds_use)

    beamlist = BeamList(beams)
    assert not beamlist.check_all_azza_beams_full_sky()


def test_telescope_init_errors(beam_objs, hera_loc):
    array_xyz = np.array([hera_loc.x.value, hera_loc.y.value, hera_loc.z.value])

    beam_list = BeamList(beam_objs)
    with pytest.raises(
        ValueError,
        match="location must be an EarthLocation object or a "
        "lunarsky.MoonLocation object",
    ):
        Telescope("telescope_name", array_xyz, beam_list)

    with pytest.raises(ValueError, match="beam_list must be a BeamList object"):
        Telescope("telescope_name", hera_loc, beam_objs)


# @pytest.mark.parallel(2, timeout=5)
@pytest.mark.skipif(not has_mpi, reason="Must have mpi for this test")
def test_share_beams(beam_objs):
    from pyuvsim import mpi

    mpi.start_mpi(block_nonroot_stdout=False)

    expected_list = BeamList(beam_objs)

    if mpi.rank == 0:
        beam_list = BeamList(beam_objs)
    else:
        beam_list = BeamList([])

    assert isinstance(beam_list, BeamList)
    mpi.world_comm.Barrier()

    beam_list.share()

    assert len(expected_list) == len(beam_list)

    for expected, other in zip(expected_list, beam_list, strict=True):
        try:
            assert expected == other
        except Exception as err:
            try:
                assert expected.beam == other.beam
            except Exception as err2:
                raise err2 from err
            raise err from None
