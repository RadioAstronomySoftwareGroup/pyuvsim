# Reference simulations

The goal of a reference simulation is to provide a simulated instrument output for a set of
precisely defined inputs (sky model, antenna positions, primary beam, etc.) to serve as a point
of comparison for other simulators and later versions of pyuvsim. Ideally, these tests should
span multiple axes at once, but computational limitations in continuous regression testing
prohibit this. The version 1 reference simulations span these axes separately (e.g., time,
frequency, baselines, sources), and each serve to test for different expected behavior. We also
implement tests of UVBeam interpolation, integration with HEALPix maps, lunar simulation, and
using multiple analytic beams.

We detail the set of reference simulations and the high-level outcomes of
comparisons with PRISim below, for more detailed writeups see the Memos folder.


## First Reference Simulations
 | Name (beam) | Purpose |
 |:-----|:------|
 |1.1 (short dipole) | Test imaging and source orientation.|
 |1.2 (uniform) | Check that sources behave appropriately near the horizon and rise/set.|
 |1.3 (gaussian) | Check that visibilities have sensible frequency evolution. Get observable fringes and delay transform.|
 |1.4 (airy) | Simulate a realistic sky model with many sources.|
 |1.5 (MWA UVBeam) | Test interpolation of UVBeam.|
 |1.6 (airy) | Test interface with HEALPix.|
 |1.7 (short dipole,uniform,gaussian,airy) | Test use of multiple analytic beams at once.|
 |1.8 (uniform) | Test simulating and imaging on the moon.|


### Details
|              Obsparam File               |              Catalog               | Ntimes  | Nfreqs  |       Layout       |        Beam       |       Results Filename       |
|:----------------------------------------:|:----------------------------------:|:-------:|:-------:|:------------------:|:-----------------:|:----------------------------:|
|     obsparam_ref_1.1_baseline_number.yaml|              RASG.txt              |    1    |    1    | MWA Phase I (128T) |    short dipole   | ref_1.1_baseline_number.uvh5 |
|     obsparam_ref_1.2_time_axis.yaml      | two_points_on_opposite_horizon.txt |  3600   |    1    |  Baseline Lite 4x  |      uniform      |    ref_1.2_time_axis.uvh5    |
|   obsparam_ref_1.3_frequency_axis.yaml   |                R.txt               |    1    |  10000  |    Baseline Lite   | 23° FWHM gaussian |  ref_1.3_frequency_axis.uvh5 |
|     obsparam_ref_1.4_source_axis.yaml    |              gleam.vot             |    1    |    1    |    5 km Triangle   |       airy        |   ref_1.4_source_axis.uvh5   | <!-- TODO: document diameter 14 -->
|       obsparam_ref_1.5_uvbeam.yaml       |                R.txt               |    2    |    2    | MWA Phase I (128T) |     MWA UVBeam    |      ref_1.5_uvbeam.uvh5     |
|      obsparam_ref_1.6_healpix.yaml       |     gsm16_nside128_100mhz.skyh5    |    1    |    1    |    Baseline Lite   |       airy        |     ref_1.6_healpix.uvh5     | <!-- TODO: document diameter 14 -->
|     obsparam_ref_1.7_multi_beam.yaml     |                R.txt               |   100   |   100   |    Baseline Lite   |   All 4 Analytic  |    ref_1.7_multi_beam.uvh5   |
|       obsparam_ref_1.8_lunar.yaml        |              MOON.txt              |    1    |    1    | MWA Phase I (128T) |      uniform      |      ref_1.8_lunar.uvh5      |



The total number of data points is chosen to be sufficient to perform relevant analysis (seen in
the Memo folder) but still lightweight enough to run quickly on a single core
for regression testing in CI using Github Actions.

The catalogs may be found in the **../catalog_files** folder, and the beam/layout files are in the **telescope_config** folder. UVBeam models and large catalog files will need to
be downloaded using `../download_data_files.py` and placed appropriately.

[second_generation](https://github.com/RadioAstronomySoftwareGroup/pyuvsim/tree/main/reference_simulations/second_generation).

For a full description of how antenna layouts, instrument configuration, and catalogs are all
written into parameter files, please see the
[parameter file section of the pyuvsim docs](https://pyuvsim.readthedocs.io/en/latest/parameter_files.html).
As a quick summary: antenna layout is specified by a csv file, overall array configuration and
primary beam assignments are defined by a yaml file, and catalogs are defined either by VOTable
files or csv files. The "obsparam" yaml files define the simulations themselves, including the
catalog, telescope configuration, array layout, and time/frequency array structures, as well as
output filing information and any additional UVData parameters that may be desired.


### Catalogs

RASG.txt:

   - This is a set of point sources near zenith at JD 2460000.0 for an observer at the MWA
   location. They spell out the word "RASG" from <!-- TODO: check the orientation --> east to
   west across the sky with the tops of the letters to the north.

two_points_on_opposite_horizon.txt:

   - This is two points near opposite horizons at JD 2460000.0 for an observer at the MWA location.
   This is chosen so that one source comes on the horizon at 15 minutes, while the other leaves
   the opposite horizon at 30 minutes.

R.txt

   - This is a set of point sources ~30° RA away from zenith at JD 2460000.0 for an observer at
   the MWA location. They spell out the letter "R" from <!-- TODO: check the orientation --> east
   to west across the sky with the tops of the letters to the north.

MOON.txt

   - This is a set of point sources near zenith at JD 2460000.0 for an observer on the moon at
   selenodetic coordinates (0.6875, 24.433, 0). They spell out the word "MOON".

gleam.vot:

   - The GLEAM catalog. It's too large to fit on github, so it's not included in the data
   directory. You can download it using `download_data_files.py`.

gsm16_nside128_100mhz.skyh5:

  -  This is a HEALPix map created using pygdsm GlobalSkyModel16 and scaled down to nside 128
  using healpy, then saved as a skyh5 file using pyradiosky. You can download it using
  `download_data_files.py`.


### Antenna Layouts

128T:

   - This is the full MWA128 Phase I layout with 128 antennas.
   - The layout is in `telescope_config/mwa_128T_layout.csv`.

baseline lite:

   - This consists of a right triangle of antennas with an additional antenna sqrt(2) meters off
   of the center of the hypotenuse. This provides a perfectly N-S and E-W and diagonal baselines,
   as well as some that don't perfectly fit the symmetry.
   - The layout is in `telescope_config/baseline_lite.csv`.

baseline lite 4x:

   - This is just baseline lite but every baseline is increased by 4x for better sensitivity --
   e.g. (50,0) --> (200,0).
   - The layout is in `telescope_config/baseline_lite_4x.csv`.

baseline lite multi beam:

   - This is just baseline lite but each of the 4 antennas has a different analytic beam (short
   dipole,uniform,gaussian,airy).
   - The layout is in `telescope_config/baseline_lite_multi_beam.csv`.

5km triangle:

   - An isosceles triangle consisting of two 5km baselines.
   - The layout is in `telescope_config/5km_triangle_layout.csv`.

See `Memos/` for plots of the array layouts.


### Beams

Four types of AnalyticBeam objects were used as primary beams in the reference simulations: short
dipole, uniform, gaussian, and airy. One UVBeam object was used as a primary beam: the MWA UVBeam
found [here](https://github.com/MWATelescope/mwa_pb).<!-- TODO: discussion on beam behavior -->


### Other design choices

All simulations have the MWA site as their telescope_location for simplicity with the exception
of the lunar simulation. This is lat/lon/alt (-26.70331941, 116.6708152, 377.827). The lunar
simulation has telescope location at selenodetic coordinates (0.6875, 24.433, 0). All simulations
start at Julian Date 2460000.0 (2023-02-24 12:00:00 UTC). Unless specified otherwise, simulations
have a 10 second integration time, a frequency channel width of 100 KHz, and a frequency of 100
MHz. For further simulation specification, see `Memos/` and `list`.


## Old data

Old reference simulation data is stored on the
[Brown Digital Repository (BDR)](https://repository.library.brown.edu/studio/collections/bdr:wte2qah8/)
and ideally updated regularly. The latest first generation reference simulation data can be
downloaded with `../download_ref_sims.py` for comparison, though comparison with latest first
generation reference simulations is also handled with pytest. All reference simulations uploaded
to the BDR are available to download via https -- through script or manual page navigation -- and
the BDR has a solid API.


## Comparison with PRISim

This comparison was performed using pyuvsim v1.0.0 in 2019, and the full comparison writeup can
be found in `Memos/`. The PRISim simulator was designed with the intent of simulating wide field,
high bandwidth interferometers specifically targeting 21 cm instruments. It has mainly been used
in the context of a delay-spectrum style analysis, though images have been made and are known to
be roughly correct. The data for this comparison can be found on the NRAO computing cluster at
`/lustre/aoc/projects/hera/djacobs/prisim_ref/`.

- In the frequency domain, PRISim simulations agree at 10^-5 with pyuvsim, with quantization type
errors at 10^-6.
- The difference shows a floor in the delay spectrum at 1e-6 which is above the expected Blackman
Harris floor at 1e-10
- In the image domain a phase offset at roughly the psf size and can be traced to differences in
uvw calculation.  Antenna positions agree to tolerance which suggests uvw difference is related
to phasing.  PRISim and pyvusim use different codes for phasing and time calculation. The PRISim
phasing code is not covered with analytic tests.
