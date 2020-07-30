# Reference simulations

The goal of a reference simulation is to provide a simulated instrument output for a set of precisely defined inputs(sky model, antenna positions, primary beam, etc.) to serve as a point of comparison for other simulators and later versions of pyuvsim. Ideally, these tests should span multiple axes at once, but computational limitations in the earlier versions of pyuvsim prohibited this. Therefore, the version 1 reference simulations span these axes separately(e.g., time, frequency, baselines, sources), and each serve to test for different expected behavior.

We detail the set of reference simulations below.


# Second Generation Reference Simulations
 | Name(beam) | Purpose |
 |: ----- | : ------|
 |2.1 (airy) | Replace the first generation reference simulations by covering multiple axes.|
 |2.2 (UVBeam) | Check that visibilities are sensible with UVBeam.|
 |2.3 (airy | Check that visibilities are sensible with a known power spectrum diffuse map. |
 |2.3 (airy) | Check that visibilities are sensible with the healpix map. |
 |2.4 | Check phasing polarized response(not done yet). |


# Details
|         Obsparam File | Catalog | Ntimes | Nfreqs | Layout | Beam | Results Filename |
|: -----------------------------: | : ------------------------------------------: | : ------: | : ------: | : -------------: | : -----------------: | : ----------------------: |
|     obsparam_ref_2.1_airy.yaml | gleam.vot | 60 | 128 | MWA_nocore | airy | ref_2.1_airy.uvh5 |
|  obsparam_ref_2.2_uvbeam.yaml | gleam.vot | 1 | 128 | MWA_nocore | HERA_UVBeam | ref_2.2_uvbeam_gleam.uvh5 |
| obsparam_ref_2.3_airy_flatPS.yaml | healpix_flat_EoR_spectrum_noise.hdf5 | 1 | 128 | MWA_nocore | airy | ref_2.3_airy_flatPS.uvh5 |
|  obsparam_ref_2.3_airy_gsm.yaml | healpix_gsm_shell.hdf5 | 1 | 128 | MWA_nocore | airy | ref_2.3_airy_gsm.uvh5

 |
The total number of data points is constrained by the current level of optimization and availability of computing resources. pyuvsim is currently running on the Oscar cluster at Brown University, given limitations of memory and processor availability.

The catalogs may be found in the ** catalog_files ** folder, and the beam / layout files are in the ** telescope_config ** folder or on Google Drive if they are too large to fit on GitHub.

For a full description of how antenna layouts, instrument configuration, and catalogs are all written into parameter files, please see the documentation at https: // pyuvsim.readthedocs.io / en / latest / parameter_files.html. As a quick summary: antenna layout is specified by a csv file, overall array configuration and primary beam assignments are defined by a yaml file, and catalogs are defined either by VOTable files or csv files. The "obsparam" yaml files define the simulations themselves, including the catalog, telescope configuration, array layout, and time / frequency array structures, as well as output filing information and any additional UVData parameters that may be desired.

# Catalogs

gleam.vot

   - The GLEAM catalog. It's too large to fit on github, so it's not included in the data directory, but it's on Google Drive.

healpix_flat_EoR_spectrum_noise.hdf5

   - The healpix map with a white power spectrum. It's too large to fit on github, so it's not included in the data directory, but it's on Google Drive.

healpix_gsm_shell.hdf5

   - The GSM healpix map. It's too large to fit on github, so it's not included in the data directory, but it's on Google Drive.


# Antenna Layouts


MWA_nocore:

   - This is the MWA128 Phase I layout with the core 40 antennas removed (88 antennas remain).
   - The layout is written in mwa_nocore_layout.csv and the telescope configuration is in mwa88_nocore_config.yaml. The configuration specifies that all antennas have the unphysical "uniform" beam.

![mwa88_layout.png](.. / first_generation / Memos / figures / mwa88_layout.png "MWA-88 layout")


# Beams

Only two types of primary beams were used in these simulations: The airy beam, which is an AnalyticBeam object, and HERA beam, which is a pyuvdata UVBeam object.


# Other design choices

All simulations chose the HERA site as their telescope_location, for simplicity. This is lat / lon / alt(-30.72153°, 21.42831°, 1073.0 m). This includes the MWA128 - like array(MWA_nocore).
