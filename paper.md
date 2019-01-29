---
title: 'pyuvsim: A comprehensive simulation package for radio interferometers in python.'
tags:
 - radio astronomy
 - simulation
 - pyuvdata

authors:
 - name: Adam E. Lanman
   orcid: 0000-0003-2116-3573
   affiliation: 4, 7
 - name: Bryna J. Hazelton
   orcid: 0000-0001-7532-645X
   affiliation: 1, 2, 7
 - name: Daniel C. Jacobs
   orcid: 0000-0002-0917-2269
   affiliation: 3, 7
 - name: Matthew Kolopanis
   orcid: 0000-0002-2950-2974
   affiliation: 3, 7
 - name: Jonathan C. Pober
   orcid: 0000-0002-3492-0433
   affiliation: 4, 7
 - name: James E. Aguirre
   orcid: 0000-0002-4810-666X
   affiliation: 6, 7
 - name: Nithyanandan Thyagarajan
   orcid: 0000-0003-1602-7868
   affiliation: 5, 7
affiliations:
 - name: University of Washington, eScience Institute
   index: 1
 - name: University of Washington, Physics Department
   index: 2
 - name: Arizona State University, School of Earth and Space Exploration
   index: 3
 - name: Brown University, Physics Department
   index: 4
 - name: National Radio Astronomy Observatory
   index: 5
 - name: University of Pennsylvania, Physics Department
   index: 6
 - name: Radio Astronomy Software Group
   index: 7
date: 29 December 2019
bibliography: paper.bib
---

# Summary

Data simulations are essential to the progress of upcoming low-frequency radio telescope arrays such as HERA, MWA, and the SKA. Simulated datasets are used to verify analysis pipelines, to provide models for sky-based calibration, and to test the effects of design choices and environmental factors. Most simulators make simplifying assumptions to reduce the computational demand of evaluating the measurement equation, compromising accuracy for speed. This can lead to unexpected numerical effects, which can be hard to distinguish from real effects seen in data. ``pyuvsim`` is a simulator designed to be accurate and verifiable up front, with strict control over any approximations being made and including all effects found to be important to 21cm cosmological experiments. The default behavior is to add contributions from each point source above the horizon to each baseline in a full Jones-matrix formalism, with floating precision source positions and interpolated E-field primary beam values. Results are tested against analytic calculations whenever possible, and compared with other simulators and data for consistency, including PRISim [@nithyanandan_2019_2548117] and FHD [@fhd]. Data handling and primary beam modeling are done using the ``pyuvdata`` package [@j_hazelton_pyuvdata:_2017], which supports a variety of output data types and ensures accurate phasing methods. Source motions and coherency rotations are calculated using ``astropy`` transformations [@astropy:2013], which takes into account higher order corrections to sky motion.

Currently, ``pyuvsim`` supports simulating point sources from the GLEAM catalog [@hurley-walker_galactic_2017] and several mock source catalogs. It includes utilities for making parameter files based on the structure of any data file readable by ``pyuvdata``, allowing users to quickly set up simulations to compare with an observation. It is parallelized using the Message Passing Interface (MPI) [@dalcin2005mpi]. Performance improvements are made by identifying and correcting bottlenecks using built-in line profiling tools. Features to be added include support for spectral models for sources, diffuse source models, and ionospheric effects. These are not expected to change the API. Users can rely on ``pyuvsim`` for precise results with steady, long-term improvements in resource usage and runtime.

# Acknowledgements

This work was supported in part by NASA Award #80NSSC18K0389.

