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
 - name: James Aguirre
   affiliation: 6, 7
 - name: Zachary Martinot
   affiliation: 6, 7
 - name: Nithyanandan Thyagarajan
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
 - name: National Radio Astronomical Observatory
   index: 5
 - name: University of Pennsylvania, Physics Department
   index: 6
 - name: Radio Astronomy Software Group
   index: 7
date: 10 December 2018
bibliography: paper.bib
---

# Summary

Data simulations is essential to the progress of upcoming low-frequency radio telescope arrays such as HERA, MWA, and the SKA. Simulated datasets are used to verify analysis pipelines, provide models for sky-based calibration, and test the effects of design choices and environmental factors. Due to the computational complexity of the interferometer measurement equation, most simulators make simplifying assumptions such as a flat sky, non-polarized response, analytic primary beams, or ignoring Earth precession and nutation or stellar aberration on source positions. pyuvsim is a simulator designed to be accurate and verifiable up front, with strict control over any approximations being made, and including all effects found to be important to 21cm cosmologicaly experiments. The default behavior is to add contributions from each point source above the horizon to each baseline in a full Jones-matrix formalism, with floating precision source positions and interpolated primary beam values. Results are tested against analytic calculations whenever possible, and compared with other simulators for consistency. Data handling and primary beam modeling are done using ``pyuvdata``, which supports a variety of output data types and ensures accurate phasing methods, and source motions and coherency rotations are calculated using ``astropy`` transformations, which takes into account higher order corrections to sky motion.

Currently, pyuvsim supports simulating point sources from the GLEAM catalog and several mock source catalogs. It includes utilities for making parameter files based on ``uvfits`` format data files, allowing users to more easily set up simulations that match an observation. It is parallelized using the Message Passing Interface (MPI). Performance improvements are made by identifying and correcting bottlenecks, which are identified through built-in line profiling tools. Features to be added include support for complicated spectral models for sources, diffuse source models, and ionospheric effects. These are not expected to change the API. Users can rely on pyuvsim for precise results with steady, long-term improvements in resource usage and runtime.
