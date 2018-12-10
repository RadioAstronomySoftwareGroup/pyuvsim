---
title: 'pyuvsim: A comprehensive simulation package for radio interferometers in python.'
tags:
 - radio astronomy
 - simulation
 - pyuvdata

authors:
 - name: Adam E. Lanman
   orcid: 0000-0003-2116-3573
   affiliation: 4
 - name: Bryna J. Hazelton
   orcid: 0000-0001-7532-645X
   affiliation: 1, 2
 - name: Daniel C. Jacobs
   orcid: 0000-0002-0917-2269
   affiliation: 3
 - name: Matthew Kolopanis
   affiliation: 3
 - name: Jonathan C. Pober
   orcid: 0000-0002-3492-0433
   affiliation: 4
 - name: James Aguirre
   affiliation: 6
 - name: Zachary Martinot
   affiliation: 6
 - name: Nithyanandan Thyagarajan
   affiliation: 5
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
date: 10 December 2018
bibliography: paper.bib
---

# Summary

pyuvsim is a precise and robust point source simulation tool for radio interferometry. Simulating data from radio interferometers is a computationally demanding task, but is essential to the progress of upcoming low-frequency experiments such as HERA, MWA, and the SKA. Simulated datasets are used to verify analysis pipelines, provide models for sky-based calibration, and test the effects of design choices and environmental factors. Existing simulators typically make simplifying assumptions to reduce the complexity of the task, and so it can be nearly impossible to determine if an unusual data artifact is a result of an approximation or an actual physical phenomenon. pyuvsim was designed to be accurate, massively-parallelized, and verifiable, with calculations tested against analytic results whenever possible and compared with other simulators for consistency. Data handling and primary beam modeling are done using ``pyuvdata``, which supports a variety of output data types and ensures accurate phasing methods.

Currently, pyuvsim supports simulating point sources from the GLEAM catalog and several mock source catalogs.
