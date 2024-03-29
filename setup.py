# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import glob
import io
import sys

from setuptools import setup

# add pyuvsim to our path in order to use the branch_scheme function
sys.path.append("pyuvsim")
from branch_scheme import branch_scheme  # noqa

with io.open("README.md", "r", encoding="utf-8") as readme_file:
    readme = readme_file.read()

sim_reqs = ["mpi4py>=3.1.1"]
moon_reqs = ["lunarsky>=0.2.2"]
healpix_reqs = ["astropy_healpix>=1.0.2"]
casa_reqs = ["python-casacore>=3.5.2"]
test_reqs = (
    sim_reqs
    + moon_reqs
    + casa_reqs
    + healpix_reqs
    + ["coverage", "line_profiler", "pre-commit", "pytest", "pytest-cov"]
)
doc_reqs = ["sphinx", "pypandoc"]

setup_args = {
    "name": "pyuvsim",
    "author": "Radio Astronomy Software Group",
    "url": "https://github.com/RadioAstronomySoftwareGroup/pyuvsim",
    "license": "BSD",
    "description": "A comprehensive simulation package for radio interferometers in python",
    "long_description": readme,
    "long_description_content_type": "text/markdown",
    "package_dir": {"pyuvsim": "pyuvsim"},
    "packages": ["pyuvsim", "pyuvsim.tests"],
    "scripts": glob.glob("scripts/*"),
    "use_scm_version": {"local_scheme": branch_scheme},
    "include_package_data": True,
    "python_requires": ">=3.10",
    "install_requires": [
        "astropy>=6.0",
        "numpy>=1.23",
        "psutil",
        "pyradiosky>=0.2.0",
        "pyuvdata>=2.4.3",
        "pyyaml>=5.4.1",
        "scipy>=1.7.3",
        "setuptools>=61",
        "setuptools_scm>=7.0.3",
    ],
    "extras_require": {
        "sim": sim_reqs,
        "moon": moon_reqs,
        "casa": casa_reqs,
        "healpix": healpix_reqs,
        "all": sim_reqs + moon_reqs + casa_reqs + healpix_reqs,
        "test": test_reqs,
        "doc": doc_reqs,
        "dev": test_reqs + doc_reqs,
    },
    "keywords": "radio astronomy interferometry",
    "classifiers": [
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Astronomy",
    ],
}

if __name__ == "__main__":
    setup(**setup_args)
