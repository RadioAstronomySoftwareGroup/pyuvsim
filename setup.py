# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

import glob
import io
import json
import os
import sys

from setuptools import setup

sys.path.append("pyuvsim")
import version  # noqa

data = [version.git_origin, version.git_hash, version.git_description, version.git_branch]
with open(os.path.join('pyuvsim', 'GIT_INFO'), 'w') as outfile:
    json.dump(data, outfile)

with io.open('README.md', 'r', encoding='utf-8') as readme_file:
    readme = readme_file.read()

setup_args = {
    'name': 'pyuvsim',
    'author': 'Radio Astronomy Software Group',
    'url': 'https://github.com/RadioAstronomySoftwareGroup/pyuvsim',
    'license': 'BSD',
    'description': 'A comprehensive simulation package for radio interferometers in python',
    'long_description': readme,
    'long_description_content_type': 'text/markdown',
    'package_dir': {'pyuvsim': 'pyuvsim'},
    'packages': ['pyuvsim', 'pyuvsim.tests'],
    'scripts': glob.glob('scripts/*'),
    'version': version.version,
    'include_package_data': True,
    'install_requires': ['numpy>=1.15', 'scipy', 'astropy>=4.0', 'pyyaml', 'pyuvdata'],
    'test_requires': ['pytest'],
    'classifiers': ['Development Status :: 5 - Production/Stable',
                    'Intended Audience :: Science/Research',
                    'License :: OSI Approved :: BSD License',
                    'Programming Language :: Python :: 3.6',
                    'Topic :: Scientific/Engineering :: Astronomy'],
    'keywords': 'radio astronomy interferometry',
    'extras_require': {
        'sim': ['mpi4py>=3.0.0', 'psutil'],
        'all': ['mpi4py>=3.0.0', 'psutil', 'line_profiler', 'h5py'],
        'dev': ['mpi4py>=3.0.0', 'psutil', 'line_profiler', 'h5py', 'pypandoc',
                'pytest', 'pytest-cov', 'sphinx', 'pre-commit']
    }
}

if __name__ == '__main__':
    setup(**setup_args)
