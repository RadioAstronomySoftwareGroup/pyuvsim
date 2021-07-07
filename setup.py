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
    'use_scm_version': {'local_scheme': branch_scheme},
    'include_package_data': True,
    'install_requires': ['numpy>=1.15', 'scipy', 'astropy>=4.0', 'pyyaml',
                         'pyradiosky>=0.1.2', 'pyuvdata>=2.1.3', 'setuptools_scm'],
    'test_requires': ['pytest'],
    'classifiers': ['Development Status :: 5 - Production/Stable',
                    'Intended Audience :: Science/Research',
                    'License :: OSI Approved :: BSD License',
                    'Programming Language :: Python :: 3.6',
                    'Topic :: Scientific/Engineering :: Astronomy'],
    'keywords': 'radio astronomy interferometry',
    'extras_require': {
        'sim': ['mpi4py>=3.0.0', 'psutil'],
        'all': ['mpi4py>=3.0.0', 'psutil', 'line_profiler', 'lunarsky'],
        'moon': ['mpi4py>=3.0.0', 'psutil', 'lunarsky'],
        'dev': ['mpi4py>=3.0.0', 'psutil', 'line_profiler', 'pypandoc',
                'pytest', 'pytest-cov', 'sphinx', 'pre-commit', 'lunarsky']
    }
}

if __name__ == '__main__':
    setup(**setup_args)
