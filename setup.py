# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

from setuptools import setup
import glob
import os
import io
import json
from pyuvsim import version

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
    'test_suite': 'nose.collector',
    'test_requires': ['nose', 'h5py'],
    'install_requires': ['numpy>=1.15', 'astropy>=2.0', 'scipy>1.0.1', 'mpi4py>=3.0.0',
                         'pyyaml>=5.1', 'six>=1.11', 'pyuvdata>=1.3.7'],
    'classifiers': ['Development Status :: 3 - Alpha',
                    'Intended Audience :: Science/Research',
                    'License :: OSI Approved :: BSD License',
                    'Programming Language :: Python :: 2.7',
                    'Programming Language :: Python :: 3.6',
                    'Topic :: Scientific/Engineering :: Astronomy'],
    'keywords': 'radio astronomy interferometry'
}

if __name__ == '__main__':
    setup(**setup_args)
