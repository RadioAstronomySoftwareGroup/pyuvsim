from setuptools import setup
import glob
import os.path as op
from os import listdir
import json

# data = [version.git_origin, version.git_hash, version.git_description, version.git_branch]
# with open(op.join('pyuvsim', 'GIT_INFO'), 'w') as outfile:
    # json.dump(data, outfile)

setup_args = {
    'name': 'pyuvsim',
    'author': '#datasim',
    'url': 'https://github.com/HERA-Team/pyuvsim',
    'license': 'BSD',
    'description': 'A radio interferometer simulator',
    'package_dir': {'pyuvsim': 'pyuvsim'},
    'packages': ['pyuvsim', 'pyuvsim.tests'],
    'scripts': glob.glob('scripts/*'),
    'version': '0.0.1',
    'include_package_data': True,
    'install_requires': ['numpy>=1.10', 'astropy>=1.2',  'pyuvdata'],
    'classifiers': ['Development Status :: 2 - Pre-Alpha',
                    'Intended Audience :: Science/Research',
                    'License :: OSI Approved :: BSD License',
                    'Programming Language :: Python :: 2.7',
                    'Topic :: Scientific/Engineering :: Astronomy'],
    'keywords': 'radio astronomy interferometry'
}

if __name__ == '__main__':
    apply(setup, (), setup_args)
