# -*- mode: python; coding: utf-8 -*
# Copyright (c) 2018 Radio Astronomy Software Group
# Licensed under the 3-clause BSD License

"""Tests for version.py.

"""
from __future__ import absolute_import, division, print_function

import nose.tools as nt
import sys
import os
import six
from six import StringIO
import subprocess
import json
import pyuvsim
from pyuvsim.data import DATA_PATH


def test_get_gitinfo_file():
    pyuvsim_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

    git_file = os.path.join(pyuvsim_dir, 'GIT_INFO')
    if not os.path.exists(git_file):
        # write a file to read in
        temp_git_file = os.path.join(DATA_PATH, 'GIT_INFO')
        version_info = pyuvsim.version.construct_version_info()
        data = [version_info['git_origin'], version_info['git_origin'],
                version_info['git_origin'], version_info['git_origin']]
        with open(temp_git_file, 'w') as outfile:
            json.dump(data, outfile)
        git_file = temp_git_file

    with open(git_file) as data_file:
        data = [pyuvsim.version._unicode_to_str(x) for x in json.loads(data_file.read().strip())]
        git_origin = data[0]
        git_hash = data[1]
        git_description = data[2]
        git_branch = data[3]

    test_file_info = {'git_origin': git_origin, 'git_hash': git_hash,
                      'git_description': git_description, 'git_branch': git_branch}

    if 'temp_git_file' in locals():
        file_info = pyuvsim.version._get_gitinfo_file(git_file=temp_git_file)
        os.remove(temp_git_file)
    else:
        file_info = pyuvsim.version._get_gitinfo_file()

    nt.assert_equal(file_info, test_file_info)


def test_construct_version_info():
    # this test is a bit silly because it uses the nearly the same code as the original,
    # but it will detect accidental changes that could cause problems.
    # It does test that the __version__ attribute is set on pyuvsim.

    # this line is modified from the main implementation since we're in pyuvsim/tests/
    pyuvsim_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

    def get_git_output(args, capture_stderr=False):
        """Get output from Git, ensuring that it is of the ``str`` type,
        not bytes."""

        argv = ['git', '-C', pyuvsim_dir] + args

        if capture_stderr:
            data = subprocess.check_output(argv, stderr=subprocess.STDOUT)
        else:
            data = subprocess.check_output(argv)

        data = data.strip()

        if six.PY2:
            return data
        return data.decode('utf8')

    def unicode_to_str(u):
        if six.PY2:
            return u.encode('utf8')
        return u

    try:
        git_origin = get_git_output(['config', '--get', 'remote.origin.url'], capture_stderr=True)
        git_hash = get_git_output(['rev-parse', 'HEAD'], capture_stderr=True)
        git_description = get_git_output(['describe', '--dirty', '--tag', '--always'])
        git_branch = get_git_output(['rev-parse', '--abbrev-ref', 'HEAD'], capture_stderr=True)
    except subprocess.CalledProcessError:
        try:
            # Check if a GIT_INFO file was created when installing package
            git_file = os.path.join(pyuvsim_dir, 'GIT_INFO')
            with open(git_file) as data_file:
                data = [unicode_to_str(x) for x in json.loads(data_file.read().strip())]
                git_origin = data[0]
                git_hash = data[1]
                git_description = data[2]
                git_branch = data[3]
        except (IOError, OSError):
            git_origin = ''
            git_hash = ''
            git_description = ''
            git_branch = ''
    test_version_info = {'version': pyuvsim.__version__, 'git_origin': git_origin,
                         'git_hash': git_hash, 'git_description': git_description,
                         'git_branch': git_branch}

    nt.assert_equal(pyuvsim.version.construct_version_info(), test_version_info)


def test_main():
    version_info = pyuvsim.version.construct_version_info()

    saved_stdout = sys.stdout
    try:
        out = StringIO()
        sys.stdout = out
        pyuvsim.version.main()
        output = out.getvalue()
        nt.assert_equal(output, 'Version = {v}\ngit origin = {o}\n'
                        'git branch = {b}\ngit description = {d}\n'
                        .format(v=version_info['version'],
                                o=version_info['git_origin'],
                                b=version_info['git_branch'],
                                d=version_info['git_description']))
    finally:
        sys.stdout = saved_stdout
