# -*- coding: utf-8 -*-

"""
Format the readme.md file into the sphinx index.rst file.

"""
from __future__ import absolute_import, division, print_function

import codecs
import os
import inspect
import re
import pypandoc
from astropy.time import Time


def write_index_rst(readme_file=None, write_file=None):
    t = Time.now()
    t.out_subfmt = 'date'
    out = ('.. pyuvsim documentation master file, created by\n'
           '   make_index.py on {date}\n\n').format(date=t.iso)

    print(readme_file)
    if readme_file is None:
        main_path = os.path.dirname(os.path.dirname(os.path.abspath(inspect.stack()[0][1])))
        readme_file = os.path.join(main_path, 'README.md')

    readme_md = pypandoc.convert_file(readme_file, 'md')

    travis_str = 'https://travis-ci.org/RadioAstronomySoftwareGroup/pyuvsim.svg'
    regex_travis = re.compile(travis_str)
    loc_travis_start = re.search(regex_travis, readme_md).start()
    loc_travis_end = re.search(regex_travis, readme_md).end()
    end_branch_str = r'\)\]'
    regex_end = re.compile(end_branch_str)
    loc_branch_end = re.search(regex_end, readme_md).start()
    branch_str = readme_md[loc_travis_end:loc_branch_end]

    cover_str = 'https://coveralls.io/repos/github/RadioAstronomySoftwareGroup/pyuvsim/badge.svg'
    regex_cover = re.compile(cover_str)
    loc_cover_start = re.search(regex_cover, readme_md).start()
    loc_cover_end = re.search(regex_cover, readme_md).end()

    readme_text = pypandoc.convert_file(readme_file, 'rst')

    out += readme_text

    out += ('\n\nFurther Documentation\n====================================\n'
            '.. toctree::\n'
            '   :maxdepth: 2\n\n'
            '   index\n'
            '   comparison\n'
            '   usage\n'
            '   parameter_files\n'
            '   classes\n')

    out.replace(u"\u2018", "'").replace(u"\u2019", "'").replace(u"\xa0", " ")

    if write_file is None:
        write_path = os.path.dirname(os.path.abspath(inspect.stack()[0][1]))
        write_file = os.path.join(write_path, 'index.rst')
    F = codecs.open(write_file, 'w', 'utf-8')
    F.write(out)
    print("wrote " + write_file)
