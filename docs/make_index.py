# -*- coding: utf-8 -*-

"""
Format the readme.md file into the sphinx index.rst file.

"""
from __future__ import absolute_import, division, print_function

import codecs
import inspect
import os

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

    pypandoc.convert_file(readme_file, 'md')
    readme_text = pypandoc.convert_file(readme_file, 'rst')

    out += readme_text

    out += ('\n\nFurther Documentation\n====================================\n'
            '.. toctree::\n'
            '   :maxdepth: 2\n\n'
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
