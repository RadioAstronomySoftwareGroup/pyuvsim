#!/bin/env python

import requests
import os

# The google drive ID for the folder containing these files:
# 14hH-zBhHGddVacc0ncqRWq7ofhGLWfND

version = 'v1'

target_dir = os.path.join('latest_ref_data', version)

if not os.path.exists(target_dir):
    os.makedirs(target_dir)

urlbase = 'https://drive.google.com/uc?export=download'

fileids = 'gdrive_file_ids.dat'

dat = []

with open(fileids, 'r') as dfile:
    for line in dfile:
        line = line.strip()
        # fileid name type size size_unit date time
        dat.append(line.split())

for d in dat:
    fname = d[1]
    fid = d[0]
    r = requests.get(urlbase, params={'id': fid})
    print(fname)
    fname = os.path.join(target_dir, fname)
    with open(fname, 'wb') as ofile:
        ofile.write(r.content)
