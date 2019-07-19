#! /bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $DIR/..

cd pyuvsim/tests
python -m pytest --cov=pyuvsim --cov-config=../../.coveragerc\
       --cov-report term --cov-report html:cover \
       "$@"
