#! /bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $DIR/..

python setup.py install

cd pyuvsim/tests
nosetests --nologcapture --with-coverage --cover-erase --cover-package=pyuvsim --cover-html "$@"
