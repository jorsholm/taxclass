#!/usr/bin/env bash

python3.9 -m venv /projappl/project_2005718/crest4

source /projappl/project_2005718/crest4/bin/activate

pip3 install crest4==4.3.7

git clone --depth=1 https://github.com/brendanf/crest4_utils.git

cd crest4_utils

git fetch --depth=1 origin e7dcbdde52ff89524b064504ec6a41f02f2b73b9

git checkout FETCH_HEAD

cd ..

