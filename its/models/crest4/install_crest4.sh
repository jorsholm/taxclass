#!/usr/bin/env bash

python3.9 -m venv /projappl/project_2005718/crest4

source /projappl/project_2005718/crest4/bin/activate

pip3 install crest4==4.3.7

git clone --depth=1 https://github.com/xapple/crest4_utils.git

cd crest4_utils

git fetch --depth=1 origin 89d6f909b4c3a41909bc6b39dbff9b33f6876886

git checkout FETCH_HEAD

cd ..

