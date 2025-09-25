#!/usr/bin/env bash
set -e -x

# download protax training code
curl -RLO https://github.com/psomervuo/PROTAX/raw/refs/heads/main/protax.zip
unzip protax.zip

# I tried to compile USEARCH but the executable gave segfault
# just download precompiled binary
INSTALL_PREFIX="/projappl/project_2005718"
PATH="$INSTALL_PREFIX/bin:$PATH"

USEARCH_VERSION="12.0-beta1"
wget -O usearch https://github.com/rcedgar/usearch12/releases/download/v$USEARCH_VERSION/usearch_linux_x86_${USEARCH_VERSION%1}
chmod +x usearch
cp usearch $INSTALL_PREFIX/bin/


