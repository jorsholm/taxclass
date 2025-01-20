#!/usr/bin/env bash

curl -JLO "https://sourceforge.net/projects/rdp-classifier/files/latest/download"
unzip rdp_classifier_2.14.zip "rdp_classifier_2.14/dist/*"
mv rdp_classifier_2.14/dist dist
rm -rf rdp_classifier_2.14