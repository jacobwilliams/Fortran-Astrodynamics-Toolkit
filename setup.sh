#!/bin/bash

#
# use pixi to install gawk and gfortran.
# then run the script to get third-party files.
#

rm -rf ./env
mkdir env
cd env
pixi init .
pixi add gawk gfortran
cd ../
pixi run --manifest-path ./env/pixi.toml ./get_third_party.sh