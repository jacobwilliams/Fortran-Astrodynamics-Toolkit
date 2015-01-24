#!/bin/bash

#
#  Build the Fortran Astrodynamics Toolkit unit test program
#
#  This is using the FoBiS, the Fortran Building System for Poor Men
#   https://github.com/szaghi/FoBiS
#
#  Also, build the documentation using ROBODoc
#   http://rfsber.home.xs4all.nl/Robo/
#

#
#  TO DO: modify this to build the FAT library.
#

echo ""
echo "build test program..."
echo ""

./FoBiS.py build -s ./src -compiler gnu  -cflags "-c -O2" -o ./bin/test

echo ""
echo "build documentation..."
echo ""

robodoc --rc ./robodoc.rc  --src ./src --doc ./doc --multidoc --html --ignore_case_when_linking --syntaxcolors --source_line_numbers --index --tabsize 4 --documenttitle Fortran-Astrodynamics-Toolkit --sections
