#!/bin/bash

#
#  Build the Fortran Astrodynamics Toolkit library and example programs
#
#  This is using the FoBiS, the Fortran Building System for Poor Men
#   https://github.com/szaghi/FoBiS
#   Install using: pip install FoBiS.py
#
#  Also, build the documentation using ROBODoc
#   Download from: http://rfsber.home.xs4all.nl/Robo/
#

PROJECTNAME='Fortran-Astrodynamics-Toolkit'        # project name for robodoc
DOCDIR='./doc/'                                    # build directory for documentation
SRCDIR='./src/'                                    # library source directory
UNITTESTDIR='./tests/unit_tests/'                  # tests source directory
TESTDIR='./tests/'                                 # tests source directory
BINDIR='./bin/'                                    # build directory for example
LIBDIR='./lib/'                                    # build directory for library
MODCODE='fortran_astrodynamics_toolkit.f90'        # FAT module file name
TESTCODE='test.f90'                                # unit test program file name
EXAMPLECODE1='gravity_test.f90'                    # example program file name
LIBOUT='libfat.a'                                  # name of FAT library

FCOMPILER='gnu'
FCOMPILERFLAGS='-c -O2'
# "-c -O2"
# "-cflags "-c -O2 -Wall -Wtabs"

#build the stand-alone library:
echo "build library..."
FoBiS.py build -compiler ${FCOMPILER} -cflags "${FCOMPILERFLAGS}" -dbld ${LIBDIR} -s ${SRCDIR} -dmod ./ -dobj ./ -t ${MODCODE} -o ${LIBOUT} -mklib static -colors 

echo "build test program : test"
FoBiS.py build -compiler ${FCOMPILER} -cflags "${FCOMPILERFLAGS}" -dbld ${BINDIR} -s ${UNITTESTDIR} -i ${LIBDIR} -dmod ./ -dobj ./ -t ${TESTCODE} -libs ${LIBDIR}${LIBOUT} -o test -colors --log

echo "build test program : gravity"
FoBiS.py build -compiler ${FCOMPILER} -cflags "${FCOMPILERFLAGS}" -dbld ${BINDIR} -s ${TESTDIR}gravity -i ${LIBDIR} -dmod ./ -dobj ./ -t ${EXAMPLECODE1} -libs ${LIBDIR}${LIBOUT} -o gravity -colors --log

echo "build documentation..."
robodoc --rc ./robodoc.rc  --src ${SRCDIR} --doc ${DOCDIR} --multidoc --html --ignore_case_when_linking --syntaxcolors --source_line_numbers --index --tabsize 4 --documenttitle ${PROJECTNAME} --sections

