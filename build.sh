#!/bin/bash

#
#  Build the Fortran Astrodynamics Toolkit library and example programs
#
#  Requirements:
#  
#  * FoBiS - https://github.com/szaghi/FoBiS
#    Install using: sudo pip install FoBiS.py
#
#  * FORD - https://github.com/cmacmackin/ford
#    Install using: sudo pip install ford
#

DOCDIR='./doc/'                             # build directory for documentation
SRCDIR='./src/'                             # library source directory
TESTDIR='./tests/'                          # tests source directory
BINDIR='./bin/'                             # build directory for example
LIBDIR='./lib/'                             # build directory for library
MODCODE='fortran_astrodynamics_toolkit.f90' # FAT module file name
LIBOUT='libfat.a'                           # name of FAT library
FORDMD='fortran-astrodynamics-toolkit.md'   # FORD MD config file
FCOMPILER='gnu'                             # Fortran compiler flag for FoBiS
FCOMPILERFLAGS='-c -O2'                     # Fortran compiler settings

#build source using FoBiS:

if hash FoBiS.py 2>/dev/null; then
    
    echo "build library..."
    FoBiS.py build -compiler ${FCOMPILER} -cflags "${FCOMPILERFLAGS}" -dbld ${LIBDIR} -s ${SRCDIR} -dmod ./ -dobj ./ -t ${MODCODE} -o ${LIBOUT} -mklib static -colors 
    
    echo "build test programs"
    FoBiS.py build -compiler ${FCOMPILER} -cflags "${FCOMPILERFLAGS}" -dbld ${BINDIR} -s ${TESTDIR} -i ${LIBDIR} -dmod ./ -dobj ./ -libs ${LIBDIR}${LIBOUT} -colors

else
	echo "FoBiS.py not found! Cannot build library. Install using: sudo pip install FoBiS.py"
fi

# build the documentation using FORD:

if hash ford 2>/dev/null; then

	echo "Building documentation..."
    ford ${FORDMD}
        
else
	echo "Ford not found! Cannot build documentation. Install using: sudo pip install ford"
fi
