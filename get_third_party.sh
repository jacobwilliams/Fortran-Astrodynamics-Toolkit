#!/bin/bash

# Get some third-party files

### Ephemeris files from JPL:

mkdir tmp
cd tmp
wget ftp://ssd.jpl.nasa.gov/pub/eph/planets/fortran/*
wget ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/de405/*
wget ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/de421/*
#edit asc2eph.f file to set NRECL = 4:
sed -i '_original' '/^C.*PARAMETER ( NRECL = 4 )/s/^C//' asc2eph.f
gfortran asc2eph.f -o asc2eph
mkdir ../eph
cat header.405 ascp*.405 | ./asc2eph
mv JPLEPH ../eph/JPLEPH.405
cat header.421 ascp*.421 | ./asc2eph
mv JPLEPH ../eph/JPLEPH.421

### Geopotential file for Earth:

wget ftp://ftp.csr.utexas.edu/pub/grace/GGM03/GGM03_Archive.zip
unzip GGM03_Archive.zip
mkdir ../grav
cp GGM03_Archive/GGM03C.GEO ../grav
