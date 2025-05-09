#!/bin/bash

# Get some third-party files
# Axel Z. : replace wget with curl to prevent 'Unsupported scheme' error

### Ephemeris files from JPL:

rm -rf ./tmp

mkdir tmp
cd tmp
while read d
do
curl -s ftp://ssd.jpl.nasa.gov/${d}/ | gawk -e '//{print $9}' | while read f
do
  echo "=== $d/$f"
  curl -s -o ./${f} ftp://ssd.jpl.nasa.gov/${d}/${f}
done
done<<EOF
pub/eph/planets/fortran
pub/eph/planets/ascii/de405
pub/eph/planets/ascii/de421
EOF

# edit asc2eph.f file to set NRECL = 4:
# seems that -i only works on mac?
if [[ $OSTYPE == 'darwin'* ]]
then
  sed -i '_original' '/^C.*PARAMETER ( NRECL = 4 )/s/^C//' asc2eph.f
else
  sed --in-place='_original' '/^C.*PARAMETER ( NRECL = 4 )/s/^C//' asc2eph.f
fi

gfortran asc2eph.f -o asc2eph
mkdir ../eph
cat header.405 ascp*.405 | ./asc2eph
mv JPLEPH ../eph/JPLEPH.405
cat header.421 ascp*.421 | ./asc2eph
mv JPLEPH ../eph/JPLEPH.421

### Geopotential file for Earth:

curl -kL -o ./GGM03_Archive.zip https://download.csr.utexas.edu/pub/grace/GGM03/GGM03_Archive.zip
unzip ./GGM03_Archive.zip
mkdir ../grav
cp GGM03_Archive/GGM03C.GEO ../grav