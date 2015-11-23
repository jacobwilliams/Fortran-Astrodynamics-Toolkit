Fortran-Astrodynamics-Toolkit
=============================

![Image](https://raw.githubusercontent.com/jacobwilliams/Fortran-Astrodynamics-Toolkit/master/tests/pork_chop/pork_chop.png)

## Overview

The goal is to produce a comprehensive library, written in modern Fortran (Fortran 2003/2008), of all the standard orbital mechanics algorithms.  This is a work in progress, and is currently in a preliminary state.  Currently-implemented and proposed capabilities include:

 * Lambert solvers
  - [x] Gooding
  - [x] Izzo
 * Kepler propagators
  - [x] Gooding
  - [ ] Goodyear
  - [ ] Shepperd
 * ODE solvers
  - [x] Runge-Kutta
  - [ ] Nystrom
  - [ ] Adams
 * Force models
  - [ ] point mass gravity
  - [x] geopotential
  - [ ] solar radiation pressure
  - [ ] drag
 * Reference frames
  - [x] IAU_EARTH
 * Alternate equations of motion
  - [x] Circular restricted three-body problem
  - [x] Clohessy-Wiltshire
 * Misc
  - [ ] orbital element conversions
  - [ ] targeting and optimization
  - [ ] spacecraft engine models

## Examples

Note that the example plots require the [pyplot-fortran](https://github.com/jacobwilliams/pyplot-fortran) module.

<a href="url"><img src="https://github.com/jacobwilliams/Fortran-Astrodynamics-Toolkit/blob/master/tests/crtbp/crtbp_test.png" align="center" height="300"></a><a href="url"><img src="https://github.com/jacobwilliams/Fortran-Astrodynamics-Toolkit/blob/master/tests/gravity/trajectory.png" align="center" height="400"></a>

## Third-Party Requirements

To use the ephemeris_module, a copy of one of the JPL binary ephemeris files must be present in the ```eph``` directory.  This can be built from the instructions at: ftp://ssd.jpl.nasa.gov/pub/eph/planets/fortran/userguide.txt.  For example (on Linux):
```bash
wget ftp://ssd.jpl.nasa.gov/pub/eph/planets/fortran/*
wget ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/de405/*
#edit asc2eph.f file to set NRECL = 4:
sed -i '_original' '/^C.*PARAMETER ( NRECL = 4 )/s/^C//' asc2eph.f
gfortran asc2eph.f -o asc2eph
cat header.405 ascp*.405 | ./asc2eph
mkdir Fortran-Astrodynamics-Toolkit/eph
mv JPLEPH Fortran-Astrodynamics-Toolkit/eph/JPLEPH.405
```
Note that some of the examples require the file to be named ```JPLEPH_2000-2100.405```.

To use the geopotential_module, you need a geopotential model file (for example ```GGM03C.GEO``` from ftp://ftp.csr.utexas.edu/pub/grace/GGM03/GGM03_Archive.zip). This should be placed in the ```grav``` directory.  For example:
```bash
wget ftp://ftp.csr.utexas.edu/pub/grace/GGM03/GGM03_Archive.zip
unzip GGM03_Archive.zip
mkdir Fortran-Astrodynamics-Toolkit/grav
cp GGM03_Archive/GGM03C.GEO Fortran-Astrodynamics-Toolkit/grav
```

## See also

 * [SPICE](http://naif.jpl.nasa.gov/naif/toolkit.html)
 * [NOVAS](http://aa.usno.navy.mil/software/novas/novas_info.php)
 * [SOFA](http://www.iausofa.org)
