Fortran-Astrodynamics-Toolkit
=============================

![Image](https://raw.githubusercontent.com/jacobwilliams/Fortran-Astrodynamics-Toolkit/master/tests/pork_chop/pork_chop.png)

A modern Fortran library for astrodynamics.

Overview
---------------

The goal is to produce a comprehensive library, written in modern Fortran (Fortran 2003/2008), of all the standard orbital mechanics algorithms such as:

 * Lambert solvers
 * Kepler propagators
 * ODE solvers (Runge-Kutta, Nystrom, Adams)
 * Orbital element conversions
 * Force models (point mass gravity, geopotential, solar radiation pressure, drag)
 * Reference frames
 * Targeting and optimization
 * Engine models

This is a work in progress, and is currently in a very preliminary state.  

Third-Party Requirements
---------------

To use the ephemeris_module, a copy of one of the JPL binary ephemeris files must be present in the ```eph``` directory.  This can be built from the instructions at: ftp://ssd.jpl.nasa.gov/pub/eph/planets/fortran/userguide.txt.  For example:
```bash
wget ftp://ssd.jpl.nasa.gov/pub/eph/planets/fortran/*
wget ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/de405/*
#edit asc2eph.f file to set NRECL = 4
gfortran asc2eph.f -o asc2eph
cat header.405 ascp*.405 | ./asc2eph
mkdir ~/Fortran-Astrodynamics-Toolkit/eph
mv ./JPLEPH ~/Fortran-Astrodynamics-Toolkit/eph/JPLEPH.405
```

See also
---------------
 * [SPICE](http://naif.jpl.nasa.gov/naif/toolkit.html)
 * [NOVAS](http://aa.usno.navy.mil/software/novas/novas_info.php)
 * [SOFA](http://www.iausofa.org)
