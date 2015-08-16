project: fortran-astrodynamics-toolkit
project_dir: ./src
output_dir: ./doc
project_github: https://github.com/jacobwilliams/Fortran-Astrodynamics-Toolkit
summary: A modern Fortran library for astrodynamics
author: Jacob Williams
github: https://github.com/jacobwilliams
predocmark_alt: >
predocmark: <
docmark_alt:
docmark: !
display: public
display: private
display: protected
display: none
source: true
exclude: brent_module.f90
         rk_module.f90

Brief description
---------------

A modern Fortran library for astrodynamics.

The goal is to produce a comprehensive library, written in modern Fortran (Fortran 2003/2008), 
of all the standard orbital mechanics algorithms such as:

 * Lambert solvers
 * Kepler propagators
 * ODE solvers (Runge-Kutta, Nystrom, Adams)
 * Orbital element conversions
 * Force models (point mass gravity, geopotential, solar radiation pressure, drag)
 * Reference frames
 * Targeting and optimization
 * Engine models

This is a work in progress, and is currently in a very preliminary state.  


License
---------------

The Fortran Astrodynamics Toolkit is released under a [BSD License](https://raw.githubusercontent.com/jacobwilliams/Fortran-Astrodynamics-Toolkit/master/LICENSE).
