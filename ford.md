project: fortran-astrodynamics-toolkit
src_dir: ./src
output_dir: ./doc
media_dir: ./media
project_github: https://github.com/jacobwilliams/Fortran-Astrodynamics-Toolkit
github: https://github.com/jacobwilliams
website: https://degenerateconic.com
twitter: https://twitter.com/degenerateconic
summary: ![Fortran Astrodynamics Toolkit](|media|/logo.png)<br>
         A modern Fortran library for Astrodynamics
author: Jacob Williams
github: https://github.com/jacobwilliams
predocmark_alt: >
predocmark: <
docmark_alt:
docmark: !
display: public
         protected
         private
source: true
graph: true
exclude: pyplot_module.f90
exclude_dir: ./tests
extra_mods: pyplot_module:https://github.com/jacobwilliams/pyplot-fortran
            iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html

{!README.md!}