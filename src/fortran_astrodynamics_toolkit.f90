!*****************************************************************************************
!> author: Jacob Williams
!  license: BSD
!
!  The main module that uses all the other modules.

    module fortran_astrodynamics_toolkit

    use brent_module
    use complex_step_module
    use conversion_module
    use ephemeris_module
    use geodesy_module
    use geometry_module
    use gooding_module
    use geopotential_module
    use iau_orientation_module
    use kind_module
    use lambert_module
    use math_module
    use minpack_module
    use numbers_module
    use relative_motion_module
    use rk_module
    use time_module
    use vector_module
    
    implicit none
    
    public
               
    end module fortran_astrodynamics_toolkit
!*****************************************************************************************