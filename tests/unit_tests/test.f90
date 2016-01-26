!*****************************************************************************************
!> author: Jacob Williams
!
!  Unit test program for the Fortran Astrodynamics Toolkit

    program test

    use fortran_astrodynamics_toolkit

    implicit none

    call brent_test()
    call ephemeris_test()
    call vector_test()
    call lambert_test()
    call rk_test()
    call complex_step_test()
    call time_module_test()
    call geopotential_module_test()
    call iau_test()
    call relative_motion_test()
    call crtbp_test()
    call bplane_test()
    call direct_inverse_test()
    call modified_equinoctial_test()

    end program test
!*****************************************************************************************
