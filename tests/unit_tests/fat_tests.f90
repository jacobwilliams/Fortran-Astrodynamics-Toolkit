!*****************************************************************************************
!> author: Jacob Williams
!
!  Unit test program for the Fortran Astrodynamics Toolkit

    program fat_tests

    use fortran_astrodynamics_toolkit

    implicit none

    call brent_test()
    call ephemeris_test()
    call vector_test()
    call lambert_test()
    call complex_step_test()
    call time_module_test()
    call geopotential_module_test()
    call iau_test()
    call relative_motion_test()
    call crtbp_test()
    call bplane_test()
    call direct_inverse_test()
    call modified_equinoctial_test()
    call transformation_module_test()
    call step_size_test()
    call rk_test()
    call rk_test_variable_step()
    call standish_module_test()
    call halo_orbit_test()
    call eispack_test()

    end program fat_tests
!*****************************************************************************************
