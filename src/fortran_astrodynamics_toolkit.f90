!*****************************************************************************************
!> author: Jacob Williams
!  license: BSD
!
!  The main module that uses all the other modules.
!  Allows for a single `use fortran_astrodynamics_toolkit`
!  to access the entire library.

    module fortran_astrodynamics_toolkit

    use analytical_ephemeris_module
    use base_class_module
    use bplane_module
    use brent_module
    use celestial_body_module
    use complex_step_module
    use conversion_module
    use crtbp_module
    use drag_module
    use eispack_module
    use ephemeris_module
    use geodesy_module
    use geometry_module
    use geopotential_module
    use gooding_module
    use gravity_module
    use halo_orbit_module
    use iau_orientation_module
    use jpl_ephemeris_module
    use kepler_module
    use kind_module
    use lambert_module
    use math_module
    use matrix_module
    use minpack_module
    use modified_equinoctial_module
    use numbers_module
    use obliquity_module
    use orbital_mechanics_module
    use random_module
    use relative_motion_module
    use rk_module
    use rk_module_variable_step
    use standish_module
    use string_module
    use time_module
    use transformation_module
    use vector_module

    use c_interface_module

    implicit none

    public

    private :: wp

    integer,parameter,public :: fat_wp = wp  !! default real kind

    end module fortran_astrodynamics_toolkit
!*****************************************************************************************
