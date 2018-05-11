!*****************************************************************************************
!> author: Jacob Williams
!
!  Defines some numeric parameters.

    module numbers_module

    use kind_module,  only: wp

    private

    real(wp),parameter,public :: zero       = 0.0_wp
    real(wp),parameter,public :: one        = 1.0_wp
    real(wp),parameter,public :: two        = 2.0_wp
    real(wp),parameter,public :: three      = 3.0_wp
    real(wp),parameter,public :: four       = 4.0_wp
    real(wp),parameter,public :: five       = 5.0_wp
    real(wp),parameter,public :: six        = 6.0_wp
    real(wp),parameter,public :: seven      = 7.0_wp
    real(wp),parameter,public :: eight      = 8.0_wp
    real(wp),parameter,public :: nine       = 9.0_wp
    real(wp),parameter,public :: ten        = 10.0_wp

    real(wp),parameter,public :: pi         = acos(-one)
    real(wp),parameter,public :: twopi      = two*pi
    real(wp),parameter,public :: fourpi     = four*pi

    real(wp),parameter,public :: universal_grav_constant = 6.67408e-20_wp !! CODATA-recommended universal gravitational
                                                                          !! constant \( km^3/kg-s^2  \)

    !> 3x3 identity matrix:
    real(wp),dimension(3,3),parameter,public :: identity_3x3 = reshape(&
                                                    [[one,zero,zero],&
                                                     [zero,one,zero],&
                                                     [zero,zero,one]],[3,3])

    !> 6x6 identity matrix:
    real(wp),dimension(6,6),parameter,public :: identity_6x6= reshape(&
                                                    [[one,zero,zero,zero,zero,zero],&
                                                     [zero,one,zero,zero,zero,zero],&
                                                     [zero,zero,one,zero,zero,zero],&
                                                     [zero,zero,zero,one,zero,zero],&
                                                     [zero,zero,zero,zero,one,zero],&
                                                     [zero,zero,zero,zero,zero,one] ],[6,6])

    end module numbers_module
!*****************************************************************************************
