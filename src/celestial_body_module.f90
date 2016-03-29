!*****************************************************************************************
!> author: Jacob Williams
!
!  Celestial body definitions

    module celestial_body_module

    use kind_module
    use numbers_module
    use base_class_module

    implicit none

    type,extends(base_class),public :: celestial_body
        !! A celestial body (Planet, moon, etc.)
        !! The `ID` from the [[base_class]] is the NAIF SPICE ID code for the body
        real(wp) :: mu = zero   !! gravitational parameter \( \mu \) [\(km^3/s^2\)]
    end type celestial_body

    !define some bodies:
    type(celestial_body),parameter,public :: body_earth = celestial_body(3,'Earth',3.9860043543609593E+05_wp) !! the earth
    type(celestial_body),parameter,public :: body_moon  = celestial_body(10,'Moon',4.9028000661637961E+03_wp) !! the moon

!*****************************************************************************************
    end module celestial_body_module
!*****************************************************************************************
