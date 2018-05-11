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
    ! MU values from: https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de431.tpc

    type(celestial_body),parameter,public :: body_sun = &
                    celestial_body(10, 'Sun',1.3271244004193938E+11_wp )
    type(celestial_body),parameter,public :: body_mercury = &
                    celestial_body(199,'Mercury',2.2031780000000021E+04_wp )
    type(celestial_body),parameter,public :: body_venus = &
                    celestial_body(299,'Venus',3.2485859200000006E+05_wp )
    type(celestial_body),parameter,public :: body_earth = &
                    celestial_body(399,'Earth',3.9860043543609598E+05_wp )
    type(celestial_body),parameter,public :: body_earth_moon_barycenter = &
                    celestial_body(3,'Earth-Moon Barycenter',4.0350323550225981E+05_wp )
    type(celestial_body),parameter,public :: body_moon = &
                    celestial_body(301,'Moon',4.9028000661637961E+03_wp )
    type(celestial_body),parameter,public :: body_mars = &
                    celestial_body(499,'Mars',4.282837362069909E+04_wp  )
    type(celestial_body),parameter,public :: body_jupiter = &
                    celestial_body(599,'Jupiter',1.266865349218008E+08_wp  )
    type(celestial_body),parameter,public :: body_saturn = &
                    celestial_body(699,'Saturn',3.793120749865224E+07_wp  )
    type(celestial_body),parameter,public :: body_uranus = &
                    celestial_body(799,'Uranus',5.793951322279009E+06_wp  )
    type(celestial_body),parameter,public :: body_neptune = &
                    celestial_body(899,'Neptune',6.835099502439672E+06_wp  )
    type(celestial_body),parameter,public :: body_pluto = &
                    celestial_body(999,'Pluto',8.696138177608748E+02_wp  )

!*****************************************************************************************
    end module celestial_body_module
!*****************************************************************************************
