!*****************************************************************************************
!> author: Jacob Williams
!
!  Conversion factors.
!
!# See also
!
!  1. A. Thompson and B. N. Taylor, "NIST Special Publication 811:
!     Guide for the use of the International System of Units".
!     http://www.nist.gov/pml/pubs/sp811/

    module conversion_module

    use kind_module,       only: wp
    use numbers_module,    only: one,pi

    implicit none

    public

    !metric/imperial:
    real(wp),parameter :: lbm2kg  = 0.45359237_wp         !! exact
    real(wp),parameter :: lbf2N   = 4.4482216152605_wp    !! exact
    real(wp),parameter :: ft2m    = 0.3048_wp             !! exact
    real(wp),parameter :: mile2km = 1.609344_wp           !! exact
    real(wp),parameter :: nmi2km  = 1.852_wp              !! exact
    real(wp),parameter :: slug2kg = lbf2N/ft2m     !! approximately 14.593902937206362
    real(wp),parameter :: kg2lbm  = one/lbm2kg     !! approximately 2.2046226218487757
    real(wp),parameter :: N2lbf   = one/lbf2N      !! approximately 0.2248089430997105
    real(wp),parameter :: m2ft    = one/ft2m       !! approximately 3.280839895013123
    real(wp),parameter :: km2mile = one/mile2km    !! approximately 0.621371192237334
    real(wp),parameter :: km2nmi  = one/nmi2km     !! approximately 0.5399568034557235
    real(wp),parameter :: kg2slug = ft2m/lbf2N     !! approximately 0.06852176585679176

    !angles:
    real(wp),parameter :: deg2rad = pi/180.0_wp
    real(wp),parameter :: rad2deg = 180.0_wp/pi

    !metric:
    real(wp),parameter :: km2m = 1000.0_wp
    real(wp),parameter :: m2km = one/km2m
    real(wp),parameter :: au2m = 149597870700.0_wp !! IAU 2012 defined value

    !time:
    real(wp),parameter :: min2sec      = 60.0_wp
    real(wp),parameter :: hr2min       = 60.0_wp
    real(wp),parameter :: day2hr       = 24.0_wp
    real(wp),parameter :: year2day     = 365.25_wp         !! julian year
    real(wp),parameter :: century2day  = year2day*100.0_wp !! julian century
    real(wp),parameter :: deg2arcmin   = 60.0_wp
    real(wp),parameter :: deg2arcsec   = 3600.0_wp
    real(wp),parameter :: hr2sec       = hr2min*min2sec
    real(wp),parameter :: day2min      = day2hr*hr2min
    real(wp),parameter :: day2sec      = day2min*min2sec
    real(wp),parameter :: century2sec  = century2day*day2sec
    real(wp),parameter :: day2year     = one/year2day
    real(wp),parameter :: day2century  = one/century2day
    real(wp),parameter :: hr2day       = one/day2hr
    real(wp),parameter :: sec2hr       = one/hr2sec
    real(wp),parameter :: sec2day      = one/day2sec
    real(wp),parameter :: sec2century  = one/century2sec
    real(wp),parameter :: arcmin2deg   = one/deg2arcmin
    real(wp),parameter :: arcsec2deg   = one/deg2arcsec

    end module conversion_module
!*****************************************************************************************
