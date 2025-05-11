!*****************************************************************************************
!> author: Jacob Williams
!
!  For computing the obliquity of the ecliptic.
!
!### Notes
!  * Mean obliquity does not include the short-period terms
!    due to nutation. True obliquity includes these terms.

    module obliquity_module

    use conversion_module
    use kind_module

    implicit none

    private

    abstract interface
        pure function mean_obliquity_func(et) result(e)
        !! a function for computing the mean obliquity of the ecliptic.
        import :: wp
        implicit none
        real(wp),intent(in) :: et  !! ephemeris time (sec)
        real(wp)            :: e   !! obliquity of ecliptic (deg)
        end function mean_obliquity_func
    end interface

    public :: mean_ecliptic_to_equatorial_rotmat
    public :: equatorial_to_mean_ecliptic_rotmat
    public :: mean_obliquity_of_ecliptic_iau1980
    public :: mean_obliquity_of_ecliptic_iau2006

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Rotation matrix from Mean Ecliptic to J2000.
!
!### Reference
! * https://en.wikipedia.org/wiki/Ecliptic_coordinate_system

    pure function mean_ecliptic_to_equatorial_rotmat(obliquity_func) result(rot)

    implicit none

    real(wp),dimension(3,3) :: rot  !! rotation matrix
    procedure(mean_obliquity_func),optional :: obliquity_func !! optional function to compute
                                                              !! the mean obliquity. If not
                                                              !! present, then
                                                              !! [[mean_obliquity_of_ecliptic_iau1980]]
                                                              !! is used.

    real(wp) :: e !! mean obliquity at J2000 (rad)
    real(wp) :: s,c

    if (present(obliquity_func)) then
        e = obliquity_func(0.0_wp)
    else
        e = mean_obliquity_of_ecliptic_iau1980(0.0_wp)
    end if
    e = e * deg2rad
    s = sin(e)
    c = cos(e)

    rot(:,1) = [1.0_wp, 0.0_wp, 0.0_wp]
    rot(:,2) = [0.0_wp, c, s]
    rot(:,3) = [0.0_wp, -s, c]

    end function mean_ecliptic_to_equatorial_rotmat
!*****************************************************************************************

!*****************************************************************************************
!>
!  Rotation matrix from J2000 to Mean Ecliptic.

    pure function equatorial_to_mean_ecliptic_rotmat(obliquity_func) result(rot)

    implicit none

    real(wp),dimension(3,3) :: rot  !! rotation matrix
    procedure(mean_obliquity_func),optional :: obliquity_func !! optional function to compute
                                                              !! the mean obliquity. If not
                                                              !! present, then
                                                              !! [[mean_obliquity_of_ecliptic_iau1980]]
                                                              !! is used.

    rot = transpose(mean_ecliptic_to_equatorial_rotmat(obliquity_func))

    end function equatorial_to_mean_ecliptic_rotmat
!*****************************************************************************************

!*****************************************************************************************
!>
!  Mean obliquity of the ecliptic, IAU 2006 formula.

    pure function mean_obliquity_of_ecliptic_iau2006(et) result(e)

    implicit none

    real(wp),intent(in) :: et  !! ephemeris time (sec)
    real(wp)            :: e   !! obliquity of ecliptic (deg)

    real(wp) :: t  !! time in centuries from the J2000 epoch

    real(wp),parameter,dimension(6) :: c = [84381.406_wp,&
                                            -46.836769_wp,&
                                            -0.0001831_wp,&
                                            0.00200340_wp,&
                                            -0.000000576_wp,&
                                            -0.0000000434_wp] !! coefficients

    ! convert input time to centuries:
    t = et*sec2day*day2century

    ! use horner's rule:
    e = (c(1)+t*(c(2)+t*(c(3)+t*(c(4)+t*(c(5)+t*c(6))))))*arcsec2deg

    end function mean_obliquity_of_ecliptic_iau2006
!*****************************************************************************************

!*****************************************************************************************
!>
!  Mean obliquity of the ecliptic, IAU 1980 formula.
!
!@note This equation is consistent with the one from the SPICE `zzmobliq` routine.

    pure function mean_obliquity_of_ecliptic_iau1980(et) result(e)

    implicit none

    real(wp),intent(in) :: et  !! ephemeris time (sec)
    real(wp)            :: e   !! obliquity of ecliptic (deg)

    real(wp) :: t  !! time in centuries from the J2000 epoch

    real(wp),dimension(0:3),parameter :: c = [84381.448_wp,&
                                              -46.8150_wp,&
                                              -0.00059_wp,&
                                              +0.001813_wp] !! coefficients

    ! convert input time to centuries:
    t = et*sec2day*day2century

    ! use horner's rule:
    e = (c(0)+t*(c(1)+t*(c(2)+t*c(3))))*arcsec2deg

    end function mean_obliquity_of_ecliptic_iau1980
!*****************************************************************************************

!*****************************************************************************************
    end module obliquity_module
!*****************************************************************************************
