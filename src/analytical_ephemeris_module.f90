!*****************************************************************************************
!>
!  Analytical ephemeris routines for solar system bodies.

    module analytical_ephemeris_module

    use numbers_module
    use kind_module
    use conversion_module

    implicit none

    private

    public :: simpson_lunar_ephemeris

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  A simple analytical lunar ephemeris model.
!  Returns Lunar cartesian coordinates
!  (mean equator and equinox of epoch J2000).
!
!### Reference
!  * "An alternative lunar ephemeris model for on-board flight software use",
!    D. G. Simpson, Proceedings of the 1999 NASA/GSFC Flight Mechanics Symposium,
!    p. 175-184).
!  * [David G. Simpson Personal Web Site](http://www.davidgsimpson.com/software.html)
!
!### Reference code
!
!```fortran
!  T = (JD - 2451545.0D0)/36525.0D0
!  X =  383.0D3 * SIN( 8399.685D0 * T + 5.381D0)   &
!      + 31.5D3 * SIN(   70.990D0 * T + 6.169D0)   &
!      + 10.6D3 * SIN(16728.377D0 * T + 1.453D0)   &
!      +  6.2D3 * SIN( 1185.622D0 * T + 0.481D0)   &
!      +  3.2D3 * SIN( 7143.070D0 * T + 5.017D0)   &
!      +  2.3D3 * SIN(15613.745D0 * T + 0.857D0)   &
!      +  0.8D3 * SIN( 8467.263D0 * T + 1.010D0)
!  Y =  351.0D3 * SIN( 8399.687D0 * T + 3.811D0)   &
!      + 28.9D3 * SIN(   70.997D0 * T + 4.596D0)   &
!      + 13.7D3 * SIN( 8433.466D0 * T + 4.766D0)   &
!      +  9.7D3 * SIN(16728.380D0 * T + 6.165D0)   &
!      +  5.7D3 * SIN( 1185.667D0 * T + 5.164D0)   &
!      +  2.9D3 * SIN( 7143.058D0 * T + 0.300D0)   &
!      +  2.1D3 * SIN(15613.755D0 * T + 5.565D0)
!  Z =  153.2D3 * SIN( 8399.672D0 * T + 3.807D0)   &
!      + 31.5D3 * SIN( 8433.464D0 * T + 1.629D0)   &
!      + 12.5D3 * SIN(   70.996D0 * T + 4.595D0)   &
!      +  4.2D3 * SIN(16728.364D0 * T + 6.162D0)   &
!      +  2.5D3 * SIN( 1185.645D0 * T + 5.167D0)   &
!      +  3.0D3 * SIN(  104.881D0 * T + 2.555D0)   &
!      +  1.8D3 * SIN( 8399.116D0 * T + 6.248D0)
!```
!
!@note Also added velocity output, which is not present in reference code.

    subroutine simpson_lunar_ephemeris(jd,r_moon,v_moon)

    implicit none

    real(wp),intent(in)                        :: jd     !! Julian date
    real(wp),dimension(3),intent(out)          :: r_moon !! Moon position (km)
    real(wp),dimension(3),intent(out),optional :: v_moon !! Moon velocity (km/s)

    real(wp) :: t  !! time in Julian centuries from J2000
    real(wp),dimension(7) :: xterms,yterms,zterms

    !coefficients:
    real(wp),dimension(7),parameter :: xcoeffs = [  383.0e3_wp, 31.5e3_wp, &
                                                     10.6e3_wp, 6.2e3_wp, &
                                                      3.2e3_wp, 2.3e3_wp, &
                                                      0.8e3_wp]
    real(wp),dimension(7),parameter :: ycoeffs = [  351.0e3_wp, 28.9e3_wp, &
                                                     13.7e3_wp, 9.7e3_wp, &
                                                      5.7e3_wp, 2.9e3_wp, &
                                                      2.1e3_wp]
    real(wp),dimension(7),parameter :: zcoeffs = [  153.2e3_wp, 31.5e3_wp, &
                                                     12.5e3_wp, 4.2e3_wp, &
                                                      2.5e3_wp, 3.0e3_wp, &
                                                      1.8e3_wp]
    real(wp),dimension(7),parameter :: xa = [ 8399.685_wp, 70.990_wp, &
                                             16728.377_wp, 1185.622_wp, &
                                              7143.070_wp, 15613.745_wp, &
                                              8467.263_wp]
    real(wp),dimension(7),parameter :: xp = [ 5.381_wp, 6.169_wp, &
                                              1.453_wp, 0.481_wp, &
                                              5.017_wp, 0.857_wp, &
                                              1.010_wp]
    real(wp),dimension(7),parameter :: ya = [ 8399.687_wp, 70.997_wp, &
                                              8433.466_wp, 16728.380_wp, &
                                              1185.667_wp, 7143.058_wp, &
                                             15613.755_wp]
    real(wp),dimension(7),parameter :: yp = [3.811_wp, 4.596_wp, &
                                             4.766_wp, 6.165_wp, &
                                             5.164_wp, 0.300_wp, &
                                             5.565_wp]
    real(wp),dimension(7),parameter :: za = [ 8399.672_wp, 8433.464_wp, &
                                                70.996_wp, 16728.364_wp, &
                                              1185.645_wp, 104.881_wp, &
                                              8399.116_wp]
    real(wp),dimension(7),parameter :: zp = [3.807_wp, 1.629_wp, &
                                             4.595_wp, 6.162_wp, &
                                             5.167_wp, 2.555_wp, &
                                             6.248_wp]
    real(wp),dimension(7),parameter :: vxcoeffs = xcoeffs*xa
    real(wp),dimension(7),parameter :: vycoeffs = ycoeffs*ya
    real(wp),dimension(7),parameter :: vzcoeffs = zcoeffs*za

    t = (jd - 2451545.0_wp)*day2century

    xterms = xa * t + xp
    yterms = ya * t + yp
    zterms = za * t + zp

    r_moon(1) = dot_product(xcoeffs, sin(xterms))
    r_moon(2) = dot_product(ycoeffs, sin(yterms))
    r_moon(3) = dot_product(zcoeffs, sin(zterms))

    !v_moon is just d(r_moon)/dt:
    ! [convert units to km/s]
    if (present(v_moon)) then
        v_moon(1) =  dot_product(vxcoeffs, cos(xterms)) / (century2day * day2sec)
        v_moon(2) =  dot_product(vycoeffs, cos(yterms)) / (century2day * day2sec)
        v_moon(3) =  dot_product(vzcoeffs, cos(zterms)) / (century2day * day2sec)
    end if

    end subroutine simpson_lunar_ephemeris
!*****************************************************************************************

!*****************************************************************************************
    end module analytical_ephemeris_module
!*****************************************************************************************
