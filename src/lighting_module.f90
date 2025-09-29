!*****************************************************************************************
!> author: Jacob Williams
!
!  Routines for computing solar fraction, lighting, eclipses, etc.

    module lighting_module

    use kind_module,      only: wp
    use numbers_module,   only: pi, zero, one
    use vector_module,    only: unit

    implicit none

    private

    public :: solar_fraction

    contains 
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute the solar fraction visible due to an eclipse by another body.
!
!### Reference
!  * J. Wertz, "Spacecraft Attitude Determination and Control", 1978. 
!    See Chapter 3 and Appendix A.

    function solar_fraction(rp, rs, d_s, d_p) result(fraction)

    real(wp),intent(in) :: rp !! radius of the planet
    real(wp),intent(in) :: rs !! radius of the Sun
    real(wp),dimension(3),intent(in)  :: d_s !! vector from the spacecraft to the Sun
    real(wp),dimension(3),intent(in)  :: d_p !! vector from the spacecraft to the planet
    real(wp) :: fraction !! fraction of the Sun visible [0=total eclipse, 1=no eclipse]

    real(wp) :: s !! distance from the planet to the Sun
    real(wp) :: c !! distance from the center of the planet to the apex of the shadow cone
    real(wp) :: rho_c !! angular radius of the shadow cone
    real(wp) :: rho_s !! angular radius of the Sun
    real(wp) :: rho_p !! angular radius of the planet
    real(wp) :: theta !! angular separation of the sun and planet as viewed by the spacecraft
    real(wp) :: ds !! distance from the spacecraft to the Sun
    real(wp) :: dp !! distance from the spacecraft to the planet
    real(wp) :: drho !! difference in angular radii of the planet and Sun
    real(wp) :: crp, crs, srp, srs, cth, sth, t1, t2, t3 !! temp variables

    if (rp<=zero) then ! no eclipse possible if the planet has no radius
        fraction = one
        return
    end if

    ds = norm2(d_s) 
    dp = norm2(d_p)

    if (ds<=rs) then ! inside the Sun
        fraction = one
        return
    else if (dp<=rp) then ! inside the planet
        fraction = zero
        return
    end if

    s     = norm2(d_s - d_p)
    c     = (rp*s) / (rs - rp)
    rho_c = asin((rs - rp) / s)  ! appx = asin(rs/s)
    rho_s = asin(rs/ds)
    rho_p = asin(rp/dp)
    theta = acos(dot_product(unit(d_s), unit(d_p)))
    drho  = rho_p - rho_s
    crp   = cos(rho_p)
    crs   = cos(rho_s)
    srp   = sin(rho_p)
    srs   = sin(rho_s)
    cth   = cos(theta)
    sth   = sin(theta)

    if ( (ds>s) .and. (rho_p+rho_s>theta) .and. (theta>abs(drho)) ) then
        ! partial eclipse
        t1 = pi - crs * acos( (crp-crs*cth)/(srs*sth) )
        t2 =     -crp * acos( (crs-crp*cth)/(srp*sth) )
        t3 =           -acos( (cth-crs*crp)/(srs*srp) )
        fraction = (t1 + t2 + t3) / (pi*(one-crs)) 
    else if ( (s<ds) .and. (ds<s+c) .and. (drho>theta) ) then 
        ! total eclipse 
        fraction = 0.0_wp
    else if ( (s+c<ds) .and. (drho>theta) ) then 
        ! annular eclipse
        fraction = (one-crp) / (one-crs)
    else 
        ! no eclipse 
        fraction = 1.0_wp
    end if

    end function solar_fraction
!*****************************************************************************************

    end module lighting_module
!*****************************************************************************************