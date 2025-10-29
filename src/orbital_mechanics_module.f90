!*****************************************************************************************
!> author: Jacob Williams
!
!  Basic orbital mechanics routines.

    module orbital_mechanics_module

    use kind_module
    use numbers_module
    use vector_module
    use math_module, only: wrap_angle

    implicit none

    private

    public :: rv_to_orbital_elements
    public :: orbital_elements_to_rv
    public :: orbit_period
    public :: orbit_energy
    public :: periapsis_apoapsis
    public :: sphere_of_influence
    public :: sphere_of_influence_earth_moon

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Convert position and velocity vectors to orbital elements.
!
!### See also
!  * The poliastro routine `rv2coe`.

    pure subroutine rv_to_orbital_elements(mu, r, v, p, ecc, inc, raan, aop, tru)

    implicit none

    real(wp),intent(in)              :: mu   !! gravitational parameter [\(km^{3}/s^{2}\)]
    real(wp),dimension(3),intent(in) :: r    !! position vector [km]
    real(wp),dimension(3),intent(in) :: v    !! velocity vector [km/s]
    real(wp),intent(out)             :: p    !! semiparameter \(a(1-e^{2})\) [km]
    real(wp),intent(out)             :: ecc  !! eccentricity [--]
    real(wp),intent(out)             :: inc  !! inclination [rad]
    real(wp),intent(out)             :: raan !! raan [rad]
    real(wp),intent(out)             :: aop  !! argument of peripsis [rad]
    real(wp),intent(out)             :: tru  !! true anomaly [rad]

    real(wp),dimension(3) :: h,n,e
    logical :: circular,equatorial
    reaL(wp) :: hmag,rmag,vmag

    rmag = norm2(r)
    vmag = norm2(v)
    h    = cross(r,v)
    hmag = norm2(h)
    n    = cross([zero,zero,one],h)/hmag
    e    = ((vmag**2-mu/rmag)*r-dot_product(r,v)*v)/mu
    ecc  = norm2(e)
    p    = hmag**2/mu
    inc  = atan2(norm2(h(1:2)),h(3))

    call orbit_check(ecc,inc,circular,equatorial)

    if (equatorial .and. .not. circular) then
        raan = zero
        aop  = wrap_angle(atan2(e(2),e(1))) ! Longitude of periapsis
        tru  = wrap_angle(atan2(dot_product(h,cross(e,r))/hmag,dot_product(r,e)))
    elseif ( .not. equatorial .and. circular) then
        raan = wrap_angle(atan2(n(2),n(1)))
        aop  = zero
        tru  = wrap_angle(atan2(dot_product(r,cross(h,n))/hmag,dot_product(r,n))) ! Argument of latitude
    elseif (equatorial .and. circular) then
        raan = zero
        aop  = zero
        tru  = wrap_angle(atan2(r(2),r(1))) ! True longitude
    else
        raan = wrap_angle(atan2(n(2),n(1)))
        aop  = wrap_angle(atan2(dot_product(e,cross(h,n))/hmag,dot_product(e,n)))
        tru  = wrap_angle(atan2(dot_product(r,cross(h,e))/hmag,dot_product(r,e)))
    endif

    end subroutine rv_to_orbital_elements
!*****************************************************************************************

!*****************************************************************************************
!>
!  Convert orbital elements to position and velocity vectors.

    pure subroutine orbital_elements_to_rv(mu, p, ecc, inc, raan, aop, tru, r, v)

    implicit none

    real(wp),intent(in)               :: mu   !! gravitational parameter [\(km^{3}/s^{2}\)]
    real(wp),intent(in)               :: p    !! semiparameter \(a(1-e^{2})\) [km]
    real(wp),intent(in)               :: ecc  !! eccentricity [--]
    real(wp),intent(in)               :: inc  !! inclination [rad]
    real(wp),intent(in)               :: raan !! raan [rad]
    real(wp),intent(in)               :: aop  !! argument of peripsis [rad]
    real(wp),intent(in)               :: tru  !! true anomaly [rad]
    real(wp),dimension(3),intent(out) :: r    !! position vector [km]
    real(wp),dimension(3),intent(out) :: v    !! velocity vector [km/s]

    real(wp),dimension(3,2) :: rotmat
    real(wp),dimension(2)   :: r_pqw,v_pqw
    logical                 :: circular,equatorial
    real(wp)                :: ctru,stru,sr,cr,si,ci,sa,ca,raan_tmp,aop_tmp

    call orbit_check(ecc,inc,circular,equatorial)

    if (circular) then   ! periapsis undefined
        aop_tmp = zero
    else
        aop_tmp = aop
    end if

    if (equatorial) then   ! node undefined
        raan_tmp = zero
    else
        raan_tmp = raan
    end if

    ! perifocal position and velocity:
    ctru   = cos(tru)
    stru   = sin(tru)
    r_pqw  = [ctru, stru] * p/(one+ecc*ctru)
    v_pqw  = [-stru, (ecc+ctru)] * sqrt(mu/p)

    ! perifocal to cartesian:
    sr          = sin(raan_tmp)
    cr          = cos(raan_tmp)
    si          = sin(inc)
    ci          = cos(inc)
    sa          = sin(aop_tmp)
    ca          = cos(aop_tmp)
    rotmat(1,:) = [cr*ca-sr*sa*ci, -cr*sa-sr*ca*ci]
    rotmat(2,:) = [sr*ca+cr*sa*ci, -sr*sa+cr*ca*ci]
    rotmat(3,:) = [sa*si, ca*si]

    ! transform:
    r = matmul(rotmat,r_pqw)
    v = matmul(rotmat,v_pqw)

    end subroutine orbital_elements_to_rv
!*****************************************************************************************

!*****************************************************************************************
!>
!  Check the orbit for singularities.

    pure subroutine orbit_check(ecc,inc,circular,equatorial)

    implicit none

    real(wp),intent(in) :: ecc        !! eccentricity
    real(wp),intent(in) :: inc        !! inclination [rad]
    logical,intent(out) :: circular   !! is the orbit circular?
    logical,intent(out) :: equatorial !! is the orbit equatorial?

    real(wp),parameter :: tol = 1.0e-10_wp !! tolerance for circular & equatorial checks

    circular   = ecc < tol
    equatorial = (one - abs(cos(inc))) < tol  ! 0 or 180 deg

    end subroutine orbit_check
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute the two-body orbital period.

    pure function orbit_period(mu,a) result(period)

    implicit none

    real(wp),intent(in) :: mu !! gravitational parameter [\(km^{3}/s^{2}\)]
    real(wp),intent(in) :: a  !! semimajor axis [km]
    real(wp) :: period        !! two-body orbital period [s]

    period = twopi/sqrt(mu/a**3)

    end function orbit_period
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute the two-body orbital energy.

    pure function orbit_energy(mu,rv) result(energy)

    implicit none

    real(wp),intent(in) :: mu               !! gravitational parameter [\(km^{3}/s^{2}\)]
    real(wp),dimension(6),intent(in) :: rv  !! position and velocity vector [km,km/s]
    real(wp) :: energy                      !! two-body orbital energy [\(km^{2}/s^{2}\)]

    real(wp) :: rmag  !! position vector magnitude [km]
    real(wp) :: vmag  !! velocity vector magnitude [km]

    rmag = norm2(rv(1:3))
    vmag = norm2(rv(4:6))

    energy = (vmag**2 / two) - (mu / rmag)

    end function orbit_energy
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute the periapsis and apoapsis position and velocity magnitudes.

    pure subroutine periapsis_apoapsis(mu,a,e,rp,ra,vp,va)

    implicit none

    real(wp),intent(in)  :: mu  !! gravitational parameter [\(km^{3}/s^{2}\)]
    real(wp),intent(in)  :: a   !! semimajor axis [km]
    real(wp),intent(in)  :: e   !! eccentricity [--]
    real(wp),intent(out) :: rp  !! periapsis position magnitude [km]
    real(wp),intent(out) :: ra  !! apoapsis position magnitude [km]
    real(wp),intent(out) :: vp  !! periapsis velocity magnitude [km/s]
    real(wp),intent(out) :: va  !! apoapsis velocity magnitude [km/s]

    real(wp) :: rarp   !! \(r_a + r_p \)
    real(wp) :: twomu  !! \( 2 \mu \)

    twomu = two * mu

    rp   = a*(one-e)
    ra   = a*(one+e)
    rarp = ra + rp
    vp   = sqrt(twomu*ra/(rp*rarp))
    va   = sqrt(twomu*rp/(ra*rarp))

    end subroutine periapsis_apoapsis
!*****************************************************************************************

!*****************************************************************************************
!>
!  Computes the sphere-of-influence radius of the secondary body.
!
!### See also
!  * R.H. Battin, "An Introduction to the Mathematics and
!    Methods of Astrodynamics, Revised Edition", AIAA, 1999.
!  * This is the approximate formula (8.74 from Battin).

    pure function sphere_of_influence(mu_primary, mu_secondary, r_ps) result(r_soi)

    implicit none

    real(wp),intent(in)   :: mu_primary   !! primary body gravitational parameter [km^3/s^2]
    real(wp),intent(in)   :: mu_secondary !! secondary body gravitational parameter [km^3/s^2]
    real(wp),intent(in)   :: r_ps         !! distance between the two bodies [km]
    real(wp)              :: r_soi        !! sphere of influence radius [km]

    real(wp),parameter :: two_fifths = two/five

    if (mu_primary>zero .and. mu_secondary>zero .and. r_ps>zero) then
        r_soi = r_ps * (mu_secondary/mu_primary)**two_fifths
    else
        r_soi = zero
    end if

    end function sphere_of_influence
!*****************************************************************************************

!*****************************************************************************************
!>
!  Computes the sphere-of-influence radius of the secondary body.
!
!### Notes
!  * `r` and `r_sp` should be in the same inertial frame.
!  * The mass of the spacecraft is neglected.
!
!### See also
!  * R.H. Battin, "An Introduction to the Mathematics and
!    Methods of Astrodynamics, Revised Edition", AIAA, 1999.
!  * This is the more complex formula, on p. 397 of Battin,
!    which is better for the Earth/Moon system.

    pure function sphere_of_influence_earth_moon(mu_primary, mu_secondary, r, r_sp) result(r_soi)

    implicit none

    real(wp),intent(in)              :: mu_primary   !! primary body gravitational parameter [km^3/s^2]
    real(wp),intent(in)              :: mu_secondary !! secondary body gravitational parameter [km^3/s^2]
    real(wp),dimension(3),intent(in) :: r            !! vector from the secondary body to the spacecraft [km]
    real(wp),dimension(3),intent(in) :: r_sp         !! vector from the secondary to the primary body [km]
    real(wp)                         :: r_soi        !! sphere of influence radius of the secondary body [km]

    real(wp) :: r_mag, r_sp_mag, alpha, ca, ca2, denom

    r_mag    = norm2(r)
    r_sp_mag = norm2(r_sp)

    if (mu_primary>zero .and. mu_secondary>zero .and. r_mag>zero .and. r_sp_mag>zero) then

        alpha = angle_between_vectors(r,r_sp)
        ca    = cos(alpha)
        ca2   = ca * ca

        denom = (mu_secondary**2/mu_primary**2)**(one/five)*(one+three*ca2)**(one/ten) + &
                (two/five)*ca*((one+six*ca2)/(one+three*ca2))

        r_soi = r_sp_mag / denom

    else
        r_soi = zero
    end if

    end function sphere_of_influence_earth_moon
!*****************************************************************************************

!*****************************************************************************************
    end module orbital_mechanics_module
!*****************************************************************************************
