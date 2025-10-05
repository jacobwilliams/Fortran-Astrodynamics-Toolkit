!*****************************************************************************************
!> author: Jacob Williams
!
!  Routines for computing solar fraction, lighting, eclipses, etc.

    module lighting_module

    use kind_module,           only: wp
    use numbers_module,        only: pi, zero, one, two
    use vector_module,         only: unit, cross, axis_angle_rotation
    use ephemeris_module,      only: ephemeris_class
    use transformation_module, only: icrf_frame
    use celestial_body_module, only: celestial_body, body_sun, body_ssb
    use conversion_module,     only: deg2rad, rad2deg
    use math_module,           only: wrap_angle

    implicit none

    private

    real(wp),parameter :: c_light = 299792.458_wp !! speed of light in km/s

    public :: from_j2000body_to_j2000ssb
    public :: apparent_position
    public :: get_sun_fraction   ! high-level routine
    public :: solar_fraction      ! low-level routine
    public :: solar_fraction_alt  ! low-level routine
    public :: solar_fraction_alt2 ! low-level routine
    public :: cubic_shadow_model  ! low-level routine

    public :: lighting_module_test

    contains
!*****************************************************************************************

!********************************************************************************
!>
!  Compute the "sun fraction" using the selected shadow model.

    function get_sun_fraction(b, rad_body, rad_sun, eph, et, rv, model, rbubble, use_geometric, info) result (phi)

    type(celestial_body),intent(in)      :: b           !! eclipsing body
    real(wp),intent(in)                  :: rad_body    !! radius of the eclipsing body [km]
    real(wp),intent(in)                  :: rad_sun     !! radius of the Sun [km]
    class(ephemeris_class),intent(inout) :: eph         !! the ephemeris to use for sun and ssb (if necessary)
    real(wp),intent(in)                  :: et          !! observer ephemeris time (sec)
    real(wp),dimension(6),intent(in)     :: rv          !! state of the spacecraft (j2000-body frame)
    integer,intent(in)                   :: model       !! 1=circular cubic shadow model
                                                        !! 2=solar fraction model
    real(wp),intent(in)                  :: rbubble     !! eclipse bubble [km]. see the reference.
                                                        !! if rbubble=0, then no bubble is used.
                                                        !! only used if model=1
    logical,intent(in),optional          :: use_geometric  !! if true, use geometric positions
                                                           !! (no light time or stellar aberration correction)
                                                           !! default = false
    character(len=:),allocatable,intent(out),optional :: info !! info string
    real(wp) :: phi !! if `model=1`, circular cubic sun frac value:
                    !!
                    !!  * >0 no eclipse,
                    !!  * <0 eclipse,
                    !!  * =0 on the eclipse line
                    !!
                    !! if `model=2`, true solar fraction value [0=total eclipse, 1=no eclipse],
                    !! with model of umbra/penumbra/antumbra (Wertz, 1978)
                    !! if `model=3`, alternate version of solar fraction (Montenbruck and Gill)
                    !! if `model=4`, alternate version of solar fraction (nyxspace)

    logical :: status_ok !! true if no problems
    real(wp),dimension(3) :: r_sun !! apparent state of the sun (j2000-ssb frame)
    real(wp),dimension(3) :: r_body !! apparent state of the eclipsing body (j2000-ssb frame)
    real(wp),dimension(6) :: rv_ssb !! state of the spacecraft !! (j2000-ssb frame)
    logical :: use_apparent

    if (present(use_geometric)) then
        use_apparent = .not. use_geometric
    else
        use_apparent = .true.
    end if

    if (use_apparent) then
        ! apparent position of sun and body wrt to the spacecraft
        call from_j2000body_to_j2000ssb(b, eph, et, rv, rv_ssb) ! state of spacecraft in j2000-ssb
        call apparent_position(eph, body_sun, et, rv_ssb, r_sun,  status_ok) ! apparent position of sun in j2000
        if (.not. status_ok) error stop 'error getting apparent sun position'
        call apparent_position(eph, b,        et, rv_ssb, r_body, status_ok) ! apparent position of body in j2000
        if (.not. status_ok) error stop 'error getting apparent body position'
    else
        ! use geometric positions
        r_body = -rv(1:3) ! geometric position of body wrt spacecraft in j2000
        call eph%get_r(et, body_sun, b, r_sun, status_ok) ! geometric position of sun wrt body in j2000
        if (.not. status_ok) error stop 'error getting geometric sun position'
        r_sun = r_body + r_sun ! geometric position of sun wrt spacecraft in j2000
    end if

    ! compute sun fraction value
    select case(model)
    case(1); call cubic_shadow_model( r_sun, rad_sun, r_body, rad_body, phi, rbubble)
    case(2); call solar_fraction(     r_sun, rad_sun, r_body, rad_body, phi, info)
    case(3); call solar_fraction_alt( r_sun, rad_sun, r_body, rad_body, phi, info)
    case(4); call solar_fraction_alt2(r_sun, rad_sun, r_body, rad_body, phi)
    case default
        error stop 'invalid sun fraction model'
    end select

    end function get_sun_fraction
!********************************************************************************

!*****************************************************************************************
!>
!  Compute the solar fraction visible due to an eclipse by another body.
!
!### Reference
!  * J. Wertz, "Spacecraft Attitude Determination and Control", 1978.
!    See Chapter 3 and Appendix A. Note that the implementation here corrects
!    a typo in this reference, and also protects for a division by zero.

    subroutine solar_fraction(d_s, rs, d_p, rp, fraction, info)

    real(wp),dimension(3),intent(in) :: d_s !! vector from the spacecraft to the Sun
    real(wp),intent(in)              :: rs  !! radius of the Sun
    real(wp),dimension(3),intent(in) :: d_p !! vector from the spacecraft to the planet
    real(wp),intent(in)              :: rp  !! radius of the planet
    real(wp),intent(out)             :: fraction !! fraction of the Sun visible [0=total eclipse, 1=no eclipse]
    character(len=:),allocatable,intent(out),optional :: info !! info string

    real(wp) :: s !! distance from the planet to the Sun
    real(wp) :: c !! distance from the center of the planet to the apex of the shadow cone
    real(wp) :: rho_c !! angular radius of the shadow cone
    real(wp) :: rho_s !! angular radius of the Sun
    real(wp) :: rho_p !! angular radius of the planet
    real(wp) :: theta !! angular separation of the sun and planet as viewed by the spacecraft
    real(wp) :: ds !! distance from the spacecraft to the Sun
    real(wp) :: dp !! distance from the spacecraft to the planet
    real(wp) :: drho !! difference in angular radii of the planet and Sun
    real(wp) :: crp, crs, srp, srs, cth, sth, t1, t2, t3, delr !! temp variables

    if (rp<=zero) then ! no eclipse possible if the planet has no radius
        if (present(info)) info = 'no eclipse: planet radius <= 0'
        fraction = one
        return
    end if

    ds = norm2(d_s)
    dp = norm2(d_p)

    if (ds<=rs) then ! inside the Sun
        if (present(info)) info = 'inside Sun'
        fraction = one
        return
    else if (dp<=rp) then ! inside the planet
        if (present(info)) info = 'inside Planet'
        fraction = zero
        return
    end if

    s    = norm2(d_s - d_p)
    delr = rs - rp
    if (delr==zero) then
        ! special case when the bodies are the same size,
        ! to avoid division by zero
        c = huge(1.0_wp)
    else
        c = (rp*s) / delr
    end if
    rho_c = asin(delr / s)  ! appx = asin(rs/s)
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
        if (present(info)) info = 'penumbra'
        t1 = pi - crs * acos( (crp-crs*cth)/(srs*sth) )
        t2 =     -crp * acos( (crs-crp*cth)/(srp*sth) )
        t3 =           -acos( (cth-crs*crp)/(srs*srp) )
        fraction = one - (t1 + t2 + t3) / (pi*(one-crs))
    else if ( (s<ds) .and. (ds<s+c) .and. (drho>theta) ) then
        ! total eclipse
        if (present(info)) info = 'umbra'
        fraction = zero
    else if ( (s+c<ds) .and. (drho<theta) ) then   ! JW : typo in original reference
        ! annular eclipse
        if (present(info)) info = 'antumbra'
        fraction = one - (one-crp) / (one-crs)
    else
        ! no eclipse
        if (present(info)) info = 'full sun'
        fraction = one
    end if

    end subroutine solar_fraction
!*****************************************************************************************

!********************************************************************************
!>
!  convert from a j2000-body frame to a j2000-ssb frame.

    subroutine from_j2000body_to_j2000ssb(b, eph, et, rv, rv_ssb)

    type(celestial_body),intent(in) :: b !! eclipsing body
    class(ephemeris_class),intent(inout) :: eph !! the ephemeris to use for body and ssb
    real(wp),intent(in) :: et !! ephemeris time (sec)
    real(wp),dimension(6),intent(in) :: rv !! j2000-body state (km, km/s)
    real(wp),dimension(6),intent(out) :: rv_ssb !! j2000-ssb state (km, km/s)

    type(icrf_frame) :: f1, f2
    logical :: status_ok

    f1 = icrf_frame(b=b)
    f2 = icrf_frame(b=body_ssb)
    call f1%transform(rv,f2,et,eph,rv_ssb,status_ok) ! from f1 to f2
    if (.not. status_ok) error stop 'transformation error in from_j2000body_to_j2000ssb'

    end subroutine from_j2000body_to_j2000ssb
!********************************************************************************

!********************************************************************************
!>
!  Return the position of a target body relative to an observer,
!  corrected for light time and stellar aberration.
!
!  see the SPICELIB routine `spkapo` (with 'lt+s')

    subroutine apparent_position(eph, b_target, et, rv_obs_ssb, r_target, status_ok)

    class(ephemeris_class),intent(inout) :: eph !! the ephemeris
    type(celestial_body),intent(in) :: b_target !! target body
    real(wp),dimension(6),intent(in) :: rv_obs_ssb !! state of the observer
                                                   !! (j2000 frame w.r.t. solar system barycenter)
    real(wp),intent(in) :: et !! observer ephemeris time (sec)
    real(wp),dimension(3),intent(out) :: r_target !! apparant state of the target (j2000 frame)
                                                  !! Corrected for one-way light time and stellar aberration
    logical,intent(out) :: status_ok !! true if no problems

    real(wp),dimension(3) :: r_targ_ssb !! target body r wrt. ssb
    real(wp) :: lt !! one-way light time [sec]

    ! Find the geometric position of the target body with respect to the
    ! solar system barycenter. Subtract the position of the observer
    ! to get the relative position. Use this to compute the one-way
    ! light time.
    call eph%get_r(et,b_target,body_ssb,r_targ_ssb,status_ok)
    if (.not. status_ok) return
    r_targ_ssb = r_targ_ssb - rv_obs_ssb(1:3) ! relative pos of target
    lt = norm2(r_targ_ssb) / c_light ! light time

    ! To correct for light time, find the position of the target body
    ! at the current epoch minus the one-way light time. Note that
    ! the observer remains where he is.
    call eph%get_r(et-lt,b_target,body_ssb,r_targ_ssb,status_ok)
    if (.not. status_ok) return
    r_targ_ssb = r_targ_ssb - rv_obs_ssb(1:3)

    ! At this point, r_targ_ssb contains the geometric or light-time
    ! corrected position of the target relative to the observer

    ! stellar aberration correction
    r_target = stellar_aberration(r_targ_ssb,rv_obs_ssb(4:6))

    contains

        function stellar_aberration ( pobj, vobs ) result(appobj)
            !!  Correct the apparent position of an object for stellar aberration.
            !!  see SPICELIB routine `STELAB`

            real(wp),dimension(3),intent(in) :: pobj
            real(wp),dimension(3),intent(in) :: vobs
            real(wp),dimension(3) :: appobj

            real(wp),dimension(3) :: u, vbyc,h
            real(wp) :: lensqr, sinphi, phi
            real(wp),parameter :: zero_tol = epsilon(1.0_wp) !! tolerance for zero

            u = unit(pobj)
            vbyc = vobs / c_light
            lensqr = dot_product ( vbyc, vbyc )
            if ( lensqr >= 1.0_wp) error stop 'velocity > speed of light'
            h = cross(u, vbyc)
            sinphi  = norm2 ( h )
            if ( abs(sinphi) > zero_tol ) then  ! if (sinphi /= 0)
                ! rotate the position of the object by phi
                ! radians about h to obtain the apparent position.
                phi = asin ( sinphi )
                call axis_angle_rotation ( pobj, h, phi, appobj )
            else
                ! observer is moving along the line of sight to the object,
                ! and no correction is required
                appobj = pobj
            end if

        end function stellar_aberration

    end subroutine apparent_position
!********************************************************************************

!********************************************************************************
!>
!  The "circular cubic" shadow model.
!
!### Reference
!  * J. Williams, et. al, "A new eclipse algorithm for use in
!    spacecraft trajectory optimization", 2023, AAS 23-243

    subroutine cubic_shadow_model(rsun, radsun, rplanet, radplanet, sunfrac, rbubble)

    real(wp),dimension(3),intent(in)   :: rsun      !! apparent position vector of sun wrt spacecraft [km]
    real(wp), intent(in)               :: radsun    !! radius of sun [km]
    real(wp),dimension(3), intent(in)  :: rplanet   !! apparent position vector of eclipsing body wrt spacecraft [km]
    real(wp), intent(in)               :: radplanet !! radius of the eclipsing body [km]
    real(wp), intent(out)              :: sunfrac   !! value of the function (>0 no eclipse,
                                                    !! <0 eclipse, =0 on the shadow line)
    real(wp),intent(in),optional       :: rbubble   !! eclipse bubble radius. if present, then `sunfrac` is
                                                    !! the value along an arc length of `rbubble`
                                                    !! in the direction of the max eclipse line.

    real(wp),dimension(3) :: r   !! radius vector from eclipsing body to spacecraft
    real(wp),dimension(3) :: rsb !! radius vector from the sun to the eclipsing body
    real(wp) :: tmp              !! temp value
    real(wp) :: alpha            !! [deg]
    real(wp) :: alpha0           !! [deg]
    real(wp) :: sin_alpha0       !! `sin(alpha0)`
    real(wp) :: rsbmag           !! magnitude of radius vector from the sun to the eclipsing body
    real(wp) :: rmag             !! magnitude of `r`
    logical :: compute_bubble    !! use the `rbubble` inputs to adjust `alpha`

    compute_bubble = present(rbubble)
    if (compute_bubble) compute_bubble = rbubble > zero

    r      = -rplanet
    rmag   = norm2(r)
    if (rmag<radplanet) then
        ! if inside the body, just return value from the surface
        r    = radplanet * unit(r)
        rmag = radplanet
    end if
    rsb        = rplanet - rsun
    alpha      = safe_acosd(dot_product(unit(r),unit(rsb)))
    if (compute_bubble) alpha = rad2deg*abs(wrap_angle(alpha*deg2rad-abs(rbubble)/rmag))
    rsbmag     = norm2(rsb)
    tmp        = (radsun + radplanet) / rsbmag
    sin_alpha0 = (one/rmag)*(radplanet*sqrt((one/tmp)**2-one)+sqrt(rmag**2-radplanet**2))*tmp
    alpha0     = safe_asind(sin_alpha0)
    sunfrac    = (alpha**2/(alpha0**2-alpha0**3/270.0_wp))*(one-alpha/270.0_wp)-one

    contains

        pure real(wp) function safe_asind(x)
            !! `asind` with range checking
            real(wp),intent(in) :: x
            safe_asind = asind(min(one,max(-one,x)))
        end function safe_asind

        pure real(wp) function safe_acosd(x)
            !! `acosd` with range checking
            real(wp),intent(in) :: x
            safe_acosd = acosd(min(one,max(-one,x)))
        end function safe_acosd

    end subroutine cubic_shadow_model
!*****************************************************************************************

!*****************************************************************************************
!>
!  Another eclipse model, using circular area assumptions.
!
!### References
!  * Montenbruck and Gill, "Satellite Orbits".
!  * The GMAT routine `ShadowState::FindShadowState`.

    subroutine solar_fraction_alt(d_s, rs, d_p, rp, percentsun, info)

    real(wp),dimension(3),intent(in) :: d_s !! vector from the spacecraft to the Sun
    real(wp),intent(in)              :: rs  !! radius of the Sun
    real(wp),dimension(3),intent(in) :: d_p !! vector from the spacecraft to the planet
    real(wp),intent(in)              :: rp  !! radius of the planet
    real(wp),intent(out)             :: percentsun !! fraction of the Sun visible [0=total eclipse, 1=no eclipse]
    character(len=:),allocatable,intent(out),optional :: info !! info string

    real(wp),dimension(3) :: unitsun, d_s_hat, d_p_hat
    real(wp) :: rho_s, rho_p, theta, rdotsun, d_s_mag, d_p_mag, c, a2, b2, x, y, area

    !              [sc]
    !            /     \
    !        d_s       d_p
    !      /              \
    !  [sun] ---------- [body]

    unitsun = unit(d_s - d_p) ! body to sun unit vector
    rdotsun = dot_product(-d_p,unitsun)

    if (rdotsun > zero) then
        ! sunny side of central body is always fully lit
        ! [the assumption here is the sun is always bigger than the body?]
        if (present(info)) info = 'full sun'
        percentsun = one
    else

        d_s_mag = norm2(d_s)
        d_p_mag = norm2(d_p)
        if (rs >= d_s_mag) then ! inside the Sun
            if (present(info)) info = 'inside Sun'
            percentsun = one
        else if (rp >= d_p_mag) then ! inside the planet
            if (present(info)) info = 'inside Planet'
            percentsun = zero
        else
            rho_s   = asin(rs/d_s_mag)
            rho_p   = asin(rp/d_p_mag)
            d_p_hat = unit(d_p)
            d_s_hat = unit(d_s)
            theta   = acos(dot_product(d_p_hat,d_s_hat)) ! apparant distance from sun to body

            if (rho_s + rho_p <= theta) then ! full sunlight
                if (present(info)) info = 'full sunlight'
                percentsun = one
            else if (theta <= rho_p-rho_s) then ! umbra
                if (present(info)) info = 'umbra'
                percentsun = zero
            else if ( (abs(rho_s-rho_p)<theta) .and. (theta < rho_s + rho_p) ) then ! penumbra
                if (present(info)) info = 'penumbra'
                ! see montenbruck and gill, eq. 3.87-3.94
                c          = acos(dot_product(d_p_hat,d_s_hat))
                a2         = rho_s*rho_s
                b2         = rho_p*rho_p
                x          = (c*c + a2 - b2) / (two * c)
                y          = sqrt(a2 - x*x)
                area       = a2*acos(x/rho_s) + b2*acos((c-x)/rho_p) - c*y
                percentsun = one - area / (pi * a2)
            else ! antumbra
                if (present(info)) info = 'antumbra'
                percentsun =  one - rho_p*rho_p/(rho_s*rho_s)
            end if
        end if
    end if

    end subroutine solar_fraction_alt
!*****************************************************************************************

!*****************************************************************************************
!>
!  Another eclipse model, using circular area assumptions,
!  coded up based on the nixspace documentation.
!  The results are very similar to `solar_fraction_alt`.
!
!### References
!  * https://nyxspace.com/nyxspace/MathSpec/celestial/eclipse/#nomenclature

    subroutine solar_fraction_alt2(r_l, Rl, r_e, Re, percentsun, info)

    real(wp),dimension(3),intent(in) :: r_l !! vector from the spacecraft to the Sun
    real(wp),intent(in)              :: Rl  !! radius of the Sun
    real(wp),dimension(3),intent(in) :: r_e !! vector from the spacecraft to the planet
    real(wp),intent(in)              :: Re  !! radius of the planet
    real(wp),intent(out)             :: percentsun !! fraction of the Sun visible [0=total eclipse, 1=no eclipse]
    character(len=:),allocatable,intent(out),optional :: info !! info string

    real(wp) :: rlp, rep, dp, r_l_mag, r_e_mag, &
                d1, d2, dp2, rlp2, rep2, At, Astar

    ! this check isn't mentioned in the reference, but needed
    ! for sc -- sun -- body case
    if (dot_product(-r_e,unit(r_l-r_e)) > zero) then
        ! sunny side of body is always fully lit
        ! [the assumption here is the sun is always bigger than the body?]
        if (present(info)) info = 'full sun'
        percentsun = one
        return
    end if

    r_l_mag = norm2(r_l)
    r_e_mag = norm2(r_e)

    ! these checks also aren't in the writeup:
    if (Rl >= r_l_mag) then ! inside the Sun
        if (present(info)) info = 'inside Sun'
        percentsun = one; return
    else if (Re >= r_e_mag) then ! inside the planet
        if (present(info)) info = 'inside Planet'
        percentsun = zero; return
    end if

    rlp = asin(Rl/r_l_mag)
    rep = asin(Re/r_e_mag)
    dp  = acos(dot_product(r_l,r_e)/(r_l_mag*r_e_mag))

    ! modified this check:
    !if (dp-rlp<rep) then  ! original
    if (rlp+rep<=dp) then  ! corrected
        if (present(info)) info = 'full sun'
        percentsun = one ! full sun
    else if (rep>dp+rlp) then
        if (present(info)) info = 'umbra'
        percentsun = zero ! umbra
    else if (rlp-rep>=dp .or. dp>=rlp+rep) then ! antumbra
        if (present(info)) info = 'antumbra'
        percentsun = one - rep*rep/(rlp*rlp)
    else ! penumbra
        if (present(info)) info = 'penumbra'
        dp2   = dp*dp
        rlp2  = rlp*rlp
        rep2  = rep*rep
        d1    = (dp2 - rlp2 + rep2) / (two*dp)
        d2    = (dp2 + rlp2 - rep2) / (two*dp)
        At    = A(rep, rep2, d1) + A(rlp, rlp2, d2)
        Astar = pi*rlp2
        percentsun = (Astar - At) / Astar
    end if

    contains
        pure real(wp) function A(r,r2,d)
            real(wp),intent(in) :: r, r2, d
            A = r2*acos(d/r) - d*sqrt(r2 - d*d)
        end function A
    end subroutine solar_fraction_alt2
!*****************************************************************************************

!*****************************************************************************************
!>
!  Unit tests for the listing module.

    subroutine lighting_module_test()

    real(wp) :: rs, rp
    real(wp),dimension(3) :: d_s, d_p

    rs = 1.0_wp ! sun radius
    rp = 1.0_wp ! planet radius

    ! sun -- body -- sc  -> 0.0
    d_s = [-100.0_wp, 0.0_wp, 0.0_wp]
    d_p = [-10.0_wp, 0.0_wp, 0.0_wp]
    call go()

    ! sc -- sun -- body  -> 1.0
    d_s = [10.0_wp, 0.0_wp, 0.0_wp]
    d_p = [100.0_wp, 0.0_wp, 0.0_wp]
    call go()

    ! sc -- body -- sun  -> 0.0
    d_s = [100.0_wp, 0.0_wp, 0.0_wp]
    d_p = [10.0_wp, 0.0_wp, 0.0_wp]
    call go()

    ! sc -- body -- sun  -> penumbra
    d_s = [100.0_wp, 0.0_wp, 0.0_wp]
    d_p = [10.0_wp, 1.0_wp, 0.0_wp]
    call go()

    ! body -- sc -- sun
    d_s = [-100.0_wp, 0.0_wp, 0.0_wp]
    d_p = [100.0_wp, 0.0_wp, 0.0_wp]
    call go()

    !....................................
    ! sc -- body -- sun  -> antumbra
    rs = 100.0_wp
    d_s = [20000.0_wp, 0.0_wp, 0.0_wp]
    d_p = [400.0_wp,  0.0_wp, 0.0_wp]
    call go()

    rs = 100.0_wp  ! umbra
    d_s = [20000.0_wp, 0.0_wp, 0.0_wp]
    d_p = [100.0_wp,  0.0_wp, 0.0_wp]
    call go()

    ! realistic sun/earth case:
    !  sun -- earth -- sc
    rs = 696000.0_wp
    rp = 6378.0_wp
    d_s = [-149597870.7_wp, 0.0_wp, 0.0_wp]
    d_p = [-6778.0_wp, 6400.0_wp, 0.0_wp]
    call go()

    ! ! an edge case, a very small sun very close to the body on x-axis,
    ! ! sc on y-axis very close to body    .. i don't think any properly handle this .. .double check...
    ! rs = 0.0001_wp ! sun radius
    ! rp = 10.0_wp ! planet radius
    ! d_p = [0.0001_wp, -rp-0.01_wp, 0.0_wp]
    ! d_s = d_p + [-rp-0.01_wp, 0.0_wp, 0.0_wp]
    ! call go()

    contains

        subroutine go()
            real(wp) :: phi1, phi2, phi3
            character(len=:),allocatable :: info1, info2, info3
            print*, '----------------------------------'
            write(*,*) ''
            call solar_fraction(     d_s, rs, d_p, rp, phi1, info1)
            call solar_fraction_alt( d_s, rs, d_p, rp, phi2, info2)
            call solar_fraction_alt2(d_s, rs, d_p, rp, phi3, info3)
            write(*,*) 'phi1 = ', phi1, info1
            write(*,*) 'phi2 = ', phi2, info2
            write(*,*) 'phi3 = ', phi3, info3
            write(*,*) 'diff 1= ', abs(phi1-phi2) ! spherical vs circular
            write(*,*) 'diff 2= ', abs(phi2-phi3) ! two circular models
            if (abs(phi1-phi2)>1.0e-4_wp) error stop 'WARNING: large difference between models'
            print*, ''
        end subroutine go
    end subroutine lighting_module_test
!*****************************************************************************************

!*****************************************************************************************
    end module lighting_module
!*****************************************************************************************