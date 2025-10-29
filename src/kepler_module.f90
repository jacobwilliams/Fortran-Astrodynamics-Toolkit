!*******************************************************************************
!>
!  Kepler propagation routines.

    module kepler_module

    use numbers_module
    use kind_module,     only: wp
    use iso_fortran_env, only: error_unit

    implicit none

    private

    public :: kepler_shepperd
    public :: kepler_goodyear_stienon_klumpp
    public :: kepler_classical

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  Classical Kepler propagator for elliptical and hyperbolic orbits.
!  Uses Lagrange formulations from Battin & Newton's method.
!
!### See also
!  * Dario Izzo: pykep/src/core_functions/propagate_lagrangian.h

    subroutine kepler_classical(x0, dt, mu, xf)

    use newton_module, only: newton

    implicit none

    real(wp),dimension(6),intent(in)  :: x0  !! initial position,velocity vector
    real(wp),intent(in)               :: dt  !! propagation time
    real(wp),intent(in)               :: mu  !! central body gravitational parameter
    real(wp),dimension(6),intent(out) :: xf  !! final position,velocity vector

    real(wp),dimension(3) :: r0  !! initial position vector
    real(wp),dimension(3) :: v0  !! initial velocity vector
    real(wp) :: de     !! eccentric anomaly difference
    real(wp) :: dh     !! hyperbolic anomaly difference
    integer  :: iflag  !! newton status flag
    real(wp) :: r,rmag,vmag,energy,a,sqrta,f,&
                g,ft,gt,sigma0,dm,dn,xs,fx,b,p,x,z,term

    real(wp),parameter :: parabolic_tol = 1.0e-12_wp !! zero tol for parabolic orbits (energy)
    real(wp),parameter :: ftol = 1.0e-12_wp          !! function tol for root finding
    real(wp),parameter :: xtol = 1.0e-12_wp          !! indep variable tol for root finding
    integer,parameter  :: max_iter = 1000            !! maximum number of iterations in newton

    ! check trivial case:
    if (dt==zero) then
        xf = x0
        return
    end if

    r0     = x0(1:3)
    v0     = x0(4:6)
    rmag   = norm2(r0)
    vmag   = norm2(v0)
    energy = (vmag * vmag / two - mu / rmag)
    sigma0 = dot_product(r0,v0) / sqrt(mu)

    ! if not parabolic, then compute semimajor axis
    if (abs(energy) > parabolic_tol) then
        a = -mu / two / energy
    end if

    if (energy < -parabolic_tol) then ! elliptical case

        sqrta = sqrt(a)
        dm = sqrt(mu / a**3) * dt
        de = dm

        call newton(de,kepde_,d_kepde_,ftol,xtol,max_iter,xs,fx,iflag)
        if (iflag<0) then
            write(error_unit,'(A)') 'Error in kepler_classical [elliptical]: newton did not converge'
            write(*,*) xs,fx,iflag
        end if

        de = xs
        r = a + (rmag - a) * cos(de) + sigma0 * sqrta * sin(de) ! eqn 4.42

        ! lagrange coefficients (Battin eqn 4.41)
        f  = one - a / rmag * (one - cos(de))
        g  = a * sigma0 / sqrt(mu) * (one - cos(de)) + rmag * sqrt(a / mu) * sin(de)
        ft = -sqrt(mu * a) / (r * rmag) * sin(de)
        gt = one - a / r * (one - cos(de))

    else if (energy > parabolic_tol) then ! hyperbolic case

        sqrta = sqrt(-a)
        dn = sqrt(-mu / a**3) * dt
        dh = sign(one,dt)      ! todo: need a better initial guess

        call newton(dh,kepdh_,d_kepdh_,ftol,xtol,max_iter,xs,fx,iflag)
        if (iflag<0) then
            write(error_unit,'(A)') 'Error in kepler_classical [hyperbola]: newton did not converge'
            write(*,*) xs,fx,iflag
        end if

        dh = xs
        r = a + (rmag - a) * cosh(dh) + sigma0 * sqrta * sinh(dh)

        ! lagrange coefficients (Battin eqn 4.62)
        f  = one - a / rmag * (one - cosh(dh))
        g  = a * sigma0 / sqrt(mu) * (one - cosh(dh)) + rmag * sqrt(-a / mu) * sinh(dh)
        ft = -sqrt(-mu * a) / (r * rmag) * sinh(dh)
        gt = one - a / r * (one - cosh(dh))

    else  ! parbolic case

        ! See Battin Section 4.2

        p    = two * rmag - sigma0**2
        B    = one / p**1.5_wp * ( sigma0 * (rmag + p) + three * sqrt(mu) * dt )
        term = B + sqrt(one+B*B)
        z    = (term)**(one/three) - (term)**(-one/three) ! Battin eqn 4.12 where z = tan(f/2)
        x    = sqrt(p) * z - sigma0
        r    = rmag + sigma0 * x + one/two * x**2

        f  = one - x**2 / two / rmag
        g  = x / two / sqrt(mu) * ( two * rmag + sigma0 * x )
        ft = - sqrt(mu) * x / rmag / r
        gt = one - x**2 / two / r

    end if

    ! results:
    xf(1:3) = f  * r0 + g  * v0
    xf(4:6) = ft * r0 + gt * v0

    contains

        ! function wrappers for newtons method:

        subroutine kepde_(x,f)
        implicit none
        real(wp),intent(in)  :: x  !! de
        real(wp),intent(out) :: f
        f = kepde(x, dm, sigma0, sqrta, a, rmag)
        end subroutine kepde_

        subroutine d_kepde_(x,f)
        implicit none
        real(wp),intent(in)  :: x  !! de
        real(wp),intent(out) :: f
        f = d_kepde(x, sigma0, sqrta, a, rmag)
        end subroutine d_kepde_

        subroutine kepdh_(x,f)
        implicit none
        real(wp),intent(in)  :: x  !! dh
        real(wp),intent(out) :: f
        f = kepdh(x, dn, sigma0, sqrta, a, rmag)
        end subroutine kepdh_

        subroutine d_kepdh_(x,f)
        implicit none
        real(wp),intent(in)  :: x  !! dh
        real(wp),intent(out) :: f
        f = d_kepdh(x, sigma0, sqrta, a, rmag)
        end subroutine d_kepdh_

    end subroutine kepler_classical
!*******************************************************************************

!*******************************************************************************
!>
!  Elliptic Kepler's equation

    pure function kepe(e, m, ecc)

    implicit none

    real(wp),intent(in) :: e   !! eccentric anomaly
    real(wp),intent(in) :: m   !! mean anomaly
    real(wp),intent(in) :: ecc !! eccentricity
    real(wp) :: kepe

    kepe = e - ecc * sin(e) - m

    end function kepe
!*******************************************************************************

!*******************************************************************************
!>
!  Derivative of [[kepe]] w.r.t. `e`

    pure function d_kepe(e, ecc)

    implicit none

    real(wp),intent(in) :: e   !! eccentric anomaly
    real(wp),intent(in) :: ecc !! eccentricity
    real(wp) :: d_kepe

    d_kepe = one - ecc * cos(e)

    end function d_kepe
!*******************************************************************************

!*******************************************************************************
!>
!  Elliptic Kepler's equation written in terms of the
!  eccentric anomaly difference.  See Battin, eqn 4.43.

    pure function kepde(de, dm, sigma0, sqrta, a, r)

    implicit none

    real(wp),intent(in) :: de      !! eccentric anomaly difference
    real(wp),intent(in) :: dm      !! mean anomaly difference
    real(wp),intent(in) :: sigma0
    real(wp),intent(in) :: sqrta
    real(wp),intent(in) :: a
    real(wp),intent(in) :: r
    real(wp) :: kepde

    kepde = -dm + de + sigma0 / sqrta * (one - cos(de)) - &
            (one - r / a) * sin(de)

    end function kepde
!*******************************************************************************

!*******************************************************************************
!>
!  Derivative of [[kepde]] w.r.t `de`.

    pure function d_kepde(de, sigma0, sqrta, a, r)

    implicit none

    real(wp),intent(in) :: de      !! eccentric anomaly difference
    real(wp),intent(in) :: sigma0
    real(wp),intent(in) :: sqrta
    real(wp),intent(in) :: a
    real(wp),intent(in) :: r
    real(wp) :: d_kepde

    d_kepde = one + sigma0 / sqrta * sin(de) - (one - r / a) * cos(de)

    end function d_kepde
!*******************************************************************************

!*******************************************************************************
!>
!  Battin, eqn. 4.64.

    pure function kepdh(dh, dn, sigma0, sqrta, a, r)

    implicit none

    real(wp),intent(in) :: dh      !! hyperbolic anomaly difference
    real(wp),intent(in) :: dn
    real(wp),intent(in) :: sigma0
    real(wp),intent(in) :: sqrta
    real(wp),intent(in) :: a
    real(wp),intent(in) :: r
    real(wp) :: kepdh

    kepdh = -dn - dh + sigma0 / sqrta * (cosh(dh) - one) + &
            (one - r / a) * sinh(dh)

    end function kepdh
!*******************************************************************************

!*******************************************************************************
!>
!  Derivative of [[kepdh]] w.r.t `dh`.

    pure function d_kepdh(dh, sigma0, sqrta, a, r)

    implicit none

    real(wp),intent(in) :: dh      !! hyperbolic anomaly difference
    real(wp),intent(in) :: sigma0
    real(wp),intent(in) :: sqrta
    real(wp),intent(in) :: a
    real(wp),intent(in) :: r
    real(wp) :: d_kepdh

    d_kepdh = -one + sigma0 / sqrta * sinh(dh) + (one - r / a) * cosh(dh)

    end function d_kepdh
!*******************************************************************************

!*******************************************************************************
!>
!  Barker time of flight equation

    pure function barker(r1, r2, mu)

    implicit none

    real(wp),dimension(3),intent(in) :: r1
    real(wp),dimension(3),intent(in) :: r2
    real(wp),intent(in) :: mu
    real(wp) :: barker

    real(wp) :: x, r1mag, r2mag, r21mag, sigma, r1pr2mag
    real(wp),dimension(3) :: r21

    r1mag    = norm2(r1)
    r2mag    = norm2(r2)
    r21      = r2 - r1
    r21mag   = norm2(r21)
    x        = r1(1) * r2(2) - r1(2) * r2(1)
    sigma    = sign(one,x)
    r1pr2mag = r1mag + r2mag

    barker = (r1pr2mag + r21mag)**1.5_wp - &
             sigma * (r1pr2mag - r21mag)**1.5_wp / (six * sqrt(mu))

    end function barker
!*******************************************************************************

!*******************************************************************************
!>
!
    pure function kepds(ds, dt, r0, vr0, alpha, mu)

    implicit none

    real(wp),intent(in) :: ds  !! universal anomaly difference
    real(wp),intent(in) :: dt
    real(wp),intent(in) :: r0
    real(wp),intent(in) :: vr0
    real(wp),intent(in) :: alpha
    real(wp),intent(in) :: mu
    real(wp) :: kepds

    real(wp) :: c,s,ads2,ds2

    ds2  = ds * ds
    ads2 = alpha * ds2
    s    = stumpff_s(ads2)
    c    = stumpff_c(ads2)

    kepds = -sqrt(mu) * dt + r0 * vr0 * ds2 * c / sqrt(mu) + &
            (one - alpha * r0) * ds2 * ds * s + r0 * ds

    end function kepds
!*******************************************************************************

!*******************************************************************************
!>
!
    pure function d_kepds(ds, r0, vr0, alpha, mu)

    implicit none

    real(wp),intent(in) :: ds
    real(wp),intent(in) :: r0
    real(wp),intent(in) :: vr0
    real(wp),intent(in) :: alpha
    real(wp),intent(in) :: mu
    real(wp) :: d_kepds

    real(wp) :: c,s,ads2,ds2

    ds2  = ds * ds
    ads2 = alpha * ds2
    s    = stumpff_s(ads2)
    c    = stumpff_c(ads2)

    d_kepds = r0 * vr0 / sqrt(mu) * ds * (one - ads2 * s) + &
              (one - alpha * r0) * ds2 * c + r0

    end function d_kepds
!*******************************************************************************

!*******************************************************************************
!>
!  Stumpff function S(z)

    pure function stumpff_s(z) result(s)

    implicit none

    real(wp),intent(in) :: z
    real(wp) :: s

    if (z > zero) then
        s = (sqrt(z) - sin(sqrt(z))) / sqrt(z)**3
    else if (z < zero) then
        s = (sinh(sqrt(-z)) - sqrt(-z)) / (-z)**(three/two)
    else
        s = one/six
    end if

    end function stumpff_s
!*******************************************************************************

!*******************************************************************************
!>
!  Stumpff function C(z)

    pure function stumpff_c(z) result(c)

    implicit none

    real(wp),intent(in) :: z
    real(wp) :: c

    if (z > zero) then
        c = (one - cos(sqrt(z))) / z
    else if (z < zero) then
        c = (cosh(sqrt(-z)) - one) / (-z)
    else
        c = 0.5_wp
    end if

    end function stumpff_c
!*******************************************************************************

!*******************************************************************************
!>
!  Kepler propagation using Shepperd's method.
!
!### Reference
!  * S.W. Shepperd, "Universal Keplerian State Transition Matrix".
!    Celestial Mechanics 35(1985) p. 129-144.
!  * Rody P.S. Oldenhuis, `progress_orbitM` Matlab function (BSD license).
!    http://www.mathworks.com/matlabcentral/fileexchange/26349-kepler-state-transition-matrix-mex
!  * C.D. Eagle, Orbital Mechanics with MATLAB, `twobody2.m` Matlab function (BSD license).
!    http://www.mathworks.com/matlabcentral/fileexchange/48723-matlab-functions-for-two-body-orbit-propagation

    subroutine kepler_shepperd(mu,rv1,dt,rv2,istat)

    implicit none

    real(wp),intent(in)               :: mu    !! gravitational parameter
    real(wp),dimension(6),intent(in)  :: rv1   !! initial position,velocity vector
    real(wp),intent(in)               :: dt    !! time step
    real(wp),dimension(6),intent(out) :: rv2   !! final position,velocity vector
    integer,intent(out),optional      :: istat !! status flag (if not present, warnings are printed):
                                               !! Linear combination of :
                                               !!
                                               !! `0` : all is well
                                               !! `-10` : failed to converge in time loop
                                               !! `-100` : failed to converge in `g` loop

    ! note: some of these could also be (optional) input:
    real(wp),parameter :: min_dt     = epsilon(one) !! time step considered as zero (sec)
    integer,parameter  :: max_iter   = 1000         !! max iterations for time
    integer,parameter  :: max_iter_g = 1000         !! max iterations for g
    real(wp),parameter :: ttol       = 1.0e-8_wp    !! tolerance for time (sec)
    real(wp),parameter :: gtol       = 1.0e-12_wp   !! tolerance for g
    real(wp),parameter :: zero_tol   = epsilon(one) !! tolerance for beta=0 (parabola)
    logical,parameter  :: use_halley = .true.       !! use the Halley update
                                                    !! rather than the original Newton

    real(wp),dimension(3) :: r1  !! initial position vector
    real(wp),dimension(3) :: v1  !! initial velocity vector
    real(wp) :: r1mag,nu0,beta,p,deltau,u,t,deltat,bu,q,&
                a,b,gprev,u0w2,u1w2,uu,u0,u1,u2,u3,r,&
                ff,f,gg,g,umin,umax,du,abs_beta
    integer :: n,iter,k,d,l,iterg

    if (present(istat)) istat = 0

    if (abs(dt) <= min_dt) then

        rv2 = rv1

    else

        r1 = rv1(1:3)
        v1 = rv1(4:6)
        r1mag = norm2(r1)
        nu0 = dot_product(r1,v1)
        beta = two*mu/r1mag - dot_product(v1,v1)
        abs_beta = abs(beta)

        if (abs_beta>zero_tol) then
            umax = one/sqrt(abs_beta)
        else
            umax = huge(one)
        end if
        umin = -umax

        if (beta > zero) then
            p = twopi*mu*beta**(-three/two)
            n = floor((dt + p/two - two*nu0/beta)/p)
            deltau = twopi*n*beta**(-five/two)
        else
            deltau = zero
        end if

        u = zero
        t = zero
        iter = 0
        deltat = dt-t
        do while (abs(deltat) > ttol)
            iter = iter + 1
            bu = beta*u*u
            q = min(0.5_wp,bu/(one+bu))  ! avoid q > 0.5 due to numerical issues
            a = one
            b = one
            g = one
            n = 0
            k = -9
            d = 15
            l = 3
            gprev = huge(one)
            iterg = 0
            do while (abs(g-gprev) > gtol)
                iterg = iterg + 1
                k = -k
                l = l + 2
                d = d + 4*l
                n = n + (1+k)*l
                a = d/(d - n*a*q)
                b = (a-one)*b
                gprev = g
                g = g + b
                if (iterg==max_iter_g) then
                    if (present(istat)) then
                        istat = istat - 100
                    else
                        write(*,*) 'Warning: kepler_shepperd failed to converge in g iteration'
                    end if
                    exit
                end if
            end do

            u0w2 = one - two*q
            u1w2 = two*(one-q)*u
            uu = 16.0_wp/15.0_wp*u1w2**5*g + deltau
            u0 = two*u0w2**2-one
            u1 = two*u0w2*u1w2
            u2 = two*u1w2**2
            u3 = beta*uu + u1*u2/three
            r = r1mag*u0 + nu0*u1 + mu*u2
            t = r1mag*u1 + nu0*u2 + mu*u3
            deltat = dt - t

            if (use_halley) then
                ! Halley version from Oldenhuis routine
                du = deltat/((one-q)*(four*r + deltat*beta*u))
            else
                ! original Newton
                du = deltat/four/(one-q)/r
            end if

            ! this logic is from the Eagle routine
            if (du<zero) then
                umax = u
                u = u + du
                if (u<umin) u = (umin + umax)/two
            else
                umin = u
                u = u + du
                if (u>umax) u = (umin + umax)/two
            end if

            if (iter==max_iter) then
                if (abs(deltat) > ttol) then
                    ! it still hasn't converged
                    if (present(istat)) then
                        istat = istat - 10
                    else
                        write(*,*) 'Warning: kepler_shepperd failed to converge in time iteration'
                    end if
                    exit
                end if
            end if

        end do

        ff = one - mu/r1mag*u2
        f  = -mu*u1/r/r1mag
        gg = r1mag*u1 + nu0*u2
        g  = one - mu/r*u2

        rv2 = [r1*ff + v1*gg, r1*f + v1*g]

    end if

    end subroutine kepler_shepperd
!*******************************************************************************

!*******************************************************************************
!>
!  Kepler propagator based on the Goodyear code with
!  modifications by Stienon and Klumpp.
!
!### See also
!  * W. H. Goodyear, "Completely General Closed-Form Solution for Coordinates
!    and Partial Derivatives of the Two-Body Problem", Astronomical Journal,
!    Vol. 70, No. 3, April 1965.
!    [pdf](http://adsabs.harvard.edu/full/1965AJ.....70..189G)
!  * W. H. Goodyear, "A General Method for the Computation of Cartesian
!    Coordinates and Partial Derivatives of the Two-Body Problem",
!    NASA CR-522, Sep 1, 1966.
!    [pdf](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19660027556.pdf)
!  * A. Klumpp, "Performance Comparision of Lambert and Kepler Algorithms",
!    JPL Interoffice Memorandum, 314.1-0426-ARK, Jan 2, 1991.
!   (See the `KEPGSK` Fortran 77 routine)

    subroutine kepler_goodyear_stienon_klumpp(rv0,tau,mu,accy,rvf)

    implicit none

    real(wp),dimension(6),intent(in)  :: rv0   !! state vector at reference
                                               !! time T0 [km,km/s]
    real(wp),intent(in)               :: tau   !! Time interval T-T0 [sec]
    real(wp),intent(in)               :: mu    !! Central body gravitational
                                               !! constant [km^3/s^2]
    real(wp),intent(in)               :: accy  !! Fractional accuracy required
                                               !! [0->0.1]
    real(wp),dimension(6),intent(out) :: rvf   !! state vector at solution time
                                               !! T [km,km/s]

    real(wp),parameter :: onethird  = 1.0_wp / 3.0_wp !! 0.333333 in original code

    integer  :: iloop             !! loop counter defining number of iterations
    real(wp) :: r                 !! radius at solution time `T`
    real(wp) :: r0                !! radius at solution time `T0`
    integer  :: i                 !! index
    integer  :: j                 !! index
    integer  :: m                 !! counter for normalizing lambda
    real(wp) :: psi               !! independent variable in Kepler's equation
    real(wp) :: psin              !! lower limit to `psi`
    real(wp) :: psip              !! upper limit to `psi`
    real(wp) :: dlim              !! tau/periapsis distance = abs(limit to `psi`)
    real(wp) :: dtau              !! Kepler's equation residual
    real(wp) :: dtaun             !! lower limit to `dtau`
    real(wp) :: dtaup             !! upper limit to `dtau`
    real(wp) :: accrcy            !! accuracy = midval(zero, `accy`, 0.1)
    real(wp) :: sig0              !! dot product: `r` with `v`
    real(wp) :: alpha             !! twice the energy per unit mass
    real(wp) :: h2                !! square of angular momentum per unit mass
    real(wp) :: a                 !! reduced argument of `s_0`, `s_1`, etc.
    real(wp) :: ap                !! actual argument of `s_0`, `s_1`, etc.
    real(wp) :: c0                !! c(n) is the nth Goodyear polynomial `s_n`
                                  !! divided by its leading term, so that
                                  !!
                                  !! * `c0 = s_0 = cosh(z)`
                                  !! * `c1 = s_1/psi = sinh(z)/z`
                                  !! * `c(n+2) = (n+1)(n+2)[c(n)-1]/z^2`
                                  !!
                                  !! where `z = psi*sqrt(alpha)`. `ap=norm2(z)`
    real(wp) :: c1
    real(wp) :: c2
    real(wp) :: c3
    real(wp) :: c4
    real(wp) :: c5x3
    real(wp) :: s1                !! `s_n` is the nth Goodyear
                                  !! polynomial, given by
                                  !! `s_n = psi^n sim{z^(2k)/(2k)!} k >=0 `
    real(wp) :: s2
    real(wp) :: s3
    real(wp) :: g                 !! the G-function
    real(wp) :: gdm1              !! time derivative of the G-function minus one
    real(wp) :: fm1               !! the F-function minus one
    real(wp) :: fd                !! time derivative of the F-function
    real(wp),dimension(4) :: c    !! terms in F(psi)

    accrcy = min(max(zero, accy), 0.1_wp)

    r0    = dot_product(rv0(1:3), rv0(1:3)) ! r dot r
    sig0  = dot_product(rv0(1:3), rv0(4:6)) ! r dot v
    alpha = dot_product(rv0(4:6), rv0(4:6)) ! v dot v
    h2    = max(r0*alpha - sig0*sig0, zero)
    r0    = sqrt(r0)
    alpha = alpha - two*mu/r0

    !compute initial limits for phi and dtau:
    c0 = sqrt(max(mu*mu + h2*alpha, zero))    ! gm * eccentricity
    dlim = div(1.1_wp * tau * (c0 + mu), h2)  ! tau/periapsis distance = max(a).
                                              ! Arbitrary factor 1.1 ensures
                                              ! psin <= psi <= psip
                                              ! when converged, despite
                                              ! numerical imprecision.

    if (tau < zero) then
        psin  = dlim
        psip  = zero
        dtaun = psin
        dtaup = -tau
    else
        psin  = zero
        psip  = dlim
        dtaun = -tau
        dtaup = psip
    end if

    !compute initial value of psi:
    psi = tau/r0

    if (alpha >= zero) then
        c2 = one
        if (tau < zero) c2 = -one
        c3 = abs(tau)
        if (alpha > zero) then
            c0 = h2*sqrt(h2)/(mu + c0)**2
            if (c3 > c0) then
                c1 = sqrt(alpha)
                psi = log(two*alpha*c1*c3/(r0*alpha + c2*c1*sig0 + mu))
                psi = max(psi,one)*c2/c1
            end if
        else    !parabola
            c0 = r0*sqrt(six*r0/mu)
            if (c3 > c0) psi = c2*(six*c3/mu)**onethird
        end if
    end if

    ! begin loop for solving Kepler's equation:

    m = 0
    iloop = 0

    do
        iloop = iloop + 1

        ! begin series summation:

        !compute argument a in reduced series obtained by factoring out psi's:
        a = alpha*psi*psi
        if (abs(a)>one) then
            !save a in ap and mod a if a exceeds unity in magnitude
            ap = a
            do
                m = m + 1
                a = a*0.25_wp
                if (abs(a)<=one) exit
            end do
        end if

        !sum series c5x3=3*s5/psi**5 and c4=s4/psi**4
        c4 = one
        c5x3 = one
        do i=9,3,-1
            j = 2*i
            c5x3 = one + a*c5x3/(j*(j+1))
            c4   = one + a*c4  /(j*(j-1))
        end do

        c5x3 = c5x3/40.0_wp
        c4   = c4  /24.0_wp

        !compute series c3=s3/psi**3,c2=s2/psi**2,c1=s1/psi,c0=rv0
        c3 = (0.5_wp + a*c5x3)/ three
        c2 = 0.5_wp + a*c4
        c1 = one + a*c3
        c0 = one + a*c2
        if (m>0) then

            !demod series c0 and c1 if necessary with double angle formulas:
            do
                c1 = c1*c0
                c0 = two*c0*c0 - one
                m = m - 1
                if (m<=0) exit
            end do

            !compute c2,c3,c4,c5x3 from c0,c1,ap if demod required:
            c2 = (c0 - one)/ap
            c3 = (c1 - one)/ap
            c4 = (c2 - 0.5_wp)/ap
            c5x3 = (three*c3 - 0.5_wp)/ap

        end if

        !compute series s1,s2,s3 from c1,c2,c3:
        s1 = c1*psi
        s2 = c2*psi*psi
        s3 = c3*psi*psi*psi

        ! compute slope r and residuals for kepler's equation:

        c(1) = r0*s1
        c(2) = sig0*s2
        c(3) = mu*s3
        c(4) = tau
        g    = c(4) - c(3)
        dtau = c(1) + c(2) - g
        r    = abs(r0*c0 + (sig0*s1 + mu*s2))

        ! compute next psi:
        do

            !method = 0
            if (dtau < zero) then
                psin  = psi
                dtaun = dtau
                psi   = psi - dtau/r
                if (psi < psip) exit    ! < (not <=) to avoid false convergence
            else
                psip  = psi
                dtaup = dtau
                psi   = psi - dtau/r
                if (psi > psin) exit    ! > (not >=) to avoid false convergence
            end if

            !reset psi within bounds psin and psip:

            !try incrementing bound with dtau nearest zero by the ratio 4*dtau/tau
            !--method = 1
            if (abs(dtaun) < abs(dtaup)) psi = psin*(one-div(four*dtaun,tau))
            if (abs(dtaup) < abs(dtaun)) psi = psip*(one-div(four*dtaup,tau))
            if (psi_in(psi)) exit

            !try doubling bound closest to zero:
            !--method = 2
            if (tau > zero) psi = psin + psin
            if (tau < zero) psi = psip + psip
            if (psi_in(psi)) exit

            !try interpolating between bounds:
            !--method = 3
            psi = psin + (psip - psin) * div(-dtaun, dtaup - dtaun)
            if (psi_in(psi)) exit

            !try halving between bounds:
            !--method = 4
            psi = psin + (psip - psin) * 0.5_wp
            exit

        end do

        ! test for convergence
        i = 1
        do j=2,4
            if (abs(c(j)) > abs(c(i))) i = j
        end do

        if (abs(dtau)>abs(c(i))*accrcy .and. &
            psip-psin>abs(psi)*accrcy .and. &
            psi/=psin .and. psi/=psip ) then
            cycle
        else
            exit
        end if

    end do

    !compute remaining three of four functions: g, gdm1, fm1, fd
    gdm1 = -mu*s2/r
    fm1  = -mu*s2/r0
    fd   = -mu*s1/r0/r

    ! compute state at time T = T0 + TAU
    rvf(1:3) = rv0(1:3) + fm1*rv0(1:3) + g*rv0(4:6)   ! terminal positions
    rvf(4:6) = rv0(4:6) + fd*rv0(1:3) + gdm1*rv0(4:6) ! terminal velocities

    contains
!*******************************************************************************

    !***************************************************************************
    !>
    !  Returns a nonoverflowing quotient

        pure function div(num,den)

        implicit none

        real(wp)            :: div
        real(wp),intent(in) :: num
        real(wp),intent(in) :: den

        real(wp),parameter :: fsmall  = 1.0e-38_wp  !! small number above underflow limit
        real(wp),parameter :: flarge  = 1.0e+38_wp  !! large number below underflow limit

        div = num / ( sign(one,den) * max(abs(den),abs(num)/flarge, fsmall) )

        end function div
    !***************************************************************************

    !***************************************************************************
    !>
    !  Returns true if `phi` is within and not on bounds

        pure function psi_in(psi)

        implicit none

        logical :: psi_in
        real(wp),intent(in) :: psi

        ! > & < (not >= & <=) to avoid false convergence
        psi_in = (psi > psin .and. psi < psip)

        end function psi_in
    !***************************************************************************

    end subroutine kepler_goodyear_stienon_klumpp
!*******************************************************************************

!*******************************************************************************
    end module kepler_module
!*******************************************************************************
