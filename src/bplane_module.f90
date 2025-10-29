!*****************************************************************************************
!> author: Jacob Williams
!
!  B-plane and hyperbolic routines.

    module bplane_module

    use kind_module, only: wp
    use numbers_module
    use vector_module

    implicit none

    private

    public :: bplane
    public :: calculate_bplane_data
    public :: hyperbolic_turning_angle
    public :: vinf_to_energy

    !unit tests:
    public :: bplane_test

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute B-plane parameters from position and velocity.
!
!# References
!  1. W. Kizner, "A method of describing miss distances for lunar and interplanetary trajectories",
!     Planetary and Space Science, Volume 7, July 1961.
!  2. W. Kizner, "[Some orbital elements useful in space trajectory calculations](http://www.dtic.mil/dtic/tr/fulltext/u2/263968.pdf)",
!     JPL Technical Release No. 34-84, July 25, 1960.
!  3. A. B. Sergeyevsky, G. C. Snyder, R. A. Cunniff,
!     "[Interplanetary Mission Design Handbook, Volume I, Part 2](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19840010158.pdf)",
!     JPL Publication 82-43, September 15, 1983.
!  4. [B-Plane Targeting](http://www.agi.com/resources/help/online/stk/10.1/index.html?page=source%2Fextfile%2Fgator%2Feq-bplane.htm)

    subroutine bplane(mu,rv,vinfvec,bmag,theta,BdotT,BdotR,status_ok)

    implicit none

    real(wp),intent(in)               :: mu        !! central body grav parameter \( (km^3/s^2) \)
    real(wp),dimension(6),intent(in)  :: rv        !! state vector (km,km/s)
    real(wp),dimension(3),intent(out) :: vinfvec   !! incoming V-infinity vector (km/s)
    real(wp),intent(out)              :: bmag      !! magnitude of B vector (km)
    real(wp),intent(out)              :: theta     !! aim point orientation [rad]
    real(wp),intent(out)              :: BdotT     !! \( \mathbf{B} \cdot \mathbf{T} \) (km)
    real(wp),intent(out)              :: BdotR     !! \( \mathbf{B} \cdot \mathbf{R} \) (km)
    logical,intent(out)               :: status_ok !! false if there were errors (non-hyperbolic or degenerate state)

    !local variables:
    real(wp),dimension(3) :: r,v,h,evec,bvec,tvec,ehat,hhat,hehat,Shat,That,Rhat,Bhat
    real(wp) :: rmag,vmag,vinf2,rdv,vmag2,a,alpha,sd2,cd2,ct,st,vinf,e

    status_ok = .true.

    r         = rv(1:3)
    v         = rv(4:6)
    h         = cross(r,v)
    hhat      = unit(h)

    if (all(hhat==zero)) then
        write(*,*) 'error: degenerate state.'
        status_ok = .false.
        return
    end if

    rdv       = dot_product(r,v)
    rmag      = norm2(r)
    vmag      = norm2(v)
    vmag2     = vmag*vmag
    vinf2     = vmag2 - two*mu/rmag
    vinf      = sqrt(vinf2)                     ! magnitude of v-infinity vector
    evec      = cross(v,h)/mu - r/rmag          ! eccentricity vector
    e         = norm2(evec)                     ! eccentricity

    if (e<=one) then
        write(*,*) 'error: state is not hyperbolic.'
        status_ok = .false.
        return
    end if

    ehat      = evec/e                          ! eccentricity unit vector
    a         = one / (two/rmag - vmag2/mu)     ! semi-major axis
    hehat     = cross(hhat,ehat)                ! h x e unit vector
    sd2       = one/e                           ! sin(delta/2)
    cd2       = sqrt(one-sd2*sd2)               ! cos(delta/2)
    Shat      = cd2*hehat + sd2*ehat            ! incoming vinf unit vector
    !Shat     = cd2*hehat - sd2*ehat            ! outgoing vinf unit vector
    That      = ucross(Shat,[zero, zero, one])  ! here we define Tvec relative to the Z-axis of
                                                ! the frame in which the state is defined

    if (all(That==zero)) then
        write(*,*) 'error: vinf vector is parallel to z-axis.'
        status_ok = .false.
        return
    end if

    Rhat      = cross(Shat,That)             ! Eqn 1 in [1]
    Bhat      = ucross(Shat,h)
    ct        = dot_product(Bhat,That)       ! cos(theta)
    st        = dot_product(Bhat,Rhat)       ! sin(theta)

    !outputs:
    Bmag      = abs(a)*sqrt( e*e - one )     ! magnitude of B vector
    theta     = atan2(st,ct)                 ! aim point orientation
    vinfvec   = vinf * Shat                  ! incoming vinf vector
    BdotT     = bmag*ct                      ! B dot T
    BdotR     = bmag*st                      ! B dot R

    end subroutine bplane
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute B-plane parameters from position and velocity -- alternate version.
!
!### See also
!  * [[bplane]]
!
!### Reference
!  * This one is based on the algorithm in GMAT: GmatCalcUtil::CalculateBPlaneData

    subroutine calculate_bplane_data(mu,state,bdotr,bdott,bmag,theta,istat)

    implicit none

    real(wp),intent(in)              :: mu     !! central body grav parameter \( (km^3/s^2) \)
    real(wp),dimension(6),intent(in) :: state  !! state vector (km,km/s)
    real(wp),intent(out)             :: bdotr  !! \( \mathbf{B} \cdot \mathbf{R} \) (km)
    real(wp),intent(out)             :: bdott  !! \( \mathbf{B} \cdot \mathbf{T} \) (km)
    real(wp),intent(out)             :: bmag   !! magnitude of B vector (km)
    real(wp),intent(out)             :: theta  !! aim point orientation [rad]
    integer,intent(out)              :: istat  !! status flag:
                                               !!
                                               !! * 0 if no errors.
                                               !! * -1 if state is not hyperbolic

    real(wp),dimension(3) :: r        !! position vector
    real(wp),dimension(3) :: v        !! velocity vector
    real(wp),dimension(3) :: evec     !! eccentricity vector
    real(wp),dimension(3) :: hvec     !! angular momentum vector
    real(wp),dimension(3) :: nvec     !! orbit normal vector
    real(wp),dimension(3) :: svec     !! incoming asymptote
    real(wp),dimension(3) :: bvec     !! B vector
    real(wp),dimension(3) :: tvec     !! T vector
    real(wp),dimension(3) :: rvec     !! R vector
    real(wp)              :: rmag     !! magnitude of `r`
    real(wp)              :: vmag     !! magnitude of `v`
    real(wp)              :: e        !! eccentricity
    real(wp)              :: hmag     !! magnitude of `hvec` vector
    real(wp)              :: b        !! semiminor axis
    real(wp)              :: oneovere !! 1/e
    real(wp)              :: temp

    r    = state(1:3)
    v    = state(4:6)
    hvec = cross(r, v)
    rmag = norm2(r)
    vmag = norm2(v)
    evec = cross(v,hvec)/mu - r/rmag
    e    = norm2(evec)

    if (e <= 1.0) then ! not hyperbolic
        write(*,*) 'error: state is not hyperbolic.'
        istat    = -1
        bdotr    = 0.0_wp
        bdott    = 0.0_wp
        bmag     = 0.0_wp
        theta    = 0.0_wp
    else
        istat    = 0
        evec     = unit(evec)
        hmag     = norm2(hvec)
        hvec     = unit(hvec)
        nvec     = cross(hvec, evec)
        b        = (hmag*hmag) / (mu * sqrt(e*e - 1.0_wp))
        oneovere = 1.0_wp/e
        temp     = sqrt(1.0_wp - oneovere*oneovere)
        svec     = (evec/e) + (temp*nvec)
        bvec     = b * (temp * evec - oneovere*nvec)
        tvec     = [svec(2), -svec(1), 0.0_wp] / sqrt(svec(1)*svec(1) + svec(2)*svec(2))
        rvec     = cross(svec, tvec)
        bdott    = dot_product(bvec, tvec)
        bdotr    = dot_product(bvec, rvec)
        bmag     = sqrt(bdott*bdott + bdotr*bdotr)
        theta    = atan(bdotr, bdott)
    end if

    end subroutine calculate_bplane_data
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute the incoming and/or outgoing v-infinity vectors, given
!  the position and velocity of a hyperbola.
!
!@note This is for testing the other routines.

    subroutine compute_vinf_vectors(mu,rv,vinfin,vinfout)

    implicit none

    real(wp),intent(in)                        :: mu      !! central body gravitational parameter
    real(wp),dimension(6),intent(in)           :: rv      !! position,velocity vector
    real(wp),dimension(3),intent(out),optional :: vinfin  !! incoming v-infinity vector
    real(wp),dimension(3),intent(out),optional :: vinfout !! outgoing v-infinity vector

    real(wp),dimension(3) :: h,e,p,q,r,v
    real(wp) :: rmag,vmag,vinfmag,emag,qmag,cbeta,sbeta

    if (present(vinfin) .or. present(vinfout)) then
        r       = rv(1:3)
        v       = rv(4:6)
        rmag    = norm2(r)
        vmag    = norm2(v)
        h       = cross(r,v)
        q       = cross(v,h)
        vinfmag = sqrt(vmag*vmag-two*mu/rmag)
        e       = q/mu-r/rmag
        emag    = norm2(e)
        q       = cross(h,e)
        qmag    = norm2(q)
        cbeta   = one/emag
        sbeta   = sqrt(one-cbeta*cbeta)
        p       = e/emag
        q       = q/qmag
        if (present(vinfin))  vinfin  = vinfmag*( cbeta*p+sbeta*q)
        if (present(vinfout)) vinfout = vinfmag*(-cbeta*p+sbeta*q)
    end if

    end subroutine compute_vinf_vectors
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 3/6/2016
!
!  Compute the hyperbolic turning angle from the eccentricity.

    pure function hyperbolic_turning_angle(e) result(delta)

    implicit none

    real(wp),intent(in) :: e       !! eccentricity [--]
    real(wp)            :: delta   !! turning angle [rad]

    delta = two*asin(one/e)

    end function hyperbolic_turning_angle
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 3/6/2016
!
!  Convert V-infinity magnitude to energy.

    pure function vinf_to_energy(vinfmag) result(energy)

    implicit none

    real(wp),intent(in) :: vinfmag  !! \(v^{\infty} \) vector magnitude [km/s]
    real(wp)            :: energy   !! two-body orbital energy [km^2/s^2]

    energy = (vinfmag**2)/two

    end function vinf_to_energy
!*****************************************************************************************

!*****************************************************************************************
!>
!  Unit test for [[bplane_module]].

    subroutine bplane_test()

    implicit none

    real(wp),dimension(3) :: vinfvec
    real(wp) :: bmag,theta,BdotT,BdotR,bmag2,theta2,BdotT2,BdotR2
    logical :: status_ok
    integer :: istat

    real(wp),parameter :: mu = 0.398600436233000e+06_wp  !! grav. param. for earth \( (km^2/s^2) \)
    real(wp),dimension(6),parameter :: rv = [ -1.518170076605391E+04_wp, &
                                               1.518170076605392E+04_wp, &
                                               2.147016712324346E+04_wp, &
                                              -4.613927243557805E+00_wp, &
                                              -2.352686048227026E+00_wp, &
                                               1.598938983116768E+00_wp  ] !! example hyperbolic state

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' bplane_test'
    write(*,*) '---------------'
    write(*,*) ''

    call compute_vinf_vectors(mu,rv,vinfin=vinfvec)

    write(*,*) 'vinfvec(1) =' , vinfvec(1)
    write(*,*) 'vinfvec(2) =' , vinfvec(2)
    write(*,*) 'vinfvec(3) =' , vinfvec(3)


    call bplane(mu,rv,vinfvec,bmag,theta,BdotT,BdotR,status_ok)

    write(*,*) ''
    write(*,*) 'vinfvec(1) =' , vinfvec(1)
    write(*,*) 'vinfvec(2) =' , vinfvec(2)
    write(*,*) 'vinfvec(3) =' , vinfvec(3)
    write(*,*) 'bmag       =' , bmag
    write(*,*) 'theta      =' , theta
    write(*,*) 'BdotT      =' , BdotT
    write(*,*) 'BdotR      =' , BdotR


    call calculate_bplane_data(mu,rv,BdotR2,BdotT2,bmag2,theta2,istat)
    write(*,*) ''
    write(*,*) 'Alternate version:'
    write(*,*) ''
    write(*,*) 'bmag       =' , bmag2
    write(*,*) 'theta      =' , theta2
    write(*,*) 'BdotT      =' , BdotT2
    write(*,*) 'BdotR      =' , BdotR2

    write(*,*) ''
    write(*,*) 'Difference:'
    write(*,*) ''
    write(*,*) 'bmag       =' , bmag2 - bmag
    write(*,*) 'theta      =' , theta2 - theta
    write(*,*) 'BdotT      =' , BdotT2 - BdotT
    write(*,*) 'BdotR      =' , BdotR2 - BdotR

    end subroutine bplane_test
!*****************************************************************************************

!*****************************************************************************************
    end module bplane_module
!*****************************************************************************************
