!*****************************************************************************************
!> author: Jacob Williams
!  date: 11/27/2015
!
!  Modified equinoctial elements routines.
!
!  The modified equinoctial elements are applicable to all orbits
!  and have non-singular equations of motion (except for a singularity
!  at \( i = \pi \) ). They are defined as:
!
!  $$ \begin{array}{rcl}
!      p & = & a (1 - e^2) \\
!      f & = & e \cos (\omega + \Omega) \\
!      g & = & e \sin (\omega + \Omega)  \\
!      h & = & \tan (i/2) \cos \Omega \\
!      k & = & \tan (i/2) \sin \Omega \\
!      L & = & \Omega + \omega + \nu \\
!  \end{array} $$
!
!  Where \(L\) is the true longitude, \(p\) is the semi-latus rectum,
!  and \(\nu\) is the true anomaly.
!
!### References
!  * Broucke, R. A. & Cefola, P. J., "On the Equinoctial Orbit Elements"
!    Celestial Mechanics, Volume 5, Issue 3, p 303-310. (1972)
!  * M. J. H. Walker, B. Ireland, Joyce Owens,
!    "A Set of Modified Equinoctial Orbit Elements"
!    Celestial Mechanics, Volume 36, Issue 4, p 409-419. (1985)
!  * Walker, M. J. H, "Erratum - a Set of Modified Equinoctial Orbit Elements"
!    Celestial Mechanics, Volume 38, Issue 4, p 391-392. (1986)

    module modified_equinoctial_module

    use kind_module,      only: wp
    use numbers_module
    use vector_module,    only: unit,cross

    implicit none

    private

    public :: cartesian_to_equinoctial
    public :: equinoctial_to_cartesian
    public :: modified_equinoctial_derivs

    public :: modified_equinoctial_test ! for testing

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Convert Cartesian coordinates to modified equinoctial elements (posigrade formulation).

    subroutine cartesian_to_equinoctial(mu,rv,evec)

    implicit none

    real(wp),intent(in)               :: mu    !! central body gravitational parameter (\(km^3/s^2\))
    real(wp),dimension(6),intent(in)  :: rv    !! Cartesian state vector
    real(wp),dimension(6),intent(out) :: evec  !! Modified equinoctial element vector

    real(wp),dimension(3) :: r,v,hvec,hhat,ecc,fhat,ghat,rhat,vhat
    real(wp) :: hmag,rmag,p,f,g,h,k,L,kk,hh,s2,tkh,rdv

    r = rv(1:3)
    v = rv(4:6)

    rdv      = dot_product(r,v)
    rhat     = unit(r)
    rmag     = norm2(r)
    hvec     = cross(r,v)
    hmag     = norm2(hvec)
    hhat     = unit(hvec)
    vhat     = (rmag*v - rdv*rhat)/hmag
    p        = hmag*hmag / mu
    k        = hhat(1)/(one + hhat(3))
    h        = -hhat(2)/(one + hhat(3))
    kk       = k*k
    hh       = h*h
    s2       = one+hh+kk
    tkh      = two*k*h
    ecc      = cross(v,hvec)/mu - rhat
    fhat(1)  = one-kk+hh
    fhat(2)  = tkh
    fhat(3)  = -two*k
    ghat(1)  = tkh
    ghat(2)  = one+kk-hh
    ghat(3)  = two*h
    fhat     = fhat/s2
    ghat     = ghat/s2
    f        = dot_product(ecc,fhat)
    g        = dot_product(ecc,ghat)
    L        = atan2(rhat(2)-vhat(1),rhat(1)+vhat(2))

    evec = [p,f,g,h,k,L]

    end subroutine cartesian_to_equinoctial
!*****************************************************************************************

!*****************************************************************************************
!>
!  Convert modified equinoctial elements (posigrade formulation) to Cartesian coordinates.

    subroutine equinoctial_to_cartesian(mu,evec,rv)

    implicit none

    real(wp),intent(in)               :: mu   !! central body gravitational parameter (\(km^3/s^2\))
    real(wp),dimension(6),intent(in)  :: evec !! Modified equinoctial element vector
    real(wp),dimension(6),intent(out) :: rv   !! Cartesian state vector

    real(wp) :: p,f,g,h,k,L,s2,r,w,cL,sL,smp,hh,kk,tkh
    real(wp) :: x,y,xdot,ydot
    real(wp),dimension(3) :: fhat,ghat

    p       = evec(1)
    f       = evec(2)
    g       = evec(3)
    h       = evec(4)
    k       = evec(5)
    L       = evec(6)

    kk      = k*k
    hh      = h*h
    tkh     = two*k*h
    s2      = one + hh + kk
    cL      = cos(L)
    sL      = sin(L)
    w       = one + f*cL + g*sL
    r       = p/w
    smp     = sqrt(mu/p)
    fhat(1) = one-kk+hh
    fhat(2) = tkh
    fhat(3) = -two*k
    ghat(1) = tkh
    ghat(2) = one+kk-hh
    ghat(3) = two*h
    fhat    = fhat/s2
    ghat    = ghat/s2
    x       = r*cL
    y       = r*sL
    xdot    = -smp * (g + sL)
    ydot    =  smp * (f + cL)

    rv(1:3) = x*fhat + y*ghat
    rv(4:6) = xdot*fhat + ydot*ghat

    end subroutine equinoctial_to_cartesian
!*****************************************************************************************

!*****************************************************************************************
!>
!  Modified equinoctial elements (posigrade formulation) equations of motion.

    subroutine modified_equinoctial_derivs(mu,evec,scn,evecd)

    implicit none

    real(wp),intent(in)               :: mu    !! central body gravitational parameter (\(km^3/s^2\))
    real(wp),dimension(6),intent(in)  :: evec  !! modified equinoctial element vector
    real(wp),dimension(3),intent(in)  :: scn   !! Perturbation (in the RSW frame)
    real(wp),dimension(6),intent(out) :: evecd !! derivative of `evec`

    real(wp) :: p,f,g,h,k,L,c,s,n,sqrtpm,sl,cl,s2no2w,s2,w

    p      = evec(1)
    f      = evec(2)
    g      = evec(3)
    h      = evec(4)
    k      = evec(5)
    L      = evec(6)
    s      = scn(1)
    c      = scn(2)
    n      = scn(3)
    sqrtpm = sqrt(p/mu)
    sl     = sin(L)
    cl     = cos(L)
    s2     = one + h*h + k*k
    w      = one + f*cl + g*sl
    s2no2w = s2*n / (two*w)

    evecd(1) = (two*p*c/w) * sqrtpm
    evecd(2) = sqrtpm * ( s*sl + ((w+one)*cl+f)*c/w - g*n*(h*sl-k*cl)/w )
    evecd(3) = sqrtpm * (-s*cl + ((w+one)*sl+g)*c/w + f*n*(h*sl-k*cl)/w )
    evecd(4) = sqrtpm * s2no2w * cl
    evecd(5) = sqrtpm * s2no2w * sl
    evecd(6) = sqrt(mu*p)*(w/p)**2 + sqrtpm * ((h*sl-k*cl)*n)/w

    end subroutine modified_equinoctial_derivs
!*****************************************************************************************

!*****************************************************************************************
!>
!  Unit tests for the modified_equinoctial_module.

    subroutine modified_equinoctial_test()

    implicit none

    real(wp),parameter :: mu = 398600.436233_wp  !! central body grav. parameter
    real(wp),dimension(6),parameter :: x0 = [-2301.67224489839_wp, &
                                             -5371.07610250925_wp, &
                                             -3421.14671530212_wp, &
                                              6.1338624555516_wp,  &
                                              .306265184163608_wp, &
                                             -4.59713439017524_wp  ] !! state vector
    real(wp),dimension(3),parameter :: scn = [zero,zero,zero] !! perturbation

    real(wp),dimension(6) :: e,ed,x2

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' modified_equinoctial_test'
    write(*,*) '---------------'
    write(*,*) ''

    call cartesian_to_equinoctial(mu,x0,e)
    call equinoctial_to_cartesian(mu,e,x2)
    call modified_equinoctial_derivs(mu,e,scn,ed)

    write(*,'(A,1X,*(F20.6,1X))' ) 'x:       ', x0
    write(*,'(A,1X,*(F20.6,1X))' ) 'x2:      ', x2
    write(*,'(A,1X,*(F20.6,1X))' ) 'x error: ', x2-x0
    write(*,'(A,1X,*(F20.6,1X))' ) 'e:       ', e
    write(*,'(A,1X,*(F20.6,1X))' ) 'ed:      ', ed

    end subroutine modified_equinoctial_test
!*****************************************************************************************

!*****************************************************************************************
    end module modified_equinoctial_module
!*****************************************************************************************
