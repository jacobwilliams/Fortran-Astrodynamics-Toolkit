!*****************************************************************************************
!> author: Jacob Williams
!
!  This module contains various routines related to the
!  Circular Restricted Three-Body Problem (CRTBP).

    module crtbp_module

    use kind_module,    only: wp
    use numbers_module

    implicit none

    private

    public :: compute_crtpb_parameter
    public :: compute_jacobi_constant
    public :: crtbp_derivs
    public :: normalize_variables, unnormalize_variables
    public :: crtbp_test

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  Compute \( \mu \), the normalized CRTBP parameter.
!  It is equal to \( M_2 / (M_1 + M_2) \).
!
!@note The inputs can either be mu's or masses,
!      as long are they are both the same units.

    pure function compute_crtpb_parameter(mu1,mu2) result(mu)

    implicit none

    real(wp),intent(in) :: mu1  !! grav param for body 1 \([km^3/s^2]\)
    real(wp),intent(in) :: mu2  !! grav param for body 2 \([km^3/s^2]\)
    real(wp)            :: mu   !! CRTBP parameter \(\mu\)

    mu = mu2 / (mu1 + mu2)

    end function compute_crtpb_parameter
!*******************************************************************************

!*******************************************************************************
!>
!  Convert state in km, km/s units to normalized CRTBP state.

    subroutine normalize_variables(mu1,mu2,d12,x,m,t,x_crtbp,m_crtbp,t_crtbp)

    implicit none

    real(wp),intent(in)                        :: mu1      !! grav. param. of body 1 \( [\mathrm{km}^3/\mathrm{s}^2] \)
    real(wp),intent(in)                        :: mu2      !! grav. param. of body 2 \( [\mathrm{km}^3/\mathrm{s}^2] \)
    real(wp),intent(in)                        :: d12      !! distance between body 1 and body 2 [km]
    real(wp),dimension(6),intent(in),optional  :: x        !! unnormalized state w.r.t. barycenter [km,km/s]
    real(wp),intent(in),optional               :: m        !! unnormalized mass [kg]
    real(wp),intent(in),optional               :: t        !! unnormalized time [sec]
    real(wp),dimension(6),intent(out),optional :: x_crtbp  !! CRTBP normalized state
    real(wp),intent(out),optional              :: m_crtbp  !! CRTBP normalized mass
    real(wp),intent(out),optional              :: t_crtbp  !! CRTBP normalized time

    real(wp) :: tp,n,tu,du,mu

    n = sqrt((mu1+mu2)/d12**3)  !mean motion (rad/sec)
    tu = one/n  !time unit
    du = d12    !distance unit

    if (present(x) .and. present(x_crtbp)) then
        x_crtbp(1:3) = x(1:3) / du       !scale distance
        x_crtbp(4:6) = x(4:6) / (du/tu)  !scale velocity
    end if

    if (present(m) .and. present(m_crtbp)) then
        mu = (mu1 + mu2)/universal_grav_constant    !mass unit
        m_crtbp = m / mu  !scale mass
    end if

    if (present(t) .and. present(t_crtbp)) t_crtbp = t / tu  !scale time

    end subroutine normalize_variables
!*******************************************************************************

!*******************************************************************************
!>
!  Convert normalized CRTBP state to km, km/s units.
!
!# Notes:
!
! See also: http://www.spaceatdia.org/uploads/mariano/ss1/2012SSLecture3.pdf
!
!   m1 = one - mu   ! mass of body 1
!   m2 = mu         ! mass of body 2
!
!   x1 = -mu        ! location of body 1
!   x2 = one - mu   ! location of body 2
!
!   normalized mass     : 1 MU -> (m1 + m2) kg
!   normalized position : 1 DU -> d12 km
!   normalized time     : 1 TU -> 1/n sec
!                         (2pi TU -> 1 rev -> 2pi/n)

    subroutine unnormalize_variables(mu1,mu2,d12,x_crtbp,m_crtbp,t_crtbp,x,m,t)

    implicit none

    real(wp),intent(in)                        :: mu1      !! grav. param. of body 1 \( [\mathrm{km}^3/\mathrm{s}^2] \)
    real(wp),intent(in)                        :: mu2      !! grav. param. of body 2 \( [\mathrm{km}^3/\mathrm{s}^2] \)
    real(wp),intent(in)                        :: d12      !! distance between body 1 and body 2 [km]
    real(wp),dimension(6),intent(in),optional  :: x_crtbp  !! CRTBP normalized state
    real(wp),intent(in),optional               :: m_crtbp  !! CRTBP normalized mass
    real(wp),intent(in),optional               :: t_crtbp  !! CRTBP normalized time
    real(wp),dimension(6),intent(out),optional :: x        !! unnormalized state w.r.t. barycenter [km,km/s]
    real(wp),intent(out),optional              :: m        !! unnormalized mass [kg]
    real(wp),intent(out),optional              :: t        !! unnormalized time [sec]

    real(wp) :: tp,n,tu,du,mu

    n = sqrt((mu1+mu2)/d12**3)  !mean motion (rad/sec)
    tu = one/n  !time unit
    du = d12    !distance unit

    if (present(x) .and. present(x_crtbp)) then
        x(1:3) = x_crtbp(1:3) * du       !unscale distance
        x(4:6) = x_crtbp(4:6) * (du/tu)  !unscale velocity
    end if

    if (present(m) .and. present(m_crtbp)) then
        mu = (mu1 + mu2)/universal_grav_constant    !mass unit
        m = m_crtbp * mu  !unscale mass
    end if

    if (present(t) .and. present(t_crtbp)) t = t_crtbp * tu  !unscale time

    end subroutine unnormalize_variables
!*******************************************************************************

!*******************************************************************************
!>
!  Compute the CRTBP Jacobi constant, given the state.

    pure function compute_jacobi_constant(mu,x) result(c)

    implicit none

    real(wp),intent(in)              :: mu   !! CRTBP parameter (See [[compute_crtpb_parameter]])
    real(wp),dimension(6),intent(in) :: x    !! normalized state vector
    real(wp)                         :: c    !! Jacobi constant

    !local variables:
    real,dimension(3) :: r,v,rb1,rb2
    real(wp) :: omm,r1,r2

    !extract variables from x vector:
    r = x(1:3)    ! position
    v = x(4:6)    ! velocity

    !other parameters:
    omm = one - mu
    rb1 = [-mu,zero,zero]   ! location of body 1
    rb2 = [omm,zero,zero]   ! location of body 2
    r1  = norm2(r - rb1)    ! body1 -> sc distance
    r2  = norm2(r - rb2)    ! body2 -> sc distance

    !compute Jacobi integral:
    ! [ See: http://cosweb1.fau.edu/~jmirelesjames/hw4Notes.pdf ]
    if (r1==zero .or. r2==zero) then
        c = huge(one)   ! a large value
    else
        c = r(1)**2 + r(2)**2 + &
            two*omm/r1 + two*mu/r2 - &
            (v(1)**2 + v(2)**2 + v(3)**2)
    end if

    end function compute_jacobi_constant
!*******************************************************************************

!*******************************************************************************
!>
!  CRTBP derivatives.

    subroutine crtbp_derivs(mu,x,dx)

    implicit none

    real(wp),intent(in)               :: mu   !! CRTBP parameter (See [[compute_crtpb_parameter]])
    real(wp),dimension(6),intent(in)  :: x    !! normalized state \([\mathbf{r},\mathbf{v}]\)
    real(wp),dimension(6),intent(out) :: dx   !! normalized state derivative \([\dot{\mathbf{r}},\dot{\mathbf{v}}]\)

    !local variables:
    real,dimension(3) :: r1,r2,rb1,rb2,r,v,g
    real(wp) :: r13,r23,omm,c1,c2

    !extract variables from x vector:
    r = x(1:3) ! position
    v = x(4:6) ! velocity

    !other parameters:
    omm = one - mu
    rb1 = [-mu,zero,zero] ! location of body 1
    rb2 = [omm,zero,zero] ! location of body 2
    r1  = r - rb1         ! body1 -> sc vector
    r2  = r - rb2         ! body2 -> sc vector
    r13 = norm2(r1)**3
    r23 = norm2(r2)**3
    c1  = omm/r13
    c2  = mu/r23

    !normalized gravity from both bodies:
    g(1) = -c1*(r(1) + mu) - c2*(r(1)-one+mu)
    g(2) = -c1*r(2)        - c2*r(2)
    g(3) = -c1*r(3)        - c2*r(3)

    ! derivative of x:
    dx(1:3) = v                       ! rdot
    dx(4)   =  two*v(2) + r(1) + g(1) ! vdot
    dx(5)   = -two*v(1) + r(2) + g(2) !
    dx(6)   =                    g(3) !

    end subroutine crtbp_derivs
!*******************************************************************************

!*******************************************************************************
!>
!  Unit tests for CRTBP routines.

    subroutine crtbp_test()

    implicit none

    real(wp),parameter :: mu_earth = 398600.435608_wp    !! \( \mu_{Earth} ~ (\mathrm{km}^3/\mathrm{s}^2) \)
    real(wp),parameter :: mu_moon  = 4902.799108_wp      !! \( \mu_{Moon}  ~ (\mathrm{km}^3/\mathrm{s}^2) \)
    real(wp),parameter :: mu_sun   = 132712440017.987_wp !! \( \mu_{Sun}   ~ (\mathrm{km}^3/\mathrm{s}^2) \)

    !< sample state (normalized)
    !< see: [Celestial Mechanics Notes Set 4: The Circular Restricted Three Body Problem](http://cosweb1.fau.edu/~jmirelesjames/hw4Notes.pdf), p.40.
    real(wp),dimension(6),parameter :: x = [  0.30910452642073_wp, &
                                              0.07738174525518_wp, &
                                              0.0_wp,              &
                                             -0.72560796964234_wp, &
                                              1.55464233412773_wp, &
                                              0.0_wp               ]

    real(wp)              :: mu  !! CRTPB parameter
    real(wp)              :: c   !! Jacobi constant
    real(wp),dimension(6) :: xd  !! derivative vector

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' crtbp_test'
    write(*,*) '---------------'
    write(*,*) ''

    mu = compute_crtpb_parameter(mu_earth,mu_moon)
    c = compute_jacobi_constant(mu,x)
    call crtbp_derivs(mu,x,xd)

    write(*,'(A,1X,*(F30.16,1X))') 'Earth-Moon mu:   ', mu
    write(*,'(A,1X,*(F12.6,1X))' ) 'Sample point x:  ', x
    write(*,'(A,1X,*(F12.6,1X))' ) 'Sample point c:  ', c
    write(*,'(A,1X,*(F12.6,1X))' ) 'Sample point xd: ', xd
    write(*,*) ''

    end subroutine crtbp_test
!*******************************************************************************

!*******************************************************************************
    end module crtbp_module
!*******************************************************************************
