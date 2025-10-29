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
    public :: crtbp_derivs_with_stm
    public :: normalize_variables, unnormalize_variables
    public :: compute_libration_points,compute_libration_points_v2
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
!  Compute the coordinates of the libration points (L1,L2,L3,L4,L5).
!  L1-L3 are computed using Newton's method. L4-L5 are known analytically.
!
!@note The coordinate are w.r.t. the barycenter of the system.

    subroutine compute_libration_points(mu,r1,r2,r3,r4,r5)

    use math_module, only: cube_root

    implicit none

    real(wp),intent(in)                        :: mu  !! CRTBP parameter
    real(wp),intent(out),optional              :: r1  !! L1 x coordinate
    real(wp),intent(out),optional              :: r2  !! L2 x coordinate
    real(wp),intent(out),optional              :: r3  !! L3 x coordinate
    real(wp),dimension(2),intent(out),optional :: r4  !! L4 [x,y] coordinates
    real(wp),dimension(2),intent(out),optional :: r5  !! L5 [x,y] coordinates

    integer  :: i  !! counter
    real(wp) :: f  !! quintic function value (to be driven to zero)
    real(wp) :: fp !! derivative of quintic function
    real(wp) :: x  !! indep. variable in the quintic functions

    integer,parameter  :: maxiter = 100    !! maximum number of
                                           !! iterations for newton's method
    real(wp),parameter :: tol = 1.0e-12_wp !! convergence tolerance for
                                           !! newton's method

    !L1, L2, and L3 are solved using iterative Newton method:

    if (present(r1)) then
        x = cube_root(mu / (3.0_wp - 3.0_wp * mu)) !initial guess
        do i = 1,maxiter
           f = x**5 - &
               (3.0_wp - mu) * x**4 + &
               (3.0_wp - 2.0_wp * mu) * x**3 - &
               mu * x**2 + &
               2.0_wp * mu * x - &
               mu
           if (abs(f) <= tol) exit
           fp = 5.0_wp * x**4 - &
                4.0_wp * x**3 * (3.0_wp - mu) + &
                3.0_wp * x**2 * (3.0_wp - 2.0_wp * mu) - &
                2.0_wp * x * mu + &
                2.0_wp * mu
           x = x - f / fp
        end do
        r1 = 1.0_wp - x  ! wrt primary body
        r1 = r1 - mu     ! wrt barycenter
    end if

    if (present(r2)) then
        x = cube_root(mu / (3.0_wp - 3.0_wp * mu)) !initial guess
        do i = 1,maxiter
           f = x**5 + &
               (3.0_wp - mu) * x**4 + &
               (3.0_wp - 2.0_wp * mu) * x**3 - &
               mu * x**2 - &
               2.0_wp * mu * x - &
               mu
           if (abs(f) <= tol) exit
           fp = 5.0_wp * x**4 + &
                4.0_wp * x**3 * (3.0_wp - mu) + &
                3.0_wp * x**2 * (3.0_wp - 2.0_wp * mu) - &
                2.0_wp * x * mu - &
                2.0_wp * mu
           x = x - f / fp
        end do
        r2 = 1.0_wp + x  ! wrt primary body
        r2 = r2 - mu     ! wrt barycenter
    end if

    if (present(r3)) then
        x = -(7.0_wp / 12.0_wp) * mu !initial guess
        do i = 1,maxiter
           f = x**5 + &
               (7.0_wp + mu) * x**4 + &
               (19.0_wp + 6.0_wp * mu) * x**3 + &
               (24.0_wp + 13.0_wp * mu) * x**2 + &
               2.0_wp * (6.0_wp + 7.0_wp * mu) * x + &
               7.0_wp * mu
           if (abs(f) <= tol) exit
           fp = 5.0_wp * x**4 + &
                4.0_wp * x**3 * (7.0_wp + mu) + &
                3.0_wp * x**2 * (19.0_wp + 6.0_wp * mu) + &
                2.0_wp * x * (24.0_wp + 13.0_wp * mu) + &
                2.0_wp * (6.0_wp + 7.0_wp * mu)
           x = x - f / fp
        end do
        r3 = -(x + 1.0_wp)  ! wrt primary body
        r3 = r3 - mu        ! wrt barycenter
    end if

    ! L4 and L5 are analytic:
    if (present(r4)) r4 = [0.5_wp - mu,  sqrt(3.0_wp)/2.0_wp]
    if (present(r5)) r5 = [0.5_wp - mu, -sqrt(3.0_wp)/2.0_wp]

    end subroutine compute_libration_points
!*******************************************************************************

!*******************************************************************************
!>
!  Compute the coordinates of the libration points (L1,L2,L3,L4,L5).
!
!  This is just an alternate version of [[compute_libration_points]].
!
!### Reference
!  * J.S. Parker, R.L. Anderson, "Low-Energy Lunar Trajectory Design", 2014.
!   (Appendix A.5)
!
!@note The coordinate are w.r.t. the barycenter of the system.

    subroutine compute_libration_points_v2(mu,r1,r2,r3,r4,r5)

    use math_module, only: cube_root

    implicit none

    real(wp),intent(in)                        :: mu  !! CRTBP parameter
    real(wp),intent(out),optional              :: r1  !! L1 x coordinate
    real(wp),intent(out),optional              :: r2  !! L2 x coordinate
    real(wp),intent(out),optional              :: r3  !! L3 x coordinate
    real(wp),dimension(2),intent(out),optional :: r4  !! L4 [x,y] coordinates
    real(wp),dimension(2),intent(out),optional :: r5  !! L5 [x,y] coordinates

    integer,parameter  :: maxiter = 100    !! maximum number of iterations
    real(wp),parameter :: tol = 1.0e-12_wp !! convergence tolerance

    real(wp) :: gamma, gamma0
    integer :: i !! counter

    if (present(r1)) then
        gamma0 = cube_root(mu * (one-mu) / three)
        gamma = gamma0 + one
        do i = 1,maxiter
            if (abs(gamma-gamma0)<=tol) exit
            gamma0 = gamma
            gamma = cube_root((mu*(gamma0-one)**2)/(three-two*mu-gamma0*(three-mu-gamma0)))
        end do
        r1 = one - mu - gamma
    end if

    if (present(r2)) then
        gamma0 = cube_root(mu * (one-mu) / three)
        gamma = gamma0 + one
        do i = 1,maxiter
            if (abs(gamma-gamma0)<=tol) exit
            gamma0 = gamma
            gamma = cube_root((mu*(gamma0+one)**2)/(three-two*mu+gamma0*(three-mu+gamma0)))
        end do
        r2 = one - mu + gamma
    end if

    if (present(r3)) then
        gamma0 = cube_root(mu * (one-mu) / three)
        gamma = gamma0 + one
        do i = 1,maxiter
            if (abs(gamma-gamma0)<=tol) exit
            gamma0 = gamma
            gamma = cube_root((one-mu)*(gamma0+one)**2/(one+two*mu+gamma0*(two+mu+gamma0)))
        end do
        r3 = - mu - gamma
    end if

    call compute_libration_points(mu,r4=r4,r5=r5)

    end subroutine compute_libration_points_v2
!*******************************************************************************

!*******************************************************************************
!>
!  CRTBP derivatives: state only.

    subroutine crtbp_derivs(mu,x,dx)

    implicit none

    real(wp),intent(in)               :: mu   !! CRTBP parameter (See [[compute_crtpb_parameter]])
    real(wp),dimension(6),intent(in)  :: x    !! normalized state \([\mathbf{r},\mathbf{v}]\)
    real(wp),dimension(6),intent(out) :: dx   !! normalized state derivative \([\dot{\mathbf{r}},\dot{\mathbf{v}}]\)

    !local variables:
    real(wp),dimension(3) :: r1,r2,rb1,rb2,r,v,g
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
!  CRTBP derivatives: state + state transition matrix.

    subroutine crtbp_derivs_with_stm(mu,x,dx)

    implicit none

    real(wp),intent(in)                :: mu   !! CRTBP parameter (See [[compute_crtpb_parameter]])
    real(wp),dimension(42),intent(in)  :: x    !! normalized state and STM \([\mathbf{r},\mathbf{v},\mathbf{\Phi}]\)
    real(wp),dimension(42),intent(out) :: dx   !! normalized state and STM derivative \([\dot{\mathbf{r}},\dot{\mathbf{v}},\dot{\mathbf{\Phi}}]\)

    !local variables:
    real,dimension(3) :: rb1,rb2,r,v,g
    real(wp),dimension(6,6) :: A, phi, phi_dot
    real(wp) :: r1,r2,r13,r23,r15,r25,omm,tmu,tomm,c1,c2
    real(wp) :: Uxx,Uxy,Uxz,Uyx,Uyy,Uyz,Uzx,Uzy,Uzz

    !extract variables from x vector:
    r = x(1:3) ! position
    v = x(4:6) ! velocity

    !other parameters:
    omm  = one - mu
    rb1  = [-mu,zero,zero] ! location of body 1
    rb2  = [omm,zero,zero] ! location of body 2
    r1   = norm2(r - rb1)  ! body1 -> sc distance
    r2   = norm2(r - rb2)  ! body2 -> sc distance
    r13  = r1**3
    r23  = r2**3
    r15  = r1**5
    r25  = r2**5
    c1   = omm/r13
    c2   = mu/r23
    tmu  = three*mu
    tomm = three*omm

    !normalized gravity from both bodies:
    g(1) = -c1*(r(1) + mu) - c2*(r(1)-one+mu)
    g(2) = -c1*r(2)        - c2*r(2)
    g(3) = -c1*r(3)        - c2*r(3)

    !STM terms:
    Uxx = one - c1 - c2 + tomm*(r(1)+mu)**2/r15 + tmu*(r(1)-one+mu)**2/r25
    Uyy = one - c1 - c2 + tomm*r(2)**2/r15 + tmu*r(2)**2/r25
    Uzz =       c1 - c2 + tomm*r(3)**2/r15 + tmu*r(3)**2/r25
    Uxy = tomm*(r(1)+mu)*r(2)/r15 + tmu*(r(1)-one+mu)*r(2)/r25
    Uxz = tomm*(r(1)+mu)*r(3)/r15 + tmu*(r(1)-one+mu)*r(3)/r25
    Uyz = tomm* r(2)*r(3)/r15 + tmu*r(2)*r(3)/r25
    Uyx = Uxy
    Uzx = Uxz
    Uzy = Uyz

    !columns of A matrix:
    A(:,1) = [zero,zero,zero,Uxx,Uyx,Uzx]
    A(:,2) = [zero,zero,zero,Uxy,Uyy,Uzy]
    A(:,3) = [zero,zero,zero,Uxz,Uyz,Uzz]
    A(:,4) = [one,zero,zero,zero,-two,zero]
    A(:,5) = [zero,one,zero,two,zero,zero]
    A(:,6) = [zero,zero,one,zero,zero,zero]

    !unpack phi into matrix:
    phi = reshape(x(7:42), shape=[6,6])

    !derivative of phi matrix:
    phi_dot = matmul(A,phi)

    !derivative of x vector:
    dx(1:3)  = v                           ! r_dot
    dx(4)    =  two*v(2) + r(1) + g(1)     ! v_dot
    dx(5)    = -two*v(1) + r(2) + g(2)     !
    dx(6)    =                    g(3)     !
    dx(7:42) = pack (phi_dot, mask=.true.) ! phi_dot

    end subroutine crtbp_derivs_with_stm
!*******************************************************************************

!*******************************************************************************
!>
!  Unit tests for CRTBP routines.

    subroutine crtbp_test()

    use celestial_body_module

    implicit none

    real(wp),parameter :: mu_earth = body_earth%mu !! \( \mu_{Earth} ~ (\mathrm{km}^3/\mathrm{s}^2) \)
    real(wp),parameter :: mu_moon  = body_moon%mu  !! \( \mu_{Moon}  ~ (\mathrm{km}^3/\mathrm{s}^2) \)
    real(wp),parameter :: mu_sun   = body_sun%mu   !! \( \mu_{Sun}   ~ (\mathrm{km}^3/\mathrm{s}^2) \)

    !< sample state (normalized)
    !< see: [Celestial Mechanics Notes Set 4: The Circular Restricted
    !< Three Body Problem](http://cosweb1.fau.edu/~jmirelesjames/hw4Notes.pdf), p.40.
    real(wp),dimension(6),parameter :: x = [  0.30910452642073_wp, &
                                              0.07738174525518_wp, &
                                              0.0_wp,              &
                                             -0.72560796964234_wp, &
                                              1.55464233412773_wp, &
                                              0.0_wp               ]

    integer                 :: i       !! counter
    real(wp)                :: mu      !! CRTPB parameter
    real(wp)                :: mu1     !! primary body mu
    real(wp)                :: mu2     !! secondary body mu
    real(wp)                :: c       !! Jacobi constant
    real(wp),dimension(6)   :: xd      !! derivative vector: state
    real(wp),dimension(42)  :: x_phi   !! initial state + phi (identity)
    real(wp),dimension(42)  :: x_phi_d !! derivative vector: state + phi
    real(wp),dimension(6,6) :: eye     !! 6x6 identity matrix
    real(wp)                :: r1      !! L1 x coordinate (normalized)
    real(wp)                :: r2      !! L2 x coordinate (normalized)
    real(wp)                :: r3      !! L3 x coordinate (normalized)
    real(wp),dimension(2)   :: r4      !! L4 x coordinate (normalized)
    real(wp),dimension(2)   :: r5      !! L5 x coordinate (normalized)

    !create an identity matrix for stm initial condition:
    eye = zero
    do i = 1, 6
        eye(i,i) = one
    end do
    x_phi = [x, pack(eye,mask=.true.)]

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' crtbp_test'
    write(*,*) '---------------'
    write(*,*) ''

    do i=1,3

        select case (i)
        case(1)
            mu1 = mu_earth
            mu2 = mu_moon
        case(2)
            mu1 = mu_earth
            mu2 = mu_earth
        case(3)
            mu1 = mu_earth + mu_moon/four
            mu2 = mu_earth
        end select

        write(*,*) ''
        mu = compute_crtpb_parameter(mu1,mu2)
        c = compute_jacobi_constant(mu,x)
        call crtbp_derivs(mu,x,xd)
        call crtbp_derivs_with_stm(mu,x_phi,x_phi_d)
        call compute_libration_points(mu,r1,r2,r3,r4,r5)

        write(*,'(A,1X,*(F30.16,1X))')  'mu:         ', mu
        write(*,'(A,1X,*(F30.16,1X))' ) 'L1 x:       ', r1
        write(*,'(A,1X,*(F30.16,1X))' ) 'L2 x:       ', r2
        write(*,'(A,1X,*(F30.16,1X))' ) 'L3 x:       ', r3
        write(*,'(A,1X,*(F30.16,1X))' ) 'L4 x:       ', r4
        write(*,'(A,1X,*(F30.16,1X))' ) 'L5 x:       ', r5
        write(*,'(A,1X,*(F30.16,1X))' ) 'x:          ', x
        write(*,'(A,1X,*(F30.16,1X))' ) 'c:          ', c
        write(*,'(A,1X,*(F30.16,1X))' ) 'xd:         ', xd
        write(*,'(A,1X,*(F30.16,1X))' ) 'x+phi:      ', x_phi
        write(*,'(A,1X,*(F30.16,1X))' ) 'xd+phi_dot: ', x_phi_d
        write(*,*) ''

        call compute_libration_points_v2(mu,r1,r2,r3,r4,r5)
        write(*,*) ''
        write(*,*) 'alternate formulation:'
        write(*,'(A,1X,*(F30.16,1X))' ) 'L1 x:       ', r1
        write(*,'(A,1X,*(F30.16,1X))' ) 'L2 x:       ', r2
        write(*,'(A,1X,*(F30.16,1X))' ) 'L3 x:       ', r3
        write(*,'(A,1X,*(F30.16,1X))' ) 'L4 x:       ', r4
        write(*,'(A,1X,*(F30.16,1X))' ) 'L5 x:       ', r5
        write(*,*) ''

    end do

    end subroutine crtbp_test
!*******************************************************************************

!*******************************************************************************
    end module crtbp_module
!*******************************************************************************
