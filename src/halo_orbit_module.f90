!*******************************************************************************
!>
!  Halo orbit routines.
!
!### References
!  * D.L. Richardson, "Analytic Construction of Periodic
!    Orbits About the Collinear Points", Celestial Mechanics 22 (1980)
!
!@todo Add differentially-corrected option using the STM derivatives.

    module halo_orbit_module

    use iso_fortran_env, only: wp => real64, error_unit
    use numbers_module
    use crtbp_module

    implicit none

    private

    public :: halo_to_rv
    public :: halo_to_rv_diffcorr
    public :: compute_halo_monodromy_matrix
    public :: compute_monodromy_matrix_eigenvalues

    public :: halo_orbit_test ! test routine

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  Compute the state vector for a halo orbit.
!  This uses the approximation, which is retargeted in the
!  real CR3BP system to produce a periodic orbit.
!
!@todo use a variable-step size integrator

    subroutine halo_to_rv_diffcorr(libpoint,mu1,mu2,dist,A_z,n,t1,rv,info,period)

    use rk_module
    use minpack_module
    use math_module,     only: wrap_angle

    implicit none

    integer,intent(in)  :: libpoint         !! Libration point number: [1,2,3]
    real(wp),intent(in) :: mu1              !! grav param for primary body [km3/s2]
    real(wp),intent(in) :: mu2              !! grav param for secondary body [km3/s2]
    real(wp),intent(in) :: dist             !! distance between bodies [km]
    real(wp),intent(in) :: A_z              !! halo z amplitude [km]
    integer,intent(in)  :: n                !! halo family: 1, 3
    real(wp),intent(in) :: t1               !! tau1 [rad]
    real(wp),dimension(6),intent(out) :: rv !! cr3bp normalized state vector
                                            !! [wrt barycenter]
    integer,intent(out) :: info             !! status code (1=no errors)
    real(wp),intent(out),optional :: period !! period of halo (normalized time units)

    integer,parameter  :: n_state_vars    = 6         !! number of state variables in the equations of motion
    integer,parameter  :: n_opt_vars      = 2         !! number of variables in the targeting problem
    real(wp),parameter :: t0              = 0.0_wp    !! initial time (normalized) (epoch doesn't matter for cr3bp)
    real(wp),parameter :: tol             = 1.0e-8_wp !! tolerance for event finding
    real(wp),parameter :: xtol            = 1.0e-6_wp !! tolerance for [[hybrd]]
    integer,parameter  :: maxfev          = 1000      !! max number of function evaluations for [[hybrd]]
    integer,parameter  :: n_steps_per_rev = 100       !! number of integration steps per orbit rev

    type(rk8_10_class)               :: prop   !! integrator
    real(wp),dimension(n_opt_vars)   :: x_vy0  !! variables in the targeting problem (x0 and vy0)
    real(wp),dimension(n_opt_vars)   :: vx_vzf !! constraints in the targeting problem (vxf and vzf)
    real(wp),dimension(n_state_vars) :: x0     !! halo initial guess from richardson approximation
    real(wp),dimension(n_state_vars) :: xf     !! state after 1/2 rev (to get the period)
    real(wp) :: tf_actual     !! 1/2 period for retargeted orbit (normalized time)
    real(wp) :: actual_period !! actual halo orbit period for retargeted orbit (normalized time)
    real(wp) :: approx_period !! period approximation (normalized time)
    real(wp) :: dt_to_t1      !! time from `t1=0` to input `t1`
    real(wp) :: gf            !! function value after 1/2 rev (y-coordinate)
    real(wp) :: tau           !! `t1` wrapped from \(- \pi \) to \( \pi \)
    real(wp) :: dt            !! time step (normalized)
    real(wp) :: tmax          !! max final time for event finding integration
    real(wp) :: mu            !! CRTBP parameter

    ! compute the CRTBP mu parameter:
    mu = compute_crtpb_parameter(mu1,mu2)

    ! first we get the halo state approximation at tau1=0:
    call halo_to_rv(libpoint,mu1,mu2,dist,A_z,n,zero,x0,approx_period)

    ! for now, fixed number of integration steps per period:
    dt = approx_period / real(n_steps_per_rev,wp)
    tmax = two * approx_period ! should be enough to find the x-z crossing

    ! initialize the integrator:
    call prop%initialize(n_state_vars,func,g=xz_plane_crossing)

    ! now, solve for a halo:
    x_vy0 = [x0(1),x0(5)]  ! x0 and vy0
    call hybrd1(halo_fcn,2,x_vy0,vx_vzf,tol=xtol,info=info)
    if (info==1) then ! solution converged

        ! now have the solution at t1=0:
        x0(1) = x_vy0(1)
        x0(5) = x_vy0(2)

        ! this is the t1 we want:
        tau = wrap_angle(t1)

        ! if we need the period:
        if (present(period) .or. tau/=zero) then
            ! integrate to the first x-axis crossings (one half rev):
            ! [need to check output...]
            call prop%integrate_to_event(t0,x0,dt,tmax,tol,tf_actual,xf,gf)
            actual_period = two * tf_actual ! normalized period
        end if

        ! now we want to propagate to the input tau1
        if (tau==zero) then
            ! already have the solution
            rv = x0
        else
            ! now, integrate from t1=0 to input t1 to get rv:
            dt_to_t1 = actual_period * (tau / twopi)
            call prop%integrate(t0,x0,dt,dt_to_t1,rv)
        end if

        if (present(period)) period = actual_period

    else ! there was an error
        write(error_unit,'(A)') 'Error: the halo targeting problem did not converge.'
        rv = x0
    end if

    !call prop%destroy()

    contains
!*******************************************************************************

    !***************************************************************************
        subroutine halo_fcn(n,xvec,fvec,iflag)
        !! Halo function for [[hybrd1]]

        implicit none

        integer,intent(in)                :: n      !! `n=2` in this case
        real(wp),dimension(n),intent(in)  :: xvec   !! x_vy0
        real(wp),dimension(n),intent(out) :: fvec   !! [vxf,vzf]
        integer,intent(inout)             :: iflag  !! status flag (set negative
                                                    !! to terminate solver)

        real(wp) :: gf
        real(wp),dimension(6) :: x,x1,xf

        x    = x0      ! initial guess state (z is held fixed)
        x(1) = xvec(1) ! x0
        x(5) = xvec(2) ! vy0

        !integrate to the next x-z-plane crossing:
        call prop%integrate_to_event(t0,x,dt,tmax,tol,tf_actual,xf,gf)

        !want x and z-velocity at the x-z-plane crossing to be zero:
        fvec = [xf(4),xf(6)]

        end subroutine halo_fcn
    !***************************************************************************

    !***************************************************************************
        subroutine func(me,t,x,xdot)
        !! CRTBP derivative function
        implicit none
        class(rk_class),intent(inout)        :: me
        real(wp),intent(in)                  :: t
        real(wp),dimension(me%n),intent(in)  :: x
        real(wp),dimension(me%n),intent(out) :: xdot

        call crtbp_derivs(mu,x,xdot)

        end subroutine func
    !**************************************************************************

    !***************************************************************************
        subroutine xz_plane_crossing(me,t,x,g)
        !! x-z-plane crossing event function
        implicit none
        class(rk_class),intent(inout)        :: me
        real(wp),intent(in)                  :: t
        real(wp),dimension(me%n),intent(in)  :: x
        real(wp),intent(out)                 :: g

        g = x(2)  ! y = 0 at x-z-plane crossing

        end subroutine xz_plane_crossing
    !***************************************************************************

    end subroutine halo_to_rv_diffcorr
!*******************************************************************************

!*******************************************************************************
!>
!  Compute the state vector from the halo orbit approximation.
!  This will be an approximation of a halo orbit in the CR3BP system,
!  and will need to be corrected to produce a real halo orbit.

    subroutine halo_to_rv(libpoint,mu1,mu2,dist,A_z,n,t1,rv,period)

    implicit none

    integer,intent(in)  :: libpoint         !! Libration point number: [1,2,3]
    real(wp),intent(in) :: mu1              !! grav param for primary body [km3/s2]
    real(wp),intent(in) :: mu2              !! grav param for secondary body [km3/s2]
    real(wp),intent(in) :: dist             !! distance between bodies [km]
    real(wp),intent(in) :: A_z              !! halo z amplitude [km]
    integer,intent(in)  :: n                !! halo family: 1, 3
    real(wp),intent(in) :: t1               !! tau1 [rad]
    real(wp),dimension(6),intent(out) :: rv !! cr3bp normalized state vector
                                            !! [wrt barycenter]
    real(wp),intent(out),optional :: period !! normalized halo orbit period

    real(wp) :: mu          !! CRTBP parameter
    integer  :: delta_n     !! 2 - n
    real(wp) :: gamma_l     !! dimensionless quantity from reference
    real(wp) :: lambda      !! linearized frequency
    real(wp) :: Ax,Ay,Az,Ax2,Az2
    real(wp) :: x,y,z,vx,vy,vz
    real(wp) :: a1,a2,a21,a22,a23,a24,a31,a32,b21,&
                b22,b31,b32,c2,c3,c4,delta,w,&
                d1,d2,d21,d3,d31,d32,k,l1,l2,s1,s2,term
    real(wp),dimension(3) :: x_libpoint  !! x-coordinates of the libration
                                         !! point (wrt barycenter, normalized)
    logical :: ok

    ! error check:
    if (n/=1 .and. n/=3) then
        error stop 'invalid n input to halo_to_rv'
    end if

    ! compute all the intermediate parameters:
    mu = compute_crtpb_parameter(mu1,mu2)

    ! lib point x-coordinate: wrt to barycenter - normalized
    select case (libpoint)
    case(1)
        call compute_libration_points(mu,r1=x_libpoint(1))
        gamma_l = (1.0_wp - mu) - x_libpoint(1)
    case(2)
        call compute_libration_points(mu,r2=x_libpoint(2))
        gamma_l = x_libpoint(2) - ( 1.0_wp - mu )
    case(3)
        call compute_libration_points(mu,r3=x_libpoint(3))
        gamma_l = - (x_libpoint(3) + mu)
    case default
        error stop 'invalid libration point input to halo_to_rv'
    end select

    Az      = A_z / (dist*gamma_l) ! normalized z-Amplitude
    c2      = c_n(libpoint, 2, mu, gamma_l)
    c3      = c_n(libpoint, 3, mu, gamma_l)
    c4      = c_n(libpoint, 4, mu, gamma_l)
    lambda  = sqrt((-c2+two+sqrt(nine*c2**2-eight*c2))/two)  ! root of quartic eqn
    k       = two*lambda/(lambda**2+one-c2)
    delta   = lambda**2-c2
    d1      = 16.0_wp*lambda**4+four*lambda**2*(c2-two)-two*c2**2+c2+one
    d2      = 81.0_wp*lambda**4+nine*lambda**2*(c2-two)-two*c2**2+c2+one
    d3      = two*lambda*(lambda*(one+k**2)-two*k)
    a21     =  three*c3*(k**2-two)/four/(one+two*c2)
    a23     = -three*lambda*c3*(three*k**3*lambda-six*k*(k-lambda)+four)/four/k/d1
    b21     = -three*c3*lambda*(three*lambda*k-four)/two/d1
    s1      = ((three/two)*c3*(two*a21*(k**2-two)-a23*(k**2+two)-&
              two*k*b21)-(three/eight)*c4*(three*k**4-eight*k**2+eight))/d3
    a22     = three*c3/four/(one+two*c2)
    a24     = -three*c3*lambda*(two+three*lambda*k)/four/k/d1
    b22     = three*lambda*c3/d1
    d21     = -c3/two/lambda**2
    s2      = ((three/two)*c3*(two*a22*(k**2-two)+&
              a24*(k**2+two)+two*k*b22+five*d21) + &
              (three/eight)*c4*(12.0_wp-k**2))/d3
    a1      = -(three/two)*c3*(two*a21+a23+five*d21) - &
              (three/eight)*c4*(12.0_wp-k**2)
    a2      = (three/two)*c3*(a24-two*a22)+(nine/eight)*c4
    l1      = two*s1*lambda**2+a1
    l2      = two*s2*lambda**2+a2
    Az2     = Az**2

    ! check if this Az is feasible
    ok = (l1/=zero)
    if (ok) then
        term = (-delta-l2*Az2)/l1
        ok = (term>=zero)
    end if
    if (.not. ok) then
        rv = zero
        if (present(period)) period = zero
        write(error_unit,'(A)') 'Error: infeasible input.'
        return
    end if

    Ax      = sqrt(term) ! equation 18
    Ax2     = Ax**2
    Ay      = k*Ax
    w       = one+s1*Ax2+s2*Az2  ! frequency correction
    a31     = -nine*lambda*(c3*(k*a23-b21)+k*c4*(one+(one/four)*k**2))/d2 + &
              (nine*lambda**2+one-c2)*(three*c3*(two*a23-k*b21)+ &
              c4*(two+three*k**2))/two/d2
    a32     = -nine*lambda*(four*c3*(k*a24-b22)+k*c4)/four/d2 - &
              three*(nine*lambda**2+one-c2)*(c3*(k*b22+d21-two*a24)-c4)/two/d2
    b31     = (three*lambda*(three*c3*(k*b21-two*a23)-c4*(two+three*k**2)) + &
              (nine*lambda**2+1+two*c2)*(12.0_wp*c3*(k*a23-b21)+&
              three*k*c4*(four+k**2))/eight)/d2
    b32     = (three*lambda*(three*c3*(k*b22+d21-two*a24)-three*c4) + &
              (nine*lambda**2+one+two*c2)*(12.0_wp*c3*(k*a24-b22)+three*c4*k)/eight)/d2
    d31     = three*(four*c3*a24+c4)/64.0_wp/lambda**2
    d32     = three*(four*c3*(a23-d21)+c4*(four+k**2))/64.0_wp/lambda**2
    delta_n = 2 - n  ! equation 21

    if (present(period)) period = twopi/(lambda*w)

    ! Equations 20a, 20b, 20c (and their derivatives):
    x  = a21*Ax2+a22*Az2-Ax*cos(t1)+&
         (a23*Ax2-a24*Az2)*cos(two*t1)+&
         (a31*Ax**3-a32*Ax*Az2)*cos(three*t1)
    y  = k*Ax*sin(t1)+(b21*Ax2-b22*Az2)*sin(two*t1)+&
         (b31*Ax**3-b32*Ax*Az2)*sin(three*t1)
    z  = delta_n*(Az*cos(t1)+d21*Ax*Az*(cos(two*t1)-three)+&
         (d32*Az*Ax2-d31*Az**3)*cos(three*t1))
    vx = Ax*sin(t1)-(a23*Ax2-a24*Az2)*sin(two*t1)*two-&
         (a31*Ax**3-a32*Ax*Az2)*sin(three*t1)*three
    vy = k*Ax*cos(t1)+(b21*Ax2-b22*Az2)*cos(two*t1)*two+&
         (b31*Ax**3-b32*Ax*Az2)*cos(three*t1)*three
    vz = delta_n*(-Az*sin(t1)+d21*Ax*Az*(-sin(two*t1)*two)-&
         (d32*Az*Ax2-d31*Az**3)*sin(three*t1)*three)

    rv = [x,y,z,vx,vy,vz]

    ! convert from richardson scale, libration point centered to
    ! standard normalized coordinates wrt barycenter:
    rv(1:3) = rv(1:3) * gamma_l
    rv(4:6) = rv(4:6) * gamma_l * (lambda*w)
    rv(1)   = rv(1) + x_libpoint(libpoint)

    contains
!*******************************************************************************

    !***************************************************************************
        pure function c_n(lib,n,mu,gl) result(cn)

        !! Equations 8a, 8b in the reference.

        implicit none

        integer,intent(in)   :: lib  !! libration point (1,2,3)
        integer,intent(in)   :: n    !! the n in cn
        real(wp),intent(in)  :: mu   !! cr3bp normalized grav parameter
        real(wp),intent(in)  :: gl   !! \( \gamma_l \)
        real(wp)             :: cn   !! result

        ! Equation 8a and 8b:
        select case(lib)
        case(1); cn = (mu+(-1)**n*(one-mu)*(gl/(one-gl))**(n+1)  )/gl**3
        case(2); cn = ((-1)**n*(mu+(one-mu)*(gl/(one+gl))**(n+1)))/gl**3
        case(3); cn = (one-mu+mu*(gl/(one+gl))**(n+1)            )/gl**3
        end select

        end function c_n
    !***************************************************************************

    end subroutine halo_to_rv
!*******************************************************************************

!*******************************************************************************
!>
!  Compute the halo orbit monodromy matrix
!  (which is the state transition matrix propagated for one period)
!  The input should be the result from the [[halo_to_rv_diffcorr]] routine.

    subroutine compute_halo_monodromy_matrix(mu,rv,period,phi)

    use rk_module

    implicit none

    real(wp),intent(in)                 :: mu     !! CRTBP parameter
    real(wp),dimension(6),intent(in)    :: rv     !! halo orbit state vector
                                                  !! (normalized)
    real(wp),intent(in)                 :: period !! halo orbit period
                                                  !! (normalized)
    real(wp),dimension(6,6),intent(out) :: phi    !! monodromy matrix

    real(wp),parameter :: t0 = zero  !! initial time (normalized)
                                     !! (epoch doesn't matter for cr3bp)
    integer,parameter  :: n_steps_per_rev = 100   !! number of integration steps
                                                  !! per orbit rev

    real(wp),dimension(42)  :: x0     !! initial normalized state and STM
    real(wp),dimension(42)  :: xf     !! final normalized state and STM
    real(wp),dimension(6,6) :: phi0   !! initial STM
    integer                 :: i      !! counter
    type(rk8_10_class)      :: prop   !! integrator
    real(wp)                :: dt     !! integration time step (normalized)

    ! initial state:
    x0(1:6) = rv

    ! initial stm is the identity matrix:
    phi0 = zero
    do i = 1, 6
       phi0(i,i) = one
    end do
    x0(7:42) = pack (phi0, mask=.true.)

    ! for now, use a fixed time step:
    ! (same as was used in [[halo_to_rv_diffcorr]])
    dt = period / real(n_steps_per_rev,wp)

    ! initialize the integrator:
    call prop%initialize(42,func)

    ! propagate for one period:
    call prop%integrate(t0,x0,dt,period,xf)

    ! extract the STM:
    phi = reshape(xf(7:42), shape=[6,6])

    contains
!*******************************************************************************

    !***************************************************************************
        subroutine func(me,t,x,xdot)
        !! CRTBP derivative function (with STM)
        implicit none
        class(rk_class),intent(inout)        :: me
        real(wp),intent(in)                  :: t
        real(wp),dimension(me%n),intent(in)  :: x
        real(wp),dimension(me%n),intent(out) :: xdot

        call crtbp_derivs_with_stm(mu,x,xdot)

        end subroutine func
    !**************************************************************************

    end subroutine compute_halo_monodromy_matrix
!*******************************************************************************

!*******************************************************************************
!>
!  Compute the eigenvalues of the monodromy matrix.
!
!### Reference
!  * J.S. Parker, R.L. Anderson, "Low-Energy Lunar Trajectory Design",
!    2014. (p 79)

    subroutine compute_monodromy_matrix_eigenvalues(phi,lambda)

    use matrix_module, only: matrix_trace

    implicit none

    real(wp),dimension(6,6),intent(in)   :: phi    !! monodromy matrix
    complex(wp),dimension(6),intent(out) :: lambda !! eigenvalues of `phi`

    real(wp) :: alpha, beta, alpha2
    complex(wp) :: p, q, a, b, c

    alpha  = two - matrix_trace(6,phi)
    alpha2 = alpha*alpha
    beta   = (alpha2 - matrix_trace(6,matmul(phi,phi)))/two + one
    a      = sqrt(alpha2 - four*beta + eight)
    p      = (alpha + a) / two
    q      = (alpha - a) / two
    b      = sqrt(p*p - four)
    c      = sqrt(q*q - four)

    ! eigenvalues:
    lambda(1) = (-p + b) / two
    lambda(2) = (-p - b) / two
    lambda(3) = (-q + c) / two
    lambda(4) = (-q - c) / two
    lambda(5) = (one, zero)
    lambda(6) = (one, zero)

    end subroutine compute_monodromy_matrix_eigenvalues
!*******************************************************************************

!*******************************************************************************
!>
!  Unit test for the halo orbit module.

    subroutine halo_orbit_test()

    use celestial_body_module

    implicit none

    integer               :: libpoint  !! Libration point number: [1,2,3]
    real(wp)              :: mu1       !! grav param for primary body [km3/s2]
    real(wp)              :: mu2       !! grav param for secondary body [km3/s2]
    real(wp)              :: dist      !! distance between bodies [km]
    real(wp)              :: Az        !! halo z amplitude [km]
    integer               :: n         !! halo family: 1, 3
    real(wp)              :: t1        !! tau1 [rad]
    real(wp),dimension(6) :: rv        !! normalized state [wrt barycenter]

    write(*,*) ''
    write(*,*) '----------------'
    write(*,*) ' halo_orbit_test'
    write(*,*) '----------------'
    write(*,*) ''

    libpoint = 2
    mu1  = body_earth%mu
    mu2  = body_moon%mu
    dist = 384400.0_wp
    Az   = 10000.0_wp
    n    = 1
    t1   = 0.0_wp
    call halo_to_rv(libpoint,mu1,mu2,dist,Az,n,t1,rv)

    write(*,*) ''
    write(*,*) 'halo orbit state:'
    write(*,*) rv
    write(*,*) ''

    end subroutine halo_orbit_test
!*******************************************************************************

!*******************************************************************************
    end module halo_orbit_module
!*******************************************************************************
