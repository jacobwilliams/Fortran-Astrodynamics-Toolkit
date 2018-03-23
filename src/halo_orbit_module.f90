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

    use iso_fortran_env, only: wp => real64
    use numbers_module
    use crtbp_module, only: compute_libration_points,&
                            unnormalize_variables,&
                            compute_crtpb_parameter,&
                            normalize_variables

    implicit none

    private

    public :: halo_to_rv
    public :: halo_orbit_test ! test routine

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  Compute the state vector from the halo orbit approximation.
!  This will be an approximation of a halo orbit in the CR3BP system,
!  and will need to be corrected to produce a real halo orbit.

    subroutine halo_to_rv(libpoint,mu1,mu2,dist,A_z,n,t1,rv)

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

    real(wp) :: mu          !! CRTBP parameter
    integer  :: delta_n     !! 2 - n
    real(wp) :: gamma_l     !! dimensionless quantity from reference
    real(wp) :: lambda      !! linearized frequency
    real(wp) :: Ax,Ay,Az,Ax2,Az2
    real(wp) :: x,y,z,vx,vy,vz
    real(wp) :: a1,a2,a21,a22,a23,a24,a31,a32,b21,&
                b22,b31,b32,c2,c3,c4,delta,w,&
                d1,d2,d21,d3,d31,d32,k,l1,l2,s1,s2
    real(wp),dimension(3) :: x_libpoint  !! x-coordinates of the libration
                                         !! point (wrt barycenter, normalized)

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
    Ax      = sqrt((-delta-l2*Az2)/l1) ! equation 18 -- NOTE: we should check if this is feasible
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

    ! write(*,*) ''
    ! write(*,*) 'mu      ', mu
    ! write(*,*) 'gamma_l ', gamma_l
    ! write(*,*) 'w       ', w
    ! write(*,*) 'c2      ', c2
    ! write(*,*) 'c3      ', c3
    ! write(*,*) 'c4      ', c4
    ! write(*,*) 'lambda  ', lambda
    ! write(*,*) 'k       ', k
    ! write(*,*) 'delta   ', delta
    ! write(*,*) 'd1      ', d1
    ! write(*,*) 'd2      ', d2
    ! write(*,*) 'd3      ', d3
    ! write(*,*) 'a21     ', a21
    ! write(*,*) 'a23     ', a23
    ! write(*,*) 'b21     ', b21
    ! write(*,*) 's1      ', s1
    ! write(*,*) 'a22     ', a22
    ! write(*,*) 'a24     ', a24
    ! write(*,*) 'b22     ', b22
    ! write(*,*) 'd21     ', d21
    ! write(*,*) 's2      ', s2
    ! write(*,*) 'a1      ', a1
    ! write(*,*) 'a2      ', a2
    ! write(*,*) 'l1      ', l1
    ! write(*,*) 'l2      ', l2
    ! write(*,*) 'Ax2     ', Ax2
    ! write(*,*) 'Az2     ', Az2
    ! write(*,*) 'Ax      ', Ax
    ! write(*,*) 'Ay      ', Ay
    ! write(*,*) 'Az      ', Az
    ! write(*,*) 'Ax      ', Ax*dist*gamma_l, 'km'
    ! write(*,*) 'Ay      ', Ay*dist*gamma_l, 'km'
    ! write(*,*) 'Az      ', Az*dist*gamma_l, 'km'
    ! write(*,*) 'a31     ', a31
    ! write(*,*) 'a32     ', a32
    ! write(*,*) 'b31     ', b31
    ! write(*,*) 'b32     ', b32
    ! write(*,*) 'd31     ', d31
    ! write(*,*) 'd32     ', d32
    ! write(*,*) 'delta_n ', delta_n
    ! write(*,*) ''
    ! write(*,*) 'rich x:   ',x
    ! write(*,*) 'rich y:   ',y
    ! write(*,*) 'rich z:   ',z
    ! write(*,*) 'rich vx:  ',vx
    ! write(*,*) 'rich vy:  ',vy
    ! write(*,*) 'rich vz:  ',vz
    ! write(*,*) ''
    ! write(*,*) 'bary: x:  ', rv(1)
    ! write(*,*) 'bary: y:  ', rv(2)
    ! write(*,*) 'bary: z:  ', rv(3)
    ! write(*,*) 'bary: vx: ', rv(4)
    ! write(*,*) 'bary: vy: ', rv(5)
    ! write(*,*) 'bary: vz: ', rv(6)
    ! write(*,*) ''

    end subroutine halo_to_rv
!*******************************************************************************

!*******************************************************************************
!>
!  Equations 8a, 8b in the reference.

    pure function c_n(lib,n,mu,gl) result(cn)

    implicit none

    integer,intent(in)   :: lib  !! libration point (1,2,3)
    integer,intent(in)   :: n    !! the n in cn
    real(wp),intent(in)  :: mu   !! cr3bp normalized grav parameter
    real(wp),intent(in)  :: gl   !! \( \gamma_l \)
    real(wp)             :: cn

    ! Equation 8a and 8b:
    select case(lib)
    case(1); cn = (mu+(-1)**n*(one-mu)*(gl/(one-gl))**(n+1)  )/gl**3
    case(2); cn = ((-1)**n*(mu+(one-mu)*(gl/(one+gl))**(n+1)))/gl**3
    case(3); cn = (one-mu+mu*(gl/(one+gl))**(n+1)            )/gl**3
    end select

    end function c_n
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
