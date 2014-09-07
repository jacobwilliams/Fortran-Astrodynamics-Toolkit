!*****************************************************************************************
    module rk_module
!*****************************************************************************************
!****h* FAT/rk_module
!
!  NAME
!    rk_module
!
!  DESCRIPTION
!    Runge-Kutta integration.
!
!*****************************************************************************************
    
    use kind_module,       only: wp
    use numbers_module,    only: zero
    
    implicit none
    
    private
       
    !main integration class:
    type,abstract,public :: rk_class
        private
        integer :: n = 0                               !user specified number of variables
        procedure(deriv_func),pointer :: f => null()   !user-specified derivative function
        procedure(report_func),pointer :: report => null() !user-specified report function
        contains
        procedure,non_overridable,public :: integrate  !main integration routine
        procedure(step_func),deferred :: step          !the step routine for the rk method
    end type rk_class
    
    !extend the abstract class to create an RK4 method:
    ! [all we need to do is set the step function]
    type,extends(rk_class),public :: rk4_class
        contains
        procedure :: step => rk4
    end type rk4_class   
    type,extends(rk_class),public :: rk8_10_class
        contains
        procedure :: step => rk8_10
    end type rk8_10_class       
    
    interface
    
        subroutine deriv_func(me,t,x,xdot)  !derivative function
        import :: rk_class,wp
        implicit none
            class(rk_class),intent(inout)        :: me
            real(wp),intent(in)                  :: t    
            real(wp),dimension(me%n),intent(in)  :: x    
            real(wp),dimension(me%n),intent(out) :: xdot    
        end subroutine deriv_func
    
        subroutine report_func(me,t,x)  !report function
        import :: rk_class,wp
        implicit none
            class(rk_class),intent(inout)        :: me
            real(wp),intent(in)                  :: t    
            real(wp),dimension(me%n),intent(in)  :: x    
        end subroutine report_func
        
        subroutine step_func(me,t,x,h,xf)   !rk step function
        import :: rk_class,wp
        implicit none
            class(rk_class),intent(inout)        :: me
            real(wp),intent(in)                  :: t
            real(wp),dimension(me%n),intent(in)  :: x
            real(wp),intent(in)                  :: h
            real(wp),dimension(me%n),intent(out) :: xf
        end subroutine step_func    
            
    end interface
    
    public :: rk_test !for testing
    
    contains
!*****************************************************************************************
    
!*****************************************************************************************
!****f* rk_module/integrate
! 
!  NAME
!    integrate
!
!  DESCRIPTION
!    Main integration routine for the rk_class.
!
!  SOURCE

    subroutine integrate(me,t0,x0,h,tf,xf) 
    
    implicit none
    
    class(rk_class),intent(inout)        :: me    
    real(wp),intent(in)                  :: t0    !initial time
    real(wp),dimension(me%n),intent(in)  :: x0    !initial state
    real(wp),intent(in)                  :: h     !time step
    real(wp),intent(in)                  :: tf    !final time
    real(wp),dimension(me%n),intent(out) :: xf    !final state
    
    real(wp) :: t,dt,t2
    real(wp),dimension(me%n) :: x
    logical :: last,export
    
    export = associated(me%report)
    
    if (export) call me%report(t0,x0)  !first point
   
    if (h==zero) then
        xf = x0
    else
    
        t = t0
        x = x0
        dt = h
        do
            t2 = t + dt
            last = ((dt>=zero .and. t2>=tf) .or. &  !adjust last time step
                    (dt<zero .and. t2<=tf))         !
            if (last) dt = tf-t                     !
            call me%step(t,x,dt,xf)
            if (last) exit
            if (export) call me%report(t2,xf)   !intermediate point
            x = xf
            t = t2
        end do
        
    end if
    
    if (export) call me%report(t2,xf)   !last point
    
    end subroutine integrate
!*****************************************************************************************
    
!*****************************************************************************************
!****f* rk_module/rk4
! 
!  NAME
!    rk4
!
!  DESCRIPTION
!    Take one Runge Kutta 4 integration step: t -> t+h (x -> xf)
!
!  SOURCE

    subroutine rk4(me,t,x,h,xf) 
     
    implicit none
    
    class(rk4_class),intent(inout)       :: me
    real(wp),intent(in)                  :: t   !initial time
    real(wp),dimension(me%n),intent(in)  :: x   !initial state
    real(wp),intent(in)                  :: h   !time step
    real(wp),dimension(me%n),intent(out) :: xf  !state at time t+h
    
    !local variables:
    real(wp),dimension(me%n) :: f1,f2,f3,f4
    real(wp) :: h2
    
    !parameters:
    real(wp),parameter :: half = 0.5_wp
    real(wp),parameter :: six  = 6.0_wp
    
    h2 = half*h
    
    call me%f(t,x,f1)
    call me%f(t+h2,x+h2*f1,f2)
    call me%f(t+h2,x+h2*f2,f3)
    call me%f(t+h,x+h*f3,f4)
    
    xf = x + h*(f1+f2+f2+f3+f3+f4)/six
    
    end subroutine rk4
!*****************************************************************************************

!*****************************************************************************************
!****f* rk_module/rk8_10
! 
!  NAME
!    rk8_10
!
!  DESCRIPTION
!    Take one Runge Kutta 8 integration step: t -> t+h (x -> xf)
!    This is Formula (8-10) from Reference [1].
!
!  SEE ALSO
!    [1] E. B. Shanks, "Higher Order Approximations of Runge-Kutta Type",
!        NASA Technical Note, NASA TN D-2920, Sept. 1965.
!
!  SOURCE

	subroutine rk8_10(me,t,x,h,xf)

	implicit none
	
    class(rk8_10_class),intent(inout)    :: me
    real(wp),intent(in)                  :: t   !initial time
    real(wp),dimension(me%n),intent(in)  :: x   !initial state
    real(wp),intent(in)                  :: h   !time step
    real(wp),dimension(me%n),intent(out) :: xf  !state at time t+h
    
	!local variables:
	real(wp),dimension(me%n) :: f0,f1,f2,f3,f4,f5,f6,f7,f8,f9
	
	!parameters:
	real(wp),parameter :: a1  = 4.0_wp/27.0_wp
	real(wp),parameter :: a2  = 2.0_wp/9.0_wp
	real(wp),parameter :: a3  = 1.0_wp/3.0_wp
	real(wp),parameter :: a4  = 1.0_wp/2.0_wp
	real(wp),parameter :: a5  = 2.0_wp/3.0_wp
	real(wp),parameter :: a6  = 1.0_wp/6.0_wp
	real(wp),parameter :: a8  = 5.0_wp/6.0_wp
	real(wp),parameter :: c   = 1.0_wp/840.0_wp
	real(wp),parameter :: c0  = 41.0_wp
	real(wp),parameter :: c3  = 27.0_wp
	real(wp),parameter :: c4  = 272.0_wp
	real(wp),parameter :: c5  = 27.0_wp
	real(wp),parameter :: c6  = 216.0_wp
	real(wp),parameter :: c8  = 216.0_wp
	real(wp),parameter :: c9  = 41.0_wp	
	real(wp),parameter :: aa1 = 4.0_wp/27.0_wp
	real(wp),parameter :: aa2 = 1.0_wp/18.0_wp
	real(wp),parameter :: aa3 = 1.0_wp/12.0_wp
	real(wp),parameter :: aa4 = 1.0_wp/8.0_wp
	real(wp),parameter :: aa5 = 1.0_wp/54.0_wp
	real(wp),parameter :: aa6 = 1.0_wp/4320.0_wp
	real(wp),parameter :: aa7 = 1.0_wp/20.0_wp
	real(wp),parameter :: aa8 = 1.0_wp/288.0_wp
	real(wp),parameter :: aa9 = 1.0_wp/820.0_wp
	real(wp),parameter :: b21 = 3.0_wp
	real(wp),parameter :: b32 = 3.0_wp
	real(wp),parameter :: b43 = 3.0_wp
	real(wp),parameter :: b50 = 13.0_wp
	real(wp),parameter :: b52 = -27.0_wp
	real(wp),parameter :: b53 = 42.0_wp
	real(wp),parameter :: b54 = 8.0_wp	
	real(wp),parameter :: b60 = 389.0_wp
	real(wp),parameter :: b62 = -54.0_wp
	real(wp),parameter :: b63 = 966.0_wp
	real(wp),parameter :: b64 = -824.0_wp
	real(wp),parameter :: b65 = 243.0_wp
	real(wp),parameter :: b70 = -231.0_wp
	real(wp),parameter :: b72 = 81.0_wp
	real(wp),parameter :: b73 = -1164.0_wp
	real(wp),parameter :: b74 = 656.0_wp
	real(wp),parameter :: b75 = -122.0_wp
	real(wp),parameter :: b76 = 800.0_wp
	real(wp),parameter :: b80 = -127.0_wp
	real(wp),parameter :: b82 = 18.0_wp
	real(wp),parameter :: b83 = -678.0_wp
	real(wp),parameter :: b84 = 456.0_wp
	real(wp),parameter :: b85 = -9.0_wp
	real(wp),parameter :: b86 = 576.0_wp
	real(wp),parameter :: b87 = 4.0_wp	
	real(wp),parameter :: b90 = 1481.0_wp
	real(wp),parameter :: b92 = -81.0_wp
	real(wp),parameter :: b93 = 7104.0_wp
	real(wp),parameter :: b94 = -3376.0_wp
	real(wp),parameter :: b95 = 72.0_wp
	real(wp),parameter :: b96 = -5040.0_wp
	real(wp),parameter :: b97 = -60.0_wp
	real(wp),parameter :: b98 = 720.0_wp
		
	call me%f(t,x,f0)
	call me%f(t+a1*h,x+aa1*h*f0,f1)
	call me%f(t+a2*h,x+aa2*h*(f0+b21*f1),f2)
	call me%f(t+a3*h,x+aa3*h*(f0+b32*f2),f3)
	call me%f(t+a4*h,x+aa4*h*(f0+b43*f3),f4)
	call me%f(t+a5*h,x+aa5*h*(b50*f0+b52*f2+b53*f3+b54*f4),f5)
	call me%f(t+a6*h,x+aa6*h*(b60*f0+b62*f2+b63*f3+b64*f4+b65*f5),f6)
	call me%f(t+h,x+aa7*h*(b70*f0+b72*f2+b73*f3+b74*f4+b75*f5+b76*f6),f7)
	call me%f(t+a8*h,x+aa8*h*(b80*f0+b82*f2+b83*f3+b84*f4+b85*f5+b86*f6+b87*f7),f8)
	call me%f(t+h,x+aa9*h*(b90*f0+b92*f2+b93*f3+b94*f4+b95*f5+b96*f6+b97*f7+b98*f8),f9)
	
	xf = x + h*c*(c0*f0+c3*f3+c4*f4+c5*f5+c6*f6+c8*f8+c9*f9)
	
	end subroutine rk8_10
!*****************************************************************************************

!*****************************************************************************************
!****f* rk_module/rk_test
! 
!  NAME
!    rk_test
!
!  DESCRIPTION
!    Unit test of the rk_module.
!    Integrate a two-body orbit around the Earth.
!
!  SOURCE

    subroutine rk_test()
    ! unit test for the integration routines
    
    !spacecraft propagation type:
    ! extend the rk class to include data used in the deriv routine
    type,extends(rk4_class) :: spacecraft
        real(wp) :: mu = zero      !central body gravitational parameter (km3/s2)
        integer :: fevals = 0      !number of function evaluations
        logical :: first = .true.  !first point is being exported
    end type spacecraft
    
    integer,parameter :: n=6    !number of state variables
    type(spacecraft) :: s
    real(wp) :: t0,tf,x0(n),dt,xf(n),x02(n)
    
    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' rk_test'
    write(*,*) '---------------'
    write(*,*) ''
        
    !constructor (main body is Earth):
    s = spacecraft(n=n,f=twobody,mu=398600.436233_wp,report=twobody_report)
    
    !initial conditions:
    x0 = [10000.0_wp,10000.0_wp,10000.0_wp,&   !initial state [r,v] (km,km/s)
            1.0_wp,2.0_wp,3.0_wp]
    t0 = zero       !initial time (sec)
    dt = 10.0_wp    !time step (sec)
    tf = 1000.0_wp  !final time (sec)
        
    s%fevals = 0
    s%first = .true.
    call s%integrate(t0,x0,dt,tf,xf)    !forward
    write(*,*) ''
    write(*,'(A/,*(F15.6/))') 'Final state:',xf
    
    s%fevals = 0
    s%report => null()    !disable reporting
    call s%integrate(tf,xf,-dt,t0,x02)  !backwards
    
    write(*,'(A/,*(E20.12/))') 'Error:',x02-x0
    write(*,'(A,I5)') 'Function evaluations:', s%fevals
    write(*,*) ''
    
!    contains: twobody
!*****************************************************************************************
    contains
    
    !*********************************************************
        subroutine twobody(me,t,x,xdot)
        ! derivative routine for two-body orbit propagation
        
        implicit none
        
        class(rk_class),intent(inout)        :: me
        real(wp),intent(in)                  :: t    
        real(wp),dimension(me%n),intent(in)  :: x    
        real(wp),dimension(me%n),intent(out) :: xdot    
    
        real(wp),dimension(3) :: r,v,a_grav
        real(wp) :: rmag
        
        select type (me)
        class is (spacecraft)
        
            r = x(1:3)
            v = x(4:6)
            rmag = norm2(r)
            a_grav = -me%mu/rmag**3 * r !acceleration due to gravity
            
            xdot(1:3) = v
            xdot(4:6) = a_grav
            
            me%fevals = me%fevals + 1
            
        end select
        
        end subroutine twobody
    !*********************************************************
        
    !*********************************************************
        subroutine twobody_report(me,t,x)
        !report function - write time,state to console

        implicit none
        
        class(rk_class),intent(inout)        :: me
        real(wp),intent(in)                  :: t
        real(wp),dimension(me%n),intent(in)  :: x
        
        select type (me)
        class is (spacecraft)
            if (me%first) then  !print header
                write(*,*) ''
                write(*,'(*(A15,1X))')  'time (sec)','x (km)','y (km)','z (km)',&
                                        'vx (km/s)','vy (km/s)','vz (km/s)'
                me%first = .false.
            end if
        end select

        write(*,'(*(F15.6,1X))') t,x    
                    
        end subroutine twobody_report
    !*********************************************************
        
    end subroutine rk_test
!*****************************************************************************************

!*****************************************************************************************
    end module rk_module
!*****************************************************************************************