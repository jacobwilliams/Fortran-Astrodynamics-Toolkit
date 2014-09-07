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