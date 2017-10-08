!*****************************************************************************************
!> author: Jacob Williams
!  date: April 1, 2016
!
!  High-order variable step size Runge-Kutta integration methods.
!
!  Currently have four methods:
!   * Fehlberg 7(8)
!   * Fehlberg 8(9)
!   * Verner 8(9)
!   * Feagin 8(10)
!
!@warning This is a work in progress.

    module rk_module_variable_step

    use kind_module, only: wp
    use numbers_module

    implicit none

    private

    type,public :: stepsize_class

        !! Algorithms for adjusting the step size for variable-step
        !! Runge-Kutta integrators.

        private

        real(wp) :: hmax           = huge(one)        !! maximum allowed step size
        real(wp) :: hmin           = two*epsilon(one) !! minimum allowed step size
        real(wp) :: hfactor_reject = 1.0e-3_wp        !! minimum allowed factor for decreasing step size after refected step
        real(wp) :: hfactor_accept = 100.0_wp         !! maximum allowed factor for increasing step size after accepted step
        integer  :: accept_mode    = 1                !! method to determine if step is accepted [1,2]
        integer  :: max_attempts   = 100              !! maximum number of attempts to decrease step size before giving up

        procedure(norm_func),nopass,pointer :: norm => maxval_func
            !! routine for computing the norm of the state

        procedure(compute_h_factor_func),pointer :: compute_h_factor => compute_h_factor_v1
            !! routine for computing the h factor when computing a new step size

        contains

        private

        procedure,public :: initialize => stepsize_class_constructor
        procedure,public :: compute_stepsize
        procedure,public :: destroy => destroy_stepsize_class

    end type stepsize_class

    type,abstract,public :: rk_variable_step_class

        !! Main integration class for variable step size Runge-Kutta methods

        private

        integer :: n = 0  !! user specified number of variables
        procedure(deriv_func),pointer  :: f      => null()  !! user-specified derivative function
        procedure(report_func),pointer :: report => null()  !! user-specified report function
        procedure(event_func),pointer  :: g      => null()  !! event function (stop when this is zero)

        class(stepsize_class),allocatable :: stepsize_method  !! the method for varying the step size

        real(wp),dimension(:),allocatable :: rtol  !! relative tolerance (`size(n)`)
        real(wp),dimension(:),allocatable :: atol  !! absolute tolerance (`size(n)`)

        integer :: p = 0 !! order of the method

        integer :: hinit_method = 1 !! if automatically computing the inital step size, which
                                    !! method to use. 1 = `hstart`, 2 = `hinit`.

        integer :: num_rejected_steps = 0 !! number of rejected steps

        contains

        private

        procedure,public                 :: initialize         !! initialize the class (set n,f, and report)
        procedure,public                 :: destroy            !! destructor
        procedure,non_overridable,public :: integrate          !! main integration routine
        procedure,non_overridable,public :: integrate_to_event !! integration with event finding
        procedure(step_func),deferred    :: step               !! the step routine for the rk method
        procedure(order_func),deferred   :: order              !! returns `p`, the order of the method
        procedure :: hstart  !! for automatically computing the initial step size [this is from DDEABM]
        procedure :: hinit   !! for automatically computing the initial step size [this is from DOP853]

    end type rk_variable_step_class

    type,extends(rk_variable_step_class),public :: rkf78_class
        !! Runga-Kutta Fehlberg 7(8) method.
        contains
        procedure :: step  => rkf78
        procedure :: order => rkf78_order
    end type rkf78_class
    type,extends(rk_variable_step_class),public :: rkf89_class
        !! Runga-Kutta Fehlberg 8(9) method.
        contains
        procedure :: step  => rkf89
        procedure :: order => rkf89_order
    end type rkf89_class
    type,extends(rk_variable_step_class),public :: rkv89_class
        !! Runga-Kutta Verner 8(9) method.
        contains
        procedure :: step  => rkv89
        procedure :: order => rkv89_order
    end type rkv89_class
    type,extends(rk_variable_step_class),public :: rkf108_class
        !! Runga-Kutta Feagin 8(10) method.
        contains
        procedure :: step  => rkf108
        procedure :: order => rkf108_order
    end type rkf108_class

    abstract interface

        pure function norm_func(x) result(xmag)
        !! Vector norm function. Must return a value \( \ge 0 \).
        import :: wp
        implicit none
            real(wp),dimension(:),intent(in) :: x    !! a vector
            real(wp)                         :: xmag !! the magnitude of the vector
        end function norm_func

        pure function order_func(me) result(p)
        import :: rk_variable_step_class
        implicit none
            class(rk_variable_step_class),intent(in) :: me
            integer                    :: p !! order of the method
        end function order_func

        subroutine deriv_func(me,t,x,xdot)  !! derivative function
        import :: rk_variable_step_class,wp
        implicit none
            class(rk_variable_step_class),intent(inout)     :: me
            real(wp),intent(in)               :: t    !! time
            real(wp),dimension(:),intent(in)  :: x    !! state vector
            real(wp),dimension(:),intent(out) :: xdot !! derivative of state vector
        end subroutine deriv_func

        subroutine event_func(me,t,x,g)  !! event function
        import :: rk_variable_step_class,wp
        implicit none
            class(rk_variable_step_class),intent(inout)     :: me
            real(wp),intent(in)               :: t !! time
            real(wp),dimension(:),intent(in)  :: x !! state vector
            real(wp),intent(out)              :: g !! g(t,x). The goal is to stop the integration when g=0.
        end subroutine event_func

        subroutine report_func(me,t,x)  !! report function
        import :: rk_variable_step_class,wp
        implicit none
            class(rk_variable_step_class),intent(inout)    :: me
            real(wp),intent(in)              :: t !! time
            real(wp),dimension(:),intent(in) :: x !! state vector
        end subroutine report_func

        subroutine step_func(me,t,x,h,xf,terr)   !! rk step function
        import :: rk_variable_step_class,wp
        implicit none
            class(rk_variable_step_class),intent(inout)        :: me
            real(wp),intent(in)                  :: t    !! initial time
            real(wp),dimension(me%n),intent(in)  :: x    !! initial state vector
            real(wp),intent(in)                  :: h    !! time step \( |\Delta t| \)
            real(wp),dimension(me%n),intent(out) :: xf   !! final state vector
            real(wp),dimension(me%n),intent(out) :: terr !! truncation error estimate
        end subroutine step_func

        pure function compute_h_factor_func(me,h,tol,err,p) result(hfactor)
            import :: wp, stepsize_class
            implicit none
            class(stepsize_class),intent(in) :: me
            real(wp),intent(in)              :: h        !! current step size
            real(wp),intent(in)              :: tol      !! abs error tolerance
            real(wp),intent(in)              :: err      !! truncation error estimate
            integer,intent(in)               :: p        !! order of the method
            real(wp)                         :: hfactor  !! new step size factor (`hnew` = `h` * `hfactor`)
        end function compute_h_factor_func

    end interface

    ! public routines:
    public :: norm2_func,maxval_func

    ! for testing:
    public :: step_size_test
    public :: rk_test_variable_step

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Use intrinsic `norm2(x)` for computing the vector norm.

    pure function norm2_func(x) result(xmag)

    implicit none

    real(wp),dimension(:),intent(in) :: x
    real(wp) :: xmag

    xmag = norm2(x)

    end function norm2_func
!*****************************************************************************************

!*****************************************************************************************
!>
!  Use `maxval(abs(x))` for computing the vector norm.

    pure function maxval_func(x) result(xmag)

    implicit none

    real(wp),dimension(:),intent(in) :: x
    real(wp) :: xmag

    xmag = maxval(abs(x))

    end function maxval_func
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[stepsize_class]].
!
!@warning The `norm` and `compute_h_factor` options aren't public in the module.
!         Need to fix this.

    pure subroutine stepsize_class_constructor(me,hmin,hmax,hfactor_reject,&
                        hfactor_accept,norm,accept_mode,compute_h_factor,max_attempts)


    implicit none

    class(stepsize_class),intent(inout)       :: me
    real(wp),intent(in),optional              :: hmin             !! minimum allowed step size (>0)
    real(wp),intent(in),optional              :: hmax             !! maximum allowed step size (>0)
    real(wp),intent(in),optional              :: hfactor_reject   !! minimum allowed factor for decreasing step size after refected step (>0)
    real(wp),intent(in),optional              :: hfactor_accept   !! maximum allowed factor for decreasing step size after accepted step (>0)
    procedure(norm_func),optional             :: norm             !! the user-specified \( ||x|| \) function
    integer,intent(in),optional               :: accept_mode      !! method to determine if step is accepted [1,2]
    procedure(compute_h_factor_func),optional :: compute_h_factor !! routine for compute the `h` factor
    integer,intent(in),optional               :: max_attempts     !! max step size change attempts after rejected step

    if (present(hmin))             me%hmin             = abs(hmin)
    if (present(hmax))             me%hmax             = abs(hmax)
    if (present(hfactor_reject))   me%hfactor_reject   = abs(hfactor_reject)
    if (present(hfactor_accept))   me%hfactor_accept   = abs(hfactor_accept)
    if (present(norm))             me%norm             => norm
    if (present(compute_h_factor)) me%compute_h_factor => compute_h_factor
    if (present(accept_mode))      me%accept_mode      = accept_mode
    if (present(max_attempts))     me%max_attempts     = max_attempts

    end subroutine stepsize_class_constructor
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for [[stepsize_class]].

    subroutine destroy_stepsize_class(me)

    implicit none

    class(stepsize_class),intent(out) :: me

    end subroutine destroy_stepsize_class
!*****************************************************************************************

!*****************************************************************************************
    pure subroutine compute_stepsize(me,h,tol,err,p,hnew,accept)

    !! Compute the new step size using the specific method.

    implicit none

    class(stepsize_class),intent(in) :: me
    real(wp),intent(in)              :: h      !! current step size (<>0)
    real(wp),intent(in)              :: tol    !! abs error tolerance (>0)
    real(wp),intent(in)              :: err    !! truncation error estimate (>0)
    integer,intent(in)               :: p      !! order of the method
    real(wp),intent(out)             :: hnew   !! new step size (<>0)
    logical,intent(out)              :: accept !! if the step is accepted

    real(wp) :: hfactor  !! step size factor (>0)
    real(wp),parameter :: small = ten * epsilon(one) !! small error value

    if (err<=small) then ! the error is extremely small

        hfactor = me%hfactor_accept
        accept = .true.

    else

        ! compute base factor based on the selected formula:
        hfactor = abs(me%compute_h_factor(h,tol,err,p))

        ! if the step is to be accepted:
        select case (me%accept_mode)
        case(1) !algorithm 17.12
            accept = (hfactor>=one)
        case(2) !algorithm 17.13
            accept = (err<=tol)
        end select

        !...notes:
        ! see: L. Shampine "Some Practical Runge-Kutta Formulas",
        !      Mathematics of Computation, 46(173), Jan 1986.
        ! different conditions for satisfying error conditions:
        !  ||err|| <= tol   -- Error per setp (EPS)
        !  ||err|| <= h*tol -- Error per unit step (EPUS)

        !compute the actual hfactor based on the limits:
        if (accept) then
            hfactor = min(me%hfactor_accept, hfactor)
        else
            hfactor = max(me%hfactor_reject, hfactor)
        end if

    end if

    ! compute the new step size (enforce min/max bounds & add sign):
    hnew = sign(max(me%hmin,min(me%hmax,abs(h)*hfactor)),h)

    end subroutine compute_stepsize
!*****************************************************************************************



    ! hfactor = 0.9_wp * abs(tol/err)  **(one/real(p,wp))   ! v1
    ! hfactor = 0.9_wp * abs(tol/err)  **(one/real(p+1,wp)) ! v2
    ! hfactor =          abs(tol/err)  **(one/real(p+1,wp)) ! v3
    ! hfactor = 0.9_wp * abs(tol*h/err)**(one/real(p,wp))   ! v4
    ! hfactor =          abs(tol*h/err)**(one/real(p,wp))   ! v5
    !
    !
    ! ... could be reduced to these parameters ...
    !
    ! USE_TOL_TIMES_H = true or false
    ! HFACTOR_SCALE   = 0.9_wp
    ! EXPONENT        = p or p+1
    !
    ! if (USE_TOL_TIMES_H) then
    !     hfactor = HFACTOR_SCALE * abs(tol*h/err)**(one/real(EXPONENT,wp))
    ! else
    !     hfactor = HFACTOR_SCALE * abs(tol/err)**(one/real(EXPONENT,wp))
    ! end if



!*****************************************************************************************
    pure function compute_h_factor_v1(me,h,tol,err,p) result(hfactor)

    !! Hull, Enright, Fellen and Sedgwick method

    implicit none

    class(stepsize_class),intent(in) :: me
    real(wp),intent(in)              :: h        !! current step size
    real(wp),intent(in)              :: tol      !! abs error tolerance
    real(wp),intent(in)              :: err      !! truncation error estimate
    integer,intent(in)               :: p        !! order of the method
    real(wp)                         :: hfactor  !! step size factor

    hfactor = 0.9_wp * abs(tol/err)**(one/real(p,wp))

    !..note: algorithm 17.13 has p+1 in the exponent...

    end function compute_h_factor_v1
!*****************************************************************************************

!*****************************************************************************************
    pure function compute_h_factor_v2(me,h,tol,err,p) result(hfactor)

    !! Stoer and Bulirsch method

    implicit none

    class(stepsize_class),intent(in) :: me
    real(wp),intent(in)              :: h        !! current step size
    real(wp),intent(in)              :: tol      !! abs error tolerance
    real(wp),intent(in)              :: err      !! truncation error estimate
    integer,intent(in)               :: p        !! order of the method
    real(wp)                         :: hfactor  !! step size factor

    hfactor = abs(tol/err)**(one/real(p+1,wp))

    end function compute_h_factor_v2
!*****************************************************************************************

!*****************************************************************************************
    pure function compute_h_factor_v3(me,h,tol,err,p) result(hfactor)

    !! Stoer and Bulirsch method version 2

    implicit none

    class(stepsize_class),intent(in) :: me
    real(wp),intent(in)              :: h        !! current step size
    real(wp),intent(in)              :: tol      !! abs error tolerance
    real(wp),intent(in)              :: err      !! truncation error estimate
    integer,intent(in)               :: p        !! order of the method
    real(wp)                         :: hfactor  !! step size factor

    hfactor = 0.9_wp * abs(tol*h/err)**(one/real(p,wp))

    end function compute_h_factor_v3
!*****************************************************************************************

!*****************************************************************************************
    pure function compute_h_factor_v4(me,h,tol,err,p) result(hfactor)

    !! Algorithm 17.12 method

    implicit none

    class(stepsize_class),intent(in) :: me
    real(wp),intent(in)              :: h        !! current step size
    real(wp),intent(in)              :: tol      !! abs error tolerance
    real(wp),intent(in)              :: err      !! truncation error estimate
    integer,intent(in)               :: p        !! order of the method
    real(wp)                         :: hfactor  !! step size factor

    hfactor = abs(tol*h/err)**(one/real(p,wp))

    end function compute_h_factor_v4
!*****************************************************************************************

!*****************************************************************************************
    pure function compute_h_factor_v5(me,h,tol,err,p) result(hfactor)

    !! Algorithm 17.13 method
    !!
    !! See: L.F. Shampine, H. A. Watts, "Global Error Estimation for Ordinary Differential Equations",
    !!      ACM Transactions on Mathematical Software (TOMS),
    !!      Volume 2 Issue 2, June 1976

    implicit none

    class(stepsize_class),intent(in) :: me
    real(wp),intent(in)              :: h        !! current step size
    real(wp),intent(in)              :: tol      !! abs error tolerance
    real(wp),intent(in)              :: err      !! truncation error estimate
    integer,intent(in)               :: p        !! order of the method
    real(wp)                         :: hfactor  !! step size factor

    hfactor = 0.9_wp * abs(tol/err)**(one/real(p+1,wp))

    end function compute_h_factor_v5
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize the [[rk_variable_step_class]].

    subroutine initialize(me,n,f,rtol,atol,stepsize_method,hinit_method,report,g)

    implicit none

    class(rk_variable_step_class),intent(inout) :: me
    integer,intent(in)                          :: n               !! number of equations
    procedure(deriv_func)                       :: f               !! derivative function
    real(wp),dimension(:),intent(in)            :: rtol            !! relative tolerance (if size=1, then same tol used for all equations)
    real(wp),dimension(:),intent(in)            :: atol            !! absolute tolerance (if size=1, then same tol used for all equations)
    class(stepsize_class),intent(in)            :: stepsize_method !! method for varying the step size
    integer,intent(in),optional                 :: hinit_method    !! which method to use for automatic initial step size computation.
                                                                   !! 1 = use `hstart`, 2 = use `hinit`.
    procedure(report_func),optional             :: report          !! for reporting the steps
    procedure(event_func),optional              :: g               !! for stopping at an event

    call me%destroy()

    me%n = n
    me%f => f

    allocate(me%rtol(n))
    allocate(me%atol(n))
    if (size(rtol)==1) then
        me%rtol = rtol(1) !use this for all equations
    else if (size(rtol)==n) then
        me%rtol = rtol
    else
        error stop 'invalid size for rtol array.'
    end if
    if (size(atol)==1) then
        me%atol = atol(1) !use this for all equations
    else if (size(atol)==n) then
        me%atol = atol
    else
        error stop 'invalid size for atol array.'
    end if

    if (present(hinit_method)) me%hinit_method = hinit_method

    if (present(report)) me%report => report
    if (present(g))      me%g      => g

    allocate(me%stepsize_method, source=stepsize_method)

    me%num_rejected_steps = 0

    end subroutine initialize
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for [[rk_variable_step_class]].

    subroutine destroy(me)

    implicit none

    class(rk_variable_step_class),intent(out) :: me

    end subroutine destroy
!*****************************************************************************************

!*****************************************************************************************
!>
!  Main integration routine for the [[rk_variable_step_class]].

    subroutine integrate(me,t0,x0,h,tf,xf)

    implicit none

    class(rk_variable_step_class),intent(inout)     :: me
    real(wp),intent(in)               :: t0    !! initial time
    real(wp),dimension(:),intent(in)  :: x0    !! initial state
    real(wp),intent(in)               :: h     !! initial abs(time step)
    real(wp),intent(in)               :: tf    !! final time
    real(wp),dimension(:),intent(out) :: xf    !! final state

    real(wp) :: t,dt,t2,err,tol,dt_new
    real(wp),dimension(me%n) :: x,terr,etol,xp0
    logical :: last,export,accept
    integer :: i,p

    if (.not. associated(me%f)) error stop 'Error in integrate: f is not associated.'

    me%num_rejected_steps = 0
    export = associated(me%report)

    if (export) call me%report(t0,x0)  !first point

    if (t0==tf) then
        xf = x0
    else

        t = t0
        x = x0

        if (h==zero) then
            ! compute an appropriate initial step size:
            ! WARNING: this may not be working in all cases .....
            etol = me%rtol * me%stepsize_method%norm(x0) + me%atol
            call me%f(t0,x0,xp0)  ! get initial dx/dt
            select case (me%hinit_method)
            case(1)
                call me%hstart(t0,tf,x0,xp0,etol,dt)
            case(2)
                dt = me%hinit(t0,x0,sign(one,tf-t0),xp0,me%stepsize_method%hmax,me%atol,me%rtol)
            case default
                error stop 'invalid hinit_method selection'
            end select
            !write(*,*) 'inital step size: ',dt
        else
            ! user-specified initial step size:
            dt = sign(h,tf-t0)  ! (correct sign)
        end if

        p = me%order()     !order of the method
        do
            t2 = t + dt
            last = ((dt>=zero .and. t2>=tf) .or. &  !adjust last time step
                    (dt<zero .and. t2<=tf))         !
            if (last) dt = tf-t                     !

            do i=0,me%stepsize_method%max_attempts

                ! take a step:
                call me%step(t,x,dt,xf,terr)

                ! evaluate error and compute new step size:
                err = me%stepsize_method%norm(terr)
                tol = me%stepsize_method%norm( me%rtol * xf + me%atol )
                call me%stepsize_method%compute_stepsize(dt,tol,err,p,dt_new,accept)
                dt = dt_new

                if (accept) then
                    !accept this step
                    exit
                else
                    !step is rejected, repeat step with new dt
                    me%num_rejected_steps = me%num_rejected_steps + 1

                    !note: if we have reached the min step size, and the error
                    !is still too large, we can't proceed.
                    if (i>=me%stepsize_method%max_attempts) then
                        error stop 'error: too many attempts to reduce step size.'
                    end if
                    if (abs(dt) <= abs(me%stepsize_method%hmin)) then
                        error stop 'warning: min step size.'
                    end if

                    !......
                    !... if we have two rejected steps and the step size hasn't changed..
                    !    then we need to abort, since no progress is being made...
                    !......

                    last = ((dt>=zero .and. t2>=tf) .or. &  !adjust last time step
                            (dt<zero .and. t2<=tf))         !
                    if (last) dt = tf-t                     !
                    t2 = t + dt

                end if

            end do

            if (last) exit
            if (export) call me%report(t2,xf)   !intermediate point
            x = xf
            t = t2
        end do

    end if

    if (export) call me%report(tf,xf)   !last point

    end subroutine integrate
!*****************************************************************************************

!*****************************************************************************************
!>
!  Event-finding integration routine for the [[rk_variable_step_class]].
!  Integrates until g(t,x)=0, or until t=tf (whichever happens first).
!
!@note There are some efficiency improvements that could be made here.
!      This is a work in progress.

    subroutine integrate_to_event(me,t0,x0,h,tmax,tol,tf,xf,gf)

    use brent_module

    implicit none

    class(rk_variable_step_class),intent(inout) :: me
    real(wp),intent(in)                  :: t0      !! initial time
    real(wp),dimension(me%n),intent(in)  :: x0      !! initial state
    real(wp),intent(in)                  :: h       !! abs(time step)
    real(wp),intent(in)                  :: tmax    !! max final time if event not located
    real(wp),intent(in)                  :: tol     !! function tolerance for root finding
    real(wp),intent(out)                 :: tf      !! actual final time reached
    real(wp),dimension(me%n),intent(out) :: xf      !! final state (at tf)
    real(wp),intent(out)                 :: gf      !! g value at tf

    real(wp),dimension(me%n) :: etol,xp0
    real(wp),dimension(me%n) :: x,g_xf
    real(wp),dimension(me%n) :: terr !! truncation error estimate
    integer :: i,p,iflag
    real(wp) :: t,dt,t2,ga,gb,dt_root,dum,err,dt_new,stol
    logical :: first,last,export,accept
    procedure(report_func),pointer :: report
    type(brent_class) :: solver

    if (.not. associated(me%f)) error stop 'Error in integrate_to_event: f is not associated.'
    if (.not. associated(me%g)) error stop 'Error in integrate_to_event: g is not associated.'

    me%num_rejected_steps = 0
    export = associated(me%report)

    if (export) call me%report(t0,x0)  !first point

    if (t0==tmax) then
        xf = x0
        tf = t0
        call me%g(t0,x0,gf)
    else

        first = .true.
        t = t0
        x = x0
        call me%g(t,x,ga)     !evaluate event function

        if (h==zero) then
            ! compute an appropriate initial step size:
            ! WARNING: this may not be working in all cases .....
            etol = me%rtol * me%stepsize_method%norm(x0) + me%atol
            call me%f(t0,x0,xp0)  ! get initial dx/dt
            select case (me%hinit_method)
            case(1)
                call me%hstart(t0,tmax,x0,xp0,etol,dt)
            case(2)
                dt = me%hinit(t0,x0,sign(one,tmax-t0),xp0,me%stepsize_method%hmax,me%atol,me%rtol)
            case default
                error stop 'invalid hinit_method selection'
            end select
        else
            ! user-specified initial step size:
            dt = sign(h,tmax-t0)  ! (correct sign)
        end if

        p = me%order()     !order of the method
        do

            t2 = t + dt
            last = ((dt>=zero .and. t2>=tmax) .or. &  !adjust last time step
                    (dt<zero .and. t2<=tmax))         !
            if (last) then
                dt = tmax-t
                t2 = tmax
            end if

            do i=0,me%stepsize_method%max_attempts

                ! take a step:
                call me%step(t,x,dt,xf,terr)

                ! evaluate error and compute new step size:
                err = me%stepsize_method%norm(terr)
                stol = me%stepsize_method%norm( me%rtol * xf + me%atol )
                call me%stepsize_method%compute_stepsize(dt,stol,err,p,dt_new,accept)
                dt = dt_new

                if (accept) then
                    !accept this step
                    exit
                else
                    !step is rejected, repeat step with new dt
                    me%num_rejected_steps = me%num_rejected_steps + 1

                    !note: if we have reached the min step size, and the error
                    !is still too large, we can't proceed.
                    if (i>=me%stepsize_method%max_attempts) then
                        error stop 'error: too many attempts to reduce step size.'
                    end if
                    if (abs(dt) <= abs(me%stepsize_method%hmin)) then
                        error stop 'warning: min step size.'
                    end if

                    !......
                    !... if we have two rejected steps and the step size hasn't changed..
                    !    then we need to abort, since no progress is being made...
                    !......

                    last = ((dt>=zero .and. t2>=tmax) .or. &  !adjust last time step
                            (dt<zero .and. t2<=tmax))         !
                    if (last) then
                        dt = tmax-t
                        t2 = tmax
                    else
                        t2 = t + dt
                    end if

                end if

            end do

            call me%g(t2,xf,gb)     !evaluate event function

            if (first .and. abs(ga)<=tol) then

                !we ignore a root at t0 after the first step
                if (abs(gb)<=tol) then !check this one since it could have landed on a root
                    gf = gb
                    tf = t2
                    exit
                else
                    if (last) then  !exiting without having found a root
                        tf = t2
                        gf = gb
                        exit
                    end if
                    if (export) call me%report(t2,xf)   !intermediate point
                    x = xf
                    t = t2
                    ga = gb
                end if

            elseif (ga*gb<=zero) then !there is a root somewhere on [t,t+dt]

                !find the root:
                call solver%set_function(solver_func)
                call solver%find_zero(zero,dt,tol,dt_root,dum,iflag,ga,gb)
                t2 = t + dt_root
                gf = solver_func(solver,dt_root)
                tf = t2
                xf = g_xf !computed in the solver function
                exit

            else  !no root yet, continue

                if (last) then  !exiting without having found a root
                    tf = t2
                    gf = gb
                    exit
                end if
                if (export) call me%report(t2,xf)   !intermediate point
                x = xf
                t = t2
                ga = gb

            end if

            if (first) first = .false.
            if (last) exit
            x = xf
            t = t2
        end do

    end if

    if (export) call me%report(tf,xf)   !last point

    contains

        function solver_func(this,delt) result(g)

        !! root solver function. The input is the dt offset from time t.

        implicit none

        class(brent_class),intent(inout) :: this
        real(wp),intent(in) :: delt  !! from [0 to dt]
        real(wp) :: g

        real(wp),dimension(me%n) :: terr !! truncation error estimate

        !take a step from t to t+delt and evaluate g function:
        ! [we don't check the error because we are within a
        !  step that was already accepted, so it should be ok]
        call me%step(t,x,delt,g_xf,terr)
        call me%g(t+delt,g_xf,g)

        end function solver_func

    end subroutine integrate_to_event
!*****************************************************************************************

!*****************************************************************************************
!>
!  Fehlberg's 7(8) algorithm.
!
!### Reference
!  * E. Fehlberg, "Classical Fifth-, Sixth-, Seventh-, and Eighth-Order
!    Runge-Kutta Formulas with Stepsize Control",
!   [NASA TR R-2870](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19680027281_1968027281.pdf).

    subroutine rkf78(me,t,x,h,xf,terr)

    implicit none

    class(rkf78_class),intent(inout)     :: me
    real(wp),intent(in)                  :: t
    real(wp),dimension(me%n),intent(in)  :: x
    real(wp),intent(in)                  :: h
    real(wp),dimension(me%n),intent(out) :: xf
    real(wp),dimension(me%n),intent(out) :: terr

    real(wp),parameter :: a1  = 2.0_wp/27.0_wp
    real(wp),parameter :: a2  = 1.0_wp/9.0_wp
    real(wp),parameter :: a3  = 1.0_wp/6.0_wp
    real(wp),parameter :: a4  = 5.0_wp/12.0_wp
    real(wp),parameter :: a5  = 1.0_wp/2.0_wp
    real(wp),parameter :: a6  = 5.0_wp/6.0_wp
    real(wp),parameter :: a7  = 1.0_wp/6.0_wp
    real(wp),parameter :: a8  = 2.0_wp/3.0_wp
    real(wp),parameter :: a9  = 1.0_wp/3.0_wp
    !real(wp),parameter :: a10 = 1.0_wp
    !real(wp),parameter :: a12 = 1.0_wp

    real(wp),parameter :: b10  = 2.0_wp/27.0_wp
    real(wp),parameter :: b20  = 1.0_wp/36.0_wp
    real(wp),parameter :: b21  = 1.0_wp/12.0_wp
    real(wp),parameter :: b30  = 1.0_wp/24.0_wp
    real(wp),parameter :: b32  = 1.0_wp/8.0_wp
    real(wp),parameter :: b40  = 5.0_wp/12.0_wp
    real(wp),parameter :: b42  = -25.0_wp/16.0_wp
    real(wp),parameter :: b43  = 25.0_wp/16.0_wp
    real(wp),parameter :: b50  = 1.0_wp/20.0_wp
    real(wp),parameter :: b53  = 1.0_wp/4.0_wp
    real(wp),parameter :: b54  = 1.0_wp/5.0_wp
    real(wp),parameter :: b60  = -25.0_wp/108.0_wp
    real(wp),parameter :: b63  = 125.0_wp/108.0_wp
    real(wp),parameter :: b64  = -65.0_wp/27.0_wp
    real(wp),parameter :: b65  = 125.0_wp/54.0_wp
    real(wp),parameter :: b70  = 31.0_wp/300.0_wp
    real(wp),parameter :: b74  = 61.0_wp/225.0_wp
    real(wp),parameter :: b75  = -2.0_wp/9.0_wp
    real(wp),parameter :: b76  = 13.0_wp/900.0_wp
    real(wp),parameter :: b80  = 2.0_wp
    real(wp),parameter :: b83  = -53.0_wp/6.0_wp
    real(wp),parameter :: b84  = 704.0_wp/45.0_wp
    real(wp),parameter :: b85  = -107.0_wp/9.0_wp
    real(wp),parameter :: b86  = 67.0_wp/90.0_wp
    real(wp),parameter :: b87  = 3.0_wp
    real(wp),parameter :: b90  = -91.0_wp/108.0_wp
    real(wp),parameter :: b93  = 23.0_wp/108.0_wp
    real(wp),parameter :: b94  = -976.0_wp/135.0_wp
    real(wp),parameter :: b95  = 311.0_wp/54.0_wp
    real(wp),parameter :: b96  = -19.0_wp/60.0_wp
    real(wp),parameter :: b97  = 17.0_wp/6.0_wp
    real(wp),parameter :: b98  = -1.0_wp/12.0_wp
    real(wp),parameter :: b100 = 2383.0_wp/4100.0_wp
    real(wp),parameter :: b103 = -341.0_wp/164.0_wp
    real(wp),parameter :: b104 = 4496.0_wp/1025.0_wp
    real(wp),parameter :: b105 = -301.0_wp/82.0_wp
    real(wp),parameter :: b106 = 2133.0_wp/4100.0_wp
    real(wp),parameter :: b107 = 45.0_wp/82.0_wp
    real(wp),parameter :: b108 = 45.0_wp/164.0_wp
    real(wp),parameter :: b109 = 18.0_wp/41.0_wp
    real(wp),parameter :: b110 = 3.0_wp/205.0_wp
    real(wp),parameter :: b115 = -6.0_wp/41.0_wp
    real(wp),parameter :: b116 = -3.0_wp/205.0_wp
    real(wp),parameter :: b117 = -3.0_wp/41.0_wp
    real(wp),parameter :: b118 = 3.0_wp/41.0_wp
    real(wp),parameter :: b119 = 6.0_wp/41.0_wp
    real(wp),parameter :: b120 = -1777.0_wp/4100.0_wp
    real(wp),parameter :: b123 = -341.0_wp/164.0_wp
    real(wp),parameter :: b124 = 4496.0_wp/1025.0_wp
    real(wp),parameter :: b125 = -289.0_wp/82.0_wp
    real(wp),parameter :: b126 = 2193.0_wp/4100.0_wp
    real(wp),parameter :: b127 = 51.0_wp/82.0_wp
    real(wp),parameter :: b128 = 33.0_wp/164.0_wp
    real(wp),parameter :: b129 = 12.0_wp/41.0_wp
    !real(wp),parameter :: b1211 = 1.0_wp

    real(wp),parameter :: c5  = 34.0_wp/105.0_wp
    real(wp),parameter :: c6  = 9.0_wp/35.0_wp
    real(wp),parameter :: c7  = 9.0_wp/35.0_wp
    real(wp),parameter :: c8  = 9.0_wp/280.0_wp
    real(wp),parameter :: c9  = 9.0_wp/280.0_wp
    real(wp),parameter :: c11 = 41.0_wp/840.0_wp
    real(wp),parameter :: c12 = 41.0_wp/840.0_wp

    real(wp),dimension(me%n) :: f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12

    if (h==zero) then
        xf = x
        terr = zero
        return
    end if

    call me%f(t,x,f0)
    call me%f(t+h*a1,x+f0*b10*h,f1)
    call me%f(t+h*a2,x+(f0*b20+f1*b21)*h,f2)
    call me%f(t+h*a3,x+(f0*b30+f2*b32)*h,f3)
    call me%f(t+h*a4,x+(f0*b40+f2*b42+f3*b43)*h,f4)
    call me%f(t+h*a5,x+(f0*b50+f3*b53+f4*b54)*h,f5)
    call me%f(t+h*a6,x+(f0*b60+f3*b63+f4*b64+f5*b65)*h,f6)
    call me%f(t+h*a7,x+(f0*b70+f4*b74+f5*b75+f6*b76)*h,f7)
    call me%f(t+h*a8,x+(f0*b80+f3*b83+f4*b84+f5*b85+f6*b86+&
                f7*b87)*h,f8)
    call me%f(t+h*a9,x+(f0*b90+f3*b93+f4*b94+f5*b95+f6*b96+&
                f7*b97+f8*b98)*h,f9)
    call me%f(t+h,x+(f0*b100+f3*b103+f4*b104+f5*b105+&
                f6*b106+f7*b107+f8*b108+f9*b109)*h,f10)
    call me%f(t,x+(f0*b110+f5*b115+f6*b116+f7*b117+f8*b118+&
                f9*b119)*h,f11)
    call me%f(t+h,x+(f0*b120+f3*b123+f4*b124+f5*b125+f6*b126+&
                f7*b127+f8*b128+f9*b129+f11)*h,f12)

    xf = x + h*(f5*c5+f6*c6+f7*c7+f8*c8+f9*c9+f11*c11+f12*c12)

    terr = (41.0_wp/840.0_wp)*(f0+f10-f11-f12)

    end subroutine rkf78
!*****************************************************************************************

!*****************************************************************************************
!>
!  Fehlberg 8(9) method.
!
!### Reference
!  * E. Fehlberg, "Classical Fifth-, Sixth-, Seventh-, and Eighth-Order
!    Runge-Kutta Formulas with Stepsize Control",
!   [NASA TR R-2870](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19680027281_1968027281.pdf).

    subroutine rkf89(me,t,x,h,xf,terr)

    implicit none

    class(rkf89_class),intent(inout)     :: me
    real(wp),intent(in)                  :: t    !! initial time
    real(wp),dimension(me%n),intent(in)  :: x    !! initial state
    real(wp),intent(in)                  :: h    !! time step
    real(wp),dimension(me%n),intent(out) :: xf   !! state at time `t+h`
    real(wp),dimension(me%n),intent(out) :: terr !! truncation error estimate

    real(wp),parameter :: a1  = 0.44368940376498183109599404281370_wp
    real(wp),parameter :: a2  = 0.66553410564747274664399106422055_wp
    real(wp),parameter :: a3  = 0.99830115847120911996598659633083_wp
    real(wp),parameter :: a4  = 0.3155_wp
    real(wp),parameter :: a5  = 0.50544100948169068626516126737384_wp
    real(wp),parameter :: a6  = 0.17142857142857142857142857142857_wp
    real(wp),parameter :: a7  = 0.82857142857142857142857142857143_wp
    real(wp),parameter :: a8  = 0.66543966121011562534953769255586_wp
    real(wp),parameter :: a9  = 0.24878317968062652069722274560771_wp
    real(wp),parameter :: a10 = 0.1090_wp
    real(wp),parameter :: a11 = 0.8910_wp
    real(wp),parameter :: a12 = 0.3995_wp
    real(wp),parameter :: a13 = 0.6005_wp
    real(wp),parameter :: a14 = 1.0_wp
    real(wp),parameter :: a16 = 1.0_wp

    real(wp),parameter :: b1    = 0.44368940376498183109599404281370_wp
    real(wp),parameter :: b20   = 0.16638352641186818666099776605514_wp
    real(wp),parameter :: b21   = 0.49915057923560455998299329816541_wp
    real(wp),parameter :: b30   = 0.24957528961780227999149664908271_wp
    real(wp),parameter :: b32   = 0.74872586885340683997448994724812_wp
    real(wp),parameter :: b40   = 0.20661891163400602426556710393185_wp
    real(wp),parameter :: b42   = 0.17707880377986347040380997288319_wp
    real(wp),parameter :: b43   = -0.68197715413869494669377076815048e-1_wp
    real(wp),parameter :: b50   = 0.10927823152666408227903890926157_wp
    real(wp),parameter :: b53   = 0.40215962642367995421990563690087e-2_wp
    real(wp),parameter :: b54   = 0.39214118169078980444392330174325_wp
    real(wp),parameter :: b60   = 0.98899281409164665304844765434355e-1_wp
    real(wp),parameter :: b63   = 0.35138370227963966951204487356703e-2_wp
    real(wp),parameter :: b64   = 0.12476099983160016621520625872489_wp
    real(wp),parameter :: b65   = -0.55745546834989799643742901466348e-1_wp
    real(wp),parameter :: b70   = -0.36806865286242203724153101080691_wp
    real(wp),parameter :: b74   = -0.22273897469476007645024020944166e+1_wp
    real(wp),parameter :: b75   = 0.13742908256702910729565691245744e+1_wp
    real(wp),parameter :: b76   = 0.20497390027111603002159354092206e+1_wp
    real(wp),parameter :: b80   = 0.45467962641347150077351950603349e-1_wp
    real(wp),parameter :: b85   = 0.32542131701589147114677469648853_wp
    real(wp),parameter :: b86   = 0.28476660138527908888182420573687_wp
    real(wp),parameter :: b87   = 0.97837801675979152435868397271099e-2_wp
    real(wp),parameter :: b90   = 0.60842071062622057051094145205182e-1_wp
    real(wp),parameter :: b95   = -0.21184565744037007526325275251206e-1_wp
    real(wp),parameter :: b96   = 0.19596557266170831957464490662983_wp
    real(wp),parameter :: b97   = -0.42742640364817603675144835342899e-2_wp
    real(wp),parameter :: b98   = 0.17434365736814911965323452558189e-1_wp
    real(wp),parameter :: b100  = 0.54059783296931917365785724111182e-1_wp
    real(wp),parameter :: b106  = 0.11029825597828926530283127648228_wp
    real(wp),parameter :: b107  = -0.12565008520072556414147763782250e-2_wp
    real(wp),parameter :: b108  = 0.36790043477581460136384043566339e-2_wp
    real(wp),parameter :: b109  = -0.57780542770972073040840628571866e-1_wp
    real(wp),parameter :: b110  = 0.12732477068667114646645181799160_wp
    real(wp),parameter :: b117  = 0.11448805006396105323658875721817_wp
    real(wp),parameter :: b118  = 0.28773020709697992776202201849198_wp
    real(wp),parameter :: b119  = 0.50945379459611363153735885079465_wp
    real(wp),parameter :: b1110 = -0.14799682244372575900242144449640_wp
    real(wp),parameter :: b120  = -0.36526793876616740535848544394333e-2_wp
    real(wp),parameter :: b125  = 0.81629896012318919777819421247030e-1_wp
    real(wp),parameter :: b126  = -0.38607735635693506490517694343215_wp
    real(wp),parameter :: b127  = 0.30862242924605106450474166025206e-1_wp
    real(wp),parameter :: b128  = -0.58077254528320602815829374733518e-1_wp
    real(wp),parameter :: b129  = 0.33598659328884971493143451362322_wp
    real(wp),parameter :: b1210 = 0.41066880401949958613549622786417_wp
    real(wp),parameter :: b1211 = -0.11840245972355985520633156154536e-1_wp
    real(wp),parameter :: b130  =  -0.12375357921245143254979096135669e+1_wp
    real(wp),parameter :: b135  =  -0.24430768551354785358734861366763e+2_wp
    real(wp),parameter :: b136  =  0.54779568932778656050436528991173_wp
    real(wp),parameter :: b137  =  -0.44413863533413246374959896569346e+1_wp
    real(wp),parameter :: b138  =  0.10013104813713266094792617851022e+2_wp
    real(wp),parameter :: b139  =  -0.14995773102051758447170985073142e+2_wp
    real(wp),parameter :: b1310 =  0.58946948523217013620824539651427e+1_wp
    real(wp),parameter :: b1311 =  0.17380377503428984877616857440542e+1_wp
    real(wp),parameter :: b1312 =  0.27512330693166730263758622860276e+2_wp
    real(wp),parameter :: b140  = -0.35260859388334522700502958875588_wp
    real(wp),parameter :: b145  = -0.18396103144848270375044198988231_wp
    real(wp),parameter :: b146  = -0.65570189449741645138006879985251_wp
    real(wp),parameter :: b147  = -0.39086144880439863435025520241310_wp
    real(wp),parameter :: b148  = 0.26794646712850022936584423271209_wp
    real(wp),parameter :: b149  = -0.10383022991382490865769858507427e+1_wp
    real(wp),parameter :: b1410 = 0.16672327324258671664727346168501e+1_wp
    real(wp),parameter :: b1411 = 0.49551925855315977067732967071441_wp
    real(wp),parameter :: b1412 = 0.11394001132397063228586738141784e+1_wp
    real(wp),parameter :: b1413 = 0.51336696424658613688199097191534e-1_wp
    real(wp),parameter :: b150  = 0.10464847340614810391873002406755e-2_wp
    real(wp),parameter :: b158  = -0.67163886844990282237778446178020e-2_wp
    real(wp),parameter :: b159  = 0.81828762189425021265330065248999e-2_wp
    real(wp),parameter :: b1510 = -0.42640342864483347277142138087561e-2_wp
    real(wp),parameter :: b1511 = 0.28009029474168936545976331153703e-3_wp
    real(wp),parameter :: b1512 = -0.87835333876238676639057813145633e-2_wp
    real(wp),parameter :: b1513 = 0.10254505110825558084217769664009e-1_wp
    real(wp),parameter :: b160  = -0.13536550786174067080442168889966e+1_wp
    real(wp),parameter :: b165  = -0.18396103144848270375044198988231_wp
    real(wp),parameter :: b166  = -0.65570189449741645138006879985251_wp
    real(wp),parameter :: b167  = -0.39086144880439863435025520241310_wp
    real(wp),parameter :: b168  =  0.27466285581299925758962207732989_wp
    real(wp),parameter :: b169  = -0.10464851753571915887035188572676e+1_wp
    real(wp),parameter :: b1610 =  0.16714967667123155012004488306588e+1_wp
    real(wp),parameter :: b1611 =  0.49523916825841808131186990740287_wp
    real(wp),parameter :: b1612 =  0.11481836466273301905225795954930e+1_wp
    real(wp),parameter :: b1613 =  0.41082191313833055603981327527525e-1_wp
    real(wp),parameter :: b1615 =  1.0_wp

    real(wp),parameter :: c0  = 0.32256083500216249913612900960247e-1_wp
    real(wp),parameter :: c8  = 0.25983725283715403018887023171963_wp
    real(wp),parameter :: c9  = 0.92847805996577027788063714302190e-1_wp
    real(wp),parameter :: c10 = 0.16452339514764342891647731842800_wp
    real(wp),parameter :: c11 = 0.17665951637860074367084298397547_wp
    real(wp),parameter :: c12 = 0.23920102320352759374108933320941_wp
    real(wp),parameter :: c13 = 0.39484274604202853746752118829325e-2_wp
    real(wp),parameter :: c14 = 0.30726495475860640406368305522124e-1_wp

    real(wp),dimension(me%n) :: f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16

    if (h==zero) then
        xf = x
        terr = zero
        return
    end if

    call me%f(t,x,f0)
    call me%f(t+h*a1,x+f0*b1*h,f1)
    call me%f(t+h*a2,x+(f0*b20+f1*b21)*h,f2)
    call me%f(t+h*a3,x+(f0*b30+f2*b32)*h,f3)
    call me%f(t+h*a4,x+(f0*b40+f2*b42+f3*b43)*h,f4)
    call me%f(t+h*a5,x+(f0*b50+f3*b53+f4*b54)*h,f5)
    call me%f(t+h*a6,x+(f0*b60+f3*b63+f4*b64+f5*b65)*h,f6)
    call me%f(t+h*a7,x+(f0*b70+f4*b74+f5*b75+f6*b76)*h,f7)
    call me%f(t+h*a8,x+(f0*b80+f5*b85+f6*b86+f7*b87)*h,f8)
    call me%f(t+h*a9,x+(f0*b90+f5*b95+f6*b96+f7*b97+f8*b98)*h,f9)
    call me%f(t+h*a10,x+(f0*b100+f6*b106+f7*b107+f8*b108+&
            f9*b109)*h,f10)
    call me%f(t+h*a11,x+(f0*b110+f7*b117+f8*b118+f9*b119+&
            f10*b1110)*h,f11)
    call me%f(t+h*a12,x+(f0*b120+f5*b125+f6*b126+f7*b127+&
            f8*b128+f9*b129+f10*b1210+f11*b1211)*h,f12)
    call me%f(t+h*a13,x+(f0*b130+f5*b135+f6*b136+f7*b137+&
            f8*b138+f9*b139+f10*b1310+f11*b1311+f12*b1312)*h,f13)
    call me%f(t+h*a14,x+(f0*b140+f5*b145+f6*b146+f7*b147+f8*b148+&
            f9*b149+f10*b1410+f11*b1411+f12*b1412+f13*b1413)*h,f14)
    call me%f(t,x+(f0*b150+f8*b158+f9*b159+f10*b1510+f11*b1511+&
            f12*b1512+f13*b1513)*h,f15)
    call me%f(t+h*a16,x+(f0*b160+f5*b165+f6*b166+f7*b167+f8*b168+&
            f9*b169+f10*b1610+f11*b1611+f12*b1612+f13*b1613+&
            f15*b1615)*h,f16)

    xf = x+h*(f0*c0+f8*c8+f9*c9+f10*c10+f11*c11+f12*c12+f13*c13+f14*c14)

    terr = c14*h*(f0+f14-f15-f16)

    end subroutine rkf89
!*****************************************************************************************

!*****************************************************************************************
!>
!  Runge Kutta Verner 8(9)
!
!### Reference
!  * J. H. Verner, "Explicit Runge–Kutta Methods with Estimates of the
!    Local Truncation Error", SIAM Journal on Numerical Analysis,
!   15(4), 772–790, 1978.

    subroutine rkv89(me,t,x,h,xf,terr)

    implicit none

    class(rkv89_class),intent(inout)     :: me
    real(wp),intent(in)                  :: t     !! initial time
    real(wp),dimension(me%n),intent(in)  :: x     !! initial state
    real(wp),intent(in)                  :: h     !! time step
    real(wp),dimension(me%n),intent(out) :: xf    !! state at time `t+h`
    real(wp),dimension(me%n),intent(out) :: terr  !! truncation error estimate

    real(wp),parameter :: s6 = sqrt(6.0_wp)

    real(wp),parameter :: a2  = 1.0_wp/12.0_wp
    real(wp),parameter :: a3  = 1.0_wp/9.0_wp
    real(wp),parameter :: a4  = 1.0_wp/6.0_wp
    real(wp),parameter :: a5  = 2.0_wp*(1.0_wp+s6)/15.0_wp
    real(wp),parameter :: a6  = (6.0_wp+s6)/15.0_wp
    real(wp),parameter :: a7  = (6.0_wp-s6)/15.0_wp
    real(wp),parameter :: a8  = 2.0_wp/3.0_wp
    real(wp),parameter :: a9  = 1.0_wp/2.0_wp
    real(wp),parameter :: a10 = 1.0_wp/3.0_wp
    real(wp),parameter :: a11 = 1.0_wp/4.0_wp
    real(wp),parameter :: a12 = 4.0_wp/3.0_wp
    real(wp),parameter :: a13 = 5.0_wp/6.0_wp
    real(wp),parameter :: a15 = 1.0_wp/6.0_wp

    real(wp),parameter :: b31   = 1.0_wp/27.0_wp
    real(wp),parameter :: b32   = 2.0_wp/27.0_wp
    real(wp),parameter :: b41   = 1.0_wp/24.0_wp
    real(wp),parameter :: b43   = 3.0_wp/24.0_wp
    real(wp),parameter :: b51   = (4.0_wp+94.0_wp*s6)/375.0_wp
    real(wp),parameter :: b53   = -(282.0_wp+252.0_wp*s6)/375.0_wp
    real(wp),parameter :: b54   = (328.0_wp+208.0_wp*s6)/375.0_wp
    real(wp),parameter :: b61   = (9.0_wp-s6)/150.0_wp
    real(wp),parameter :: b64   = (312.0_wp+32.0_wp*s6)/1425.0_wp
    real(wp),parameter :: b65   = (69.0_wp+29.0_wp*s6)/570.0_wp
    real(wp),parameter :: b71   = (927.0_wp-347.0_wp*s6)/1250.0_wp
    real(wp),parameter :: b74   = (-16248.0_wp+7328.0_wp*s6)/9375.0_wp
    real(wp),parameter :: b75   = (-489.0_wp+179.0_wp*s6)/3750.0_wp
    real(wp),parameter :: b76   = (14268.0_wp-5798.0_wp*s6)/9375.0_wp
    real(wp),parameter :: b81   = 4.0_wp/54.0_wp
    real(wp),parameter :: b86   = (16.0_wp-s6)/54.0_wp
    real(wp),parameter :: b87   = (16.0_wp+s6)/54.0_wp
    real(wp),parameter :: b91   = 38.0_wp/512.0_wp
    real(wp),parameter :: b96   = (118.0_wp-23.0_wp*s6)/512.0_wp
    real(wp),parameter :: b97   = (118.0_wp+23.0_wp*s6)/512.0_wp
    real(wp),parameter :: b98   = -18.0_wp/512.0_wp
    real(wp),parameter :: b101  = 11.0_wp/144.0_wp
    real(wp),parameter :: b106  = (266.0_wp-s6)/864.0_wp
    real(wp),parameter :: b107  = (266.0_wp+s6)/864.0_wp
    real(wp),parameter :: b108  = -1.0_wp/16.0_wp
    real(wp),parameter :: b109  = -8.0_wp/27.0_wp
    real(wp),parameter :: b111  = (5034.0_wp-271.0_wp*s6)/61440.0_wp
    real(wp),parameter :: b117  = (7859.0_wp-1626.0_wp*s6)/10240.0_wp
    real(wp),parameter :: b118  = (-2232.0_wp+813.0_wp*s6)/20480.0_wp
    real(wp),parameter :: b119  = (-594.0_wp+271.0_wp*s6)/960.0_wp
    real(wp),parameter :: b1110 = (657.0_wp-813.0_wp*s6)/5120.0_wp
    real(wp),parameter :: b121  = (5996.0_wp-3794.0_wp*s6)/405.0_wp
    real(wp),parameter :: b126  = (-4342.0_wp-338.0_wp*s6)/9.0_wp
    real(wp),parameter :: b127  = (154922.0_wp-40458.0_wp*s6)/135.0_wp
    real(wp),parameter :: b128  = (-4176.0_wp+3794.0_wp*s6)/45.0_wp
    real(wp),parameter :: b129  = (-340864.0_wp+242816.0_wp*s6)/405.0_wp
    real(wp),parameter :: b1210 = (26304.0_wp-15176.0_wp*s6)/45.0_wp
    real(wp),parameter :: b1211 = -26624.0_wp/81.0_wp
    real(wp),parameter :: b131  = (3793.0_wp+2168.0_wp*s6)/103680.0_wp
    real(wp),parameter :: b136  = (4042.0_wp+2263.0_wp*s6)/13824.0_wp
    real(wp),parameter :: b137  = (-231278.0_wp+40717.0_wp*s6)/69120.0_wp
    real(wp),parameter :: b138  = (7947.0_wp-2168.0_wp*s6)/11520.0_wp
    real(wp),parameter :: b139  = (1048.0_wp-542.0_wp*s6)/405.0_wp
    real(wp),parameter :: b1310 = (-1383.0_wp+542.0_wp*s6)/720.0_wp
    real(wp),parameter :: b1311 = 2624.0_wp/1053.0_wp
    real(wp),parameter :: b1312 = 3.0_wp/1664.0_wp
    real(wp),parameter :: b141  = -137.0_wp/1296.0_wp
    real(wp),parameter :: b146  = (5642.0_wp-337.0_wp*s6)/864.0_wp
    real(wp),parameter :: b147  = (5642.0_wp+337.0_wp*s6)/864.0_wp
    real(wp),parameter :: b148  = -299.0_wp/48.0_wp
    real(wp),parameter :: b149  = 184.0_wp/81.0_wp
    real(wp),parameter :: b1410 = -44.0_wp/9.0_wp
    real(wp),parameter :: b1411 = -5120.0_wp/1053.0_wp
    real(wp),parameter :: b1412 = -11.0_wp/468.0_wp
    real(wp),parameter :: b1413 = 16.0_wp/9.0_wp
    real(wp),parameter :: b151  = (33617.0_wp-2168.0_wp*s6)/518400.0_wp
    real(wp),parameter :: b156  = (-3846.0_wp+31.0_wp*s6)/13824.0_wp
    real(wp),parameter :: b157  = (155338.0_wp-52807.0_wp*s6)/345600.0_wp
    real(wp),parameter :: b158  = (-12537.0_wp+2168.0_wp*s6)/57600.0_wp
    real(wp),parameter :: b159  = (92.0_wp+542.0_wp*s6)/2025.0_wp
    real(wp),parameter :: b1510 = (-1797.0_wp-542.0_wp*s6)/3600.0_wp
    real(wp),parameter :: b1511 = 320.0_wp/567.0_wp
    real(wp),parameter :: b1512 = -1.0_wp/1920.0_wp
    real(wp),parameter :: b1513 = 4.0_wp/105.0_wp
    real(wp),parameter :: b161  = (-36487.0_wp-30352.0_wp*s6)/279600.0_wp
    real(wp),parameter :: b166  = (-29666.0_wp-4499.0_wp*s6)/7456.0_wp
    real(wp),parameter :: b167  = (2779182.0_wp-615973.0_wp*s6)/186400.0_wp
    real(wp),parameter :: b168  = (-94329.0_wp+91056.0_wp*s6)/93200.0_wp
    real(wp),parameter :: b169  = (-232192.0_wp+121408.0_wp*s6)/17475.0_wp
    real(wp),parameter :: b1610 = (101226.0_wp-22764.0_wp*s6)/5825.0_wp
    real(wp),parameter :: b1611 = -169984.0_wp/9087.0_wp
    real(wp),parameter :: b1612 = -87.0_wp/30290.0_wp
    real(wp),parameter :: b1613 = 492.0_wp/1165.0_wp
    real(wp),parameter :: b1615 = 1260.0_wp/233.0_wp

    real(wp),parameter :: c1  = 103.0_wp/1680.0_wp
    real(wp),parameter :: c8  = -27.0_wp/140.0_wp
    real(wp),parameter :: c9  = 76.0_wp/105.0_wp
    real(wp),parameter :: c10 = -201.0_wp/280.0_wp
    real(wp),parameter :: c11 = 1024.0_wp/1365.0_wp
    real(wp),parameter :: c12 = 3.0_wp/7280.0_wp
    real(wp),parameter :: c13 = 12.0_wp/35.0_wp
    real(wp),parameter :: c14 = 9.0_wp/280.0_wp

    real(wp),parameter :: e1  = -1911.0_wp/109200.0_wp
    real(wp),parameter :: e8  = 34398.0_wp/109200.0_wp
    real(wp),parameter :: e9  = -61152.0_wp/109200.0_wp
    real(wp),parameter :: e10 = 114660.0_wp/109200.0_wp
    real(wp),parameter :: e11 = -114688.0_wp/109200.0_wp
    real(wp),parameter :: e12 = -63.0_wp/109200.0_wp
    real(wp),parameter :: e13 = -13104.0_wp/109200.0_wp
    real(wp),parameter :: e14 = -3510.0_wp/109200.0_wp
    real(wp),parameter :: e15 = 39312.0_wp/109200.0_wp
    real(wp),parameter :: e16 = 6058.0_wp/109200.0_wp

    real(wp),dimension(me%n) :: f1,f2,f3,f4,f5,f6,f7,f8,f9,&
                                f10,f11,f12,f13,f14,f15,f16

    call me%f(t,x,f1)
    call me%f(t+a2*h,x+h*(a2*f1),f2)
    call me%f(t+a3*h,x+h*(b31*f1+b32*f2),f3)
    call me%f(t+a4*h,x+h*(b41*f1+b43*f3),f4)
    call me%f(t+a5*h,x+h*(b51*f1+b53*f3+b54*f4),f5)
    call me%f(t+a6*h,x+h*(b61*f1+b64*f4+b65*f5),f6)
    call me%f(t+a7*h,x+h*(b71*f1+b74*f4+b75*f5+b76*f6),f7)
    call me%f(t+a8*h,x+h*(b81*f1+b86*f6+b87*f7),f8)
    call me%f(t+a9*h,x+h*(b91*f1+b96*f6+b97*f7+b98*f8),f9)
    call me%f(t+a10*h,x+h*(b101*f1+b106*f6+b107*f7+b108*f8+b109*f9),f10)
    call me%f(t+a11*h,x+h*(b111*f1+b117*f7+b118*f8+b119*f9+b1110*f10),f11)
    call me%f(t+a12*h,x+h*(b121*f1+b126*f6+b127*f7+b128*f8+b129*f9+&
                b1210*f10+b1211*f11),f12)
    call me%f(t+a13*h,x+h*(b131*f1+b136*f6+b137*f7+b138*f8+b139*f9+&
                b1310*f10+b1311*f11+b1312*f12),f13)
    call me%f(t+h,x+h*(b141*f1+b146*f6+b147*f7+b148*f8+b149*f9+b1410*f10+&
                b1411*f11+b1412*f12+b1413*f13),f14)
    call me%f(t+a15*h,x+h*(b151*f1+b156*f6+b157*f7+b158*f8+b159*f9+b1510*f10+&
                b1511*f11+b1512*f12+b1513*f13),f15)
    call me%f(t+h,x+h*(b161*f1+b166*f6+b167*f7+b168*f8+b169*f9+b1610*f10+&
                b1611*f11+b1612*f12+b1613*f13+b1615*f15),f16)

    xf = x+h*(c1*f1+c8*f8+c9*f9+c10*f10+c11*f11+c12*f12+c13*f13+c14*f14)

    terr = e1*f1+e8*f8+e9*f9+e10*f10+e11*f11+e12*f12+e13*f13+e14*f14+e15*f15+e16*f16

    end subroutine rkv89
!*****************************************************************************************

!*****************************************************************************************
!>
!  Feagin's RK8(10) method.
!
!### Reference
!  * T. Feagin, "[A Tenth-Order Runge-Kutta Method with Error Estimate]
!    (http://sce.uhcl.edu/feagin/courses/rk10.pdf)",

    subroutine rkf108(me,t,x,h,xf,terr)

    implicit none

    class(rkf108_class),intent(inout)    :: me
    real(wp),intent(in)                  :: t     !! initial time
    real(wp),dimension(me%n),intent(in)  :: x     !! initial state
    real(wp),intent(in)                  :: h     !! time step
    real(wp),dimension(me%n),intent(out) :: xf    !! state at time `t+h`
    real(wp),dimension(me%n),intent(out) :: terr  !! truncation error estimate

    real(wp),parameter :: a1  = 0.100000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: a2  = 0.539357840802981787532485197881302436857273449701009015505500_wp
    real(wp),parameter :: a3  = 0.809036761204472681298727796821953655285910174551513523258250_wp
    real(wp),parameter :: a4  = 0.309036761204472681298727796821953655285910174551513523258250_wp
    real(wp),parameter :: a5  = 0.981074190219795268254879548310562080489056746118724882027805_wp
    real(wp),parameter :: a6  = 0.833333333333333333333333333333333333333333333333333333333333_wp
    real(wp),parameter :: a7  = 0.354017365856802376329264185948796742115824053807373968324184_wp
    real(wp),parameter :: a8  = 0.882527661964732346425501486979669075182867844268052119663791_wp
    real(wp),parameter :: a9  = 0.642615758240322548157075497020439535959501736363212695909875_wp
    real(wp),parameter :: a10 = 0.357384241759677451842924502979560464040498263636787304090125_wp
    real(wp),parameter :: a11 = 0.117472338035267653574498513020330924817132155731947880336209_wp
    real(wp),parameter :: a12 = 0.833333333333333333333333333333333333333333333333333333333333_wp
    real(wp),parameter :: a13 = 0.309036761204472681298727796821953655285910174551513523258250_wp
    real(wp),parameter :: a14 = 0.539357840802981787532485197881302436857273449701009015505500_wp
    real(wp),parameter :: a15 = 0.100000000000000000000000000000000000000000000000000000000000_wp
    !real(wp),parameter :: a16 = 1.00000000000000000000000000000000000000000000000000000000000_wp

    real(wp),parameter :: c0  = 0.0333333333333333333333333333333333333333333333333333333333333_wp
    real(wp),parameter :: c1  = 0.0250000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c2  = 0.0333333333333333333333333333333333333333333333333333333333333_wp
    real(wp),parameter :: c4  = 0.0500000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c6  = 0.0400000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c8  = 0.189237478148923490158306404106012326238162346948625830327194_wp
    real(wp),parameter :: c9  = 0.277429188517743176508360262560654340428504319718040836339472_wp
    real(wp),parameter :: c10 = 0.277429188517743176508360262560654340428504319718040836339472_wp
    real(wp),parameter :: c11 = 0.189237478148923490158306404106012326238162346948625830327194_wp
    real(wp),parameter :: c12 = -0.0400000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c13 = -0.0500000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c14 = -0.0333333333333333333333333333333333333333333333333333333333333_wp
    real(wp),parameter :: c15 = -0.0250000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c16 = 0.0333333333333333333333333333333333333333333333333333333333333_wp

    real(wp),parameter :: b10 = 0.100000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: b20 = -0.915176561375291440520015019275342154318951387664369720564660_wp
    real(wp),parameter :: b21 = 1.45453440217827322805250021715664459117622483736537873607016_wp
    real(wp),parameter :: b30 = 0.202259190301118170324681949205488413821477543637878380814562_wp
    real(wp),parameter :: b32 = 0.606777570903354510974045847616465241464432630913635142443687_wp
    real(wp),parameter :: b40 = 0.184024714708643575149100693471120664216774047979591417844635_wp
    real(wp),parameter :: b42 = 0.197966831227192369068141770510388793370637287463360401555746_wp
    real(wp),parameter :: b43 = -0.0729547847313632629185146671595558023015011608914382961421311_wp
    real(wp),parameter :: b50 = 0.0879007340206681337319777094132125475918886824944548534041378_wp
    real(wp),parameter :: b53 = 0.410459702520260645318174895920453426088035325902848695210406_wp
    real(wp),parameter :: b54 = 0.482713753678866489204726942976896106809132737721421333413261_wp
    real(wp),parameter :: b60 = 0.0859700504902460302188480225945808401411132615636600222593880_wp
    real(wp),parameter :: b63 = 0.330885963040722183948884057658753173648240154838402033448632_wp
    real(wp),parameter :: b64 = 0.489662957309450192844507011135898201178015478433790097210790_wp
    real(wp),parameter :: b65 = -0.0731856375070850736789057580558988816340355615025188195854775_wp
    real(wp),parameter :: b70 = 0.120930449125333720660378854927668953958938996999703678812621_wp
    real(wp),parameter :: b74 = 0.260124675758295622809007617838335174368108756484693361887839_wp
    real(wp),parameter :: b75 = 0.0325402621549091330158899334391231259332716675992700000776101_wp
    real(wp),parameter :: b76 = -0.0595780211817361001560122202563305121444953672762930724538856_wp
    real(wp),parameter :: b80 = 0.110854379580391483508936171010218441909425780168656559807038_wp
    real(wp),parameter :: b85 = -0.0605761488255005587620924953655516875526344415354339234619466_wp
    real(wp),parameter :: b86 = 0.321763705601778390100898799049878904081404368603077129251110_wp
    real(wp),parameter :: b87 = 0.510485725608063031577759012285123416744672137031752354067590_wp
    real(wp),parameter :: b90 = 0.112054414752879004829715002761802363003717611158172229329393_wp
    real(wp),parameter :: b95 = -0.144942775902865915672349828340980777181668499748506838876185_wp
    real(wp),parameter :: b96 = -0.333269719096256706589705211415746871709467423992115497968724_wp
    real(wp),parameter :: b97 = 0.499269229556880061353316843969978567860276816592673201240332_wp
    real(wp),parameter :: b98 = 0.509504608929686104236098690045386253986643232352989602185060_wp
    real(wp),parameter :: b100 = 0.113976783964185986138004186736901163890724752541486831640341_wp
    real(wp),parameter :: b105 = -0.0768813364203356938586214289120895270821349023390922987406384_wp
    real(wp),parameter :: b106 = 0.239527360324390649107711455271882373019741311201004119339563_wp
    real(wp),parameter :: b107 = 0.397774662368094639047830462488952104564716416343454639902613_wp
    real(wp),parameter :: b108 = 0.0107558956873607455550609147441477450257136782823280838547024_wp
    real(wp),parameter :: b109 = -0.327769124164018874147061087350233395378262992392394071906457_wp
    real(wp),parameter :: b110 = 0.0798314528280196046351426864486400322758737630423413945356284_wp
    real(wp),parameter :: b115 = -0.0520329686800603076514949887612959068721311443881683526937298_wp
    real(wp),parameter :: b116 = -0.0576954146168548881732784355283433509066159287152968723021864_wp
    real(wp),parameter :: b117 = 0.194781915712104164976306262147382871156142921354409364738090_wp
    real(wp),parameter :: b118 = 0.145384923188325069727524825977071194859203467568236523866582_wp
    real(wp),parameter :: b119 = -0.0782942710351670777553986729725692447252077047239160551335016_wp
    real(wp),parameter :: b1110 = -0.114503299361098912184303164290554670970133218405658122674674_wp
    real(wp),parameter :: b120 = 0.985115610164857280120041500306517278413646677314195559520529_wp
    real(wp),parameter :: b123 = 0.330885963040722183948884057658753173648240154838402033448632_wp
    real(wp),parameter :: b124 = 0.489662957309450192844507011135898201178015478433790097210790_wp
    real(wp),parameter :: b125 = -1.37896486574843567582112720930751902353904327148559471526397_wp
    real(wp),parameter :: b126 = -0.861164195027635666673916999665534573351026060987427093314412_wp
    real(wp),parameter :: b127 = 5.78428813637537220022999785486578436006872789689499172601856_wp
    real(wp),parameter :: b128 = 3.28807761985103566890460615937314805477268252903342356581925_wp
    real(wp),parameter :: b129 = -2.38633905093136384013422325215527866148401465975954104585807_wp
    real(wp),parameter :: b1210 = -3.25479342483643918654589367587788726747711504674780680269911_wp
    real(wp),parameter :: b1211 = -2.16343541686422982353954211300054820889678036420109999154887_wp
    real(wp),parameter :: b130 = 0.895080295771632891049613132336585138148156279241561345991710_wp
    real(wp),parameter :: b132 = 0.197966831227192369068141770510388793370637287463360401555746_wp
    real(wp),parameter :: b133 = -0.0729547847313632629185146671595558023015011608914382961421311_wp
    real(wp),parameter :: b135 = -0.851236239662007619739049371445966793289359722875702227166105_wp
    real(wp),parameter :: b136 = 0.398320112318533301719718614174373643336480918103773904231856_wp
    real(wp),parameter :: b137 = 3.63937263181035606029412920047090044132027387893977804176229_wp
    real(wp),parameter :: b138 = 1.54822877039830322365301663075174564919981736348973496313065_wp
    real(wp),parameter :: b139 = -2.12221714704053716026062427460427261025318461146260124401561_wp
    real(wp),parameter :: b1310 = -1.58350398545326172713384349625753212757269188934434237975291_wp
    real(wp),parameter :: b1311 = -1.71561608285936264922031819751349098912615880827551992973034_wp
    real(wp),parameter :: b1312 = -0.0244036405750127452135415444412216875465593598370910566069132_wp
    real(wp),parameter :: b140 = -0.915176561375291440520015019275342154318951387664369720564660_wp
    real(wp),parameter :: b141 = 1.45453440217827322805250021715664459117622483736537873607016_wp
    real(wp),parameter :: b144 = -0.777333643644968233538931228575302137803351053629547286334469_wp
    real(wp),parameter :: b146 = -0.0910895662155176069593203555807484200111889091770101799647985_wp
    real(wp),parameter :: b1412 = 0.0910895662155176069593203555807484200111889091770101799647985_wp
    real(wp),parameter :: b1413 = 0.777333643644968233538931228575302137803351053629547286334469_wp
    real(wp),parameter :: b150 = 0.100000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: b152 = -0.157178665799771163367058998273128921867183754126709419409654_wp
    real(wp),parameter :: b1514 = 0.157178665799771163367058998273128921867183754126709419409654_wp
    real(wp),parameter :: b160 = 0.181781300700095283888472062582262379650443831463199521664945_wp
    real(wp),parameter :: b161 = 0.675000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: b162 = 0.342758159847189839942220553413850871742338734703958919937260_wp
    real(wp),parameter :: b164 = 0.259111214548322744512977076191767379267783684543182428778156_wp
    real(wp),parameter :: b165 = -0.358278966717952089048961276721979397739750634673268802484271_wp
    real(wp),parameter :: b166 = -1.04594895940883306095050068756409905131588123172378489286080_wp
    real(wp),parameter :: b167 = 0.930327845415626983292300564432428777137601651182965794680397_wp
    real(wp),parameter :: b168 = 1.77950959431708102446142106794824453926275743243327790536000_wp
    real(wp),parameter :: b169 = 0.100000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: b1610 = -0.282547569539044081612477785222287276408489375976211189952877_wp
    real(wp),parameter :: b1611 = -0.159327350119972549169261984373485859278031542127551931461821_wp
    real(wp),parameter :: b1612 = -0.145515894647001510860991961081084111308650130578626404945571_wp
    real(wp),parameter :: b1613 = -0.259111214548322744512977076191767379267783684543182428778156_wp
    real(wp),parameter :: b1614 = -0.342758159847189839942220553413850871742338734703958919937260_wp
    real(wp),parameter :: b1615 = -0.675000000000000000000000000000000000000000000000000000000000_wp

    real(wp),dimension(me%n) :: f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16

    if (h==zero) then
        xf = x
        terr = zero
        return
    end if

    call me%f(t,x,f0)
    call me%f(t+a1*h,x+h*(b10*f0),f1)
    call me%f(t+a2*h,x+h*(b20*f0+b21*f1),f2)
    call me%f(t+a3*h,x+h*(b30*f0+b32*f2),f3)
    call me%f(t+a4*h,x+h*(b40*f0+b42*f2+b43*f3),f4)
    call me%f(t+a5*h,x+h*(b50*f0+b53*f3+b54*f4),f5)
    call me%f(t+a6*h,x+h*(b60*f0+b63*f3+b64*f4+b65*f5),f6)
    call me%f(t+a7*h,x+h*(b70*f0+b74*f4+b75*f5+b76*f6),f7)
    call me%f(t+a8*h,x+h*(b80*f0+b85*f5+b86*f6+b87*f7),f8)
    call me%f(t+a9*h,x+h*(b90*f0+b95*f5+b96*f6+b97*f7+b98*f8),f9)
    call me%f(t+a10*h,x+h*(b100*f0+b105*f5+b106*f6+b107*f7+b108*f8+&
                b109*f9),f10)
    call me%f(t+a11*h,x+h*(b110*f0+b115*f5+b116*f6+b117*f7+b118*f8+&
                b119*f9+b1110*f10),f11)
    call me%f(t+a12*h,x+h*(b120*f0+b123*f3+b124*f4+b125*f5+b126*f6+&
                b127*f7+b128*f8+b129*f9+b1210*f10+b1211*f11),f12)
    call me%f(t+a13*h,x+h*(b130*f0+b132*f2+b133*f3+b135*f5+b136*f6+&
                b137*f7+b138*f8+b139*f9+b1310*f10+b1311*f11+&
                b1312*f12),f13)
    call me%f(t+a14*h,x+h*(b140*f0+b141*f1+b144*f4+b146*f6+b1412*f12+&
                b1413*f13),f14)
    call me%f(t+a15*h,x+h*(b150*f0+b152*f2+b1514*f14),f15)
    call me%f(t+h,x+h*(b160*f0+b161*f1+b162*f2+b164*f4+b165*f5+&
                b166*f6+b167*f7+b168*f8+b169*f9+b1610*f10+b1611*f11+&
                b1612*f12+b1613*f13+b1614*f14+b1615*f15),f16)

    xf = x+h*(c0*f0+c1*f1+c2*f2+c4*f4+c6*f6+c8*f8+c9*f9+&
              c10*f10+c11*f11+c12*f12+c13*f13+c14*f14+c15*f15+c16*f16)

    terr = (1.0_wp/360.0_wp)*h*(f1-f15)

    end subroutine rkf108
!*****************************************************************************************

!*****************************************************************************************
!>
!  Returns the order of the rkf78 method.

    pure function rkf78_order(me) result(p)

    implicit none

    class(rkf78_class),intent(in) :: me
    integer                       :: p    !! order of the method

    p = 7

    end function rkf78_order
!*****************************************************************************************

!*****************************************************************************************
!>
!  Returns the order of the rkf89 method.

    pure function rkf89_order(me) result(p)

    implicit none

    class(rkf89_class),intent(in) :: me
    integer                       :: p    !! order of the method

    p = 8

    end function rkf89_order
!*****************************************************************************************

!*****************************************************************************************
!>
!  Returns the order of the rkv89 method.

    pure function rkv89_order(me) result(p)

    implicit none

    class(rkv89_class),intent(in) :: me
    integer                       :: p    !! order of the method

    p = 8

    end function rkv89_order
!*****************************************************************************************

!*****************************************************************************************
!>
!  Returns the order of the [[rkf108]] method.

    pure function rkf108_order(me) result(p)

    implicit none

    class(rkf108_class),intent(in) :: me
    integer                        :: p    !! order of the method

    p = 10

    end function rkf108_order
!*****************************************************************************************

!*****************************************************************************************
!>
!  Computes a starting step size to be used in solving initial
!  value problems in ordinary differential equations.
!
!  It is based on an estimate of the local lipschitz constant for the
!  differential equation (lower bound on a norm of the jacobian) ,
!  a bound on the differential equation (first derivative), and
!  a bound on the partial derivative of the equation with respect to
!  the independent variable. (all approximated near the initial point a)
!
!@note Subroutine hstart also uses the `me%stepsize_method%norm`
!      function for computing vector norms
!
!@note This routine is from [DDEABM](https://github.com/jacobwilliams/ddeabm).
!
!# History
!   * 820301  date written -- watts, h. a., (snla)
!   * 890531  changed all specific intrinsics to generic.  (wrb)
!   * 890831  modified array declarations.  (wrb)
!   * 890911  removed unnecessary intrinsics.  (wrb)
!   * 891024  changed references from dvnorm to dhvnrm.  (wrb)
!   * 891214  prologue converted to version 4.0 format.  (bab)
!   * 900328  added type section.  (wrb)
!   * 910722  updated author section.  (als)
!   * December, 2015 : Refactored this routine (jw)
!   * April 2016 : Some modifications for the variable-step RK module (jw)

    subroutine hstart(me,a,b,y,yprime,etol,h)

    implicit none

    class(rk_variable_step_class),intent(inout) :: me
    real(wp),intent(in)                :: a         !! the initial point of integration.
    real(wp),intent(in)                :: b         !! a value of the independent variable used to define
                                                    !! the direction of integration. a reasonable choice is to
                                                    !! set `b` to the first point at which a solution is desired.
                                                    !! you can also use `b`, if necessary, to restrict the length
                                                    !! of the first integration step because the algorithm will
                                                    !! not compute a starting step length which is bigger than
                                                    !! `abs(b-a)`, unless `b` has been chosen too close to `a`.
                                                    !! (it is presumed that hstart has been called with `b`
                                                    !! different from `a` on the machine being used. also see the
                                                    !! discussion about the parameter `small`.)
    real(wp),dimension(me%n),intent(in) :: y        !! the vector of initial values of the `neq` solution
                                                    !! components at the initial point `a`.
    real(wp),dimension(me%n),intent(in) :: yprime   !! the vector of derivatives of the `neq`
                                                    !! solution components at the initial point `a`.
                                                    !! (defined by the differential equations in subroutine `me%f`)
    real(wp),dimension(me%n),intent(in) :: etol     !! the vector of error tolerances corresponding to
                                                    !! the `neq` solution components. it is assumed that all
                                                    !! elements are positive. following the first integration
                                                    !! step, the tolerances are expected to be used by the
                                                    !! integrator in an error test which roughly requires that
                                                    !! `abs(local error) <= etol` for each vector component.
    real(wp),intent(out) :: h       !! appropriate starting step size to be attempted by the
                                    !! differential equation method.

    real(wp),dimension(me%n) :: spy  !! work array which provide the routine with needed storage space.
    real(wp),dimension(me%n) :: pv   !! work array which provide the routine with needed storage space.
    real(wp),dimension(me%n) :: yp   !! work array which provide the routine with needed storage space.
    real(wp),dimension(me%n) :: sf   !! work array which provide the routine with needed storage space.

    real(wp),parameter :: small  = epsilon(one)
    real(wp),parameter :: big    = huge(one)
    real(wp),parameter :: relper = small**0.375_wp

    integer :: j, k, lk
    real(wp) :: absdx, da, delf, dely,&
                dfdub, dfdxb,&
                dx, dy, fbnd,&
                srydpb, tolexp, tolmin, tolp, tolsum, ydpb
    integer :: morder    !! the order of the formula which will be used by
                         !! the initial value method for taking the first integration
                         !! step.

    morder = me%p
    dx = b - a
    absdx = abs(dx)

    ! compute an approximate bound (dfdxb) on the partial
    ! derivative of the equation with respect to the
    ! independent variable. protect against an overflow.
    ! also compute a bound (fbnd) on the first derivative
    ! locally.

    da = sign(max(min(relper*abs(a),absdx),100.0_wp*small*abs(a)),dx)
    if (da == zero) da = relper*dx
    call me%f(a+da,y,sf)
    yp = sf - yprime
    delf = me%stepsize_method%norm(yp)
    dfdxb = big
    if (delf < big*abs(da)) dfdxb = delf/abs(da)
    fbnd = me%stepsize_method%norm(sf)

    ! compute an estimate (dfdub) of the local lipschitz
    ! constant for the system of differential equations. this
    ! also represents an estimate of the norm of the jacobian
    ! locally.  three iterations (two when neq=1) are used to
    ! estimate the lipschitz constant by numerical differences.
    ! the first perturbation vector is based on the initial
    ! derivatives and direction of integration. the second
    ! perturbation vector is formed using another evaluation of
    ! the differential equation.  the third perturbation vector
    ! is formed using perturbations based only on the initial
    ! values. components that are zero are always changed to
    ! non-zero values (except on the first iteration). when
    ! information is available, care is taken to ensure that
    ! components of the perturbation vector have signs which are
    ! consistent with the slopes of local solution curves.
    ! also choose the largest bound (fbnd) for the first
    ! derivative.
    !
    ! perturbation vector size is held
    ! constant for all iterations. compute
    ! this change from the
    ! size of the vector of initial
    ! values.

    dely = relper*me%stepsize_method%norm(y)
    if (dely == zero) dely = relper
    dely = sign(dely,dx)
    delf = me%stepsize_method%norm(yprime)
    fbnd = max(fbnd,delf)
    if (delf == zero) then
        ! cannot have a null perturbation vector
        spy  = zero
        yp   = one
        delf = me%stepsize_method%norm(yp)
    else
        ! use initial derivatives for first perturbation
        spy = yprime
        yp  = yprime
    end if

    dfdub = zero
    lk = min(me%n+1,3)
    do k = 1, lk
        ! define perturbed vector of initial values
        pv = y + yp * (dely/delf)
        if (k == 2) then
            ! use a shifted value of the independent variable
            ! in computing one estimate
            call me%f(a+da,pv,yp)
            pv = yp - sf
        else
            ! evaluate derivatives associated with perturbed
            ! vector and compute corresponding differences
            call me%f(a,pv,yp)
            pv = yp - yprime
        end if
        ! choose largest bounds on the first derivative
        ! and a local lipschitz constant
        fbnd = max(fbnd,me%stepsize_method%norm(yp))
        delf = me%stepsize_method%norm(pv)
        if (delf >= big*abs(dely)) then
            ! protect against an overflow
            dfdub = big
            exit
        end if
        dfdub = max(dfdub,delf/abs(dely))
        if (k == lk) exit

        ! choose next perturbation vector
        if (delf == zero) delf = one
        do j = 1, me%n
            if (k == 2) then
                dy = y(j)
                if (dy == zero) dy = dely/relper
            else
                dy = abs(pv(j))
                if (dy == zero) dy = delf
            end if
            if (spy(j) == zero) spy(j) = yp(j)
            if (spy(j) /= zero) dy = sign(dy,spy(j))
            yp(j) = dy
        end do
        delf = me%stepsize_method%norm(yp)
    end do

    ! compute a bound (ydpb) on the norm of the second derivative

    ydpb = dfdxb + dfdub*fbnd

    ! define the tolerance parameter upon which the starting step
    ! size is to be based.  a value in the middle of the error
    ! tolerance range is selected.

    tolmin = big
    tolsum = zero
    do k = 1, me%n
        tolexp = log10(etol(k))
        tolmin = min(tolmin,tolexp)
        tolsum = tolsum + tolexp
    end do
    tolp = 10.0_wp**(0.5_wp*(tolsum/me%n + tolmin)/(morder+1))

    ! compute a starting step size based on the above first and
    ! second derivative information
    !
    ! restrict the step length to be not bigger
    ! than abs(b-a). (unless b is too close to a)

    h = absdx

    if (ydpb == zero .and. fbnd == zero) then
        ! both first derivative term (fbnd) and second
        ! derivative term (ydpb) are zero
        if (tolp < one) h = absdx*tolp
    elseif (ydpb == zero) then
        ! only second derivative term (ydpb) is zero
        if (tolp < fbnd*absdx) h = tolp/fbnd
    else
        ! second derivative term (ydpb) is non-zero
        srydpb = sqrt(0.5_wp*ydpb)
        if (tolp < srydpb*absdx) h = tolp/srydpb
    end if

    ! further restrict the step length to be not bigger than  1/dfdub
    if (h*dfdub > one) h = one/dfdub

    ! finally, restrict the step length to be not
    ! smaller than 100*small*abs(a). however, if
    ! a=0. and the computed h underflowed to zero,
    ! the algorithm returns small*abs(b) for the step length.
    h = max(h,100.0_wp*small*abs(a))
    if (h == zero) h = small*abs(b)

    ! now set direction of integration
    h = sign(h,dx)

    end subroutine hstart
!*****************************************************************************************

!*****************************************************************************************
!>
!  computation of an initial step size guess
!
!@note This routine is from dop853. It was modified for this module.

    function hinit(me,x,y,posneg,f0,hmax,atol,rtol)

    implicit none

    class(rk_variable_step_class),intent(inout) :: me
    real(wp),intent(in)               :: x
    real(wp),dimension(:),intent(in)  :: y       !! dimension(n)
    real(wp),intent(in)               :: posneg  !! posneg = sign(1.0_wp,xend-x)
    real(wp),dimension(:),intent(in)  :: f0      !! dimension(n)
    real(wp),intent(in)               :: hmax
    real(wp),dimension(:),intent(in)  :: atol
    real(wp),dimension(:),intent(in)  :: rtol

    real(wp) :: der12,der2,dnf,dny,h,h1,hinit,sk
    integer :: i
    integer :: iord  !! order of the method
    real(wp),dimension(me%n) :: f1,y1

    iord = me%p

    ! compute a first guess for explicit euler as
    !   h = 0.01 * norm (y0) / norm (f0)
    ! the increment for explicit euler is small
    ! compared to the solution
    dnf = zero
    dny = zero
    do i = 1 , me%n
        sk = atol(i) + rtol(i)*abs(y(i))
        dnf = dnf + (f0(i)/sk)**2
        dny = dny + (y(i)/sk)**2
    end do
    if ( dnf<=1.0e-10_wp .or. dny<=1.0e-10_wp ) then
        h = 1.0e-6_wp
    else
        h = sqrt(dny/dnf)*0.01_wp
    end if
    h = min(h,hmax)
    h = sign(h,posneg)
    ! perform an explicit euler step
    do i = 1 , me%n
        y1(i) = y(i) + h*f0(i)
    end do
    call me%f(x+h,y1,f1)
    ! estimate the second derivative of the solution
    der2 = zero
    do i = 1 , me%n
        sk = atol(i) + rtol(i)*abs(y(i))
        der2 = der2 + ((f1(i)-f0(i))/sk)**2
    end do
    der2 = sqrt(der2)/h
    ! step size is computed such that
    !  h**iord * max ( norm (f0), norm (der2)) = 0.01
    der12 = max(abs(der2),sqrt(dnf))
    if ( der12<=1.0e-15_wp ) then
        h1 = max(1.0e-6_wp,abs(h)*1.0e-3_wp)
    else
        h1 = (0.01_wp/der12)**(1.0_wp/iord)
    end if

    h = min(100.0_wp*abs(h),h1,hmax)
    hinit = sign(h,posneg)

    end function hinit
!*****************************************************************************************

!*****************************************************************************************
!>
!  Unit tests for step size adjustment routines.

    subroutine step_size_test()

    implicit none

    type(stepsize_class)   :: s1     !! for testing the different methods
    type(stepsize_class)   :: s2     !! for testing the different methods
    type(stepsize_class)   :: s3     !! for testing the different methods
    real(wp)               :: h      !! current step size
    real(wp)               :: tol    !! abs error tolerance
    real(wp)               :: err    !! truncation error estimate
    integer                :: p      !! order of the method
    real(wp)               :: hnew   !! new step size
    logical                :: accept !! if the step is accepted

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' step_size_test'
    write(*,*) '---------------'
    write(*,*) ''

    h   = 10.0_wp
    tol = 1.0e-9_wp
    err = 1.0e-7_wp
    p   = 4

    call s1%initialize()
    call s1%compute_stepsize(h,tol,err,p,hnew,accept)
    write(*,*) 'stepsize_hull    : hnew = ', hnew

    call s2%initialize()
    call s2%compute_stepsize(h,tol,err,p,hnew,accept)
    write(*,*) 'stepsize_stoer_1 : hnew = ', hnew

    call s3%initialize()
    call s3%compute_stepsize(h,tol,err,p,hnew,accept)
    write(*,*) 'stepsize_stoer_2 : hnew = ', hnew

    end subroutine step_size_test
!*****************************************************************************************

!*****************************************************************************************
!>
!  Unit test of the [[rk_module]].
!  Integrate a two-body orbit around the Earth.

    subroutine rk_test_variable_step()

    implicit none

    type,extends(rkf78_class) :: spacecraft
    !type,extends(rkf89_class) :: spacecraft
    !type,extends(rkv89_class) :: spacecraft
    !type,extends(rkf108_class) :: spacecraft
        !! spacecraft propagation type.
        real(wp) :: mu     = zero     !! central body gravitational parameter (km3/s2)
        integer  :: fevals = 0        !! number of function evaluations
        logical  :: first  = .true.   !! first point is being exported
    end type spacecraft

    integer,parameter  :: n = 6            !! number of state variables
    real(wp),parameter :: tol = 1.0e-12_wp !! event location tolerance

    type(spacecraft) :: s, s2
    real(wp) :: t0,tf,x0(n),dt,xf(n),x02(n),gf,tf_actual,rtol,atol

    type(stepsize_class) :: sz

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' rk_variable_step_test'
    write(*,*) '---------------'
    write(*,*) ''

    !***************************************************************************

    !step size constructor:
    !call sz%initialize(hmin=1.0e-6_wp,hmax=1.0e6_wp)
    call sz%initialize( hmin           = 1.0e-6_wp,&
                        hmax           = 1.0e+6_wp,&
                        hfactor_reject = 0.5_wp,   &
                        hfactor_accept = 2.0_wp,   &
                        max_attempts   = 5,        &
                        accept_mode    = 1,        &
                        norm = maxval_func, &
                        compute_h_factor = compute_h_factor_v5 )

    !integrator constructor:
    call s%initialize(n=n,f=twobody,rtol=[1.0e-12_wp],atol=[1.0e-12_wp],&
                        stepsize_method=sz,report=twobody_report)

    !initial conditions:
    x0   = [10000.0_wp,10000.0_wp,10000.0_wp,&   !initial state [r,v] (km,km/s)
            1.0_wp,2.0_wp,3.0_wp]
    t0   = zero              ! initial time (sec)
    !dt   = 0.0_wp           ! automatically compute an initial time step (sec)
    dt   = 10.0_wp           ! initial time step (sec)
    tf   = 1000.0_wp         ! final time (sec)
    s%mu = 398600.436233_wp  ! main body is Earth

    s%fevals = 0
    s%first = .true.
    call s%integrate(t0,x0,dt,tf,xf)    !forward
    write(*,*) ''
    write(*,'(A/,*(F15.6/))') 'Final state:',xf
    write(*,'(A,I5)') 'Function evaluations:', s%fevals
    write(*,'(A,I5)') 'Number of rejected steps:',s%num_rejected_steps

    s%fevals = 0
    s%report => null()    !disable reporting
    call s%integrate(tf,xf,-dt,t0,x02)  !backwards
    s%num_rejected_steps = 0

    write(*,'(A/,*(E20.12/))') 'Error:',x02-x0
    write(*,'(A,I5)') 'Function evaluations:', s%fevals
    write(*,'(A,I5)') 'Number of rejected steps:',s%num_rejected_steps
    write(*,*) ''

   !***************************************************************************
   !event finding test:

   write(*,*) ' Event test - integrate until z = 12,000'

   ! NOTE: the following causes an ICE in gfortran 7.1, but works with ifort:
   ! s2 = spacecraft(n=n,f=twobody,g=twobody_event,mu=398600.436233_wp,&
   !                  rtol=[1.0e-12_wp],atol=[1.0e-12_wp],&
   !                  stepsize_method=sz,report=twobody_report)
   ! do it this way instead:
   call s2%initialize(n=n,f=twobody,g=twobody_event,&
                      rtol=[1.0e-12_wp],atol=[1.0e-12_wp],&
                      stepsize_method=sz,&
                      report=twobody_report)
   s2%mu = 398600.436233_wp

   s2%fevals = 0
   s2%first = .true.
   x0 = [10000.0_wp,10000.0_wp,10000.0_wp,&   !initial state [r,v] (km,km/s)
           1.0_wp,2.0_wp,3.0_wp]
   t0 = zero       !initial time (sec)
   dt = 10.0_wp    !time step (sec)
   tf = 1000.0_wp  !final time (sec)

   call s2%integrate_to_event(t0,x0,dt,tf,tol,tf_actual,xf,gf)
   write(*,*) ''
   write(*,'(A/,*(F15.6/))') 'Final time: ',tf_actual
   write(*,'(A/,*(F15.6/))') 'Final state:',xf
   write(*,'(A/,*(F15.6/))') 'Event func :',gf
   write(*,'(A,I5)') 'Function evaluations:', s2%fevals
   write(*,'(A,I5)') 'Number of rejected steps:',s2%num_rejected_steps

    contains
!*****************************************************************************************

    !*********************************************************
        subroutine twobody(me,t,x,xdot)

        !! derivative routine for two-body orbit propagation

        implicit none

        class(rk_variable_step_class),intent(inout)     :: me
        real(wp),intent(in)               :: t
        real(wp),dimension(:),intent(in)  :: x
        real(wp),dimension(:),intent(out) :: xdot

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

        !! report function - write time,state to console

        implicit none

        class(rk_variable_step_class),intent(inout)    :: me
        real(wp),intent(in)              :: t
        real(wp),dimension(:),intent(in) :: x

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

    !*********************************************************
        subroutine twobody_event(me,t,x,g)

        !! event function (z = 12,000)

        implicit none

        class(rk_variable_step_class),intent(inout)    :: me
        real(wp),intent(in)              :: t
        real(wp),dimension(:),intent(in) :: x
        real(wp),intent(out)             :: g

        g = x(3) - 12000.0_wp

        end subroutine twobody_event
    !*********************************************************

    end subroutine rk_test_variable_step
!*****************************************************************************************

!*****************************************************************************************
    end module rk_module_variable_step
!*****************************************************************************************
