!*****************************************************************************************
    module brent_module
!*****************************************************************************************
!****h* FAT/brent_module
!
!  NAME
!    brent_module
!
!  DESCRIPTION
!    Brent algorithms for minimization without derivatives.
!
!  SEE ALSO
!    [1] R. Brent, "Algorithms for Minimization Without Derivatives",
!        Prentice-Hall, Inc.m 1973.
!
!*****************************************************************************************    
    
    use kind_module,  only: wp
    use numbers_module

    implicit none
    
    private
    
    !the main class:
    type,public :: brent_class
        procedure(func),pointer :: f => null()  !function to be minimized
    contains
        procedure :: set_function               !set f
        procedure :: minimize => fmin           !call the brent algorithm
    end type brent_class
    
    !interface to the function to be minimized:
    abstract interface
        function func(me,x) result(f)
            import :: wp,brent_class
            implicit none
            class(brent_class),intent(inout) :: me
            real(wp),intent(in) :: x
            real(wp) :: f
        end function func
    end interface
        
    public :: brent_example
    
    contains
!*****************************************************************************************    
      
!*****************************************************************************************    
!****f* brent_module/set_function
!
!  NAME
!    set_function
!
!  DESCRIPTION
!    Set the function to be minimized.
!
!  AUTHOR
!    Jacob Williams, 7/19/2014
!
!  SOURCE

    subroutine set_function(me,f)
    
    implicit none
    
    class(brent_class),intent(inout) :: me
    procedure(func)                  :: f
    
    me%f => f
    
    end subroutine set_function
!*****************************************************************************************    
    
!*****************************************************************************************    
!****f* brent_module/fmin
!
!  NAME
!    fmin
!
!  DESCRIPTION
!    An approximation x to the point where f attains a minimum on
!    the interval (ax,bx) is determined.
!
!    the method used is a combination of golden section search and
!    successive parabolic interpolation. convergence is never much slower
!    than that for a fibonacci search. if f has a continuous second
!    derivative which is positive at the minimum (which is not at ax or
!    bx), then convergence is superlinear, and usually of the order of
!    about 1.324.
!
!    the function f is never evaluated at two points closer together
!    than eps*abs(fmin) + (tol/3), where eps is approximately the square
!    root of the relative machine precision. if f is a unimodal
!    function and the computed values of f are always unimodal when
!    separated by at least eps*abs(x) + (tol/3), then fmin approximates
!    the abcissa of the global minimum of f on the interval ax,bx with
!    an error less than 3*eps*abs(fmin) + tol. if f is not unimodal,
!    then fmin may approximate a local, but perhaps non-global, minimum to
!    the same accuracy.
!
!    this function subprogram is a slightly modified version of the
!    algol 60 procedure localmin given in richard brent, algorithms for
!    minimization without derivatives, prentice - hall, inc. (1973).
!
!  INPUTS
!
!    ax    left endpoint of initial interval
!    bx    right endpoint of initial interval
!    f     function subprogram which evaluates f(x) for any x in the interval (ax,bx)
!    tol   desired length of the interval of uncertainty of the final result ( >= zero)
!
!  OUTPUT
!
!    fmin abcissa approximating the point where f attains a minimum
!
!  SEE ALSO
!    [1] http://www.netlib.no/netlib/fmm/fmin.f
!
!  SOURCE

    real(wp) function fmin(me,ax,bx,tol)

    implicit none

    class(brent_class),intent(inout) :: me
    real(wp),intent(in) :: ax
    real(wp),intent(in) :: bx
    real(wp),intent(in) :: tol

    real(wp) :: a,b,d,e,xm,p,q,r,tol1,tol2,u,v,w
    real(wp) :: fu,fv,fw,fx,x
    real(wp) :: abs,sqrt,sign

    real(wp),parameter :: c = (three-sqrt(five))/two    !squared inverse of golden ratio
    real(wp),parameter :: half = 0.5_wp
    real(wp),parameter :: eps = sqrt(epsilon(one))

    !initialization

    a = ax
    b = bx
    v = a + c*(b - a)
    w = v
    x = v
    e = zero
    fx = me%f(x)
    fv = fx
    fw = fx

    do    !  main loop starts here

        xm = half*(a + b)
        tol1 = eps*abs(x) + tol/three
        tol2 = two*tol1

        !  check stopping criterion

        if (abs(x - xm) <= (tol2 - half*(b - a))) exit

        ! is golden-section necessary

        if (abs(e) <= tol1) go to 40

        !  fit parabola

        r = (x - w)*(fx - fv)
        q = (x - v)*(fx - fw)
        p = (x - v)*q - (x - w)*r
        q = two*(q - r)
        if (q > zero) p = -p
        q =  abs(q)
        r = e
        e = d

        !  is parabola acceptable

        if (abs(p) >= abs(half*q*r)) go to 40
        if (p <= q*(a - x)) go to 40
        if (p >= q*(b - x)) go to 40

        !  a parabolic interpolation step

        d = p/q
        u = x + d

        !  f must not be evaluated too close to ax or bx

        if ((u - a) < tol2) d = sign(tol1, xm - x)
        if ((b - u) < tol2) d = sign(tol1, xm - x)
        go to 50

        !  a golden-section step

    40  if (x >= xm) e = a - x
        if (x < xm) e = b - x
        d = c*e

        !  f must not be evaluated too close to x

    50  if (abs(d) >= tol1) u = x + d
        if (abs(d) < tol1) u = x + sign(tol1, d)
        fu = me%f(u)

        !  update a, b, v, w, and x

        if (fu <= fx) then
            if (u >= x) a = x
            if (u < x) b = x
            v = w
            fv = fw
            w = x
            fw = fx
            x = u
            fx = fu
            cycle
        end if

        if (u < x) a = u
        if (u >= x) b = u

        if (fu <= fw .or. w == x) then
            v = w
            fv = fw
            w = u
            fw = fu
            cycle      
        end if

        if (fu <= fv .or. v == x .or. v == w ) then
            v = u
            fv = fu
            cycle
        end if

    end do    !  end of main loop

    fmin = x

    end function fmin
!*****************************************************************************************

!*****************************************************************************************    
!****f* brent_module/brent_example
!
!  NAME
!    fmin_test
!
!  DESCRIPTION
!    Test of the fmin function
!
!  AUTHOR
!    Jacob Williams, 7/16/2014
!    
!  SOURCE

    subroutine brent_example()
    
    implicit none
    
    real(wp) :: x,f
    
    real(wp),parameter :: ax = zero
    real(wp),parameter :: bx = two*pi
    real(wp),parameter :: tol = 1.0e-6_wp
    
    type,extends(brent_class) :: myfunc_type
        integer :: i = 0    !function counter
    end type myfunc_type   
    type(myfunc_type) :: myfunc
    
    call myfunc%set_function(sin_func)    !set the function
            
    !call fmin:
    ! [the minimum is at 270 deg]
    write(*,*) 'minimum of sin(x) at: ', &
                myfunc%minimize(ax,bx,tol)*180.0_wp/pi,' deg'
    write(*,*) 'number of function calls: ', myfunc%i
    
    contains
 
        function sin_func(me,x) result(f)
        ! Example function to minimize: sin(x)

        implicit none
        
        class(brent_class),intent(inout) :: me
        real(wp),intent(in) :: x
        real(wp) :: f
        
        f = sin(x)
        
        select type (me)
        class is (myfunc_type)
            me%i = me%i + 1 !number of function calls
        end select
        
        end function sin_func
    
    end subroutine brent_example
 !*****************************************************************************************
  
!*****************************************************************************************    
    end module brent_module
!*****************************************************************************************    