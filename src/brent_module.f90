!*****************************************************************************************
    module brent_module
!*****************************************************************************************
!****h* FAT/brent_module
!
!  NAME
!    brent_module
!
!  DESCRIPTION
!    Brent algorithms for minimization and root solving without derivatives.
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
        procedure :: set_function
        procedure :: minimize => fmin
        procedure :: find_zero => zeroin
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
    
    !unit test routine:
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
!****f* brent_module/zeroin
!
!  NAME
!    zeroin
!
!  DESCRIPTION
!    a zero of the function f(x) is computed in the interval ax,bx.
!
! INPUTS
!
!    ax : left endpoint of initial interval
!    bx : right endpoint of initial interval
!    f : function subprogram which evaluates f(x) for any x in the interval ax,bx
!    tol : desired length of the interval of uncertainty of the final result (>=0.)
!
!  OUTPUT
!    zeroin : abscissa approximating a zero of f in the interval ax,bx
!
!  NOTES
!    it is assumed that f(ax) and f(bx) have opposite signs
!    this is checked, and an error message is printed if this is not
!    satisfied. zeroin returns a zero x in the given interval
!    ax,bx to within a tolerance 4*macheps*abs(x)+tol, where macheps is
!    the relative machine precision defined as the smallest representable
!    number such that 1.+macheps > 1.
!
!    this function subprogram is a slightly modified translation of
!    the algol 60 procedure zero given in richard brent, algorithms for
!    minimization without derivatives, prentice-hall, inc. (1973).
!
!  SEE ALSO
!    [1] http://www.netlib.org/go/zeroin.f
!
! SOURCE

    real(wp) function zeroin(me,ax,bx,tol)

    implicit none

    class(brent_class),intent(inout) :: me
    real(wp),intent(in) :: ax
    real(wp),intent(in) :: bx
    real(wp),intent(in) :: tol

    real(wp) :: a,b,c,d,e,fa,fb,fc,tol1,xm,p,q,r,s

    real(wp),parameter :: eps = epsilon(one)    !original code had d1mach(4) 
    
    tol1 = eps+one

    a=ax
    b=bx
    fa=me%f(a)
    fb=me%f(b)

    !check that f(ax) and f(bx) have different signs
    if (fa == zero .or. fb == zero) go to 20
    if (fa * (fb/abs(fb)) <= zero) go to 20

    write(*,'(A)') 'Error: f(ax) and f(bx) do not have different signs:'//&
    ' zeroin is aborting'
    return

20  c=a
    fc=fa
    d=b-a
    e=d

30  if (abs(fc)<abs(fb)) then
        a=b
        b=c
        c=a
        fa=fb
        fb=fc
        fc=fa
    end if

40  tol1=two*eps*abs(b)+0.5_wp*tol
    xm = 0.5_wp*(c-b)
    if ((abs(xm)<=tol1).or.(fb==zero)) go to 150

    ! see if a bisection is forced

    if ((abs(e)>=tol1).and.(abs(fa)>abs(fb))) go to 50
    d=xm
    e=d
    go to 110
50  s=fb/fa
    if (a/=c) go to 60

    ! linear interpolation

    p=two*xm*s
    q=one-s
    go to 70

    ! inverse quadratic interpolation

60  q=fa/fc
    r=fb/fc
    p=s*(two*xm*q*(q-r)-(b-a)*(r-one))
    q=(q-one)*(r-one)*(s-one)
70  if (p<=zero) go to 80
    q=-q
    go to 90

80  p=-p
90  s=e
    e=d
    if (((two*p)>=(three*xm*q-abs(tol1*q))).or.(p>=&
        abs(0.5_wp*s*q))) go to 100
    d=p/q
    go to 110

100 d=xm
    e=d
110 a=b
    fa=fb
    if (abs(d)<=tol1) go to 120
    b=b+d
    go to 140

120 if (xm<=zero) go to 130
    b=b+tol1
    go to 140

130 b=b-tol1
140 fb=me%f(b)
    if ((fb*(fc/abs(fc)))>zero) go to 20
    go to 30

150 zeroin=b

    end function zeroin
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
    
    !call zeroin:
    ! [the root is at pi]
    myfunc%i = 0
    write(*,*) 'root of sin(x) at: ', &
                myfunc%minimize(ax+0.0001_wp,bx/two,tol)*180.0_wp/pi,' deg'
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