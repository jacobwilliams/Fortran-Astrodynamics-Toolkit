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
    
    abstract interface
        function func(x) result(f)
            import :: wp
            implicit none
            real(wp),intent(in) :: x
            real(wp) :: f
        end function func
    end interface
    
    public :: fmin
    
    contains
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
!  INPUTS
!
!    ax    left endpoint of initial interval
!    bx    right endpoint of initial interval
!    f     function subprogram which evaluates f(x) for any x in the interval (ax,bx)
!    tol   desired length of the interval of uncertainty of the final result ( >= zero)
!
!  OUTPUT
!
!	fmin abcissa approximating the point where f attains a minimum
!
!	the method used is a combination of golden section search and
!	successive parabolic interpolation. convergence is never much slower
!	than that for a fibonacci search. if f has a continuous second
!	derivative which is positive at the minimum (which is not at ax or
!	bx), then convergence is superlinear, and usually of the order of
!	about 1.324.
!
!	the function f is never evaluated at two points closer together
!	than eps*abs(fmin) + (tol/3), where eps is approximately the square
!	root of the relative machine precision. if f is a unimodal
!	function and the computed values of f are always unimodal when
!	separated by at least eps*abs(x) + (tol/3), then fmin approximates
!	the abcissa of the global minimum of f on the interval ax,bx with
!	an error less than 3*eps*abs(fmin) + tol. if f is not unimodal,
!	then fmin may approximate a local, but perhaps non-global, minimum to
!	the same accuracy.
!
!	this function subprogram is a slightly modified version of the
!	algol 60 procedure localmin given in richard brent, algorithms for
!	minimization without derivatives, prentice - hall, inc. (1973).
!
!  SEE ALSO
!    [1] http://www.netlib.no/netlib/fmm/fmin.f
!
!  SOURCE

    real(wp) function fmin(ax,bx,f,tol)

    implicit none
    
    real(wp),intent(in) :: ax,bx,tol
    procedure(func) :: f
    
    real(wp) :: a,b,d,e,xm,p,q,r,tol1,tol2,u,v,w
    real(wp) :: fu,fv,fw,fx,x
    real(wp) :: abs,sqrt,sign

    real(wp),parameter :: c = (three-sqrt(five))/two    !squared inverse of golden ratio
    real(wp),parameter :: half = 0.5_wp
    real(wp),parameter :: eps = sqrt(epsilon(one))

!
!  initialization
!
      a = ax
      b = bx
      v = a + c*(b - a)
      w = v
      x = v
      e = zero
      fx = f(x)
      fv = fx
      fw = fx
!
!  main loop starts here
!
   20 xm = half*(a + b)
      tol1 = eps*abs(x) + tol/three
      tol2 = two*tol1
!
!  check stopping criterion
!
      if (abs(x - xm) <= (tol2 - half*(b - a))) go to 90
!
! is golden-section necessary
!
      if (abs(e) <= tol1) go to 40
!
!  fit parabola
!
      r = (x - w)*(fx - fv)
      q = (x - v)*(fx - fw)
      p = (x - v)*q - (x - w)*r
      q = two*(q - r)
      if (q > zero) p = -p
      q =  abs(q)
      r = e
      e = d
!
!  is parabola acceptable
!
   30 if (abs(p) >= abs(half*q*r)) go to 40
      if (p <= q*(a - x)) go to 40
      if (p >= q*(b - x)) go to 40
!
!  a parabolic interpolation step
!
      d = p/q
      u = x + d
!
!  f must not be evaluated too close to ax or bx
!
      if ((u - a) < tol2) d = sign(tol1, xm - x)
      if ((b - u) < tol2) d = sign(tol1, xm - x)
      go to 50
!
!  a golden-section step
!
   40 if (x >= xm) e = a - x
      if (x < xm) e = b - x
      d = c*e
!
!  f must not be evaluated too close to x
!
   50 if (abs(d) >= tol1) u = x + d
      if (abs(d) < tol1) u = x + sign(tol1, d)
      fu = f(u)
!
!  update  a, b, v, w, and x
!
      if (fu > fx) go to 60
      if (u >= x) a = x
      if (u < x) b = x
      v = w
      fv = fw
      w = x
      fw = fx
      x = u
      fx = fu
      go to 20
   60 if (u < x) a = u
      if (u >= x) b = u
      if (fu <= fw) go to 70
      if (w == x) go to 70
      if (fu <= fv) go to 80
      if (v == x) go to 80
      if (v == w) go to 80
      go to 20
   70 v = w
      fv = fw
      w = u
      fw = fu
      go to 20
   80 v = u
      fv = fu
      go to 20
!
!  end of main loop
!
   90 fmin = x

    end function fmin
!*****************************************************************************************    
   
!*****************************************************************************************    
    end module brent_module
!*****************************************************************************************    