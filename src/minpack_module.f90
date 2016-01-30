!*****************************************************************************************
!>
!  Minpack routines for solving a set of nonlinear equations.
!  The two main routines here are [[hybrj]] (user-provided Jacobian) and
!  [[hybrd]] (estimates the Jacobian using finite differences).
!
!### License
!
!  *** Original Minpack License ***
!
!     Minpack Copyright Notice (1999) University of Chicago.  All rights reserved
!
!     Redistribution and use in source and binary forms, with or
!     without modification, are permitted provided that the
!     following conditions are met:
!
!     1. Redistributions of source code must retain the above
!     copyright notice, this list of conditions and the following
!     disclaimer.
!
!     2. Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials
!     provided with the distribution.
!
!     3. The end-user documentation included with the
!     redistribution, if any, must include the following
!     acknowledgment:
!
!        "This product includes software developed by the
!        University of Chicago, as Operator of Argonne National
!        Laboratory.
!
!     Alternately, this acknowledgment may appear in the software
!     itself, if and wherever such third-party acknowledgments
!     normally appear.
!
!     4. WARRANTY DISCLAIMER. THE SOFTWARE IS SUPPLIED "AS IS"
!     WITHOUT WARRANTY OF ANY KIND. THE COPYRIGHT HOLDER, THE
!     UNITED STATES, THE UNITED STATES DEPARTMENT OF ENERGY, AND
!     THEIR EMPLOYEES: (1) DISCLAIM ANY WARRANTIES, EXPRESS OR
!     IMPLIED, INCLUDING BUT NOT LIMITED TO ANY IMPLIED WARRANTIES
!     OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE
!     OR NON-INFRINGEMENT, (2) DO NOT ASSUME ANY LEGAL LIABILITY
!     OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR
!     USEFULNESS OF THE SOFTWARE, (3) DO NOT REPRESENT THAT USE OF
!     THE SOFTWARE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS, (4)
!     DO NOT WARRANT THAT THE SOFTWARE WILL FUNCTION
!     UNINTERRUPTED, THAT IT IS ERROR-FREE OR THAT ANY ERRORS WILL
!     BE CORRECTED.
!
!     5. LIMITATION OF LIABILITY. IN NO EVENT WILL THE COPYRIGHT
!     HOLDER, THE UNITED STATES, THE UNITED STATES DEPARTMENT OF
!     ENERGY, OR THEIR EMPLOYEES: BE LIABLE FOR ANY INDIRECT,
!     INCIDENTAL, CONSEQUENTIAL, SPECIAL OR PUNITIVE DAMAGES OF
!     ANY KIND OR NATURE, INCLUDING BUT NOT LIMITED TO LOSS OF
!     PROFITS OR LOSS OF DATA, FOR ANY REASON WHATSOEVER, WHETHER
!     SUCH LIABILITY IS ASSERTED ON THE BASIS OF CONTRACT, TORT
!     (INCLUDING NEGLIGENCE OR STRICT LIABILITY), OR OTHERWISE,
!     EVEN IF ANY OF SAID PARTIES HAS BEEN WARNED OF THE
!     POSSIBILITY OF SUCH LOSS OR DAMAGES.
!
!  *** Modifications ***
!
!  Modifications for the Fortran Astrodynamics Toolkit are covered
!  under the [following license](https://github.com/jacobwilliams/Fortran-Astrodynamics-Toolkit/blob/master/LICENSE).
!
!### History
!  * Argonne National Laboratory. minpack project. march 1980.
!    burton s. garbow, kenneth e. hillstrom, jorge j. more, john l. nazareth
!  * Jacob Williams, Jan 2016, extensive refactoring into modern Fortran.

    module minpack_module

    use kind_module,    only: wp
    use numbers_module

    implicit none

    abstract interface
        subroutine fcn_hybrd(n,x,fvec,iflag)
            !! function for [[hybrd]].
            !! calculate the functions at `x` and
            !! return this vector in `fvec`.
            import :: wp
            implicit none
            integer,intent(in)                :: n
            real(wp),dimension(n),intent(in)  :: x
            real(wp),dimension(n),intent(out) :: fvec
            integer,intent(inout)             :: iflag  !! the value of `iflag` should not be changed by fcn unless
                                                        !! the user wants to terminate execution of [[hybrd]].
                                                        !! in this case set `iflag` to a negative integer.
        end subroutine fcn_hybrd
        subroutine fcn_hybrj(n,x,fvec,fjac,ldfjac,iflag)
            !! function for [[hybrj]]
            import :: wp
            implicit none
            integer,intent(in)                       :: n
            real(wp),dimension(n),intent(in)         :: x
            integer,intent(in)                       :: ldfjac
            real(wp),dimension(n),intent(out)        :: fvec
            real(wp),dimension(ldfjac,n),intent(out) :: fjac
            integer,intent(inout)                    :: iflag
        end subroutine fcn_hybrj
    end interface

    public :: hybrd,hybrd1
    public :: hybrj,hybrj1

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!     given an m by n matrix a, an n by n nonsingular diagonal
!     matrix d, an m-vector b, and a positive number delta, the
!     problem is to determine the convex combination x of the
!     gauss-newton and scaled gradient directions that minimizes
!     (a*x - b) in the least squares sense, subject to the
!     restriction that the euclidean norm of d*x be at most delta.
!
!     this subroutine completes the solution of the problem
!     if it is provided with the necessary information from the
!     qr factorization of a. that is, if a = q*r, where q has
!     orthogonal columns and r is an upper triangular matrix,
!     then dogleg expects the full upper triangle of r and
!     the first n components of (q transpose)*b.
!
!     the subroutine statement is
!
!       subroutine dogleg(n,r,lr,diag,qtb,delta,x,wa1,wa2)
!
!     where
!
!       n is a positive integer input variable set to the order of r.
!
!       r is an input array of length lr which must contain the upper
!         triangular matrix r stored by rows.
!
!       lr is a positive integer input variable not less than
!         (n*(n+1))/2.
!
!       diag is an input array of length n which must contain the
!         diagonal elements of the matrix d.
!
!       qtb is an input array of length n which must contain the first
!         n elements of the vector (q transpose)*b.
!
!       delta is a positive input variable which specifies an upper
!         bound on the euclidean norm of d*x.
!
!       x is an output array of length n which contains the desired
!         convex combination of the gauss-newton direction and the
!         scaled gradient direction.
!
!       wa1 and wa2 are work arrays of length n.

    subroutine dogleg(n,r,lr,diag,qtb,delta,x,wa1,wa2)

    implicit none

    integer n , lr
    real(wp) delta
    real(wp) r(lr) , diag(n) , qtb(n) , x(n) , wa1(n) , wa2(n)

      integer i , j , jj , jp1 , k , l
      real(wp) alpha , bnorm , epsmch , gnorm , qnorm , sgnorm , sum , temp

      epsmch = dpmpar(1)  ! the machine precision
!
!     first, calculate the gauss-newton direction.
!
      jj = (n*(n+1))/2 + 1
      do k = 1 , n
         j = n - k + 1
         jp1 = j + 1
         jj = jj - k
         l = jj + 1
         sum = zero
         if ( n>=jp1 ) then
            do i = jp1 , n
               sum = sum + r(l)*x(i)
               l = l + 1
            enddo
         endif
         temp = r(jj)
         if ( temp==zero ) then
            l = j
            do i = 1 , j
               temp = max(temp,abs(r(l)))
               l = l + n - i
            enddo
            temp = epsmch*temp
            if ( temp==zero ) temp = epsmch
         endif
         x(j) = (qtb(j)-sum)/temp
      enddo
!
!     test whether the gauss-newton direction is acceptable.
!
      do j = 1 , n
         wa1(j) = zero
         wa2(j) = diag(j)*x(j)
      enddo
      qnorm = enorm(n,wa2)
      if ( qnorm>delta ) then
!
!     the gauss-newton direction is not acceptable.
!     next, calculate the scaled gradient direction.
!
         l = 1
         do j = 1 , n
            temp = qtb(j)
            do i = j , n
               wa1(i) = wa1(i) + r(l)*temp
               l = l + 1
            enddo
            wa1(j) = wa1(j)/diag(j)
         enddo
!
!     calculate the norm of the scaled gradient and test for
!     the special case in which the scaled gradient is zero.
!
         gnorm = enorm(n,wa1)
         sgnorm = zero
         alpha = delta/qnorm
         if ( gnorm/=zero ) then
!
!     calculate the point along the scaled gradient
!     at which the quadratic is minimized.
!
            do j = 1 , n
               wa1(j) = (wa1(j)/gnorm)/diag(j)
            enddo
            l = 1
            do j = 1 , n
               sum = zero
               do i = j , n
                  sum = sum + r(l)*wa1(i)
                  l = l + 1
               enddo
               wa2(j) = sum
            enddo
            temp = enorm(n,wa2)
            sgnorm = (gnorm/temp)/temp
!
!     test whether the scaled gradient direction is acceptable.
!
            alpha = zero
            if ( sgnorm<delta ) then
!
!     the scaled gradient direction is not acceptable.
!     finally, calculate the point along the dogleg
!     at which the quadratic is minimized.
!
               bnorm = enorm(n,qtb)
               temp = (bnorm/gnorm)*(bnorm/qnorm)*(sgnorm/delta)
               temp = temp - (delta/qnorm)*(sgnorm/delta)**2 + &
                      sqrt((temp-(delta/qnorm))**2+&
                      (one-(delta/qnorm)**2)*(one-(sgnorm/delta)**2))
               alpha = ((delta/qnorm)*(one-(sgnorm/delta)**2))/temp
            endif
         endif
!
!     form appropriate convex combination of the gauss-newton
!     direction and the scaled gradient direction.
!
         temp = (one-alpha)*min(sgnorm,delta)
         do j = 1 , n
            x(j) = temp*wa1(j) + alpha*x(j)
         enddo
      endif

    end subroutine dogleg
!*****************************************************************************************

!*****************************************************************************************
!>
!  Replacement for the original Minpack routine.

    real(wp) function dpmpar(i)
    implicit none

    integer,intent(in) :: i

    real(wp),dimension(3),parameter :: dmach = [epsilon(1.0_wp),&
                                                  tiny(1.0_wp),&
                                                  huge(1.0_wp)]

    dpmpar = dmach(i)

    end function dpmpar
!*****************************************************************************************

!*****************************************************************************************
!>
!     given an n-vector x, this function calculates the
!     euclidean norm of x.
!
!     the euclidean norm is computed by accumulating the sum of
!     squares in three different sums. the sums of squares for the
!     small and large components are scaled so that no overflows
!     occur. non-destructive underflows are permitted. underflows
!     and overflows do not occur in the computation of the unscaled
!     sum of squares for the intermediate components.
!     the definitions of small, intermediate and large components
!     depend on two constants, rdwarf and rgiant. the main
!     restrictions on these constants are that rdwarf**2 not
!     underflow and rgiant**2 not overflow. the constants
!     given here are suitable for every known computer.

    real(wp) function enorm(n,x)

    implicit none

    integer,intent(in) :: n                 !! size of `x`
    real(wp),dimension(n),intent(in) :: x   !! input array

      integer i
      real(wp) agiant , floatn , s1 , s2 , s3 , xabs , x1max , x3max
      real(wp),parameter :: rdwarf = 3.834e-20_wp
      real(wp),parameter :: rgiant = 1.304e19_wp

      s1 = zero
      s2 = zero
      s3 = zero
      x1max = zero
      x3max = zero
      floatn = n
      agiant = rgiant/floatn
      do i = 1 , n
         xabs = abs(x(i))
         if ( xabs>rdwarf .and. xabs<agiant ) then
!
!           sum for intermediate components.
!
            s2 = s2 + xabs**2
         elseif ( xabs<=rdwarf ) then
!
!              sum for small components.
!
            if ( xabs<=x3max ) then
               if ( xabs/=zero ) s3 = s3 + (xabs/x3max)**2
            else
               s3 = one + s3*(x3max/xabs)**2
               x3max = xabs
            endif
!
!              sum for large components.
!
         elseif ( xabs<=x1max ) then
            s1 = s1 + (xabs/x1max)**2
         else
            s1 = one + s1*(x1max/xabs)**2
            x1max = xabs
         endif
      enddo
!
!     calculation of norm.
!
      if ( s1/=zero ) then
         enorm = x1max*sqrt(s1+(s2/x1max)/x1max)
      elseif ( s2==zero ) then
         enorm = x3max*sqrt(s3)
      else
         if ( s2>=x3max ) enorm = sqrt(s2*(one+(x3max/s2)*(x3max*s3)))
         if ( s2<x3max ) enorm = sqrt(x3max*((s2/x3max)+(x3max*s3)))
      endif

    end function enorm
!*****************************************************************************************

!*****************************************************************************************
!>
!     this subroutine computes a forward-difference approximation
!     to the n by n jacobian matrix associated with a specified
!     problem of n functions in n variables. if the jacobian has
!     a banded form, then function evaluations are saved by only
!     approximating the nonzero terms.
!
!     the subroutine statement is
!
!       subroutine fdjac1(fcn,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,wa1,wa2)
!
!     where
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions. fcn must be declared
!         in an external statement in the user calling
!         program, and should be written as follows.
!
!         subroutine fcn(n,x,fvec,iflag)
!         integer n,iflag
!         real(wp) x(n),fvec(n)
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ----------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of fdjac1.
!         in this case set iflag to a negative integer.
!
!       n is a positive integer input variable set to the number
!         of functions and variables.
!
!       x is an input array of length n.
!
!       fvec is an input array of length n which must contain the
!         functions evaluated at x.
!
!       fjac is an output n by n array which contains the
!         approximation to the jacobian matrix evaluated at x.
!
!       ldfjac is a positive integer input variable not less than n
!         which specifies the leading dimension of the array fjac.
!
!       iflag is an integer variable which can be used to terminate
!         the execution of fdjac1. see description of fcn.
!
!       ml is a nonnegative integer input variable which specifies
!         the number of subdiagonals within the band of the
!         jacobian matrix. if the jacobian is not banded, set
!         ml to at least n - 1.
!
!       epsfcn is an input variable used in determining a suitable
!         step length for the forward-difference approximation. this
!         approximation assumes that the relative errors in the
!         functions are of the order of epsfcn. if epsfcn is less
!         than the machine precision, it is assumed that the relative
!         errors in the functions are of the order of the machine
!         precision.
!
!       mu is a nonnegative integer input variable which specifies
!         the number of superdiagonals within the band of the
!         jacobian matrix. if the jacobian is not banded, set
!         mu to at least n - 1.
!
!       wa1 and wa2 are work arrays of length n. if ml + mu + 1 is at
!         least n, then the jacobian is considered dense, and wa2 is
!         not referenced.

    subroutine fdjac1(fcn,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,wa1,wa2)

    implicit none

    integer n , ldfjac , iflag , ml , mu
    real(wp) epsfcn
    real(wp) x(n) , fvec(n) , fjac(ldfjac,n) , wa1(n) , wa2(n)
    procedure(fcn_hybrd) :: fcn

      integer i , j , k , msum
      real(wp) eps , epsmch , h , temp

      epsmch = dpmpar(1) ! the machine precision
!
      eps = sqrt(max(epsfcn,epsmch))
      msum = ml + mu + 1
      if ( msum<n ) then
!
!        computation of banded approximate jacobian.
!
         do k = 1 , msum
            do j = k , n , msum
               wa2(j) = x(j)
               h = eps*abs(wa2(j))
               if ( h==zero ) h = eps
               x(j) = wa2(j) + h
            enddo
            call fcn(n,x,wa1,iflag)
            if ( iflag<0 ) return
            do j = k , n , msum
               x(j) = wa2(j)
               h = eps*abs(wa2(j))
               if ( h==zero ) h = eps
               do i = 1 , n
                  fjac(i,j) = zero
                  if ( i>=j-mu .and. i<=j+ml ) fjac(i,j) = (wa1(i)-fvec(i))/h
               enddo
            enddo
         enddo
      else
!
!        computation of dense approximate jacobian.
!
         do j = 1 , n
            temp = x(j)
            h = eps*abs(temp)
            if ( h==zero ) h = eps
            x(j) = temp + h
            call fcn(n,x,wa1,iflag)
            if ( iflag<0 ) return
            x(j) = temp
            do i = 1 , n
               fjac(i,j) = (wa1(i)-fvec(i))/h
            enddo
         enddo
      endif

    end subroutine fdjac1
!*****************************************************************************************

!*****************************************************************************************
!>
!  The purpose of hybrd is to find a zero of a system of
!  n nonlinear functions in n variables by a modification
!  of the powell hybrid method. the user must provide a
!  subroutine which calculates the functions. the jacobian is
!  then calculated by a forward-difference approximation.
!
!### Characteristics of the algorithm.
!  HYBRD is a modification of the Powell hybrid method.  Two of its
!  main characteristics involve the choice of the correction as a
!  convex combination of the Newton and scaled gradient directions
!  and the updating of the Jacobian by the rank-1 method of Broy-
!  den.  The choice of the correction guarantees (under reasonable
!  conditions) global convergence for starting points far from the
!  solution and a fast rate of convergence.  The Jacobian is
!  approximated by forward differences at the starting point, but
!  forward differences are not used again until the rank-1 method
!  fails to produce satisfactory progress.
!
!### References
!  * M. J. D. Powell, A Hybrid Method for Nonlinear Equations.
!    Numerical Methods for Nonlinear Algebraic Equations,
!    P. Rabinowitz, editor. Gordon and Breach, 1970.

    subroutine hybrd(fcn,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,diag,mode, &
                     factor,nprint,info,nfev,fjac,ldfjac,r,lr,qtf,wa1,&
                     wa2,wa3,wa4)

    implicit none

    procedure(fcn_hybrd) :: fcn             !! user-supplied subroutine which calculates the functions
    integer,intent(in) :: n                 !! a positive integer input variable set to the number
                                            !! of functions and variables.
    integer,intent(in) :: maxfev            !! a positive integer input variable. termination
                                            !! occurs when the number of calls to `fcn` is at least `maxfev`
                                            !! by the end of an iteration.
    integer,intent(in) :: ml                !! a nonnegative integer input variable which specifies
                                            !! the number of subdiagonals within the band of the
                                            !! jacobian matrix. if the jacobian is not banded, set
                                            !! `ml` to at least `n - 1`.
    integer,intent(in) :: mu                !! a nonnegative integer input variable which specifies
                                            !! the number of superdiagonals within the band of the
                                            !! jacobian matrix. if the jacobian is not banded, set
                                            !! `mu` to at least` n - 1`.
    integer,intent(in) :: mode              !! if `mode = 1`, the
                                            !! variables will be scaled internally. if `mode = 2`,
                                            !! the scaling is specified by the input `diag`. other
                                            !! values of `mode` are equivalent to `mode = 1`.
    integer,intent(in)  :: nprint           !! an integer input variable that enables controlled
                                            !! printing of iterates if it is positive. in this case,
                                            !! `fcn` is called with `iflag = 0` at the beginning of the first
                                            !! iteration and every `nprint` iterations thereafter and
                                            !! immediately prior to return, with `x` and `fvec` available
                                            !! for printing. if `nprint` is not positive, no special calls
                                            !! of `fcn` with `iflag = 0` are made.
    integer,intent(out) :: info             !! an integer output variable. if the user has
                                            !! terminated execution, `info` is set to the (negative)
                                            !! value of `iflag`. see description of `fcn`. otherwise,
                                            !! `info` is set as follows:
                                            !!  * ***info = 0*** improper input parameters.
                                            !!  * ***info = 1*** relative error between two consecutive iterates
                                            !!    is at most `xtol`.
                                            !!  * ***info = 2*** number of calls to `fcn` has reached or exceeded
                                            !!    `maxfev`.
                                            !!  * ***info = 3*** `xtol` is too small. no further improvement in
                                            !!    the approximate solution `x` is possible.
                                            !!  * ***info = 4*** iteration is not making good progress, as
                                            !!    measured by the improvement from the last
                                            !!    five jacobian evaluations.
                                            !!  * ***info = 5*** iteration is not making good progress, as
                                            !!    measured by the improvement from the last
                                            !!    ten iterations.
    integer,intent(out) :: nfev             !! output variable set to the number of calls to `fcn`.
    integer,intent(in):: ldfjac             !! a positive integer input variable not less than `n`
                                            !! which specifies the leading dimension of the array `fjac`.
    integer,intent(in) :: lr                !! a positive integer input variable not less than `(n*(n+1))/2`.
    real(wp),intent(in) :: xtol             !! a nonnegative input variable. termination
                                            !! occurs when the relative error between two consecutive
                                            !! iterates is at most `xtol`.
    real(wp),intent(in) :: epsfcn           !! an input variable used in determining a suitable
                                            !! step length for the forward-difference approximation. this
                                            !! approximation assumes that the relative errors in the
                                            !! functions are of the order of `epsfcn`. if `epsfcn` is less
                                            !! than the machine precision, it is assumed that the relative
                                            !! errors in the functions are of the order of the machine
                                            !! precision.
    real(wp),intent(in) :: factor           !! a positive input variable used in determining the
                                            !! initial step bound. this bound is set to the product of
                                            !! `factor` and the euclidean norm of `diag*x` if nonzero, or else
                                            !! to `factor` itself. in most cases factor should lie in the
                                            !! interval (.1,100.). 100. is a generally recommended value.
    real(wp),intent(inout) :: x(n)          !! array of length n. on input `x` must contain
                                            !! an initial estimate of the solution vector. on output `x`
                                            !! contains the final estimate of the solution vector.
    real(wp),intent(out) :: fvec(n)         !! an output array of length `n` which contains
                                            !! the functions evaluated at the output `x`.
    real(wp),intent(inout) :: diag(n)       !! an array of length `n`. if `mode = 1` (see
                                            !! below), `diag` is internally set. if `mode = 2`, `diag`
                                            !! must contain positive entries that serve as
                                            !! multiplicative scale factors for the variables.
    real(wp),intent(out) :: fjac(ldfjac,n)  !! array which contains the
                                            !! orthogonal matrix `q` produced by the QR factorization
                                            !! of the final approximate jacobian.
    real(wp),intent(out) :: r(lr)           !! an output array which contains the
                                            !! upper triangular matrix produced by the QR factorization
                                            !! of the final approximate jacobian, stored rowwise.
    real(wp),intent(out) :: qtf(n)          !! an output array of length `n` which contains
                                            !! the vector `(q transpose)*fvec`.
    real(wp),intent(inout) :: wa1(n)  !! work array
    real(wp),intent(inout) :: wa2(n)  !! work array
    real(wp),intent(inout) :: wa3(n)  !! work array
    real(wp),intent(inout) :: wa4(n)  !! work array

    integer :: i , iflag , iter , j , jm1 , l , msum , ncfail , ncsuc , nslow1 , nslow2
    integer :: iwa(1)
    logical :: jeval , sing
    real(wp) :: actred , delta , epsmch , fnorm , fnorm1 , &
                  pnorm , prered , ratio ,&
                  sum , temp , xnorm

    real(wp),parameter :: p1    = 1.0e-1_wp
    real(wp),parameter :: p5    = 5.0e-1_wp
    real(wp),parameter :: p001  = 1.0e-3_wp
    real(wp),parameter :: p0001 = 1.0e-4_wp

    epsmch = dpmpar(1)  ! the machine precision

    info = 0
    iflag = 0
    nfev = 0
    !
    !     check the input parameters for errors.
    !
    if ( n<=0 .or. xtol<zero .or. maxfev<=0 .or. ml<0 .or. mu<0 .or.  &
         factor<=zero .or. ldfjac<n .or. lr<(n*(n+1))/2 ) goto 300
    if ( mode==2 ) then
       do j = 1 , n
          if ( diag(j)<=zero ) goto 300
       enddo
    endif
    !
    !     evaluate the function at the starting point
    !     and calculate its norm.
    !
    iflag = 1
    call fcn(n,x,fvec,iflag)
    nfev = 1
    if ( iflag<0 ) goto 300
    fnorm = enorm(n,fvec)
    !
    !     determine the number of calls to fcn needed to compute
    !     the jacobian matrix.
    !
    msum = min0(ml+mu+1,n)
    !
    !     initialize iteration counter and monitors.
    !
    iter = 1
    ncsuc = 0
    ncfail = 0
    nslow1 = 0
    nslow2 = 0
    !
    !     beginning of the outer loop.
    !
    100  jeval = .true.
    !
    !        calculate the jacobian matrix.
    !
    iflag = 2
    call fdjac1(fcn,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,wa1,wa2)
    nfev = nfev + msum
    if ( iflag<0 ) goto 300
    !
    !        compute the qr factorization of the jacobian.
    !
    call qrfac(n,n,fjac,ldfjac,.false.,iwa,1,wa1,wa2,wa3)
    !
    !        on the first iteration and if mode is 1, scale according
    !        to the norms of the columns of the initial jacobian.
    !
    if ( iter==1 ) then
       if ( mode/=2 ) then
          do j = 1 , n
             diag(j) = wa2(j)
             if ( wa2(j)==zero ) diag(j) = one
          enddo
       endif
    !
    !        on the first iteration, calculate the norm of the scaled x
    !        and initialize the step bound delta.
    !
       do j = 1 , n
          wa3(j) = diag(j)*x(j)
       enddo
       xnorm = enorm(n,wa3)
       delta = factor*xnorm
       if ( delta==zero ) delta = factor
    endif
    !
    !        form (q transpose)*fvec and store in qtf.
    !
    do i = 1 , n
       qtf(i) = fvec(i)
    enddo
    do j = 1 , n
       if ( fjac(j,j)/=zero ) then
          sum = zero
          do i = j , n
             sum = sum + fjac(i,j)*qtf(i)
          enddo
          temp = -sum/fjac(j,j)
          do i = j , n
             qtf(i) = qtf(i) + fjac(i,j)*temp
          enddo
       endif
    enddo
    !
    !        copy the triangular factor of the qr factorization into r.
    !
    sing = .false.
    do j = 1 , n
       l = j
       jm1 = j - 1
       if ( jm1>=1 ) then
          do i = 1 , jm1
             r(l) = fjac(i,j)
             l = l + n - i
          enddo
       endif
       r(l) = wa1(j)
       if ( wa1(j)==zero ) sing = .true.
    enddo
    !
    !        accumulate the orthogonal factor in fjac.
    !
    call qform(n,n,fjac,ldfjac,wa1)
    !
    !        rescale if necessary.
    !
    if ( mode/=2 ) then
       do j = 1 , n
          diag(j) = dmax1(diag(j),wa2(j))
       enddo
    endif
    !
    !        beginning of the inner loop.
    !
    !
    !           if requested, call fcn to enable printing of iterates.
    !
    200  if ( nprint>0 ) then
       iflag = 0
       if ( mod(iter-1,nprint)==0 ) call fcn(n,x,fvec,iflag)
       if ( iflag<0 ) goto 300
    endif
    !
    !           determine the direction p.
    !
    call dogleg(n,r,lr,diag,qtf,delta,wa1,wa2,wa3)
    !
    !           store the direction p and x + p. calculate the norm of p.
    !
    do j = 1 , n
       wa1(j) = -wa1(j)
       wa2(j) = x(j) + wa1(j)
       wa3(j) = diag(j)*wa1(j)
    enddo
    pnorm = enorm(n,wa3)
    !
    !           on the first iteration, adjust the initial step bound.
    !
    if ( iter==1 ) delta = dmin1(delta,pnorm)
    !
    !           evaluate the function at x + p and calculate its norm.
    !
    iflag = 1
    call fcn(n,wa2,wa4,iflag)
    nfev = nfev + 1
    if ( iflag>=0 ) then
       fnorm1 = enorm(n,wa4)
    !
    !           compute the scaled actual reduction.
    !
       actred = -one
       if ( fnorm1<fnorm ) actred = one - (fnorm1/fnorm)**2
    !
    !           compute the scaled predicted reduction.
    !
       l = 1
       do i = 1 , n
          sum = zero
          do j = i , n
             sum = sum + r(l)*wa1(j)
             l = l + 1
          enddo
          wa3(i) = qtf(i) + sum
       enddo
       temp = enorm(n,wa3)
       prered = zero
       if ( temp<fnorm ) prered = one - (temp/fnorm)**2
    !
    !           compute the ratio of the actual to the predicted
    !           reduction.
    !
       ratio = zero
       if ( prered>zero ) ratio = actred/prered
    !
    !           update the step bound.
    !
       if ( ratio>=p1 ) then
          ncfail = 0
          ncsuc = ncsuc + 1
          if ( ratio>=p5 .or. ncsuc>1 ) delta = dmax1(delta,pnorm/p5)
          if ( dabs(ratio-one)<=p1 ) delta = pnorm/p5
       else
          ncsuc = 0
          ncfail = ncfail + 1
          delta = p5*delta
       endif
    !
    !           test for successful iteration.
    !
       if ( ratio>=p0001 ) then
    !
    !           successful iteration. update x, fvec, and their norms.
    !
          do j = 1 , n
             x(j) = wa2(j)
             wa2(j) = diag(j)*x(j)
             fvec(j) = wa4(j)
          enddo
          xnorm = enorm(n,wa2)
          fnorm = fnorm1
          iter = iter + 1
       endif
    !
    !           determine the progress of the iteration.
    !
       nslow1 = nslow1 + 1
       if ( actred>=p001 ) nslow1 = 0
       if ( jeval ) nslow2 = nslow2 + 1
       if ( actred>=p1 ) nslow2 = 0
    !
    !           test for convergence.
    !
       if ( delta<=xtol*xnorm .or. fnorm==zero ) info = 1
       if ( info==0 ) then
    !
    !           tests for termination and stringent tolerances.
    !
          if ( nfev>=maxfev ) info = 2
          if ( p1*dmax1(p1*delta,pnorm)<=epsmch*xnorm ) info = 3
          if ( nslow2==5 ) info = 4
          if ( nslow1==10 ) info = 5
          if ( info==0 ) then
    !
    !           criterion for recalculating jacobian approximation
    !           by forward differences.
    !
             if ( ncfail==2 ) goto 100
    !
    !           calculate the rank one modification to the jacobian
    !           and update qtf if necessary.
    !
             do j = 1 , n
                sum = zero
                do i = 1 , n
                   sum = sum + fjac(i,j)*wa4(i)
                enddo
                wa2(j) = (sum-wa3(j))/pnorm
                wa1(j) = diag(j)*((diag(j)*wa1(j))/pnorm)
                if ( ratio>=p0001 ) qtf(j) = sum
             enddo
    !
    !           compute the qr factorization of the updated jacobian.
    !
             call r1updt(n,n,r,lr,wa1,wa2,wa3,sing)
             call r1mpyq(n,n,fjac,ldfjac,wa2,wa3)
             call r1mpyq(1,n,qtf,1,wa2,wa3)
    !
    !           end of the inner loop.
    !
             jeval = .false.
    !
    !        end of the outer loop.
    !
             goto 200
          endif
       endif
    endif
    !
    !     termination, either normal or user imposed.
    !
    300  if ( iflag<0 ) info = iflag
    iflag = 0
    if ( nprint>0 ) call fcn(n,x,fvec,iflag)

    end subroutine hybrd
!*****************************************************************************************

!*****************************************************************************************
!>
!  the purpose of `hybrd1` is to find a zero of a system of
!  `n` nonlinear functions in `n` variables by a modification
!  of the powell hybrid method. this is done by using the
!  more general nonlinear equation solver [[hybrd]]. the user
!  must provide a subroutine which calculates the functions.
!  the jacobian is then calculated by a forward-difference
!  approximation.

    subroutine hybrd1(fcn,n,x,fvec,tol,info)

    implicit none

    procedure(fcn_hybrd)                :: fcn      !! user-supplied subroutine which calculates the functions
    integer,intent(in)                  :: n        !! a positive integer input variable set to the number
                                                    !! of functions and variables.
    integer,intent(out)                 :: info     !! an integer output variable. if the user has
                                                    !! terminated execution, info is set to the (negative)
                                                    !! value of `iflag`. see description of `fcn`. otherwise,
                                                    !! `info` is set as follows:
                                                    !!  ***info = 0*** improper input parameters.
                                                    !!  ***info = 1*** algorithm estimates that the relative error
                                                    !!  between `x` and the solution is at most `tol`.
                                                    !!  ***info = 2*** number of calls to `fcn` has reached or exceeded
                                                    !!  `200*(n+1)`.
                                                    !!  ***info = 3*** `tol` is too small. no further improvement in
                                                    !!  the approximate solution `x` is possible.
                                                    !!  ***info = 4*** iteration is not making good progress.
    real(wp),intent(in)                 :: tol      !! a nonnegative input variable. termination occurs
                                                    !! when the algorithm estimates that the relative error
                                                    !! between `x` and the solution is at most `tol`.
    real(wp),dimension(n),intent(inout) :: x        !! an array of length `n`. on input `x` must contain
                                                    !! an initial estimate of the solution vector. on output `x`
                                                    !! contains the final estimate of the solution vector.
    real(wp),dimension(n),intent(out)   :: fvec     !! an output array of length `n` which contains
                                                    !! the functions evaluated at the output `x`.

    integer :: lwa                           !! length of `wa` work array
    real(wp),dimension(:),allocatable :: wa  !! work array
    integer :: index , j , lr , maxfev , ml , mode , mu , nfev , nprint
    real(wp) :: epsfcn , xtol
    real(wp),dimension(n) :: diag

    real(wp),parameter :: factor = 100.0_wp

    info = 0

    ! check the input parameters for errors.
    if ( n>0 .and. tol>=zero ) then

        !set up inputs:
        lwa     = (n*(3*n+13))/2    ! the work array was formerly an input
        allocate(wa(lwa))           !
        wa      = 0.0_wp            !
        maxfev  = 200*(n+1)
        xtol    = tol
        ml      = n - 1
        mu      = n - 1
        epsfcn  = zero
        mode    = 2
        diag    = one
        nprint  = 0
        lr      = (n*(n+1))/2
        index   = 6*n + lr

        ! call hybrd:
        call hybrd(fcn,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,diag,mode,&
                   factor,nprint,info,nfev,wa(index+1),n,wa(6*n+1),lr,&
                   wa(n+1),wa(2*n+1),wa(3*n+1),wa(4*n+1),wa(5*n+1))
        if ( info==5 ) info = 4

        deallocate(wa)

    endif

    end subroutine hybrd1
!*****************************************************************************************

!*****************************************************************************************
!>
!  the purpose of hybrj is to find a zero of a system of
!  n nonlinear functions in n variables by a modification
!  of the powell hybrid method. the user must provide a
!  subroutine which calculates the functions and the jacobian.
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions and the jacobian. fcn must
!         be declared in an external statement in the user
!         calling program, and should be written as follows.
!
!         subroutine fcn(n,x,fvec,fjac,ldfjac,iflag)
!         integer n,ldfjac,iflag
!         real(wp) x(n),fvec(n),fjac(ldfjac,n)
!         ----------
!         if iflag = 1 calculate the functions at x and
!         return this vector in fvec. do not alter fjac.
!         if iflag = 2 calculate the jacobian at x and
!         return this matrix in fjac. do not alter fvec.
!         ---------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of hybrj.
!         in this case set iflag to a negative integer.
!
!       n is a positive integer input variable set to the number
!         of functions and variables.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length n which contains
!         the functions evaluated at the output x.
!
!       fjac is an output n by n array which contains the
!         orthogonal matrix q produced by the qr factorization
!         of the final approximate jacobian.
!
!       ldfjac is a positive integer input variable not less than n
!         which specifies the leading dimension of the array fjac.
!
!       xtol is a nonnegative input variable. termination
!         occurs when the relative error between two consecutive
!         iterates is at most xtol.
!
!       maxfev is a positive integer input variable. termination
!         occurs when the number of calls to fcn with iflag = 1
!         has reached maxfev.
!
!       diag is an array of length n. if mode = 1 (see
!         below), diag is internally set. if mode = 2, diag
!         must contain positive entries that serve as
!         multiplicative scale factors for the variables.
!
!       mode is an integer input variable. if mode = 1, the
!         variables will be scaled internally. if mode = 2,
!         the scaling is specified by the input diag. other
!         values of mode are equivalent to mode = 1.
!
!       factor is a positive input variable used in determining the
!         initial step bound. this bound is set to the product of
!         factor and the euclidean norm of diag*x if nonzero, or else
!         to factor itself. in most cases factor should lie in the
!         interval (.1,100.). 100. is a generally recommended value.
!
!       nprint is an integer input variable that enables controlled
!         printing of iterates if it is positive. in this case,
!         fcn is called with iflag = 0 at the beginning of the first
!         iteration and every nprint iterations thereafter and
!         immediately prior to return, with x and fvec available
!         for printing. fvec and fjac should not be altered.
!         if nprint is not positive, no special calls of fcn
!         with iflag = 0 are made.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0   improper input parameters.
!
!         info = 1   relative error between two consecutive iterates
!                    is at most xtol.
!
!         info = 2   number of calls to fcn with iflag = 1 has
!                    reached maxfev.
!
!         info = 3   xtol is too small. no further improvement in
!                    the approximate solution x is possible.
!
!         info = 4   iteration is not making good progress, as
!                    measured by the improvement from the last
!                    five jacobian evaluations.
!
!         info = 5   iteration is not making good progress, as
!                    measured by the improvement from the last
!                    ten iterations.
!
!       nfev is an integer output variable set to the number of
!         calls to fcn with iflag = 1.
!
!       njev is an integer output variable set to the number of
!         calls to fcn with iflag = 2.
!
!       r is an output array of length lr which contains the
!         upper triangular matrix produced by the qr factorization
!         of the final approximate jacobian, stored rowwise.
!
!       lr is a positive integer input variable not less than
!         (n*(n+1))/2.
!
!       qtf is an output array of length n which contains
!         the vector (q transpose)*fvec.
!
!       wa1, wa2, wa3, and wa4 are work arrays of length n.
!
!### Characteristics of the algorithm
!  HYBRJ is a modification of the Powell hybrid method.  Two of its
!  main characteristics involve the choice of the correction as a
!  convex combination of the Newton and scaled gradient directions
!  and the updating of the Jacobian by the rank-1 method of Broy-
!  den.  The choice of the correction guarantees (under reasonable
!  conditions) global convergence for starting points far from the
!  solution and a fast rate of convergence.  The Jacobian is calcu-
!  lated at the starting point, but it is not recalculated until
!  the rank-1 method fails to produce satisfactory progress.
!
!### References.
!  * M. J. D. Powell, A Hybrid Method for Nonlinear Equations.
!   Numerical Methods for Nonlinear Algebraic Equations,
!   P. Rabinowitz, editor. Gordon and Breach, 1970.

    subroutine hybrj(fcn,n,x,fvec,fjac,ldfjac,xtol,maxfev,diag,mode,  &
                   & factor,nprint,info,nfev,njev,r,lr,qtf,wa1,wa2,   &
                   & wa3,wa4)

    implicit none

    procedure(fcn_hybrj) :: fcn
    integer n , ldfjac , maxfev , mode , nprint , info , nfev , njev , lr
    real(wp) xtol , factor
    real(wp) x(n) , fvec(n) , fjac(ldfjac,n) , diag(n) , &
             r(lr) , qtf(n) , wa1(n) , wa2(n) , wa3(n) , &
             wa4(n)

      integer i , iflag , iter , j , jm1 , l , ncfail , ncsuc , nslow1 , nslow2
      integer iwa(1)
      logical jeval , sing
      real(wp) actred , delta , epsmch , fnorm , fnorm1 , &
               pnorm , prered , ratio ,&
               sum , temp , xnorm

      real(wp),parameter :: p1    = 1.0e-1_wp
      real(wp),parameter :: p5    = 5.0e-1_wp
      real(wp),parameter :: p001  = 1.0e-3_wp
      real(wp),parameter :: p0001 = 1.0e-4_wp

      epsmch = dpmpar(1)  ! the machine precision
!
      info = 0
      iflag = 0
      nfev = 0
      njev = 0
!
!     check the input parameters for errors.
!
      if ( n<=0 .or. ldfjac<n .or. xtol<zero .or. maxfev<=0 .or. &
           factor<=zero .or. lr<(n*(n+1))/2 ) goto 300
      if ( mode==2 ) then
         do j = 1 , n
            if ( diag(j)<=zero ) goto 300
         enddo
      endif
!
!     evaluate the function at the starting point
!     and calculate its norm.
!
      iflag = 1
      call fcn(n,x,fvec,fjac,ldfjac,iflag)
      nfev = 1
      if ( iflag<0 ) goto 300
      fnorm = enorm(n,fvec)
!
!     initialize iteration counter and monitors.
!
      iter = 1
      ncsuc = 0
      ncfail = 0
      nslow1 = 0
      nslow2 = 0
!
!     beginning of the outer loop.
!
 100  jeval = .true.
!
!        calculate the jacobian matrix.
!
      iflag = 2
      call fcn(n,x,fvec,fjac,ldfjac,iflag)
      njev = njev + 1
      if ( iflag<0 ) goto 300
!
!        compute the qr factorization of the jacobian.
!
      call qrfac(n,n,fjac,ldfjac,.false.,iwa,1,wa1,wa2,wa3)
!
!        on the first iteration and if mode is 1, scale according
!        to the norms of the columns of the initial jacobian.
!
      if ( iter==1 ) then
         if ( mode/=2 ) then
            do j = 1 , n
               diag(j) = wa2(j)
               if ( wa2(j)==zero ) diag(j) = one
            enddo
         endif
!
!        on the first iteration, calculate the norm of the scaled x
!        and initialize the step bound delta.
!
         do j = 1 , n
            wa3(j) = diag(j)*x(j)
         enddo
         xnorm = enorm(n,wa3)
         delta = factor*xnorm
         if ( delta==zero ) delta = factor
      endif
!
!        form (q transpose)*fvec and store in qtf.
!
      do i = 1 , n
         qtf(i) = fvec(i)
      enddo
      do j = 1 , n
         if ( fjac(j,j)/=zero ) then
            sum = zero
            do i = j , n
               sum = sum + fjac(i,j)*qtf(i)
            enddo
            temp = -sum/fjac(j,j)
            do i = j , n
               qtf(i) = qtf(i) + fjac(i,j)*temp
            enddo
         endif
      enddo
!
!        copy the triangular factor of the qr factorization into r.
!
      sing = .false.
      do j = 1 , n
         l = j
         jm1 = j - 1
         if ( jm1>=1 ) then
            do i = 1 , jm1
               r(l) = fjac(i,j)
               l = l + n - i
            enddo
         endif
         r(l) = wa1(j)
         if ( wa1(j)==zero ) sing = .true.
      enddo
!
!        accumulate the orthogonal factor in fjac.
!
      call qform(n,n,fjac,ldfjac,wa1)
!
!        rescale if necessary.
!
      if ( mode/=2 ) then
         do j = 1 , n
            diag(j) = max(diag(j),wa2(j))
         enddo
      endif
!
!        beginning of the inner loop.
!
!
!           if requested, call fcn to enable printing of iterates.
!
 200  if ( nprint>0 ) then
         iflag = 0
         if ( mod(iter-1,nprint)==0 )                                   &
            & call fcn(n,x,fvec,fjac,ldfjac,iflag)
         if ( iflag<0 ) goto 300
      endif
!
!           determine the direction p.
!
      call dogleg(n,r,lr,diag,qtf,delta,wa1,wa2,wa3)
!
!           store the direction p and x + p. calculate the norm of p.
!
      do j = 1 , n
         wa1(j) = -wa1(j)
         wa2(j) = x(j) + wa1(j)
         wa3(j) = diag(j)*wa1(j)
      enddo
      pnorm = enorm(n,wa3)
!
!           on the first iteration, adjust the initial step bound.
!
      if ( iter==1 ) delta = min(delta,pnorm)
!
!           evaluate the function at x + p and calculate its norm.
!
      iflag = 1
      call fcn(n,wa2,wa4,fjac,ldfjac,iflag)
      nfev = nfev + 1
      if ( iflag>=0 ) then
         fnorm1 = enorm(n,wa4)
!
!           compute the scaled actual reduction.
!
         actred = -one
         if ( fnorm1<fnorm ) actred = one - (fnorm1/fnorm)**2
!
!           compute the scaled predicted reduction.
!
         l = 1
         do i = 1 , n
            sum = zero
            do j = i , n
               sum = sum + r(l)*wa1(j)
               l = l + 1
            enddo
            wa3(i) = qtf(i) + sum
         enddo
         temp = enorm(n,wa3)
         prered = zero
         if ( temp<fnorm ) prered = one - (temp/fnorm)**2
!
!           compute the ratio of the actual to the predicted
!           reduction.
!
         ratio = zero
         if ( prered>zero ) ratio = actred/prered
!
!           update the step bound.
!
         if ( ratio>=p1 ) then
            ncfail = 0
            ncsuc = ncsuc + 1
            if ( ratio>=p5 .or. ncsuc>1 ) delta = max(delta,pnorm/p5)
            if ( abs(ratio-one)<=p1 ) delta = pnorm/p5
         else
            ncsuc = 0
            ncfail = ncfail + 1
            delta = p5*delta
         endif
!
!           test for successful iteration.
!
         if ( ratio>=p0001 ) then
!
!           successful iteration. update x, fvec, and their norms.
!
            do j = 1 , n
               x(j) = wa2(j)
               wa2(j) = diag(j)*x(j)
               fvec(j) = wa4(j)
            enddo
            xnorm = enorm(n,wa2)
            fnorm = fnorm1
            iter = iter + 1
         endif
!
!           determine the progress of the iteration.
!
         nslow1 = nslow1 + 1
         if ( actred>=p001 ) nslow1 = 0
         if ( jeval ) nslow2 = nslow2 + 1
         if ( actred>=p1 ) nslow2 = 0
!
!           test for convergence.
!
         if ( delta<=xtol*xnorm .or. fnorm==zero ) info = 1
         if ( info==0 ) then
!
!           tests for termination and stringent tolerances.
!
            if ( nfev>=maxfev ) info = 2
            if ( p1*max(p1*delta,pnorm)<=epsmch*xnorm ) info = 3
            if ( nslow2==5 ) info = 4
            if ( nslow1==10 ) info = 5
            if ( info==0 ) then
!
!           criterion for recalculating jacobian.
!
               if ( ncfail==2 ) goto 100
!
!           calculate the rank one modification to the jacobian
!           and update qtf if necessary.
!
               do j = 1 , n
                  sum = zero
                  do i = 1 , n
                     sum = sum + fjac(i,j)*wa4(i)
                  enddo
                  wa2(j) = (sum-wa3(j))/pnorm
                  wa1(j) = diag(j)*((diag(j)*wa1(j))/pnorm)
                  if ( ratio>=p0001 ) qtf(j) = sum
               enddo
!
!           compute the qr factorization of the updated jacobian.
!
               call r1updt(n,n,r,lr,wa1,wa2,wa3,sing)
               call r1mpyq(n,n,fjac,ldfjac,wa2,wa3)
               call r1mpyq(1,n,qtf,1,wa2,wa3)
!
!           end of the inner loop.
!
               jeval = .false.
!
!        end of the outer loop.
!
               goto 200
            endif
         endif
      endif
!
!     termination, either normal or user imposed.
!
 300  if ( iflag<0 ) info = iflag
      iflag = 0
      if ( nprint>0 ) call fcn(n,x,fvec,fjac,ldfjac,iflag)

      end subroutine hybrj
!*****************************************************************************************

!*****************************************************************************************
!>
!     the purpose of hybrj1 is to find a zero of a system of
!     n nonlinear functions in n variables by a modification
!     of the powell hybrid method. this is done by using the
!     more general nonlinear equation solver hybrj. the user
!     must provide a subroutine which calculates the functions
!     and the jacobian.
!
!     the subroutine statement is
!
!       subroutine hybrj1(fcn,n,x,fvec,fjac,ldfjac,tol,info,wa,lwa)
!
!     where
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions and the jacobian. fcn must
!         be declared in an external statement in the user
!         calling program, and should be written as follows.
!
!         subroutine fcn(n,x,fvec,fjac,ldfjac,iflag)
!         integer n,ldfjac,iflag
!         real(wp) x(n),fvec(n),fjac(ldfjac,n)
!         ----------
!         if iflag = 1 calculate the functions at x and
!         return this vector in fvec. do not alter fjac.
!         if iflag = 2 calculate the jacobian at x and
!         return this matrix in fjac. do not alter fvec.
!         ---------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of hybrj1.
!         in this case set iflag to a negative integer.
!
!       n is a positive integer input variable set to the number
!         of functions and variables.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length n which contains
!         the functions evaluated at the output x.
!
!       fjac is an output n by n array which contains the
!         orthogonal matrix q produced by the qr factorization
!         of the final approximate jacobian.
!
!       ldfjac is a positive integer input variable not less than n
!         which specifies the leading dimension of the array fjac.
!
!       tol is a nonnegative input variable. termination occurs
!         when the algorithm estimates that the relative error
!         between x and the solution is at most tol.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0   improper input parameters.
!
!         info = 1   algorithm estimates that the relative error
!                    between x and the solution is at most tol.
!
!         info = 2   number of calls to fcn with iflag = 1 has
!                    reached 100*(n+1).
!
!         info = 3   tol is too small. no further improvement in
!                    the approximate solution x is possible.
!
!         info = 4   iteration is not making good progress.
!
!       wa is a work array of length lwa.
!
!       lwa is a positive integer input variable not less than
!         (n*(n+13))/2.

    subroutine hybrj1(fcn,n,x,fvec,fjac,ldfjac,tol,info,wa,lwa)

    implicit none

    procedure(fcn_hybrj) :: fcn
    integer n , ldfjac , info , lwa
    real(wp) tol
    real(wp) x(n) , fvec(n) , fjac(ldfjac,n) , wa(lwa)

      integer j , lr , maxfev , mode , nfev , njev , nprint
      real(wp) xtol

      real(wp),parameter :: factor = 100.0_wp

      info = 0

      ! check the input parameters for errors.

      if ( n>0 .and. ldfjac>=n .and. tol>=zero .and. lwa>=(n*(n+13))/2 ) then

         ! call hybrj.

         maxfev = 100*(n+1)
         xtol = tol
         mode = 2
         do j = 1 , n
            wa(j) = one
         enddo
         nprint = 0
         lr = (n*(n+1))/2
         call hybrj(fcn,n,x,fvec,fjac,ldfjac,xtol,maxfev,wa(1),mode,    &
                    factor,nprint,info,nfev,njev,wa(6*n+1),lr,wa(n+1),  &
                    wa(2*n+1),wa(3*n+1),wa(4*n+1),wa(5*n+1))
         if ( info==5 ) info = 4

      endif

      end subroutine hybrj1
!*****************************************************************************************

!*****************************************************************************************
!>
!     this subroutine proceeds from the computed qr factorization of
!     an m by n matrix a to accumulate the m by m orthogonal matrix
!     q from its factored form.
!
!     the subroutine statement is
!
!       subroutine qform(m,n,q,ldq,wa)
!
!     where
!
!       m is a positive integer input variable set to the number
!         of rows of a and the order of q.
!
!       n is a positive integer input variable set to the number
!         of columns of a.
!
!       q is an m by m array. on input the full lower trapezoid in
!         the first min(m,n) columns of q contains the factored form.
!         on output q has been accumulated into a square matrix.
!
!       ldq is a positive integer input variable not less than m
!         which specifies the leading dimension of the array q.
!
!       wa is a work array of length m.

    subroutine qform(m,n,q,ldq,wa)

    implicit none

    integer m , n , ldq
    real(wp) q(ldq,m) , wa(m)

      integer i , j , jm1 , k , l , minmn , np1
      real(wp) sum , temp
!
!     zero out upper triangle of q in the first min(m,n) columns.
!
      minmn = min(m,n)
      if ( minmn>=2 ) then
         do j = 2 , minmn
            jm1 = j - 1
            do i = 1 , jm1
               q(i,j) = zero
            enddo
         enddo
      endif
!
!     initialize remaining columns to those of the identity matrix.
!
      np1 = n + 1
      if ( m>=np1 ) then
         do j = np1 , m
            do i = 1 , m
               q(i,j) = zero
            enddo
            q(j,j) = one
         enddo
      endif
!
!     accumulate q from its factored form.
!
      do l = 1 , minmn
         k = minmn - l + 1
         do i = k , m
            wa(i) = q(i,k)
            q(i,k) = zero
         enddo
         q(k,k) = one
         if ( wa(k)/=zero ) then
            do j = k , m
               sum = zero
               do i = k , m
                  sum = sum + q(i,j)*wa(i)
               enddo
               temp = sum/wa(k)
               do i = k , m
                  q(i,j) = q(i,j) - temp*wa(i)
               enddo
            enddo
         endif
      enddo

      end subroutine qform
!*****************************************************************************************

!*****************************************************************************************
!>
!     this subroutine uses householder transformations with column
!     pivoting (optional) to compute a qr factorization of the
!     m by n matrix a. that is, qrfac determines an orthogonal
!     matrix q, a permutation matrix p, and an upper trapezoidal
!     matrix r with diagonal elements of nonincreasing magnitude,
!     such that a*p = q*r. the householder transformation for
!     column k, k = 1,2,...,min(m,n), is of the form
!
!                           t
!           i - (1/u(k))*u*u
!
!     where u has zeros in the first k-1 positions. the form of
!     this transformation and the method of pivoting first
!     appeared in the corresponding linpack subroutine.
!
!     the subroutine statement is
!
!       subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
!
!     where
!
!       m is a positive integer input variable set to the number
!         of rows of a.
!
!       n is a positive integer input variable set to the number
!         of columns of a.
!
!       a is an m by n array. on input a contains the matrix for
!         which the qr factorization is to be computed. on output
!         the strict upper trapezoidal part of a contains the strict
!         upper trapezoidal part of r, and the lower trapezoidal
!         part of a contains a factored form of q (the non-trivial
!         elements of the u vectors described above).
!
!       lda is a positive integer input variable not less than m
!         which specifies the leading dimension of the array a.
!
!       pivot is a logical input variable. if pivot is set true,
!         then column pivoting is enforced. if pivot is set false,
!         then no column pivoting is done.
!
!       ipvt is an integer output array of length lipvt. ipvt
!         defines the permutation matrix p such that a*p = q*r.
!         column j of p is column ipvt(j) of the identity matrix.
!         if pivot is false, ipvt is not referenced.
!
!       lipvt is a positive integer input variable. if pivot is false,
!         then lipvt may be as small as 1. if pivot is true, then
!         lipvt must be at least n.
!
!       rdiag is an output array of length n which contains the
!         diagonal elements of r.
!
!       acnorm is an output array of length n which contains the
!         norms of the corresponding columns of the input matrix a.
!         if this information is not needed, then acnorm can coincide
!         with rdiag.
!
!       wa is a work array of length n. if pivot is false, then wa
!         can coincide with rdiag.

    subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)

    implicit none

    integer m , n , lda , lipvt
    integer ipvt(lipvt)
    logical pivot
    real(wp) a(lda,n) , rdiag(n) , acnorm(n) , wa(n)

      integer i , j , jp1 , k , kmax , minmn
      real(wp) ajnorm , epsmch , sum , temp

      real(wp),parameter :: p05 = 5.0e-2_wp

      epsmch = dpmpar(1)  ! the machine precision
!
!     compute the initial column norms and initialize several arrays.
!
      do j = 1 , n
         acnorm(j) = enorm(m,a(1,j))
         rdiag(j) = acnorm(j)
         wa(j) = rdiag(j)
         if ( pivot ) ipvt(j) = j
      enddo
!
!     reduce a to r with householder transformations.
!
      minmn = min(m,n)
      do j = 1 , minmn
         if ( pivot ) then
!
!        bring the column of largest norm into the pivot position.
!
            kmax = j
            do k = j , n
               if ( rdiag(k)>rdiag(kmax) ) kmax = k
            enddo
            if ( kmax/=j ) then
               do i = 1 , m
                  temp = a(i,j)
                  a(i,j) = a(i,kmax)
                  a(i,kmax) = temp
               enddo
               rdiag(kmax) = rdiag(j)
               wa(kmax) = wa(j)
               k = ipvt(j)
               ipvt(j) = ipvt(kmax)
               ipvt(kmax) = k
            endif
         endif
!
!        compute the householder transformation to reduce the
!        j-th column of a to a multiple of the j-th unit vector.
!
         ajnorm = enorm(m-j+1,a(j,j))
         if ( ajnorm/=zero ) then
            if ( a(j,j)<zero ) ajnorm = -ajnorm
            do i = j , m
               a(i,j) = a(i,j)/ajnorm
            enddo
            a(j,j) = a(j,j) + one
!
!        apply the transformation to the remaining columns
!        and update the norms.
!
            jp1 = j + 1
            if ( n>=jp1 ) then
               do k = jp1 , n
                  sum = zero
                  do i = j , m
                     sum = sum + a(i,j)*a(i,k)
                  enddo
                  temp = sum/a(j,j)
                  do i = j , m
                     a(i,k) = a(i,k) - temp*a(i,j)
                  enddo
                  if ( .not.(.not.pivot .or. rdiag(k)==zero) ) then
                     temp = a(j,k)/rdiag(k)
                     rdiag(k) = rdiag(k)*sqrt(max(zero,one-temp**2))
                     if ( p05*(rdiag(k)/wa(k))**2<=epsmch ) then
                        rdiag(k) = enorm(m-j,a(jp1,k))
                        wa(k) = rdiag(k)
                     endif
                  endif
               enddo
            endif
         endif
         rdiag(j) = -ajnorm
      enddo

      end subroutine qrfac
!*****************************************************************************************

!*****************************************************************************************
!>
!     given an m by n matrix a, this subroutine computes a*q where
!     q is the product of 2*(n - 1) transformations
!
!           gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
!
!     and gv(i), gw(i) are givens rotations in the (i,n) plane which
!     eliminate elements in the i-th and n-th planes, respectively.
!     q itself is not given, rather the information to recover the
!     gv, gw rotations is supplied.
!
!     the subroutine statement is
!
!       subroutine r1mpyq(m,n,a,lda,v,w)
!
!     where
!
!       m is a positive integer input variable set to the number
!         of rows of a.
!
!       n is a positive integer input variable set to the number
!         of columns of a.
!
!       a is an m by n array. on input a must contain the matrix
!         to be postmultiplied by the orthogonal matrix q
!         described above. on output a*q has replaced a.
!
!       lda is a positive integer input variable not less than m
!         which specifies the leading dimension of the array a.
!
!       v is an input array of length n. v(i) must contain the
!         information necessary to recover the givens rotation gv(i)
!         described above.
!
!       w is an input array of length n. w(i) must contain the
!         information necessary to recover the givens rotation gw(i)
!         described above.

    subroutine r1mpyq(m,n,a,lda,v,w)

    implicit none

    integer m , n , lda
    real(wp) a(lda,n) , v(n) , w(n)

      integer i , j , nmj , nm1
      real(wp) cos , sin , temp
!
!     apply the first set of givens rotations to a.
!
      nm1 = n - 1
      if ( nm1>=1 ) then
         do nmj = 1 , nm1
            j = n - nmj
            if ( abs(v(j))>one ) cos = one/v(j)
            if ( abs(v(j))>one ) sin = sqrt(one-cos**2)
            if ( abs(v(j))<=one ) sin = v(j)
            if ( abs(v(j))<=one ) cos = sqrt(one-sin**2)
            do i = 1 , m
               temp = cos*a(i,j) - sin*a(i,n)
               a(i,n) = sin*a(i,j) + cos*a(i,n)
               a(i,j) = temp
            enddo
         enddo
!
!     apply the second set of givens rotations to a.
!
         do j = 1 , nm1
            if ( abs(w(j))>one ) cos = one/w(j)
            if ( abs(w(j))>one ) sin = sqrt(one-cos**2)
            if ( abs(w(j))<=one ) sin = w(j)
            if ( abs(w(j))<=one ) cos = sqrt(one-sin**2)
            do i = 1 , m
               temp = cos*a(i,j) + sin*a(i,n)
               a(i,n) = -sin*a(i,j) + cos*a(i,n)
               a(i,j) = temp
            enddo
         enddo
      endif

      end subroutine r1mpyq
!*****************************************************************************************

!*****************************************************************************************
!>
!     given an m by n lower trapezoidal matrix s, an m-vector u,
!     and an n-vector v, the problem is to determine an
!     orthogonal matrix q such that
!
!                   t
!           (s + u*v )*q
!
!     is again lower trapezoidal.
!
!     this subroutine determines q as the product of 2*(n - 1)
!     transformations
!
!           gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
!
!     where gv(i), gw(i) are givens rotations in the (i,n) plane
!     which eliminate elements in the i-th and n-th planes,
!     respectively. q itself is not accumulated, rather the
!     information to recover the gv, gw rotations is returned.
!
!     the subroutine statement is
!
!       subroutine r1updt(m,n,s,ls,u,v,w,sing)
!
!     where
!
!       m is a positive integer input variable set to the number
!         of rows of s.
!
!       n is a positive integer input variable set to the number
!         of columns of s. n must not exceed m.
!
!       s is an array of length ls. on input s must contain the lower
!         trapezoidal matrix s stored by columns. on output s contains
!         the lower trapezoidal matrix produced as described above.
!
!       ls is a positive integer input variable not less than
!         (n*(2*m-n+1))/2.
!
!       u is an input array of length m which must contain the
!         vector u.
!
!       v is an array of length n. on input v must contain the vector
!         v. on output v(i) contains the information necessary to
!         recover the givens rotation gv(i) described above.
!
!       w is an output array of length m. w(i) contains information
!         necessary to recover the givens rotation gw(i) described
!         above.
!
!       sing is a logical output variable. sing is set true if any
!         of the diagonal elements of the output s are zero. otherwise
!         sing is set false.

    subroutine r1updt(m,n,s,ls,u,v,w,sing)

    implicit none

    integer m , n , ls
    logical sing
    real(wp) s(ls) , u(m) , v(n) , w(m)

      integer i , j , jj , l , nmj , nm1
      real(wp) cos , cotan , giant , sin ,     &
                     & tan , tau , temp

      real(wp),parameter :: p5  = 5.0e-1_wp
      real(wp),parameter :: p25 = 2.5e-1_wp
!
!     giant is the largest magnitude.
!
      giant = dpmpar(3)
!
!     initialize the diagonal element pointer.
!
      jj = (n*(2*m-n+1))/2 - (m-n)
!
!     move the nontrivial part of the last column of s into w.
!
      l = jj
      do i = n , m
         w(i) = s(l)
         l = l + 1
      enddo
!
!     rotate the vector v into a multiple of the n-th unit vector
!     in such a way that a spike is introduced into w.
!
      nm1 = n - 1
      if ( nm1>=1 ) then
         do nmj = 1 , nm1
            j = n - nmj
            jj = jj - (m-j+1)
            w(j) = zero
            if ( v(j)/=zero ) then
!
!        determine a givens rotation which eliminates the
!        j-th element of v.
!
               if ( abs(v(n))>=abs(v(j)) ) then
                  tan = v(j)/v(n)
                  cos = p5/sqrt(p25+p25*tan**2)
                  sin = cos*tan
                  tau = sin
               else
                  cotan = v(n)/v(j)
                  sin = p5/sqrt(p25+p25*cotan**2)
                  cos = sin*cotan
                  tau = one
                  if ( abs(cos)*giant>one ) tau = one/cos
               endif
!
!        apply the transformation to v and store the information
!        necessary to recover the givens rotation.
!
               v(n) = sin*v(j) + cos*v(n)
               v(j) = tau
!
!        apply the transformation to s and extend the spike in w.
!
               l = jj
               do i = j , m
                  temp = cos*s(l) - sin*w(i)
                  w(i) = sin*s(l) + cos*w(i)
                  s(l) = temp
                  l = l + 1
               enddo
            endif
         enddo
      endif
!
!     add the spike from the rank 1 update to w.
!
      do i = 1 , m
         w(i) = w(i) + v(n)*u(i)
      enddo
!
!     eliminate the spike.
!
      sing = .false.
      if ( nm1>=1 ) then
         do j = 1 , nm1
            if ( w(j)/=zero ) then
!
!        determine a givens rotation which eliminates the
!        j-th element of the spike.
!
               if ( abs(s(jj))>=abs(w(j)) ) then
                  tan = w(j)/s(jj)
                  cos = p5/sqrt(p25+p25*tan**2)
                  sin = cos*tan
                  tau = sin
               else
                  cotan = s(jj)/w(j)
                  sin = p5/sqrt(p25+p25*cotan**2)
                  cos = sin*cotan
                  tau = one
                  if ( abs(cos)*giant>one ) tau = one/cos
               endif
!
!        apply the transformation to s and reduce the spike in w.
!
               l = jj
               do i = j , m
                  temp = cos*s(l) + sin*w(i)
                  w(i) = -sin*s(l) + cos*w(i)
                  s(l) = temp
                  l = l + 1
               enddo
!
!        store the information necessary to recover the
!        givens rotation.
!
               w(j) = tau
            endif
!
!        test for zero diagonal elements in the output s.
!
            if ( s(jj)==zero ) sing = .true.
            jj = jj + (m-j+1)
         enddo
      endif
!
!     move w back into the last column of the output s.
!
      l = jj
      do i = n , m
         s(l) = w(i)
         l = l + 1
      enddo
      if ( s(jj)==zero ) sing = .true.

      end subroutine r1updt
!*****************************************************************************************

!*****************************************************************************************
    end module minpack_module
!*****************************************************************************************
