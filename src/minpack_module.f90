!*****************************************************************************************
    module minpack_module
!*****************************************************************************************
!****h* FAT/minpack_module
!
!  NAME
!    minpack_module
!
!  DESCRIPTION
!    Minpack includes software for solving nonlinear equations and
!    nonlinear least squares problems.  Five algorithmic paths each include
!    a core subroutine and an easy-to-use driver.  The algorithms proceed
!    either from an analytic specification of the Jacobian matrix or
!    directly from the problem functions.  The paths include facilities for
!    systems of equations with a banded Jacobian matrix, for least squares
!    problems with a large amount of data, and for checking the consistency
!    of the Jacobian matrix with the functions.
!
!*****************************************************************************************    
    
    implicit none
    
    contains
!*****************************************************************************************    

      subroutine chkder(m,n,x,fvec,fjac,ldfjac,xp,fvecp,mode,err)
      integer m,n,ldfjac,mode
      double precision x(n),fvec(m),fjac(ldfjac,n),xp(n),fvecp(m),&
                       err(m)
!     **********
!
!     subroutine chkder
!
!     this subroutine checks the gradients of m nonlinear functions
!     in n variables, evaluated at a point x, for consistency with
!     the functions themselves. the user must call chkder twice,
!     first with mode = 1 and then with mode = 2.
!
!     mode = 1. on input, x must contain the point of evaluation.
!               on output, xp is set to a neighboring point.
!
!     mode = 2. on input, fvec must contain the functions and the
!                         rows of fjac must contain the gradients
!                         of the respective functions each evaluated
!                         at x, and fvecp must contain the functions
!                         evaluated at xp.
!               on output, err contains measures of correctness of
!                          the respective gradients.
!
!     the subroutine does not perform reliably if cancellation or
!     rounding errors cause a severe loss of significance in the
!     evaluation of a function. therefore, none of the components
!     of x should be unusually small (in particular, zero) or any
!     other value which may cause loss of significance.
!
!     the subroutine statement is
!
!       subroutine chkder(m,n,x,fvec,fjac,ldfjac,xp,fvecp,mode,err)
!
!     where
!
!       m is a positive integer input variable set to the number
!         of functions.
!
!       n is a positive integer input variable set to the number
!         of variables.
!
!       x is an input array of length n.
!
!       fvec is an array of length m. on input when mode = 2,
!         fvec must contain the functions evaluated at x.
!
!       fjac is an m by n array. on input when mode = 2,
!         the rows of fjac must contain the gradients of
!         the respective functions evaluated at x.
!
!       ldfjac is a positive integer input parameter not less than m
!         which specifies the leading dimension of the array fjac.
!
!       xp is an array of length n. on output when mode = 1,
!         xp is set to a neighboring point of x.
!
!       fvecp is an array of length m. on input when mode = 2,
!         fvecp must contain the functions evaluated at xp.
!
!       mode is an integer input variable set to 1 on the first call
!         and 2 on the second. other values of mode are equivalent
!         to mode = 1.
!
!       err is an array of length m. on output when mode = 2,
!         err contains measures of correctness of the respective
!         gradients. if there is no severe loss of significance,
!         then if err(i) is 1.0 the i-th gradient is correct,
!         while if err(i) is 0.0 the i-th gradient is incorrect.
!         for values of err between 0.0 and 1.0, the categorization
!         is less certain. in general, a value of err(i) greater
!         than 0.5 indicates that the i-th gradient is probably
!         correct, while a value of err(i) less than 0.5 indicates
!         that the i-th gradient is probably incorrect.
!
!     subprograms called
!
!       minpack supplied ... dpmpar
!
!       fortran supplied ... dabs,dlog10,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      integer i,j
      double precision eps,epsf,epslog,epsmch,factor,one,temp,zero
      data factor,one,zero /1.0d2,1.0d0,0.0d0/
!
!     epsmch is the machine precision.
!
      epsmch = dpmpar(1)
!
      eps = dsqrt(epsmch)
!
      if (mode == 2) go to 20
!
!        mode = 1.
!
         do 10 j = 1, n
            temp = eps*dabs(x(j))
            if (temp == zero) temp = eps
            xp(j) = x(j) + temp
   10       continue
         go to 70
   20 continue
!
!        mode = 2.
!
         epsf = factor*epsmch
         epslog = dlog10(eps)
         do 30 i = 1, m
            err(i) = zero
   30       continue
         do 50 j = 1, n
            temp = dabs(x(j))
            if (temp == zero) temp = one
            do 40 i = 1, m
               err(i) = err(i) + temp*fjac(i,j)
   40          continue
   50       continue
         do 60 i = 1, m
            temp = one
            if (fvec(i) /= zero .and. fvecp(i) /= zero&
                .and. dabs(fvecp(i)-fvec(i)) >= epsf*dabs(fvec(i)))&
               temp = eps*dabs((fvecp(i)-fvec(i))/eps-err(i))&
                      /(dabs(fvec(i)) + dabs(fvecp(i)))
            err(i) = one
            if (temp > epsmch .and. temp < eps)&
               err(i) = (dlog10(temp) - epslog)/epslog
            if (temp >= eps) err(i) = zero
   60       continue
   70 continue
!
      return
!
!     last card of subroutine chkder.
!
      end subroutine chkder
      
      subroutine dogleg(n,r,lr,diag,qtb,delta,x,wa1,wa2)
      integer n,lr
      double precision delta
      double precision r(lr),diag(n),qtb(n),x(n),wa1(n),wa2(n)
!     **********
!
!     subroutine dogleg
!
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
!
!     subprograms called
!
!       minpack-supplied ... dpmpar,enorm
!
!       fortran-supplied ... dabs,dmax1,dmin1,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      integer i,j,jj,jp1,k,l
      double precision alpha,bnorm,epsmch,gnorm,one,qnorm,sgnorm,sum,&
                       temp,zero
      data one,zero /1.0d0,0.0d0/
!
!     epsmch is the machine precision.
!
      epsmch = dpmpar(1)
!
!     first, calculate the gauss-newton direction.
!
      jj = (n*(n + 1))/2 + 1
      do 50 k = 1, n
         j = n - k + 1
         jp1 = j + 1
         jj = jj - k
         l = jj + 1
         sum = zero
         if (n < jp1) go to 20
         do 10 i = jp1, n
            sum = sum + r(l)*x(i)
            l = l + 1
   10       continue
   20    continue
         temp = r(jj)
         if (temp /= zero) go to 40
         l = j
         do 30 i = 1, j
            temp = dmax1(temp,dabs(r(l)))
            l = l + n - i
   30       continue
         temp = epsmch*temp
         if (temp == zero) temp = epsmch
   40    continue
         x(j) = (qtb(j) - sum)/temp
   50    continue
!
!     test whether the gauss-newton direction is acceptable.
!
      do 60 j = 1, n
         wa1(j) = zero
         wa2(j) = diag(j)*x(j)
   60    continue
      qnorm = enorm(n,wa2)
      if (qnorm <= delta) go to 140
!
!     the gauss-newton direction is not acceptable.
!     next, calculate the scaled gradient direction.
!
      l = 1
      do 80 j = 1, n
         temp = qtb(j)
         do 70 i = j, n
            wa1(i) = wa1(i) + r(l)*temp
            l = l + 1
   70       continue
         wa1(j) = wa1(j)/diag(j)
   80    continue
!
!     calculate the norm of the scaled gradient and test for
!     the special case in which the scaled gradient is zero.
!
      gnorm = enorm(n,wa1)
      sgnorm = zero
      alpha = delta/qnorm
      if (gnorm == zero) go to 120
!
!     calculate the point along the scaled gradient
!     at which the quadratic is minimized.
!
      do 90 j = 1, n
         wa1(j) = (wa1(j)/gnorm)/diag(j)
   90    continue
      l = 1
      do 110 j = 1, n
         sum = zero
         do 100 i = j, n
            sum = sum + r(l)*wa1(i)
            l = l + 1
  100       continue
         wa2(j) = sum
  110    continue
      temp = enorm(n,wa2)
      sgnorm = (gnorm/temp)/temp
!
!     test whether the scaled gradient direction is acceptable.
!
      alpha = zero
      if (sgnorm >= delta) go to 120
!
!     the scaled gradient direction is not acceptable.
!     finally, calculate the point along the dogleg
!     at which the quadratic is minimized.
!
      bnorm = enorm(n,qtb)
      temp = (bnorm/gnorm)*(bnorm/qnorm)*(sgnorm/delta)
      temp = temp - (delta/qnorm)*(sgnorm/delta)**2&
             + dsqrt((temp-(delta/qnorm))**2&
                     +(one-(delta/qnorm)**2)*(one-(sgnorm/delta)**2))
      alpha = ((delta/qnorm)*(one - (sgnorm/delta)**2))/temp
  120 continue
!
!     form appropriate convex combination of the gauss-newton
!     direction and the scaled gradient direction.
!
      temp = (one - alpha)*dmin1(sgnorm,delta)
      do 130 j = 1, n
         x(j) = temp*wa1(j) + alpha*x(j)
  130    continue
  140 continue
      return
!
!     last card of subroutine dogleg.
!
      end subroutine dogleg
      
    double precision function dpmpar(i)
    !
    !  Replacement for the original Minpack routine.
    !
    implicit none
    
    integer,intent(in) :: i
    
    double precision,dimension(3),parameter :: dmach = [epsilon(1.0d0),&
                                                        tiny(1.0d0),&
                                                        huge(1.0d0)]

    dpmpar = dmach(i)

    end function dpmpar
      
      double precision function enorm(n,x)
      integer n
      double precision x(n)
!     **********
!
!     function enorm
!
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
!
!     the function statement is
!
!       double precision function enorm(n,x)
!
!     where
!
!       n is a positive integer input variable.
!
!       x is an input array of length n.
!
!     subprograms called
!
!       fortran-supplied ... dabs,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      integer i
      double precision agiant,floatn,one,rdwarf,rgiant,s1,s2,s3,xabs,&
                       x1max,x3max,zero
      data one,zero,rdwarf,rgiant /1.0d0,0.0d0,3.834d-20,1.304d19/
      s1 = zero
      s2 = zero
      s3 = zero
      x1max = zero
      x3max = zero
      floatn = n
      agiant = rgiant/floatn
      do 90 i = 1, n
         xabs = dabs(x(i))
         if (xabs > rdwarf .and. xabs < agiant) go to 70
            if (xabs <= rdwarf) go to 30
!
!              sum for large components.
!
               if (xabs <= x1max) go to 10
                  s1 = one + s1*(x1max/xabs)**2
                  x1max = xabs
                  go to 20
   10          continue
                  s1 = s1 + (xabs/x1max)**2
   20          continue
               go to 60
   30       continue
!
!              sum for small components.
!
               if (xabs <= x3max) go to 40
                  s3 = one + s3*(x3max/xabs)**2
                  x3max = xabs
                  go to 50
   40          continue
                  if (xabs /= zero) s3 = s3 + (xabs/x3max)**2
   50          continue
   60       continue
            go to 80
   70    continue
!
!           sum for intermediate components.
!
            s2 = s2 + xabs**2
   80    continue
   90    continue
!
!     calculation of norm.
!
      if (s1 == zero) go to 100
         enorm = x1max*dsqrt(s1+(s2/x1max)/x1max)
         go to 130
  100 continue
         if (s2 == zero) go to 110
            if (s2 >= x3max)&
               enorm = dsqrt(s2*(one+(x3max/s2)*(x3max*s3)))
            if (s2 < x3max)&
               enorm = dsqrt(x3max*((s2/x3max)+(x3max*s3)))
            go to 120
  110    continue
            enorm = x3max*dsqrt(s3)
  120    continue
  130 continue
      return
!
!     last card of function enorm.
!
      end function enorm
      
      subroutine fdjac1(fcn,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,&
                        wa1,wa2)
      integer n,ldfjac,iflag,ml,mu
      double precision epsfcn
      double precision x(n),fvec(n),fjac(ldfjac,n),wa1(n),wa2(n)
!     **********
!
!     subroutine fdjac1
!
!     this subroutine computes a forward-difference approximation
!     to the n by n jacobian matrix associated with a specified
!     problem of n functions in n variables. if the jacobian has
!     a banded form, then function evaluations are saved by only
!     approximating the nonzero terms.
!
!     the subroutine statement is
!
!       subroutine fdjac1(fcn,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,
!                         wa1,wa2)
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
!         double precision x(n),fvec(n)
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
!
!     subprograms called
!
!       minpack-supplied ... dpmpar
!
!       fortran-supplied ... dabs,dmax1,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      integer i,j,k,msum
      double precision eps,epsmch,h,temp,zero
      data zero /0.0d0/
!
!     epsmch is the machine precision.
!
      epsmch = dpmpar(1)
!
      eps = dsqrt(dmax1(epsfcn,epsmch))
      msum = ml + mu + 1
      if (msum < n) go to 40
!
!        computation of dense approximate jacobian.
!
         do 20 j = 1, n
            temp = x(j)
            h = eps*dabs(temp)
            if (h == zero) h = eps
            x(j) = temp + h
            call fcn(n,x,wa1,iflag)
            if (iflag < 0) go to 30
            x(j) = temp
            do 10 i = 1, n
               fjac(i,j) = (wa1(i) - fvec(i))/h
   10          continue
   20       continue
   30    continue
         go to 110
   40 continue
!
!        computation of banded approximate jacobian.
!
         do 90 k = 1, msum
            do 60 j = k, n, msum
               wa2(j) = x(j)
               h = eps*dabs(wa2(j))
               if (h == zero) h = eps
               x(j) = wa2(j) + h
   60          continue
            call fcn(n,x,wa1,iflag)
            if (iflag < 0) go to 100
            do 80 j = k, n, msum
               x(j) = wa2(j)
               h = eps*dabs(wa2(j))
               if (h == zero) h = eps
               do 70 i = 1, n
                  fjac(i,j) = zero
                  if (i >= j - mu .and. i <= j + ml)&
                     fjac(i,j) = (wa1(i) - fvec(i))/h
   70             continue
   80          continue
   90       continue
  100    continue
  110 continue
      return
!
!     last card of subroutine fdjac1.
!
      end subroutine fdjac1

      subroutine fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa)
      integer m,n,ldfjac,iflag
      double precision epsfcn
      double precision x(n),fvec(m),fjac(ldfjac,n),wa(m)
!     **********
!
!     subroutine fdjac2
!
!     this subroutine computes a forward-difference approximation
!     to the m by n jacobian matrix associated with a specified
!     problem of m functions in n variables.
!
!     the subroutine statement is
!
!       subroutine fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa)
!
!     where
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions. fcn must be declared
!         in an external statement in the user calling
!         program, and should be written as follows.
!
!         subroutine fcn(m,n,x,fvec,iflag)
!         integer m,n,iflag
!         double precision x(n),fvec(m)
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ----------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of fdjac2.
!         in this case set iflag to a negative integer.
!
!       m is a positive integer input variable set to the number
!         of functions.
!
!       n is a positive integer input variable set to the number
!         of variables. n must not exceed m.
!
!       x is an input array of length n.
!
!       fvec is an input array of length m which must contain the
!         functions evaluated at x.
!
!       fjac is an output m by n array which contains the
!         approximation to the jacobian matrix evaluated at x.
!
!       ldfjac is a positive integer input variable not less than m
!         which specifies the leading dimension of the array fjac.
!
!       iflag is an integer variable which can be used to terminate
!         the execution of fdjac2. see description of fcn.
!
!       epsfcn is an input variable used in determining a suitable
!         step length for the forward-difference approximation. this
!         approximation assumes that the relative errors in the
!         functions are of the order of epsfcn. if epsfcn is less
!         than the machine precision, it is assumed that the relative
!         errors in the functions are of the order of the machine
!         precision.
!
!       wa is a work array of length m.
!
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... dpmpar
!
!       fortran-supplied ... dabs,dmax1,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      integer i,j
      double precision eps,epsmch,h,temp,zero
      data zero /0.0d0/
!
!     epsmch is the machine precision.
!
      epsmch = dpmpar(1)
!
      eps = dsqrt(dmax1(epsfcn,epsmch))
      do 20 j = 1, n
         temp = x(j)
         h = eps*dabs(temp)
         if (h == zero) h = eps
         x(j) = temp + h
         call fcn(m,n,x,wa,iflag)
         if (iflag < 0) go to 30
         x(j) = temp
         do 10 i = 1, m
            fjac(i,j) = (wa(i) - fvec(i))/h
   10       continue
   20    continue
   30 continue
      return
!
!     last card of subroutine fdjac2.
!
      end subroutine fdjac2
      
      subroutine hybrd(fcn,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,diag,&
                       mode,factor,nprint,info,nfev,fjac,ldfjac,r,lr,&
                       qtf,wa1,wa2,wa3,wa4)
      integer n,maxfev,ml,mu,mode,nprint,info,nfev,ldfjac,lr
      double precision xtol,epsfcn,factor
      double precision x(n),fvec(n),diag(n),fjac(ldfjac,n),r(lr),&
                       qtf(n),wa1(n),wa2(n),wa3(n),wa4(n)
      external fcn
!     **********
!
!     subroutine hybrd
!
!     the purpose of hybrd is to find a zero of a system of
!     n nonlinear functions in n variables by a modification
!     of the powell hybrid method. the user must provide a
!     subroutine which calculates the functions. the jacobian is
!     then calculated by a forward-difference approximation.
!
!     the subroutine statement is
!
!       subroutine hybrd(fcn,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,
!                        diag,mode,factor,nprint,info,nfev,fjac,
!                        ldfjac,r,lr,qtf,wa1,wa2,wa3,wa4)
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
!         double precision x(n),fvec(n)
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ---------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of hybrd.
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
!       xtol is a nonnegative input variable. termination
!         occurs when the relative error between two consecutive
!         iterates is at most xtol.
!
!       maxfev is a positive integer input variable. termination
!         occurs when the number of calls to fcn is at least maxfev
!         by the end of an iteration.
!
!       ml is a nonnegative integer input variable which specifies
!         the number of subdiagonals within the band of the
!         jacobian matrix. if the jacobian is not banded, set
!         ml to at least n - 1.
!
!       mu is a nonnegative integer input variable which specifies
!         the number of superdiagonals within the band of the
!         jacobian matrix. if the jacobian is not banded, set
!         mu to at least n - 1.
!
!       epsfcn is an input variable used in determining a suitable
!         step length for the forward-difference approximation. this
!         approximation assumes that the relative errors in the
!         functions are of the order of epsfcn. if epsfcn is less
!         than the machine precision, it is assumed that the relative
!         errors in the functions are of the order of the machine
!         precision.
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
!         for printing. if nprint is not positive, no special calls
!         of fcn with iflag = 0 are made.
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
!         info = 2   number of calls to fcn has reached or exceeded
!                    maxfev.
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
!         calls to fcn.
!
!       fjac is an output n by n array which contains the
!         orthogonal matrix q produced by the qr factorization
!         of the final approximate jacobian.
!
!       ldfjac is a positive integer input variable not less than n
!         which specifies the leading dimension of the array fjac.
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
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... dogleg,dpmpar,enorm,fdjac1,
!                            qform,qrfac,r1mpyq,r1updt
!
!       fortran-supplied ... dabs,dmax1,dmin1,min0,mod
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      integer i,iflag,iter,j,jm1,l,msum,ncfail,ncsuc,nslow1,nslow2
      integer iwa(1)
      logical jeval,sing
      double precision actred,delta,epsmch,fnorm,fnorm1,one,pnorm,&
                       prered,p1,p5,p001,p0001,ratio,sum,temp,xnorm,&
                       zero
      data one,p1,p5,p001,p0001,zero&
           /1.0d0,1.0d-1,5.0d-1,1.0d-3,1.0d-4,0.0d0/
!
!     epsmch is the machine precision.
!
      epsmch = dpmpar(1)
!
      info = 0
      iflag = 0
      nfev = 0
!
!     check the input parameters for errors.
!
      if (n <= 0 .or. xtol < zero .or. maxfev <= 0&
          .or. ml < 0 .or. mu < 0 .or. factor <= zero&
          .or. ldfjac < n .or. lr < (n*(n + 1))/2) go to 300
      if (mode /= 2) go to 20
      do 10 j = 1, n
         if (diag(j) <= zero) go to 300
   10    continue
   20 continue
!
!     evaluate the function at the starting point
!     and calculate its norm.
!
      iflag = 1
      call fcn(n,x,fvec,iflag)
      nfev = 1
      if (iflag < 0) go to 300
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
   30 continue
         jeval = .true.
!
!        calculate the jacobian matrix.
!
         iflag = 2
         call fdjac1(fcn,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,wa1,&
                     wa2)
         nfev = nfev + msum
         if (iflag < 0) go to 300
!
!        compute the qr factorization of the jacobian.
!
         call qrfac(n,n,fjac,ldfjac,.false.,iwa,1,wa1,wa2,wa3)
!
!        on the first iteration and if mode is 1, scale according
!        to the norms of the columns of the initial jacobian.
!
         if (iter /= 1) go to 70
         if (mode == 2) go to 50
         do 40 j = 1, n
            diag(j) = wa2(j)
            if (wa2(j) == zero) diag(j) = one
   40       continue
   50    continue
!
!        on the first iteration, calculate the norm of the scaled x
!        and initialize the step bound delta.
!
         do 60 j = 1, n
            wa3(j) = diag(j)*x(j)
   60       continue
         xnorm = enorm(n,wa3)
         delta = factor*xnorm
         if (delta == zero) delta = factor
   70    continue
!
!        form (q transpose)*fvec and store in qtf.
!
         do 80 i = 1, n
            qtf(i) = fvec(i)
   80       continue
         do 120 j = 1, n
            if (fjac(j,j) == zero) go to 110
            sum = zero
            do 90 i = j, n
               sum = sum + fjac(i,j)*qtf(i)
   90          continue
            temp = -sum/fjac(j,j)
            do 100 i = j, n
               qtf(i) = qtf(i) + fjac(i,j)*temp
  100          continue
  110       continue
  120       continue
!
!        copy the triangular factor of the qr factorization into r.
!
         sing = .false.
         do 150 j = 1, n
            l = j
            jm1 = j - 1
            if (jm1 < 1) go to 140
            do 130 i = 1, jm1
               r(l) = fjac(i,j)
               l = l + n - i
  130          continue
  140       continue
            r(l) = wa1(j)
            if (wa1(j) == zero) sing = .true.
  150       continue
!
!        accumulate the orthogonal factor in fjac.
!
         call qform(n,n,fjac,ldfjac,wa1)
!
!        rescale if necessary.
!
         if (mode == 2) go to 170
         do 160 j = 1, n
            diag(j) = dmax1(diag(j),wa2(j))
  160       continue
  170    continue
!
!        beginning of the inner loop.
!
  180    continue
!
!           if requested, call fcn to enable printing of iterates.
!
            if (nprint <= 0) go to 190
            iflag = 0
            if (mod(iter-1,nprint) == 0) call fcn(n,x,fvec,iflag)
            if (iflag < 0) go to 300
  190       continue
!
!           determine the direction p.
!
            call dogleg(n,r,lr,diag,qtf,delta,wa1,wa2,wa3)
!
!           store the direction p and x + p. calculate the norm of p.
!
            do 200 j = 1, n
               wa1(j) = -wa1(j)
               wa2(j) = x(j) + wa1(j)
               wa3(j) = diag(j)*wa1(j)
  200          continue
            pnorm = enorm(n,wa3)
!
!           on the first iteration, adjust the initial step bound.
!
            if (iter == 1) delta = dmin1(delta,pnorm)
!
!           evaluate the function at x + p and calculate its norm.
!
            iflag = 1
            call fcn(n,wa2,wa4,iflag)
            nfev = nfev + 1
            if (iflag < 0) go to 300
            fnorm1 = enorm(n,wa4)
!
!           compute the scaled actual reduction.
!
            actred = -one
            if (fnorm1 < fnorm) actred = one - (fnorm1/fnorm)**2
!
!           compute the scaled predicted reduction.
!
            l = 1
            do 220 i = 1, n
               sum = zero
               do 210 j = i, n
                  sum = sum + r(l)*wa1(j)
                  l = l + 1
  210             continue
               wa3(i) = qtf(i) + sum
  220          continue
            temp = enorm(n,wa3)
            prered = zero
            if (temp < fnorm) prered = one - (temp/fnorm)**2
!
!           compute the ratio of the actual to the predicted
!           reduction.
!
            ratio = zero
            if (prered > zero) ratio = actred/prered
!
!           update the step bound.
!
            if (ratio >= p1) go to 230
               ncsuc = 0
               ncfail = ncfail + 1
               delta = p5*delta
               go to 240
  230       continue
               ncfail = 0
               ncsuc = ncsuc + 1
               if (ratio >= p5 .or. ncsuc > 1)&
                  delta = dmax1(delta,pnorm/p5)
               if (dabs(ratio-one) <= p1) delta = pnorm/p5
  240       continue
!
!           test for successful iteration.
!
            if (ratio < p0001) go to 260
!
!           successful iteration. update x, fvec, and their norms.
!
            do 250 j = 1, n
               x(j) = wa2(j)
               wa2(j) = diag(j)*x(j)
               fvec(j) = wa4(j)
  250          continue
            xnorm = enorm(n,wa2)
            fnorm = fnorm1
            iter = iter + 1
  260       continue
!
!           determine the progress of the iteration.
!
            nslow1 = nslow1 + 1
            if (actred >= p001) nslow1 = 0
            if (jeval) nslow2 = nslow2 + 1
            if (actred >= p1) nslow2 = 0
!
!           test for convergence.
!
            if (delta <= xtol*xnorm .or. fnorm == zero) info = 1
            if (info /= 0) go to 300
!
!           tests for termination and stringent tolerances.
!
            if (nfev >= maxfev) info = 2
            if (p1*dmax1(p1*delta,pnorm) <= epsmch*xnorm) info = 3
            if (nslow2 == 5) info = 4
            if (nslow1 == 10) info = 5
            if (info /= 0) go to 300
!
!           criterion for recalculating jacobian approximation
!           by forward differences.
!
            if (ncfail == 2) go to 290
!
!           calculate the rank one modification to the jacobian
!           and update qtf if necessary.
!
            do 280 j = 1, n
               sum = zero
               do 270 i = 1, n
                  sum = sum + fjac(i,j)*wa4(i)
  270             continue
               wa2(j) = (sum - wa3(j))/pnorm
               wa1(j) = diag(j)*((diag(j)*wa1(j))/pnorm)
               if (ratio >= p0001) qtf(j) = sum
  280          continue
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
            go to 180
  290    continue
!
!        end of the outer loop.
!
         go to 30
  300 continue
!
!     termination, either normal or user imposed.
!
      if (iflag < 0) info = iflag
      iflag = 0
      if (nprint > 0) call fcn(n,x,fvec,iflag)
      return
!
!     last card of subroutine hybrd.
!
      end subroutine hybrd
      
      subroutine hybrd1(fcn,n,x,fvec,tol,info,wa,lwa)
      integer n,info,lwa
      double precision tol
      double precision x(n),fvec(n),wa(lwa)
      external fcn
!     **********
!
!     subroutine hybrd1
!
!     the purpose of hybrd1 is to find a zero of a system of
!     n nonlinear functions in n variables by a modification
!     of the powell hybrid method. this is done by using the
!     more general nonlinear equation solver hybrd. the user
!     must provide a subroutine which calculates the functions.
!     the jacobian is then calculated by a forward-difference
!     approximation.
!
!     the subroutine statement is
!
!       subroutine hybrd1(fcn,n,x,fvec,tol,info,wa,lwa)
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
!         double precision x(n),fvec(n)
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ---------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of hybrd1.
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
!         info = 2   number of calls to fcn has reached or exceeded
!                    200*(n+1).
!
!         info = 3   tol is too small. no further improvement in
!                    the approximate solution x is possible.
!
!         info = 4   iteration is not making good progress.
!
!       wa is a work array of length lwa.
!
!       lwa is a positive integer input variable not less than
!         (n*(3*n+13))/2.
!
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... hybrd
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      integer index,j,lr,maxfev,ml,mode,mu,nfev,nprint
      double precision epsfcn,factor,one,xtol,zero
      data factor,one,zero /1.0d2,1.0d0,0.0d0/
      info = 0
!
!     check the input parameters for errors.
!
      if (n <= 0 .or. tol < zero .or. lwa < (n*(3*n + 13))/2)&
         go to 20
!
!     call hybrd.
!
      maxfev = 200*(n + 1)
      xtol = tol
      ml = n - 1
      mu = n - 1
      epsfcn = zero
      mode = 2
      do 10 j = 1, n
         wa(j) = one
   10    continue
      nprint = 0
      lr = (n*(n + 1))/2
      index = 6*n + lr
      call hybrd(fcn,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,wa(1),mode,&
                 factor,nprint,info,nfev,wa(index+1),n,wa(6*n+1),lr,&
                 wa(n+1),wa(2*n+1),wa(3*n+1),wa(4*n+1),wa(5*n+1))
      if (info == 5) info = 4
   20 continue
      return
!
!     last card of subroutine hybrd1.
!
      end subroutine hybrd1
      
      subroutine hybrj(fcn,n,x,fvec,fjac,ldfjac,xtol,maxfev,diag,mode,&
                       factor,nprint,info,nfev,njev,r,lr,qtf,wa1,wa2,&
                       wa3,wa4)
      integer n,ldfjac,maxfev,mode,nprint,info,nfev,njev,lr
      double precision xtol,factor
      double precision x(n),fvec(n),fjac(ldfjac,n),diag(n),r(lr),&
                       qtf(n),wa1(n),wa2(n),wa3(n),wa4(n)
!     **********
!
!     subroutine hybrj
!
!     the purpose of hybrj is to find a zero of a system of
!     n nonlinear functions in n variables by a modification
!     of the powell hybrid method. the user must provide a
!     subroutine which calculates the functions and the jacobian.
!
!     the subroutine statement is
!
!       subroutine hybrj(fcn,n,x,fvec,fjac,ldfjac,xtol,maxfev,diag,
!                        mode,factor,nprint,info,nfev,njev,r,lr,qtf,
!                        wa1,wa2,wa3,wa4)
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
!         double precision x(n),fvec(n),fjac(ldfjac,n)
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
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... dogleg,dpmpar,enorm,
!                            qform,qrfac,r1mpyq,r1updt
!
!       fortran-supplied ... dabs,dmax1,dmin1,mod
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      integer i,iflag,iter,j,jm1,l,ncfail,ncsuc,nslow1,nslow2
      integer iwa(1)
      logical jeval,sing
      double precision actred,delta,epsmch,fnorm,fnorm1,one,pnorm,&
                       prered,p1,p5,p001,p0001,ratio,sum,temp,xnorm,&
                       zero
      data one,p1,p5,p001,p0001,zero&
           /1.0d0,1.0d-1,5.0d-1,1.0d-3,1.0d-4,0.0d0/
!
!     epsmch is the machine precision.
!
      epsmch = dpmpar(1)
!
      info = 0
      iflag = 0
      nfev = 0
      njev = 0
!
!     check the input parameters for errors.
!
      if (n <= 0 .or. ldfjac < n .or. xtol < zero&
          .or. maxfev <= 0 .or. factor <= zero&
          .or. lr < (n*(n + 1))/2) go to 300
      if (mode /= 2) go to 20
      do 10 j = 1, n
         if (diag(j) <= zero) go to 300
   10    continue
   20 continue
!
!     evaluate the function at the starting point
!     and calculate its norm.
!
      iflag = 1
      call fcn(n,x,fvec,fjac,ldfjac,iflag)
      nfev = 1
      if (iflag < 0) go to 300
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
   30 continue
         jeval = .true.
!
!        calculate the jacobian matrix.
!
         iflag = 2
         call fcn(n,x,fvec,fjac,ldfjac,iflag)
         njev = njev + 1
         if (iflag < 0) go to 300
!
!        compute the qr factorization of the jacobian.
!
         call qrfac(n,n,fjac,ldfjac,.false.,iwa,1,wa1,wa2,wa3)
!
!        on the first iteration and if mode is 1, scale according
!        to the norms of the columns of the initial jacobian.
!
         if (iter /= 1) go to 70
         if (mode == 2) go to 50
         do 40 j = 1, n
            diag(j) = wa2(j)
            if (wa2(j) == zero) diag(j) = one
   40       continue
   50    continue
!
!        on the first iteration, calculate the norm of the scaled x
!        and initialize the step bound delta.
!
         do 60 j = 1, n
            wa3(j) = diag(j)*x(j)
   60       continue
         xnorm = enorm(n,wa3)
         delta = factor*xnorm
         if (delta == zero) delta = factor
   70    continue
!
!        form (q transpose)*fvec and store in qtf.
!
         do 80 i = 1, n
            qtf(i) = fvec(i)
   80       continue
         do 120 j = 1, n
            if (fjac(j,j) == zero) go to 110
            sum = zero
            do 90 i = j, n
               sum = sum + fjac(i,j)*qtf(i)
   90          continue
            temp = -sum/fjac(j,j)
            do 100 i = j, n
               qtf(i) = qtf(i) + fjac(i,j)*temp
  100          continue
  110       continue
  120       continue
!
!        copy the triangular factor of the qr factorization into r.
!
         sing = .false.
         do 150 j = 1, n
            l = j
            jm1 = j - 1
            if (jm1 < 1) go to 140
            do 130 i = 1, jm1
               r(l) = fjac(i,j)
               l = l + n - i
  130          continue
  140       continue
            r(l) = wa1(j)
            if (wa1(j) == zero) sing = .true.
  150       continue
!
!        accumulate the orthogonal factor in fjac.
!
         call qform(n,n,fjac,ldfjac,wa1)
!
!        rescale if necessary.
!
         if (mode == 2) go to 170
         do 160 j = 1, n
            diag(j) = dmax1(diag(j),wa2(j))
  160       continue
  170    continue
!
!        beginning of the inner loop.
!
  180    continue
!
!           if requested, call fcn to enable printing of iterates.
!
            if (nprint <= 0) go to 190
            iflag = 0
            if (mod(iter-1,nprint) == 0)&
               call fcn(n,x,fvec,fjac,ldfjac,iflag)
            if (iflag < 0) go to 300
  190       continue
!
!           determine the direction p.
!
            call dogleg(n,r,lr,diag,qtf,delta,wa1,wa2,wa3)
!
!           store the direction p and x + p. calculate the norm of p.
!
            do 200 j = 1, n
               wa1(j) = -wa1(j)
               wa2(j) = x(j) + wa1(j)
               wa3(j) = diag(j)*wa1(j)
  200          continue
            pnorm = enorm(n,wa3)
!
!           on the first iteration, adjust the initial step bound.
!
            if (iter == 1) delta = dmin1(delta,pnorm)
!
!           evaluate the function at x + p and calculate its norm.
!
            iflag = 1
            call fcn(n,wa2,wa4,fjac,ldfjac,iflag)
            nfev = nfev + 1
            if (iflag < 0) go to 300
            fnorm1 = enorm(n,wa4)
!
!           compute the scaled actual reduction.
!
            actred = -one
            if (fnorm1 < fnorm) actred = one - (fnorm1/fnorm)**2
!
!           compute the scaled predicted reduction.
!
            l = 1
            do 220 i = 1, n
               sum = zero
               do 210 j = i, n
                  sum = sum + r(l)*wa1(j)
                  l = l + 1
  210             continue
               wa3(i) = qtf(i) + sum
  220          continue
            temp = enorm(n,wa3)
            prered = zero
            if (temp < fnorm) prered = one - (temp/fnorm)**2
!
!           compute the ratio of the actual to the predicted
!           reduction.
!
            ratio = zero
            if (prered > zero) ratio = actred/prered
!
!           update the step bound.
!
            if (ratio >= p1) go to 230
               ncsuc = 0
               ncfail = ncfail + 1
               delta = p5*delta
               go to 240
  230       continue
               ncfail = 0
               ncsuc = ncsuc + 1
               if (ratio >= p5 .or. ncsuc > 1)&
                  delta = dmax1(delta,pnorm/p5)
               if (dabs(ratio-one) <= p1) delta = pnorm/p5
  240       continue
!
!           test for successful iteration.
!
            if (ratio < p0001) go to 260
!
!           successful iteration. update x, fvec, and their norms.
!
            do 250 j = 1, n
               x(j) = wa2(j)
               wa2(j) = diag(j)*x(j)
               fvec(j) = wa4(j)
  250          continue
            xnorm = enorm(n,wa2)
            fnorm = fnorm1
            iter = iter + 1
  260       continue
!
!           determine the progress of the iteration.
!
            nslow1 = nslow1 + 1
            if (actred >= p001) nslow1 = 0
            if (jeval) nslow2 = nslow2 + 1
            if (actred >= p1) nslow2 = 0
!
!           test for convergence.
!
            if (delta <= xtol*xnorm .or. fnorm == zero) info = 1
            if (info /= 0) go to 300
!
!           tests for termination and stringent tolerances.
!
            if (nfev >= maxfev) info = 2
            if (p1*dmax1(p1*delta,pnorm) <= epsmch*xnorm) info = 3
            if (nslow2 == 5) info = 4
            if (nslow1 == 10) info = 5
            if (info /= 0) go to 300
!
!           criterion for recalculating jacobian.
!
            if (ncfail == 2) go to 290
!
!           calculate the rank one modification to the jacobian
!           and update qtf if necessary.
!
            do 280 j = 1, n
               sum = zero
               do 270 i = 1, n
                  sum = sum + fjac(i,j)*wa4(i)
  270             continue
               wa2(j) = (sum - wa3(j))/pnorm
               wa1(j) = diag(j)*((diag(j)*wa1(j))/pnorm)
               if (ratio >= p0001) qtf(j) = sum
  280          continue
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
            go to 180
  290    continue
!
!        end of the outer loop.
!
         go to 30
  300 continue
!
!     termination, either normal or user imposed.
!
      if (iflag < 0) info = iflag
      iflag = 0
      if (nprint > 0) call fcn(n,x,fvec,fjac,ldfjac,iflag)
      return
!
!     last card of subroutine hybrj.
!
      end subroutine hybrj
      
      subroutine hybrj1(fcn,n,x,fvec,fjac,ldfjac,tol,info,wa,lwa)
      integer n,ldfjac,info,lwa
      double precision tol
      double precision x(n),fvec(n),fjac(ldfjac,n),wa(lwa)
      external fcn
!     **********
!
!     subroutine hybrj1
!
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
!         double precision x(n),fvec(n),fjac(ldfjac,n)
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
!
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... hybrj
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      integer j,lr,maxfev,mode,nfev,njev,nprint
      double precision factor,one,xtol,zero
      data factor,one,zero /1.0d2,1.0d0,0.0d0/
      info = 0
!
!     check the input parameters for errors.
!
      if (n <= 0 .or. ldfjac < n .or. tol < zero&
          .or. lwa < (n*(n + 13))/2) go to 20
!
!     call hybrj.
!
      maxfev = 100*(n + 1)
      xtol = tol
      mode = 2
      do 10 j = 1, n
         wa(j) = one
   10    continue
      nprint = 0
      lr = (n*(n + 1))/2
      call hybrj(fcn,n,x,fvec,fjac,ldfjac,xtol,maxfev,wa(1),mode,&
                 factor,nprint,info,nfev,njev,wa(6*n+1),lr,wa(n+1),&
                 wa(2*n+1),wa(3*n+1),wa(4*n+1),wa(5*n+1))
      if (info == 5) info = 4
   20 continue
      return
!
!     last card of subroutine hybrj1.
!
      end subroutine hybrj1
      
      subroutine lmder(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,&
                       maxfev,diag,mode,factor,nprint,info,nfev,njev,&
                       ipvt,qtf,wa1,wa2,wa3,wa4)
      integer m,n,ldfjac,maxfev,mode,nprint,info,nfev,njev
      integer ipvt(n)
      double precision ftol,xtol,gtol,factor
      double precision x(n),fvec(m),fjac(ldfjac,n),diag(n),qtf(n),&
                       wa1(n),wa2(n),wa3(n),wa4(m)
!     **********
!
!     subroutine lmder
!
!     the purpose of lmder is to minimize the sum of the squares of
!     m nonlinear functions in n variables by a modification of
!     the levenberg-marquardt algorithm. the user must provide a
!     subroutine which calculates the functions and the jacobian.
!
!     the subroutine statement is
!
!       subroutine lmder(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,
!                        maxfev,diag,mode,factor,nprint,info,nfev,
!                        njev,ipvt,qtf,wa1,wa2,wa3,wa4)
!
!     where
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions and the jacobian. fcn must
!         be declared in an external statement in the user
!         calling program, and should be written as follows.
!
!         subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
!         integer m,n,ldfjac,iflag
!         double precision x(n),fvec(m),fjac(ldfjac,n)
!         ----------
!         if iflag = 1 calculate the functions at x and
!         return this vector in fvec. do not alter fjac.
!         if iflag = 2 calculate the jacobian at x and
!         return this matrix in fjac. do not alter fvec.
!         ----------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of lmder.
!         in this case set iflag to a negative integer.
!
!       m is a positive integer input variable set to the number
!         of functions.
!
!       n is a positive integer input variable set to the number
!         of variables. n must not exceed m.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length m which contains
!         the functions evaluated at the output x.
!
!       fjac is an output m by n array. the upper n by n submatrix
!         of fjac contains an upper triangular matrix r with
!         diagonal elements of nonincreasing magnitude such that
!
!                t     t           t
!               p *(jac *jac)*p = r *r,
!
!         where p is a permutation matrix and jac is the final
!         calculated jacobian. column j of p is column ipvt(j)
!         (see below) of the identity matrix. the lower trapezoidal
!         part of fjac contains information generated during
!         the computation of r.
!
!       ldfjac is a positive integer input variable not less than m
!         which specifies the leading dimension of the array fjac.
!
!       ftol is a nonnegative input variable. termination
!         occurs when both the actual and predicted relative
!         reductions in the sum of squares are at most ftol.
!         therefore, ftol measures the relative error desired
!         in the sum of squares.
!
!       xtol is a nonnegative input variable. termination
!         occurs when the relative error between two consecutive
!         iterates is at most xtol. therefore, xtol measures the
!         relative error desired in the approximate solution.
!
!       gtol is a nonnegative input variable. termination
!         occurs when the cosine of the angle between fvec and
!         any column of the jacobian is at most gtol in absolute
!         value. therefore, gtol measures the orthogonality
!         desired between the function vector and the columns
!         of the jacobian.
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
!         interval (.1,100.).100. is a generally recommended value.
!
!       nprint is an integer input variable that enables controlled
!         printing of iterates if it is positive. in this case,
!         fcn is called with iflag = 0 at the beginning of the first
!         iteration and every nprint iterations thereafter and
!         immediately prior to return, with x, fvec, and fjac
!         available for printing. fvec and fjac should not be
!         altered. if nprint is not positive, no special calls
!         of fcn with iflag = 0 are made.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0  improper input parameters.
!
!         info = 1  both actual and predicted relative reductions
!                   in the sum of squares are at most ftol.
!
!         info = 2  relative error between two consecutive iterates
!                   is at most xtol.
!
!         info = 3  conditions for info = 1 and info = 2 both hold.
!
!         info = 4  the cosine of the angle between fvec and any
!                   column of the jacobian is at most gtol in
!                   absolute value.
!
!         info = 5  number of calls to fcn with iflag = 1 has
!                   reached maxfev.
!
!         info = 6  ftol is too small. no further reduction in
!                   the sum of squares is possible.
!
!         info = 7  xtol is too small. no further improvement in
!                   the approximate solution x is possible.
!
!         info = 8  gtol is too small. fvec is orthogonal to the
!                   columns of the jacobian to machine precision.
!
!       nfev is an integer output variable set to the number of
!         calls to fcn with iflag = 1.
!
!       njev is an integer output variable set to the number of
!         calls to fcn with iflag = 2.
!
!       ipvt is an integer output array of length n. ipvt
!         defines a permutation matrix p such that jac*p = q*r,
!         where jac is the final calculated jacobian, q is
!         orthogonal (not stored), and r is upper triangular
!         with diagonal elements of nonincreasing magnitude.
!         column j of p is column ipvt(j) of the identity matrix.
!
!       qtf is an output array of length n which contains
!         the first n elements of the vector (q transpose)*fvec.
!
!       wa1, wa2, and wa3 are work arrays of length n.
!
!       wa4 is a work array of length m.
!
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... dpmpar,enorm,lmpar,qrfac
!
!       fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      integer i,iflag,iter,j,l
      double precision actred,delta,dirder,epsmch,fnorm,fnorm1,gnorm,&
                       one,par,pnorm,prered,p1,p5,p25,p75,p0001,ratio,&
                       sum,temp,temp1,temp2,xnorm,zero
      data one,p1,p5,p25,p75,p0001,zero&
           /1.0d0,1.0d-1,5.0d-1,2.5d-1,7.5d-1,1.0d-4,0.0d0/
!
!     epsmch is the machine precision.
!
      epsmch = dpmpar(1)
!
      info = 0
      iflag = 0
      nfev = 0
      njev = 0
!
!     check the input parameters for errors.
!
      if (n <= 0 .or. m < n .or. ldfjac < m&
          .or. ftol < zero .or. xtol < zero .or. gtol < zero&
          .or. maxfev <= 0 .or. factor <= zero) go to 300
      if (mode /= 2) go to 20
      do 10 j = 1, n
         if (diag(j) <= zero) go to 300
   10    continue
   20 continue
!
!     evaluate the function at the starting point
!     and calculate its norm.
!
      iflag = 1
      call fcn(m,n,x,fvec,fjac,ldfjac,iflag)
      nfev = 1
      if (iflag < 0) go to 300
      fnorm = enorm(m,fvec)
!
!     initialize levenberg-marquardt parameter and iteration counter.
!
      par = zero
      iter = 1
!
!     beginning of the outer loop.
!
   30 continue
!
!        calculate the jacobian matrix.
!
         iflag = 2
         call fcn(m,n,x,fvec,fjac,ldfjac,iflag)
         njev = njev + 1
         if (iflag < 0) go to 300
!
!        if requested, call fcn to enable printing of iterates.
!
         if (nprint <= 0) go to 40
         iflag = 0
         if (mod(iter-1,nprint) == 0)&
            call fcn(m,n,x,fvec,fjac,ldfjac,iflag)
         if (iflag < 0) go to 300
   40    continue
!
!        compute the qr factorization of the jacobian.
!
         call qrfac(m,n,fjac,ldfjac,.true.,ipvt,n,wa1,wa2,wa3)
!
!        on the first iteration and if mode is 1, scale according
!        to the norms of the columns of the initial jacobian.
!
         if (iter /= 1) go to 80
         if (mode == 2) go to 60
         do 50 j = 1, n
            diag(j) = wa2(j)
            if (wa2(j) == zero) diag(j) = one
   50       continue
   60    continue
!
!        on the first iteration, calculate the norm of the scaled x
!        and initialize the step bound delta.
!
         do 70 j = 1, n
            wa3(j) = diag(j)*x(j)
   70       continue
         xnorm = enorm(n,wa3)
         delta = factor*xnorm
         if (delta == zero) delta = factor
   80    continue
!
!        form (q transpose)*fvec and store the first n components in
!        qtf.
!
         do 90 i = 1, m
            wa4(i) = fvec(i)
   90       continue
         do 130 j = 1, n
            if (fjac(j,j) == zero) go to 120
            sum = zero
            do 100 i = j, m
               sum = sum + fjac(i,j)*wa4(i)
  100          continue
            temp = -sum/fjac(j,j)
            do 110 i = j, m
               wa4(i) = wa4(i) + fjac(i,j)*temp
  110          continue
  120       continue
            fjac(j,j) = wa1(j)
            qtf(j) = wa4(j)
  130       continue
!
!        compute the norm of the scaled gradient.
!
         gnorm = zero
         if (fnorm == zero) go to 170
         do 160 j = 1, n
            l = ipvt(j)
            if (wa2(l) == zero) go to 150
            sum = zero
            do 140 i = 1, j
               sum = sum + fjac(i,j)*(qtf(i)/fnorm)
  140          continue
            gnorm = dmax1(gnorm,dabs(sum/wa2(l)))
  150       continue
  160       continue
  170    continue
!
!        test for convergence of the gradient norm.
!
         if (gnorm <= gtol) info = 4
         if (info /= 0) go to 300
!
!        rescale if necessary.
!
         if (mode == 2) go to 190
         do 180 j = 1, n
            diag(j) = dmax1(diag(j),wa2(j))
  180       continue
  190    continue
!
!        beginning of the inner loop.
!
  200    continue
!
!           determine the levenberg-marquardt parameter.
!
            call lmpar(n,fjac,ldfjac,ipvt,diag,qtf,delta,par,wa1,wa2,&
                       wa3,wa4)
!
!           store the direction p and x + p. calculate the norm of p.
!
            do 210 j = 1, n
               wa1(j) = -wa1(j)
               wa2(j) = x(j) + wa1(j)
               wa3(j) = diag(j)*wa1(j)
  210          continue
            pnorm = enorm(n,wa3)
!
!           on the first iteration, adjust the initial step bound.
!
            if (iter == 1) delta = dmin1(delta,pnorm)
!
!           evaluate the function at x + p and calculate its norm.
!
            iflag = 1
            call fcn(m,n,wa2,wa4,fjac,ldfjac,iflag)
            nfev = nfev + 1
            if (iflag < 0) go to 300
            fnorm1 = enorm(m,wa4)
!
!           compute the scaled actual reduction.
!
            actred = -one
            if (p1*fnorm1 < fnorm) actred = one - (fnorm1/fnorm)**2
!
!           compute the scaled predicted reduction and
!           the scaled directional derivative.
!
            do 230 j = 1, n
               wa3(j) = zero
               l = ipvt(j)
               temp = wa1(l)
               do 220 i = 1, j
                  wa3(i) = wa3(i) + fjac(i,j)*temp
  220             continue
  230          continue
            temp1 = enorm(n,wa3)/fnorm
            temp2 = (dsqrt(par)*pnorm)/fnorm
            prered = temp1**2 + temp2**2/p5
            dirder = -(temp1**2 + temp2**2)
!
!           compute the ratio of the actual to the predicted
!           reduction.
!
            ratio = zero
            if (prered /= zero) ratio = actred/prered
!
!           update the step bound.
!
            if (ratio > p25) go to 240
               if (actred >= zero) temp = p5
               if (actred < zero)&
                  temp = p5*dirder/(dirder + p5*actred)
               if (p1*fnorm1 >= fnorm .or. temp < p1) temp = p1
               delta = temp*dmin1(delta,pnorm/p1)
               par = par/temp
               go to 260
  240       continue
               if (par /= zero .and. ratio < p75) go to 250
               delta = pnorm/p5
               par = p5*par
  250          continue
  260       continue
!
!           test for successful iteration.
!
            if (ratio < p0001) go to 290
!
!           successful iteration. update x, fvec, and their norms.
!
            do 270 j = 1, n
               x(j) = wa2(j)
               wa2(j) = diag(j)*x(j)
  270          continue
            do 280 i = 1, m
               fvec(i) = wa4(i)
  280          continue
            xnorm = enorm(n,wa2)
            fnorm = fnorm1
            iter = iter + 1
  290       continue
!
!           tests for convergence.
!
            if (dabs(actred) <= ftol .and. prered <= ftol&
                .and. p5*ratio <= one) info = 1
            if (delta <= xtol*xnorm) info = 2
            if (dabs(actred) <= ftol .and. prered <= ftol&
                .and. p5*ratio <= one .and. info == 2) info = 3
            if (info /= 0) go to 300
!
!           tests for termination and stringent tolerances.
!
            if (nfev >= maxfev) info = 5
            if (dabs(actred) <= epsmch .and. prered <= epsmch&
                .and. p5*ratio <= one) info = 6
            if (delta <= epsmch*xnorm) info = 7
            if (gnorm <= epsmch) info = 8
            if (info /= 0) go to 300
!
!           end of the inner loop. repeat if iteration unsuccessful.
!
            if (ratio < p0001) go to 200
!
!        end of the outer loop.
!
         go to 30
  300 continue
!
!     termination, either normal or user imposed.
!
      if (iflag < 0) info = iflag
      iflag = 0
      if (nprint > 0) call fcn(m,n,x,fvec,fjac,ldfjac,iflag)
      return
!
!     last card of subroutine lmder.
!
      end subroutine lmder
      
      subroutine lmder1(fcn,m,n,x,fvec,fjac,ldfjac,tol,info,ipvt,wa,&
                        lwa)
      integer m,n,ldfjac,info,lwa
      integer ipvt(n)
      double precision tol
      double precision x(n),fvec(m),fjac(ldfjac,n),wa(lwa)
      external fcn
!     **********
!
!     subroutine lmder1
!
!     the purpose of lmder1 is to minimize the sum of the squares of
!     m nonlinear functions in n variables by a modification of the
!     levenberg-marquardt algorithm. this is done by using the more
!     general least-squares solver lmder. the user must provide a
!     subroutine which calculates the functions and the jacobian.
!
!     the subroutine statement is
!
!       subroutine lmder1(fcn,m,n,x,fvec,fjac,ldfjac,tol,info,
!                         ipvt,wa,lwa)
!
!     where
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions and the jacobian. fcn must
!         be declared in an external statement in the user
!         calling program, and should be written as follows.
!
!         subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
!         integer m,n,ldfjac,iflag
!         double precision x(n),fvec(m),fjac(ldfjac,n)
!         ----------
!         if iflag = 1 calculate the functions at x and
!         return this vector in fvec. do not alter fjac.
!         if iflag = 2 calculate the jacobian at x and
!         return this matrix in fjac. do not alter fvec.
!         ----------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of lmder1.
!         in this case set iflag to a negative integer.
!
!       m is a positive integer input variable set to the number
!         of functions.
!
!       n is a positive integer input variable set to the number
!         of variables. n must not exceed m.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length m which contains
!         the functions evaluated at the output x.
!
!       fjac is an output m by n array. the upper n by n submatrix
!         of fjac contains an upper triangular matrix r with
!         diagonal elements of nonincreasing magnitude such that
!
!                t     t           t
!               p *(jac *jac)*p = r *r,
!
!         where p is a permutation matrix and jac is the final
!         calculated jacobian. column j of p is column ipvt(j)
!         (see below) of the identity matrix. the lower trapezoidal
!         part of fjac contains information generated during
!         the computation of r.
!
!       ldfjac is a positive integer input variable not less than m
!         which specifies the leading dimension of the array fjac.
!
!       tol is a nonnegative input variable. termination occurs
!         when the algorithm estimates either that the relative
!         error in the sum of squares is at most tol or that
!         the relative error between x and the solution is at
!         most tol.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0  improper input parameters.
!
!         info = 1  algorithm estimates that the relative error
!                   in the sum of squares is at most tol.
!
!         info = 2  algorithm estimates that the relative error
!                   between x and the solution is at most tol.
!
!         info = 3  conditions for info = 1 and info = 2 both hold.
!
!         info = 4  fvec is orthogonal to the columns of the
!                   jacobian to machine precision.
!
!         info = 5  number of calls to fcn with iflag = 1 has
!                   reached 100*(n+1).
!
!         info = 6  tol is too small. no further reduction in
!                   the sum of squares is possible.
!
!         info = 7  tol is too small. no further improvement in
!                   the approximate solution x is possible.
!
!       ipvt is an integer output array of length n. ipvt
!         defines a permutation matrix p such that jac*p = q*r,
!         where jac is the final calculated jacobian, q is
!         orthogonal (not stored), and r is upper triangular
!         with diagonal elements of nonincreasing magnitude.
!         column j of p is column ipvt(j) of the identity matrix.
!
!       wa is a work array of length lwa.
!
!       lwa is a positive integer input variable not less than 5*n+m.
!
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... lmder
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      integer maxfev,mode,nfev,njev,nprint
      double precision factor,ftol,gtol,xtol,zero
      data factor,zero /1.0d2,0.0d0/
      info = 0
!
!     check the input parameters for errors.
!
      if (n <= 0 .or. m < n .or. ldfjac < m .or. tol < zero&
          .or. lwa < 5*n + m) go to 10
!
!     call lmder.
!
      maxfev = 100*(n + 1)
      ftol = tol
      xtol = tol
      gtol = zero
      mode = 1
      nprint = 0
      call lmder(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,maxfev,&
                 wa(1),mode,factor,nprint,info,nfev,njev,ipvt,wa(n+1),&
                 wa(2*n+1),wa(3*n+1),wa(4*n+1),wa(5*n+1))
      if (info == 8) info = 4
   10 continue
      return
!
!     last card of subroutine lmder1.
!
      end subroutine lmder1
      
      subroutine lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,&
                       diag,mode,factor,nprint,info,nfev,fjac,ldfjac,&
                       ipvt,qtf,wa1,wa2,wa3,wa4)
      integer m,n,maxfev,mode,nprint,info,nfev,ldfjac
      integer ipvt(n)
      double precision ftol,xtol,gtol,epsfcn,factor
      double precision x(n),fvec(m),diag(n),fjac(ldfjac,n),qtf(n),&
                       wa1(n),wa2(n),wa3(n),wa4(m)
      external fcn
!     **********
!
!     subroutine lmdif
!
!     the purpose of lmdif is to minimize the sum of the squares of
!     m nonlinear functions in n variables by a modification of
!     the levenberg-marquardt algorithm. the user must provide a
!     subroutine which calculates the functions. the jacobian is
!     then calculated by a forward-difference approximation.
!
!     the subroutine statement is
!
!       subroutine lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,
!                        diag,mode,factor,nprint,info,nfev,fjac,
!                        ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
!
!     where
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions. fcn must be declared
!         in an external statement in the user calling
!         program, and should be written as follows.
!
!         subroutine fcn(m,n,x,fvec,iflag)
!         integer m,n,iflag
!         double precision x(n),fvec(m)
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ----------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of lmdif.
!         in this case set iflag to a negative integer.
!
!       m is a positive integer input variable set to the number
!         of functions.
!
!       n is a positive integer input variable set to the number
!         of variables. n must not exceed m.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length m which contains
!         the functions evaluated at the output x.
!
!       ftol is a nonnegative input variable. termination
!         occurs when both the actual and predicted relative
!         reductions in the sum of squares are at most ftol.
!         therefore, ftol measures the relative error desired
!         in the sum of squares.
!
!       xtol is a nonnegative input variable. termination
!         occurs when the relative error between two consecutive
!         iterates is at most xtol. therefore, xtol measures the
!         relative error desired in the approximate solution.
!
!       gtol is a nonnegative input variable. termination
!         occurs when the cosine of the angle between fvec and
!         any column of the jacobian is at most gtol in absolute
!         value. therefore, gtol measures the orthogonality
!         desired between the function vector and the columns
!         of the jacobian.
!
!       maxfev is a positive integer input variable. termination
!         occurs when the number of calls to fcn is at least
!         maxfev by the end of an iteration.
!
!       epsfcn is an input variable used in determining a suitable
!         step length for the forward-difference approximation. this
!         approximation assumes that the relative errors in the
!         functions are of the order of epsfcn. if epsfcn is less
!         than the machine precision, it is assumed that the relative
!         errors in the functions are of the order of the machine
!         precision.
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
!         for printing. if nprint is not positive, no special calls
!         of fcn with iflag = 0 are made.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0  improper input parameters.
!
!         info = 1  both actual and predicted relative reductions
!                   in the sum of squares are at most ftol.
!
!         info = 2  relative error between two consecutive iterates
!                   is at most xtol.
!
!         info = 3  conditions for info = 1 and info = 2 both hold.
!
!         info = 4  the cosine of the angle between fvec and any
!                   column of the jacobian is at most gtol in
!                   absolute value.
!
!         info = 5  number of calls to fcn has reached or
!                   exceeded maxfev.
!
!         info = 6  ftol is too small. no further reduction in
!                   the sum of squares is possible.
!
!         info = 7  xtol is too small. no further improvement in
!                   the approximate solution x is possible.
!
!         info = 8  gtol is too small. fvec is orthogonal to the
!                   columns of the jacobian to machine precision.
!
!       nfev is an integer output variable set to the number of
!         calls to fcn.
!
!       fjac is an output m by n array. the upper n by n submatrix
!         of fjac contains an upper triangular matrix r with
!         diagonal elements of nonincreasing magnitude such that
!
!                t     t           t
!               p *(jac *jac)*p = r *r,
!
!         where p is a permutation matrix and jac is the final
!         calculated jacobian. column j of p is column ipvt(j)
!         (see below) of the identity matrix. the lower trapezoidal
!         part of fjac contains information generated during
!         the computation of r.
!
!       ldfjac is a positive integer input variable not less than m
!         which specifies the leading dimension of the array fjac.
!
!       ipvt is an integer output array of length n. ipvt
!         defines a permutation matrix p such that jac*p = q*r,
!         where jac is the final calculated jacobian, q is
!         orthogonal (not stored), and r is upper triangular
!         with diagonal elements of nonincreasing magnitude.
!         column j of p is column ipvt(j) of the identity matrix.
!
!       qtf is an output array of length n which contains
!         the first n elements of the vector (q transpose)*fvec.
!
!       wa1, wa2, and wa3 are work arrays of length n.
!
!       wa4 is a work array of length m.
!
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... dpmpar,enorm,fdjac2,lmpar,qrfac
!
!       fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      integer i,iflag,iter,j,l
      double precision actred,delta,dirder,epsmch,fnorm,fnorm1,gnorm,&
                       one,par,pnorm,prered,p1,p5,p25,p75,p0001,ratio,&
                       sum,temp,temp1,temp2,xnorm,zero
      data one,p1,p5,p25,p75,p0001,zero&
           /1.0d0,1.0d-1,5.0d-1,2.5d-1,7.5d-1,1.0d-4,0.0d0/
!
!     epsmch is the machine precision.
!
      epsmch = dpmpar(1)
!
      info = 0
      iflag = 0
      nfev = 0
!
!     check the input parameters for errors.
!
      if (n <= 0 .or. m < n .or. ldfjac < m&
          .or. ftol < zero .or. xtol < zero .or. gtol < zero&
          .or. maxfev <= 0 .or. factor <= zero) go to 300
      if (mode /= 2) go to 20
      do 10 j = 1, n
         if (diag(j) <= zero) go to 300
   10    continue
   20 continue
!
!     evaluate the function at the starting point
!     and calculate its norm.
!
      iflag = 1
      call fcn(m,n,x,fvec,iflag)
      nfev = 1
      if (iflag < 0) go to 300
      fnorm = enorm(m,fvec)
!
!     initialize levenberg-marquardt parameter and iteration counter.
!
      par = zero
      iter = 1
!
!     beginning of the outer loop.
!
   30 continue
!
!        calculate the jacobian matrix.
!
         iflag = 2
         call fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa4)
         nfev = nfev + n
         if (iflag < 0) go to 300
!
!        if requested, call fcn to enable printing of iterates.
!
         if (nprint <= 0) go to 40
         iflag = 0
         if (mod(iter-1,nprint) == 0) call fcn(m,n,x,fvec,iflag)
         if (iflag < 0) go to 300
   40    continue
!
!        compute the qr factorization of the jacobian.
!
         call qrfac(m,n,fjac,ldfjac,.true.,ipvt,n,wa1,wa2,wa3)
!
!        on the first iteration and if mode is 1, scale according
!        to the norms of the columns of the initial jacobian.
!
         if (iter /= 1) go to 80
         if (mode == 2) go to 60
         do 50 j = 1, n
            diag(j) = wa2(j)
            if (wa2(j) == zero) diag(j) = one
   50       continue
   60    continue
!
!        on the first iteration, calculate the norm of the scaled x
!        and initialize the step bound delta.
!
         do 70 j = 1, n
            wa3(j) = diag(j)*x(j)
   70       continue
         xnorm = enorm(n,wa3)
         delta = factor*xnorm
         if (delta == zero) delta = factor
   80    continue
!
!        form (q transpose)*fvec and store the first n components in
!        qtf.
!
         do 90 i = 1, m
            wa4(i) = fvec(i)
   90       continue
         do 130 j = 1, n
            if (fjac(j,j) == zero) go to 120
            sum = zero
            do 100 i = j, m
               sum = sum + fjac(i,j)*wa4(i)
  100          continue
            temp = -sum/fjac(j,j)
            do 110 i = j, m
               wa4(i) = wa4(i) + fjac(i,j)*temp
  110          continue
  120       continue
            fjac(j,j) = wa1(j)
            qtf(j) = wa4(j)
  130       continue
!
!        compute the norm of the scaled gradient.
!
         gnorm = zero
         if (fnorm == zero) go to 170
         do 160 j = 1, n
            l = ipvt(j)
            if (wa2(l) == zero) go to 150
            sum = zero
            do 140 i = 1, j
               sum = sum + fjac(i,j)*(qtf(i)/fnorm)
  140          continue
            gnorm = dmax1(gnorm,dabs(sum/wa2(l)))
  150       continue
  160       continue
  170    continue
!
!        test for convergence of the gradient norm.
!
         if (gnorm <= gtol) info = 4
         if (info /= 0) go to 300
!
!        rescale if necessary.
!
         if (mode == 2) go to 190
         do 180 j = 1, n
            diag(j) = dmax1(diag(j),wa2(j))
  180       continue
  190    continue
!
!        beginning of the inner loop.
!
  200    continue
!
!           determine the levenberg-marquardt parameter.
!
            call lmpar(n,fjac,ldfjac,ipvt,diag,qtf,delta,par,wa1,wa2,&
                       wa3,wa4)
!
!           store the direction p and x + p. calculate the norm of p.
!
            do 210 j = 1, n
               wa1(j) = -wa1(j)
               wa2(j) = x(j) + wa1(j)
               wa3(j) = diag(j)*wa1(j)
  210          continue
            pnorm = enorm(n,wa3)
!
!           on the first iteration, adjust the initial step bound.
!
            if (iter == 1) delta = dmin1(delta,pnorm)
!
!           evaluate the function at x + p and calculate its norm.
!
            iflag = 1
            call fcn(m,n,wa2,wa4,iflag)
            nfev = nfev + 1
            if (iflag < 0) go to 300
            fnorm1 = enorm(m,wa4)
!
!           compute the scaled actual reduction.
!
            actred = -one
            if (p1*fnorm1 < fnorm) actred = one - (fnorm1/fnorm)**2
!
!           compute the scaled predicted reduction and
!           the scaled directional derivative.
!
            do 230 j = 1, n
               wa3(j) = zero
               l = ipvt(j)
               temp = wa1(l)
               do 220 i = 1, j
                  wa3(i) = wa3(i) + fjac(i,j)*temp
  220             continue
  230          continue
            temp1 = enorm(n,wa3)/fnorm
            temp2 = (dsqrt(par)*pnorm)/fnorm
            prered = temp1**2 + temp2**2/p5
            dirder = -(temp1**2 + temp2**2)
!
!           compute the ratio of the actual to the predicted
!           reduction.
!
            ratio = zero
            if (prered /= zero) ratio = actred/prered
!
!           update the step bound.
!
            if (ratio > p25) go to 240
               if (actred >= zero) temp = p5
               if (actred < zero)&
                  temp = p5*dirder/(dirder + p5*actred)
               if (p1*fnorm1 >= fnorm .or. temp < p1) temp = p1
               delta = temp*dmin1(delta,pnorm/p1)
               par = par/temp
               go to 260
  240       continue
               if (par /= zero .and. ratio < p75) go to 250
               delta = pnorm/p5
               par = p5*par
  250          continue
  260       continue
!
!           test for successful iteration.
!
            if (ratio < p0001) go to 290
!
!           successful iteration. update x, fvec, and their norms.
!
            do 270 j = 1, n
               x(j) = wa2(j)
               wa2(j) = diag(j)*x(j)
  270          continue
            do 280 i = 1, m
               fvec(i) = wa4(i)
  280          continue
            xnorm = enorm(n,wa2)
            fnorm = fnorm1
            iter = iter + 1
  290       continue
!
!           tests for convergence.
!
            if (dabs(actred) <= ftol .and. prered <= ftol&
                .and. p5*ratio <= one) info = 1
            if (delta <= xtol*xnorm) info = 2
            if (dabs(actred) <= ftol .and. prered <= ftol&
                .and. p5*ratio <= one .and. info == 2) info = 3
            if (info /= 0) go to 300
!
!           tests for termination and stringent tolerances.
!
            if (nfev >= maxfev) info = 5
            if (dabs(actred) <= epsmch .and. prered <= epsmch&
                .and. p5*ratio <= one) info = 6
            if (delta <= epsmch*xnorm) info = 7
            if (gnorm <= epsmch) info = 8
            if (info /= 0) go to 300
!
!           end of the inner loop. repeat if iteration unsuccessful.
!
            if (ratio < p0001) go to 200
!
!        end of the outer loop.
!
         go to 30
  300 continue
!
!     termination, either normal or user imposed.
!
      if (iflag < 0) info = iflag
      iflag = 0
      if (nprint > 0) call fcn(m,n,x,fvec,iflag)
      return
!
!     last card of subroutine lmdif.
!
      end subroutine lmdif
      
      subroutine lmdif1(fcn,m,n,x,fvec,tol,info,iwa,wa,lwa)
      integer m,n,info,lwa
      integer iwa(n)
      double precision tol
      double precision x(n),fvec(m),wa(lwa)
      external fcn
!     **********
!
!     subroutine lmdif1
!
!     the purpose of lmdif1 is to minimize the sum of the squares of
!     m nonlinear functions in n variables by a modification of the
!     levenberg-marquardt algorithm. this is done by using the more
!     general least-squares solver lmdif. the user must provide a
!     subroutine which calculates the functions. the jacobian is
!     then calculated by a forward-difference approximation.
!
!     the subroutine statement is
!
!       subroutine lmdif1(fcn,m,n,x,fvec,tol,info,iwa,wa,lwa)
!
!     where
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions. fcn must be declared
!         in an external statement in the user calling
!         program, and should be written as follows.
!
!         subroutine fcn(m,n,x,fvec,iflag)
!         integer m,n,iflag
!         double precision x(n),fvec(m)
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ----------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of lmdif1.
!         in this case set iflag to a negative integer.
!
!       m is a positive integer input variable set to the number
!         of functions.
!
!       n is a positive integer input variable set to the number
!         of variables. n must not exceed m.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length m which contains
!         the functions evaluated at the output x.
!
!       tol is a nonnegative input variable. termination occurs
!         when the algorithm estimates either that the relative
!         error in the sum of squares is at most tol or that
!         the relative error between x and the solution is at
!         most tol.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0  improper input parameters.
!
!         info = 1  algorithm estimates that the relative error
!                   in the sum of squares is at most tol.
!
!         info = 2  algorithm estimates that the relative error
!                   between x and the solution is at most tol.
!
!         info = 3  conditions for info = 1 and info = 2 both hold.
!
!         info = 4  fvec is orthogonal to the columns of the
!                   jacobian to machine precision.
!
!         info = 5  number of calls to fcn has reached or
!                   exceeded 200*(n+1).
!
!         info = 6  tol is too small. no further reduction in
!                   the sum of squares is possible.
!
!         info = 7  tol is too small. no further improvement in
!                   the approximate solution x is possible.
!
!       iwa is an integer work array of length n.
!
!       wa is a work array of length lwa.
!
!       lwa is a positive integer input variable not less than
!         m*n+5*n+m.
!
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... lmdif
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      integer maxfev,mode,mp5n,nfev,nprint
      double precision epsfcn,factor,ftol,gtol,xtol,zero
      data factor,zero /1.0d2,0.0d0/
      info = 0
!
!     check the input parameters for errors.
!
      if (n <= 0 .or. m < n .or. tol < zero&
          .or. lwa < m*n + 5*n + m) go to 10
!
!     call lmdif.
!
      maxfev = 200*(n + 1)
      ftol = tol
      xtol = tol
      gtol = zero
      epsfcn = zero
      mode = 1
      nprint = 0
      mp5n = m + 5*n
      call lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,wa(1),&
                 mode,factor,nprint,info,nfev,wa(mp5n+1),m,iwa,&
                 wa(n+1),wa(2*n+1),wa(3*n+1),wa(4*n+1),wa(5*n+1))
      if (info == 8) info = 4
   10 continue
      return
!
!     last card of subroutine lmdif1.
!
      end subroutine lmdif1
      
      subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,wa1,&
                       wa2)
      integer n,ldr
      integer ipvt(n)
      double precision delta,par
      double precision r(ldr,n),diag(n),qtb(n),x(n),sdiag(n),wa1(n),&
                       wa2(n)
!     **********
!
!     subroutine lmpar
!
!     given an m by n matrix a, an n by n nonsingular diagonal
!     matrix d, an m-vector b, and a positive number delta,
!     the problem is to determine a value for the parameter
!     par such that if x solves the system
!
!           a*x = b ,     sqrt(par)*d*x = 0 ,
!
!     in the least squares sense, and dxnorm is the euclidean
!     norm of d*x, then either par is zero and
!
!           (dxnorm-delta) <= 0.1*delta ,
!
!     or par is positive and
!
!           abs(dxnorm-delta) <= 0.1*delta .
!
!     this subroutine completes the solution of the problem
!     if it is provided with the necessary information from the
!     qr factorization, with column pivoting, of a. that is, if
!     a*p = q*r, where p is a permutation matrix, q has orthogonal
!     columns, and r is an upper triangular matrix with diagonal
!     elements of nonincreasing magnitude, then lmpar expects
!     the full upper triangle of r, the permutation matrix p,
!     and the first n components of (q transpose)*b. on output
!     lmpar also provides an upper triangular matrix s such that
!
!            t   t                   t
!           p *(a *a + par*d*d)*p = s *s .
!
!     s is employed within lmpar and may be of separate interest.
!
!     only a few iterations are generally needed for convergence
!     of the algorithm. if, however, the limit of 10 iterations
!     is reached, then the output par will contain the best
!     value obtained so far.
!
!     the subroutine statement is
!
!       subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,
!                        wa1,wa2)
!
!     where
!
!       n is a positive integer input variable set to the order of r.
!
!       r is an n by n array. on input the full upper triangle
!         must contain the full upper triangle of the matrix r.
!         on output the full upper triangle is unaltered, and the
!         strict lower triangle contains the strict upper triangle
!         (transposed) of the upper triangular matrix s.
!
!       ldr is a positive integer input variable not less than n
!         which specifies the leading dimension of the array r.
!
!       ipvt is an integer input array of length n which defines the
!         permutation matrix p such that a*p = q*r. column j of p
!         is column ipvt(j) of the identity matrix.
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
!       par is a nonnegative variable. on input par contains an
!         initial estimate of the levenberg-marquardt parameter.
!         on output par contains the final estimate.
!
!       x is an output array of length n which contains the least
!         squares solution of the system a*x = b, sqrt(par)*d*x = 0,
!         for the output par.
!
!       sdiag is an output array of length n which contains the
!         diagonal elements of the upper triangular matrix s.
!
!       wa1 and wa2 are work arrays of length n.
!
!     subprograms called
!
!       minpack-supplied ... dpmpar,enorm,qrsolv
!
!       fortran-supplied ... dabs,dmax1,dmin1,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      integer i,iter,j,jm1,jp1,k,l,nsing
      double precision dxnorm,dwarf,fp,gnorm,parc,parl,paru,p1,p001,&
                       sum,temp,zero
      data p1,p001,zero /1.0d-1,1.0d-3,0.0d0/
!
!     dwarf is the smallest positive magnitude.
!
      dwarf = dpmpar(2)
!
!     compute and store in x the gauss-newton direction. if the
!     jacobian is rank-deficient, obtain a least squares solution.
!
      nsing = n
      do 10 j = 1, n
         wa1(j) = qtb(j)
         if (r(j,j) == zero .and. nsing == n) nsing = j - 1
         if (nsing < n) wa1(j) = zero
   10    continue
      if (nsing < 1) go to 50
      do 40 k = 1, nsing
         j = nsing - k + 1
         wa1(j) = wa1(j)/r(j,j)
         temp = wa1(j)
         jm1 = j - 1
         if (jm1 < 1) go to 30
         do 20 i = 1, jm1
            wa1(i) = wa1(i) - r(i,j)*temp
   20       continue
   30    continue
   40    continue
   50 continue
      do 60 j = 1, n
         l = ipvt(j)
         x(l) = wa1(j)
   60    continue
!
!     initialize the iteration counter.
!     evaluate the function at the origin, and test
!     for acceptance of the gauss-newton direction.
!
      iter = 0
      do 70 j = 1, n
         wa2(j) = diag(j)*x(j)
   70    continue
      dxnorm = enorm(n,wa2)
      fp = dxnorm - delta
      if (fp <= p1*delta) go to 220
!
!     if the jacobian is not rank deficient, the newton
!     step provides a lower bound, parl, for the zero of
!     the function. otherwise set this bound to zero.
!
      parl = zero
      if (nsing < n) go to 120
      do 80 j = 1, n
         l = ipvt(j)
         wa1(j) = diag(l)*(wa2(l)/dxnorm)
   80    continue
      do 110 j = 1, n
         sum = zero
         jm1 = j - 1
         if (jm1 < 1) go to 100
         do 90 i = 1, jm1
            sum = sum + r(i,j)*wa1(i)
   90       continue
  100    continue
         wa1(j) = (wa1(j) - sum)/r(j,j)
  110    continue
      temp = enorm(n,wa1)
      parl = ((fp/delta)/temp)/temp
  120 continue
!
!     calculate an upper bound, paru, for the zero of the function.
!
      do 140 j = 1, n
         sum = zero
         do 130 i = 1, j
            sum = sum + r(i,j)*qtb(i)
  130       continue
         l = ipvt(j)
         wa1(j) = sum/diag(l)
  140    continue
      gnorm = enorm(n,wa1)
      paru = gnorm/delta
      if (paru == zero) paru = dwarf/dmin1(delta,p1)
!
!     if the input par lies outside of the interval (parl,paru),
!     set par to the closer endpoint.
!
      par = dmax1(par,parl)
      par = dmin1(par,paru)
      if (par == zero) par = gnorm/dxnorm
!
!     beginning of an iteration.
!
  150 continue
         iter = iter + 1
!
!        evaluate the function at the current value of par.
!
         if (par == zero) par = dmax1(dwarf,p001*paru)
         temp = dsqrt(par)
         do 160 j = 1, n
            wa1(j) = temp*diag(j)
  160       continue
         call qrsolv(n,r,ldr,ipvt,wa1,qtb,x,sdiag,wa2)
         do 170 j = 1, n
            wa2(j) = diag(j)*x(j)
  170       continue
         dxnorm = enorm(n,wa2)
         temp = fp
         fp = dxnorm - delta
!
!        if the function is small enough, accept the current value
!        of par. also test for the exceptional cases where parl
!        is zero or the number of iterations has reached 10.
!
         if (dabs(fp) <= p1*delta&
             .or. parl == zero .and. fp <= temp&
                  .and. temp < zero .or. iter == 10) go to 220
!
!        compute the newton correction.
!
         do 180 j = 1, n
            l = ipvt(j)
            wa1(j) = diag(l)*(wa2(l)/dxnorm)
  180       continue
         do 210 j = 1, n
            wa1(j) = wa1(j)/sdiag(j)
            temp = wa1(j)
            jp1 = j + 1
            if (n < jp1) go to 200
            do 190 i = jp1, n
               wa1(i) = wa1(i) - r(i,j)*temp
  190          continue
  200       continue
  210       continue
         temp = enorm(n,wa1)
         parc = ((fp/delta)/temp)/temp
!
!        depending on the sign of the function, update parl or paru.
!
         if (fp > zero) parl = dmax1(parl,par)
         if (fp < zero) paru = dmin1(paru,par)
!
!        compute an improved estimate for par.
!
         par = dmax1(parl,par+parc)
!
!        end of an iteration.
!
         go to 150
  220 continue
!
!     termination.
!
      if (iter == 0) par = zero
      return
!
!     last card of subroutine lmpar.
!
      end subroutine lmpar
      
      subroutine lmstr(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,&
                       maxfev,diag,mode,factor,nprint,info,nfev,njev,&
                       ipvt,qtf,wa1,wa2,wa3,wa4)
      integer m,n,ldfjac,maxfev,mode,nprint,info,nfev,njev
      integer ipvt(n)
      logical sing
      double precision ftol,xtol,gtol,factor
      double precision x(n),fvec(m),fjac(ldfjac,n),diag(n),qtf(n),&
                       wa1(n),wa2(n),wa3(n),wa4(m)
!     **********
!
!     subroutine lmstr
!
!     the purpose of lmstr is to minimize the sum of the squares of
!     m nonlinear functions in n variables by a modification of
!     the levenberg-marquardt algorithm which uses minimal storage.
!     the user must provide a subroutine which calculates the
!     functions and the rows of the jacobian.
!
!     the subroutine statement is
!
!       subroutine lmstr(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,
!                        maxfev,diag,mode,factor,nprint,info,nfev,
!                        njev,ipvt,qtf,wa1,wa2,wa3,wa4)
!
!     where
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions and the rows of the jacobian.
!         fcn must be declared in an external statement in the
!         user calling program, and should be written as follows.
!
!         subroutine fcn(m,n,x,fvec,fjrow,iflag)
!         integer m,n,iflag
!         double precision x(n),fvec(m),fjrow(n)
!         ----------
!         if iflag = 1 calculate the functions at x and
!         return this vector in fvec.
!         if iflag = i calculate the (i-1)-st row of the
!         jacobian at x and return this vector in fjrow.
!         ----------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of lmstr.
!         in this case set iflag to a negative integer.
!
!       m is a positive integer input variable set to the number
!         of functions.
!
!       n is a positive integer input variable set to the number
!         of variables. n must not exceed m.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length m which contains
!         the functions evaluated at the output x.
!
!       fjac is an output n by n array. the upper triangle of fjac
!         contains an upper triangular matrix r such that
!
!                t     t           t
!               p *(jac *jac)*p = r *r,
!
!         where p is a permutation matrix and jac is the final
!         calculated jacobian. column j of p is column ipvt(j)
!         (see below) of the identity matrix. the lower triangular
!         part of fjac contains information generated during
!         the computation of r.
!
!       ldfjac is a positive integer input variable not less than n
!         which specifies the leading dimension of the array fjac.
!
!       ftol is a nonnegative input variable. termination
!         occurs when both the actual and predicted relative
!         reductions in the sum of squares are at most ftol.
!         therefore, ftol measures the relative error desired
!         in the sum of squares.
!
!       xtol is a nonnegative input variable. termination
!         occurs when the relative error between two consecutive
!         iterates is at most xtol. therefore, xtol measures the
!         relative error desired in the approximate solution.
!
!       gtol is a nonnegative input variable. termination
!         occurs when the cosine of the angle between fvec and
!         any column of the jacobian is at most gtol in absolute
!         value. therefore, gtol measures the orthogonality
!         desired between the function vector and the columns
!         of the jacobian.
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
!         for printing. if nprint is not positive, no special calls
!         of fcn with iflag = 0 are made.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0  improper input parameters.
!
!         info = 1  both actual and predicted relative reductions
!                   in the sum of squares are at most ftol.
!
!         info = 2  relative error between two consecutive iterates
!                   is at most xtol.
!
!         info = 3  conditions for info = 1 and info = 2 both hold.
!
!         info = 4  the cosine of the angle between fvec and any
!                   column of the jacobian is at most gtol in
!                   absolute value.
!
!         info = 5  number of calls to fcn with iflag = 1 has
!                   reached maxfev.
!
!         info = 6  ftol is too small. no further reduction in
!                   the sum of squares is possible.
!
!         info = 7  xtol is too small. no further improvement in
!                   the approximate solution x is possible.
!
!         info = 8  gtol is too small. fvec is orthogonal to the
!                   columns of the jacobian to machine precision.
!
!       nfev is an integer output variable set to the number of
!         calls to fcn with iflag = 1.
!
!       njev is an integer output variable set to the number of
!         calls to fcn with iflag = 2.
!
!       ipvt is an integer output array of length n. ipvt
!         defines a permutation matrix p such that jac*p = q*r,
!         where jac is the final calculated jacobian, q is
!         orthogonal (not stored), and r is upper triangular.
!         column j of p is column ipvt(j) of the identity matrix.
!
!       qtf is an output array of length n which contains
!         the first n elements of the vector (q transpose)*fvec.
!
!       wa1, wa2, and wa3 are work arrays of length n.
!
!       wa4 is a work array of length m.
!
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... dpmpar,enorm,lmpar,qrfac,rwupdt
!
!       fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, dudley v. goetschel, kenneth e. hillstrom,
!     jorge j. more
!
!     **********
      integer i,iflag,iter,j,l
      double precision actred,delta,dirder,epsmch,fnorm,fnorm1,gnorm,&
                       one,par,pnorm,prered,p1,p5,p25,p75,p0001,ratio,&
                       sum,temp,temp1,temp2,xnorm,zero
      data one,p1,p5,p25,p75,p0001,zero&
           /1.0d0,1.0d-1,5.0d-1,2.5d-1,7.5d-1,1.0d-4,0.0d0/
!
!     epsmch is the machine precision.
!
      epsmch = dpmpar(1)
!
      info = 0
      iflag = 0
      nfev = 0
      njev = 0
!
!     check the input parameters for errors.
!
      if (n <= 0 .or. m < n .or. ldfjac < n&
          .or. ftol < zero .or. xtol < zero .or. gtol < zero&
          .or. maxfev <= 0 .or. factor <= zero) go to 340
      if (mode /= 2) go to 20
      do 10 j = 1, n
         if (diag(j) <= zero) go to 340
   10    continue
   20 continue
!
!     evaluate the function at the starting point
!     and calculate its norm.
!
      iflag = 1
      call fcn(m,n,x,fvec,wa3,iflag)
      nfev = 1
      if (iflag < 0) go to 340
      fnorm = enorm(m,fvec)
!
!     initialize levenberg-marquardt parameter and iteration counter.
!
      par = zero
      iter = 1
!
!     beginning of the outer loop.
!
   30 continue
!
!        if requested, call fcn to enable printing of iterates.
!
         if (nprint <= 0) go to 40
         iflag = 0
         if (mod(iter-1,nprint) == 0) call fcn(m,n,x,fvec,wa3,iflag)
         if (iflag < 0) go to 340
   40    continue
!
!        compute the qr factorization of the jacobian matrix
!        calculated one row at a time, while simultaneously
!        forming (q transpose)*fvec and storing the first
!        n components in qtf.
!
         do 60 j = 1, n
            qtf(j) = zero
            do 50 i = 1, n
               fjac(i,j) = zero
   50          continue
   60       continue
         iflag = 2
         do 70 i = 1, m
            call fcn(m,n,x,fvec,wa3,iflag)
            if (iflag < 0) go to 340
            temp = fvec(i)
            call rwupdt(n,fjac,ldfjac,wa3,qtf,temp,wa1,wa2)
            iflag = iflag + 1
   70       continue
         njev = njev + 1
!
!        if the jacobian is rank deficient, call qrfac to
!        reorder its columns and update the components of qtf.
!
         sing = .false.
         do 80 j = 1, n
            if (fjac(j,j) == zero) sing = .true.
            ipvt(j) = j
            wa2(j) = enorm(j,fjac(1,j))
   80       continue
         if (.not.sing) go to 130
         call qrfac(n,n,fjac,ldfjac,.true.,ipvt,n,wa1,wa2,wa3)
         do 120 j = 1, n
            if (fjac(j,j) == zero) go to 110
            sum = zero
            do 90 i = j, n
               sum = sum + fjac(i,j)*qtf(i)
   90          continue
            temp = -sum/fjac(j,j)
            do 100 i = j, n
               qtf(i) = qtf(i) + fjac(i,j)*temp
  100          continue
  110       continue
            fjac(j,j) = wa1(j)
  120       continue
  130    continue
!
!        on the first iteration and if mode is 1, scale according
!        to the norms of the columns of the initial jacobian.
!
         if (iter /= 1) go to 170
         if (mode == 2) go to 150
         do 140 j = 1, n
            diag(j) = wa2(j)
            if (wa2(j) == zero) diag(j) = one
  140       continue
  150    continue
!
!        on the first iteration, calculate the norm of the scaled x
!        and initialize the step bound delta.
!
         do 160 j = 1, n
            wa3(j) = diag(j)*x(j)
  160       continue
         xnorm = enorm(n,wa3)
         delta = factor*xnorm
         if (delta == zero) delta = factor
  170    continue
!
!        compute the norm of the scaled gradient.
!
         gnorm = zero
         if (fnorm == zero) go to 210
         do 200 j = 1, n
            l = ipvt(j)
            if (wa2(l) == zero) go to 190
            sum = zero
            do 180 i = 1, j
               sum = sum + fjac(i,j)*(qtf(i)/fnorm)
  180          continue
            gnorm = dmax1(gnorm,dabs(sum/wa2(l)))
  190       continue
  200       continue
  210    continue
!
!        test for convergence of the gradient norm.
!
         if (gnorm <= gtol) info = 4
         if (info /= 0) go to 340
!
!        rescale if necessary.
!
         if (mode == 2) go to 230
         do 220 j = 1, n
            diag(j) = dmax1(diag(j),wa2(j))
  220       continue
  230    continue
!
!        beginning of the inner loop.
!
  240    continue
!
!           determine the levenberg-marquardt parameter.
!
            call lmpar(n,fjac,ldfjac,ipvt,diag,qtf,delta,par,wa1,wa2,&
                       wa3,wa4)
!
!           store the direction p and x + p. calculate the norm of p.
!
            do 250 j = 1, n
               wa1(j) = -wa1(j)
               wa2(j) = x(j) + wa1(j)
               wa3(j) = diag(j)*wa1(j)
  250          continue
            pnorm = enorm(n,wa3)
!
!           on the first iteration, adjust the initial step bound.
!
            if (iter == 1) delta = dmin1(delta,pnorm)
!
!           evaluate the function at x + p and calculate its norm.
!
            iflag = 1
            call fcn(m,n,wa2,wa4,wa3,iflag)
            nfev = nfev + 1
            if (iflag < 0) go to 340
            fnorm1 = enorm(m,wa4)
!
!           compute the scaled actual reduction.
!
            actred = -one
            if (p1*fnorm1 < fnorm) actred = one - (fnorm1/fnorm)**2
!
!           compute the scaled predicted reduction and
!           the scaled directional derivative.
!
            do 270 j = 1, n
               wa3(j) = zero
               l = ipvt(j)
               temp = wa1(l)
               do 260 i = 1, j
                  wa3(i) = wa3(i) + fjac(i,j)*temp
  260             continue
  270          continue
            temp1 = enorm(n,wa3)/fnorm
            temp2 = (dsqrt(par)*pnorm)/fnorm
            prered = temp1**2 + temp2**2/p5
            dirder = -(temp1**2 + temp2**2)
!
!           compute the ratio of the actual to the predicted
!           reduction.
!
            ratio = zero
            if (prered /= zero) ratio = actred/prered
!
!           update the step bound.
!
            if (ratio > p25) go to 280
               if (actred >= zero) temp = p5
               if (actred < zero)&
                  temp = p5*dirder/(dirder + p5*actred)
               if (p1*fnorm1 >= fnorm .or. temp < p1) temp = p1
               delta = temp*dmin1(delta,pnorm/p1)
               par = par/temp
               go to 300
  280       continue
               if (par /= zero .and. ratio < p75) go to 290
               delta = pnorm/p5
               par = p5*par
  290          continue
  300       continue
!
!           test for successful iteration.
!
            if (ratio < p0001) go to 330
!
!           successful iteration. update x, fvec, and their norms.
!
            do 310 j = 1, n
               x(j) = wa2(j)
               wa2(j) = diag(j)*x(j)
  310          continue
            do 320 i = 1, m
               fvec(i) = wa4(i)
  320          continue
            xnorm = enorm(n,wa2)
            fnorm = fnorm1
            iter = iter + 1
  330       continue
!
!           tests for convergence.
!
            if (dabs(actred) <= ftol .and. prered <= ftol&
                .and. p5*ratio <= one) info = 1
            if (delta <= xtol*xnorm) info = 2
            if (dabs(actred) <= ftol .and. prered <= ftol&
                .and. p5*ratio <= one .and. info == 2) info = 3
            if (info /= 0) go to 340
!
!           tests for termination and stringent tolerances.
!
            if (nfev >= maxfev) info = 5
            if (dabs(actred) <= epsmch .and. prered <= epsmch&
                .and. p5*ratio <= one) info = 6
            if (delta <= epsmch*xnorm) info = 7
            if (gnorm <= epsmch) info = 8
            if (info /= 0) go to 340
!
!           end of the inner loop. repeat if iteration unsuccessful.
!
            if (ratio < p0001) go to 240
!
!        end of the outer loop.
!
         go to 30
  340 continue
!
!     termination, either normal or user imposed.
!
      if (iflag < 0) info = iflag
      iflag = 0
      if (nprint > 0) call fcn(m,n,x,fvec,wa3,iflag)
      return
!
!     last card of subroutine lmstr.
!
      end subroutine lmstr
      
      subroutine lmstr1(fcn,m,n,x,fvec,fjac,ldfjac,tol,info,ipvt,wa,&
                        lwa)
      integer m,n,ldfjac,info,lwa
      integer ipvt(n)
      double precision tol
      double precision x(n),fvec(m),fjac(ldfjac,n),wa(lwa)
      external fcn
!     **********
!
!     subroutine lmstr1
!
!     the purpose of lmstr1 is to minimize the sum of the squares of
!     m nonlinear functions in n variables by a modification of
!     the levenberg-marquardt algorithm which uses minimal storage.
!     this is done by using the more general least-squares solver
!     lmstr. the user must provide a subroutine which calculates
!     the functions and the rows of the jacobian.
!
!     the subroutine statement is
!
!       subroutine lmstr1(fcn,m,n,x,fvec,fjac,ldfjac,tol,info,
!                         ipvt,wa,lwa)
!
!     where
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions and the rows of the jacobian.
!         fcn must be declared in an external statement in the
!         user calling program, and should be written as follows.
!
!         subroutine fcn(m,n,x,fvec,fjrow,iflag)
!         integer m,n,iflag
!         double precision x(n),fvec(m),fjrow(n)
!         ----------
!         if iflag = 1 calculate the functions at x and
!         return this vector in fvec.
!         if iflag = i calculate the (i-1)-st row of the
!         jacobian at x and return this vector in fjrow.
!         ----------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of lmstr1.
!         in this case set iflag to a negative integer.
!
!       m is a positive integer input variable set to the number
!         of functions.
!
!       n is a positive integer input variable set to the number
!         of variables. n must not exceed m.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length m which contains
!         the functions evaluated at the output x.
!
!       fjac is an output n by n array. the upper triangle of fjac
!         contains an upper triangular matrix r such that
!
!                t     t           t
!               p *(jac *jac)*p = r *r,
!
!         where p is a permutation matrix and jac is the final
!         calculated jacobian. column j of p is column ipvt(j)
!         (see below) of the identity matrix. the lower triangular
!         part of fjac contains information generated during
!         the computation of r.
!
!       ldfjac is a positive integer input variable not less than n
!         which specifies the leading dimension of the array fjac.
!
!       tol is a nonnegative input variable. termination occurs
!         when the algorithm estimates either that the relative
!         error in the sum of squares is at most tol or that
!         the relative error between x and the solution is at
!         most tol.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0  improper input parameters.
!
!         info = 1  algorithm estimates that the relative error
!                   in the sum of squares is at most tol.
!
!         info = 2  algorithm estimates that the relative error
!                   between x and the solution is at most tol.
!
!         info = 3  conditions for info = 1 and info = 2 both hold.
!
!         info = 4  fvec is orthogonal to the columns of the
!                   jacobian to machine precision.
!
!         info = 5  number of calls to fcn with iflag = 1 has
!                   reached 100*(n+1).
!
!         info = 6  tol is too small. no further reduction in
!                   the sum of squares is possible.
!
!         info = 7  tol is too small. no further improvement in
!                   the approximate solution x is possible.
!
!       ipvt is an integer output array of length n. ipvt
!         defines a permutation matrix p such that jac*p = q*r,
!         where jac is the final calculated jacobian, q is
!         orthogonal (not stored), and r is upper triangular.
!         column j of p is column ipvt(j) of the identity matrix.
!
!       wa is a work array of length lwa.
!
!       lwa is a positive integer input variable not less than 5*n+m.
!
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... lmstr
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, dudley v. goetschel, kenneth e. hillstrom,
!     jorge j. more
!
!     **********
      integer maxfev,mode,nfev,njev,nprint
      double precision factor,ftol,gtol,xtol,zero
      data factor,zero /1.0d2,0.0d0/
      info = 0
!
!     check the input parameters for errors.
!
      if (n <= 0 .or. m < n .or. ldfjac < n .or. tol < zero&
          .or. lwa < 5*n + m) go to 10
!
!     call lmstr.
!
      maxfev = 100*(n + 1)
      ftol = tol
      xtol = tol
      gtol = zero
      mode = 1
      nprint = 0
      call lmstr(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,maxfev,&
                 wa(1),mode,factor,nprint,info,nfev,njev,ipvt,wa(n+1),&
                 wa(2*n+1),wa(3*n+1),wa(4*n+1),wa(5*n+1))
      if (info == 8) info = 4
   10 continue
      return
!
!     last card of subroutine lmstr1.
!
      end subroutine lmstr1
      
      subroutine qform(m,n,q,ldq,wa)
      integer m,n,ldq
      double precision q(ldq,m),wa(m)
!     **********
!
!     subroutine qform
!
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
!
!     subprograms called
!
!       fortran-supplied ... min0
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      integer i,j,jm1,k,l,minmn,np1
      double precision one,sum,temp,zero
      data one,zero /1.0d0,0.0d0/
!
!     zero out upper triangle of q in the first min(m,n) columns.
!
      minmn = min0(m,n)
      if (minmn < 2) go to 30
      do 20 j = 2, minmn
         jm1 = j - 1
         do 10 i = 1, jm1
            q(i,j) = zero
   10       continue
   20    continue
   30 continue
!
!     initialize remaining columns to those of the identity matrix.
!
      np1 = n + 1
      if (m < np1) go to 60
      do 50 j = np1, m
         do 40 i = 1, m
            q(i,j) = zero
   40       continue
         q(j,j) = one
   50    continue
   60 continue
!
!     accumulate q from its factored form.
!
      do 120 l = 1, minmn
         k = minmn - l + 1
         do 70 i = k, m
            wa(i) = q(i,k)
            q(i,k) = zero
   70       continue
         q(k,k) = one
         if (wa(k) == zero) go to 110
         do 100 j = k, m
            sum = zero
            do 80 i = k, m
               sum = sum + q(i,j)*wa(i)
   80          continue
            temp = sum/wa(k)
            do 90 i = k, m
               q(i,j) = q(i,j) - temp*wa(i)
   90          continue
  100       continue
  110    continue
  120    continue
      return
!
!     last card of subroutine qform.
!
      end subroutine qform
      
      subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
      integer m,n,lda,lipvt
      integer ipvt(lipvt)
      logical pivot
      double precision a(lda,n),rdiag(n),acnorm(n),wa(n)
!     **********
!
!     subroutine qrfac
!
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
!
!     subprograms called
!
!       minpack-supplied ... dpmpar,enorm
!
!       fortran-supplied ... dmax1,dsqrt,min0
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      integer i,j,jp1,k,kmax,minmn
      double precision ajnorm,epsmch,one,p05,sum,temp,zero
      data one,p05,zero /1.0d0,5.0d-2,0.0d0/
!
!     epsmch is the machine precision.
!
      epsmch = dpmpar(1)
!
!     compute the initial column norms and initialize several arrays.
!
      do 10 j = 1, n
         acnorm(j) = enorm(m,a(1,j))
         rdiag(j) = acnorm(j)
         wa(j) = rdiag(j)
         if (pivot) ipvt(j) = j
   10    continue
!
!     reduce a to r with householder transformations.
!
      minmn = min0(m,n)
      do 110 j = 1, minmn
         if (.not.pivot) go to 40
!
!        bring the column of largest norm into the pivot position.
!
         kmax = j
         do 20 k = j, n
            if (rdiag(k) > rdiag(kmax)) kmax = k
   20       continue
         if (kmax == j) go to 40
         do 30 i = 1, m
            temp = a(i,j)
            a(i,j) = a(i,kmax)
            a(i,kmax) = temp
   30       continue
         rdiag(kmax) = rdiag(j)
         wa(kmax) = wa(j)
         k = ipvt(j)
         ipvt(j) = ipvt(kmax)
         ipvt(kmax) = k
   40    continue
!
!        compute the householder transformation to reduce the
!        j-th column of a to a multiple of the j-th unit vector.
!
         ajnorm = enorm(m-j+1,a(j,j))
         if (ajnorm == zero) go to 100
         if (a(j,j) < zero) ajnorm = -ajnorm
         do 50 i = j, m
            a(i,j) = a(i,j)/ajnorm
   50       continue
         a(j,j) = a(j,j) + one
!
!        apply the transformation to the remaining columns
!        and update the norms.
!
         jp1 = j + 1
         if (n < jp1) go to 100
         do 90 k = jp1, n
            sum = zero
            do 60 i = j, m
               sum = sum + a(i,j)*a(i,k)
   60          continue
            temp = sum/a(j,j)
            do 70 i = j, m
               a(i,k) = a(i,k) - temp*a(i,j)
   70          continue
            if (.not.pivot .or. rdiag(k) == zero) go to 80
            temp = a(j,k)/rdiag(k)
            rdiag(k) = rdiag(k)*dsqrt(dmax1(zero,one-temp**2))
            if (p05*(rdiag(k)/wa(k))**2 > epsmch) go to 80
            rdiag(k) = enorm(m-j,a(jp1,k))
            wa(k) = rdiag(k)
   80       continue
   90       continue
  100    continue
         rdiag(j) = -ajnorm
  110    continue
      return
!
!     last card of subroutine qrfac.
!
      end subroutine qrfac
      
      subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)
      integer n,ldr
      integer ipvt(n)
      double precision r(ldr,n),diag(n),qtb(n),x(n),sdiag(n),wa(n)
!     **********
!
!     subroutine qrsolv
!
!     given an m by n matrix a, an n by n diagonal matrix d,
!     and an m-vector b, the problem is to determine an x which
!     solves the system
!
!           a*x = b ,     d*x = 0 ,
!
!     in the least squares sense.
!
!     this subroutine completes the solution of the problem
!     if it is provided with the necessary information from the
!     qr factorization, with column pivoting, of a. that is, if
!     a*p = q*r, where p is a permutation matrix, q has orthogonal
!     columns, and r is an upper triangular matrix with diagonal
!     elements of nonincreasing magnitude, then qrsolv expects
!     the full upper triangle of r, the permutation matrix p,
!     and the first n components of (q transpose)*b. the system
!     a*x = b, d*x = 0, is then equivalent to
!
!                  t       t
!           r*z = q *b ,  p *d*p*z = 0 ,
!
!     where x = p*z. if this system does not have full rank,
!     then a least squares solution is obtained. on output qrsolv
!     also provides an upper triangular matrix s such that
!
!            t   t               t
!           p *(a *a + d*d)*p = s *s .
!
!     s is computed within qrsolv and may be of separate interest.
!
!     the subroutine statement is
!
!       subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)
!
!     where
!
!       n is a positive integer input variable set to the order of r.
!
!       r is an n by n array. on input the full upper triangle
!         must contain the full upper triangle of the matrix r.
!         on output the full upper triangle is unaltered, and the
!         strict lower triangle contains the strict upper triangle
!         (transposed) of the upper triangular matrix s.
!
!       ldr is a positive integer input variable not less than n
!         which specifies the leading dimension of the array r.
!
!       ipvt is an integer input array of length n which defines the
!         permutation matrix p such that a*p = q*r. column j of p
!         is column ipvt(j) of the identity matrix.
!
!       diag is an input array of length n which must contain the
!         diagonal elements of the matrix d.
!
!       qtb is an input array of length n which must contain the first
!         n elements of the vector (q transpose)*b.
!
!       x is an output array of length n which contains the least
!         squares solution of the system a*x = b, d*x = 0.
!
!       sdiag is an output array of length n which contains the
!         diagonal elements of the upper triangular matrix s.
!
!       wa is a work array of length n.
!
!     subprograms called
!
!       fortran-supplied ... dabs,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      integer i,j,jp1,k,kp1,l,nsing
      double precision cos,cotan,p5,p25,qtbpj,sin,sum,tan,temp,zero
      data p5,p25,zero /5.0d-1,2.5d-1,0.0d0/
!
!     copy r and (q transpose)*b to preserve input and initialize s.
!     in particular, save the diagonal elements of r in x.
!
      do 20 j = 1, n
         do 10 i = j, n
            r(i,j) = r(j,i)
   10       continue
         x(j) = r(j,j)
         wa(j) = qtb(j)
   20    continue
!
!     eliminate the diagonal matrix d using a givens rotation.
!
      do 100 j = 1, n
!
!        prepare the row of d to be eliminated, locating the
!        diagonal element using p from the qr factorization.
!
         l = ipvt(j)
         if (diag(l) == zero) go to 90
         do 30 k = j, n
            sdiag(k) = zero
   30       continue
         sdiag(j) = diag(l)
!
!        the transformations to eliminate the row of d
!        modify only a single element of (q transpose)*b
!        beyond the first n, which is initially zero.
!
         qtbpj = zero
         do 80 k = j, n
!
!           determine a givens rotation which eliminates the
!           appropriate element in the current row of d.
!
            if (sdiag(k) == zero) go to 70
            if (dabs(r(k,k)) >= dabs(sdiag(k))) go to 40
               cotan = r(k,k)/sdiag(k)
               sin = p5/dsqrt(p25+p25*cotan**2)
               cos = sin*cotan
               go to 50
   40       continue
               tan = sdiag(k)/r(k,k)
               cos = p5/dsqrt(p25+p25*tan**2)
               sin = cos*tan
   50       continue
!
!           compute the modified diagonal element of r and
!           the modified element of ((q transpose)*b,0).
!
            r(k,k) = cos*r(k,k) + sin*sdiag(k)
            temp = cos*wa(k) + sin*qtbpj
            qtbpj = -sin*wa(k) + cos*qtbpj
            wa(k) = temp
!
!           accumulate the tranformation in the row of s.
!
            kp1 = k + 1
            if (n < kp1) go to 70
            do 60 i = kp1, n
               temp = cos*r(i,k) + sin*sdiag(i)
               sdiag(i) = -sin*r(i,k) + cos*sdiag(i)
               r(i,k) = temp
   60          continue
   70       continue
   80       continue
   90    continue
!
!        store the diagonal element of s and restore
!        the corresponding diagonal element of r.
!
         sdiag(j) = r(j,j)
         r(j,j) = x(j)
  100    continue
!
!     solve the triangular system for z. if the system is
!     singular, then obtain a least squares solution.
!
      nsing = n
      do 110 j = 1, n
         if (sdiag(j) == zero .and. nsing == n) nsing = j - 1
         if (nsing < n) wa(j) = zero
  110    continue
      if (nsing < 1) go to 150
      do 140 k = 1, nsing
         j = nsing - k + 1
         sum = zero
         jp1 = j + 1
         if (nsing < jp1) go to 130
         do 120 i = jp1, nsing
            sum = sum + r(i,j)*wa(i)
  120       continue
  130    continue
         wa(j) = (wa(j) - sum)/sdiag(j)
  140    continue
  150 continue
!
!     permute the components of z back to components of x.
!
      do 160 j = 1, n
         l = ipvt(j)
         x(l) = wa(j)
  160    continue
      return
!
!     last card of subroutine qrsolv.
!
      end subroutine qrsolv
      
      subroutine r1mpyq(m,n,a,lda,v,w)
      integer m,n,lda
      double precision a(lda,n),v(n),w(n)
!     **********
!
!     subroutine r1mpyq
!
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
!
!     subroutines called
!
!       fortran-supplied ... dabs,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      integer i,j,nmj,nm1
      double precision cos,one,sin,temp
      data one /1.0d0/
!
!     apply the first set of givens rotations to a.
!
      nm1 = n - 1
      if (nm1 < 1) go to 50
      do 20 nmj = 1, nm1
         j = n - nmj
         if (dabs(v(j)) > one) cos = one/v(j)
         if (dabs(v(j)) > one) sin = dsqrt(one-cos**2)
         if (dabs(v(j)) <= one) sin = v(j)
         if (dabs(v(j)) <= one) cos = dsqrt(one-sin**2)
         do 10 i = 1, m
            temp = cos*a(i,j) - sin*a(i,n)
            a(i,n) = sin*a(i,j) + cos*a(i,n)
            a(i,j) = temp
   10       continue
   20    continue
!
!     apply the second set of givens rotations to a.
!
      do 40 j = 1, nm1
         if (dabs(w(j)) > one) cos = one/w(j)
         if (dabs(w(j)) > one) sin = dsqrt(one-cos**2)
         if (dabs(w(j)) <= one) sin = w(j)
         if (dabs(w(j)) <= one) cos = dsqrt(one-sin**2)
         do 30 i = 1, m
            temp = cos*a(i,j) + sin*a(i,n)
            a(i,n) = -sin*a(i,j) + cos*a(i,n)
            a(i,j) = temp
   30       continue
   40    continue
   50 continue
      return
!
!     last card of subroutine r1mpyq.
!
      end subroutine r1mpyq
      
      subroutine r1updt(m,n,s,ls,u,v,w,sing)
      integer m,n,ls
      logical sing
      double precision s(ls),u(m),v(n),w(m)
!     **********
!
!     subroutine r1updt
!
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
!
!     subprograms called
!
!       minpack-supplied ... dpmpar
!
!       fortran-supplied ... dabs,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more,
!     john l. nazareth
!
!     **********
      integer i,j,jj,l,nmj,nm1
      double precision cos,cotan,giant,one,p5,p25,sin,tan,tau,temp,&
                       zero
      data one,p5,p25,zero /1.0d0,5.0d-1,2.5d-1,0.0d0/
!
!     giant is the largest magnitude.
!
      giant = dpmpar(3)
!
!     initialize the diagonal element pointer.
!
      jj = (n*(2*m - n + 1))/2 - (m - n)
!
!     move the nontrivial part of the last column of s into w.
!
      l = jj
      do 10 i = n, m
         w(i) = s(l)
         l = l + 1
   10    continue
!
!     rotate the vector v into a multiple of the n-th unit vector
!     in such a way that a spike is introduced into w.
!
      nm1 = n - 1
      if (nm1 < 1) go to 70
      do 60 nmj = 1, nm1
         j = n - nmj
         jj = jj - (m - j + 1)
         w(j) = zero
         if (v(j) == zero) go to 50
!
!        determine a givens rotation which eliminates the
!        j-th element of v.
!
         if (dabs(v(n)) >= dabs(v(j))) go to 20
            cotan = v(n)/v(j)
            sin = p5/dsqrt(p25+p25*cotan**2)
            cos = sin*cotan
            tau = one
            if (dabs(cos)*giant > one) tau = one/cos
            go to 30
   20    continue
            tan = v(j)/v(n)
            cos = p5/dsqrt(p25+p25*tan**2)
            sin = cos*tan
            tau = sin
   30    continue
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
         do 40 i = j, m
            temp = cos*s(l) - sin*w(i)
            w(i) = sin*s(l) + cos*w(i)
            s(l) = temp
            l = l + 1
   40       continue
   50    continue
   60    continue
   70 continue
!
!     add the spike from the rank 1 update to w.
!
      do 80 i = 1, m
         w(i) = w(i) + v(n)*u(i)
   80    continue
!
!     eliminate the spike.
!
      sing = .false.
      if (nm1 < 1) go to 140
      do 130 j = 1, nm1
         if (w(j) == zero) go to 120
!
!        determine a givens rotation which eliminates the
!        j-th element of the spike.
!
         if (dabs(s(jj)) >= dabs(w(j))) go to 90
            cotan = s(jj)/w(j)
            sin = p5/dsqrt(p25+p25*cotan**2)
            cos = sin*cotan
            tau = one
            if (dabs(cos)*giant > one) tau = one/cos
            go to 100
   90    continue
            tan = w(j)/s(jj)
            cos = p5/dsqrt(p25+p25*tan**2)
            sin = cos*tan
            tau = sin
  100    continue
!
!        apply the transformation to s and reduce the spike in w.
!
         l = jj
         do 110 i = j, m
            temp = cos*s(l) + sin*w(i)
            w(i) = -sin*s(l) + cos*w(i)
            s(l) = temp
            l = l + 1
  110       continue
!
!        store the information necessary to recover the
!        givens rotation.
!
         w(j) = tau
  120    continue
!
!        test for zero diagonal elements in the output s.
!
         if (s(jj) == zero) sing = .true.
         jj = jj + (m - j + 1)
  130    continue
  140 continue
!
!     move w back into the last column of the output s.
!
      l = jj
      do 150 i = n, m
         s(l) = w(i)
         l = l + 1
  150    continue
      if (s(jj) == zero) sing = .true.
      return
!
!     last card of subroutine r1updt.
!
      end subroutine r1updt
      
      subroutine rwupdt(n,r,ldr,w,b,alpha,cos,sin)
      integer n,ldr
      double precision alpha
      double precision r(ldr,n),w(n),b(n),cos(n),sin(n)
!     **********
!
!     subroutine rwupdt
!
!     given an n by n upper triangular matrix r, this subroutine
!     computes the qr decomposition of the matrix formed when a row
!     is added to r. if the row is specified by the vector w, then
!     rwupdt determines an orthogonal matrix q such that when the
!     n+1 by n matrix composed of r augmented by w is premultiplied
!     by (q transpose), the resulting matrix is upper trapezoidal.
!     the matrix (q transpose) is the product of n transformations
!
!           g(n)*g(n-1)* ... *g(1)
!
!     where g(i) is a givens rotation in the (i,n+1) plane which
!     eliminates elements in the (n+1)-st plane. rwupdt also
!     computes the product (q transpose)*c where c is the
!     (n+1)-vector (b,alpha). q itself is not accumulated, rather
!     the information to recover the g rotations is supplied.
!
!     the subroutine statement is
!
!       subroutine rwupdt(n,r,ldr,w,b,alpha,cos,sin)
!
!     where
!
!       n is a positive integer input variable set to the order of r.
!
!       r is an n by n array. on input the upper triangular part of
!         r must contain the matrix to be updated. on output r
!         contains the updated triangular matrix.
!
!       ldr is a positive integer input variable not less than n
!         which specifies the leading dimension of the array r.
!
!       w is an input array of length n which must contain the row
!         vector to be added to r.
!
!       b is an array of length n. on input b must contain the
!         first n elements of the vector c. on output b contains
!         the first n elements of the vector (q transpose)*c.
!
!       alpha is a variable. on input alpha must contain the
!         (n+1)-st element of the vector c. on output alpha contains
!         the (n+1)-st element of the vector (q transpose)*c.
!
!       cos is an output array of length n which contains the
!         cosines of the transforming givens rotations.
!
!       sin is an output array of length n which contains the
!         sines of the transforming givens rotations.
!
!     subprograms called
!
!       fortran-supplied ... dabs,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, dudley v. goetschel, kenneth e. hillstrom,
!     jorge j. more
!
!     **********
      integer i,j,jm1
      double precision cotan,one,p5,p25,rowj,tan,temp,zero
      data one,p5,p25,zero /1.0d0,5.0d-1,2.5d-1,0.0d0/
!
      do 60 j = 1, n
         rowj = w(j)
         jm1 = j - 1
!
!        apply the previous transformations to
!        r(i,j), i=1,2,...,j-1, and to w(j).
!
         if (jm1 < 1) go to 20
         do 10 i = 1, jm1
            temp = cos(i)*r(i,j) + sin(i)*rowj
            rowj = -sin(i)*r(i,j) + cos(i)*rowj
            r(i,j) = temp
   10       continue
   20    continue
!
!        determine a givens rotation which eliminates w(j).
!
         cos(j) = one
         sin(j) = zero
         if (rowj == zero) go to 50
         if (dabs(r(j,j)) >= dabs(rowj)) go to 30
            cotan = r(j,j)/rowj
            sin(j) = p5/dsqrt(p25+p25*cotan**2)
            cos(j) = sin(j)*cotan
            go to 40
   30    continue
            tan = rowj/r(j,j)
            cos(j) = p5/dsqrt(p25+p25*tan**2)
            sin(j) = cos(j)*tan
   40    continue
!
!        apply the current transformation to r(j,j), b(j), and alpha.
!
         r(j,j) = cos(j)*r(j,j) + sin(j)*rowj
         temp = cos(j)*b(j) + sin(j)*alpha
         alpha = -sin(j)*b(j) + cos(j)*alpha
         b(j) = temp
   50    continue
   60    continue
      return
!
!     last card of subroutine rwupdt.
!
      end subroutine rwupdt
      
    end module minpack_module


!    Minpack Copyright Notice (1999) University of Chicago.  All rights reserved
!    
!    Redistribution and use in source and binary forms, with or
!    without modification, are permitted provided that the
!    following conditions are met:
!    
!    1. Redistributions of source code must retain the above
!    copyright notice, this list of conditions and the following
!    disclaimer.
!    
!    2. Redistributions in binary form must reproduce the above
!    copyright notice, this list of conditions and the following
!    disclaimer in the documentation and/or other materials
!    provided with the distribution.
!    
!    3. The end-user documentation included with the
!    redistribution, if any, must include the following
!    acknowledgment:
!    
!       "This product includes software developed by the
!       University of Chicago, as Operator of Argonne National
!       Laboratory.
!    
!    Alternately, this acknowledgment may appear in the software
!    itself, if and wherever such third-party acknowledgments
!    normally appear.
!    
!    4. WARRANTY DISCLAIMER. THE SOFTWARE IS SUPPLIED "AS IS"
!    WITHOUT WARRANTY OF ANY KIND. THE COPYRIGHT HOLDER, THE
!    UNITED STATES, THE UNITED STATES DEPARTMENT OF ENERGY, AND
!    THEIR EMPLOYEES: (1) DISCLAIM ANY WARRANTIES, EXPRESS OR
!    IMPLIED, INCLUDING BUT NOT LIMITED TO ANY IMPLIED WARRANTIES
!    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE
!    OR NON-INFRINGEMENT, (2) DO NOT ASSUME ANY LEGAL LIABILITY
!    OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR
!    USEFULNESS OF THE SOFTWARE, (3) DO NOT REPRESENT THAT USE OF
!    THE SOFTWARE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS, (4)
!    DO NOT WARRANT THAT THE SOFTWARE WILL FUNCTION
!    UNINTERRUPTED, THAT IT IS ERROR-FREE OR THAT ANY ERRORS WILL
!    BE CORRECTED.
!    
!    5. LIMITATION OF LIABILITY. IN NO EVENT WILL THE COPYRIGHT
!    HOLDER, THE UNITED STATES, THE UNITED STATES DEPARTMENT OF
!    ENERGY, OR THEIR EMPLOYEES: BE LIABLE FOR ANY INDIRECT,
!    INCIDENTAL, CONSEQUENTIAL, SPECIAL OR PUNITIVE DAMAGES OF
!    ANY KIND OR NATURE, INCLUDING BUT NOT LIMITED TO LOSS OF
!    PROFITS OR LOSS OF DATA, FOR ANY REASON WHATSOEVER, WHETHER
!    SUCH LIABILITY IS ASSERTED ON THE BASIS OF CONTRACT, TORT
!    (INCLUDING NEGLIGENCE OR STRICT LIABILITY), OR OTHERWISE,
!    EVEN IF ANY OF SAID PARTIES HAS BEEN WARNED OF THE
!    POSSIBILITY OF SUCH LOSS OR DAMAGES.

