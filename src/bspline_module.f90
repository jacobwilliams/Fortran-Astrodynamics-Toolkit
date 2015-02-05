!*****************************************************************************************
    module bspline_module
!*****************************************************************************************
!****h* FAT/bspline_module
!
!  NAME
!    bspline_module
!
!  DESCRIPTION
!    B-Spline routines.
!
!    Currently this module contains the routines from [1], 
!    with minimal conversion to free-format source code.
!
!    More refactoring will be done...
!
!  SEE ALSO
!    [1] DBSPLIN from the NIST Core Math Library (CMLIB)
!        http://www.nist.gov/itl/math/mcsd-software.cfm
!
!*****************************************************************************************
    
    implicit none
    
    contains
    
!*****************************************************************************************

      subroutine dbfqad(f,t,bcoef,n,k,id,x1,x2,tol,quad,ierr,work)
!***begin prologue  dbfqad
!***date written   800901   (yymmdd)
!***revision date  820801   (yymmdd)
!***revision history  (yymmdd)
!   000330  modified array declarations.  (jec)
!
!***category no.  h2a2a1,e3,k6
!***keywords  b-spline,data fitting,double precision,interpolation,
!             quadrature,spline
!***author  amos, d. e., (snla)
!***purpose  computes the integral on (x1,x2) of a product of a function
!            f and the id-th derivative of a k-th order b-spline
!            (b-representation).
!***description
!
!     written by d. e. amos, june, 1979.
!
!     reference
!         sand-79-1825
!
!     abstract    **** a double precision routine ****
!
!         dbfqad computes the integral on (x1,x2) of a product of a
!         function f and the id-th derivative of a k-th order b-spline,
!         using the b-representation (t,bcoef,n,k).  (x1,x2) must be a
!         subinterval of t(k) <= x <= t(n+1).  an integration rou-
!         tine, dbsgq8 (a modification of gaus8), integrates the product
!         on subintervals of (x1,x2) formed by included (distinct) knots
!
!         the maximum number of significant digits obtainable in
!         dbsqad is the smaller of 18 and the number of digits
!         carried in double precision arithmetic.
!
!         dbfqad calls dintrv, dbvalu, dbsgq8, d1mach, xerror
!
!     description of arguments
!         input      f,t,bcoef,x1,x2,tol are double precision
!           f      - external function of one argument for the
!                    integrand bf(x)=f(x)*dbvalu(t,bcoef,n,k,id,x,inbv,
!                    work)
!           t      - knot array of length n+k
!           bcoef  - coefficient array of length n
!           n      - length of coefficient array
!           k      - order of b-spline, k >= 1
!           id     - order of the spline derivative, 0 <= id <= k-1
!                    id=0 gives the spline function
!           x1,x2  - end points of quadrature interval in
!                    t(k) <= x <= t(n+1)
!           tol    - desired accuracy for the quadrature, suggest
!                    10.*dtol < tol <= .1 where dtol is the maximum
!                    of 1.0d-18 and double precision unit roundoff for
!                    the machine = d1mach(4)
!
!         output     quad,work are double precision
!           quad   - integral of bf(x) on (x1,x2)
!           ierr   - a status code
!                    ierr=1  normal return
!                         2  some quadrature on (x1,x2) does not meet
!                            the requested tolerance.
!           work   - work vector of length 3*k
!
!     error conditions
!         improper input is a fatal error
!         some quadrature fails to meet the requested tolerance
!***references  d.e. amos, *quadrature subroutines for splines and
!                 b-splines*, sand79-1825, sandia laboratories,
!                 december 1979.
!***routines called  d1mach,dbsgq8,dintrv,xerror
!***end prologue  dbfqad
!
!
      integer id, ierr, iflg, ilo, il1, il2, k, left, mflag, n, npk, np1
      double precision a,aa,ans,b,bb,bcoef,q,quad,t,ta,tb,tol,work,wtol,&
       x1, x2
      double precision f
      dimension t(*), bcoef(*), work(*)
      integer inbv
      external f
!***first executable statement  dbfqad
      ierr = 1
      quad = 0.0d0
      if(k<1) go to 100
      if(n<k) go to 105
      if(id<0 .or. id>=k) go to 110
      wtol = d1mach(4)
      wtol = dmax1(wtol,1.d-18)
      if (tol<wtol .or. tol>0.1d0) go to 30
      aa = dmin1(x1,x2)
      bb = dmax1(x1,x2)
      if (aa<t(k)) go to 20
      np1 = n + 1
      if (bb>t(np1)) go to 20
      if (aa==bb) return
      npk = n + k
!
      ilo = 1
      call dintrv(t, npk, aa, ilo, il1, mflag)
      call dintrv(t, npk, bb, ilo, il2, mflag)
      if (il2>=np1) il2 = n
      inbv = 1
      q = 0.0d0
      do 10 left=il1,il2
        ta = t(left)
        tb = t(left+1)
        if (ta==tb) go to 10
        a = dmax1(aa,ta)
        b = dmin1(bb,tb)
        call dbsgq8(f,t,bcoef,n,k,id,a,b,inbv,tol,ans,iflg,work)
        if (iflg>1) ierr = 2
        q = q + ans
   10 continue
      if (x1>x2) q = -q
      quad = q
      return
!
!
   20 continue
      call xerror( ' dbfqad,  x1 or x2 or both do not satisfy t(k)<=x.&
      le.t(n+1)',  61, 2, 1)
      return
   30 continue
      call xerror( ' dbfqad,  tol is less dtol or greater than 0.1',&
       46, 2, 1)
      return
  100 continue
      call xerror( ' dbfqad,  k does not satisfy k>=1', 35, 2, 1)
      return
  105 continue
      call xerror( ' dbfqad,  n does not satisfy n>=k', 35, 2, 1)
      return
  110 continue
      call xerror( ' dbfqad,  id does not satisfy 0<=id<k',&
       42, 2, 1)
      return
      end subroutine dbfqad
      
      subroutine dbint4(x,y,ndata,ibcl,ibcr,fbcl,fbcr,kntopt,t,bcoef,n,&
         k,w)
!***begin prologue  dbint4
!***date written   800901   (yymmdd)
!***revision date  820801   (yymmdd)
!***revision history  (yymmdd)
!   000330  modified array declarations.  (jec)
!
!***category no.  e1a
!***keywords  b-spline,data fitting,double precision,interpolation,
!             spline
!***author  amos, d. e., (snla)
!***purpose  computes the b-representation of a cubic spline
!            which interpolates data (x(i),y(i)),i=1,ndata.
!***description
!
!     written by d. e. amos, august, 1979.
!
!     references
!         sand-78-1968
!
!         a practical guide to splines by c. de boor, applied
!         mathematics series 27, springer, 1979.
!
!         siam j. numerical analysis, 14, no. 3, june, 1977, pp.441-472.
!
!     abstract    **** a double precision routine ****
!
!         dbint4 computes the b representation (t,bcoef,n,k) of a
!         cubic spline (k=4) which interpolates data (x(i),y(i)),
!         i=1,ndata.  parameters ibcl, ibcr, fbcl, fbcr allow the
!         specification of the spline first or second derivative at
!         both x(1) and x(ndata).  when this data is not specified
!         by the problem, it is common practice to use a natural
!         spline by setting second derivatives at x(1) and x(ndata)
!         to zero (ibcl=ibcr=2,fbcl=fbcr=0.0).  the spline is defined
!         on t(4) <= x <= t(n+1) with (ordered) interior knots at
!         x(i) values where n=ndata+2.  the knots t(1),t(2),t(3) lie to
!         the left of t(4)=x(1) and the knots t(n+2), t(n+3), t(n+4)
!         lie to the right of t(n+1)=x(ndata) in increasing order.  if
!         no extrapolation outside (x(1),x(ndata)) is anticipated, the
!         knots t(1)=t(2)=t(3)=t(4)=x(1) and t(n+2)=t(n+3)=t(n+4)=
!         t(n+1)=x(ndata) can be specified by kntopt=1.  kntopt=2
!         selects a knot placement for t(1), t(2), t(3) to make the
!         first 7 knots symmetric about t(4)=x(1) and similarly for
!         t(n+2), t(n+3), t(n+4) about t(n+1)=x(ndata).  kntopt=3
!         allows the user to make his own selection, in increasing
!         order, for t(1), t(2), t(3) to the left of x(1) and t(n+2),
!         t(n+3), t(n+4) to the right of x(ndata) in the work array
!         w(1) through w(6).  in any case, the interpolation on
!         t(4) <= x <= t(n+1) by using function dbvalu is unique
!         for given boundary conditions.
!
!         dbint4 calls dbspvd, dbnfac, dbnslv, d1mach, xerror
!
!     description of arguments
!
!         input      x,y,fbcl,fbcr,w are double precision
!           x      - x vector of abscissae of length ndata, distinct
!                    and in increasing order
!           y      - y vector of ordinates of length ndata
!           ndata  - number of data points, ndata >= 2
!           ibcl   - selection parameter for left boundary condition
!                    ibcl = 1 constrain the first derivative at
!                             x(1) to fbcl
!                         = 2 constrain the second derivative at
!                             x(1) to fbcl
!           ibcr   - selection parameter for right boundary condition
!                    ibcr = 1 constrain first derivative at
!                             x(ndata) to fbcr
!                    ibcr = 2 constrain second derivative at
!                             x(ndata) to fbcr
!           fbcl   - left boundary values governed by ibcl
!           fbcr   - right boundary values governed by ibcr
!           kntopt - knot selection parameter
!                    kntopt = 1 sets knot multiplicity at t(4) and
!                               t(n+1) to 4
!                           = 2 sets a symmetric placement of knots
!                               about t(4) and t(n+1)
!                           = 3 sets t(i)=w(i) and t(n+1+i)=w(3+i),i=1,3
!                               where w(i),i=1,6 is supplied by the user
!           w      - work array of dimension at least 5*(ndata+2)
!                    if kntopt=3, then w(1),w(2),w(3) are knot values to
!                    the left of x(1) and w(4),w(5),w(6) are knot
!                    values to the right of x(ndata) in increasing
!                    order to be supplied by the user
!
!         output     t,bcoef are double precision
!           t      - knot array of length n+4
!           bcoef  - b spline coefficient array of length n
!           n      - number of coefficients, n=ndata+2
!           k      - order of spline, k=4
!
!     error conditions
!         improper  input is a fatal error
!         singular system of equations is a fatal error
!***references  d.e. amos, *computation with splines and b-splines*,
!                 sand78-1968, sandia laboratories, march 1979.
!               c. de boor, *package for calculating with b-splines*,
!                 siam journal on numerical analysis, volume 14, no. 3,
!                 june 1977, pp. 441-472.
!               c. de boor, *a practical guide to splines*, applied
!                 mathematics series 27, springer, 1979.
!***routines called  d1mach,dbnfac,dbnslv,dbspvd,xerror
!***end prologue  dbint4
!
!
      integer i, ibcl, ibcr, iflag, ilb, ileft, it, iub, iw, iwp, j,&
       jw, k, kntopt, n, ndata, ndm, np, nwrow
      double precision bcoef,fbcl,fbcr,t,tol,txn,tx1,vnikx,w,wdtol,&
       work,x,xl,y
      dimension x(*), y(*), t(*), bcoef(*), w(5,*), vnikx(4,4), work(15)
!***first executable statement  dbint4
      wdtol = d1mach(4)
      tol = dsqrt(wdtol)
      if (ndata<2) go to 200
      ndm = ndata - 1
      do 10 i=1,ndm
        if (x(i)>=x(i+1)) go to 210
   10 continue
      if (ibcl<1 .or. ibcl>2) go to 220
      if (ibcr<1 .or. ibcr>2) go to 230
      if (kntopt<1 .or. kntopt>3) go to 240
      k = 4
      n = ndata + 2
      np = n + 1
      do 20 i=1,ndata
        t(i+3) = x(i)
   20 continue
      go to (30, 50, 90), kntopt
!     set up knot array with multiplicity 4 at x(1) and x(ndata)
   30 continue
      do 40 i=1,3
        t(4-i) = x(1)
        t(np+i) = x(ndata)
   40 continue
      go to 110
!     set up knot array with symmetric placement about end points
   50 continue
      if (ndata>3) go to 70
      xl = (x(ndata)-x(1))/3.0d0
      do 60 i=1,3
        t(4-i) = t(5-i) - xl
        t(np+i) = t(np+i-1) + xl
   60 continue
      go to 110
   70 continue
      tx1 = x(1) + x(1)
      txn = x(ndata) + x(ndata)
      do 80 i=1,3
        t(4-i) = tx1 - x(i+1)
        t(np+i) = txn - x(ndata-i)
   80 continue
      go to 110
!     set up knot array less than x(1) and greater than x(ndata) to be
!     supplied by user in work locations w(1) through w(6) when kntopt=3
   90 continue
      do 100 i=1,3
        t(4-i) = w(4-i,1)
        jw = max0(1,i-1)
        iw = mod(i+2,5)+1
        t(np+i) = w(iw,jw)
        if (t(4-i)>t(5-i)) go to 250
        if (t(np+i)<t(np+i-1)) go to 250
  100 continue
  110 continue
!
      do 130 i=1,5
        do 120 j=1,n
          w(i,j) = 0.0d0
  120   continue
  130 continue
!     set up left interpolation point and left boundary condition for
!     right limits
      it = ibcl + 1
      call dbspvd(t, k, it, x(1), k, 4, vnikx, work)
      iw = 0
      if (dabs(vnikx(3,1))<tol) iw = 1
      do 140 j=1,3
        w(j+1,4-j) = vnikx(4-j,it)
        w(j,4-j) = vnikx(4-j,1)
  140 continue
      bcoef(1) = y(1)
      bcoef(2) = fbcl
!     set up interpolation equations for points i=2 to i=ndata-1
      ileft = 4
      if (ndm<2) go to 170
      do 160 i=2,ndm
        ileft = ileft + 1
        call dbspvd(t, k, 1, x(i), ileft, 4, vnikx, work)
        do 150 j=1,3
          w(j+1,3+i-j) = vnikx(4-j,1)
  150   continue
        bcoef(i+1) = y(i)
  160 continue
!     set up right interpolation point and right boundary condition for
!     left limits(ileft is associated with t(n)=x(ndata-1))
  170 continue
      it = ibcr + 1
      call dbspvd(t, k, it, x(ndata), ileft, 4, vnikx, work)
      jw = 0
      if (dabs(vnikx(2,1))<tol) jw = 1
      do 180 j=1,3
        w(j+1,3+ndata-j) = vnikx(5-j,it)
        w(j+2,3+ndata-j) = vnikx(5-j,1)
  180 continue
      bcoef(n-1) = fbcr
      bcoef(n) = y(ndata)
!     solve system of equations
      ilb = 2 - jw
      iub = 2 - iw
      nwrow = 5
      iwp = iw + 1
      call dbnfac(w(iwp,1), nwrow, n, ilb, iub, iflag)
      if (iflag==2) go to 190
      call dbnslv(w(iwp,1), nwrow, n, ilb, iub, bcoef)
      return
!
!
  190 continue
      call xerror( ' dbint4,  the system of equations is singular',&
       45, 2, 1)
      return
  200 continue
      call xerror( ' dbint4,  ndata is less than 2', 30, 2, 1)
      return
  210 continue
      call xerror( ' dbint4,  x values are not distinct or not ordered',&
       50, 2, 1)
      return
  220 continue
      call xerror( ' dbint4,  ibcl is not 1 or 2', 28, 2, 1)
      return
  230 continue
      call xerror( ' dbint4,  ibcr is not 1 or 2', 28, 2, 1)
      return
  240 continue
      call xerror( ' dbint4,  kntopt is not 1, 2, or 3', 34, 2, 1)
      return
  250 continue
      call xerror(  ' dbint4,  knot input through w array is not ordered&
       properly',  60, 2, 1)
      return
      end subroutine dbint4
      
      subroutine dbintk(x,y,t,n,k,bcoef,q,work)
!***begin prologue  dbintk
!***date written   800901   (yymmdd)
!***revision date  820801   (yymmdd)
!***revision history  (yymmdd)
!   000330  modified array declarations.  (jec)
!
!***category no.  e1a
!***keywords  b-spline,data fitting,double precision,interpolation,
!             spline
!***author  amos, d. e., (snla)
!***purpose  produces the b-spline coefficients, bcoef, of the
!            b-spline of order k with knots t(i), i=1,...,n+k, which
!            takes on the value y(i) at x(i), i=1,...,n.
!***description
!
!     written by carl de boor and modified by d. e. amos
!
!     references
!
!         a practical guide to splines by c. de boor, applied
!         mathematics series 27, springer, 1979.
!
!     abstract    **** a double precision routine ****
!
!         dbintk is the splint routine of the reference.
!
!         dbintk produces the b-spline coefficients, bcoef, of the
!         b-spline of order k with knots t(i), i=1,...,n+k, which
!         takes on the value y(i) at x(i), i=1,...,n.  the spline or
!         any of its derivatives can be evaluated by calls to dbvalu.
!
!         the i-th equation of the linear system a*bcoef = b for the
!         coefficients of the interpolant enforces interpolation at
!         x(i), i=1,...,n.  hence, b(i) = y(i), for all i, and a is
!         a band matrix with 2k-1 bands if a is invertible.  the matrix
!         a is generated row by row and stored, diagonal by diagonal,
!         in the rows of q, with the main diagonal going into row k.
!         the banded system is then solved by a call to dbnfac (which
!         constructs the triangular factorization for a and stores it
!         again in q), followed by a call to dbnslv (which then
!         obtains the solution bcoef by substitution).  dbnfac does no
!         pivoting, since the total positivity of the matrix a makes
!         this unnecessary.  the linear system to be solved is
!         (theoretically) invertible if and only if
!                 t(i) < x(i) < t(i+k),        for all i.
!         equality is permitted on the left for i=1 and on the right
!         for i=n when k knots are used at x(1) or x(n).  otherwise,
!         violation of this condition is certain to lead to an error.
!
!         dbintk calls dbspvn, dbnfac, dbnslv, xerror
!
!     description of arguments
!
!         input       x,y,t are double precision
!           x       - vector of length n containing data point abscissa
!                     in strictly increasing order.
!           y       - corresponding vector of length n containing data
!                     point ordinates.
!           t       - knot vector of length n+k
!                     since t(1),..,t(k) <= x(1) and t(n+1),..,t(n+k)
!                     >= x(n), this leaves only n-k knots (not nec-
!                     essarily x(i) values) interior to (x(1),x(n))
!           n       - number of data points, n >= k
!           k       - order of the spline, k >= 1
!
!         output      bcoef,q,work are double precision
!           bcoef   - a vector of length n containing the b-spline
!                     coefficients
!           q       - a work vector of length (2*k-1)*n, containing
!                     the triangular factorization of the coefficient
!                     matrix of the linear system being solved.  the
!                     coefficients for the interpolant of an
!                     additional data set (x(i),yy(i)), i=1,...,n
!                     with the same abscissa can be obtained by loading
!                     yy into bcoef and then executing
!                         call dbnslv(q,2k-1,n,k-1,k-1,bcoef)
!           work    - work vector of length 2*k
!
!     error conditions
!         improper input is a fatal error
!         singular system of equations is a fatal error
!***references  c. de boor, *a practical guide to splines*, applied
!                 mathematics series 27, springer, 1979.
!               d.e. amos, *computation with splines and b-splines*,
!                 sand78-1968,sandia laboratories,march,1979.
!***routines called  dbnfac,dbnslv,dbspvn,xerror
!***end prologue  dbintk
!
!
      integer iflag, iwork, k, n, i, ilp1mx, j, jj, km1, kpkm2, left,&
       lenq, np1
      double precision bcoef(n), y(n), q(*), t(*), x(n), xi, work(*)
!     dimension q(2*k-1,n), t(n+k)
!***first executable statement  dbintk
      if(k<1) go to 100
      if(n<k) go to 105
      jj = n - 1
      if(jj==0) go to 6
      do 5 i=1,jj
      if(x(i)>=x(i+1)) go to 110
    5 continue
    6 continue
      np1 = n + 1
      km1 = k - 1
      kpkm2 = 2*km1
      left = k
!                zero out all entries of q
      lenq = n*(k+km1)
      do 10 i=1,lenq
        q(i) = 0.0d0
   10 continue
!
!  ***   loop over i to construct the  n  interpolation equations
      do 50 i=1,n
        xi = x(i)
        ilp1mx = min0(i+k,np1)
!        *** find  left  in the closed interval (i,i+k-1) such that
!                t(left) <= x(i) < t(left+1)
!        matrix is singular if this is not possible
        left = max0(left,i)
        if (xi<t(left)) go to 80
   20   if (xi<t(left+1)) go to 30
        left = left + 1
        if (left<ilp1mx) go to 20
        left = left - 1
        if (xi>t(left+1)) go to 80
!        *** the i-th equation enforces interpolation at xi, hence
!        a(i,j) = b(j,k,t)(xi), all j. only the  k  entries with  j =
!        left-k+1,...,left actually might be nonzero. these  k  numbers
!        are returned, in  bcoef (used for temp.storage here), by the
!        following
   30   call dbspvn(t, k, k, 1, xi, left, bcoef, work, iwork)
!        we therefore want  bcoef(j) = b(left-k+j)(xi) to go into
!        a(i,left-k+j), i.e., into  q(i-(left+j)+2*k,(left+j)-k) since
!        a(i+j,j)  is to go into  q(i+k,j), all i,j,  if we consider  q
!        as a two-dim. array , with  2*k-1  rows (see comments in
!        dbnfac). in the present program, we treat  q  as an equivalent
!        one-dimensional array (because of fortran restrictions on
!        dimension statements) . we therefore want  bcoef(j) to go into
!        entry
!            i -(left+j) + 2*k + ((left+j) - k-1)*(2*k-1)
!                   =  i-left+1 + (left -k)*(2*k-1) + (2*k-2)*j
!        of  q .
        jj = i - left + 1 + (left-k)*(k+km1)
        do 40 j=1,k
          jj = jj + kpkm2
          q(jj) = bcoef(j)
   40   continue
   50 continue
!
!     ***obtain factorization of  a  , stored again in  q.
      call dbnfac(q, k+km1, n, km1, km1, iflag)
      go to (60, 90), iflag
!     *** solve  a*bcoef = y  by backsubstitution
   60 do 70 i=1,n
        bcoef(i) = y(i)
   70 continue
      call dbnslv(q, k+km1, n, km1, km1, bcoef)
      return
!
!
   80 continue
      call xerror( ' dbintk,  some abscissa was not in the support of th&
      e corresponding basis function and the system is singular.',109,2,&
      1)
      return
   90 continue
      call xerror( ' dbintk,  the system of solver detects a singular sy&
      stem although the theoretical conditions for a solution were satis&
      fied.',123,8,1)
      return
  100 continue
      call xerror( ' dbintk,  k does not satisfy k>=1', 35, 2, 1)
      return
  105 continue
      call xerror( ' dbintk,  n does not satisfy n>=k', 35, 2, 1)
      return
  110 continue
      call xerror( ' dbintk,  x(i) does not satisfy x(i)<x(i+1) for s&
      ome i', 57, 2, 1)
      return
      end subroutine dbintk
      
      subroutine dbnfac(w,nroww,nrow,nbandl,nbandu,iflag)
!***begin prologue  dbnfac
!***refer to  dbint4,dbintk
!
!  dbnfac is the banfac routine from
!        * a practical guide to splines *  by c. de boor
!
!  dbnfac is a double precision routine
!
!  returns in  w  the lu-factorization (without pivoting) of the banded
!  matrix  a  of order  nrow  with  (nbandl + 1 + nbandu) bands or diag-
!  onals in the work array  w .
!
! *****  i n p u t  ****** w is double precision
!  w.....work array of size  (nroww,nrow)  containing the interesting
!        part of a banded matrix  a , with the diagonals or bands of  a
!        stored in the rows of  w , while columns of  a  correspond to
!        columns of  w . this is the storage mode used in  linpack  and
!        results in efficient innermost loops.
!           explicitly,  a  has  nbandl  bands below the diagonal
!                            +     1     (main) diagonal
!                            +   nbandu  bands above the diagonal
!        and thus, with    middle = nbandu + 1,
!          a(i+j,j)  is in  w(i+middle,j)  for i=-nbandu,...,nbandl
!                                              j=1,...,nrow .
!        for example, the interesting entries of a (1,2)-banded matrix
!        of order  9  would appear in the first  1+1+2 = 4  rows of  w
!        as follows.
!                          13 24 35 46 57 68 79
!                       12 23 34 45 56 67 78 89
!                    11 22 33 44 55 66 77 88 99
!                    21 32 43 54 65 76 87 98
!
!        all other entries of  w  not identified in this way with an en-
!        try of  a  are never referenced .
!  nroww.....row dimension of the work array  w .
!        must be  >=  nbandl + 1 + nbandu  .
!  nbandl.....number of bands of  a  below the main diagonal
!  nbandu.....number of bands of  a  above the main diagonal .
!
! *****  o u t p u t  ****** w is double precision
!  iflag.....integer indicating success( = 1) or failure ( = 2) .
!     if  iflag = 1, then
!  w.....contains the lu-factorization of  a  into a unit lower triangu-
!        lar matrix  l  and an upper triangular matrix  u (both banded)
!        and stored in customary fashion over the corresponding entries
!        of  a . this makes it possible to solve any particular linear
!        system  a*x = b  for  x  by a
!              call dbnslv ( w, nroww, nrow, nbandl, nbandu, b )
!        with the solution x  contained in  b  on return .
!     if  iflag = 2, then
!        one of  nrow-1, nbandl,nbandu failed to be nonnegative, or else
!        one of the potential pivots was found to be zero indicating
!        that  a  does not have an lu-factorization. this implies that
!        a  is singular in case it is totally positive .
!
! *****  m e t h o d  ******
!     gauss elimination  w i t h o u t  pivoting is used. the routine is
!  intended for use with matrices  a  which do not require row inter-
!  changes during factorization, especially for the  t o t a l l y
!  p o s i t i v e  matrices which occur in spline calculations.
!     the routine should not be used for an arbitrary banded matrix.
!***routines called  (none)
!***end prologue  dbnfac
!
      integer iflag, nbandl, nbandu, nrow, nroww, i, ipk, j, jmax, k,&
       kmax, middle, midmk, nrowm1
      double precision w(nroww,nrow), factor, pivot
!
!***first executable statement  dbnfac
      iflag = 1
      middle = nbandu + 1
!                         w(middle,.) contains the main diagonal of  a .
      nrowm1 = nrow - 1
      if (nrowm1) 120, 110, 10
   10 if (nbandl>0) go to 30
!                a is upper triangular. check that diagonal is nonzero .
      do 20 i=1,nrowm1
        if (w(middle,i)==0.0d0) go to 120
   20 continue
      go to 110
   30 if (nbandu>0) go to 60
!              a is lower triangular. check that diagonal is nonzero and
!                 divide each column by its diagonal .
      do 50 i=1,nrowm1
        pivot = w(middle,i)
        if (pivot==0.0d0) go to 120
        jmax = min0(nbandl,nrow-i)
        do 40 j=1,jmax
          w(middle+j,i) = w(middle+j,i)/pivot
   40   continue
   50 continue
      return
!
!        a  is not just a triangular matrix. construct lu factorization
   60 do 100 i=1,nrowm1
!                                  w(middle,i)  is pivot for i-th step .
        pivot = w(middle,i)
        if (pivot==0.0d0) go to 120
!                 jmax  is the number of (nonzero) entries in column  i
!                     below the diagonal .
        jmax = min0(nbandl,nrow-i)
!              divide each entry in column  i  below diagonal by pivot .
        do 70 j=1,jmax
          w(middle+j,i) = w(middle+j,i)/pivot
   70   continue
!                 kmax  is the number of (nonzero) entries in row  i  to
!                     the right of the diagonal .
        kmax = min0(nbandu,nrow-i)
!                  subtract  a(i,i+k)*(i-th column) from (i+k)-th column
!                  (below row  i ) .
        do 90 k=1,kmax
          ipk = i + k
          midmk = middle - k
          factor = w(midmk,ipk)
          do 80 j=1,jmax
            w(midmk+j,ipk) = w(midmk+j,ipk) - w(middle+j,i)*factor
   80     continue
   90   continue
  100 continue
!                                       check the last diagonal entry .
  110 if (w(middle,nrow)/=0.0d0) return
  120 iflag = 2
      return
      end subroutine dbnfac
      
      subroutine dbnslv(w,nroww,nrow,nbandl,nbandu,b)
!***begin prologue  dbnslv
!***refer to  dbint4,dbintk
!
!  dbnslv is the banslv routine from
!        * a practical guide to splines *  by c. de boor
!
!  dbnslv is a double precision routine
!
!  companion routine to  dbnfac . it returns the solution  x  of the
!  linear system  a*x = b  in place of  b , given the lu-factorization
!  for  a  in the work array  w from dbnfac.
!
! *****  i n p u t  ****** w,b are double precision
!  w, nroww,nrow,nbandl,nbandu.....describe the lu-factorization of a
!        banded matrix  a  of order  nrow  as constructed in  dbnfac .
!        for details, see  dbnfac .
!  b.....right side of the system to be solved .
!
! *****  o u t p u t  ****** b is double precision
!  b.....contains the solution  x , of order  nrow .
!
! *****  m e t h o d  ******
!     (with  a = l*u, as stored in  w,) the unit lower triangular system
!  l(u*x) = b  is solved for  y = u*x, and  y  stored in  b . then the
!  upper triangular system  u*x = y  is solved for  x  . the calcul-
!  ations are so arranged that the innermost loops stay within columns.
!***routines called  (none)
!***end prologue  dbnslv
!
      integer nbandl, nbandu, nrow, nroww, i, j, jmax, middle, nrowm1
      double precision w(nroww,nrow), b(nrow)
!***first executable statement  dbnslv
      middle = nbandu + 1
      if (nrow==1) go to 80
      nrowm1 = nrow - 1
      if (nbandl==0) go to 30
!                                 forward pass
!            for i=1,2,...,nrow-1, subtract  right side(i)*(i-th column
!            of  l )  from right side  (below i-th row) .
      do 20 i=1,nrowm1
        jmax = min0(nbandl,nrow-i)
        do 10 j=1,jmax
          b(i+j) = b(i+j) - b(i)*w(middle+j,i)
   10   continue
   20 continue
!                                 backward pass
!            for i=nrow,nrow-1,...,1, divide right side(i) by i-th diag-
!            onal entry of  u, then subtract  right side(i)*(i-th column
!            of  u)  from right side  (above i-th row).
   30 if (nbandu>0) go to 50
!                                a  is lower triangular .
      do 40 i=1,nrow
        b(i) = b(i)/w(1,i)
   40 continue
      return
   50 i = nrow
   60 b(i) = b(i)/w(middle,i)
      jmax = min0(nbandu,i-1)
      do 70 j=1,jmax
        b(i-j) = b(i-j) - b(i)*w(middle-j,i)
   70 continue
      i = i - 1
      if (i>1) go to 60
   80 b(1) = b(1)/w(middle,1)
      return
      end subroutine dbnslv
      
      subroutine dbsgq8(fun,xt,bc,n,kk,id,a,b,inbv,err,ans,ierr,work)
!***begin prologue  dbsgq8
!***revision history  (yymmdd)
!   000330  modified array declarations.  (jec)
!
!***refer to  dbfqad
!
!     written by r.e. jones and modified by d.e. amos
!
!     abstract    **** a double precision routine ****
!
!        dbsgq8, a modification of gaus8, integrates the
!        product of fun(x) by the id-th derivative of a spline
!        dbvalu(xt,bc,n,kk,id,x,inbv,work)  between limits a and b.
!
!        dbsgq8 calls dbvalu, dintrv, i1mach, d1mach, xerror
!
!     description of arguments
!
!        input-- fun,xt,bc,a,b,err are double precision
!        fun - name of external function of one argument which
!              multiplies dbvalu.
!        xt  - knot array for dbvalu
!        bc  - b-coefficient array for dbvalu
!        n   - number of b-coefficients for dbvalu
!        kk  - order of the spline, kk>=1
!        id  - order of the spline derivative, 0<=id<=kk-1
!        a   - lower limit of integral
!        b   - upper limit of integral (may be less than a)
!        inbv- initialization parameter for dbvalu
!        err - is a requested pseudorelative error tolerance.  normally
!              pick a value of dabs(err)<1d-3.  ans will normally
!              have no more error than dabs(err) times the integral of
!              the absolute value of fun(x)*dbvalu(xt,bc,n,kk,x,id,
!              inbv,work).
!
!
!        output-- err,ans,work are double precision
!        err - will be an estimate of the absolute error in ans if the
!              input value of err was negative.  (err is unchanged if
!              the input value of err was nonnegative.)  the estimated
!              error is solely for information to the user and should
!              not be used as a correction to the computed integral.
!        ans - computed value of integral
!        ierr- a status code
!            --normal codes
!               1 ans most likely meets requested error tolerance,
!                 or a=b.
!              -1 a and b are too nearly equal to allow normal
!                 integration.  ans is set to zero.
!            --abnormal code
!               2 ans probably does not meet requested error tolerance.
!        work- work vector of length 3*k for dbvalu
!***routines called  d1mach,dbvalu,i1mach,xerror
!***end prologue  dbsgq8
!
      integer icall,id,ierr,inbv,k, kk, kml, kmx, l, lmn, lmx, lr, mxl,&
       n, nbits, nib, nlmn, nlmx
      double precision a,aa,ae,anib,ans,area,b,bc,c,ce,ee,ef,eps,err,&
       est,gl,glr,gr,hh,sq2,tol,vl,vr,work,w1, w2, w3, w4, xt, x1,&
       x2, x3, x4, x, h
      double precision,external :: fun
      dimension xt(*), bc(*), work(*)
      dimension aa(60), hh(60), lr(60), vl(60), gr(60)
      data x1, x2, x3, x4/&
           1.83434642495649805d-01,     5.25532409916328986d-01,&
           7.96666477413626740d-01,     9.60289856497536232d-01/
      data w1, w2, w3, w4/&
           3.62683783378361983d-01,     3.13706645877887287d-01,&
           2.22381034453374471d-01,     1.01228536290376259d-01/
      data icall  /  0  /
      data sq2/1.41421356d0/
      data nlmn/1/,kmx/5000/,kml/6/
!
!     initialize
!
!***first executable statement  dbsgq8
      if (icall/=0)  call xerror(  'dbsgq8- dbsgq8 called recursively.&
        recursive calls are illegal in fortran.', 75, 7, 2)
      icall = 1
      k = DIGITS(1.0d0)  !i1mach(14)
      anib = d1mach(5)*dble(float(k))/0.30102000d0
      nbits = int(sngl(anib))
      nlmx = min0((nbits*5)/8,60)
      ans = 0.0d0
      ierr = 1
      ce = 0.0d0
      if (a==b) go to 140
      lmx = nlmx
      lmn = nlmn
      if (b==0.0d0) go to 10
      if (dsign(1.0d0,b)*a<=0.0d0) go to 10
      c = dabs(1.0d0-a/b)
      if (c>0.1d0) go to 10
      if (c<=0.0d0) go to 140
      anib = 0.5d0 - dlog(c)/0.69314718d0
      nib = int(sngl(anib))
      lmx = min0(nlmx,nbits-nib-7)
      if (lmx<1) go to 130
      lmn = min0(lmn,lmx)
   10 tol = dmax1(dabs(err),2.0d0**(5-nbits))/2.0d0
      if (err==0.0d0) tol = dsqrt(d1mach(4))
      eps = tol
      hh(1) = (b-a)/4.0d0
      aa(1) = a
      lr(1) = 1
      l = 1
      est = g8(aa(l)+2.0d0*hh(l),2.0d0*hh(l),fun)
      k = 8
      area = dabs(est)
      ef = 0.5d0
      mxl = 0
!
!     compute refined estimates, estimate the error, etc.
!
   20 gl = g8(aa(l)+hh(l),hh(l),fun)
      gr(l) = g8(aa(l)+3.0d0*hh(l),hh(l),fun)
      k = k + 16
      area = area + (dabs(gl)+dabs(gr(l))-dabs(est))
      glr = gl + gr(l)
      ee = dabs(est-glr)*ef
      ae = dmax1(eps*area,tol*dabs(glr))
      if (ee-ae) 40, 40, 50
   30 mxl = 1
   40 ce = ce + (est-glr)
      if (lr(l)) 60, 60, 80
!
!     consider the left half of this level
!
   50 if (k>kmx) lmx = kml
      if (l>=lmx) go to 30
      l = l + 1
      eps = eps*0.5d0
      ef = ef/sq2
      hh(l) = hh(l-1)*0.5d0
      lr(l) = -1
      aa(l) = aa(l-1)
      est = gl
      go to 20
!
!     proceed to right half at this level
!
   60 vl(l) = glr
   70 est = gr(l-1)
      lr(l) = 1
      aa(l) = aa(l) + 4.0d0*hh(l)
      go to 20
!
!     return one level
!
   80 vr = glr
   90 if (l<=1) go to 120
      l = l - 1
      eps = eps*2.0d0
      ef = ef*sq2
      if (lr(l)) 100, 100, 110
  100 vl(l) = vl(l+1) + vr
      go to 70
  110 vr = vl(l+1) + vr
      go to 90
!
!      exit
!
  120 ans = vr
      if ((mxl==0) .or. (dabs(ce)<=2.0d0*tol*area)) go to 140
      ierr = 2
      call xerror( 'dbsgq8- ans is probably insufficiently accurate.',&
       48, 3, 1)
      go to 140
  130 ierr = -1
      call xerror( 'dbsgq8- the following temporary diagnostic will appe&
      ar only once.  a and b are too nearly equal to allow normal integr&
      ation.  ans is set to zero, and ierr=-1.', 158, 1, -1)
  140 icall = 0
      if (err<0.0d0) err = ce
      return
	  
	  contains
	  
	  function g8(x,h,fun)  !replaced statement function in original code
	  
	  implicit none
	  
	  double precision,intent(in) :: x
	  double precision,intent(in) :: h
	  double precision :: g8
	  double precision,external :: fun
	  
		  g8 = h*((w1*(fun(x-x1*h)*dbvalu(xt,bc,n,kk,id,x-x1*h,inbv,work)&
					  + fun(x+x1*h)*dbvalu(xt,bc,n,kk,id,x+x1*h,inbv,work))&
					+w2*(fun(x-x2*h)*dbvalu(xt,bc,n,kk,id,x-x2*h,inbv,work)+&
						fun(x+x2*h)*dbvalu(xt,bc,n,kk,id,x+x2*h,inbv,work)))&
				   +(w3*(fun(x-x3*h)*dbvalu(xt,bc,n,kk,id,x-x3*h,inbv,work)+&
						 fun(x+x3*h)*dbvalu(xt,bc,n,kk,id,x+x3*h,inbv,work))&
					+w4*(fun(x-x4*h)*dbvalu(xt,bc,n,kk,id,x-x4*h,inbv,work)+&
					   fun(x+x4*h)*dbvalu(xt,bc,n,kk,id,x+x4*h,inbv,work))))

	  end function g8
	  
      end subroutine dbsgq8
      
      subroutine dbspdr(t,a,n,k,nderiv,ad)
!***begin prologue  dbspdr
!***date written   800901   (yymmdd)
!***revision date  820801   (yymmdd)
!***revision history  (yymmdd)
!   000330  modified array declarations.  (jec)
!
!***category no.  e3,k6
!***keywords  b-spline,data fitting,double precision,interpolation,
!             spline
!***author  amos, d. e., (snla)
!***purpose  uses the b-representation to construct a divided difference
!            table ad preparatory to a (right) derivative calculation
!            in bspev.
!***description
!
!     written by carl de boor and modified by d. e. amos
!
!     reference
!         siam j. numerical analysis, 14, no. 3, june, 1977, pp.441-472.
!
!     abstract     **** a double precision routine ****
!         dbspdr is the bspldr routine of the reference.
!
!         dbspdr uses the b-representation (t,a,n,k) to construct a
!         divided difference table adif preparatory to a (right)
!         derivative calculation in dbspev.  the lower triangular matrix
!         adif is stored in vector ad by columns.  the arrays are
!         related by
!
!           adif(i,j) = ad(i-j+1 + (2*n-j+2)*(j-1)/2)
!
!         i = j,n   ,   j=1,nderiv.
!
!     description of arguments
!
!         input      t,a are double precision
!          t       - knot vector of length n+k
!          a       - b-spline coefficient vector of length n
!          n       - number of b-spline coefficients
!                    n = sum of knot multiplicities-k
!          k       - order of the spline, k >= 1
!          nderiv  - number of derivatives, 1 <= nderiv <= k.
!                    nderiv=1 gives the zero-th derivative =
!                    function value
!
!         output     ad is double precision
!          ad      - table of differences in a vector of length
!                    (2*n-nderiv+1)*nderiv/2 for input to dbspev
!
!     error conditions
!         improper input is a fatal error
!***references  c. de boor, *package for calculating with b-splines*,
!                 siam journal on numerical analysis, volume 14, no. 3,
!                 june 1977, pp. 441-472.
!***routines called  xerror
!***end prologue  dbspdr
!
!
      integer i, id, ii, ipkmid, jj, jm, k, kmid, n, nderiv
      double precision a, ad, diff, fkmid, t
!     dimension t(n+k), ad((2*n-nderiv+1)*nderiv/2)
      dimension t(*), a(n), ad(*)
!***first executable statement  dbspdr
      if(k<1) go to 100
      if(n<k) go to 105
      if(nderiv<1 .or. nderiv>k) go to 110
      do 10 i=1,n
        ad(i) = a(i)
   10 continue
      if (nderiv==1) return
      kmid = k
      jj = n
      jm = 0
      do 30 id=2,nderiv
        kmid = kmid - 1
        fkmid = dble(float(kmid))
        ii = 1
        do 20 i=id,n
          ipkmid = i + kmid
          diff = t(ipkmid) - t(i)
          if (diff/=0.0d0) ad(ii+jj) = (ad(ii+jm+1)-ad(ii+jm))/&
           diff*fkmid
          ii = ii + 1
   20   continue
        jm = jj
        jj = jj + n - id + 1
   30 continue
      return
!
!
  100 continue
      call xerror( ' dbspdr,  k does not satisfy k>=1', 35, 2, 1)
      return
  105 continue
      call xerror( ' dbspdr,  n does not satisfy n>=k', 35, 2, 1)
      return
  110 continue
      call xerror( ' dbspdr,  nderiv does not satisfy 1<=nderiv<=k',&
       50, 2, 1)
      return
      end subroutine dbspdr
      
      subroutine dbspev(t,ad,n,k,nderiv,x,inev,svalue,work)
!***begin prologue  dbspev
!***date written   800901   (yymmdd)
!***revision date  840425   (yymmdd)
!***revision history  (yymmdd)
!   000330  modified array declarations.  (jec)
!
!***category no.  e3,k6
!***keywords  b-spline,data fitting,double precision,interpolation,
!             spline
!***author  amos, d. e., (snla)
!***purpose  calculates the value of the spline and its derivatives at x
!            from the b-representation .
!***description
!
!     written by carl de boor and modified by d. e. amos
!
!     reference
!         siam j. numerical analysis, 14, no. 3, june, 1977, pp.441-472.
!
!     abstract    **** a double precision routine ****
!         dbspev is the bsplev routine of the reference.
!
!         dbspev calculates the value of the spline and its derivatives
!         at x from the b-representation (t,a,n,k) and returns them in
!         svalue(i),i=1,nderiv, t(k) <= x <= t(n+1).  ad(i) can be
!         the b-spline coefficients a(i), i=1,n) if nderiv=1.  otherwise
!         ad must be computed before hand by a call to dbspdr (t,a,n,k,
!         nderiv,ad).  if x=t(i),i=k,n), right limiting values are
!         obtained.
!
!         to compute left derivatives or left limiting values at a
!         knot t(i), replace n by i-1 and set x=t(i), i=k+1,n+1.
!
!         dbspev calls dintrv, dbspvn
!
!     description of arguments
!
!         input      t,ad,x, are double precision
!          t       - knot vector of length n+k
!          ad      - vector of length (2*n-nderiv+1)*nderiv/2 containing
!                    the difference table from dbspdr.
!          n       - number of b-spline coefficients
!                    n = sum of knot multiplicities-k
!          k       - order of the b-spline, k >= 1
!          nderiv  - number of derivatives, 1 <= nderiv <= k.
!                    nderiv=1 gives the zero-th derivative =
!                    function value
!          x       - argument, t(k) <= x <= t(n+1)
!          inev    - an initialization parameter which must be set
!                    to 1 the first time dbspev is called.
!
!         output     svalue,work are double precision
!          inev    - inev contains information for efficient process-
!                    ing after the initial call and inev must not
!                    be changed by the user.  distinct splines require
!                    distinct inev parameters.
!          svalue  - vector of length nderiv containing the spline
!                    value in svalue(1) and the nderiv-1 derivatives
!                    in the remaining components.
!          work    - work vector of length 3*k
!
!     error conditions
!         improper input is a fatal error.
!***references  c. de boor, *package for calculating with b-splines*,
!                 siam journal on numerical analysis, volume 14, no. 3,
!                 june 1977, pp. 441-472.
!***routines called  dbspvn,dintrv,xerror
!***end prologue  dbspev
!
!
      integer i,id,inev,iwork,jj,k,kp1,kp1mn,l,left,ll,mflag,&
       n, nderiv
      double precision ad, svalue, sum, t, work, x
!     dimension t(n+k)
      dimension t(*), ad(*), svalue(nderiv), work(*)
!***first executable statement  dbspev
      if(k<1) go to 100
      if(n<k) go to 105
      if(nderiv<1 .or. nderiv>k) go to 115
      id = nderiv
      call dintrv(t, n+1, x, inev, i, mflag)
      if (x<t(k)) go to 110
      if (mflag==0) go to 30
      if (x>t(i)) go to 110
   20 if (i==k) go to 120
      i = i - 1
      if (x==t(i)) go to 20
!
! *i* has been found in (k,n) so that t(i) <= x < t(i+1)
!     (or <= t(i+1), if t(i) < t(i+1) = t(n+1) ).
   30 kp1mn = k + 1 - id
      kp1 = k + 1
      call dbspvn(t, kp1mn, k, 1, x, i, work(1),work(kp1),iwork)
      jj = (n+n-id+2)*(id-1)/2
!     adif(leftpl,id) = ad(leftpl-id+1 + (2*n-id+2)*(id-1)/2)
!     leftpl = left + l
   40 left = i - kp1mn
      sum = 0.0d0
      ll = left + jj + 2 - id
      do 50 l=1,kp1mn
        sum = sum + work(l)*ad(ll)
        ll = ll + 1
   50 continue
      svalue(id) = sum
      id = id - 1
      if (id==0) go to 60
      jj = jj-(n-id+1)
      kp1mn = kp1mn + 1
      call dbspvn(t, kp1mn, k, 2, x, i, work(1), work(kp1),iwork)
      go to 40
!
   60 return
!
!
  100 continue
      call xerror( ' dbspev,  k does not satisfy k>=1',35,2,1)
      return
  105 continue
      call xerror( ' dbspev,  n does not satisfy n>=k',35,2,1)
      return
  110 continue
      call xerror( ' dbspev,  x is not in t(k)<=x<=t(n+1)',41,2,1)
      return
  115 continue
      call xerror( ' dbspev,  nderiv does not satisfy 1<=nderiv<=k',&
       50, 2, 1)
      return
  120 continue
      call xerror( ' dbspev, a left limiting value cannot be obtained at&
       t(k)',57,2,1)
      return
      end subroutine dbspev
      
      subroutine dbsppp(t,a,n,k,ldc,c,xi,lxi,work)
!***begin prologue  dbsppp
!***date written   800901   (yymmdd)
!***revision date  820801   (yymmdd)
!***revision history  (yymmdd)
!   000330  modified array declarations.  (jec)
!
!***category no.  e3,k6
!***keywords  b-spline,data fitting,double precision,interpolation,
!             spline
!***author  amos, d. e., (snla)
!***purpose  converts the b-representation to the piecewise
!            polynomial (pp) form for use with ppval.
!***description
!
!     written by carl de boor and modified by d. e. amos
!
!     reference
!         siam j. numerical analysis, 14, no. 3, june, 1977, pp.441-472.
!
!     abstract    **** a double precision routine ****
!         dbsppp is the bsplpp routine of the reference.
!
!         dbsppp converts the b-representation (t,a,n,k) to the
!         piecewise polynomial (pp) form (c,xi,lxi,k) for use with
!         dppval.  here xi(*), the break point array of length lxi, is
!         the knot array t(*) with multiplicities removed.  the columns
!         of the matrix c(i,j) contain the right taylor derivatives
!         for the polynomial expansion about xi(j) for the intervals
!         xi(j) <= x <= xi(j+1), i=1,k, j=1,lxi.  function dppval
!         makes this evaluation at a specified point x in
!         xi(1) <= x <= xi(lxi+1)
!
!         dbsppp calls dbspdr, dbspev, dintrv, dbspvn
!
!     description of arguments
!
!         input      t,a are double precision
!          t       - knot vector of length n+k
!          a       - b-spline coefficient vector of length n
!          n       - number of b-spline coefficients
!                    n = sum of knot multiplicities-k
!          k       - order of the b-spline, k >= 1
!          ldc     - leading dimension of c, ldc >= k
!
!         output     c,xi,work are double precision
!          c       - matrix of dimension at least (k,lxi) containing
!                    right derivatives at break points
!          xi      - xi break point vector of length lxi+1
!          lxi     - number of break points, lxi <= n-k+1
!          work    - work vector of length k*(n+3)
!
!     error conditions
!         improper input is a fatal error
!***references  c. de boor, *package for calculating with b-splines*,
!                 siam journal on numerical analysis, volume 14, no. 3,
!                 june 1977, pp. 441-472.
!***routines called  dbspdr,dbspev,xerror
!***end prologue  dbsppp
!
!
      integer ileft, inev, k, ldc, lxi, n, nk
      double precision a, c, t, work, xi
!     dimension t(n+k),xi(lxi+1),c(ldc,*)
!     here, * = the final value of the output parameter lxi.
      dimension t(*), a(n), work(*), xi(*), c(ldc,*)
!***first executable statement  dbsppp
      if(k<1) go to 100
      if(n<k) go to 105
      if(ldc<k) go to 110
      call dbspdr(t, a, n, k, k, work)
      lxi = 0
      xi(1) = t(k)
      inev = 1
      nk = n*k + 1
      do 10 ileft=k,n
        if (t(ileft+1)==t(ileft)) go to 10
        lxi = lxi + 1
        xi(lxi+1) = t(ileft+1)
        call dbspev(t,work(1),n,k, k,xi(lxi),inev,c(1,lxi),work(nk))
   10 continue
      return
  100 continue
      call xerror( ' dbsppp,  k does not satisfy k>=1', 35, 2, 1)
      return
  105 continue
      call xerror( ' dbsppp,  n does not satisfy n>=k', 35, 2, 1)
      return
  110 continue
      call xerror( ' dbsppp,  ldc does not satisfy ldc>=k', 39, 2, 1)
      return
      end subroutine dbsppp
      
      subroutine dbspvd(t,k,nderiv,x,ileft,ldvnik,vnikx,work)
!***begin prologue  dbspvd
!***date written   800901   (yymmdd)
!***revision date  820801   (yymmdd)
!***revision history  (yymmdd)
!   000330  modified array declarations.  (jec)
!
!***category no.  e3,k6
!***keywords  b-spline,data fitting,double precision,interpolation,
!             spline
!***author  amos, d. e., (snla)
!***purpose  calculates the value and all derivatives of order less than
!            nderiv of all basis functions which do not vanish at x.
!***description
!
!     written by carl de boor and modified by d. e. amos
!
!     reference
!         siam j. numerical analysis, 14, no. 3, june, 1977, pp.441-472.
!
!     abstract    **** a double precision routine ****
!
!         dbspvd is the bsplvd routine of the reference.
!
!         dbspvd calculates the value and all derivatives of order
!         less than nderiv of all basis functions which do not
!         (possibly) vanish at x.  ileft is input such that
!         t(ileft) <= x < t(ileft+1).  a call to intrv(t,n+1,x,
!         ilo,ileft,mflag) will produce the proper ileft.  the output of
!         dbspvd is a matrix vnikx(i,j) of dimension at least (k,nderiv)
!         whose columns contain the k nonzero basis functions and
!         their nderiv-1 right derivatives at x, i=1,k, j=1,nderiv.
!         these basis functions have indices ileft-k+i, i=1,k,
!         k <= ileft <= n.  the nonzero part of the i-th basis
!         function lies in (t(i),t(i+k)), i=1,n).
!
!         if x=t(ileft+1) then vnikx contains left limiting values
!         (left derivatives) at t(ileft+1).  in particular, ileft = n
!         produces left limiting values at the right end point
!         x=t(n+1).  to obtain left limiting values at t(i), i=k+1,n+1,
!         set x= next lower distinct knot, call intrv to get ileft,
!         set x=t(i), and then call dbspvd.
!
!         dbspvd calls dbspvn
!
!     description of arguments
!         input      t,x are double precision
!          t       - knot vector of length n+k, where
!                    n = number of b-spline basis functions
!                    n = sum of knot multiplicities-k
!          k       - order of the b-spline, k >= 1
!          nderiv  - number of derivatives = nderiv-1,
!                    1 <= nderiv <= k
!          x       - argument of basis functions,
!                    t(k) <= x <= t(n+1)
!          ileft   - largest integer such that
!                    t(ileft) <= x <  t(ileft+1)
!          ldvnik  - leading dimension of matrix vnikx
!
!         output     vnikx,work are double precision
!          vnikx   - matrix of dimension at least (k,nderiv) contain-
!                    ing the nonzero basis functions at x and their
!                    derivatives columnwise.
!          work    - a work vector of length (k+1)*(k+2)/2
!
!     error conditions
!         improper input is a fatal error
!***references  c. de boor, *package for calculating with b-splines*,
!                 siam journal on numerical analysis, volume 14, no. 3,
!                 june 1977, pp. 441-472.
!***routines called  dbspvn,xerror
!***end prologue  dbspvd
!
!
      integer i,ideriv,ileft,ipkmd,j,jj,jlow,jm,jp1mid,k,kmd, kp1, l,&
       ldummy, m, mhigh, nderiv, ldvnik, iwork
      double precision factor, fkmd, t, v, vnikx, work, x
!     dimension t(ileft+k), work((k+1)*(k+2)/2)
!     a(i,j) = work(i+j*(j+1)/2),  i=1,j+1  j=1,k-1
!     a(i,k) = w0rk(i+k*(k-1)/2)  i=1.k
!     work(1) and work((k+1)*(k+2)/2) are not used.
      dimension t(*), vnikx(ldvnik,nderiv), work(*)
!***first executable statement  dbspvd
      if(k<1) go to 200
      if(nderiv<1 .or. nderiv>k) go to 205
      if(ldvnik<k) go to 210
      ideriv = nderiv
      kp1 = k + 1
      jj = kp1 - ideriv
      call dbspvn(t, jj, k, 1, x, ileft, vnikx, work, iwork)
      if (ideriv==1) go to 100
      mhigh = ideriv
      do 20 m=2,mhigh
        jp1mid = 1
        do 10 j=ideriv,k
          vnikx(j,ideriv) = vnikx(jp1mid,1)
          jp1mid = jp1mid + 1
   10   continue
        ideriv = ideriv - 1
        jj = kp1 - ideriv
        call dbspvn(t, jj, k, 2, x, ileft, vnikx, work, iwork)
   20 continue
!
      jm = kp1*(kp1+1)/2
      do 30 l = 1,jm
        work(l) = 0.0d0
   30 continue
!     a(i,i) = work(i*(i+3)/2) = 1.0       i = 1,k
      l = 2
      j = 0
      do 40 i = 1,k
        j = j + l
        work(j) = 1.0d0
        l = l + 1
   40 continue
      kmd = k
      do 90 m=2,mhigh
        kmd = kmd - 1
        fkmd = float(kmd)
        i = ileft
        j = k
        jj = j*(j+1)/2
        jm = jj - j
        do 60 ldummy=1,kmd
          ipkmd = i + kmd
          factor = fkmd/(t(ipkmd)-t(i))
          do 50 l=1,j
            work(l+jj) = (work(l+jj)-work(l+jm))*factor
   50     continue
          i = i - 1
          j = j - 1
          jj = jm
          jm = jm - j
   60   continue
!
        do 80 i=1,k
          v = 0.0d0
          jlow = max0(i,m)
          jj = jlow*(jlow+1)/2
          do 70 j=jlow,k
            v = work(i+jj)*vnikx(j,m) + v
            jj = jj + j + 1
   70     continue
          vnikx(i,m) = v
   80   continue
   90 continue
  100 return
!
!
  200 continue
      call xerror( ' dbspvd,  k does not satisfy k>=1', 35, 2, 1)
      return
  205 continue
      call xerror( ' dbspvd,  nderiv does not satisfy 1<=nderiv<=k',&
       50, 2, 1)
      return
  210 continue
      call xerror( ' dbspvd,  ldvnik does not satisfy ldvnik>=k',45,&
       2, 1)
      return
      end subroutine dbspvd
      
      subroutine dbspvn(t,jhigh,k,index,x,ileft,vnikx,work,iwork)
!***begin prologue  dbspvn
!***date written   800901   (yymmdd)
!***revision date  820801   (yymmdd)
!***revision history  (yymmdd)
!   000330  modified array declarations.  (jec)
!
!***category no.  e3,k6
!***keywords  b-spline,data fitting,double precision,interpolation,
!             spline
!***author  amos, d. e., (snla)
!***purpose  calculates the value of all (possibly) nonzero basis
!            functions at x.
!***description
!
!     written by carl de boor and modified by d. e. amos
!
!     reference
!         siam j. numerical analysis, 14, no. 3, june, 1977, pp.441-472.
!
!     abstract    **** a double precision routine ****
!         dbspvn is the bsplvn routine of the reference.
!
!         dbspvn calculates the value of all (possibly) nonzero basis
!         functions at x of order max(jhigh,(j+1)*(index-1)), where t(k)
!         <= x <= t(n+1) and j=iwork is set inside the routine on
!         the first call when index=1.  ileft is such that t(ileft) <=
!         x < t(ileft+1).  a call to dintrv(t,n+1,x,ilo,ileft,mflag)
!         produces the proper ileft.  dbspvn calculates using the basic
!         algorithm needed in dbspvd.  if only basis functions are
!         desired, setting jhigh=k and index=1 can be faster than
!         calling dbspvd, but extra coding is required for derivatives
!         (index=2) and dbspvd is set up for this purpose.
!
!         left limiting values are set up as described in dbspvd.
!
!     description of arguments
!
!         input      t,x are double precision
!          t       - knot vector of length n+k, where
!                    n = number of b-spline basis functions
!                    n = sum of knot multiplicities-k
!          jhigh   - order of b-spline, 1 <= jhigh <= k
!          k       - highest possible order
!          index   - index = 1 gives basis functions of order jhigh
!                          = 2 denotes previous entry with work, iwork
!                              values saved for subsequent calls to
!                              dbspvn.
!          x       - argument of basis functions,
!                    t(k) <= x <= t(n+1)
!          ileft   - largest integer such that
!                    t(ileft) <= x <  t(ileft+1)
!
!         output     vnikx, work are double precision
!          vnikx   - vector of length k for spline values.
!          work    - a work vector of length 2*k
!          iwork   - a work parameter.  both work and iwork contain
!                    information necessary to continue for index = 2.
!                    when index = 1 exclusively, these are scratch
!                    variables and can be used for other purposes.
!
!     error conditions
!         improper input is a fatal error.
!***references  c. de boor, *package for calculating with b-splines*,
!                 siam journal on numerical analysis, volume 14, no. 3,
!                 june 1977, pp. 441-472.
!***routines called  xerror
!***end prologue  dbspvn
!
!
      integer ileft, imjp1, index, ipj, iwork, jhigh, jp1, jp1ml, k, l
      double precision t, vm, vmprev, vnikx, work, x
!     dimension t(ileft+jhigh)
      dimension t(*), vnikx(k), work(*)
!     content of j, deltam, deltap is expected unchanged between calls.
!     work(i) = deltap(i), work(k+i) = deltam(i), i = 1,k
!***first executable statement  dbspvn
      if(k<1) go to 90
      if(jhigh>k .or. jhigh<1) go to 100
      if(index<1 .or. index>2) go to 105
      if(x<t(ileft) .or. x>t(ileft+1)) go to 110
      go to (10, 20), index
   10 iwork = 1
      vnikx(1) = 1.0d0
      if (iwork>=jhigh) go to 40
!
   20 ipj = ileft + iwork
      work(iwork) = t(ipj) - x
      imjp1 = ileft - iwork + 1
      work(k+iwork) = x - t(imjp1)
      vmprev = 0.0d0
      jp1 = iwork + 1
      do 30 l=1,iwork
        jp1ml = jp1 - l
        vm = vnikx(l)/(work(l)+work(k+jp1ml))
        vnikx(l) = vm*work(l) + vmprev
        vmprev = vm*work(k+jp1ml)
   30 continue
      vnikx(jp1) = vmprev
      iwork = jp1
      if (iwork<jhigh) go to 20
!
   40 return
!
!
   90 continue
      call xerror( ' dbspvn,  k does not satisfy k>=1', 35, 2, 1)
      return
  100 continue
      call xerror( ' dbspvn,  jhigh does not satisfy 1<=jhigh<=k',&
       48, 2, 1)
      return
  105 continue
      call xerror( ' dbspvn,  index is not 1 or 2',29,2,1)
      return
  110 continue
      call xerror( ' dbspvn,  x does not satisfy t(ileft)<=x<=t(ilef&
      t+1)', 56, 2, 1)
      return
      end subroutine dbspvn
      
      subroutine dbsqad(t,bcoef,n,k,x1,x2,bquad,work)
!***begin prologue  dbsqad
!***date written   800901   (yymmdd)
!***revision date  820801   (yymmdd)
!***revision history  (yymmdd)
!   000330  modified array declarations.  (jec)
!
!***category no.  h2a2a1,e3,k6
!***keywords  b-spline,data fitting,double precision,interpolation,
!             quadrature,spline
!***author  amos, d. e., (snla)
!***purpose  computes the integral on (x1,x2) of a k-th order
!            b-spline using the b-representation.
!***description
!
!     written by d. e. amos, june, 1979.
!
!      reference
!         sand-79-1825
!
!     abstract    **** a double precision routine ****
!
!         dbsqad computes the integral on (x1,x2) of a k-th order
!         b-spline using the b-representation (t,bcoef,n,k).  orders
!         k as high as 20 are permitted by applying a 2, 6, or 10
!         point gauss formula on subintervals of (x1,x2) which are
!         formed by included (distinct) knots.
!
!         if orders k greater than 20 are needed, use dbfqad with
!         f(x) = 1.
!
!         the maximum number of significant digits obtainable in
!         dbsqad is the smaller of 18 and the number of digits
!         carried in double precision arithmetic.
!
!         dbsqad calls dintrv, dbvalu, xerror
!
!     description of arguments
!         input      t,bcoef,x1,x2 are double precision
!           t      - knot array of length n+k
!           bcoef  - b-spline coefficient array of length n
!           n      - length of coefficient array
!           k      - order of b-spline, 1 <= k <= 20
!           x1,x2  - end points of quadrature interval in
!                    t(k) <= x <= t(n+1)
!
!         output     bquad,work are double precision
!           bquad  - integral of the b-spline over (x1,x2)
!           work   - work vector of length 3*k
!
!     error conditions
!         improper input is a fatal error
!***references  d.e. amos, *quadrature subroutines for splines and
!                 b-splines*, sand79-1825, sandia laboratories,
!                 december 1979.
!***routines called  dbvalu,dintrv,xerror
!***end prologue  dbsqad
!
!
      integer i,il1,il2,ilo,inbv, jf,k,left,m,mf,mflag,n, npk, np1
      double precision a,aa,b,bb,bcoef,bma,bpa,bquad,c1,gpts,gwts,gx,q,&
       sum, t, ta, tb, work, x1, x2, y1, y2
      dimension t(*), bcoef(*), gpts(9), gwts(9), sum(5), work(*)
!
      data gpts(1), gpts(2), gpts(3), gpts(4), gpts(5), gpts(6),&
           gpts(7), gpts(8), gpts(9)/&
           5.77350269189625764d-01,     2.38619186083196909d-01,&
           6.61209386466264514d-01,     9.32469514203152028d-01,&
           1.48874338981631211d-01,     4.33395394129247191d-01,&
           6.79409568299024406d-01,     8.65063366688984511d-01,&
           9.73906528517171720d-01/
      data gwts(1), gwts(2), gwts(3), gwts(4), gwts(5), gwts(6),&
           gwts(7), gwts(8), gwts(9)/&
           1.00000000000000000d+00,     4.67913934572691047d-01,&
           3.60761573048138608d-01,     1.71324492379170345d-01,&
           2.95524224714752870d-01,     2.69266719309996355d-01,&
           2.19086362515982044d-01,     1.49451349150580593d-01,&
           6.66713443086881376d-02/
!
!***first executable statement  dbsqad
      bquad = 0.0d0
      if(k<1 .or. k>20) go to 65
      if(n<k) go to 70
      aa = dmin1(x1,x2)
      bb = dmax1(x1,x2)
      if (aa<t(k)) go to 60
      np1 = n + 1
      if (bb>t(np1)) go to 60
      if (aa==bb) return
      npk = n + k
!     selection of 2, 6, or 10 point gauss formula
      jf = 0
      mf = 1
      if (k<=4) go to 10
      jf = 1
      mf = 3
      if (k<=12) go to 10
      jf = 4
      mf = 5
   10 continue
!
      do 20 i=1,mf
        sum(i) = 0.0d0
   20 continue
      ilo = 1
      inbv = 1
      call dintrv(t, npk, aa, ilo, il1, mflag)
      call dintrv(t, npk, bb, ilo, il2, mflag)
      if (il2>=np1) il2 = n
      do 40 left=il1,il2
        ta = t(left)
        tb = t(left+1)
        if (ta==tb) go to 40
        a = dmax1(aa,ta)
        b = dmin1(bb,tb)
        bma = 0.5d0*(b-a)
        bpa = 0.5d0*(b+a)
        do 30 m=1,mf
          c1 = bma*gpts(jf+m)
          gx = -c1 + bpa
          y2 = dbvalu(t,bcoef,n,k,0,gx,inbv,work)
          gx = c1 + bpa
          y1 = dbvalu(t,bcoef,n,k,0,gx,inbv,work)
          sum(m) = sum(m) + (y1+y2)*bma
   30   continue
   40 continue
      q = 0.0d0
      do 50 m=1,mf
        q = q + gwts(jf+m)*sum(m)
   50 continue
      if (x1>x2) q = -q
      bquad = q
      return
!
!
   60 continue
      call xerror( ' dbsqad,  x1 or x2 or both do not satisfy t(k)<=x.&
      le.t(n+1)',  61, 2, 1)
      return
   65 continue
      call xerror( ' dbsqad,  k does not satisfy 1<=k<=20',&
       41, 2, 1)
      return
   70 continue
      call xerror( ' dbsqad,  n does not satisfy n>=k', 35, 2, 1)
      return
      end subroutine dbsqad
      
      double precision function dbvalu(t,a,n,k,ideriv,x,inbv,work)
!***begin prologue  dbvalu
!***date written   800901   (yymmdd)
!***revision date  820801   (yymmdd)
!***revision history  (yymmdd)
!   000330  modified array declarations.  (jec)
!
!***category no.  e3,k6
!***keywords  b-spline,data fitting,double precision,interpolation,
!             spline
!***author  amos, d. e., (snla)
!***purpose  evaluates the b-representation of a b-spline at x for the
!            function value or any of its derivatives.
!***description
!
!     written by carl de boor and modified by d. e. amos
!
!     reference
!         siam j. numerical analysis, 14, no. 3, june, 1977, pp.441-472.
!
!     abstract   **** a double precision routine ****
!         dbvalu is the bvalue function of the reference.
!
!         dbvalu evaluates the b-representation (t,a,n,k) of a b-spline
!         at x for the function value on ideriv=0 or any of its
!         derivatives on ideriv=1,2,...,k-1.  right limiting values
!         (right derivatives) are returned except at the right end
!         point x=t(n+1) where left limiting values are computed.  the
!         spline is defined on t(k) <= x <= t(n+1).  dbvalu returns
!         a fatal error message when x is outside of this interval.
!
!         to compute left derivatives or left limiting values at a
!         knot t(i), replace n by i-1 and set x=t(i), i=k+1,n+1.
!
!         dbvalu calls dintrv
!
!     description of arguments
!
!         input      t,a,x are double precision
!          t       - knot vector of length n+k
!          a       - b-spline coefficient vector of length n
!          n       - number of b-spline coefficients
!                    n = sum of knot multiplicities-k
!          k       - order of the b-spline, k >= 1
!          ideriv  - order of the derivative, 0 <= ideriv <= k-1
!                    ideriv = 0 returns the b-spline value
!          x       - argument, t(k) <= x <= t(n+1)
!          inbv    - an initialization parameter which must be set
!                    to 1 the first time dbvalu is called.
!
!         output     work,dbvalu are double precision
!          inbv    - inbv contains information for efficient process-
!                    ing after the initial call and inbv must not
!                    be changed by the user.  distinct splines require
!                    distinct inbv parameters.
!          work    - work vector of length 3*k.
!          dbvalu  - value of the ideriv-th derivative at x
!
!     error conditions
!         an improper input is a fatal error
!***references  c. de boor, *package for calculating with b-splines*,
!                 siam journal on numerical analysis, volume 14, no. 3,
!                 june 1977, pp. 441-472.
!***routines called  dintrv,xerror
!***end prologue  dbvalu
!
!
      integer i,ideriv,iderp1,ihi,ihmkmj,ilo,imk,imkpj, inbv, ipj,&
       ip1, ip1mj, j, jj, j1, j2, k, kmider, kmj, km1, kpk, mflag, n
      double precision a, fkmj, t, work, x
      dimension t(*), a(n), work(*)
!***first executable statement  dbvalu
      dbvalu = 0.0d0
      if(k<1) go to 102
      if(n<k) go to 101
      if(ideriv<0 .or. ideriv>=k) go to 110
      kmider = k - ideriv
!
! *** find *i* in (k,n) such that t(i) <= x < t(i+1)
!     (or, <= t(i+1) if t(i) < t(i+1) = t(n+1)).
      km1 = k - 1
      call dintrv(t, n+1, x, inbv, i, mflag)
      if (x<t(k)) go to 120
      if (mflag==0) go to 20
      if (x>t(i)) go to 130
   10 if (i==k) go to 140
      i = i - 1
      if (x==t(i)) go to 10
!
! *** difference the coefficients *ideriv* times
!     work(i) = aj(i), work(k+i) = dp(i), work(k+k+i) = dm(i), i=1.k
!
   20 imk = i - k
      do 30 j=1,k
        imkpj = imk + j
        work(j) = a(imkpj)
   30 continue
      if (ideriv==0) go to 60
      do 50 j=1,ideriv
        kmj = k - j
        fkmj = dble(float(kmj))
        do 40 jj=1,kmj
          ihi = i + jj
          ihmkmj = ihi - kmj
          work(jj) = (work(jj+1)-work(jj))/(t(ihi)-t(ihmkmj))*fkmj
   40   continue
   50 continue
!
! *** compute value at *x* in (t(i),(t(i+1)) of ideriv-th derivative,
!     given its relevant b-spline coeff. in aj(1),...,aj(k-ideriv).
   60 if (ideriv==km1) go to 100
      ip1 = i + 1
      kpk = k + k
      j1 = k + 1
      j2 = kpk + 1
      do 70 j=1,kmider
        ipj = i + j
        work(j1) = t(ipj) - x
        ip1mj = ip1 - j
        work(j2) = x - t(ip1mj)
        j1 = j1 + 1
        j2 = j2 + 1
   70 continue
      iderp1 = ideriv + 1
      do 90 j=iderp1,km1
        kmj = k - j
        ilo = kmj
        do 80 jj=1,kmj
          work(jj) = (work(jj+1)*work(kpk+ilo)+work(jj)&
                    *work(k+jj))/(work(kpk+ilo)+work(k+jj))
          ilo = ilo - 1
   80   continue
   90 continue
  100 dbvalu = work(1)
      return
!
!
  101 continue
      call xerror( ' dbvalu,  n does not satisfy n>=k',35,2,1)
      return
  102 continue
      call xerror( ' dbvalu,  k does not satisfy k>=1',35,2,1)
      return
  110 continue
      call xerror( ' dbvalu,  ideriv does not satisfy 0<=ideriv<k',&
       50, 2, 1)
      return
  120 continue
      call xerror( ' dbvalu,  x is n0t greater than or equal to t(k)',&
       48, 2, 1)
      return
  130 continue
      call xerror( ' dbvalu,  x is not less than or equal to t(n+1)',&
       47, 2, 1)
      return
  140 continue
      call xerror( ' dbvalu,  a left limiting value cann0t be obtained a&
      t t(k)',    58, 2, 1)
      return
      end function dbvalu
      
      subroutine dintrv(xt,lxt,x,ilo,ileft,mflag)
!***begin prologue  dintrv
!***date written   800901   (yymmdd)
!***revision date  820801   (yymmdd)
!***category no.  e3,k6
!***keywords  b-spline,data fitting,double precision,interpolation,
!             spline
!***author  amos, d. e., (snla)
!***purpose  computes the largest integer ileft in 1<=ileft<=lxt
!            such that xt(ileft)<=x where xt(*) is a subdivision of
!            the x interval.
!***description
!
!     written by carl de boor and modified by d. e. amos
!
!     reference
!         siam j.  numerical analysis, 14, no. 3, june 1977, pp.441-472.
!
!     abstract    **** a double precision routine ****
!         dintrv is the interv routine of the reference.
!
!         dintrv computes the largest integer ileft in 1 <= ileft <=
!         lxt such that xt(ileft) <= x where xt(*) is a subdivision of
!         the x interval.  precisely,
!
!                      x < xt(1)                1         -1
!         if  xt(i) <= x < xt(i+1)  then  ileft=i  , mflag=0
!           xt(lxt) <= x                         lxt        1,
!
!         that is, when multiplicities are present in the break point
!         to the left of x, the largest index is taken for ileft.
!
!     description of arguments
!
!         input      xt,x are double precision
!          xt      - xt is a knot or break point vector of length lxt
!          lxt     - length of the xt vector
!          x       - argument
!          ilo     - an initialization parameter which must be set
!                    to 1 the first time the spline array xt is
!                    processed by dintrv.
!
!         output
!          ilo     - ilo contains information for efficient process-
!                    ing after the initial call and ilo must not be
!                    changed by the user.  distinct splines require
!                    distinct ilo parameters.
!          ileft   - largest integer satisfying xt(ileft) <= x
!          mflag   - signals when x lies out of bounds
!
!     error conditions
!         none
!***references  c. de boor, *package for calculating with b-splines*,
!                 siam journal on numerical analysis, volume 14, no. 3,
!                 june 1977, pp. 441-472.
!***routines called  (none)
!***end prologue  dintrv
!
!
      integer ihi, ileft, ilo, istep, lxt, mflag, middle
      double precision x, xt
      dimension xt(lxt)
!***first executable statement  dintrv
      ihi = ilo + 1
      if (ihi<lxt) go to 10
      if (x>=xt(lxt)) go to 110
      if (lxt<=1) go to 90
      ilo = lxt - 1
      ihi = lxt
!
   10 if (x>=xt(ihi)) go to 40
      if (x>=xt(ilo)) go to 100
!
! *** now x < xt(ihi) . find lower bound
      istep = 1
   20 ihi = ilo
      ilo = ihi - istep
      if (ilo<=1) go to 30
      if (x>=xt(ilo)) go to 70
      istep = istep*2
      go to 20
   30 ilo = 1
      if (x<xt(1)) go to 90
      go to 70
! *** now x >= xt(ilo) . find upper bound
   40 istep = 1
   50 ilo = ihi
      ihi = ilo + istep
      if (ihi>=lxt) go to 60
      if (x<xt(ihi)) go to 70
      istep = istep*2
      go to 50
   60 if (x>=xt(lxt)) go to 110
      ihi = lxt
!
! *** now xt(ilo) <= x < xt(ihi) . narrow the interval
   70 middle = (ilo+ihi)/2
      if (middle==ilo) go to 100
!     note. it is assumed that middle = ilo in case ihi = ilo+1
      if (x<xt(middle)) go to 80
      ilo = middle
      go to 70
   80 ihi = middle
      go to 70
! *** set output and return
   90 mflag = -1
      ileft = 1
      return
  100 mflag = 0
      ileft = ilo
      return
  110 mflag = 1
      ileft = lxt
      return
      end subroutine dintrv
      
      subroutine dpfqad(f,ldc,c,xi,lxi,k,id,x1,x2,tol,quad,ierr)
!***begin prologue  dpfqad
!***date written   800901   (yymmdd)
!***revision date  820801   (yymmdd)
!***revision history  (yymmdd)
!   000330  modified array declarations.  (jec)
!
!***category no.  h2a2a1,e3,k6
!***keywords  b-spline,data fitting,double precision,interpolation,
!             quadrature,spline
!***author  amos, d. e., (snla)
!***purpose  computes the integral on (x1,x2) of a product of a
!            function f and the id-th derivative of a b-spline,
!            (pp-representation).
!***description
!
!     written by d. e. amos, june, 1979.
!
!     reference sand-79-1825
!
!     abstract    **** a double precision routine ****
!
!         dpfqad computes the integral on (x1,x2) of a product of a
!         function f and the id-th derivative of a b-spline, using the
!         pp-representation (c,xi,lxi,k).  (x1,x2) is normally a sub-
!         interval of xi(1) <= x <= xi(lxi+1).  an integration
!         routine, dppgq8 (a modification of gaus8), integrates the
!         product on subintervals of (x1,x2) formed by the included
!         break points.  integration outside of (xi(1),xi(lxi+1)) is
!         permitted provided f is defined.
!
!         the maximum number of significant digits obtainable in
!         dbsqad is the smaller of 18 and the number of digits
!         carried in double precision arithmetic.
!
!         dpfqad calls dintrv, dppval, dppgq8, d1mach, xerror
!
!     description of arguments
!         input      f,c,xi,x1,x2,tol are double precision
!           f      - external function of one argument for the
!                    integrand pf(x)=f(x)*dppval(ldc,c,xi,lxi,k,id,x,
!                    inppv)
!           ldc    - leading dimension of matrix c, ldc >= k
!           c(i,j) - right taylor derivatives at xi(j), i=1,k , j=1,lxi
!           xi(*)  - break point array of length lxi+1
!           lxi    - number of polynomial pieces
!           k      - order of b-spline, k >= 1
!           id     - order of the spline derivative, 0 <= id <= k-1
!                    id=0 gives the spline function
!           x1,x2  - end points of quadrature interval, normally in
!                    xi(1) <= x <= xi(lxi+1)
!           tol    - desired accuracy for the quadrature, suggest
!                    10.*dtol < tol <= 0.1 where dtol is the
!                    maximum of 1.0d-18 and double precision unit
!                    roundoff for the machine = d1mach(4)
!
!         output     quad is double precision
!           quad   - integral of pf(x) on (x1,x2)
!           ierr   - a status code
!                    ierr=1 normal return
!                         2 some quadrature does not meet the
!                           requested tolerance
!
!     error conditions
!         improper input is a fatal error.
!         some quadrature does not meet the requested tolerance.
!***references  d.e. amos, *quadrature subroutines for splines and
!                 b-splines*, sand79-1825, sandia laboratories,
!                 december 1979.
!***routines called  d1mach,dintrv,dppgq8,xerror
!***end prologue  dpfqad
!
!
      integer id,ierr,iflg,ilo,il1,il2,inppv,k,ldc,left,lxi,mf1,mf2
      double precision a,aa,ans,b,bb,c,q,quad,ta,tb,tol,wtol,xi,x1,x2
      double precision f
      dimension xi(*), c(ldc,*)
      external f
!
!***first executable statement  dpfqad
      ierr = 1
      quad = 0.0d0
      if(k<1) go to 100
      if(ldc<k) go to 105
      if(id<0 .or. id>=k) go to 110
      if(lxi<1) go to 115
      wtol = d1mach(4)
      wtol = dmax1(wtol,1.0d-18)
      if (tol<wtol .or. tol>0.1d0) go to 20
      aa = dmin1(x1,x2)
      bb = dmax1(x1,x2)
      if (aa==bb) return
      ilo = 1
      call dintrv(xi, lxi, aa, ilo, il1, mf1)
      call dintrv(xi, lxi, bb, ilo, il2, mf2)
      q = 0.0d0
      inppv = 1
      do 10 left=il1,il2
        ta = xi(left)
        a = dmax1(aa,ta)
        if (left==1) a = aa
        tb = bb
        if (left<lxi) tb = xi(left+1)
        b = dmin1(bb,tb)
        call dppgq8(f,ldc,c,xi,lxi,k,id,a,b,inppv,tol,ans,iflg)
        if (iflg>1) ierr = 2
        q = q + ans
   10 continue
      if (x1>x2) q = -q
      quad = q
      return
!
   20 continue
      call xerror( ' dpfqad,  tol is less dtol or greater than 0.1',&
       46, 2, 1)
      return
  100 continue
      call xerror( ' dpfqad,  k does not satisfy k>=1', 35, 2, 1)
      return
  105 continue
      call xerror( ' dpfqad,  ldc does not satisfy ldc>=k', 39, 2, 1)
      return
  110 continue
      call xerror( ' dpfqad,  id does not satisfy 0<=id<k', 42,&
       2, 1)
      return
  115 continue
      call xerror( ' dpfqad,  lxi does not satisfy lxi>=1', 39, 2, 1)
      return
      end subroutine dpfqad
      
      subroutine dppgq8(fun,ldc,c,xi,lxi,kk,id,a,b,inppv,err,ans,ierr)
!***begin prologue  dppgq8
!***revision history  (yymmdd)
!   000330  modified array declarations.  (jec)
!
!***refer to  dpfqad
!
!     written by r.e. jones and modified by d.e. amos
!
!     abstract    **** a double precision routine ****
!
!        dppgq8, a modification of gaus8, integrates the
!        product of fun(x) by the id-th derivative of a spline
!        dppval(ldc,c,xi,lxi,kk,id,x,inppv)  between limits a and b.
!
!        dppgq8 calls dppval, dintrv, i1mach, d1mach, xerror
!
!     description of arguments
!
!      input-- fun,c,xi,a,b,err are double precision
!        fun - name of external function of one argument which
!              multiplies dppval.
!        ldc - leading dimension of matrix c, ldc >= kk
!        c   - matrix of tayor derivatives of dimension at least
!              (k,lxi)
!        xi  - breakpoint vector of length lxi+1
!        lxi - number of polynomial pieces
!        kk  - order of the spline, kk >= 1
!        id  - order of the spline derivative, 0 <= id <= kk-1
!        a   - lower limit of integral
!        b   - upper limit of integral (may be less than a)
!        inppv- initialization parameter for dppval
!        err - is a requested pseudorelative error tolerance.  normally
!              pick a value of dabs(err) < 1d-3.  ans will normally
!              have no more error than dabs(err) times the integral of
!              the absolute value of fun(x)*dppval(ldc,c,xi,lxi,kk,id,x,
!              inppv).
!
!
!      output-- err,ans are double precision
!        err - will be an estimate of the absolute error in ans if the
!              input value of err was negative.  (err is unchanged if
!              the input value of err was nonnegative.)  the estimated
!              error is solely for information to the user and should
!              not be used as a correction to the computed integral.
!        ans - computed value of integral
!        ierr- a status code
!            --normal codes
!               1 ans most likely meets requested error tolerance,
!                 or a=b.
!              -1 a and b are too nearly equal to allow normal
!                 integration.  ans is set to zero.
!            --abnormal code
!               2 ans probably does not meet requested error tolerance.
!***routines called  d1mach,dppval,i1mach,xerror
!***end prologue  dppgq8
!
      integer icall,id,ierr,inppv,k,kk,kml,kmx,l,ldc,lmn,lmx,lr,lxi,mxl,&
       nbits, nib, nlmn, nlmx
      double precision a,aa,ae,anib,ans,area,b,be,c,cc,ee,ef,eps,err,&
       est,gl,glr,gr,hh,sq2,tol,vl,vr,w1, w2, w3, w4, xi, x1,&
       x2, x3, x4, x, h
      double precision,external :: fun
      dimension xi(*), c(ldc,*)
      dimension aa(60), hh(60), lr(60), vl(60), gr(60)
      data x1, x2, x3, x4/&
           1.83434642495649805d-01,     5.25532409916328986d-01,&
           7.96666477413626740d-01,     9.60289856497536232d-01/
      data w1, w2, w3, w4/&
           3.62683783378361983d-01,     3.13706645877887287d-01,&
           2.22381034453374471d-01,     1.01228536290376259d-01/
      data icall  /  0  /
      data sq2/1.41421356d0/
      data nlmn/1/,kmx/5000/,kml/6/
!
!     initialize
!
!***first executable statement  dppgq8
      if (icall/=0) call xerror(   'dppgq8- dppgq8 called recursively.&
        recursive calls are illegal in fortran.', 75, 7, 2)
      icall = 1
      k = DIGITS(1.0d0)  !i1mach(14)
      anib = d1mach(5)*dble(float(k))/0.30102000d0
      nbits = int(sngl(anib))
      nlmx = min0((nbits*5)/8,60)
      ans = 0.0d0
      ierr = 1
      be = 0.0d0
      if (a==b) go to 140
      lmx = nlmx
      lmn = nlmn
      if (b==0.0d0) go to 10
      if (dsign(1.0d0,b)*a<=0.0d0) go to 10
      cc = dabs(1.0d0-a/b)
      if (cc>0.1d0) go to 10
      if (cc<=0.0d0) go to 140
      anib = 0.5d0 - dlog(cc)/0.69314718d0
      nib = int(sngl(anib))
      lmx = min0(nlmx,nbits-nib-7)
      if (lmx<1) go to 130
      lmn = min0(lmn,lmx)
   10 tol = dmax1(dabs(err),2.0d0**(5-nbits))/2.0d0
      if (err==0.0d0) tol = dsqrt(d1mach(4))
      eps = tol
      hh(1) = (b-a)/4.0d0
      aa(1) = a
      lr(1) = 1
      l = 1
      est = g8(aa(l)+2.0d0*hh(l),2.0d0*hh(l),fun)
      k = 8
      area = dabs(est)
      ef = 0.5d0
      mxl = 0
!
!     compute refined estimates, estimate the error, etc.
!
   20 gl = g8(aa(l)+hh(l),hh(l),fun)
      gr(l) = g8(aa(l)+3.0d0*hh(l),hh(l),fun)
      k = k + 16
      area = area + (dabs(gl)+dabs(gr(l))-dabs(est))
      glr = gl + gr(l)
      ee = dabs(est-glr)*ef
      ae = dmax1(eps*area,tol*dabs(glr))
      if (ee-ae) 40, 40, 50
   30 mxl = 1
   40 be = be + (est-glr)
      if (lr(l)) 60, 60, 80
!
!     consider the left half of this level
!
   50 if (k>kmx) lmx = kml
      if (l>=lmx) go to 30
      l = l + 1
      eps = eps*0.5d0
      ef = ef/sq2
      hh(l) = hh(l-1)*0.5d0
      lr(l) = -1
      aa(l) = aa(l-1)
      est = gl
      go to 20
!
!     proceed to right half at this level
!
   60 vl(l) = glr
   70 est = gr(l-1)
      lr(l) = 1
      aa(l) = aa(l) + 4.0d0*hh(l)
      go to 20
!
!     return one level
!
   80 vr = glr
   90 if (l<=1) go to 120
      l = l - 1
      eps = eps*2.0d0
      ef = ef*sq2
      if (lr(l)) 100, 100, 110
  100 vl(l) = vl(l+1) + vr
      go to 70
  110 vr = vl(l+1) + vr
      go to 90
!
!      exit
!
  120 ans = vr
      if ((mxl==0) .or. (dabs(be)<=2.0d0*tol*area)) go to 140
      ierr = 2
      call xerror( 'dppgq8- ans is probably insufficiently accurate.',&
       48, 3, 1)
      go to 140
  130 ierr = -1
      call xerror( 'dppgq8- the following temporary diagnostic will appe&
      ar only once.  a and b are too nearly equal to allow normal integr&
      ation.  ans is set to zero, and ierr=-1.', 158, 1, -1)
  140 icall = 0
      if (err<0.0d0) err = be
      return
	  
	  contains
	  
	  function g8(x,h,fun)  !replaces the statement function in the original code
	  
	  implicit none
	  
	  double precision,intent(in) :: x
	  double precision,intent(in) :: h
	  double precision :: g8
	  double precision,external :: fun
	  
	  g8 = h*((w1*(fun(x-x1*h)*dppval(ldc,c,xi,lxi,kk,id,x-x1*h,inppv&
                 ) +fun(x+x1*h)*dppval(ldc,c,xi,lxi,kk,id,x+x1*h,inppv))&
                +w2*(fun(x-x2*h)*dppval(ldc,c,xi,lxi,kk,id,x-x2*h,inppv)&
                  +fun(x+x2*h)*dppval(ldc,c,xi,lxi,kk,id,x+x2*h,inppv)))&
              +(w3*(fun(x-x3*h)*dppval(ldc,c,xi,lxi,kk,id,x-x3*h,inppv)&
                   +fun(x+x3*h)*dppval(ldc,c,xi,lxi,kk,id,x+x3*h,inppv))&
               +w4*(fun(x-x4*h)*dppval(ldc,c,xi,lxi,kk,id,x-x4*h,inppv)&
                +fun(x+x4*h)*dppval(ldc,c,xi,lxi,kk,id,x+x4*h,inppv))))

	  end function g8
	  
      end subroutine dppgq8
      
      subroutine dppqad(ldc,c,xi,lxi,k,x1,x2,pquad)
!***begin prologue  dppqad
!***date written   800901   (yymmdd)
!***revision date  820801   (yymmdd)
!***revision history  (yymmdd)
!   000330  modified array declarations.  (jec)
!
!***category no.  h2a2a1,e3,k6
!***keywords  b-spline,data fitting,double precision,interpolation,
!             quadrature,spline
!***author  amos, d. e., (snla)
!***purpose  computes the integral on (x1,x2) of a k-th order b-spline
!            using the piecewise polynomial representation.
!***description
!
!     written by d. e. amos, june, 1979.
!
!     reference sand-79-1825
!
!     abstract    **** a double precision routine ****
!
!         dppqad computes the integral on (x1,x2) of a k-th order
!         b-spline using the piecewise polynomial representation
!         (c,xi,lxi,k).  here the taylor expansion about the left
!         end point xi(j) of the j-th interval is integrated and
!         evaluated on subintervals of (x1,x2) which are formed by
!         included break points.  integration outside (xi(1),xi(lxi+1))
!         is permitted.
!
!         dppqad calls dintrv
!
!     description of arguments
!         input      c,xi,x1,x2 are double precision
!           ldc    - leading dimension of matrix c, ldc >= k
!           c(i,j) - right taylor derivatives at xi(j), i=1,k , j=1,lxi
!           xi(*)  - break point array of length lxi+1
!           lxi    - number of polynomial pieces
!           k      - order of b-spline, k >= 1
!           x1,x2  - end points of quadrature interval, normally in
!                    xi(1) <= x <= xi(lxi+1)
!
!         output     pquad is double precision
!           pquad  - integral of the pp representation over (x1,x2)
!
!     error conditions
!         improper input is a fatal error
!***references  d.e. amos, *quadrature subroutines for splines and
!                 b-splines*, sand79-1825, sandia laboratories,
!                 december 1979.
!***routines called  dintrv,xerror
!***end prologue  dppqad
!
!
      integer i, ii, il, ilo, il1, il2, im, k, ldc, left, lxi, mf1, mf2
      double precision a,aa,bb,c,dx,flk,pquad,q,s,ss,ta,tb,x,xi,x1,x2
      dimension xi(*), c(ldc,*), ss(2)
!
!***first executable statement  dppqad
      pquad = 0.0d0
      if(k<1) go to 100
      if(lxi<1) go to 105
      if(ldc<k) go to 110
      aa = dmin1(x1,x2)
      bb = dmax1(x1,x2)
      if (aa==bb) return
      ilo = 1
      call dintrv(xi, lxi, aa, ilo, il1, mf1)
      call dintrv(xi, lxi, bb, ilo, il2, mf2)
      q = 0.0d0
      do 40 left=il1,il2
        ta = xi(left)
        a = dmax1(aa,ta)
        if (left==1) a = aa
        tb = bb
        if (left<lxi) tb = xi(left+1)
        x = dmin1(bb,tb)
        do 30 ii=1,2
          ss(ii) = 0.0d0
          dx = x - xi(left)
          if (dx==0.0d0) go to 20
          s = c(k,left)
          flk = dble(float(k))
          im = k - 1
          il = im
          do 10 i=1,il
            s = s*dx/flk + c(im,left)
            im = im - 1
            flk = flk - 1.0d0
   10     continue
          ss(ii) = s*dx
   20     continue
          x = a
   30   continue
        q = q + (ss(1)-ss(2))
   40 continue
      if (x1>x2) q = -q
      pquad = q
      return
!
!
  100 continue
      call xerror( ' dppqad,  k does not satisfy k>=1', 35, 2, 1)
      return
  105 continue
      call xerror( ' dppqad,  lxi does not satisfy lxi>=1', 39, 2, 1)
      return
  110 continue
      call xerror( ' dppqad,  ldc does not satisfy ldc>=k', 39, 2, 1)
      return
      end subroutine dppqad
      
      double precision function dppval(ldc,c,xi,lxi,k,ideriv,x,inppv)
!***begin prologue  dppval
!***date written   800901   (yymmdd)
!***revision date  820801   (yymmdd)
!***revision history  (yymmdd)
!   000330  modified array declarations.  (jec)
!
!***category no.  e3,k6
!***keywords  b-spline,data fitting,double precision,interpolation,
!             spline
!***author  amos, d. e., (snla)
!***purpose  calculates (at x) the value of the ideriv-th derivative
!            of the b-spline from the pp-representation.
!***description
!
!     written by carl de boor and modified by d. e. amos
!
!     reference
!         siam j. numerical analysis, 14, no. 3, june, 1977, pp.441-472.
!
!     abstract    **** a double precision routine ****
!         dppval is the ppvalu function of the reference.
!
!         dppval calculates (at x) the value of the ideriv-th
!         derivative of the b-spline from the pp-representation
!         (c,xi,lxi,k).  the taylor expansion about xi(j) for x in
!         the interval xi(j) <= x < xi(j+1) is evaluated, j=1,lxi.
!         right limiting values at x=xi(j) are obtained.  dppval will
!         extrapolate beyond xi(1) and xi(lxi+1).
!
!         to obtain left limiting values (left derivatives) at xi(j)
!         replace lxi by j-1 and set x=xi(j),j=2,lxi+1.
!
!         dppval calls dintrv
!
!     description of arguments
!
!         input      c,xi,x are double precision
!          ldc     - leading dimension of c matrix, ldc >= k
!          c       - matrix of dimension at least (k,lxi) containing
!                    right derivatives at break points xi(*).
!          xi      - break point vector of length lxi+1
!          lxi     - number of polynomial pieces
!          k       - order of b-spline, k >= 1
!          ideriv  - order of the derivative, 0 <= ideriv <= k-1
!                    ideriv=0 gives the b-spline value
!          x       - argument, xi(1) <= x <= xi(lxi+1)
!          inppv   - an initialization parameter which must be set
!                    to 1 the first time dppval is called.
!
!         output     dppval is double precision
!          inppv   - inppv contains information for efficient process-
!                    ing after the initial call and inppv must not
!                    be changed by the user.  distinct splines require
!                    distinct inppv parameters.
!          dppval  - value of the ideriv-th derivative at x
!
!     error conditions
!         improper input is a fatal error
!***references  c. de boor, *package for calculating with b-splines*,
!                 siam journal on numerical analysis, volume 14, no. 3,
!                 june 1977, pp. 441-472.
!***routines called  dintrv,xerror
!***end prologue  dppval
!
!
      integer i, ideriv, inppv, j, k, ldc, lxi, ndummy
      double precision c, dx, fltk, x, xi
      dimension xi(*), c(ldc,lxi)
!***first executable statement  dppval
      dppval = 0.0d0
      if(k<1) go to 90
      if(ldc<k) go to 80
      if(lxi<1) go to 85
      if(ideriv<0 .or. ideriv>=k) go to 95
      i = k - ideriv
      fltk = dble(float(i))
      call dintrv(xi, lxi, x, inppv, i, ndummy)
      dx = x - xi(i)
      j = k
   10 dppval = (dppval/fltk)*dx + c(j,i)
      j = j - 1
      fltk = fltk - 1.0d0
      if (fltk>0.0d0) go to 10
      return
!
!
   80 continue
      call xerror( ' dppval,  ldc does not satisfy ldc>=k', 39, 2, 1)
      return
   85 continue
      call xerror( ' dppval,  lxi does not satisfy lxi>=1', 39, 2, 1)
      return
   90 continue
      call xerror( ' dppval,  k does not satisfy k>=1', 35, 2, 1)
      return
   95 continue
      call xerror( ' dppval,  ideriv does not satisfy 0<=ideriv<k',&
       50, 2, 1)
      return
      end function dppval

!*****************************************************************************************
    pure function d1mach (i) result(d)
    
    implicit none
    
    double precision :: d
    integer,intent(in) :: i
    double precision :: b, x
    
    x = 1.0d0
    b = radix(x)
    select case (i)
    case (1)
      d = b**(minexponent(x)-1) ! smallest positive magnitude
    case (2) 
      d = huge(x)               ! largest magnitude
    case (3) 
      d = b**(-digits(x))       ! smallest relative spacing
    case (4) 
      d = b**(1-digits(x))      ! largest relative spacing
    case (5)
      d = log10(b)
    end select
    
    end function d1mach     
!*****************************************************************************************
    
!*****************************************************************************************
    subroutine xerror(messg,nmessg,nerr,level)

    implicit none

    character(len=*),intent(in) :: messg
    integer,intent(in) :: nmessg,nerr,level
    
    write(*,'(A)') trim(messg)

    end subroutine xerror
!*****************************************************************************************
            
!*****************************************************************************************
    end module bspline_module
!*****************************************************************************************