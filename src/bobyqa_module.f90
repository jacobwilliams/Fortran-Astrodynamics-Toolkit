!*****************************************************************************************
    module bobyqa_module
!*****************************************************************************************
!****h* FAT/bobyqa_module
!
!  NAME
!    bobyqa_module
!
!  DESCRIPTION
!    The BOBYQA algorithm by M.J.D. Powell.
!    Its purpose is to seek the least value of a function of 
!      several variables, when derivatives are not available.  
!    The name BOBYQA denotes Bound Approximation BY Quadratic
!      Approximation, the constraints being lower and upper bounds on every
!      variable, which can be set to huge values for unconstrained variables.
!
!  SEE ALSO
!    http://www.damtp.cam.ac.uk/user/na/NA_papers/NA2009_06.pdf
!
!*****************************************************************************************    
    
    use kind_module,  only: wp
    
    private
    
    public :: bobyqa
    public :: bobyqa_test
    
    abstract interface
        subroutine func(n,x,f)
            import :: wp
            implicit none
            integer,intent(in) :: n
            real(wp),dimension(n),intent(in) :: x
            real(wp),intent(out) :: f
        end subroutine func
    end interface
    
    procedure(func),pointer :: calfun => null()
    
    contains
!*****************************************************************************************

!*****************************************************************************************
! test problem for bobyqa, the objective function being the sum of
! the reciprocals of all pairwise distances between the points p_i,
! i=1,2,...,m in two dimensions, where m=n/2 and where the components
! of p_i are x(2*i-1) and x(2*i). thus each vector x of n variables
! defines the m points p_i. the initial x gives equally spaced points
! on a circle. four different choices of the pairs (n,npt) are tried,
! namely (10,16), (10,21), (20,26) and (20,41). convergence to a local
! minimum that is not global occurs in both the n=10 cases. the details
! of the results are highly sensitive to computer rounding errors. the
! choice iprint=2 provides the current x and optimal f so far whenever
! rho is reduced. the bound constraints of the problem require every
! component of x to be in the interval [-1,1].

    subroutine bobyqa_test()

      implicit real(wp) (a-h,o-z)
      dimension x(100),xl(100),xu(100)
    
      calfun => calfun_test        !set the function
    
      twopi=8.0d0*atan(1.0d0)
      bdl=-1.0d0
      bdu=1.0d0
      iprint=2
      maxfun=500000
      rhobeg=1.0d-1
      rhoend=1.0d-6
      m=5
      do
          n=2*m
          do i=1,n
            xl(i)=bdl
            xu(i)=bdu
          end do
          do jcase=1,2
              npt=n+6
              if (jcase == 2) npt=2*n+1
              print 30, m,n,npt
    30        format (//5x,'2d output with m =',i4,',  n =',i4,'  and  npt =',i4)
              do j=1,m
                temp=dble(j)*twopi/dble(m)
                x(2*j-1)=cos(temp)
                x(2*j)=sin(temp)
              end do
              call bobyqa (n,npt,x,xl,xu,rhobeg,rhoend,iprint,maxfun)
          end do
          m=m+m
          if (m > 10) exit
      end do
      
    end subroutine bobyqa_test
!*****************************************************************************************

!*****************************************************************************************
      subroutine calfun_test (n,x,f)

      implicit none
      integer,intent(in) :: n
      real(wp),dimension(n),intent(in) :: x
      real(wp),intent(out) :: f
      
      integer :: i,j
      real(wp) :: temp

      f=0.0d0
      do i=4,n,2
        do j=2,i-2,2
          temp=(x(i-1)-x(j-1))**2+(x(i)-x(j))**2
          temp=max(temp,1.0d-6)
          f=f+1.0d0/sqrt(temp)
        end do
      end do
      return
      end subroutine calfun_test
!*****************************************************************************************


!*****************************************************************************************
      subroutine bobyqa (n,npt,x,xl,xu,rhobeg,rhoend,iprint,maxfun)
      implicit real(wp) (a-h,o-z)
      dimension x(*),xl(*),xu(*)
      
      real(wp),dimension((npt+5)*(npt+n)+3*n*(n+5)/2) :: w

!
!     this subroutine seeks the least value of a function of many variables,
!     by applying a trust region method that forms quadratic models by
!     interpolation. there is usually some freedom in the interpolation
!     conditions, which is taken up by minimizing the frobenius norm of
!     the change to the second derivative of the model, beginning with the
!     zero matrix. the values of the variables are constrained by upper and
!     lower bounds. the arguments of the subroutine are as follows.
!
!     n must be set to the number of variables and must be at least two.
!     npt is the number of interpolation conditions. its value must be in
!       the interval [n+2,(n+1)(n+2)/2]. choices that exceed 2*n+1 are not
!       recommended.
!     initial values of the variables must be set in x(1),x(2),...,x(n). they
!       will be changed to the values that give the least calculated f.
!     for i=1,2,...,n, xl(i) and xu(i) must provide the lower and upper
!       bounds, respectively, on x(i). the construction of quadratic models
!       requires xl(i) to be strictly less than xu(i) for each i. further,
!       the contribution to a model from changes to the i-th variable is
!       damaged severely by rounding errors if xu(i)-xl(i) is too small.
!     rhobeg and rhoend must be set to the initial and final values of a trust
!       region radius, so both must be positive with rhoend no greater than
!       rhobeg. typically, rhobeg should be about one tenth of the greatest
!       expected change to a variable, while rhoend should indicate the
!       accuracy that is required in the final values of the variables. an
!       error return occurs if any of the differences xu(i)-xl(i), i=1,...,n,
!       is less than 2*rhobeg.
!     the value of iprint should be set to 0, 1, 2 or 3, which controls the
!       amount of printing. specifically, there is no output if iprint=0 and
!       there is output only at the return if iprint=1. otherwise, each new
!       value of rho is printed, with the best vector of variables so far and
!       the corresponding value of the objective function. further, each new
!       value of f with its variables are output if iprint=3.
!     maxfun must be set to an upper bound on the number of calls of calfun.
!     the array w will be used for working space. its length must be at least
!       (npt+5)*(npt+n)+3*n*(n+5)/2.
!
!     subroutine calfun (n,x,f) has to be provided by the user. it must set
!     f to the value of the objective function for the current values of the
!     variables x(1),x(2),...,x(n), which are generated automatically in a
!     way that satisfies the bounds given in xl and xu.
!
!     return if the value of npt is unacceptable.
!
      np=n+1
      if (npt < n+2 .or. npt > ((n+2)*np)/2) then
          print 10
   10     format (/4x,'return from bobyqa because npt is not in',&
            ' the required interval')
          go to 40
      end if
!
!     partition the working space array, so that different parts of it can
!     be treated separately during the calculation of bobyqb. the partition
!     requires the first (npt+2)*(npt+n)+3*n*(n+5)/2 elements of w plus the
!     space that is taken by the last array in the argument list of bobyqb.
!
      ndim=npt+n
      ixb=1
      ixp=ixb+n
      ifv=ixp+n*npt
      ixo=ifv+npt
      igo=ixo+n
      ihq=igo+n
      ipq=ihq+(n*np)/2
      ibmat=ipq+npt
      izmat=ibmat+ndim*n
      isl=izmat+npt*(npt-np)
      isu=isl+n
      ixn=isu+n
      ixa=ixn+n
      id=ixa+n
      ivl=id+n
      iw=ivl+ndim
!
!     return if there is insufficient space between the bounds. modify the
!     initial x if necessary in order to avoid conflicts between the bounds
!     and the construction of the first quadratic model. the lower and upper
!     bounds on moves from the updated x are set now, in the isl and isu
!     partitions of w, in order to provide useful and exact information about
!     components of x that become within distance rhobeg from their bounds.
!
      zero=0.0d0
      do 30 j=1,n
      temp=xu(j)-xl(j)
      if (temp < rhobeg+rhobeg) then
          print 20
   20     format (/4x,'return from bobyqa because one of the',&
            ' differences xu(i)-xl(i)'/6x,' is less than 2*rhobeg.')
          go to 40
      end if
      jsl=isl+j-1
      jsu=jsl+n
      w(jsl)=xl(j)-x(j)
      w(jsu)=xu(j)-x(j)
      if (w(jsl) >= -rhobeg) then
          if (w(jsl) >= zero) then
              x(j)=xl(j)
              w(jsl)=zero
              w(jsu)=temp
          else
              x(j)=xl(j)+rhobeg
              w(jsl)=-rhobeg
              w(jsu)=max(xu(j)-x(j),rhobeg)
          end if
      else if (w(jsu) <= rhobeg) then
          if (w(jsu) <= zero) then
              x(j)=xu(j)
              w(jsl)=-temp
              w(jsu)=zero
          else
              x(j)=xu(j)-rhobeg
              w(jsl)=min(xl(j)-x(j),-rhobeg)
              w(jsu)=rhobeg
          end if
      end if
   30 continue
!
!     make the call of bobyqb.
!
      call bobyqb (n,npt,x,xl,xu,rhobeg,rhoend,iprint,maxfun,w(ixb),&
        w(ixp),w(ifv),w(ixo),w(igo),w(ihq),w(ipq),w(ibmat),w(izmat),&
        ndim,w(isl),w(isu),w(ixn),w(ixa),w(id),w(ivl),w(iw))
   40 return
      end subroutine bobyqa
!*****************************************************************************************


!*****************************************************************************************
      subroutine bobyqb (n,npt,x,xl,xu,rhobeg,rhoend,iprint,&
        maxfun,xbase,xpt,fval,xopt,gopt,hq,pq,bmat,zmat,ndim,&
        sl,su,xnew,xalt,d,vlag,w)
      implicit real(wp) (a-h,o-z)
      dimension x(*),xl(*),xu(*),xbase(*),xpt(npt,*),fval(*),&
        xopt(*),gopt(*),hq(*),pq(*),bmat(ndim,*),zmat(npt,*),&
        sl(*),su(*),xnew(*),xalt(*),d(*),vlag(*),w(*)
!
!     the arguments n, npt, x, xl, xu, rhobeg, rhoend, iprint and maxfun
!       are identical to the corresponding arguments in subroutine bobyqa.
!     xbase holds a shift of origin that should reduce the contributions
!       from rounding errors to values of the model and lagrange functions.
!     xpt is a two-dimensional array that holds the coordinates of the
!       interpolation points relative to xbase.
!     fval holds the values of f at the interpolation points.
!     xopt is set to the displacement from xbase of the trust region centre.
!     gopt holds the gradient of the quadratic model at xbase+xopt.
!     hq holds the explicit second derivatives of the quadratic model.
!     pq contains the parameters of the implicit second derivatives of the
!       quadratic model.
!     bmat holds the last n columns of h.
!     zmat holds the factorization of the leading npt by npt submatrix of h,
!       this factorization being zmat times zmat^t, which provides both the
!       correct rank and positive semi-definiteness.
!     ndim is the first dimension of bmat and has the value npt+n.
!     sl and su hold the differences xl-xbase and xu-xbase, respectively.
!       all the components of every xopt are going to satisfy the bounds
!       sl(i) .leq. xopt(i) .leq. su(i), with appropriate equalities when
!       xopt is on a constraint boundary.
!     xnew is chosen by subroutine trsbox or altmov. usually xbase+xnew is the
!       vector of variables for the next call of calfun. xnew also satisfies
!       the sl and su constraints in the way that has just been mentioned.
!     xalt is an alternative to xnew, chosen by altmov, that may replace xnew
!       in order to increase the denominator in the updating of update.
!     d is reserved for a trial step from xopt, which is usually xnew-xopt.
!     vlag contains the values of the lagrange functions at a new point x.
!       they are part of a product that requires vlag to be of length ndim.
!     w is a one-dimensional array that is used for working space. its length
!       must be at least 3*ndim = 3*(npt+n).
!
!     set some constants.
!
      half=0.5d0
      one=1.0d0
      ten=10.0d0
      tenth=0.1d0
      two=2.0d0
      zero=0.0d0
      np=n+1
      nptm=npt-np
      nh=(n*np)/2
!
!     the call of prelim sets the elements of xbase, xpt, fval, gopt, hq, pq,
!     bmat and zmat for the first iteration, with the corresponding values of
!     of nf and kopt, which are the number of calls of calfun so far and the
!     index of the interpolation point at the trust region centre. then the
!     initial xopt is set too. the branch to label 720 occurs if maxfun is
!     less than npt. gopt will be updated if kopt is different from kbase.
!
      call prelim (n,npt,x,xl,xu,rhobeg,iprint,maxfun,xbase,xpt,&
        fval,gopt,hq,pq,bmat,zmat,ndim,sl,su,nf,kopt)
      xoptsq=zero
      do 10 i=1,n
      xopt(i)=xpt(kopt,i)
   10 xoptsq=xoptsq+xopt(i)**2
      fsave=fval(1)
      if (nf < npt) then
          if (iprint > 0) print 390
          goto 720
      end if
      kbase=1
!
!     complete the settings that are required for the iterative procedure.
!
      rho=rhobeg
      delta=rho
      nresc=nf
      ntrits=0
      diffa=zero
      diffb=zero
      itest=0
      nfsav=nf
!
!     update gopt if necessary before the first iteration and after each
!     call of rescue that makes a call of calfun.
!
   20 if (kopt /= kbase) then
          ih=0
          do 30 j=1,n
          do 30 i=1,j
          ih=ih+1
          if (i < j) gopt(j)=gopt(j)+hq(ih)*xopt(i)
   30     gopt(i)=gopt(i)+hq(ih)*xopt(j)
          if (nf > npt) then
              do 50 k=1,npt
              temp=zero
              do 40 j=1,n
   40         temp=temp+xpt(k,j)*xopt(j)
              temp=pq(k)*temp
              do 50 i=1,n
   50         gopt(i)=gopt(i)+temp*xpt(k,i)
          end if
      end if
!
!     generate the next point in the trust region that provides a small value
!     of the quadratic model subject to the constraints on the variables.
!     the integer ntrits is set to the number "trust region" iterations that
!     have occurred since the last "alternative" iteration. if the length
!     of xnew-xopt is less than half*rho, however, then there is a branch to
!     label 650 or 680 with ntrits=-1, instead of calculating f at xnew.
!
   60 call trsbox (n,npt,xpt,xopt,gopt,hq,pq,sl,su,delta,xnew,d,&
        w,w(np),w(np+n),w(np+2*n),w(np+3*n),dsq,crvmin)
      dnorm=min(delta,sqrt(dsq))
      if (dnorm < half*rho) then
          ntrits=-1
          distsq=(ten*rho)**2
          if (nf <= nfsav+2) goto 650
!
!     the following choice between labels 650 and 680 depends on whether or
!     not our work with the current rho seems to be complete. either rho is
!     decreased or termination occurs if the errors in the quadratic model at
!     the last three interpolation points compare favourably with predictions
!     of likely improvements to the model within distance half*rho of xopt.
!
          errbig=max(diffa,diffb,diffc)
          frhosq=0.125d0*rho*rho
          if (crvmin > zero .and. errbig > frhosq*crvmin)&
             goto 650
          bdtol=errbig/rho
          do 80 j=1,n
          bdtest=bdtol
          if (xnew(j) == sl(j)) bdtest=w(j)
          if (xnew(j) == su(j)) bdtest=-w(j)
          if (bdtest < bdtol) then
              curv=hq((j+j*j)/2)
              do 70 k=1,npt
   70         curv=curv+pq(k)*xpt(k,j)**2
              bdtest=bdtest+half*curv*rho
              if (bdtest < bdtol) goto 650
          end if
   80     continue
          goto 680
      end if
      ntrits=ntrits+1
!
!     severe cancellation is likely to occur if xopt is too far from xbase.
!     if the following test holds, then xbase is shifted so that xopt becomes
!     zero. the appropriate changes are made to bmat and to the second
!     derivatives of the current model, beginning with the changes to bmat
!     that do not depend on zmat. vlag is used temporarily for working space.
!
   90 if (dsq <= 1.0d-3*xoptsq) then
          fracsq=0.25d0*xoptsq
          sumpq=zero
          do 110 k=1,npt
          sumpq=sumpq+pq(k)
          sum=-half*xoptsq
          do 100 i=1,n
  100     sum=sum+xpt(k,i)*xopt(i)
          w(npt+k)=sum
          temp=fracsq-half*sum
          do 110 i=1,n
          w(i)=bmat(k,i)
          vlag(i)=sum*xpt(k,i)+temp*xopt(i)
          ip=npt+i
          do 110 j=1,i
  110     bmat(ip,j)=bmat(ip,j)+w(i)*vlag(j)+vlag(i)*w(j)
!
!     then the revisions of bmat that depend on zmat are calculated.
!
          do 150 jj=1,nptm
          sumz=zero
          sumw=zero
          do 120 k=1,npt
          sumz=sumz+zmat(k,jj)
          vlag(k)=w(npt+k)*zmat(k,jj)
  120     sumw=sumw+vlag(k)
          do 140 j=1,n
          sum=(fracsq*sumz-half*sumw)*xopt(j)
          do 130 k=1,npt
  130     sum=sum+vlag(k)*xpt(k,j)
          w(j)=sum
          do 140 k=1,npt
  140     bmat(k,j)=bmat(k,j)+sum*zmat(k,jj)
          do 150 i=1,n
          ip=i+npt
          temp=w(i)
          do 150 j=1,i
  150     bmat(ip,j)=bmat(ip,j)+temp*w(j)
!
!     the following instructions complete the shift, including the changes
!     to the second derivative parameters of the quadratic model.
!
          ih=0
          do 170 j=1,n
          w(j)=-half*sumpq*xopt(j)
          do 160 k=1,npt
          w(j)=w(j)+pq(k)*xpt(k,j)
  160     xpt(k,j)=xpt(k,j)-xopt(j)
          do 170 i=1,j
          ih=ih+1
          hq(ih)=hq(ih)+w(i)*xopt(j)+xopt(i)*w(j)
  170     bmat(npt+i,j)=bmat(npt+j,i)
          do 180 i=1,n
          xbase(i)=xbase(i)+xopt(i)
          xnew(i)=xnew(i)-xopt(i)
          sl(i)=sl(i)-xopt(i)
          su(i)=su(i)-xopt(i)
  180     xopt(i)=zero
          xoptsq=zero
      end if
      if (ntrits == 0) goto 210
      goto 230
!
!     xbase is also moved to xopt by a call of rescue. this calculation is
!     more expensive than the previous shift, because new matrices bmat and
!     zmat are generated from scratch, which may include the replacement of
!     interpolation points whose positions seem to be causing near linear
!     dependence in the interpolation conditions. therefore rescue is called
!     only if rounding errors have reduced by at least a factor of two the
!     denominator of the formula for updating the h matrix. it provides a
!     useful safeguard, but is not invoked in most applications of bobyqa.
!
  190 nfsav=nf
      kbase=kopt
      call rescue (n,npt,xl,xu,iprint,maxfun,xbase,xpt,fval,&
        xopt,gopt,hq,pq,bmat,zmat,ndim,sl,su,nf,delta,kopt,&
        vlag,w,w(n+np),w(ndim+np))
!
!     xopt is updated now in case the branch below to label 720 is taken.
!     any updating of gopt occurs after the branch below to label 20, which
!     leads to a trust region iteration as does the branch to label 60.
!
      xoptsq=zero
      if (kopt /= kbase) then
          do 200 i=1,n
          xopt(i)=xpt(kopt,i)
  200     xoptsq=xoptsq+xopt(i)**2
      end if
      if (nf < 0) then
          nf=maxfun
          if (iprint > 0) print 390
          goto 720
      end if
      nresc=nf
      if (nfsav < nf) then
          nfsav=nf
          goto 20
      end if
      if (ntrits > 0) goto 60
!
!     pick two alternative vectors of variables, relative to xbase, that
!     are suitable as new positions of the knew-th interpolation point.
!     firstly, xnew is set to the point on a line through xopt and another
!     interpolation point that minimizes the predicted value of the next
!     denominator, subject to ||xnew - xopt|| .leq. adelt and to the sl
!     and su bounds. secondly, xalt is set to the best feasible point on
!     a constrained version of the cauchy step of the knew-th lagrange
!     function, the corresponding value of the square of this function
!     being returned in cauchy. the choice between these alternatives is
!     going to be made when the denominator is calculated.
!
  210 call altmov (n,npt,xpt,xopt,bmat,zmat,ndim,sl,su,kopt,&
        knew,adelt,xnew,xalt,alpha,cauchy,w,w(np),w(ndim+1))
      do 220 i=1,n
  220 d(i)=xnew(i)-xopt(i)
!
!     calculate vlag and beta for the current choice of d. the scalar
!     product of d with xpt(k,.) is going to be held in w(npt+k) for
!     use when vquad is calculated.
!
  230 do 250 k=1,npt
      suma=zero
      sumb=zero
      sum=zero
      do 240 j=1,n
      suma=suma+xpt(k,j)*d(j)
      sumb=sumb+xpt(k,j)*xopt(j)
  240 sum=sum+bmat(k,j)*d(j)
      w(k)=suma*(half*suma+sumb)
      vlag(k)=sum
  250 w(npt+k)=suma
      beta=zero
      do 270 jj=1,nptm
      sum=zero
      do 260 k=1,npt
  260 sum=sum+zmat(k,jj)*w(k)
      beta=beta-sum*sum
      do 270 k=1,npt
  270 vlag(k)=vlag(k)+sum*zmat(k,jj)
      dsq=zero
      bsum=zero
      dx=zero
      do 300 j=1,n
      dsq=dsq+d(j)**2
      sum=zero
      do 280 k=1,npt
  280 sum=sum+w(k)*bmat(k,j)
      bsum=bsum+sum*d(j)
      jp=npt+j
      do 290 i=1,n
  290 sum=sum+bmat(jp,i)*d(i)
      vlag(jp)=sum
      bsum=bsum+sum*d(j)
  300 dx=dx+d(j)*xopt(j)
      beta=dx*dx+dsq*(xoptsq+dx+dx+half*dsq)+beta-bsum
      vlag(kopt)=vlag(kopt)+one
!
!     if ntrits is zero, the denominator may be increased by replacing
!     the step d of altmov by a cauchy step. then rescue may be called if
!     rounding errors have damaged the chosen denominator.
!
      if (ntrits == 0) then
          denom=vlag(knew)**2+alpha*beta
          if (denom < cauchy .and. cauchy > zero) then
              do 310 i=1,n
              xnew(i)=xalt(i)
  310         d(i)=xnew(i)-xopt(i)
              cauchy=zero
              go to 230
          end if
          if (denom <= half*vlag(knew)**2) then
              if (nf > nresc) goto 190
              if (iprint > 0) print 320
  320         format (/5x,'return from bobyqa because of much',&
                ' cancellation in a denominator.')
              goto 720
          end if
!
!     alternatively, if ntrits is positive, then set knew to the index of
!     the next interpolation point to be deleted to make room for a trust
!     region step. again rescue may be called if rounding errors have damaged
!     the chosen denominator, which is the reason for attempting to select
!     knew before calculating the next value of the objective function.
!
      else
          delsq=delta*delta
          scaden=zero
          biglsq=zero
          knew=0
          do 350 k=1,npt
          if (k == kopt) goto 350
          hdiag=zero
          do 330 jj=1,nptm
  330     hdiag=hdiag+zmat(k,jj)**2
          den=beta*hdiag+vlag(k)**2
          distsq=zero
          do 340 j=1,n
  340     distsq=distsq+(xpt(k,j)-xopt(j))**2
          temp=max(one,(distsq/delsq)**2)
          if (temp*den > scaden) then
              scaden=temp*den
              knew=k
              denom=den
          end if
          biglsq=max(biglsq,temp*vlag(k)**2)
  350     continue
          if (scaden <= half*biglsq) then
              if (nf > nresc) goto 190
              if (iprint > 0) print 320
              goto 720
          end if
      end if
!
!     put the variables for the next calculation of the objective function
!       in xnew, with any adjustments for the bounds.
!
!
!     calculate the value of the objective function at xbase+xnew, unless
!       the limit on the number of calculations of f has been reached.
!
  360 do 380 i=1,n
      x(i)=min(max(xl(i),xbase(i)+xnew(i)),xu(i))
      if (xnew(i) == sl(i)) x(i)=xl(i)
      if (xnew(i) == su(i)) x(i)=xu(i)
  380 continue
      if (nf >= maxfun) then
          if (iprint > 0) print 390
  390     format (/4x,'return from bobyqa because calfun has been',&
            ' called maxfun times.')
          goto 720
      end if
      nf=nf+1
      call calfun (n,x,f)
      if (iprint == 3) then
          print 400, nf,f,(x(i),i=1,n)
  400      format (/4x,'function number',i6,'    f =',1pd18.10,&
             '    the corresponding x is:'/(2x,5d15.6))
      end if
      if (ntrits == -1) then
          fsave=f
          goto 720
      end if
!
!     use the quadratic model to predict the change in f due to the step d,
!       and set diff to the error of this prediction.
!
      fopt=fval(kopt)
      vquad=zero
      ih=0
      do 410 j=1,n
      vquad=vquad+d(j)*gopt(j)
      do 410 i=1,j
      ih=ih+1
      temp=d(i)*d(j)
      if (i == j) temp=half*temp
  410 vquad=vquad+hq(ih)*temp
      do 420 k=1,npt
  420 vquad=vquad+half*pq(k)*w(npt+k)**2
      diff=f-fopt-vquad
      diffc=diffb
      diffb=diffa
      diffa=abs(diff)
      if (dnorm > rho) nfsav=nf
!
!     pick the next value of delta after a trust region step.
!
      if (ntrits > 0) then
          if (vquad >= zero) then
              if (iprint > 0) print 430
  430         format (/4x,'return from bobyqa because a trust',&
                ' region step has failed to reduce q.')
              goto 720
          end if
          ratio=(f-fopt)/vquad
          if (ratio <= tenth) then
              delta=min(half*delta,dnorm)
          else if (ratio <= 0.7d0) then
              delta=max(half*delta,dnorm)
          else
              delta=max(half*delta,dnorm+dnorm)
          end if
          if (delta <= 1.5d0*rho) delta=rho
!
!     recalculate knew and denom if the new f is less than fopt.
!
          if (f < fopt) then
              ksav=knew
              densav=denom
              delsq=delta*delta
              scaden=zero
              biglsq=zero
              knew=0
              do 460 k=1,npt
              hdiag=zero
              do 440 jj=1,nptm
  440         hdiag=hdiag+zmat(k,jj)**2
              den=beta*hdiag+vlag(k)**2
              distsq=zero
              do 450 j=1,n
  450         distsq=distsq+(xpt(k,j)-xnew(j))**2
              temp=max(one,(distsq/delsq)**2)
              if (temp*den > scaden) then
                  scaden=temp*den
                  knew=k
                  denom=den
              end if
  460         biglsq=max(biglsq,temp*vlag(k)**2)
              if (scaden <= half*biglsq) then
                  knew=ksav
                  denom=densav
              end if
          end if
      end if
!
!     update bmat and zmat, so that the knew-th interpolation point can be
!     moved. also update the second derivative terms of the model.
!
      call update (n,npt,bmat,zmat,ndim,vlag,beta,denom,knew,w)
      ih=0
      pqold=pq(knew)
      pq(knew)=zero
      do 470 i=1,n
      temp=pqold*xpt(knew,i)
      do 470 j=1,i
      ih=ih+1
  470 hq(ih)=hq(ih)+temp*xpt(knew,j)
      do 480 jj=1,nptm
      temp=diff*zmat(knew,jj)
      do 480 k=1,npt
  480 pq(k)=pq(k)+temp*zmat(k,jj)
!
!     include the new interpolation point, and make the changes to gopt at
!     the old xopt that are caused by the updating of the quadratic model.
!
      fval(knew)=f
      do 490 i=1,n
      xpt(knew,i)=xnew(i)
  490 w(i)=bmat(knew,i)
      do 520 k=1,npt
      suma=zero
      do 500 jj=1,nptm
  500 suma=suma+zmat(knew,jj)*zmat(k,jj)
      sumb=zero
      do 510 j=1,n
  510 sumb=sumb+xpt(k,j)*xopt(j)
      temp=suma*sumb
      do 520 i=1,n
  520 w(i)=w(i)+temp*xpt(k,i)
      do 530 i=1,n
  530 gopt(i)=gopt(i)+diff*w(i)
!
!     update xopt, gopt and kopt if the new calculated f is less than fopt.
!
      if (f < fopt) then
          kopt=knew
          xoptsq=zero
          ih=0
          do 540 j=1,n
          xopt(j)=xnew(j)
          xoptsq=xoptsq+xopt(j)**2
          do 540 i=1,j
          ih=ih+1
          if (i < j) gopt(j)=gopt(j)+hq(ih)*d(i)
  540     gopt(i)=gopt(i)+hq(ih)*d(j)
          do 560 k=1,npt
          temp=zero
          do 550 j=1,n
  550     temp=temp+xpt(k,j)*d(j)
          temp=pq(k)*temp
          do 560 i=1,n
  560     gopt(i)=gopt(i)+temp*xpt(k,i)
      end if
!
!     calculate the parameters of the least frobenius norm interpolant to
!     the current data, the gradient of this interpolant at xopt being put
!     into vlag(npt+i), i=1,2,...,n.
!
      if (ntrits > 0) then
          do 570 k=1,npt
          vlag(k)=fval(k)-fval(kopt)
  570     w(k)=zero
          do 590 j=1,nptm
          sum=zero
          do 580 k=1,npt
  580     sum=sum+zmat(k,j)*vlag(k)
          do 590 k=1,npt
  590     w(k)=w(k)+sum*zmat(k,j)
          do 610 k=1,npt
          sum=zero
          do 600 j=1,n
  600     sum=sum+xpt(k,j)*xopt(j)
          w(k+npt)=w(k)
  610     w(k)=sum*w(k)
          gqsq=zero
          gisq=zero
          do 630 i=1,n
          sum=zero
          do 620 k=1,npt
  620     sum=sum+bmat(k,i)*vlag(k)+xpt(k,i)*w(k)
          if (xopt(i) == sl(i)) then
              gqsq=gqsq+min(zero,gopt(i))**2
              gisq=gisq+min(zero,sum)**2
          else if (xopt(i) == su(i)) then
              gqsq=gqsq+max(zero,gopt(i))**2
              gisq=gisq+max(zero,sum)**2
          else
              gqsq=gqsq+gopt(i)**2
              gisq=gisq+sum*sum
          end if
  630     vlag(npt+i)=sum
!
!     test whether to replace the new quadratic model by the least frobenius
!     norm interpolant, making the replacement if the test is satisfied.
!
          itest=itest+1
          if (gqsq < ten*gisq) itest=0
          if (itest >= 3) then
              do 640 i=1,max0(npt,nh)
              if (i <= n) gopt(i)=vlag(npt+i)
              if (i <= npt) pq(i)=w(npt+i)
              if (i <= nh) hq(i)=zero
              itest=0
  640         continue
          end if
      end if
!
!     if a trust region step has provided a sufficient decrease in f, then
!     branch for another trust region calculation. the case ntrits=0 occurs
!     when the new interpolation point was reached by an alternative step.
!
      if (ntrits == 0) goto 60
      if (f <= fopt+tenth*vquad) goto 60
!
!     alternatively, find out if the interpolation points are close enough
!       to the best point so far.
!
      distsq=max((two*delta)**2,(ten*rho)**2)
  650 knew=0
      do 670 k=1,npt
      sum=zero
      do 660 j=1,n
  660 sum=sum+(xpt(k,j)-xopt(j))**2
      if (sum > distsq) then
          knew=k
          distsq=sum
      end if
  670 continue
!
!     if knew is positive, then altmov finds alternative new positions for
!     the knew-th interpolation point within distance adelt of xopt. it is
!     reached via label 90. otherwise, there is a branch to label 60 for
!     another trust region iteration, unless the calculations with the
!     current rho are complete.
!
      if (knew > 0) then
          dist=sqrt(distsq)
          if (ntrits == -1) then
              delta=min(tenth*delta,half*dist)
              if (delta <= 1.5d0*rho) delta=rho
          end if
          ntrits=0
          adelt=max(min(tenth*dist,delta),rho)
          dsq=adelt*adelt
          goto 90
      end if
      if (ntrits == -1) goto 680
      if (ratio > zero) goto 60
      if (max(delta,dnorm) > rho) goto 60
!
!     the calculations with the current value of rho are complete. pick the
!       next values of rho and delta.
!
  680 if (rho > rhoend) then
          delta=half*rho
          ratio=rho/rhoend
          if (ratio <= 16.0d0) then
              rho=rhoend
          else if (ratio <= 250.0d0) then
              rho=sqrt(ratio)*rhoend
          else
              rho=tenth*rho
          end if
          delta=max(delta,rho)
          if (iprint >= 2) then
              if (iprint >= 3) print 690
  690         format (5x)
              print 700, rho,nf
  700         format (/4x,'new rho =',1pd11.4,5x,'number of',&
                ' function values =',i6)
              print 710, fval(kopt),(xbase(i)+xopt(i),i=1,n)
  710         format (4x,'least value of f =',1pd23.15,9x,&
                'the corresponding x is:'/(2x,5d15.6))
          end if
          ntrits=0
          nfsav=nf
          goto 60
      end if
!
!     return from the calculation, after another newton-raphson step, if
!       it is too short to have been tried before.
!
      if (ntrits == -1) goto 360
  720 if (fval(kopt) <= fsave) then
          do 730 i=1,n
          x(i)=min(max(xl(i),xbase(i)+xopt(i)),xu(i))
          if (xopt(i) == sl(i)) x(i)=xl(i)
          if (xopt(i) == su(i)) x(i)=xu(i)
  730     continue
          f=fval(kopt)
      end if
      if (iprint >= 1) then
          print 740, nf
  740     format (/4x,'at the return from bobyqa',5x,&
            'number of function values =',i6)
          print 710, f,(x(i),i=1,n)
      end if
      return
      end subroutine bobyqb
!*****************************************************************************************


!*****************************************************************************************
      subroutine altmov (n,npt,xpt,xopt,bmat,zmat,ndim,sl,su,kopt,&
        knew,adelt,xnew,xalt,alpha,cauchy,glag,hcol,w)
      implicit real(wp) (a-h,o-z)
      dimension xpt(npt,*),xopt(*),bmat(ndim,*),zmat(npt,*),sl(*),&
        su(*),xnew(*),xalt(*),glag(*),hcol(*),w(*)
!
!     the arguments n, npt, xpt, xopt, bmat, zmat, ndim, sl and su all have
!       the same meanings as the corresponding arguments of bobyqb.
!     kopt is the index of the optimal interpolation point.
!     knew is the index of the interpolation point that is going to be moved.
!     adelt is the current trust region bound.
!     xnew will be set to a suitable new position for the interpolation point
!       xpt(knew,.). specifically, it satisfies the sl, su and trust region
!       bounds and it should provide a large denominator in the next call of
!       update. the step xnew-xopt from xopt is restricted to moves along the
!       straight lines through xopt and another interpolation point.
!     xalt also provides a large value of the modulus of the knew-th lagrange
!       function subject to the constraints that have been mentioned, its main
!       difference from xnew being that xalt-xopt is a constrained version of
!       the cauchy step within the trust region. an exception is that xalt is
!       not calculated if all components of glag (see below) are zero.
!     alpha will be set to the knew-th diagonal element of the h matrix.
!     cauchy will be set to the square of the knew-th lagrange function at
!       the step xalt-xopt from xopt for the vector xalt that is returned,
!       except that cauchy is set to zero if xalt is not calculated.
!     glag is a working space vector of length n for the gradient of the
!       knew-th lagrange function at xopt.
!     hcol is a working space vector of length npt for the second derivative
!       coefficients of the knew-th lagrange function.
!     w is a working space vector of length 2n that is going to hold the
!       constrained cauchy step from xopt of the lagrange function, followed
!       by the downhill version of xalt when the uphill step is calculated.
!
!     set the first npt components of w to the leading elements of the
!     knew-th column of the h matrix.
!
      half=0.5d0
      one=1.0d0
      zero=0.0d0
      const=one+sqrt(2.0d0)
      do 10 k=1,npt
   10 hcol(k)=zero
      do 20 j=1,npt-n-1
      temp=zmat(knew,j)
      do 20 k=1,npt
   20 hcol(k)=hcol(k)+temp*zmat(k,j)
      alpha=hcol(knew)
      ha=half*alpha
!
!     calculate the gradient of the knew-th lagrange function at xopt.
!
      do 30 i=1,n
   30 glag(i)=bmat(knew,i)
      do 50 k=1,npt
      temp=zero
      do 40 j=1,n
   40 temp=temp+xpt(k,j)*xopt(j)
      temp=hcol(k)*temp
      do 50 i=1,n
   50 glag(i)=glag(i)+temp*xpt(k,i)
!
!     search for a large denominator along the straight lines through xopt
!     and another interpolation point. slbd and subd will be lower and upper
!     bounds on the step along each of these lines in turn. predsq will be
!     set to the square of the predicted denominator for each line. presav
!     will be set to the largest admissible value of predsq that occurs.
!
      presav=zero
      do 80 k=1,npt
      if (k == kopt) goto 80
      dderiv=zero
      distsq=zero
      do 60 i=1,n
      temp=xpt(k,i)-xopt(i)
      dderiv=dderiv+glag(i)*temp
   60 distsq=distsq+temp*temp
      subd=adelt/sqrt(distsq)
      slbd=-subd
      ilbd=0
      iubd=0
      sumin=min(one,subd)
!
!     revise slbd and subd if necessary because of the bounds in sl and su.
!
      do 70 i=1,n
      temp=xpt(k,i)-xopt(i)
      if (temp > zero) then
          if (slbd*temp < sl(i)-xopt(i)) then
              slbd=(sl(i)-xopt(i))/temp
              ilbd=-i
          end if
          if (subd*temp > su(i)-xopt(i)) then
              subd=max(sumin,(su(i)-xopt(i))/temp)
              iubd=i
          end if
      else if (temp < zero) then
          if (slbd*temp > su(i)-xopt(i)) then
              slbd=(su(i)-xopt(i))/temp
              ilbd=i
          end if
          if (subd*temp < sl(i)-xopt(i)) then
              subd=max(sumin,(sl(i)-xopt(i))/temp)
              iubd=-i
          end if
      end if
   70 continue
!
!     seek a large modulus of the knew-th lagrange function when the index
!     of the other interpolation point on the line through xopt is knew.
!
      if (k == knew) then
          diff=dderiv-one
          step=slbd
          vlag=slbd*(dderiv-slbd*diff)
          isbd=ilbd
          temp=subd*(dderiv-subd*diff)
          if (abs(temp) > abs(vlag)) then
              step=subd
              vlag=temp
              isbd=iubd
          end if
          tempd=half*dderiv
          tempa=tempd-diff*slbd
          tempb=tempd-diff*subd
          if (tempa*tempb < zero) then
              temp=tempd*tempd/diff
              if (abs(temp) > abs(vlag)) then
                  step=tempd/diff
                  vlag=temp
                  isbd=0
              end if
          end if
!
!     search along each of the other lines through xopt and another point.
!
      else
          step=slbd
          vlag=slbd*(one-slbd)
          isbd=ilbd
          temp=subd*(one-subd)
          if (abs(temp) > abs(vlag)) then
              step=subd
              vlag=temp
              isbd=iubd
          end if
          if (subd > half) then
              if (abs(vlag) < 0.25d0) then
                  step=half
                  vlag=0.25d0
                  isbd=0
              end if
          end if
          vlag=vlag*dderiv
      end if
!
!     calculate predsq for the current line search and maintain presav.
!
      temp=step*(one-step)*distsq
      predsq=vlag*vlag*(vlag*vlag+ha*temp*temp)
      if (predsq > presav) then
          presav=predsq
          ksav=k
          stpsav=step
          ibdsav=isbd
      end if
   80 continue
!
!     construct xnew in a way that satisfies the bound constraints exactly.
!
      do 90 i=1,n
      temp=xopt(i)+stpsav*(xpt(ksav,i)-xopt(i))
   90 xnew(i)=max(sl(i),min(su(i),temp))
      if (ibdsav < 0) xnew(-ibdsav)=sl(-ibdsav)
      if (ibdsav > 0) xnew(ibdsav)=su(ibdsav)
!
!     prepare for the iterative method that assembles the constrained cauchy
!     step in w. the sum of squares of the fixed components of w is formed in
!     wfixsq, and the free components of w are set to bigstp.
!
      bigstp=adelt+adelt
      iflag=0
  100 wfixsq=zero
      ggfree=zero
      do 110 i=1,n
      w(i)=zero
      tempa=min(xopt(i)-sl(i),glag(i))
      tempb=max(xopt(i)-su(i),glag(i))
      if (tempa > zero .or. tempb < zero) then
          w(i)=bigstp
          ggfree=ggfree+glag(i)**2
      end if
  110 continue
      if (ggfree == zero) then
          cauchy=zero
          goto 200
      end if
!
!     investigate whether more components of w can be fixed.
!
  120 temp=adelt*adelt-wfixsq
      if (temp > zero) then
          wsqsav=wfixsq
          step=sqrt(temp/ggfree)
          ggfree=zero
          do 130 i=1,n
          if (w(i) == bigstp) then
              temp=xopt(i)-step*glag(i)
              if (temp <= sl(i)) then
                  w(i)=sl(i)-xopt(i)
                  wfixsq=wfixsq+w(i)**2
              else if (temp >= su(i)) then
                  w(i)=su(i)-xopt(i)
                  wfixsq=wfixsq+w(i)**2
              else
                  ggfree=ggfree+glag(i)**2
              end if
          end if
  130     continue
          if (wfixsq > wsqsav .and. ggfree > zero) goto 120
      end if
!
!     set the remaining free components of w and all components of xalt,
!     except that w may be scaled later.
!
      gw=zero
      do 140 i=1,n
      if (w(i) == bigstp) then
          w(i)=-step*glag(i)
          xalt(i)=max(sl(i),min(su(i),xopt(i)+w(i)))
      else if (w(i) == zero) then
          xalt(i)=xopt(i)
      else if (glag(i) > zero) then
          xalt(i)=sl(i)
      else
          xalt(i)=su(i)
      end if
  140 gw=gw+glag(i)*w(i)
!
!     set curv to the curvature of the knew-th lagrange function along w.
!     scale w by a factor less than one if that can reduce the modulus of
!     the lagrange function at xopt+w. set cauchy to the final value of
!     the square of this function.
!
      curv=zero
      do 160 k=1,npt
      temp=zero
      do 150 j=1,n
  150 temp=temp+xpt(k,j)*w(j)
  160 curv=curv+hcol(k)*temp*temp
      if (iflag == 1) curv=-curv
      if (curv > -gw .and. curv < -const*gw) then
          scale=-gw/curv
          do 170 i=1,n
          temp=xopt(i)+scale*w(i)
  170     xalt(i)=max(sl(i),min(su(i),temp))
          cauchy=(half*gw*scale)**2
      else
          cauchy=(gw+half*curv)**2
      end if
!
!     if iflag is zero, then xalt is calculated as before after reversing
!     the sign of glag. thus two xalt vectors become available. the one that
!     is chosen is the one that gives the larger value of cauchy.
!
      if (iflag == 0) then
          do 180 i=1,n
          glag(i)=-glag(i)
  180     w(n+i)=xalt(i)
          csave=cauchy
          iflag=1
          goto 100
      end if
      if (csave > cauchy) then
          do 190 i=1,n
  190     xalt(i)=w(n+i)
          cauchy=csave
      end if
  200 return
      end subroutine altmov
!*****************************************************************************************


!*****************************************************************************************
      subroutine prelim (n,npt,x,xl,xu,rhobeg,iprint,maxfun,xbase,&
        xpt,fval,gopt,hq,pq,bmat,zmat,ndim,sl,su,nf,kopt)
      implicit real(wp) (a-h,o-z)
      dimension x(*),xl(*),xu(*),xbase(*),xpt(npt,*),fval(*),gopt(*),&
        hq(*),pq(*),bmat(ndim,*),zmat(npt,*),sl(*),su(*)
!
!     the arguments n, npt, x, xl, xu, rhobeg, iprint and maxfun are the
!       same as the corresponding arguments in subroutine bobyqa.
!     the arguments xbase, xpt, fval, hq, pq, bmat, zmat, ndim, sl and su
!       are the same as the corresponding arguments in bobyqb, the elements
!       of sl and su being set in bobyqa.
!     gopt is usually the gradient of the quadratic model at xopt+xbase, but
!       it is set by prelim to the gradient of the quadratic model at xbase.
!       if xopt is nonzero, bobyqb will change it to its usual value later.
!     nf is maintaned as the number of calls of calfun so far.
!     kopt will be such that the least calculated value of f so far is at
!       the point xpt(kopt,.)+xbase in the space of the variables.
!
!     subroutine prelim sets the elements of xbase, xpt, fval, gopt, hq, pq,
!     bmat and zmat for the first iteration, and it maintains the values of
!     nf and kopt. the vector x is also changed by prelim.
!
!     set some constants.
!
      half=0.5d0
      one=1.0d0
      two=2.0d0
      zero=0.0d0
      rhosq=rhobeg*rhobeg
      recip=one/rhosq
      np=n+1
!
!     set xbase to the initial vector of variables, and set the initial
!     elements of xpt, bmat, hq, pq and zmat to zero.
!
      do 20 j=1,n
      xbase(j)=x(j)
      do 10 k=1,npt
   10 xpt(k,j)=zero
      do 20 i=1,ndim
   20 bmat(i,j)=zero
      do 30 ih=1,(n*np)/2
   30 hq(ih)=zero
      do 40 k=1,npt
      pq(k)=zero
      do 40 j=1,npt-np
   40 zmat(k,j)=zero
!
!     begin the initialization procedure. nf becomes one more than the number
!     of function values so far. the coordinates of the displacement of the
!     next initial interpolation point from xbase are set in xpt(nf+1,.).
!
      nf=0
   50 nfm=nf
      nfx=nf-n
      nf=nf+1
      if (nfm <= 2*n) then
          if (nfm >= 1 .and. nfm <= n) then
              stepa=rhobeg
              if (su(nfm) == zero) stepa=-stepa
              xpt(nf,nfm)=stepa
          else if (nfm > n) then
              stepa=xpt(nf-n,nfx)
              stepb=-rhobeg
              if (sl(nfx) == zero) stepb=min(two*rhobeg,su(nfx))
              if (su(nfx) == zero) stepb=max(-two*rhobeg,sl(nfx))
              xpt(nf,nfx)=stepb
          end if
      else
          itemp=(nfm-np)/n
          jpt=nfm-itemp*n-n
          ipt=jpt+itemp
          if (ipt > n) then
              itemp=jpt
              jpt=ipt-n
              ipt=itemp
          end if
          xpt(nf,ipt)=xpt(ipt+1,ipt)
          xpt(nf,jpt)=xpt(jpt+1,jpt)
      end if
!
!     calculate the next value of f. the least function value so far and
!     its index are required.
!
      do 60 j=1,n
      x(j)=min(max(xl(j),xbase(j)+xpt(nf,j)),xu(j))
      if (xpt(nf,j) == sl(j)) x(j)=xl(j)
      if (xpt(nf,j) == su(j)) x(j)=xu(j)
   60 continue
      call calfun (n,x,f)
      if (iprint == 3) then
          print 70, nf,f,(x(i),i=1,n)
   70      format (/4x,'function number',i6,'    f =',1pd18.10,&
             '    the corresponding x is:'/(2x,5d15.6))
      end if
      fval(nf)=f
      if (nf == 1) then
          fbeg=f
          kopt=1
      else if (f < fval(kopt)) then
          kopt=nf
      end if
!
!     set the nonzero initial elements of bmat and the quadratic model in the
!     cases when nf is at most 2*n+1. if nf exceeds n+1, then the positions
!     of the nf-th and (nf-n)-th interpolation points may be switched, in
!     order that the function value at the first of them contributes to the
!     off-diagonal second derivative terms of the initial quadratic model.
!
      if (nf <= 2*n+1) then
          if (nf >= 2 .and. nf <= n+1) then
              gopt(nfm)=(f-fbeg)/stepa
              if (npt < nf+n) then
                  bmat(1,nfm)=-one/stepa
                  bmat(nf,nfm)=one/stepa
                  bmat(npt+nfm,nfm)=-half*rhosq
              end if
          else if (nf >= n+2) then
              ih=(nfx*(nfx+1))/2
              temp=(f-fbeg)/stepb
              diff=stepb-stepa
              hq(ih)=two*(temp-gopt(nfx))/diff
              gopt(nfx)=(gopt(nfx)*stepb-temp*stepa)/diff
              if (stepa*stepb < zero) then
                  if (f < fval(nf-n)) then
                      fval(nf)=fval(nf-n)
                      fval(nf-n)=f
                      if (kopt == nf) kopt=nf-n
                      xpt(nf-n,nfx)=stepb
                      xpt(nf,nfx)=stepa
                  end if
              end if
              bmat(1,nfx)=-(stepa+stepb)/(stepa*stepb)
              bmat(nf,nfx)=-half/xpt(nf-n,nfx)
              bmat(nf-n,nfx)=-bmat(1,nfx)-bmat(nf,nfx)
              zmat(1,nfx)=sqrt(two)/(stepa*stepb)
              zmat(nf,nfx)=sqrt(half)/rhosq
              zmat(nf-n,nfx)=-zmat(1,nfx)-zmat(nf,nfx)
          end if
!
!     set the off-diagonal second derivatives of the lagrange functions and
!     the initial quadratic model.
!
      else
          ih=(ipt*(ipt-1))/2+jpt
          zmat(1,nfx)=recip
          zmat(nf,nfx)=recip
          zmat(ipt+1,nfx)=-recip
          zmat(jpt+1,nfx)=-recip
          temp=xpt(nf,ipt)*xpt(nf,jpt)
          hq(ih)=(fbeg-fval(ipt+1)-fval(jpt+1)+f)/temp
      end if
      if (nf < npt .and. nf < maxfun) goto 50
      return
      end subroutine prelim
!*****************************************************************************************


!*****************************************************************************************
      subroutine rescue (n,npt,xl,xu,iprint,maxfun,xbase,xpt,&
        fval,xopt,gopt,hq,pq,bmat,zmat,ndim,sl,su,nf,delta,&
        kopt,vlag,ptsaux,ptsid,w)
      implicit real(wp) (a-h,o-z)
      dimension xl(*),xu(*),xbase(*),xpt(npt,*),fval(*),xopt(*),&
        gopt(*),hq(*),pq(*),bmat(ndim,*),zmat(npt,*),sl(*),su(*),&
        vlag(*),ptsaux(2,*),ptsid(*),w(*)
!
!     the arguments n, npt, xl, xu, iprint, maxfun, xbase, xpt, fval, xopt,
!       gopt, hq, pq, bmat, zmat, ndim, sl and su have the same meanings as
!       the corresponding arguments of bobyqb on the entry to rescue.
!     nf is maintained as the number of calls of calfun so far, except that
!       nf is set to -1 if the value of maxfun prevents further progress.
!     kopt is maintained so that fval(kopt) is the least calculated function
!       value. its correct value must be given on entry. it is updated if a
!       new least function value is found, but the corresponding changes to
!       xopt and gopt have to be made later by the calling program.
!     delta is the current trust region radius.
!     vlag is a working space vector that will be used for the values of the
!       provisional lagrange functions at each of the interpolation points.
!       they are part of a product that requires vlag to be of length ndim.
!     ptsaux is also a working space array. for j=1,2,...,n, ptsaux(1,j) and
!       ptsaux(2,j) specify the two positions of provisional interpolation
!       points when a nonzero step is taken along e_j (the j-th coordinate
!       direction) through xbase+xopt, as specified below. usually these
!       steps have length delta, but other lengths are chosen if necessary
!       in order to satisfy the given bounds on the variables.
!     ptsid is also a working space array. it has npt components that denote
!       provisional new positions of the original interpolation points, in
!       case changes are needed to restore the linear independence of the
!       interpolation conditions. the k-th point is a candidate for change
!       if and only if ptsid(k) is nonzero. in this case let p and q be the
!       integer parts of ptsid(k) and (ptsid(k)-p) multiplied by n+1. if p
!       and q are both positive, the step from xbase+xopt to the new k-th
!       interpolation point is ptsaux(1,p)*e_p + ptsaux(1,q)*e_q. otherwise
!       the step is ptsaux(1,p)*e_p or ptsaux(2,q)*e_q in the cases q=0 or
!       p=0, respectively.
!     the first ndim+npt elements of the array w are used for working space.
!     the final elements of bmat and zmat are set in a well-conditioned way
!       to the values that are appropriate for the new interpolation points.
!     the elements of gopt, hq and pq are also revised to the values that are
!       appropriate to the final quadratic model.
!
!     set some constants.
!
      half=0.5d0
      one=1.0d0
      zero=0.0d0
      np=n+1
      sfrac=half/dble(np)
      nptm=npt-np
!
!     shift the interpolation points so that xopt becomes the origin, and set
!     the elements of zmat to zero. the value of sumpq is required in the
!     updating of hq below. the squares of the distances from xopt to the
!     other interpolation points are set at the end of w. increments of winc
!     may be added later to these squares to balance the consideration of
!     the choice of point that is going to become current.
!
      sumpq=zero
      winc=zero
      do 20 k=1,npt
      distsq=zero
      do 10 j=1,n
      xpt(k,j)=xpt(k,j)-xopt(j)
   10 distsq=distsq+xpt(k,j)**2
      sumpq=sumpq+pq(k)
      w(ndim+k)=distsq
      winc=max(winc,distsq)
      do 20 j=1,nptm
   20 zmat(k,j)=zero
!
!     update hq so that hq and pq define the second derivatives of the model
!     after xbase has been shifted to the trust region centre.
!
      ih=0
      do 40 j=1,n
      w(j)=half*sumpq*xopt(j)
      do 30 k=1,npt
   30 w(j)=w(j)+pq(k)*xpt(k,j)
      do 40 i=1,j
      ih=ih+1
   40 hq(ih)=hq(ih)+w(i)*xopt(j)+w(j)*xopt(i)
!
!     shift xbase, sl, su and xopt. set the elements of bmat to zero, and
!     also set the elements of ptsaux.
!
      do 50 j=1,n
      xbase(j)=xbase(j)+xopt(j)
      sl(j)=sl(j)-xopt(j)
      su(j)=su(j)-xopt(j)
      xopt(j)=zero
      ptsaux(1,j)=min(delta,su(j))
      ptsaux(2,j)=max(-delta,sl(j))
      if (ptsaux(1,j)+ptsaux(2,j) < zero) then
          temp=ptsaux(1,j)
          ptsaux(1,j)=ptsaux(2,j)
          ptsaux(2,j)=temp
      end if
      if (abs(ptsaux(2,j)) < half*abs(ptsaux(1,j))) then
          ptsaux(2,j)=half*ptsaux(1,j)
      end if
      do 50 i=1,ndim
   50 bmat(i,j)=zero
      fbase=fval(kopt)
!
!     set the identifiers of the artificial interpolation points that are
!     along a coordinate direction from xopt, and set the corresponding
!     nonzero elements of bmat and zmat.
!
      ptsid(1)=sfrac
      do 60 j=1,n
      jp=j+1
      jpn=jp+n
      ptsid(jp)=dble(j)+sfrac
      if (jpn <= npt) then
          ptsid(jpn)=dble(j)/dble(np)+sfrac
          temp=one/(ptsaux(1,j)-ptsaux(2,j))
          bmat(jp,j)=-temp+one/ptsaux(1,j)
          bmat(jpn,j)=temp+one/ptsaux(2,j)
          bmat(1,j)=-bmat(jp,j)-bmat(jpn,j)
          zmat(1,j)=sqrt(2.0d0)/abs(ptsaux(1,j)*ptsaux(2,j))
          zmat(jp,j)=zmat(1,j)*ptsaux(2,j)*temp
          zmat(jpn,j)=-zmat(1,j)*ptsaux(1,j)*temp
      else
          bmat(1,j)=-one/ptsaux(1,j)
          bmat(jp,j)=one/ptsaux(1,j)
          bmat(j+npt,j)=-half*ptsaux(1,j)**2
      end if
   60 continue
!
!     set any remaining identifiers with their nonzero elements of zmat.
!
      if (npt >= n+np) then
          do 70 k=2*np,npt
          iw=(dble(k-np)-half)/dble(n)
          ip=k-np-iw*n
          iq=ip+iw
          if (iq > n) iq=iq-n
          ptsid(k)=dble(ip)+dble(iq)/dble(np)+sfrac
          temp=one/(ptsaux(1,ip)*ptsaux(1,iq))
          zmat(1,k-np)=temp
          zmat(ip+1,k-np)=-temp
          zmat(iq+1,k-np)=-temp
   70     zmat(k,k-np)=temp
      end if
      nrem=npt
      kold=1
      knew=kopt
!
!     reorder the provisional points in the way that exchanges ptsid(kold)
!     with ptsid(knew).
!
   80 do 90 j=1,n
      temp=bmat(kold,j)
      bmat(kold,j)=bmat(knew,j)
   90 bmat(knew,j)=temp
      do 100 j=1,nptm
      temp=zmat(kold,j)
      zmat(kold,j)=zmat(knew,j)
  100 zmat(knew,j)=temp
      ptsid(kold)=ptsid(knew)
      ptsid(knew)=zero
      w(ndim+knew)=zero
      nrem=nrem-1
      if (knew /= kopt) then
          temp=vlag(kold)
          vlag(kold)=vlag(knew)
          vlag(knew)=temp
!
!     update the bmat and zmat matrices so that the status of the knew-th
!     interpolation point can be changed from provisional to original. the
!     branch to label 350 occurs if all the original points are reinstated.
!     the nonnegative values of w(ndim+k) are required in the search below.
!
          call update (n,npt,bmat,zmat,ndim,vlag,beta,denom,knew,w)
          if (nrem == 0) goto 350
          do 110 k=1,npt
  110     w(ndim+k)=abs(w(ndim+k))
      end if
!
!     pick the index knew of an original interpolation point that has not
!     yet replaced one of the provisional interpolation points, giving
!     attention to the closeness to xopt and to previous tries with knew.
!
  120 dsqmin=zero
      do 130 k=1,npt
      if (w(ndim+k) > zero) then
          if (dsqmin == zero .or. w(ndim+k) < dsqmin) then
              knew=k
              dsqmin=w(ndim+k)
          end if
      end if
  130 continue
      if (dsqmin == zero) goto 260
!
!     form the w-vector of the chosen original interpolation point.
!
      do 140 j=1,n
  140 w(npt+j)=xpt(knew,j)
      do 160 k=1,npt
      sum=zero
      if (k == kopt) then
          continue
      else if (ptsid(k) == zero) then
          do 150 j=1,n
  150     sum=sum+w(npt+j)*xpt(k,j)
      else
          ip=ptsid(k)
          if (ip > 0) sum=w(npt+ip)*ptsaux(1,ip)
          iq=dble(np)*ptsid(k)-dble(ip*np)
          if (iq > 0) then
              iw=1
              if (ip == 0) iw=2
              sum=sum+w(npt+iq)*ptsaux(iw,iq)
          end if
      end if
  160 w(k)=half*sum*sum
!
!     calculate vlag and beta for the required updating of the h matrix if
!     xpt(knew,.) is reinstated in the set of interpolation points.
!
      do 180 k=1,npt
      sum=zero
      do 170 j=1,n
  170 sum=sum+bmat(k,j)*w(npt+j)
  180 vlag(k)=sum
      beta=zero
      do 200 j=1,nptm
      sum=zero
      do 190 k=1,npt
  190 sum=sum+zmat(k,j)*w(k)
      beta=beta-sum*sum
      do 200 k=1,npt
  200 vlag(k)=vlag(k)+sum*zmat(k,j)
      bsum=zero
      distsq=zero
      do 230 j=1,n
      sum=zero
      do 210 k=1,npt
  210 sum=sum+bmat(k,j)*w(k)
      jp=j+npt
      bsum=bsum+sum*w(jp)
      do 220 ip=npt+1,ndim
  220 sum=sum+bmat(ip,j)*w(ip)
      bsum=bsum+sum*w(jp)
      vlag(jp)=sum
  230 distsq=distsq+xpt(knew,j)**2
      beta=half*distsq*distsq+beta-bsum
      vlag(kopt)=vlag(kopt)+one
!
!     kold is set to the index of the provisional interpolation point that is
!     going to be deleted to make way for the knew-th original interpolation
!     point. the choice of kold is governed by the avoidance of a small value
!     of the denominator in the updating calculation of update.
!
      denom=zero
      vlmxsq=zero
      do 250 k=1,npt
      if (ptsid(k) /= zero) then
          hdiag=zero
          do 240 j=1,nptm
  240     hdiag=hdiag+zmat(k,j)**2
          den=beta*hdiag+vlag(k)**2
          if (den > denom) then
              kold=k
              denom=den
          end if
      end if
  250 vlmxsq=max(vlmxsq,vlag(k)**2)
      if (denom <= 1.0d-2*vlmxsq) then
          w(ndim+knew)=-w(ndim+knew)-winc
          goto 120
      end if
      goto 80
!
!     when label 260 is reached, all the final positions of the interpolation
!     points have been chosen although any changes have not been included yet
!     in xpt. also the final bmat and zmat matrices are complete, but, apart
!     from the shift of xbase, the updating of the quadratic model remains to
!     be done. the following cycle through the new interpolation points begins
!     by putting the new point in xpt(kpt,.) and by setting pq(kpt) to zero,
!     except that a return occurs if maxfun prohibits another value of f.
!
  260 do 340 kpt=1,npt
      if (ptsid(kpt) == zero) goto 340
      if (nf >= maxfun) then
          nf=-1
          goto 350
      end if
      ih=0
      do 270 j=1,n
      w(j)=xpt(kpt,j)
      xpt(kpt,j)=zero
      temp=pq(kpt)*w(j)
      do 270 i=1,j
      ih=ih+1
  270 hq(ih)=hq(ih)+temp*w(i)
      pq(kpt)=zero
      ip=ptsid(kpt)
      iq=dble(np)*ptsid(kpt)-dble(ip*np)
      if (ip > 0) then
          xp=ptsaux(1,ip)
          xpt(kpt,ip)=xp
      end if
      if (iq > 0) then
          xq=ptsaux(1,iq)
          if (ip == 0) xq=ptsaux(2,iq)
          xpt(kpt,iq)=xq
      end if
!
!     set vquad to the value of the current model at the new point.
!
      vquad=fbase
      if (ip > 0) then
          ihp=(ip+ip*ip)/2
          vquad=vquad+xp*(gopt(ip)+half*xp*hq(ihp))
      end if
      if (iq > 0) then
          ihq=(iq+iq*iq)/2
          vquad=vquad+xq*(gopt(iq)+half*xq*hq(ihq))
          if (ip > 0) then
              iw=max0(ihp,ihq)-iabs(ip-iq)
              vquad=vquad+xp*xq*hq(iw)
          end if
      end if
      do 280 k=1,npt
      temp=zero
      if (ip > 0) temp=temp+xp*xpt(k,ip)
      if (iq > 0) temp=temp+xq*xpt(k,iq)
  280 vquad=vquad+half*pq(k)*temp*temp
!
!     calculate f at the new interpolation point, and set diff to the factor
!     that is going to multiply the kpt-th lagrange function when the model
!     is updated to provide interpolation to the new function value.
!
      do 290 i=1,n
      w(i)=min(max(xl(i),xbase(i)+xpt(kpt,i)),xu(i))
      if (xpt(kpt,i) == sl(i)) w(i)=xl(i)
      if (xpt(kpt,i) == su(i)) w(i)=xu(i)
  290 continue
      nf=nf+1
      call calfun (n,w,f)
      if (iprint == 3) then
          print 300, nf,f,(w(i),i=1,n)
  300     format (/4x,'function number',i6,'    f =',1pd18.10,&
            '    the corresponding x is:'/(2x,5d15.6))
      end if
      fval(kpt)=f
      if (f < fval(kopt)) kopt=kpt
      diff=f-vquad
!
!     update the quadratic model. the return from the subroutine occurs when
!     all the new interpolation points are included in the model.
!
      do 310 i=1,n
  310 gopt(i)=gopt(i)+diff*bmat(kpt,i)
      do 330 k=1,npt
      sum=zero
      do 320 j=1,nptm
  320 sum=sum+zmat(k,j)*zmat(kpt,j)
      temp=diff*sum
      if (ptsid(k) == zero) then
          pq(k)=pq(k)+temp
      else
          ip=ptsid(k)
          iq=dble(np)*ptsid(k)-dble(ip*np)
          ihq=(iq*iq+iq)/2
          if (ip == 0) then
              hq(ihq)=hq(ihq)+temp*ptsaux(2,iq)**2
          else
              ihp=(ip*ip+ip)/2
              hq(ihp)=hq(ihp)+temp*ptsaux(1,ip)**2
              if (iq > 0) then
                  hq(ihq)=hq(ihq)+temp*ptsaux(1,iq)**2
                  iw=max0(ihp,ihq)-iabs(iq-ip)
                  hq(iw)=hq(iw)+temp*ptsaux(1,ip)*ptsaux(1,iq)
              end if
          end if
      end if
  330 continue
      ptsid(kpt)=zero
  340 continue
  350 return
      end subroutine rescue
!*****************************************************************************************


!*****************************************************************************************
      subroutine trsbox (n,npt,xpt,xopt,gopt,hq,pq,sl,su,delta,&
        xnew,d,gnew,xbdi,s,hs,hred,dsq,crvmin)
      implicit real(wp) (a-h,o-z)
      dimension xpt(npt,*),xopt(*),gopt(*),hq(*),pq(*),sl(*),su(*),&
        xnew(*),d(*),gnew(*),xbdi(*),s(*),hs(*),hred(*)
!
!     the arguments n, npt, xpt, xopt, gopt, hq, pq, sl and su have the same
!       meanings as the corresponding arguments of bobyqb.
!     delta is the trust region radius for the present calculation, which
!       seeks a small value of the quadratic model within distance delta of
!       xopt subject to the bounds on the variables.
!     xnew will be set to a new vector of variables that is approximately
!       the one that minimizes the quadratic model within the trust region
!       subject to the sl and su constraints on the variables. it satisfies
!       as equations the bounds that become active during the calculation.
!     d is the calculated trial step from xopt, generated iteratively from an
!       initial value of zero. thus xnew is xopt+d after the final iteration.
!     gnew holds the gradient of the quadratic model at xopt+d. it is updated
!       when d is updated.
!     xbdi is a working space vector. for i=1,2,...,n, the element xbdi(i) is
!       set to -1.0, 0.0, or 1.0, the value being nonzero if and only if the
!       i-th variable has become fixed at a bound, the bound being sl(i) or
!       su(i) in the case xbdi(i)=-1.0 or xbdi(i)=1.0, respectively. this
!       information is accumulated during the construction of xnew.
!     the arrays s, hs and hred are also used for working space. they hold the
!       current search direction, and the changes in the gradient of q along s
!       and the reduced d, respectively, where the reduced d is the same as d,
!       except that the components of the fixed variables are zero.
!     dsq will be set to the square of the length of xnew-xopt.
!     crvmin is set to zero if d reaches the trust region boundary. otherwise
!       it is set to the least curvature of h that occurs in the conjugate
!       gradient searches that are not restricted by any constraints. the
!       value crvmin=-1.0d0 is set, however, if all of these searches are
!       constrained.
!
!     a version of the truncated conjugate gradient is applied. if a line
!     search is restricted by a constraint, then the procedure is restarted,
!     the values of the variables that are at their bounds being fixed. if
!     the trust region boundary is reached, then further changes may be made
!     to d, each one being in the two dimensional space that is spanned
!     by the current d and the gradient of q at xopt+d, staying on the trust
!     region boundary. termination occurs when the reduction in q seems to
!     be close to the greatest reduction that can be achieved.
!
!     set some constants.
!
      half=0.5d0
      one=1.0d0
      onemin=-1.0d0
      zero=0.0d0
!
!     the sign of gopt(i) gives the sign of the change to the i-th variable
!     that will reduce q from its value at xopt. thus xbdi(i) shows whether
!     or not to fix the i-th variable at one of its bounds initially, with
!     nact being set to the number of fixed variables. d and gnew are also
!     set for the first iteration. delsq is the upper bound on the sum of
!     squares of the free variables. qred is the reduction in q so far.
!
      iterc=0
      nact=0
      sqstp=zero
      do 10 i=1,n
      xbdi(i)=zero
      if (xopt(i) <= sl(i)) then
          if (gopt(i) >= zero) xbdi(i)=onemin
      else if (xopt(i) >= su(i)) then
          if (gopt(i) <= zero) xbdi(i)=one
      end if
      if (xbdi(i) /= zero) nact=nact+1
      d(i)=zero
   10 gnew(i)=gopt(i)
      delsq=delta*delta
      qred=zero
      crvmin=onemin
!
!     set the next search direction of the conjugate gradient method. it is
!     the steepest descent direction initially and when the iterations are
!     restarted because a variable has just been fixed by a bound, and of
!     course the components of the fixed variables are zero. itermax is an
!     upper bound on the indices of the conjugate gradient iterations.
!
   20 beta=zero
   30 stepsq=zero
      do 40 i=1,n
      if (xbdi(i) /= zero) then
          s(i)=zero
      else if (beta == zero) then
          s(i)=-gnew(i)
      else
          s(i)=beta*s(i)-gnew(i)
      end if
   40 stepsq=stepsq+s(i)**2
      if (stepsq == zero) goto 190
      if (beta == zero) then
          gredsq=stepsq
          itermax=iterc+n-nact
      end if
      if (gredsq*delsq <= 1.0d-4*qred*qred) go to 190
!
!     multiply the search direction by the second derivative matrix of q and
!     calculate some scalars for the choice of steplength. then set blen to
!     the length of the the step to the trust region boundary and stplen to
!     the steplength, ignoring the simple bounds.
!
      goto 210
   50 resid=delsq
      ds=zero
      shs=zero
      do 60 i=1,n
      if (xbdi(i) == zero) then
          resid=resid-d(i)**2
          ds=ds+s(i)*d(i)
          shs=shs+s(i)*hs(i)
      end if
   60 continue
      if (resid <= zero) goto 90
      temp=sqrt(stepsq*resid+ds*ds)
      if (ds < zero) then
          blen=(temp-ds)/stepsq
      else
          blen=resid/(temp+ds)
      end if
      stplen=blen
      if (shs > zero) then
          stplen=min(blen,gredsq/shs)
      end if

!
!     reduce stplen if necessary in order to preserve the simple bounds,
!     letting iact be the index of the new constrained variable.
!
      iact=0
      do 70 i=1,n
      if (s(i) /= zero) then
          xsum=xopt(i)+d(i)
          if (s(i) > zero) then
              temp=(su(i)-xsum)/s(i)
          else
              temp=(sl(i)-xsum)/s(i)
          end if
          if (temp < stplen) then
              stplen=temp
              iact=i
          end if
      end if
   70 continue
!
!     update crvmin, gnew and d. set sdec to the decrease that occurs in q.
!
      sdec=zero
      if (stplen > zero) then
          iterc=iterc+1
          temp=shs/stepsq
          if (iact == 0 .and. temp > zero) then
              crvmin=min(crvmin,temp)
              if (crvmin == onemin) crvmin=temp
          end if
          ggsav=gredsq
          gredsq=zero
          do 80 i=1,n
          gnew(i)=gnew(i)+stplen*hs(i)
          if (xbdi(i) == zero) gredsq=gredsq+gnew(i)**2
   80     d(i)=d(i)+stplen*s(i)
          sdec=max(stplen*(ggsav-half*stplen*shs),zero)
          qred=qred+sdec
      end if
!
!     restart the conjugate gradient method if it has hit a new bound.
!
      if (iact > 0) then
          nact=nact+1
          xbdi(iact)=one
          if (s(iact) < zero) xbdi(iact)=onemin
          delsq=delsq-d(iact)**2
          if (delsq <= zero) goto 90
          goto 20
      end if
!
!     if stplen is less than blen, then either apply another conjugate
!     gradient iteration or return.
!
      if (stplen < blen) then
          if (iterc == itermax) goto 190
          if (sdec <= 0.01d0*qred) goto 190
          beta=gredsq/ggsav
          goto 30
      end if
   90 crvmin=zero
!
!     prepare for the alternative iteration by calculating some scalars
!     and by multiplying the reduced d by the second derivative matrix of
!     q, where s holds the reduced d in the call of ggmult.
!
  100 if (nact >= n-1) goto 190
      dredsq=zero
      dredg=zero
      gredsq=zero
      do 110 i=1,n
      if (xbdi(i) == zero) then
          dredsq=dredsq+d(i)**2
          dredg=dredg+d(i)*gnew(i)
          gredsq=gredsq+gnew(i)**2
          s(i)=d(i)
      else
          s(i)=zero
      end if
  110 continue
      itcsav=iterc
      goto 210
!
!     let the search direction s be a linear combination of the reduced d
!     and the reduced g that is orthogonal to the reduced d.
!
  120 iterc=iterc+1
      temp=gredsq*dredsq-dredg*dredg
      if (temp <= 1.0d-4*qred*qred) goto 190
      temp=sqrt(temp)
      do 130 i=1,n
      if (xbdi(i) == zero) then
          s(i)=(dredg*d(i)-dredsq*gnew(i))/temp
      else
          s(i)=zero
      end if
  130 continue
      sredg=-temp
!
!     by considering the simple bounds on the variables, calculate an upper
!     bound on the tangent of half the angle of the alternative iteration,
!     namely angbd, except that, if already a free variable has reached a
!     bound, there is a branch back to label 100 after fixing that variable.
!
      angbd=one
      iact=0
      do 140 i=1,n
      if (xbdi(i) == zero) then
          tempa=xopt(i)+d(i)-sl(i)
          tempb=su(i)-xopt(i)-d(i)
          if (tempa <= zero) then
              nact=nact+1
              xbdi(i)=onemin
              goto 100
          else if (tempb <= zero) then
              nact=nact+1
              xbdi(i)=one
              goto 100
          end if
          ratio=one
          ssq=d(i)**2+s(i)**2
          temp=ssq-(xopt(i)-sl(i))**2
          if (temp > zero) then
              temp=sqrt(temp)-s(i)
              if (angbd*temp > tempa) then
                  angbd=tempa/temp
                  iact=i
                  xsav=onemin
              end if
          end if
          temp=ssq-(su(i)-xopt(i))**2
          if (temp > zero) then
              temp=sqrt(temp)+s(i)
              if (angbd*temp > tempb) then
                  angbd=tempb/temp
                  iact=i
                  xsav=one
              end if
          end if
      end if
  140 continue
!
!     calculate hhd and some curvatures for the alternative iteration.
!
      goto 210
  150 shs=zero
      dhs=zero
      dhd=zero
      do 160 i=1,n
      if (xbdi(i) == zero) then
          shs=shs+s(i)*hs(i)
          dhs=dhs+d(i)*hs(i)
          dhd=dhd+d(i)*hred(i)
      end if
  160 continue
!
!     seek the greatest reduction in q for a range of equally spaced values
!     of angt in [0,angbd], where angt is the tangent of half the angle of
!     the alternative iteration.
!
      redmax=zero
      isav=0
      redsav=zero
      iu=17.0d0*angbd+3.1d0
      do 170 i=1,iu
      angt=angbd*dble(i)/dble(iu)
      sth=(angt+angt)/(one+angt*angt)
      temp=shs+angt*(angt*dhd-dhs-dhs)
      rednew=sth*(angt*dredg-sredg-half*sth*temp)
      if (rednew > redmax) then
          redmax=rednew
          isav=i
          rdprev=redsav
      else if (i == isav+1) then
          rdnext=rednew
      end if
  170 redsav=rednew
!
!     return if the reduction is zero. otherwise, set the sine and cosine
!     of the angle of the alternative iteration, and calculate sdec.
!
      if (isav == 0) goto 190
      if (isav < iu) then
          temp=(rdnext-rdprev)/(redmax+redmax-rdprev-rdnext)
          angt=angbd*(dble(isav)+half*temp)/dble(iu)
      end if
      cth=(one-angt*angt)/(one+angt*angt)
      sth=(angt+angt)/(one+angt*angt)
      temp=shs+angt*(angt*dhd-dhs-dhs)
      sdec=sth*(angt*dredg-sredg-half*sth*temp)
      if (sdec <= zero) goto 190
!
!     update gnew, d and hred. if the angle of the alternative iteration
!     is restricted by a bound on a free variable, that variable is fixed
!     at the bound.
!
      dredg=zero
      gredsq=zero
      do 180 i=1,n
      gnew(i)=gnew(i)+(cth-one)*hred(i)+sth*hs(i)
      if (xbdi(i) == zero) then
          d(i)=cth*d(i)+sth*s(i)
          dredg=dredg+d(i)*gnew(i)
          gredsq=gredsq+gnew(i)**2
      end if
  180 hred(i)=cth*hred(i)+sth*hs(i)
      qred=qred+sdec
      if (iact > 0 .and. isav == iu) then
          nact=nact+1
          xbdi(iact)=xsav
          goto 100
      end if
!
!     if sdec is sufficiently small, then return after setting xnew to
!     xopt+d, giving careful attention to the bounds.
!
      if (sdec > 0.01d0*qred) goto 120
  190 dsq=zero
      do 200 i=1,n
      xnew(i)=max(min(xopt(i)+d(i),su(i)),sl(i))
      if (xbdi(i) == onemin) xnew(i)=sl(i)
      if (xbdi(i) == one) xnew(i)=su(i)
      d(i)=xnew(i)-xopt(i)
  200 dsq=dsq+d(i)**2
      return

!     the following instructions multiply the current s-vector by the second
!     derivative matrix of the quadratic model, putting the product in hs.
!     they are reached from three different parts of the software above and
!     they can be regarded as an external subroutine.
!
  210 ih=0
      do 220 j=1,n
      hs(j)=zero
      do 220 i=1,j
      ih=ih+1
      if (i < j) hs(j)=hs(j)+hq(ih)*s(i)
  220 hs(i)=hs(i)+hq(ih)*s(j)
      do 250 k=1,npt
      if (pq(k) /= zero) then
          temp=zero
          do 230 j=1,n
  230     temp=temp+xpt(k,j)*s(j)
          temp=temp*pq(k)
          do 240 i=1,n
  240     hs(i)=hs(i)+temp*xpt(k,i)
      end if
  250 continue
      if (crvmin /= zero) goto 50
      if (iterc > itcsav) goto 150
      do 260 i=1,n
  260 hred(i)=hs(i)
      goto 120
      end subroutine trsbox
!*****************************************************************************************


!*****************************************************************************************
      subroutine update (n,npt,bmat,zmat,ndim,vlag,beta,denom,&
        knew,w)
      implicit real(wp) (a-h,o-z)
      dimension bmat(ndim,*),zmat(npt,*),vlag(*),w(*)
!
!     the arrays bmat and zmat are updated, as required by the new position
!     of the interpolation point that has the index knew. the vector vlag has
!     n+npt components, set on entry to the first npt and last n components
!     of the product hw in equation (4.11) of the powell (2006) paper on
!     newuoa. further, beta is set on entry to the value of the parameter
!     with that name, and denom is set to the denominator of the updating
!     formula. elements of zmat may be treated as zero if their moduli are
!     at most ztest. the first ndim elements of w are used for working space.
!
!     set some constants.
!
      one=1.0d0
      zero=0.0d0
      nptm=npt-n-1
      ztest=zero
      do 10 k=1,npt
      do 10 j=1,nptm
   10 ztest=max(ztest,abs(zmat(k,j)))
      ztest=1.0d-20*ztest
!
!     apply the rotations that put zeros in the knew-th row of zmat.
!
      jl=1
      do 30 j=2,nptm
      if (abs(zmat(knew,j)) > ztest) then
          temp=sqrt(zmat(knew,1)**2+zmat(knew,j)**2)
          tempa=zmat(knew,1)/temp
          tempb=zmat(knew,j)/temp
          do 20 i=1,npt
          temp=tempa*zmat(i,1)+tempb*zmat(i,j)
          zmat(i,j)=tempa*zmat(i,j)-tempb*zmat(i,1)
   20     zmat(i,1)=temp
      end if
      zmat(knew,j)=zero
   30 continue
!
!     put the first npt components of the knew-th column of hlag into w,
!     and calculate the parameters of the updating formula.
!
      do 40 i=1,npt
      w(i)=zmat(knew,1)*zmat(i,1)
   40 continue
      alpha=w(knew)
      tau=vlag(knew)
      vlag(knew)=vlag(knew)-one
!
!     complete the updating of zmat.
!
      temp=sqrt(denom)
      tempb=zmat(knew,1)/temp
      tempa=tau/temp
      do 50 i=1,npt
   50 zmat(i,1)=tempa*zmat(i,1)-tempb*vlag(i)
!
!     finally, update the matrix bmat.
!
      do 60 j=1,n
      jp=npt+j
      w(jp)=bmat(knew,j)
      tempa=(alpha*vlag(jp)-tau*w(jp))/denom
      tempb=(-beta*w(jp)-tau*vlag(jp))/denom
      do 60 i=1,jp
      bmat(i,j)=bmat(i,j)+tempa*vlag(i)+tempb*w(i)
      if (i > npt) bmat(jp,i-npt)=bmat(i,j)
   60 continue
      return
      end subroutine update
!*****************************************************************************************
 
!*****************************************************************************************
    end module bobyqa_module
!*****************************************************************************************