!*****************************************************************************************
!> author: Jacob Williams
!
!  This module contains the Izzo and Gooding algorithms for solving Lambert's problem.

    module lambert_module

    use kind_module,      only: wp
    use numbers_module
    use vector_module,    only: cross, unit, ucross

    implicit none

    private

    !constants:
    real(wp),parameter :: log2       = log(two)
    real(wp),parameter :: two_third  = two/three
    real(wp),parameter :: four_third = four/three
    real(wp),parameter :: five_half  = five/two
    real(wp),parameter :: three_half = three/two

    abstract interface
        function func(t) result(f)  !! interface to the [[zeroin]] input function
        import :: wp
        implicit none
        real(wp),intent(in)  :: t  !! Independant variable for the function.
        real(wp)             :: f  !! The function evaluated at `t`.
        end function func
    end interface

    !public routines:
    public :: solve_lambert_izzo
    public :: solve_lambert_gooding
    public :: solve_lambert_arorarussell

    public :: lambert_test

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Solve Lambert's problem using Izzo's method.
!
!# References
!
!  1. D. Izzo, [Revisiting Lambert's Problem](http://arxiv-web3.library.cornell.edu/abs/1403.2705)
!     [v2] Tue, 24 Jun 2014 13:08:37 GMT (606kb,D)
!  2. [PyKEP](https://github.com/esa/pykep)
!  3. R. A. Battin, "An Introduction to the Mathematics and Methods of
!     Astrodynamics (Revised Edition)", AIAA Education Series, 1999.

    subroutine solve_lambert_izzo(r1,r2,tof,mu,long_way,multi_revs,v1,v2,status_ok)

    implicit none

    real(wp),dimension(3),intent(in)                :: r1            !! first cartesian position [km]
    real(wp),dimension(3),intent(in)                :: r2            !! second cartesian position [km]
    real(wp),intent(in)                             :: tof           !! time of flight [sec]
    real(wp),intent(in)                             :: mu            !! gravity parameter [km^3/s^2]
    logical,intent(in)                              :: long_way      !! when true, do "long way" (>pi) transfers
    integer,intent(in)                              :: multi_revs    !! maximum number of multi-rev solutions to compute
    real(wp),dimension(:,:),allocatable,intent(out) :: v1            !! vector containing 3d arrays with the cartesian components of the velocities at r1
    real(wp),dimension(:,:),allocatable,intent(out) :: v2            !! vector containing 3d arrays with the cartesian components of the velocities at r2
    logical,intent(out)                             :: status_ok     !! true if everything is OK

    !local variables:
    real(wp),dimension(:),allocatable :: x
    real(wp),dimension(3) :: r1_hat,r2_hat,h_hat,it1,it2,c
    real(wp) :: s,cmag,lambda2,lambda3,lambda5,t,t00,t0,t1,r1mag,r2mag,&
                d3t,d2t,dt,err,t_min,x_old,x_new,term,lambda,&
                gamma,rho,sigma,vr1,vt1,vr2,vt2,y,vt,ly
    integer :: n_solutions,it,m_nmax,i,iter

    !tolerances are from [2]
    integer,parameter   :: max_halley_iters = 12         !! for halley iterations
    real(wp),parameter  :: halley_tol       = 1e-13_wp   !! for halley iterations
    real(wp),parameter  :: htol_singlerev   = 1e-5_wp    !! for householder iterations
    real(wp),parameter  :: htol_multirev    = 1e-8_wp    !! for householder iterations

    !======= Begin Algorithm 1 in [1] =======

    r1mag = norm2(r1)
    r2mag = norm2(r2)

    !check for valid inputs:
    if (tof<=zero .or. mu<=zero .or. r1mag==zero .or. r2mag==zero) then
        write(*,*) 'Error in solve_lambert_izzo: invalid input'
        status_ok = .false.
        return
    end if

    status_ok = .true.

    c       = r2 - r1
    cmag    = norm2(c)
    s       = (cmag + r1mag + r2mag) / two
    t       = sqrt(two*mu/(s*s*s)) * tof
    r1_hat  = unit(r1)
    r2_hat  = unit(r2)
    h_hat   = ucross(r1_hat,r2_hat)
    lambda2 = one - cmag/s
    lambda  = sqrt(lambda2)

    if ( all(h_hat == zero) ) then
        write(*,*) 'Warning: pi transfer in solve_lambert_izzo'
        !arbitrarily choose the transfer plane:
        h_hat = [zero,zero,one]
    end if

    it1 = ucross(h_hat,r1_hat)
    it2 = ucross(h_hat,r2_hat)
    if (long_way) then
        lambda = -lambda
        it1 = -it1
        it2 = -it2
    end if

    lambda3 = lambda*lambda2
    lambda5 = lambda2*lambda3
    t1 = two_third * (one - lambda3)

    !======= Begin Algorithm 2 in [1] =======
    ![xlist, ylist] = findxy(lambda, tof)

    ! maximum number of revolutions for which a solution exists:
    m_nmax = floor(t/pi)

    t00 = acos(lambda) + lambda*sqrt(one-lambda2)
    t0 = t00 + m_nmax*pi

    if (t < t0 .and. m_nmax > 0) then    ! Compute xm and tm using Halley

        dt = zero
        d2t = zero
        d3t = zero
        it = 0
        err = one
        t_min = t0
        x_old = zero
        x_new = zero

        do

            call dtdx(dt,d2t,d3t,x_old,t_min)
            if (dt /= zero) x_new = x_old - dt*d2t/(d2t * d2t - dt * d3t / two)
            err = abs(x_old-x_new)
            if ( (err<halley_tol) .or. (it>max_halley_iters) ) exit
            call compute_tof(x_new,m_nmax,t_min)
            x_old = x_new
            it = it + 1

        end do

        if (t_min > t) m_nmax = m_nmax - 1

    end if
    !======= End Algorithm 2 =======

    !mmax is the maximum number of revolutions.
    !Truncate to user-input multi_revs value if it is larger.
    m_nmax = min(multi_revs,m_nmax)

    !the number of solutions to the problem:
    n_solutions = m_nmax*2 + 1

    !allocate output arrays:
    allocate(v1(3,n_solutions))
    allocate(v2(3,n_solutions))
    allocate(x(n_solutions))

    ! Find the x value for each solution:

    ! initial guess for 0 rev solution:

    if (t>=t00) then

        x(1) = -(t-t00)/(t-t00+four)                        !from [2]

    elseif (t<=t1) then

        x(1) = five_half * (t1*(t1-t))/(t*(one-lambda5)) + one

    else

        x(1) = (t/t00) ** ( log2 / log(t1/t00) ) - one      !from [2]

    end if

    ! 0 rev solution:

    iter = householder(t, x(1), 0, htol_singlerev)

    ! multi-rev solutions:

    do i = 1,m_nmax

        !Eqn 31:

        ! left solution:
        term     = ((i*pi+pi)/(eight*t)) ** two_third
        x(2*i)   = (term-one)/(term+one)
        iter     = householder(t, x(2*i), i, htol_multirev)

        ! right solution:
        term     = ((eight*t)/(i*pi)) ** two_third
        x(2*i+1) = (term-one)/(term+one)
        iter     = householder(t, x(2*i+1), i, htol_multirev)

    end do

    ! construct terminal velocity vectors using each x:

    gamma = sqrt(mu*s/two)
    rho   = (r1mag-r2mag) / cmag
    sigma = sqrt(one-rho*rho)

    do i=1,n_solutions

        y   = sqrt(one - lambda2 + lambda2*x(i)*x(i))
        ly  = lambda*y
        vr1 = gamma*((ly-x(i))-rho*(ly+x(i)))/r1mag
        vr2 = -gamma*((ly-x(i))+rho*(ly+x(i)))/r2mag
        vt  = gamma*sigma*(y+lambda*x(i))
        vt1 = vt/r1mag
        vt2 = vt/r2mag

        v1(:,i) = vr1*r1_hat + vt1*it1   !terminal velocity vectors
        v2(:,i) = vr2*r2_hat + vt2*it2   !

    end do

    deallocate(x)

    contains
!*****************************************************************************************

    !*************************************************************************************
        function householder(t,x,n,eps) result(it)

        !! Householder root solver for x.

        implicit none

        integer                 :: it
        real(wp),intent(in)     :: t
        real(wp),intent(inout)  :: x    !! input is initial guess
        integer,intent(in)      :: n
        real(wp),intent(in)     :: eps

        real(wp) :: xnew,tof,delta,dt,d2t,d3t,dt2,term

        integer,parameter :: max_iters = 15

        do it = 1, max_iters

            call compute_tof(x,n,tof)
            call dtdx(dt,d2t,d3t,x,tof)

            delta = tof-t
            dt2   = dt*dt
            term  = delta*(dt2-delta*d2t/two)/&
                    (dt*(dt2-delta*d2t)+d3t*delta*delta/six)
            xnew  = x - term    ! Ref. [1], p. 12.
            x     = xnew

            if (abs(term)<=eps) exit

        end do

        end function householder
    !*************************************************************************************

    !*************************************************************************************
        subroutine dtdx(dt,d2t,d3t,x,t)

        !! Compute 1st-3rd derivatives for the Householder iterations.

        implicit none

        real(wp),intent(out)  :: dt
        real(wp),intent(out)  :: d2t
        real(wp),intent(out)  :: d3t
        real(wp),intent(in)   :: x
        real(wp),intent(in)   :: t

        real(wp) :: umx2,y,y2,y3,y5,umx2_inv

        umx2     = one-x*x
        umx2_inv = one/umx2
        y        = sqrt(one-lambda2*umx2)    !Ref [1], p. 6
        y2       = y*y
        y3       = y2*y
        y5       = y3*y2

        !Ref [1], Eqn. 22:

        dt       = umx2_inv * (three*t*x - two + two*lambda3*x/y)
        d2t      = umx2_inv * (three*t + five*x*dt + two*(one-lambda2)*lambda3/y3)
        d3t      = umx2_inv * (seven*x*d2t + eight*dt - six*(one-lambda2)*lambda5*x/y5)

        end subroutine dtdx
    !*************************************************************************************

    !*************************************************************************************
        subroutine compute_tof(x,n,tof)

        !!  Compute time of flight from x

        implicit none

        real(wp),intent(in)   :: x
        integer,intent(in)    :: n
        real(wp),intent(out)  :: tof

        real(wp),parameter :: battin   = 0.01_wp
        real(wp),parameter :: lagrange = 0.2_wp

        real(wp) :: dist,k,e,rho,z,eta,s1,q,y,g,d,l,f,a,alpha,beta

        dist = abs(x-one)

        if (dist < lagrange .and. dist > battin) then    !use lagrange tof expression

            !See Ref. [1], Eqn. 9

            a = one / (one-x*x)

            if (a>zero) then   !ellipse

                alpha = two * acos(x)
                beta = two * asin(sqrt(lambda2/a))
                if (lambda<zero) beta = -beta
                tof = ((a*sqrt(a)*((alpha-sin(alpha))-(beta-sin(beta))+two*pi*n))/two)

            else   !hyperbola

                alpha = two * acosh(x)
                beta = two * asinh(sqrt(-lambda2/a))
                if (lambda<zero) beta = -beta
                tof = (-a*sqrt(-a)*((beta-sinh(beta))-(alpha-sinh(alpha)))/two)

            end if

        else

            k    = lambda2
            e    = x*x - one
            rho  = abs(e)
            z    = sqrt(one+k*e)

            if (dist < battin) then  ! use battin series tof expression

                !Equation 20 in [1]:
                !
                !See also: Ref. [3], Eqn. 7.30, p. 304.

                eta = z - lambda*x
                s1  = (one - lambda - x*eta)/two
                q   = four_third*hypergeo(s1)
                tof = (eta*eta*eta*q + four*lambda*eta)/two + n*pi / (rho**three_half)

            else  ! use lancaster tof expresion

                y = sqrt(rho)
                g = x*z - lambda*e
                d = zero
                if (e < zero) then
                    l = acos(g)
                    d = n*pi+l
                else
                    f = y*(z-lambda*x)
                    d = log(f+g)
                end if
                tof = (x-lambda*z-d/y)/e

            end if

        end if

        end subroutine compute_tof
    !*************************************************************************************

    !*************************************************************************************
        pure function hypergeo(x) result(f)

        !!  Evaluate the Gaussian (or ordinary) hypergeometric function: F(3,1,5/2,x)
        !!  See Ref. [3], p. 34.

        implicit none

        real(wp)             :: f
        real(wp),intent(in)  :: x

        real(wp) :: term
        integer  :: i

        real(wp),parameter :: tol = 1e-11_wp
        integer,parameter  :: max_iters = 10000

        !initialize:
        f    = one
        term = one

        !compute the series until the last term is within convergence tolerance:
        do i = 0, max_iters

            term = term*(three+i)*(one+i) / (five_half+i)*x / (i+one)
            f = f + term
            if (abs(term)<=tol) exit

        end do

        end function hypergeo
    !*************************************************************************************

    end subroutine solve_lambert_izzo
!*****************************************************************************************

!*****************************************************************************************
!>
!  Solve Lambert's problem using Gooding's method.
!
!# References
!
!  1. R. H, Gooding. "[A procedure for the solution of Lambert's orbital
!     boundary-value problem](http://adsabs.harvard.edu/abs/1990CeMDA..48..145G)"
!     Celestial Mechanics and Dynamical Astronomy,
!     vol. 48, no. 2, 1990, p. 145-165.
!  2. A. Klumpp, "Performance Comparision of Lambert and Kepler Algorithms",
!     JPL Interoffice Memorandum, 314.1-0426-ARK, Jan 2, 1991.
!     [Zip](http://derastrodynamics.com/docs/lambert_papers_v1.zip)

    subroutine solve_lambert_gooding(r1,r2,tof,mu,long_way,multi_revs,v1,v2,status_ok)

    implicit none

    real(wp),dimension(3),intent(in)                :: r1         !! first cartesian position [km]
    real(wp),dimension(3),intent(in)                :: r2         !! second cartesian position [km]
    real(wp),intent(in)                             :: tof        !! time of flight [sec]
    real(wp),intent(in)                             :: mu         !! gravity parameter [km^3/s^2]
    logical,intent(in)                              :: long_way   !! when true, do "long way" (>pi) transfers
    integer,intent(in)                              :: multi_revs !! maximum number of multi-rev solutions to compute
    real(wp),dimension(:,:),allocatable,intent(out) :: v1         !! vector containing 3d arrays with the cartesian components of the velocities at r1
    real(wp),dimension(:,:),allocatable,intent(out) :: v2         !! vector containing 3d arrays with the cartesian components of the velocities at r2
    logical,intent(out)                             :: status_ok  !! true if everything is OK

    integer                 :: i,j,k,n,n_solutions
    real(wp)                :: num_revs,pa,ta,r1mag,r2mag,dr,r1r2
    real(wp),dimension(3,2) :: vt1,vt2
    real(wp),dimension(3)   :: r1hat,r2hat,r1xr2,rho,r1xr2_hat,etai,etaf
    real(wp),dimension(2)   :: vri,vti,vrf,vtf

    !temp arrays to hold all the solutions:
    ! they will be packed into the output arrays
    logical,dimension(2*multi_revs+1) :: solution_exists
    real(wp),dimension(3,1+2*multi_revs) :: all_vt1, all_vt2

    r1mag = norm2(r1)
    r2mag = norm2(r2)

    if ( r1mag==0.0_wp .or. r2mag==0.0_wp .or. mu<=0.0_wp .or. tof<=0.0_wp ) then
        write(*,*) 'Error in solve_lambert_gooding: invalid input'
        status_ok = .false.
        return
    end if

    !initialize:
    solution_exists = .false.
    status_ok = .true.

    dr       = r1mag - r2mag
    r1r2     = r1mag*r2mag
    r1hat    = r1/r1mag
    r2hat    = r2/r2mag
    r1xr2    = cross(r1,r2)
    if (all(r1xr2==0.0_wp)) then    !the vectors are parallel,
                                    ! so the transfer plane is undefined
        write(*,*) 'Warning: pi transfer in solve_lambert_gooding'
        r1xr2 = [0.0_wp,0.0_wp,1.0_wp]    !degenerate conic...choose the x-y plane
    end if
    r1xr2_hat = unit(r1xr2)

    !a trick to make sure argument is between [-1 and 1]:
    pa = acos(max(-1.0_wp,min(1.0_wp,dot_product(r1hat,r2hat))))

    do i=0,multi_revs

        num_revs = real(i,wp)    !number of complete revs for this case

        !transfer angle and normal vector:
        if (long_way) then ! greater than pi
            ta    =  num_revs * twopi + (twopi - pa)
            rho   = -r1xr2_hat
        else ! less than pi
            ta    = num_revs * twopi + pa
            rho   = r1xr2_hat
        end if

        etai = cross(rho,r1hat)
        etaf = cross(rho,r2hat)

        !Gooding routine:
        call vlamb(mu,r1mag,r2mag,ta,tof,n,vri,vti,vrf,vtf)

        select case (n)    !number of solutions

        case(1)

            vt1(:,1) = vri(1)*r1hat + vti(1)*etai
            vt2(:,1) = vrf(1)*r2hat + vtf(1)*etaf

        case(2)

            vt1(:,1) = vri(1)*r1hat + vti(1)*etai
            vt2(:,1) = vrf(1)*r2hat + vtf(1)*etaf

            vt1(:,2) = vri(2)*r1hat + vti(2)*etai
            vt2(:,2) = vrf(2)*r2hat + vtf(2)*etaf

        end select

        if (i==0 .and. n==1) then    !there can be only one solution
            all_vt1(:,1) = vt1(:,1)
            all_vt2(:,1) = vt2(:,1)
            solution_exists(1) = .true.
        else
            select case(n)
            case(1)
                all_vt1(:,2*i)         = vt1(:,1)
                all_vt2(:,2*i)         = vt2(:,1)
                solution_exists(2*i)   = .true.
            case(2)
                all_vt1(:,2*i)         = vt1(:,1)
                all_vt2(:,2*i)         = vt2(:,1)
                solution_exists(2*i)   = .true.
                all_vt1(:,2*i+1)       = vt1(:,2)
                all_vt2(:,2*i+1)       = vt2(:,2)
                solution_exists(2*i+1) = .true.
            end select
        end if

    end do

    !return all the solutions:
    n_solutions = count(solution_exists)

    allocate(v1(3,n_solutions))
    allocate(v2(3,n_solutions))

    k=0
    do i=1,size(solution_exists)
        if (solution_exists(i)) then
            k=k+1
            v1(:,k) = all_vt1(:,i)
            v2(:,k) = all_vt2(:,i)
        end if
    end do

    contains
!*****************************************************************************************

    !*************************************************************************************
        subroutine vlamb(gm,r1,r2,th,tdelt,n,vri,vti,vrf,vtf)

        !!  Gooding support routine
        !!  Note: this contains the modification from [2]

        implicit none

        real(wp),intent(in) :: gm
        real(wp),intent(in) :: r1
        real(wp),intent(in) :: r2
        real(wp),intent(in) :: th
        real(wp),intent(in) :: tdelt
        integer,intent(out) :: n
        real(wp),dimension(2),intent(out) :: vri
        real(wp),dimension(2),intent(out) :: vti
        real(wp),dimension(2),intent(out) :: vrf
        real(wp),dimension(2),intent(out) :: vtf

        integer :: m,i
        real(wp) :: thr2,r1r2th,csq,c,s,gms,qsqfm1,q,rho,sig,t,x1,x2,x,unused,&
                    qzminx,qzplx,zplqx,vt2,vr1,vt1,vr2

        !the following yields m = 0 when th = 2 pi exactly
        ! neither this nor the original code works for th < 0.0
        thr2 = th
        m = 0
        do while (thr2 > twopi)
            thr2 = thr2 - twopi
            m = m + 1
        end do
        thr2   = thr2 / 2.0_wp

        !note: dr and r1r2 are computed in the calling routine

        r1r2th = 4.0_wp*r1r2*sin(thr2)**2
        csq    = dr*dr + r1r2th
        c      = sqrt(csq)
        s      = (r1 + r2 + c)/2.0_wp
        gms    = sqrt(gm*s/2.0_wp)
        qsqfm1 = c/s
        q      = sqrt(r1r2)*cos(thr2)/s

        if (c/=0.0_wp) then
            rho = dr/c
            sig = r1r2th/csq
        else
            rho = 0.0_wp
            sig = 1.0_wp
        end if

        t = 4.0_wp*gms*tdelt/s**2

        call xlamb(m,q,qsqfm1,t,n,x1,x2)

        !proceed for single solution, or a pair

        do i=1,n

            if (i==1) then
                x = x1
            else
                x = x2
            end if

            call tlamb(m,q,qsqfm1,x,-1,unused,qzminx,qzplx,zplqx)

            vt2 = gms*zplqx*sqrt(sig)
            vr1 = gms*(qzminx - qzplx*rho)/r1
            vt1 = vt2/r1
            vr2 = -gms*(qzminx + qzplx*rho)/r2
            vt2 = vt2/r2

            vri(i) = vr1
            vti(i) = vt1
            vrf(i) = vr2
            vtf(i) = vt2

        end do

        end subroutine vlamb
    !*************************************************************************************

    !*************************************************************************************
        subroutine tlamb(m,q,qsqfm1,x,n,t,dt,d2t,d3t)

        !!  Gooding support routine

          implicit none

          real(wp),intent(in)  :: q
          real(wp),intent(in)  :: qsqfm1
          real(wp),intent(in)  :: x
          integer,intent(in)   :: n
          real(wp),intent(out) :: t
          real(wp),intent(out) :: dt
          real(wp),intent(out) :: d2t
          real(wp),intent(out) :: d3t

          integer   :: m,i
          real(wp)  :: qsq,xsq,u,y,z,&
                        qx,a,b,aa,bb,g,f,fg1,term,fg1sq,twoi1,&
                        told,qz,qz2,u0i,u1i,u2i,u3i,tq,tqsum,&
                        ttmold,p,tterm,tqterm
          logical   :: lm1, l1, l2, l3

          real(wp), parameter :: sw  = 0.4_wp

          lm1 = n==-1
          l1 = n>=1
          l2 = n>=2
          l3 = n==3
          qsq = q*q
          xsq = x*x
          u = (1.0_wp - x)*(1.0_wp + x)
          if (.not.lm1) then
            !(needed if series, and otherwise useful when z = 0)
            dt = 0.0_wp
            d2t = 0.0_wp
            d3t = 0.0_wp
          end if
          if (lm1 .or. m>0 .or. x<0.0_wp .or. abs(u)>sw) then
            !direct computation (not series)
            y = sqrt(abs(u))
            z = sqrt(qsqfm1 + qsq*xsq)
            qx = q*x
            if (qx<=0.0_wp) then
              a = z - qx
              b = q*z - x
            end if
            if (qx<0.0_wp .and. lm1) then
              aa = qsqfm1/a
              bb = qsqfm1*(qsq*u - xsq)/b
            end if
            if (qx==0.0_wp.and.lm1 .or. qx>0.0_wp) then
              aa = z + qx
              bb = q*z + x
            end if
            if (qx>0.0_wp) then
              a = qsqfm1/aa
              b = qsqfm1*(qsq*u - xsq)/bb
            end if
            if (.not.lm1) then
              if (qx*u>=0.0_wp) then
                g = x*z + q*u
               else
                g = (xsq - qsq*u)/(x*z - q*u)
              end if
              f = a*y
              if (x<=1.0_wp) then
                t = m*pi + atan2(f, g)
               else
                if (f>sw) then
                  t = log(f + g)
                 else
                  fg1 = f/(g + 1.0_wp)
                  term = 2.0_wp*fg1
                  fg1sq = fg1*fg1
                  t = term
                  twoi1 = 1.0_wp
                  do
                      twoi1 = twoi1 + 2.0_wp
                      term = term*fg1sq
                      told = t
                      t = t + term/twoi1
                      if (t/=told) cycle
                      exit
                  end do    !(continue looping for inverse tanh)
                end if
              end if
              t = 2.0_wp*(t/y + b)/u
              if (l1 .and. z/=0.0_wp) then
                qz = q/z
                qz2 = qz*qz
                qz = qz*qz2
                dt = (3.0_wp*x*t - 4.0_wp*(a + qx*qsqfm1)/z)/u
                if (l2) d2t = (3.0_wp*t + 5.0_wp*x*dt + 4.0_wp*qz*qsqfm1)/u
                if (l3) d3t = (8.0_wp*dt + 7.0_wp*x*d2t - 12.0_wp*qz*qz2*x*qsqfm1)/u
              end if
             else
              dt = b
              d2t = bb
              d3t = aa
            end if
           else
            !compute by series
            u0i = 1.0_wp
            if (l1) u1i = 1.0_wp
            if (l2) u2i = 1.0_wp
            if (l3) u3i = 1.0_wp
            term = 4.0_wp
            tq = q*qsqfm1
            i = 0
            if (q<0.5_wp) tqsum = 1.0_wp - q*qsq
            if (q>=0.5_wp) tqsum = (1.0_wp/(1.0_wp + q) + q)*qsqfm1
            ttmold = term/3.0_wp
            t = ttmold*tqsum

            do
                i = i + 1
                p = i
                u0i = u0i*u
                if (l1 .and. i>1) u1i = u1i*u
                if (l2 .and. i>2) u2i = u2i*u
                if (l3 .and. i>3) u3i = u3i*u
                term = term*(p - 0.5_wp)/p
                tq = tq*qsq
                tqsum = tqsum + tq
                told = t
                tterm = term/(2.0_wp*p + 3.0_wp)
                tqterm = tterm*tqsum
                t = t - u0i*((1.5_wp*p + 0.25_wp)*tqterm/(p*p - 0.25_wp)-ttmold*tq)
                ttmold = tterm
                tqterm = tqterm*p
                if (l1) dt = dt + tqterm*u1i
                if (l2) d2t = d2t + tqterm*u2i*(p - 1.0_wp)
                if (l3) d3t = d3t + tqterm*u3i*(p - 1.0_wp)*(p - 2.0_wp)
                if (i<n .or. t/=told) cycle
                exit
            end do

            if (l3) d3t = 8.0_wp*x*(1.5_wp*d2t - xsq*d3t)
            if (l2) d2t = 2.0_wp*(2.0_wp*xsq*d2t - dt)
            if (l1) dt = -2.0_wp*x*dt
            t = t/xsq
          end if

        end subroutine tlamb
    !*************************************************************************************

    !*************************************************************************************
        pure function d8rt(x)

        !!  8th root function, used by xlamb

        implicit none

        real(wp) :: d8rt
        real(wp),intent(in) :: x

        d8rt = sqrt(sqrt(sqrt(x)))

        end function d8rt
    !*************************************************************************************

    !*************************************************************************************
        subroutine xlamb(m,q,qsqfm1,tin,n,x,xpl)

        !!  Gooding support routine

        implicit none

        integer,intent(in)   :: m
        real(wp),intent(in)  :: q
        real(wp),intent(in)  :: qsqfm1
        real(wp),intent(in)  :: tin
        integer,intent(out)  :: n
        real(wp),intent(out) :: x
        real(wp),intent(out) :: xpl

        integer  :: i,ij
        real(wp) :: thr2,t0,dt,d2t,d3t,tdiff,w,xm,tmin,&
                    xmold,xtest,tdiffm,d2t2,t,tdiff0

        real(wp),parameter :: tol = 3.0e-7_wp
        real(wp),parameter :: c0  = 1.7_wp
        real(wp),parameter :: c1  = 0.5_wp
        real(wp),parameter :: c2  = 0.03_wp
        real(wp),parameter :: c3  = 0.15_wp
        real(wp),parameter :: c41 = 1.0_wp
        real(wp),parameter :: c42 = 0.24_wp

        thr2 = atan2(qsqfm1, 2.0_wp*q)/pi

        if (m==0) then

            !single-rev starter from t (at x = 0) & bilinear (usually)

            n = 1
            call tlamb(m,q,qsqfm1,0.0_wp,0,t0,dt,d2t,d3t)
            tdiff = tin - t0
            if (tdiff<=0.0_wp) then
                x = t0*tdiff/(-4.0_wp*tin)
                !(-4 is the value of dt, for x = 0)
            else
                x = -tdiff/(tdiff + 4.0_wp)
                w = x + c0*sqrt(2.0_wp*(1.0_wp - thr2))
                if (w<0.0_wp) x = x - sqrt(d8rt(-w))*(x + sqrt(tdiff/(tdiff + 1.5_wp*t0)))
                w = 4.0_wp/(4.0_wp + tdiff)
                x = x*(1.0_wp + x*(c1*w - c2*x*sqrt(w)))
            end if

        else

            !with multirevs, first get t(min) as basis for starter

            xm = 1.0_wp/(1.5_wp*(m + 0.5_wp)*pi)
            if (thr2<0.5_wp) xm = d8rt(2.0_wp*thr2)*xm
            if (thr2>0.5_wp) xm = (2.0_wp - d8rt(2.0_wp - 2.0_wp*thr2))*xm
            !(starter for tmin)

            do i=1,12
                call tlamb(m,q,qsqfm1,xm,3,tmin,dt,d2t,d3t)
                if (d2t==0.0_wp) exit
                xmold = xm
                xm = xm - dt*d2t/(d2t*d2t - dt*d3t/2.0_wp)
                xtest = abs(xmold/xm - 1.0_wp)
                if (xtest<=tol) exit
            end do

            if (i>12) then
                !(break off & exit if tmin not located - should never happen)
                !now proceed from t(min) to full starter
                n = -1
                return
            end if

            tdiffm = tin - tmin
            if (tdiffm<0.0_wp) then

                n = 0
                return
                !(exit if no solution with this m)

            else if (tdiffm==0.0_wp) then

                x = xm
                n = 1
                return
                !(exit if unique solution already from x(tmin))

            else

                n = 3
                if (d2t==0.0_wp) d2t = 6.0_wp*m*pi
                x = sqrt(tdiffm/(d2t/2.0_wp + tdiffm/(1.0_wp - xm)**2))
                w = xm + x
                w = w*4.0_wp/(4.0_wp + tdiffm) + (1.0_wp - w)**2
                x = x*(1.0_wp - (1.0_wp + m + c41*(thr2 - 0.5_wp))/&
                    (1.0_wp + c3*m)*x*(c1*w + c2*x*sqrt(w))) + xm
                d2t2 = d2t/2.0_wp
                if (x>=1.0_wp) then
                    n = 1
                    goto 3
                end if
                !(no finite solution with x > xm)

            end if

        end if

    !(now have a starter, so proceed by halley)

    5    continue

            do i=1,3
                call tlamb(m,q,qsqfm1,x,2,t,dt,d2t,d3t)
                t = tin - t
                if (dt/=0.0_wp) x = x + t*dt/(dt*dt + t*d2t/2.0_wp)
            end do
            if (n/=3) return
            !(exit if only one solution, normally when m = 0)

            n = 2
            xpl = x
            !(second multi-rev starter)
        3   call tlamb(m,q,qsqfm1,0.0_wp,0,t0,dt,d2t,d3t)
            tdiff0 = t0 - tmin
            tdiff = tin - t0
            if (tdiff<=0) then
                x = xm - sqrt(tdiffm/(d2t2 - tdiffm*(d2t2/tdiff0 - 1.0_wp/xm**2)))
            else
                x = -tdiff/(tdiff + 4.0_wp)
                ij = 200
                w = x + c0*sqrt(2.0_wp*(1.0_wp - thr2))
                if (w<0.0_wp) x = x - sqrt(d8rt(-w))*(x + sqrt(tdiff/(tdiff+1.5_wp*t0)))
                w = 4.0_wp/(4.0_wp + tdiff)
                x = x*(1.0_wp + (1.0_wp + m + c42*(thr2 - 0.5_wp))/&
                    (1.0_wp + c3*m)*x*(c1*w - c2*x*sqrt(w)))
                if (x<=-1.0_wp) then
                    n = n - 1
                    !(no finite solution with x < xm)
                    if (n==1) x = xpl
                end if
            end if

        goto 5

        end subroutine xlamb
    !*************************************************************************************

    end subroutine solve_lambert_gooding
!*****************************************************************************************

!*****************************************************************************************
!>
!  Solve Lambert's problem using the Arora/Russell method.
!
!### Reference
! * Arora and Russell, "A fast and robust multiple revolution lambert algorithm
!   using a cosine transformation", AAS Hilton Head 2013, AAS 13-728.
!   [see also the journal article]
! * `Lambert_AroraRussell.cpp` from [EMTG](https://sourceforge.net/projects/emtg/).

    subroutine solve_lambert_arorarussell(r1,r2,tofin,mu,lu,nrev,long_way,&
                                            shortperiod,tolerance,&
                                            max_iterations,v1,v2)

    implicit none

    real(wp),dimension(3),intent(in)  :: r1
    real(wp),dimension(3),intent(in)  :: r2
    real(wp),intent(in)               :: tofin
    real(wp),intent(in)               :: mu
    real(wp),intent(in)               :: lu  !! scale factor in km
    integer,intent(in)                :: nrev
    logical,intent(in)                :: long_way
    logical,intent(in)                :: shortperiod
    real(wp),intent(in)               :: tolerance
    integer,intent(in)                :: max_iterations
    real(wp),dimension(3),intent(out) :: v1
    real(wp),dimension(3),intent(out) :: v2

    !integer,parameter  :: max_iterations_bisection = 50         ! original
    !real(wp),parameter :: bisection_error_tolerance = 1e-8_wp
    integer,parameter  :: max_iterations_bisection = 100         !! modified from original - better match for izzo and gooding
    real(wp),parameter :: bisection_error_tolerance = 1e-14_wp   !! modified from original - better match for izzo and gooding
    real(wp),parameter :: sq2 = sqrt(two)

    logical,parameter :: use_zeroin = .false.   !! jw test

    ! precomputed values for for dele_bi0 for first 20 revs
    ! (deltae_b point where tau crosses zero from arora eqn 55)
    real(wp),dimension(0:19),parameter :: dele_bi0 = [ &
        2.848574_wp, 2.969742_wp, 3.019580_wp, 3.046927_wp,&
        3.064234_wp, 3.076182_wp, 3.084929_wp, 3.091610_wp,&
        3.096880_wp, 3.101145_wp, 3.104666_wp, 3.107623_wp,&
        3.110142_wp, 3.112312_wp, 3.114203_wp, 3.115864_wp,&
        3.117335_wp, 3.118646_wp, 3.119824_wp, 3.120886_wp ]

    real(wp) :: cos_transfer_angle,error
    real(wp) :: k,k_n,k_m,deltak,r1mag,r2mag,tu,r1mag_n,r2mag_n,tof,mu_n
    real(wp) :: theta,s,eps,d,tau,r_buff,t_parabolic,m_k,w_k
    real(wp) :: k_i, z, alpha, f_0, f_1, f_i, f_star
    real(wp) :: tof20,tof100,m_i,w_ki,x_star,w,t0,t1,sgn_tau
    real(wp) :: var2,var1,dele_biguess,var3,sgn_var3,k_biguess
    real(wp) :: m_kbiguess,w_kbiguess,t_biguess,k_bi,t_bi,m_bi,w_bi
    real(wp) :: w_km1,w_kmp5,w_k0,w_kp5,w_k1,tof_km1,tof_kmp5,tof_k0,tof_kp5,tof_k1
    real(wp) :: m_ki,tof_ki,w_km1p41,w_km1p38,w_k1dsq2
    real(wp) :: tof_km1p41,tof_km1p38,tof_k1dsq2
    real(wp) :: c_1,c_2,c_3,c_4,f_n,gamma_1,gamma_2,gamma_3
    real(wp) :: k_initialguess,m
    real(wp) :: dw, ddw
    real(wp) :: v,vv,v3,v4,v5,v6,v7
    real(wp) :: tofc,c,sqrtctau,dtofc,ddtofc
    real(wp) :: low_k,high_k,middle_k,tof_low_k,tof_high_k,tof_middle_k,l_low_k,l_high_k
    real(wp) :: l_middle_k
    real(wp) :: f,g,gdot
    real(wp) :: lu3,luptu
    real(wp) :: r1mag_n_plus_r2mag_n,r1mag_n_times_r2mag_n,omkt
    integer :: sgn_l_middle_k,sgn_l_low_k
    integer :: iterations,iterations_bisection,iflag
    real(wp),dimension(3) :: r1_n,r2_n,v1_n,v2_n

    lu3     = lu*lu*lu
    k       = -one                  ! iteration variable, initialized to -1. it gets overwritten soon but this prevents an annoying compiler warning.
    k_n     = -sq2                  ! lower bound for iteration variable k (will be overwritten based on transfer type)
    k_m     = huge(one)             ! upper bound for iteration variable k (will be overwritten based on transfer type)
    deltak  = huge(one)             ! error in current solution for k
    r1mag   = norm2(r1)             ! magnitude of initial position
    r2mag   = norm2(r2)             ! magnitude of initial position
    tu      = sqrt((one / mu)*lu3)  ! time unit set so that mu = 1
    luptu   = lu / tu

    ! normalize r1, r2, and mu
    r1_n    = r1 / lu
    r2_n    = r2 / lu
    r1mag_n = r1mag / lu
    r2mag_n = r2mag / lu
    tof     = tofin / tu
    mu_n    = mu*(tu*tu) / lu3

    r1mag_n_plus_r2mag_n  = r1mag_n + r2mag_n
    r1mag_n_times_r2mag_n = r1mag_n * r2mag_n

    ! define transfer angle based on long_way flag
    cos_transfer_angle = dot_product(r1_n, r2_n) / r1mag_n_times_r2mag_n !cosine of the transfer angle
    theta = safe_acos(cos_transfer_angle) !transfer angle
    if (long_way) then  !long_way == 1
        theta = twopi - theta
    else
        theta = theta
    end if

    ! calculate s and tau
    s = sqrt(r1mag_n_plus_r2mag_n**3 / mu_n) !semi-perimeter
    iterations = 0 !iterations counter
    eps = 2.0e-2_wp !to prevent singularity at k = sqrt(2)

    !d = theta <= pi ? one : -one
    if (theta<=pi) then
        d = one
    else
        d = -one
    end if

    tau = d * sqrt(r1mag_n_times_r2mag_n*(one + cos_transfer_angle)) / &
                   r1mag_n_plus_r2mag_n !lambert geometry parameter
    r_buff = 0.2_wp  !user-defined parameter to determine when to skip k_bi root solve

    !step 1: generate appropriate initial guess
    !declare some variables that will be used in the initial guess

    !step 1.1 compare the desired time of flight to the parabolic time of flight
    t_parabolic = s * sqrt(one - sq2*tau) * (tau + sq2) / three

    if (tof <= t_parabolic) then !a hyperbolic trajectory is desired

        tof20  = s * sqrt(one - 20.0_wp * tau)  * (tau + 0.04940968903_wp * (one - 20.0_wp * tau))
        tof100 = s * sqrt(one - 100.0_wp * tau) * (tau + 0.00999209404_wp * (one - 100.0_wp * tau))
        if (d > zero) then
            k_n    = sq2
            k_m    = one / tau
            k_i    = (k_n + k_m) / two
            z      = one / sq2
            alpha  = 0.5_wp
            f_0    = t_parabolic
            f_1    = zero
            m_i    = two - k_i * k_i
            w_ki   = compute_w(k_i, m_i, nrev)
            f_i    = s * sqrt(one - k_i*tau) * (tau + (one - k_i*tau) * w_ki)
            f_star = tof
            x_star = pow((z * (f_0 - f_star)*(f_1 - f_i)) / &
                     ((f_i - f_star)*(f_1 - f_0)*z + (f_0 - f_i)*(f_1 - f_star)), &
                     one / alpha)
            k = k_n + (k_m - k_n) * x_star
        else if (tof > tof20) then !arora's "h1" region
            k_n    = sq2
            k_m    = 20.0_wp
            k_i    = (two * k_n + k_m) / three
            z      = one / three
            alpha  = one
            f_0    = t_parabolic
            f_1    = tof20
            m_i    = two - k_i * k_i
            w      = compute_w(k_i, m_i, nrev)
            f_i    = s * sqrt(one - k_i*tau) * (tau + (one - k_i*tau) * w)
            f_star = tof
            x_star = pow((z * (f_0 - f_star)*(f_1 - f_i)) / &
                     ((f_i - f_star)*(f_1 - f_0)*z + (f_0 - f_i)*(f_1 - f_star)), &
                     one / alpha)
            k = k_n + (k_m - k_n) * x_star
        else !arora's "h2" region
            k_n = 20.0_wp
            k_m = huge(1.0_wp)
            t0  = tof20
            t1  = tof100
            k = ((t1*(t0 - tof)*10.0_wp - &
                  t0*sqrt(20.0_wp)*(t1 - tof)) / (tof*(t0 - t1))) * &
                  ((t1*(t0 - tof)*10.0_wp - &
                  t0*sqrt(20.0_wp)*(t1 - tof)) / (tof*(t0 - t1)))
        end if

    else if (nrev >= 1) then !a multi-revolution elliptical orbit is desired

        m_k = two - k*k
        w_k = compute_w(k, m_k, nrev)
        error = tof - compute_tof(k, s, tau, w_k)

        ! calculate estimate for k_bi (k_biguess)
        !sgn_tau = tau >= 0 ? one : -one
        sgn_tau = sign(one,tau)

        var2 = dele_bi0(nrev - 1) ! dummy variable 2
        var1 = eight*abs(tau) / (var2*(sq2 - two*abs(tau))) ! dummy variable 1
        dele_biguess = var2*(one - sgn_tau) + var2*sgn_tau*pow((one / (one + var1)), 0.25_wp)
        var3 = pi - dele_biguess
        !sgn_var3 = var3 >= 0 ? one : -one
        sgn_var3 = sign(one,var3)

        k_biguess = sgn_var3*sqrt(cos(dele_biguess) + one)

        ! calculate t_biguess
        m_kbiguess = two - k_biguess*k_biguess ! m based on k_biguess
        w_kbiguess = compute_w(k_biguess, m_kbiguess, nrev)
        t_biguess  = s*sqrt(one - k_biguess*tau)*(tau + (one - k_biguess*tau)*w_kbiguess)

        ! root solve to find k_bi and t_bi
        if ( (abs(tof - t_biguess) > (r_buff*tof)) .and. (tof > t_biguess)) then
            !do not need to root solve
            k_bi = k_biguess
            t_bi = t_biguess
        else ! find k_bi using newton raphson
            k_bi = compute_kb(k_biguess, tau, s, nrev, &
                              tolerance, max_iterations, sq2, eps)  ! solve via newton raphson
            m_bi = two - k_bi*k_bi
            w_bi = compute_w(k_bi, m_bi, nrev)
            t_bi = s*sqrt(one - k_bi*tau)*(tau + (one - k_bi*tau)*w_bi)
        end if

        if (tof < t_bi) then
            !return - no solution for this nrev
            write(*,*) 'no solution for this nrev'
            return
        end if

        w_km1    = 5.71238898_wp + two * pi*nrev
        w_kmp5   = 1.95494660_wp + 2.71408094_wp*nrev
        w_k0     = sq2 / four*(pi + two * pi*nrev)
        w_kp5    = 0.75913433_wp + 2.71408094_wp*nrev
        w_k1     = 0.57079632_wp + two * pi*nrev
        tof_km1  = compute_tof(-one,    s, tau, w_km1)   ! (k,s,tau,w)
        tof_kmp5 = compute_tof(-0.5_wp, s, tau, w_kmp5)  ! (k,s,tau,w)
        tof_k0   = compute_tof(zero,    s, tau, w_k0)    ! (k,s,tau,w)
        tof_kp5  = compute_tof(0.5_wp,  s, tau, w_kp5)   ! (k,s,tau,w)
        tof_k1   = compute_tof(one,     s, tau, w_k1)    ! (k,s,tau,w)

        ! generate initial guess for k (k_star) for long period and short
        ! period solutions for a given nrev
        if (k_bi >= zero) then

            if (.not. shortperiod) then ! long period solution

                if ((tof >= tof_k1) .and. (k_bi >= one)) then ! use first row of table 5

                    ! compute intial guess for k
                    k_n    = k_bi
                    k_m    = sq2
                    k_i    = (k_bi + sq2)*0.5_wp
                    z      = 0.25_wp
                    alpha  = two
                    f_0    = one / t_bi
                    f_1    = zero
                    m_ki   = two - k_i*k_i
                    w_ki   = compute_w(k_i, m_ki, nrev)
                    tof_ki = compute_tof(k_i, s, tau, w_ki)
                    f_i    = one / tof_ki
                    f_star = one / tof
                    x_star = pow((z * (f_0 - f_star)*(f_1 - f_i)) / &
                                ((f_i - f_star)*(f_1 - f_0)*z + (f_0 - f_i)*(f_1 - f_star)), &
                                one / alpha)
                    k = k_n + (k_m - k_n) * x_star

                else if ((tof >= tof_k1) .and. (k_bi <= one)) then ! use second row of table 5

                    ! compute intial guess for k
                    k_n    = one
                    k_m    = sq2
                    k_i    = (one + two*sq2) / three
                    z      = 4.0_wp / 9.0_wp
                    alpha  = two
                    f_0    = one / tof_k1
                    f_1    = zero
                    m_ki   = two - k_i*k_i
                    w_ki   = compute_w(k_i, m_ki, nrev)
                    tof_ki = compute_tof(k_i, s, tau, w_ki)
                    f_i    = one / tof_ki
                    f_star = one / tof
                    x_star = pow((z * (f_0 - f_star)*(f_1 - f_i)) / &
                                ((f_i - f_star)*(f_1 - f_0)*z + (f_0 - f_i)*(f_1 - f_star)), &
                                one / alpha)
                    k = k_n + (k_m - k_n) * x_star

                else !((tof < tof_k1) .and. (k_bi <= one)) ! use third row of table 5

                    ! compute intial guess for k
                    k_n    = k_bi
                    k_m    = one
                    k_i    = (one + k_bi)*0.5_wp
                    z      = 0.25_wp
                    alpha  = two
                    f_0    = t_bi
                    f_1    = tof_k1
                    m_ki   = two - k_i*k_i
                    w_ki   = compute_w(k_i, m_ki, nrev)
                    tof_ki = compute_tof(k_i, s, tau, w_ki)
                    f_i    = tof_ki
                    f_star = tof
                    x_star = pow((z * (f_0 - f_star)*(f_1 - f_i)) / &
                                ((f_i - f_star)*(f_1 - f_0)*z + &
                                (f_0 - f_i)*(f_1 - f_star)), &
                                one / alpha)
                    k = k_n + (k_m - k_n) * x_star

                end if

            else ! long period solution

                if (tof < tof_k0) then ! use fourth row of table 5

                    ! compute intial guess for k
                    k_n    = zero
                    k_m    = k_bi
                    k_i    = k_bi*0.5_wp
                    z      = 0.5_wp ** ( 6.0_wp / 5.0_wp)
                    alpha  = 6.0_wp / 5.0_wp
                    f_0    = tof_k0
                    f_1    = t_bi
                    m_ki   = two - k_i*k_i
                    w_ki   = compute_w(k_i, m_ki, nrev)
                    tof_ki = compute_tof(k_i, s, tau, w_ki)
                    f_i    = tof_ki
                    f_star = tof
                    x_star = pow((z * (f_0 - f_star)*(f_1 - f_i)) / &
                                ((f_i - f_star)*(f_1 - f_0)*z + (f_0 - f_i)*(f_1 - f_star)), &
                                one / alpha)
                    k = k_n + (k_m - k_n) * x_star

                else if ((tof > tof_k0) .and. (tof < tof_km1)) then ! use fifth row of table 5

                    ! compute intial guess for k
                    k_n    = -one
                    k_m    = zero
                    k_i    = -0.5_wp
                    z      = 0.5_wp
                    alpha  = one
                    f_0    = tof_km1
                    f_1    = tof_k0
                    f_i    = tof_kmp5
                    f_star = tof
                    x_star = pow((z * (f_0 - f_star)*(f_1 - f_i)) / &
                                ((f_i - f_star)*(f_1 - f_0)*z + (f_0 - f_i)*(f_1 - f_star)), &
                                one / alpha)
                    k = k_n + (k_m - k_n) * x_star

                else !(tof > tof_km1) ! use sixth row of table 5

                    ! compute intial guess for k
                    k_n    = -one
                    k_m    = -one*sq2
                    k_i    = (-one - two * sq2) / three
                    z      = 4.0_wp / 9.0_wp
                    alpha  = two
                    f_0    = one / tof_km1
                    f_1    = zero
                    m_ki   = two - k_i*k_i
                    w_ki   = compute_w(k_i, m_ki, nrev)
                    tof_ki = compute_tof(k_i, s, tau, w_ki)
                    f_i    = one / tof_ki
                    f_star = one / tof
                    x_star = pow((z * (f_0 - f_star)*(f_1 - f_i)) / &
                            ((f_i - f_star)*(f_1 - f_0)*z + (f_0 - f_i)*(f_1 - f_star)), &
                            one / alpha)
                    k = k_n + (k_m - k_n) * x_star

                end if

            end if

        else ! k_bi < 0

            if (shortperiod) then ! short period solution

                if ((tof >= tof_km1) .and. (k_bi <= -one)) then ! use first row of table 6

                    ! compute intial guess for k
                    k_n    = k_bi
                    k_m    = -one*sq2
                    k_i    = (k_bi - sq2)*0.5_wp
                    z      = 0.25_wp
                    alpha  = two
                    f_0    = one / t_bi
                    f_1    = zero
                    m_ki   = two - k_i*k_i
                    w_ki   = compute_w(k_i, m_ki, nrev)
                    tof_ki = compute_tof(k_i, s, tau, w_ki)
                    f_i    = one / tof_ki
                    f_star = one / tof
                    x_star = pow((z * (f_0 - f_star)*(f_1 - f_i)) / &
                            ((f_i - f_star)*(f_1 - f_0)*z + (f_0 - f_i)*(f_1 - f_star)), &
                            one / alpha)
                    k = k_n + (k_m - k_n) * x_star

                else if ((tof >= tof_km1) .and. (k_bi >= -one)) then ! use second row of table 6

                    ! compute intial guess for k
                    k_n    = -one
                    k_m    = -sq2
                    k_i    = (-one - two*sq2) / three
                    z      = 4.0_wp / 9.0_wp
                    alpha  = two
                    f_0    = one / tof_km1
                    f_1    = zero
                    m_ki   = two - k_i*k_i
                    w_ki   = compute_w(k_i, m_ki, nrev)
                    tof_ki = compute_tof(k_i, s, tau, w_ki)
                    f_i    = one / tof_ki
                    f_star = one / tof
                    x_star = pow((z * (f_0 - f_star)*(f_1 - f_i)) / &
                            ((f_i - f_star)*(f_1 - f_0)*z + (f_0 - f_i)*(f_1 - f_star)), &
                            one / alpha)
                    k = k_n + (k_m - k_n) * x_star

                else ! ((tof <= tof_km1) .and. (k_bi >= -one)) ! use third row of table 6

                    ! compute intial guess for k
                    k_n    = k_bi
                    k_m    = -one
                    k_i    = (-one + k_bi)*0.5_wp
                    z      = 0.25_wp
                    alpha  = two
                    f_0    = t_bi
                    f_1    = tof_km1
                    m_ki   = two - k_i*k_i
                    w_ki   = compute_w(k_i, m_ki, nrev)
                    tof_ki = compute_tof(k_i, s, tau, w_ki)
                    f_i    = tof_ki
                    f_star = tof
                    x_star = pow((z * (f_0 - f_star)*(f_1 - f_i)) / &
                            ((f_i - f_star)*(f_1 - f_0)*z + (f_0 - f_i)*(f_1 - f_star)), &
                            one / alpha)
                    k = k_n + (k_m - k_n) * x_star

                end if

            else ! long period solution

                if (tof <= tof_k0) then ! use fourth row of table 6

                    ! compute intial guess for k
                    k_n    = k_bi
                    k_m    = zero
                    k_i    = k_bi*0.5_wp
                    z      = 0.25_wp
                    alpha  = two
                    f_0    = t_bi
                    f_1    = tof_k0
                    m_ki   = two - k_i*k_i
                    w_ki   = compute_w(k_i, m_ki, nrev)
                    tof_ki = compute_tof(k_i, s, tau, w_ki)
                    f_i    = tof_ki
                    f_star = tof
                    x_star = pow((z * (f_0 - f_star)*(f_1 - f_i)) / &
                            ((f_i - f_star)*(f_1 - f_0)*z + (f_0 - f_i)*(f_1 - f_star)), &
                            one / alpha)
                    k = k_n + (k_m - k_n) * x_star

                else if ((tof > tof_k0) .and. (tof < tof_k1)) then ! use fifth row of table 6

                    ! compute intial guess for k
                    k_n    = zero
                    k_m    = one
                    k_i    = 0.5_wp
                    z      = 0.5_wp ** (6.0_wp / 5.0_wp)
                    alpha  = 6.0_wp / 5.0_wp
                    f_0    = tof_k0
                    f_1    = tof_k1
                    f_i    = tof_kp5
                    f_star = tof
                    x_star = pow((z * (f_0 - f_star)*(f_1 - f_i)) / &
                            ((f_i - f_star)*(f_1 - f_0)*z + (f_0 - f_i)*(f_1 - f_star)), &
                            one / alpha)
                    k = k_n + (k_m - k_n) * x_star

                else !(tof > tof_k1) ! use sixth row of table 6

                    ! compute intial guess for k
                    k_n    = one
                    k_m    = sq2
                    k_i    = (one + two*sq2) / three
                    z      = 4.0_wp / 9.0_wp
                    alpha  = two
                    f_0    = one / tof_k1
                    f_1    = zero
                    m_ki   = two - k_i*k_i
                    w_ki   = compute_w(k_i, m_ki, nrev)
                    tof_ki = compute_tof(k_i, s, tau, w_ki)
                    f_i    = one / tof_ki
                    f_star = one / tof
                    x_star = pow((z * (f_0 - f_star)*(f_1 - f_i)) / &
                            ((f_i - f_star)*(f_1 - f_0)*z + (f_0 - f_i)*(f_1 - f_star)), &
                            one / alpha)
                    k = k_n + (k_m - k_n) * x_star

                end if

            end if

        end if

    else if (nrev < 1) then ! a single-revolution elliptical orbit is desired

        !determine elliptical region by comparing actual tof (i.e., tstar) with tof(k) function value
        w_km1      = 5.712388981_wp    ! w calculated at k=-1
        w_km1p41   = 4839.684497246_wp ! w calculated at k=-1.41
        w_km1p38   = 212.087279879_wp  ! w calculated at k=-1.38
        w_kmp5     = 1.954946607_wp    ! w calculated at k=-0.5
        w_k1dsq2   = 0.6686397730_wp   ! w calculate at k=a/sqrt(2)
        tof_km1p41 = s*sqrt(one + 1.41_wp*tau)*(tau + (one + 1.41_wp*tau)*w_km1p41)
        tof_km1p38 = s*sqrt(one + 1.38_wp*tau)*(tau + (one + 1.38_wp*tau)*w_km1p38)
        tof_km1    = s*sqrt(one + tau)*(tau + (one + tau)*w_km1)
        tof_k0     = s*(sq2 / 4.0_wp * pi + tau)
        tof_k1dsq2 = s*sqrt(one - one / sq2*tau)*(tau + (one - one / sq2*tau)*w_k1dsq2)
        tof_kmp5   = s*sqrt(one + 0.5_wp*tau)*(tau + (one + 0.5_wp*tau)*w_kmp5)

        if (tof >= tof_km1p38) then ! use region e4 for guess

            k_n     = -1.38_wp
            k_m     = -sq2
            k_i     = -1.41_wp
            c_1     = 49267.0_wp / 27059.0_wp
            c_2     = 67286.0_wp / 17897.0_wp
            c_3     = 2813.0_wp / 287443.0_wp
            c_4     = 4439.0_wp / 3156.0_wp
            alpha   = 243.0_wp
            f_n     = 1.0_wp / tof_km1p38
            f_i     = 1.0_wp / tof_km1p41
            f_star  = 1.0_wp / tof
            gamma_1 = f_i*(f_star - f_n)
            gamma_2 = f_star*(f_n - f_i)
            gamma_3 = f_n*(f_star - f_i)
            k = -c_4*pow(((gamma_1*c_1 - c_3*gamma_3)*c_2 + c_3*c_1*gamma_2) / &
                (gamma_3*c_1 - gamma_1*c_3 - gamma_2*c_2), (one / alpha))

        else if ((tof_km1p38 >= tof) .and. (tof >= tof_km1)) then ! use region e3 for guess

            k_n     = -one
            k_m     = -one*sq2
            k_i     = -1.38_wp
            c_1     = 540649.0_wp / 3125.0_wp
            c_2     = 256.0_wp
            c_3     = one
            c_4     = one
            alpha   = 16.0_wp
            f_n     = 1.0_wp / tof_km1
            f_i     = 1.0_wp / tof_km1p38
            f_star  = 1.0_wp / tof
            gamma_1 = f_i*(f_star - f_n)
            gamma_2 = f_star*(f_n - f_i)
            gamma_3 = f_n*(f_star - f_i)
            k = -c_4*pow(((gamma_1*c_1 - c_3*gamma_3)*c_2 + c_3*c_1*gamma_2) / &
                (gamma_3*c_1 - gamma_1*c_3 - gamma_2*c_2), &
                (one / alpha))

        else if ((tof_km1 >= tof) .and. (tof >= tof_k0)) then ! use region e2 for guess

            k_n    = zero
            k_m    = -one
            k_i    = -0.5_wp
            z      = 0.5_wp
            alpha  = one
            f_0    = tof_k0
            f_1    = tof_km1
            f_i    = tof_kmp5
            f_star = tof
            x_star = (z * (f_0 - f_star)*(f_1 - f_i)) / &
                    ((f_i - f_star)*(f_1 - f_0)*z + &
                    (f_0 - f_i)*(f_1 - f_star))
            k = k_n + (k_m - k_n) * x_star

        else ! use region e1 for guess

            k_n    = zero
            k_m    = sq2
            k_i    = 1.0_wp / sq2
            z      = 0.5_wp
            alpha  = one
            f_0    = tof_k0
            f_1    = t_parabolic
            f_i    = tof_k1dsq2
            f_star = tof
            x_star = (z * (f_0 - f_star)*(f_1 - f_i)) / &
                    ((f_i - f_star)*(f_1 - f_0)*z + &
                    (f_0 - f_i)*(f_1 - f_star))
            k = k_n + (k_m - k_n) * x_star

        end if

    end if

    !step 2: iterate to find k
    k_initialguess = k
    do while (abs(deltak) > tolerance .and. iterations < max_iterations)

        !step 2.1 increment the iterations counter
        iterations = iterations + 1

        !step 2.2 compute w, dw, ddw
        m = two - k*k
        !double sgnk = k >= 0 ? one : -one  ! JW : unnecessary
        w = compute_w(k, m, nrev)

        !dw and ddw
        if ((k < sq2 - eps) .or. (k > sq2 + eps)) then
            dw = (-two + three * w * k) / m
            ddw = (five * dw * k + three * w) / m
        else
            v = k - sq2
            vv = v*v
            v3 = v*vv
            v4 = v3*v
            v5 = v4*v
            v6 = v5*v
            v7 = v6*v
            dw = -one / five &
                + sq2 * v * (four / 35.0_wp) - vv * (6.0_wp / 63.0_wp) &
                + sq2 * v3 * (eight / 231.0_wp) - v4 * (10.0_wp / 429.0_wp) &
                + sq2 * v5 * (48.0_wp / 6435.0_wp) - v6 * (56.0_wp / 12155.0_wp) &
                + sq2 * v7 * (64.0_wp / 46189.0_wp)
            ddw = sq2 * (four / 35.0_wp) - v * (12.0_wp / 63.0_wp) &
                + sq2 * vv * (24.0_wp / 231.0_wp) - v3 * (40.0_wp / 429.0_wp) &
                + sq2 * v4 * (240.0_wp / 6435.0_wp) - v5 * (336.0_wp / 12155.0_wp) &
                + sq2 * v6 * (448.0_wp / 46189.0_wp)
        end if

        omkt = one - k*tau

        !step 2.3 compute tofc, dtofc, ddtofc
        tofc = s * sqrt(omkt) * (tau + (omkt) * w)
        if (abs(tof - tofc) < tolerance) exit

        c = (one - k * tau) / tau
        sqrtctau = sqrt(omkt)
        dtofc = -tofc / (two * c) + s * tau * sqrtctau * (dw * c - w)
        ddtofc = -tofc / (four * c*c) + s * tau * sqrtctau * (w / c + c*ddw - three*dw)

        !step 2.4 compute deltak
        deltak = -(tofc - tof) / (dtofc - (tofc - tof) * ddtofc / (two * dtofc))

        !step 2.5 update k from deltak
        k = k + deltak

        ! step 2.6 bound k
        ! ensure k is not less than -sqrt(2)
        if (k < -sq2) k = -sq2 + 1e-12_wp

        ! ensure k doesn't bleed into to elliptic area when hyperbolic
        if ((tof < t_parabolic) .and. (k < sq2)) k = sq2 + 1e-12_wp

        ! ensure k doesn't bleed into to hyperbolic area when elliptic
        if ((tof > t_parabolic) .and. (k > sq2)) k = sq2 - 1e-12_wp

        ! ensure tof doesn't become indeterminate when d=1
        if ((tof < t_parabolic) .and. (d > zero) .and. &
           ((one - tau*k) < zero)) k = one / tau - 1e-5_wp

    end do

    m_k = two - k*k
    w_k = compute_w(k, m_k, nrev)
    error = tof - compute_tof(k, s, tau, w_k)

    ! step 3: if error is too large from halley's method, try bisection method
    if (abs(error) > 1.0e-4_wp) then

        iterations_bisection = 0

        ! set initial low and high bounds for bisection method based on k_n and k_m
        low_k = k_n + 1.0e-14_wp
        high_k = k_m - 1.0e-14_wp

        ! calculate initial low value of tof(k) and l(k)
        m = two - low_k*low_k
        w = compute_w(low_k, m, nrev)
        tof_low_k = s * sqrt(one - low_k*tau) * (tau + (one - low_k*tau) * w)
        l_low_k = tof_low_k - tof

        ! calculate initial high value of tof(k) and l(k)
        m = two - high_k*high_k
        w = compute_w(high_k, m, nrev)
        tof_high_k = s * sqrt(one - high_k*tau) * (tau + (one - high_k*tau) * w)
        l_high_k = tof_high_k - tof

        if (use_zeroin) then !!!!!!!! test: use ZEROIN for the rootfinder rather than bisection

            call zeroin(l_func,low_k,high_k,l_low_k,l_high_k,&
                        bisection_error_tolerance,k,l_middle_k,iflag)
            error = abs(l_middle_k)

        else !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ORIGINAL bisection

            ! iterate using bisection method until within tolerance or maximum number iterations
            do while (iterations_bisection < max_iterations_bisection)
                middle_k = (low_k + high_k)/two
                m = two - middle_k*middle_k
                w = compute_w(middle_k, m, nrev)
                tof_middle_k = s * sqrt(one - middle_k*tau) * (tau + (one - middle_k*tau) * w)
                l_middle_k = tof_middle_k - tof  !! function for root solve: l(k) = tof(k) - t*

                if (abs(l_middle_k) < bisection_error_tolerance) exit

                !sgn_l_middle_k = l_middle_k >= 0 ? 1 : -1
                !sgn_l_low_k = l_low_k >= 0 ? 1 : -1
                sgn_l_middle_k = sign(one,l_middle_k)
                sgn_l_low_k    = sign(one,l_low_k)
                if (sgn_l_middle_k == sgn_l_low_k) then
                    low_k = middle_k
                    tof_low_k = tof_middle_k
                    l_low_k = l_middle_k
                else
                    high_k = middle_k
                    tof_high_k = tof_middle_k
                    l_high_k = l_middle_k
                end if
                iterations_bisection = iterations_bisection + 1
            end do

            k = middle_k
            error = abs(l_middle_k)
            iterations = iterations + iterations_bisection

        end if

    end if ! end bisection method error catch

    ! step 3: calculate f & g (typos in arora-russell 2013 aas paper, fixed in journal paper)
    omkt = one - k*tau
    f    = one - omkt*r1mag_n_plus_r2mag_n/r1mag_n
    g    = s*tau*sqrt(omkt*mu_n)
    gdot = one - omkt*r1mag_n_plus_r2mag_n/r2mag_n

    ! step 4: calculate v1 and v2
    v1_n = (r2_n - f * r1_n ) / g
    v2_n = (gdot*r2_n - r1_n) / g

    !step 5: transform to input length and time units (final output)
    v1 = v1_n * luptu
    v2 = v2_n * luptu

    contains

    !***************************************************************************
        function l_func(k) result(lk)

        !! function for root solve: `l(k) = tof(k) - t*`
        !! interface to the [[zeroin]] input function

        implicit none

        real(wp),intent(in)  :: k
        real(wp)             :: lk

        real(wp) :: tof_k,omkt

        omkt = one - k*tau

        m     = two - k*k
        w     = compute_w(k, m, nrev)
        tof_k = s * sqrt(omkt) * (tau + omkt * w)
        lk    = tof_k - tof

        end function l_func
    !***************************************************************************

    !***************************************************************************
        pure elemental real(wp) function acoshar(b)

            !! fast computation of acosh

            implicit none

            real(wp),intent(in) :: b

            acoshar = log(b + sqrt(b*b - one))

        end function acoshar
    !***************************************************************************

    !***************************************************************************
        pure function compute_kb(k_bguess,tau,s,nrev,tolerance,max_iterations,sq2,eps) result(k)

            !! newton's method to find k_bi given tau, s, nrev, and initial guess for k_bi

            implicit none

            real(wp),intent(in) :: k_bguess
            real(wp),intent(in) :: tau
            real(wp),intent(in) :: s
            integer,intent(in)  :: nrev
            real(wp),intent(in) :: tolerance
            integer,intent(in)  :: max_iterations  !note: was double in original (typo?)
            real(wp),intent(in) :: sq2
            real(wp),intent(in) :: eps
            real(wp)            :: k

            integer :: iterations
            real(wp) :: deltak,m,w,dw,ddw,tofc,c,sqrtctau,dtofc,ddtofc,omkt !,sgnk

            ! initialize
            iterations = 0
            deltak = 1.0e+10_wp
            k = k_bguess

            ! perform iteration loop
            do while (iterations < max_iterations)

                !step 1.1 increment the iterations counter
                iterations = iterations + 1

                ! step 1.2 bound k
                ! ensure k is not less than -sqrt(2)
                if (k < -sq2) k = -sq2 + 0.00001_wp
                if (k > sq2)  k = sq2 - 0.00001_wp

                !step 1.3 compute w, dw, ddw
                m = two - k*k
                !sgnk = k >= 0 ? one : -one  ! JW : unnecessary

                w   = compute_w(k, m, nrev)
                dw  = (-two + three * w * k) / m
                ddw = (five * dw * k + three * w) / m

                !step 2.3 compute tofc, dtofc, ddtofc
                omkt = one - k*tau
                tofc = s * sqrt(omkt) * (tau + omkt * w)
                c = (one - k * tau) / tau
                sqrtctau = sqrt(omkt)
                dtofc = -tofc / (two * c) + s * tau * sqrtctau * (dw * c - w)

                ! check for convergence
                if (abs(dtofc) < tolerance) exit

                ddtofc = -tofc / (four * c*c) + s * tau * sqrtctau * (w / c + c*ddw - three*dw)

                !step 2.4 compute deltak
                deltak = -one*dtofc / ddtofc

                !step 2.5 update k from deltak
                k = k + deltak

            end do

        end function compute_kb
    !***************************************************************************

    !***************************************************************************
        pure function compute_tof(k,s,tau,w) result(tofc)

            !! calculate tof function for a given k, tau,

            implicit none

            real(wp),intent(in) :: k
            real(wp),intent(in) :: s
            real(wp),intent(in) :: tau
            real(wp),intent(in) :: w
            real(wp)            :: tofc

            real(wp) :: omkt

            omkt = one - k*tau
            tofc = s*sqrt(omkt)*(tau + omkt*w)

        end function compute_tof
    !***************************************************************************

    !***************************************************************************
        pure elemental real(wp) function acosar(x)

            !! fast, rough computation of acos()

            implicit none

            real(wp),intent(in) :: x

            real(wp) :: coeff,fx,sgnx

            fx = abs(x)
            !sgnx = x >= zero ? one : -one
            if (x >= zero) then
                sgnx = one
            else
                sgnx = -one
            end if

            if (fx <= 0.6_wp) then
                coeff = (0.000014773722_wp + (1.1782782_wp - 0.52020038_wp * fx) * fx) / &
                        (1.1793469_wp + (-0.53277664_wp - 0.14454764_wp * fx) * fx)
            else if (fx <= 0.97_wp) then
                coeff = (0.011101554_wp + (8.9810074_wp + (-14.816468_wp + 5.9249913_wp * fx) * fx) * fx) / &
                        (9.2299851_wp + (-16.001036_wp + 6.8381053_wp * fx) * fx)
            else if (fx <= 0.99_wp) then
                coeff = (-35.750586_wp + (107.24325_wp - 70.780244_wp * fx) * fx) / &
                        (27.105764_wp - 26.638535_wp * fx)
            else
                coeff = safe_asin(fx)
            end if

            acosar = pi/two - sgnx * coeff

        end function acosar
    !***************************************************************************

    !***************************************************************************
        pure function compute_w(k,m,nrev) result(w)
            !! function to compute the parameter w

            implicit none

            real(wp),intent(in) :: k
            real(wp),intent(in) :: m
            integer,intent(in)  :: nrev
            real(wp)            :: w

            real(wp),parameter :: sq2 = sqrt(two)
            real(wp),parameter :: eps = 2.0e-2_wp

            real(wp) :: sgnk,v,v2,v3,v4,v5,v6,v7,v8

            !sgnk = k < zero ? -one : one
            sgnk = sign(one,k)

            if (-sq2 <= k .and. k < sq2 - eps) then !elliptical orbit case
                w = ((one - sgnk) * pi + sgnk*acos(one - m) + two * pi*nrev) / sqrt(m*m*m) - k / m
            else if (k < sq2 .and. nrev > 0) then
                w = ((one - sgnk) * pi + sgnk*acos(one - m) + two * pi*nrev) / sqrt(m*m*m) - k / m
            else if (k > sq2 + eps) then !hyperbolic orbits
                w = -one*acoshar(one - m) / sqrt(-m*m*m) - k / m
            else if (sq2 - eps <= k .and. k <= sq2 + eps) then !nrev = 0 case
                v = k - sq2
                v2 = v*v
                v3 = v*v2
                v4 = v3*v
                v5 = v4*v
                v6 = v5*v
                v7 = v6*v
                v8 = v7*v
                w = sq2 / three - v / five &
                    + sq2 * v2 * (two / 35.0_wp) - v3 * (two / 63.0_wp) &
                    + sq2 * v4 * (two / 231.0_wp) - v5 * (two / 429.0_wp) &
                    + sq2 * v6 * (eight / 6435.0_wp) - v7 * (eight / 12155.0_wp) &
                    + sq2 * v8 * (eight / 46189.0_wp)
            else
                !write(*,*) 'error on w compute *************************', k    !!! this needs to set status_ok = .false.
                w = huge(1.0_wp)
            end if

        end function compute_w
    !***************************************************************************

    !***************************************************************************
        pure elemental real(wp) function safe_acos(x)

            !! return x > one ? zero : (x < -one ? pi : acos(x))

            implicit none

            real(wp),intent(in) :: x

            if (x>one) then
                safe_acos = zero
            else
                if (x<-one) then
                    safe_acos = pi
                else
                    safe_acos = acos(x)
                end if
            end if

        end function safe_acos
    !***************************************************************************

    !***************************************************************************
        pure elemental real(wp) function safe_asin(x)

            !! safe_asin = x > one ? -piover2 : (x < -one ? piover2 : asin(x))

            implicit none

            real(wp),intent(in) :: x

            real(wp),parameter :: piover2 = pi/two

            if (x > one) then
                safe_asin = -piover2
            else
                if (x < -one) then
                    safe_asin = piover2
                else
                    safe_asin = asin(x)
                end if
            end if

        end function safe_asin
    !***************************************************************************

    !***************************************************************************
        pure elemental real(wp) function pow(x,y)

            !! return x**y

            implicit none

            real(wp),intent(in) :: x
            real(wp),intent(in) :: y

            pow = x ** y

        end function pow
    !***************************************************************************

    !*****************************************************************************************
    !>
    !  Find a zero of the function \( f(x) \) in the given interval
    !  \( [a_x,b_x] \) to within a tolerance \( 4 \epsilon |x| + tol \),
    !  where \( \epsilon \) is the relative machine precision defined as
    !  the smallest representable number such that \( 1.0 + \epsilon > 1.0 \).
    !
    !  It is assumed that \( f(a_x) \) and \( f(b_x) \) have opposite signs.
    !
    !#References
    !  * R. P. Brent, "[An algorithm with guaranteed convergence for
    !    finding a zero of a function](http://maths-people.anu.edu.au/~brent/pd/rpb005.pdf)",
    !    The Computer Journal, Vol 14, No. 4., 1971.
    !  * R. P. Brent, "[Algorithms for minimization without derivatives](http://maths-people.anu.edu.au/~brent/pub/pub011.html)",
    !    Prentice-Hall, Inc., 1973.
    !
    !# See also
    !  1. [zeroin.f](http://www.netlib.org/go/zeroin.f) from Netlib

        subroutine zeroin(f,ax,bx,fax,fbx,tol,xzero,fzero,iflag)

        use iso_fortran_env, only: error_unit

        implicit none

        procedure(func)       :: f       !! the function to find the root of
        real(wp),intent(in)   :: ax      !! left endpoint of initial interval
        real(wp),intent(in)   :: bx      !! right endpoint of initial interval
        real(wp),intent(in)   :: fax     !! `f(ax)` for endpoint
        real(wp),intent(in)   :: fbx     !! `f(ax)` for endpoint
        real(wp),intent(in)   :: tol     !! desired length of the interval of
                                         !! uncertainty of the final result (>=0)
        real(wp),intent(out)  :: xzero   !! abscissa approximating a zero of `f`
                                         !! in the interval `ax`,`bx`
        real(wp),intent(out)  :: fzero   !! value of `f` at the root (`f(xzero)`)
        integer,intent(out)   :: iflag   !! status flag (`-1`=error, `0`=root found)

        real(wp) :: a,b,c,d,e,fa,fb,fc,tol1,xm,p,q,r,s

        real(wp),parameter :: eps = epsilon(one)

        tol1 = eps+one

        a=ax
        b=bx
        fa = fax
        fb = fbx

        !check trivial cases first:
        if (fa==zero) then

            iflag = 0
            xzero = a
            fzero = fa

        elseif (fb==zero) then

            iflag = 0
            xzero = b
            fzero = fb

        elseif (fa*(fb/abs(fb))<zero) then  ! check that f(ax) and f(bx) have different signs

            c=a
            fc=fa
            d=b-a
            e=d

            do

                if (abs(fc)<abs(fb)) then
                    a=b
                    b=c
                    c=a
                    fa=fb
                    fb=fc
                    fc=fa
                end if

                tol1=two*eps*abs(b)+0.5_wp*tol
                xm = 0.5_wp*(c-b)
                if ((abs(xm)<=tol1).or.(fb==zero)) exit

                ! see if a bisection is forced
                if ((abs(e)>=tol1).and.(abs(fa)>abs(fb))) then
                    s=fb/fa
                    if (a/=c) then
                        ! inverse quadratic interpolation
                        q=fa/fc
                        r=fb/fc
                        p=s*(two*xm*q*(q-r)-(b-a)*(r-one))
                        q=(q-one)*(r-one)*(s-one)
                    else
                        ! linear interpolation
                        p=two*xm*s
                        q=one-s
                    end if
                    if (p<=zero) then
                        p=-p
                    else
                        q=-q
                    end if
                    s=e
                    e=d
                    if (((two*p)>=(three*xm*q-abs(tol1*q))) .or. &
                        (p>=abs(0.5_wp*s*q))) then
                        d=xm
                        e=d
                    else
                        d=p/q
                    end if
                else
                    d=xm
                    e=d
                end if

                a=b
                fa=fb
                if (abs(d)<=tol1) then
                    if (xm<=zero) then
                        b=b-tol1
                    else
                        b=b+tol1
                    end if
                else
                    b=b+d
                end if
                fb=f(b)
                if ((fb*(fc/abs(fc)))>zero) then
                    c=a
                    fc=fa
                    d=b-a
                    e=d
                end if

            end do

            iflag = 0
            xzero = b
            fzero = fb

        else

            iflag = -1
            write(error_unit,'(A)')&
                'Error in zeroin: f(ax) and f(bx) do not have different signs.'

        end if

        end subroutine zeroin
    !*****************************************************************************************

    end subroutine solve_lambert_arorarussell
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!
!  Compare the Lambert routines.

    subroutine lambert_test()

    use gooding_module,    only: pv3els
    use random_module,     only: get_random_number

    implicit none

    real(wp),parameter :: tol = 1.0e-11_wp

    real(wp),dimension(:,:),allocatable :: izzo_v1,izzo_v2
    real(wp),dimension(:,:),allocatable :: gooding_v1,gooding_v2
    real(wp),dimension(6)   :: rv1,rv2,pv,e
    integer                 :: imeth,icases,i,multi_revs,iway,n_cases
    real(wp)                :: tof,err1,err2
    logical                 :: long_way,status_ok
    character(len=100)      :: str
    real                    :: tstart, tend
    real(wp),dimension(3)   :: arora_v1,arora_v2

    real(wp),dimension(6),parameter :: rv1_base = &
                                        [1e6_wp,1e6_wp,1e6_wp,10.0_wp,10.0_wp,10.0_wp]
    real(wp),dimension(6),parameter :: rv2_base = &
                                        [1e6_wp,1e6_wp,1e6_wp,10.0_wp,10.0_wp,10.0_wp]
    real(wp),parameter :: tof_base = 86400.0_wp  !sec
    real(wp),parameter :: mu = 3.986004362330e+05_wp

    !real(wp),parameter :: auora_tol      = 1.0e-14_wp   !  122946 cases/sec
    real(wp),parameter :: auora_tol       = 1.0e-9_wp    !! 203025 cases/sec
    integer,parameter  :: max_iterations = 100
    logical,parameter  :: shortperiod    = .true.   !! "solution 1" for mult-rev case.
    real(wp),parameter :: lu             = 1.0e5_wp

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' lambert_test'
    write(*,*) '---------------'
    write(*,*) ''

    write(*,*) ''
    write(*,*) '   Test 1'
    write(*,*) ''

    n_cases = 10
    multi_revs = 0  !1

    !reseed the random number generator:
    call random_seed()

    do icases=1,n_cases

        write(*,*) ''
        do i=1,6
            rv1(i) = get_random_number(-rv1_base(i),rv1_base(i))
        end do
        do i=1,6
            rv2(i) = get_random_number(-rv2_base(i),rv2_base(i))
        end do
        tof = get_random_number(1000.0_wp, tof_base)

        do iway=1,2    !short and long way

            long_way = (iway==1)

            do imeth=1,3    ![izzo, gooding, auora]

                !if (icases==1 .and. imeth==1 .and. iway==1) &
                !        write(*,'(*(A30,1X))') &
                !        'case',&
                !        'x1','y1','z1','vx1','vy1','vz1',&
                !        'x2','y2','z2','vx2','vy2','vz2','tof'
                !if (imeth==1) write(*,'(I30,1X,*(F30.6,1X))') icases, rv1, rv2, tof

                select case (imeth)
                case(1)
                    call solve_lambert_izzo(   rv1(1:3),rv2(1:3),tof,mu,long_way,&
                                                multi_revs,izzo_v1,izzo_v2,&
                                                status_ok)
                case(2)
                    call solve_lambert_gooding(rv1(1:3),rv2(1:3),tof,mu,long_way,&
                                                multi_revs,gooding_v1,gooding_v2,&
                                                status_ok)
                case(3)
                    call solve_lambert_arorarussell(rv1(1:3),rv2(1:3),tof,mu,lu,multi_revs,long_way,&
                                                            shortperiod,auora_tol,&
                                                            max_iterations,arora_v1,arora_v2)
                end select

            end do

            !results:
            if (size(izzo_v1,2)==size(gooding_v1,2)) then

                do i = 1, size(izzo_v1,2)

                    !orb. elements of transfer orbit:
                    pv = [rv1(1:3),gooding_v1(:,i)]
                    call pv3els(mu, pv, e)

                    err1 = norm2( izzo_v1(:,i) - gooding_v1(:,i) )
                    err2 = norm2( izzo_v2(:,i) - gooding_v2(:,i) )

                    if (err1>tol) then
                        str = '*****IZZO-GOODING ERROR*****'
                    else
                        str = '     IZZO-GOODING v1'
                    end if
                    write(*,'(I5,1X,E25.16,1X,E25.16,1X,E25.16,1X,A)')&
                                 icases, e(1:2), err1, trim(str)

                    if (err2>tol) then
                        str = '*****IZZO-GOODING ERROR*****'
                    else
                        str = '     IZZO-GOODING v2'
                    end if
                    write(*,'(I5,1X,E25.16,1X,E25.16,1X,E25.16,1X,A)')&
                                 icases, e(1:2), err2, trim(str)

                     err1 = norm2( arora_v1 - gooding_v1(:,i) )
                     err2 = norm2( arora_v2 - gooding_v2(:,i) )

                     if (err1>tol) then
                         str = '*****ARORA-GOODING ERROR*****'
                     else
                         str = '     ARORA-GOODING v1'
                     end if
                     write(*,'(I5,1X,E25.16,1X,E25.16,1X,E25.16,1X,A)')&
                                  icases, e(1:2), err1, trim(str)

                     if (err2>tol) then
                         str = '*****ARORA-GOODING ERROR*****'
                     else
                         str = '     ARORA-GOODING v2'
                     end if
                     write(*,'(I5,1X,E25.16,1X,E25.16,1X,E25.16,1X,A)')&
                                  icases, e(1:2), err2, trim(str)

                end do

            else
                write(*,*) icases, 'Error: arrays sizes are not the same'
                stop
            end if

        end do

    end do

    write(*,*) ''
    write(*,*) '   Test 2: Speed test'
    write(*,*) ''

    n_cases = 1000000

    do imeth=1,3    ![izzo, gooding, auora]

        !reseed the random number generator:
        call random_seed()

        call cpu_time(tstart)

        do icases=1,n_cases

            do i=1,6
                rv1(i) = get_random_number(-rv1_base(i),rv1_base(i))
            end do
            do i=1,6
                rv2(i) = get_random_number(-rv2_base(i),rv2_base(i))
            end do
            tof = get_random_number(1000.0_wp, tof_base)

            do iway=1,2    !short and long way

                long_way = (iway==1)

                select case (imeth)
                case(1)
                    call solve_lambert_izzo(   rv1(1:3),rv2(1:3),tof,mu,long_way,&
                                                multi_revs,izzo_v1,izzo_v2,&
                                                status_ok)
                case(2)
                    call solve_lambert_gooding(rv1(1:3),rv2(1:3),tof,mu,long_way,&
                                                multi_revs,gooding_v1,gooding_v2,&
                                                status_ok)
                case(3)
                    call solve_lambert_arorarussell(rv1(1:3),rv2(1:3),tof,mu,lu,multi_revs,long_way,&
                                                            shortperiod,auora_tol,&
                                                            max_iterations,arora_v1,arora_v2)
                end select

            end do

        end do

        call cpu_time(tend)
        select case (imeth)
        case(1)
            write(*,*) 'run time for Izzo   : ', tend-tstart, ' sec.  ', dble(n_cases)/(tend-tstart),' cases/sec'
        case(2)
            write(*,*) 'run time for Gooding: ', tend-tstart, ' sec.  ', dble(n_cases)/(tend-tstart),' cases/sec'
        case(3)
            write(*,*) 'run time for Arora  : ', tend-tstart, ' sec.  ', dble(n_cases)/(tend-tstart),' cases/sec'
        end select
        write(*,*) ''

    end do

    end subroutine lambert_test
!*****************************************************************************************

    end module lambert_module
!*****************************************************************************************
