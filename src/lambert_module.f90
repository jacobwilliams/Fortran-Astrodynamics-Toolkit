!*****************************************************************************************
    module lambert_module
!*****************************************************************************************
!****h* FAT/lambert_module
!
!  NAME
!    lambert_module
!
!  DESCRIPTION
!    This module contains the Izzo and Gooding algorithms for solving Lambert's problem.
!
!  AUTHOR
!    Jacob Williams
!
!*****************************************************************************************    
    
    use kind_module,      only: wp
    use vector_module,    only: cross, unit, ucross

    implicit none

    private

    !constants:
    real(wp),parameter :: zero       = 0.0_wp
    real(wp),parameter :: one        = 1.0_wp
    real(wp),parameter :: two        = 2.0_wp
    real(wp),parameter :: three      = 3.0_wp
    real(wp),parameter :: four       = 4.0_wp
    real(wp),parameter :: five       = 5.0_wp
    real(wp),parameter :: six        = 6.0_wp
    real(wp),parameter :: seven      = 7.0_wp
    real(wp),parameter :: eight      = 8.0_wp
    real(wp),parameter :: pi         = acos(-one)
    real(wp),parameter :: twopi      = two*pi
    real(wp),parameter :: log2       = log(two)
    real(wp),parameter :: two_third  = two/three
    real(wp),parameter :: four_third = four/three
    real(wp),parameter :: five_half  = five/two
    real(wp),parameter :: three_half = three/two
    
    !public routines:
    public :: solve_lambert_izzo
    public :: solve_lambert_gooding
    
    contains

!*****************************************************************************************

!*****************************************************************************************
    subroutine solve_lambert_izzo(r1,r2,tof,mu,long_way,multi_revs,v1,v2,status_ok)
!*****************************************************************************************
!****f* vector_module/lambert_module
!
!  NAME
!    solve_lambert_izzo
!
!  DESCRIPTION
!    Solve Lambert's problem using Izzo's method.
!
!  SEE ALSO
!    [1] D. Izzo, "Revisiting Lambert's Problem"
!        http://arxiv-web3.library.cornell.edu/abs/1403.2705
!        [v2] Tue, 24 Jun 2014 13:08:37 GMT (606kb,D)
!    [2] https://github.com/esa/pykep
!    [3] R. A. Battin, "An Introduction to the Mathematics and Methods of 
!        Astrodynamics (Revised Edition)", AIAA Education Series, 1999.
!
!  SOURCE

    implicit none

    real(wp),dimension(3),intent(in)                :: r1            ! first cartesian position [km]
    real(wp),dimension(3),intent(in)                :: r2            ! second cartesian position [km]
    real(wp),intent(in)                             :: tof           ! time of flight [sec]
    real(wp),intent(in)                             :: mu            ! gravity parameter [km^3/s^2]
    logical,intent(in)                              :: long_way      ! when true, do "long way" (>pi) transfers
    integer,intent(in)                              :: multi_revs    ! maximum number of multi-rev solutions to compute
    real(wp),dimension(:,:),allocatable,intent(out) :: v1            ! vector containing 3d arrays with the cartesian components of the velocities at r1
    real(wp),dimension(:,:),allocatable,intent(out) :: v2            ! vector containing 3d arrays with the cartesian components of the velocities at r2
    logical,intent(out)                             :: status_ok     ! true if everything is OK

    !local variables:
    real(wp),dimension(:),allocatable :: x
    real(wp),dimension(3) :: r1_hat,r2_hat,h_hat,it1,it2,c
    real(wp) :: s,cmag,lambda2,lambda3,lambda5,t,t00,t0,t1,r1mag,r2mag,&
                d3t,d2t,dt,err,t_min,x_old,x_new,term,lambda,&
                gamma,rho,sigma,vr1,vt1,vr2,vt2,y,vt,ly
    integer :: n_solutions,it,m_nmax,i,iter

    !tolerances are from [2]
    integer,parameter   :: max_halley_iters = 12         !for halley iterations
    real(wp),parameter  :: halley_tol       = 1e-13_wp   !
    real(wp),parameter  :: htol_singlerev   = 1e-5_wp    !for householder iterations
    real(wp),parameter  :: htol_multirev    = 1e-8_wp    !

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
    !*************************************************************************************
    !    Householder root solver for x.
    !*************************************************************************************

        implicit none

        integer                 :: it
        real(wp),intent(in)     :: t
        real(wp),intent(inout)  :: x    !input is initial guess
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

    !*************************************************************************************
        end function householder
    !*************************************************************************************

    !*************************************************************************************
        subroutine dtdx(dt,d2t,d3t,x,t)
    !*************************************************************************************
    !    Compute 1st-3rd derivatives for the Householder iterations.
    !*************************************************************************************

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

    !*************************************************************************************
        end subroutine dtdx
    !*************************************************************************************

    !*************************************************************************************
        subroutine compute_tof(x,n,tof)
    !*************************************************************************************
    !    Compute time of flight from x
    !*************************************************************************************

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

    !*************************************************************************************
        end subroutine compute_tof
    !*************************************************************************************

    !*************************************************************************************
        pure function hypergeo(x) result(f)
    !*************************************************************************************
    !    Evaluate the Gaussian (or ordinary) hypergeometric function: F(3,1,5/2,x)
    !    See Ref. [3], p. 34.
    !*************************************************************************************

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

    !*************************************************************************************
        end function hypergeo
    !*************************************************************************************

!*****************************************************************************************
    end subroutine solve_lambert_izzo
!*****************************************************************************************

!*****************************************************************************************
    subroutine solve_lambert_gooding(r1,r2,tof,mu,long_way,multi_revs,v1,v2,status_ok)
!*****************************************************************************************
!****f* vector_module/lambert_module
!
!  NAME
!    solve_lambert_gooding
!
!  DESCRIPTION
!    Solve Lambert's problem using Gooding's method.
!
!  SEE ALSO
!    [1] R. H, Gooding. "A procedure for the solution of Lambert's orbital 
!        boundary-value problem" Celestial Mechanics and Dynamical Astronomy,
!        vol. 48, no. 2, 1990, p. 145-165.
!        http://adsabs.harvard.edu/abs/1990CeMDA..48..145G
!    [2] A. Klumpp, "Performance Comparision of Lambert and Kepler Algorithms",
!        JPL Interoffice Memorandum, 314.1-0426-ARK, Jan 2, 1991.
!        http://derastrodynamics.com/docs/lambert_papers_v1.zip
!
!  SOURCE

    implicit none

    real(wp),dimension(3),intent(in)                :: r1         ! first cartesian position [km]
    real(wp),dimension(3),intent(in)                :: r2         ! second cartesian position [km]
    real(wp),intent(in)                             :: tof        ! time of flight [sec]
    real(wp),intent(in)                             :: mu         ! gravity parameter [km^3/s^2]
    logical,intent(in)                              :: long_way   ! when true, do "long way" (>pi) transfers
    integer,intent(in)                              :: multi_revs ! maximum number of multi-rev solutions to compute
    real(wp),dimension(:,:),allocatable,intent(out) :: v1         ! vector containing 3d arrays with the cartesian components of the velocities at r1
    real(wp),dimension(:,:),allocatable,intent(out) :: v2         ! vector containing 3d arrays with the cartesian components of the velocities at r2
    logical,intent(out)                             :: status_ok  ! true if everything is OK

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
        if (long_way) then !>pi
            ta    =  num_revs * twopi + (twopi - pa)
            rho   = -r1xr2_hat    
        else !<pi
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
    !*************************************************************************************
    !    Gooding support routine
    !    Note: this contains the modification from [2]
    !*************************************************************************************
    
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

    !*************************************************************************************
        end subroutine vlamb
    !*************************************************************************************

    !*************************************************************************************
        subroutine tlamb(m,q,qsqfm1,x,n,t,dt,d2t,d3t)
    !*************************************************************************************
    !    Gooding support routine
    !*************************************************************************************

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

    !*************************************************************************************
        end subroutine tlamb
    !*************************************************************************************

    !*************************************************************************************
        pure function d8rt(x)
    !*************************************************************************************
    !    8th root function, used by xlamb
    !*************************************************************************************

        implicit none
    
        real(wp) :: d8rt
        real(wp),intent(in) :: x
    
        d8rt = sqrt(sqrt(sqrt(x)))
    
    !*************************************************************************************
        end function d8rt
    !*************************************************************************************

    !*************************************************************************************
        subroutine xlamb(m,q,qsqfm1,tin,n,x,xpl)
    !*************************************************************************************
    !    Gooding support routine
    !*************************************************************************************

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
      
    !*************************************************************************************
        end subroutine xlamb
    !*************************************************************************************
    
!*****************************************************************************************
    end subroutine solve_lambert_gooding
!*****************************************************************************************

!*****************************************************************************************
    end module lambert_module
!*****************************************************************************************