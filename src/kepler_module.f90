!*******************************************************************************
!>
!  Kepler propagation routines.

    module kepler_module

    use kind_module,    only: wp
    use numbers_module

    implicit none

    public :: kepler_shepperd

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  Kepler propagation using Shepperd's method.
!
!### Reference
!  * S.W. Shepperd, "Universal Keplerian State Transition Matrix".
!    Celestial Mechanics 35(1985) p. 129-144.
!  * Rody P.S. Oldenhuis, `progress_orbitM` Matlab function (BSD license).

    subroutine kepler_shepperd(mu,rv1,dt,rv2,istat)

    implicit none

    real(wp),intent(in)               :: mu    !! gravitational parameter
    real(wp),dimension(6),intent(in)  :: rv1   !! initial position,velocity vector
    real(wp),intent(in)               :: dt    !! time step
    real(wp),dimension(6),intent(out) :: rv2   !! final position,velocity vector
    integer,intent(out),optional      :: istat !! status flag (if not present, warnings are printed):
                                               !!
                                               !! `0` : all is well
                                               !! `-1` : failed to converge in time loop
                                               !! `-2` : failed to converge in `g` loop

    ! note: these could also be (optional) input:
    real(wp),parameter :: min_dt     = epsilon(one) !! time step considred as zero (sec)
    real(wp),parameter :: ttol       = 1.0e-8_wp    !! tolerance for time (sec)
    integer,parameter  :: max_iter   = 1000         !! max iterations for time
    real(wp),parameter :: gtol       = 1.0e-14_wp   !! tolerance for g
    integer,parameter  :: max_iter_g = 1000         !! max iterations for g

    real(wp),dimension(3) :: r1  !! initial position vector
    real(wp),dimension(3) :: v1  !! initial velocity vector
    real(wp) :: r1mag,nu0,beta,p,deltau,u,t,deltat,bu,q,&
                a,b,gprev,u0w2,u1w2,uu,u0,u1,u2,u3,r,&
                ff,f,gg,g
    integer :: n,iter,k,d,l,iterg

    if (present(istat)) istat = 0

    if (abs(dt) <= min_dt) then

        rv2 = rv1

    else

        r1 = rv1(1:3)
        v1 = rv1(4:6)
        r1mag = norm2(r1)
        nu0 = dot_product(r1,v1)
        beta = two*mu/r1mag - dot_product(v1,v1)

        if (beta > zero) then
            p = twopi*mu*beta**(-three/two)
            n = floor((dt + p/two - two*nu0/beta)/p)
            deltau = twopi*n*beta**(-five/two)
        else
            deltau = zero
        end if

        u = zero
        t = zero
        iter = 0
        deltat = t-dt
        do while (abs(deltat) > ttol)
            iter = iter + 1
            bu = beta*u*u
            q  = bu/(one + bu)
            a = one
            b = one
            g = one
            n = 0
            k = -9
            d = 15
            l = 3
            gprev = huge(one)
            iterg = 0
            do while (abs(g-gprev) > gtol)
                iterg = iterg + 1
                k = -k
                l = l + 2
                d = d + 4*l
                n = n + (1+k)*l
                a = d/(d - n*a*q)
                b = (a-one)*b
                gprev = g
                g = g + b
                if (iterg==max_iter_g) then
                    if (present(istat)) then
                        istat = -2
                    else
                        write(*,*) 'Warning: kepler_shepperd failed to converge in g iteration'
                    end if
                    exit
                end if
            end do

            u0w2 = one - two*q
            u1w2 = two*(one-q)*u
            uu = 16.0_wp/15.0_wp*u1w2**5*g + deltau
            u0 = two*u0w2**2-one
            u1 = two*u0w2*u1w2
            u2 = two*u1w2**2
            u3 = beta*uu + u1*u2/three
            r = r1mag*u0 + nu0*u1 + mu*u2
            t = r1mag*u1 + nu0*u2 + mu*u3
            deltat = t - dt

            ! u = u - deltat/4/(1-q)/r    ! original
            u = u - deltat/((one-q)*(four*r + deltat*beta*u)) ! Halley version from oldenhuis routine

            if (iter==max_iter) then
                if (present(istat)) then
                    istat = -1
                else
                    write(*,*) 'Warning: kepler_shepperd failed to converge in time iteration'
                end if
                exit
            end if

        end do

        ff = one - mu/r1mag*u2
        f  = -mu*u1/r/r1mag
        gg = r1mag*u1 + nu0*u2
        g  = one - mu/r*u2

        rv2 = [r1*ff + v1*gg, r1*f + v1*g]

    end if

    end subroutine kepler_shepperd
!*******************************************************************************

!*******************************************************************************
    end module kepler_module
!*******************************************************************************
