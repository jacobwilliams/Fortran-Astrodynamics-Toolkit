!*****************************************************************************************
!> author: Jacob Williams
!
!  Propagate a state using the various kepler methods and compare.

    program kepler_tests

    use kind_module,    only: wp
    use gooding_module
    use kepler_module
    use conversion_module
    use random_module
    use numbers_module

    implicit none

    real(wp),parameter :: dt1 = 0.74060940549140621E+06_wp !! time step (sec)
    real(wp),dimension(6),parameter :: rv1 = [  0.14624621707935020E+06_wp,&
                                                -0.20413097635909508E+04_wp,&
                                                -0.13806223844513239E+06_wp,&
                                                0.82912050920300651E+00_wp,&
                                                0.10284105194164532E+01_wp,&
                                                -0.19717889495193228E+01_wp ] !! initial state [km,km/s]

    logical,parameter :: run_all = .true.  !! run the full set of tests
    real(wp),parameter :: mu = 398600.0_wp !! Earth grav param.
    integer,parameter :: n_cases = 10000   !! number of cases to run
    integer,parameter :: n_methods = 4 !! number of kepler methods to test
    character(len=*),dimension(n_methods),parameter :: methods = ['gooding  ',&
                                                                  'shepperd ',&
                                                                  'goodyear ',&
                                                                  'classical'] !! method names

    integer :: i  !! counter
    real(wp) :: alpha,rp,inc,raan,w,tau,sma,per,dt
    real(wp),dimension(6) :: rv  !! cartesian state
    real(wp),dimension(6) :: e   !! universal elements
    real(wp),dimension(6) :: err !! state error estimate
    real(wp),dimension(n_methods,n_cases) :: pos_err,vel_err
    integer,dimension(n_methods,n_cases) :: istats
    integer :: imethod,istat

    pos_err = zero
    vel_err = zero
    istats  = 0

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' kepler_tests'
    write(*,*) '---------------'
    write(*,*) ''

    write(*,*) ''
    write(*,*) ' Unusual case ...'
    write(*,*) ''
    do imethod = 1,n_methods
        call prop_error(rv1,dt1,imethod,err,istat)
        write(*,'(A,1X,6(E30.17,1X),I5)') trim(methods(imethod))//' err:  ', err, istat
    end do
    write(*,*) ''

    if (run_all) then

        write(*,*) ''
        write(*,*) ' All cases ...'
        write(*,*) ''

        do i = 1, n_cases

            ! get random elements:
            rp = get_random_number(7000.0_wp,100000.0_wp)
            if (get_random_number(0.0_wp,1.0_wp)>0.5_wp) then
                ! ellipse
                sma = get_random_number(rp, 1000000.0_wp)
                per = twopi*sqrt(sma**3/mu)
                tau = get_random_number(0.0_wp,per)
                dt = get_random_number(0.0_wp,per)
            else
                ! hyperbola
                sma = get_random_number(-1000000.0_wp,-rp)
                tau = get_random_number(0.0_wp,10.0_wp) * day2sec
                dt = get_random_number(0.0_wp,10.0_wp) * day2sec
            end if
            alpha = mu / sma
            ! add a specific test for a parabola (a=0)
            if (i==n_cases) then  ! the last one
                alpha = zero
            end if
            inc   = get_random_number(0.0_wp, pi)
            raan  = get_random_number(0.0_wp, twopi)
            w     = get_random_number(0.0_wp, twopi)

            e = [alpha,rp,inc,raan,w,tau]

            call els3pv (mu, e, rv)

            do imethod = 1,n_methods
                call prop_error(rv,dt,imethod,err,istats(imethod,i))
                pos_err(imethod,i) = norm2(err(1:3))
                vel_err(imethod,i) = norm2(err(4:6))

                if (istats(imethod,i) /= 0 .or. pos_err(imethod,i)>one .or. vel_err(imethod,i)>one) then
                    write(*,*) ''
                    write(*,*) 'error for method: '//trim(methods(imethod))
                    write(*,*) 'i: ', i
                    write(*,'(A,1X,I6)')           'istat   =', istats(imethod,i)
                    write(*,'(A,1X,*(E30.16,1X))') 'e       =',e
                    write(*,'(A,1X,*(E30.16,1X))') 'pos err =',err(1:3)
                    write(*,'(A,1X,*(E30.16,1X))') 'vel err =',err(4:6)
                    write(*,*) ''
                    write(*,'(A,1X,*(E30.17,1X))') 'rv   (km,km/s) = ',rv
                    write(*,'(A,1X,*(E30.17))')    'dt (sec)       = ', dt
                    if (e(1)/=zero) write(*,'(A,1X,*(E30.17))')    'sma  (km)      = ', mu / e(1)
                    write(*,'(A,1X,*(E30.17))')    'rp   (km)      = ', e(2)
                    write(*,'(A,1X,*(E30.17))')    'inc  (deg)     = ', e(3) * rad2deg
                    write(*,'(A,1X,*(E30.17))')    'raan (deg)     = ', e(4) * rad2deg
                    write(*,'(A,1X,*(E30.17))')    'w    (deg)     = ', e(5) * rad2deg
                    write(*,'(A,1X,*(E30.17))')    'tau  (days)    = ', e(6) * sec2day
                    write(*,*) ''
                    !stop
                end if

            end do

        end do

        ! results:
        write(*,*) ''
        write(*,*) '-------------------'
        do imethod = 1,n_methods
            do i = 1, n_cases
                if (istats(imethod,i)/=0) &
                    write(*,'(A,E30.16,1X,I5)') trim(methods(imethod))//' pos error: ', pos_err(imethod,i), istats(imethod,i)
            end do
            write(*,*) ''
        end do
        do imethod = 1,n_methods
            do i = 1, n_cases
                if (istats(imethod,i)/=0) &
                    write(*,'(A,E30.16,1X,I5)') trim(methods(imethod))//' vel error: ', vel_err(imethod,i), istats(imethod,i)
            end do
            write(*,*) ''
        end do
        write(*,*) '-------------------'
        write(*,*) ''

        write(*,*) ''
        write(*,*) 'summary'
        do imethod = 1,n_methods
            write(*,*) '-------------------'
            write(*,'(A,E30.16)') trim(methods(imethod))//' max pos error: ', maxval(abs(pos_err(imethod,:)))
            write(*,'(A,E30.16)') trim(methods(imethod))//' max vel error: ', maxval(abs(vel_err(imethod,:)))
            write(*,*) '-------------------'
        end do
        write(*,*) ''

    end if

    contains

        subroutine prop_error(x0,dt,imethod,xerr,istat)

        !! propagate the state forward the backwards, and return
        !! the error vector.

        implicit none

        real(wp),dimension(6),intent(in)  :: x0
        real(wp),intent(in)               :: dt
        integer,intent(in)                :: imethod
        real(wp),dimension(6),intent(out) :: xerr
        integer,intent(out)               :: istat

        real(wp),dimension(6) :: rv2,rv3,rv4

        real(wp),parameter :: tol = 1.0e-15_wp !! tolerance for goodyear

        istat = 0

        select case (imethod)
        case(1)
            call propagate(mu,x0,dt,rv2)
            call propagate(mu,rv2,-dt,rv3)
        case(2)
            call kepler_shepperd(mu,x0,dt,rv2,istat)
            call kepler_shepperd(mu,rv2,-dt,rv3,istat)
        case(3)
            call kepler_goodyear_stienon_klumpp(x0,dt,mu,tol,rv2)
            call kepler_goodyear_stienon_klumpp(rv2,-dt,mu,tol,rv3)
        case(4)
            call kepler_classical(x0,dt,mu,rv2)
            call kepler_classical(rv2,-dt,mu,rv3)
        case default
            error stop 'invalid method'
        end select

        xerr = rv3 - x0

        end subroutine prop_error

    end program kepler_tests
!*****************************************************************************************
