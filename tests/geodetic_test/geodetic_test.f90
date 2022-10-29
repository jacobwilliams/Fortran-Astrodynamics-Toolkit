!*****************************************************************************************
!>
!  Unit test for the [[geodesy_module]] geodetic routines.

program geodetic_test

    use kind_module,    only: wp
    use random_module,  only: get_random_number
    use math_module,    only: wrap_angle
    use geodesy_module
    use numbers_module
    use conversion_module
    use vector_module

    implicit none

    integer,parameter :: n_repeat = 1000  !! number of times to repeat speed test
    real(wp),parameter :: ax  = 6378137.0_wp    !! Ellipsoid : WGS84 Equatorial radius
    real(wp),parameter :: ay  = 6378137.0_wp    !! Ellipsoid : WGS84 Equatorial radius (same)
    real(wp),parameter :: b   = 6356752.3142_wp !! Ellipsoid : WGS84 Polar radius
    real(wp),parameter :: tol = 1.0e-12_wp !! tolerance

    real(wp),parameter :: test_tol = 1.0e-9_wp !! tolerance for a failed test

    real(wp) :: h, tmp, phi, lambda, xi, yi, zi, phi_, lambda_, h_
    real(wp) :: phi_error, lambda_error, h_error
    integer :: i !! counter
    real(wp),dimension(3) :: err !! error vector
    real(wp),dimension(3) :: r

    real(wp) :: tstart, tstop

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' geodetic_test'
    write(*,*) '---------------'
    write(*,*) ''

    write(*,*) ''
    write(*,*) ' 1. Triaxial'
    write(*,*) ''

    tmp = zero
    call cpu_time(tstart)
    do i = 1, n_repeat

        h      = get_random_number(6778.0_wp, 10000.0_wp)
        lambda = get_random_number(0.0_wp, 360.0_wp) * deg2rad
        phi    = get_random_number(-90.0_wp, 90.0_wp) * deg2rad

        call geodetic_to_cartesian_triaxial(ax, ay, b, phi, lambda, h, xi, yi, zi)
        call cartesian_to_geodetic_triaxial(ax, ay, b, xi, yi, zi, tol, phi_, lambda_, h_)

        phi_error    = rel_error(phi_ ,    phi,    .true.)
        lambda_error = rel_error(lambda_ , lambda, .true.)
        h_error      = rel_error(h_ ,      h,      .false.)
        err = [phi_error,lambda_error,h_error]

        if (any(err > test_tol)) then
        write(*,*) phi_error, lambda_error, h_error
        error stop 'FAILURE'
        end if
        tmp = tmp + norm2(err) !compute something so loop isn't optimized away

    end do
    call cpu_time(tstop)
    write(*,'(A30,1X,E30.16,1X,F17.1,1X,A)') 'Consistency', tmp, n_repeat/(tstop-tstart), 'cases/sec'

    !----------------
    tmp = zero
    call cpu_time(tstart)
    do i = 1, n_repeat

        h      = get_random_number(6778.0_wp, 10000.0_wp)
        lambda = get_random_number(0.0_wp, 360.0_wp) * deg2rad
        phi    = get_random_number(-90.0_wp, 90.0_wp) * deg2rad

        call geodetic_to_cartesian_triaxial(ax, ay, b, phi, lambda, h, xi, yi, zi)
        tmp = tmp + xi + yi + zi !compute something so loop isn't optimized away

    end do
    call cpu_time(tstop)
    write(*,'(A30,1X,E30.16,1X,F17.1,1X,A)') 'geodetic_to_cartesian_triaxial : Speed', tmp, n_repeat/(tstop-tstart), 'cases/sec'

    !----------------
    tmp = zero
    call cpu_time(tstart)
    do i = 1, n_repeat

        !  xi = get_random_number(-10000.0_wp, 10000.0_wp)
        !  yi = get_random_number(-10000.0_wp, 10000.0_wp)
        !  zi = get_random_number(-10000.0_wp, 10000.0_wp)
        ! point on or outside the body:
        r = get_random_number(6378137.0_wp, 7378137.0_wp) * &
                unit([get_random_number(-10000.0_wp, 10000.0_wp),&
                      get_random_number(-10000.0_wp, 10000.0_wp),&
                      get_random_number(-10000.0_wp, 10000.0_wp) ])
        xi = r(1); yi = r(2); zi = r(3)

        call cartesian_to_geodetic_triaxial(ax, ay, b, xi, yi, zi, tol, phi_, lambda_, h_)
        tmp = tmp + phi_ + lambda_ + h_ !compute something so loop isn't optimized away

    end do
    call cpu_time(tstop)
    write(*,'(A30,1X,E30.16,1X,F17.1,1X,A)') 'cartesian_to_geodetic_triaxial : Speed', tmp, n_repeat/(tstop-tstart), 'cases/sec'


    write(*,*) ''
    write(*,*) ' 1. Oblate Spheroid'
    write(*,*) ''

    !----------------
    tmp = zero
    call cpu_time(tstart)
    do i = 1, n_repeat

        h      = get_random_number(6778.0_wp, 10000.0_wp)
        lambda = get_random_number(0.0_wp, 360.0_wp) * deg2rad
        phi    = get_random_number(-90.0_wp, 90.0_wp) * deg2rad

        call geodetic_to_cartesian(ax,b,phi,lambda,h,r)
        tmp = tmp + sum(r) !compute something so loop isn't optimized away

    end do
    call cpu_time(tstop)
    write(*,'(A30,1X,E30.16,1X,F17.1,1X,A)') 'geodetic_to_cartesian : Speed', tmp, n_repeat/(tstop-tstart), 'cases/sec'

    !----------------
    tmp = zero
    call cpu_time(tstart)
    do i = 1, n_repeat

        ! point on or outside the body:
        r = get_random_number(6378137.0_wp, 7378137.0_wp) * &
            unit([get_random_number(-10000.0_wp, 10000.0_wp),&
                    get_random_number(-10000.0_wp, 10000.0_wp),&
                    get_random_number(-10000.0_wp, 10000.0_wp) ])

        call heikkinen(r, ax, b, h_, lambda_, phi_)
        tmp = tmp + phi_ + lambda_ + h_ !compute something so loop isn't optimized away

    end do
    call cpu_time(tstop)
    write(*,'(A30,1X,E30.16,1X,F17.1,1X,A)') 'heikkinen : Speed', tmp, n_repeat/(tstop-tstart), 'cases/sec'

    write(*,*) ''
    write(*,*) ' 3. Comparison'
    write(*,*) ''

    tmp = zero
    call cpu_time(tstart)
    do i = 1, n_repeat

        h      = get_random_number(6778.0_wp, 10000.0_wp)
        lambda = get_random_number(0.0_wp, 360.0_wp) * deg2rad
        phi    = get_random_number(-90.0_wp, 90.0_wp) * deg2rad

        call geodetic_to_cartesian_triaxial(ax, ay, b, phi, lambda, h, xi, yi, zi)
        r = [xi, yi, zi]
        call cartesian_to_geodetic_triaxial(ax, ay, b, xi, yi, zi, tol, phi_, lambda_, h_)
        call heikkinen(r, ax, b, h, lambda, phi)

        phi_error    = rel_error(phi_ ,    phi,    .true.)
        lambda_error = rel_error(lambda_ , lambda, .true.)
        h_error      = rel_error(h_ ,      h,      .false.)
        err = [phi_error,lambda_error,h_error]

        if (any(err > test_tol)) then
            write(*,*) phi_error, lambda_error, h_error
            error stop 'FAILURE'
        end if
        !write(*,*) 'err: ', err
        tmp = tmp + norm2(err) !compute something so loop isn't optimized away

    end do
    call cpu_time(tstop)
    write(*,'(A30,1X,E30.16,1X,F17.1,1X,A)') 'Comparison', tmp, n_repeat/(tstop-tstart), 'cases/sec'

    contains

    pure real(wp) function rel_error(r1,r2,angle) result(e)
        !! compute the relative error
        real(wp),intent(in) :: r1, r2
        logical,intent(in) :: angle

        real(wp) :: diff, denom

        denom = one
        if (angle) then
            diff = wrap_angle(r2 - r1)
            if (r1 /= zero) denom = wrap_angle(r1)
        else
            diff = r2 - r1
            if (r1 /= zero) denom = r1
        end if

        e = abs(diff / denom)

    end function rel_error

end program geodetic_test
!*****************************************************************************************