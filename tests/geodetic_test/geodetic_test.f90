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
    real(wp) :: ax  = 6378137.0_wp    !! Ellipsoid : WGS84 Equatorial radius
    real(wp) :: ay  = 6378137.0_wp    !! Ellipsoid : WGS84 Equatorial radius (same)
    real(wp) :: b   = 6356752.3142_wp !! Ellipsoid : WGS84 Polar radius
    real(wp),parameter :: tol = 1.0e-13_wp !! tolerance
    real(wp),parameter :: test_tol = 1.0e-6_wp !! tolerance for a failed test

    real(wp) :: h, tmp, phi, lambda, phi_, lambda_, h_
    real(wp) :: phi_error, lambda_error, h_error
    integer :: i !! counter
    real(wp),dimension(3) :: err !! error vector
    real(wp),dimension(3) :: r !! cartesian position vector
    integer :: isize !! size if `iseed`
    integer,dimension(:),allocatable :: iseed !! for random number generator

    real(wp) :: tstart, tstop

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' geodetic_test'
    write(*,*) '---------------'
    write(*,*) ''

    write(*,*) ''
    write(*,*) ' 1. Triaxial'
    write(*,*) ''

    call random_seed(size=isize)
    allocate(iseed(isize)); iseed = 42
    call random_seed(put=iseed)

    tmp = zero
    call cpu_time(tstart)
    do i = 1, n_repeat

        h      = get_random_number(1.0_wp, 10000.0_wp)
        lambda = get_random_number(0.0_wp, 360.0_wp) * deg2rad
        phi    = get_random_number(-90.0_wp, 90.0_wp) * deg2rad

        call geodetic_to_cartesian_triaxial(ax, ay, b, phi, lambda, h, r)
        call cartesian_to_geodetic_triaxial(ax, ay, b, r(1), r(2), r(3), tol, phi_, lambda_, h_)

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

        h      = get_random_number(1.0_wp, 10000.0_wp)
        lambda = get_random_number(0.0_wp, 360.0_wp) * deg2rad
        phi    = get_random_number(-90.0_wp, 90.0_wp) * deg2rad

        call geodetic_to_cartesian_triaxial(ax, ay, b, phi, lambda, h, r)
        tmp = tmp + sum(r) !compute something so loop isn't optimized away

    end do
    call cpu_time(tstop)
    write(*,'(A30,1X,E30.16,1X,F17.1,1X,A)') 'geodetic_to_cartesian_triaxial : Speed', tmp, n_repeat/(tstop-tstart), 'cases/sec'

    !----------------
    tmp = zero
    call cpu_time(tstart)
    do i = 1, n_repeat

        ! point on or outside the body:
        r = get_random_number(6378137.0_wp, 7378137.0_wp) * &
                unit([get_random_number(-10000.0_wp, 10000.0_wp),&
                      get_random_number(-10000.0_wp, 10000.0_wp),&
                      get_random_number(-10000.0_wp, 10000.0_wp) ])

        call cartesian_to_geodetic_triaxial(ax, ay, b, r(1), r(2), r(3), tol, phi_, lambda_, h_)
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

        h      = get_random_number(1.0_wp, 10000.0_wp)
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
    write(*,*) ' 3. Comparison: heikkinen vs. triaxial'
    write(*,*) ''

    tmp = zero
    call cpu_time(tstart)
    do i = 1, n_repeat

        h      = get_random_number(1.0_wp, 10000.0_wp)
        lambda = get_random_number(0.0_wp, 360.0_wp) * deg2rad
        phi    = get_random_number(-90.0_wp, 90.0_wp) * deg2rad

        call geodetic_to_cartesian_triaxial(ax, ay, b, phi, lambda, h, r)
        call cartesian_to_geodetic_triaxial(ax, ay, b, r(1), r(2), r(3), tol, phi_, lambda_, h_)
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

    !----------------------------------------
    ay  = 6378147.0_wp  ! make it a triaxial ellipsoid

    write(*,*) ''
    write(*,*) ' 4. Comparison: triaxial vs. CartesianIntoGeodeticI'
    write(*,*) ''

    tmp = zero
    call cpu_time(tstart)
    do i = 1, n_repeat

        !  what is meant by "octant" for the CartesianIntoGeodeticI inputs?
        !  the lat and long values returned are not correct for general inputs
        !  do we have to do some transformation of the inputs ???
        !  for now, just generate points from 0 -> 90 deg

        h      = get_random_number(1.0_wp, 10000.0_wp)
        lambda = get_random_number(0.0_wp, 90.0_wp) * deg2rad
        phi    = get_random_number(0.0_wp, 90.0_wp) * deg2rad

        call geodetic_to_cartesian_triaxial(ax, ay, b, phi, lambda, h, r)
        call cartesian_to_geodetic_triaxial(ax, ay, b, r(1), r(2), r(3), tol, phi_, lambda_, h_)
        call CartesianIntoGeodeticI(ax, ay, b, r(1), r(2), r(3), phi_, lambda_, h_, error=tol)

        phi_error    = rel_error(phi_ ,    phi,    .true.)
        lambda_error = rel_error(lambda_ , lambda, .true.)
        h_error      = rel_error(h_ ,      h,      .false.)
        err = [phi_error,lambda_error,h_error]

        if (any(err > test_tol)) then
            write(*,*) 'x,y,z:', r
            write(*,*) 'var   ', 'CartesianIntoGeodeticI      ', 'cartesian_to_geodetic_triaxial'
            write(*,*) 'lat:  ', phi_*rad2deg,     phi*rad2deg
            write(*,*) 'lon:  ', lambda_*rad2deg , lambda*rad2deg
            write(*,*) 'alt:  ', h_ ,      h
            write(*,*) 'errors: ', phi_error, lambda_error, h_error
            error stop 'FAILURE'
        end if
        tmp = tmp + norm2(err) !compute something so loop isn't optimized away

    end do
    call cpu_time(tstop)
    write(*,'(A30,1X,E30.16,1X,F17.1,1X,A)') 'Comparison', tmp, n_repeat/(tstop-tstart), 'cases/sec'

    write(*,*) ''
    write(*,*) ' 5. Comparison: triaxial vs. CartesianIntoGeodeticII'
    write(*,*) ''

    tmp = zero
    call cpu_time(tstart)
    do i = 1, n_repeat

        !  what is meant by "octant" for the CartesianIntoGeodeticI inputs?
        !  the lat and long values returned are not correct for general inputs
        !  do we have to do some transformation of the inputs ???
        !  for now, just generate points from 0 -> 90 deg

        h      = get_random_number(1.0_wp, 10000.0_wp)
        lambda = get_random_number(0.0_wp, 90.0_wp) * deg2rad
        phi    = get_random_number(0.0_wp, 90.0_wp) * deg2rad

        call geodetic_to_cartesian_triaxial(ax, ay, b, phi, lambda, h, r)
        call cartesian_to_geodetic_triaxial(ax, ay, b, r(1), r(2), r(3), tol, phi_, lambda_, h_)
        call CartesianIntoGeodeticII(ax, ay, b, r(1), r(2), r(3), phi_, lambda_, h_, error=tol)

        phi_error    = rel_error(phi_ ,    phi,    .true.)
        lambda_error = rel_error(lambda_ , lambda, .true.)
        h_error      = rel_error(h_ ,      h,      .false.)
        err = [phi_error,lambda_error,h_error]

        if (any(err > test_tol)) then
            write(*,*) 'x,y,z:', r
            write(*,*) 'var   ', 'CartesianIntoGeodeticII     ', 'cartesian_to_geodetic_triaxial'
            write(*,*) 'lat:  ', phi_*rad2deg,     phi*rad2deg
            write(*,*) 'lon:  ', lambda_*rad2deg , lambda*rad2deg
            write(*,*) 'alt:  ', h_ ,      h
            write(*,*) 'errors: ', phi_error, lambda_error, h_error
            error stop 'FAILURE'
        end if
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