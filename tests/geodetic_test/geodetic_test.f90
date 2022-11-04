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
    use ieee_arithmetic, only: ieee_is_nan

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
    integer :: icase
    character(len=:),allocatable :: method

    real(wp) :: tstart, tstop

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' geodetic_test'
    write(*,*) '---------------'
    write(*,*) ''

    write(*,*) ''
    write(*,*) ' 0. Test case 1 : Triaxial'
    write(*,*) ''

    write(*,'(A26,1X,*(A26,1X))') 'Method', 'Lat (rad)', 'Long (rad)', 'Alt (m)'
    ax = 6378137.0_wp
    ay = 6378102.0_wp
    b  = 6356752.0_wp
    r  = [-5734871.6046899008_wp, -2808462.3459780114_wp, 2937.9300139431466_wp]
    call cartesian_to_geodetic_triaxial(ax, ay, b, r(1), r(2), r(3), tol, phi_, lambda_, h_)
    write(*,'(A26,1X,*(E26.16,1X))') 'Panou', phi_, lambda_, h_
    call cartesian_to_geodetic_triaxial_2(ax, ay, b, r, tol, phi_, lambda_, h_)
    write(*,'(A26,1X,*(E26.16,1X))') 'Bektas', phi_, lambda_, h_
    call CartesianIntoGeodeticI(ax, ay, b, r(1), r(2), r(3), phi_, lambda_, h_, error=tol)
    write(*,'(A26,1X,*(E26.16,1X))') 'CartesianIntoGeodeticI', phi_, lambda_, h_
    call CartesianIntoGeodeticII(ax, ay, b, r(1), r(2), r(3), phi_, lambda_, h_, error=tol)
    write(*,'(A26,1X,*(E26.16,1X))') 'CartesianIntoGeodeticII', phi_, lambda_, h_

    ! From C++ code:
   ! Latitude (rad) = 0.000463181
   ! Longitude (rad) = -2.686201153
   !  Altitude = 7495.954101011

    write(*,*) ''
    write(*,*) ' 1. Triaxial'
    write(*,*) ''

    ax  = 6378137.0_wp    ! make it a triaxial ellipsoid
    ay  = 6378102.0_wp    ! see: https://link.springer.com/article/10.1007/s10291-020-01033-7/tables/1
    b   = 6356752.0_wp

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

        if (any(err > test_tol) .or. any(ieee_is_nan(err))) then
            write(*,*) 'r      = ', r
            write(*,*) 'h      = ', h, h_
            write(*,*) 'lambda = ', lambda, lambda_
            write(*,*) 'phi    = ', phi, phi_
            error stop 'FAILURE'
        end if
        tmp = tmp + norm2(err) !compute something so loop isn't optimized away

    end do
    call cpu_time(tstop)
    write(*,'(A40,1X,E30.16,1X,F17.1,1X,A)') 'Panou Consistency', tmp, n_repeat/(tstop-tstart), 'cases/sec'

    !----------------
    call random_seed(put=iseed)
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
    write(*,'(A40,1X,E30.16,1X,F17.1,1X,A)') 'geodetic_to_cartesian_triaxial : Speed', tmp, n_repeat/(tstop-tstart), 'cases/sec'

    !----------------
    call random_seed(put=iseed)
    tmp = zero
    call cpu_time(tstart)
    do i = 1, n_repeat

        h      = get_random_number(1.0_wp, 10000.0_wp)
        lambda = get_random_number(0.0_wp, 360.0_wp) * deg2rad
        phi    = get_random_number(-90.0_wp, 90.0_wp) * deg2rad

        call geodetic_to_cartesian_triaxial_2(ax, ay, b, phi, lambda, h, r)
        tmp = tmp + sum(r) !compute something so loop isn't optimized away

    end do
    call cpu_time(tstop)
    write(*,'(A40,1X,E30.16,1X,F17.1,1X,A)') 'geodetic_to_cartesian_triaxial_2 : Speed', tmp, n_repeat/(tstop-tstart), 'cases/sec'

    !..... triaxial speed tests ......
    do icase = 1, 4

        select case (icase)
        case(1); method = 'Panou'
        case(2); method = 'Bektas'
        case(3); method = 'CartesianIntoGeodeticI'
        case(4); method = 'CartesianIntoGeodeticII'
        end select

        call random_seed(put=iseed)
        tmp = zero
        call cpu_time(tstart)
        do i = 1, n_repeat

            ! point on or outside the body:
            r = get_random_number(6378137.0_wp, 7378137.0_wp) * &
                    unit([get_random_number(-10000.0_wp, 10000.0_wp),&
                        get_random_number(-10000.0_wp, 10000.0_wp),&
                        get_random_number(-10000.0_wp, 10000.0_wp) ])

            select case (icase)
            case(1); call cartesian_to_geodetic_triaxial(ax, ay, b, r(1), r(2), r(3), tol, phi_, lambda_, h_)
            case(2); call cartesian_to_geodetic_triaxial_2(ax, ay, b, r, tol, phi_, lambda_, h_)
            case(3); call CartesianIntoGeodeticI(ax, ay, b, r(1), r(2), r(3), phi_, lambda_, h_, error=tol)
            case(4); call CartesianIntoGeodeticII(ax, ay, b, r(1), r(2), r(3), phi_, lambda_, h_, error=tol)
            end select
            tmp = sin(tmp) + phi_ + lambda_ + h_ ! compute something so loop isn't optimized away

        end do
        call cpu_time(tstop)
        write(*,'(A40,1X,E30.16,1X,F17.1,1X,A)') method//' : Speed', tmp, n_repeat/(tstop-tstart), 'cases/sec'

    end do

    write(*,*) ''
    write(*,*) ' 1. Oblate Spheroid'
    write(*,*) ''

    ay  = ax ! make it an Oblate Spheroid

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
    write(*,*) ' 3. Comparison: heikkinen vs. Panou'
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
    ay  = 6378127.0_wp  ! make it a triaxial ellipsoid

    write(*,*) ''
    write(*,*) ' 4. Comparison: Panou vs. CartesianIntoGeodeticI'
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
            write(*,*) 'var   ', 'CartesianIntoGeodeticI      ', 'Panou             '
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
    write(*,*) ' 5. Comparison: Panou vs. CartesianIntoGeodeticII'
    write(*,*) ''

    tmp = zero
    call cpu_time(tstart)
    do i = 1, n_repeat

        !  what is meant by "octant" for the CartesianIntoGeodeticI inputs?
        !  the lat and long values returned are not correct for general inputs
        !  do we have to do some transformation of the inputs ???
        !  for now, just generate points from 0 -> 90 deg

        h      = get_random_number(500.0_wp, 10000.0_wp)
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
            write(*,*) 'var   ', 'CartesianIntoGeodeticII     ', 'Panou        '
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


    !.....................
    write(*,*) ''
    write(*,*) ' 6. Comparison: Panou vs. Bektas'
    write(*,*) ''

    tmp = zero
    call cpu_time(tstart)
    do i = 1, n_repeat

        h      = get_random_number(1.0_wp, 10000.0_wp)
        lambda = get_random_number(-180.0_wp, 180.0_wp) * deg2rad
        phi    = get_random_number(-90.0_wp, 90.0_wp) * deg2rad

        call geodetic_to_cartesian_triaxial(ax, ay, b, phi, lambda, h, r)
        call cartesian_to_geodetic_triaxial(ax, ay, b, r(1), r(2), r(3), tol, phi_, lambda_, h_)
        call cartesian_to_geodetic_triaxial_2(ax, ay, b, r, tol,phi_, lambda_, h_)

        phi_error    = rel_error(phi_ ,    phi,    .true.)
        lambda_error = rel_error(lambda_ , lambda, .true.)
        h_error      = rel_error(h_ ,      h,      .false.)
        err = [phi_error,lambda_error,h_error]

        if (any(err > test_tol)) then
            write(*,*) 'x,y,z:', r
            write(*,*) 'var   ', 'Bektas           ', 'Panou          '
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