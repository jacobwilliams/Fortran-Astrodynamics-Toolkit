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

    implicit none

    integer,parameter :: n_repeat = 1000  !! number of times to repeat speed test
    real(wp),parameter :: ax  = 6378137.0_wp    !! Ellipsoid : WGS84 Equatorial radius
    real(wp),parameter :: ay  = 6378137.0_wp    !! Ellipsoid : WGS84 Equatorial radius (same)
    real(wp),parameter :: b   = 6356752.3142_wp !! Ellipsoid : WGS84 Polar radius
    real(wp),parameter :: tol = 1.0e-12_wp !! tolerance

    real(wp) :: h, tmp, phi, lambda, xi, yi, zi, phi_, lambda_, h_
    integer :: i !! counter
    real(wp),dimension(3) :: err !! error vector

    real(wp) :: tstart, tstop

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' geodetic_test'
    write(*,*) '---------------'
    write(*,*) ''

    write(*,*) ''
    write(*,*) ' 1. oblate spheroid test'
    write(*,*) ''

    tmp = zero
    call cpu_time(tstart)
    do i = 1, n_repeat

         h   = get_random_number(6778.0_wp, 10000.0_wp)
         lambda = get_random_number(0.0_wp, 360.0_wp) * deg2rad
         phi = get_random_number(-90.0_wp, 90.0_wp) * deg2rad

         call geodetic_to_cartesian_triaxial(ax, ay, b, phi, lambda, h, xi, yi, zi)
         call cartesian_to_geodetic_triaxial(ax, ay, b, xi, yi, zi, tol, phi_, lambda_, h_)

         err = [wrap_angle(phi_ - phi), &
                wrap_angle(lambda_ - lambda), &
                h_ - h]
         !write(*,*) err

         tmp = tmp + norm2(err) !compute something so loop isn't optimized away

    end do
    call cpu_time(tstop)

    write(*,'(A10,1X,E30.16,1X,F13.6,1X,A)') 'cartesian_to_geodetic_triaxial', tmp, tstop-tstart, 'sec'

end program geodetic_test
!*****************************************************************************************