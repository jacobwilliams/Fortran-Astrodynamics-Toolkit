!*****************************************************************************************
!> author: Jacob Williams
!  date: 11/11/2014
!
!  Earth-Mars porkchop routines.
!
!  ```f2py -c porkchop.f90 -L../../lib/ -lfat -I../../lib/ -m porkchop only: generate_porkchop```
!
!@note To run requires the ```JPLEPH_2000-2100.405``` ephemeris file.

    module porkchop

    use fortran_astrodynamics_toolkit, wp => fat_wp

    implicit none

    private

    public :: generate_porkchop

    contains
!*****************************************************************************************

    !************************************************************
    !>
    !  Generate the data for an Earth-Mars porkchop plot.
    !  This routine is meant to be called from Python.

        subroutine generate_porkchop(   n_t0,n_dt,&
                                        y,m,d,&
                                        initial_t0,&
                                        delta_t0,&
                                        final_t0,&
                                        initial_dt,&
                                        delta_dt,&
                                        final_dt,&
                                        t0_out,tf_out,&
                                        c30_lt_pi_out,c30_gt_pi_out,&
                                        c3f_lt_pi_out,c3f_gt_pi_out)

        implicit none

        !subroutine arguments:
        integer,intent(in) :: n_t0, n_dt      ! size of the output arrays:
                                              !   n_t0 = size([(i,i=t0_days,tf_days,delta_t_days)])
                                              !   n_dt = size([(i,i=dt0_days,dtf_days,dt_step_days)])
                                              ! [these are so we can call routine from Python]
        integer,intent(in) :: y,m,d           ! yyyy,mm,dd of base epoch
        integer,intent(in) :: initial_t0      !! initial departure epoch [days from base epoch]
        integer,intent(in) :: delta_t0        !! departure epoch step size [days]
        integer,intent(in) :: final_t0        !! final departure epoch [days from base epoch]
        integer,intent(in) :: initial_dt      !! initial time of flight [days]
        integer,intent(in) :: delta_dt        !! time of flight step [days]
        integer,intent(in) :: final_dt        !! final time of flight [days]
        real(8),dimension(n_t0,n_dt),intent(out) :: t0_out,tf_out    !! departure and arrival [python day number counts]
        real(8),dimension(n_t0,n_dt),intent(out) :: c30_lt_pi_out    !! initial c3 (<pi solution) [km/s]
        real(8),dimension(n_t0,n_dt),intent(out) :: c30_gt_pi_out    !! initial c3 (>pi solution) [km/s]
        real(8),dimension(n_t0,n_dt),intent(out) :: c3f_lt_pi_out    !! final c3 (<pi solution) [km/s]
        real(8),dimension(n_t0,n_dt),intent(out) :: c3f_gt_pi_out    !! final c3 (>pi solution) [km/s]

        !parameters:
        real(wp),parameter :: day2sec       = one * 24.0d0 * 3600.0d0  !! conversion factor
        real(wp),parameter :: sec2day       = one / day2sec            !! conversion factor
        real(wp),parameter :: mu            = 132712440018.0d0         !! Sun gravity parameter [km^3/s^2]
        real(wp),parameter :: j2000_datenum = 730120.5d0               !! J2000 epoch - python date number (matplotlib.dates.date2num)
        real(wp),parameter :: j2000_jd      = 2451545.0d0              !! J2000 epoch - julian date
        character(len=*),parameter :: ephemeris_file = '../../eph/JPLEPH.405'    !! ephemeris file

        !local variables:
        real(8),dimension(3) :: r1                   !! first cartesian position [km]
        real(8),dimension(3) :: r2                   !! second cartesian position [km]
        real(8) :: tof                               !! time of flight [sec]
        logical :: long_way                          !! when true, do "long way" (>pi) transfers
        real(8),dimension(:,:),allocatable :: v1     !! vector containing 3d arrays with the cartesian components of the velocities at r1
        real(8),dimension(:,:),allocatable :: v2     !! vector containing 3d arrays with the cartesian components of the velocities at r2
        logical :: status_ok                         !! true if everything is OK
        real(8),dimension(6) :: rv_earth_departure   !! earth position vector at departure [km]
        real(8),dimension(6) :: rv_mars_arrival      !! mars position vector at arrival [km]
        real(8) :: jd0                               !! initial julian date epoch [days]
        real(8) :: jd_departure                      !! julian date of earth departure [days]
        real(8) :: jd_arrival                        !! julian date of mars arrival [days]
        real(8) :: dt                                !! time of flight [days]
        real(8) :: dv0                               !! earth departure delta-v [km/s]
        real(8) :: dvf                               !! mars arrival delta-v [km/s]
        real(8) :: c30_lt_pi                         !! earth departure C3 (<pi transfer) [km^2/s^2]
        real(8) :: c3f_lt_pi                         !! mars arrival C3 (<pi transfer) [km^2/s^2]
        real(8) :: c30_gt_pi                         !! earth departure C3 (>pi transfer) [km^2/s^2]
        real(8) :: c3f_gt_pi                         !! mars arrival C3 (>pi transfer) [km^2/s^2]
        integer :: i                                 !! counter
        integer :: j                                 !! counter
        integer :: n                                 !! lambert multi-rev input
        real(8) :: dep                               !! temp variable
        integer :: t0_days,delta_t_days,tf_days      ! departure time ranges [days]
        integer :: dt0_days,dt_step_days,dtf_days    ! flight time ranges [days]
        integer :: ii,jj                             ! counters
        type(jpl_ephemeris) :: eph405      !! the ephemeris
        !test:
        t0_out = 0.0d0
        tf_out = 0.0d0
        c30_lt_pi_out = 0.0d0
        c30_gt_pi_out = 0.0d0
        c3f_lt_pi_out = 0.0d0
        c3f_gt_pi_out = 0.0d0

        !get inputs:
        t0_days         = initial_t0
        delta_t_days    = delta_t0
        tf_days         = final_t0
        dt0_days        = initial_dt
        dt_step_days    = delta_dt
        dtf_days        = final_dt

        call eph405%initialize(filename=ephemeris_file,status_ok=status_ok)
        if (.not. status_ok) error stop 'error initializing ephemeris'

        jd0 = jd(y,m,d)        ! initial julian date epoch [days]:
        n = 0                  ! only consider <1 revolution solutions
        ii = 0

        do i = t0_days, tf_days, delta_t_days            ! departure epoch loop

            ii = ii + 1
            jj = 0
            jd_departure = jd0 + real(i,8)                          ! jd of departure [days]
            call eph405%get_state(jd_departure,3,11,rv_earth_departure,status_ok)    ! state of earth at departure
            if (.not. status_ok) error stop 'error getting state'
            r1 = rv_earth_departure(1:3)                            ! initial r vector for lambert [km]

            do j = dt0_days, dtf_days, dt_step_days      ! flight time loop

                jj = jj + 1
                dt = real(j,8)                                      ! flight time [days]
                tof = dt*day2sec                                    ! flight time [sec]
                jd_arrival = jd_departure + dt                      ! jd of arrival [days]
                call eph405%get_state(jd_arrival,4,11,rv_mars_arrival,status_ok)     ! state of mars at arrival
                if (.not. status_ok) error stop 'error getting state'
                r2 = rv_mars_arrival(1:3)                           ! final r vector for lambert [km]

                !Earth-Mars transfers:

                !---short way:
                long_way=.false.
                call solve_lambert_gooding(r1,r2,tof,mu,long_way,n,v1,v2,status_ok)
                dv0 = norm2( v1(:,1) - rv_earth_departure(4:6) )
                dvf = norm2( rv_mars_arrival(4:6) - v2(:,1) )
                c30_lt_pi = dv0*dv0
                c3f_lt_pi = dvf*dvf

                !---long way:
                long_way=.true.
                call solve_lambert_gooding(r1,r2,tof,mu,long_way,n,v1,v2,status_ok)
                dv0 = norm2( v1(:,1) - rv_earth_departure(4:6) )
                dvf = norm2( rv_mars_arrival(4:6) - v2(:,1) )
                c30_gt_pi = dv0*dv0
                c3f_gt_pi = dvf*dvf

                !departure and arrival (python day counts):
                dep = jd_departure - j2000_jd      !departure date in days since j2000:
                t0_out(ii,jj) = j2000_datenum + dep        !departure epoch [python date - days]
                tf_out(ii,jj) = j2000_datenum + dep + dt   !arrival epoch [python date - days]

                !add to output arrays:
                c30_lt_pi_out(ii,jj) = c30_lt_pi
                c30_gt_pi_out(ii,jj) = c30_gt_pi
                c3f_lt_pi_out(ii,jj) = c3f_lt_pi
                c3f_gt_pi_out(ii,jj) = c3f_gt_pi

            end do

        end do

        !close files:
        call eph405%close()

        end subroutine generate_porkchop
    !************************************************************

    !************************************************************
    !>
    !  Returns JD, the Julian date at Greenwich noon on the
    !  specified YEAR, MONTH, and DAY.
    !
    !  Valid for any Gregorian calendar date producing a
    !  Julian date greater than zero.
    !
    !# See also
    !   * http://aa.usno.navy.mil/faq/docs/JD_Formula.php

        pure integer function jd(y,m,d)

        implicit none

        integer,intent(in) :: y   !! year (YYYY)
        integer,intent(in) :: m   !! month (MM)
        integer,intent(in) :: d   !! day (DD)

        jd = d-32075+1461*(y+4800+(m-14)/12)/4+367*&
             (m-2-(m-14)/12*12)/12-3*((y+4900+(m-14)/12)/100)/4

        end function jd
    !************************************************************

    end module porkchop
!*****************************************************************************************
