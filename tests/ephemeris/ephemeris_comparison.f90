!*****************************************************************************************
!>
!  Comparison of the JPL ephemeris with the simple analytic version.

    program ephemeris_comparison

    use fortran_astrodynamics_toolkit, wp => fat_wp

    implicit none

    real(wp)              :: jd
    real(wp)              :: jd0
    real(wp),dimension(6) :: rv
    real(wp),dimension(3) :: r_s,v_s
    real(wp)              :: max_err_s,err_s
    integer               :: ntarg,nctr
    type(jpl_ephemeris)   :: eph405
    logical               :: status_ok
    integer               :: i
    integer               :: jd_start, jd_end

    character(len=*),parameter :: ephemeris_file_405 = '../eph/JPLEPH.405' !! JPL DE405 ephemeris file

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' ephemeris_comparison'
    write(*,*) '---------------'
    write(*,*) ''

    !initialize:
    call eph405%initialize(filename=ephemeris_file_405,status_ok=status_ok)

    if (status_ok) then

        ntarg = 10     ! moon
        nctr  = 3      ! earth

        max_err_s = zero
        jd_start  = julian_day(2016,1,1)
        jd_end    = julian_day(2050,1,1)

        do i = jd_start, jd_end

            jd = real(i,wp)

            call eph405%get_state(jd, ntarg, nctr, rv, status_ok)
            if (.not.status_ok) error stop 'error computing state.'

            call simpson_lunar_ephemeris(jd, r_s, v_s)

            if (i==jd_end) then
                write(*,'(A,1X,F15.2)') 'Position of Moon wrt Earth at jd=',jd
                write(*,'(A,1X,*(F25.16,1X))') 'DE405', rv(1:3), rv(4:6)
                write(*,'(A,1X,*(F25.16,1X))') 'SIMPS', r_s, v_s
                write(*,'(A,1X,*(F25.16,1X))') 'DE405-SIMPS ', &
                            norm2(rv(1:3) - r_s), norm2(rv(4:6) - v_s)
            end if

            err_s = norm2(rv(1:3) - r_s)

            if (err_s > max_err_s) max_err_s = err_s

        end do
        call eph405%close()

        write(*,*) 'max errors:'
        write(*,'(A,1X,*(F25.16,1X))') 'DE405-SIMPS ', max_err_s

    else
        write(*,*) 'Error opening DE405 ephemeris file'
    end if

!*****************************************************************************************
    end program ephemeris_comparison
!*****************************************************************************************
