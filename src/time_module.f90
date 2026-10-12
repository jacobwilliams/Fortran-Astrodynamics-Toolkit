!*****************************************************************************************
!> author: Jacob Williams
!
!  Time conversion routines.

    module time_module

    use kind_module

    implicit none

    private

    !parameters:
    real(wp),parameter :: jd_j2000 = 2451545.0_wp  !! julian date of J2000 epoch

    interface julian_date
        !! calendar date to julian date
        module procedure :: julian_date_realsec, &
                            julian_date_intsec
    end interface

    interface julian_date_to_calendar_date
        module procedure :: calendar_date_realsec
    end interface

    !public routines:
    public :: julian_day
    public :: julian_date
    public :: et_to_jd
    public :: jd_to_et
    public :: jd_to_mjd
    public :: mjd_to_jd
    public :: julian_date_to_calendar_date
    public :: calendar_date_to_et
    public :: is_leap_year

    !test routine:
    public :: time_module_test

    contains
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 2/3/2015
!
!  Convert ephemeris time (seconds from J2000 epoch) to Julian date.

    pure function et_to_jd(et) result(jd)

    use conversion_module, only: sec2day

    implicit none

    real(wp),intent(in) :: et   !! ephemeris time [sec from J2000 epoch]
    real(wp)            :: jd   !! Julian date [days]

    jd = jd_j2000 + et*sec2day

    end function et_to_jd
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 3/19/2016
!
!  Convert Julian date to ephemeris time (seconds from J2000 epoch).

    pure function jd_to_et(jd) result(et)

    use conversion_module, only: day2sec

    implicit none

    real(wp),intent(in) :: jd   !! Julian date [days]
    real(wp)            :: et   !! ephemeris time [sec from J2000 epoch]

    et = (jd - jd_j2000) * day2sec

    end function jd_to_et
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 3/15/2015
!
!  Converts Julian date to Modified Julian date.
!
!### Reference
!   * [USNO](http://tycho.usno.navy.mil/mjd.html)

    pure function jd_to_mjd(jd) result(mjd)

    implicit none

    real(wp)            :: mjd  !! modified julian date
    real(wp),intent(in) :: jd   !! julian date

    mjd = jd - 2400000.5_wp

    end function jd_to_mjd
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 3/15/2015
!
!  Converts Modified Julian date to Julian date.
!
!### Reference
!   * [USNO](http://tycho.usno.navy.mil/mjd.html)

    pure function mjd_to_jd(mjd) result(jd)

    implicit none

    real(wp)            :: jd   !! julian date
    real(wp),intent(in) :: mjd  !! modified julian date

    jd = mjd + 2400000.5_wp

    end function mjd_to_jd
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!
!  Returns the Julian day number (i.e., the Julian date at Greenwich noon)
!  on the specified YEAR, MONTH, and DAY.
!
!  Valid for any Gregorian calendar date producing a
!  Julian date greater than zero.
!
!### Reference
!   * [USNO](http://aa.usno.navy.mil/faq/docs/JD_Formula.php)

    pure integer function julian_day(y,m,d)

    implicit none

    integer,intent(in) :: y   !! year (YYYY)
    integer,intent(in) :: m   !! month (MM)
    integer,intent(in) :: d   !! day (DD)

    julian_day = d-32075+1461*(y+4800+(m-14)/12)/4+367*&
                 (m-2-(m-14)/12*12)/12-3*((y+4900+(m-14)/12)/100)/4

    end function julian_day
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 1/21/2015
!
!  Returns the Julian date for the specified YEAR, MONTH, DAY, HR, MIN, SEC.
!
!  Valid for any Gregorian calendar date producing a
!  Julian date greater than zero.
!
!### History
!  * JW : 10/4/2017 : moved main code to [[julian_date_realsec]] routine.

    pure function julian_date_intsec(y,m,d,hour,minute,second) result(julian_date)

    implicit none

    real(wp)           :: julian_date
    integer,intent(in) :: y
    integer,intent(in) :: m
    integer,intent(in) :: d
    integer,intent(in) :: hour
    integer,intent(in) :: minute
    integer,intent(in) :: second

    ! call the other routine:
    julian_date = julian_date_realsec(y,m,d,hour,minute,real(second,wp))

    end function julian_date_intsec
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 1/21/2015
!
!  Returns the Julian date for the specified YEAR, MONTH, DAY, HR, MIN, SEC.
!
!  Valid for any Gregorian calendar date producing a
!  Julian date greater than zero.
!
!### History
!  * JW : 10/4/2017 : made `second` a real value & renamed routine.

    pure function julian_date_realsec(y,m,d,hour,minute,second) result(julian_date)

    implicit none

    real(wp)            :: julian_date
    integer,intent(in)  :: y
    integer,intent(in)  :: m
    integer,intent(in)  :: d
    integer,intent(in)  :: hour
    integer,intent(in)  :: minute
    real(wp),intent(in) :: second

    integer :: julian_day_number

    julian_day_number = julian_day(y,m,d)

    julian_date = real(julian_day_number,wp) + &
                    (hour-12.0_wp)/24.0_wp + &
                    minute/1440.0_wp + &
                    second/86400.0_wp

    end function julian_date_realsec
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!
!  Returns the year, month, day, hr, min, sec for the specified Julian date.
!
!### See also
!  * https://aa.usno.navy.mil/faq/JD_formula.html
!  * http://www.davidgsimpson.com/software/jd2greg_f90.txt

    pure subroutine calendar_date_realsec(julian_date,year,month,day,hrs,min,sec)

    implicit none

    real(wp),intent(in)  :: julian_date !! julian date
    integer,intent(out)  :: year
    integer,intent(out)  :: month
    integer,intent(out)  :: day
    integer,intent(out)  :: hrs
    integer,intent(out)  :: min
    real(wp),intent(out) :: sec

    integer :: i,j,k,l,n,jd
    real(wp) :: frac_day

    jd = int(julian_date)

    l = jd+68569
    n = 4*l/146097
    l = l-(146097*n+3)/4
    i = 4000*(l+1)/1461001
    l = l-1461*i/4+31
    j = 80*l/2447
    k = l-2447*j/80
    l = j/11
    j = j+2-12*l
    i = 100*(n-49)+i+l

    year  = i
    month = j
    day   = k

    frac_day = julian_date - real(jd,wp) + 0.5_wp

    hrs = int(frac_day*24.0_wp)
    min = int((frac_day - hrs/24.0_wp) * 1440.0_wp)
    sec = (frac_day - hrs/24.0_wp - min/1440.0_wp) * 86400.0_wp

    if (sec == 60.0_wp) then
        sec = 0.0_wp
        min = min + 1
     end if

     if (min == 60) then
        min = 0
        hrs = hrs + 1
     end if

    end subroutine calendar_date_realsec
!*****************************************************************************************

!*****************************************************************************************
!>
!  Return true if the specified year is a leap year.

    pure logical function is_leap_year(y)
    integer, intent(in) :: y !! year
    is_leap_year = (mod(y, 4) == 0 .and. (mod(y, 100) /= 0 .or. mod(y, 400) == 0))
    end function is_leap_year
!*****************************************************************************************

!*****************************************************************************************
!>
!  Directly converts a calendar date to seconds since the J2000 epoch.

    function calendar_date_to_et(year, month, day, hour, minute, second) result(et)

    use conversion_module, only: hr2sec, min2sec, day2sec

    integer, intent(in) :: year, month, day, hour, minute
    real(wp), intent(in) :: second
    real(wp) :: et

    ! Constants for the J2000 epoch: January 1, 2000, 12:00:00
    integer, parameter :: j2000_year = 2000
    integer, parameter :: j2000_month = 1
    integer, parameter :: j2000_day = 1
    integer, parameter :: j2000_hour = 12
    integer, parameter :: j2000_minute = 0
    real(wp), parameter :: j2000_second = 0.0_wp

    ! Days in each month (non-leap year)
    integer, dimension(12), parameter :: days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    integer :: y, m
    integer :: total_days

    ! first do time of day:
    et = real(hour - j2000_hour, wp) * hr2sec + &
         real(minute - j2000_minute, wp) * min2sec + &
         (second - j2000_second)

    ! now, calculate the total number of days from J2000 to the given date
    total_days = 0
    do y = j2000_year, year - 1
        total_days = total_days + 365 + merge(1, 0, is_leap_year(y))
    end do
    do m = 1, month - 1
        if (m == 2) then
            total_days = total_days + days_in_month(m) + merge(1, 0, is_leap_year(year))
        else
            total_days = total_days + days_in_month(m)
        end if
    end do

    ! add days:
    et = et + real(total_days + (day - j2000_day), wp) * day2sec

    end function calendar_date_to_et
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 1/21/2015
!
!  Test routine for the Julian date routines.

    subroutine time_module_test()

    implicit none

    real(wp) :: jd, sec
    integer :: year,month,day,hrs,min
    real(wp) :: et1, et2

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' time_module_test'
    write(*,*) '---------------'
    write(*,*) ''

    ! JD = 2451545.0
    jd = julian_date(2000,1,1,12,0,0)

    call calendar_date_realsec(jd,year,month,day,hrs,min,sec)

    write(*,*) 'jd    ', jd
    write(*,*) 'year  ', year
    write(*,*) 'month ', month
    write(*,*) 'day   ', day
    write(*,*) 'hrs   ', hrs
    write(*,*) 'min   ', min
    write(*,*) 'sec   ', sec

    if (year/=2000)  error stop 'error: incorrect year'
    if (month/=1)    error stop 'error: incorrect month'
    if (day/=1)      error stop 'error: incorrect day'
    if (hrs/=12)     error stop 'error: incorrect hrs'
    if (min/=0)      error stop 'error: incorrect min'
    if (sec/=0.0_wp) error stop 'error: incorrect sec'

    ! compare two ways to convert calendar date to ephemeris time:
    ! et1 should be more accurate since it doesn't go through jd
    write(*,*) ''
    et1 = calendar_date_to_et(2026,3,4,8,9,10.12345678_wp)
    et2 = jd_to_et(julian_date(2026,3,4,8,9,10.12345678_wp))
    write(*,*) 'et1   ', et1
    write(*,*) 'et2   ', et2
    write(*,*) 'et diff: ', abs(et1-et2)

    write(*,*) ''
    write(*,*) 'PASSED'

    end subroutine time_module_test
!*****************************************************************************************

!*****************************************************************************************
    end module time_module
!*****************************************************************************************
