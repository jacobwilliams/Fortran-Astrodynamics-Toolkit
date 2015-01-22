!*****************************************************************************************
    module time_module
!*****************************************************************************************
!****h* FAT/time_module
!
!  NAME
!    time_module
!
!  DESCRIPTION
!    Time conversion routines.
!
!*****************************************************************************************    
    
    use kind_module

    implicit none
    
    private
    
    public :: julian_day
    public :: julian_date
    
    public :: time_module_test

    contains
!*****************************************************************************************

!*****************************************************************************************
!****f* time_module/julian_day
!
!  NAME    
!    julian_day
!
!  DESCRIPTION
!    Returns the Julian day number (i.e., the Julian date at Greenwich noon) 
!        on the specified YEAR, MONTH, and DAY.
!    Valid for any Gregorian calendar date producing a 
!        Julian date greater than zero.
!
!  NOTES
!    http://aa.usno.navy.mil/faq/docs/JD_Formula.php
!
!  SOURCE

    pure integer function julian_day(y,m,d)
    
    implicit none
    
    integer,intent(in) :: y        ! year (YYYY)
    integer,intent(in) :: m        ! month (MM)
    integer,intent(in) :: d        ! day (DD)
    
    julian_day = d-32075+1461*(y+4800+(m-14)/12)/4+367*&
                 (m-2-(m-14)/12*12)/12-3*((y+4900+(m-14)/12)/100)/4
    
    end function julian_day
!*****************************************************************************************
    
!*****************************************************************************************
!****f* time_module/julian_date
!
!  NAME
!    julian_date
!
!  DESCRIPTION
!    Returns the Julian date for the specified YEAR, MONTH, DAY, HR, MIN, SEC.
!    Valid for any Gregorian calendar date producing a 
!        Julian date greater than zero.
!
!  AUTHOR
!    Jacob Williams : 1/21/2015
!
!  SOURCE

    function julian_date(y,m,d,hour,minute,second)
   
    implicit none
    
    real(wp) :: julian_date
    integer,intent(in) :: y,m,d,hour,minute,second
    
    integer :: julian_day_number
    
    julian_day_number = julian_day(y,m,d)
    
    julian_date = real(julian_day_number,wp) + &
                    (hour-12.0_wp)/24.0_wp + &
                    minute/1440.0_wp + &
                    second/86400.0_wp
                                
    end function julian_date
!*****************************************************************************************
    
!*****************************************************************************************
!****f* time_module/time_module_test
!
!  NAME
!    time_module_test
!
!  DESCRIPTION
!    Test routine for the Julian date routines.
!
!  AUTHOR
!    Jacob Williams : 1/21/2015
!
!  SOURCE

    subroutine time_module_test()

    implicit none
    
    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' time_module_test'
    write(*,*) '---------------'
    write(*,*) ''

    ! JD = 2451545.0
    write(*,*) julian_date(2000,1,1,12,0,0)

    end subroutine time_module_test
!*****************************************************************************************
       
!*****************************************************************************************
    end module time_module
!*****************************************************************************************