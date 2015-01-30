!*****************************************************************************************
    module random_module
!*****************************************************************************************
!****h* FAT/random_module
!
!  NAME
!    random_module
!
!  DESCRIPTION
!    Random number generation.
!
!*****************************************************************************************    
    
    use kind_module,      only: wp
    
    implicit none
    
    private
    
    !public routines:
    public :: get_random_number
    
    contains
!*****************************************************************************************

!*****************************************************************************************
!****f* random_module/get_random_number
!
!  NAME
!    get_random_number
!
!  DESCRIPTION
!    Returns a uniform random number x, such that: a <= x < b.
!
!  SOURCE

    function get_random_number(a,b) result(x)

    implicit none

    real(wp)            :: x
    real(wp),intent(in) :: a
    real(wp),intent(in) :: b

    call random_number(x)
    
    x = a + (b-a)*x
    
    end function get_random_number
!*****************************************************************************************

!*****************************************************************************************
    end module random_module
!*****************************************************************************************