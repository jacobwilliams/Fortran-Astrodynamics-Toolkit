!*****************************************************************************************
!> author: Jacob Williams
!
!  Random number generation.

    module random_module

    use kind_module,      only: wp

    implicit none

    private

    !public routines:
    public :: get_random_number

    contains
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!
!  Returns a uniform random number `x`, such that: `a <= x < b`.

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
