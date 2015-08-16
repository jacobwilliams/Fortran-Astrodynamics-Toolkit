!*****************************************************************************************
!> author: Jacob Williams
!
!  General math routines

    module math_module

    use kind_module,    only: wp
    
    implicit none
    
    private
    
    public :: cube_root
    
    contains
!*****************************************************************************************
    
!*****************************************************************************************
!> author: Jacob Williams
!
!  Cube root of a number (real solution only).

    pure elemental function cube_root(x) result(y)
    
    use numbers_module, only: one,three
    
    implicit none

    real(wp),intent(in) :: x
    real(wp)            :: y

    real(wp),parameter :: one_third = one/three

    y = sign( abs(x) ** one_third , x )
    
    end function cube_root
!*****************************************************************************************

!*****************************************************************************************
    end module math_module
!*****************************************************************************************