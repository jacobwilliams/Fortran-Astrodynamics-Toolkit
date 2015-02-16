!*****************************************************************************************
    module math_module
!*****************************************************************************************
!****h* FAT/math_module
!
!  NAME
!    math_module
!
!  DESCRIPTION
!    General math routines
!
!*****************************************************************************************

    use kind_module,    only: wp
    
    implicit none
    
    private
    
    public :: cube_root
    
    contains
!*****************************************************************************************
    
!*****************************************************************************************
!****f* math_module/cube_root
!
!  NAME
!    cube_root
!
!  DESCRIPTION
!    Cube root of a number (real solution only).
!
!  AUTHOR
!    Jacob Williams
!
!  SOURCE

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