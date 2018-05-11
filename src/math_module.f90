!*****************************************************************************************
!> author: Jacob Williams
!
!  General math routines

    module math_module

    use kind_module,    only: wp
    use numbers_module, only: pi,twopi

    implicit none

    private

    public :: cube_root
    public :: wrap_angle

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Wrap an angle (in rad) from -pi to pi.

    pure elemental function wrap_angle(a) result(r)

    implicit none

    real(wp),intent(in) :: a
    real(wp)            :: r

    r = mod(a,twopi)
    if (abs(r)>pi) r = r - sign(twopi,r)

    end function wrap_angle
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
