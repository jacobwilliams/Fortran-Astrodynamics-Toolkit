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
    public :: magnitude

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Returns a positive number the same magnitude as the input,
!  with only one significant digit.
!
!  If `mina` is present, then `max(mina,mag(a))` is returned
!
!  Examples:
!```
!     mag(1234.56)  -> 1000.0
!     mag(-999.99)  -> 900.0
!     mag(1.456e-4) -> 0.0001
!```

    pure elemental function magnitude(a,mina) result(m)

    implicit none

    real(wp),intent(in) :: a
    real(wp),intent(in),optional :: mina
    real(wp) :: m

    real(wp) :: x,tmp

    x = abs(a)

    if (x==0.0_wp) then
        if (.not. present(mina)) then
            m = 1.0_wp
        else
            m = mina
        end if
    else
        tmp = 10.0_wp ** floor(log10(x))
        m = tmp * floor(x/tmp)
        if (present(mina)) m = max(mina,m)
    end if

    end function magnitude
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
