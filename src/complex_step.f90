!*****************************************************************************************
!> author: Jacob Williams
!
!  Complex-step differentiation routines.
!
!# See also
!  1. J.R.R.A. Martins, P. Sturdza, J.J. Alonso,
!     "The Complex-Step Derivative Approximation",
!     ACM Transactions on Mathematical Software,
!     Vol. 29, No. 3, September 2003, Pages 245262.

    module complex_step_module

    use kind_module,    only: wp

    implicit none

    private

    interface
        function func(x) result(f)
        import :: wp
        implicit none
        complex(wp),intent(in) :: x
        complex(wp) :: f
        end function func
    end interface

    public :: complex_step_derivative
    public :: complex_step_test  !for testing

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute the first derivative using the complex-step method.
!  This is Equation 6 from Reference [1].

    subroutine complex_step_derivative(f,x,h,dfdx)

    implicit none

    procedure(func)        :: f
    complex(wp),intent(in) :: x
    real(wp),intent(in)    :: h
    real(wp),intent(out)   :: dfdx

    dfdx = AIMAG(f(cmplx(real(x,wp),AIMAG(x)+h,wp)))/h

    end subroutine complex_step_derivative
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute the first derivative using a forward difference.
!  This is Equation 1 from Reference [1].

    subroutine forward_diff(f,x,h,dfdx)

    implicit none

    procedure(func)        :: f
    complex(wp),intent(in) :: x
    real(wp),intent(in)    :: h
    real(wp),intent(out)   :: dfdx

    dfdx = (f(x+h) - f(x)) / h

    end subroutine forward_diff
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute the first derivative using a 2-point central difference [-h,h].

    subroutine central_diff(f,x,h,dfdx)

    implicit none

    procedure(func)        :: f
    complex(wp),intent(in) :: x
    real(wp),intent(in)    :: h
    real(wp),intent(out)   :: dfdx

    dfdx = (f(x+h) - f(x-h)) / (2.0_wp*h)

    end subroutine central_diff
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute the first derivative using a 4-point central difference [-2h,-h,h,2h].

    subroutine central_diff_4(f,x,h,dfdx)

    implicit none

    procedure(func)        :: f
    complex(wp),intent(in) :: x
    real(wp),intent(in)    :: h
    real(wp),intent(out)   :: dfdx

    real(wp) :: h2

    h2 = 2.0_wp * h

    dfdx = (f(x-h2) - 8.0_wp*f(x-h) + 8.0_wp*f(x+h) - f(x+h2)) / (12.0_wp*h)

    end subroutine central_diff_4
!*****************************************************************************************

!*****************************************************************************************
!>
!  Unit test for the complex_step module.

    subroutine complex_step_test()

    implicit none

    integer     :: i
    complex(wp) :: x
    real(wp)    :: h,dfdx,dfdx2,dfdx3,dfdx4,err,err2,err3,err4

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' complex_step_test'
    write(*,*) '---------------'
    write(*,*) ''

    x = cmplx(2.0_wp,0.0_wp,wp)
    h = 1.0e-10_wp

    call complex_step_derivative(test_func,x,h,dfdx)

    !write(*,*) ''
    !write(*,*) 'x     :',x
    !write(*,*) 'dfdx  :',dfdx
    !write(*,*) 'error :',real(test_deriv(x),wp) - dfdx
    !write(*,*) ''

    write(*,'(*(A30))') 'h', 'forward diff err', 'central diff err', 'central diff 4 err', 'complex step err'

    do i=1,200
        h = 10.0_wp**(-i/10.0_wp)
        call complex_step_derivative(test_func,x,h,dfdx)
        call forward_diff(test_func,x,h,dfdx2)
        call central_diff(test_func,x,h,dfdx3)
        call central_diff_4(test_func,x,h,dfdx4)
        err  = real(test_deriv(x),wp) - dfdx
        err2 = real(test_deriv(x),wp) - dfdx2
        err3 = real(test_deriv(x),wp) - dfdx3
        err4 = real(test_deriv(x),wp) - dfdx4
        write(*,'(*(E30.16,1X))') h, err2, err3, err4, err
    end do

    contains
!*****************************************************************************************

    !****************************************
        function test_func(x) result(f)

        implicit none

        complex(wp),intent(in) :: x
        complex(wp) :: f

        f = exp(x) + sin(x)

        end function test_func
    !****************************************

    !****************************************
        function test_deriv(x) result(f)

        implicit none

        complex(wp),intent(in) :: x
        complex(wp) :: f

        f = exp(x) + cos(x)

        end function test_deriv
    !****************************************

    end subroutine complex_step_test
!*****************************************************************************************

!*****************************************************************************************
    end module complex_step_module
!*****************************************************************************************
