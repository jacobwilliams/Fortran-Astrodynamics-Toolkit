!*****************************************************************************************
    module complex_step_module
!*****************************************************************************************
!****h* FAT/complex_step_module
!
!  NAME
!    complex_step_module
!
!  DESCRIPTION
!    Complex-step differentiation routines.
!
!  SEE ALSO
!    [1] J.R.R.A. Martins, P. Sturdza, J.J. Alonso, 
!        "The Complex-Step Derivative Approximation",
!        ACM Transactions on Mathematical Software, 
!        Vol. 29, No. 3, September 2003, Pages 245–262.
!
!*****************************************************************************************
    
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
!****f* complex_step_module/complex_step_derivative
! 
!  NAME
!    complex_step_derivative
!
!  DESCRIPTION
!    Compute the first derivative using the complex-step method.
!    This is Equation 6 from Reference [1].
!
!  SOURCE

    subroutine complex_step_derivative(f,x,h,dfdx)
    
    implicit none
    
    procedure(func)        :: f
    complex(wp),intent(in) :: x
    real(wp),intent(in)    :: h
    real(wp),intent(out)   :: dfdx
    
    dfdx = imag(f(cmplx(real(x,wp),imag(x)+h,wp)))/h
    
    end subroutine complex_step_derivative
!*****************************************************************************************

!*****************************************************************************************
!****f* complex_step_module/forward_diff
! 
!  NAME
!    forward_diff
!
!  DESCRIPTION
!    Compute the first derivative using a forward difference.
!    This is Equation 1 from Reference [1].
!
!  SOURCE

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
!****f* complex_step_module/complex_step_test
! 
!  NAME
!    complex_step_test
!
!  DESCRIPTION
!    Unit test for the complex_step module.
!
!  SOURCE

    subroutine complex_step_test()
    
    implicit none
    
    integer     :: i
    complex(wp) :: x
    real(wp)    :: dfdx,h,err,dfdx2,err2

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
    
    write(*,'(*(A30))') 'h', 'forward diff err', 'complex step err'
    
    do i=1,50
        h = 10.0_wp**(-i)
        call complex_step_derivative(test_func,x,h,dfdx)
        call forward_diff(test_func,x,h,dfdx2)
        err  = real(test_deriv(x),wp) - dfdx
        err2 = real(test_deriv(x),wp) - dfdx2
        write(*,'(*(E30.16))') h, err2, err
    end do
    
    contains
    
!   contains: test_func, test_deriv
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