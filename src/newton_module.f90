!*******************************************************************************
!>
!  Basic Newton solver.

    module newton_module

    use kind_module,    only: wp
    use numbers_module

    implicit none

    private

    public :: newton

    abstract interface
        subroutine func(x,f)
        !! interface for the function and derivative
        import :: wp
        implicit none
        real(wp),intent(in)  :: x
        real(wp),intent(out) :: f
        end subroutine func
    end interface

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  Newton's method for root finding of scalar function f(x)

    subroutine newton(x,f,dfdx,ftol,xtol,max_iter,xs,fx,iflag)

    implicit none

    real(wp),intent(in)    :: x         !! initial point (initial guess)
    procedure(func)        :: f         !! function f(x)
    procedure(func)        :: dfdx      !! first derivative function f'(x)
    real(wp),intent(in)    :: ftol      !! convergence tolerance for f(x)
    real(wp),intent(in)    :: xtol      !! convergence tolerance for x
    integer,intent(in)     :: max_iter  !! the maximum number of iterations
    real(wp),intent(out)   :: xs        !! the value where f(x) is zero
    real(wp),intent(out)   :: fx        !! the value of f(x) at the root xs
    integer,intent(out)    :: iflag     !! status flag:
                                        !!  0 : absolute convergence in f
                                        !!  1 : relative convergence in x
                                        !! -1 : Error: derivative is zero
                                        !! -2 : Error: max iterations exceeded

    real(wp)  :: f1,df1,x1,x1_prev
    integer   :: iter
    real(wp),parameter :: alpha = one !! step factor

    x1 = x

    do iter = 1,max_iter

        call f(x1, f1)
        if (abs(f1)<=ftol) then
            iflag = 0
            exit
        end if

        call dfdx(x1, df1)
        if (df1==zero) then
            iflag = -1
            exit
        end if

        x1_prev = x1
        x1 = x1 - alpha * ( f1 / df1 )

        if (abs(x1_prev-x1)<=xtol) then
            iflag = 1
            exit
        end if

    end do

    ! max iterations exceeded
    if (iter>max_iter) iflag = -2

    !return results:
    xs  = x1
    fx  = f1

    end subroutine newton
!*******************************************************************************

!*******************************************************************************
    end module newton_module
!*******************************************************************************
