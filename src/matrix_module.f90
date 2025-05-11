!*****************************************************************************************
!> author: Jacob Williams
!
!  Various matrix routines

    module matrix_module

    use kind_module,    only: wp
    use numbers_module, only: zero,one

    implicit none

    private

    public :: print_matrix
    public :: matrix_trace
    public :: matrix_determinant
    public :: matrix_cofactor

    contains
!*****************************************************************************************

!*******************************************************************************
!>
!  Matrix determinant of an \( n \times n \) matrix (recursive formulation).
!
!### Reference
!  * https://rosettacode.org/wiki/Matrix_arithmetic#Fortran

    pure recursive function matrix_determinant(n,a) result(det)

    implicit none

    integer,intent(in)                  :: n   !! size of `a` matrix
    real(wp),dimension(n,n),intent(in)  :: a   !! the matrix
    real(wp)                            :: det !! the determinant of `a` matrix

    integer :: i   !! counter
    integer :: sgn !! `(-1)**(i-1)` term
    real(wp), dimension(n-1, n-1) :: b  !! temp matrix

    if (n == 1) then
        det = a(1,1)
    else
        det = zero
        sgn = 1
        do i = 1, n
            b(:, :(i-1)) = a(2:, :i-1)
            b(:, i:) = a(2:, i+1:)
            det = det + sgn * a(1, i) * matrix_determinant(n-1,b)
            sgn = -sgn
        end do
    end if

    end function matrix_determinant
!*******************************************************************************

!*******************************************************************************
!>
!  Compute the cofactors matrix (the transpose of the adjugate matrix).
!
!### References
!  * https://warwick.ac.uk/fac/sci/physics/research/condensedmatt/imr_cdt/students/david_goodwin/teaching/cis008-2/determinant_algorithm_cis008-2_lec_21.pdf
!  * https://groups.google.com/forum/#!topic/comp.lang.fortran/Y6jCv-QdDhc

    pure function matrix_cofactor (n,a) result(c)

    implicit none

    integer,intent(in)                  :: n   !! size of `a` matrix
    real(wp),dimension(n,n),intent(in)  :: a   !! the matrix
    real(wp),dimension(n,n)             :: c   !! the cofactors of `a` matrix

    integer :: i !! counter
    integer :: j !! counter
    real(wp),dimension(n-1,n-1) :: c_temp  !! temp matrix
    logical,dimension(n,n) :: m  !! mask for row/col removal

    if (n==1) then
        c(1,1) = one
    else
        do i=1,n
            do j=1,n
                ! remove the ith row and jth col from a:
                m = .true.
                m(:,j) = .false.
                m(i,:) = .false.
                c_temp = reshape(pack(a, m),[n-1,n-1])
                c(i,j) = ( (-1)**(i+j) ) * matrix_determinant(n-1,c_temp)
            end do
        end do
    end if

    end function matrix_cofactor
!*******************************************************************************

!*****************************************************************************************
!>
!  Print a matrix to the console.

    subroutine print_matrix(mat,unit)

    use iso_fortran_env, only: output_unit

    implicit none

    real(wp),dimension(:,:),intent(in) :: mat !! the matrix to print
    integer,intent(in),optional :: unit  !! unit number (assumed to be an open file).
                                         !! if not present, then the standard output is used.

    integer :: i      !! counter
    integer :: n      !! number of rows in the matrix
    integer :: iunit  !! the file unit to print to
    integer :: istat  !! `iostat` flag for write statement

    character(len=*),parameter :: fmt = 'E26.16'  !! real number format statement

    if (present(unit)) then
        iunit = unit
    else
        iunit = output_unit
    end if

    n = size(mat,1)

    do i=1,n
        write(iunit,fmt='(*('//fmt//',1X))',iostat=istat) mat(i,:)
    end do

    end subroutine print_matrix
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute the matrix trace (sum of the diagonal elements).

    pure function matrix_trace(n,mat) result(trace)

    implicit none

    integer,intent(in)                 :: n     !! size of the matrix
    real(wp),dimension(n,n),intent(in) :: mat   !! the matrix
    real(wp)                           :: trace !! the matrix trace

    integer :: i !! counter

    trace = zero
    do i = 1, n
        trace = trace + mat(i,i)
    end do

    end function matrix_trace
!*****************************************************************************************

!*****************************************************************************************
    end module matrix_module
!*****************************************************************************************
