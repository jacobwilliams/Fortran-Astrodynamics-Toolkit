!*****************************************************************************************
!> author: Jacob Williams
!
!  Various matrix routines

    module matrix_module

    use kind_module, only: wp

    implicit none

    private

    public :: print_matrix

    contains
!*****************************************************************************************

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
    end module matrix_module
!*****************************************************************************************
