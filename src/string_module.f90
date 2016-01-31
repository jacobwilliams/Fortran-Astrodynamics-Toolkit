!*****************************************************************************************
!> author: Jacob Williams
!
!  Module for string manipulation.

    module string_module

    implicit none

    private

    character(len=*),parameter :: lower = 'abcdefghijklmnopqrstuvwxyz' !! lowercase characters
    character(len=*),parameter :: upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' !! uppercase characters

    public :: lowercase,uppercase

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Convert the string to uppercase.

    subroutine uppercase(str)

    implicit none

    character(len=*),intent(inout) :: str

    integer :: i,idx

    do i=1,len_trim(str)
        idx = index(str(i:i),lower)
        if (idx>0) str(i:i) = upper(idx:idx)
    end do

    end subroutine uppercase
!*****************************************************************************************

!*****************************************************************************************
!>
!  Convert the string to lowercase.

    subroutine lowercase(str)

    implicit none

    character(len=*),intent(inout) :: str

    integer :: i,idx

    do i=1,len_trim(str)
        idx = index(str(i:i),upper)
        if (idx>0) str(i:i) = lower(idx:idx)
    end do

    end subroutine lowercase
!*****************************************************************************************

    end module string_module
!*****************************************************************************************
