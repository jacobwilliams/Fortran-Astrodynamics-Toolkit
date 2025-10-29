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
    public :: replace_char
    public :: reverse
    public :: lchop, rchop
    public :: strip

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Convert the string to uppercase.

    pure subroutine uppercase(str)

    character(len=*),intent(inout) :: str

    integer :: i,idx

    do i=1,len_trim(str)
        idx = index(lower,str(i:i))
        if (idx>0) str(i:i) = upper(idx:idx)
    end do

    end subroutine uppercase
!*****************************************************************************************

!*****************************************************************************************
!>
!  Convert the string to lowercase.

    pure subroutine lowercase(str)

    character(len=*),intent(inout) :: str

    integer :: i,idx

    do i=1,len_trim(str)
        idx = index(upper,str(i:i))
        if (idx>0) str(i:i) = lower(idx:idx)
    end do

    end subroutine lowercase
!*****************************************************************************************

!*****************************************************************************************
!>
!  Replace all occurrences of a single character `s1` in `str` with `s2`.

    pure function replace_char(str, s1, s2) result(newstr)

    character(len=*),intent(in) :: str  !! original string
    character(len=1),intent(in) :: s1   !! find all occurrences of this character
    character(len=1),intent(in) :: s2   !! replace with this character
    character(len=:),allocatable :: newstr !! new string

    integer :: i !! counter

    newstr = str
    do i = 1, len(newstr)
        if (newstr(i:i) == s1) newstr(i:i) = s2
    end do

    end function replace_char
!*****************************************************************************************

!*****************************************************************************************
!>
!  Chop leading `chars` string from `str`.
!  Note that trailing spaces are not ignored in either string.

    pure function lchop(str, chars) result(newstr)

    character(len=*),intent(in) :: str  !! original string
    character(len=*),intent(in) :: chars !! characters to strip
    character(len=:),allocatable :: newstr !! new string

    ! this logic here is to account for trailing spaces, which we preserve
    if (len(chars)>len(str)) then
        newstr = str  ! not possible to chop
    else
        if (str==chars) then
            if (len(str)>len(chars)) then
                newstr = str(len(chars)+1:)  ! only trailing spaces remain
            else
                newstr = ''  ! string is now empty
            end if
        else
            if (index(str,chars) == 1) then
                newstr = str(len(chars)+1:)  ! chop leading chars, keep rest of string
            else
                newstr = str  ! original string, noting to chop
            end if
        end if
    end if

    end function lchop
!*****************************************************************************************

!*****************************************************************************************
!>
!  Chop trailing `chars` string from `str`.
!  Note that trailing spaces are not ignored in either string.

    pure function rchop(str, chars) result(newstr)

    character(len=*),intent(in) :: str  !! original string
    character(len=*),intent(in) :: chars !! characters to strip
    character(len=:),allocatable :: newstr !! new string

    newstr = reverse(lchop(reverse(str), reverse(chars)))

    end function rchop
!*****************************************************************************************

!*****************************************************************************************
!>
!  Reverse a string.

    pure function reverse(str) result(newstr)

    character(len=*),intent(in) :: str  !! original string
    character(len=:),allocatable :: newstr !! new string
    integer :: i, j !! counter
    integer :: n !! length of `str`

    n = len(str)
    allocate(character(len=n) :: newstr)
    if (n==0) then
        newstr = ''
    else
        j = 0
        do i = n, 1, -1
            j = j + 1
            newstr(j:j) = str(i:i)
        end do
    end if

    end function reverse
!*****************************************************************************************

!*****************************************************************************************
!>
!  Strip all occurances of chars from the beginning and end of the string.

    pure function strip(str, chars) result(newstr)

    character(len=*),intent(in) :: str  !! original string
    character(len=*),intent(in),optional :: chars !! characters to strip
    character(len=:),allocatable :: newstr !! new string

    integer :: i !! counter
    integer :: n !! length of str
    integer :: idx !! for using scan
    integer :: start_idx, end_idx  !! substring to keep

    if (present(chars)) then
        if (chars /= '') then
            ! have to step through manually from beginning and end
            n = len(str)
            start_idx = 1
            end_idx = n
            ! forward search
            do i = 1, n
                idx = scan(str(i:i), chars)
                if (idx > 0) then
                    start_idx = start_idx + 1
                else
                    exit
                end if
            end do
            ! backward search
            do i = n, 1, -1
                idx = scan(str(i:i), chars)
                if (idx > 0) then
                    end_idx = end_idx - 1
                else
                    exit
                end if
            end do
            if (end_idx <= start_idx) then
                newstr = ''  ! all stripped
            else
                newstr = str(start_idx:end_idx)  ! return substring
            end if
            return ! done
        end if
    end if

    ! in this case assuming it's a space, so use intrinsic functions
    newstr = trim(adjustl(str))

    end function strip
!*****************************************************************************************

!*****************************************************************************************
    end module string_module
!*****************************************************************************************