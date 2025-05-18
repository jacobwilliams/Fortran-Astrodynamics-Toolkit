!*****************************************************************************************
!> author: Jacob Williams
!
!  A base class for defining other classes.

    module base_class_module

    implicit none

    integer,parameter :: name_len = 100  !! length of name strings

    type,abstract,public :: base_class
        !! A base class for defining other classes.
        integer :: id = 0  !! a unique ID code that distinguishes a
                           !! variable from other variables of the same type.
        character(len=name_len) :: name = '' !! the variable name
    contains
        generic,public :: operator(==) => base_class_equal
        generic,public :: operator(/=) => base_class_not_equal
        procedure,private :: base_class_equal
        procedure,private :: base_class_not_equal
    end type base_class

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  `==` operator for [[base_class]] variables.
!  To be equal, they must be the same type and have the same `ID`.

    pure elemental function base_class_equal(b1,b2) result(is_equal)

    implicit none

    class(base_class),intent(in) :: b1
    class(base_class),intent(in) :: b2
    logical :: is_equal

    is_equal = same_type_as(b1,b2) .and. (b1%id == b2%id)

    end function base_class_equal
!*****************************************************************************************

!*****************************************************************************************
!>
!  `/=` operator for [[base_class]] variables.
!  To be equal, they must be the same type and have the same `ID`.

    pure elemental function base_class_not_equal(b1,b2) result(not_equal)

    implicit none

    class(base_class),intent(in) :: b1
    class(base_class),intent(in) :: b2
    logical :: not_equal

    not_equal = .not. (b1%id == b2%id)

    end function base_class_not_equal
!*****************************************************************************************

!*****************************************************************************************
    end module base_class_module
!*****************************************************************************************
