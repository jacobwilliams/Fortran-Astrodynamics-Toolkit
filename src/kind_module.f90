!*****************************************************************************************
!> author: Jacob Williams
!
!  Define the numeric kinds.

    module kind_module

    use, intrinsic :: iso_fortran_env,    only: real32,real64,real128

    implicit none

    private

    !integer,parameter,public :: wp = real32     !! single precision reals
    integer,parameter,public :: wp = real64      !! double precision reals
    !integer,parameter,public :: wp = real128    !! quad precision reals

    end module kind_module
!*****************************************************************************************
