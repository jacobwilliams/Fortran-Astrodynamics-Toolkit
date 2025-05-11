!*****************************************************************************************
!> author: Jacob Williams
!
!  Define the numeric kinds.

    module kind_module

    use, intrinsic :: iso_fortran_env,    only: real32,real64,real128

    implicit none

    private

#ifdef REAL32
    integer,parameter,public :: wp = real32   !! real kind used by this module [4 bytes]
#elif REAL64
    integer,parameter,public :: wp = real64   !! real kind used by this module [8 bytes]
#elif REAL128
    integer,parameter,public :: wp = real128  !! real kind used by this module [16 bytes]
#else
    integer,parameter,public :: wp = real64   !! real kind used by this module [8 bytes]
#endif

    end module kind_module
!*****************************************************************************************
