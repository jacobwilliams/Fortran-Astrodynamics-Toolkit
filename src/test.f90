!*****************************************************************************************
!****h* FAT/test
!
!  NAME
!    test
!
!  DESCRIPTION
!    Unit test program for the Fortran Astrodynamics Toolkit
!
!  SOURCE

    program test
    
    use fortran_astrodynamics_toolkit
    
    implicit none
        
    call bobyqa_test()        
    call brent_test()
    call ephemeris_test()
    
    end program test
!*****************************************************************************************