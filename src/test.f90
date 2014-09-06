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
    call vector_test()
    call lambert_test()
    call rk_test()
    call complex_step_test()
    
    end program test
!*****************************************************************************************