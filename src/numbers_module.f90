!*****************************************************************************************
    module numbers_module
!*****************************************************************************************
!****h* FAT/numbers_module
!
!  NAME
!    numbers_module
!
!  DESCRIPTION
!    Defines some numeric parameters.
!
!*****************************************************************************************    
    
    use kind_module,  only: wp
    
    private
    
    real(wp),parameter,public :: zero       = 0.0_wp
    real(wp),parameter,public :: one        = 1.0_wp
    real(wp),parameter,public :: two        = 2.0_wp
    real(wp),parameter,public :: three      = 3.0_wp
    real(wp),parameter,public :: four       = 4.0_wp
    real(wp),parameter,public :: five       = 5.0_wp
    real(wp),parameter,public :: six        = 6.0_wp
    real(wp),parameter,public :: seven      = 7.0_wp
    real(wp),parameter,public :: eight      = 8.0_wp
    real(wp),parameter,public :: nine       = 9.0_wp
    real(wp),parameter,public :: ten        = 10.0_wp
   
    real(wp),parameter,public :: pi         = acos(-one)
    real(wp),parameter,public :: twopi      = two*pi    
    real(wp),parameter,public :: fourpi     = four*pi
    
!*****************************************************************************************    
    end module numbers_module
!*****************************************************************************************    