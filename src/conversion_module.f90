!*****************************************************************************************
    module conversion_module
!*****************************************************************************************
!****h* FAT/conversion_module
!
!  NAME
!    conversion_module
!
!  DESCRIPTION
!    Conversion factors.
!
!  SEE ALSO
!    [1] A. Thompson and B. N. Taylor, "NIST Special Publication 811: 
!            Guide for the use of the International System of Units".
!            http://www.nist.gov/pml/pubs/sp811/
!
!*****************************************************************************************

    use kind_module,       only: wp
    use numbers_module,    only: one
 
    implicit none
  
    real(wp),parameter :: lbm2kg  = 0.45359237_wp         !these are exact
    real(wp),parameter :: lbf2N   = 4.4482216152605_wp    !
    real(wp),parameter :: ft2m    = 0.3048_wp             !
    real(wp),parameter :: mile2km = 1.609344_wp           !
    real(wp),parameter :: nmi2km  = 1.852_wp              !
    real(wp),parameter :: slug2kg = lbf2N/ft2m     ! approximately 14.593902937206362
 
    real(wp),parameter :: kg2lbm  = one/lbm2kg     ! approximately 2.2046226218487757
    real(wp),parameter :: N2lbf   = one/lbf2N      ! approximately 0.2248089430997105
    real(wp),parameter :: m2ft    = one/ft2m       ! approximately 3.280839895013123
    real(wp),parameter :: km2mile = one/mile2km    ! approximately 0.621371192237334
    real(wp),parameter :: km2nmi  = one/nmi2km     ! approximately 0.5399568034557235
    real(wp),parameter :: kg2slug = one/slug2kg    ! approximately 0.06852176585679176
 
!*****************************************************************************************
    end module conversion_module
!*****************************************************************************************