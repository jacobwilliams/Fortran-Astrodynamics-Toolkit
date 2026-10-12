!*****************************************************************************************
!> author: Jacob Williams
!
!  IAU orientation models for the Earth and planets.
!
!# See also
!  1. NAIF SPICE pck.html documentation
!  2. "Report of the IAU Working Group on Cartographic Coordinates
!     and Rotational Elements: 2009", Celestial Mechanics and
!     Dynamical Astronomy, February 2011, Vol 109, Issue 2, p 101-135.
!
!@warning Not Finished.
!
!@todo Use a class to provide access to the different models

    module iau_orientation_module

    use kind_module
    use numbers_module
    use conversion_module

    implicit none

    private

    public :: iau_rotation_matrix  !base routine
    public :: icrf_to_iau_earth
    public :: icrf_to_iau_moon

    !test routines:
    public :: iau_test

    contains
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 2/3/2015
!
!  Returns the rotation matrix for a coordinate transformation from
!  the International Celestial Reference Frame (ICRF) frame to the
!  IAU rotating frame associated with a body.
!  The IAU orientation models use three Euler angles to describe
!  the pole and prime meridian location (ra, dec, and w).

    pure function iau_rotation_matrix(w,dec,ra) result(rotmat)

    use vector_module, only: rotation_matrix,x_axis,y_axis,z_axis

    implicit none

    real(wp),intent(in)     :: w         !! right ascension of the pole [rad]
    real(wp),intent(in)     :: dec       !! declination of the pole [rad]
    real(wp),intent(in)     :: ra        !! prime meridian location [rad]
    real(wp),dimension(3,3) :: rotmat    !! the rotation matrix

    real(wp),parameter :: pi2 = pi / two

    real(wp),dimension(3,3) :: tmp

    !it is a 3-1-3 rotation:
    tmp    = matmul( rotation_matrix(x_axis,pi2-dec), rotation_matrix(z_axis,pi2+ra) )
    rotmat = matmul( rotation_matrix(z_axis,w), tmp )

    end function iau_rotation_matrix
!*****************************************************************************************

!
! TO DO:
!... also need to add computation of rotmatdot
! see: ftp://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/tisbod.html
!

!*****************************************************************************************
!> author: Jacob Williams
!  date: 2/3/2015
!
!  Rotation matrix from ICRF to IAU_EARTH.

    pure function icrf_to_iau_earth(et) result(rotmat)

    implicit none

    real(wp),intent(in)     :: et        !! ephemeris time [sec]
    real(wp),dimension(3,3) :: rotmat    !! the rotation matrix

    real(wp) :: w,dec,ra,d,t

    d = et * sec2day     !interval in days from the J2000 epoch
    t = et * sec2century !interval in Julian centuries (36525 days) from the J2000 epoch

    ra  = (         - 0.641_wp * t          ) * deg2rad
    dec = ( 90.0_wp - 0.557_wp * t          ) * deg2rad
    w   = ( 190.147_wp + 360.9856235_wp * d ) * deg2rad

    rotmat = iau_rotation_matrix(w,dec,ra)

    end function icrf_to_iau_earth
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 3/5/2016
!
!  Rotation matrix from ICRF to IAU_MOON.

    pure function icrf_to_iau_moon(et) result(rotmat)

    implicit none

    real(wp),intent(in)     :: et        !! ephemeris time [sec]
    real(wp),dimension(3,3) :: rotmat    !! the rotation matrix

    real(wp) :: w,dec,ra,d,t
    real(wp) :: e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13

    d = et * sec2day     !interval in days from the J2000 epoch
    t = et * sec2century !interval in Julian centuries (36525 days) from the J2000 epoch

    ! Nutation/precession angles:
    e1  = ( 125.045_wp - 0.0529921_wp  * d ) * deg2rad
    e2  = ( 250.089_wp - 0.1059842_wp  * d ) * deg2rad
    e3  = ( 260.008_wp + 13.0120009_wp * d ) * deg2rad
    e4  = ( 176.625_wp + 13.3407154_wp * d ) * deg2rad
    e5  = ( 357.529_wp + 0.9856003_wp  * d ) * deg2rad
    e6  = ( 311.589_wp + 26.4057084_wp * d ) * deg2rad
    e7  = ( 134.963_wp + 13.0649930_wp * d ) * deg2rad
    e8  = ( 276.617_wp + 0.3287146_wp  * d ) * deg2rad
    e9  = ( 34.226_wp  + 1.7484877_wp  * d ) * deg2rad
    e10 = ( 15.134_wp  - 0.1589763_wp  * d ) * deg2rad
    e11 = ( 119.743_wp + 0.0036096_wp  * d ) * deg2rad
    e12 = ( 239.961_wp + 0.1643573_wp  * d ) * deg2rad
    e13 = ( 25.053_wp  + 12.9590088_wp * d ) * deg2rad

    ra  = ( 269.9949_wp + 0.0031_wp * t &
            - 3.8787_wp * sin(E1 ) &
            - 0.1204_wp * sin(E2 ) &
            + 0.0700_wp * sin(E3 ) &
            - 0.0172_wp * sin(E4 ) &
            + 0.0072_wp * sin(E6 ) &
            - 0.0052_wp * sin(E10) &
            + 0.0043_wp * sin(E13) ) * deg2rad

    dec = ( 66.5392_wp + 0.0130_wp * t &
            + 1.5419_wp * cos(E1 ) &
            + 0.0239_wp * cos(E2 ) &
            - 0.0278_wp * cos(E3 ) &
            + 0.0068_wp * cos(E4 ) &
            - 0.0029_wp * cos(E6 ) &
            + 0.0009_wp * cos(E7 ) &
            + 0.0008_wp * cos(E10) &
            - 0.0009_wp * cos(E13) ) * deg2rad

    w = ( 38.3213_wp + 13.17635815_wp * d - 1.4e-12_wp * d**2 &
          + 3.5610_wp * sin(E1 ) &
          + 0.1208_wp * sin(E2 ) &
          - 0.0642_wp * sin(E3 ) &
          + 0.0158_wp * sin(E4 ) &
          + 0.0252_wp * sin(E5 ) &
          - 0.0066_wp * sin(E6 ) &
          - 0.0047_wp * sin(E7 ) &
          - 0.0046_wp * sin(E8 ) &
          + 0.0028_wp * sin(E9 ) &
          + 0.0052_wp * sin(E10) &
          + 0.0040_wp * sin(E11) &
          + 0.0019_wp * sin(E12) &
          - 0.0044_wp * sin(E13) ) * deg2rad

    rotmat = iau_rotation_matrix(w,dec,ra)

    end function icrf_to_iau_moon
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 2/3/2015
!
!  Unit test routine for iau_module.

    subroutine iau_test()

    implicit none

    real(wp)                :: et
    real(wp),dimension(3,3) :: rotmat

    integer :: i

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' iau_test'
    write(*,*) '---------------'
    write(*,*) ''

    et     = 1000.0_wp               !epoch [sec from J2000]
    rotmat = icrf_to_iau_earth(et)   !rotation matrix from J2000 to IAU_EARTH

    write(*,*) '---icrf_to_iau_earth---'
    write(*,*) 'et    =',et
    write(*,*) 'rotmat='
    do i=1,3
        write(*,*) rotmat(i,:)
    end do

    ! ---icrf_to_iau_earth---
    ! et    =   1000.0000000000000
    ! rotmat=
    !  0.24742305587604752      -0.96890754534215406       -7.6219971891330182E-010
    !  0.96890754534215406       0.24742305587604752       -2.9847705387722482E-009
    !  3.0805524797727912E-009  -1.0920940632167532E-017   1.0000000000000000
    !
    ! Compare to SPICE:
    ! call PXFORM( 'J2000', 'IAU_EARTH', ET, ROTMAT )
    ! rotmat=
    !  0.24742305587604752      -0.96890754534215406       -7.6219971891330182E-010
    !  0.96890754534215406       0.24742305587604752       -2.9847705387722482E-009
    !  3.0805524797727912E-009  -1.0920940632167532E-017    1.0000000000000000

    et     = 1000.0_wp               !epoch [sec from J2000]
    rotmat = icrf_to_iau_moon(et)    !rotation matrix from J2000 to IAU_MOON

    write(*,*) '---icrf_to_iau_moon---'
    write(*,*) 'et    =',et
    write(*,*) 'rotmat='
    do i=1,3
        write(*,*) rotmat(i,:)
    end do

    ! ---icrf_to_iau_moon---
    ! et    =   1000.0000000000000
    ! rotmat=
    !  0.78257369718829173       0.55976292119480831       0.27247730276942850
    ! -0.62214729926548507       0.71906872303263469       0.30963350847878057
    ! -2.2608548951909870E-002 -0.41183205753276536       0.91097925876642116
    !
    ! Compare to SPICE:
    ! call PXFORM( 'J2000', 'IAU_MOON', ET, ROTMAT )
    ! rotmat=
    !  0.78257369718829173       0.55976292119480819       0.27247730276942861
    ! -0.62214729926548507       0.71906872303263458       0.30963350847878074
    ! -2.2608548951909880E-002 -0.41183205753276558       0.91097925876642105


    end subroutine iau_test
!*****************************************************************************************

!*****************************************************************************************
    end module iau_orientation_module
!*****************************************************************************************
