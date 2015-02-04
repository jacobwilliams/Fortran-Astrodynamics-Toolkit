!*****************************************************************************************
    module iau_orientation_module
!*****************************************************************************************
!****h* FAT/iau_orientation_module
!
!  NAME
!    iau_orientation_module
!
!  DESCRIPTION
!    IAU orientation models for the Earth and planets.
!
!    !!!! NOT FINISHED !!!!
!
!  SEE ALSO
!    [1] NAIF SPICE pck.html documentation
!    [2] "Report of the IAU Working Group on Cartographic Coordinates 
!        and Rotational Elements: 2009", Celestial Mechanics and 
!        Dynamical Astronomy, February 2011, Vol 109, Issue 2, p 101-135.
!
!*****************************************************************************************    
    
    use kind_module
    use numbers_module
    use conversion_module

    implicit none
    
    private
    
    !
    !  TO DO: use a class to provide access to the different models...
    !

    public :: iau_rotation_matrix  !base routine
    public :: j2000_to_iau_earth   
    
    !test routines:
    public :: iau_test

    contains
!*****************************************************************************************

!*****************************************************************************************
!****f* iau_orientation_module/iau_rotation_matrix
!
!  NAME
!    iau_rotation_matrix
!
!  DESCRIPTION
!    Returns the rotation matrix for a coordinate transformation from 
!    the J2000 frame to the IAU rotating frame associated with a body.
!    The IAU orientation models use three Euler angles to describe 
!    the pole and prime meridian location (ra, dec, and w). 
!
!  AUTHOR
!    Jacob Williams, 2/3/2015
!
!  SOURCE

    pure function iau_rotation_matrix(w,dec,ra) result(rotmat)
    
    use vector_module, only: rotation_matrix,x_axis,y_axis,z_axis
    
    implicit none
    
    real(wp),intent(in)     :: w         !right ascension of the pole [rad]
    real(wp),intent(in)     :: dec       !declination of the pole [rad]
    real(wp),intent(in)     :: ra        !prime meridian location [rad]
    real(wp),dimension(3,3) :: rotmat    !the rotation matrix
    
    real(wp),parameter :: pi2 = pi / two
    
    real(wp),dimension(3,3) :: tmp
    
    !it is a 3-1-3 rotation:
    tmp    = matmul( rotation_matrix(x_axis,pi2-dec), rotation_matrix(z_axis,pi2+ra) )
    rotmat = matmul( rotation_matrix(z_axis,w), tmp )
    
    end function iau_rotation_matrix
!*****************************************************************************************

!*****************************************************************************************
!****f* iau_orientation_module/j2000_to_iau_earth
!
!  NAME
!    j2000_to_iau_earth
!
!  DESCRIPTION
!    Rotation matrix from J2000 to IAU_EARTH.
!
!  AUTHOR
!    Jacob Williams, 2/3/2015
!
!  SOURCE

    pure function j2000_to_iau_earth(et) result(rotmat)
    
    implicit none
 
    real(wp),intent(in)     :: et        !ephemeris time [sec]
    real(wp),dimension(3,3) :: rotmat    !the rotation matrix
    
    real(wp) :: w,dec,ra,d,t
    
    d = et * sec2day     !interval in days from the J2000 epoch
    t = et * sec2century !interval in Julian centuries (36525 days) from the J2000 epoch
    
    ra  = (         - 0.641_wp * t          ) * deg2rad
    dec = ( 90.0_wp - 0.557_wp * t          ) * deg2rad
    w   = ( 190.147_wp + 360.9856235_wp * d ) * deg2rad 
    
    rotmat = iau_rotation_matrix(w,dec,ra)
    
    end function j2000_to_iau_earth
!*****************************************************************************************

!*****************************************************************************************
!****f* iau_orientation_module/iau_test
!
!  NAME
!    iau_test
!
!  DESCRIPTION
!    Unit test routine for iau_module.
!
!  AUTHOR
!    Jacob Williams, 2/3/2015
!
!  SOURCE

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
    rotmat = j2000_to_iau_earth(et)  !rotation matrix from J2000 to IAU_EARTH
    
    write(*,*) 'et    =',et
    write(*,*) 'rotmat='
    do i=1,3
        write(*,*) rotmat(i,:)
    end do
    ! Note: agrees with SPICE PXFORM( 'J2000', 'IAU_EARTH', ET, ROTMAT )
    ! rotmat=
    !  0.24742305587604752      -0.96890754534215406       -7.6219971891330182E-010
    !  0.96890754534215406       0.24742305587604752       -2.9847705387722482E-009
    !  3.0805524797727912E-009  -1.0920940632167532E-017    1.0000000000000000     

    end subroutine iau_test
!*****************************************************************************************

!*****************************************************************************************
    end module iau_orientation_module
!*****************************************************************************************