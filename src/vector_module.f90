!*****************************************************************************************
    module vector_module
!*****************************************************************************************
!****h* FAT/vector_module
!
!  NAME
!    vector_module
!
!  DESCRIPTION
!    Routines for the manipulation of vectors.
!
!*****************************************************************************************

    use kind_module,    only: wp
    use numbers_module, only: one,zero,pi
    
    implicit none
    
    private
    
    public :: cross
    public :: unit
    public :: ucross
    public :: axis_angle_rotation
    public :: cross_matrix
    public :: axis_angle_rotation_to_rotation_matrix
    
    public :: vector_test
    
    contains
!*****************************************************************************************
    
!*****************************************************************************************
!****f* vector_module/cross
!
!  NAME
!    cross
!
!  DESCRIPTION
!    Cross product of two vectors
!
!  AUTHOR
!    Jacob Williams
!
!  SOURCE

    pure function cross(r,v) result(rxv)

    implicit none

    real(wp),dimension(3),intent(in) :: r
    real(wp),dimension(3),intent(in) :: v
    real(wp),dimension(3)            :: rxv

    rxv = [r(2)*v(3) - v(2)*r(3), &
           r(3)*v(1) - v(3)*r(1), &
           r(1)*v(2) - v(1)*r(2) ]

    end function cross
!*****************************************************************************************

!*****************************************************************************************
!****f* vector_module/unit
!
!  NAME
!    unit
!
!  DESCRIPTION
!    Unit vector
!
!  AUTHOR
!    Jacob Williams
!
!  SOURCE

    pure function unit(r) result(u)

    implicit none

    real(wp),dimension(3),intent(in) :: r
    real(wp),dimension(3)            :: u

    real(wp) :: rmag

    rmag = norm2(r)

    if (rmag==0.0_wp) then
        u = 0.0_wp
    else
        u = r / rmag
    end if

    end function unit
!*****************************************************************************************
    
!*****************************************************************************************
!****f* vector_module/ucross
!
!  NAME
!    ucross
!
!  DESCRIPTION
!    Unit vector of the cross product of two vectors
!
!  AUTHOR
!    Jacob Williams
!
!  SOURCE

    pure function ucross(v1,v2) result(u)

    implicit none

    real(wp),dimension(3),intent(in) :: v1
    real(wp),dimension(3),intent(in) :: v2
    real(wp),dimension(3)            :: u

    u = unit(cross(v1,v2))

    end function ucross
!*****************************************************************************************
    
!*****************************************************************************************
!****f* vector_module/axis_angle_rotation
!
!  NAME
!    axis_angle_rotation
!
!  DESCRIPTION
!    Rotate a vector in space, given an axis and angle of rotation. 
!
!  SEE ALSO
!    http://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
!
!  AUTHOR
!    Jacob Williams, 7/20/2014
!
!  SOURCE

    subroutine axis_angle_rotation(v,k,theta,vrot)
    
    implicit none
    
    real(wp),dimension(3),intent(in)  :: v      !vector to rotate
    real(wp),dimension(3),intent(in)  :: k      !rotation axis
    real(wp),intent(in)               :: theta  !rotation angle [rad]
    real(wp),dimension(3),intent(out) :: vrot   !result
    
    real(wp),dimension(3) :: khat
    real(wp) :: ct,st
    
    ct = cos(theta)
    st = sin(theta)      
    khat = unit(k)   !rotation axis unit vector
    
    vrot = v*ct + cross(khat,v)*st + khat*dot_product(khat,v)*(one-ct)
     
    end subroutine axis_angle_rotation
!*****************************************************************************************

!*****************************************************************************************
!****f* vector_module/cross_matrix
!
!  NAME
!    cross_matrix
!
!  DESCRIPTION
!    Computes the cross product matrix, where:
!      cross(a,b) == matmul(cross_matrix(a),b)
!
!  AUTHOR
!    Jacob Williams, 7/20/2014
!
!  SOURCE

    function cross_matrix(r) result(rcross)
    
    implicit none
    
    real(wp),dimension(3),intent(in) :: r 
    real(wp),dimension(3,3)          :: rcross
    
    rcross(:,1) = [zero,r(3),-r(2)]
    rcross(:,2) = [-r(3),zero,r(1)]
    rcross(:,3) = [r(2),-r(1),zero]
    
    end function cross_matrix
!*****************************************************************************************

!*****************************************************************************************
!****f* vector_module/axis_angle_rotation_to_rotation_matrix
!
!  NAME
!    axis_angle_rotation_to_rotation_matrix
!
!  DESCRIPTION
!    Computes the rotation matrix that corresponds to a 
!      rotation about the axis k by an angle theta.
!
!  AUTHOR
!    Jacob Williams, 7/20/2014
!
!  SOURCE

    subroutine axis_angle_rotation_to_rotation_matrix(k,theta,rotmat)
    
    implicit none
    
    real(wp),dimension(3),intent(in)    :: k        !rotation axis
    real(wp),intent(in)                 :: theta    !rotation angle [rad]
    real(wp),dimension(3,3),intent(out) :: rotmat   !rotation matrix
    
    !3x3 identity matrix:
    real(wp),dimension(3,3),parameter :: I = &
            reshape([one,zero,zero,zero,one,zero,zero,zero,one],[3,3])

    real(wp),dimension(3,3) :: w
    real(wp),dimension(3) :: khat
    real(wp) :: ct,st

    ct = cos(theta)
    st = sin(theta)  
    khat = unit(k)
    w  = cross_matrix(khat)
    
    rotmat = I + w*st + matmul(w,w)*(one-ct)
       
    end subroutine axis_angle_rotation_to_rotation_matrix
!*****************************************************************************************

!*****************************************************************************************
!****f* vector_module/vector_test
!
!  NAME
!    vector_test
!
!  DESCRIPTION
!    Unit test routine for the vector module.
!
!  AUTHOR
!    Jacob Williams, 7/20/2014
!
!  SOURCE

    subroutine vector_test()
    
    implicit none
    
    integer :: i
    real(wp) :: theta
    real(wp),dimension(3) :: v,k,v2,v3
    real(wp),dimension(3,3) :: rotmat

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' vector_test'
    write(*,*) '---------------'
    write(*,*) ''
    
    v = [1.2_wp, 3.0_wp, -5.0_wp]
    k = [-0.1_wp, 16.2_wp, 2.1_wp]
    theta = 0.123_wp
    
    call axis_angle_rotation(v,k,theta,v2)
    
    call axis_angle_rotation_to_rotation_matrix(k,theta,rotmat)
    v3 = matmul(rotmat,v)
    
    write(*,*) 'Single test:'
    write(*,*) ''
    write(*,*) '  v1   :', v
    write(*,*) '  v2   :', v2
    write(*,*) '  v3   :', v3    
    write(*,*) '  Error:', v3-v2
    
    write(*,*) ''
    write(*,*) '0-360 test:'
    write(*,*) ''
    do i=0,360,10
    
        theta = i * 180.0_wp/pi
        
        call axis_angle_rotation(v,k,theta,v2)
    
        call axis_angle_rotation_to_rotation_matrix(k,theta,rotmat)
        v3 = matmul(rotmat,v)
            
        write(*,*) 'Error:', norm2(v3-v2)
        
    end do
    
    end subroutine vector_test
 !*****************************************************************************************
   
!*****************************************************************************************
    end module vector_module
!*****************************************************************************************