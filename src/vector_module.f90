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
    use numbers_module, only: one
    
    implicit none
    
    private
    
    public :: cross
    public :: unit
    public :: ucross
    public :: rodrigues_rotation
    
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
!****f* vector_module/rodrigues_rotation
!
!  NAME
!    rodrigues_rotation
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

    subroutine rodrigues_rotation(v,k,theta,vrot)
    
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
     
    end subroutine rodrigues_rotation
!*****************************************************************************************
    
!*****************************************************************************************
    end module vector_module
!*****************************************************************************************