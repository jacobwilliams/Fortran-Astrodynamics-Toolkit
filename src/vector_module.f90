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
    
    implicit none
    
    private
    
    public :: cross
    public :: unit
    public :: ucross
    
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
    end module vector_module
!*****************************************************************************************    