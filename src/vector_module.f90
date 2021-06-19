!*****************************************************************************************
!> author: Jacob Williams
!
!  Routines for the manipulation of vectors.

    module vector_module

    use kind_module,    only: wp
    use numbers_module, only: one,zero,pi

    implicit none

    private

    integer,parameter,public :: x_axis = 1
    integer,parameter,public :: y_axis = 2
    integer,parameter,public :: z_axis = 3

    real(wp),dimension(3),parameter,public :: x_unit = [one,zero,zero] !! x-axis unit vector
    real(wp),dimension(3),parameter,public :: y_unit = [zero,one,zero] !! y-axis unit vector
    real(wp),dimension(3),parameter,public :: z_unit = [zero,zero,one] !! z-axis unit vector

    public :: cross
    public :: unit
    public :: uhat_dot
    public :: ucross
    public :: axis_angle_rotation
    public :: cross_matrix
    public :: outer_product
    public :: box_product
    public :: vector_projection
    public :: vector_projection_on_plane
    public :: axis_angle_rotation_to_rotation_matrix
    public :: spherical_to_cartesian
    public :: cartesian_to_spherical
    public :: rotation_matrix
    public :: rotation_matrix_dot
    public :: angle_between_vectors

    !test routine:
    public :: vector_test

    contains
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!
!  Cross product of two 3x1 vectors

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
!> author: Jacob Williams
!
!  Unit vector

    pure function unit(r) result(u)

    implicit none

    real(wp),dimension(:),intent(in) :: r
    real(wp),dimension(size(r))      :: u

    real(wp) :: rmag

    rmag = norm2(r)

    if (rmag==zero) then
        u = zero
    else
        u = r / rmag
    end if

    end function unit
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!
!  Time derivative of a unit vector.

    pure function uhat_dot(u,udot) result(uhatd)

    implicit none

    real(wp),dimension(3),intent(in) :: u      !! vector [`u`]
    real(wp),dimension(3),intent(in) :: udot   !! derivative of vector [`du/dt`]
    real(wp),dimension(3)            :: uhatd  !! derivative of unit vector [`d(uhat)/dt`]

    real(wp)              :: umag  !! vector magnitude
    real(wp),dimension(3) :: uhat  !! unit vector

    umag = norm2(u)

    if (umag == zero) then  !singularity
        uhatd = zero
    else
        uhat = u / umag
        uhatd = ( udot - dot_product(uhat,udot)*uhat ) / umag
    end if

    end function uhat_dot
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!
!  Unit vector of the cross product of two 3x1 vectors

    pure function ucross(v1,v2) result(u)

    implicit none

    real(wp),dimension(3),intent(in) :: v1
    real(wp),dimension(3),intent(in) :: v2
    real(wp),dimension(3)            :: u

    u = unit(cross(v1,v2))

    end function ucross
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 7/20/2014
!
!  Rotate a 3x1 vector in space, given an axis and angle of rotation.
!
!# Reference
!   * [Wikipedia](http://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula)

    pure subroutine axis_angle_rotation(v,k,theta,vrot)

    implicit none

    real(wp),dimension(3),intent(in)  :: v      !! vector to rotate
    real(wp),dimension(3),intent(in)  :: k      !! rotation axis
    real(wp),intent(in)               :: theta  !! rotation angle [rad]
    real(wp),dimension(3),intent(out) :: vrot   !! result

    real(wp),dimension(3) :: khat
    real(wp) :: ct,st

    ct = cos(theta)
    st = sin(theta)
    khat = unit(k)   !rotation axis unit vector

    vrot = v*ct + cross(khat,v)*st + khat*dot_product(khat,v)*(one-ct)

    end subroutine axis_angle_rotation
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 7/20/2014
!
!  Computes the cross product matrix, where:
!  ``cross(a,b) == matmul(cross_matrix(a),b)``

    pure function cross_matrix(r) result(rcross)

    implicit none

    real(wp),dimension(3),intent(in) :: r
    real(wp),dimension(3,3)          :: rcross

    rcross(:,1) = [zero,r(3),-r(2)]
    rcross(:,2) = [-r(3),zero,r(1)]
    rcross(:,3) = [r(2),-r(1),zero]

    end function cross_matrix
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 7/21/2014
!
!  Computes the outer product of the two vectors.

    pure function outer_product(a,b) result(c)

    implicit none

    real(wp),dimension(:),intent(in)    :: a
    real(wp),dimension(:),intent(in)    :: b
    real(wp),dimension(size(a),size(b)) :: c

    integer :: i

    do i=1,size(b)
        c(:,i) = a*b(i)
    end do

    end function outer_product
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 7/21/2014
!
!  Computes the box product (scalar triple product) of the three vectors.

    pure function box_product(a,b,c) result(d)

    implicit none

    real(wp),dimension(:),intent(in) :: a
    real(wp),dimension(:),intent(in) :: b
    real(wp),dimension(:),intent(in) :: c
    real(wp) :: d

    d = dot_product(a,cross(b,c))

    end function box_product
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 7/21/2014
!
!  The projection of one vector onto another vector.
!
!# Reference
!   * [Wikipedia](http://en.wikipedia.org/wiki/Gram-Schmidt_process)

    pure function vector_projection(a,b) result(c)

    implicit none

    real(wp),dimension(:),intent(in)       :: a  !! the original vector
    real(wp),dimension(size(a)),intent(in) :: b  !! the vector to project on to
    real(wp),dimension(size(a))            :: c  !! the projection of a onto b

    real(wp) :: bmag2

    bmag2 = dot_product(b,b)

    if (bmag2==zero) then
        c = zero
    else
        c = b * dot_product(a,b) / bmag2
    end if

    end function vector_projection
!*****************************************************************************************

!*****************************************************************************************
!>
!  Project a vector onto a plane.
!
!# Reference
!   * [Projection of a Vector onto a Plane](http://www.maplesoft.com/support/help/Maple/view.aspx?path=MathApps/ProjectionOfVectorOntoPlane)

    pure subroutine vector_projection_on_plane(a,b,c)

    implicit none

    real(wp),dimension(3),intent(in)  :: a !! the original vector
    real(wp),dimension(3),intent(in)  :: b !! the plane to project on to (a normal vector)
    real(wp),dimension(3),intent(out) :: c !! the projection of a onto the b plane

    c = a - vector_projection(a,b)

    end subroutine vector_projection_on_plane
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 7/20/2014
!
!  Computes the rotation matrix that corresponds to a
!  rotation about the axis `k` by an angle `theta`.

    pure subroutine axis_angle_rotation_to_rotation_matrix(k,theta,rotmat)

    implicit none

    real(wp),dimension(3),intent(in)    :: k        !! rotation axis
    real(wp),intent(in)                 :: theta    !! rotation angle [rad]
    real(wp),dimension(3,3),intent(out) :: rotmat   !! rotation matrix

    real(wp),dimension(3,3),parameter :: I = &
            reshape([one,zero,zero,zero,one,zero,zero,zero,one],[3,3]) !! 3x3 identity matrix

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
!> author: Jacob Williams
!  date: 9/24/2014
!
!  Convert spherical (r,alpha,beta) to Cartesian (x,y,z).

    pure function spherical_to_cartesian(r,alpha,beta) result(rvec)

    implicit none

    real(wp),intent(in)   :: r        !! magnitude
    real(wp),intent(in)   :: alpha    !! right ascension [rad]
    real(wp),intent(in)   :: beta     !! declination [rad]
    real(wp),dimension(3) :: rvec     !! [x,y,z] vector

    rvec(1) = r * cos(alpha) * cos(beta)
    rvec(2) = r * sin(alpha) * cos(beta)
    rvec(3) = r * sin(beta)

    end function spherical_to_cartesian
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 3/6/2016
!
!  Convert Cartesian (x,y,z) to spherical (r,alpha,beta).

    pure subroutine cartesian_to_spherical(rvec,r,alpha,beta)

    implicit none

    real(wp),dimension(3),intent(in) :: rvec     !! [x,y,z] vector
    real(wp),intent(out)             :: r        !! magnitude
    real(wp),intent(out)             :: alpha    !! right ascension [rad]
    real(wp),intent(out)             :: beta     !! declination [rad]

    real(wp) :: r1

    r1 = rvec(1)*rvec(1)+rvec(2)*rvec(2)
    r  = sqrt(r1+rvec(3)*rvec(3))

    if (r/=zero) then
        beta = atan2(rvec(3),sqrt(r1))
        if (r1/=zero) then
            alpha = atan2(rvec(2),rvec(1))
        else
            alpha = zero
        end if
    else
        alpha = zero
        beta = zero
    end if

    end subroutine cartesian_to_spherical
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 2/3/2015
!
!  The 3x3 rotation matrix for a rotation about the x, y, or z-axis.
!
!  EXAMPLE
!```Fortran
!    real(wp),dimension(3,3) :: rotmat
!    real(wp),dimension(3) :: vec,vec2
!    real(wp) :: ang
!    ang = pi / 4.0_wp
!    vec = [1.414_wp, 0.0_wp, 0.0_wp]
!    rotmat = rotation_matrix(z_axis,ang)
!    vec2 = matmul(rotmat,vec)
!```

    pure function rotation_matrix(axis,angle) result(rotmat)

    implicit none

    real(wp),dimension(3,3) :: rotmat   !! the rotation matrix
    integer,intent(in)      :: axis     !! x_axis, y_axis, or z_axis
    real(wp),intent(in)     :: angle    !! angle in radians

    real(wp) :: c,s

    !precompute these:
    c = cos(angle)
    s = sin(angle)

    select case (axis)
    case(x_axis); rotmat = reshape([one, zero, zero, zero, c, -s, zero, s, c],[3,3])
    case(y_axis); rotmat = reshape([c, zero, s, zero, one, zero, -s, zero, c],[3,3])
    case(z_axis); rotmat = reshape([c, -s, zero, s, c, zero, zero, zero, one],[3,3])
    case default; rotmat = zero
    end select

    end function rotation_matrix
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 3/19/2016
!
!  Time derivative of the 3x3 rotation matrix
!  for a rotation about the x, y, or z-axis.

    pure function rotation_matrix_dot(axis,angle,angledot) result(rotmatdot)

    implicit none

    real(wp),dimension(3,3) :: rotmatdot   !! the rotation matrix derivative \( d \mathbf{C} / d t \)
    integer,intent(in)      :: axis        !! x_axis, y_axis, or z_axis
    real(wp),intent(in)     :: angle       !! angle in radians
    real(wp),intent(in)     :: angledot    !! time derivative of angle in radians/sec

    real(wp) :: c,s

    !precompute these:
    c = cos(angle)
    s = sin(angle)

    !first compute d[C]/da (time derivate w.r.t. the angle):
    select case (axis)
    case(x_axis); rotmatdot = reshape([zero, zero, zero, zero, -s, -c, zero, c, -s],[3,3])
    case(y_axis); rotmatdot = reshape([-s, zero, c, zero, zero, zero, -c, zero, -s],[3,3])
    case(z_axis); rotmatdot = reshape([-s, -c, zero, c, -s, zero, zero, zero, zero],[3,3])
    case default
        rotmatdot = zero
        return
    end select

    rotmatdot = rotmatdot * angledot    ! d[C]/dt = d[C]/da * da/dt

    end function rotation_matrix_dot
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 3/13/2015
!
!  The angle between two vectors (in radians).

    pure function angle_between_vectors(v1,v2) result(ang)

    implicit none

    real(wp)                         :: ang  !! [rad]
    real(wp),dimension(3),intent(in) :: v1
    real(wp),dimension(3),intent(in) :: v2

    real(wp) :: d,c

    d   = dot_product(v1,v2)
    c   = norm2(cross(v1,v2))
    ang = atan2(c,d)

    end function angle_between_vectors
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 7/20/2014
!
!  Unit test routine for the [[vector_module]].

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

    !z-axis rotation test:
    theta = pi / 4.0_wp
    v = [one/cos(theta), 0.0_wp, 0.0_wp]
    rotmat = rotation_matrix(z_axis,theta)
    v2 = matmul(rotmat,v)
    write(*,*) v2    !should be [1, -1, 0]


    end subroutine vector_test
 !****************************************************************************************

!*****************************************************************************************
    end module vector_module
!*****************************************************************************************
