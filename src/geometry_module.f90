!*****************************************************************************************
!> author: Jacob Williams
!  date: 2012
!
!  Geometry routines.

    module geometry_module

    use kind_module,      only: wp
    use vector_module,    only: cross

    implicit none

    private

    !public routines:
    public :: locpt
    public :: distance_from_point_to_line
    public :: distance_from_point_to_line_segment
    public :: distance_from_point_to_path

    !unit test routine:
    public :: geometry_unit_test

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  given a polygonal line connecting the vertices (x(i),y(i))
!  (i = 1,...,n) taken in this order. it is assumed that the
!  polygonal path is a loop, where (x(n),y(n)) = (x(1),y(1))
!  or there is an arc from (x(n),y(n)) to (x(1),y(1)).
!
!  (x0,y0) is an arbitrary point and l and m are variables.
!  l and m are assigned the following values:
!
!     l = -1   if (x0,y0) is outside the polygonal path
!     l =  0   if (x0,y0) lies on the polygonal path
!     l =  1   if (x0,y0) is inside the polygonal path
!
!  m = 0 if (x0,y0) is on or outside the path. if (x0,y0)
!  is inside the path then m is the winding number of the
!  path around the point (x0,y0).
!
!# History
!  * Original version from the NSWC Library
!  * Modified by J. Williams : 08/04/2012 : refactored to modern Fortran

    pure subroutine locpt (x0, y0, x, y, n, l, m)

    implicit none

    !arguments:
    integer,intent(in)               :: n
    real(wp),intent(in)              :: x0
    real(wp),intent(in)              :: y0
    real(wp),dimension(n),intent(in) :: x
    real(wp),dimension(n),intent(in) :: y
    integer,intent(out)              :: l
    integer,intent(out)              :: m

    !constants:
    real(wp),parameter :: eps = epsilon(1.0_wp)
    real(wp),parameter :: pi  = atan2(0.0_wp, -1.0_wp)
    real(wp),parameter :: pi2 = 2.0_wp*pi
    real(wp),parameter :: tol = 4.0_wp*eps*pi

    !local variables:
    integer  :: i,n0
    real(wp) :: u,v,theta1,sum,theta,angle,thetai

    n0 = n
    if (x(1) == x(n) .and. y(1) == y(n)) n0 = n - 1
    l = -1
    m = 0
    u = x(1) - x0
    v = y(1) - y0

    if (u == 0.0_wp .and. v == 0.0_wp) then

        l = 0  !(x0, y0) is on the boundary of the path

    else

        if (n0 >= 2) then

            theta1 = atan2(v, u)
            sum = 0.0_wp
            theta = theta1

            do i = 2,n0

                u = x(i) - x0
                v = y(i) - y0

                if (u == 0.0_wp .and. v == 0.0_wp) then
                    l = 0  !(x0, y0) is on the boundary of the path
                    exit
                end if

                thetai = atan2(v, u)
                angle = abs(thetai - theta)
                if (abs(angle - pi) < tol) then
                    l = 0  !(x0, y0) is on the boundary of the path
                    exit
                end if

                if (angle > pi) angle = angle - pi2
                if (theta > thetai) angle = -angle
                sum = sum + angle
                theta = thetai

            end do

            if (l/=0) then

                angle = abs(theta1 - theta)
                if (abs(angle - pi) < tol) then

                    l = 0  !(x0, y0) is on the boundary of the path

                else

                    if (angle > pi) angle = angle - pi2
                    if (theta > theta1) angle = -angle
                    sum = sum + angle    !sum = 2*pi*m where m is the winding number
                    m = abs(sum)/pi2 + 0.2_wp
                    if (m /= 0) then
                        l = 1
                        if (sum < 0.0_wp) m = -m
                    end if

                end if

            end if

        end if

    end if

    end subroutine locpt
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date:8/2012
!
!  Compute the distance between the point X and the line defined
!  by the two points X1 and X2.
!
!# References
!  1. http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html

    pure function distance_from_point_to_line (x1, x2, x) result(d)

    implicit none

    real(wp),dimension(3),intent(in) :: x1
    real(wp),dimension(3),intent(in) :: x2
    real(wp),dimension(3),intent(in) :: x
    real(wp)                         :: d

    real(wp),dimension(3) :: x21
    real(wp) :: x21_mag

    x21 = x2 - x1
    x21_mag = norm2(x21)

    if (x21_mag/=0.0_wp) then

        d = norm2( cross( x21, x1 - x ) ) / x21_mag

    else

        d = norm2(x1 - x)

    end if

    end function distance_from_point_to_line
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date:8/2012
!
!  Compute the distance between a line segment and a point.
!
!# References
!  1. http://forums.codeguru.com/showthread.php?194400-Distance-between-point-and-line-segment
!
!@note x,x1,x2 should all be the same length

    pure function distance_from_point_to_line_segment(x1, x2, x) result(d)

    implicit none

    real(wp),dimension(:),intent(in) :: x1
    real(wp),dimension(:),intent(in) :: x2
    real(wp),dimension(:),intent(in) :: x
    real(wp)                         :: d

    real(wp),dimension(size(x1)) :: x12
    real(wp) :: s
    real(wp) :: x12_mag

    x12 = x1 - x2
    x12_mag = norm2(x12)

    if (x12_mag==0.0_wp) then

        d = norm2(x1 - x)

    else

        s = dot_product( x1-x, x12 ) / dot_product( x12, x12 )

        !if the projection is not on the segment,
        ! then use the closest end point:
        if (s<0.0_wp) then
            s = 0.0_wp
        else if (s>1.0_wp) then
            s = 1.0_wp
        end if

        d = norm2( x - x1 + s*x12 )

    end if

    end function distance_from_point_to_line_segment
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date:8/2012
!
!  Compute the distance between a point and a polygonal path.
!  Given a point (x0,y0), and a path (x(n),y(n)), the distance
!  to the path is the distance to the closest line segment (x(i),y(i)).

    function distance_from_point_to_path(x0, y0, x, y, n) result(d)

    implicit none

    !arguments:
    integer,intent(in)                  :: n
    real(wp),intent(in)                 :: x0
    real(wp),intent(in)                 :: y0
    real(wp),dimension(n),intent(in)    :: x
    real(wp),dimension(n),intent(in)    :: y
    real(wp)                            :: d

    integer  :: l,m,i
    real(wp) :: dmin

    !is the point inside, outside, or on the path:
    call locpt (x0, y0, x, y, n, l, m)

    select case(l)

    case(1,-1)    !point is not on the path

        if (n==1) then

            !only one point in the path:
            d = norm2([x0-x(1), y0-y(1)])

        else

            do i=1,n    !loop through all line segments in the path
                        !the distance to the path is the distance from the closest line segment

                if (i<n) then
                    d = distance_from_point_to_line_segment( &
                            [x(i),y(i)], [x(i+1),y(i+1)], [x0,y0])
                else
                    !note: if 1st /= nth point, then have to check that line segment also.
                    if (x(1) /= x(n) .or. y(1) /= y(n)) then
                        d = distance_from_point_to_line_segment( &
                            [x(n),y(n)], [x(1),y(1)], [x0,y0])
                    end if
                end if

                !get lowest value:
                if (d<dmin .or. i==1) dmin = d

            end do

        end if

        !set sign of d:
        d = dmin * l   ! <0 if outside the path
                       ! >0 if inside the path

    case default    !point is on the path

        d = 0.0_wp

    end select

    end function distance_from_point_to_path
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date:8/2012
!
!  Unit test routine
!
!# Output
!
!    x0,y0=  0.59999999999999998       0.59999999999999998
!    l=           1
!    m=          -1
!    dist to path=  0.40000000000000002
!
!    x0,y0=   1.5000000000000000       0.50000000000000000
!    l=          -1
!    m=           0
!    dist to path= -0.50000000000000000
!
!    x0,y0=   1.0000000000000000        0.0000000000000000
!    l=           0
!    m=           0
!    dist to path=   0.0000000000000000

    subroutine geometry_unit_test()

    implicit none

    integer   :: l,m
    real(wp)  :: x0,y0

    !a 1x1 square:
    integer,parameter :: n = 4
    real(wp),dimension(n),parameter :: x = [0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp]
    real(wp),dimension(n),parameter :: y = [0.0_wp, 1.0_wp, 1.0_wp, 0.0_wp]

    x0 = 0.6_wp    !inside the path
    y0 = 0.6_wp
    call go()

    x0 = 1.5_wp    !outside the path
    y0 = 0.5_wp
    call go()

    x0 = 10.0_wp   !outside the path
    y0 = 0.0_wp
    call go()

    x0 = 1.0_wp    !on the path
    y0 = 0.0_wp
    call go()

    contains

        subroutine go()    !call locpt for the x0,y0 point, and print results
          implicit none

          real(wp) :: d

          call locpt (x0, y0, x, y, n, l, m)

          write(*,*) ''
          write(*,*) 'x0,y0=',x0,y0
          write(*,*) 'l=',l
          write(*,*) 'm=',m

          d = distance_from_point_to_path(x0, y0, x, y, n)
          write(*,*) 'dist to path=',d


        end subroutine go

    end subroutine geometry_unit_test
!*****************************************************************************************

!*****************************************************************************************
    end module geometry_module
!*****************************************************************************************
