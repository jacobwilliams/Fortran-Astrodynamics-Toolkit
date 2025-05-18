!*****************************************************************************************
!>
!  Drag model

    module drag_module

    use kind_module

    implicit none

    private

    public :: drag_acceleration

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Acceleration due to atmospheric drag.

    pure subroutine drag_acceleration(vrel,cd,area,m,rho,acc)

    implicit none

    real(wp),dimension(3),intent(in)  :: vrel   !! velocity relative to the atmosphere [km/s]
    real(wp),intent(in)               :: cd     !! spacecraft drag coefficient [--]
    real(wp),intent(in)               :: area   !! cross-section area [km^2]
    real(wp),intent(in)               :: m      !! spacecraft mass [kg]
    real(wp),intent(in)               :: rho    !! atmospheric density [kg/km^3]
    real(wp),dimension(3),intent(out) :: acc    !! drag acceleration vector [km/s^2]

    real(wp) :: vrel_mag  !! magnitude of the relative velocity [km/s]

    vrel_mag = norm2(vrel)

    acc = - (0.5_wp * rho * cd * area * vrel_mag / m) * vrel

    end subroutine drag_acceleration
!*****************************************************************************************

!*****************************************************************************************
    end module drag_module
!*****************************************************************************************
