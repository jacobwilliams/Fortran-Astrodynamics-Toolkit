!*****************************************************************************************
!>
!  Gravity model
!
!### See also
!  * [[geopotential_module]]

    module gravity_module

    use kind_module
    use numbers_module
    use iso_fortran_env, only: output_unit

    implicit none

    private

    public :: third_body_gravity
    public :: gravity_j2_j3_j4

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Third-body (pointmass) gravitational acceleration.

    subroutine third_body_gravity(r,rb,mu,acc)

    implicit none

    real(wp),dimension(3),intent(in)  :: r   !! satellite position vector [km]
    real(wp),dimension(3),intent(in)  :: rb  !! third-body position vector [km]
    real(wp),intent(in)               :: mu  !! third-body gravitational parameter [km^3/s^2]
    real(wp),dimension(3),intent(out) :: acc !! gravity acceleration vector [km/s^2]

    real(wp),dimension(3) :: r_sc_b  !! vector from third-body to spacecraft [km]
    real(wp) :: rb_mag               !! distance between origin and third-body [km]
    real(wp) :: r_sc_b_mag           !! distance between spacecraft and third-body [km]

    r_sc_b     = rb - r
    r_sc_b_mag = norm2(r_sc_b)
    rb_mag     = norm2(rb)
    acc        = (mu/r_sc_b_mag**3)*r_sc_b - (mu/rb_mag**3)*rb

    end subroutine third_body_gravity
!*****************************************************************************************

!*****************************************************************************************
!>
!  Gravitational acceleration due to simplified spherical harmonic
!  expansion (only the J2-J4 terms are used).
!
!### Reference
!  * http://www.ni.com/pdf/manuals/370762a.pdf

    subroutine gravity_j2_j3_j4(r,mu,req,j2,j3,j4,acc)

    implicit none

    real(wp),dimension(3),intent(in)  :: r   !! satellite position vector [km]
    real(wp),intent(in)               :: mu  !! central body gravitational parameter [km^3/s^2]
    real(wp),intent(in)               :: req !! body equatorial radius [km]
    real(wp),intent(in)               :: j2  !! j2 coefficient
    real(wp),intent(in)               :: j3  !! j3 coefficient
    real(wp),intent(in)               :: j4  !! j4 coefficient
    real(wp),dimension(3),intent(out) :: acc !! gravity acceleration vector [km/s^2]

    real(wp) :: mmor3,reqor,reqor2,reqor3,reqor4,&
                rmag,rmag2,rmag3,rzor,rzor2,rzor3,rzor4,c,d

    rmag2 = dot_product(r,r)
    rmag  = sqrt(rmag2)

    if (rmag==zero) then

        write(output_unit,'(A)') 'Error in gravity_j2_j3_j4: spacecraft at center of body.'
        acc = zero

    else

        rmag3  = rmag*rmag2
        mmor3  = -mu/rmag3
        reqor  = req/rmag
        reqor2 = reqor*reqor
        reqor3 = reqor2*reqor
        reqor4 = reqor3*reqor
        rzor   = r(3)/rmag
        rzor2  = rzor*rzor
        rzor3  = rzor2*rzor
        rzor4  = rzor3*rzor

        c = mmor3 * (1.0_wp - 1.5_wp*J2*reqor2*(5.0_wp*rzor2-1.0_wp) + &
                     2.5_wp*J3*reqor3*(3.0_wp*rzor-7.0_wp*rzor3) - &
                     0.625_wp*J4*reqor4*(3.0_wp-42.0_wp*rzor2+63.0_wp*rzor4))
        d = mmor3 * (r(3) + 1.5_wp*J2*reqor2*(3.0_wp-5.0_wp*rzor2)*r(3) + &
                   0.5_wp*J3*reqor3*(30.0_wp*rzor*r(3)-35.0_wp*rzor3*r(3)-3.0_wp*rmag) - &
                   0.625_wp*J4*reqor4*(15.0_wp-70.0_wp*rzor2+63.0_wp*rzor4)*r(3))

        acc(1) = c * r(1)
        acc(2) = c * r(2)
        acc(3) = d

    end if

    end subroutine gravity_J2_J3_J4
!*****************************************************************************************

!*****************************************************************************************
    end module gravity_module
!*****************************************************************************************
