!*****************************************************************************************
!> author: Jacob Williams
!
!  Routines for relative motion.

    module relative_motion_module

    use kind_module,    only: wp
    use numbers_module
    
    implicit none
    
    private
    
    public :: cw_equations
    public :: from_ijk_to_lvlh, from_lvlh_to_ijk
    public :: from_ijk_to_rsw, from_rsv_to_ijk
    
    public :: relative_motion_test  !test routine
    
    contains
!*****************************************************************************************
    
!*****************************************************************************************
!> author: Jacob Williams
!  date: 6/14/2015
!
!  Clohessy-Wiltshire equations for relative motion (RSW frame).
!
!  The frame is centered at the target spacecraft and is defined by:
!
!   * x-axis : unit(r)
!   * z-axis : unit(h)
!   * y-axis : completes the right handed system
!     (for a perfectly-circular orbit, the y-axis is unit(v))
!
!# References
!   * [The Clohessy Wiltshire Model](http://courses.ae.utexas.edu/ase366k/cw_equations.pdf)

    pure function cw_equations(x0,dt,n) result(x)

    implicit none

    real(wp),dimension(6),intent(in) :: x0    !! initial state of chaser (at t0) [km, km/s]
    real(wp),intent(in)              :: dt    !! elapsed time from t0 [sec]
    real(wp),intent(in)              :: n     !! mean motion of target orbit (`sqrt(mu/a**3)`) [1/sec]
    real(wp),dimension(6)            :: x     !! final state of chaser [km, km/s]
   
    real(wp) :: nt,cnt,snt
    
    if (dt==zero) then
        x = x0
    else
    
        if (n==zero) then
            error stop 'Error: Target orbit mean motion must be non-zero.'
        else
        
            nt = n*dt
            cnt = cos(nt)
            snt = sin(nt)
   
            x(1) = (four-three*cnt)*x0(1) + (snt/n)*x0(4) + (two/n)*(one-cnt)*x0(5)
            x(2) = six*(snt-nt)*x0(1) + x0(2) - (two/n)*(one-cnt)*x0(4) + one/n*(four*snt-three*nt)*x0(5)
            x(3) = cnt*x0(3) + (snt/n)*x0(6)
            x(4) = three*n*snt*x0(1) + cnt*x0(4) + two*snt*x0(5)
            x(5) = -(six*n*(one-cnt))*x0(1) - two*snt*x0(4) + (four*cnt-three)*x0(5)
            x(6) = -n*snt*x0(3) + cnt*x0(6)
        
        end if
    
    end if

    end function cw_equations
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 4/19/2014
!
!  Compute the transformation matrices to convert IJK to LVLH.
!
!  The LVLH frame is defined by:
!
!   * z-axis : -unit(r)
!   * y-axis : -unit(cross(r,v))
!   * x-axis completes the right handed system
!     (for a perfectly-circular orbit, the x-axis is unit(v))
!
!# See also
!   * [LVLH Transformations](http://degenerateconic.com/wp-content/uploads/2015/03/lvlh.pdf)

    subroutine from_ijk_to_lvlh(r,v,a,c,cdot) 

    use vector_module, only: unit, cross, uhat_dot

    implicit none

    real(wp),dimension(3),intent(in)             :: r      !! position vector of target [km]
    real(wp),dimension(3),intent(in)             :: v      !! velocity vector of target [km/s]
    real(wp),dimension(3),intent(in),optional    :: a      !! acceleration vector of target [km/s^2] 
                                                           !! (if not present, then a torque-free force model is assumed)
    real(wp),dimension(3,3),intent(out)          :: c      !! C transformation matrix
    real(wp),dimension(3,3),intent(out),optional :: cdot   !! CDOT transformation matrix

    real(wp),dimension(3) :: ex_hat,ex_hat_dot
    real(wp),dimension(3) :: ey_hat,ey_hat_dot
    real(wp),dimension(3) :: ez_hat,ez_hat_dot
    real(wp),dimension(3) :: h,h_hat,h_dot

    h      = cross(r,v)
    h_hat  = unit(h)
    ez_hat = -unit(r)
    ey_hat = -h_hat
    ex_hat = cross(ey_hat,ez_hat)
    
    c(1,:) = ex_hat
    c(2,:) = ey_hat
    c(3,:) = ez_hat
    
    if (present(cdot)) then
    
        ez_hat_dot = -uhat_dot(r,v)
 
        if (present(a)) then
            h_dot      = cross(r,a)
            ey_hat_dot = -uhat_dot(h, h_dot) 
            ex_hat_dot = cross(ey_hat_dot,ez_hat) + cross(ey_hat,ez_hat_dot)
        else !assume no external torque
            ey_hat_dot = zero
            ex_hat_dot = cross(ey_hat,ez_hat_dot)
        end if

        cdot(1,:) = ex_hat_dot
        cdot(2,:) = ey_hat_dot
        cdot(3,:) = ez_hat_dot
    
    end if

    end subroutine from_ijk_to_lvlh      
!*****************************************************************************************    
    
!*****************************************************************************************
!> author: Jacob Williams
!  date: 4/19/2014
!
!  Compute the transformation matrices to convert LVLH to IJK.
!
!# See also
!   * [LVLH Transformations](http://degenerateconic.com/wp-content/uploads/2015/03/lvlh.pdf)

    subroutine from_lvlh_to_ijk(r,v,a,c,cdot)

    use vector_module, only: unit, cross
 
    implicit none

    real(wp),dimension(3),intent(in)             :: r     !! position vector of target [km]
    real(wp),dimension(3),intent(in)             :: v     !! velocity vector of target [km/s]
    real(wp),dimension(3),intent(in),optional    :: a     !! acceleration vector of target [km/s^2] 
                                                          !! (if not present, then a torque-free force model is assumed)
    real(wp),dimension(3,3),intent(out)          :: c     !! C transformation matrix
    real(wp),dimension(3,3),intent(out),optional :: cdot  !! CDOT transformation matrix

    call from_ijk_to_lvlh(r,v,a,c,cdot)

    c = transpose(c)
    
    if (present(cdot)) cdot = transpose(cdot)

    end subroutine from_lvlh_to_ijk
!***************************************************************************************** 

!*****************************************************************************************
!> author: Jacob Williams
!  date: 4/19/2014
!
!  Compute the transformation matrices to convert IJK to RSW.
!
!  The RSW frame is defined by:
!
!   * x-axis : unit(r)
!   * z-axis : unit(h)
!   * y-axis : completes the right handed system
!     (for a perfectly-circular orbit, the y-axis is unit(v))

    subroutine from_ijk_to_rsw(r,v,a,c,cdot) 

    use vector_module, only: unit, cross, uhat_dot

    implicit none

    real(wp),dimension(3),intent(in)             :: r      !! position vector of target [km]
    real(wp),dimension(3),intent(in)             :: v      !! velocity vector of target [km/s]
    real(wp),dimension(3),intent(in),optional    :: a      !! acceleration vector of target [km/s^2] 
                                                           !! (if not present, then a torque-free force model is assumed)
    real(wp),dimension(3,3),intent(out)          :: c      !! C transformation matrix
    real(wp),dimension(3,3),intent(out),optional :: cdot   !! CDOT transformation matrix

    real(wp),dimension(3) :: ex_hat,ex_hat_dot
    real(wp),dimension(3) :: ey_hat,ey_hat_dot
    real(wp),dimension(3) :: ez_hat,ez_hat_dot
    real(wp),dimension(3) :: h,h_hat,h_dot

    h      = cross(r,v)
    h_hat  = unit(h)
    ex_hat = unit(r)
    ez_hat = h_hat
    ey_hat = cross(ez_hat,ex_hat)
    
    c(1,:) = ex_hat
    c(2,:) = ey_hat
    c(3,:) = ez_hat

    !... need to verify the following ...
 
    if (present(cdot)) then
    
        ex_hat_dot = uhat_dot(r,v)
 
        if (present(a)) then
            h_dot      = cross(r,a)
            ez_hat_dot = uhat_dot(h, h_dot) 
            ey_hat_dot = cross(ez_hat_dot,ex_hat) + cross(ez_hat,ex_hat_dot)
        else !assume no external torque
            ez_hat_dot = zero
            ey_hat_dot = cross(ez_hat,ex_hat_dot)
        end if

        cdot(1,:) = ex_hat_dot
        cdot(2,:) = ey_hat_dot
        cdot(3,:) = ez_hat_dot
    
    end if
 
    end subroutine from_ijk_to_rsw
!*****************************************************************************************  

!*****************************************************************************************
!> author: Jacob Williams
!  date: 4/19/2014
!
!  Compute the transformation matrices to convert LVLH to RSW.

    subroutine from_rsv_to_ijk(r,v,a,c,cdot)

    use vector_module, only: unit, cross
 
    implicit none

    real(wp),dimension(3),intent(in)             :: r     !! position vector of target [km]
    real(wp),dimension(3),intent(in)             :: v     !! velocity vector of target [km/s]
    real(wp),dimension(3),intent(in),optional    :: a     !! acceleration vector of target [km/s^2] 
                                                          !! (if not present, then a torque-free force model is assumed)
    real(wp),dimension(3,3),intent(out)          :: c     !! C transformation matrix
    real(wp),dimension(3,3),intent(out),optional :: cdot  !! CDOT transformation matrix

    call from_ijk_to_rsw(r,v,a,c,cdot)

    c = transpose(c)
    
    if (present(cdot)) cdot = transpose(cdot)

    end subroutine from_rsv_to_ijk
!***************************************************************************************** 

!*****************************************************************************************  
!> author: Jacob Williams
!  date: 8/22/2015
!
!  Unit tests for the [[relative_motion_module]].

    subroutine relative_motion_test()
    
    implicit none

    real(wp),dimension(6),parameter :: target_eci_state = [-2301672.24489839_wp, &
                                                           -5371076.10250925_wp, &
                                                           -3421146.71530212_wp, &
                                                            6133.8624555516_wp,  &
                                                            306.265184163608_wp, &
                                                           -4597.13439017524_wp  ]

    real(wp),dimension(6),parameter :: chaser_eci_state = [-2255213.51862763_wp, &
                                                           -5366553.94133467_wp, &
                                                           -3453871.15040494_wp, &
                                                            6156.89588163809_wp, &
                                                            356.79933181917_wp,  &
                                                           -4565.88915429063_wp  ]

    real(wp),dimension(3) :: r_12_I, r1_I, r2_I
    real(wp),dimension(3) :: v_12_I, v1_I, v2_I
    real(wp),dimension(3) :: r_12_R, v_12_R
    real(wp),dimension(3,3) :: c,cdot
    
    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' relative_motion_test'
    write(*,*) '---------------'
    write(*,*) ''
    
    r1_I = target_eci_state(1:3)
    v1_I = target_eci_state(4:6)
    r2_I = chaser_eci_state(1:3)
    v2_I = chaser_eci_state(4:6)
    r_12_I = r2_I - r1_I
    v_12_I = v2_I - v1_I

    call from_ijk_to_lvlh(r1_I,v1_I,c=c,cdot=cdot)     

    r_12_R = matmul( c, r_12_I )
    v_12_R = matmul( cdot, r_12_I ) + matmul( c, v_12_I )

    write(*,'(A,*(D30.16,1X))') 'r_12_LVLH : ', r_12_R
    write(*,'(A,*(D30.16,1X))') 'v_12_LVLH : ', v_12_R
    write(*,*) ''

    call from_ijk_to_rsw(r1_I,v1_I,c=c,cdot=cdot)     

    r_12_R = matmul( c, r_12_I )
    v_12_R = matmul( cdot, r_12_I ) + matmul( c, v_12_I )

    write(*,'(A,*(D30.16,1X))') 'r_12_RSW : ', r_12_R
    write(*,'(A,*(D30.16,1X))') 'v_12_RSW : ', v_12_R
    write(*,*) ''

    end subroutine relative_motion_test
!*****************************************************************************************  

!*****************************************************************************************
    end module relative_motion_module
!*****************************************************************************************