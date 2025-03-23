!*****************************************************************************************
!> author: Jacob Williams
!
!  This module contains various routines related to relative motion.
!
!## Axis systems
!
!  Three different axis systems are used here.  They are:
!
!  * The **IJK** frame
!
!  * The **LVLH** frame, defined by:
!
!     * x-axis : completes the right handed system
!       (for a perfectly-circular orbit, the x-axis is \( \hat{\mathbf{v}} \))
!     * y-axis : \( -\hat{\mathbf{h}} \)
!     * z-axis : \( -\hat{\mathbf{r}} \)
!
!  * The **RSW** frame, defined by:
!
!     * x-axis : \( \hat{\mathbf{r}} \)
!     * y-axis : completes the right handed system
!       (for a perfectly-circular orbit, the y-axis is \( \hat{\mathbf{v}} \))
!     * z-axis : \( \hat{\mathbf{h}} \)

    module relative_motion_module

    use kind_module,    only: wp
    use numbers_module

    implicit none

    private

    abstract interface
        subroutine report_func(t,rv)  !! for reporting the points in the [[cw_propagator]].
        import :: wp
        implicit none
        real(wp),intent(in)              :: t   !! time [sec]
        real(wp),dimension(6),intent(in) :: rv  !! state [km,km/s]
        end subroutine report_func
    end interface

    interface from_ijk_to_lvlh  !! Conversion from IJK to LVLH
        procedure :: from_ijk_to_lvlh_mat  !! just returns matrices
        procedure :: from_ijk_to_lvlh_rv   !! transforms the r,v vectors
    end interface from_ijk_to_lvlh
    public :: from_ijk_to_lvlh

    interface from_lvlh_to_ijk  !! Conversion from LVLH to IJK
        procedure :: from_lvlh_to_ijk_mat  !! just returns matrices
        procedure :: from_lvlh_to_ijk_rv   !! transforms the r,v vectors
    end interface from_lvlh_to_ijk
    public :: from_lvlh_to_ijk

    interface from_ijk_to_rsw !! Conversion from IJK to RSW
        procedure :: from_ijk_to_rsw_mat  !! just returns matrices
        procedure :: from_ijk_to_rsw_rv   !! transforms the r,v vectors
    end interface from_ijk_to_rsw
    public :: from_ijk_to_rsw

    interface from_rsw_to_ijk  !! Conversion from RSW to IJK
        procedure :: from_rsw_to_ijk_mat  !! just returns matrices
        procedure :: from_rsw_to_ijk_rv   !! transforms the r,v vectors
    end interface from_rsw_to_ijk
    public :: from_rsw_to_ijk

    interface from_lvlh_to_rsw  !! Conversion from LVLH to RSW
        procedure :: from_lvlh_to_rsw_rv   !! transforms the r,v vectors
    end interface from_lvlh_to_rsw
    public :: from_lvlh_to_rsw

    interface from_rsw_to_lvlh !! Conversion from RSW to LVLH
        procedure :: from_rsw_to_lvlh_rv   !! transforms the r,v vectors
    end interface from_rsw_to_lvlh
    public :: from_rsw_to_lvlh

    public :: cw_equations
    public :: cw_propagator

    public :: relative_motion_test  !test routine

    contains
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 6/14/2015
!
!  Clohessy-Wiltshire equations for relative motion.
!
!  These apply to an RSW frame centered at the target spacecraft.
!
!# References
!   * [The Clohessy Wiltshire Model](http://courses.ae.utexas.edu/ase366k/cw_equations.pdf)

    function cw_equations(x0,dt,n) result(x)

    implicit none

    real(wp),dimension(6),intent(in) :: x0    !! initial state [r,v] of chaser (at t0) [km, km/s]
    real(wp),intent(in)              :: dt    !! elapsed time from t0 [sec]
    real(wp),intent(in)              :: n     !! mean motion of target orbit (`sqrt(mu/a**3)`) [1/sec]
    real(wp),dimension(6)            :: x     !! final state [r,v] of chaser [km, km/s]

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
!  date: 8/23/2015
!
!  Clohessy-Wiltshire propagation routine.
!
!# See also
!   * [[rk_module]]

    subroutine cw_propagator(t0,x0,h,n,tf,xf,report)

    implicit none

    real(wp),intent(in)               :: t0      !! initialize time [sec]
    real(wp),dimension(6),intent(in)  :: x0      !! initial state in RST coordinates [km,km/s]
    real(wp),intent(in)               :: h       !! abs(time step) [sec]
    real(wp),intent(in)               :: n       !! mean motion of target orbit (`sqrt(mu/a**3)`) [1/sec]
    real(wp),intent(in)               :: tf      !! final time [sec]
    real(wp),dimension(6),intent(out) :: xf      !! final state in RST coordinates [km,km/s]
    procedure(report_func),optional   :: report  !! to report each point

    real(wp) :: t,dt,t2
    real(wp),dimension(6) :: x
    logical :: last,export

    export = present(report)

    if (export) call report(t0,x0)  !first point

    if (h==zero) then
        xf = x0
    else

        t = t0
        x = x0
        dt = sign(h,tf-t0)  !time step  (correct sign)
        do
            t2 = t + dt
            last = ((dt>=zero .and. t2>=tf) .or. &  !adjust last time step
                    (dt<zero .and. t2<=tf))         !
            if (last) dt = tf-t                     !
            xf = cw_equations(x,dt,n)  ! propagate
            if (last) exit
            if (export) call report(t2,xf)   !intermediate point
            x = xf
            t = t2
        end do

    end if

    if (export) call report(tf,xf)   !last point

    end subroutine cw_propagator
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 4/19/2014
!
!  Compute the transformation matrices to convert IJK to LVLH.
!
!# See also
!   * [LVLH Transformations](http://degenerateconic.com/wp-content/uploads/2015/03/lvlh.pdf)

    subroutine from_ijk_to_lvlh_mat(r,v,a,c,cdot)

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

    end subroutine from_ijk_to_lvlh_mat
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 8/23/2014
!
!  Transform a position (and optionally velocity) vector from IJK to LVLH.

    subroutine from_ijk_to_lvlh_rv(rt_ijk,vt_ijk,r_ijk,v_ijk,dr_lvlh,dv_lvlh)

    implicit none

    real(wp),dimension(3),intent(in)           :: rt_ijk  !! Target IJK absolute position vector [km]
    real(wp),dimension(3),intent(in)           :: vt_ijk  !! Target IJK absolute position vector [km]
    real(wp),dimension(3),intent(in)           :: r_ijk   !! Chaser IJK absolute position vector [km]
    real(wp),dimension(3),intent(in)           :: v_ijk   !! Chaser IJK absolute position vector [km]
    real(wp),dimension(3),intent(out)          :: dr_lvlh !! Chaser LVLH position vector relative to target [km]
    real(wp),dimension(3),intent(out),optional :: dv_lvlh !! Chaser LVLH position vector relative to target [km]

    real(wp),dimension(3,3) :: c
    real(wp),dimension(3,3) :: cdot
    real(wp),dimension(3) :: dr_ijk, dv_ijk

    !IJK state of chaser relative to target:
    dr_ijk = r_ijk - rt_ijk    ! [target + delta = chaser]

    if (present(dv_lvlh)) then

        dv_ijk = v_ijk - vt_ijk ! [target + delta = chaser]

        call from_ijk_to_lvlh(rt_ijk,vt_ijk,c=c,cdot=cdot)

        dr_lvlh = matmul( c, dr_ijk )
        dv_lvlh = matmul( cdot, dr_ijk ) + matmul( c, dv_ijk )

    else

        call from_ijk_to_lvlh(r_ijk,v_ijk,c=c)

        dr_lvlh = matmul( c, dr_ijk )

    end if

    end subroutine from_ijk_to_lvlh_rv
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 4/19/2014
!
!  Compute the transformation matrices to convert LVLH to IJK.
!
!# See also
!   * [LVLH Transformations](http://degenerateconic.com/wp-content/uploads/2015/03/lvlh.pdf)

    subroutine from_lvlh_to_ijk_mat(r,v,a,c,cdot)

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

    end subroutine from_lvlh_to_ijk_mat
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 8/23/2014
!
!  Transform a position (and optionally velocity) vector from LVLH to IJK.

    subroutine from_lvlh_to_ijk_rv(rt_ijk,vt_ijk,dr_lvlh,dv_lvlh,r_ijk,v_ijk)

    implicit none

    real(wp),dimension(3),intent(in)            :: rt_ijk   !! Target IJK absolute position vector [km]
    real(wp),dimension(3),intent(in)            :: vt_ijk   !! Target IJK absolute position vector [km]
    real(wp),dimension(3),intent(in)            :: dr_lvlh  !! Chaser LVLH position vector relative to target [km]
    real(wp),dimension(3),intent(in)            :: dv_lvlh  !! Chaser LVLH position vector relative to target [km]
    real(wp),dimension(3),intent(out)           :: r_ijk    !! Chaser IJK absolute position vector [km]
    real(wp),dimension(3),intent(out),optional  :: v_ijk    !! Chaser IJK absolute position vector [km]

    real(wp),dimension(3,3) :: c
    real(wp),dimension(3,3) :: cdot

    if (present(v_ijk)) then

        call from_lvlh_to_ijk(rt_ijk,vt_ijk,c=c,cdot=cdot)

        !chaser = target + delta:
        r_ijk = rt_ijk + matmul( c, dr_lvlh )
        v_ijk = vt_ijk + matmul( cdot, dr_lvlh ) + matmul( c, dv_lvlh )

    else

        call from_lvlh_to_ijk(rt_ijk,vt_ijk,c=c)

        r_ijk = rt_ijk + matmul( c, dr_lvlh )

    end if

    end subroutine from_lvlh_to_ijk_rv
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 4/19/2014
!
!  Compute the transformation matrices to convert IJK to RSW.

    subroutine from_ijk_to_rsw_mat(r,v,a,c,cdot)

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

    end subroutine from_ijk_to_rsw_mat
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 8/23/2014
!
!  Transform a position (and optionally velocity) vector from IJK to RSW.

    subroutine from_ijk_to_rsw_rv(rt_ijk,vt_ijk,r_ijk,v_ijk,dr_rsw,dv_rsw)

    implicit none

    real(wp),dimension(3),intent(in)           :: rt_ijk   !! Target IJK absolute position vector [km]
    real(wp),dimension(3),intent(in)           :: vt_ijk   !! Target IJK absolute position vector [km]
    real(wp),dimension(3),intent(in)           :: r_ijk    !! Chaser IJK absolute position vector [km]
    real(wp),dimension(3),intent(in)           :: v_ijk    !! Chaser IJK absolute position vector [km]
    real(wp),dimension(3),intent(out)          :: dr_rsw   !! Chaser RSW position vector relative to target [km]
    real(wp),dimension(3),intent(out),optional :: dv_rsw   !! Chaser RSW position vector relative to target [km]

    real(wp),dimension(3,3) :: c
    real(wp),dimension(3,3) :: cdot
    real(wp),dimension(3) :: dr_ijk, dv_ijk

    dr_ijk = r_ijk - rt_ijk  ! delta = chaser - target

    if (present(dv_rsw)) then

        dv_ijk = v_ijk - vt_ijk  ! delta = chaser - target

        call from_ijk_to_rsw(rt_ijk,vt_ijk,c=c,cdot=cdot)

        dr_rsw = matmul( c, dr_ijk )
        dv_rsw = matmul( cdot, dr_ijk ) + matmul( c, dv_ijk )

    else

        call from_ijk_to_rsw(rt_ijk,vt_ijk,c=c)

        dr_rsw = matmul( c, dr_ijk )

    end if

    end subroutine from_ijk_to_rsw_rv
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 4/19/2014
!
!  Compute the transformation matrices to convert RSW to IJK.

    subroutine from_rsw_to_ijk_mat(r,v,a,c,cdot)

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

    end subroutine from_rsw_to_ijk_mat
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 8/23/2014
!
!  Transform a position (and optionally velocity) vector from RSW to IJK.

    subroutine from_rsw_to_ijk_rv(rt_ijk,vt_ijk,dr_rsw,dv_rsw,r_ijk,v_ijk)

    implicit none

    real(wp),dimension(3),intent(in)            :: rt_ijk   !! Target IJK absolute position vector [km]
    real(wp),dimension(3),intent(in)            :: vt_ijk   !! Target IJK absolute position vector [km]
    real(wp),dimension(3),intent(in)            :: dr_rsw   !! Chaser RSW position vector [km]
    real(wp),dimension(3),intent(in)            :: dv_rsw   !! Chaser RSW velocity vector [km/s]
    real(wp),dimension(3),intent(out)           :: r_ijk    !! Chaser IJK absolute position vector [km]
    real(wp),dimension(3),intent(out),optional  :: v_ijk    !! Chaser IJK absolute velocity vector [km/s]

    real(wp),dimension(3,3) :: c
    real(wp),dimension(3,3) :: cdot

    if (present(v_ijk)) then

        call from_rsw_to_ijk(rt_ijk,vt_ijk,c=c,cdot=cdot)

        r_ijk = rt_ijk + matmul( c, dr_rsw )
        v_ijk = vt_ijk + matmul( cdot, dr_rsw ) + matmul( c, dv_rsw )

    else

        call from_rsw_to_ijk(rt_ijk,vt_ijk,c=c)

        r_ijk = rt_ijk + matmul( c, dr_rsw )

    end if

    end subroutine from_rsw_to_ijk_rv
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 8/23/2014
!
!  Transform a position (and optionally velocity) vector from RSW to LVLH.

    subroutine from_rsw_to_lvlh_rv(dr_rsw,dv_rsw,dr_lvlh,dv_lvlh)

    implicit none

    real(wp),dimension(3),intent(in)           :: dr_rsw  !! Chaser RSW position vector relative to target [km]
    real(wp),dimension(3),intent(in)           :: dv_rsw  !! Chaser RSW velocity vector relative to target [km/s]
    real(wp),dimension(3),intent(out)          :: dr_lvlh !! Chaser LVLH position vector relative to target [km]
    real(wp),dimension(3),intent(out),optional :: dv_lvlh !! Chaser LVLH position vector relative to target [km]

    dr_lvlh(1) =  dr_rsw(2)
    dr_lvlh(2) = -dr_rsw(3)
    dr_lvlh(3) = -dr_rsw(1)

    if (present(dv_lvlh)) then
        dv_lvlh(1) =  dv_rsw(2)
        dv_lvlh(2) = -dv_rsw(3)
        dv_lvlh(3) = -dv_rsw(1)
    end if

    end subroutine from_rsw_to_lvlh_rv
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 8/23/2014
!
!  Transform a position (and optionally velocity) vector from LVLH to RSW.

    subroutine from_lvlh_to_rsw_rv(dr_lvlh,dv_lvlh,dr_rsw,dv_rsw)

    implicit none

    real(wp),dimension(3),intent(in)           :: dr_lvlh !! Chaser LVLH position vector relative to target [km]
    real(wp),dimension(3),intent(in)           :: dv_lvlh !! Chaser LVLH position vector relative to target [km]
    real(wp),dimension(3),intent(out)          :: dr_rsw  !! Chaser RSW position vector relative to target [km]
    real(wp),dimension(3),intent(out),optional :: dv_rsw  !! Chaser RSW velocity vector relative to target [km/s]

    dr_rsw(2) =  dr_lvlh(1)
    dr_rsw(3) = -dr_lvlh(2)
    dr_rsw(1) = -dr_lvlh(3)

    if (present(dv_rsw)) then
        dv_rsw(2) =  dv_lvlh(1)
        dv_rsw(3) = -dv_lvlh(2)
        dv_rsw(1) = -dv_lvlh(3)
    end if

    end subroutine from_lvlh_to_rsw_rv
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
