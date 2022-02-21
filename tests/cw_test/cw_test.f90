!*****************************************************************************************
!> author: Jacob Williams
!  date: 8/23/2015
!
!  Test program for the CW-equation propagator.
!
!  Propagate a relative state with an RK method and also using the CW equations.
!  Note that this test requires [pyplot-fortran](https://github.com/jacobwilliams/pyplot-fortran).
!
!# See also
!   * [[gravity_test]]

    program cw_test

    use fortran_astrodynamics_toolkit, wp => fat_wp
    use pyplot_module

    implicit none

    type,extends(rk8_10_class) :: spacecraft
        !! spacecraft propagation type:
        !! extend the rk class to include data used in the deriv routine
        type(geopotential_model_mueller) :: grav   !! central body geopotential model
        logical :: first = .true.  !! first point is being exported
        logical :: use_geopotential = .true. !! if true, use the geopotential model.
                                             !! if false, use pointmass.
    end type spacecraft

    integer,parameter               :: n             = 6                        !! number of state variables
    real(wp),parameter              :: mu_earth      = 398600.4418_wp           !! grav. param of earth [km3/s2]
    real(wp),parameter              :: orbit_inc     = 45.0_wp*deg2rad          !! orbit inclination [rad]
    real(wp),parameter              :: orbit_sma     = 6778.0_wp                !! orbit semimajor axis [km]
    character(len=*),parameter      :: gravfile      = '../grav/GGM03C.GEO'     !! grav coefficient file
    integer,parameter               :: grav_n        = 10                       !! max degree
    integer,parameter               :: grav_m        = 10                       !! max order
    real(wp),dimension(3),parameter :: r_chaser_lvlh = [1.0_wp,1.0_wp,1.0_wp]   !! chaser r relative to target (LVLH)
    real(wp),dimension(3),parameter :: v_chaser_lvlh = [0.001_wp,0.001_wp,0.001_wp] !! chaser v relative to target (LVLH)

    real(wp),dimension(:),allocatable :: x_lvlh_int, y_lvlh_int, z_lvlh_int
    real(wp),dimension(:),allocatable :: x_lvlh_int_pm, y_lvlh_int_pm, z_lvlh_int_pm
    real(wp),dimension(:),allocatable :: x_lvlh_cw, y_lvlh_cw, z_lvlh_cw
    real(wp),dimension(n)   :: rv_target_eci_t0, rv_target_eci_tf
    real(wp),dimension(n)   :: rv_chaser_rsw_t0, rv_chaser_rsw_tf
    real(wp),dimension(n)   :: rv_chaser_eci_t0, rv_chaser_eci_tf
    real(wp),dimension(n+n) :: x0               !! initial [target, chaser] states
    real(wp),dimension(n+n) :: xf               !! final [target, chaser] states
    type(pyplot)            :: plt              !! for making the plot
    real(wp)                :: mean_motion      !! [1/sec]
    type(spacecraft)        :: s
    real(wp)                :: t0,tf,dt
    logical                 :: status_ok
    integer                 :: iunit,istat,iunit_cw,iunit_pm
    logical                 :: cw_header_written
    integer                 :: icase

    cw_header_written = .false.

    !open output files:
    open(newunit=iunit,file   ='cw_test_integration_geopot.txt',   status='REPLACE',iostat=istat)
    open(newunit=iunit_cw,file='cw_test_cw_prop.txt',              status='REPLACE',iostat=istat)
    open(newunit=iunit_pm,file='cw_test_integration_pointmass.txt',status='REPLACE',iostat=istat)

    !constructors:
    call s%initialize(n=n+n,f=deriv,report=integrator_report) ! integrating target and chaser
    call s%grav%initialize(gravfile,grav_n,grav_m,status_ok)  ! gravitational field
    if (.not. status_ok) stop 'Error'

    !get an initial cartesian state vector:
    ! [circular, 45 deg inclination, radius of 6778 km]
    call els3pv(mu_earth,[mu_earth/orbit_sma,orbit_sma,orbit_inc,zero,zero,zero],rv_target_eci_t0)
    mean_motion = sqrt(mu_earth / orbit_sma**3) ! mean motion of target orbit

    !convert chaser state to eci frame (for integration):
    call from_lvlh_to_ijk(rt_ijk  = rv_target_eci_t0(1:3),&
                          vt_ijk  = rv_target_eci_t0(4:6),&
                          dr_lvlh = r_chaser_lvlh,&
                          dv_lvlh = v_chaser_lvlh,&
                          r_ijk   = rv_chaser_eci_t0(1:3),&  ! chaser in eci-Earth frame.
                          v_ijk   = rv_chaser_eci_t0(4:6))   !

    !also convert chaser state to RST (for CW propagation)
    call from_lvlh_to_rsw(r_chaser_lvlh,v_chaser_lvlh,rv_chaser_rsw_t0(1:3),rv_chaser_rsw_t0(4:6))

    write(*,*) ''
    write(*,'(A/,*(E30.16/))') 'rv_target_eci_t0:',rv_target_eci_t0
    write(*,'(A/,*(E30.16/))') 'rv_chaser_rsw_t0:',rv_chaser_rsw_t0
    write(*,'(A/,*(E30.16/))') 'rv_chaser_eci_t0:',rv_chaser_eci_t0
    write(*,'(A/,*(E30.16/))') 'r_chaser_lvlh_t0:',r_chaser_lvlh,v_chaser_lvlh
    write(*,*) ''

    do icase=1,2

        !initial conditions:
        t0 = zero           !initial ephemeris time (sec) -> this is 1/1/2000 12:00:00
        dt = 10.0_wp        !time step (sec)
        tf = 4.0_wp*hr2sec  !final time (sec)

        !integrate (both states):

        s%first = .true.
        s%use_geopotential = icase==1

        x0 = [rv_target_eci_t0,rv_chaser_eci_t0]
        call s%integrate(t0,x0,dt,tf,xf)
        rv_target_eci_tf = xf(1:6)
        rv_chaser_eci_tf = xf(7:12)

        if (s%use_geopotential) then
            write(*,'(A/,*(E30.16/))') 'rv_target_eci_tf (integrated, GEOPOT):',rv_target_eci_tf
        else
            write(*,'(A/,*(E30.16/))') 'rv_target_eci_tf (integrated, POINTMASS):',rv_target_eci_tf
        end if

    end do

    !also propagate with the CW equations:
    call cw_propagator(t0,rv_chaser_rsw_t0,dt,mean_motion,tf,rv_chaser_rsw_tf,cw_report)
    write(*,'(A/,*(E30.16/))') 'rv_chaser_eci_tf (integrated):',rv_chaser_eci_tf

    !plot both trajectories:

    write(*,*) 'generate plots...'

    call plt%initialize(grid=.true.,xlabel='LVLH x [km]',ylabel='LVLH y [km]',&
                            title='Integration vs CW Propagation',legend=.true.,figsize=[10,5])
    call plt%add_plot([zero],[zero],label='Target', &
                            linestyle='bo',markersize=5,istat=istat)
    call plt%add_plot(x_lvlh_int,y_lvlh_int, &
                            label='Integrated (10x10 gravity)', &
                            linestyle='b-',markersize=5,linewidth=2,istat=istat)
    call plt%add_plot(x_lvlh_int_pm,y_lvlh_int_pm, &
                            label='Integrated (pointmass gravity)', &
                            linestyle='g-',markersize=5,linewidth=2,istat=istat)
    call plt%add_plot(x_lvlh_cw,y_lvlh_cw, &
                            label='CW Equations', linestyle='r-',markersize=5,linewidth=2,istat=istat)
    call plt%savefig('cw_test_xy.png',istat=istat)
    call plt%destroy()

    call plt%initialize(grid=.true.,xlabel='LVLH x [km]',ylabel='LVLH z [km]',&
                            title='Integration vs CW Propagation',legend=.true.,figsize=[10,5])
    call plt%add_plot([zero],[zero],label='Target',linestyle='bo',markersize=5,istat=istat)
    call plt%add_plot(x_lvlh_int,z_lvlh_int, &
                            label='Integrated (10x10 gravity)', &
                            linestyle='b-',markersize=5,linewidth=2,istat=istat)
    call plt%add_plot(x_lvlh_int_pm,z_lvlh_int_pm, &
                            label='Integrated (pointmass gravity)', &
                            linestyle='g-',markersize=5,linewidth=2,istat=istat)
    call plt%add_plot(x_lvlh_cw,z_lvlh_cw, &
                            label='CW Equations', &
                            linestyle='r-',markersize=5,linewidth=2,istat=istat)
    call plt%savefig('cw_test_xz.png',istat=istat)
    call plt%destroy()

    call plt%initialize(grid=.true.,xlabel='LVLH y [km]',ylabel='LVLH z [km]',&
                            title='Integration vs CW Propagation',legend=.true.,figsize=[10,5])
    call plt%add_plot([zero],[zero],label='Target',linestyle='bo',markersize=5,istat=istat)
    call plt%add_plot(y_lvlh_int,z_lvlh_int, &
                            label='Integrated (10x10 gravity)', &
                            linestyle='b-',markersize=5,linewidth=2,istat=istat)
    call plt%add_plot(y_lvlh_int_pm,z_lvlh_int_pm,&
                            label='Integrated (pointmass gravity)', &
                            linestyle='g-',markersize=5,linewidth=2,istat=istat)
    call plt%add_plot(y_lvlh_cw,z_lvlh_cw, &
                            label='CW Equations', &
                            linestyle='r-',markersize=5,linewidth=2,istat=istat)
    call plt%savefig('cw_test_yz.png',istat=istat)
    call plt%destroy()

    ! 3d plot:
    call plt%initialize(grid=.true.,xlabel='LVLH x [km]',ylabel='LVLH y [km]',zlabel='LVLH z[km]',&
                            title='Integration vs CW Propagation',legend=.true.,mplot3d=.true.,figsize=[10,5])
    call plt%add_3d_plot([zero],[zero],[zero],&
                            label='Target',linestyle='bo',markersize=5,istat=istat)
    call plt%add_3d_plot(x_lvlh_int,y_lvlh_int,z_lvlh_int,&
                            label='Integrated (10x10 gravity)', &
                            linestyle='b-',markersize=5,linewidth=2,istat=istat)
    call plt%add_3d_plot(x_lvlh_int_pm,y_lvlh_int_pm,z_lvlh_int_pm,&
                            label='Integrated (pointmass gravity)', &
                            linestyle='g-',markersize=5,linewidth=2,istat=istat)
    call plt%add_3d_plot(x_lvlh_cw,y_lvlh_cw,z_lvlh_cw, &
                            label='CW Equations', &
                            linestyle='r-',markersize=5,linewidth=2,istat=istat)
    call plt%savefig('cw_test_xyz.png',istat=istat)
    call plt%destroy()

    !cleanup:
    call s%grav%destroy()
    close(iunit,   iostat=istat)
    close(iunit_cw,iostat=istat)
    close(iunit_pm,iostat=istat)

    contains
!*****************************************************************************************

    !*********************************************************
        subroutine deriv(me,et,x,xdot)

        !! derivative routine for target + chaser
        !!
        !! For the spherical harmonic mode, the gravity frame is
        !! assumed to be the IAU_EARTH frame.

        implicit none

        class(rk_class),intent(inout)         :: me
        real(wp),intent(in)                   :: et    !! ephemeris time [sec]
        real(wp),dimension(me%n),intent(in)   :: x     !! [target,chaser] state vectors
        real(wp),dimension(me%n),intent(out)  :: xdot  !! [target,chaser] state dot vectors

        real(wp),dimension(3)   :: r,rb,v,a_grav
        real(wp),dimension(3,3) :: rotmat
        real(wp) :: rmag
        integer :: i

        select type (me)
        class is (spacecraft)

            !rotation matrix from inertial to body-fixed:
            rotmat = icrf_to_iau_earth(et)

            do i=1,2

                !get state:
                select case (i)
                case(1)    !target
                    r = x(1:3)
                    v = x(4:6)
                case(2) !chaser
                    r = x(7:9)
                    v = x(10:12)
                end select

                if (me%use_geopotential) then  !! spherical harmonics

                    !input state is inertial frame, have to convert
                    ! to body-fixed Earth frame:
                    rb = matmul(rotmat,r)    !r in body-fixed frame

                    !get the acceleration due to the geopotential:
                    call me%grav%get_acc(rb,grav_n,grav_m,a_grav)

                    !convert acceleration back to inertial frame:
                    a_grav = matmul(transpose(rotmat),a_grav)

                else

                    rmag = norm2(r)
                    a_grav = -mu_earth/rmag**3 * r ! pointmass gravity

                end if

                !derivative vector:
                 select case (i)
                case(1)    !target
                    xdot(1:3) = v
                    xdot(4:6) = a_grav
                case(2) !chaser
                    xdot(7:9) = v
                    xdot(10:12) = a_grav
                end select

            end do

        end select

        end subroutine deriv
    !*********************************************************

    !*********************************************************
        subroutine integrator_report(me,t,x)

        !! report function for the integrator

        implicit none

        class(rk_class),intent(inout)        :: me
        real(wp),intent(in)                  :: t  !! time [sec]
        real(wp),dimension(me%n),intent(in)  :: x  !! [target,chaser] state vectors

        real(wp),dimension(6) :: rv_chaser_lvlh  !! chaser [r,v] state vector in LVLH relative to target
        integer :: iu

        select type (me)
        class is (spacecraft)

            if (me%use_geopotential) then
                iu = iunit     ! geopotential model output file
            else
                iu = iunit_pm  ! pointmass model output file
            end if

            if (me%first) then  !print header
                write(iu,*) ''
                write(iu,'(*(A25,1X))')  't (sec)',&
                                         'x_TARGET_ECI (km)',    'y_TARGET_ECI (km)',    'z_TARGET_ECI (km)',   &
                                         'vx_TARGET_ECI (km/s)', 'vy_TARGET_ECI (km/s)', 'vz_TARGET_ECI (km/s)',&
                                         'x_CHASER_ECI (km)',    'y_CHASER_ECI (km)',    'z_CHASER_ECI (km)',   &
                                         'vx_CHASER_ECI (km/s)', 'vy_CHASER_ECI (km/s)', 'vz_CHASER_ECI (km/s)',&
                                         'x_CHASER_LVLH (km)',   'y_CHASER_LVLH (km)',   'z_CHASER_LVLH (km)',  &
                                         'vx_CHASER_LVLH (km/s)','vy_CHASER_LVLH (km/s)','vz_CHASER_LVLH (km/s)'
                me%first = .false.
            end if

            !chaser in LVLH relative to target:

            call from_ijk_to_lvlh(rt_ijk  = x(1:3),&
                                  vt_ijk  = x(4:6),&
                                  r_ijk   = x(7:9),&
                                  v_ijk   = x(10:12),&
                                  dr_lvlh = rv_chaser_lvlh(1:3),&
                                  dv_lvlh = rv_chaser_lvlh(4:6) )

            !output states:
            write(iu,'(*(F25.6,1X))') t,x,rv_chaser_lvlh

            !also accumulate the LVLH vectors for plotting:
            if (me%use_geopotential) then
                if (allocated(x_lvlh_int)) then
                    x_lvlh_int = [x_lvlh_int, rv_chaser_lvlh(1)]
                    y_lvlh_int = [y_lvlh_int, rv_chaser_lvlh(2)]
                    z_lvlh_int = [z_lvlh_int, rv_chaser_lvlh(3)]
                else
                    x_lvlh_int = [rv_chaser_lvlh(1)]
                    y_lvlh_int = [rv_chaser_lvlh(2)]
                    z_lvlh_int = [rv_chaser_lvlh(3)]
                end if
            else
                if (allocated(x_lvlh_int_pm)) then
                    x_lvlh_int_pm = [x_lvlh_int_pm, rv_chaser_lvlh(1)]
                    y_lvlh_int_pm = [y_lvlh_int_pm, rv_chaser_lvlh(2)]
                    z_lvlh_int_pm = [z_lvlh_int_pm, rv_chaser_lvlh(3)]
                else
                    x_lvlh_int_pm = [rv_chaser_lvlh(1)]
                    y_lvlh_int_pm = [rv_chaser_lvlh(2)]
                    z_lvlh_int_pm = [rv_chaser_lvlh(3)]
                end if
            end if
        end select

        end subroutine integrator_report
    !*********************************************************

    !*********************************************************
        subroutine cw_report(t,x)

        !! report function for the CW propagator.

        implicit none

        real(wp),intent(in)               :: t   !! time [sec]
        real(wp),dimension(6),intent(in)  :: x   !! chaser RSW state vector [km,km/s]

        real(wp),dimension(6) :: rv_chaser_lvlh  !! chaser [r,v] state vector in LVLH relative to target

        if (.not. cw_header_written) then  !print header
            write(iunit_cw,*) ''
            write(iunit_cw,'(*(A25,1X))') 't (sec)',&
                                          'x_CHASER_LVLH (km)',   'y_CHASER_LVLH (km)',   'z_CHASER_LVLH (km)',  &
                                          'vx_CHASER_LVLH (km/s)','vy_CHASER_LVLH (km/s)','vz_CHASER_LVLH (km/s)'
            cw_header_written = .true.
        end if

        !convert state to LVLH:
        call from_rsw_to_lvlh(x(1:3),x(4:6),rv_chaser_lvlh(1:3),rv_chaser_lvlh(4:6))

        !output states:
        write(iunit_cw,'(*(F25.6,1X))') t,rv_chaser_lvlh

        !also accumulate the LVLH vectors for plotting:
        if (allocated(x_lvlh_cw)) then
            x_lvlh_cw = [x_lvlh_cw, rv_chaser_lvlh(1)]
            y_lvlh_cw = [y_lvlh_cw, rv_chaser_lvlh(2)]
            z_lvlh_cw = [z_lvlh_cw, rv_chaser_lvlh(3)]
        else
            x_lvlh_cw = [rv_chaser_lvlh(1)]
            y_lvlh_cw = [rv_chaser_lvlh(2)]
            z_lvlh_cw = [rv_chaser_lvlh(3)]
        end if

        end subroutine cw_report
    !*********************************************************

    end program cw_test
!*****************************************************************************************
