!*****************************************************************************************
!>
!  Test of Fortran Astrodynamics Toolkit - integrate a trajectory near the Earth,
!  using the GGM03C geopotential model.

    program gravity_test

    use fortran_astrodynamics_toolkit, wp => fat_wp

    implicit none

    !spacecraft propagation type:
    ! extend the rk class to include data used in the deriv routine
    type,extends(rk8_10_class) :: spacecraft
        type(geopotential_model_mueller) :: grav   !central body geopotential model
        logical :: first = .true.  !first point is being exported
    end type spacecraft

    integer,parameter  :: n=6    !number of state variables
    real(wp),parameter :: mu_earth  = 398600.4418_wp   !grav. param of earth [km3/s2]
    real(wp),parameter :: orbit_inc = 45.0_wp*deg2rad  !orbit inclination [rad]
    real(wp),parameter :: orbit_sma = 6778.0_wp        !orbit semimajor axis [km]

    !the coefficient file:
    character(len=*),parameter :: gravfile = '../grav/GGM03C.GEO'
    integer,parameter :: grav_n = 10    !max degree
    integer,parameter :: grav_m = 10    !max order

    type(spacecraft) :: s
    real(wp) :: t0,tf,x0(n),dt,xf(n),x02(n)
    logical :: status_ok
    integer :: iunit,istat

    !open output file:
    open(newunit=iunit,file='gravity_test.txt',status='REPLACE',iostat=istat)

    !constructor (main body is Earth):
    call s%initialize(n=n,f=deriv,report=twobody_report)
    call s%grav%initialize(gravfile,grav_n,grav_m,status_ok)
    if (.not. status_ok) stop 'Error'

    !get an initial cartesian state vector:
    ! [circular, 45 deg inclination, radius of 6778 km]
    call els3pv(mu_earth,[mu_earth/orbit_sma,orbit_sma,orbit_inc,zero,zero,zero],x0)

    write(iunit,'(A/,*(E30.16/))') 'Initial state:',x0

    !initial conditions:
    t0 = zero         !initial time (sec)
    dt = 100.0_wp     !time step (sec)
    tf = two*day2sec  !final time (sec)

    !integrate:
    s%first = .true.
    call s%integrate(t0,x0,dt,tf,xf)
    write(iunit,*) ''
    write(iunit,'(A/,*(E30.16/))') 'Final state:',xf

    !cleanup:
    call s%grav%destroy()
    close(iunit,iostat=istat)

!*****************************************************************************************
    contains

    !*********************************************************
        subroutine deriv(me,et,x,xdot)

        !! derivative routine

        implicit none

        class(rk_class),intent(inout)         :: me
        real(wp),intent(in)                   :: et    !! this is ephemeris time
        real(wp),dimension(me%n),intent(in)   :: x
        real(wp),dimension(me%n),intent(out)  :: xdot

        real(wp),dimension(3) :: r,rb,v,a_grav
        real(wp),dimension(3,3) :: rotmat

        select type (me)
        class is (spacecraft)

            !get state:
            r = x(1:3)
            v = x(4:6)

            !input state is inertial frame, have to convert
            ! to body-fixed Earth frame:
            rotmat = icrf_to_iau_earth(et)  !rotation matrix
            rb = matmul(rotmat,r)    !r in body-fixed frame

            !get the acceleration due to the geopotential:
            call me%grav%get_acc(rb,grav_n,grav_m,a_grav)

            !convert acceleration back to inertial frame:
            a_grav = matmul(transpose(rotmat),a_grav)

            !derivative vector:
            xdot(1:3) = v
            xdot(4:6) = a_grav

        end select

        end subroutine deriv
    !*********************************************************

    !*********************************************************
        subroutine twobody_report(me,t,x)

        !! report function - write time,state to console

        implicit none

        class(rk_class),intent(inout)        :: me
        real(wp),intent(in)                  :: t
        real(wp),dimension(me%n),intent(in)  :: x

        select type (me)
        class is (spacecraft)
            if (me%first) then  !print header
                write(iunit,*) ''
                write(iunit,'(*(A15,1X))')  't (sec)','x (km)','y (km)','z (km)',&
                                        'vx (km/s)','vy (km/s)','vz (km/s)'
                me%first = .false.
            end if
        end select

        write(iunit,'(*(F15.6,1X))') t,x

        end subroutine twobody_report
    !*********************************************************

    end program gravity_test
!*****************************************************************************************