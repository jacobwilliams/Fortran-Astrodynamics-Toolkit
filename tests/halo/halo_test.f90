!*****************************************************************************************
!> author: Jacob Williams
!
!  Halo Orbit generation example (Earth-Moon system).
!
!@note Requires pyplot-fortran module.

    program halo_test

    use fortran_astrodynamics_toolkit, wp => fat_wp
    use pyplot_module
    use halo_orbit_module

    implicit none

    real(wp),parameter :: mu_earth = 398600.436233_wp    !! \( \mu_{Earth} ~ (\mathrm{km}^3/\mathrm{s}^2) \)
    real(wp),parameter :: mu_moon  = 4902.800076_wp      !! \( \mu_{Moon}  ~ (\mathrm{km}^3/\mathrm{s}^2) \)
    real(wp),parameter :: em_distance = 384400.0_wp      !! earth-moon distance [km]

    real(wp),parameter :: A_z = 1000.0_wp   !! halo z-amplitude
    integer,parameter  :: sol_family = 1    !! 1 or 3 (family)
    real(wp),parameter :: t1 = 0.0_wp       !! tau for halo orbit

    integer,parameter  :: n  = 6            !! number of state variables
    real(wp),parameter :: t0 = 0.0_wp       !! initial time (normalized)
    real(wp),parameter :: tmax = 1000.0_wp  !! max final time (normalized)
    real(wp),parameter :: dt = 0.001_wp     !! time step (normalized)

    real(wp),parameter :: tol    = 1.0e-8_wp !! tolerance for event finding
    real(wp),parameter :: xtol   = 1.0e-6_wp !! tolerance for [[hybrd]]
    integer,parameter  :: maxfev = 1000      !! max number of function evaluations for [[hybrd]]

    !< initial state (normalized)
    real(wp),dimension(:),allocatable   :: x_crtbp,y_crtbp,z_crtbp
    real(wp),dimension(:),allocatable   :: x_vec
    real(wp),dimension(:),allocatable   :: y_vec
    real(wp),dimension(:,:),allocatable :: c_mat
    real(wp),dimension(2)   :: xy_L4,xy_L5
    real(wp),dimension(2)   :: x_vy0,vx_vzf
    real(wp),dimension(n)   :: x0
    real(wp)                :: mu    !! CRTPB parameter
    real(wp)                :: c     !! CRTPB Jacobi constant
    type(rk8_10_class)      :: prop  !! integrator
    real(wp),dimension(n)   :: xf    !! final state
    type(pyplot)            :: plt   !! for making the plot
    real(wp),dimension(6)   :: x0_km,x1
    real(wp),dimension(6)   :: x_zvc
    integer                 :: i,j,nx,ny
    real(wp)                :: x,y,gf,tf_actual
    real(wp)                :: x_L1,x_L2,x_L3
    integer                 :: info
    real(wp)                :: Az
    integer                 :: libpoint !! libration point
    real(wp)                :: period   !! halo orbit period
    integer                 :: istat

    character(len=1),dimension(4),parameter :: colors = ['r','g','b','k']  !! line colors for plots

    !for zero-velocity countours:
    real(wp),parameter  :: xmin  = -2.0_wp
    real(wp),parameter  :: xmax  = 2.0_wp
    real(wp),parameter  :: xstep = 0.01_wp
    real(wp),parameter  :: ymin  = -2.0_wp
    real(wp),parameter  :: ymax  = 2.0_wp
    real(wp),parameter  :: ystep = 0.01_wp

    ! get the halo orbit initial guess using the analytic approximation:
    libpoint = 2 ! L2 point
    call halo_to_rv(libpoint,mu_earth,mu_moon,em_distance,A_z,sol_family,t1,x0)

    call unnormalize_variables(mu_earth,mu_moon,em_distance,x_crtbp=x0,x=x0_km)

    write(*,*) ''
    write(*,*) ' initial guess in km,km/s: '
    write(*,'(*(F30.16,1X))') x0_km
    write(*,*) ''

    !compute the CRTBP parameter & Jacobi constant:
    mu = compute_crtpb_parameter(mu_earth,mu_moon)
    c  = compute_jacobi_constant(mu,x0)

    ! !*************************************************
    ! !compute the zero velocity surface for this c

    !Just get the number of elements by counting...
    i=0
    nx=0
    do
        i=i+1
        x = xmin + (i-1)*xstep
        if (x>xmax) exit
        nx = nx + 1
    end do
    j=0
    ny=0
    do
        j=j+1
        y = ymin + (j-1)*ystep
        if (y>ymax) exit
        ny = ny + 1
    end do

    !size the arrays:
    allocate(x_vec(nx))
    allocate(y_vec(ny))
    allocate(c_mat(nx,ny))

    !compute the jacobi constant for each combo:
    do i = 1, nx
        x = xmin + (i-1)*xstep
        x_vec(i) = x
        do j = 1, ny
            y = ymin + (j-1)*ystep
            if (i==1) y_vec(j) = y
            x_zvc(1:3) = [x,y,zero]
            x_zvc(4:6) = zero
            c_mat(i,j) = compute_jacobi_constant(mu,x_zvc)
        end do
    end do

    !*************************************************
    !compute the libration point locations:

    call compute_libration_points(mu,x_L1,x_L2,x_L3,xy_L4,xy_L5)

    write(*,*) ''
    write(*,*) 'mu:  ', mu
    write(*,*) 'Jacobi constant:  ', c
    write(*,*) 'x_L1:', x_L1
    write(*,*) 'x_L2:', x_L2
    write(*,*) 'x_L3:', x_L3
    write(*,*) ''
    write(*,*) 'x_moon :', 1.0_wp - mu
    write(*,*) 'x_earth:', -mu

    !*************************************************
    !integrate initial guess:

    write(*,*) ''
    write(*,*) ' Halo Initial Guess'
    write(*,*) ''
    write(*,*) 't0=',t0
    write(*,*) 'r0=',x0(1:3)
    write(*,*) 'v0=',x0(4:6)

    call prop%initialize(n,func,report,x_axis_crossing)
    call prop%integrate_to_event(t0,x0,dt,tmax,tol,tf_actual,xf,gf)

    write(*,*) ''
    write(*,*) 'tf=',tf_actual
    write(*,*) 'rf=',xf(1:3)
    write(*,*) 'vf=',xf(4:6)

    !plot the 2D trajectory, zero-velocity curves, and libration point locations:
    call plt%initialize(grid=.true.,xlabel='x',ylabel='y',&
                            title='Halo Initial Guess',legend=.false.,figsize=[10,10],&
                            use_numpy=.true.,axis_equal=.true.,mplot3d=.true.)

    !trajectory:
    call plt%add_3d_plot(x_crtbp,y_crtbp,z_crtbp,label='trajectory',linestyle='b-',linewidth=2,istat=istat)
    !bodies:
    call plt%add_3d_plot([-mu],[zero],[zero], label='B1',linestyle='bo',markersize=5,linewidth=5,istat=istat)
    call plt%add_3d_plot([1.0_wp-mu],[zero],[zero], label='B2',linestyle='ko',markersize=5,linewidth=5,istat=istat)
    !libration point locations:
    call plt%add_3d_plot([x_L1],[zero],[zero], label='L1',linestyle='rx',markersize=3,linewidth=3,istat=istat)
    call plt%add_3d_plot([x_L2],[zero],[zero], label='L2',linestyle='rx',markersize=3,linewidth=3,istat=istat)
    call plt%add_3d_plot([x_L3],[zero],[zero], label='L3',linestyle='rx',markersize=3,linewidth=3,istat=istat)
    !zero-velocity curve (for this jacobi constant):
    call plt%add_contour(x_vec, y_vec, c_mat, &
                           linestyle='-', linewidth=2, &
                           levels=[c], color='r',istat=istat)
    call plt%savefig('halo_guess.png',pyfile='halo_guess.py',istat=istat)
    call plt%destroy()

    !-------------------------------------------
    !
    !  solve for a Halo using HYBRD
    !

    write(*,*) ''
    write(*,*) 'solve for halo...'

    !turn off integration reporting for iterations:
    call prop%destroy()
    call prop%initialize(n,func,g=x_axis_crossing)

    !initialize the plot:
    call plt%initialize(grid=.true.,xlabel='x',ylabel='y',&
                            title='Earth-Moon Halo Orbits',&
                            legend=.false.,figsize=[10,10],&
                            use_numpy=.true.,axis_equal=.true.,&
                            mplot3d=.true.)

    ! see also: http://ccar.colorado.edu/asen5050/projects/projects_2012/bezrouk/Mission%20Analysis.html

    do libpoint=1,2 ! libration point L1, L2, L3
        do i=1,10

            Az = A_z+(i-1)*1000.0_wp  ! in km
            !Az = A_z+(i-1)*10000.0_wp  ! in km  - for L3 cases

            call prop%destroy()
            call prop%initialize(n,func,g=x_axis_crossing)

            call halo_to_rv(libpoint,mu_earth,mu_moon,em_distance,&
                            Az,sol_family,t1,x0,period=period)
            if (period==zero) cycle

            write(*,*) ''
            write(*,*) ' r0=',x0(1:3)
            write(*,*) ' v0=',x0(4:6)

            if (allocated(x_crtbp)) then ! clear these for new trajectory
                deallocate(x_crtbp)
                deallocate(y_crtbp)
                deallocate(z_crtbp)
            end if

            !now, solve for a halo:
            ! [note: we could also use the STM for better derivatives here....]
            x_vy0 = [x0(1),x0(5)]
            call hybrd1(halo_fcn,2,x_vy0,vx_vzf,tol=xtol,info=info)
            write(*,*) ''
            write(*,*) ' info=',info
            write(*,*) ' Az =',Az
            write(*,*) ' x0 =',x_vy0(1)
            write(*,*) ' vy0=',x_vy0(2)
            write(*,*) ' vxf=',vx_vzf(1)
            write(*,*) ' vzf=',vx_vzf(2)

            !now plot solution:
            call prop%initialize(n,func,report,g=x_axis_crossing)
            x1    = x0
            x1(1) = x_vy0(1)  !solution from hybrd
            x1(5) = x_vy0(2)  !solution from hybrd
            !integrate one full rev (two x-axis crossings):
            call prop%integrate_to_event(t0,x1,dt,tmax,tol,tf_actual,xf,gf)  !1/2 rev
            call prop%integrate_to_event(t0,xf,dt,tmax,tol,tf_actual,x1,gf)  !1 rev

            write(*,*) 'period: ', 2.0_wp * tf_actual

            !plot the 3D trajectory:
            call plt%add_3d_plot(x_crtbp,y_crtbp,z_crtbp,label='solution',&
                                linestyle=colors(1+mod(i,size(colors)))//'-',linewidth=2,istat=istat)

        end do
    end do

    !also plot the libration point locations:
    call plt%add_3d_plot([x_L1],[zero],[zero], label='L1',linestyle='rx',markersize=4,linewidth=3,istat=istat)
    call plt%add_3d_plot([x_L2],[zero],[zero], label='L2',linestyle='rx',markersize=4,linewidth=3,istat=istat)
    !call plt%add_3d_plot([x_L3],[zero],[zero], label='L3',linestyle='rx',markersize=3,linewidth=3,istat=istat)
    !call plt%add_3d_plot([-mu],[zero],[zero], label='B1',linestyle='bo',markersize=5,linewidth=5,istat=istat)
    call plt%add_3d_plot([1.0_wp-mu],[zero],[zero], label='B2',linestyle='ko',markersize=5,linewidth=5,istat=istat)
    call plt%savefig('halos.png',pyfile='halos.py',istat=istat)
    call plt%destroy()

    contains

!**************************************************************************
    subroutine halo_fcn(n,xvec,fvec,iflag)  !! Halo function for [[hybrd1]]

    implicit none

    integer,intent(in)                :: n      !! `n=2` in this case
    real(wp),dimension(n),intent(in)  :: xvec   !! x_vy0
    real(wp),dimension(n),intent(out) :: fvec   !! [vxf,vzf]
    integer,intent(inout)             :: iflag

    real(wp) :: gf,t0
    real(wp),dimension(6) :: x,x1,xf

    real(wp),parameter :: tol = 1.0e-8_wp !! event finding tolerance

    t0   = zero    ! epoch doesn't matter for crtbp
    x    = x0      ! initial guess state
    x(1) = xvec(1) ! x0
    x(5) = xvec(2) ! vy0

    !integrate to the next x-axis crossing:
    call prop%integrate_to_event(t0,x,dt,tmax,tol,tf_actual,xf,gf)

    !want x and z-velocity at the x-axis crossing to be zero:
    fvec = [xf(4),xf(6)]

    end subroutine halo_fcn
!**************************************************************************

!**************************************************************************
    subroutine func(me,t,x,xdot)  !! CRTBP derivative function
    implicit none
    class(rk_class),intent(inout)        :: me
    real(wp),intent(in)                  :: t
    real(wp),dimension(me%n),intent(in)  :: x
    real(wp),dimension(me%n),intent(out) :: xdot

    call crtbp_derivs(mu,x,xdot)

    end subroutine func
!**************************************************************************

!**************************************************************************
    subroutine x_axis_crossing(me,t,x,g)  !! x-axis crossing event function
    implicit none
    class(rk_class),intent(inout)        :: me
    real(wp),intent(in)                  :: t
    real(wp),dimension(me%n),intent(in)  :: x
    real(wp),intent(out)                 :: g

    g = x(2)  ! y = 0 at x-z-plane crossing

    end subroutine x_axis_crossing
!**************************************************************************

!**************************************************************************
    subroutine report(me,t,x)  !! report function

    implicit none

    class(rk_class),intent(inout)       :: me
    real(wp),intent(in)                 :: t
    real(wp),dimension(me%n),intent(in) :: x

    !write(*,'(*(F30.16,1X))') t, x, compute_jacobi_constant(mu,x)

    if (allocated(x_crtbp)) then   ! uses Fortran 2008 auto LHS reallocations
        x_crtbp  = [x_crtbp, x(1)]
        y_crtbp  = [y_crtbp, x(2)]
        z_crtbp  = [z_crtbp, x(3)]
    else
        x_crtbp  = [x(1)]
        y_crtbp  = [x(2)]
        z_crtbp  = [x(3)]
    end if

    end subroutine report
!**************************************************************************

    end program halo_test
!*****************************************************************************************
