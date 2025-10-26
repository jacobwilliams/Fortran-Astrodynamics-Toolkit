!*****************************************************************************************
!> author: Jacob Williams
!  date: 1/29/2016
!
!  Distant Retrograde Orbit generation example (Earth-Moon system).
!
!@note Requires pyplot-fortran module.

    program dro_test

    use fortran_astrodynamics_toolkit, wp => fat_wp
    use pyplot_module

    implicit none

    real(wp),parameter :: mu_earth = 398600.436233_wp    !! \( \mu_{Earth} ~ (\mathrm{km}^3/\mathrm{s}^2) \)
    real(wp),parameter :: mu_moon  = 4902.800076_wp      !! \( \mu_{Moon}  ~ (\mathrm{km}^3/\mathrm{s}^2) \)

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
    real(wp),dimension(1)   :: vy0,vxf
    real(wp),dimension(n)   :: x0
    real(wp)                :: mu    !! CRTPB parameter
    real(wp)                :: c     !! CRTPB Jacobi constant
    type(rk8_10_class)      :: prop  !! integrator
    real(wp),dimension(n)   :: xf    !! final state
    type(pyplot)            :: plt   !! for making the plot
    real(wp),dimension(6)   :: x0_km,x1
    real(wp),dimension(6)   :: x_zvc
    integer                 :: i,j,nx,ny
    real(wp)                :: x,y,gf,tf_actual,t1
    real(wp)                :: x_L1,x_L2,x_L3
    integer                 :: info
    integer                 :: istat

    character(len=1),dimension(4),parameter :: colors = ['r','g','b','k']  !! line colors for plots

    !for zero-velocity countours:
    real(wp),parameter  :: xmin  = -2.0_wp
    real(wp),parameter  :: xmax  = 2.0_wp
    real(wp),parameter  :: xstep = 0.01_wp
    real(wp),parameter  :: ymin  = -2.0_wp
    real(wp),parameter  :: ymax  = 2.0_wp
    real(wp),parameter  :: ystep = 0.01_wp

    !initial guess for DRO state (normalized):
    x0 = [ 1.2_wp, &
           0.0_wp, &
           0.0_wp, &
           0.0_wp, &
           -0.7_wp, &
           0.0_wp ]

    call unnormalize_variables(mu_earth,mu_moon,384400.0_wp,x_crtbp=x0,x=x0_km)

    write(*,*) ''
    write(*,*) ' initial guess in km,km/s: '
    write(*,'(*(F30.16,1X))') x0_km
    write(*,*) ''

    !compute the CRTBP parameter & Jacobi constant:
    mu = compute_crtpb_parameter(mu_earth,mu_moon)
    c  = compute_jacobi_constant(mu,x0)

    !*************************************************
    !compute the zero velocity surface for this c

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

    !*************************************************
    !integrate initial guess:

    write(*,*) ''
    write(*,*) ' DRO Initial Guess'
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
                            title='DRO Initial Guess',legend=.false.,figsize=[10,10],&
                            use_numpy=.true.,axis_equal=.true.)

    !trajectory:
    call plt%add_plot(x_crtbp,y_crtbp,label='trajectory',linestyle='b-',linewidth=2,istat=istat)
    !libration point locations:
    call plt%add_plot([x_L1],[zero],label='L1',linestyle='rx',markersize=3,linewidth=3,istat=istat)
    call plt%add_plot([x_L2],[zero],label='L2',linestyle='rx',markersize=3,linewidth=3,istat=istat)
    !zero-velocity curve (for this jacobi constant):
    call plt%add_contour(x_vec, y_vec, c_mat, &
                            linestyle='-', linewidth=2, &
                            levels=[c], color='r',istat=istat)
    call plt%savefig('dro_guess.png',istat=istat)
    call plt%destroy()

    !-------------------------------------------
    !
    !  solve for a DRO using HYBRD
    !

    write(*,*) ''
    write(*,*) 'solve for DRO...'

    !turn off integration reporting for iterations:
    call prop%destroy()
    call prop%initialize(n,func,g=x_axis_crossing)

    !initialize the plot:
    call plt%initialize(grid=.true.,xlabel='x',ylabel='y',&
                            title='Earth-Moon Distant Retrograde Orbits',&
                            legend=.false.,figsize=[10,10],&
                            use_numpy=.true.,axis_equal=.true.)

    vy0(1) = -0.7_wp

    do i=1,50,3

        call prop%destroy()
        call prop%initialize(n,func,g=x_axis_crossing)

        x0 = [ 1.2_wp + real(i,wp)/100.0_wp, &
               0.0_wp, &
               0.0_wp, &
               0.0_wp, &
               vy0(1), &  !use previous solution as initial guess
               0.0_wp ]

        write(*,*) ''
        write(*,*) ' r0=',x0(1:3)
        write(*,*) ' v0=',x0(4:6)

        if (allocated(x_crtbp)) then ! clear these for new trajectory
            deallocate(x_crtbp)
            deallocate(y_crtbp)
            deallocate(z_crtbp)
        end if

        !now, solve for a DRO:
        vy0 = x0(5)
        call hybrd1(dro_fcn,1,vy0,vxf,tol=1.0e-4_wp,info=info)
        write(*,*) ''
        write(*,*) ' info=',info
        write(*,*) ' vy0=',vy0
        write(*,*) ' vxf=',vxf

        !now plot solution:
        call prop%initialize(n,func,report,g=x_axis_crossing)
        x1    = x0
        x1(5) = vy0(1)  !solution from hybrd
        !integrate one full rev (two x-axis crossings):
        call prop%integrate_to_event(t0,x1,dt,tmax,tol,tf_actual,xf,gf)  !1/2 rev
        call prop%integrate_to_event(t0,xf,dt,tmax,tol,tf_actual,x1,gf)  !1 rev

        !plot the 2D trajectory:
        call plt%add_plot(x_crtbp,y_crtbp,label='solution',&
                            linestyle=colors(1+mod(i,size(colors)))//'-',linewidth=2,istat=istat)

    end do

    !also plot the libration point locations:
    call plt%add_plot([x_L1],[zero], label='L1',linestyle='rx',markersize=4,linewidth=3,istat=istat)
    call plt%add_plot([x_L2],[zero], label='L2',linestyle='rx',markersize=4,linewidth=3,istat=istat)
    call plt%savefig('dros.png',istat=istat)
    call plt%destroy()

    contains

!**************************************************************************
    subroutine dro_fcn(n,xvec,fvec,iflag)  !! DRO function for [[hybrd1]]

    implicit none

    integer,intent(in)                :: n      !! `n=1` in this case
    real(wp),dimension(n),intent(in)  :: xvec   !! vy0
    real(wp),dimension(n),intent(out) :: fvec   !! vxf
    integer,intent(inout)             :: iflag

    real(wp) :: gf,t0
    real(wp),dimension(6) :: x,x1,xf

    real(wp),parameter :: tol = 1.0e-8_wp !! event finding tolerance

    t0   = zero    ! epoch doesn't matter for crtbp
    x    = x0      ! initial guess state
    x(5) = xvec(1) ! vy0 (only optimization variable)

    !integrate to the next x-axis crossing:
    call prop%integrate_to_event(t0,x,dt,tmax,tol,tf_actual,xf,gf)

    !want x-velocity at the x-axis crossing to be zero:
    fvec = xf(4)

    end subroutine dro_fcn
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

    g = x(2)  ! y = 0 at x-axis crossing

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

    end program dro_test
!*****************************************************************************************
