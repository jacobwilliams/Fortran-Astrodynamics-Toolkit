!*****************************************************************************************
!> author: Jacob Williams
!  date: 11/21/2015
!
!  Test program for the CRTBP routines.

    program crtbp_propagation_test

    use fortran_astrodynamics_toolkit, wp => fat_wp
    use pyplot_module, only: pyplot

    implicit none

    real(wp),parameter :: mu_earth = 398600.436233_wp    !! \( \mu_{Earth} ~ (\mathrm{km}^3/\mathrm{s}^2) \)
    real(wp),parameter :: mu_moon  = 4902.800076_wp      !! \( \mu_{Moon}  ~ (\mathrm{km}^3/\mathrm{s}^2) \)
    real(wp),parameter :: mu_sun   = 132712440017.987_wp !! \( \mu_{Sun}   ~ (\mathrm{km}^3/\mathrm{s}^2) \)


    integer,parameter  :: n  = 6         !! number of state variables
    real(wp),parameter :: t0 = 0.0_wp    !! initial time (normalized)
    real(wp),parameter :: tf = 50.0_wp   !! final time (normalized)
    real(wp),parameter :: dt = 0.01_wp   !! time step (normalized)

    !< initial state (normalized)
    !< see: [Celestial Mechanics Notes Set 4: The Circular Restricted Three Body Problem](http://cosweb1.fau.edu/~jmirelesjames/hw4Notes.pdf), p.40.
    real(wp),dimension(n),parameter :: x0 = [ 0.30910452642073_wp, &
                                              0.07738174525518_wp, &
                                              0.0_wp,              &
                                             -0.72560796964234_wp, &
                                              1.55464233412773_wp, &
                                              0.0_wp               ]

    real(wp),dimension(:),allocatable :: x_crtbp,y_crtbp,z_crtbp
    real(wp)              :: mu    !! CRTPB parameter
    real(wp)              :: c     !! CRTPB Jacobi constant
    type(rk8_10_class)    :: prop  !! integrator
    real(wp),dimension(n) :: xf    !! final state
    type(pyplot)          :: plt   !! for making the plot
    real(wp),dimension(6) :: x0_km
    real(wp),dimension(6) :: x_zvc
    integer :: i,j,nx,ny
    real(wp) :: x,y
    real(wp),dimension(:),allocatable   :: x_vec
    real(wp),dimension(:),allocatable   :: y_vec
    real(wp),dimension(:,:),allocatable :: c_mat
    real(wp) :: x_L1,x_L2,x_L3
    real(wp),dimension(2) :: xy_L4,xy_L5
    integer :: istat

    !for zero-velocity countours:
    real(wp),parameter  :: xmin  = -2.0_wp
    real(wp),parameter  :: xmax  = 2.0_wp
    real(wp),parameter  :: xstep = 0.01_wp
    real(wp),parameter  :: ymin  = -2.0_wp
    real(wp),parameter  :: ymax  = 2.0_wp
    real(wp),parameter  :: ystep = 0.01_wp

    call unnormalize_variables(mu_earth,mu_moon,384400.0_wp,x_crtbp=x0,x=x0_km)

    write(*,*) ''
    write(*,*) ' initial state in km,km/s: '
    write(*,'(*(F30.16,1X))') x0_km
    write(*,*) ''

    !compute the CRTBP parameter & Jacobi constant:
    mu = compute_crtpb_parameter(mu_earth,mu_moon)
    c  = compute_jacobi_constant(mu,x0)

    !*************************************************
    !compute the zero velocity surface for the c

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
    write(*,*) 'x_L1:', x_L1
    write(*,*) 'x_L2:', x_L2
    write(*,*) 'x_L3:', x_L3
    write(*,*) ''

    !*************************************************
    !integrate:
    call prop%initialize(n,func,report)
    call prop%integrate(t0,x0,dt,tf,xf)

    !plot the 2D trajectory, zero-velocity curves, and libration point locations:
    call plt%initialize(grid=.true.,xlabel='x [km]',ylabel='y [km]',&
                            title='CRTBP Example',legend=.false.,figsize=[10,10],&
                            use_numpy=.true.,axis_equal=.true.)

    !trajectory:
    call plt%add_plot(x_crtbp,y_crtbp,label='trajectory',linestyle='b-',linewidth=2,istat=istat)
    !libration point locations:
    call plt%add_plot([x_L1],[zero],        label='L1',linestyle='rx',markersize=3,linewidth=3,istat=istat)
    call plt%add_plot([x_L2],[zero],        label='L2',linestyle='rx',markersize=3,linewidth=3,istat=istat)
    call plt%add_plot([x_L3],[zero],        label='L3',linestyle='rx',markersize=3,linewidth=3,istat=istat)
    call plt%add_plot([xy_L4(1)],[xy_L4(2)],label='L4',linestyle='rx',markersize=3,linewidth=3,istat=istat)
    call plt%add_plot([xy_L5(1)],[xy_L5(2)],label='L5',linestyle='rx',markersize=3,linewidth=3,istat=istat)
    !zero-velocity curve (for this jacobi constant):
    call plt%add_contour(x_vec, y_vec, c_mat, &
                            linestyle='-', linewidth=2, levels=[c], color='r',istat=istat)
    call plt%savefig('crtbp_test.png',istat=istat)
    call plt%destroy()

    contains

    subroutine func(me,t,x,xdot)  !! CRTBP derivative function
    implicit none
    class(rk_class),intent(inout)        :: me
    real(wp),intent(in)                  :: t
    real(wp),dimension(me%n),intent(in)  :: x
    real(wp),dimension(me%n),intent(out) :: xdot

    call crtbp_derivs(mu,x,xdot)

    end subroutine func

    subroutine report(me,t,x)  !! report function
    implicit none
    class(rk_class),intent(inout)       :: me
    real(wp),intent(in)                 :: t
    real(wp),dimension(me%n),intent(in) :: x

    !write(*,'(*(F30.16,1X))') t, x, compute_jacobi_constant(mu,x)

    if (allocated(x_crtbp)) then
        x_crtbp  = [x_crtbp, x(1)]
        y_crtbp  = [y_crtbp, x(2)]
        z_crtbp  = [z_crtbp, x(3)]
    else
        x_crtbp  = [x(1)]
        y_crtbp  = [x(2)]
        z_crtbp  = [x(3)]
    end if

    end subroutine report

    end program crtbp_propagation_test
!*****************************************************************************************
