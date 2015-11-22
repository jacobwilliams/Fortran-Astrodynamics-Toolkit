!*****************************************************************************************
!> author: Jacob Williams
!  date: 11/21/2015
!
!  Test program for the CRTBP routines.

    program crtbp_propagation_test
    
    use fortran_astrodynamics_toolkit
    use pyplot_module
        
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
    
   ! real(wp),dimension(n),parameter :: x0 = [ 0.5_wp, &
   !                                           1.0_wp, &
   !                                           0.0_wp,              &
   !                                           0.5_wp,              &
   !                                           0.5_wp,              &
   !                                           0.0_wp               ]
    
    real(wp),dimension(:),allocatable :: x_crtbp,y_crtbp,z_crtbp
    real(wp)              :: mu    !! CRTPB parameter
    type(rk8_10_class)    :: prop  !! integrator
    real(wp),dimension(n) :: xf    !! final state
    type(pyplot)          :: plt   !! for making the plot
    real(wp),dimension(6) :: x0_km
    
    call unnormalize_state(x0,mu_earth,mu_moon,384400.0_wp,x0_km)
    
    write(*,*) ''
    write(*,*) ' initial state in km,km/s: '
    write(*,'(*(F30.16,1X))') x0_km
    write(*,*) ''
    
    !compute the CRTBP parameter:
    mu = compute_crtpb_parameter(mu_earth,mu_moon)
    
    !integrate:
    call prop%initialize(n,func,report) 
    call prop%integrate(t0,x0,dt,tf,xf)
    
    !plot the trajectory (2D):
    call plt%initialize(grid=.true.,xlabel='x [km]',ylabel='y [km]',&
                            title='CRTBP',legend=.false.,figsize=[10,5])
    call plt%add_plot(x_crtbp,y_crtbp,label='trajectory',linestyle='b-',linewidth=2)
    call plt%savefig('crtbp_test.png')
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
    
    !write(*,'(*(F30.16,1X))') t, x
    
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