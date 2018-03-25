!*****************************************************************************************
!> author: Jacob Williams
!
!  Halo Orbit generation example (Earth-Moon system).
!
!@note Requires pyplot-fortran module.

    program halo_test_2

    use fortran_astrodynamics_toolkit
    use pyplot_module
    use halo_orbit_module

    implicit none

    real(wp),parameter :: mu_earth = 398600.436233_wp    !! \( \mu_{Earth} ~ (\mathrm{km}^3/\mathrm{s}^2) \)
    real(wp),parameter :: mu_moon  = 4902.800076_wp      !! \( \mu_{Moon}  ~ (\mathrm{km}^3/\mathrm{s}^2) \)
    real(wp),parameter :: em_distance = 384400.0_wp      !! earth-moon distance [km]
    real(wp),parameter :: A_z = 1000.0_wp    !! initial halo z-amplitude
    real(wp),parameter :: t1 = 0.0_wp        !! tau for halo orbit
    integer,parameter  :: n  = 6             !! number of state variables
    real(wp),parameter :: Az_step =1000.0_wp !! Az step size (km)

    real(wp),dimension(n) :: rv  !! halo orbit initial state
    integer  :: libpoint         !! libration point (1,2,3)
    integer  :: i                !! counter
    real(wp) :: Az               !! halo z-amplitude (km)
    integer  :: info             !! minpack status code
    integer  :: sol_family       !! 1 or 3 (family)
    real(wp) :: period           !! halo period (normalized time)

    do libpoint = 1, 3
        write(*,*) ''
        write(*,*) '---------------------------'
        write(*,*) ' libpoint: ', libpoint
        write(*,*) '---------------------------'
        write(*,*) ''
        do sol_family = 1, 3, 2
            write(*,*) ''
            write(*,*) ' family: ', sol_family
            write(*,*) ''
            do i=1,10
                Az = A_z+(i-1)*Az_step
                call halo_to_rv_diffcorr(libpoint,mu_earth,mu_moon,em_distance,&
                                            Az,sol_family,t1,rv,info,period)
                if (info/=1) error stop 'error'
                write(*,*) ''
                write(*,*) az,' [', rv(1),',',&
                                    rv(2),',',&
                                    rv(3),',',&
                                    rv(4),',',&
                                    rv(5),',',&
                                    rv(6),']', period
            end do
        end do
    end do

    end program halo_test_2
!*****************************************************************************************
