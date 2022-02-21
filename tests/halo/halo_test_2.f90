!*****************************************************************************************
!> author: Jacob Williams
!
!  Halo Orbit generation example (Earth-Moon system).
!
!@note Requires pyplot-fortran module.

    program halo_test_2

    use fortran_astrodynamics_toolkit, wp => fat_wp
    use pyplot_module

    implicit none

    real(wp),parameter :: mu_earth = 398600.436233_wp    !! \( \mu_{Earth} ~ (\mathrm{km}^3/\mathrm{s}^2) \)
    real(wp),parameter :: mu_moon  = 4902.800076_wp      !! \( \mu_{Moon}  ~ (\mathrm{km}^3/\mathrm{s}^2) \)
    real(wp),parameter :: em_distance = 384400.0_wp      !! earth-moon distance [km]
    real(wp),parameter :: A_z = 1000.0_wp    !! initial halo z-amplitude
    real(wp),parameter :: t1 = 0.0_wp        !! tau for halo orbit
    integer,parameter  :: n  = 6             !! number of state variables
    real(wp),parameter :: Az_step =1000.0_wp !! Az step size (km)

    real(wp),dimension(n)    :: rv     !! halo orbit initial state
    real(wp),dimension(n,n)  :: phi    !! monodromy matrix
    real(wp),dimension(n,2)  :: w      !! real and and imaginary parts of the eigenvalues of `phi`
    real(wp),dimension(n,n)  :: z      !! real and imaginary parts of the eigenvectors of `phi`
    complex(wp),dimension(n) :: lambda !! eigenvalues of `phi`
    integer  :: libpoint         !! libration point (1,2,3)
    integer  :: i                !! counter
    integer  :: j                !! counter
    real(wp) :: Az               !! halo z-amplitude (km)
    integer  :: info             !! minpack status code
    integer  :: sol_family       !! 1 or 3 (family)
    real(wp) :: period           !! halo period (normalized time)
    integer  :: ierr             !! output flag from [[compute_eigenvalues_and_eigenvectors]]
    real(wp) :: mu               !! CRTBP parameter
    real(wp),dimension(:),allocatable   :: e !! real eigenvalues
    real(wp),dimension(:,:),allocatable :: v !! normalized eigenvectors associated with the real eigenvalues
    integer :: n_results !! number of real eigenvalues

    mu = compute_crtpb_parameter(mu_earth,mu_moon)

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

                if (libpoint==1 .and. sol_family==1 .and. i==1) then

                    call compute_halo_monodromy_matrix(mu,rv,period,phi)
                    write(*,*) ''
                    write(*,*) 'monodromy matrix:'
                    write(*,*) ''
                    call print_matrix(phi)
                    write(*,*) ''

                    ! general method:
                    call compute_eigenvalues_and_eigenvectors(6, phi, w, z, ierr)
                    write(*,*) ''
                    write(*,*) 'eigenvalues:'
                    write(*,*) ''
                    call print_matrix(w)
                    write(*,*) ''
                    write(*,*) 'eigenvectors:'
                    write(*,*) ''
                    call print_matrix(z)
                    write(*,*) ''

                    ! quicker method for just the eigenvalues:
                    call compute_monodromy_matrix_eigenvalues(phi,lambda)
                    write(*,*) ''
                    write(*,*) 'lambda:'
                    write(*,*) lambda
                    write(*,*) ''

                    ! normalized values:
                    call compute_real_eigenvalues_and_normalized_eigenvectors(6, phi, e, v, n_results, ierr)
                    write(*,*) ''
                    write(*,*) 'real eigenvalues:', e
                    write(*,*) ''
                    write(*,*) 'normalized eigenvectors:'
                    write(*,*) ''
                    call print_matrix(v)
                    write(*,*) ''

                    ! the determinant should be 1.0 since matrix is sympletic
                    write(*,*) ''
                    write(*,*) '-----'
                    write(*,*) 'determinate of phi: ', matrix_determinant(6,phi)
                    write(*,*) '-----'
                    write(*,*) ''

                end if
            end do
        end do
    end do

    end program halo_test_2
!*****************************************************************************************
