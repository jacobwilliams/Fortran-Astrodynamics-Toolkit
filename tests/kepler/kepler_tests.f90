!*****************************************************************************************
!> author: Jacob Williams
!
!  Propagate a state using the various kepler methods and compare.

    program kepler_tests

    use kind_module,    only: wp
    use gooding_module
    use kepler_module
    use conversion_module

    implicit none

    real(wp),dimension(6),parameter :: rv1 = [20.0e3_wp,20.0e3_wp,100.0_wp,1.0_wp,2.0_wp,3.0_wp]  !! initial state
    real(wp),parameter :: dt = 1.0_wp * day2sec  !! time step (sec)
    real(wp),parameter :: mu = 398600.0_wp       !! Earth grav param.

    real(wp),dimension(6) :: rv2_gooding
    real(wp),dimension(6) :: rv2_shepperd

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' kepler_tests'
    write(*,*) '---------------'
    write(*,*) ''

    call propagate(mu,rv1,dt,rv2_gooding)
    call kepler_shepperd(mu,rv1,dt,rv2_shepperd)

    write(*,*) 'gooding-shepperd diff: ', rv2_shepperd-rv2_gooding
    write(*,*) ''

    end program kepler_tests
!*****************************************************************************************
