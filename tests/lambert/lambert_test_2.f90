!*****************************************************************************************
!>
!  Just an example using the lambert module.

    program lambert_test_2

    use fortran_astrodynamics_toolkit, wp => fat_wp

    implicit none

    real(wp),dimension(3)               :: r1         !! first cartesian position [km]
    real(wp),dimension(3)               :: r2         !! second cartesian position [km]
    real(wp)                            :: tof        !! time of flight [sec]
    real(wp)                            :: mu         !! gravity parameter [km^3/s^2]
    logical                             :: long_way   !! when true, do "long way" (>pi) transfers
    integer                             :: multi_revs !! maximum number of multi-rev solutions to compute
    real(wp),dimension(:,:),allocatable :: v1         !! vector containing 3d arrays with the cartesian components of the velocities at r1
    real(wp),dimension(:,:),allocatable :: v2         !! vector containing 3d arrays with the cartesian components of the velocities at r2
    logical                             :: status_ok  !! true if everything is OK

    integer                 :: i,j,n,nsol,iunit_short,iunit_long
    real(wp),dimension(6)   :: rv,els ![al, q, ei, bom, om, tau]
    real(wp),dimension(20)  :: alpha
    real(wp)                :: err

    open(newunit=iunit_short,file='test2_short.csv',status='REPLACE')
    open(newunit=iunit_long, file='test2_long.csv', status='REPLACE')

    r1 = [20.0d3, 20.0d3, zero]  ! initial position
    r2 = [-20.0d3, 10.0d3, zero] ! final position
    mu = 398600.0_wp             ! earth
    multi_revs = 5               ! try up to 5 revs
    n = 2*day2sec                ! 2 days

    do i=5000,n,30  !transfer time loop [sec]

        tof = dble(i)    !sec

        !long way solutions:
        nsol = 0
        long_way=.true.
        call solve_lambert_gooding(r1,r2,tof,mu,long_way,multi_revs,v1,v2,status_ok)
        !call solve_lambert_izzo(r1,r2,tof,mu,long_way,multi_revs,v1,v2,status_ok)
        if (.not. status_ok) write(*,*) 'error!!'
        do j=1,size(v1,2)
            nsol = nsol + 1
            call pv3els (mu, [r1,v1(:,j)], els)  !convert transfer orbit state to orb elements
            alpha(nsol) = els(1)

            !verify solution:
            call propagate(mu, [r1,v1(:,j)], tof, rv)
            err = norm2(rv(1:3)-r2)
            if (err>1.0e-8_wp) write(*,*) 'Error:',err

        end do
        write(iunit_short,'(I3,A,1X,F10.0,A,1X,*(F20.6,1X,1H,))') nsol, ',',tof, ',',alpha(1:nsol)

        !short way solutions:
        nsol = 0
        long_way=.false.
        call solve_lambert_gooding(r1,r2,tof,mu,long_way,multi_revs,v1,v2,status_ok)
        !call solve_lambert_izzo(r1,r2,tof,mu,long_way,multi_revs,v1,v2,status_ok)
        if (.not. status_ok) write(*,*) 'error!!'
        do j=1,size(v1,2)
            nsol = nsol + 1
            call pv3els (mu, [r1,v1(:,j)], els)  !convert transfer orbit state to orb elements
            alpha(nsol) = els(1)

            !verify solution:
            call propagate(mu, [r1,v1(:,j)], tof, rv)
            err = norm2(rv(1:3)-r2)
            if (err>1.0e-8_wp) write(*,*) 'Error:',err

        end do
        write(iunit_long,'(I3,A,1X,F10.0,A,1X,*(F20.6,1X,1H,))') nsol, ',',tof, ',',alpha(1:nsol)

    end do

    close(iunit_short)
    close(iunit_long)

    end program lambert_test_2
!*****************************************************************************************
