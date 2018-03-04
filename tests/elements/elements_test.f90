!*****************************************************************************************
!>
!  Orbital elements conversion tests.

    program elements_test

    use fortran_astrodynamics_toolkit

    implicit none

    real(wp),dimension(3) :: r,v,r2,v2
    real(wp) :: mu
    real(wp) :: a,p,ecc,inc,raan,aop,tru
    real(wp) :: a2,p2,ecc2,inc2,raan2,aop2,tru2

    mu = 3.9860043543609593E+05_wp ! for earth

    ! general:
    a    = 8000.0_wp ! km
    ecc  = 0.1_wp
    inc  = 45.0_wp * deg2rad
    raan = 45.0_wp * deg2rad
    aop  = 45.0_wp * deg2rad
    tru  = 45.0_wp * deg2rad
    p = a * (one - ecc**2)

    call orbital_elements_to_rv(mu, p, ecc, inc, raan, aop, tru, r, v)
    call rv_to_orbital_elements(mu, r, v, p2, ecc2, inc2, raan2, aop2, tru2)
    call orbital_elements_to_rv(mu, p2, ecc2, inc2, raan2, aop2, tru2, r2, v2)
    write(*,*) ''
    write(*,*) 'r diff: ', r2 - r
    write(*,*) 'v diff: ', v2 - v

    ! circular:
    a    = 8000.0_wp ! km
    ecc  = 0.00001_wp
    inc  = 45.0_wp * deg2rad
    raan = 45.0_wp * deg2rad
    aop  = 45.0_wp * deg2rad
    tru  = 45.0_wp * deg2rad
    p = a * (one - ecc**2)

    call orbital_elements_to_rv(mu, p, ecc, inc, raan, aop, tru, r, v)
    call rv_to_orbital_elements(mu, r, v, p2, ecc2, inc2, raan2, aop2, tru2)
    call orbital_elements_to_rv(mu, p2, ecc2, inc2, raan2, aop2, tru2, r2, v2)
    write(*,*) ''
    write(*,*) 'r diff: ', r2 - r
    write(*,*) 'v diff: ', v2 - v

    ! equatorial:
    a    = 8000.0_wp ! km
    ecc  = 0.1_wp
    inc  = 0.00001_wp * deg2rad
    raan = 45.0_wp * deg2rad
    aop  = 45.0_wp * deg2rad
    tru  = 45.0_wp * deg2rad
    p = a * (one - ecc**2)

    call orbital_elements_to_rv(mu, p, ecc, inc, raan, aop, tru, r, v)
    call rv_to_orbital_elements(mu, r, v, p2, ecc2, inc2, raan2, aop2, tru2)
    call orbital_elements_to_rv(mu, p2, ecc2, inc2, raan2, aop2, tru2, r2, v2)
    write(*,*) ''
    write(*,*) 'r diff: ', r2 - r
    write(*,*) 'v diff: ', v2 - v

    ! circular and equatorial:
    a    = 8000.0_wp ! km
    ecc  = 0.00001_wp
    inc  = 0.00001_wp * deg2rad
    raan = 45.0_wp * deg2rad
    aop  = 45.0_wp * deg2rad
    tru  = 45.0_wp * deg2rad
    p = a * (one - ecc**2)

    call orbital_elements_to_rv(mu, p, ecc, inc, raan, aop, tru, r, v)
    call rv_to_orbital_elements(mu, r, v, p2, ecc2, inc2, raan2, aop2, tru2)
    call orbital_elements_to_rv(mu, p2, ecc2, inc2, raan2, aop2, tru2, r2, v2)
    write(*,*) ''
    write(*,*) 'r diff: ', r2 - r
    write(*,*) 'v diff: ', v2 - v

!*****************************************************************************************
    end program elements_test
!*****************************************************************************************
