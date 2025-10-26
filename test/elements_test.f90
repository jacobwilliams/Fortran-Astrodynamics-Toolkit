!*****************************************************************************************
!>
!  Orbital elements conversion tests.
!
!  This tests three types of orbital elements:
!
!  * classical
!  * modified equinoctial
!  * gooding

    program elements_test

    use fortran_astrodynamics_toolkit, wp => fat_wp

    implicit none

    real(wp),dimension(3) :: r,v,r2,v2
    real(wp),dimension(6) :: rv,rv2,rv3,evec,gvec
    real(wp) :: mu
    real(wp) :: a,p,ecc,inc,raan,aop,tru
    real(wp) :: a2,p2,ecc2,inc2,raan2,aop2,tru2
    integer  :: i  !! counter

    mu = 3.9860043543609593E+05_wp ! for earth

    do i = 1, 8

        write(*,*) ''
        write(*,*) '---------------------------'

        select case (i)
        case (1)
            write(*,*) 'general elliptical:'
            a    = 8000.0_wp ! km
            ecc  = 0.1_wp
            inc  = 45.0_wp * deg2rad
            raan = 45.0_wp * deg2rad
            aop  = 45.0_wp * deg2rad
            tru  = 45.0_wp * deg2rad
            p = a * (one - ecc**2)
        case (2)
            write(*,*) 'circular:'
            a    = 8000.0_wp ! km
            ecc  = 0.00001_wp
            inc  = 45.0_wp * deg2rad
            raan = 45.0_wp * deg2rad
            aop  = 45.0_wp * deg2rad
            tru  = 45.0_wp * deg2rad
            p = a * (one - ecc**2)
        case (3)
            write(*,*) 'equatorial elliptical:'
            a    = 8000.0_wp ! km
            ecc  = 0.1_wp
            inc  = 0.00001_wp * deg2rad
            raan = 45.0_wp * deg2rad
            aop  = 45.0_wp * deg2rad
            tru  = 45.0_wp * deg2rad
            p = a * (one - ecc**2)
        case (4)
            write(*,*) 'circular and equatorial:'
            a    = 8000.0_wp ! km
            ecc  = 0.00001_wp
            inc  = 0.00001_wp * deg2rad
            raan = 45.0_wp * deg2rad
            aop  = 45.0_wp * deg2rad
            tru  = 45.0_wp * deg2rad
            p = a * (one - ecc**2)
        case (5)
            write(*,*) 'general hyperbolic:'
            a    = -80000.0_wp ! km
            ecc  = 1.1_wp
            inc  = 45.0_wp * deg2rad
            raan = 45.0_wp * deg2rad
            aop  = 45.0_wp * deg2rad
            tru  = 45.0_wp * deg2rad
            p = a * (one - ecc**2)
        case (6)
            write(*,*) 'equatorial hyperbolic:'
            a    = -80000.0_wp ! km
            ecc  = 1.1_wp
            inc  = 0.00001_wp * deg2rad
            raan = 45.0_wp * deg2rad
            aop  = 45.0_wp * deg2rad
            tru  = 45.0_wp * deg2rad
            p = a * (one - ecc**2)
        case (7)
            write(*,*) 'general parabolic:'
            p    = 80000.0_wp
            ecc  = 1.0_wp
            inc  = 45.0_wp * deg2rad
            raan = 45.0_wp * deg2rad
            aop  = 45.0_wp * deg2rad
            tru  = 45.0_wp * deg2rad
        case (8)
            write(*,*) 'equatorial parabolic:'
            p    = 80000.0_wp
            ecc  = 1.0_wp
            inc  = 0.00001_wp * deg2rad
            raan = 45.0_wp * deg2rad
            aop  = 45.0_wp * deg2rad
            tru  = 45.0_wp * deg2rad
        case default
            error stop 'invalid case number'
        end select

        write(*,*) '---------------------------'

        call orbital_elements_to_rv(mu, p, ecc, inc, raan, aop, tru, r, v)
        rv = [r,v]

        ! classical orbital elements:
        call rv_to_orbital_elements(mu, r, v, p2, ecc2, inc2, raan2, aop2, tru2)
        call orbital_elements_to_rv(mu, p2, ecc2, inc2, raan2, aop2, tru2, r2, v2)

        ! gooding elements:
        call pv3els(mu, rv, gvec)
        call els3pv(mu, gvec, rv2)

        ! modified equinoctial elements:
        call cartesian_to_equinoctial(mu, rv, evec)
        call equinoctial_to_cartesian(mu, evec, rv3)

        write(*,*) ''
        write(*,'(A22,*(E26.16,1X))') 'cartesian:', rv
        write(*,'(A22,*(E26.16,1X))') 'classical elements:', p, ecc, inc, raan, aop, tru
        write(*,'(A22,*(E26.16,1X))') 'gooding elements:', gvec
        write(*,'(A22,*(E26.16,1X))') 'equinoctial elements:', evec

        write(*,*) ''
        write(*,'(A22,*(E26.16,1X))') 'elements r diff:', r2 - r
        write(*,'(A22,*(E26.16,1X))') 'elements v diff:', v2 - v
        write(*,*) ''
        write(*,'(A22,*(E26.16,1X))') 'gooding r diff:', rv2(1:3) - rv(1:3)
        write(*,'(A22,*(E26.16,1X))') 'gooding v diff:', rv2(4:6) - rv(4:6)
        write(*,*) ''
        write(*,'(A22,*(E26.16,1X))') 'equinoctial r diff:', rv3(1:3) - rv(1:3)
        write(*,'(A22,*(E26.16,1X))') 'equinoctial v diff:', rv3(4:6) - rv(4:6)

    end do

    write(*,*) ''

!*****************************************************************************************
    end program elements_test
!*****************************************************************************************
