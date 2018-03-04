!*****************************************************************************************
!>
!  Approximate positions of the major planets.
!
!### See also
!  * [[analytical_ephemeris_module]] -- analytical ephemeris for Earth's Moon.
!
!### Reference
!  * E.M. Standish, Solar System Dynamics Group, JPL/Caltech,
!    "[Keplerian Elements for Approximate Positions of the Major Planets](http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf)"
!
!### History
!  * Original version copyright 2018 [https://github.com/CumuloEpsilon](CumuloEpsilon)
!    [BSD License](https://github.com/CumuloEpsilon/Standish-Ephemeris).
!  * Jacob Williams, extensive refactoring with some modifications,
!    and integration into the FAT ephemeris module.
!
!### Original license
!
!```
!   Copyright (c) 2018, CumuloEpsilon
!   All rights reserved.
!
!   Redistribution and use in source and binary forms, with or without
!   modification, are permitted provided that the following conditions are met:
!
!   * Redistributions of source code must retain the above copyright notice, this
!     list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above copyright notice,
!     this list of conditions and the following disclaimer in the documentation
!     and/or other materials provided with the distribution.
!
!   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
!   FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!   DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!   SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!   CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!   OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!```

    module standish_module

    use celestial_body_module
    use kind_module,            only: wp
    use conversion_module,      only: deg2rad,rad2deg,year2day,day2sec,au2m,m2km
    use numbers_module,         only: zero,one,two,twopi
    use ephemeris_module,       only: ephemeris_class
    use time_module,            only: et_to_jd
    use base_class_module,      only: base_class

    implicit none

    private

    type,extends(ephemeris_class),public :: standish_ephemeris
    contains
        procedure,public :: get_rv => standish_rv_func
    end type standish_ephemeris

    ! constants

    real(wp),parameter :: obliquity = 23.43928_wp         !! obliquity at J2000 [deg]
    real(wp),parameter :: s_sobl = sin(obliquity*deg2rad) !! sin of j2000 obliquity
    real(wp),parameter :: s_cobl = cos(obliquity*deg2rad) !! cos of j2000 obliquity
    real(wp),parameter :: s_dpc = 36525.0_wp              !! julian days per century
    real(wp),parameter :: mu_sun = 39.47692641_wp         !! au^3/yr^2
    real(wp),parameter :: epoch = 2451545.0_wp            !! Julian date of J2000 epoch

    type,extends(base_class) :: ephem

        !! an ephemeris defined using a date range
        !! and a set of elements from the reference.
        !!
        !! There are two that can be used.
        !!
        !!@note This should probably be merged into [[standish_ephemeris]]

        real(wp),dimension(2) :: jd_range = zero  !! valid julian date range
        real(wp),dimension (16, 9) :: o = zero    !! keplerian elements terms

    end type ephem

    !>
    !  The first ephemeris data table (this is standish's table 1).
    !  keplerian elements valid 1800 ad - 2050 ad
    type(ephem),parameter :: eph1 = ephem(1,&       ! id
        'Keplerian Elements Valid 1800AD-2050AD',&  ! name
        [2378497.0_wp, 2470172.0_wp],&              ! jd range
        reshape([ &   ! element table (in au and radians). perturbations are zero.
            0.38709927, 0.20563594, 0.12225995, 4.4025989, 1.3518935, &
            0.84353095, 3.70000009e-07, 1.90600003e-05, - 1.03803286e-04, &
            2608.7903, 2.80085020e-03, - 2.18760967e-03, 0.0000000, &
            0.0000000, 0.0000000, 0.0000000, 0.72333568, 6.77671982e-03, &
            5.92482723e-02, 3.1761343, 2.2968962, 1.3383157, 3.90000014e-06, &
            - 4.10700013e-05, - 1.37689030e-05, 1021.3286, 4.68322469e-05, - &
            4.84667765e-03, 0.0000000, 0.0000000, 0.0000000, 0.0000000, &
            1.0000026, 1.67112295e-02, - 2.67209913e-07, 1.7534375, &
            1.7966015, 0.0000000, 5.62000014e-06, - 4.39200012e-05, - &
            2.25962198e-04, 628.30756, 5.64218918e-03, 0.0000000, 0.0000000, &
            0.0000000, 0.0000000, 0.0000000, 1.5237104, 9.33941007e-02, &
            3.22832055e-02, - 7.94723779e-02, - 0.41789517, 0.86497712, &
            1.84700002e-05, 7.88199977e-05, - 1.41918135e-04, 334.06131, &
            7.75643345e-03, - 5.10636950e-03, 0.0000000, 0.0000000, &
            0.0000000, 0.0000000, 5.2028871, 4.83862385e-02, 2.27660220e-02, &
            0.60033119, 0.25706047, 1.7536005, - 1.16069998e-04, - &
            1.32529996e-04, - 3.20641411e-05, 52.966312, 3.70929041e-03, &
            3.57253314e-03, 0.0000000, 0.0000000, 0.0000000, 0.0000000, &
            9.5366764, 5.38617894e-02, 4.33887430e-02, 0.87186599, &
            1.6161553, 1.9837835, - 1.25059998e-03, - 5.09909994e-04, &
            3.37911442e-05, 21.336540, - 7.31244357e-03, - 5.03838016e-03, &
            0.0000000, 0.0000000, 0.0000000, 0.0000000, 19.189165, &
            4.72574383e-02, 1.34850740e-02, 5.4670362, 2.9837148, 1.2918390, &
            - 1.96175999e-03, - 4.39700016e-05, - 4.24008576e-05, 7.4784222, &
            7.12186471e-03, 7.40122399e-04, 0.0000000, 0.0000000, 0.0000000, &
            0.0000000, 30.069923, 8.59048031e-03, 3.08930874e-02, - &
            0.96202600, 0.78478318, 2.3000686, 2.62910005e-04, &
            5.10499995e-05, 6.17357864e-06, 3.8128369, - 5.62719675e-03, - &
            8.87786155e-05, 0.0000000, 0.0000000, 0.0000000, 0.0000000, &
            39.482117, 0.24882729, 0.29914966, 4.1700983, 3.9107401, &
            1.9251670, - 3.15960002e-04, 5.17000008e-05, 8.40899645e-07, &
            2.5343544, - 7.09117157e-04, - 2.06556579e-04, 0.0000000, &
            0.0000000, 0.0000000, 0.0000000 ], [16,9])&
        )

    !>
    !  The first ephemeris data table (this is standish's table 2).
    !  keplerian elements valid 3000 bc - 3000 ad
    type(ephem),parameter :: eph2 = ephem(2,&           ! id
        'Keplerian Elements Valid 3000BC-3000AD',&  ! name
        [625674.0_wp, 2816788.0_wp],&               ! jd range
        reshape([ &   ! element table (in au and radians). perturbations are not zero.
            0.38709843, 0.20563661, 0.12227069, 4.4026222, 1.3518922, &
            0.84368551, 0.0000000, 2.12300001e-05, - 1.03002007e-04, &
            2608.7903, 2.78205727e-03, - 2.13177688e-03, 0.0000000, &
            0.0000000, 0.0000000, 0.0000000, 0.72332102, 6.76399004e-03, &
            5.93023673e-02, 3.1761451, 2.2997777, 1.3381896, - &
            2.60000007e-07, - 5.10700011e-05, 7.59113527e-06, 1021.3286, &
            9.91285546e-04, - 4.76024114e-03, 0.0000000, 0.0000000, &
            0.0000000, 0.0000000, 1.0000002, 1.67316291e-02, - &
            9.48516663e-06, 1.7534785, 1.7964685, - 8.92317668e-02, - &
            2.99999989e-08, - 3.66099994e-05, - 2.33381579e-04, 628.30762, &
            5.54932002e-03, - 4.21040738e-03, 0.0000000, 0.0000000, &
            0.0000000, 0.0000000, 1.5237124, 9.33651105e-02, 3.23203318e-02, &
            - 7.97289312e-02, - 0.41743821, 0.86765921, 9.69999974e-07, &
            9.14900011e-05, - 1.26493964e-04, 334.06125, 7.89301097e-03, - &
            4.68663359e-03, 0.0000000, 0.0000000, 0.0000000, 0.0000000, &
            5.2024803, 4.85358983e-02, 2.26650927e-02, 0.59925520, &
            0.24914493, 1.7504400, - 2.86400009e-05, 1.80260002e-04, - &
            5.63216017e-05, 52.969063, 3.17635899e-03, 2.27322499e-03, - &
            2.17328397e-06, 1.05837814e-03, - 6.21955749e-03, 0.66935557, &
            9.5414991, 5.55082485e-02, 4.35327180e-02, 0.87398607, &
            1.6207365, 1.9833919, - 3.06500006e-05, - 3.20440013e-04, &
            7.88834659e-05, 21.329931, 9.45610274e-03, - 4.36594151e-03, &
            4.52022823e-06, - 2.34475732e-03, 1.52402408e-02, 0.66935557, &
            19.187979, 4.68574017e-02, 1.34910680e-02, 5.4838729, 3.0095420, &
            1.2908891, - 2.04550000e-04, - 1.54999998e-05, - 3.14429781e-05, &
            7.4786506, 1.61739404e-03, 1.00176642e-03, 1.01806800e-05, - &
            1.70574244e-02, 3.08735552e-03, 0.13387112, 30.069527, &
            8.95438995e-03, 3.08932904e-02, 5.3096914, 0.81474739, &
            2.3001058, 6.44699976e-05, 8.17999990e-06, 3.90953755e-06, &
            3.8129361, 1.76267436e-04, - 1.05819658e-04, - 7.21658762e-06, &
            1.19286822e-02, - 1.77369907e-03, 0.13387112, 39.486862, &
            0.24885239, 0.29916763, 4.1707320, 3.9112310, 1.9251275, &
            4.49750992e-03, 6.01600004e-05, 8.74410020e-08, 2.5338767, - &
            1.69092222e-04, - 1.41368364e-04, - 2.20386923e-04, 0.0000000, &
            0.0000000, 0.0000000 ], [16,9])&
        )

    type(ephem),dimension(2),parameter :: eph_set = [eph1,eph2]  !! the set of [[ephem]] options
                                                                 !! that are available.

    public :: standish_module_test  ! unit test

    contains
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 3/4/2018
!
!  Convert the NAIF SPICE ID code to the one used by the standish ephemeris.
!  Returns `0` if the body was not found.
!
!### See also
!  * [[spice_id_to_old_id]]

    pure function spice_id_to_standish_id(spice_id) result(old_id)

    implicit none

    integer,intent(in) :: spice_id !! the ID code used by SPICE
    integer            :: old_id   !! the ID code used by this module

    integer :: i !! counter

    !>
    !  The index of this array is the ID code. The value is the NAIF code.
    !  See: [NAIF Integer ID codes](http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/naif_ids.html)
    integer,parameter,dimension(9) :: new_ids = &
        [   199,&  ! mercury
            299,&  ! venus
            3,  &  ! earth-moon barycenter
            499,&  ! mars
            599,&  ! jupiter
            699,&  ! saturn
            799,&  ! uranus
            899,&  ! neptune
            999]   ! pluto

    !just a simple search of the list:
    ! [small enough that bisection search probably not worth it]
    do i=1,size(new_ids)
        if (new_ids(i)==spice_id) then
            old_id = i
            return
        end if
    end do

    !not found:
    old_id = 0

    end function spice_id_to_standish_id
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 3/4/2018
!
!  Return the state of the `targ` body relative to
!  the `obs` body, in the inertial frame [ICRF].
!
!  This is the function that can be used in the [[ephemeris_class]].

    subroutine standish_rv_func(me,et,targ,obs,rv,status_ok)

    implicit none

    class(standish_ephemeris),intent(inout) :: me
    real(wp),intent(in)                     :: et         !! ephemeris time [sec]
    type(celestial_body),intent(in)         :: targ       !! target body
    type(celestial_body),intent(in)         :: obs        !! observer body
    real(wp),dimension(6),intent(out)       :: rv         !! state of targ w.r.t. obs [km, km/s]
    logical,intent(out)                     :: status_ok  !! true if there were no problems

    real(wp) :: jd
    integer :: itbl_targ, itbl_obs
    integer :: itarg
    integer :: iobs
    real(wp), dimension(6) :: targ_rv_au
    real(wp), dimension(6) :: obs_rv_au
    real(wp),dimension(6) :: rv_ecliptic

    integer,parameter :: naif_id_sun = 10 !! NAIF ID for the Sun

    if (targ==obs) then
        rv = zero
        return
    end if

    jd = et_to_jd(et)

    if (targ%id==naif_id_sun) then
        itarg = -1 ! dummy
    else
        itarg = spice_id_to_standish_id(targ%id)
    end if

    if (obs%id==naif_id_sun) then
        iobs = -1 ! dummy
    else
        iobs = spice_id_to_standish_id(obs%id)
    end if

    if (itarg==0 .or. iobs==0) then

        ! if the bodies were not found in this ephemeris
        rv = zero
        status_ok = .false.

    else

        if (targ%id/=naif_id_sun) then
            ! targ w.r.t sun [j2000 ecliptic]
            call helio (itarg, jd, targ_rv_au, itbl_targ)
        else
            targ_rv_au = zero
            itbl_targ = 3 ! dummy value
        end if

        if (obs%id/=naif_id_sun) then
            ! obs w.r.t sun [j2000 ecliptic]
            call helio (iobs, jd, obs_rv_au, itbl_obs )
        else
            obs_rv_au = zero
            itbl_obs = 3 ! dummy value
        end if

        if (itbl_targ>0 .and. itbl_obs>0) then

            ! vector from obs to targ [j2000 ecliptic]
            rv_ecliptic = targ_rv_au - obs_rv_au

            ! convert to ICRF:
            call ec2eq (rv_ecliptic, rv)

            ! Convert to km, km/s:
            rv = rv * au2m * m2km
            rv(4:6) = rv(4:6) / (year2day*day2sec)

            status_ok = .true.

        else
            ! out of range of the ephemeris:
            rv = zero
            status_ok = .false.
        end if

    end if

    end subroutine standish_rv_func
!*****************************************************************************************

!*****************************************************************************************
!>
!  For planet np and julian date jd and using using table `itbl`,
!  return j2000 ecliptic position (au) and velocity (au/yr).
!  in cartesian coordinates (p(1)-p(6)).

    pure subroutine helio (np, jd, p, itbl)

    implicit none

    integer, intent (in) :: np      !! planet 1-9
    real(wp), intent (in) :: jd     !! julian date
    real(wp), dimension(6), intent(out) :: p !! position (au)/velocity (au/yr)
    integer, intent (out) :: itbl   !! table used or error if zero

    real(wp),dimension(8) :: z  !! elements [a, e, i, l, w, o, ma, ea]
    real(wp),dimension(8) :: po

    z = zero
    po = zero
    itbl = tbl (jd)
    if (itbl > 0) then
        call calcelements (np, jd, itbl, z)
        call el2op (z, po)
        call op2ec (z, po, p)
    end if

    end subroutine helio
!*****************************************************************************************

!*****************************************************************************************
!>
!  solve kepler's equation ma = ea + ec*sin(ea)

    pure real(wp) function kepler (ma, ec)

    implicit none

    real(wp), intent (in) :: ma !! eccentricity
    real(wp), intent (in) :: ec !! mean anomaly in rad

    real(wp) :: r, ea
    integer :: i

    ! acceptable accuracy for this calculation
    real(wp),parameter :: tol = 1.0e-08_wp  !! max error in eccentric anomaly `ea` in rad
    integer,parameter :: maxit = 12         !! max iterations (1-4 typical for `ec<0.3`)

    ea = ma + ec * sin (ma)! starting value
    do i = 1, maxit ! newton(-raphson) iterations
        r = (ma-ea+ec*sin(ea)) / (one-ec*cos(ea))
        ea = ea + r
        if (abs(r) <= tol) exit
    end do
    kepler = modulo (ea, twopi) ! eccentric anomaly adjusted 0-2pi

    end function kepler
!*****************************************************************************************

!*****************************************************************************************
!>
!  Determine which data set to use for highest accuracy for the given julian date.
!
!@note There was a typo in the original routine.
!@note Assumes that [[eph_set]] has only two elements.

    pure function tbl (jd) result(itbl)

    implicit none

    real(wp), intent (in) :: jd  !! julian date (eg 2451545.0)
    integer :: itbl !! Which data set to use:
                    !!
                    !! * itbl=1 jd in range of table 1
                    !!   (1800ad-2050ad) - highest accuracy
                    !! * itbl=2 jd outside range of table 1
                    !!   but in range of table 2 (3000bc-3000ad)
                    !! * itbl=0 3000bc<jd or jd>3000ad  julian
                    !!   date out of range for ephemeris.

    if ((jd > eph_set(1)%jd_range(1)) .and. (jd < eph_set(1)%jd_range(2))) then
        itbl = 1
    else
        if ((jd > eph_set(2)%jd_range(1)) .and. (jd < eph_set(2)%jd_range(2))) then
            itbl = 2
        else
            itbl = 0
        end if
    end if

    end function tbl
!*****************************************************************************************

!*****************************************************************************************
!>
!  Calculate current elements `z(jd)` for planet `j` from jpl data

    pure subroutine calcelements (np, jd, itbl, z)

    implicit none

    integer, intent (in) :: np        !! planet number (1-9)
    integer, intent (in) :: itbl      !! which table to use (1-2)
    real(wp), intent (in) :: jd       !! julian date
    real(wp), dimension(8), intent(out) :: z  !! elements for `jd`
                                              !!
                                              !!  * a = semi major axis (au)
                                              !!  * e = eccentricity (rad)
                                              !!  * i = inclination (rad)
                                              !!  * l = mean longitude (rad)
                                              !!  * w = longitude of perihelion (rad)
                                              !!  * o = longitude of ascending mode (rad)
                                              !!  * ma = mean anomaly (rad)
                                              !!  * ea = eccentric anomaly (rad)

    select case (itbl)
    case(1)
        z = eph(eph_set(1)%o)
    case(2)
        z = eph(eph_set(2)%o)
    case default
        error stop 'invalid value of itbl in calcelements'
    end select

    contains

        pure function eph(o) result(z)

        implicit none

        real(wp),dimension(:,:),intent(in) :: o !! data table to use
        real(wp), dimension(8) :: z  !! result

        integer :: i   !! counter
        real(wp) :: t  !! centuries since epoch
        real(wp) :: tz !! perturbation term

        t = (jd-epoch) / s_dpc

        do i = 1, 6      ! a,e,i,l,w,o
            z (i) = o(i, np) + o(i+6, np) * t
            ! if (i>2) z(i) = modulo(z(i), twopi)  !optional scaling
        end do

        if (itbl==2) then
            ! perturbation term nonzero for planets 5-9 if table 2 used
            tz =  o(13, np) * t ** 2 + &
                  o(14, np) * cos(o(16, np)*t) + &
                  o(15, np) * sin(o(16, np)*t)
        else
            tz = zero
        end if

        z (7) = modulo ((z(4)-z(5)+tz), twopi)  ! mean anomaly
        z (8) = kepler (z(7), z(2))             ! eccentric anomaly

        end function eph

    end subroutine calcelements
!*****************************************************************************************

!*****************************************************************************************
!>
!  heliocentric coordinates for orbital plane from elements

    pure subroutine el2op (z, po)

    implicit none

    real(wp), dimension(8), intent (in) :: z   !! elements [a,e,i,l,w,o,ma,ea]
    real(wp), dimension(6), intent (out) :: po !! coordinates and velocities

    real(wp) :: v, s1, c1, s2

    ! heliocentric orbital plane
    po = zero
    s1 = sin (z(8))
    c1 = cos (z(8))
    s2 = sqrt (one-z(2)*z(2))
    v = twopi / (sqrt(z(1))*(one-z(2)*c1))  ! velocity au/yr

    po (1) = z (1) * (c1-z(2)) ! xp (plane of orbit)
    po (2) = z (1) * s1 * s2   ! yp
    po (4) = - v * s1          ! vxp
    po (5) = v * c1 * s2       ! vyp

    end subroutine el2op
!*****************************************************************************************

!*****************************************************************************************
!>
!  heliocentric coordinates j2000 ecliptic plane from orbital plane

    pure subroutine op2ec (z, po, pe)

    implicit none

    real(wp),dimension(8),intent(in)  :: z  !! elements a,e,i,l,w,o,ma,ea
    real(wp),dimension(6),intent(in)  :: po !! orbital plane coordinates
    real(wp),dimension(6),intent(out) :: pe !! j2000 ecliptic plane coordinates

    real(wp) :: s1, s2, s3, c1, c2, c3

    ! heliocentric au, au/yr
    s1 = sin (z(5)-z(6))
    s2 = sin (z(3))
    s3 = sin (z(6))
    c1 = cos (z(5)-z(6))
    c2 = cos (z(3))
    c3 = cos (z(6))

    pe (1) = (c1*c3-s1*s3*c2) * po (1) - (s1*c3+c1*s3*c2) * po (2) ! xec
    pe (2) = (c1*s3+s1*c3*c2) * po (1) - (s1*s3-c1*c3*c2) * po (2) ! yec
    pe (3) = s1 * s2 * po (1) + c1 * s2 * po (2)                   ! zec
    pe (4) = (c1*c3-s1*s3*c2) * po (4) - (s1*c3+c1*s3*c2) * po (5) ! vxec
    pe (5) = (c1*s3+s1*c3*c2) * po (4) - (s1*s3-c1*c3*c2) * po (5) ! vyec
    pe (6) = s1 * s2 * po (4) + c1 * s2 * po (5)                   ! vzec

    end subroutine op2ec
!*****************************************************************************************

!*****************************************************************************************
!>
!  converts cartesian heliocentric j2000 ecliptic to equatorial

    pure subroutine ec2eq (pe, pq)

    implicit none

    real(wp), dimension(6), intent (in) :: pe  !! ecliptic
    real(wp), dimension(6), intent (out) :: pq !! equatorial

    pq (1) = pe (1)                             ! xeq same as xec
    pq (2) = s_cobl * pe (2) - s_sobl * pe (3)  ! yeq
    pq (3) = s_sobl * pe (2) + s_cobl * pe (3)  ! zeq
    pq (4) = pe (4)                             ! vxeq same as vxec
    pq (5) = s_cobl * pe (5) - s_sobl * pe (6)  ! vyeq
    pq (6) = s_sobl * pe (5) + s_cobl * pe (6)  ! vzeq

    end subroutine ec2eq
!*****************************************************************************************

! the following two aren't used:

!*****************************************************************************************
!>
!  converts cartesian heliocentric equatorial to ecliptic

    pure subroutine eq2ec (pq, pe)

    implicit none

    real(wp), dimension(6), intent (out) :: pe  !! ecliptic
    real(wp), dimension(6), intent (in) :: pq   !! equatorial

    pe (1) = pq (1)                              ! xec same as xeq
    pe (2) = s_cobl * pq (2) + s_sobl * pq (3)   ! yec
    pe (3) = - s_sobl * pq (2) + s_cobl * pq (3) ! zec
    pe (4) = pq (4)                              ! vxec same as vxeq
    pe (5) = s_cobl * pq (5) + s_sobl * pq (6)   ! vyec
    pe (6) = - s_sobl * pq (5) + s_cobl * pq (6) ! vzec

    end subroutine eq2ec
!*****************************************************************************************

!*****************************************************************************************
!>
!  cartesian to spherical coordinates (angles in radians)

    pure subroutine sphere (x, y, z, rho, theta, phi)

    implicit none

    real(wp), intent (in) :: x       !! x = r cos(phi) cos (theta)
    real(wp), intent (in) :: y       !! y = r cos(phi) sin(theta)
    real(wp), intent (in) :: z       !! z = r sin(phi)
    real(wp), intent (out) :: rho    !! distance
    real(wp), intent (out) :: theta  !! longitude
    real(wp), intent (out) :: phi    !! latitude

    real(wp) :: r

    theta = zero
    phi = zero
    rho = sqrt (x*x+y*y+z*z)
    r = sqrt (x*x+y*y)
    if (r /= zero) then
        theta = modulo (atan2(y, x), twopi)
        phi = atan2 (z, r)
    end if

    end subroutine sphere
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 3/4/2018
!
!  Test routine for the standish_module routines.

    subroutine standish_module_test()

    implicit none

    type(standish_ephemeris) :: eph
    real(wp)              :: et         !! ephemeris time (sec)
    real(wp)              :: jd         !! julian date
    real(wp),dimension(6) :: rv         !! state of targ w.r.t. obs
    logical               :: status_ok  !! true if there were no problems
    integer :: itbl

    !>
    !  State of Earth w.r.t. Sun from SPICE
    !  See: http://wgc.jpl.nasa.gov:8080/webgeocalc/#NewCalculation
    real(wp),dimension(6),parameter :: rv_from_spice = [-26502576.96907434_wp,&
                                                        132754176.58943012_wp,&
                                                        57555793.70952077_wp,&
                                                        -29.78644078_wp,&
                                                        -5.02614568_wp,&
                                                        -2.17905509_wp]

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' standish_module_test'
    write(*,*) '---------------'
    write(*,*) ''

    et = 0.0_wp ! J2000

    ! get earth w.r.t. sun in J2000 frame:
    call eph%get_rv(et,body_earth_moon_barycenter,body_sun,rv,status_ok)

    if (status_ok) then
        write(*,*) ''
        write(*,*) 'State of Earth wrt Sun @ J2000'
        write(*,'(A,1X,*(E26.16,1X))') 'standish:', rv
        write(*,'(A,1X,*(E26.16,1X))') 'SPICE:   ', rv_from_spice
        write(*,'(A,1X,*(E26.16,1X))') 'diff:    ', rv - rv_from_spice
        write(*,*) ''
    else
        error stop 'error calling standish_ephemeris'
    end if

    ! low-level routine tests:
    write(*,*) 'helio:'
    jd = et_to_jd(et)
    call helio (3, jd, rv, itbl)
    write(*,*) ''
    write(*,*) 'State of Earth wrt Sun [ecliptic]'
    write(*,'(*(E26.16,1X))') rv
    write(*,*) ''

    end subroutine standish_module_test
!*****************************************************************************************

!*****************************************************************************************
    end module standish_module
!*****************************************************************************************
