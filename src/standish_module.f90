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
    use conversion_module,      only: deg2rad,rad2deg,year2day,&
                                      day2sec,au2m,m2km,day2century
    use numbers_module,         only: zero,one,two,twopi
    use ephemeris_module,       only: ephemeris_class
    use time_module,            only: et_to_jd
    use base_class_module,      only: base_class

    implicit none

    private

    type,extends(ephemeris_class),public :: standish_ephemeris
        !! Standish ephemeris class for computing the
        !! approximate positions of the major planets.
    contains
        procedure,public :: get_rv => standish_rv_func
    end type standish_ephemeris

    ! constants

    real(wp),parameter :: obliquity = 23.43928_wp         !! obliquity at J2000 [deg]
    real(wp),parameter :: s_sobl = sin(obliquity*deg2rad) !! sin of j2000 obliquity
    real(wp),parameter :: s_cobl = cos(obliquity*deg2rad) !! cos of j2000 obliquity
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
            0.38709927_wp, 0.20563594_wp, 0.12225995_wp, 4.4025989_wp, 1.3518935_wp, &
            0.84353095_wp, 3.70000009e-07_wp, 1.90600003e-05_wp, - 1.03803286e-04_wp, &
            2608.7903_wp, 2.80085020e-03_wp, - 2.18760967e-03_wp, 0.0_wp, &
            0.0_wp, 0.0_wp, 0.0_wp, 0.72333568_wp, 6.77671982e-03_wp, &
            5.92482723e-02_wp, 3.1761343_wp, 2.2968962_wp, 1.3383157_wp, 3.90000014e-06_wp, &
            - 4.10700013e-05_wp, - 1.37689030e-05_wp, 1021.3286_wp, 4.68322469e-05_wp, - &
            4.84667765e-03_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
            1.0000026_wp, 1.67112295e-02_wp, - 2.67209913e-07_wp, 1.7534375_wp, &
            1.7966015_wp, 0.0_wp, 5.62000014e-06_wp, - 4.39200012e-05_wp, - &
            2.25962198e-04_wp, 628.30756_wp, 5.64218918e-03_wp, 0.0_wp, 0.0_wp, &
            0.0_wp, 0.0_wp, 0.0_wp, 1.5237104_wp, 9.33941007e-02_wp, &
            3.22832055e-02_wp, - 7.94723779e-02_wp, - 0.41789517_wp, 0.86497712_wp, &
            1.84700002e-05_wp, 7.88199977e-05_wp, - 1.41918135e-04_wp, 334.06131_wp, &
            7.75643345e-03_wp, - 5.10636950e-03_wp, 0.0_wp, 0.0_wp, &
            0.0_wp, 0.0_wp, 5.2028871_wp, 4.83862385e-02_wp, 2.27660220e-02_wp, &
            0.60033119_wp, 0.25706047_wp, 1.7536005_wp, - 1.16069998e-04_wp, - &
            1.32529996e-04_wp, - 3.20641411e-05_wp, 52.966312_wp, 3.70929041e-03_wp, &
            3.57253314e-03_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
            9.5366764_wp, 5.38617894e-02_wp, 4.33887430e-02_wp, 0.87186599_wp, &
            1.6161553_wp, 1.9837835_wp, - 1.25059998e-03_wp, - 5.09909994e-04_wp, &
            3.37911442e-05_wp, 21.336540_wp, - 7.31244357e-03_wp, - 5.03838016e-03_wp, &
            0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 19.189165_wp, &
            4.72574383e-02_wp, 1.34850740e-02_wp, 5.4670362_wp, 2.9837148_wp, 1.2918390_wp, &
            - 1.96175999e-03_wp, - 4.39700016e-05_wp, - 4.24008576e-05_wp, 7.4784222_wp, &
            7.12186471e-03_wp, 7.40122399e-04_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
            0.0_wp, 30.069923_wp, 8.59048031e-03_wp, 3.08930874e-02_wp, - &
            0.96202600_wp, 0.78478318_wp, 2.3000686_wp, 2.62910005e-04_wp, &
            5.10499995e-05_wp, 6.17357864e-06_wp, 3.8128369_wp, - 5.62719675e-03_wp, - &
            8.87786155e-05_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
            39.482117_wp, 0.24882729_wp, 0.29914966_wp, 4.1700983_wp, 3.9107401_wp, &
            1.9251670_wp, - 3.15960002e-04_wp, 5.17000008e-05_wp, 8.40899645e-07_wp, &
            2.5343544_wp, - 7.09117157e-04_wp, - 2.06556579e-04_wp, 0.0_wp, &
            0.0_wp, 0.0_wp, 0.0_wp ], [16,9])&
        )

    !>
    !  The first ephemeris data table (this is standish's table 2).
    !  keplerian elements valid 3000 bc - 3000 ad
    type(ephem),parameter :: eph2 = ephem(2,&           ! id
        'Keplerian Elements Valid 3000BC-3000AD',&  ! name
        [625674.0_wp, 2816788.0_wp],&               ! jd range
        reshape([ &   ! element table (in au and radians). perturbations are not zero.
            0.38709843_wp, 0.20563661_wp, 0.12227069_wp, 4.4026222_wp, 1.3518922_wp, &
            0.84368551_wp, 0.0_wp, 2.12300001e-05_wp, - 1.03002007e-04_wp, &
            2608.7903_wp, 2.78205727e-03_wp, - 2.13177688e-03_wp, 0.0_wp, &
            0.0_wp, 0.0_wp, 0.0_wp, 0.72332102_wp, 6.76399004e-03_wp, &
            5.93023673e-02_wp, 3.1761451_wp, 2.2997777_wp, 1.3381896_wp, - &
            2.60000007e-07_wp, - 5.10700011e-05_wp, 7.59113527e-06_wp, 1021.3286_wp, &
            9.91285546e-04_wp, - 4.76024114e-03_wp, 0.0_wp, 0.0_wp, &
            0.0_wp, 0.0_wp, 1.0000002_wp, 1.67316291e-02_wp, - &
            9.48516663e-06_wp, 1.7534785_wp, 1.7964685_wp, - 8.92317668e-02_wp, - &
            2.99999989e-08_wp, - 3.66099994e-05_wp, - 2.33381579e-04_wp, 628.30762_wp, &
            5.54932002e-03_wp, - 4.21040738e-03_wp, 0.0_wp, 0.0_wp, &
            0.0_wp, 0.0_wp, 1.5237124_wp, 9.33651105e-02_wp, 3.23203318e-02_wp, &
            - 7.97289312e-02_wp, - 0.41743821_wp, 0.86765921_wp, 9.69999974e-07_wp, &
            9.14900011e-05_wp, - 1.26493964e-04_wp, 334.06125_wp, 7.89301097e-03_wp, - &
            4.68663359e-03_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
            5.2024803_wp, 4.85358983e-02_wp, 2.26650927e-02_wp, 0.59925520_wp, &
            0.24914493_wp, 1.7504400_wp, - 2.86400009e-05_wp, 1.80260002e-04_wp, - &
            5.63216017e-05_wp, 52.969063_wp, 3.17635899e-03_wp, 2.27322499e-03_wp, - &
            2.17328397e-06_wp, 1.05837814e-03_wp, - 6.21955749e-03_wp, 0.66935557_wp, &
            9.5414991_wp, 5.55082485e-02_wp, 4.35327180e-02_wp, 0.87398607_wp, &
            1.6207365_wp, 1.9833919_wp, - 3.06500006e-05_wp, - 3.20440013e-04_wp, &
            7.88834659e-05_wp, 21.329931_wp, 9.45610274e-03_wp, - 4.36594151e-03_wp, &
            4.52022823e-06_wp, - 2.34475732e-03_wp, 1.52402408e-02_wp, 0.66935557_wp, &
            19.187979_wp, 4.68574017e-02_wp, 1.34910680e-02_wp, 5.4838729_wp, 3.0095420_wp, &
            1.2908891_wp, - 2.04550000e-04_wp, - 1.54999998e-05_wp, - 3.14429781e-05_wp, &
            7.4786506_wp, 1.61739404e-03_wp, 1.00176642e-03_wp, 1.01806800e-05_wp, - &
            1.70574244e-02_wp, 3.08735552e-03_wp, 0.13387112_wp, 30.069527_wp, &
            8.95438995e-03_wp, 3.08932904e-02_wp, 5.3096914_wp, 0.81474739_wp, &
            2.3001058_wp, 6.44699976e-05_wp, 8.17999990e-06_wp, 3.90953755e-06_wp, &
            3.8129361_wp, 1.76267436e-04_wp, - 1.05819658e-04_wp, - 7.21658762e-06_wp, &
            1.19286822e-02_wp, - 1.77369907e-03_wp, 0.13387112_wp, 39.486862_wp, &
            0.24885239_wp, 0.29916763_wp, 4.1707320_wp, 3.9112310_wp, 1.9251275_wp, &
            4.49750992e-03_wp, 6.01600004e-05_wp, 8.74410020e-08_wp, 2.5338767_wp, - &
            1.69092222e-04_wp, - 1.41368364e-04_wp, - 2.20386923e-04_wp, 0.0_wp, &
            0.0_wp, 0.0_wp ], [16,9])&
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
!
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

    integer :: i   !! counter
    real(wp) :: t  !! centuries since epoch
    real(wp) :: tz !! perturbation term

    t = (jd-epoch) * day2century

    do i = 1, 6      ! a,e,i,l,w,o
        z (i) = eph_set(itbl)%o(i, np) + eph_set(itbl)%o(i+6, np) * t
        ! if (i>2) z(i) = modulo(z(i), twopi)  !optional scaling
    end do

    if (itbl==2) then
        ! perturbation term nonzero for planets 5-9 if table 2 used
        tz =  eph_set(itbl)%o(13, np) * t ** 2 + &
              eph_set(itbl)%o(14, np) * cos(eph_set(itbl)%o(16, np)*t) + &
              eph_set(itbl)%o(15, np) * sin(eph_set(itbl)%o(16, np)*t)
    else
        tz = zero
    end if

    z (7) = modulo ((z(4)-z(5)+tz), twopi)  ! mean anomaly
    z (8) = kepler (z(7), z(2))             ! eccentric anomaly

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
