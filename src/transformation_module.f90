!*****************************************************************************************
!> author: Jacob Williams
!
!  Coordinate transformations.

    module transformation_module

    use numbers_module
    use kind_module
    use vector_module
    use ephemeris_module
    use iau_orientation_module
    use time_module
    use celestial_body_module
    use iso_fortran_env, only: error_unit

    implicit none

    private

    !frame center options for two_body_rotating_frame:
    integer,parameter,public :: center_at_primary_body   = 1
    integer,parameter,public :: center_at_secondary_body = 2
    integer,parameter,public :: center_at_barycenter     = 3

    !abstract frame classes:

    type,abstract,public :: reference_frame
        !! a reference frame defines an orientation at
        !! a specified frame center at a specified epoch.
        !! Usually, the center is the `primary_body` of the
        !! frame, but can be otherwise for a
        !! [[two_body_rotating_frame]].
        private
        type(celestial_body) :: primary_body = body_earth  !! the primary body of the frame
        real(wp) :: et = zero  !! epoch at which the frame is defined [sec]
    contains
        private
        procedure,pass(from),public :: transform        !! coordinate transformation routine
        procedure(c_cdot_func),deferred :: get_c_cdot   !! to get the rotating matrix for the frame orientation
    end type reference_frame

    type,abstract,extends(reference_frame),public :: inertial_frame_class
        !! a non-rotating frame (a frame where the orientation axes are time invariant)
    end type inertial_frame_class

    type,abstract,extends(reference_frame),public :: rotating_frame_class
        !! a rotating frame (a frame where the orientation axes vary with time)
    end type rotating_frame_class

    type,abstract,extends(rotating_frame_class),public :: iau_rotating_frame_class
        !! frame defined by the orientation of a celestial body using the IAU models.
    end type iau_rotating_frame_class

    !concrete rotating frames:

    type,extends(rotating_frame_class),public :: two_body_rotating_frame
        !! The two-body rotating frame is constructed from the states
        !! of two celestial bodies.
        type(celestial_body) :: secondary_body = body_moon  !! the secondary body used to construct the frame
        integer :: center = center_at_barycenter  !! the frame center (can be primary_body,secondary_body, or barycenter)
        real(wp),dimension(6) :: rv12 = zero  !! [r,v] of secondary body w.r.t. primary body
    contains
        procedure :: from_primary_to_center
        procedure :: get_c_cdot => get_c_cdot_two_body_rotating
    end type two_body_rotating_frame
    interface two_body_rotating_frame
        module procedure :: two_body_rotating_frame_constructor
    end interface two_body_rotating_frame

    type,extends(two_body_rotating_frame),public :: two_body_rotating_pulsating_frame
        !! This frame is an extension of the two-body rotating frame, where
        !! a scale factor is used to scale the position and velocity of the state
        !! based on the distance between the primary and secondary bodies.
        real(wp) :: scale = zero !! scale factor
    contains
        procedure :: get_c_cdot => get_c_cdot_two_body_rotating_pulsating
    end type two_body_rotating_pulsating_frame
    interface two_body_rotating_pulsating_frame
        module procedure :: two_body_rotating_pulsating_frame_constructor
    end interface two_body_rotating_pulsating_frame

    !... note: could also have a two_body_rotating_pulsating_crtbp_frame
    !          which scales r,v,t based on the crtbp assumptions...

    !IAU rotating frames [not finished...]:
    type,extends(iau_rotating_frame_class),public :: iau_earth_rotating_frame
        !! IAU Earth frame
    contains
        procedure :: get_c_cdot => get_c_cdot_iau_earth
    end type iau_earth_rotating_frame
    interface iau_earth_rotating_frame
        module procedure :: iau_earth_rotating_frame_constructor
    end interface iau_earth_rotating_frame

    type,extends(iau_rotating_frame_class),public :: iau_moon_rotating_frame
        !! IAU Moon frame
    contains
        procedure :: get_c_cdot => get_c_cdot_iau_moon
    end type iau_moon_rotating_frame
    interface iau_moon_rotating_frame
        module procedure :: iau_moon_rotating_frame_constructor
    end interface iau_moon_rotating_frame

    !inertial frame definitions:

    type,extends(inertial_frame_class),public :: icrf_frame
        !! the fundamental inertial frame
        !! for the ephemeris (i.e., J2000).
    contains
        procedure :: get_c_cdot => get_c_cdot_icrf
    end type icrf_frame
    interface icrf_frame
        module procedure :: icrf_frame_constructor
    end interface icrf_frame

    type,extends(inertial_frame_class),public :: ecliptic_frame
        !! Mean ecliptic frame.
    contains
        procedure :: get_c_cdot => get_c_cdot_ecliptic
    end type ecliptic_frame
    interface ecliptic_frame
        module procedure :: ecliptic_frame_constructor
    end interface ecliptic_frame

    abstract interface
        subroutine c_cdot_func(me,eph,to_icrf,c,cdot,status_ok)
            import :: wp,reference_frame,ephemeris_class
            implicit none
            class(reference_frame),intent(inout)         :: me
            class(ephemeris_class),intent(inout)         :: eph     !! for ephemeris computations
                                                                    !! (assumed to have already
                                                                    !! been initialized)
            logical,intent(in)                           :: to_icrf
            real(wp),dimension(3,3),intent(out)          :: c
            real(wp),dimension(3,3),intent(out),optional :: cdot
            logical,intent(out)                          :: status_ok
        end subroutine c_cdot_func
    end interface

    !test routine:
    public :: transformation_module_test

    contains
!********************************************************************************

!********************************************************************************
!>
!  Constructor for a [[two_body_rotating_frame]]

    pure function two_body_rotating_frame_constructor(&
                primary_body,secondary_body,center,et) result(f)

    implicit none

    type(two_body_rotating_frame) :: f
    type(celestial_body),intent(in) :: primary_body   !! the primary body of
                                                      !! the frame
    type(celestial_body),intent(in) :: secondary_body !! the secondary body used
                                                      !! to construct the frame
    integer,intent(in)              :: center         !! the frame center (can
                                                      !! be `primary_body`,
                                                      !! `secondary_body`, or
                                                      !! `barycenter`)
    real(wp),intent(in)             :: et             !! epoch at which the
                                                      !! frame is defined [sec]

    f%primary_body   = primary_body
    f%secondary_body = secondary_body
    f%center         = center
    f%et             = et

    end function two_body_rotating_frame_constructor
!********************************************************************************

!********************************************************************************
!>
!  Constructor for a [[two_body_rotating_pulsating_frame]]

    pure function two_body_rotating_pulsating_frame_constructor(&
                primary_body,secondary_body,center,scale,et) result(f)

    implicit none

    type(two_body_rotating_pulsating_frame) :: f
    type(celestial_body),intent(in) :: primary_body   !! the primary body of
                                                      !! the frame
    type(celestial_body),intent(in) :: secondary_body !! the secondary body used
                                                      !! to construct the frame
    integer,intent(in)              :: center         !! the frame center (can
                                                      !! be `primary_body`,
                                                      !! `secondary_body`, or
                                                      !! `barycenter`)
    real(wp),intent(in)             :: scale          !! scale factor
    real(wp),intent(in)             :: et             !! epoch at which the
                                                      !! frame is defined [sec]

    f%primary_body   = primary_body
    f%secondary_body = secondary_body
    f%center         = center
    f%scale          = scale
    f%et             = et

    end function two_body_rotating_pulsating_frame_constructor
!********************************************************************************

!********************************************************************************
!>
!  Constructor for a [[iau_earth_rotating_frame]]

    pure function iau_earth_rotating_frame_constructor(b,et) result(f)

    implicit none

    type(iau_earth_rotating_frame)  :: f
    type(celestial_body),intent(in) :: b   !! the central body
    real(wp),intent(in)             :: et  !! epoch at which the
                                           !! frame is defined [sec]

    f%primary_body   = b
    f%et             = et

    end function iau_earth_rotating_frame_constructor
!********************************************************************************

!********************************************************************************
!>
!  Constructor for a [[iau_moon_rotating_frame]]

    pure function iau_moon_rotating_frame_constructor(b,et) result(f)

    implicit none

    type(iau_earth_rotating_frame)  :: f
    type(celestial_body),intent(in) :: b   !! the central body
    real(wp),intent(in)             :: et  !! epoch at which the
                                           !! frame is defined [sec]

    f%primary_body   = b
    f%et             = et

    end function iau_moon_rotating_frame_constructor
!********************************************************************************

!********************************************************************************
!>
!  Constructor for a [[icrf_frame]]
!
!@note the `et` doesn't matter for inertial frames

    pure function icrf_frame_constructor(b) result(f)

    implicit none

    type(icrf_frame)  :: f
    type(celestial_body),intent(in) :: b   !! the central body

    f%primary_body = b

    end function icrf_frame_constructor
!********************************************************************************

!********************************************************************************
!>
!  Constructor for a [[ecliptic_frame]]
!
!@note the `et` doesn't matter for inertial frames

    pure function ecliptic_frame_constructor(b) result(f)

    implicit none

    type(ecliptic_frame)  :: f
    type(celestial_body),intent(in) :: b   !! the central body

    f%primary_body = b

    end function ecliptic_frame_constructor
!********************************************************************************

!********************************************************************************
!>
!  Transform a Cartesian state from one reference frame to another at
!  a specified epoch. The `from` and `to` [[reference_frame]]s may each
!  be defined at a different epoch. The `et` ephemeris time is the time
!  the transformation is to be done, and accounts for the motion of the two
!  frame centers from `from%et` and `to%et` to `et`.

    subroutine transform(from,rv,to,et,eph,rv_out,status_ok)

    implicit none

    class(reference_frame),intent(inout) :: from
    real(wp),dimension(6),intent(in)     :: rv
    class(reference_frame),intent(inout) :: to
    real(wp),intent(in)                  :: et      !! the time of the transformation [sec]
    class(ephemeris_class),intent(inout) :: eph     !! for ephemeris computations (assumed to have already been initialized)
    real(wp),dimension(6),intent(out)    :: rv_out
    logical,intent(out)                  :: status_ok

    ! local variables
    real(wp), dimension(3,3) :: c, cdot
    real(wp), dimension(3,3) :: rot1, rotd1
    real(wp), dimension(3,3) :: rot2, rotd2
    real(wp), dimension(3)   :: rc21_out, vc21_out
    real(wp), dimension(3)   :: rc21_icrf, vc21_icrf

    ! rotation matrix: input -> inertial
    call from%get_c_cdot(eph=eph,to_icrf=.true.,c=rot1,cdot=rotd1,status_ok=status_ok)
    if (.not. status_ok) then
        write(error_unit,'(A)') 'Error in transform: '//&
                        'Could not compute rotation matrix from FROM frame to inertial.'
        rv_out = zero
        return
    end if
    ! rotation matrix: inertial -> output
    call to%get_c_cdot(eph=eph,to_icrf=.false.,c=rot2,cdot=rotd2,status_ok=status_ok)
    if (.not. status_ok) then
        write(error_unit,'(A)') 'Error in transform: '//&
                        'Could not compute rotation matrix from inertial to TO frame.'
        rv_out = zero
        return
    end if

    ! rotation matrix: input -> output
    c    = matmul(rot2, rot1)
    cdot = matmul(rotd2, rot1) + matmul(rot2, rotd1)

    ! get the state of the `from` frame center w.r.t. the `to` frame center,
    ! at the transformation time:
    call rvcto_rvcfrom_icrf(from, to, eph, et, &
                rc21_icrf, vc21_icrf, status_ok)

    if (status_ok) then

        ! to->from frame center state: inertial -> output frame
        rc21_out = matmul(rot2, rc21_icrf)
        vc21_out = matmul(rotd2, rc21_icrf) + matmul(rot2, vc21_icrf)

        ! rotation + translation:
        rv_out(1:3) = matmul(c, rv(1:3)) + rc21_out
        rv_out(4:6) = matmul(c, rv(4:6)) + matmul(cdot,rv(1:3)) + vc21_out

    else !error
        rv_out = zero
    end if

    end subroutine transform
!********************************************************************************

!********************************************************************************
    subroutine rvcto_rvcfrom_icrf(from, to, eph, et, rc21, vc21, status_ok)

    !! Returns the state of the `from` frame center w.r.t. the `to` frame center,
    !! at the specified ephemeris time `et`.

    implicit none

    class(reference_frame),intent(in)    :: from
    class(reference_frame),intent(in)    :: to
    class(ephemeris_class),intent(inout) :: eph       !! for ephemeris computations (assumed to have already been initialized)
    real(wp),intent(in)                  :: et        !! ephemeris time [sec]
    real(wp),dimension(3),intent(out)    :: rc21      !! position of `from` frame center w.r.t. `to` frame center
    real(wp),dimension(3),intent(out)    :: vc21      !! velocity of `from` frame center w.r.t. `to` frame center
    logical,intent(out)                  :: status_ok !! true if there were no errors

    ! local variables
    real(wp),dimension(6) :: rvc1  !! inertial state of `from` frame center w.r.t. `from` primary body
    real(wp),dimension(6) :: rvc2  !! inertial state of `to`   frame center w.r.t. `to`   primary body
    real(wp),dimension(6) :: rvb21 !! inertial state of `from` primary body w.r.t. `to`   primary body
    real(wp),dimension(6) :: rvc21 !! inertial state of `from` frame center w.r.t. `to`   frame center

    ! get TO primary body -> FROM primary body state [inertial]
    call eph%get_rv(et,from%primary_body,to%primary_body,rvb21,status_ok)
    if (status_ok) then

        ! currently, only the two-body rotating frames may be
        ! centered somewhere other than the primary body.
        select type (from)
        class is (two_body_rotating_frame)
            call from%from_primary_to_center(eph,et,rvc1,status_ok)
        class default
            rvc1 = zero
        end select

        if (status_ok) then

            select type (to)
            class is (two_body_rotating_frame)
                call to%from_primary_to_center(eph,et,rvc2,status_ok)
            class default
                rvc2 = zero
            end select

            if (status_ok) then
                rvc21 = rvb21 + rvc1 - rvc2
                rc21  = rvc21(1:3)
                vc21  = rvc21(4:6)
                return
            else
                write(error_unit,'(A)') 'Error in rvcto_rvcfrom_icrf: '//&
                                         'Could not compute center of TO frame.'
            end if

        else
            write(error_unit,'(A)') 'Error in rvcto_rvcfrom_icrf: '//&
                                     'Could not compute center of FROM frame.'
        end if

    else !error
        write(error_unit,'(A)') 'Error in rvcto_rvcfrom_icrf: '//&
                                 'Could not compute translation.'
    end if

    ! we end up here if there was an error:
    rc21  = zero
    vc21 = zero

    end subroutine rvcto_rvcfrom_icrf
!********************************************************************************

!********************************************************************************
    subroutine from_primary_to_center(me,eph,et,rc,status_ok)

    !! returns the state of the frame center w.r.t. the frame primary body.

    implicit none

    class(two_body_rotating_frame),intent(in) :: me
    class(ephemeris_class),intent(inout)      :: eph        !! for ephemeris computations (assumed to have already been initialized)
    real(wp),intent(in)                       :: et         !! ephemeris time [sec]
    real(wp),dimension(6),intent(out)         :: rc         !! state of frame center w.r.t. primary body [inertial]
    logical,intent(out)                       :: status_ok  !! true if no errors.

    real(wp) :: mu1 !! gravitational parameter of primary body
    real(wp) :: mu2 !! gravitational parameter of secondary body

    if (me%center == center_at_primary_body) then
        rc = zero
    else

        ! primary body -> secondary body state [inertial]:
        call eph%get_rv(et,me%secondary_body,me%primary_body,rc,status_ok)
        if (status_ok) then
            if (me%center == center_at_barycenter) then
                mu1 = me%primary_body%mu
                mu2 = me%secondary_body%mu
                rc = rc * ( mu2/(mu1 + mu2) )
            elseif (me%center == center_at_secondary_body) then
                !frame center is secondary body (rc already computed)
            else
                error stop 'invalid rotating frame center selection.'
            end if
        else
            write(error_unit,'(A)') 'Error in from_primary_to_center: '//&
                                    'Could not compute primary to secondary body state.'
        end if
    end if

    end subroutine from_primary_to_center
!********************************************************************************

!********************************************************************************
    subroutine get_c_cdot_two_body_rotating(me,eph,to_icrf,c,cdot,status_ok)

    !! rotation matrix for ROTATING <-> ICRF

    implicit none

    class(two_body_rotating_frame),intent(inout) :: me
    class(ephemeris_class),intent(inout)         :: eph     !! for ephemeris computations (assumed to have already been initialized)
    logical,intent(in)                           :: to_icrf
    real(wp),dimension(3,3),intent(out)          :: c
    real(wp),dimension(3,3),intent(out),optional :: cdot
    logical,intent(out)                          :: status_ok

    real(wp),dimension(3) :: r         !! position of secondary body w.r.t. primary body [inertial frame]
    real(wp),dimension(3) :: v         !! velocity of secondary body w.r.t. primary body [inertial frame]
    real(wp),dimension(3) :: h         !! angular momentum vector
    real(wp),dimension(3) :: w         !! angular velocity of frame
    logical               :: need_cdot !! if we need to compute `cdot`

    need_cdot = present(cdot)

    ! get position & velocity of secondary body w.r.t. primary body, in the inertial frame
    call eph%get_rv(me%et,me%secondary_body,me%primary_body,me%rv12,status_ok)

    if (status_ok) then

        r      = me%rv12(1:3)
        v      = me%rv12(4:6)
        h      = cross(r,v)
        c(1,:) = unit(r)
        c(3,:) = unit(h)
        c(2,:) = cross(c(3,:),c(1,:))

        if (need_cdot) then
            w = h / dot_product(r,r)            ! see: https://en.wikipedia.org/wiki/Angular_velocity
            cdot = -matmul(c,cross_matrix(w))   ! see: http://arxiv.org/pdf/1311.6010.pdf
        end if

        if (to_icrf) then
            c = transpose(c)
            if (need_cdot) cdot = transpose(cdot)
        end if

    else
        write(error_unit,'(A)') 'Error in get_c_cdot_two_body_rotating: '//&
                                'Could not compute rotation matrix.'
        c = zero
        if (need_cdot) cdot = zero
    end if

    end subroutine get_c_cdot_two_body_rotating
!********************************************************************************

!********************************************************************************
    subroutine get_c_cdot_two_body_rotating_pulsating(me,eph,to_icrf,c,cdot,status_ok)

    !! rotation matrix for ROTATING_PULSATING <-> ICRF

    implicit none

    class(two_body_rotating_pulsating_frame),intent(inout) :: me
    class(ephemeris_class),intent(inout)                   :: eph     !! for ephemeris computations (assumed to have already been initialized)
    logical,intent(in)                                     :: to_icrf
    real(wp),dimension(3,3),intent(out)                    :: c
    real(wp),dimension(3,3),intent(out),optional           :: cdot
    logical,intent(out)                                    :: status_ok

    ! local variables
    real(wp),dimension(3,3) :: cr, crdot
    real(wp),dimension(3) :: r12,v12
    real(wp) :: r12mag,factor
    logical :: need_cdot

    need_cdot = present(cdot)

    ! rotating frame transformation matrices:
    if (need_cdot) then
        call me%two_body_rotating_frame%get_c_cdot(eph,to_icrf,cr,crdot,status_ok=status_ok)
    else
        call me%two_body_rotating_frame%get_c_cdot(eph,to_icrf,cr,status_ok=status_ok)
    end if

    if (status_ok) then

        r12    = me%rv12(1:3) ! was computed in get_c_cdot_two_body_rotating
        v12    = me%rv12(4:6)
        r12mag = norm2(r12)

        if (to_icrf) then
            factor = r12mag/me%scale
            c = factor*cr
            if (need_cdot) cdot = factor*(dot_product(v12,r12)*cr/(r12mag**2)+crdot)
        else
            factor = me%scale/r12mag
            c = factor*cr
            if (need_cdot) cdot = factor*(-dot_product(v12,r12)*cr/(r12mag**2)+crdot)
        end if

    else  !error
        write(error_unit,'(A)') 'Error in get_c_cdot_two_body_rotating_pulsating: '//&
                                'Could not compute rotation matrix.'
        c = zero
        if (need_cdot) cdot = zero
    end if

    end subroutine get_c_cdot_two_body_rotating_pulsating
!********************************************************************************

!********************************************************************************
    subroutine get_c_cdot_icrf(me,eph,to_icrf,c,cdot,status_ok)

    !! rotation matrix for ICRF <-> ICRF

    implicit none

    class(icrf_frame),intent(inout)              :: me
    class(ephemeris_class),intent(inout)         :: eph     !! for ephemeris computations (assumed to have already been initialized)
    logical,intent(in)                           :: to_icrf
    real(wp),dimension(3,3),intent(out)          :: c
    real(wp),dimension(3,3),intent(out),optional :: cdot
    logical,intent(out)                          :: status_ok

    c = identity_3x3
    if (present(cdot)) cdot = zero
    status_ok = .true.

    end subroutine get_c_cdot_icrf
!********************************************************************************

!********************************************************************************
    subroutine get_c_cdot_ecliptic(me,eph,to_icrf,c,cdot,status_ok)

    !! rotation matrix for ICRF <-> Mean Ecliptic

    use obliquity_module

    implicit none

    class(ecliptic_frame),intent(inout)          :: me
    class(ephemeris_class),intent(inout)         :: eph
    logical,intent(in)                           :: to_icrf
    real(wp),dimension(3,3),intent(out)          :: c
    real(wp),dimension(3,3),intent(out),optional :: cdot
    logical,intent(out)                          :: status_ok

    if (to_icrf) then
        c = mean_ecliptic_to_equatorial_rotmat()
    else
        c = equatorial_to_mean_ecliptic_rotmat()
    end if

    if (present(cdot)) cdot = zero
    status_ok = .true.

    end subroutine get_c_cdot_ecliptic
!********************************************************************************

!********************************************************************************
    subroutine get_c_cdot_iau_earth(me,eph,to_icrf,c,cdot,status_ok)

    !! rotation matrix for IAU_EARTH <-> ICRF

    implicit none

    class(iau_earth_rotating_frame),intent(inout) :: me
    class(ephemeris_class),intent(inout)          :: eph     !! for ephemeris computations (assumed to have already been initialized)
    logical,intent(in)                            :: to_icrf
    real(wp),dimension(3,3),intent(out)           :: c
    real(wp),dimension(3,3),intent(out),optional  :: cdot
    logical,intent(out)                           :: status_ok

    c = icrf_to_iau_earth(me%et)

    !... don't have the cdot code yet... need to refactor iau code ...
    ! see also: ftp://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/tisbod.html

    !if (present(cdot)) cdot =
    error stop 'not yet supported'

    status_ok = .true.

    end subroutine get_c_cdot_iau_earth
!********************************************************************************

!********************************************************************************
    subroutine get_c_cdot_iau_moon(me,eph,to_icrf,c,cdot,status_ok)

    !! rotation matrix for IAU_MOON <-> ICRF

    implicit none

    class(iau_moon_rotating_frame),intent(inout) :: me
    class(ephemeris_class),intent(inout)         :: eph     !! for ephemeris computations (assumed to have already been initialized)
    logical,intent(in)                           :: to_icrf
    real(wp),dimension(3,3),intent(out)          :: c
    real(wp),dimension(3,3),intent(out),optional :: cdot
    logical,intent(out)                          :: status_ok

    c = icrf_to_iau_moon(me%et)

    !... don't have the cdot code yet... need to refactor iau code ...
    ! see also: ftp://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/tisbod.html

    !if (present(cdot)) cdot =
    error stop 'not yet supported'

    status_ok = .true.

    end subroutine get_c_cdot_iau_moon
!********************************************************************************

!********************************************************************************
!>
!  Transformation units test

    subroutine transformation_module_test()

    use jpl_ephemeris_module, only: jpl_ephemeris

    implicit none

    real(wp),dimension(6),parameter :: initial_state = [10000.0_wp,&
                                                        0.0_wp,&
                                                        0.0_wp,&
                                                        1.0_wp,&
                                                        2.0_wp,&
                                                        3.0_wp] !! km, km/s
    real(wp),parameter :: et = zero           !! ephemeris time [sec]
    real(wp),parameter :: scale = 384400.0_wp !! scale factor [km]
    character(len=*),parameter :: ephemeris_file_421 = '../eph/JPLEPH.421' !! JPL DE421 ephemeris file

    type(icrf_frame) :: from
    type(two_body_rotating_pulsating_frame) :: to
    type(jpl_ephemeris) :: eph421
    logical :: status_ok
    real(wp),dimension(6) :: rv_out
    type(ecliptic_frame) :: ecliptic_f
    real(wp),dimension(3,3) :: c
    integer :: i !! counter

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' transformation_module_test'
    write(*,*) '---------------'
    write(*,*) ''

    !open the ephemeris:
    call eph421%initialize(filename=ephemeris_file_421,status_ok=status_ok)
    if (.not. status_ok) error stop 'Error initializing DE421 ephemeris.'

    !initialize frames:
    !  [NOTE: need to make constructors for the frames]
    ! note: default is Earth-Moon-barycenter:
    from = icrf_frame(et=et)
    to   = two_body_rotating_pulsating_frame(et=et,scale=scale)

    call from%transform(initial_state,to,et,eph421,rv_out,status_ok)
    if (.not. status_ok) error stop 'Error in state transformation.'

    !results:
    write(*,*) ''
    write(*,'(A/,*(E30.16/))') 'initial state (J2000-Earth):', initial_state
    write(*,'(A/,*(E30.16/))') 'final state (Earth-Moon rotating, centered at barycenter, scale=384400):', rv_out
    write(*,*) ''

    ! ecliptic frame:
    ecliptic_f = ecliptic_frame(b=body_earth)
    call ecliptic_f%get_c_cdot(eph421,to_icrf=.true.,c=c,status_ok=status_ok)
    write(*,*) ''
    write(*,*) 'ecliptic to j2000:'
    do i=1,3
        write(*,*) c(i,:)
    end do

    ! from SPICE:
    ! rot = [1.0000000000000000E+00, 0.0000000000000000E+00, 0.0000000000000000E+00,
    !        0.0000000000000000E+00, 9.1748206206918181E-01, -3.9777715593191371E-01,
    !        0.0000000000000000E+00, 3.9777715593191371E-01, 9.1748206206918181E-01]
    !
    ! from FAT:
    !        1.0000000000000000        0.0000000000000000        0.0000000000000000
    !        0.0000000000000000       0.91748206206918181      -0.39777715593191371
    !        0.0000000000000000       0.39777715593191371       0.91748206206918181

    !close the ephemeris:
    call eph421%close()

    end subroutine transformation_module_test
!********************************************************************************

!********************************************************************************
    end module transformation_module
!********************************************************************************
