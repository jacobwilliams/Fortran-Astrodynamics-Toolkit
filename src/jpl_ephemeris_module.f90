!*****************************************************************************************
!>
!  For reading the JPL planetary and lunar ephemerides.
!  This is an extensively modified version of the original FORTRAN 77 code from JPL.
!
!### Ephemeris files
!  Note that this module uses the JPL binary ephemeris files, which
!  can be obtained using the instructions
!  [here](ftp://ssd.jpl.nasa.gov/pub/eph/planets/fortran/userguide.txt).
!  See also the comments in [[ephemeris_test]] for more details.
!
!### License
!
!  ***Original JPL License***
!
!  THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
!  CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
!  GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
!  ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
!  PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS"
!  TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY
!  WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A
!  PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
!  SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
!  SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
!
!  IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA
!  BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT
!  LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND,
!  INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS,
!  REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE
!  REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!
!  RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF
!  THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY
!  CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE
!  ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
!
!  ***Modifications***
!
!  Modifications for the Fortran Astrodynamics Toolkit are covered
!  under the [following license](https://github.com/jacobwilliams/Fortran-Astrodynamics-Toolkit/blob/master/LICENSE).
!
!### History
!  * [Original code from JPL](ftp://ssd.jpl.nasa.gov/pub/eph/planets/fortran/) Version : March 25, 2013
!  * Extensive modifications by Jacob Williams for the Fortran Astrodynamics Toolkit.

    module jpl_ephemeris_module

    use, intrinsic :: iso_fortran_env, only: real64,error_unit
    use ephemeris_module

    implicit none

    private

    integer,parameter :: wp = real64  !! using double precision reals since that is how the ephemeris is stored.

    integer,parameter,public :: nmax   = 1000   !! Current maximum number of ephemeris constants
                                                !! used in the integration and listed
                                                !! in the `header.xxx` file for the ephemeris.
    integer,parameter        :: oldmax = 400    !! For earlier versions of the code, the maximum
                                                !! number of ephemeris constants used in the
                                                !! integration and listed in the `header.xxx`
                                                !! file for the ephemeris.

    integer,parameter :: nrecl = 4  !! `nrecl=1` if `recl` in the open statement is the record length in s.p. words
                                    !! `nrecl=4` if `recl` in the open statement is the record length in bytes

    type,extends(ephemeris_class),public :: jpl_ephemeris

        !! Main class for accessing a JPL ephemeris file.

        character(len=:),allocatable :: namfil  !! name of the binary ephemeris file

        integer :: ksize = 2036  !! ksize must be set by the user according to the ephemeris to be read:
                                 !!   for ***de200***, set `ksize=1652`,
                                 !!   for ***de405***, set `ksize=2036`,
                                 !!   for ***de406***, set `ksize=1456`,
                                 !!   for ***de414***, set `ksize=2036`,
                                 !!   for ***de418***, set `ksize=2036`,
                                 !!   for ***de421***, set `ksize=2036`,
                                 !!   for ***de422***, set `ksize=2036`,
                                 !!   for ***de423***, set `ksize=2036`,
                                 !!   for ***de424***, set `ksize=2036`,
                                 !!   for ***de430***, set `ksize=2036`.

        !ephhdr
        integer,dimension(3,13)  :: ipt   = 0       !! ipt(39)
        real(wp),dimension(nmax) :: cval  = 0.0_wp
        real(wp),dimension(3)    :: ss    = 0.0_wp
        real(wp)                 :: au    = 0.0_wp
        real(wp)                 :: emrat = 0.0_wp
        integer                  :: ncon  = 0
        integer                  :: numde = 0

        !stcomx
        real(wp),dimension(6) :: pvsun = 0.0_wp   !! dp 6-word array containing the barycentric position and
                                                  !! velocity of the sun.
        logical               :: km    = .true.   !! logical flag defining physical units of the output states.
                                                  !!   km = .true.  : km and km/sec
                                                  !!   km = .false. : au and au/day
                                                  !! for nutations and librations.  angle unit is always radians.
        logical               :: bary  = .false.  !! logical flag defining output center.
                                                  !! only the 9 planets are affected.
                                                  !!   bary = .true.  : center is solar-system barycenter
                                                  !!   bary = .false. : center is sun

        !chrhdr
        character(len=6),dimension(14,3) :: ttl  = ''
        character(len=6),dimension(nmax) :: cnam = ''

        logical :: initialized = .false.  !! is the ephemeris initialized?
        integer :: nrfile      = 0        !! file unit for the ephemeris file
        integer :: nrl         = -1       !! this was formerly in state
        integer :: ncoeffs     = 0

        ! formerly in state:
        real(wp),dimension(1500) :: buf = 0.0_wp

        ! formerly in interp:
        real(wp),dimension(18) :: pc = 0.0_wp
        real(wp),dimension(18) :: vc = 0.0_wp
        integer                :: np = 2
        integer                :: nv = 3
        real(wp)               :: twot = 0.0_wp

    contains

        procedure,public :: get_rv => get_rv_from_jpl_ephemeris

        procedure,public :: initialize => initialize_ephemeris
        procedure,public :: get_state
        procedure,public :: get_constants
        procedure,public :: close => close_ephemeris

        procedure :: interp
        procedure :: state

    end type jpl_ephemeris

    !>
    !  These correspond to the numbering convention for 'ntarg' and 'ncent' in [[get_state]]:
    character(len=*),dimension(15),parameter :: list_of_bodies = [ &
                                                        'mercury                ',&
                                                        'venus                  ',&
                                                        'earth                  ',&
                                                        'mars                   ',&
                                                        'jupiter                ',&
                                                        'saturn                 ',&
                                                        'uranus                 ',&
                                                        'neptune                ',&
                                                        'pluto                  ',&
                                                        'moon                   ',&
                                                        'sun                    ',&
                                                        'solar-system barycenter',&
                                                        'earth-moon barycenter  ',&
                                                        'nutations              ',&
                                                        'librations             ']

    !public routines:
    public :: ephemeris_test

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Interface for the [[ephemeris_module]].

    subroutine get_rv_from_jpl_ephemeris(me,et,targ,obs,rv,status_ok)

    use time_module,           only: et_to_jd
    use conversion_module,     only: day2sec
    use celestial_body_module, only: celestial_body
    use numbers_module,        only: zero

    implicit none

    class(jpl_ephemeris),intent(inout) :: me
    real(wp),intent(in)                :: et         !! ephemeris time [sec]
    type(celestial_body),intent(in)    :: targ       !! target body
    type(celestial_body),intent(in)    :: obs        !! observer body
    real(wp),dimension(6),intent(out)  :: rv         !! state of targ w.r.t. obs [km,km/s] in ICRF frame
    logical,intent(out)                :: status_ok  !! true if there were no problems

    real(wp) :: jd     !! julian date for input to [[get_state]].
    integer :: ntarg   !! id code for target body
    integer :: ncent   !! id code for observer body

    if (targ==obs) then
        !don't bother if target and observer are the same body
        rv = zero
        status_ok = .true.
    else

        !convert to expected inputs:
        jd    = et_to_jd(et)
        ntarg = spice_id_to_old_id(targ%id)
        ncent = spice_id_to_old_id(obs%id)

        if (ntarg>0 .and. ncent>0) then
            call me%get_state(jd,ntarg,ncent,rv,status_ok)
            if (status_ok) then
                if (.not. me%km) then
                    !we must return in units of km/s
                    !so, convert from AU, AU/day to km, km/s
                    rv      = rv * me%au
                    rv(4:6) = rv(4:6) / day2sec
                end if
            else
                write(error_unit,'(A)') 'Error in get_rv_from_jpl_ephemeris: '//&
                            'Error calling ephemeris.'
            end if
        else
            write(error_unit,'(A)') 'Error in get_rv_from_jpl_ephemeris: '//&
                        'No ephemeris for this body.'
            status_ok = .false.
        end if

    end if

    end subroutine get_rv_from_jpl_ephemeris
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 3/20/2016
!
!  Convert the NAIF SPICE ID code to the old one used by the JPL ephemeris.
!  Returns `0` if the body was not found.

    pure function spice_id_to_old_id(spice_id) result(old_id)

    implicit none

    integer,intent(in) :: spice_id !! the ID code used by SPICE
    integer            :: old_id   !! the ID code used by this module (old JPL ephemeris code)

    integer :: i !! counter

    !>
    !  The index of this array is the old ID code. The value is the new code.
    !  See: [NAIF Integer ID codes](http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/naif_ids.html)
    integer,parameter,dimension(13) :: new_ids = &
        [   199,&  ! mercury
            299,&  ! venus
            399,&  ! earth
            499,&  ! mars
            599,&  ! jupiter
            699,&  ! saturn
            799,&  ! uranus
            899,&  ! neptune
            999,&  ! pluto
            301,&  ! moon
            10, &  ! sun
            0,  &  ! solar-system barycenter
            3   ]  ! earth-moon barycenter

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

    end function spice_id_to_old_id
!*****************************************************************************************

!*****************************************************************************************
!>
!  This subroutine reads the JPL planetary ephemeris
!  and gives the position and velocity of the point `ntarg`
!  with respect to `ncent`.
!
!### Notes
!
!  The numbering convention for `ntarg` and `ncent` is:
!
!    1 = mercury      8 = neptune
!    2 = venus        9 = pluto
!    3 = earth       10 = moon
!    4 = mars        11 = sun
!    5 = jupiter     12 = solar-system barycenter
!    6 = saturn      13 = earth-moon barycenter
!    7 = uranus      14 = nutations (longitude and obliq)
!                    15 = librations, if on eph file
!
!  (if nutations are wanted, set ntarg = 14.
!  for librations, set ntarg = 15. set ncent=0.)

    subroutine get_state(me,jd,ntarg,ncent,rrd,status_ok)

    implicit none

    class(jpl_ephemeris),intent(inout) :: me
    real(wp),intent(in)                :: jd         !! d.p. Julian ephemeris date at which interpolation is wanted.
    integer,intent(in)                 :: ntarg      !! integer number of 'target' point.
    integer,intent(in)                 :: ncent      !! integer number of 'center' point.
    real(wp),dimension(6),intent(out)  :: rrd        !! output 6-word d.p. array containing position and velocity
                                                     !! of point `ntarg` relative to `ncent`.
                                                     !! the units are AU and AU/day (or km and km/sec if `me%km=.true.`).
                                                     !! For librations the units are radians and radians
                                                     !! per day. In the case of nutations the first four words of
                                                     !! `rrd` will be set to nutations and rates, having units of
                                                     !! radians and radians/day.
    logical,intent(out)                :: status_ok  !! true if there were no problems

    real(wp),dimension(2)    :: et2
    real(wp),dimension(6,13) :: pv
    real(wp),dimension(6,11) :: pvst
    real(wp),dimension(4)    :: pnut
    integer,dimension(12)    :: list
    integer                  :: i,j,k
    logical                  :: bsave

    status_ok = .false.

    if (me%initialized) then

        ! initialize et2 for 'state' and set up component count
        et2(1) = jd
        et2(2) = 0.0_wp
        list   = 0
        rrd    = 0.0_wp

        if (ntarg /= ncent) then

            select case (ntarg)

            case (14)    !nutation

                if (me%ipt(2,12)>0) then            !me%ipt(35)
                    list(11) = 2
                    call me%state(et2,list,pvst,pnut,status_ok)
                    if (.not. status_ok) then
                        write(error_unit,'(A)') 'error interpolating nutations in get_state.'
                        return
                    end if
                    rrd(1:4) = pnut(1:4)
                    rrd(5) = 0.0_wp
                    rrd(6) = 0.0_wp
                    return
                else
                    rrd(1:4) = 0.0_wp
                    write(error_unit,'(A)') 'error in get_state: the ephemeris file does not contain nutations.'
                    return
                endif

            case (15)    !librations

                if (me%ipt(2,13)>0) then            !me%ipt(38)
                    list(12) = 2
                    call me%state(et2,list,pvst,pnut,status_ok)
                    if (.not. status_ok) then
                        write(error_unit,'(A)') 'error interpolating librations in get_state.'
                        return
                    end if
                    rrd = pvst(:,11)
                    return
                else
                    write(error_unit,'(A)') 'error in get_state: the ephemeris file does not contain librations.'
                    return
                endif

               case default

                ! force barycentric output by 'state'

                bsave = me%bary
                me%bary  = .true.

                ! set up proper entries in 'list' array for state call

                do i=1,2
                    k=ntarg
                    if (i == 2)  k = ncent
                    if (k <= 10) list(k)  = 2
                    if (k == 10) list(3)  = 2
                    if (k == 3)  list(10) = 2
                    if (k == 13) list(3)  = 2
                enddo

                ! make call to state

                call me%state(et2,list,pvst,pnut,status_ok)
                if (.not. status_ok) then
                    write(error_unit,'(A)') 'error interpolating state in get_state.'
                    return
                end if
                do i=1,10
                    do j = 1,6
                        pv(j,i) = pvst(j,i)
                    end do
                enddo

                if (ntarg == 11 .or. ncent == 11) pv(:,11) = me%pvsun
                if (ntarg == 12 .or. ncent == 12) pv(:,12) = 0.0_wp
                if (ntarg == 13 .or. ncent == 13) pv(:,13) = pvst(:,3)

                if (ntarg*ncent == 30 .and. ntarg+ncent == 13) then
                    pv(:,3) = 0.0_wp
                else
                    if (list(3) == 2)  pv(:,3)  = pvst(:,3) - pvst(:,10)/(1.0_wp+me%emrat)
                    if (list(10) == 2) pv(:,10) = pv(:,3) + pvst(:,10)
                end if

                rrd = pv(:,ntarg) - pv(:,ncent)

                me%bary = bsave

            end select

        end if

        status_ok = .true.

    else
        write(error_unit,'(A)') 'error in get_state: the ephemeris is not initialized.'
    end if

    end subroutine get_state
!*****************************************************************************************

!*****************************************************************************************
!>
!  this subroutine differentiates and interpolates a
!  set of chebyshev coefficients to give position and velocity.

    subroutine interp(me,buf,t,ncf,ncm,na,ifl,pv)

    implicit none

    class(jpl_ephemeris),intent(inout) :: me
    integer,intent(in)                 :: ncf  !! # of coefficients per component
    integer,intent(in)                 :: ncm  !! # of components per set of coefficients
    real(wp),dimension(ncf,ncm,*)      :: buf  !! 1st location of array of d.p. chebyshev coefficients of position
    real(wp),dimension(2),intent(in)   :: t
    integer,intent(in)                 :: na   !! # of sets of coefficients in full array
                                               !! (i.e., # of sub-intervals in full interval)
    integer,intent(in)                 :: ifl  !! integer flag
                                               !! = 1 for positions only
                                               !! = 2 for pos and vel
    real(wp),dimension(ncm,*)          :: pv   !! interpolated quantities requested.  dimension
                                               !! expected is pv(ncm,ifl), dp.

    real(wp) :: dna,dt1,temp,vfac,tc
    integer :: l,i,j

    ! entry point. get correct sub-interval number for this set
    ! of coefficients and then get normalized chebyshev time
    ! within that subinterval.

    dna  = dble(na)
    dt1  = int(t(1))
    temp = dna*t(1)
    l    = int(temp-dt1)+1

    ! tc is the normalized chebyshev time (-1 <= tc <= 1)

    tc=2.0_wp*(mod(temp,1.0_wp)+dt1)-1.0_wp

    ! check to see whether chebyshev time has changed,
    ! and compute new polynomial values if it has.
    ! (the element pc(2) is the value of t1(tc) and hence
    ! contains the value of tc on the previous call.)

    if (tc/=me%pc(2)) then
        me%np    = 2
        me%nv    = 3
        me%pc(2) = tc
        me%twot  = tc+tc
    endif

    ! be sure that at least 'ncf' polynomials have been evaluated
    ! and are stored in the array 'pc'.

    if (me%np<ncf) then
        do i=me%np+1,ncf
            me%pc(i)=me%twot*me%pc(i-1)-me%pc(i-2)
        end do
        me%np=ncf
    endif

    ! interpolate to get position for each component

    do i=1,ncm
        pv(i,1)=0.0_wp
        do j=ncf,1,-1
            pv(i,1)=pv(i,1)+me%pc(j)*buf(j,i,l)
        end do
    end do
    if (ifl<=1) return

    ! if velocity interpolation is wanted, be sure enough
    ! derivative polynomials have been generated and stored.

    vfac=(dna+dna)/t(2)
    me%vc(3)=me%twot+me%twot
    if (me%nv<ncf) then
        do i=me%nv+1,ncf
            me%vc(i)=me%twot*me%vc(i-1)+me%pc(i-1)+me%pc(i-1)-me%vc(i-2)
        end do
        me%nv=ncf
    endif

    ! interpolate to get velocity for each component

    do i=1,ncm
        pv(i,2)=0.0_wp
        do j=ncf,2,-1
            pv(i,2)=pv(i,2)+me%vc(j)*buf(j,i,l)
        end do
        pv(i,2)=pv(i,2)*vfac
    end do

    end subroutine interp
!*****************************************************************************************

!*****************************************************************************************
!>
!  this subroutine breaks a d.p. number into a d.p. integer
!  and a d.p. fractional part.

    subroutine split(tt,fr)

    implicit none

    real(wp),intent(in)               :: tt  !! d.p. input number
    real(wp),dimension(2),intent(out) :: fr  !! d.p. 2-word output array.
                                             !! fr(1) contains integer part
                                             !! fr(2) contains fractional part
                                             !! for negative input numbers, fr(1) contains the next
                                             !! more negative integer; fr(2) contains a positive fraction.

    ! get integer and fractional parts

    fr(1) = int(tt)
    fr(2) = tt-fr(1)

    if (tt>=0.0_wp .or. fr(2)==0.0_wp) return

    ! make adjustments for negative input number

    fr(1) = fr(1) - 1.0_wp
    fr(2) = fr(2) + 1.0_wp

    end subroutine split
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize the ephemeris.
!  This routine may be called to load a different ephemeris file.
!  Otherwise, it is called on the first call to get_state, and loads
!  the file specified in the module header.
!
!# Note
!  * Based on code formerly in [[state]].

    subroutine initialize_ephemeris(me,filename,ksize,km,bary,status_ok)

    implicit none

    class(jpl_ephemeris),intent(inout) :: me
    character(len=*),intent(in)        :: filename   !! ephemeris file name
    integer,intent(in),optional        :: ksize      !! corresponding `ksize`
    logical,intent(in),optional        :: km         !! defining physical units of the output states.
                                                     !! `km = .true.`  : km and km/sec [default],
                                                     !! `km = .false.` : au and au/day.
    logical,intent(in),optional        :: bary       !! logical flag defining output center.
                                                     !! only the 9 planets are affected.
                                                     !! `bary = .true.`  : center is solar-system barycenter,
                                                     !! `bary = .false.` : center is sun [default].
    logical,intent(out) :: status_ok                 !! true if there were not problems.

    !local variables:
    integer :: irecsz,istat,i,j,k,l

    !just in case it was already open:
    call me%close()  ! clears everything in the class

    !ephemeris file name:
    me%namfil = trim(filename)

    !optional inputs:
    if (present(ksize)) me%ksize  = ksize
    if (present(km))    me%km     = km
    if (present(bary))  me%bary   = bary

    irecsz     = nrecl*me%ksize
    me%ncoeffs = me%ksize/2

    open(newunit = me%nrfile,     &
         file    = me%namfil,     &
         access  = 'DIRECT',      &
         form    = 'UNFORMATTED', &
         action  = 'READ',        &   !JW added
         recl    = irecsz,        &
         iostat  = istat,         &
         status  = 'OLD'          )

         !write(*,*) "istat=",istat

    status_ok = (istat==0)  !if there were no problems opening the file

    if (status_ok) then

        read(me%nrfile,rec=1,iostat=istat) &
                 me%ttl,(me%cnam(k),k=1,oldmax),me%ss,me%ncon,me%au,me%emrat,&
                 ((me%ipt(i,j),i=1,3),j=1,12),me%numde,(me%ipt(i,13),i=1,3), &
                 (me%cnam(l),l=oldmax+1,me%ncon)

        if (istat==0) then
            if (me%ncon <= oldmax) then
                read(me%nrfile,rec=2,iostat=istat) (me%cval(i),i=1,oldmax)
            else
                read(me%nrfile,rec=2,iostat=istat) (me%cval(i),i=1,me%ncon)
            endif
            if (istat==0) then
                me%nrl = 0
                me%initialized = .true.
            end if
        end if

        ! check if the reads went OK:
        status_ok = me%initialized
        if (status_ok) then
            !initialize some of the class variables:
            ! [note: this was formerly done in the interp routine]
            me%pc(1) = 1.0_wp
            me%vc(2) = 1.0_wp
        else
            write(error_unit,'(A)') 'Error reading ephemeris file: '//trim(me%namfil)
        end if

    else
        write(error_unit,'(A)') 'Error opening ephemeris file: '//trim(me%namfil)
    end if

    end subroutine initialize_ephemeris
!*****************************************************************************************

!*****************************************************************************************
!>
!  Close the ephemeris.

    subroutine close_ephemeris(me)

    implicit none

    class(jpl_ephemeris),intent(inout) :: me

    logical :: is_open
    integer :: istat

    if (me%initialized) then

        !close the file:
        inquire(unit=me%nrfile,opened=is_open,iostat=istat)
        if (is_open) close(unit=me%nrfile,iostat=istat)

        !initialize all class variables to defaults:
        call clear(me)

    end if

    contains

        subroutine clear(eph)
            !! clear all the variables in the [[jpl_ephemeris]] class.
            implicit none
            class(jpl_ephemeris),intent(out) :: eph
        end subroutine clear

    end subroutine close_ephemeris
!*****************************************************************************************

!*****************************************************************************************
!>
!  This subroutine reads and interpolates the JPL planetary ephemeris file.
!
!@note The ephemeris is assumed to have been initialized.

    subroutine state(me,et2,list,pv,pnut,status_ok)

    implicit none

    class(jpl_ephemeris),intent(inout)   :: me
    real(wp),dimension(2),intent(in)     :: et2     !! 2-word Julian ephemeris epoch at which interpolation
                                                    !! is wanted.  any combination of `et2(1)+et2(2)` which falls
                                                    !! within the time span on the file is a permissible epoch.
                                                    !!
                                                    !! ***a.*** for ease in programming, the user may put the
                                                    !!    entire epoch in `et2(1)` and set `et2(2)=0`.
                                                    !! ***b.*** for maximum interpolation accuracy, set `et2(1)` =
                                                    !!    the most recent midnight at or before interpolation
                                                    !!    epoch and set `et2(2)` = fractional part of a day
                                                    !!    elapsed between `et2(1)` and epoch.
                                                    !! ***c.*** as an alternative, it may prove convenient to set
                                                    !!    `et2(1)` = some fixed epoch, such as start of integration,
                                                    !!    and `et2(2)` = elapsed interval between then and epoch.

    integer,dimension(12),intent(in)     :: list    !! 12-word integer array specifying what interpolation
                                                    !! is wanted for each of the bodies on the file.
                                                    !!    ***list(i) = 0*** : no interpolation for body i,
                                                    !!    ***list(i) = 1*** : position only,
                                                    !!    ***list(i) = 2*** : position and velocity.
                                                    !!
                                                    !! The designation of the astronomical bodies by i is:
                                                    !!    ***i = 1*** : mercury,
                                                    !!    ***i = 2*** : venus,
                                                    !!    ***i = 3*** : earth-moon barycenter,
                                                    !!    ***i = 4*** : mars,
                                                    !!    ***i = 5*** : jupiter,
                                                    !!    ***i = 6*** : saturn,
                                                    !!    ***i = 7*** : uranus,
                                                    !!    ***i = 8*** : neptune,
                                                    !!    ***i = 9*** : pluto,
                                                    !!    ***i =10*** : geocentric moon,
                                                    !!    ***i =11*** : nutations in longitude and obliquity,
                                                    !!    ***i =12*** : lunar librations (if on file).

    real(wp),dimension(6,11),intent(out) :: pv      !! dp 6 x 11 array that will contain requested interpolated
                                                    !! quantities (other than nutation, stored in `pnut`).
                                                    !! the body specified by `list(i)` will have its
                                                    !! state in the array starting at `pv(1,i)`.
                                                    !! (on any given call, only those words in `pv` which are
                                                    !! affected by the first 10 `list` entries, and by `list(12)`
                                                    !! if librations are on the file, are set.
                                                    !! the rest of the `pv` array is untouched.)
                                                    !! the order of components starting in `pv(1,i)` is:
                                                    !! `x`,`y`,`z`,`dx`,`dy`,`dz`.
                                                    !!
                                                    !! All output vectors are referenced to the earth mean
                                                    !! equator and equinox of J2000 if the DE number is 200 or
                                                    !! greater; of B1950 if the DE number is less than 200.
                                                    !!
                                                    !! The moon state is always geocentric; the other nine states
                                                    !! are either heliocentric or solar-system barycentric,
                                                    !! depending on the setting of the `bary` variable in the class.
                                                    !!
                                                    !! Lunar librations, if on file, are put into `pv(k,11)` if
                                                    !! `list(12)` is `1` or `2`.

    real(wp),dimension(4),intent(out)    :: pnut    !! dp 4-word array that will contain nutations and rates,
                                                    !! depending on the setting of `list(11)`.  the order of
                                                    !! quantities in `pnut` is:
                                                    !!
                                                    !! * `d psi`  (nutation in longitude),
                                                    !! * `d epsilon` (nutation in obliquity),
                                                    !! * `d psi dot`,
                                                    !! * `d epsilon dot`.
    logical,intent(out) :: status_ok !! true if there were no problems

    real(wp),dimension(2)    :: t
    real(wp),dimension(4)    :: pjd
    real(wp) :: aufac,s,tmp1,tmp2
    integer :: istat,i,j,k,nr

    status_ok = .true.

    s = et2(1)-0.5_wp
    call split(s,pjd(1:2))
    call split(et2(2),pjd(3:4))
    pjd(1) = pjd(1)+pjd(3)+0.5_wp
    pjd(2) = pjd(2)+pjd(4)
    call split(pjd(2),pjd(3:4))
    pjd(1) = pjd(1)+pjd(3)

    ! error return for epoch out of range

    if (pjd(1)+pjd(4)<me%ss(1) .or. pjd(1)+pjd(4)>me%ss(2)) then
        write(error_unit,'(A,F12.2,A,2F12.2)') &
                'Error: requested jed,',&
                et2(1)+et2(2),&
                ' not within ephemeris limits,',&
                me%ss(1),me%ss(2)
        status_ok = .false.
        return
    end if

    ! calculate record # and relative time in interval

    nr = int((pjd(1)-me%ss(1))/me%ss(3))+3
    if (pjd(1)==me%ss(2)) nr = nr-1

    tmp1 = dble(nr-3)*me%ss(3) + me%ss(1)
    tmp2 = pjd(1) - tmp1
    t(1) = (tmp2 + pjd(4))/me%ss(3)

    ! read correct record if not in core

    if (nr/=me%nrl) then
        me%nrl = nr
        read(me%nrfile,rec=nr,iostat=istat) (me%buf(k),k=1,me%ncoeffs)
        if (istat/=0) then
            write(error_unit,'(2F12.2,A80)') et2,'Error return in state'
            status_ok = .false.
            return
        end if
    endif

    if (me%km) then
        t(2) = me%ss(3)*86400.0_wp
        aufac = 1.0_wp
    else
        t(2) = me%ss(3)
        aufac = 1.0_wp/me%au
    endif

    ! interpolate ssbary sun

    call me%interp(me%buf(me%ipt(1,11)),t,me%ipt(2,11),3,me%ipt(3,11),2,me%pvsun)

    me%pvsun = me%pvsun * aufac

    ! check and interpolate whichever bodies are requested

    do i=1,10
        if (list(i)==0) cycle
        call me%interp(me%buf(me%ipt(1,i)),t,me%ipt(2,i),3,me%ipt(3,i),list(i),pv(1,i))
        do j=1,6
            if (i<=9 .and. .not.me%bary) then
                pv(j,i) = pv(j,i)*aufac - me%pvsun(j)
            else
                pv(j,i) = pv(j,i)*aufac
            endif
        enddo
    end do

    ! do nutations if requested (and if on file)

    if (list(11)>0 .and. me%ipt(2,12)>0) &
        call me%interp(me%buf(me%ipt(1,12)),t,me%ipt(2,12),2,me%ipt(3,12),list(11),pnut)

    ! get librations if requested (and if on file)

    if (list(12)>0 .and. me%ipt(2,13)>0) &
        call me%interp(me%buf(me%ipt(1,13)),t,me%ipt(2,13),3,me%ipt(3,13),list(12),pv(1,11))

    end subroutine state
!*****************************************************************************************

!*****************************************************************************************
!>
!  Obtain the constants from the ephemeris file.

    subroutine get_constants(me,nam,val,sss,n)

    implicit none

    class(jpl_ephemeris),intent(inout)        :: me
    character(len=6),dimension(:),intent(out) :: nam  !! array of constant names
    real(wp),dimension(:),intent(out)         :: val  !! array of values of constants
    real(wp),dimension(3),intent(out)         :: sss  !! jd start, jd stop, step of ephemeris
    integer,intent(out)                       :: n    !! number of entries in `nam` and `val` arrays

    integer :: i

    if (me%initialized) then

        n   = me%ncon
        sss = me%ss

        do i=1,n
            nam(i) = me%cnam(i)
            val(i) = me%cval(i)
        enddo

    else
        write(error_unit,'(A)') 'error in get_constants: the ephemeris is not initialized.'
    end if

    end subroutine get_constants
!*****************************************************************************************

!*****************************************************************************************
!>
!  Ephemeris test routine.
!
!### Ephemeris files
!  This routine requires the DE405 and DE421 JPL binary ephemeris files
!  to be present in the `../eph` directory.
!  These can be built by using the instructions
!  [here](ftp://ssd.jpl.nasa.gov/pub/eph/planets/fortran/userguide.txt).
!
!  For example (on Linux):
!```bash
!  wget ftp://ssd.jpl.nasa.gov/pub/eph/planets/fortran/*
!  wget ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/de405/*
!  #edit asc2eph.f file to set NRECL = 4:
!  sed -i '_original' '/^C.*PARAMETER ( NRECL = 4 )/s/^C//' asc2eph.f
!  gfortran asc2eph.f -o asc2eph
!  cat header.405 ascp*.405 | ./asc2eph
!  mkdir Fortran-Astrodynamics-Toolkit/eph
!  mv JPLEPH Fortran-Astrodynamics-Toolkit/eph/JPLEPH.405
!```

    subroutine ephemeris_test()

    implicit none

    character(len=6),dimension(nmax) :: nams
    real(wp) :: jd
    real(wp),dimension(6) :: rv,rv1,rv2,diffrv
    real(wp),dimension(3) :: ss
    real(wp),dimension(nmax) :: vals
    integer :: nvs,ntarg,nctr,i,j
    type(jpl_ephemeris) :: eph405, eph421
    logical :: status_ok_405,status_ok_421

    character(len=*),parameter :: ephemeris_file_405 = '../eph/JPLEPH.405' !! JPL DE405 ephemeris file
    character(len=*),parameter :: ephemeris_file_421 = '../eph/JPLEPH.421' !! JPL DE421 ephemeris file

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' ephemeris_test'
    write(*,*) '---------------'
    write(*,*) ''

    !initialize:
    call eph405%initialize(filename=ephemeris_file_405,status_ok=status_ok_405)
    call eph421%initialize(filename=ephemeris_file_421,status_ok=status_ok_421)

    if (status_ok_405) then

        !get some constants from the file:
        call eph405%get_constants(nams,vals,ss,nvs)

        write(*,'(A)') ''
        write(*,'(A)') 'Ephemeris initialized'
        write(*,'(A,1X,F15.3,1X,A,1X,F15.3)') 'JD range: ',ss(1), 'to ', ss(2)
        write(*,'(A)') ''
        do i=1,nvs
            write(*,'(A,1X,D25.16)') nams(i), vals(i)
        end do

        jd = 2451536.5d0       ! julian date

        if (jd < ss(1) .or. jd > ss(2)) then

            write(*,'(A)') ''
            write(*,*) 'error: jed out of bounds.'
            write(*,*) 'jed   = ', jd
            write(*,*) 'ss(1) = ', ss(1)
            write(*,*) 'ss(2) = ', ss(2)

        else

            !test DE405:
            do j=1,2

                if (j==1) then
                    ntarg = 3      !earth
                    nctr  = 11     !sun
                else
                    ntarg = 10     !moon
                    nctr  = 3      !earth
                end if

                write(*,*) ''
                write(*,*) 'DE405'
                write(*,*) 'state of "'//trim(list_of_bodies(ntarg))//&
                            '" wrt "'//trim(list_of_bodies(nctr))//'"'

                do i=1,10
                    call eph405%get_state( jd, ntarg, nctr, rv, status_ok_405)
                    write(*,'(F15.2,1X,*(E25.16,1X))') jd, norm2(rv(1:3)), rv
                    jd = jd + 10.0_wp
                end do

            end do
        end if

    else
        write(*,*) 'Error opening DE405 ephemeris file'
    end if

    if (status_ok_405 .and. status_ok_421) then

        !compare DE405 with DE421
        do j=1,2

            if (j==1) then
                ntarg = 3      !earth
                nctr  = 11     !sun
            else
                ntarg = 10     !moon
                nctr  = 3      !earth
            end if

            write(*,*) ''
            write(*,*) 'DE421 - DE405 Difference'
            write(*,*) 'state of "'//trim(list_of_bodies(ntarg))//&
                        '" wrt "'//trim(list_of_bodies(nctr))//'"'

            do i=1,10
                call eph405%get_state( jd, ntarg, nctr, rv1, status_ok_405)
                call eph421%get_state( jd, ntarg, nctr, rv2, status_ok_421)
                diffrv = rv2 - rv1
                write(*,'(F15.2,1X,*(E25.16,1X))') jd, norm2(diffrv(1:3)), norm2(diffrv(4:6))
                jd = jd + 10.0_wp
            end do

        end do

    else
        write(*,*) 'Error opening DE421 ephemeris file'
    end if

    !cleanup:
    call eph405%close()
    call eph421%close()

    end subroutine ephemeris_test
!*****************************************************************************************

!*****************************************************************************************
    end module jpl_ephemeris_module
!*****************************************************************************************
