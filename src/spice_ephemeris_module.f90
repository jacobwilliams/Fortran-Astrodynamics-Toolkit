!*****************************************************************************************
!>
!  Interface to the SPICE ephemeris library.
!
!  Not a standard part of FAT. If used, it requires linking with the Fortran
!  [SPICELIB SPICE Toolkit](https://naif.jpl.nasa.gov/naif/toolkit.html).
!
!@note Haven't validated this yet.

    module spice_ephemeris_module

    use, intrinsic :: iso_fortran_env, only: error_unit
    use ephemeris_module
    use kind_module

    implicit none

    private

    type,extends(ephemeris_class),public :: spice_ephemeris

        !! Main class for accessing the SPICE ephemeris system.
        !!
        !! Note: SPICE is not object-oriented or threadsafe. So,
        !! while this class provides an object-oriented like interface
        !! to SPICE, it should really be treated as a singleton.

        character(len=:),dimension(:),allocatable :: kernels  !! the list of kernels

    contains

        procedure,public :: get_rv     => get_rv_from_spice_ephemeris
        procedure,public :: initialize => initialize_spice_ephemeris
        procedure,public :: close      => close_spice_ephemeris

    end type spice_ephemeris

    !these routines are in the SPICELIB:
    interface
        subroutine trcoff()
            implicit none
        end subroutine trcoff
        function failed()
            !! see: ftp://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/failed.html
            implicit none
            logical :: failed
        end function failed
        subroutine furnsh(file)
            !! see: ftp://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/furnsh.html
            implicit none
            character(len=*),intent(in) :: file
        end subroutine furnsh
        subroutine unload(file)
            !! see: ftp://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/unload.html
            implicit none
            character(len=*),intent(in) :: file
        end subroutine unload
        subroutine kclear()
            !! see: ftp://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/kclear.html
            implicit none
        end subroutine kclear
        subroutine spkgeo ( targ, et, ref, obs, state, lt )
            !! see: ftp://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/spkgeo.html
            import :: wp
            implicit none
            integer               :: targ
            real(wp)              :: et
            character(len=*)      :: ref
            integer               :: obs
            real(wp),dimension(6) :: state
            real(wp)              :: lt
        end subroutine spkgeo
    end interface

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Close the SPICE ephemeris and unload all the kernels.

    subroutine close_spice_ephemeris(me)

    implicit none

    class(spice_ephemeris),intent(inout) :: me

    integer :: i  !! counter

    !unload all the kernels:
    if (allocated(me%kernels)) then
        do i=1,size(me%kernels)
            call unload(trim(me%kernels(i)))
        end do
        deallocate(me%kernels)
    end if

    !clear the system:
    call kclear()

    end subroutine close_spice_ephemeris
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize a SPICE ephemeris by loading the specified kernels.

    subroutine initialize_spice_ephemeris(me,kernels)

    implicit none

    class(spice_ephemeris),intent(inout) :: me
    character(len=*),dimension(:),intent(in) :: kernels  !! list of kernels to load

    integer :: i !! counter

    ! disable the SPICE traceback system to speed it up.
    call trcoff()

    call me%close() ! just in case

    !save the kernel list for unloading later:
    allocate( character(len=len(kernels)) :: me%kernels(size(kernels)))

    !load all the kernels:
    do i=1,size(kernels)
        call furnsh(trim(kernels(i)))
    end do

    end subroutine initialize_spice_ephemeris
!*****************************************************************************************

!*****************************************************************************************
!>
!  Interface for the [[ephemeris_module]].
!
!  Return the Cartesian state of `targ` relative to `obs` in the `J2000` frame.

    subroutine get_rv_from_spice_ephemeris(me,et,targ,obs,rv,status_ok)

    use celestial_body_module, only: celestial_body
    use numbers_module,        only: zero

    implicit none

    class(spice_ephemeris),intent(inout) :: me
    real(wp),intent(in)                :: et         !! ephemeris time [sec]
    type(celestial_body),intent(in)    :: targ       !! target body
    type(celestial_body),intent(in)    :: obs        !! observer body
    real(wp),dimension(6),intent(out)  :: rv         !! state of targ w.r.t. obs [km,km/s] in ICRF frame
    logical,intent(out)                :: status_ok  !! true if there were no problems

    real(wp) :: lt  !! light time output from spkgeo

    call spkgeo ( targ%id, et, 'J2000', obs%id, rv, lt )

    status_ok = .not. failed()

    end subroutine get_rv_from_spice_ephemeris
!*****************************************************************************************

!*****************************************************************************************
    end module spice_ephemeris_module
!*****************************************************************************************
