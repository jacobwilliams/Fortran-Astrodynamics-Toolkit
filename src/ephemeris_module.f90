!*****************************************************************************************
!>
!  Abstract class for celestial body ephemerides.

    module ephemeris_module

    use kind_module
    use celestial_body_module

    implicit none

    private

    type,abstract,public :: ephemeris_class
        !! abstract class for interfacing with ephemeris systems.
        private
    contains
        private
        procedure(rv_func),deferred,public :: get_rv  !! get the state of one body w.r.t. another body.
    end type ephemeris_class

    abstract interface
        subroutine rv_func(me,et,targ,obs,rv,status_ok)
        !! function to return the state of the `targ` body relative to
        !! the `obs` body, in the inertial frame [ICRF].
        import :: wp,ephemeris_class,celestial_body
        implicit none
        class(ephemeris_class),intent(inout) :: me
        real(wp),intent(in)                  :: et         !! ephemeris time [sec]
        type(celestial_body),intent(in)      :: targ       !! target body
        type(celestial_body),intent(in)      :: obs        !! observer body
        real(wp),dimension(6),intent(out)    :: rv         !! state of targ w.r.t. obs
        logical,intent(out)                  :: status_ok  !! true if there were no problems
        end subroutine rv_func
    end interface

!*****************************************************************************************
    end module ephemeris_module
!*****************************************************************************************
