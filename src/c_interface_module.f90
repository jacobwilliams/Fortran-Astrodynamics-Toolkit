!*****************************************************************************************
!> author: Jacob Williams
!
!  C interfaces to some of the routines.
!  This is an experiment to be able to call them from Python.
!  (see the `python_test.py` file in `tests`)

    module c_interface_module

    use iso_c_binding
    use geopotential_module

    implicit none

    private

    type :: container
        !! a container for data that is
        !! to be passed to C. We include
        !! it here so that we can use `c_loc()`
        private
        class(*),pointer :: data
    end type container

    interface
        function strlen(str) result(isize) bind(C, name='strlen')
            !! C string length
            import
            type(c_ptr),value :: str
            integer(c_int) :: isize
        end function strlen
    end interface

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!
!@note This is just a wapper for `initialize` in [[geopotential_model]].

    function initialize_geopotential_model(itype,gravfile,n,m) &
        result(cp) bind(c,name='initialize_geopotential_model')

    implicit none

    integer(c_int),intent(in),value :: itype  !! mode :
                                              !!
                                              !! * 1 (Mueller) is only mode
                                              !! currently supported
    type(c_ptr),intent(in),value :: gravfile !! gravity coefficient file name
    integer(c_int),intent(in),value :: n    !! degree
    integer(c_int),intent(in),value :: m    !! order
    type(c_ptr) :: cp   !! pointer to a [[container]]
                        !! containing a [[geopotential_model]]

    type(container),pointer :: grav_container  !! Fortran version of `cp`
    class(geopotential_model),pointer :: grav  !! the data in the container
    logical :: status_ok !! initialization status flag
    character(len=:),allocatable :: gravfile_f  !! Fortran version of `gravfile`

    allocate(grav_container)

    select case (itype)
    case(1) !! mueller method
        allocate(geopotential_model_mueller :: grav_container%data)
        select type (g => grav_container%data)
        class is (geopotential_model_mueller)

            ! get the gravity file name:
            call c_ptr_to_f_string(gravfile,gravfile_f)

            call g%initialize(gravfile_f,n,m,status_ok)
            if (.not. status_ok) then
                write(*,*) 'error in initialize!'
                call g%destroy()
                cp = c_null_ptr
            else
                cp = c_loc(grav_container)
            end if

        end select

    case default
        error stop 'error: invalid itype input'
    end select

    ! if there was an error:
    if (c_associated(cp,c_null_ptr)) then
        deallocate(grav_container)
    end if

    end function initialize_geopotential_model
!*****************************************************************************************

!*****************************************************************************************
!>
!
!@note This is just a wapper for `destroy` in [[geopotential_model]].

    subroutine destroy_geopotential_model(cp) bind(c,name='destroy_geopotential_model')

    implicit none

    type(c_ptr),intent(in),value :: cp  !! pointer to a [[container]]
                                        !! containing a [[geopotential_model]]

    type(container),pointer :: grav_container !! Fortran version of `cp`

    ! convert cp to fortran:
    call c_f_pointer(cp,grav_container)

    if (associated(grav_container)) then
        select type (g => grav_container%data)
        class is (geopotential_model)
            call g%destroy()
            !cp = c_null_ptr  ! should we do this too (make inout ?)
        end select
        deallocate(grav_container)
    else
        error stop 'error: pointer is not associated'
    end if

    end subroutine destroy_geopotential_model
!*****************************************************************************************

!*****************************************************************************************
!>
!
!@note This is just a wapper for `get_acc` in [[geopotential_model]].

    subroutine get_acceleration(cp,n,m,rvec,acc) bind(c,name='get_acceleration')

    implicit none

    type(c_ptr),intent(in),value :: cp  !! pointer to a [[container]]
                                        !! containing a [[geopotential_model]]
    integer(c_int),intent(in),value :: n !! degree
    integer(c_int),intent(in),value :: m !! order
    real(c_double),dimension(3),intent(in) :: rvec !! position vector
    real(c_double),dimension(3),intent(out) :: acc !! acceleration vector

    type(container),pointer :: grav_container  !! Fortran version of `cp`

    ! convert cp to fortran:
    call c_f_pointer(cp,grav_container)

    if (associated(grav_container)) then
        select type (g => grav_container%data)
        class is (geopotential_model)
            call g%get_acc(rvec,n,m,acc)
        end select
    else
        error stop 'error: pointer is not associated'
    end if

    end subroutine get_acceleration
!*****************************************************************************************

!*****************************************************************************************
!>
!  Convert a `c_ptr` to a string into a Fortran string.

    subroutine c_ptr_to_f_string(cp,fstr)

    implicit none

    type(c_ptr),intent(in) :: cp
    character(len=:),allocatable,intent(out) :: fstr

    integer :: ilen  !! string length

    ilen = strlen(cp)
    block
        !convert the C string to a Fortran string
        character(kind=c_char,len=ilen+1),pointer :: s
        call c_f_pointer(cp,s)
        fstr = s(1:ilen)
        nullify(s)
    end block

    end subroutine c_ptr_to_f_string
!*****************************************************************************************

!*****************************************************************************************
!>
!  Convert a Fortran string to a `c_ptr` to a string.
!  (the C string must already have been allocated to a fixed size)

    subroutine f_string_to_c_ptr(fstr,buffer)

    implicit none

    character(len=*),intent(in) :: fstr
    type(c_ptr),intent(inout) :: buffer   !! a preallocated string buffer

    integer :: ilen !! string length of buffer

    ilen = strlen(buffer)

    block
        character(kind=c_char,len=ilen+1),pointer :: s
        call c_f_pointer(buffer,s)
        s(1:min(len(fstr),ilen)) = fstr(1:min(len(fstr),ilen))
        buffer = c_loc(s)
    end block

    end subroutine f_string_to_c_ptr
!*****************************************************************************************

!*****************************************************************************************
!>
!  Just a test of [[f_string_to_c_ptr]].

    subroutine return_a_string(ival,buffer) bind(c,name='return_a_string')

    implicit none

    integer(c_int),intent(in),value :: ival
    type(c_ptr),intent(inout)       :: buffer  !! a preallocated string buffer

    call f_string_to_c_ptr(repeat('*',ival),buffer)

    end subroutine return_a_string
!*****************************************************************************************

!*****************************************************************************************
    end module c_interface_module
!*****************************************************************************************
