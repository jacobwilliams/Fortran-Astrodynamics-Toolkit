!*****************************************************************************************
!> author: Jacob Williams
!
!  C interfaces to some of the routines.
!  This is an experiment to be able to call them from Python.
!  (see the `python_test.py` file in `tests`)

    module c_interface_module

    use iso_c_binding

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
    function initialize_geopotential_model(itype,gravfile,n,m) &
        result(cp) bind(c,name='initialize_geopotential_model')

    use geopotential_module

    implicit none

    integer(c_int),intent(in),value :: itype
    type(c_ptr),intent(in),value :: gravfile
    integer(c_int),intent(in),value :: n
    integer(c_int),intent(in),value :: m
    type(c_ptr) :: cp   !! pointer to a [[container]] containing a [[geopotential_model]]

    type(container),pointer :: grav_container
    class(geopotential_model),pointer :: grav
    logical :: status_ok
    character(len=:),allocatable :: gravfile_f
    integer :: ilen

    allocate(grav_container)

    select case (itype)
    case(1)
        allocate(geopotential_model_mueller :: grav_container%data)
        select type (g => grav_container%data)
        class is (geopotential_model_mueller)

            ! get the gravity file name:
            ilen = strlen(gravfile) ! get string length
            block
                !convert the C string to a Fortran string
                character(kind=c_char,len=ilen+1),pointer :: s
                call c_f_pointer(cptr=gravfile,fptr=s)
                gravfile_f = s(1:ilen)
                nullify(s)
            end block

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

    end function initialize_geopotential_model
!*****************************************************************************************

!*****************************************************************************************
    subroutine destroy_geopotential_model(cp) bind(c,name='destroy_geopotential_model')

    use geopotential_module

    implicit none

    type(c_ptr),intent(in),value :: cp  !! pointer to a [[container]] containing a [[geopotential_model]]

    type(container),pointer :: grav_container

    ! convert cp to fortran:
    call c_f_pointer(cp,grav_container)

    if (associated(grav_container)) then
        select type (g => grav_container%data)
        class is (geopotential_model)
            call g%destroy()
            !cp = c_null_ptr  ! should we do this too (make inout ?)
        end select
    else
        error stop 'error: pointer is not associated'
    end if

    end subroutine destroy_geopotential_model
!*****************************************************************************************

!*****************************************************************************************
    subroutine get_acceleration(cp,n,m,rvec,acc) bind(c,name='get_acceleration')

    use geopotential_module

    implicit none

    type(c_ptr),intent(in),value :: cp  !! pointer to a [[container]] containing a [[geopotential_model]]
    integer(c_int),intent(in),value :: n
    integer(c_int),intent(in),value :: m
    real(c_double),dimension(3),intent(in) :: rvec
    real(c_double),dimension(3),intent(out) :: acc

    type(container),pointer :: grav_container

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
    end module c_interface_module
!*****************************************************************************************
