!*****************************************************************************************
    module geopotential_module
!*****************************************************************************************
!****h* FAT/geopotential_module
!
!  NAME
!    geopotential_module
!
!  DESCRIPTION
!    Gravity models for computing gravitational acceleration due to geopotential.
!    The routines return the acceleration in the body-fixed frame.
!
!  EXAMPLE
!    type(geopotential_model_mueller),target :: g
!    call g%initialize(gravfile,n,m,status_ok)   
!    call g%get_acc(rvec,n,m,acc)
!
!  NOTES
!    Need to update to make sure they all work when N /= M
!
!  AUTHOR
!    Jacob Williams
!
!*****************************************************************************************
    use kind_module,    only: wp
    use numbers_module
   
    implicit none
    
    private
    
    !
    ! The base abstract class for the various geopotential models
    !
    type,abstract,public :: geopotential_model
    
        character(len=:),allocatable :: name      !model name
        character(len=:),allocatable :: filename  !model file name
        
        integer  :: nmax = 0     ! degree of the model
        integer  :: mmax = 0     ! order of the model
        real(wp) :: re   = zero  ! body radius [km]
        real(wp) :: mu   = zero  ! body grav. parameter [km3/s2]
        
        contains
        
        procedure,public :: initialize  => read_geopotential_file
        procedure,public :: destroy     => destroy_geopotential_model
        
        procedure(acc_function),deferred,public :: get_acc 
        
    end type geopotential_model
    
    !
    ! The models where the C,S coefficients are stored in vectors
    !
    type,extends(geopotential_model),abstract,public :: geopotential_model_vector_coeff
        real(wp),dimension(:),allocatable :: c 
        real(wp),dimension(:),allocatable :: s    
    end type geopotential_model_vector_coeff
        
    !
    ! The models where the C,S coefficients are stored in matrices
    !
    type,extends(geopotential_model),abstract,public :: geopotential_model_matrix_coeff
        real(wp),dimension(:,:),allocatable :: cnm
        real(wp),dimension(:,:),allocatable :: snm
    end type geopotential_model_matrix_coeff
        
    !
    !  Mueller method
    !
    type,extends(geopotential_model_vector_coeff),public :: geopotential_model_mueller
    contains
        procedure,public :: get_acc => compute_gravity_acceleration_mueller
    end type geopotential_model_mueller

    !
    !  Lear method
    !
    type,extends(geopotential_model_matrix_coeff),public :: geopotential_model_lear
    contains
        procedure,public :: get_acc => compute_gravity_acceleration_lear
    end type geopotential_model_lear
    
    !
    !  Pines method
    !
    type,extends(geopotential_model_matrix_coeff),public :: geopotential_model_pines
    contains
        procedure,public :: get_acc => compute_gravity_acceleration_pines
    end type geopotential_model_pines
    
    !
    ! Interface to the acceleration function for the different methods:
    !
    abstract interface
        subroutine acc_function(me,r,n,m,a)
            import
            implicit none
            class(geopotential_model),intent(inout)  :: me
            real(wp),dimension(3),intent(in)    :: r
            integer,intent(in)                  :: n
            integer,intent(in)                  :: m
            real(wp),dimension(3),intent(out)   :: a
        end subroutine acc_function
    end interface
    
    public :: geopotential_module_test   !for testing
    
    contains
!*****************************************************************************************

!*****************************************************************************************
!****f* geopotential_module/destroy_geopotential_model
!
!  NAME
!    destroy_geopotential_model
!
!  DESCRIPTION
!    Destroy a gravity model.
!
!  SOURCE

    subroutine destroy_geopotential_model(me)
    
    implicit none
    
    class(geopotential_model),intent(inout)  :: me
    
    !common to all:
    if (allocated(me%name))     deallocate(me%name)
    if (allocated(me%filename)) deallocate(me%filename)
    me%nmax = 0
    me%mmax = 0
    me%re   = zero
    me%mu   = zero    
   
    select type (me)
    
    class is (geopotential_model_vector_coeff)
    
        if (allocated(me%c)) deallocate(me%c)
        if (allocated(me%s)) deallocate(me%s)

    class is (geopotential_model_matrix_coeff)
        
        if (allocated(me%cnm)) deallocate(me%cnm)
        if (allocated(me%snm)) deallocate(me%snm)

    end select
    
    end subroutine destroy_geopotential_model
!*****************************************************************************************

!*****************************************************************************************
!****f* geopotential_module/compute_gravity_acceleration_mueller
!
!  NAME
!    compute_gravity_acceleration_mueller
!
!  DESCRIPTION
!    Wrapper for Mueller method.
!
!  SOURCE

    subroutine compute_gravity_acceleration_mueller(me,r,n,m,a)
    
    implicit none
    
    class(geopotential_model_mueller),intent(inout) :: me
    real(wp),dimension(3),intent(in)           :: r
    integer,intent(in)                         :: n
    integer,intent(in)                         :: m
    real(wp),dimension(3),intent(out)          :: a
    
    call geopot(r(1),r(2),r(3),n,m+1,me%re,me%mu,me%c,me%s,a(1),a(2),a(3))    
    
    end subroutine compute_gravity_acceleration_mueller
!*****************************************************************************************

!*****************************************************************************************
!****f* geopotential_module/compute_gravity_acceleration_pines
!
!  NAME
!    compute_gravity_acceleration_pines
!
!  DESCRIPTION
!    Wrapper for Pines method.
!
!  SOURCE

    subroutine compute_gravity_acceleration_pines(me,r,n,m,a)
    
    implicit none
    
    class(geopotential_model_pines),intent(inout) :: me
    real(wp),dimension(3),intent(in)         :: r
    integer,intent(in)                       :: n
    integer,intent(in)                       :: m
    real(wp),dimension(3),intent(out)        :: a
                
    call gravpot(r,n,me%re,me%mu,me%cnm,me%snm,a)
     
    end subroutine compute_gravity_acceleration_pines
!*****************************************************************************************

!*****************************************************************************************
!****f* geopotential_module/compute_gravity_acceleration_lear
!
!  NAME
!    compute_gravity_acceleration_lear
!
!  DESCRIPTION
!    Wrapper for Lear method.
!
!  SOURCE

    subroutine compute_gravity_acceleration_lear(me,r,n,m,a)
    
    implicit none
    
    class(geopotential_model_lear),intent(inout) :: me
    real(wp),dimension(3),intent(in)        :: r
    integer,intent(in)                      :: n
    integer,intent(in)                      :: m
    real(wp),dimension(3),intent(out)       :: a
                
    call grav(me%mu,r,me%re,n,m,me%cnm,me%snm,a)
     
    end subroutine compute_gravity_acceleration_lear
!*****************************************************************************************

!*****************************************************************************************
!****f* geopotential_module/read_geopotential_file
!
!  NAME
!    read_geopotential_file
!
!  DESCRIPTION
!    Read the gravity coefficient file.
!    Example file: ftp://ftp.csr.utexas.edu/pub/grav/EGM96.GEO.Z
!
!  AUTHOR
!    Jacob Williams : 9/20/2014
!
!  SOURCE

    subroutine read_geopotential_file(me,filename,nmax,mmax,status_ok)

    implicit none
    
    class(geopotential_model),intent(inout) :: me
    character(len=*),intent(in)        :: filename
    integer,intent(in)                 :: nmax
    integer,intent(in)                 :: mmax    
    logical,intent(out)                :: status_ok    

    character(len=100) :: c1,c2,fmt    ! for reading the format statements
    integer  :: iunit,istat,i1,i2,nc,i,ii,n,m
    real(wp) :: d1,d2,d3,d4,f1,t1,t2
    
    call me%destroy()
    
    if (nmax>0 .and. mmax>0 .and. mmax<=nmax) then
    
        status_ok = .true.       
        nc = number_of_coefficients(nmax,mmax)    
        me%nmax = nmax
        me%mmax = mmax
        me%filename = trim(filename)
    
        open(newunit=iunit,file=me%filename,status='OLD',iostat=istat)
    
        if (istat==0) then
    
            !size the coefficient arrays:
            select type (me)
            class is (geopotential_model_vector_coeff)
            
                allocate(me%c(nc))    !they are stored compressed into arrays
                allocate(me%s(nc))    !
                me%c = zero
                me%s = zero
                
            class is (geopotential_model_matrix_coeff)
                            
                !Note: will get all the nmax x nmax coefficients, even if mmax < nmax
                
                allocate(me%cnm(nmax,0:nmax))    !probably could replace with 2:nmax
                allocate(me%snm(nmax,0:nmax))
                me%cnm = zero
                me%snm = zero
                
            class default
            
                write(*,*) 'ERROR: INVALID geopotential_model CLASS!'
                status_ok = .false.
                return
                
            end select
                
            !Read the file:
            
            ! (2A10,2E20.10)                                               J2-DOT = -26x10-12
            !EGM 96                  398600.44150E+09          6378136.30
            ! (A6,2I3,2D19.12,2D13.6,F4.0)                                                   
            !RECOEF  2  0-0.484165371736E-03 0.000000000000E+00 0.356106E-10 0.000000E+00 -1.
            !HONKCR  2  0-0.484169544736E-03 0.000000000000E+00 0.356106E-10 0.000000E+00 -1.
            !IERS    2  1-0.186987640000E-09 0.119528010000E-08 0.100000E-29 0.100000E-29 -1.
            !RECOEF  2  2 0.243914352398E-05-0.140016683654E-05 0.537392E-10 0.543533E-10 -1.
            !RECOEF  3  0 0.957254173792E-06 0.000000000000E+00 0.180942E-10 0.000000E+00 -1.
            !RECOEF  3  1 0.202998882184E-05 0.248513158716E-06 0.139652E-09 0.136459E-09 -1.
        
            read(iunit,'(A)',iostat=istat) c1    !this is the FMT statement for the next line
            call get_format_statement(c1,fmt)    !

            read(iunit,trim(fmt),iostat=istat) c1,c2,d1,d2
                        
            me%name = trim(c1)//trim(c2)
            me%mu = d1/1000.0_wp**3    !km3/s2
            me%re = d2/1000.0_wp       !km
        
            read(iunit,'(A)',iostat=istat) c1    !this is the FMT statement for the next lines
            call get_format_statement(c1,fmt)    !
                        
            do i=1,nc    !...until end of file or all the coefficients have been read...
                
                read(iunit,trim(fmt),iostat=istat) c1,i1,i2,d1,d2,d3,d4,f1
                    
                if (istat > 0)  then
                          
                    write(*,*) 'Error reading file:' //trim(filename)
                    call me%destroy()
                    status_ok = .false.
                    exit
                    
                else if (istat < 0) then
                
                    ! end of file:
                    if (i>nc) then
                        write(*,*) 'Error: not enough coefficients in file.'
                        call me%destroy()
                        status_ok = .false.
                        exit
                    end if
                  
                end if
                
                select type (me)
                class is (geopotential_model_vector_coeff)            
                    me%c(i) = d1
                    me%s(i) = d2
                class is (geopotential_model_matrix_coeff)
                    me%cnm(i1,i2) = d1
                    me%snm(i1,i2) = d2
                end select
                
               end do
           
            close(iunit,iostat=istat)
        
        else
            write(*,*) 'Error reading file: '//trim(filename)
            call me%destroy()
            status_ok = .false.
        end if
        
    else
        write(*,*) 'Error: invalid n,m values: ',nmax,mmax
        call me%destroy()
        status_ok = .false.
    end if
    
    !unnormalize the coefficients:
    if (status_ok) then    
    
        select type (me)
        class is (geopotential_model_vector_coeff)
        
            !for this one, the coefficients are stored in arrays
            
            ii = 0    !counter    
            do n = 2, me%nmax    ! based on Lear's CONVERT routine
                t1 = 2*n+1
                ii=ii+1    
                me%c(ii) = sqrt(t1)*me%c(ii)
                me%s(ii) = zero
                do m = 1, n
                    ii=ii+1    
                    t2 = sqrt(FL(n-m)*t1*two / FL(n+m))
                    me%c(ii) = t2*me%c(ii)
                    me%s(ii) = t2*me%s(ii)
                end do
            end do
                 
        class is (geopotential_model_matrix_coeff)
        
            !for this one, the coefficients are stored in matrices
            
            call convert(me%nmax,me%cnm,me%snm)
                        
        end select
        
    end if
            
    end subroutine read_geopotential_file
!*****************************************************************************************
    
!*****************************************************************************************
!****f* geopotential_module/get_format_statement
!
!  NAME
!    get_format_statement
!
!  DESCRIPTION
!    Returns the format statement from a line 
!     in a .GEO gravity coefficient file.
!
!  AUTHOR
!    Jacob Williams : 1/24/2015
!
!  SOURCE

    subroutine get_format_statement(str,fmt)
    
    implicit none
    
    character(len=*),intent(in)  :: str
    character(len=*),intent(out) :: fmt
    
    integer :: i1,i2
    
    !note: other text after the format statement is ignored
    
    i1 = index(str,'(')
    i2 = index(str,')')
    
    if (i1/=0 .and. i2/=0 .and. i2>i1) then
        fmt = str(i1:i2)
    else
        write(*,*) 'ERROR: THE STRING DOES NOT CONTAIN A FORMAT STATEMENT: '//trim(str)
        fmt = ''
    end if
    
    end subroutine get_format_statement
!*****************************************************************************************
   
!*****************************************************************************************
!****f* geopotential_module/number_of_coefficients
!
!  NAME
!    number_of_coefficients
!
!  DESCRIPTION
!    Number of (c,s) coefficients for n x m geopotential model
!    Starting with n=2,m=0.
!
!  AUTHOR
!    Jacob Williams : 9/20/2014
!
!  SOURCE

    pure function number_of_coefficients(n,m) result(np)
    
    implicit none
    
    integer            :: np  !number of coefficients
    integer,intent(in) :: n   !degree
    integer,intent(in) :: m   !order
    
    integer :: i  !counter
    
    if (n>=m .and. n>1) then        
        np = m - 1 + sum( [ (i, i=n,2,-1) ] )
    else
        np = -999    !error
    end if
    
    end function number_of_coefficients
!*****************************************************************************************

!*****************************************************************************************
!****f* geopotential_module/FL
!
!  NAME
!    FL
!
!  DESCRIPTION
!    The FL factorial function from [1].
!
!  SEE ALSO
!    [1] W. M. Lear, "The Programs TRAJ1 and TRAJ2", 
!        JSC Mission Planning and Analysis Division,
!        JSC-22512, 87-FM-4, April 1987
!
!  HISTORY
!    Jacob Williams : 9/20/2014 : coded from [1] with some modifications.
!
!  SOURCE

    function FL(n)
    
    implicit none
    
    real(wp) :: FL
    integer,intent(in) :: n
    
    integer :: i
    
    FL = one
    if (n==0 .or. n==1) return
    
    do i=2,n
        FL = FL*i
    end do
    
    end function FL
!*****************************************************************************************

!*****************************************************************************************
!****f* geopotential_module/convert
!
!  NAME
!    convert
!
!  DESCRIPTION
!    Based on the CONVERT subroutine from [1].
!    Unnormalizes the C,S coefficients.
!
!  SEE ALSO
!    [1] W. M. Lear, "The Programs TRAJ1 and TRAJ2", 
!        JSC Mission Planning and Analysis Division,
!        JSC-22512, 87-FM-4, April 1987
!
!  HISTORY
!    Jacob Williams : 9/20/2014 : coded from [1]. some modifications.
!
!  SOURCE

    subroutine convert(nmodel,cnm,snm)
    
    implicit none
    
    integer,intent(in) :: nmodel
    real(wp),dimension(nmodel,0:nmodel),intent(inout) :: cnm
    real(wp),dimension(nmodel,0:nmodel),intent(inout) :: snm

    integer :: n,m
    real(wp) :: t1,t2
    
    do n = 1, nmodel        !JW : this could be 2,nmodel ...
        t1 = 2*n+1
        cnm(n,0) = sqrt(t1)*cnm(n,0)
        do m = 1, n
            t2 = sqrt(FL(n-m)*t1*two / FL(n+m))
            cnm(n,m) = t2*cnm(n,m)
            snm(n,m) = t2*snm(n,m)
        end do
    end do

    end subroutine convert
!*****************************************************************************************

!*****************************************************************************************
!****f* geopotential_module/gravpot
!
!  NAME
!    gravpot
!
!  DESCRIPTION
!    Spencer's implementation of the Pines algorithms from [1]
!
!  SEE ALSO
!    [1] J.L. Spencer, "Pines' nonsingular gravitational potential
!        derivation, description, and implementation", 
!        NASA-CR-147478, MDC-W0013, Feb 9, 1976.
!        http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19760011100.pdf
!
!  HISTORY
!    Jacob Williams, 1/25/2014 : updated and fixed bugs in the original code.
!
!  SOURCE

    subroutine gravpot(r,nmax,re,mu,c,s,fg)
        
    implicit none
    
    real(wp),dimension(3),intent(in)           :: r     !position vector
    integer,intent(in)                         :: nmax  !degree/order
    real(wp),intent(in)                        :: re    !body radius
    real(wp),intent(in)                        :: mu    !grav constant
    real(wp),dimension(nmax,0:nmax),intent(in) :: c     !coefficients
    real(wp),dimension(nmax,0:nmax),intent(in) :: s     !    
    real(wp),dimension(3),intent(out)          :: fg    !grav acceleration
    
    !local variables:
    real(wp),dimension(nmax+1) :: creal, cimag, rho
    real(wp),dimension(nmax+1,nmax+1) :: a,d,e,f
    integer :: nax0,i,j,k,l,n,m
    real(wp) :: rinv,ess,t,u,r0,rhozero,a1,a2,a3,a4,fac1,fac2,fac3,fac4
    real(wp) :: ci1,si1,temp
    
    !JW : not done in original paper,
    !     but seems to be necessary
    !     (probably assumed the compiler 
    !      did it automatically)
    a = zero
    d = zero
    e = zero
    f = zero
    
    !get the direction cosines ess, t and u:
    
    nax0 = nmax + 1
    rinv = one/norm2(r)
    ess  = r(1) * rinv
    t    = r(2) * rinv
    u    = r(3) * rinv
    
    !generate the functions creal, cimag, a, d, e, f and rho:
    
    r0       = re*rinv 
    rhozero  = mu*rinv    !JW: typo in original paper (see p.18)
    rho(1)   = r0*rhozero    
    creal(1) = ess
    cimag(1) = t
    d(1,1)   = zero
    e(1,1)   = zero
    f(1,1)   = zero
    a(1,1)   = one
        
    main_loop: do i=2,nax0
    
        if (i/=nax0) then !JW : to prevent access of c,s outside bounds
            ci1 = c(i,1)
            si1 = s(i,1)
        else
            ci1 = zero
            si1 = zero
        end if
        
        rho(i)   = r0*rho(i-1)
        creal(i) = ess*creal(i-1) - t*cimag(i-1)
        cimag(i) = ess*cimag(i-1) + t*creal(i-1)
        d(i,1)   = ess*ci1 + t*si1
        e(i,1)   = ci1
        f(i,1)   = si1                        
        a(i,i)   = (2*i-1)*a(i-1,i-1)
        a(i,i-1) = u*a(i,i)
        
        do k=2,i
                    
            if (i/=nax0) then        
                d(i,k) = c(i,k)*creal(k)   + s(i,k)*cimag(k)
                e(i,k) = c(i,k)*creal(k-1) + s(i,k)*cimag(k-1)
                f(i,k) = s(i,k)*creal(k-1) - c(i,k)*cimag(k-1)    
            end if
 
            !JW : typo in original paper 
            !  (should be GOTO 1, rather than GOTO 10)
            if (i/=2) then     
                L = i-2
                do j=1,L            
                    a(i,i-j-1) = (u*a(i,i-j)-a(i-1,i-j))/(j+1)                
                end do
            end if
            
        end do
        
    end do main_loop
  
    !compute auxiliary quantities a1, a2, a3, a4
    
    a1 = zero
    a2 = zero
    a3 = zero
    a4 = rhozero*rinv
         
    do n=2,nmax
    
        fac1 = zero
        fac2 = zero
        fac3 = a(n,1)  *c(n,0)
        fac4 = a(n+1,1)*c(n,0)
        
        do m=1,n
            temp = m*a(n,m)
            fac1 = fac1 + temp      *e(n,m)
            fac2 = fac2 + temp      *f(n,m)
            fac3 = fac3 + a(n,m+1)  *d(n,m)    
            fac4 = fac4 + a(n+1,m+1)*d(n,m)
        end do
        
        temp = rinv*rho(n)
        a1 = a1 + temp*fac1
        a2 = a2 + temp*fac2
        a3 = a3 + temp*fac3
        a4 = a4 + temp*fac4
        
    end do
    
    fg(1) = a1 - ess*a4
    fg(2) = a2 - t*a4
    fg(3) = a3 - u*a4
          
    end subroutine gravpot
!*****************************************************************************************
     
!*****************************************************************************************
!****f* geopotential_module/geopot
!
!  NAME
!    geopot
!
!  DESCRIPTION
!    Compute the gravitational acceleration vector using the Mueller method.
!
!    !!!!! WARNING: this one crashes if nmax/=mmax !!!!!
!
!  SEE ALSO
!    [1] Alan C. Mueller, "A Fast Recursive Algorithm for Calculating the 
!            Forces due to the Geopotential (Program GEOPOT)", 
!            JSC Internal Note 75-FM-42, June 9, 1975.
!
!  SOURCE
    
    subroutine geopot(x,y,z,nmax,mmax,re,ksq,c,s,fx,fy,fz)
        
    implicit none
    
    !subroutine arguments:
    real(wp),intent(in)              :: x,y,z        ! position vector
    integer,intent(in)               :: nmax         ! degree of model
    integer,intent(in)               :: mmax         ! order+1 of model
    real(wp),intent(in)              :: re           ! body radius
    real(wp),intent(in)              :: ksq          ! body GM
    real(wp),dimension(:),intent(in) :: c            ! C,S coefficients
    real(wp),dimension(:),intent(in) :: s            !
    real(wp),intent(out)             :: fx,fy,fz     ! gravitational acceleration
    
    !local variables:
    real(wp) :: r,ri,reor,reorn,ksqor2,xor,yor,zor,rdedx,rdedy,rdedz,&
                sum1,sum2,sum3,sum4,temp1,temp2,temp3,temp4,fact,dcstld,temp
    integer :: i,j,k,im1,l,jm1,jp1,kk
    real(wp),dimension(0:mmax) :: p0,p1,p2,ctil,stil
    
    !write(*,'(A,1x,*(e20.5,1x/))') 'c=',c
    !write(*,'(A,1x,*(e20.5,1x/))') 's=',s
    
    !abbreviations:
    
    r      = sqrt(x*x + y*y + z*z)
    ri     = one/r
    reor   = re*ri
    reorn  = reor
    ksqor2 = ksq*ri*ri
    zor    = z*ri
    xor    = x*ri
    yor    = y*ri
    
    !the derivatives of the argument of the legendre polynomial - zor
    
    rdedz = zor*zor - one
    rdedx = zor*xor
    rdedy = zor*yor
    
    !initialization:
    
    k = 0    
    do i=1,mmax
        p0(i) = zero
        p1(i) = zero
    end do
    p0(0)   = one
    p1(0)   = zor
    p1(1)   = one
    ctil(0) = one
    stil(0) = zero
    ctil(1) = xor
    stil(1) = yor
    sum1    = zero
    !sum2   = zero       !original
    sum2    = one        !JW : include central body term
    sum3    = zero
    sum4    = zero
    
    !computation of forces:
    
    do i = 2,nmax
        
        reorn = reorn*reor
        fact  = 2*i - 1
        im1   = i-1
        l     = 1
        
        !recursion formulas for legendre polynomial - p2(0)
        
        p2(0) = (fact*zor*p1(0)-im1*p0(0))/i
        k     = k + 1
        p2(1) = p0(1)+fact*p1(0)
        temp1 = p2(1)*c(k)
        temp2 = p2(0)*c(k)*(i+1)
        
        if (i < mmax) then
            
            !recursive formulas for: 
            !    'ctilda' - ctil
            !    'stilda' - stil
            
            ctil(i) = ctil(1)*ctil(im1) - stil(1)*stil(im1)
            stil(i) = stil(1)*ctil(im1) + ctil(1)*stil(im1)
            temp3   = zero
            temp4   = zero
            
            do j=1,i
                
                jm1 = j-1
                jp1 = j+1
                
                !recursive formula for derivative of legendre polynomial - p2(j)
                
                p2(jp1) = p0(jp1) + fact*p1(j)
                kk      = k + j
                dcstld  = j*p2(j)
                temp    = (c(kk)*ctil(j)+s(kk)*stil(j))
                temp1   = temp1+p2(jp1)*temp
                temp2   = temp2+(i+jp1)*p2(j)*temp
                temp3   = temp3+dcstld*(c(kk)*ctil(jm1)+s(kk)*stil(jm1))
                temp4   = temp4-dcstld*(c(kk)*stil(jm1)-s(kk)*ctil(jm1))
                
            end do
            
            l = i
            sum3 = sum3+reorn*temp3
            sum4 = sum4+reorn*temp4
        
        end if
        
        sum1 = sum1+reorn*temp1
        sum2 = sum2+reorn*temp2
        k = k + i
        
        !shift indices:
        
        do j = 0, l
            p0(j) = p1(j)
            p1(j) = p2(j)    
        end do
        
    end do
    
    fx = -ksqor2*(sum1*rdedx + sum2*xor - sum3 )
    fy = -ksqor2*(sum1*rdedy + sum2*yor - sum4 )
    fz = -ksqor2*(sum1*rdedz + sum2*zor        )
                
    end subroutine geopot
!*****************************************************************************************
    
!*****************************************************************************************
!****f* geopotential_module/grav
!
!  NAME
!    grav
!
!  DESCRIPTION
!    Based on the GRAV subroutine from [1].
!
!  SEE ALSO
!    [1] W. M. Lear, "The Programs TRAJ1 and TRAJ2", 
!        JSC Internal Note 87-FM-4, April 1987.
!    [2] W. M. Lear, "The Gravitational Acceleration Equations", 
!        JSC Internal Note 86-FM-15, April 1986.
!
!  AUTHOR
!    Jacob Williams : 9/20/2014 : updated
!
!  SOURCE

    subroutine grav(mu,rgr,rbar,nmodel,mmodel,cnm,snm,agr)
    
    implicit none
    
    real(wp),dimension(3),intent(in)                :: rgr      ! position vector [body-fixed coordinates]
    real(wp),intent(in)                             :: mu       ! gravitational constant
    real(wp),intent(in)                             :: rbar     ! gravitational scaling radius (generally the equatorial radius)
    integer,intent(in)                              :: nmodel   ! the degree of the gravity model (>=2)
    integer,intent(in)                              :: mmodel   ! the order of the gravity model (>=0, <=nmodel)
    real(wp),dimension(nmodel,0:nmodel),intent(in)  :: cnm      ! gravity coefficients
    real(wp),dimension(nmodel,0:nmodel),intent(in)  :: snm      !
    real(wp),dimension(3),intent(out)               :: agr      ! gravitational acceleration vector [body-fixed coordinates]
    
    !local variables:
    real(wp),dimension(nmodel,nmodel) :: pnm,ppnm
    real(wp),dimension(nmodel) :: cm,sm,pn,rb,ppn
    real(wp),dimension(3) :: asph
     real(wp) :: e1,e2,e3,e4,e5,r1,r2,t1,t3,absr,sphi,cphi,tcm,tsm,tsnm,tcnm,tpnm
    integer :: n,nm1,nm2,m
        
    do n=2,nmodel
        pnm(n-1,n) = zero
    end do
        
    e1 = rgr(1)**2 + rgr(2)**2
    r2 = e1 + rgr(3)**2
    absr = sqrt(r2)
    r1 = sqrt(e1)
    sphi = rgr(3)/absr
    cphi = r1/absr
    if (r1==zero) then
        sm(1) = zero
        cm(1) = one    
    else
        sm(1) = rgr(2)/r1
        cm(1) = rgr(1)/r1
    end if
    rb(1) = rbar/absr
    rb(2) = rb(1)**2
    sm(2) = two*cm(1)*sm(1)
    cm(2) = two*cm(1)**2 - one
    pn(1) = sphi
    pn(2) = (three*sphi**2 - one)/two
    ppn(1) = one
    ppn(2) = three*sphi
    pnm(1,1) = one
    pnm(2,2) = three*cphi
    pnm(2,1) = ppn(2)
    ppnm(1,1) = -sphi
    ppnm(2,2) = -six*sphi*cphi
    ppnm(2,1) = three-six*sphi**2
    
    if (nmodel>=3) then
    
        do n = 3, nmodel    
            nm1 = n-1
            nm2 = n-2
            rb(n) = rb(nm1)*rb(1)
            sm(n) = two*cm(1)*sm(nm1)-sm(nm2)
            cm(n) = two*cm(1)*cm(nm1)-cm(nm2)
            e1 = 2*n-1
            pn(n) = (e1*sphi*pn(nm1)-nm1*pn(nm2))/n
            ppn(n) = sphi*ppn(nm1)+n*pn(nm1)
            pnm(n,n) = e1*cphi*pnm(nm1,nm1)
            ppnm(n,n) = -n*sphi*pnm(n,n)
        end do
        
        do n = 3, nmodel
            nm1 = n-1
            e1 = (2*n-1)*sphi
            e2 = -n*sphi
            do m = 1, nm1
                e3 = pnm(nm1,m)
                e4 = n+m
                e5 = (e1*e3-(e4-one)*pnm(n-2,m))/(n-m)
                pnm(n,m) = e5
                ppnm(n,m) = e2*e5 + e4*e3
            end do
        end do
    
    end if
    
    asph(1) = -one        ![NOTE: set to zero to only output the harmonic terms]
    asph(3) = zero
    
    do n = 2,nmodel
        e1 = cnm(n,0)*rb(n)
        asph(1) = asph(1)-(n+1)*e1*pn(n)
        asph(3) = asph(3) + e1*ppn(n)
    end do
    asph(3) = cphi*asph(3)
    t1 = zero
    t3 = zero
    asph(2) = zero
    
    do n = 2, nmodel
        e1 = zero
        e2 = zero
        e3 = zero
        !do m = 1, n                !original
        do m = 1, min(n,mmodel)     !JW - allow for specifying order !!!!!!
            tsnm = snm(n,m)
            tcnm = cnm(n,m)
            tsm = sm(m)
            tcm = cm(m)
            tpnm = pnm(n,m)
            e4 = tsnm*tsm+tcnm*tcm
            e1 = e1+e4*tpnm
            e2 = e2+m*(tsnm*tcm-tcnm*tsm)*tpnm
            e3 = e3+e4*ppnm(n,m)
        end do
        t1 = t1 + (n+1)*rb(n)*e1
        asph(2) = asph(2) + rb(n)*e2
        t3 = t3+rb(n)*e3
    end do
    
    e4 = mu/r2
    asph(1) = e4*(asph(1) - cphi*t1)
    asph(2) = e4*asph(2)
    asph(3) = e4*(asph(3) + t3)
    
    e5 = asph(1)*cphi - asph(3)*sphi
    
    agr(1) = e5*cm(1) - asph(2)*sm(1)
    agr(2) = e5*sm(1) + asph(2)*cm(1)
    agr(3) = asph(1)*sphi + asph(3)*cphi
    
    end subroutine grav
!*****************************************************************************************

!*****************************************************************************************
!****f* geopotential_module/geopotential_module_test
!
!  NAME
!    geopotential_module_test
!
!  DESCRIPTION
!    Unit test routine for geopotential_module
!
!  AUTHOR
!    Jacob Williams : 9/20/2014
!
!  SOURCE

    subroutine geopotential_module_test()
    
    use conversion_module
    use vector_module,  only: spherical_to_cartesian
    use random_module,  only: get_random_number
    
    implicit none

    !the coefficient file:    
    character(len=*),parameter :: gravfile = '../grav/GGM03C.GEO'
    
    class(geopotential_model),pointer :: g
    type(geopotential_model_mueller),target :: g_mueller
    type(geopotential_model_lear)   ,target :: g_lear
    type(geopotential_model_pines)  ,target :: g_pines
    
    real(wp),dimension(3) :: a1,a2,a3,rvec
    logical :: status_ok
    integer :: lat,lon,i,j,nmax,mmax
    real(wp) :: h,err1,err2,rlon,rlat,tmp
    character(len=20) :: name
    real :: tstart, tstop
    
    !test case:
    real(wp),dimension(3),parameter :: r = [0.1275627320e+05_wp, &
                                            0.1275627320e+05_wp, &
                                            0.1275627320e+05_wp ]   !km
    
    integer,parameter :: n=5    !degree
    integer,parameter :: m=5    !order
    
    integer,parameter :: n_repeat = 1000000  !number of time to repeat speed test
       
    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' geopotential_module_test'
    write(*,*) '---------------'
    write(*,*) ''

    write(*,*) ''
    write(*,*) '*** number_of_coefficients routine ***'
    write(*,*) ''
    
    !tests:
    write(*,*) ''
    write(*,*) 6,6,number_of_coefficients(6,6)    ! 25
    write(*,*) 3,2,number_of_coefficients(3,2)    ! 6
    write(*,*) 4,0,number_of_coefficients(4,0)    ! 8
    write(*,*) 4,1,number_of_coefficients(4,1)    ! 9
    write(*,*) 0,0,number_of_coefficients(0,0)    ! -99 (error)
    write(*,*) ''
    
    write(*,*) ''
    write(*,*) '*** test case for the three methods ***'
    write(*,*) ''
    
    write(*,*) ''
    write(*,*) 'reading file: '//gravfile
    write(*,'(A,*(E30.16,1X))') 'r =',r

    call g_mueller%initialize(gravfile,n,m,status_ok)   
    if (.not. status_ok) stop 'Error'
    call g_mueller%get_acc(r,n,m,a1)

    call g_lear%initialize(gravfile,n,m,status_ok)   
    if (.not. status_ok) stop 'Error'
    call g_lear%get_acc(r,n,m,a2)

    call g_pines%initialize(gravfile,n,m,status_ok)   
    if (.not. status_ok) stop 'Error'
    call g_pines%get_acc(r,n,m,a3)
    
    write(*,*) ''
    write(*,'(A,*(E30.16,1X))') 'mueller =',a1
    write(*,'(A,*(E30.16,1X))') 'lear    =',a2
    write(*,'(A,*(E30.16,1X))') 'pines   =',a3      
    write(*,*) ''
    write(*,'(A,*(E30.16,1X))') 'mueller-lear difference  =',norm2(a1-a2)
    write(*,'(A,*(E30.16,1X))') 'mueller-pines difference =',norm2(a1-a3)
    write(*,'(A,*(E30.16,1X))') 'lear-pines difference    =',norm2(a2-a3)
    write(*,*) ''
        
    write(*,*) ''
    write(*,*) '*** accuracy test ***'
    write(*,*) ''
    
    h = 6778.0_wp    !radius magnitude
    status_ok = .true.
    do lat = -90, 90
        do lon = 0, 360
    
            rvec = spherical_to_cartesian(h,lon*deg2rad,lat*deg2rad)
                       
            call g_mueller%get_acc(rvec,n,m,a1)     ! mueller
            call g_lear%get_acc(rvec,n,m,a2)        ! lear
            call g_pines%get_acc(rvec,n,m,a3)       ! pines
 
            err1 = norm2( a2 - a1 )
            err2 = norm2( a3 - a2 )
            if (err2>1.0e-15_wp .or. err1>1.0e-15_wp) then
                write(*,*) lat,lon,norm2(a1),norm2(a2),norm2(a2),err1,err2
                status_ok = .false.
            end if
            if (abs(lat)==90) exit    !only do poles once
        
        end do
    end do
    call g_mueller%destroy()
    call g_lear%destroy()
    call g_pines%destroy()
    if (status_ok) write(*,*) 'All tests passed.'
     
    write(*,*) ''
    write(*,*) '*** speed test ***'
    write(*,*) ''
    
    nmax = 10
    mmax = 10
    
    do i=1,3
    
        select case(i)
        case(1)
            g => g_mueller
            name = 'Mueller'
        case(2)
            g => g_lear
            name = 'Lear'
        case(3)
            g => g_pines
            name = 'Pines'
        end select
        call g%initialize(gravfile,nmax,mmax,status_ok)
        if (.not. status_ok) stop 'Error'
       
        call random_seed()
        tmp = zero
        call cpu_time(tstart)
        do j=1,n_repeat
        
             h    = get_random_number(6778.0_wp, 10000.0_wp)
             rlon = get_random_number(0.0_wp, 360.0_wp)
             rlat = get_random_number(-90.0_wp, 90.0_wp)
             
             rvec = spherical_to_cartesian(h,rlon*deg2rad,rlat*deg2rad)        
             call g%get_acc(rvec,nmax,mmax,a1)
             tmp = tmp + norm2(a1)
    
        end do   
        call cpu_time(tstop)
       
        call g%destroy()
        
        write(*,'(A10,1X,E30.16,1X,F13.6,1X,A)') trim(name), tmp, tstop-tstart, 'sec'
    
    end do
    
    end subroutine geopotential_module_test
!*****************************************************************************************

!*****************************************************************************************
    end module geopotential_module
!*****************************************************************************************