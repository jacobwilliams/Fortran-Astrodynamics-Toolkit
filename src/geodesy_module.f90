!*****************************************************************************************    
    module geodesy_module
!*****************************************************************************************
!****h* FAT/geodesy_module
!
!  NAME
!    geodesy_module
!
!  DESCRIPTION
!    Geodesy routines.
!
!*****************************************************************************************    
    
    use kind_module,    only: wp
    
    implicit none
    
    private
    
    public :: heikkinen
    public :: olson
    public :: direct
    public :: geodetic_to_cartesian
    public :: great_circle_distance
    
    contains
!*****************************************************************************************    
     
!*****************************************************************************************
!****f* geodesy_module/heikkinen
!
!  NAME
!    heikkinen
!
!  DESCRIPTION
!    Heikkinen routine for cartesian to geodetic transformation
!
!  SEE ALSO
!    [1] M. Heikkinen, "Geschlossene formeln zur berechnung raumlicher 
!        geodatischer koordinaten aus rechtwinkligen Koordinaten". 
!        Z. Ermess., 107 (1982), 207-211 (in German).
!    [2] E. D. Kaplan, "Understanding GPS: Principles and Applications", 
!        Artech House, 1996.
!
!  AUTHOR
!    Jacob Williams
!
!  SOURCE

    pure subroutine heikkinen(rvec, a, b, h, lon, lat)

    implicit none
    
    real(wp),dimension(3),intent(in) :: rvec     !position vector [km]
    real(wp),intent(in)  :: a                    !geoid semimajor axis [km]
    real(wp),intent(in)  :: b                    !geoid semiminor axis [km]
    real(wp),intent(out) :: h                    !geodetic altitude [km]
    real(wp),intent(out) :: lon                  !longitude [rad]
    real(wp),intent(out) :: lat                  !geodetic latitude [rad]

    real(wp) :: f,e_2,ep,r,e2,ff,g,c,s,pp,q,r0,u,v,z0,x,y,z,z2,r2,tmp,a2,b2
            
    x   = rvec(1)
    y   = rvec(2)
    z   = rvec(3)       
    a2  = a*a
    b2  = b*b
    f   = (a-b)/a
    e_2 = (2.0_wp*f-f*f)
    ep  = sqrt(a2/b2 - 1.0_wp)
    z2  = z*z
    r   = sqrt(x**2 + y**2)
    r2  = r*r
    e2  = a2 - b2
    ff  = 54.0_wp * b2 * z2
    g   = r2 + (1.0_wp - e_2)*z2 - e_2*e2
    c   = e_2**2 * ff * r2 / g**3
    s   = (1.0_wp + c + sqrt(c**2 + 2.0_wp*c))**(1.0_wp/3.0_wp)
    pp  = ff / ( 3.0_wp*(s + 1.0_wp/s + 1.0_wp)**2 * g**2 )
    q   = sqrt( 1.0_wp + 2.0_wp*e_2**2 * pp )
    r0  = -pp*e_2*r/(1.0_wp+q) + &
            sqrt( max(0.0_wp, 1.0_wp/2.0_wp * a2 * (1.0_wp + 1.0_wp/q) - &
                ( pp*(1.0_wp-e_2)*z2 )/(q*(1.0_wp+q)) - &
                1.0_wp/2.0_wp * pp * r2) )
    u   = sqrt( (r - e_2*r0)**2 + z2 )
    v   = sqrt( (r - e_2*r0)**2 + (1.0_wp - e_2)*z2 )
    z0  = b**2 * z / (a*v)
    
    h   = u*(1.0_wp - b2/(a*v) )
    lat = atan2( (z + ep**2*z0), r )
    lon = atan2( y, x )
        
    end subroutine heikkinen
!*****************************************************************************************

!*****************************************************************************************
!****f* geodesy_module/olson
!
!  NAME
!    olson
!
!  DESCRIPTION
!    Olson routine for cartesian to geodetic transformation.
!
!  SEE ALSO
!    [1] Olson, D. K., Converting Earth-Centered, Earth-Fixed Coordinates to
!       Geodetic Coordinates, IEEE Transactions on Aerospace and Electronic
!       Systems, 32 (1996) 473-476.
!
!  AUTHOR
!    Jacob Williams
!
!  SOURCE

    pure subroutine olson(rvec, a, b, h, long, lat)
    
    implicit none
    
    real(wp),dimension(3),intent(in) :: rvec     !position vector [km]
    real(wp),intent(in)  :: a                    !geoid semimajor axis [km]
    real(wp),intent(in)  :: b                    !geoid semiminor axis [km]
    real(wp),intent(out) :: h                    !geodetic altitude [km]
    real(wp),intent(out) :: long                 !longitude [rad]
    real(wp),intent(out) :: lat                  !geodetic latitude [rad]
    
    real(wp) :: f,x,y,z,e2,a1,a2,a3,a4,a5,a6,w,zp,&
                w2,r2,r,s2,c2,u,v,s,ss,c,g,rg,rf,m,p,z2
    
    x  = rvec(1)                
    y  = rvec(2)                
    z  = rvec(3)                
    f  = (a-b)/a
    e2 = f * (2.0_wp - f)
    a1 = a * e2
    a2 = a1 * a1
    a3 = a1 * e2 / 2.0_wp
    a4 = 2.5_wp * a2
    a5 = a1 + a3
    a6 = 1.0_wp - e2 
    zp = abs(z)
    w2 = x*x + y*y
    w  = sqrt(w2)
    z2 = z * z
    r2 = z2 + w2
    r  = sqrt(r2)
    
    if (r < 100.0_wp) then
    
        lat = 0.0_wp
        long = 0.0_wp
        h = -1.0e7_wp
    
    else
    
        s2 = z2 / r2
        c2 = w2 / r2
        u  = a2 / r
        v  = a3 - a4 / r
    
        if (c2 > 0.3_wp) then
            s = (zp / r) * (1.0_wp + c2 * (a1 + u + s2 * v) / r)
            lat = asin(s)
            ss = s * s
            c = sqrt(1.0_wp - ss)
        else
            c = (w / r) * (1.0_wp - s2 * (a5 - u - c2 * v) / r)
            lat = acos(c)
            ss = 1.0_wp - c * c
            s = sqrt(ss)
        end if
    
        g   = 1.0_wp - e2 * ss
        rg  = a / sqrt(g)
        rf  = a6 * rg
        u   = w - rg * c
        v   = zp - rf * s
        f   = c * u + s * v
        m   = c * v - s * u
        p   = m / (rf / g + f)
        lat = lat + p
        if (z < 0.0_wp) lat = -lat
        h = f + m * p / 2.0_wp
        long = atan2( y, x )   
    
    end if
            
    end subroutine olson    
!*****************************************************************************************
 
!*****************************************************************************************
!****f* geodesy_module/direct
!
!  NAME
!    direct
!
!  DESCRIPTION
!    Solve the "direct" geodetic problem: given the latitude and longitude of one 
!    point and the azimuth and distance to a second point, determine the latitude 
!    and longitude of that second point.  The solution is obtained using the 
!    algorithm by Vincenty.
!
!  SEE ALSO
!    [1] T. Vincenty, "Direct and Inverse Solutions of Geodesics on the 
!        Ellipsoid with Application of Nested Equations", 
!        Survey Review XXII. 176, April 1975.
!        http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf
!    [2] http://www.ngs.noaa.gov/PC_PROD/Inv_Fwd/
!
!  AUTHOR
!    Jacob Williams
!
!  SOURCE

    subroutine direct(a,f,glat1,glon1,glat2,glon2,faz,baz,s)

    implicit none

    real(wp),intent(in)  :: a        !semimajor axis of ellipsoid [m]
    real(wp),intent(in)  :: f        !flattening of ellipsoid [-]            
    real(wp),intent(in)  :: glat1    !latitude of 1 [rad]
    real(wp),intent(in)  :: glon1    !longitude of 1 [rad]
    real(wp),intent(in)  :: faz      !forward azimuth 1->2 [rad]
    real(wp),intent(in)  :: s        !distance from 1->2 [m]
    real(wp),intent(out) :: glat2    !latitude of 2 [rad]
    real(wp),intent(out) :: glon2    !longitude of 2 [rad]
    real(wp),intent(out) :: baz      !back azimuth 2->1 [rad]

    real(wp) :: r,tu,sf,cf,cu,su,sa,csa,c2a,x,c,d,y,sy,cy,cz,e

    real(wp),parameter :: pi  = acos(-1.0_wp)
    real(wp),parameter :: eps = 0.5e-13_wp

    r = 1.0_wp-f
    tu = r*sin(glat1)/cos(glat1)
    sf = sin(faz)
    cf = cos(faz)
    baz = 0.0_wp
    if (cf/=0.0_wp) baz = atan2(tu,cf)*2.0_wp
    cu = 1.0_wp/sqrt(tu*tu+1.0_wp)
    su = tu*cu
    sa = cu*sf
    c2a = -sa*sa+1.0_wp
    x = sqrt((1.0_wp/r/r-1.0_wp)*c2a+1.0_wp)+1.0_wp
    x = (x-2.0_wp)/x
    c = 1.0_wp-x
    c = (x*x/4.0_wp+1.0_wp)/c
    d = (0.375_wp*x*x-1.0_wp)*x
    tu = s/r/a/c
    y = tu
    do
        sy = sin(y)
        cy = cos(y)
        cz = cos(baz+y)
        e = cz*cz*2.0_wp-1.0_wp
        c = y
        x = e*cy
        y = e+e-1.0_wp
        y = (((sy*sy*4.0_wp-3.0_wp)*y*cz*d/6.0_wp+x)*d/4.0_wp-cz)*sy*d+tu
        if (abs(y-c)<=eps) exit
    end do
    baz = cu*cy*cf-su*sy
    c = r*sqrt(sa*sa+baz*baz)
    d = su*cy+cu*sy*cf
    glat2 = atan2(d,c)
    c = cu*cy-su*sy*cf
    x = atan2(sy*sf,c)
    c = ((-3.0_wp*c2a+4.0_wp)*f+4.0_wp)*c2a*f/16.0_wp
    d = ((e*cy*c+cz)*sy*c+y)*sa
    glon2 = glon1+x-(1.0_wp-c)*d*f
    baz = atan2(sa,baz)+pi

    end subroutine direct
!*****************************************************************************************

!*****************************************************************************************
!****f* geodesy_module/geodetic_to_cartesian
!
!  NAME
!    geodetic_to_cartesian
!
!  DESCRIPTION
!    Geodetic latitude, longitude, and height to Cartesian position vector.
!
!  SEE ALSO
!    [1] E. D. Kaplan, "Understanding GPS: Principles and Applications", 
!        Artech House, 1996.
!
!  AUTHOR
!    Jacob Williams
!
!  SOURCE

    subroutine geodetic_to_cartesian(a,b,glat,lon,h,r)
    
    implicit none
    
    real(wp),intent(in) :: a                !geoid semimajor axis [km]
    real(wp),intent(in) :: b                !geoid semiminor axis [km]
    real(wp),intent(in) :: glat             !geodetic latitude [rad]
    real(wp),intent(in) :: lon              !longitude [rad]
    real(wp),intent(in) :: h                !geodetic altitude [km]
    real(wp),dimension(3),intent(out) :: r  !Cartesian position vector [x,y,z]

    real(wp) :: e2,slat,clat,slon,clon,tlat,ome2,d,q,aod

    slat    = sin(glat)
    clat    = cos(glat)
    tlat    = tan(glat)
    slon    = sin(lon)
    clon    = cos(lon)
    e2      = 1.0_wp - (b*b)/(a*a)
    ome2    = 1.0_wp - e2
    d       = sqrt( 1.0_wp + ome2*tlat*tlat )
    q       = sqrt( 1.0_wp - e2*slat*slat   )
    aod     = a/d
    
    r(1) = aod*clon + h*clon*clat
    r(2) = aod*slon + h*slon*clat
    r(3) = a*ome2*slat/q + h*slat
        
    end subroutine geodetic_to_cartesian
!*****************************************************************************************

!*****************************************************************************************    
!****f* geodesy_module/great_circle_distance
!
!  NAME
!    great_circle_distance
!
!  DESCRIPTION
!    Great circle distance on a spherical body, using the Vincenty algorithm.
!
!  SEE ALSO
!    [1] T. Vincenty, "Direct and Inverse Solutions of Geodesics on the 
!        Ellipsoid with Application of Nested Equations", 
!        Survey Review XXII. 176, April 1975.
!        http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf
!
!  AUTHOR
!    Jacob Williams, 7/13/2014
!
!  SOURCE

    function great_circle_distance(r,long1,lat1,long2,lat2) result(d)

    implicit none

    real(wp)            :: d        ! great circle distance from 1 to 2 [km]
    real(wp),intent(in) :: r        ! radius of the body [km]
    real(wp),intent(in) :: long1    ! longitude of first site [rad]
    real(wp),intent(in) :: lat1     ! latitude of the first site [rad]
    real(wp),intent(in) :: long2    ! longitude of the second site [rad]
    real(wp),intent(in) :: lat2     ! latitude of the second site [rad]

    real(wp) :: c1,s1,c2,s2,dlon,clon,slon

    !Compute aux variables:
    c1    = cos(lat1)
    s1    = sin(lat1)
    c2    = cos(lat2)
    s2    = sin(lat2)
    dlon  = long1-long2
    clon  = cos(dlon)
    slon  = sin(dlon)
        
    d = r*atan2( sqrt((c2*slon)**2+(c1*s2-s1*c2*clon)**2), (s1*s2+c1*c2*clon) )
                     
    end function great_circle_distance
!*****************************************************************************************

!*****************************************************************************************    
    end module geodesy_module
!*****************************************************************************************    