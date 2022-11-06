!*****************************************************************************************
!> author: Jacob Williams
!
!  Geodesy routines.

    module geodesy_module

    use kind_module,    only: wp
    use numbers_module

    implicit none

    private

    public :: heikkinen
    public :: olson
    public :: direct
    public :: inverse
    public :: cartesian_to_geodetic_triaxial
    public :: CartesianIntoGeodeticI, CartesianIntoGeodeticII
    public :: cartesian_to_geodetic_triaxial_2

    public :: geodetic_to_cartesian
    public :: geodetic_to_cartesian_triaxial
    public :: geodetic_to_cartesian_triaxial_2

    public :: great_circle_distance
    public :: geocentric_radius

    ! unit tests:
    public :: direct_inverse_test

    contains
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!
!  Heikkinen routine for cartesian to geodetic transformation
!
!# References
!  1. M. Heikkinen, "Geschlossene formeln zur berechnung raumlicher
!     geodatischer koordinaten aus rechtwinkligen Koordinaten".
!     Z. Ermess., 107 (1982), 207-211 (in German).
!  2. E. D. Kaplan, "Understanding GPS: Principles and Applications",
!     Artech House, 1996.

    pure subroutine heikkinen(rvec, a, b, h, lon, lat)

    implicit none

    real(wp),dimension(3),intent(in) :: rvec  !! position vector [km]
    real(wp),intent(in)  :: a                 !! geoid semimajor axis [km]
    real(wp),intent(in)  :: b                 !! geoid semiminor axis [km]
    real(wp),intent(out) :: h                 !! geodetic altitude [km]
    real(wp),intent(out) :: lon               !! longitude [rad]
    real(wp),intent(out) :: lat               !! geodetic latitude [rad]

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
!> author: Jacob Williams
!
!  Olson routine for cartesian to geodetic transformation.
!
!# References
!  1. Olson, D. K., Converting Earth-Centered, Earth-Fixed Coordinates to
!     Geodetic Coordinates, IEEE Transactions on Aerospace and Electronic
!     Systems, 32 (1996) 473-476.

    pure subroutine olson(rvec, a, b, h, long, lat)

    implicit none

    real(wp),dimension(3),intent(in) :: rvec !!position vector [km]
    real(wp),intent(in)  :: a                !!geoid semimajor axis [km]
    real(wp),intent(in)  :: b                !!geoid semiminor axis [km]
    real(wp),intent(out) :: h                !!geodetic altitude [km]
    real(wp),intent(out) :: long             !!longitude [rad]
    real(wp),intent(out) :: lat              !!geodetic latitude [rad]

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
!> author: Jacob Williams
!
!  Solve the "direct" geodetic problem: given the latitude and longitude of one
!  point and the azimuth and distance to a second point, determine the latitude
!  and longitude of that second point.  The solution is obtained using the
!  algorithm by Vincenty.
!
!# References
!  1. T. Vincenty, "[Direct and Inverse Solutions of Geodesics on the
!     Ellipsoid with Application of Nested Equations](http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf)",
!     Survey Review XXII. 176, April 1975.
!  2. [PC Software Download - INVERSE and FORWARD](http://www.ngs.noaa.gov/PC_PROD/Inv_Fwd/),
!     National Geodetic Survey. Version 3.0 (November, 2012).

    subroutine direct(a,f,glat1,glon1,faz,s,glat2,glon2,baz)

    implicit none

    real(wp),intent(in)  :: a        !! semimajor axis of ellipsoid [m]
    real(wp),intent(in)  :: f        !! flattening of ellipsoid [-]
    real(wp),intent(in)  :: glat1    !! latitude of 1 [rad]
    real(wp),intent(in)  :: glon1    !! longitude of 1 [rad]
    real(wp),intent(in)  :: faz      !! forward azimuth 1->2 [rad]
    real(wp),intent(in)  :: s        !! distance from 1->2 [m]
    real(wp),intent(out) :: glat2    !! latitude of 2 [rad]
    real(wp),intent(out) :: glon2    !! longitude of 2 [rad]
    real(wp),intent(out) :: baz      !! back azimuth 2->1 [rad]

    real(wp) :: r,tu,sf,cf,cu,su,sa,csa,c2a,x,c,d,y,sy,cy,cz,e

    real(wp),parameter :: pi  = acos(-1.0_wp)
    real(wp),parameter :: eps = 0.5e-13_wp

    r  = 1.0_wp - f
    tu = r*sin(glat1)/cos(glat1)
    sf = sin(faz)
    cf = cos(faz)
    if ( cf/=0.0_wp ) then
        baz = atan2(tu,cf)*2.0_wp
    else
        baz = 0.0_wp
    end if
    cu  = 1.0_wp/sqrt(tu*tu+1.0_wp)
    su  = tu*cu
    sa  = cu*sf
    c2a = -sa*sa + 1.0_wp
    x   = sqrt((1.0_wp/r/r-1.0_wp)*c2a+1.0_wp) + 1.0_wp
    x   = (x-2.0_wp)/x
    c   = 1.0_wp - x
    c   = (x*x/4.0_wp+1.0_wp)/c
    d   = (0.375_wp*x*x-1.0_wp)*x
    tu  = s/r/a/c
    y   = tu
    do
        sy = sin(y)
        cy = cos(y)
        cz = cos(baz+y)
        e  = cz*cz*2.0_wp - 1.0_wp
        c  = y
        x  = e*cy
        y  = e + e - 1.0_wp
        y  = (((sy*sy*4.0_wp-3.0_wp)*y*cz*d/6.0_wp+x)*d/4.0_wp-cz)*sy*d + tu
        if ( abs(y-c)<=eps ) exit
    end do
    baz   = cu*cy*cf - su*sy
    c     = r*sqrt(sa*sa+baz*baz)
    d     = su*cy + cu*sy*cf
    glat2 = atan2(d,c)
    c     = cu*cy - su*sy*cf
    x     = atan2(sy*sf,c)
    c     = ((-3.0_wp*c2a+4.0_wp)*f+4.0_wp)*c2a*f/16.0_wp
    d     = ((e*cy*c+cz)*sy*c+y)*sa
    glon2 = glon1 + x - (1.0_wp-c)*d*f
    baz   = atan2(sa,baz) + pi

    end subroutine direct
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!
!  Geodetic latitude, longitude, and height to Cartesian position vector.
!
!# References
!  1. E. D. Kaplan, "Understanding GPS: Principles and Applications",
!     Artech House, 1996.

    subroutine geodetic_to_cartesian(a,b,glat,lon,h,r)

    implicit none

    real(wp),intent(in) :: a                !! geoid semimajor axis [km]
    real(wp),intent(in) :: b                !! geoid semiminor axis [km]
    real(wp),intent(in) :: glat             !! geodetic latitude [rad]
    real(wp),intent(in) :: lon              !! longitude [rad]
    real(wp),intent(in) :: h                !! geodetic altitude [km]
    real(wp),dimension(3),intent(out) :: r  !! Cartesian position vector [x,y,z]

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
!> author: Jacob Williams
!  date: 7/13/2014
!
!  Great circle distance on a spherical body, using the Vincenty algorithm.
!
!# References
!  * T. Vincenty, "[Direct and Inverse Solutions of Geodesics on the Ellipsoid
!    with Application of Nested Equations](http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf)",
!    Survey Review XXII. 176, April 1975.

    pure function great_circle_distance(r,long1,lat1,long2,lat2) result(d)

    implicit none

    real(wp)            :: d        !! great circle distance from 1 to 2 [km]
    real(wp),intent(in) :: r        !! radius of the body [km]
    real(wp),intent(in) :: long1    !! longitude of first site [rad]
    real(wp),intent(in) :: lat1     !! latitude of the first site [rad]
    real(wp),intent(in) :: long2    !! longitude of the second site [rad]
    real(wp),intent(in) :: lat2     !! latitude of the second site [rad]

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
!>
!  The distance from the center of a celestial body (e.g., the Earth) to a point
!  on the spheroid surface at a specified geodetic latitude.
!
!### Reference
!  * [Geocentric radius](https://en.wikipedia.org/wiki/Earth_radius#Geocentric_radius)

    pure function geocentric_radius(a,b,lat) result(r)

    use numbers_module, only: zero

    implicit none

    real(wp),intent(in) :: a    !! equatorial radius (km)
    real(wp),intent(in) :: b    !! polar radius of point (km)
    real(wp),intent(in) :: lat  !! geodetic latitude of point (rad)
    real(wp)            :: r    !! distance from center of body to point (km)

    !local variables:
    real(wp) :: num,den,cl2,sl2,a2,b2

    if (a==zero .and. b==zero) then
        r = zero
    else
        cl2 = cos(lat)**2
        sl2 = sin(lat)**2
        a2  = a*a
        b2  = b*b
        num = cl2 * a2**2 + sl2 * b2**2
        den = cl2 * a2    + sl2 * b2
        r   = sqrt(num/den)
    end if

    end function geocentric_radius
!*****************************************************************************************

!*****************************************************************************************
!>
!  INVERSE computes the geodetic azimuth and distance between two points,
!  given their geographic positions.
!
!  Version for long-line and antipodal cases.
!  Latitudes may be 90 degrees exactly.
!
!### Reference
!  * T. Vincenty, "[Direct and Inverse Solutions of Geodesics on the Ellipsoid
!    with Application of Nested Equations](http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf)",
!    Survey Review XXII. 176, April 1975.
!  * [inverse.for](http://www.ngs.noaa.gov/PC_PROD/Inv_Fwd/source/inverse.for)
!    Version 3.0 (November, 2012).
!
!### History
!  * Original programmed by thaddeus vincenty, 1975, 1976
!  * Removed back side solution option, debugged, revised -- 2011may01 -- dgm
!    this version of code is interim -- antipodal boundary needs work
!  * Jacob Williams, 1/25/2016 : refactored into modern Fortran.

    subroutine inverse(a,rf,b1,l1,b2,l2,faz,baz,s,it,sig,lam,kind)

    implicit none

    real(wp),intent(in)  :: a     !! Equatorial semimajor axis
    real(wp),intent(in)  :: rf    !! reciprocal flattening (1/f)
    real(wp),intent(in)  :: b1    !! latitude of point 1 (rad, positive north)
    real(wp),intent(in)  :: l1    !! longitude of point 1 (rad, positive east)
    real(wp),intent(in)  :: b2    !! latitude of point 2 (rad, positive north)
    real(wp),intent(in)  :: l2    !! longitude of point 2 (rad, positive east)
    real(wp),intent(out) :: faz   !! Forward azimuth (rad, clockwise from north)
    real(wp),intent(out) :: baz   !! Back azimuth (rad, clockwise from north)
    real(wp),intent(out) :: s     !! Ellipsoidal distance
    integer,intent(out)  :: it    !! iteration count
    real(wp),intent(out) :: sig   !! spherical distance on auxiliary sphere
    real(wp),intent(out) :: lam   !! longitude difference on auxiliary sphere
    integer,intent(out)  :: kind  !! solution flag: kind=1, long-line; kind=2, antipodal

    real(wp) :: beta1,beta2,biga,bigb,bige,bigf,boa,c,cosal2,coslam,&
                cossig,costm,costm2,cosu1,cosu2,d,dsig,ep2,l,prev,&
                sinal,sinlam,sinsig,sinu1,sinu2,tem1,tem2,temp,test,z

    real(wp),parameter :: pi     = acos(-1.0_wp)
    real(wp),parameter :: two_pi = 2.0_wp * pi
    real(wp),parameter :: tol    = 1.0e-14_wp   !! convergence tolerance
    real(wp),parameter :: eps    = 1.0e-15_wp   !! tolerance for zero

    boa = 1.0_wp - 1.0_wp/rf   ! b/a

    beta1 = atan(boa*tan(b1))  ! better reduced latitude
    sinu1 = sin(beta1)
    cosu1 = cos(beta1)
    beta2 = atan(boa*tan(b2))  ! better reduced latitude
    sinu2 = sin(beta2)
    cosu2 = cos(beta2)

    l = l2 - l1  ! longitude difference [-pi,pi]
    if ( l>pi ) l = l - two_pi
    if ( l<-pi ) l = l + two_pi
    prev = l
    test = l
    it   = 0
    kind = 1
    lam  = l

    longline : do  ! long-line loop (kind=1)

        sinlam = sin(lam)
        coslam = cos(lam)
        temp   = cosu1*sinu2 - sinu1*cosu2*coslam
        sinsig = sqrt((cosu2*sinlam)**2+temp**2)
        cossig = sinu1*sinu2 + cosu1*cosu2*coslam
        sig    = atan2(sinsig,cossig)
        if ( abs(sinsig)<eps ) then
            sinal = cosu1*cosu2*sinlam/sign(eps,sinsig)
        else
            sinal = cosu1*cosu2*sinlam/sinsig
        endif
        cosal2 = -sinal**2 + 1.0_wp
        if ( abs(cosal2)<eps ) then
            costm = -2.0_wp*(sinu1*sinu2/sign(eps,cosal2)) + cossig
        else
            costm = -2.0_wp*(sinu1*sinu2/cosal2) + cossig
        endif
        costm2 = costm*costm
        c = ((-3.0_wp*cosal2+4.0_wp)/rf+4.0_wp)*cosal2/rf/16.0_wp

        antipodal : do  ! antipodal loop (kind=2)

            it = it + 1
            d = (((2.0_wp*costm2-1.0_wp)*cossig*c+costm)*sinsig*c+sig)*(1.0_wp-c)/rf
            if ( kind==1 ) then
                lam = l + d*sinal
                if ( abs(lam-test)>=tol ) then
                    if ( abs(lam)>pi ) then
                        kind   = 2
                        lam    = pi
                        if ( l<0.0_wp ) lam = -lam
                        sinal  = 0.0_wp
                        cosal2 = 1.0_wp
                        test   = 2.0_wp
                        prev   = test
                        sig    = pi - abs(atan(sinu1/cosu1)+atan(sinu2/cosu2))
                        sinsig = sin(sig)
                        cossig = cos(sig)
                        c      = ((-3.0_wp*cosal2+4.0_wp)/rf+4.0_wp)*cosal2/rf/16.0_wp
                        if ( abs(sinal-prev)<tol ) exit longline
                        if ( abs(cosal2)<eps ) then
                            costm = -2.0_wp*(sinu1*sinu2/sign(eps,cosal2)) + cossig
                        else
                            costm = -2.0_wp*(sinu1*sinu2/cosal2) + cossig
                        endif
                        costm2 = costm*costm
                        cycle antipodal
                    endif
                    if ( ((lam-test)*(test-prev))<0.0_wp .and. it>5 ) &
                         lam = (2.0_wp*lam+3.0_wp*test+prev)/6.0_wp  ! refined converge.
                    prev = test
                    test = lam
                    cycle longline
                endif
            else
                sinal  = (lam-l)/d
                if ( ((sinal-test)*(test-prev))<0.0_wp .and. it>5 ) &
                       sinal = (2.0_wp*sinal+3.0_wp*test+prev)/6.0_wp  ! refined converge.
                prev   = test
                test   = sinal
                cosal2 = -sinal**2 + 1.0_wp
                sinlam = sinal*sinsig/(cosu1*cosu2)
                coslam = -sqrt(abs(-sinlam**2+1.0_wp))
                lam    = atan2(sinlam,coslam)
                temp   = cosu1*sinu2 - sinu1*cosu2*coslam
                sinsig = sqrt((cosu2*sinlam)**2+temp**2)
                cossig = sinu1*sinu2 + cosu1*cosu2*coslam
                sig    = atan2(sinsig,cossig)
                c      = ((-3.0_wp*cosal2+4.0_wp)/rf+4.0_wp)*cosal2/rf/16.0_wp
                if ( abs(sinal-prev)>=tol ) then
                    if ( abs(cosal2)<eps ) then
                        costm = -2.0_wp*(sinu1*sinu2/sign(eps,cosal2)) + cossig
                    else
                        costm = -2.0_wp*(sinu1*sinu2/cosal2) + cossig
                    endif
                    costm2 = costm*costm
                    cycle antipodal
                endif
            endif

            exit longline  !finished

        end do antipodal

    end do longline

    ! convergence

    if ( kind==2 ) then  ! antipodal
        faz  = sinal/cosu1
        baz  = sqrt(-faz**2+1.0_wp)
        if ( temp<0.0_wp ) baz = -baz
        faz  = atan2(faz,baz)
        tem1 = -sinal
        tem2 = sinu1*sinsig - cosu1*cossig*baz
        baz  = atan2(tem1,tem2)
    else  ! long-line
        tem1 = cosu2*sinlam
        tem2 = cosu1*sinu2 - sinu1*cosu2*coslam
        faz  = atan2(tem1,tem2)
        tem1 = -cosu1*sinlam
        tem2 = sinu1*cosu2 - cosu1*sinu2*coslam
        baz  = atan2(tem1,tem2)
    endif
    if ( faz<0.0_wp ) faz = faz + two_pi
    if ( baz<0.0_wp ) baz = baz + two_pi

    ! helmert 1880 from vincenty "geodetic inverse solution between antipodal points"

    ep2  = 1.0_wp/(boa*boa) - 1.0_wp
    bige = sqrt(1.0_wp+ep2*cosal2)
    bigf = (bige-1.0_wp)/(bige+1.0_wp)
    biga = (1.0_wp+bigf*bigf/4.0_wp)/(1.0_wp-bigf)
    bigb = bigf*(1.0_wp-0.375_wp*bigf*bigf)
    z    = bigb/6.0_wp*costm*(-3.0_wp+4.0_wp*sinsig**2)*(-3.0_wp+4.0_wp*costm2)
    dsig = bigb*sinsig*(costm+bigb/4.0_wp*(cossig*(-1.0_wp+2.0_wp*costm2)-z))
    s    = (boa*a)*biga*(sig-dsig)

    end subroutine inverse
!*****************************************************************************************

!*****************************************************************************************
!>
!  Function computes the Cartesian coordinates given the
!  geodetic latitude (phi), longitude (lambda) and
!  height (h) of a point related to an ellipsoid
!  defined by its three semiaxes ax, ay and b
!
!### History
!  * Jacob Williams, 10/29/2022 : Fortran verison of this algorithm,
!    based on the Matlab (v1.0 01/03/2019) code.

subroutine geodetic_to_cartesian_triaxial(ax, ay, b, phi, lambda, h, r)

    real(wp),intent(in) :: ax !! semiaxes (0 < b <= ay <= ax)
    real(wp),intent(in) :: ay !! semiaxes (0 < b <= ay <= ax)
    real(wp),intent(in) :: b  !! semiaxes (0 < b <= ay <= ax)
    real(wp),intent(in) :: phi !! geodetic latitude (radians)
    real(wp),intent(in) :: lambda !! geodetic longitude (radians)
    real(wp),intent(in) :: h !! geodetic height
    real(wp),dimension(3),intent(out) :: r  !! Cartesian position vector [x,y,z]

    real(wp) :: ee2,ex2,N,cp,sp,cl,sl

    cp  = cos(phi)
    sp  = sin(phi)
    cl  = cos(lambda)
    sl  = sin(lambda)
    ee2 = (ax*ax-ay*ay)/(ax*ax)
    ex2 = (ax*ax-b*b)/(ax*ax)
    N   = ax/sqrt(one - ex2*sp*sp - ee2*cp*cp*sl*sl)

    r = [(N+h)*cp*cl, &
         (N*(one-ee2)+h)*cp*sl, &
         (N*(one-ex2)+h)*sp ]

end subroutine geodetic_to_cartesian_triaxial

!********************************************************************************
!>
!  Geodetic to Cartesian for Triaxial Ellipsoid.
!
!### References
!  * S. Bektas, "Geodetic Computations on Triaxial Ellipsoid",
!    International Journal of Mining Science (IJMS),
!    Volume 1, Issue 1, June 2015, p 25-34

    pure subroutine geodetic_to_cartesian_triaxial_2(a,b,c,lat,long,h,r)

    implicit none

    real(wp),intent(in) :: a    !! ellipsoid radii `a >= b >= c`
    real(wp),intent(in) :: b    !! ellipsoid radii `a >= b >= c`
    real(wp),intent(in) :: c    !! ellipsoid radii `a >= b >= c`
    real(wp),intent(in) :: lat  !! latitude (rad)
    real(wp),intent(in) :: long !! longitude (rad)
    real(wp),intent(in) :: h    !! altitude
    real(wp),dimension(3),intent(out) :: r  !! Cartesian coordinates (x,y,z)

    real(wp) :: ex2,ee2,v,a2,clat,slat,clon,slon,omee2,omex2

    a2    = a * a
    ex2   = (a2-c**2)/a2
    ee2   = (a2-b**2)/a2
    clat  = cos(lat)
    slat  = sin(lat)
    clon  = cos(long)
    slon  = sin(long)
    omee2 = 1.0_wp-ee2
    omex2 = 1.0_wp-ex2
    v     = a/sqrt(1.0_wp-ex2*slat**2-ee2*clat**2*slon**2)

    r = [(v+h)*clon*clat, &
         (v*omee2+h)*slon*clat, &
         (v*omex2+h)*slat ]

    end subroutine geodetic_to_cartesian_triaxial_2
!*****************************************************************************************

!*****************************************************************************************
!>
!  Function computes the geodetic latitude (phi), longitude (lambda) and
!  height (h) of a point related to an ellipsoid
!  defined by its three semiaxes ax, ay and b (0 < b <= ay <= ax)
!  given Cartesian coordinates Xi, Yi, Zi and tolerance (tol).
!  Latitude and longitude are returned in radians.
!
!### Reference
!  * Panou G. and Korakitis R. (2019) "Cartesian to geodetic coordinates conversion
!    on an ellipsoid using the bisection method".
!    Journal of Geodesy volume 96, Article number: 66 (2022).
!    [(link)](https://link.springer.com/article/10.1007/s00190-022-01650-9)
!
!### History
!  * Jacob Williams, 10/29/2022 : Fortran verison of this algorithm,
!    based on the Matlab (v1.0 01/03/2019) code.

subroutine cartesian_to_geodetic_triaxial(ax, ay, b, r, tol, phi, lambda, h)

    real(wp),intent(in) :: ax !! semiaxes (0 < b <= ay <= ax)
    real(wp),intent(in) :: ay !! semiaxes (0 < b <= ay <= ax)
    real(wp),intent(in) :: b  !! semiaxes (0 < b <= ay <= ax)
    real(wp),dimension(3),intent(in) :: r !! Cartesian coordinates (x,y,z)
    real(wp),intent(in) :: tol !! tolerance (may be set to zero)
    real(wp),intent(out) :: phi !! geodetic latitude (radians)
    real(wp),intent(out) :: lambda !! geodetic longitude (radians)
    real(wp),intent(out) :: h !! geodetic height

    real(wp) :: kx,ky,cx,cy,cz,XX,YY,ZZ,x,y,z,Xo,Yo,Zo,m,Mm,axax,ayay,b2
    integer :: n

    if (ax < ay .or. ay < b) error stop 'error in cartesian_to_geodetic_triaxial: invalid ax,ay,b'

    axax = ax*ax
    ayay = ay*ay
    b2   = b*b
    kx   = (axax-b2)/ax
    ky   = (ayay-b2)/ay
    cx   = (axax)/(b2)
    cy   = (ayay)/(b2)
    cz   = (axax)/(ayay)

    XX = abs(r(1))
    YY = abs(r(2))
    ZZ = abs(r(3))

    ! Compute geodetic latitude/longitude
    if (ZZ == zero) then
        if (XX == zero .and. YY == zero) then
            x = zero
            y = zero
            z = b
        elseif (ky*XX*ky*XX+kx*YY*kx*YY < kx*ky*kx*ky) then
            x = ax*XX/kx
            y = ay*YY/ky
            z = b*sqrt(one-((x*x)/(axax))-((y*y)/(ayay)))
        elseif (XX == zero) then
            x = zero
            y = ay
            z = zero
        elseif (YY == zero) then
            x = ax
            y = zero
            z = zero
        else
            Xo = XX/ax
            Yo = YY/ay
            call bisection_special_2(cz, Xo, Yo, tol, n, m, Mm)
            x = cz*XX/(cz+m)
            y = YY/(one+m)
            z = zero
        end if
    else
        if (XX == zero .and. YY == zero) then
            x = zero
            y = zero
            z = b
        else
            Xo = XX/ax
            Yo = YY/ay
            Zo = ZZ/b
            call bisection_special_3(cx, cy, Xo, Yo, Zo, tol, n, m, Mm)
            x = cx*XX/(cx+m)
            y = cy*YY/(cy+m)
            if (m < zero .and. ky*XX*ky*XX + kx*YY*kx*YY < kx*ky*kx*ky) then
                z = b*sqrt(one-((x*x)/(axax))-((y*y)/(ayay)))
            else
                z = ZZ/(one+m)
            end if
        end if
    end if

    call xyz2fl(ax, ay, b, x, y, z, phi, lambda)        ! analytic method used for initial guess
    call xyz2philambda(ax, ay, b, x, y, z, phi, lambda) ! iterative method
    call philambda_quadrant(r(1), r(2), r(3), phi, lambda)

    ! Compute geodetic height
    h = norm2([XX-x, YY-y, ZZ-z])
    if ((XX+YY+ZZ) < (x+y+z)) h = -h

end subroutine cartesian_to_geodetic_triaxial

!*****************************************************************************************
!>

subroutine bisection_special_2(cz, Xo, Yo, tol, n, m, Gm)

    real(wp),intent(in) :: cz, Xo, Yo, tol
    integer,intent(out) :: n
    real(wp),intent(out) :: m, Gm

    real(wp) :: d1,Gd1,d2,d,MM

    d1 = -one+Yo
    Gd1 = (cz*Xo*cz*Xo)/((cz+d1)*(cz+d1))
    d2 = -one+sqrt(cz*Xo*cz*Xo+Yo*Yo)
    d = (d2 - d1)/two
    n = 0
    m = -two

    do while (d > tol)
        n = n + 1
        MM = m
        m = d1 + d
        Gm = ((cz*Xo*cz*Xo)/((cz+m)*(cz+m)))+((Yo*Yo)/((one+m)**2))-one
        if (MM == m + tol .or. Gm == zero) return
        if (sign(one,Gm) == sign(one,Gd1)) then
            d1 = m
            Gd1 = Gm
        else
            d2 = m
        end if
        d = (d2 - d1)/two
    end do

    n = n + 1
    m = d1 + d
    Gm = ((cz*Xo*cz*Xo)/((cz+m)*(cz+m)))+((Yo*Yo)/((one+m)**2))-one

end subroutine bisection_special_2

!*****************************************************************************************
!>

subroutine bisection_special_3(cx, cy, Xo, Yo, Zo, tol, n, m, Hm)

    real(wp),intent(in) :: cx, cy, Xo, Yo, Zo, tol
    integer,intent(out) :: n
    real(wp),intent(out) :: m, Hm

    real(wp) :: d1,Hd1,d2,d,MM

    d1 = -one+Zo
    Hd1 = ((cx*Xo*cx*Xo)/((cx+d1)*(cx+d1)))+((cy*Yo*cy*Yo)/((cy+d1)*(cy+d1)))
    d2 = -one+sqrt(cx*Xo*cx*Xo+cy*Yo*cy*Yo+Zo*Zo)
    d = (d2 - d1)/two
    n = 0
    m = -two

    do while (d > tol)
        n = n + 1
        MM = m
        m = d1 + d
        Hm = ((cx*Xo*cx*Xo)/((cx+m)*(cx+m)))+((cy*Yo*cy*Yo)/&
             ((cy+m)*(cy+m)))+((Zo*Zo)/((one+m)**2))-one
        if (MM == m + tol .or. Hm == zero) return
        if (sign(one,Hm) == sign(one,Hd1)) then
            d1 = m
            Hd1 = Hm
        else
            d2 = m
        end if
        d = (d2 - d1)/two
    end do

    n = n + 1
    m = d1 + d
    Hm = ((cx*Xo*cx*Xo)/((cx+m)*(cx+m)))+((cy*Yo*cy*Yo)/&
         ((cy+m)*(cy+m)))+((Zo*Zo)/((one+m)**2))-one

end subroutine bisection_special_3

!*****************************************************************************************
!>

subroutine philambda_quadrant(x, y, z, phi, lambda)

    real(wp),intent(in) :: x, y, z
    real(wp),intent(inout) :: phi, lambda

    if (z < zero) then
        phi = -phi
    end if

    if (x >= zero) then
        if (y >= zero) then
            lambda = lambda
        else
            lambda = -lambda
        end if
    else
        if (y >= zero) then
            lambda = pi-lambda
        else
            lambda = lambda-pi
        end if
    end if

end subroutine philambda_quadrant

!******************************************************************************
!>
!  Determination of the geodetic latitude and longitude

    subroutine xyz2philambda(ax, ay, b, x, y, z, phi, lambda)

     real(wp),intent(in) :: ax, ay, b, x, y, z
     real(wp),intent(inout) :: phi, lambda !! input: initial guess, output: refined values

     real(wp) :: ee2,ex2,Sphi,Cphi,Slambda,Clambda,&
                 Den,NN,onemee2,onemex2,dndphi,dxdphi,&
                 dydphi,dzdphi,dndlam,dxdlam,dydlam,dzdlam
     integer :: n
     real(wp),dimension(3,2) :: J
     real(wp),dimension(2,3) :: Jt !! transpose of J
     real(wp),dimension(3,1) :: dl
     real(wp),dimension(2,2) :: Nmat, Ninv
     real(wp),dimension(2,1) :: dx
     real(wp),dimension(3) :: r0

     integer,parameter :: maxiter = 100 !! maximum number of iterations

     ee2 = (ax*ax - ay*ay)/(ax*ax) ! eqn. 5
     ex2 = (ax*ax - b*b)/(ax*ax)   !
     onemee2 = one - ee2
     onemex2 = one - ex2

     do n = 1, maxiter

        ! Design Matrix J
        Sphi = sin(phi)
        Cphi = cos(phi)
        Slambda = sin(lambda)
        Clambda = cos(lambda)

        NN  = ax/sqrt(one-ex2*Sphi*Sphi-ee2*Cphi*Cphi*Slambda*Slambda) ! eqn. 4
        Den = two*(one-ex2*Sphi**2-ee2*Cphi**2*Slambda**2)**(three/two)
        dndphi = -ax*sin(two*phi)*(ex2 - ee2*Slambda**2) / Den
        dxdphi = (dndphi*Cphi - NN*Sphi) * Clambda
        dydphi = onemee2*(dndphi*Cphi - NN*Sphi) * Slambda
        dzdphi = onemex2*(dndphi*Sphi + NN*Cphi)
        dndlam = -ax*ee2*Cphi**2*sin(two*lambda) / Den
        dxdlam = (dndlam*Clambda - NN*Slambda)*Cphi
        dydlam = onemee2*(dndlam*Slambda + NN*Clambda)*Cphi
        dzdlam = onemex2*dndlam*Sphi
        J = reshape([dxdphi,dydphi,dzdphi,dxdlam,dydlam,dzdlam],[3,2])

        ! Vector dl
        call geodetic_to_cartesian_triaxial_2(ax,ay,b,phi,lambda,0.0_wp,r0) ! just use the main one with alt=0
        dl(:,1) = [x,y,z] - r0 ! eqn. 51

        ! Solution
        Jt      = transpose(J)
        Nmat    = matmul(Jt,J) ! eqn. 53
        Ninv    = (one / (Nmat(1,1)*Nmat(2,2) - Nmat(1,2)*Nmat(2,1))) * &
                  reshape([Nmat(2,2),-Nmat(2,1),-Nmat(1,2),Nmat(1,1)], [2,2]) ! eqn. 54
        dx      = matmul(Ninv, matmul(Jt,dl)) ! eqn. 52
        phi     = phi    + dx(1,1)   ! corrections. eqn. 55
        lambda  = lambda + dx(2,1)   !

        if (all(abs(dx) <= ten*epsilon(1.0_wp))) exit

     end do

    end subroutine xyz2philambda
!*****************************************************************************************

!*****************************************************************************************
!>
!  Computes the transformation of Cartesian to geodetic coordinates on the surface of the ellipsoid
!  assuming x,y,z are all non-negative
!  Angular coordinates in radians
!
!  This is based on the [C++ version](https://www.researchgate.net/publication/353739609_PK-code)

subroutine xyz2fl(ax, ay, b, x, y, z, latitude, longitude)

    real(wp),intent(in) :: ax, ay, b, x, y, z
    real(wp),intent(out) :: latitude, longitude

    real(wp) :: nom,den,dex,xme,rot
    real(wp) :: ax2,ay2,b2,Ex2,Ee2,lex2,lee2,mex,mee

    ! note: these could be precomputed:
    ax2  = ax*ax
    ay2  = ay*ay
    b2   = b*b
    Ex2  = ax2-b2
    Ee2  = ax2-ay2
    lex2 = Ex2/ax2
    lee2 = Ee2/ax2
    mex  = one-lex2
    mee  = one-lee2

    nom = mee*z
    xme = mee*x
    dex = xme*xme+y*y
    den = mex*sqrt(dex)
    rot = sqrt(dex)

    if (den==zero)  then
        latitude = halfpi
        longitude = zero
    else
        if (nom<=den) then
            latitude = atan(nom/den)
        else
            latitude = halfpi-atan(den/nom)
        end if
        if (y<=xme) then
            den = xme+rot
            longitude = 2.0_wp*atan(y/den)
        else
            den = y+rot
            longitude = halfpi - 2.0_wp*atan(xme/den)
        end if
    end if

end subroutine xyz2fl

!****************************************************************
!>
!  Numerical solution to polynomial equation using Newton-Raphson method

pure function solve_polynomial(B, x0, error) result(x)

    real(wp),dimension(0:6),intent(in) :: B !! Polynomial `B = B(6) x^6 + B(5) x^5 + ... + B(0)`
    real(wp),intent(in) :: x0 !! Initial point
    real(wp),intent(in) :: error !! Maximum error
    real(wp) :: x !! root found after applying Newton-Raphson method to `B`
                  !! The function returns the value when the correction
                  !! is smaller than error.

    real(wp) :: f,fp,corr
    integer :: i, j !! counter

    integer,parameter :: maxiter = 100 !! maximum number of iterations

    x = x0
    do i = 1, maxiter
        f  = B(6)
        do j = 5, 0, -1
            if (j==5) then
                fp = f
            else
                fp = x*fp + f
            end if
            f  = x*f + B(j)
        end do
        if (fp==0.0_wp) exit ! singular point
        corr = f/fp
        x = x - corr
        if (abs(corr)<=error) exit
    end do

end function solve_polynomial

!****************************************************************
!>
!  Horner's method to compute `B(x-c)` in terms of `B(x)`.

pure subroutine horner(B, c, BB)

  real(wp),dimension(0:6),intent(in) :: B !! Polynomial `B = B(6) x^6 + B(5) x^5 + ... + B(0)`
  real(wp),intent(in) :: c
  real(wp),dimension(0:6),intent(out) :: BB !! Polynomial `BB` such that `B(x-c) = BB(x)`

  integer :: i,j !! counters

  BB = B

  do i = 0,6
    do j = 5,i,-1
      BB(j) = BB(j) - BB(j+1)*c
    end do
  end do

end subroutine horner

!******************************************************************
!>
!  Cartesian to Geodetic I
!
!### See also
!  * [[CartesianIntoGeodeticII]]
!
!### Reference
!  * Gema Maria Diaz-Toca, Leandro Marin, Ioana Necula,
!    "Direct transformation from Cartesian into geodetic coordinates on a triaxial ellipsoid"
!    Computers & Geosciences, Volume 142, September 2020, 104551.
!    [link](https://www.sciencedirect.com/science/article/pii/S0098300420305410?via%3Dihub),
!    [C++ code](https://data.mendeley.com/datasets/s5f6sww86x/2)

subroutine CartesianIntoGeodeticI(ax, ay, az, r, latitude, longitude, altitude, error)

    real(wp),intent(in) :: ax, ay, az !! semiaxes of the celestial body: ax>ay>az
    real(wp),dimension(3),intent(in) :: r !! cartesian coordinates of the considered point
                                          !! in the first octant: xG, yG, zG with (xG,yG,zG)<>(0,0,0)
    real(wp),intent(out) :: latitude, longitude, altitude !! geodetic coordinates of the considered point
    real(wp),intent(in) :: error !! Values smaller than error treated as 0.0

    real(wp) :: ax2,ay2,az2,ax4,ay4,az4,b5,b4,b3,b3x,b3y,b3z,&
                b2,b2x,b2y,b2z,b1,b1x,b1y,b1z,b0,b0x,b0y,b0z,eec,exc
    real(wp) :: xg2,yg2,zg2,aux,xG,yG,zG
    real(wp) :: xE,yE,zE,k,B(0:6),BB(0:6)
    logical :: done

    call special_cases(ax,ay,az,r(1),r(2),r(3),latitude,longitude,altitude,done)
    if (done) return

    ! Computations independent of xG,yG,zG. They can be precomputed, if necessary.
    ax2 = ax*ax
    ay2 = ay*ay
    az2 = az*az
    ax4 = ax2*ax2
    ay4 = ay2*ay2
    az4 = az2*az2
    b5  = 2.0_wp*(ax2+ay2+az2)
    b4  = ax4 + 4.0_wp*ax2*ay2 + ay4 + 4.0_wp*ax2*az2 + 4.0_wp*ay2*az2 + az4
    b3  = 2.0_wp*ax4*ay2 + 2.0_wp*ax2*ay4 + 2.0_wp*ax4*az2 + 8.0_wp*ax2*ay2*az2 + 2.0_wp*ay4*az2 + 2.0_wp*ax2*az4 + 2.0_wp*ay2*az4
    b3x = - 2.0_wp*ax2*ay2 - 2.0_wp*ax2*az2
    b3y = - 2.0_wp*ax2*ay2 - 2.0_wp*ay2*az2
    b3z = - 2.0_wp*ay2*az2 - 2.0_wp*ax2*az2
    b2  = 4.0_wp*ax4*ay2*az2 + 4.0_wp*ax2*ay4*az2 + ax4*az4 + 4.0_wp*ax2*ay2*az4 + ax4*ay4 + ay4*az4
    b2x = -ax2*ay4 -4.0_wp*ax2*ay2*az2 -ax2*az4
    b2y = -ax4*ay2 -4.0_wp*ax2*ay2*az2 -ay2*az4
    b2z = -ax4*az2 -4.0_wp*ax2*ay2*az2 -ay4*az2
    b1  = 2.0_wp*ax4*ay4*az2 + 2.0_wp*ax4*ay2*az4 + 2.0_wp*ax2*ay4*az4
    b1x = - 2.0_wp*ax2*ay4*az2 - 2.0_wp*ax2*ay2*az4
    b1y = - 2.0_wp*ax4*ay2*az2 - 2.0_wp*ax2*ay2*az4
    b1z = - 2.0_wp*ax4*ay2*az2 - 2.0_wp*ax2*ay4*az2
    b0  = ax4*ay4*az4
    b0x = - ax2*ay4*az4
    b0y = - ax4*ay2*az4
    b0z = - ax4*ay4*az2
    eec = (ax2-ay2)/ax2
    exc = (ax2-az2)/ax2

    ! Computations dependant of xG, yG, zG
    xG = abs(r(1))
    yG = abs(r(2))
    zG = abs(r(3))
    xg2 = xG*xG
    yg2 = yG*yG
    zg2 = zG*zG
    aux = xg2/ax2+yg2/ay2+zg2/az2

    B = [b0+b0x*xg2+b0y*yg2+b0z*zg2, &
         b1+b1x*xg2+b1y*yg2+b1z*zg2, &
         b2+b2x*xg2+b2y*yg2+b2z*zg2, &
         b3+b3x*xg2+b3y*yg2+b3z*zg2, &
         b4-(ax2*xg2+ay2*yg2+az2*zg2), &
         b5, &
         1.0_wp ]

  if (abs(aux-1.0_wp) < error) then ! The point is on the ellipsoid

    xE = xG
    yE = yG
    zE = zG

  else if (aux > 1.0_wp) then ! The point is outside the ellipsoid

    k = solve_polynomial(B,(xg2+yg2+zg2)/3.0_wp,error)
    xE = ax2*xG/(ax2+k)
    yE = ay2*yG/(ay2+k)
    zE = az2*zG/(az2+k)

  else if (zG > 0.0_wp) then ! The point  is inside the ellipsoid and zG>0

    call horner(B,az2,BB) ! B(x-az2) = BB(x)
    k = solve_polynomial(BB,(xg2+yg2+zg2)/3.0_wp+az2,error) - az2
    xE = ax2*xG/(ax2+k)
    yE = ay2*yG/(ay2+k)
    zE = az2*zG/(az2+k)

  else if (xG > 0.0_wp .and. yG > 0.0_wp) then ! The point is inside the ellipsoid and zG=0, yG > 0, xG > 0

    call horner(B,ay2,BB)
    k = solve_polynomial(BB,(xg2+yg2+zg2)/3.0_wp+ay2,error) - ay2
    xE = ax2*xG/(ax2+k)
    yE = ay2*yG/(ay2+k)
    zE = 0.0_wp

  else if (xG < error .and. yG > 0.0_wp) then

    xE = 0.0_wp
    yE = ay
    zE = 0.0_wp

  else if (xG > 0.0_wp .and. yG < error) then

    xE = ax
    yE = 0.0_wp
    zE = 0.0_wp

  end if

  ! Computing longitude
  if (xG > 0.0_wp) then
    longitude = atan(yE/((1.0_wp-eec)*xE))
  else if (yG > 0.0_wp) then
    longitude = halfpi
  else
    longitude = huge(1.0_wp)  ! undefined
  end if

  ! Computing latitude
  if (xE > 0.0_wp .or. yE > 0.0_wp) then
    latitude = atan((1.0_wp-eec)/(1.0_wp-exc)*zE/norm2([xE*(1.0_wp-eec),yE]))
  else
    latitude = halfpi
  end if

  ! Computing altitude
  if (aux>=1.0_wp) then
    altitude = norm2([xE-xG,yE-yG,zE-zG])
  else
    altitude = -norm2([xE-xG,yE-yG,zE-zG])
  end if

  call philambda_quadrant(r(1), r(2), r(3), latitude, longitude)

  end subroutine CartesianIntoGeodeticI

!******************************************************************
!>
!  Cartesian into Geodetic II
!
!### See also
!  * [[CartesianIntoGeodeticI]]

subroutine CartesianIntoGeodeticII(ax, ay, az, r, latitude, longitude, altitude, error)

    real(wp),intent(in) :: ax, ay, az !! semiaxes of the celestial body: ax>ay>az
    real(wp),dimension(3),intent(in) :: r !! cartesian coordinates of the considered point
                                          !! in the first octant: xG, yG, zG with (xG,yG,zG)<>(0,0,0)
    real(wp),intent(out) :: latitude, longitude, altitude !! geodetic coordinates of the considered point
    real(wp),intent(in) :: error !! Values smaller than error treated as 0.0

    real(wp) :: aymaz,aypaz,axmaz,axpaz,axpaz2,ax2,ay2,az2,ax4,ay4,az4,az6,&
                az8,temp0,temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,&
                temp9,tempa,az6ax2,az6ay2,tempb,maz10,excc,eecc
    real(wp) :: xg2,yg2,zg2,zgxg2,zgyg2,zg3,zg4,aux,xG,yG,zG
    real(wp) :: xE,yE,zE,k,B(0:6)
    logical :: done

    call special_cases(ax,ay,az,r(1),r(2),r(3),latitude,longitude,altitude,done)
    if (done) return

    ! Computations independent of xG,yG,zG. They can be precomputed, if necessary.
    aymaz = ay-az
    aypaz = ay+az
    axmaz = ax-az
    axpaz = ax+az
    axpaz2 = axpaz*axpaz
    ax2 = ax*ax
    ay2 = ay*ay
    az2 = az*az
    ax4 = ax2*ax2
    ay4 = ay2*ay2
    az4 = az2*az2
    az6 = az4*az2
    az8 = az4*az4
    temp0 = aymaz*aymaz*aypaz*aypaz*axmaz*axmaz*axpaz2
    temp1 = 2*az2*aymaz*aypaz*axmaz*axpaz*(ax2+ay2-2*az2)
    temp2 = -az2*(ax4*ay4-2*ax4*ay2*az2+ax4*az4-2*ax2*ay4*az2+4*ax2*ay2*az4-2*ay2*az6+az8-2*ax2*az6+ay4*az4)
    temp3 = -az2*(-ax2*ay4+2*ax2*ay2*az2-ax2*az4)
    temp4 = -az2*(-ax4*ay2+2*ax2*ay2*az2-ay2*az4)
    temp5 = -az2*(-ax4*az2-4*ax2*ay2*az2+6*ax2*az4-ay4*az2+6*ay2*az4-6*az6)
    temp6 =  -2*az4*(ax4*ay2-ax4*az2+ax2*ay4-4*ax2*ay2*az2+3*ax2*az4-ay4*az2+3*ay2*az4-2*az6)
    temp7 = -2*az4*(-ax2*ay2+ax2*az2)
    temp8 = -2*az4*(-ax2*ay2+ay2*az2)
    temp9 = -2*az4*(-ax2*az2-ay2*az2+2*az4)
    tempa = -az6*(ax4+4*ax2*ay2-6*ax2*az2+ay4-6*ay2*az2+6*az4)
    az6ax2 = az6*ax2
    az6ay2 = az6*ay2
    tempb = -2*az8*(ax2+ay2-2*az2)
    maz10 = -az6*az4
    excc = (ax2-az2)/(ax2)
    eecc = (ax2-ay2)/(ax2)

    xG = abs(r(1))
    yG = abs(r(2))
    zG = abs(r(3))
    xg2 = xG*xG
    yg2 = yG*yG
    zg2 = zG*zG
    zgxg2 = zG*xg2
    zgyg2 = zG*yg2
    zg3 = zg2*zG
    zg4 = zg2*zg2
    aux = xg2/ax2+yg2/ay2+zg2/az2

  if (abs(aux-1.0_wp) < error) then ! The point is on the ellipsoid

    xE = xG
    yE = yG
    zE = zG

  else if (zG > error) then ! The point is inside or outside the ellipsoid with zG != 0

    B(6) = temp0
    B(5) = temp1*zG
    B(4) = temp2+temp3*xg2+temp4*yg2+temp5*zg2
    B(3) = zG*temp6+temp7*zgxg2+temp8*zgyg2+temp9*zg3
    B(2) = zg2*(tempa+az6ax2*xg2+az6ay2*yg2+az8*zg2)
    B(1) = tempb*zg3
    B(0) = maz10*zg4

    k = solve_polynomial(B,az*zG/norm2([xG,yG,zG]),error)
    xE = ax2*xG*k/(ax2*k-az2*k+az2*zG)
    yE = ay2*yG*k/(ay2*k-az2*k+az2*zG)
    zE = k

  else if (yG > error) then

    B = [-ay4*ay2*zg2, &
         -2*ay4*(ax2-ay2)*zG, &
         -ay2*(ax4-2*ax2*ay2-ax2*yg2+ay4-ay2*zg2), &
         2*ay2*(ax2-ay2)*zG, &
         ax4+ay4-2*ax2*ay2, &
         0.0_wp, &
         0.0_wp ]

    k = solve_polynomial(B,ay*yG/norm2([xG,yG,zG]),error)
    xE = k*ax2*xG/(ax2*k-ay2*k+ay2*yG)
    yE = k
    zE = 0.0_wp

  else

    xE = ax
    yE = 0.0_wp
    zE = 0.0_wp

  end if

  ! Computing longitude

  if (xG > 0.0_wp) then
    longitude = atan(yE/((1.0_wp-eecc)*xE))
  else if (yG > 0.0_wp) then
    longitude = halfpi
  else
    longitude = huge(1.0_wp) ! undefined
  end if

  ! Computing latitude

  if (xE > 0.0_wp .or. yE > 0.0_wp) then
    latitude = atan((1.0_wp-eecc)/(1.0_wp-excc)*zE/norm2([xE*(1.0_wp-eecc),yE]))
  else
    latitude = halfpi
  end if

  ! Computing altitude

  if (aux>=1.0) then
    altitude = norm2([xE-xG,yE-yG,zE-zG])
  else
    altitude = -norm2([xE-xG,yE-yG,zE-zG])
  end if

  call philambda_quadrant(r(1), r(2), r(3), latitude, longitude)

end subroutine CartesianIntoGeodeticII

!********************************************************************************
!>
!  Cartesian to geodetic for Triaxial Ellipsoid.
!
!### References
!  * S. Bektas, "Geodetic Computations on Triaxial Ellipsoid",
!    International Journal of Mining Science (IJMS),
!    Volume 1, Issue 1, June 2015, p 25-34

    subroutine cartesian_to_geodetic_triaxial_2(a,b,c,r,eps,phi,lambda,h)

    implicit none

    real(wp),intent(in) :: a    !! ellipsoid radii `a >= b >= c`
    real(wp),intent(in) :: b    !! ellipsoid radii `a >= b >= c`
    real(wp),intent(in) :: c    !! ellipsoid radii `a >= b >= c`
    real(wp),dimension(3),intent(in) :: r !! Cartesian coordinates (x,y,z)
    real(wp),intent(in)  :: eps  !! convergence tolerance
    real(wp),intent(out) :: phi !! latitude (rad)
    real(wp),intent(out) :: lambda !! longitude (rad)
    real(wp),intent(out) :: h !! altitude

    integer,parameter :: maxiter = 20 !! maximum number of iterations

    integer :: i  !! iteration counter
    real(wp),dimension(3,3) :: AA
    real(wp),dimension(3) :: bvec, xvec
    real(wp) :: a2,b2,c2,x,y,z,ex2,ee2,e,f,g,xo,yo,zo,j11,j12,j21,j23,rmag
    logical :: success

    x = r(1)
    y = r(2)
    z = r(3)

    if (a<b .or. b<c) error stop 'error in cartesian_to_geodetic_triaxial_2: invalid a,b,c'
    call special_cases(a,b,c,x,y,z,phi,lambda,h,success)
    if (success) return

    rmag = norm2(r)
    a2  = a*a
    b2  = b*b
    c2  = c*c
    ex2 = (a2-c2)/a2
    ee2 = (a2-b2)/a2
    E   = 1.0_wp/a2
    F   = 1.0_wp/b2
    G   = 1.0_wp/c2
    xo  = a*x/rmag
    yo  = b*y/rmag
    zo  = c*z/rmag

    do i = 1, maxiter

        j11 = F*yo-(yo-y)*E
        j12 = (xo-x)*F-E*xo
        j21 = G*zo-(zo-z)*E
        j23 = (xo-x)*G-E*xo

        ! solve the linear system:
        AA = reshape(-[j11,j21,2.0_wp*E*xo,&
                       j12,0.0_wp,2.0_wp*F*yo,&
                       0.0_wp,j23,2.0_wp*G*zo], [3,3])
        bvec = [ (xo-x)*F*yo-(yo-y)*E*xo, &
                 (xo-x)*G*zo-(zo-z)*E*xo, &
                 E*xo**2+F*yo**2+G*zo**2-1.0_wp ]
        call linear_solver(AA,bvec,xvec,success)
        if (.not. success) then
            write(*,*) 'error in cartesian_to_geodetic_triaxial_2: matrix is singular'
            phi    = 0.0_wp
            lambda = 0.0_wp
            h      = 0.0_wp
            return
        end if
        xo = xo + xvec(1)
        yo = yo + xvec(2)
        zo = zo + xvec(3)

        if (maxval(abs(xvec))<eps) exit

    end do

    ! outputs:
    phi = atan(zo*(1.0_wp-ee2)/(1.0_wp-ex2)/sqrt((1.0_wp-ee2)**2*xo**2+yo**2))
    lambda = atan2(yo, (1.0_wp-ee2)*xo)
    h = sign(1.0_wp,z-zo)*sign(1.0_wp,zo)*sqrt((x-xo)**2+(y-yo)**2+(z-zo)**2)

    contains

    subroutine linear_solver(a,b,x,success)

        !!  Solve the 3x3 system: `A * x = b`
        !!  Reference: https://caps.gsfc.nasa.gov/simpson/software/m33inv_f90.txt

        implicit none

        real(wp),dimension(3,3),intent(in) :: a
        real(wp),dimension(3),intent(in)   :: b
        real(wp),dimension(3),intent(out)  :: x
        logical,intent(out)                :: success

        real(wp) :: det !! determinant of a
        real(wp),dimension(3,3) :: adj !! adjoint of a
        real(wp),dimension(3,3) :: ainv !! inverse of a

        det =   a(1,1)*a(2,2)*a(3,3)  &
              - a(1,1)*a(2,3)*a(3,2)  &
              - a(1,2)*a(2,1)*a(3,3)  &
              + a(1,2)*a(2,3)*a(3,1)  &
              + a(1,3)*a(2,1)*a(3,2)  &
              - a(1,3)*a(2,2)*a(3,1)

        success = abs(det) > tiny(1.0_wp) ! check for singularity

        if (success) then
            adj(:,1) = [a(2,2)*a(3,3)-a(2,3)*a(3,2),&
                        a(2,3)*a(3,1)-a(2,1)*a(3,3),&
                        a(2,1)*a(3,2)-a(2,2)*a(3,1)]
            adj(:,2) = [a(1,3)*a(3,2)-a(1,2)*a(3,3),&
                        a(1,1)*a(3,3)-a(1,3)*a(3,1),&
                        a(1,2)*a(3,1)-a(1,1)*a(3,2)]
            adj(:,3) = [a(1,2)*a(2,3)-a(1,3)*a(2,2),&
                        a(1,3)*a(2,1)-a(1,1)*a(2,3),&
                        a(1,1)*a(2,2)-a(1,2)*a(2,1)]
            ainv = adj/det
            x = matmul(ainv,b)
        else
            x = zero
        end if

        end subroutine linear_solver

    end subroutine cartesian_to_geodetic_triaxial_2
!********************************************************************************

!********************************************************************************
!>
!  Special cases for lat/lon/altitude

    subroutine special_cases(a,b,c,x,y,z,phi,lambda,h,done)

    real(wp),intent(in)  :: a      !! ellipsoid radii `a >= b >= c`
    real(wp),intent(in)  :: b      !! ellipsoid radii `a >= b >= c`
    real(wp),intent(in)  :: c      !! ellipsoid radii `a >= b >= c`
    real(wp),intent(in)  :: x      !! Cartesian x coordinate
    real(wp),intent(in)  :: y      !! Cartesian y coordinate
    real(wp),intent(in)  :: z      !! Cartesian z coordinate
    real(wp),intent(out) :: phi    !! latitude (rad)
    real(wp),intent(out) :: lambda !! longitude (rad)
    real(wp),intent(out) :: h      !! altitude
    logical,intent(out)  :: done   !! true if one of the special cases was computed

    logical :: x0, y0, z0

    real(wp),parameter :: zero_tol  = 10.0_wp * epsilon(1.0_wp)  !! zero tolerance for singularities

    x0 = abs(x) <= zero_tol
    y0 = abs(y) <= zero_tol
    z0 = abs(z) <= zero_tol

    if (x0 .and. y0 .and. z0) then ! center of the body
        phi    = zero
        lambda = zero
        h      = -c    ! just pick this value
        done = .true.
        return
    else if (x0 .and. y0) then ! (on the z-axis)
        if (z>=zero) then
            phi    = halfpi
            lambda = zero
            h      = z-c
        else
            phi    = -halfpi
            lambda = zero
            h      = -(z+c)
        end if
        done = .true.
        return
    else if (x0 .and. z0) then  ! on the y-axis
        if (y>=zero) then
            phi    = zero
            lambda = halfpi
            h      = y-b
        else
            phi    = zero
            lambda = -halfpi
            h      = -(y+b)
        end if
        done = .true.
        return
    else if (y0 .and. z0) then  ! on the x-axis
        if (x>=zero) then
            phi    = zero
            lambda = zero
            h      = x-a
        else
            phi    = zero
            lambda = pi
            h      = -(x+a)
        end if
        done = .true.
        return
    end if

    phi    = zero
    lambda = zero
    h      = zero
    done = .false.

    end subroutine special_cases
!*****************************************************************************************

!*****************************************************************************************
!>
!  Unit test for the [[direct]] and [[inverse]] geodetic routines.

    subroutine direct_inverse_test()

    implicit none

    !Ellipsoid : GRS80 / WGS84 (NAD83)
    real(wp),parameter :: a  = 6378137.0000_wp      !! Equatorial radius
    real(wp),parameter :: rf = 298.25722210088_wp   !! Inverse flattening
    real(wp),parameter :: f  = 1.0_wp / rf          !! flattening

    real(wp) :: glat1,glon1,glat2,glon2,faz,baz,s,sig,lam,glat2_,glon2_,baz_
    integer :: it, kind

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' direct_inverse_test'
    write(*,*) '---------------'
    write(*,*) ''

    !specify two points:
    glat1 = 0.523599_wp
    glon1 = 1.74533_wp
    glat2 = 0.698132_wp
    glon2 = 2.0944_wp

    call inverse(a,rf,glat1,glon1,glat2,glon2,faz,baz,s,it,sig,lam,kind)
    call direct(a,f,glat1,glon1,faz,s,glat2_,glon2_,baz_)

    write(*,*) 'lat error: ', glat2_ - glat2
    write(*,*) 'lon error: ', glon2_ - glon2
    write(*,*) 'baz error: ', baz_ - baz

    end subroutine direct_inverse_test
!*****************************************************************************************

!*****************************************************************************************
    end module geodesy_module
!*****************************************************************************************
