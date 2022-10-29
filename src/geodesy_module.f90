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

    public :: geodetic_to_cartesian
    public :: geodetic_to_cartesian_triaxial

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

subroutine cartesian_to_geodetic_triaxial(ax, ay, b, Xi, Yi, Zi, tol, phi, lambda, h)

    real(wp),intent(in) :: ax !! semiaxes (0 < b <= ay <= ax)
    real(wp),intent(in) :: ay !! semiaxes (0 < b <= ay <= ax)
    real(wp),intent(in) :: b  !! semiaxes (0 < b <= ay <= ax)
    real(wp),intent(in) :: Xi !! Cartesian coordinates
    real(wp),intent(in) :: Yi !! Cartesian coordinates
    real(wp),intent(in) :: Zi !! Cartesian coordinates
    real(wp),intent(in) :: tol !! tolerance (may be set to zero)
    real(wp),intent(out) :: phi !! geodetic latitude (radians)
    real(wp),intent(out) :: lambda !! geodetic longitude (radians)
    real(wp),intent(out) :: h !! geodetic height

    real(wp) :: kx,ky,cx,cy,cz,XX,YY,ZZ,x,y,z,Xo,Yo,Zo,m,Mm
    integer :: n
    real(wp),dimension(3) :: d

    kx = (ax*ax-b*b)/ax
    ky = (ay*ay-b*b)/ay

    cx = (ax*ax)/(b*b)
    cy = (ay*ay)/(b*b)
    cz = (ax*ax)/(ay*ay)

    XX = abs(Xi)
    YY = abs(Yi)
    ZZ = abs(Zi)

    ! Compute geodetic latitude/longitude
    if (ZZ == zero) then
        if (XX == zero .and. YY == zero) then
            x = zero
            y = zero
            z = b
        elseif (ky*XX*ky*XX+kx*YY*kx*YY < kx*ky*kx*ky) then
            x = ax*XX/kx
            y = ay*YY/ky
            z = b*sqrt(one-((x*x)/(ax*ax))-((y*y)/(ay*ay)))
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
                z = b*sqrt(one-((x*x)/(ax*ax))-((y*y)/(ay*ay)))
            else
                z = ZZ/(one+m)
            end if
        end if

    end if

    call xyz2philambda(ax, ay, b, x, y, z, phi, lambda)

    ! Compute geodetic height
    d = [XX-x, YY-y, ZZ-z]
    h = norm2(d)
    if ((XX+YY+ZZ) < (x+y+z)) h = -h

    call philambda_quadrant(Xi, Yi, Zi, phi, lambda)

end subroutine cartesian_to_geodetic_triaxial

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

subroutine geodetic_to_cartesian_triaxial(ax, ay, b, phi, lambda, h, xi, yi, zi)

    real(wp),intent(in) :: ax !! semiaxes (0 < b <= ay <= ax)
    real(wp),intent(in) :: ay !! semiaxes (0 < b <= ay <= ax)
    real(wp),intent(in) :: b  !! semiaxes (0 < b <= ay <= ax)
    real(wp),intent(in) :: phi !! geodetic latitude (radians)
    real(wp),intent(in) :: lambda !! geodetic longitude (radians)
    real(wp),intent(in) :: h !! geodetic height
    real(wp),intent(out) :: xi !! Cartesian coordinates
    real(wp),intent(out) :: yi !! Cartesian coordinates
    real(wp),intent(out) :: zi !! Cartesian coordinates

    real(wp) :: ee2,ex2,N,cp,sp,cl,sl

    cp  = cos(phi)
    sp  = sin(phi)
    cl  = cos(lambda)
    sl  = sin(lambda)
    ee2 = (ax*ax-ay*ay)/(ax*ax)
    ex2 = (ax*ax-b*b)/(ax*ax)
    N   = ax/sqrt(one - ex2*sp*sp - ee2*cp*cp*sl*sl)

    xi = (N+h)*cp*cl
    yi = (N*(one-ee2)+h)*cp*sl
    zi = (N*(one-ex2)+h)*sp

end subroutine geodetic_to_cartesian_triaxial

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

subroutine philambda2xyz(ax, ay, b, phi, lambda, x, y, z)

    real(wp),intent(in) :: ax, ay, b, phi, lambda
    real(wp),intent(out) :: x, y, z

    real(wp) :: ee2,ex2,N

    ee2 = (ax*ax-ay*ay)/(ax*ax)
    ex2 = (ax*ax-b*b)/(ax*ax)

    N = ax/sqrt(one-ex2*sin(phi)*sin(phi)-ee2*cos(phi)*cos(phi)*sin(lambda)*sin(lambda))

    x = N*cos(phi)*cos(lambda)
    y = N*(one-ee2)*cos(phi)*sin(lambda)
    z = N*(one-ex2)*sin(phi)

end subroutine philambda2xyz

!*****************************************************************************************
!>

subroutine philambda_quadrant(XX, YY, ZZ, phi, lambda)

    real(wp),intent(in) :: XX, YY, ZZ
    real(wp),intent(inout) :: phi, lambda

    if (ZZ < zero) then
        phi = -phi
    end if

    if (XX >= zero) then
        if (YY >= zero) then
            lambda = lambda
        else
            lambda = -lambda
        end if
    else
        if (YY >= zero) then
            lambda = pi-lambda
        else
            lambda = lambda-pi
        end if
    end if

end subroutine philambda_quadrant

!*****************************************************************************************
!>

subroutine xyz2philambda(ax, ay, b, x, y, z, phi, lambda)

    real(wp),intent(in) :: ax, ay, b, x, y, z
    real(wp),intent(out) :: phi, lambda

    real(wp) :: ee2,ex2,s0,SS0,Sphi,Cphi,Slambda,Clambda,&
                Den,P,L,x0,y0,z0,D,NN,onemee2,onemex2,term,&
                at,bt,onemee2x
    integer :: n
    real(wp),dimension(3,2) :: J
    real(wp),dimension(3,1) :: dl,UU
    real(wp),dimension(2,2) :: Nmat
    real(wp),dimension(2,1) :: u, dx
    real(wp),dimension(1,1) :: tmp
    real(wp),dimension(2,2) :: NmatTmp

    ee2 = (ax*ax - ay*ay)/(ax*ax)
    ex2 = (ax*ax - b*b)/(ax*ax)
    onemee2 = one - ee2
    onemex2 = one - ex2
    onemee2x = onemee2*x
    term = sqrt(onemee2*onemee2x*x + y*y)
    at = onemee2*z
    bt = onemex2*term

    ! Initial values
    if (at <= bt) then
        phi = atan(at/bt)
    else
        phi = halfpi - atan(bt/at)
    end if

    if (x == zero .and. y == zero) then
        lambda = zero
    elseif (y <= onemee2x) then
        lambda = two*atan(y/(onemee2x+term))
    else
        lambda = halfpi - two*atan(onemee2x/(y+term))
    end if

    s0 = zero
    do n = 1, 100
        SS0 = s0

        ! Design Matrix J
        Sphi = sin(phi)
        Cphi = cos(phi)
        Slambda = sin(lambda)
        Clambda = cos(lambda)

        NN = ax/sqrt(onemex2*Sphi*Sphi-ee2*Cphi*Cphi*Slambda*Slambda)

        Den = two*(onemex2*Sphi*Sphi-ee2*Cphi*Cphi*Slambda*Slambda)**(3.0_wp/two)
        P = -ax*(ex2-ee2*Slambda*Slambda)*sin(two*phi)/Den
        L = -ax*ee2*Cphi*Cphi*sin(2*lambda)/Den

        J = zero

        J(1,1) = (P*Cphi-NN*Sphi)*Clambda
        J(2,1) = onemee2*(P*Cphi-NN*Sphi)*Slambda
        J(3,1) = onemex2*(P*Sphi+NN*Cphi)

        J(1,2) = (L*Clambda-NN*Slambda)*Cphi
        J(2,2) = onemee2*(L*Slambda+NN*Clambda)*Cphi
        J(3,2) = onemex2*L*Sphi

        ! Vector dl
        call philambda2xyz(ax, ay, b, phi, lambda, x0, y0, z0)
        dl(:,1) = [x-x0, y-y0, z-z0]

        ! Solution
        Nmat    = matmul(transpose(J),J)
        u       = matmul(transpose(J),dl)
        D       = Nmat(1,1)*Nmat(2,2) - Nmat(1,2)*Nmat(2,1)
        NmatTmp = transpose(reshape([Nmat(2,2), -Nmat(1,2), -Nmat(2,1), Nmat(1,1)], [2,2]))
        dx      = matmul((NmatTmp / D) , u)
        phi     = phi + dx(1,1)
        lambda  = lambda + dx(2,1)
        UU      = matmul(J,dx) - dl
        tmp     = sqrt(matmul(transpose(UU),UU))
        s0      = tmp(1,1)

        if (s0 == SS0) exit

    end do

end subroutine xyz2philambda
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
