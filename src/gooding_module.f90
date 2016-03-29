!*****************************************************************************************
!> author: Jacob Williams
!
!  Gooding's Kepler and universal elements conversion routines.
!
!# Notes
!  The Gooding universal elements are:
!
!   * `alpha` - mu/a [km^2/s^2]
!   * `rp`    - periapsis radius [km]
!   * `inc`   - inclination [rad]
!   * `raan`  - right ascension of the ascending node [rad]
!   * `w`     - argument of periapsis [rad]
!   * `tau`   - time since last periapsis passage [sec]
!
!# References
!  1. A. W. Odell, R. H. Gooding, "Procedures for solving Kepler's equation"
!     Celestial Mechanics 38 (1986), 307-334.
!  2. R. H. Gooding, "On universal elements, and conversion procedures
!     to and from position and velocity"
!     Celestial Mechanics 44 (1988), 283-298.
!  3. R. H. Gooding, A. W. Odell. "The hyperbolic Kepler equation
!     (and the elliptic equation revisited)"
!     Celestial Mechanics 44 (1988), 267-282.

    module gooding_module

    use kind_module,    only: wp
    use numbers_module

    implicit none

    private

    !constants:
    real(wp),parameter :: ntwo   = -two
    real(wp),parameter :: pineg  = -pi
    real(wp),parameter :: half   = 0.5_wp
    real(wp),parameter :: halfpi = pi/two
    real(wp),parameter :: athird = one/three
    real(wp),parameter :: asixth = one/six

    public :: els3pv,pv3els
    public :: ekepl,ekepl1,ekepl2,emkepl,emkep,shkepl,shmkep
    public :: propagate

    contains
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!
!  Basic two-body propagator using the Gooding universal element routines.

    pure subroutine propagate(mu, rv0, dt, rvf)

    implicit none

    real(wp),intent(in)               :: mu    !! grav. parameter [km^3/s^2]
    real(wp),dimension(6),intent(in)  :: rv0   !! initial state [km, km/s]
    real(wp),intent(in)               :: dt    !! time step [sec]
    real(wp),dimension(6),intent(out) :: rvf   !! final state [km, km/s]

    real(wp),dimension(6) :: e

    !convert to elements, increment time,
    ! then convert back to cartesian:

    call pv3els(mu, rv0, e)
    e(6) = e(6) + dt
    call els3pv(mu, e, rvf)

    end subroutine propagate
!*****************************************************************************************

!*****************************************************************************************
!>
!  Kepler's equation, `em = ekepl - (1 - e1)*sin(ekepl)`,
!  with `e1` in range 1 to 0 inclusive, solved accurately
!  (based on ekepl3, but entering `e1`, not `e`)

    pure function ekepl(em, e1)

    implicit none

    real(wp) :: ekepl
    real(wp),intent(in) :: em
    real(wp),intent(in) :: e1

    real(wp) :: emr,ee,e,w,fdd,fddd,f,fd,dee
    integer :: iter

    real(wp),parameter :: sw = 0.25_wp

    !range-reduce em to lie in range -pi to pi
    emr = mod(em,twopi)
    if (emr<pineg) emr = emr + twopi
    if (emr>pi) emr = emr - twopi
    ee = emr

    if (ee/=zero) then

        if (ee<zero) ee = -ee

        !(emr is range-reduced em & ee is absolute value of emr)
        !starter by first solving cubic equation
        e = one - e1
         w = dcbsol(e,two*e1, three*ee)

        !effectively interpolate in emr (absolute value)
        ee = (ee*ee + (pi - ee)*w)/pi
        if (emr<zero) ee = -ee

        !do two iterations of halley, each followed by newton
        do iter=1,2
            fdd = e*sin(ee)
            fddd = e*cos(ee)
            if (ee*ee/six + e1 >= sw) then
                f = (ee - fdd) - emr
                fd = one - fddd
            else
                f = emkep(e1,ee) - emr
                fd = two*e*sin(half*ee)**2 + e1
            end if
            dee = f*fd/(half*f*fdd - fd*fd)
            f = f + dee*(fd + half*dee*(fdd + athird*dee*fddd))
            !to reduce the danger of underflow replace the last line by
            !    w = fd + half*dee*(fdd + athird*dee*fddd)
            fd = fd + dee*(fdd + half*dee*fddd)
            ee = ee + dee - f/fd
            !if replacing as above, then also replace the last line by
            !ee = ee - (f - dee*(fd - w))/fd
        end do

    end if

    !range-expand
    ekepl = ee + (em - emr)

    end function ekepl
!*****************************************************************************************

!*****************************************************************************************
!>
!  Solve kepler's equation, `em = ekepl - e*sin(ekepl)`,
!  with legendre-based starter and halley iterator
!  (function has also been used under the name eafkep)

    pure function ekepl1(em, e)

    implicit none

    real(wp) :: ekepl1
    real(wp),intent(in) :: em
    real(wp),intent(in) :: e

    real(wp) :: c,s,psi,xi,eta,fd,fdd,f

    real(wp),parameter :: testsq = 1.0e-8_wp

    c = e*cos(em)
    s = e*sin(em)
    psi = s/sqrt(one - c - c + e*e)

    do

         xi = cos(psi)
         eta = sin(psi)
         fd = (one - c*xi) + s*eta
         fdd = c*eta + s*xi
         f = psi - fdd
         psi = psi - f*fd/(fd*fd - half*f*fdd)
         if (f*f < testsq) exit

     end do

     ekepl1 = em + psi

    end function ekepl1
!*****************************************************************************************

!*****************************************************************************************
!>
!  Kepler's equation, `em = ekepl - e*sin(ekepl)` with
!  e in range 0 to 1 inclusive, solved accurately

    pure function ekepl2(em, e)

    implicit none

    real(wp)            :: ekepl2
    real(wp),intent(in) :: em
    real(wp),intent(in) :: e

    real(wp) :: emr, ee, w, e1, fdd, fddd, f, fd, dee
    logical :: l
    integer :: iter

    real(wp),parameter :: sw = 0.1_wp
    real(wp),parameter :: a  = (pi-one)**2/(pi+two/three)
    real(wp),parameter :: b  = two*(pi-asixth)**2/(pi+two/three)

    !range-reduce em to line in range -pi to pi
    emr = mod(em,twopi)
    if (emr<pineg) emr = emr + twopi
    if (emr>pi) emr = emr - twopi
    ee = emr

    if (ee/=zero) then

        if (ee<zero) ee = -ee

        !(emr is range-reduced em & ee is absolute value of emr)
        !started for e = 1 by cube root of bilinear function
        if (ee<asixth) then
            ee = (six*ee)**athird
        else
            w = pi - ee
            ee = pi - a*w/(b - w)
        end if
        if (emr<zero) ee = -ee

        !interpolate for e
        ee = emr + (ee - emr)*e

        !do two iterations of halley, each followed by newton
        e1 = one - e
        l = (e1 + ee*ee/six) >= sw
        do iter=1,2
            fdd = e*sin(ee)
            fddd = e*cos(ee)
            if (l) then
                f = (ee - fdd) - emr
                fd = one - fddd
            else
                f = emkepl(e,ee) - emr
                fd = e1 + two*e*sin(half*ee)**2
            end if
            dee = f*fd/(half*f*fdd - fd*fd)
            f = f + dee*(fd + half*dee*(fdd + athird*dee*fddd))
            !to reduce the danger of underflow replace the last line by
            !w = fd + half*dee*(fdd + athird*dee*fddd)
            fd = fd + dee*(fdd + half*dee*fddd)
            ee = ee + dee - f/fd
            !if replacing as above, then also replace the last line by
            !ee = ee - (f - dee*(fd - w))/fd
        end do

    end if

    !range-expand
    ekepl2 = ee + (em - emr)

    end function ekepl2
!*****************************************************************************************

!*****************************************************************************************
!>
!  Accurate computation of `ee - e*sin(ee)`
!  when (e, ee) is close to (1, 0)
!
!@note must not be used for large ee (absolute)
!      as then rounding worse not better

    pure function emkepl(e, ee)

    implicit none

    real(wp)            :: emkepl
    real(wp),intent(in) :: e
    real(wp),intent(in) :: ee

    real(wp) :: x, ee2, term, d, x0

    x    = (one - e)*sin(ee)
    ee2  = -ee*ee
    term = ee
    d    = zero

    do

        d = d + two
        term = term*ee2/(d*(d + one))
        x0 = x
        x = x - term
        if (x==x0) exit

    end do

    emkepl = x

    end function emkepl
!*****************************************************************************************

!*****************************************************************************************
!>
!  Similar to emkepl, except input is `1-e`.

    pure function emkep(e1,ee)

    implicit none

    real(wp)            :: emkep
    real(wp),intent(in) :: e1
    real(wp),intent(in) :: ee

    real(wp) :: x, ee2, term, d, x0

    x    = e1*sin(ee)
    ee2  = -ee*ee
    term = ee
    d    = zero

    do

        d = d + two
        term = term*ee2/(d*(d + one))
        x0 = x
        x = x - term
        if (x==x0) exit

    end do

    emkep = x

    end function emkep
!*****************************************************************************************

!*****************************************************************************************
!>
!  Equation `el = shkepl + (g1 - 1)*asinh(shkepl)`,
!  with g1 in range 0 to 1 inclusive, solved accurately.

    pure function shkepl (el, g1)

    implicit none

    real(wp)            :: shkepl
    real(wp),intent(in) :: el
    real(wp),intent(in) :: g1

    real(wp) :: s,g,cl,al,w,s0,s1,s2,s3,fdd,fddd,f,fd,ds,stemp
    integer  :: iter

    real(wp),parameter :: sw=half

    s = el

    if (el/=zero) then

        !started based on lagrange's theorem
        g = one - g1
        cl = sqrt(one + el**2)
        al = asinh(el)
        w = g**2*al/cl**3
        s = one - g/cl
        s = el + g*al/dcubrt(s**3 + w*el*(1.5_wp - g/0.75_wp))

        !two iterations (at most) of halley-then-newton process
        do iter=1,2
            s0 = s*s
            s1 = s0 + one
            s2 = sqrt(s1)
            s3 = s1*s2
            fdd = g*s/s3
            fddd = g*(one - two*s0)/(s1*s3)
            if (asixth*s0 + g1 >= sw) then
                f = (s - g*asinh(s)) - el
                fd = one - g/s2
            else
                f = shmkep(g1, s) - el
                fd = (s0/(s2 + one) + g1)/s2
            end if
            ds = f*fd/(half*f*fdd - fd*fd)
            stemp = s + ds
            if (stemp==s) exit
            f = f + ds*(fd + half*ds*(fdd + athird*ds*fddd))
            fd = fd + ds*(fdd + half*ds*fddd)
             s = stemp - f/fd
        end do

    end if

    shkepl = s

     end function shkepl
!*****************************************************************************************

!*****************************************************************************************
!>
!  Accurate computation of `s - (1 - g1)*asinh(s)`
!  when (g1, s) is close to (0, 0)

    pure function shmkep (g1, s)

    implicit none

    real(wp)            :: shmkep
    real(wp),intent(in) :: g1
    real(wp),intent(in) :: s

    real(wp) :: g,t,tsq,x,term,twoi1,x0

    g = one - g1
    t = s/(one + sqrt(one + s*s))
    tsq = t*t
    x = s*(g1 + g*tsq)
    term = two*g*t
    twoi1 = one

    do

        twoi1 = twoi1 + two
        term = term*tsq
        x0 = x
        x = x - term/twoi1
        if (x==x0) exit

     end do

     shmkep = x

     end function shmkep
!*****************************************************************************************

!*****************************************************************************************
!>
!  Algorithm for two-dimensional conversion
!  from orbital elements to position and velocity.

    pure subroutine els2pv (gm, al, q, om, tau, r, u, vr, vt)

    implicit none

    real(wp),intent(in)  :: gm    !! grav. parameter [km^3/s^2]
    real(wp),intent(in)  :: al    !! alpha [km^2/s^2]
    real(wp),intent(in)  :: q     !! periapsis distance [km]
    real(wp),intent(in)  :: om    !! argument of periapsis relative to assumed reference direction [rad]
    real(wp),intent(in)  :: tau   !! time from periapsis [sec]
    real(wp),intent(out) :: r     !! radial distance [km]
    real(wp),intent(out) :: u     !! angle from reference direction [rad]
    real(wp),intent(out) :: vr    !! radial velocity [km/2]
    real(wp),intent(out) :: vt    !! transverse velocity >=0 [km/s]

    real(wp) :: d,h,v,e1,e,ep1,alp,rtal,em,ee2,s2,c2,emv,s,c

    if (al==zero) then

        !(parabola - gm cannot be zero)

        d = dcbsol(half/gm, q, 1.5_wp*gm*tau)
        r = q + half*d*d/gm
        h = sqrt(two*gm*q)
        v = two*atan2(d,h)

    else

        !(ellipse or hyperbola)

        e1 = al*q
        e = gm - e1
        ep1 = gm + e
        h = sqrt(q*ep1)
        alp = abs(al)
        rtal = sqrt(alp)
        !(last 6 items could be saved if repeating gm, al & q)

        em = tau*alp*rtal
        if (al>zero) then

            !(ellipse - gm cannot be zero)
            ! make sure e1 argument to ekepl is between [0,1]
            ee2 = half*ekepl(em/gm, max(zero,min(one,e1/gm)))
            s2 = sin(ee2)
            c2 = cos(ee2)
            r = q + two*e*s2*s2/al
            d = two*e*s2*c2/rtal
            v = two*atan2(ep1*s2, h*rtal*c2)
            emv = em/gm - v
            v = v + fourpi*sign(real(int(abs(emv/fourpi) + half),wp), emv)

        else

            !(hyperbola)
            s = shkepl(em/e, -e1/e)
            s2 = s*s
            c = sqrt(one + s2)
            s2 = s2/(c+one)
            r = q - e*s2/al
            d = e*s/rtal
            v = atan2(s*h*rtal, -gm*s2 - e1)

        end if

    end if

    !(all orbits)
    u = om + v
    vr = d/r
    vt = h/r

    end subroutine els2pv
!*****************************************************************************************

!*****************************************************************************************
!>
!  Algorithm for three-dimensional conversion
!  from orbital elements to position and velocity

    pure subroutine els3pv (gm, e, pv)

    implicit none

    real(wp),intent(in)               :: gm   !! grav. parameter [km^3/sec^2]
    real(wp),dimension(6),intent(in)  :: e    !! [al, q, ei, bom, om, tau]
    real(wp),dimension(6),intent(out) :: pv   !! [x, y, z, xdot, ydot, zdot]

    real(wp) :: x,y,z,xdot,ydot,zdot,al,q,ei,bom,om,tau
    real(wp) :: r,u,vr,vt,c,s,x1,x2,y1,y2

    if (all(e==zero)) then

        pv = zero

    else

        al  = e(1)
        q   = e(2)
        ei  = e(3)
        bom = e(4)
        om  = e(5)
        tau = e(6)

        call els2pv (gm, al, q, om, tau, r, u, vr, vt)

        c = cos(u)
        s = sin(u)
        x1 = r*c
        y1 = r*s
        x2 = vr*c - vt*s
        y2 = vr*s + vt*c
        c = cos(ei)
        s = sin(ei)
        z = y1*s
        y1 = y1*c
        zdot = y2*s
        y2 = y2*c
        c = cos(bom)
        s = sin(bom)
        x = x1*c - y1*s
        y = x1*s + y1*c
        xdot = x2*c - y2*s
        ydot = x2*s + y2*c

        pv(1) = x
        pv(2) = y
        pv(3) = z
        pv(4) = xdot
        pv(5) = ydot
        pv(6) = zdot

    end if

    end subroutine els3pv
!*****************************************************************************************

!*****************************************************************************************
!>
!  Algorithm for two-dimensional conversion
!  from position and velocity to orbital elements.

    pure subroutine pv2els (gm, r, u, vr, vt, al, q, om, tau)

    implicit none

    real(wp),intent(in)  :: gm    !! grav. parameter [km^3/s^2]
    real(wp),intent(in)  :: r     !! radial distance [km]
    real(wp),intent(in)  :: u     !! angle from assumed reference direction [rad]
    real(wp),intent(in)  :: vr    !! radial velocity [km/2]
    real(wp),intent(in)  :: vt    !! transverse velocity >=0 [km/s]
    real(wp),intent(out) :: al    !! alpha: gm/a [km^2/s^2]
    real(wp),intent(out) :: q     !! periapsis distance [km]
    real(wp),intent(out) :: om    !! argument of periapsis relative to reference direction [rad]
    real(wp),intent(out) :: tau   !! time from periapsis [sec]

    real(wp) :: esq1,es,eses,ec,ecec,esq,e,v,e1
    real(wp) :: eh,em,ecesq,en,adj,vsq,rtal,d,h,p,alp

    real(wp),parameter   :: sw = 0.25_wp
    logical,parameter    :: l = .false.

    !(all orbits)
    vsq = vr*vr + vt*vt
    al = two*gm/r - vsq
    alp = abs(al)
    rtal = sqrt(alp)
    d = r*vr
    h = r*vt
    p = h*h
    esq1 = p*al
    es = d*rtal
    eses = es*es
    ec = r*vsq - gm
    ecec = ec*ec
    if (al>zero) then
        !(one esq formula superior for the ellipse)
        esq = ecec + eses
    else
        !(different formula superior for the hyperbola)
        esq = gm*gm - esq1
    end if
    e = sqrt(esq)
    q = p/(gm + e)
    if (al==zero) then
        !(parabola)
        tau = d*(two*q + r)/(three*gm)
        v = two*atan2(vr, vt)
    else if (e==zero) then
        !(circle)
        tau = zero
        v = zero
    else
        !(ellipse or hyperbola)
        e1 = al*q
        if (al>zero) then
            !(ellipse)
            eh = atan2(es, ec)
            if (gm*eh*eh/six + e1 >= gm*sw) then
                !(general case)
                em = gm*eh - es
                ecesq = gm*ec - esq
            else
                !(for e1 & eh both near zero)
                em = gm*emkep (e1/gm, eh)
                ecesq = (esq1*ecec - esq*eses)/(esq + gm*ec)
            end if
        else
            !(hyperbola)
            eh = asinh(es/e)
            if (gm*eh*eh/six - e1 >= gm*sw) then
                !(general case)
                em = es - gm*eh
                ecesq = esq - gm*ec
            else
                !(for e1 & eh both near zero)
                em = e*shmkep(-e1/e, es/e)
                ecesq = -(esq1*ecec + esq*eses)/(esq + gm*ec)
            end if
        end if
            !(ellipse or hyperbola still)
            en = alp*rtal
            tau = em/en
            v = atan2(es*h*rtal, ecesq)
    end if

    !(all orbits)
    om = u - v

    !
    !  note: the following is never executed... set l=true and test...
    !

    if (l .and. al>zero) then
        !(for ellipse, adjust revolutions if required (using l))
        adj = twopi*sign(real(int(abs(om/twopi) + half),wp), om)
        om = om - adj
        tau = tau + adj/en
    end if

    end subroutine pv2els
!*****************************************************************************************

!*****************************************************************************************
!>
!  Algorithm for three-dimensional conversion
!  from position and velocity to orbital elements.

    pure subroutine pv3els (gm, pv, e)

    implicit none

    real(wp),intent(in)               :: gm   !! grav. parameter [km^3/s^2]
    real(wp),dimension(6),intent(in)  :: pv   !! [x, y, z, xdot, ydot, zdot]
    real(wp),dimension(6),intent(out) :: e    !! [al, q, ei, bom, om, tau]

    real(wp) :: x,y,z,xdot,ydot,zdot,al,q,ei,bom,om,tau,xsqysq,&
                rsq,r,vr,hx,hy,hz,hsq,u,vt,bx,by,bz,w,h

    if (all(pv==zero)) then

        e = zero

    else

        x    = pv(1)
        y    = pv(2)
        z    = pv(3)
        xdot = pv(4)
        ydot = pv(5)
        zdot = pv(6)

        xsqysq = x*x + y*y
        rsq = xsqysq + z*z
        r = sqrt(rsq)
        vr = (x*xdot + y*ydot + z*zdot)/r
        hx = y*zdot - z*ydot
        hy = z*xdot - x*zdot
        hz = x*ydot - y*xdot
        hsq = hx*hx + hy*hy + hz*hz
        if (hsq==zero) then
            !(rectilinear orbit)
            ei = halfpi
            if (xsqysq==zero) then
                !(axial orbit)
                bom = zero
            else
                !(general rectilinear orbit)
                bom = atan2(y, x)
            end if
            u = atan2(z, sqrt(xsqysq))
            vt = zero
        else
            !(non-degenerate orbit)
            bx = hy*z - hz*y
            by = hz*x - hx*z
            bz = hx*y - hy*x
            hx = y*bz - z*by
            hy = z*bx - x*bz
            hz = x*by - y*bx
            w = hx*hx + hy*hy
            h = sqrt(w + hz*hz)
            ei = atan2(sqrt(w), hz)
            if (w==zero) then
                !(orbit in reference plane)
                bom = zero
                u = atan2(y*sign(one,hz), x)
            else
                !(general orbit)
                bom = atan2(hx, -hy)
                u = atan2(h*z, rsq*bz)
            end if
            vt = h/(r*rsq)
        end if
        call pv2els (gm, r, u, vr, vt, al, q, om, tau)

        e(1) = al
        e(2) = q
        e(3) = ei
        e(4) = bom
        e(5) = om
        e(6) = tau

    end if

    end subroutine pv3els
!*****************************************************************************************

!*****************************************************************************************
!>
!  Solution to `a*x**3 + 3*b*x - 2c = 0`, where
!  `a` and `b**3 + a*c**2` are both non-negative
!  (zero generated, in lieu of infinity, if `a = b = 0`)

    pure function dcbsol (a, b, c) result(x)

    implicit none

    real(wp)               :: x
    real(wp),intent(in)    :: a
    real(wp),intent(in)    :: b
    real(wp),intent(in)    :: c

    real(wp) :: bsq,d

    if (a==zero .and. b==zero .or. c==zero) then
        x = zero
    else
        bsq = b*b
        d = sqrt(a) * abs(c)
        d = dcubrt(d + sqrt(b*bsq + d*d))**2
        x = two * c / (d + b + bsq / d)
    end if

    end function dcbsol
!*****************************************************************************************

!*****************************************************************************************
!>
!  Cube root computed accurately, by incorporating
!  one Newton-Raphson iteration.

    pure function dcubrt(x) result(c)

    implicit none

    real(wp)            :: c
    real(wp),intent(in) :: x

    real(wp) :: y

    if (x==zero) then
        c = zero
    else
        y = abs(x)
        c = y**athird
        c = c - athird*(c - y/c**2)
        c = sign(c,x)
    end if

    end function dcubrt
!*****************************************************************************************

!*****************************************************************************************
    end module gooding_module
!*****************************************************************************************
