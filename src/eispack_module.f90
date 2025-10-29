!*****************************************************************************************
!>
!  Refactored SLATEC/EISPACK routines for computing eigenvalues and eigenvectors.

    module eispack_module

    use kind_module, only: wp

    implicit none

    private

    public :: compute_eigenvalues_and_eigenvectors
    public :: compute_real_eigenvalues_and_normalized_eigenvectors
    public :: eispack_test

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Balance a real general matrix and isolate eigenvalues
!  whenever possible.
!
!  This subroutine is a translation of the ALGOL procedure BALANCE,
!  NUM. MATH. 13, 293-304(1969) by Parlett and Reinsch.
!  HANDBOOK FOR AUTO. COMP., Vol.II-LINEAR ALGEBRA, 315-326(1971).
!
!  This subroutine balances a REAL matrix and isolates
!  eigenvalues whenever possible.
!
!### On INPUT
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameter, A, as declared in the calling program
!          dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix A.  N is an INTEGER variable.
!          N must be less than or equal to NM.
!
!        A contains the input matrix to be balanced.  A is a
!          two-dimensional REAL array, dimensioned A(NM,N).
!
!### On OUTPUT
!
!        A contains the balanced matrix.
!
!        LOW and IGH are two INTEGER variables such that A(I,J)
!          is equal to zero if
!           (1) I is greater than J and
!           (2) J=1,...,LOW-1 or I=IGH+1,...,N.
!
!        SCALE contains information determining the permutations and
!          scaling factors used.  SCALE is a one-dimensional REAL array,
!          dimensioned SCALE(N).
!
!     Suppose that the principal submatrix in rows LOW through IGH
!     has been balanced, that P(J) denotes the index interchanged
!     with J during the permutation step, and that the elements
!     of the diagonal matrix used are denoted by D(I,J).  Then
!        SCALE(J) = P(J),    for J = 1,...,LOW-1
!                 = D(J,J),      J = LOW,...,IGH
!                 = P(J)         J = IGH+1,...,N.
!     The order in which the interchanges are made is N to IGH+1,
!     then 1 TO LOW-1.
!
!     Note that 1 is returned for IGH if IGH is zero formally.
!
!  The ALGOL procedure EXC contained in BALANCE appears in
!  BALANC in line.  (Note that the ALGOL roles of identifiers
!  K,L have been reversed.)
!
!  Questions and comments should be directed to B. S. Garbow,
!  Applied Mathematics Division, ARGONNE NATIONAL LABORATORY
!
!### References
!  * B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!    Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!    system Routines - EISPACK Guide, Springer-Verlag,
!    1976.
!
!### Revision history
!  * Author: Smith, B. T., et al.
!  * 760101  DATE WRITTEN
!  * 890831  Modified array declarations.  (WRB)
!  * 890831  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 920501  Reformatted the REFERENCES section.  (WRB)

    subroutine balanc(Nm, n, a, Low, Igh, Scale)

    implicit none

    integer :: i, j, k, l, m, n, jj, Nm, Igh, Low, iexc, igo, igo1, igo2
    real(wp) :: a(Nm, *), Scale(*)
    real(wp) :: c, f, g, r, s, b2
    logical :: noconv

    real(wp),parameter :: radix = 16.0_wp

    b2 = radix*radix
    k = 1
    l = n
    igo = 1
    igo1 = 0
    igo2 = 1
    ! IN-LINE PROCEDURE FOR ROW AND
    ! COLUMN EXCHANGE
    do while (igo2 == 1)
       igo2 = 0
       if (igo1 == 1) then
          Scale(m) = j
          if (j /= m) then
             do i = 1, l
                f = a(i, j)
                a(i, j) = a(i, m)
                a(i, m) = f
             enddo
             do i = k, n
                f = a(j, i)
                a(j, i) = a(m, i)
                a(m, i) = f
             enddo
          endif
          if (iexc == 2) then
             ! SEARCH FOR COLUMNS ISOLATING AN EIGENVALUE
             ! AND PUSH THEM LEFT
             k = k + 1
             igo = 0
          else
             ! SEARCH FOR ROWS ISOLATING AN EIGENVALUE
             ! AND PUSH THEM DOWN
             if (l == 1) then
                Low = k
                Igh = l
                return
             end if
             l = l - 1
          endif
       end if
       ! FOR J=L STEP -1 UNTIL 1 DO --
       igo1 = 1
       if (igo == 1) then
          do jj = 1, l
             igo = 1
             j = l + 1 - jj
             do i = 1, l
                if (i /= j) then
                   if (a(j, i) /= 0.0_wp) then
                      igo = 0
                      exit
                   end if
                endif
             enddo
             if (igo == 0) cycle
             m = l
             iexc = 1
             igo2 = 1
             exit
          enddo
          if (igo2 == 1) cycle
       end if
       do j = k, l
          igo = 1
          do i = k, l
             if (i /= j) then
                if (a(i, j) /= 0.0_wp) then
                   igo = 0
                   exit
                end if
             endif
          enddo
          if (igo == 0) cycle
          m = k
          iexc = 2
          igo2 = 1
          exit
       enddo
       if (igo2 == 1) cycle
    end do
    ! NOW BALANCE THE SUBMATRIX IN ROWS K TO L
    do i = k, l
       Scale(i) = 1.0_wp
    enddo
    ! ITERATIVE LOOP FOR NORM REDUCTION
    noconv = .true.
    do while (noconv)
       noconv = .false.
       do i = k, l
          c = 0.0_wp
          r = 0.0_wp
          do j = k, l
             if (j /= i) then
                c = c + abs(a(j, i))
                r = r + abs(a(i, j))
             endif
          enddo
          ! GUARD AGAINST ZERO C OR R DUE TO UNDERFLOW
          if (c /= 0.0_wp .and. r /= 0.0_wp) then
             g = r/radix
             f = 1.0_wp
             s = c + r

             do while (c < g)
                f = f*radix
                c = c*b2
             end do
             g = r*radix
             do while (c >= g)
                f = f/radix
                c = c/b2
             end do
             ! NOW BALANCE
             if ((c + r)/f < 0.95_wp*s) then
                g = 1.0_wp/f
                Scale(i) = Scale(i)*f
                noconv = .true.
                do j = k, n
                   a(i, j) = a(i, j)*g
                enddo
                do j = 1, l
                   a(j, i) = a(j, i)*f
                enddo
             endif
          endif
       enddo
    end do
    Low = k
    Igh = l

 end subroutine balanc
!*****************************************************************************************

!*****************************************************************************************
!>
!  Form the eigenvectors of a real general matrix from the
!  eigenvectors of matrix output from BALANC.
!
!  This subroutine is a translation of the ALGOL procedure BALBAK,
!  NUM. MATH. 13, 293-304(1969) by Parlett and Reinsch.
!  HANDBOOK FOR AUTO. COMP., Vol.II-LINEAR ALGEBRA, 315-326(1971).
!
!  This subroutine forms the eigenvectors of a REAL GENERAL
!  matrix by back transforming those of the corresponding
!  balanced matrix determined by  BALANC.
!
!### On Input
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameter, Z, as declared in the calling program
!          dimension statement.  NM is an INTEGER variable.
!
!        N is the number of components of the vectors in matrix Z.
!          N is an INTEGER variable.  N must be less than or equal
!          to NM.
!
!        LOW and IGH are INTEGER variables determined by  BALANC.
!
!        SCALE contains information determining the permutations and
!          scaling factors used by  BALANC.  SCALE is a one-dimensional
!          REAL array, dimensioned SCALE(N).
!
!        M is the number of columns of Z to be back transformed.
!          M is an INTEGER variable.
!
!        Z contains the real and imaginary parts of the eigen-
!          vectors to be back transformed in its first M columns.
!          Z is a two-dimensional REAL array, dimensioned Z(NM,M).
!
!### On Output
!
!        Z contains the real and imaginary parts of the
!          transformed eigenvectors in its first M columns.
!
!     Questions and comments should be directed to B. S. Garbow,
!     Applied Mathematics Division, ARGONNE NATIONAL LABORATORY
!
!### References
!  * B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!    Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!    system Routines - EISPACK Guide, Springer-Verlag,
!    1976.
!
!### Revision History
!  * Author: Smith, B. T., et al.
!  * 760101  DATE WRITTEN
!  * 890831  Modified array declarations.  (WRB)
!  * 890831  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 920501  Reformatted the REFERENCES section.  (WRB)

 subroutine balbak(Nm, n, Low, Igh, Scale, m, z)

    implicit none

    integer :: i, j, k, m, n, ii, Nm, Igh, Low
    real(wp) :: Scale(*), z(Nm, *)
    real(wp) :: s

    if (m /= 0) then
       if (Igh /= Low) then

          do i = Low, Igh
             s = Scale(i)
             ! LEFT HAND EIGENVECTORS ARE BACK TRANSFORMED
             ! IF THE FOREGOING STATEMENT IS REPLACED BY
             ! S=1.0_wp/SCALE(I).
             do j = 1, m
                z(i, j) = z(i, j)*s
             enddo

          enddo
       endif
       ! FOR I=LOW-1 STEP -1 UNTIL 1,
       ! IGH+1 STEP 1 UNTIL N DO --
       do ii = 1, n
          i = ii
          if (i < Low .or. i > Igh) then
             if (i < Low) i = Low - ii
             k = Scale(i)
             if (k /= i) then

                do j = 1, m
                   s = z(i, j)
                   z(i, j) = z(k, j)
                   z(k, j) = s
                enddo
             endif
          endif

       enddo
    endif

 end subroutine balbak
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute the complex quotient of two complex numbers.
!
!  Complex division, (CR,CI) = (AR,AI)/(BR,BI)
!
!### Revision History
!  * 811101  DATE WRITTEN
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900402  Added TYPE section.  (WRB)

 subroutine cdiv(Ar, Ai, Br, Bi, Cr, Ci)

    implicit none

    real(wp) :: Ar, Ai, Br, Bi, Cr, Ci

    real(wp) :: s, ars, ais, brs, bis

    s = abs(Br) + abs(Bi)
    ars = Ar/s
    ais = Ai/s
    brs = Br/s
    bis = Bi/s
    s = brs**2 + bis**2
    Cr = (ars*brs + ais*bis)/s
    Ci = (ais*brs - ars*bis)/s

 end subroutine cdiv
!*****************************************************************************************

!*****************************************************************************************
!>
!  Reduce a real general matrix to upper Hessenberg form
!  using stabilized elementary similarity transformations.
!
!  This subroutine is a translation of the ALGOL procedure ELMHES,
!  NUM. MATH. 12, 349-368(1968) by Martin and Wilkinson.
!  HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
!
!  Given a REAL GENERAL matrix, this subroutine
!  reduces a submatrix situated in rows and columns
!  LOW through IGH to upper Hessenberg form by
!  stabilized elementary similarity transformations.
!
!### On Input
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameter, A, as declared in the calling program
!          dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix, A.  N is an INTEGER variable.
!          N must be less than or equal to NM.
!
!        LOW and IGH are two INTEGER variables determined by the
!          balancing subroutine  BALANC.  If  BALANC  has not been
!          used, set LOW=1 and IGH equal to the order of the matrix, N.
!
!        A contains the input matrix.  A is a two-dimensional REAL
!          array, dimensioned A(NM,N).
!
!### On Output
!
!        A contains the upper Hessenberg matrix.  The multipliers which
!          were used in the reduction are stored in the remaining
!          triangle under the Hessenberg matrix.
!
!        INTV contains information on the rows and columns interchanged
!          in the reduction.  Only elements LOW through IGH are used.
!          INTV is a one-dimensional INTEGER array, dimensioned INTV(IGH).
!
!     Questions and comments should be directed to B. S. Garbow,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!
!### References
!  * B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!    Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!    system Routines - EISPACK Guide, Springer-Verlag,
!    1976.
!
!### Revision History
!  * Author: Smith, B. T., et al.
!  * 760101  DATE WRITTEN
!  * 890831  Modified array declarations.  (WRB)
!  * 890831  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 920501  Reformatted the REFERENCES section.  (WRB)

 subroutine elmhes(Nm, n, Low, Igh, a, Intv)

    implicit none

    integer :: i, j, m, n, la, Nm, Igh, kp1, Low, mm1, mp1
    real(wp) :: a(Nm, *)
    real(wp) :: x, y
    integer :: Intv(*)

    la = Igh - 1
    kp1 = Low + 1
    if (la >= kp1) then

       do m = kp1, la
          mm1 = m - 1
          x = 0.0_wp
          i = m

          do j = m, Igh
             if (abs(a(j, mm1)) > abs(x)) then
                x = a(j, mm1)
                i = j
             endif
          enddo

          Intv(m) = i
          if (i /= m) then
             ! INTERCHANGE ROWS AND COLUMNS OF A
             do j = mm1, n
                y = a(i, j)
                a(i, j) = a(m, j)
                a(m, j) = y
             enddo

             do j = 1, Igh
                y = a(j, i)
                a(j, i) = a(j, m)
                a(j, m) = y
             enddo
          endif
          ! END INTERCHANGE
          if (x /= 0.0_wp) then
             mp1 = m + 1

             do i = mp1, Igh
                y = a(i, mm1)
                if (y /= 0.0_wp) then
                   y = y/x
                   a(i, mm1) = y

                   do j = m, n
                      a(i, j) = a(i, j) - y*a(m, j)
                   enddo

                   do j = 1, Igh
                      a(j, m) = a(j, m) + y*a(j, i)
                   enddo
                endif

             enddo
          endif

       enddo
    endif

 end subroutine elmhes
!*****************************************************************************************

!*****************************************************************************************
!>
!  Accumulates the stabilized elementary similarity
!  transformations used in the reduction of a real general
!  matrix to upper Hessenberg form by ELMHES.
!
!  This subroutine is a translation of the ALGOL procedure ELMTRANS,
!  NUM. MATH. 16, 181-204(1970) by Peters and Wilkinson.
!  HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
!
!  This subroutine accumulates the stabilized elementary
!  similarity transformations used in the reduction of a
!  REAL GENERAL matrix to upper Hessenberg form by  ELMHES.
!
!### On Input
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, A and Z, as declared in the calling
!          program dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix A.  N is an INTEGER variable.
!          N must be less than or equal to NM.
!
!        LOW and IGH are two INTEGER variables determined by the
!          balancing subroutine  BALANC.  If  BALANC  has not been
!          used, set LOW=1 and IGH equal to the order of the matrix, N.
!
!        A contains the multipliers which were used in the reduction
!          by  ELMHES  in its lower triangle below the subdiagonal.
!          A is a two-dimensional REAL array, dimensioned A(NM,IGH).
!
!        INT contains information on the rows and columns interchanged
!          in the reduction by  ELMHES.  Only elements LOW through IGH
!          are used.  INT is a one-dimensional INTEGER array,
!          dimensioned INT(IGH).
!
!### On Output
!
!        Z contains the transformation matrix produced in the reduction
!          by  ELMHES.  Z is a two-dimensional REAL array, dimensioned
!          Z(NM,N).
!
!     Questions and comments should be directed to B. S. Garbow,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!
!### References
!  * B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!    Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!    system Routines - EISPACK Guide, Springer-Verlag,
!    1976.
!
!### Revision History
!  * Author: Smith, B. T., et al.
!  * 760101  DATE WRITTEN
!  * 890831  Modified array declarations.  (WRB)
!  * 890831  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 920501  Reformatted the REFERENCES section.  (WRB)

 subroutine eltran(Nm, n, Low, Igh, a, Int, z)
    implicit none

    integer i, j, n, kl, mm, mp, Nm, Igh, Low, mp1
    real(wp) a(Nm, *), z(Nm, *)
    integer Int(*)

    do i = 1, n
       do j = 1, n
          z(i, j) = 0.0_wp
       enddo
       z(i, i) = 1.0_wp
    enddo

    kl = Igh - Low - 1
    if (kl >= 1) then
       ! for mp=igh-1 step -1 until low+1 do --
       do mm = 1, kl
          mp = Igh - mm
          mp1 = mp + 1
          do i = mp1, Igh
             z(i, mp) = a(i, mp - 1)
          enddo
          i = Int(mp)
          if (i /= mp) then

             do j = mp, Igh
                z(mp, j) = z(i, j)
                z(i, j) = 0.0_wp
             enddo
             z(i, mp) = 1.0_wp
          endif
       enddo
    endif

 end subroutine eltran
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute the eigenvalues of a real upper Hessenberg matrix
!  using the QR method.
!
!  This subroutine is a translation of the ALGOL procedure HQR,
!  NUM. MATH. 14, 219-231(1970) by Martin, Peters, and Wilkinson.
!  HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 359-371(1971).
!
!  This subroutine finds the eigenvalues of a REAL
!  UPPER Hessenberg matrix by the QR method.
!
!### On Input
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameter, H, as declared in the calling program
!          dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix H.  N is an INTEGER variable.
!          N must be less than or equal to NM.
!
!        LOW and IGH are two INTEGER variables determined by the
!          balancing subroutine  BALANC.  If  BALANC  has not been
!          used, set LOW=1 and IGH equal to the order of the matrix, N.
!
!        H contains the upper Hessenberg matrix.  Information about
!          the transformations used in the reduction to Hessenberg
!          form by  ELMHES  or  ORTHES, if performed, is stored
!          in the remaining triangle under the Hessenberg matrix.
!          H is a two-dimensional REAL array, dimensioned H(NM,N).
!
!### On Output
!
!        H has been destroyed.  Therefore, it must be saved before
!          calling  HQR  if subsequent calculation and back
!          transformation of eigenvectors is to be performed.
!
!        WR and WI contain the real and imaginary parts, respectively,
!          of the eigenvalues.  The eigenvalues are unordered except
!          that complex conjugate pairs of values appear consecutively
!          with the eigenvalue having the positive imaginary part first.
!          If an error exit is made, the eigenvalues should be correct
!          for indices IERR+1, IERR+2, ..., N.  WR and WI are one-
!          dimensional REAL arrays, dimensioned WR(N) and WI(N).
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          J          if the J-th eigenvalue has not been
!                     determined after a total of 30*N iterations.
!                     The eigenvalues should be correct for indices
!                     IERR+1, IERR+2, ..., N.
!
!  Questions and comments should be directed to B. S. Garbow,
!  APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!
!### References
!  * B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!    Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!    system Routines - EISPACK Guide, Springer-Verlag,
!    1976.
!
!### Revision History
!  * Author: Smith, B. T., et al.
!  * 760101  DATE WRITTEN
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 890831  Modified array declarations.  (WRB)
!  * 890831  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 920501  Reformatted the REFERENCES section.  (WRB)

 subroutine hqr(Nm, n, Low, Igh, h, Wr, Wi, Ierr)

    implicit none

    integer :: i, j, k, l, m, n, en, ll, mm, na, Nm, Igh, &
               itn, its, Low, mp2, enm2, Ierr, gt
    real(wp) :: h(Nm, *), Wr(*), Wi(*)
    real(wp) :: p, q, r, s, t, w, x, y, zz, norm, s1, s2
    logical :: notlas

    gt = 0
    Ierr = 0
    norm = 0.0_wp
    k = 1
    ! STORE ROOTS ISOLATED BY BALANC
    ! AND COMPUTE MATRIX NORM
    do i = 1, n
       do j = k, n
          norm = norm + abs(h(i, j))
       enddo
       k = i
       if (i < Low .or. i > Igh) then
          Wr(i) = h(i, i)
          Wi(i) = 0.0_wp
       endif
    enddo

    en = Igh
    t = 0.0_wp
    itn = 30*n
    ! SEARCH FOR NEXT EIGENVALUES
    do while (en >= Low)
       gt = 0
       its = 0
       na = en - 1
       enm2 = na - 1
       ! LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
       ! FOR L=EN STEP -1 UNTIL LOW DO --
       do while (.true.)
          do ll = Low, en
             l = en + Low - ll
             if (l == Low) exit
             s = abs(h(l - 1, l - 1)) + abs(h(l, l))
             if (s == 0.0_wp) s = norm
             s2 = s + abs(h(l, l - 1))
             if (s2 == s) exit
          enddo
          ! FORM SHIFT
          x = h(en, en)
          if (l == en) then
             ! ONE ROOT FOUND
             Wr(en) = x + t
             Wi(en) = 0.0_wp
             en = na
             gt = 1
             exit
          else
             y = h(na, na)
             w = h(en, na)*h(na, en)
             if (l == na) exit
             if (itn == 0) then
                ! SET ERROR -- NO CONVERGENCE TO AN
                ! EIGENVALUE AFTER 30*N ITERATIONS
                Ierr = en
                return
             else
                if (its == 10 .or. its == 20) then
                   ! FORM EXCEPTIONAL SHIFT
                   t = t + x
                   do i = Low, en
                      h(i, i) = h(i, i) - x
                   enddo
                   s = abs(h(en, na)) + abs(h(na, enm2))
                   x = 0.75_wp*s
                   y = x
                   w = -0.4375_wp*s*s
                endif
                its = its + 1
                itn = itn - 1
                ! LOOK FOR TWO CONSECUTIVE SMALL
                ! SUB-DIAGONAL ELEMENTS.
                ! FOR M=EN-2 STEP -1 UNTIL L DO --
                do mm = l, enm2
                   m = enm2 + l - mm
                   zz = h(m, m)
                   r = x - zz
                   s = y - zz
                   p = (r*s - w)/h(m + 1, m) + h(m, m + 1)
                   q = h(m + 1, m + 1) - zz - r - s
                   r = h(m + 2, m + 1)
                   s = abs(p) + abs(q) + abs(r)
                   p = p/s
                   q = q/s
                   r = r/s
                   if (m == l) exit
                   s1 = abs(p)*(abs(h(m - 1, m - 1)) + abs(zz) + abs(h(m + 1, m + 1)))
                   s2 = s1 + abs(h(m, m - 1))*(abs(q) + abs(r))
                   if (s2 == s1) exit
                enddo
                mp2 = m + 2
                do i = mp2, en
                   h(i, i - 2) = 0.0_wp
                   if (i /= mp2) h(i, i - 3) = 0.0_wp
                enddo
                ! DOUBLE QR STEP INVOLVING ROWS L TO EN AND
                ! COLUMNS M TO EN
                do k = m, na
                   notlas = k /= na
                   if (k /= m) then
                      p = h(k, k - 1)
                      q = h(k + 1, k - 1)
                      r = 0.0_wp
                      if (notlas) r = h(k + 2, k - 1)
                      x = abs(p) + abs(q) + abs(r)
                      if (x == 0.0_wp) cycle
                      p = p/x
                      q = q/x
                      r = r/x
                   endif
                   s = sign(sqrt(p*p + q*q + r*r), p)
                   if (k == m) then
                      if (l /= m) h(k, k - 1) = -h(k, k - 1)
                   else
                      h(k, k - 1) = -s*x
                   endif
                   p = p + s
                   x = p/s
                   y = q/s
                   zz = r/s
                   q = q/p
                   r = r/p
                   ! ROW MODIFICATION
                   do j = k, en
                      p = h(k, j) + q*h(k + 1, j)
                      if (notlas) then
                         p = p + r*h(k + 2, j)
                         h(k + 2, j) = h(k + 2, j) - p*zz
                      endif
                      h(k + 1, j) = h(k + 1, j) - p*y
                      h(k, j) = h(k, j) - p*x
                   enddo
                   j = min(en, k + 3)
                   ! COLUMN MODIFICATION
                   do i = l, j
                      p = x*h(i, k) + y*h(i, k + 1)
                      if (notlas) then
                         p = p + zz*h(i, k + 2)
                         h(i, k + 2) = h(i, k + 2) - p*r
                      endif
                      h(i, k + 1) = h(i, k + 1) - p*q
                      h(i, k) = h(i, k) - p
                   enddo
                enddo
             endif
          endif
       enddo
       ! TWO ROOTS FOUND
       if (gt == 0) then
          p = (y - x)/2.0_wp
          q = p*p + w
          zz = sqrt(abs(q))
          x = x + t
          if (q < 0.0_wp) then
             ! COMPLEX PAIR
             Wr(na) = x + p
             Wr(en) = x + p
             Wi(na) = zz
             Wi(en) = -zz
          else
             ! REAL PAIR
             zz = p + sign(zz, p)
             Wr(na) = x + zz
             Wr(en) = Wr(na)
             if (zz /= 0.0_wp) Wr(en) = x - w/zz
             Wi(na) = 0.0_wp
             Wi(en) = 0.0_wp
          endif
          en = enm2
       end if
    enddo

 end subroutine hqr
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute the eigenvalues and eigenvectors of a real upper
!  Hessenberg matrix using QR method.
!
!  This subroutine is a translation of the ALGOL procedure HQR2,
!  NUM. MATH. 16, 181-204(1970) by Peters and Wilkinson.
!  HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
!
!  This subroutine finds the eigenvalues and eigenvectors
!  of a REAL UPPER Hessenberg matrix by the QR method.  The
!  eigenvectors of a REAL GENERAL matrix can also be found
!  if  ELMHES  and  ELTRAN  or  ORTHES  and  ORTRAN  have
!  been used to reduce this general matrix to Hessenberg form
!  and to accumulate the similarity transformations.
!
!### On Input
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, H and Z, as declared in the calling
!          program dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix H.  N is an INTEGER variable.
!          N must be less than or equal to NM.
!
!        LOW and IGH are two INTEGER variables determined by the
!          balancing subroutine  BALANC.  If  BALANC  has not been
!          used, set LOW=1 and IGH equal to the order of the matrix, N.
!
!        H contains the upper Hessenberg matrix.  H is a two-dimensional
!          REAL array, dimensioned H(NM,N).
!
!        Z contains the transformation matrix produced by  ELTRAN
!          after the reduction by  ELMHES, or by  ORTRAN  after the
!          reduction by  ORTHES, if performed.  If the eigenvectors
!          of the Hessenberg matrix are desired, Z must contain the
!          identity matrix.  Z is a two-dimensional REAL array,
!          dimensioned Z(NM,M).
!
!### On Output
!
!        H has been destroyed.
!
!        WR and WI contain the real and imaginary parts, respectively,
!          of the eigenvalues.  The eigenvalues are unordered except
!          that complex conjugate pairs of values appear consecutively
!          with the eigenvalue having the positive imaginary part first.
!          If an error exit is made, the eigenvalues should be correct
!          for indices IERR+1, IERR+2, ..., N.  WR and WI are one-
!          dimensional REAL arrays, dimensioned WR(N) and WI(N).
!
!        Z contains the real and imaginary parts of the eigenvectors.
!          If the J-th eigenvalue is real, the J-th column of Z
!          contains its eigenvector.  If the J-th eigenvalue is complex
!          with positive imaginary part, the J-th and (J+1)-th
!          columns of Z contain the real and imaginary parts of its
!          eigenvector.  The eigenvectors are unnormalized.  If an
!          error exit is made, none of the eigenvectors has been found.
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          J          if the J-th eigenvalue has not been
!                     determined after a total of 30*N iterations.
!                     The eigenvalues should be correct for indices
!                     IERR+1, IERR+2, ..., N, but no eigenvectors are
!                     computed.
!
!     Calls CDIV for complex division.
!
!  Questions and comments should be directed to B. S. Garbow,
!  APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!
!### References
!  * B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!    Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!    system Routines - EISPACK Guide, Springer-Verlag,
!    1976.
!
!### Revision History
!  * Smith, B. T., et al.
!  * 760101  DATE WRITTEN
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 890831  Modified array declarations.  (WRB)
!  * 890831  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 920501  Reformatted the REFERENCES section.  (WRB)

 subroutine hqr2(Nm, n, Low, Igh, h, Wr, Wi, z, Ierr)
    implicit none

    integer :: i, j, k, l, m, n, en, ii, jj, ll, mm, na, Nm, &
               nn, gt
    integer :: Igh, itn, its, Low, mp2, enm2, Ierr
    real(wp) :: h(Nm, *), Wr(*), Wi(*), z(Nm, *)
    real(wp) :: p, q, r, s, t, w, x, y, ra, sa, vi, vr, zz, &
                norm, s1, s2
    logical :: notlas

    gt = 0
    Ierr = 0
    norm = 0.0_wp
    k = 1
    ! STORE ROOTS ISOLATED BY BALANC
    ! AND COMPUTE MATRIX NORM
    do i = 1, n
       do j = k, n
          norm = norm + abs(h(i, j))
       enddo
       k = i
       if (i < Low .or. i > Igh) then
          Wr(i) = h(i, i)
          Wi(i) = 0.0_wp
       endif
    enddo

    en = Igh
    t = 0.0_wp
    itn = 30*n
    ! SEARCH FOR NEXT EIGENVALUES
    do while (en >= Low)
       gt = 0
       its = 0
       na = en - 1
       enm2 = na - 1
       ! LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
       ! FOR L=EN STEP -1 UNTIL LOW DO --
       do while (.true.)
          do ll = Low, en
             l = en + Low - ll
             if (l == Low) exit
             s = abs(h(l - 1, l - 1)) + abs(h(l, l))
             if (s == 0.0_wp) s = norm
             s2 = s + abs(h(l, l - 1))
             if (s2 == s) exit
          enddo
          ! FORM SHIFT
          x = h(en, en)
          if (l == en) then
             ! ONE ROOT FOUND
             h(en, en) = x + t
             Wr(en) = h(en, en)
             Wi(en) = 0.0_wp
             en = na
             gt = 1
             exit
          else
             y = h(na, na)
             w = h(en, na)*h(na, en)
             if (l == na) exit
             if (itn == 0) then
                ! SET ERROR -- NO CONVERGENCE TO AN
                ! EIGENVALUE AFTER 30*N ITERATIONS
                Ierr = en
                return
             else
                if (its == 10 .or. its == 20) then
                   ! FORM EXCEPTIONAL SHIFT
                   t = t + x
                   do i = Low, en
                      h(i, i) = h(i, i) - x
                   enddo
                   s = abs(h(en, na)) + abs(h(na, enm2))
                   x = 0.75_wp*s
                   y = x
                   w = -0.4375_wp*s*s
                endif
                its = its + 1
                itn = itn - 1
                ! LOOK FOR TWO CONSECUTIVE SMALL
                ! SUB-DIAGONAL ELEMENTS.
                ! FOR M=EN-2 STEP -1 UNTIL L DO --
                do mm = l, enm2
                   m = enm2 + l - mm
                   zz = h(m, m)
                   r = x - zz
                   s = y - zz
                   p = (r*s - w)/h(m + 1, m) + h(m, m + 1)
                   q = h(m + 1, m + 1) - zz - r - s
                   r = h(m + 2, m + 1)
                   s = abs(p) + abs(q) + abs(r)
                   p = p/s
                   q = q/s
                   r = r/s
                   if (m == l) exit
                   s1 = abs(p)*(abs(h(m - 1, m - 1)) + abs(zz) + abs(h(m + 1, m + 1)))
                   s2 = s1 + abs(h(m, m - 1))*(abs(q) + abs(r))
                   if (s2 == s1) exit
                enddo
                mp2 = m + 2
                do i = mp2, en
                   h(i, i - 2) = 0.0_wp
                   if (i /= mp2) h(i, i - 3) = 0.0_wp
                enddo
                ! DOUBLE QR STEP INVOLVING ROWS L TO EN AND
                ! COLUMNS M TO EN
                do k = m, na
                   notlas = k /= na
                   if (k /= m) then
                      p = h(k, k - 1)
                      q = h(k + 1, k - 1)
                      r = 0.0_wp
                      if (notlas) r = h(k + 2, k - 1)
                      x = abs(p) + abs(q) + abs(r)
                      if (x == 0.0_wp) cycle
                      p = p/x
                      q = q/x
                      r = r/x
                   endif
                   s = sign(sqrt(p*p + q*q + r*r), p)
                   if (k == m) then
                      if (l /= m) h(k, k - 1) = -h(k, k - 1)
                   else
                      h(k, k - 1) = -s*x
                   endif
                   p = p + s
                   x = p/s
                   y = q/s
                   zz = r/s
                   q = q/p
                   r = r/p
                   ! ROW MODIFICATION
                   do j = k, n
                      p = h(k, j) + q*h(k + 1, j)
                      if (notlas) then
                         p = p + r*h(k + 2, j)
                         h(k + 2, j) = h(k + 2, j) - p*zz
                      endif
                      h(k + 1, j) = h(k + 1, j) - p*y
                      h(k, j) = h(k, j) - p*x
                   enddo
                   j = min(en, k + 3)
                   ! COLUMN MODIFICATION
                   do i = 1, j
                      p = x*h(i, k) + y*h(i, k + 1)
                      if (notlas) then
                         p = p + zz*h(i, k + 2)
                         h(i, k + 2) = h(i, k + 2) - p*r
                      endif
                      h(i, k + 1) = h(i, k + 1) - p*q
                      h(i, k) = h(i, k) - p
                   enddo
                   ! ACCUMULATE TRANSFORMATIONS
                   do i = Low, Igh
                      !
                      p = x*z(i, k) + y*z(i, k + 1)
                      if (notlas) then
                         p = p + zz*z(i, k + 2)
                         z(i, k + 2) = z(i, k + 2) - p*r
                      endif
                      z(i, k + 1) = z(i, k + 1) - p*q
                      z(i, k) = z(i, k) - p
                      !
                   enddo
                enddo
             endif
          endif
       enddo
       if (gt == 1) cycle
       ! TWO ROOTS FOUND
       p = (y - x)/2.0_wp
       q = p*p + w
       zz = sqrt(abs(q))
       h(en, en) = x + t
       x = h(en, en)
       h(na, na) = y + t
       if (q < 0.0_wp) then
          ! COMPLEX PAIR
          Wr(na) = x + p
          Wr(en) = x + p
          Wi(na) = zz
          Wi(en) = -zz
       else
          ! REAL PAIR
          zz = p + sign(zz, p)
          Wr(na) = x + zz
          Wr(en) = Wr(na)
          if (zz /= 0.0_wp) Wr(en) = x - w/zz
          Wi(na) = 0.0_wp
          Wi(en) = 0.0_wp
          x = h(en, na)
          s = abs(x) + abs(zz)
          p = x/s
          q = zz/s
          r = sqrt(p*p + q*q)
          p = p/r
          q = q/r
          ! ROW MODIFICATION
          do j = na, n
             zz = h(na, j)
             h(na, j) = q*zz + p*h(en, j)
             h(en, j) = q*h(en, j) - p*zz
          enddo
          ! COLUMN MODIFICATION
          do i = 1, en
             zz = h(i, na)
             h(i, na) = q*zz + p*h(i, en)
             h(i, en) = q*h(i, en) - p*zz
          enddo
          ! ACCUMULATE TRANSFORMATIONS
          do i = Low, Igh
             zz = z(i, na)
             z(i, na) = q*zz + p*z(i, en)
             z(i, en) = q*z(i, en) - p*zz
          enddo
       endif
       en = enm2
    enddo
    ! ALL ROOTS FOUND.  BACKSUBSTITUTE TO FIND
    ! VECTORS OF UPPER TRIANGULAR FORM
    if (norm /= 0.0_wp) then
       ! FOR EN=N STEP -1 UNTIL 1 DO --
       do nn = 1, n
          en = n + 1 - nn
          p = Wr(en)
          q = Wi(en)
          na = en - 1
          if (q < 0) then
             ! END COMPLEX VECTOR
             ! COMPLEX VECTOR
             m = na
             ! LAST VECTOR COMPONENT CHOSEN IMAGINARY SO THAT
             ! EIGENVECTOR MATRIX IS TRIANGULAR
             if (abs(h(en, na)) <= abs(h(na, en))) then
                call cdiv(0.0_wp, -h(na, en), h(na, na) - p, q, h(na, na), h(na, en))
             else
                h(na, na) = q/h(en, na)
                h(na, en) = -(h(en, en) - p)/h(en, na)
             endif
             h(en, na) = 0.0_wp
             h(en, en) = 1.0_wp
             enm2 = na - 1
             if (enm2 /= 0) then
                ! FOR I=EN-2 STEP -1 UNTIL 1 DO --
                do ii = 1, enm2
                   i = na - ii
                   w = h(i, i) - p
                   ra = 0.0_wp
                   sa = h(i, en)
                   !
                   do j = m, na
                      ra = ra + h(i, j)*h(j, na)
                      sa = sa + h(i, j)*h(j, en)
                   enddo
                   !
                   if (Wi(i) >= 0.0_wp) then
                      m = i
                      if (Wi(i) == 0.0_wp) then
                         call cdiv(-ra, -sa, w, q, h(i, na), h(i, en))
                      else
                         ! SOLVE COMPLEX EQUATIONS
                         x = h(i, i + 1)
                         y = h(i + 1, i)
                         vr = (Wr(i) - p)*(Wr(i) - p) + Wi(i)*Wi(i) - q*q
                         vi = (Wr(i) - p)*2.0_wp*q
                         if (vr == 0.0_wp .and. vi == 0.0_wp) then
                            s1 = norm*(abs(w) + abs(q) + abs(x) + abs(y) + abs(zz))
                            vr = s1
                            do while (.true.)
                               vr = 0.5_wp*vr
                               if (s1 + vr <= s1) exit
                            enddo
                            vr = 2.0_wp*vr
                         endif
                         call cdiv(x*r - zz*ra + q*sa, x*s - zz*sa - q*ra, vr, vi, h(i, na), &
                                   h(i, en))
                         if (abs(x) <= abs(zz) + abs(q)) then
                            call cdiv(-r - y*h(i, na), -s - y*h(i, en), zz, q, h(i + 1, na), &
                                      h(i + 1, en))
                         else
                            h(i + 1, na) = (-ra - w*h(i, na) + q*h(i, en))/x
                            h(i + 1, en) = (-sa - w*h(i, en) - q*h(i, na))/x
                         endif
                      endif
                   else
                      zz = w
                      r = ra
                      s = sa
                   endif
                enddo
             endif
          elseif (q == 0) then
             ! REAL VECTOR
             m = en
             h(en, en) = 1.0_wp
             if (na /= 0) then
                ! FOR I=EN-1 STEP -1 UNTIL 1 DO --
                do ii = 1, na
                   i = en - ii
                   w = h(i, i) - p
                   r = h(i, en)
                   if (m <= na) then
                      do j = m, na
                         r = r + h(i, j)*h(j, en)
                      enddo
                   endif
                   if (Wi(i) >= 0.0_wp) then
                      ! END REAL VECTOR
                      m = i
                      if (Wi(i) == 0.0_wp) then
                         t = w
                         if (t == 0.0_wp) then
                            t = norm
                            do while (.true.)
                               t = 0.5_wp*t
                               if (norm + t <= norm) exit
                            enddo
                            t = 2.0_wp*t
                         endif
                         h(i, en) = -r/t
                      else
                         ! SOLVE REAL EQUATIONS
                         x = h(i, i + 1)
                         y = h(i + 1, i)
                         q = (Wr(i) - p)*(Wr(i) - p) + Wi(i)*Wi(i)
                         t = (x*s - zz*r)/q
                         h(i, en) = t
                         if (abs(x) <= abs(zz)) then
                            h(i + 1, en) = (-s - y*t)/zz
                         else
                            h(i + 1, en) = (-r - w*t)/x
                         endif
                      endif
                   else
                      zz = w
                      s = r
                   endif
                enddo
             endif
          endif
       enddo
       ! END BACK SUBSTITUTION.
       ! VECTORS OF ISOLATED ROOTS
       do i = 1, n
          if (i < Low .or. i > Igh) then
             do j = i, n
                z(i, j) = h(i, j)
             enddo
          endif
       enddo
       ! MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
       ! VECTORS OF ORIGINAL FULL MATRIX.
       ! FOR J=N STEP -1 UNTIL LOW DO --
       do jj = Low, n
          j = n + Low - jj
          m = min(j, Igh)
          do i = Low, Igh
             zz = 0.0_wp
             do k = Low, m
                zz = zz + z(i, k)*h(k, j)
             enddo
             z(i, j) = zz
          enddo
       enddo
    endif
 end subroutine hqr2
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute the eigenvalues and, optionally, the eigenvectors
!  of a real general matrix.
!
!  This subroutine calls the recommended sequence of
!  subroutines from the eigensystem subroutine package (EISPACK)
!  To find the eigenvalues and eigenvectors (if desired)
!  of a REAL GENERAL matrix.
!
!### References
!  * B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!    Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!    system Routines - EISPACK Guide, Springer-Verlag,
!    1976.
!
!### Author
!  * Smith, B. T., et al.
!
!### History  (YYMMDD)
!  * 760101  DATE WRITTEN
!  * 890831  Modified array declarations.  (WRB)
!  * 890831  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 920501  Reformatted the REFERENCES section.  (WRB)
!  * 921103  Corrected description of IV1.  (DWL, FNF and WRB)
!  * Jacob Williams, refactored into modern Fortran (3/25/2018)

 subroutine rg(Nm, n, a, Wr, Wi, Matz, z, Iv1, Fv1, Ierr)

    implicit none

    integer,intent(in)  :: n   !! the order of the matrix A.
                               !! N must be less than or equal to NM.
    integer,intent(in)  :: Nm  !! must be set to the row dimension of the two-dimensional
                               !! array parameters, A and Z, as declared in the calling
                               !! program dimension statement.
    integer,intent(in)  :: Matz !! an INTEGER variable set equal to zero if only
                                !! eigenvalues are desired.  Otherwise, it is set to any
                                !! non-zero integer for both eigenvalues and eigenvectors.
    real(wp),intent(inout) :: a(Nm, *)   !! contains the real general matrix.
                                         !! dimensioned A(NM,N).
                                         !! Note: A is destroyed on output.
    integer,intent(out)  :: Ierr  !! an INTEGER flag set to:
                                  !!
                                  !! * 0 -- for normal return,
                                  !! * 10*N -- if N is greater than NM,
                                  !! * J    -- if the J-th eigenvalue has not been
                                  !!           determined after a total of 30 iterations.
                                  !!           The eigenvalues should be correct for indices
                                  !!           IERR+1, IERR+2, ..., N, but no eigenvectors are
                                  !!           computed.
    real(wp),intent(out) :: Wr(*)  !! real part of the eigenvalues.  The eigenvalues are unordered except
                                   !! that complex conjugate pairs of eigenvalues appear consecutively
                                   !! with the eigenvalue having the positive imaginary part
                                   !! first.  If an error exit is made, the eigenvalues should be
                                   !! correct for indices IERR+1, IERR+2, ..., N.  WR and WI are
                                   !! one-dimensional REAL arrays, dimensioned WR(N) and WI(N).
    real(wp),intent(out) :: Wi(*)  !! imaginary part of the eigenvalues.
    real(wp),intent(out) :: z(Nm, *) !! contains the real and imaginary parts of the eigenvectors
                                     !! if MATZ is not zero.  If the J-th eigenvalue is real, the
                                     !! J-th column of Z contains its eigenvector.  If the J-th
                                     !! eigenvalue is complex with positive imaginary part, the
                                     !! J-th and (J+1)-th columns of Z contain the real and
                                     !! imaginary parts of its eigenvector.  The conjugate of this
                                     !! vector is the eigenvector for the conjugate eigenvalue.
                                     !! Z is a two-dimensional REAL array, dimensioned Z(NM,N).
    real(wp),intent(inout) :: Fv1(*) !! one-dimensional temporary storage arrays of dimension N.
    integer,intent(inout)  :: Iv1(*) !! one-dimensional temporary storage arrays of dimension N.

    integer  :: is1
    integer  :: is2

    if (n <= Nm) then
       call balanc(Nm, n, a, is1, is2, Fv1)
       call elmhes(Nm, n, is1, is2, a, Iv1)
       if (Matz /= 0) then
          ! find both eigenvalues and eigenvectors
          call eltran(Nm, n, is1, is2, a, Iv1, z)
          call hqr2(Nm, n, is1, is2, a, Wr, Wi, z, Ierr)
          if (Ierr == 0) call balbak(Nm, n, is1, is2, Fv1, n, z)
       else
          ! find eigenvalues only
          call hqr(Nm, n, is1, is2, a, Wr, Wi, Ierr)
       endif
    else
       Ierr = 10*n
    endif

 end subroutine rg
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute the eigenvalues and, optionally, the eigenvectors
!  of a real general matrix.
!
!### See also
!  * See [[rg]] for more details. This routine is just a wrapper to that one.
!
!### Author
!  * Jacob Williams, 3/25/2018

    subroutine compute_eigenvalues_and_eigenvectors(n, a, w, z, ierr)

    use numbers_module

    implicit none

    integer,intent(in)                  :: n    !! the order of the matrix `a`
    real(wp),dimension(n,n),intent(in)  :: a    !! contains the real general matrix
    real(wp),dimension(n,2),intent(out) :: w    !! real and imaginary parts of the eigenvalues
    real(wp),dimension(n,n),intent(out) :: z    !! real and imaginary parts of the eigenvectors
    integer,intent(out)                 :: ierr !! output flag from [[rg]]

    integer,parameter :: matz = 1 !! tells [[rg]] to compute eigenvalues and eigenvectors

    integer                 :: i     !! counter
    real(wp),dimension(n,n) :: a_tmp !! copy of [[a]] matrix
    real(wp),dimension(n)   :: fv1   !! work array for [[rg]]
    integer,dimension(n)    :: iv1   !! work array for [[rg]]
    real(wp),dimension(n)   :: wr    !! real part of the eigenvalues
    real(wp),dimension(n)   :: wi    !! imaginary part of the eigenvalues

    ! temp arrays:
    a_tmp = a
    wr = zero
    wi = zero

    ! call the general routine:
    call rg(n, n, a_tmp, wr, wi, matz, z, iv1, fv1, ierr)

    ! pack outputs:
    do i=1,n
        w(i,1) = wr(i)
        w(i,2) = wi(i)
    end do

    end subroutine compute_eigenvalues_and_eigenvectors
!*****************************************************************************************

!*****************************************************************************************
!>
!  Returns only the real eigenvalues and the associated eigenvectors.
!  Wrapper for [[compute_eigenvalues_and_eigenvectors]].

    subroutine compute_real_eigenvalues_and_normalized_eigenvectors(n, a, e, v, n_results, ierr)

    use numbers_module

    implicit none

    integer,intent(in)                              :: n         !! the order of the matrix `a`
    real(wp),dimension(n,n),intent(in)              :: a         !! contains the real general matrix
    real(wp),dimension(:),allocatable,intent(out)   :: e         !! eigenvalues (size `n_results`)
    real(wp),dimension(:,:),allocatable,intent(out) :: v         !! eigenvectors (size `n,n_results`)
    integer,intent(out)                             :: n_results !! number of real eigenvalues
    integer,intent(out)                             :: ierr      !! output flag from [[rg]]

    real(wp),dimension(n,2) :: w  !! real and imaginary parts of the eigenvalues
    real(wp),dimension(n,n) :: z  !! real and imaginary parts of the eigenvectors
    integer :: i !! counter
    integer :: j !! counter

    call compute_eigenvalues_and_eigenvectors(n, a, w, z, ierr)

    if (ierr==0) then

        n_results = count(w(:,2)==0.0_wp)
        if (n_results>0) then
            allocate(e(n_results))
            allocate(v(n,n_results))
            j = 0
            do i = 1, n
                if (w(i,2)==0.0_wp) then ! real eigenvalue
                    j = j + 1
                    e(j) = w(i,1)
                    v(:,j) = z(:,i) / norm2(z(:,i))  ! normalized eigenvector
                end if
            end do
        end if

    else
        n_results = 0
    end if

    end subroutine compute_real_eigenvalues_and_normalized_eigenvectors
!*****************************************************************************************

!*****************************************************************************************
!>
!  Unit test

    subroutine eispack_test()

    implicit none

    real(wp),dimension(3,3),parameter :: a = reshape([1.0_wp,4.0_wp,-3.0_wp,&
                                                      2.0_wp,3.0_wp,-8.0_wp,&
                                                      3.0_wp,2.0_wp,1.001_wp], [3,3])

    real(wp),dimension(3,2) :: w    !! real and imaginary parts of the eigenvalues
    real(wp),dimension(3,3) :: z    !! real and imaginary parts of the eigenvectors
    integer :: ierr !! output flag
    integer :: i !! counter
    integer :: j !! counter
    complex(wp),dimension(3) :: v

    call compute_eigenvalues_and_eigenvectors(3, a, w, z, ierr)

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' eispack_test'
    write(*,*) '---------------'
    write(*,*) ''

    write(*,*) ''
    write(*,*) 'ierr = ', ierr
    write(*,*) ''
    write(*,*) 'eigenvalues:'
    do i = 1, 3
    write(*,*) w(i,:)
    end do
    write(*,*) ''
    write(*,*) 'eigenvectors (normalized):'

    do i = 1, 3

        if (w(i,2)==0.0_wp) then
            ! If the J-th eigenvalue is real, the
            ! J-th column of Z contains its eigenvector
            do j = 1, 3
                v(j) = cmplx(z(j,i), 0.0_wp, wp)
            end do
        elseif (w(i,2)>0.0_wp) then
            ! If the J-th eigenvalue is complex with positive imaginary part, the
            ! J-th and (J+1)-th columns of Z contain the real and
            ! imaginary parts of its eigenvector.
            do j = 1, 3
                v(j) = cmplx(z(j,i), z(j,i+1), wp)
            end do
        else
            do j = 1, 3
                v(j) = cmplx(z(j,i-1), -z(j,i), wp)
            end do
        end if
        v = v / sqrt(dot_product(v,v))
        do j = 1, 3
            write(*,'(F16.6,F16.6)') v(j)%re, v(j)%im
        end do
        write(*,*) ''

    end do

   ! ... results:
   ! eigenvalues:
   ! -1.89041207397761       0.000000000000000E+000
   ! 3.44570603698881        5.01584673789593
   ! 3.44570603698881       -5.01584673789593
   !
   ! eigenvectors:
   ! 0.650411095383963      -0.379058329718174      -0.373946474570154
   ! 0.605673686824173       0.640277015571315      -7.629236177444867E-002
   ! 8.565310483790509E-002 -0.395692850335186        1.34627813418190

   ! ... from numpy:
   ! eigenvalues:
   ! array([-1.89041207+0.j        ,  3.44570604+5.01584674j, 3.44570604-5.01584674j])
   !
   ! eigenvectors:
   ! array([[ 0.77377504  ,  0.03085326-0.3669729j  ,  0.03085326+0.3669729j ],
   !        [-0.4509546   , -0.25965047-0.37137631j , -0.25965047+0.37137631j],
   !        [-0.44487317  ,  0.81181293             ,  0.81181293            ] ])

   end subroutine eispack_test
!*****************************************************************************************

!*****************************************************************************************
  end module eispack_module
!*****************************************************************************************
