C**************************** D P L P *******************************
      FUNCTION DPLP(n,x,y,m,t)
C
C        Version 3.0  - 17/1/97       Author: A. SALIN
C
C        Calculation of the derivative of a function by interpolation.
C        Double precision. See comments in SUBROUTINE DSPLIN.
C        If t is outside the interval over which the function is
C        defined, the derivative is extrapolated linearly from its
C        value for the first two (last two) points of the interval.
C
      IMPLICIT REAL*8 (a-h,o-z)
c
      DIMENSION x(n),y(n)
      REAL*8 m(n)
      LOGICAL order
c
      COMMON/KODSPL/kod,k,iw,iligne
c
    1 FORMAT(' DPLP: extrapolation - ',1PE15.8,' <',E15.8)
    2 FORMAT(' DPLP: extrapolation - ',1PE15.8,' >',E15.8)
c
      order=x(2).GT.x(1)
      IF(k.LE.1.OR.k.GT.n) THEN
        klo=0
        khi=n+1
        GO TO 100
      END IF
C
      inc=1
      klo=k-1
      IF(t.GT.x(klo).EQV.order) THEN
 10     khi=klo+inc
        IF(khi.GT.n) THEN
          khi=n+1
        ELSE IF(t.GT.x(khi).EQV.order) THEN
          klo=khi
          inc=inc+inc
          GO TO 10
        END IF
      ELSE
        khi=klo
 20     klo=khi-inc
        IF(klo.LT.1) THEN
          klo=0
        ELSE IF(t.LT.x(klo).EQV.order) THEN
          khi=klo
          inc=inc+inc
          GO TO 20
        END IF
      END IF
C
 100  CONTINUE
      IF(khi-klo.GT.1) THEN
        km=(khi+klo)/2
        IF(t.GT.x(km).EQV.order) THEN
          klo=km
        ELSE
          khi=km
        END IF
        GO TO 100
      END IF
        k=klo+1
C
      IF(t.EQ.x(1)) k=2
      IF(k.LE.1) THEN
        e=x(2)-x(1)
        g=t-x(1)
        f=x(2)-t
        DPLP=(y(2)-y(1))/e+(m(2)*(2*g-f)+m(1)*(g-2*f))/6.D0
        kod=1
        IF(iw.GT.0.AND.iligne.GT.0) THEN
          IF(order) THEN
            WRITE(iw,1) t,x(1)
          ELSE
            WRITE(iw,2) t,x(1)
          END IF
          iligne=iligne-1
        END IF
      ELSE IF(k.GT.n) THEN
        e=x(n)-x(n-1)
        g=t-x(n-1)
        f=x(n)-t
        DPLP=(y(n)-y(n-1))/e+(m(n)*(2*g-f)+m(n-1)*(g-2*f))/6.D0
        kod=2
        IF(iw.GT.0.AND.iligne.GT.0) THEN
          IF(order) THEN
            WRITE(iw,2) t,x(n)
          ELSE
            WRITE(iw,1) t,x(n)
          END IF
          iligne=iligne-1
        END IF
      ELSE
        f=x(k)-t
        g=t-x(k-1)
        e=f+g
        e2=e*e
        DPLP=(-m(k-1)*(3.D0*f*f-e2)+m(k)*(3.D0*g*g-e2)
     1       +6.D0*(y(k)-y(k-1)))/(6.D0*e)
        kod=0
      END IF
      END
