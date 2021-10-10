C**************************** D P L *********************************
      FUNCTION DPL(n,x,y,m,t)
C
C        Version 3.0  - 15/1/97       Author: A. SALIN
C
C        Calculation of a function by interpolation. Double precision.
C        See comments in SUBROUTINE DSPLIN
C
      IMPLICIT REAL*8 (a-h,o-z)
c
      DIMENSION x(n),y(n)
      REAL*8 m(n)
      LOGICAL order
c
      COMMON/KODSPL/kod,k,iw,iligne
c
    1 FORMAT(' DPL: extrapolation - ',1PE15.8,' <',E15.8)
    2 FORMAT(' DPL: extrapolation - ',1PE15.8,' >',E15.8)
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
        DPL=((y(2)-y(1))/e-m(2)*e/6.D0)*(t-x(1))+y(1)
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
        DPL=((y(n)-y(n-1))/e+m(n-1)*e/6.D0)*(t-x(n))+y(n)
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
        DPL=(-g*f*(m(k-1)*(f+e)+m(k)*(g+e))+6.D0*(g*y(k)+f*y(k-1)))/
     1      (6.D0*e)
        kod=0
      END IF
      END
