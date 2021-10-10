C**************************** D P L C O F T ***************************
      SUBROUTINE DPLCOFT(n,x,t,c,dc,k)

C     Calculates coefficients for the spline interpolation of a
C     function (C) and its derivative (DC) at point T. 
C     The N nodes are stored in X by INCREASING order.
C     On output, K is such that t.GE.x(k-1) and t.LT.x(k).
C     Calculations are accelerated if K is on input close to the final
C     value 

      IMPLICIT NONE
c
      INTEGER n,k,klo,khi,inc,km
      DOUBLE PRECISION x(n),c(4),dc(4),t,f,g,e,se,se6,h
c
      IF(k.LE.1.OR.k.GT.n) THEN
        klo=0
        khi=n+1
      ELSE
C
      inc=1
      klo=k-1
      IF(t.GT.x(klo)) THEN
 10     khi=klo+inc
        IF(khi.GT.n) THEN
          khi=n+1
        ELSE IF(t.GT.x(khi)) THEN
          klo=khi
          inc=inc+inc
          GO TO 10
        END IF
      ELSE
        khi=klo
 20     klo=khi-inc
        IF(klo.LT.1) THEN
          klo=0
        ELSE IF(t.LT.x(klo)) THEN
          khi=klo
          inc=inc+inc
          GO TO 20
        END IF
      END IF
C
      END IF

100   CONTINUE
      IF(khi-klo.GT.1) THEN
        km=(khi+klo)/2
        IF(t.GT.x(km)) THEN
          klo=km
        ELSE
          khi=km
        END IF
        GO TO 100
      END IF
        k=klo+1
C
      IF(t.EQ.x(1)) k=2

        f=x(k)-t
        g=t-x(k-1)
        e=f+g
        se=1.D0/e
        se6=se/6.D0
        h=-g*f*se6
        c(1)=h*(f+e)
        c(2)=h*(g+e)
        c(3)=f*se
        c(4)=g*se
        dc(1)=-h-(f+e)*(f-g)*se6
        dc(2)=h-(g+e)*(f-g)*se6
        dc(3)=-se
        dc(4)=se

      END
