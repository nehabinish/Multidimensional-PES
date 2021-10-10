C************************ U B I C A R ********************************
      SUBROUTINE UBICAR(xa,x,n,klo,paro)

C     Determine KLO such that X is .GE.XA(KLO) and .LT.XA(KLO+1) 
C     If X is outside the interval [XA(1),XA(N)], the logical
C     variable PARO takes the value .TRUE.

      IMPLICIT NONE

      LOGICAL paro
      INTEGER k,n,klo,khi
      DOUBLE PRECISION xa(n),x

      paro=.FALSE.
      IF (x.LT.xa(1).OR.x.GT.xa(n)) THEN
         paro=.TRUE.
         RETURN
      END IF
      klo=1
      khi=n
1     CONTINUE
      IF (khi-klo.GT.1) THEN
        k=(khi+klo)/2
        IF(xa(k).GT.x)THEN
          khi=k
        ELSE
          klo=k
        END IF
        GO TO 1

      END IF

      END
