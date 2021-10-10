C************************  O R I G E N 1 0 0 *****************************
      SUBROUTINE ORIGEN100(x,y,delta,xt,yt)

C     100 surfaces.
C     Shifts the coordinate to the inside of the quadrangle with 
C     bottom left corner at the origin. Input and output are 
C     cartesian coordinates.
C     X,Y are the initial coordinates. DELTA is the nearest neighbour
C     distance on the surface and XT,YT are the calculated shifted
C     values.
C             Version 1.1              01/02/02 - 11/02/05

      IMPLICIT NONE

      INTEGER i,ind(2)
      DOUBLE PRECISION x,y,delta,xt,yt,ori(2)

      ori(1)=x/delta
      ori(2)=y/delta
      DO i=1,2
        IF (ori(i).GE.0.D0) THEN
          ind(i)=ori(i)
        ELSE
          ind(i)=ori(i)-1
        END IF
      END DO
      xt=x-ind(1)*delta
      yt=y-ind(2)*delta

      END
