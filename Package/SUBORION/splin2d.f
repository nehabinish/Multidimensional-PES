C**************************** S P L I N 2 D *************************      
      SUBROUTINE SPLIN2D(x,m,y,n,f,cx,cy,cint,idx,idy)
C
C     Two dimensional interpolation by Spline method.
C
C       Author: A. Salin     Version: 1.1    15/04/98 - 18/1/99

C     X,M: set of m abcissae
C     Y,N: set of n ordinates
C     F(i,j): value of functions at (x_i,y_j)
C     CX, CY, CINT: arrays of dimension (idx,idy) determined by
C                      subroutine SPLIN2D
C     IDX, IDY: dimension parameters in calling program
C
      IMPLICIT REAL*8 (a-h,o-z)
      PARAMETER (idn=999)
C
      DIMENSION x(idx),y(idy),f(idx,idy),cx(idx,idy),cy(idx,idy),
     1          cint(idx,idy)
      DIMENSION a(idn),aa(idn),aaa(idn),cg(idn),g(idn)
C      
      IF(idn.LT.idx.OR.idn.LT.idy) THEN
        WRITE(*,*) 'SPLIN2D: idn too small'
        STOP
      END IF
C
      DO j=1,n
        CALL DSPLIN(m,x,f(1,j),cx(1,j),0.D0,0,0.D0,0,a,aa,aaa)
      END DO

      DO i=1,m
        DO j=1,n
          g(j)=f(i,j)
        END DO
        CALL DSPLIN(n,y,g,cg,0.D0,0,0.D0,0,a,aa,aaa)
        DO j=1,n
          cy(i,j)=cg(j)
        END DO
      END DO  

      DO j=1,n
        CALL DSPLIN(m,x,cy(1,j),cint(1,j),0.D0,0,0.D0,0,a,aa,aaa)
      END DO

      END
