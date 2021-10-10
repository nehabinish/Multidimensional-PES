      BLOCK DATA BLKSPL
      COMMON/KODSPL/kod,k,iw,iligne
      DATA iw,iligne,k/6,100,2/
      END
C**************************** T R I D I A ***************************
      SUBROUTINE TRIDIA(alpha,beta,gamma,b,x,n)
C
C        Solution of tridiagonal systems
C
C     ALPHA: main diagonal (destroyed in TRIDIA).
C     BETA: subdiagonal (first element has index 1).
C     GAMMA: superdiagonal (first element has index 1).
C     B: right-hand side (destroyed in TRIDIA).
C     X: solution.
C     N: dimension of system.
C
      IMPLICIT REAL*8 (a-h,o-z)
      DIMENSION alpha(n),beta(n),gamma(n),b(n),x(n)
C
      DO i=2,n
        rap = beta(i-1) / alpha(i-1)
        alpha(i) = alpha(i) - rap*gamma(i-1)
        b(i) = b(i) - rap*b(i-1)
      END DO
      x(n) = b(n) / alpha(n)
      
      DO j=2,n
        i = n-j+1
        x(i) = ( b(i) - gamma(i)*x(i+1) ) / alpha(i)
      END DO

      END
