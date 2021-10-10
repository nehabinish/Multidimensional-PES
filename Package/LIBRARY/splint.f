C************************ S P L I N T ********************************
      SUBROUTINE SPLINT(xa,ya,cya,n,x,y,yp)  
      IMPLICIT NONE  

C     Calculates a function Y and its derivative YP at X by spline
C     interpolation. The N data are XA (nodes) and YA (function at
C     nodes). CYA is calculated by DSPLIN

      INTEGER klo,khi,k,n
      DOUBLE PRECISION xa(n),ya(n),cya(n),x,y,yp,h,a,b

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

      h=xa(khi)-xa(klo)  
      a=(xa(khi)-x)/h  
      b=(x-xa(klo))/h  
      y=a*ya(klo)+b*ya(khi)+  
     1      ((a**3-a)*cya(klo)+(b**3-b)*cya(khi))*(h*h)/6.D0  
      yp=(ya(khi)-ya(klo))/h+
     1    h*(cya(khi)*(3.D0*b*b-1.D0)-cya(klo)*(3.D0*a*a-1.D0))/6.D0

      END  
