C************************ S P L I N I I ******************************
      SUBROUTINE SPLINII(f,d2fdx,d2fdy,d4fdydx,ia,ja,vint,idx,idy)
C     ------------------------------------------------------------------ 

C     Interpolation over a grid. SPLCOFII must be called previously to 
c     determine the coefficients stored in COMMON/UBI/.
C     On output, the vector VINT contains the value of the function,
C     the derivative with respect to y (index 2) and the derivative
C     with respect to x (index 3). 

C     IA and JA indicate the position of X0 and y0 in the grid. They
C     give the index corresponding to the left/bottom corner of the 
C     rectangle in which the point (x0,y0) is located.

      IMPLICIT NONE
        
      INTEGER ia,ja,idx,idy
      DOUBLE PRECISION f(idx,idy),vint(3)
      DOUBLE PRECISION d2fdx(idx,idy),d2fdy(idx,idy),d4fdydx(idx,idy)

      DOUBLE PRECISION ax0,bx0,cx0,dx0,ay0,by0,cy0,dy0,ax0p,bx0p,cx0p,
     1   dx0p,ay0p,by0p,cy0p,dy0p
      COMMON/UBI/ax0,bx0,cx0,dx0,ay0,by0,cy0,dy0,ax0p,bx0p,cx0p,dx0p,
     1           ay0p,by0p,cy0p,dy0p
C                                                              Function
      vint(1)= ax0*( ay0*f(ia,ja)+by0*f(ia,ja+1)
     1             + cy0*d2fdy(ia,ja)+dy0*d2fdy(ia,ja+1) )
     2       + bx0*( ay0*f(ia+1,ja)+by0*f(ia+1,ja+1)
     3             + cy0*d2fdy(ia+1,ja)+dy0*d2fdy(ia+1,ja+1) )
     4       + cx0*( ay0*d2fdx(ia,ja)+by0*d2fdx(ia,ja+1)
     5             + cy0*d4fdydx(ia,ja)+dy0*d4fdydx(ia,ja+1) )
     6       + dx0*( ay0*d2fdx(ia+1,ja)+by0*d2fdx(ia+1,ja+1)
     7             + cy0*d4fdydx(ia+1,ja)+dy0*d4fdydx(ia+1,ja+1) )
C                                                  Derivative w.r. to y     
      vint(2)= ay0p*( ax0*f(ia,ja)+bx0*f(ia+1,ja)
     1              + cx0*d2fdx(ia,ja)+dx0*d2fdx(ia+1,ja) )
     2       + by0p*( ax0*f(ia,ja+1)+bx0*f(ia+1,ja+1)
     3              + cx0*d2fdx(ia,ja+1)+dx0*d2fdx(ia+1,ja+1) )
     4       + cy0p*( ax0*d2fdy(ia,ja)+bx0*d2fdy(ia+1,ja)
     5              + cx0*d4fdydx(ia,ja)+dx0*d4fdydx(ia+1,ja) )
     6       + dy0p*( ax0*d2fdy(ia,ja+1)+bx0*d2fdy(ia+1,ja+1)
     7              + cx0*d4fdydx(ia,ja+1)+dx0*d4fdydx(ia+1,ja+1) )  
C                                                  Derivative w.r. to x
      vint(3)= ax0p*( ay0*f(ia,ja)+by0*f(ia,ja+1)
     1              + cy0*d2fdy(ia,ja)+dy0*d2fdy(ia,ja+1) )
     2       + bx0p*( ay0*f(ia+1,ja)+by0*f(ia+1,ja+1)
     3              + cy0*d2fdy(ia+1,ja)+dy0*d2fdy(ia+1,ja+1) )
     4       + cx0p*( ay0*d2fdx(ia,ja)+by0*d2fdx(ia,ja+1)
     5              + cy0*d4fdydx(ia,ja)+dy0*d4fdydx(ia,ja+1) )
     6       + dx0p*( ay0*d2fdx(ia+1,ja)+by0*d2fdx(ia+1,ja+1)  
     7              + cy0*d4fdydx(ia+1,ja)+dy0*d4fdydx(ia+1,ja+1) )

      END
