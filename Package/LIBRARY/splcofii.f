C************************ S P L C O F I I ****************************

      SUBROUTINE SPLCOFII(x,y,x0,y0,ia,ja,idx,idy)

C     Calculates the coefficients stored in COMMON/UBI/ that are used 
C     by SPLINII to interpolate over x,y.
C     x,y are the coordinates of the nodes and x0,y0 the point for
C     which an interpolation will be done. 

C     This routine assumes that x0 and y0 are in the correct range.
c     If not, no message is issued and the results are wrong.
C     IA and JA indicate the position of X0 and y0 in the grid. They
C     give the index corresponding to the left/bottom corner of the 
C     rectangle in which the point (x0,y0) is located.

      IMPLICIT NONE  

      INTEGER ia,ja,idx,idy
      DOUBLE PRECISION x(idx),y(idy),x0,y0,dxia2s6,dyja2s6

      DOUBLE PRECISION ax0,bx0,cx0,dx0,ay0,by0,cy0,dy0,ax0p,bx0p,cx0p,
     1   dx0p,ay0p,by0p,cy0p,dy0p
      COMMON/UBI/ax0,bx0,cx0,dx0,ay0,by0,cy0,dy0,ax0p,bx0p,cx0p,dx0p,
     1           ay0p,by0p,cy0p,dy0p

      dxia2s6=(x(ia+1)-x(ia))**2/6.D0
      dyja2s6=(y(ja+1)-y(ja))**2/6.D0
C                                         Coefficients for the function
      ax0=(x(ia+1)-x0)/(x(ia+1)-x(ia))  
      bx0=1.d0-ax0  
      cx0=(ax0**3-ax0)*dxia2s6
      dx0=(bx0**3-bx0)*dxia2s6

      ay0=(y(ja+1)-y0)/(y(ja+1)-y(ja))  
      by0=1.d0-ay0  
      cy0=(ay0**3-ay0)*dyja2s6
      dy0=(by0**3-by0)*dyja2s6
C                                 Coefficients for derivative w.r. to x
      ax0p=-1.d0/(x(ia+1)-x(ia))  
      bx0p=-ax0p  
      cx0p=(3.d0*ax0**2*ax0p-ax0p)*dxia2s6
      dx0p=(3.d0*bx0**2*bx0p-bx0p)*dxia2s6
C                                 Coefficients for derivative w.r. to y
      ay0p=-1.d0/(y(ja+1)-y(ja))  
      by0p=-ay0p  
      cy0p=(3.d0*ay0**2*ay0p-ay0p)*dyja2s6
      dy0p=(3.d0*by0**2*by0p-by0p)*dyja2s6

      END
