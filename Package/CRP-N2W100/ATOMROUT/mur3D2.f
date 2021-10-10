C************************ M U R 3 D 2 **********************************
      SUBROUTINE MUR3D2(x,y,z,v,dvdx,dvdy,dvdz)
      IMPLICIT NONE
       
c     Calculates 1D potential and derivatives with respect to x,y,z.
C     A sum is carried out over all atoms at a
C     distance from the elementary quadrangle smaller than a
C     prescribed value (see program SNEARN3D). The crystal coordinates
C     for these atoms are given in IVEC. The first index of IVEC
C     is the atom index and the second is 1 for X, 2 for Y and 3 for Z.      
C     The actual number of nearest neighbours is NVEC.
C     The 1D potential and its derivative are given through VREP.

c     ATTENTION: ne marche que pour 4 couches (sinon il faut modifier
C     la contribution de la relaxation des couches).

C     Version 1.0  H.F. Busnengo
C     Version 2.1  A. Salin                 27/01/03       

      INTEGER nvec
      PARAMETER(nvec=37)

      INTEGER idelta
      DOUBLE PRECISION delta
      COMMON/SURF/delta,idelta

      REAL*8 relax,relaxi
      PARAMETER(relax=0.4529913D0)
C     Relaxi: relaxation des deux couche en subsurface
      PARAMETER(relaxi=0.5062107D0)

      INTEGER i,ivec(nvec,3)
      REAL*8 x,y,z,v,dvdx,dvdy,dvdz,yij,xij,r,vr,dvdr,
     #       dvdrsr,zij

       DATA ivec/0, -2, -2, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0,
     1     0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2,
     2     2, 2, 2,
     3     0, -1, 0, 1, -2, -1, -1, 0, 0, 1, 1, 2, 2, -2, -1, -1, 0, 1,
     4     1, 2, 2, -2, -1, -1, 0, 0, 1, 1, 2, 2, -1, -1, 0, 0, 1, 1,2,
     5     0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0,
     6     1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0/
 
C
      INTEGER idrep,nrep
      PARAMETER(idrep=9999)
      REAL*8 vasint,vrep,cvrep,zrep
      COMMON/POTREP/zrep(idrep),vrep(idrep),cvrep(idrep),vasint,nrep
C
C     -----------------
      v=0.D0
      dvdx=0.D0
      dvdy=0.D0
      dvdz=0.D0
    
      DO i=1,nvec
            xij=ivec(i,1)*delta+ivec(i,3)*delta/2.D0
            yij=ivec(i,2)*delta+ivec(i,3)*delta/2.D0
         IF(ivec(i,3).EQ.0) THEN
	   zij=0.D0
	 ELSE
           zij=( (ivec(i,3)-1)*relaxi +relax)*delta
	 END IF

        r=DSQRT( (x-xij)**2+(y-yij)**2+(z+zij)**2 )
C     -----------------------------------------
C     Calculating the 1D repulsive potential

      IF(r.GT.zrep(nrep)) THEN 
         vr=0.D0 
         dvdr=0.D0 
      ELSE 
         IF (r.LT.zrep(1)) THEN
            WRITE(*,*) 'Error en MUR3D: R too small'
            STOP
         ELSE
            CALL SPLINT(zrep,vrep,cvrep,nrep,r,vr,dvdr)
         ENDIF
      ENDIF
C     -------------------------------------------

        v=v+vr
        IF (r.GT.0.D0) THEN
          dvdrsr=dvdr/r
          dvdx=dvdx+dvdrsr*(x-xij)
          dvdy=dvdy+dvdrsr*(y-yij)
          dvdz=dvdz+dvdrsr*(z+zij)
        ELSE
            dvdz=dvdz+dvdr
        END IF
      END DO

      END
