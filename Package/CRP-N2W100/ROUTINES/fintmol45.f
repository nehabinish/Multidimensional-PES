C************************ F I N T M O L 4 5 ****************************

C     Calculates the 6D interpolation function for N2/W(100) using
C     17 configurations, except for the top-vertical configuration.
C     The z and r values must be the same for all configurations.

C     Version 1.1    29/10/03     Author: A. Salin

      IMPLICIT NONE

      LOGICAL fin
      INTEGER i,j,nr,nz,nr0,nz0
      DOUBLE PRECISION sq2,xa,ya,za,xb,yb,zb,z,r,v,va,vb,dvdx,dvdy,dvdz,
     1                 z0(999),r0(999)

      INTEGER idelta
      DOUBLE PRECISION delta
      COMMON/SURF/delta,idelta

      sq2=DSQRT(2.D0)
      fin=.false.

      CALL LECTATOM
C---------------------------------------------------------------------
C     TOP

C     TOP-vertical is calculated elsewhere and stored in ftop-v.dat 

C                                                     theta=pi/2, phi=0
      OPEN(10,FILE='t-pl.dat',STATUS='OLD')
      OPEN(2,FILE='ftop-pl.dat',STATUS='unknown')
      READ(10,*) nz,nr
      WRITE(2,*) nz,nr
      nr0=nr
      nz0=nz      
      DO i=1,nz
        DO j=1,nr
          READ(10,*) z,r,v
	  z0(i)=z
	  r0(j)=r
          xa=-r/2.D0
          xb=r/2.D0
          ya=0.D0
	  yb=ya
	  za=z
	  zb=z
          CALL AT3D2(xa,ya,za,va,dvdx,dvdy,dvdz,fin)
          CALL AT3D2(xb,yb,zb,vb,dvdx,dvdy,dvdz,fin)
          IF(fin) THEN
	    WRITE(*,*) 'Error'
	    STOP
	  END IF
          v=v-va-vb
	  WRITE(2,*) z,r,v
         END DO
      END DO
      CLOSE(10)
      CLOSE(2)
C                                                   theta=pi/2,phi=pi/4
      OPEN(10,FILE='t-pd.dat',STATUS='OLD')
      OPEN(2,FILE='ftop-pd.dat',STATUS='unknown')
      READ(10,*) nz,nr
      WRITE(2,*) nz,nr
      IF(nr.NE.nr0.OR.nz.NE.nz0) THEN
        WRITE(*,*) 'ERROR'
	STOP
      END IF	
      DO i=1,nz
        DO j=1,nr
          READ(10,*) z,r,v
          IF(r.NE.r0(j).OR.z.NE.z0(i)) THEN
            WRITE(*,*) 'ERROR'
	    STOP
          END IF	
          xa=-r/(2.D0*sq2)
          xb=r/(2.D0*sq2)
	  ya=xa
	  yb=xb
	  za=z
	  zb=z
          CALL AT3D2(xa,ya,za,va,dvdx,dvdy,dvdz,fin)
          CALL AT3D2(xb,yb,zb,vb,dvdx,dvdy,dvdz,fin)
          IF(fin) THEN
	    WRITE(*,*) 'Error'
	    STOP
	  END IF
          v=v-va-vb
	  WRITE(2,*) z,r,v
         END DO
      END DO
      CLOSE(10)
      CLOSE(2)
C                                                      theta=pi/4,phi=0
      OPEN(10,FILE='t-pl45.dat',STATUS='OLD')
      OPEN(2,FILE='ftop-pl45.dat',STATUS='unknown')
      READ(10,*) nz,nr
      WRITE(2,*) nz,nr
      IF(nr.NE.nr0.OR.nz.NE.nz0) THEN
        WRITE(*,*) 'ERROR'
	STOP
      END IF	
      DO i=1,nz
        DO j=1,nr
          READ(10,*) z,r,v
          IF(r.NE.r0(j).OR.z.NE.z0(i)) THEN
            WRITE(*,*) 'ERROR'
	    STOP
          END IF	
          xa=-r/(2.D0*sq2)
          xb=r/(2.D0*sq2)
          ya=0.D0
          yb=0.D0
	  za=z+xa
	  zb=z+xb
          CALL AT3D2(xa,ya,za,va,dvdx,dvdy,dvdz,fin)
          CALL AT3D2(xb,yb,zb,vb,dvdx,dvdy,dvdz,fin)
          IF(fin) THEN
	    WRITE(*,*) 'Error'
	    STOP
	  END IF
          v=v-va-vb
	  WRITE(2,*) z,r,v
         END DO
      END DO
      CLOSE(10)
      CLOSE(2)
C                                                   theta=pi/4,phi=pi/4
      OPEN(10,FILE='t-pd45.dat',STATUS='OLD')
      OPEN(2,FILE='ftop-pd45.dat',STATUS='unknown')
      READ(10,*) nz,nr
      WRITE(2,*) nz,nr
      IF(nr.NE.nr0.OR.nz.NE.nz0) THEN
        WRITE(*,*) 'ERROR'
	STOP
      END IF	
      DO i=1,nz
        DO j=1,nr
          READ(10,*) z,r,v
          IF(r.NE.r0(j).OR.z.NE.z0(i)) THEN
            WRITE(*,*) 'ERROR'
	    STOP
          END IF	
          xa=-r/(4.D0)
          xb=r/(4.D0)
          ya=xa
	  yb=xb
	  za=z-r/(2.D0*sq2)
	  zb=z+r/(2.D0*sq2)
          CALL AT3D2(xa,ya,za,va,dvdx,dvdy,dvdz,fin)
          CALL AT3D2(xb,yb,zb,vb,dvdx,dvdy,dvdz,fin)
          IF(fin) THEN
	    WRITE(*,*) 'Error'
	    STOP
	  END IF
          v=v-va-vb
	  WRITE(2,*) z,r,v
         END DO
      END DO
      CLOSE(10)
      CLOSE(2)
C---------------------------------------------------------------------
C     BRIDGE

C                                                              theta=0
      OPEN(10,FILE='b-v.dat',STATUS='OLD')
      OPEN(2,FILE='fbri-v.dat',STATUS='unknown')
      READ(10,*) nz,nr
      WRITE(2,*) nz,nr
      IF(nr.NE.nr0.OR.nz.NE.nz0) THEN
        WRITE(*,*) 'ERROR'
	STOP
      END IF	
      DO i=1,nz
        DO j=1,nr
          READ(10,*) z,r,v
          IF(r.NE.r0(j).OR.z.NE.z0(i)) THEN
            WRITE(*,*) 'ERROR'
	    STOP
          END IF	
          xa=0.5D0*delta
	  xb=xa
          ya=0.D0
	  yb=ya
          za=z-r/2.D0
          zb=z+r/2.D0
          CALL AT3D2(xa,ya,za,va,dvdx,dvdy,dvdz,fin)
          CALL AT3D2(xb,yb,zb,vb,dvdx,dvdy,dvdz,fin)
          IF(fin) THEN
	    WRITE(*,*) 'Error'
	    STOP
	  END IF
          v=v-va-vb
	  WRITE(2,*) z,r,v
         END DO
      END DO
      CLOSE(10)
      CLOSE(2)
C                                                      theta=pi/2,phi=0
      OPEN(10,FILE='b-pl.dat',STATUS='OLD')
      OPEN(2,FILE='fbri-pl.dat',STATUS='unknown')
      READ(10,*) nz,nr
      WRITE(2,*) nz,nr
      IF(nr.NE.nr0.OR.nz.NE.nz0) THEN
        WRITE(*,*) 'ERROR'
	STOP
      END IF	
      DO i=1,nz
        DO j=1,nr
          READ(10,*) z,r,v
          IF(r.NE.r0(j).OR.z.NE.z0(i)) THEN
            WRITE(*,*) 'ERROR'
	    STOP
          END IF	
          xa=0.5D0*delta-r/2.D0
          xb=0.5D0*delta+r/2.D0
          ya=0.D0
	  yb=0.D0
	  za=z
	  zb=z
          CALL AT3D2(xa,ya,za,va,dvdx,dvdy,dvdz,fin)
          CALL AT3D2(xb,yb,zb,vb,dvdx,dvdy,dvdz,fin)
          IF(fin) THEN
	    WRITE(*,*) 'Error'
	    STOP
	  END IF
          v=v-va-vb
	  WRITE(2,*) z,r,v
         END DO
      END DO
      CLOSE(10)
      CLOSE(2)
C                                                   theta=pi/2,phi=pi/2
      OPEN(10,FILE='b-pp.dat',STATUS='OLD')
      OPEN(2,FILE='fbri-pp.dat',STATUS='unknown')
      READ(10,*) nz,nr
      WRITE(2,*) nz,nr
      IF(nr.NE.nr0.OR.nz.NE.nz0) THEN
        WRITE(*,*) 'ERROR'
	STOP
      END IF	
      DO i=1,nz
        DO j=1,nr
          READ(10,*) z,r,v
          IF(r.NE.r0(j).OR.z.NE.z0(i)) THEN
            WRITE(*,*) 'ERROR'
	    STOP
          END IF	
          xa=0.5D0*delta
	  xb=xa
          ya=-r/2.D0
          yb=r/2.D0
	  za=z
	  zb=z
          CALL AT3D2(xa,ya,za,va,dvdx,dvdy,dvdz,fin)
          CALL AT3D2(xb,yb,zb,vb,dvdx,dvdy,dvdz,fin)
          IF(fin) THEN
	    WRITE(*,*) 'Error'
	    STOP
	  END IF
          v=v-va-vb
	  WRITE(2,*) z,r,v
         END DO
      END DO
      CLOSE(10)
      CLOSE(2)
C                                                     theta=pi/4,phi=0
      
      OPEN(10,FILE='b-pl45.dat',STATUS='OLD')
      OPEN(2,FILE='fbri-pl45.dat',STATUS='unknown')
      READ(10,*) nz,nr
      WRITE(2,*) nz,nr
      IF(nr.NE.nr0.OR.nz.NE.nz0) THEN
        WRITE(*,*) 'ERROR'
	STOP
      END IF	
      DO i=1,nz
        DO j=1,nr
          READ(10,*) z,r,v
          IF(r.NE.r0(j).OR.z.NE.z0(i)) THEN
            WRITE(*,*) 'ERROR'
	    STOP
          END IF	
          xa=0.5D0*delta-r/(2.D0*sq2)
          xb=0.5D0*delta+r/(2.D0*sq2)
          ya=0.D0
	  yb=0.D0
          za=z-r/(2.D0*sq2)
          zb=z+r/(2.D0*sq2)
          CALL AT3D2(xa,ya,za,va,dvdx,dvdy,dvdz,fin)
          CALL AT3D2(xb,yb,zb,vb,dvdx,dvdy,dvdz,fin)
          IF(fin) THEN
	    WRITE(*,*) 'Error'
	    STOP
	  END IF
          v=v-va-vb
	  WRITE(2,*) z,r,v
         END DO
      END DO
      CLOSE(10)
      CLOSE(2)
C                                                   theta=pi/4,phi=pi/2
      OPEN(10,FILE='b-pp45.dat',STATUS='OLD')
      OPEN(2,FILE='fbri-pp45.dat',STATUS='unknown')
      READ(10,*) nz,nr
      WRITE(2,*) nz,nr
      IF(nr.NE.nr0.OR.nz.NE.nz0) THEN
        WRITE(*,*) 'ERROR'
	STOP
      END IF	
      DO i=1,nz
        DO j=1,nr
          READ(10,*) z,r,v
          IF(r.NE.r0(j).OR.z.NE.z0(i)) THEN
            WRITE(*,*) 'ERROR'
	    STOP
          END IF	
          xa=0.5D0*delta
	  xb=xa
          ya=-r/(2.D0*sq2)
          yb=r/(2.D0*sq2)
          za=z-r/(2.D0*sq2)
          zb=z+r/(2.D0*sq2)
          CALL AT3D2(xa,ya,za,va,dvdx,dvdy,dvdz,fin)
          CALL AT3D2(xb,yb,zb,vb,dvdx,dvdy,dvdz,fin)
          IF(fin) THEN
	    WRITE(*,*) 'Error'
	    STOP
	  END IF
          v=v-va-vb
	  WRITE(2,*) z,r,v
         END DO
      END DO
      CLOSE(10)
      CLOSE(2)
C                                                   theta=pi/2,phi=pi/4
      OPEN(10,FILE='b-pd.dat',STATUS='OLD')
      OPEN(2,FILE='fbri-pd.dat',STATUS='unknown')
      READ(10,*) nz,nr
      WRITE(2,*) nz,nr
      IF(nr.NE.nr0.OR.nz.NE.nz0) THEN
        WRITE(*,*) 'ERROR'
	STOP
      END IF	
      DO i=1,nz
        DO j=1,nr
          READ(10,*) z,r,v
          IF(r.NE.r0(j).OR.z.NE.z0(i)) THEN
            WRITE(*,*) 'ERROR'
	    STOP
          END IF	
          xa=0.5D0*delta-r/(2.D0*sq2)
          xb=0.5D0*delta+r/(2.D0*sq2)
          ya=-r/(2.D0*sq2)
          yb=r/(2.D0*sq2)
	  za=z
	  zb=z
          CALL AT3D2(xa,ya,za,va,dvdx,dvdy,dvdz,fin)
          CALL AT3D2(xb,yb,zb,vb,dvdx,dvdy,dvdz,fin)
          IF(fin) THEN
	    WRITE(*,*) 'Error'
	    STOP
	  END IF
          v=v-va-vb
	  WRITE(2,*) z,r,v
         END DO
      END DO
      CLOSE(10)
      CLOSE(2)
C                                                  theta=pi/4,phi=pi/4
      OPEN(10,FILE='b-pd45.dat',STATUS='OLD')
      OPEN(2,FILE='fbri-pd45.dat',STATUS='unknown')
      READ(10,*) nz,nr
      WRITE(2,*) nz,nr
      IF(nr.NE.nr0.OR.nz.NE.nz0) THEN
        WRITE(*,*) 'ERROR'
	STOP
      END IF	
      DO i=1,nz
        DO j=1,nr
          READ(10,*) z,r,v
          IF(r.NE.r0(j).OR.z.NE.z0(i)) THEN
            WRITE(*,*) 'ERROR'
	    STOP
          END IF	
          xa=0.5D0*delta-r/4.D0
          xb=0.5D0*delta+r/4.D0
          ya=-r/4.D0
          yb=r/4.D0
          za=z-r/(2.D0*sq2)
          zb=z+r/(2.D0*sq2)
          CALL AT3D2(xa,ya,za,va,dvdx,dvdy,dvdz,fin)
          CALL AT3D2(xb,yb,zb,vb,dvdx,dvdy,dvdz,fin)
          IF(fin) THEN
	    WRITE(*,*) 'Error'
	    STOP
	  END IF
          v=v-va-vb
	  WRITE(2,*) z,r,v
         END DO
      END DO
      CLOSE(10)
      CLOSE(2)

C---------------------------------------------------------------------
C     HOLE

C                                                              theta=0
      OPEN(10,FILE='h-v.dat',STATUS='OLD')
      OPEN(2,FILE='fhol-v.dat',STATUS='unknown')
      READ(10,*) nz,nr
      WRITE(2,*) nz,nr
      IF(nr.NE.nr0.OR.nz.NE.nz0) THEN
        WRITE(*,*) 'ERROR'
	STOP
      END IF	
      DO i=1,nz
        DO j=1,nr
          READ(10,*) z,r,v
          IF(r.NE.r0(j).OR.z.NE.z0(i)) THEN
            WRITE(*,*) 'ERROR'
	    STOP
          END IF	
          xa=0.5D0*delta
	  xb=xa
	  ya=xa
	  yb=xb
          za=z-r/2.D0
          zb=z+r/2.D0
          CALL AT3D2(xa,ya,za,va,dvdx,dvdy,dvdz,fin)
          CALL AT3D2(xb,yb,zb,vb,dvdx,dvdy,dvdz,fin)
          IF(fin) THEN
	    WRITE(*,*) 'Error'
	    STOP
	  END IF
          v=v-va-vb
	  WRITE(2,*) z,r,v
         END DO
      END DO
      CLOSE(10)
      CLOSE(2)
C                                                      theta=pi/2,phi=0
      OPEN(10,FILE='h-pl.dat',STATUS='OLD')
      OPEN(2,FILE='fhol-pl.dat',STATUS='unknown')
      READ(10,*) nz,nr
      WRITE(2,*) nz,nr
      IF(nr.NE.nr0.OR.nz.NE.nz0) THEN
        WRITE(*,*) 'ERROR'
	STOP
      END IF	
      DO i=1,nz
        DO j=1,nr
          READ(10,*) z,r,v
          IF(r.NE.r0(j).OR.z.NE.z0(i)) THEN
            WRITE(*,*) 'ERROR'
	    STOP
          END IF	
          xa=0.5D0*(delta-r)
          xb=0.5D0*(delta+r)
          ya=0.5D0*delta
	  yb=ya
	  za=z
	  zb=z
          CALL AT3D2(xa,ya,za,va,dvdx,dvdy,dvdz,fin)
          CALL AT3D2(xb,yb,zb,vb,dvdx,dvdy,dvdz,fin)
          IF(fin) THEN
	    WRITE(*,*) 'Error'
	    STOP
	  END IF
          v=v-va-vb
	  WRITE(2,*) z,r,v
         END DO
      END DO
      CLOSE(10)
      CLOSE(2)
C                                                   theta=pi/2,phi=pi/4
      OPEN(10,FILE='h-pd.dat',STATUS='OLD')
      OPEN(2,FILE='fhol-pd.dat',STATUS='unknown')
      READ(10,*) nz,nr
      WRITE(2,*) nz,nr
      IF(nr.NE.nr0.OR.nz.NE.nz0) THEN
        WRITE(*,*) 'ERROR'
	STOP
      END IF	
      DO i=1,nz
        DO j=1,nr
          READ(10,*) z,r,v
          IF(r.NE.r0(j).OR.z.NE.z0(i)) THEN
            WRITE(*,*) 'ERROR'
	    STOP
          END IF	
          xa=0.5D0*(delta-r/sq2)
          xb=0.5D0*(delta+r/sq2)
	  ya=xa
	  yb=xb
	  za=z
	  zb=z
          CALL AT3D2(xa,ya,za,va,dvdx,dvdy,dvdz,fin)
          CALL AT3D2(xb,yb,zb,vb,dvdx,dvdy,dvdz,fin)
          IF(fin) THEN
	    WRITE(*,*) 'Error'
	    STOP
	  END IF
          v=v-va-vb
	  WRITE(2,*) z,r,v
         END DO
      END DO
      CLOSE(10)
      CLOSE(2)
C                                                      theta=pi/4,phi=0
      OPEN(10,FILE='h-pl45.dat',STATUS='OLD')
      OPEN(2,FILE='fhol-pl45.dat',STATUS='unknown')
      READ(10,*) nz,nr
      WRITE(2,*) nz,nr
      IF(nr.NE.nr0.OR.nz.NE.nz0) THEN
        WRITE(*,*) 'ERROR'
	STOP
      END IF	
      DO i=1,nz
        DO j=1,nr
          READ(10,*) z,r,v
          IF(r.NE.r0(j).OR.z.NE.z0(i)) THEN
            WRITE(*,*) 'ERROR'
	    STOP
          END IF	
          xa=0.5D0*delta-r/(2.D0*sq2)
          xb=0.5D0*delta+r/(2.D0*sq2)
          ya=0.5D0*delta
	  yb=ya
          za=z-r/(2.D0*sq2)
          zb=z+r/(2.D0*sq2)
          CALL AT3D2(xa,ya,za,va,dvdx,dvdy,dvdz,fin)
          CALL AT3D2(xb,yb,zb,vb,dvdx,dvdy,dvdz,fin)
          IF(fin) THEN
	    WRITE(*,*) 'Error'
	    STOP
	  END IF
          v=v-va-vb
	  WRITE(2,*) z,r,v
         END DO
      END DO
      CLOSE(10)
      CLOSE(2)
C                                                   theta=pi/4,phi=pi/4
      OPEN(10,FILE='h-pd45.dat',STATUS='OLD')
      OPEN(2,FILE='fhol-pd45.dat',STATUS='unknown')
      READ(10,*) nz,nr
      WRITE(2,*) nz,nr
      IF(nr.NE.nr0.OR.nz.NE.nz0) THEN
        WRITE(*,*) 'ERROR'
	STOP
      END IF	
      DO i=1,nz
        DO j=1,nr
          READ(10,*) z,r,v
          IF(r.NE.r0(j).OR.z.NE.z0(i)) THEN
            WRITE(*,*) 'ERROR'
	    STOP
          END IF	
          xa=0.5D0*delta-r/4.D0
          xb=0.5D0*delta+r/4.D0
	  ya=xa
	  yb=xb
          za=z-r/(2.D0*sq2)
          zb=z+r/(2.D0*sq2)
          CALL AT3D2(xa,ya,za,va,dvdx,dvdy,dvdz,fin)
          CALL AT3D2(xb,yb,zb,vb,dvdx,dvdy,dvdz,fin)
          IF(fin) THEN
	    WRITE(*,*) 'Error'
	    STOP
	  END IF
          v=v-va-vb
	  WRITE(2,*) z,r,v
         END DO
      END DO
      CLOSE(10)
      CLOSE(2)

      END
