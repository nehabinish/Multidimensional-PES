C************************ P R E P M O L 4 5***************************

C     N2/W(100) -17 configurations.

C     This program reads the interpolation for the various
C     configurations, as calculated by FINTMOL45, and writes down in 
C     TOUT45.DAT the input that is used by LECMOL45 for the calculation
C     of the PES. It must be run once before starting a calculation. 
C     Once the file TOUT45.DAT has been created, there is no point to
C     run it again. 
C     Note: the values of r and Z must be the same in all data sets.

C      Version 1.0   27/07/03      A. Salin

      IMPLICIT NONE

      INTEGER idm,idc
      PARAMETER(idm=16,idc=17)
C     idm: maximum number of node points
C     idc: maximum number of configurations
      
      INTEGER i,j,nr,nz,k,nr0,nz0
      DOUBLE PRECISION z(idm),r(idm),v(idm,idm,idc),cvr(idm,idm,idc),
     1             cvz(idm,idm,idc),cvzr(idm,idm,idc),r0(idm),z0(idm)

      nr=11
      nz=16
C----------------------------------------------------------------------
C     Potential on TOP, 

C                                                               theta=0
      OPEN(10,FILE='ftop-v.dat',STATUS='OLD')
      READ(10,*) nz,nr
      nr0=nr
      nz0=nz      
      DO i=1,nz
        DO j=1,nr
          READ(10,*) z(i),r(j),v(j,i,1)
	  z0(i)=z(i)
	  r0(j)=r(j)
        END DO
      END DO
      CLOSE(10)
      CALL SPLIN2D(r,nr,z,nz,v(1,1,1),cvr(1,1,1),
     #             cvz(1,1,1),cvzr(1,1,1),idm,idm)
C                                                      theta=pi/2,phi=0
      OPEN(10,FILE='ftop-pl.dat',STATUS='OLD')
      READ(10,*) nz,nr
      IF(nr.NE.nr0.OR.nz.NE.nz0) THEN
        WRITE(*,*) 'ERROR'
	STOP
      END IF	
      DO i=1,nz
        DO j=1,nr
          READ(10,*) z(i),r(j),v(j,i,2)
          IF(r(j).NE.r0(j).OR.z(i).NE.z0(i)) THEN
            WRITE(*,*) 'ERROR'
	    STOP
          END IF	
        END DO
      END DO
      CLOSE(10)
      CALL SPLIN2D(r,nr,z,nz,v(1,1,2),cvr(1,1,2),
     #             cvz(1,1,2),cvzr(1,1,2),idm,idm)
C                                                   theta=pi/2,phi=pi/4
      OPEN(10,FILE='ftop-pd.dat',STATUS='OLD')
      READ(10,*) nz,nr
      IF(nr.NE.nr0.OR.nz.NE.nz0) THEN
        WRITE(*,*) 'ERROR'
	STOP
      END IF	
      DO i=1,nz
        DO j=1,nr
          READ(10,*) z(i),r(j),v(j,i,3)
          IF(r(j).NE.r0(j).OR.z(i).NE.z0(i)) THEN
            WRITE(*,*) 'ERROR'
	    STOP
          END IF	
        END DO
      END DO
      CLOSE(10)
      CALL SPLIN2D(r,nr,z,nz,v(1,1,3),cvr(1,1,3),
     #             cvz(1,1,3),cvzr(1,1,3),idm,idm)
C                                                  theta=pi/4 and phi=0
      OPEN(10,FILE='ftop-pl45.dat',STATUS='OLD')
      READ(10,*) nz,nr
      IF(nr.NE.nr0.OR.nz.NE.nz0) THEN
        WRITE(*,*) 'ERROR'
	STOP
      END IF	
      DO i=1,nz
        DO j=1,nr
          READ(10,*) z(i),r(j),v(j,i,10)
          IF(r(j).NE.r0(j).OR.z(i).NE.z0(i)) THEN
            WRITE(*,*) 'ERROR'
	    STOP
          END IF	
        END DO
      END DO
      CLOSE(10)
      CALL SPLIN2D(r,nr,z,nz,v(1,1,10),cvr(1,1,10),
     #             cvz(1,1,10),cvzr(1,1,10),idm,idm)
C                                               theta=pi/4 and phi=pi/4
      OPEN(10,FILE='ftop-pd45.dat',STATUS='OLD')
      READ(10,*) nz,nr
      IF(nr.NE.nr0.OR.nz.NE.nz0) THEN
        WRITE(*,*) 'ERROR'
	STOP
      END IF	
      DO i=1,nz
        DO j=1,nr
          READ(10,*) z(i),r(j),v(j,i,11)
          IF(r(j).NE.r0(j).OR.z(i).NE.z0(i)) THEN
            WRITE(*,*) 'ERROR'
	    STOP
          END IF	
        END DO
      END DO
      CLOSE(10)
      CALL SPLIN2D(r,nr,z,nz,v(1,1,11),cvr(1,1,11),
     #             cvz(1,1,11),cvzr(1,1,11),idm,idm)
C----------------------------------------------------------------------
C     Potential on BRIDGE


C                                                               theta=0
      OPEN(10,FILE='fbri-v.dat',STATUS='OLD')
      READ(10,*) nz,nr
      IF(nr.NE.nr0.OR.nz.NE.nz0) THEN
        WRITE(*,*) 'ERROR'
	STOP
      END IF	
      DO i=1,nz
        DO j=1,nr
          READ(10,*) z(i),r(j),v(j,i,4)
          IF(r(j).NE.r0(j).OR.z(i).NE.z0(i)) THEN
            WRITE(*,*) 'ERROR'
	    STOP
          END IF	
        END DO
      END DO
      CLOSE(10)
      CALL SPLIN2D(r,nr,z,nz,v(1,1,4),cvr(1,1,4),
     #             cvz(1,1,4),cvzr(1,1,4),idm,idm)
C                                                      theta=pi/2,phi=0
      OPEN(10,FILE='fbri-pl.dat',STATUS='OLD')
      READ(10,*) nz,nr
      IF(nr.NE.nr0.OR.nz.NE.nz0) THEN
        WRITE(*,*) 'ERROR'
	STOP
      END IF	
      DO i=1,nz
        DO j=1,nr
          READ(10,*) z(i),r(j),v(j,i,5)
          IF(r(j).NE.r0(j).OR.z(i).NE.z0(i)) THEN
            WRITE(*,*) 'ERROR'
	    STOP
          END IF	
        END DO
      END DO
      CLOSE(10)
      CALL SPLIN2D(r,nr,z,nz,v(1,1,5),cvr(1,1,5),
     #             cvz(1,1,5),cvzr(1,1,5),idm,idm)
C                                                   theta=pi/2,phi=pi/2
      OPEN(10,FILE='fbri-pp.dat',STATUS='OLD')
      READ(10,*) nz,nr
      IF(nr.NE.nr0.OR.nz.NE.nz0) THEN
        WRITE(*,*) 'ERROR'
	STOP
      END IF	
      DO i=1,nz
        DO j=1,nr
          READ(10,*) z(i),r(j),v(j,i,6)
          IF(r(j).NE.r0(j).OR.z(i).NE.z0(i)) THEN
            WRITE(*,*) 'ERROR'
	    STOP
          END IF	
        END DO
      END DO
      CLOSE(10)
      CALL SPLIN2D(r,nr,z,nz,v(1,1,6),cvr(1,1,6),
     #             cvz(1,1,6),cvzr(1,1,6),idm,idm)
C                                                      theta=pi/4,phi=0
      OPEN(10,FILE='fbri-pl45.dat',STATUS='OLD')
      READ(10,*) nz,nr
      IF(nr.NE.nr0.OR.nz.NE.nz0) THEN
        WRITE(*,*) 'ERROR'
	STOP
      END IF	
      DO i=1,nz
        DO j=1,nr
          READ(10,*) z(i),r(j),v(j,i,12)
          IF(r(j).NE.r0(j).OR.z(i).NE.z0(i)) THEN
            WRITE(*,*) 'ERROR'
	    STOP
          END IF	
        END DO
      END DO
      CLOSE(10)
      CALL SPLIN2D(r,nr,z,nz,v(1,1,12),cvr(1,1,12),
     #             cvz(1,1,12),cvzr(1,1,12),idm,idm)
C                                                   theta=pi/4,phi=pi/2
      OPEN(10,FILE='fbri-pp45.dat',STATUS='OLD')
      READ(10,*) nz,nr
      IF(nr.NE.nr0.OR.nz.NE.nz0) THEN
        WRITE(*,*) 'ERROR'
	STOP
      END IF	
      DO i=1,nz
        DO j=1,nr
          READ(10,*) z(i),r(j),v(j,i,13)
          IF(r(j).NE.r0(j).OR.z(i).NE.z0(i)) THEN
            WRITE(*,*) 'ERROR'
	    STOP
          END IF	
        END DO
      END DO
      CLOSE(10)
      CALL SPLIN2D(r,nr,z,nz,v(1,1,13),cvr(1,1,13),
     #             cvz(1,1,13),cvzr(1,1,13),idm,idm)
C                                                   theta=pi/2,phi=pi/4
      OPEN(10,FILE='fbri-pd.dat',STATUS='OLD')
      READ(10,*) nz,nr
      IF(nr.NE.nr0.OR.nz.NE.nz0) THEN
        WRITE(*,*) 'ERROR'
	STOP
      END IF	
      DO i=1,nz
        DO j=1,nr
          READ(10,*) z(i),r(j),v(j,i,14)
          IF(r(j).NE.r0(j).OR.z(i).NE.z0(i)) THEN
            WRITE(*,*) 'ERROR'
	    STOP
          END IF	
        END DO
      END DO
      CLOSE(10)
      CALL SPLIN2D(r,nr,z,nz,v(1,1,14),cvr(1,1,14),
     #             cvz(1,1,14),cvzr(1,1,14),idm,idm)
C                                                   theta=pi/4,phi=pi/4
      OPEN(10,FILE='fbri-pd45.dat',STATUS='OLD')
      READ(10,*) nz,nr
      IF(nr.NE.nr0.OR.nz.NE.nz0) THEN
        WRITE(*,*) 'ERROR'
	STOP
      END IF	
      DO i=1,nz
        DO j=1,nr
          READ(10,*) z(i),r(j),v(j,i,15)
          IF(r(j).NE.r0(j).OR.z(i).NE.z0(i)) THEN
            WRITE(*,*) 'ERROR'
	    STOP
          END IF	
        END DO
      END DO
      CLOSE(10)
      CALL SPLIN2D(r,nr,z,nz,v(1,1,15),cvr(1,1,15),
     #             cvz(1,1,15),cvzr(1,1,15),idm,idm)

C----------------------------------------------------------------------
C     Potential on HOLE

C                                                               theta=0
      OPEN(10,FILE='fhol-v.dat',STATUS='OLD')
      READ(10,*) nz,nr
      IF(nr.NE.nr0.OR.nz.NE.nz0) THEN
        WRITE(*,*) 'ERROR'
	STOP
      END IF	
      DO i=1,nz
        DO j=1,nr
          READ(10,*) z(i),r(j),v(j,i,7)
          IF(r(j).NE.r0(j).OR.z(i).NE.z0(i)) THEN
            WRITE(*,*) 'ERROR'
	    STOP
          END IF	
        END DO
      END DO
      CLOSE(10)
      CALL SPLIN2D(r,nr,z,nz,v(1,1,7),cvr(1,1,7),
     #             cvz(1,1,7),cvzr(1,1,7),idm,idm)
C                                                      theta=pi/2,phi=0
      OPEN(10,FILE='fhol-pl.dat',STATUS='OLD')
      READ(10,*) nz,nr
      IF(nr.NE.nr0.OR.nz.NE.nz0) THEN
        WRITE(*,*) 'ERROR'
	STOP
      END IF	
      DO i=1,nz
        DO j=1,nr
          READ(10,*) z(i),r(j),v(j,i,8)
          IF(r(j).NE.r0(j).OR.z(i).NE.z0(i)) THEN
            WRITE(*,*) 'ERROR'
	    STOP
          END IF	
        END DO
      END DO
      CLOSE(10)
      CALL SPLIN2D(r,nr,z,nz,v(1,1,8),cvr(1,1,8),
     #             cvz(1,1,8),cvzr(1,1,8),idm,idm)
C                                                   theta=pi/2,phi=pi/4
      OPEN(10,FILE='fhol-pd.dat',STATUS='OLD')
      READ(10,*) nz,nr
      IF(nr.NE.nr0.OR.nz.NE.nz0) THEN
        WRITE(*,*) 'ERROR'
	STOP
      END IF	
      DO i=1,nz
        DO j=1,nr
          READ(10,*) z(i),r(j),v(j,i,9)
          IF(r(j).NE.r0(j).OR.z(i).NE.z0(i)) THEN
            WRITE(*,*) 'ERROR'
	    STOP
          END IF	
        END DO
      END DO
      CLOSE(10)
      CALL SPLIN2D(r,nr,z,nz,v(1,1,9),cvr(1,1,9),
     #             cvz(1,1,9),cvzr(1,1,9),idm,idm)
C                                                      theta=pi/4,phi=0
      OPEN(10,FILE='fhol-pl45.dat',STATUS='OLD')
      READ(10,*) nz,nr
      IF(nr.NE.nr0.OR.nz.NE.nz0) THEN
        WRITE(*,*) 'ERROR'
	STOP
      END IF	
      DO i=1,nz
        DO j=1,nr
          READ(10,*) z(i),r(j),v(j,i,16)
          IF(r(j).NE.r0(j).OR.z(i).NE.z0(i)) THEN
            WRITE(*,*) 'ERROR'
	    STOP
          END IF	
        END DO
      END DO
      CLOSE(10)
      CALL SPLIN2D(r,nr,z,nz,v(1,1,16),cvr(1,1,16),
     #             cvz(1,1,16),cvzr(1,1,16),idm,idm)
C                                                   theta=pi/4,phi=pi/4
      OPEN(10,FILE='fhol-pd45.dat',STATUS='OLD')
      READ(10,*) nz,nr
      IF(nr.NE.nr0.OR.nz.NE.nz0) THEN
        WRITE(*,*) 'ERROR'
	STOP
      END IF	
      DO i=1,nz
        DO j=1,nr
          READ(10,*) z(i),r(j),v(j,i,17)
          IF(r(j).NE.r0(j).OR.z(i).NE.z0(i)) THEN
            WRITE(*,*) 'ERROR'
	    STOP
          END IF	
        END DO
      END DO
      CLOSE(10)
      CALL SPLIN2D(r,nr,z,nz,v(1,1,17),cvr(1,1,17),
     #             cvz(1,1,17),cvzr(1,1,17),idm,idm)
      
      OPEN(1,FILE='TOUT45.DAT',STATUS='unknown',Form='unformatted')
      WRITE(1) nr,nz,(r(i),i=1,nr),(z(i),i=1,nz)

      DO i=1,17
        
        DO j=1,nz
          WRITE(1) (v(k,j,i),cvr(k,j,i),cvz(k,j,i),cvzr(k,j,i),k=1,nr)
        END DO

      END DO  
       
      END
