C************************ P R E P 5 D (N2W) ******************************

C     This program reads the molecular data for the various
C     configurations and writes down in TOUT5D.DAT the input that is
C     used for the calculation of the PES. Once the file TOUT5D.DAT
C     has been created, there is no point to run PREP5D again. The
C     same set of Z values is used for ALL configurations which
C     speeds up the subsequent interpolations. One uses the property
C     of some spline representations: they are invariant by addition
C     of nodes obtained by interpolation between the given nodes.
C     The only constraint in the program is that the minimum and
C     maximum value of r and Z should be the same for all sets.
C     Derivative is fixed as zero at upper bound.

C     DATA: in file prep5D.dat
C     zmin,zmax: interval of Z values for which the potential is 
C                calculated
C     nz: number of Z values

C       Version 1       A. Salin                   

      IMPLICIT NONE

      INTEGER idm,idc
      PARAMETER(idm=120,idc=18)
C     IDM: maximum number of node points
C     IDC: maximum number of configurations

      INTEGER nconf
      DATA nconf/9/
C     NCONF: actual number of configurations

      LOGICAL fin
      INTEGER nz,i,ndat,j
      DOUBLE PRECISION zmin,zmax,DPL
      DOUBLE PRECISION z(idm),v(idm,idc),cv(idm,idc),vper(idm),
     1     cvper(idm),zv(idm),a(idm),aa(idm),aaa(idm)

C     ----------------------------------------------------------------
      OPEN(1,FILE='prep5D.dat',STATUS='OLD')
        READ(1,*) zmin,zmax
        WRITE(*,*) 'zmin: ',zmin,' zmax: ',zmax     
        READ(1,*) nz
        IF((nz/2)*2.EQ.nz) nz=nz+1
        IF(nz.LE.1) nz=3
      CLOSE(1)

      DO i=1,nz
        z(i)=(i-1)*(zmax-zmin)/DFLOAT(nz-1)+zmin
      END DO
      z(nz)=zmax
C---------------------------------------------------------------------
C     Potential on top, 

C                                                              theta=0
      OPEN(10,FILE='t-p.dat',STATUS='OLD')
      READ(10,*) ndat
      DO i=1,ndat
        READ(10,*) zv(i),vper(i)
      END DO
      CLOSE(10)
      CALL DSPLIN(ndat,zv,vper,cvper,0.D0,0,0.D0,2,a,aa,aaa)

      DO i=1,nz
        v(i,1)=DPL(ndat,zv,vper,cvper,z(i))
      END DO
      CALL DSPLIN(nz,z,v(1,1),cv(1,1),0.D0,0,0.D0,2,a,aa,aaa)
C                                                    theta=pi/2, phi=0
      OPEN(10,FILE='t-pl.dat',STATUS='OLD')
      READ(10,*) ndat
      DO i=1,ndat
        READ(10,*) zv(i),vper(i)
      END DO
      CLOSE(10)
      CALL DSPLIN(ndat,zv,vper,cvper,0.D0,0,0.D0,2,a,aa,aaa)

      DO i=1,nz
        v(i,2)=DPL(ndat,zv,vper,cvper,z(i))
      END DO
      CALL DSPLIN(nz,z,v(1,2),cv(1,2),0.D0,0,0.D0,2,a,aa,aaa)
C                                                  theta=pi/2,phi=pi/4
      OPEN(10,FILE='t-pd.dat',STATUS='OLD')
      READ(10,*) ndat
      DO i=1,ndat
        READ(10,*) zv(i),vper(i)
      END DO
      CLOSE(10)
      CALL DSPLIN(ndat,zv,vper,cvper,0.D0,0,0.D0,2,a,aa,aaa)

      DO i=1,nz
        v(i,3)=DPL(ndat,zv,vper,cvper,z(i))
      END DO
      CALL DSPLIN(nz,z,v(1,3),cv(1,3),0.D0,0,0.D0,2,a,aa,aaa)
C---------------------------------------------------------------------
C     Potential on bridge

C                                                              theta=0
      OPEN(10,FILE='b-p.dat',STATUS='OLD')
      READ(10,*) ndat
      DO i=1,ndat
        READ(10,*) zv(i),vper(i)
      END DO
      CLOSE(10)
      CALL DSPLIN(ndat,zv,vper,cvper,0.D0,0,0.D0,2,a,aa,aaa)

      DO i=1,nz
        v(i,4)=DPL(ndat,zv,vper,cvper,z(i))
      END DO
      CALL DSPLIN(nz,z,v(1,4),cv(1,4),0.D0,0,0.D0,2,a,aa,aaa)
C                                                     theta=pi/2,phi=0
      OPEN(10,FILE='b-pl.dat',STATUS='OLD')
      READ(10,*) ndat
      DO i=1,ndat
        READ(10,*) zv(i),vper(i)
      END DO
      CLOSE(10)
      CALL DSPLIN(ndat,zv,vper,cvper,0.D0,0,0.D0,2,a,aa,aaa)

      DO i=1,nz
        v(i,5)=DPL(ndat,zv,vper,cvper,z(i))
      END DO
      CALL DSPLIN(nz,z,v(1,5),cv(1,5),0.D0,0,0.D0,2,a,aa,aaa)
C                                                  theta=pi/2,phi=pi/2
      OPEN(10,FILE='b-pp.dat',STATUS='OLD')
      READ(10,*) ndat
      DO i=1,ndat
        READ(10,*) zv(i),vper(i)
      END DO
      CLOSE(10)
      CALL DSPLIN(ndat,zv,vper,cvper,0.D0,0,0.D0,2,a,aa,aaa)

      DO i=1,nz
        v(i,6)=DPL(ndat,zv,vper,cvper,z(i))
      END DO
      CALL DSPLIN(nz,z,v(1,6),cv(1,6),0.D0,0,0.D0,2,a,aa,aaa)
C---------------------------------------------------------------------
C     Potential on hole

C                                                              theta=0
      OPEN(10,FILE='h-p.dat',STATUS='OLD')
      READ(10,*) ndat
      DO i=1,ndat
        READ(10,*) zv(i),vper(i)
      END DO
      CLOSE(10)
      CALL DSPLIN(ndat,zv,vper,cvper,0.D0,0,0.D0,2,a,aa,aaa)

      DO i=1,nz
        v(i,7)=DPL(ndat,zv,vper,cvper,z(i))
      END DO
      CALL DSPLIN(nz,z,v(1,7),cv(1,7),0.D0,0,0.D0,2,a,aa,aaa)
C                                                     theta=pi/2,phi=0
      OPEN(10,FILE='h-pl.dat',STATUS='OLD')
      READ(10,*) ndat
      DO i=1,ndat
        READ(10,*) zv(i),vper(i)
      END DO
      CLOSE(10)
      CALL DSPLIN(ndat,zv,vper,cvper,0.D0,0,0.D0,2,a,aa,aaa)

      DO i=1,nz
        v(i,8)=DPL(ndat,zv,vper,cvper,z(i))
      END DO
      CALL DSPLIN(nz,z,v(1,8),cv(1,8),0.D0,0,0.D0,2,a,aa,aaa)
C                                                  theta=pi/2,phi=pi/4
      OPEN(10,FILE='h-pd.dat',STATUS='OLD')
      READ(10,*) ndat
      DO i=1,ndat
        READ(10,*) zv(i),vper(i)
      END DO
      CLOSE(10)
      CALL DSPLIN(ndat,zv,vper,cvper,0.D0,0,0.D0,2,a,aa,aaa)

      DO i=1,nz
        v(i,9)=DPL(ndat,zv,vper,cvper,z(i))
      END DO
      CALL DSPLIN(nz,z,v(1,9),cv(1,9),0.D0,0,0.D0,2,a,aa,aaa)

     
      OPEN(1,FILE='TOUT5D.DAT',STATUS='unknown',form='unformatted')
      WRITE(1) nz,(z(i),i=1,nz)
      DO j=1,nconf
        WRITE(1) (v(i,j),cv(i,j),i=1,nz)
      END DO
      CLOSE(1)
c      OPEN(1,FILE='tout5D.dat',STATUS='unknown')
c      DO i=1,nz
c       WRITE(1,fmt='(1PE20.12,9E20.12)') z(i),(v(i,j),j=1,nconf)
c      END DO
c      CLOSE(1)

      END
