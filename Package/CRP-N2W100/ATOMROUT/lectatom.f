C************************ L E C T A T O M  ***************************
      SUBROUTINE LECTATOM

C     This programs reads the data files for the 3D PES and determines
C     the 3D Interpolation function. The 1D potential is read from the
C     file INTREP.DAT and calculated by interpolation in MUR3D2.

C     Warning: the Z values must be the same for all sites.

C     All distances are in angstroems and energies in eV.

      IMPLICIT NONE

      INTEGER idz,idrep
      PARAMETER(idz=9999,idrep=9999)
C     IDZ: maximum number of Z values
C     IDREP: maximum number of values for 1D potential.

C     DELTA: nearest neighbour distance on the surface (in Angstroems)
C     IDELTA=233 means delta is in Angstroms
C     IDELTA=247 means delta is in a.u.
      INTEGER idelta
      DOUBLE PRECISION delta
      COMMON/SURF/delta,idelta

      DOUBLE PRECISION alpha(idz),beta(idz),b(idz)
C     ----------------------------------------------------------------
      INTEGER nrep
      DOUBLE PRECISION zrep,vrep,cvrep,vasint
      COMMON/POTREP/zrep(idrep),vrep(idrep),cvrep(idrep),vasint,nrep 
C     ----------------------------------------------------------------
C     Top site
      INTEGER n1p
      DOUBLE PRECISION z1p,vz1p,cvz1p
      COMMON/SITIO1P/z1p(idz),vz1p(idz),cvz1p(idz),n1p
      INTEGER n1m
      DOUBLE PRECISION z1m,vz1m,cvz1m
      COMMON/SITIO1M/z1m(idz),vz1m(idz),cvz1m(idz),n1m
C     -----------------------------------------------
C     Bridge site
      INTEGER n2
      DOUBLE PRECISION z2,vz2,cvz2
      COMMON/SITIO2/z2(idz),vz2(idz),cvz2(idz),n2
C     -----------------------------------------------
C     Hole site
      INTEGER n3
      DOUBLE PRECISION z3,vz3,cvz3
      COMMON/SITIO3/z3(idz),vz3(idz),cvz3(idz),n3
C     -----------------------------------------------
C     Bridge-threefold site
      INTEGER n4
      DOUBLE PRECISION z4,vz4,cvz4
      COMMON/SITIO4/z4(idz),vz4(idz),cvz4(idz),n4
C     -----------------------------------------------
C     Top-hole
      INTEGER n5
      DOUBLE PRECISION z5,vz5,cvz5
      COMMON/SITIO5/z5(idz),vz5(idz),cvz5(idz),n5
C     -----------------------------------------------
C     Top-bridge
      INTEGER n6
      DOUBLE PRECISION z6,vz6,cvz6
      COMMON/SITIO6/z6(idz),vz6(idz),cvz6(idz),n6

      INTEGER i
      DOUBLE PRECISION pot,dz1,dz2,vm1,b0,b1,deriz,der1,der2

C     Read 1D potential:

      OPEN(10,FILE='intrep.dat',STATUS='OLD') 
      READ(10,*) delta,idelta
      READ(10,*) vasint
      READ(10,*) nrep 

      READ(10,*) zrep(1),pot,dz1 
      vrep(1)=pot-vasint 
      DO i=2,nrep-1 
        READ(10,*) zrep(i),pot 
        vrep(i)=pot-vasint 
      END DO 
      READ(10,*) zrep(nrep),pot,dz2 
      vrep(nrep)=pot-vasint 
      CLOSE(10) 

      CALL DSPLIN(nrep,zrep,vrep,cvrep,dz1,2,dz2,2,alpha,beta,b)
C     ---------------------------------------------------------- 
C     Read potential for the various sites:
 
C     Top z>0:
 
      OPEN(10,FILE='topp.dat',STATUS='OLD') 
      READ(10,*) n1p

      READ(10,*) z1p(1),pot,dz1
      CALL MUR3D2(0.D0,0.D0,z1p(1),vm1,b0,b1,deriz)
      vz1p(1)=pot-vasint-vm1
      der1=dz1-deriz
      DO  i=2,n1p-1
        READ(10,*) z1p(i),pot
        CALL MUR3D2(0.D0,0.D0,z1p(i),vm1,b0,b1,deriz)
        vz1p(i)=pot-vasint-vm1
      END DO
      READ(10,*) z1p(n1p),pot,dz2
      CALL MUR3D2(0.D0,0.D0,z1p(n1p),vm1,b0,b1,deriz)
      vz1p(n1p)=pot-vasint-vm1
      der2=dz2-deriz
      CLOSE(10)

      CALL DSPLIN(n1p,z1p,vz1p,cvz1p,der1,2,der2,2,alpha,beta,b)

C      Top z<0

      OPEN(10,FILE='topm.dat',STATUS='OLD') 
      READ(10,*) n1m

      READ(10,*) z1m(1),pot,dz1
      CALL MUR3D2(0.D0,0.D0,z1m(1),vm1,b0,b1,deriz)
      vz1m(1)=pot-vasint-vm1
      der1=dz1-deriz
      DO  i=2,n1m-1
        READ(10,*) z1m(i),pot
        CALL MUR3D2(0.D0,0.D0,z1m(i),vm1,b0,b1,deriz)
        vz1m(i)=pot-vasint-vm1
      END DO
      READ(10,*) z1m(n1m),pot,dz2
      CALL MUR3D2(0.D0,0.D0,z1m(n1m),vm1,b0,b1,deriz)
      vz1m(n1m)=pot-vasint-vm1
C     Change of sign for derivative at r=0      
      der2=dz2+deriz
      CLOSE(10)

      CALL DSPLIN(n1m,z1m,vz1m,cvz1m,der1,2,der2,2,alpha,beta,b)
C     ----------------------------------------------------------

      OPEN(10,FILE='bridge.dat',STATUS='OLD')
      READ(10,*) n2

      READ(10,*) z2(1),pot,dz1
      CALL MUR3D2(delta/2.D0,0.D0,z2(1),vm1,b0,b1,deriz)
      vz2(1)=pot-vasint-vm1
      der1=dz1-deriz
      DO i=2,n2-1
        READ(10,*) z2(i),pot
        CALL MUR3D2(delta/2.D0,0.D0,z2(i),vm1,b0,b1,deriz)
        vz2(i)=pot-vasint-vm1
      END DO
      READ(10,*) z2(n2),pot,dz2
      CALL MUR3D2(delta/2.D0,0.D0,z2(n2),vm1,b0,b1,deriz)
      vz2(n2)=pot-vasint-vm1
      der2=dz2-deriz
      CLOSE(10)

      CALL DSPLIN(n2,z2,vz2,cvz2,der1,2,der2,2,alpha,beta,b)

C     ----------------------------------------------------------

      OPEN(10,FILE='hole.dat',STATUS='OLD')
      READ(10,*) n3

      READ(10,*) z3(1),pot,dz1
      IF(z3(1).NE.z2(1)) THEN
        WRITE(*,*) 'Error in data for site 3'
        STOP
      END IF
      CALL MUR3D2(delta/2.D0,delta/2.D0,z3(1),vm1,b0,b1,deriz)
      vz3(1)=pot-vasint-vm1
      der1=dz1-deriz
      DO i=2,n3-1
        READ(10,*) z3(i),pot
        IF(z3(i).NE.z2(i)) THEN
          WRITE(*,*) 'Error in data for site 3'
          STOP
        END IF
        CALL MUR3D2(delta/2.D0,delta/2.D0,z3(i),vm1,b0,b1,deriz)
        vz3(i)=pot-vasint-vm1
      END DO
      READ(10,*) z3(n3),pot,dz2
      IF(z3(n3).NE.z2(n3)) THEN
        WRITE(*,*) 'Error in data for site 3'
        STOP
      END IF
      CALL MUR3D2(delta/2.D0,delta/2.D0,z3(n3),vm1,b0,b1,deriz)
      vz3(n3)=pot-vasint-vm1
      der2=dz2-deriz
      CLOSE(10)

      CALL DSPLIN(n3,z3,vz3,cvz3,der1,2,der2,2,alpha,beta,b)

C     ----------------------------------------------------------

      OPEN(10,FILE='brihol.dat',STATUS='OLD')
      READ(10,*) n4

      READ(10,*) z4(1),pot,dz1
      IF(z4(1).NE.z2(1)) THEN
        WRITE(*,*) 'Error in data for site 4'
        STOP
      END IF
      CALL MUR3D2(delta/2.D0,delta/4.D0,z4(1),vm1,b0,b1,deriz)
      vz4(1)=pot-vasint-vm1
      der1=dz1-deriz
      DO i=2,n4-1
        READ(10,*) z4(i),pot
        IF(z4(i).NE.z2(i)) THEN
          WRITE(*,*) 'Error in data for site 4'
          STOP
        END IF
        CALL MUR3D2(delta/2.D0,delta/4.D0,z4(i),vm1,b0,b1,deriz)
        vz4(i)=pot-vasint-vm1
      END DO
      READ(10,*) z4(n4),pot,dz2
      IF(z4(n4).NE.z2(n4)) THEN
        WRITE(*,*) 'Error in data for site 4'
        STOP
      END IF
      CALL MUR3D2(delta/2.D0,delta/4.D0,z4(n4),vm1,b0,b1,deriz)
      vz4(n4)=pot-vasint-vm1
      der2=dz2-deriz
      CLOSE(10)

      CALL DSPLIN(n4,z4,vz4,cvz4,der1,2,der2,2,alpha,beta,b)

C     ----------------------------------------------------------

      OPEN(10,FILE='tophol.dat',STATUS='OLD')
      READ(10,*) n5

      READ(10,*) z5(1),pot,dz1
      IF(z5(1).NE.z2(1)) THEN
        WRITE(*,*) 'Error in data for site 5'
        STOP
      END IF
      CALL MUR3D2(delta/4.D0,delta/4.D0,z5(1),vm1,b0,b1,deriz)
      vz5(1)=pot-vasint-vm1
      der1=dz1-deriz
      DO i=2,n5-1
        READ(10,*) z5(i),pot
        IF(z5(i).NE.z2(i)) THEN
          WRITE(*,*) 'Error in data for site 5'
          STOP
        END IF
        CALL MUR3D2(delta/4.D0,delta/4.D0,z5(i),vm1,b0,b1,deriz)
        vz5(i)=pot-vasint-vm1
      END DO
      READ(10,*) z5(n5),pot,dz2
      IF(z5(n5).NE.z2(n5)) THEN
        WRITE(*,*) 'Error in data for site 5'
        STOP
      END IF
      CALL MUR3D2(delta/4.D0,delta/4.D0,z5(n5),vm1,b0,b1,deriz)
      vz5(n5)=pot-vasint-vm1
      der2=dz2-deriz
      CLOSE(10)

      CALL DSPLIN(n5,z5,vz5,cvz5,der1,2,der2,2,alpha,beta,b)

C     ----------------------------------------------------------

      OPEN(10,FILE='topbri.dat',STATUS='OLD')
      READ(10,*) n6

      READ(10,*) z6(1),pot,dz1
      IF(z6(1).NE.z2(1)) THEN
        WRITE(*,*) 'Error in data for site 6'
        STOP
      END IF
      CALL MUR3D2(delta/4.D0,0.D0,z6(1),vm1,b0,b1,deriz)
      vz6(1)=pot-vasint-vm1
      der1=dz1-deriz
      DO i=2,n6-1
        READ(10,*) z6(i),pot
        IF(z6(i).NE.z2(i)) THEN
          WRITE(*,*) 'Error in data for site 6'
          STOP
        END IF
        CALL MUR3D2(delta/4.D0,0.D0,z6(i),vm1,b0,b1,deriz)
        vz6(i)=pot-vasint-vm1
      END DO
      READ(10,*) z6(n6),pot,dz2
      IF(z6(n6).NE.z2(n6)) THEN
        WRITE(*,*) 'Error in data for site 6'
        STOP
      END IF
      CALL MUR3D2(delta/4.D0,0.D0,z6(n6),vm1,b0,b1,deriz)
      vz6(n6)=pot-vasint-vm1
      der2=dz2-deriz
      CLOSE(10)

      CALL DSPLIN(n6,z6,vz6,cvz6,der1,2,der2,2,alpha,beta,b)

      CALL FOUR3D

      END
C************************ F O U R 3 D ******************************

C     Matrix allowing to calculate the coefficients of the Fourier 
C     interpolation for the 3D Interpolation Function.
C     This matrix depends on the set of reference sites.

      SUBROUTINE FOUR3D
      IMPLICIT NONE

      INTEGER idf3
      PARAMETER(idf3=6)
C     IDF3: number of reference sites used for the interpolation.

      DOUBLE PRECISION fouri3D
      COMMON/MAFOU3D/fouri3D(idf3,idf3)

      INTEGER i,j
      DOUBLE PRECISION tour(idf3,idf3)

      DATA tour/0.0625D0,0.125D0,0.0625D0,0.25D0,0.25D0,0.25D0,
     1          0.125D0,0.D0,-0.125D0,-0.25D0,0.D0,0.25D0,
     2          0.125D0,-0.25D0,0.125D0,0.D0,0.D0,0.D0,
     3          0.0625D0,0.125D0,0.0625D0,0.D0,-0.25D0,0.D0,
     4          0.0625D0,0.D0,-0.0625D0,0.125D0,0.D0,-0.125D0,
     5          0.03125D0,0.0625D0,0.03125D0,-0.125D0,0.125D0,-0.125D0/

      DO i=1,idf3
        DO j=1,idf3
        fouri3D(i,j)=tour(i,j)
        END DO
      END DO

      END

