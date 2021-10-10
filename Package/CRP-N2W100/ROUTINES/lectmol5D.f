C************************ L E C T M O L 5 D ***************************

C     READS the potential from TOUT5D.DAT prepared
C     by PREP5D.

      SUBROUTINE LECTMOL5D
      IMPLICIT NONE

      INTEGER idas
      PARAMETER(idas=120)
C     idas: maximum number of values of Z.

      INTEGER isatat
      DOUBLE PRECISION disatat
      COMMON/SURF/disatat,isatat
      
C     -----------------------------------------------------------------
      INTEGER i,j
      DOUBLE PRECISION a(idas),aa(idas),aaa(idas)
      
      DOUBLE PRECISION zswinf,zswsup
      COMMON/ASYSW/zswinf,zswsup
     
      INTEGER nzvib
      DOUBLE PRECISION zvib,omeg,comeg,rm,crm
      COMMON/VIBASYM/zvib(idas),omeg(idas),comeg(idas),rm(idas),
     1               crm(idas),nzvib
     
      INTEGER nzn
      DOUBLE PRECISION zn,v1,cv1z,v2,cv2z,v3,cv3z,v4,cv4z,v5,cv5z,v6,
     1     cv6z,v7,cv7z,v8,cv8z,v9,cv9z
      COMMON/DAT5D0/zn(idas),nzn
      COMMON/DAT5D1/v1(idas),cv1z(idas)
      COMMON/DAT5D2/v2(idas),cv2z(idas) 
      COMMON/DAT5D3/v3(idas),cv3z(idas)
      COMMON/DAT5D4/v4(idas),cv4z(idas) 
      COMMON/DAT5D5/v5(idas),cv5z(idas)
      COMMON/DAT5D6/v6(idas),cv6z(idas) 
      COMMON/DAT5D7/v7(idas),cv7z(idas)
      COMMON/DAT5D8/v8(idas),cv8z(idas)
      COMMON/DAT5D9/v9(idas),cv9z(idas)

C     Read Z values for switch from 6D to 5D.
      OPEN(1,FILE='zsw.dat',STATUS='old')
        READ(1,*) zswinf,zswsup
      CLOSE(1)
      	
C     Read force constant and equilibrium distance for asymptotic 
C     vibration potential
      OPEN(1,FILE='vibasym.dat',STATUS='old')
        READ(1,*) nzvib
	DO i=1,nzvib
	  READ(1,*) zvib(i),omeg(i),rm(i)
	END DO
	CALL DSPLIN(nzvib,zvib,omeg,comeg,0.D0,0,0.D0,2,a,aa,aaa)
	CALL DSPLIN(nzvib,zvib,rm,crm,0.D0,0,0.D0,2,a,aa,aaa)
      CLOSE(1)	               

C     ------------------------------------------------------------------
      OPEN(1,FILE='TOUT5D.DAT',STATUS='old',Form='unformatted')
      READ(1) nzn,(zn(i),i=1,nzn)
C
C     Here the length units of the input data is Angstroems!!!!
C     --------------------------------------------
C     Reading: Potential on top, theta=0
        READ(1) (v1(i),cv1z(i),i=1,nzn)
C     --------------------------------------------
C     Reading: Potential on top, theta=pi/2 and phi=0.
        READ(1) (v2(i),cv2z(i),i=1,nzn)
C     --------------------------------------------
C     Reading: Potential on top, theta=pi/2 and phi=pi/4.
        READ(1) (v3(i),cv3z(i),i=1,nzn)
C     --------------------------------------------
C     Reading: Potential on bridge, theta=0.
        READ(1) (v4(i),cv4z(i),i=1,nzn)
C     --------------------------------------------
C     Reading: Potential on bridge, theta=pi/2, phi=0.
        READ(1) (v5(i),cv5z(i),i=1,nzn)
C     --------------------------------------------
C     Reading: Potential on bridge, theta=pi/2, phi=pi/2.
        READ(1) (v6(i),cv6z(i),i=1,nzn)
C     --------------------------------------------
C     Reading: Potential on hole, theta=0.
        READ(1) (v7(i),cv7z(i),i=1,nzn)
C     --------------------------------------------
C     Reading: Potential on hole, theta=pi/2, phi=0.
        READ(1) (v8(i),cv8z(i),i=1,nzn)
C     --------------------------------------------
C     Reading: Potential on hole, theta=pi/2 and phi=pi/4
        READ(1) (v9(i),cv9z(i),i=1,nzn)
      CLOSE(1)

      END
