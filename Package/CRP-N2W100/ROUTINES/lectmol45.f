C************************ L E C T M O L 4 5 *************************

C     N2/W(100) -17 configurations.

C     READS the Interpolation Function from the file TOUT45.DAT 
C     prepared by PREPMOL45.

C     WARNING: this program must be called AFTER LECTATOM
C     ===================================================
c     since LECTATOM provides the crystal parameter DELTA

C      Version 1.0   27/07/03      A. Salin

      SUBROUTINE LECTMOL45
      IMPLICIT NONE


      INTEGER i,j
      
C     ------------------------------------------------------------------
      INTEGER idm,nrn,nzn
      PARAMETER(idm=16)
C     IDM: maximum number of values of Z and r.
      DOUBLE PRECISION rn,zn,v1,cv1r,cv1z,cv1zr,v2,cv2r,cv2z,cv2zr,v3,
     1 cv3r,cv3z,cv3zr,v4,cv4r,cv4z,cv4zr,v5,cv5r,cv5z,cv5zr,v6,cv6r,
     2 cv6z,cv6zr,v7,cv7r,cv7z,cv7zr,v8,cv8r,cv8z,cv8zr,v9,cv9r,cv9z,
     3 cv9zr,v10,cv10r,cv10z,cv10zr,v11,cv11r,cv11z,cv11zr,v12,cv12r,
     4 cv12z,cv12zr,v13,cv13r,cv13z,cv13zr,v14,cv14r,cv14z,cv14zr,v15,
     5 cv15r,cv15z,cv15zr,v16,cv16r,cv16z,cv16zr,v17,cv17r,cv17z,cv17zr
     
      COMMON/DAT0/rn(idm),zn(idm),nrn,nzn
      COMMON/DAT1/v1(idm,idm),cv1r(idm,idm),
     #            cv1z(idm,idm),cv1zr(idm,idm)
      COMMON/DAT2/v2(idm,idm),cv2r(idm,idm), 
     #            cv2z(idm,idm),cv2zr(idm,idm)
      COMMON/DAT3/v3(idm,idm),cv3r(idm,idm),
     #            cv3z(idm,idm),cv3zr(idm,idm)
      COMMON/DAT4/v4(idm,idm),cv4r(idm,idm), 
     #            cv4z(idm,idm),cv4zr(idm,idm)
      COMMON/DAT5/v5(idm,idm),cv5r(idm,idm),
     #            cv5z(idm,idm),cv5zr(idm,idm)
      COMMON/DAT6/v6(idm,idm),cv6r(idm,idm), 
     #            cv6z(idm,idm),cv6zr(idm,idm)
      COMMON/DAT7/v7(idm,idm),cv7r(idm,idm),
     #            cv7z(idm,idm),cv7zr(idm,idm)
      COMMON/DAT8/v8(idm,idm),cv8r(idm,idm),
     #            cv8z(idm,idm),cv8zr(idm,idm)
      COMMON/DAT9/v9(idm,idm),cv9r(idm,idm),
     #            cv9z(idm,idm),cv9zr(idm,idm)
      COMMON/DAT10/v10(idm,idm),cv10r(idm,idm), 
     #            cv10z(idm,idm),cv10zr(idm,idm)
      COMMON/DAT11/v11(idm,idm),cv11r(idm,idm),
     #            cv11z(idm,idm),cv11zr(idm,idm)
      COMMON/DAT12/v12(idm,idm),cv12r(idm,idm), 
     #            cv12z(idm,idm),cv12zr(idm,idm)
      COMMON/DAT13/v13(idm,idm),cv13r(idm,idm),
     #            cv13z(idm,idm),cv13zr(idm,idm)
      COMMON/DAT14/v14(idm,idm),cv14r(idm,idm), 
     #            cv14z(idm,idm),cv14zr(idm,idm)
      COMMON/DAT15/v15(idm,idm),cv15r(idm,idm),
     #            cv15z(idm,idm),cv15zr(idm,idm)
      COMMON/DAT16/v16(idm,idm),cv16r(idm,idm),
     #            cv16z(idm,idm),cv16zr(idm,idm)
      COMMON/DAT17/v17(idm,idm),cv17r(idm,idm),
     #            cv17z(idm,idm),cv17zr(idm,idm)

C     ------------------------------------------------------------------
      OPEN(1,FILE='TOUT45.DAT',STATUS='old',FORM='unformatted')
      READ(1) nrn,nzn,(rn(i),i=1,nrn),(zn(i),i=1,nzn)
C
C     --------------------------------------------
C     Reading: Potential on top, theta=0
      DO j=1,nzn
        READ(1) (v1(i,j),cv1r(i,j),cv1z(i,j),cv1zr(i,j),i=1,nrn)
      END DO
C     --------------------------------------------
C     Reading: Potential on top, theta=pi/2 and phi=0.
      DO j=1,nzn
        READ(1) (v2(i,j),cv2r(i,j),cv2z(i,j),cv2zr(i,j),i=1,nrn)
      END DO
C     --------------------------------------------
C     Reading: Potential on top, theta=pi/2 and phi=pi/4.
      DO j=1,nzn
        READ(1) (v3(i,j),cv3r(i,j),cv3z(i,j),cv3zr(i,j),i=1,nrn)
      END DO
C     --------------------------------------------
C     Reading: Potential on bridge, theta=0.
      DO j=1,nzn
        READ(1) (v4(i,j),cv4r(i,j),cv4z(i,j),cv4zr(i,j),i=1,nrn)
      END DO
C     --------------------------------------------
C     Reading: Potential on bridge, theta=pi/2, phi=0.
      DO j=1,nzn
        READ(1) (v5(i,j),cv5r(i,j),cv5z(i,j),cv5zr(i,j),i=1,nrn)
      END DO
C     --------------------------------------------
C     Reading: Potential on bridge, theta=pi/2, phi=pi/2.
      DO j=1,nzn
        READ(1) (v6(i,j),cv6r(i,j),cv6z(i,j),cv6zr(i,j),i=1,nrn)
      END DO
C     --------------------------------------------
C     Reading: Potential on hole, theta=0.
      DO j=1,nzn
        READ(1) (v7(i,j),cv7r(i,j),cv7z(i,j),cv7zr(i,j),i=1,nrn)
      END DO
C     --------------------------------------------
C     Reading: Potential on hole, theta=pi/2, phi=0.
      DO j=1,nzn
        READ(1) (v8(i,j),cv8r(i,j),cv8z(i,j),cv8zr(i,j),i=1,nrn)
      END DO
C     --------------------------------------------
C     Reading: Potential on hole, theta=pi/2 and phi=pi/4
      DO j=1,nzn
        READ(1) (v9(i,j),cv9r(i,j),cv9z(i,j),cv9zr(i,j),i=1,nrn)
      END DO
C     --------------------------------------------
C     Reading: Potential on top, theta=pi/4 and phi=0
      DO j=1,nzn
        READ(1) (v10(i,j),cv10r(i,j),cv10z(i,j),cv10zr(i,j),i=1,nrn)
      END DO
C     --------------------------------------------
C     Reading: Potential on top, theta=pi/4 and phi=pi/4
      DO j=1,nzn
        READ(1) (v11(i,j),cv11r(i,j),cv11z(i,j),cv11zr(i,j),i=1,nrn)
      END DO
C     --------------------------------------------
C     Reading: Potential on bridge, theta=pi/4 and phi=0.
      DO j=1,nzn
        READ(1) (v12(i,j),cv12r(i,j),cv12z(i,j),cv12zr(i,j),i=1,nrn)
      END DO
C     --------------------------------------------
C     Reading: Potential on bridge, theta=pi/4 and phi=pi/2
      DO j=1,nzn
        READ(1) (v13(i,j),cv13r(i,j),cv13z(i,j),cv13zr(i,j),i=1,nrn)
      END DO
C     --------------------------------------------
C     Reading: Potential on bridge, theta=pi/2 and phi=pi/4.
      DO j=1,nzn
        READ(1) (v14(i,j),cv14r(i,j),cv14z(i,j),cv14zr(i,j),i=1,nrn)
      END DO
C     --------------------------------------------
C     Reading: Potential on bridge, theta=pi/4 and phi=pi/4.
      DO j=1,nzn
        READ(1) (v15(i,j),cv15r(i,j),cv15z(i,j),cv15zr(i,j),i=1,nrn)
      END DO
C     --------------------------------------------
C     Reading: Potential on hole, theta=pi/4 and phi=0.
      DO j=1,nzn
        READ(1) (v16(i,j),cv16r(i,j),cv16z(i,j),cv16zr(i,j),i=1,nrn)
      END DO
C     --------------------------------------------
C     Reading: Potential on hole, theta=pi/4 and phi=pi/4
      DO j=1,nzn
        READ(1) (v17(i,j),cv17r(i,j),cv17z(i,j),cv17zr(i,j),i=1,nrn)
      END DO
      CLOSE(1)

      CALL FOUR

      END
C************************ F O U R ********************************

C     Matrix allowing to calculate the coefficients of the Fourier 
C     interpolation for the Interpolation Function.
C     This matrix depends on the set of reference sites.
C     The first index labels the site, the second the Fourier component
C     Site index: 1: top, 2: bridge1 (x=0.5, Y=0), 3: bridge2 (x=0,
C                 y=0.5), 4: hole.

      SUBROUTINE FOUR
      IMPLICIT NONE

      INTEGER idf6
      PARAMETER(idf6=4)
C     IDF6: number of reference sites used for the X,Y interpolation.

      DOUBLE PRECISION fourixy
      COMMON/MAFOUXY/fourixy(idf6,idf6)

      INTEGER i,j

      DOUBLE PRECISION tour(idf6,idf6)
      DATA tour/0.25D0,0.25D0,0.25D0,0.25D0,0.25D0,-0.25D0,0.25D0,
     1          -0.25D0,0.25D0,0.25D0,-0.25D0,-0.25D0,0.125D0,
     2          -0.125D0,-0.125D0,0.125D0/

      DO i=1,idf6
        DO j=1,idf6
        fourixy(i,j)=tour(i,j)
        END DO
      END DO

      END
