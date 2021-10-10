C     Calculation of the 6D potential for N2/W(100).
C     Routine required for trajectory calculations.

      SUBROUTINE LECTURE

      CALL LECTATOM
      CALL LECTMOL45
      CALL LECTMOL5D

      END
C************************ BLOCK DATA PESII ****************************
      BLOCK DATA PESNW
      IMPLICIT NONE
      
      CHARACTER*50 pesname
      COMMON/PESDEF/pesname
      
      DOUBLE PRECISION rmax,zccmin,ccdef
      COMMON/DATLIM/rmax,zccmin,ccdef
      
      DATA pesname/'N2-W(100)----------------------------------------.'/

C       RMAX: maximum internuclear distance for which the potential can
C             be calculated (Angstroems). .
C       ZCCMIN: minimum distance to surface of the origin O on the 
C             internuclear axis [such that OA=|r_B-r_A|*ccdef and
C             OB=|r_B-r_A|*(1-ccdef)] for which the potential can be
C             calculated (Angstroems).
C       CCDEF: defines the coordinate origin O used to produce the set
C              of molecular data, as explained above.

      DATA rmax,zccmin,ccdef/2.25D0,0.25D0,0.5D0/
      
      END
C************************ P O T 6 D ***********************************
      SUBROUTINE POT6D(rab,potv,dvdxa,dvdya,dvdza,dvdxb,dvdyb,dvdzb,
     1                 iunit,switch)

      IMPLICIT NONE

      DOUBLE PRECISION tranlg,trane
      PARAMETER(tranlg=5.29177249D-1,trane=27.2D0)      

      LOGICAL switch
      INTEGER iunit
      DOUBLE PRECISION rab(6),potv,dvdxa,dvdya,dvdza,dvdxb,dvdyb,dvdzb,
     1       xa,ya,za,xb,yb,zb,xcc,ycc,zcc,r,teta,fi
C     ------------------------------------------------------------------
C     iunit=1: CALL arguments and output in ev and Angstroems.
C     iunit=2: Call arguments and output in atomic units. 
      IF(iunit.LT.1.OR.iunit.GT.2) THEN
         WRITE(*,*) 'Incorrect IUNIT'
         STOP
      END IF
      switch=.FALSE.

      potv=0.D0
      dvdxa=0.D0
      dvdya=0.D0
      dvdza=0.D0
      dvdxb=0.D0
      dvdyb=0.D0
      dvdzb=0.D0
      
      IF(iunit.EQ.2) THEN
C     Transform RAB from atomic units into Angstroems
        za=rab(1)*tranlg
        zb=rab(4)*tranlg 
        xa=rab(2)*tranlg 
        ya=rab(3)*tranlg 
        xb=rab(5)*tranlg 
        yb=rab(6)*tranlg 
      ELSE
        za=rab(1)
        zb=rab(4)
        xa=rab(2)
        ya=rab(3)
        xb=rab(5)
        yb=rab(6)
      END IF

C     From here, everything is in Angstroems and eV.

C     ------------------------------------------------------------------
      CALL CAMBCORD(xa,ya,za,xb,yb,zb,xcc,ycc,zcc,r,teta,fi)
C     ------------------------------------------------------------------
      CALL INTMOL45(xa,ya,za,xb,yb,zb,xcc,ycc,zcc,r,teta,fi,
     1                 potv,dvdxa,dvdya,dvdza,dvdxb,dvdyb,dvdzb,switch)
      IF (switch) RETURN

C     -------------------------

      IF(iunit.EQ.2) THEN
C     Switch from Angstroems and eV to atomic units.
        potv=potv/trane 
        dvdxa=dvdxa*tranlg/trane 
        dvdya=dvdya*tranlg/trane 
        dvdza=dvdza*tranlg/trane 
        dvdxb=dvdxb*tranlg/trane 
        dvdyb=dvdyb*tranlg/trane 
        dvdzb=dvdzb*tranlg/trane
      END IF 

C     ------------------------

      END
