C************************ I N T M O L 5 D ****************************

C     Calculation of the 5D potential by direct interpolation.
C     The vibration is represented by a Z-dependent harmonic term.

      SUBROUTINE INTMOL5D(xa,ya,za,xb,yb,zb,xcc,ycc,zcc,r,teta,fi,
     1                    vinter,switch)

C     XA, YA, ZA: cartesian coordinates of atom A
C     XB, YB, ZB: cartesian coordiantes of atom B
C     XCC, YCC, ZCC: cartesian coordinates of center of charge.
C     R, TETA, PHI: spherical coordinates for internuclear distance.
C     VINTER: potential and derivatives.
C            1: potential, 2: d/dZ, 3: d/dr,
C            4: d/d(phi), 5: d/d(theta), 6: d/dX, 7: d/dY.
C     SWITCH: set to .TRUE. on output in case of problem.

C     Calculations are in Angstroems and eV.

C     ----------------------------------------------------------------
      IMPLICIT NONE

      INTEGER idf6,ids
      PARAMETER(idf6=4,ids=3)
C     IDF6: number of reference sites in (x,y) interpolation.
C     IDS: max. number of configurations for one site.

      INTEGER idas
      PARAMETER (idas=120)
           
      DOUBLE PRECISION pi,pis2
      PARAMETER (pi=3.141592653589793D0,pis2=pi/2.D0)

      LOGICAL switch
      
      INTEGER idelta
      DOUBLE PRECISION delta
      COMMON/SURF/delta,idelta

C     ----------------------------------------------------------------
      DOUBLE PRECISION vint(3,ids),psite(5,idf6),vinter(7),xcc,ycc,zcc,
     1       r,teta,fi,potv,dvdxa,dvdya,dvdza,dvdxb,dvdyb,dvdzb,
     2       xa,ya,za,xb,yb,zb,DPL,DPLP,drxa,drxb,drya,dryb,drza,drzb,
     3       dtetaxa,dtetaxb,dtetaya,dtetayb,dfixa,dfixb,dfiya,dfiyb,
     4       dtetaza,dtetazb,rho2,rho,ome,dome,rme,drme,a(idas),
     5       aa(idas),aaa(idas)
C     ----------------------------------------------------------------

C     Parameters of vibration potential read in LECTMOL5D
      INTEGER nzvib
      DOUBLE PRECISION zvib,omeg,comeg,rm,crm
      COMMON/VIBASYM/zvib(idas),omeg(idas),comeg(idas),rm(idas),
     1               crm(idas),nzvib

C     DATA read by LECTMOL5D
      INTEGER i,j,nzn
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
C     ------------------------------------------------------------------
      IF(zcc.LT.zn(1).OR.zcc.GT.zn(nzn)) THEN
        switch=.true.
	RETURN
      END IF	

C     No dependence on r.
      DO i=1,ids
        vint(3,i)=0.D0
      END DO	
C     ------------------------------------------------------------------ 
C     ------------------------------------------------------------------
C     Potential on top
C     -----------------------------------------
C     Potential on top and theta=0.  
      vint(1,1)=DPL(nzn,zn,v1,cv1z,zcc)
      vint(2,1)=DPLP(nzn,zn,v1,cv1z,zcc)
C     -----------------------------------------
C     Potential on top, theta=pi/2 and fi=0.

      vint(1,2)=DPL(nzn,zn,v2,cv2z,zcc)
      vint(2,2)=DPLP(nzn,zn,v2,cv2z,zcc)
C     -----------------------------------------
C     Potential on top, theta=pi/2 and phi=pi/4

      vint(1,3)=DPL(nzn,zn,v3,cv3z,zcc)
      vint(2,3)=DPLP(nzn,zn,v3,cv3z,zcc)
C     -----------------------------------------
        CALL VHT5D(teta,fi,vint,psite(1,1))
C     ------------------------------------------------------------------
C     Potential on bridge
C     -----------------------------------------
C     Potential on bridge and theta=0.

      vint(1,1)=DPL(nzn,zn,v4,cv4z,zcc)
      vint(2,1)=DPLP(nzn,zn,v4,cv4z,zcc)
C     -----------------------------------------
C     Potential on bridge, theta=pi/2 and phi=0

      vint(1,2)=DPL(nzn,zn,v5,cv5z,zcc)
      vint(2,2)=DPLP(nzn,zn,v5,cv5z,zcc)
C     -----------------------------------------
C     Potential on bridge, theta=pi/2 and phi=pi/2

      vint(1,3)=DPL(nzn,zn,v6,cv6z,zcc)
      vint(2,3)=DPLP(nzn,zn,v6,cv6z,zcc)
C     -----------------------------------------
C     -----------------------------------------
        CALL VBRI5D(teta,fi,vint,psite(1,2))
C     Potential at equivalent point.
        CALL VBRI5D(teta,fi+pis2,vint,psite(1,3))
C     ------------------------------------------------------------------
C     Potential on hole
C     -----------------------------------------
C     Potential on hole and theta=0.
      
      vint(1,1)=DPL(nzn,zn,v7,cv7z,zcc)
      vint(2,1)=DPLP(nzn,zn,v7,cv7z,zcc)
C     -----------------------------------------
C     Potential on hole, theta=pi/2 and phi=0

      vint(1,2)=DPL(nzn,zn,v8,cv8z,zcc)
      vint(2,2)=DPLP(nzn,zn,v8,cv8z,zcc)
C     -----------------------------------------
C     Potential on hole, theta=pi/2 and phi=pi/4

      vint(1,3)=DPL(nzn,zn,v9,cv9z,zcc)
      vint(2,3)=DPLP(nzn,zn,v9,cv9z,zcc)

        CALL VHT5D(teta,fi,vint,psite(1,4))
C     ------------------------------------------------------------------
        CALL SIS(xcc,ycc,psite,vinter)
C     ------------------------------------------------------------------

C      Add vibration potential and derivatives.
      IF(zcc.GE.zvib(nzvib)) THEN
        ome=omeg(nzvib)
	rme=rm(nzvib)
	dome=0.D0
	drme=0.D0
      ELSE
        ome=DPL(nzvib,zvib,omeg,comeg,zcc)
        dome=DPLP(nzvib,zvib,omeg,comeg,zcc)
        rme=DPL(nzvib,zvib,rm,crm,zcc)
        drme=DPLP(nzvib,zvib,rm,crm,zcc)
      END IF
        vinter(1)=vinter(1)+ome*(r-rme)**2
        vinter(2)=vinter(2)+dome*(r-rme)**2-2.D0*ome*(r-rme)*drme
        vinter(3)=2.D0*ome*(r-rme)
      END
C************************ V H T 5 D ******************************
      SUBROUTINE VHT5D(teta,fi,vint,psite)

C     VINT: first index= pot, d/dZ, d/dr.
C           second index= top site geometry (teta=0; teta=pi/2,fi=0;
C            theta=pi/2,fi=pi/4)
C     On output, PSITE contains potential and derivatives for 
C     the given value of TETA and FI.
C     The index corresponds to: 1: function, 2: d/dZ, 3: d/dr,
C     4: d/d(phi), 5: d/d(theta).

      IMPLICIT NONE

      INTEGER ids
      PARAMETER (ids=3)
               
      INTEGER i
      DOUBLE PRECISION vint(3,ids),psite(5),cos4fi,sin4fi,cos2te,sin2te,
     1       v0,v90,dv90dfi,teta,fi
      
      cos4fi=DCOS(4.D0*fi)
      sin4fi=DSIN(4.D0*fi)
      cos2te=DCOS(2.D0*teta)
      sin2te=DSIN(2.D0*teta)

      DO i=1,3
        v0=vint(i,1)
        v90=0.5D0*( (vint(i,2)+vint(i,3))
     1             +(vint(i,2)-vint(i,3))*cos4fi )

        IF(i.EQ.1) THEN 
C                                                         Potential
          psite(1)=0.5D0*( (v0+v90)+(v0-v90)*cos2te )
C                                First derivatives over fi and teta
          dv90dfi=-2.D0*(vint(1,2)-vint(1,3))*sin4fi
          psite(4)=0.5D0*dv90dfi*(1.D0-cos2te)
          psite(5)=-(v0-v90)*sin2te

        ELSE 
C                                       First derivatives over r,Z        
          psite(i)=0.5D0*( (v0+v90)+(v0-v90)*cos2te )

        END IF

      END DO

      END
C********************** V B R I 5 D ********************************
      SUBROUTINE VBRI5D(teta,fi,vint,psite)

C     VINT: first index= pot, d/dZ, d/dr.
C        second index= bridge site geometry (theta=0; theta=pi/2,fi=0;
C                      theta=pi/2,fi=pi/2)
C     On output, PSITE contains potential and derivatives for 
C     the given value of TETA and FI.
C     The index corresponds to: 1: function, 2: d/dZ, 3: d/dr,
C     4: d/d(phi), 5: d/d(theta).

      IMPLICIT NONE

      INTEGER ids
      PARAMETER (ids=3)
               
      INTEGER i      
      DOUBLE PRECISION vint(3,ids),psite(5),cos2fi,sin2fi,cos2te,sin2te,
     1       v0,v90,dv90dfi,teta,fi
      
      cos2fi=DCOS(2.D0*fi)
      sin2fi=DSIN(2.D0*fi)
      cos2te=DCOS(2.D0*teta)
      sin2te=DSIN(2.D0*teta)
 
      DO i=1,3
        v0=vint(i,1)
        v90=0.5D0*( (vint(i,2)+vint(i,3))
     1      +(vint(i,2)-vint(i,3))*cos2fi )

        IF (i.EQ.1) THEN
C                                                            Potential
          psite(i)=0.5D0*( (v0+v90)+(v0-v90)*cos2te) 
C                                  First derivatives over fi and theta
          dv90dfi=-(vint(1,2)-vint(1,3))*sin2fi
          psite(4)=0.5D0*dv90dfi*(1.D0-cos2te)
          psite(5)=-(v0-v90)*sin2te

        ELSE
C                                           First derivatives over r,Z
          psite(i)=0.5D0*( (v0+v90)+(v0-v90)*cos2te )

        END IF

      END DO

      END
