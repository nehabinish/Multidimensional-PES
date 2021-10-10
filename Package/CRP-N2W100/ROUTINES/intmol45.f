C************************ I N T M O L 4 5 *****************************

C     Calculation of the 6D potential for N2/W(100) - 17 configurations

      SUBROUTINE INTMOL45(xa,ya,za,xb,yb,zb,xcc,ycc,zcc,r,teta,fi,
     1                potv,dvdxa,dvdya,dvdza,dvdxb,dvdyb,dvdzb,switch)
C     -----------------------------------------------------------------
c     XA, YA, ZA: cartesian coordinates of atom A
C     XB, YB, ZB: cartesian coordiantes of atom B
C     XCC, YCC, ZCC: cartesian coordinates of center of charge.
C     R, THETA, PHI: spherical coordinates for internuclear distance.
C     POTV: calculated potential
C     DVDXA... : derivatives with respect to XA...
C     SWITCH: set to .TRUE. on output in case of problem.

C     Calculations are in Angstroems and eV.

      IMPLICIT NONE

      INTEGER idf6,ids
      PARAMETER(idf6=4,ids=7)
C     IDF6: number of reference sites in (x,y) interpolation.
C     IDS: max. number of configurations for one site.
           
      DOUBLE PRECISION pi,pis2
      PARAMETER (pi=3.141592653589793D0,pis2=pi/2.D0)

      INTEGER idelta
      DOUBLE PRECISION delta
      COMMON/SURF/delta,idelta

C     -----------------------------------------------------------------
      LOGICAL switch
      INTEGER i,j,irp,izp
      DOUBLE PRECISION vint(3,ids),psite(5,idf6),vinter(7),xcc,ycc,zcc,
     1       r,teta,fi,potv,dvdxa,dvdya,dvdza,dvdxb,dvdyb,dvdzb,
     2       xa,ya,za,xb,yb,zb,DPL,DPLP,drxa,drxb,drya,dryb,drza,drzb,
     3       dtetaxa,dtetaxb,dtetaya,dtetayb,dfixa,dfixb,dfiya,dfiyb,
     4       dtetaza,dtetazb,rho2,rho,dxa,dya,dza,dxb,dyb,dzb,va,vb,
     5       vint5D(7),vfin(7),fsw,dfsw,dzsw
C     -----------------------------------------------------------------

C     Z<ZSWINF: full 6D interpolation
C     ZSWINF<Z<ZSWSUP: switch from 6D to 5D interpolation
C     Z>ZSWSUP: 5D interpolation + harmonic potential for vibration.
      DOUBLE PRECISION zswinf,zswsup
      COMMON/ASYSW/zswinf,zswsup

C     Data read in LECTMOL45
      INTEGER idm,nrn,nzn
      PARAMETER (idm=16)
C     IDM: maximum number of values of Z and r.
      DOUBLE PRECISION rn,zn,v1,cv1r,cv1z,cv1zr,v2,cv2r,cv2z,cv2zr,v3,
     1 cv3r,cv3z,cv3zr,v4,cv4r,cv4z,cv4zr,v5,cv5r,cv5z,cv5zr,v6,cv6r,
     2 cv6z,cv6zr,v7,cv7r,cv7z,cv7zr,v8,cv8r,cv8z,cv8zr,v9,cv9r,cv9z,
     3 cv9zr,v10,cv10r,cv10z,cv10zr,v11,cv11r,cv11z,cv11zr,v12,cv12r,
     4 cv12z,cv12zr,v13,cv13r,cv13z,cv13zr,v14,cv14r,cv14z,cv14zr,v15,
     5 cv15r,cv15z,cv15zr,v16,cv16r,cv16z,cv16zr,v17,cv17r,cv17z,cv17zr
     
      COMMON/DAT0/rn(idm),zn(idm),nrn,nzn
      COMMON/DAT1/v1(idm,idm),cv1r(idm,idm),
     1            cv1z(idm,idm),cv1zr(idm,idm)
      COMMON/DAT2/v2(idm,idm),cv2r(idm,idm), 
     1            cv2z(idm,idm),cv2zr(idm,idm)
      COMMON/DAT3/v3(idm,idm),cv3r(idm,idm),
     1            cv3z(idm,idm),cv3zr(idm,idm)
      COMMON/DAT4/v4(idm,idm),cv4r(idm,idm), 
     1            cv4z(idm,idm),cv4zr(idm,idm)
      COMMON/DAT5/v5(idm,idm),cv5r(idm,idm),
     1            cv5z(idm,idm),cv5zr(idm,idm)
      COMMON/DAT6/v6(idm,idm),cv6r(idm,idm), 
     1            cv6z(idm,idm),cv6zr(idm,idm)
      COMMON/DAT7/v7(idm,idm),cv7r(idm,idm),
     1            cv7z(idm,idm),cv7zr(idm,idm)
      COMMON/DAT8/v8(idm,idm),cv8r(idm,idm),
     1            cv8z(idm,idm),cv8zr(idm,idm)
      COMMON/DAT9/v9(idm,idm),cv9r(idm,idm),
     1            cv9z(idm,idm),cv9zr(idm,idm)
      COMMON/DAT10/v10(idm,idm),cv10r(idm,idm), 
     1            cv10z(idm,idm),cv10zr(idm,idm)
      COMMON/DAT11/v11(idm,idm),cv11r(idm,idm),
     1            cv11z(idm,idm),cv11zr(idm,idm)
      COMMON/DAT12/v12(idm,idm),cv12r(idm,idm), 
     1            cv12z(idm,idm),cv12zr(idm,idm)
      COMMON/DAT13/v13(idm,idm),cv13r(idm,idm),
     1            cv13z(idm,idm),cv13zr(idm,idm)
      COMMON/DAT14/v14(idm,idm),cv14r(idm,idm), 
     1            cv14z(idm,idm),cv14zr(idm,idm)
      COMMON/DAT15/v15(idm,idm),cv15r(idm,idm),
     1            cv15z(idm,idm),cv15zr(idm,idm)
      COMMON/DAT16/v16(idm,idm),cv16r(idm,idm),
     1            cv16z(idm,idm),cv16zr(idm,idm)
      COMMON/DAT17/v17(idm,idm),cv17r(idm,idm),
     1            cv17z(idm,idm),cv17zr(idm,idm)

C     =================================================================

      IF(zcc.LT.zswsup) THEN
      
        IF(zcc.LT.zn(1).OR.zcc.GT.zn(nzn)) THEN
          switch=.true.
	  RETURN
        END IF	
        IF(r.LT.rn(1).OR.r.GT.rn(nrn)) THEN
          switch=.true.
	  RETURN
        END IF	

        CALL UBICAR(rn,r,nrn,irp,switch)
        IF (switch) RETURN
        CALL UBICAR(zn,zcc,nzn,izp,switch)
        IF (switch) RETURN
        CALL SPLCOFII(rn,zn,r,zcc,irp,izp,idm,idm)
C     The latter instructions implies that the first index in VINT
C     correspond to: 1- function, 2-d/dz, 3-d/dr	
C     -----------------------------------------------------------------
C     Potential on top
C     -----------------------------------------
C     Potential on top and theta=0.  
        CALL SPLINII(v1,cv1r,cv1z,cv1zr,irp,izp,vint(1,1),idm,idm)
C     -----------------------------------------
C     Potential on top, theta=pi/2 and phi=0.

        CALL SPLINII(v2,cv2r,cv2z,cv2zr,irp,izp,vint(1,2),idm,idm)
C     -----------------------------------------
C     Potential on top, theta=pi/2 and phi=pi/4

        CALL SPLINII(v3,cv3r,cv3z,cv3zr,irp,izp,vint(1,3),idm,idm)
C     -----------------------------------------
C     Potential on top, theta=pi/4 and phi=0.

        CALL SPLINII(v10,cv10r,cv10z,cv10zr,irp,izp,vint(1,4),idm,idm)
C     -----------------------------------------
C     Potential on top, theta=pi/4 and phi=pi/4

        CALL SPLINII(v11,cv11r,cv11z,cv11zr,irp,izp,vint(1,5),idm,idm)
C     -----------------------------------------
        CALL VHT45(teta,fi,vint,psite(1,1))
C     ------------------------------------------------------------------
C     Potential on bridge
C     -----------------------------------------
C     Potential on bridge and theta=0.

        CALL SPLINII(v4,cv4r,cv4z,cv4zr,irp,izp,vint(1,1),idm,idm)
C     -----------------------------------------
C     Potential on bridge, theta=pi/2 and phi=0

        CALL SPLINII(v5,cv5r,cv5z,cv5zr,irp,izp,vint(1,2),idm,idm)
C     -----------------------------------------
C     Potential on bridge, theta=pi/2 and phi=pi/2

        CALL SPLINII(v6,cv6r,cv6z,cv6zr,irp,izp,vint(1,3),idm,idm)
C     -----------------------------------------
C     Potential on bridge, theta=pi/4 and phi=0

        CALL SPLINII(v12,cv12r,cv12z,cv12zr,irp,izp,vint(1,4),idm,idm)
C     -----------------------------------------
C     Potential on bridge, theta=pi/4 and phi=pi/2

        CALL SPLINII(v13,cv13r,cv13z,cv13zr,irp,izp,vint(1,5),idm,idm)
C     -----------------------------------------
C     Potential on bridge, theta=pi/2 and phi=pi/4

        CALL SPLINII(v14,cv14r,cv14z,cv14zr,irp,izp,vint(1,6),idm,idm)
C     -----------------------------------------
C     Potential on bridge, theta=pi/4 and phi=pi/4

        CALL SPLINII(v15,cv15r,cv15z,cv15zr,irp,izp,vint(1,7),idm,idm)
C     -----------------------------------------
        CALL VBRI45(teta,fi,vint,psite(1,2))
C     Potential at equivalent point.
        CALL VBRI45(teta,fi+pis2,vint,psite(1,3))
C     ------------------------------------------------------------------
C     Potential on hole
C     -----------------------------------------
C     Potential on hole and theta=0.
      
        CALL SPLINII(v7,cv7r,cv7z,cv7zr,irp,izp,vint(1,1),idm,idm)
C     -----------------------------------------
C     Potential on hole, theta=pi/2 and phi=0

        CALL SPLINII(v8,cv8r,cv8z,cv8zr,irp,izp,vint(1,2),idm,idm)
C     -----------------------------------------
C     Potential on hole, theta=pi/2 and phi=pi/4

        CALL SPLINII(v9,cv9r,cv9z,cv9zr,irp,izp,vint(1,3),idm,idm)
C     -----------------------------------------
C     Potential on hole, theta=pi/4 and phi=0

        CALL SPLINII(v16,cv16r,cv16z,cv16zr,irp,izp,vint(1,4),idm,idm)
C     -----------------------------------------
C     Potential on hole, theta=pi/4 and phi=pi/4

        CALL SPLINII(v17,cv17r,cv17z,cv17zr,irp,izp,vint(1,5),idm,idm)

        CALL VHT45(teta,fi,vint,psite(1,4))
C     ------------------------------------------------------------------
        CALL SIS(xcc,ycc,psite,vinter)
C     ------------------------------------------------------------------

        CALL AT3D2(xa,ya,za,va,dxa,dya,dza,switch)
        CALL AT3D2(xb,yb,zb,vb,dxb,dyb,dzb,switch)

      END IF
      
C     =================================================================
      IF(zcc.GT.zswinf) THEN
        CALL INTMOL5D(xa,ya,za,xb,yb,zb,xcc,ycc,zcc,r,teta,fi,
     1                vint5D,switch)
      END IF
C     =================================================================

C     Transformation from (xcc,ycc,zcc,r,teta,fi) to atomic coordinates
        drxa=-(xb-xa)/r
        drxb=-drxa
        drya=-(yb-ya)/r
        dryb=-drya
        drza=-(zb-za)/r
        drzb=-drza
        rho2=(xb-xa)**2+(yb-ya)**2
        rho=DSQRT(rho2)
        IF (rho2.GT.0.D0) THEN
          dtetaxa=-(zb-za)*(xb-xa)/(r*r*rho)
          dtetaya=-(zb-za)*(yb-ya)/(r*r*rho)
          dfixa=(yb-ya)/rho2
          dfiya=-(xb-xa)/rho2
        ELSE
          dtetaxa=0.D0
          dtetaya=0.D0
          dfixa=0.D0
          dfiya=0.D0
        END IF
        dtetaza=rho/(r*r)
        dtetaxb=-dtetaxa
        dtetayb=-dtetaya
        dtetazb=-dtetaza
        dfixb=-dfixa
        dfiyb=-dfiya

      IF(zcc.GE.zswsup) THEN

          potv=vint5D(1)
          dvdxa=0.5D0*vint5D(6)+drxa*vint5D(3)+dtetaxa*vint5D(5)+
     1          dfixa*vint5D(4)
          dvdxb=0.5D0*vint5D(6)+drxb*vint5D(3)+dtetaxb*vint5D(5)+
     1          dfixb*vint5D(4)
          dvdya=0.5D0*vint5D(7)+drya*vint5D(3)+dtetaya*vint5D(5)+
     1          dfiya*vint5D(4)
          dvdyb=0.5D0*vint5D(7)+dryb*vint5D(3)+dtetayb*vint5D(5)+
     1          dfiyb*vint5D(4)
          dvdza=0.5D0*vint5D(2)+drza*vint5D(3)+dtetaza*vint5D(5)
          dvdzb=0.5D0*vint5D(2)+drzb*vint5D(3)+dtetazb*vint5D(5)

      ELSE 

        IF(zcc.GT.zswinf) THEN

          dzsw=(zswsup-zswinf)     
          fsw=DCOS(pi*(zcc-zswinf)/dzsw)
          dfsw=-pi*DSIN(pi*(zcc-zswinf)/dzsw)/dzsw
          vinter(1)=vinter(1)+va+vb
	  dxa=dxa*(1.D0+fsw)*0.5D0
	  dya=dya*(1.D0+fsw)*0.5D0
	  dza=dza*(1.D0+fsw)*0.5D0
	  dxb=dxb*(1.D0+fsw)*0.5D0
	  dyb=dyb*(1.D0+fsw)*0.5D0
	  dzb=dzb*(1.D0+fsw)*0.5D0
	  DO i=1,7
	    vfin(i)=0.5D0*((vinter(i)+vint5D(i))
     1              +(vinter(i)-vint5D(i))*fsw)
          END DO
	  vfin(2)=vfin(2)+(vinter(1)-vint5D(1))*dfsw*0.5D0
        
        ELSE

          vfin(1)=vinter(1)+va+vb
          DO i=2,7
	    vfin(i)=vinter(i)
	  END DO	

        END IF

        potv=vfin(1)
        dvdxa=dxa+0.5D0*vfin(6)+drxa*vfin(3)+dtetaxa*vfin(5)+
     1        dfixa*vfin(4)
        dvdxb=dxb+0.5D0*vfin(6)+drxb*vfin(3)+dtetaxb*vfin(5)+
     1        dfixb*vfin(4)
        dvdya=dya+0.5D0*vfin(7)+drya*vfin(3)+dtetaya*vfin(5)+
     1        dfiya*vfin(4)
        dvdyb=dyb+0.5D0*vfin(7)+dryb*vfin(3)+dtetayb*vfin(5)+
     1        dfiyb*vfin(4)
        dvdza=dza+0.5D0*vfin(2)+drza*vfin(3)+dtetaza*vfin(5)
        dvdzb=dzb+0.5D0*vfin(2)+drzb*vfin(3)+dtetazb*vfin(5)

      END IF
	
C     ----------------------------------------------------------------

      END
C************************ S I S  **********************************

C     Calculation of the 6D interpolation function as well as its
C     first and second derivatives using the Fourier series expansion.
C     The position the molecule center of charge is xcc,ycc. 
C     PSITE: input. Contains the function and its derivatives with
C            respect to r,Z,phi and theta (first index up to 5, see 
C            below for definition) for various sites labelled with 
C            the second index.
C     FOURIXY: comes from the solution of the linear system associated
C              with the Fourier expansion and reference sites. It is
C              calculated in FOUR which must be called once before
C              the first call to SIS
C     VINTER is the output.
C     The index corresponds to: 1: function, 2: d/dZ, 3: d/dr,
C     4: d/d(phi), 5: d/d(theta), 6: d/dX, 7: d/dY.

      SUBROUTINE SIS(xcc,ycc,psite,vinter)
      IMPLICIT NONE
      
      INTEGER idf6
      PARAMETER(idf6=4)
C     IDF6: number of sites used for the (x,y) interpolation.
      
      DOUBLE PRECISION fourixy
      COMMON/MAFOUXY/fourixy(idf6,idf6)
C     see subroutine FOUR

      INTEGER i,j,ksit      
      DOUBLE PRECISION a(idf6),psite(5,idf6),fxy(idf6),dfdx(idf6),
     1       dfdy(idf6),vinter(7),xcc,ycc

      CALL FOURXY(xcc,ycc,fxy,dfdx,dfdy)
      
C                     Calcul du potentiel et des derivees premieres/X,Y
      vinter(1)=0.D0
      vinter(6)=0.D0
      vinter(7)=0.D0
      DO j=1,idf6
        a(j)=0.D0
        DO ksit=1,idf6
          a(j)=a(j)+fourixy(ksit,j)*psite(1,ksit)
        END DO
        vinter(1)=vinter(1)+a(j)*fxy(j)
        vinter(6)=vinter(6)+a(j)*dfdx(j)
        vinter(7)=vinter(7)+a(j)*dfdy(j)
      END DO

C          First order derivatives with respect to r, Z, phi, theta

      DO i=2,5
        vinter(i)=0.D0
        DO j=1,idf6
          a(j)=0.D0
          DO ksit=1,idf6
            a(j)=a(j)+fourixy(ksit,j)*psite(i,ksit)
          END DO
          vinter(i)=vinter(i)+a(j)*fxy(j)
        END DO
      END DO

      END
C************************ F O U R X Y ******************************

C     Calculates the basis functions of the Fourier expansion and 
C     their derivatives with respect to X and Y.

      SUBROUTINE FOURXY(xcc,ycc,fxy,dfdx,dfdy)
      IMPLICIT NONE

      INTEGER idf6
      PARAMETER (idf6=4)
C     IDF6: number of sites used for the (X,Y) interpolation.

      DOUBLE PRECISION pi
      PARAMETER (pi=3.141592653589793D0)

      INTEGER idelta
      DOUBLE PRECISION delta
      COMMON/SURF/delta,idelta
C     WARNING: check carefully units for delta (interatomic distance)
C     ========

      DOUBLE PRECISION fxy(idf6),dfdx(idf6),dfdy(idf6),fac,fbigx,fbigy,
     1       xcc,ycc


      fac=2.D0*pi/delta
      fbigx=fac*xcc
      fbigy=fac*ycc
C                                    Basis functions
      fxy(1)=1.D0
      fxy(2)=DCOS(fbigx)
      fxy(3)=DCOS(fbigy)
      fxy(4)=DCOS(fbigx+fbigy)+DCOS(fbigx-fbigy)

C                    Derivatives with respect to crystal coordinates
      dfdx(1)=0.D0
      dfdx(2)=-fac*DSIN(fbigx)
      dfdx(3)=0.D0
      dfdx(4)=-fac*(DSIN(fbigx+fbigy)+DSIN(fbigx-fbigy))

      dfdy(1)=0.D0
      dfdy(2)=0.D0
      dfdy(3)=-fac*DSIN(fbigy)
      dfdy(4)=-fac*(DSIN(fbigx+fbigy)-DSIN(fbigx-fbigy))

      END
C************************ V H T 4 5 ******************************
      SUBROUTINE VHT45(teta,fi,vint,psite)

C     Interpolation over hole and top with three values of theta
C     (0, pi/4, pi/2) and two values of phi (0, pi/4)
C     VINT: first index= pot, d/dZ, d/dr.
C           second index= top site geometry (teta=0; teta=pi/2,phi=0;
C            theta=pi/2,phi=pi/4; theta=45,phi=0; theta=45,phi=pi/4)
C     On output, PSITE contains potential and derivatives for 
C     the given value of TETA and FI.
C     The index corresponds to: 1: function, 2: d/dZ, 3: d/dr,
C     4: d/d(phi), 5: d/d(theta).

      IMPLICIT NONE

      INTEGER ids
      PARAMETER (ids=5)
               
      INTEGER i
      DOUBLE PRECISION vint(3,ids),psite(5),cos4fi,sin4fi,cos2te,sin2te,
     1       v0,v45,v90,dv45dfi,dv90dfi,teta,fi,cos4te,sin4te
      
      cos4fi=DCOS(4.D0*fi)
      sin4fi=DSIN(4.D0*fi)
      cos2te=DCOS(2.D0*teta)
      sin2te=DSIN(2.D0*teta)
      cos4te=DCOS(4.D0*teta)
      sin4te=DSIN(4.D0*teta)

      DO i=1,3
        v0=vint(i,1)
        v90=0.5D0*( (vint(i,2)+vint(i,3))
     1             +(vint(i,2)-vint(i,3))*cos4fi )
        v45=0.5D0*( (vint(i,4)+vint(i,5))
     1             +(vint(i,4)-vint(i,5))*cos4fi )

        IF(i.EQ.1) THEN 
C                                                         Potential
          psite(1)=0.25D0*( (v0+v90+2.D0*v45)+2.D0*(v0-v90)*cos2te 
     1                    + (v0+v90-2.D0*v45)*cos4te )
C                                First derivatives over fi and teta
          dv90dfi=-2.D0*(vint(1,2)-vint(1,3))*sin4fi
          dv45dfi=-2.D0*(vint(1,4)-vint(1,5))*sin4fi
          psite(4)=0.25D0*(dv90dfi*(1.D0-2.D0*cos2te+cos4te)
     1             + 2.D0*dv45dfi*(1.D0-cos4te) )
          psite(5)=-(v0-v90)*sin2te -(v0+v90-2.D0*v45)*sin4te

        ELSE 
C                                       First derivatives over r,Z        
         psite(i)=0.25D0*( (v0+v90+2.D0*v45)+2.D0*(v0-v90)*cos2te 
     1                   + (v0+v90-2.D0*v45)*cos4te )

        END IF

      END DO

      END
C********************** V B R I 4 5 *******************************
      SUBROUTINE VBRI45(teta,fi,vint,psite)

C     Interpolation over bridge with three values of theta
C     (0, pi/4, pi/2) and three values of phi (0, pi/4, pi/2).
C     VINT: first index= pot, d/dZ, d/dr.
C        second index= bridge site geometry (1: theta=0; 
C     2: theta=pi/2,phi=0; 3: theta=pi/2,phi=pi/2;
C     4: theta=pi/4, phi=0; 5: theta=pi/4, phi=pi/2;
C     6: theta=pi/2, phi=pi/4; 7: theta=pi/4, phi=pi/4)
C     On output, PSITE contains potential and derivatives for 
C     the given value of TETA and FI.
C     The index corresponds to: 1: function, 2: d/dZ, 3: d/dr,
C     4: d/d(phi), 5: d/d(theta).

      IMPLICIT NONE

      INTEGER ids
      PARAMETER (ids=7)
               
      INTEGER i      
      DOUBLE PRECISION vint(3,ids),psite(5),cos2fi,sin2fi,cos2te,sin2te,
     1       v0,v45,v90,dv90dfi,dv45dfi,teta,fi,cos4fi,sin4fi,cos4te,
     2       sin4te
      
      cos2fi=DCOS(2.D0*fi)
      sin2fi=DSIN(2.D0*fi)
      cos4fi=DCOS(4.D0*fi)
      sin4fi=DSIN(4.D0*fi)
      cos2te=DCOS(2.D0*teta)
      sin2te=DSIN(2.D0*teta)
      cos4te=DCOS(4.D0*teta)
      sin4te=DSIN(4.D0*teta)
 
      DO i=1,3
        v0=vint(i,1)
        v90=0.25D0*( (vint(i,2)+vint(i,3)+2.D0*vint(i,6))
     1      +2.D0*(vint(i,2)-vint(i,3))*cos2fi 
     2      + (vint(i,2)+vint(i,3)-2.D0*vint(i,6))*cos4fi )
        v45=0.25D0*( (vint(i,4)+vint(i,5)+2.D0*vint(i,7))
     1      +2.D0*(vint(i,4)-vint(i,5))*cos2fi 
     2      + (vint(i,4)+vint(i,5)-2.D0*vint(i,7))*cos4fi )

        IF (i.EQ.1) THEN
C                                                            Potential
          psite(1)=0.25D0*( (v0+v90+2.D0*v45)+2.D0*(v0-v90)*cos2te 
     1                    + (v0+v90-2.D0*v45)*cos4te )
C                                  First derivatives over fi and theta
          dv90dfi=-(vint(1,2)-vint(1,3))*sin2fi
     1            - (vint(i,2)+vint(i,3)-2.D0*vint(i,6))*sin4fi 
          dv45dfi=-(vint(1,4)-vint(1,5))*sin2fi
     1	          - (vint(i,4)+vint(i,5)-2.D0*vint(i,7))*sin4fi 
          psite(4)=0.25d0*(dv90dfi*(1.D0-2.d0*cos2te+cos4te)
     1	           + 2.d0*dv45dfi*(1.D0-cos4te) )
          psite(5)=-(v0-v90)*sin2te - (v0+v90-2.D0*v45)*sin4te

        ELSE
C                                           First derivatives over r,Z
          psite(i)=0.25D0*( (v0+v90+2.D0*v45)+2.D0*(v0-v90)*cos2te 
     1	                    + (v0+v90-2.D0*v45)*cos4te )

        END IF

      END DO

      END
