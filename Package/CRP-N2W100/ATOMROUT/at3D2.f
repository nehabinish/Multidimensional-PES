C************************ A T 3 D 2 *******************************

C     Calculation for six sites: top, bridge, hole, bridge-hole, 
C     top-hole and top-bridge 

C     Calculation of the 3D potential as well as its
C     first derivatives using the Fourier series expansion.
C     The atom position is xa, ya, za. 

C     All distances are in angstroems and energies in eV.
   
C     FOURI3D: comes from the solution of the linear system associated
C              with the Fourier expansion and reference sites. It is
C              calculated in FOUR3D which is called within the input 
C              data reading routine because it must be called once before
C              the first call to AT3D2
C      Version 1.0           A. Salin            09/12/02
C
      SUBROUTINE AT3D2(xa,ya,za,v,dvdx,dvdy,dvdz,fin)
      IMPLICIT NONE
      
      INTEGER idf3,idz
      PARAMETER(idf3=6,idz=9999)
C     IDF3: number of reference sites used for the interpolation.
C     IDZ: maximum number of Z values in data.

      DOUBLE PRECISION fouri3D
      COMMON/MAFOU3D/fouri3D(idf3,idf3)
C     See subroutine FOUR3D.

      INTEGER k,ktop
      COMMON/CONSERV/k,ktop
C     Used to conserve the previous index K for the spline
C     interpolation.

C     DELTA: nearest neighbour distance on the surface (in Angstroems)
      INTEGER idelta
      DOUBLE PRECISION delta
      COMMON/SURF/delta,idelta

      LOGICAL fin
      INTEGER j,ksit
      DOUBLE PRECISION xa,ya,za,v,dvdx,dvdy,dvdz,x0c,y0c,vrep,dvrepdx,
     1       dvrepdy,dvrepdz
      DOUBLE PRECISION a(idf3),fxy(idf3),dfdx(idf3),dfdy(idf3),b(idf3),
     1       vs(idf3),dvsdz(idf3),c(4),dc(4)

C     ----------------------------------------------------------------
C     Data for top
      INTEGER n1p
      DOUBLE PRECISION z1p,vz1p,cvz1p
      COMMON/SITIO1P/z1p(idz),vz1p(idz),cvz1p(idz),n1p
      INTEGER n1m
      DOUBLE PRECISION z1m,vz1m,cvz1m
      COMMON/SITIO1M/z1m(idz),vz1m(idz),cvz1m(idz),n1m

C     Data for bridge
      INTEGER n2
      DOUBLE PRECISION z2,vz2,cvz2
      COMMON/SITIO2/z2(idz),vz2(idz),cvz2(idz),n2

C     Data for threefold
      INTEGER n3
      DOUBLE PRECISION z3,vz3,cvz3
      COMMON/SITIO3/z3(idz),vz3(idz),cvz3(idz),n3

C     Data for bridge-threefold site
      INTEGER n4
      DOUBLE PRECISION z4,vz4,cvz4
      COMMON/SITIO4/z4(idz),vz4(idz),cvz4(idz),n4

C     Data for top-hole site
      INTEGER n5
      DOUBLE PRECISION z5,vz5,cvz5
      COMMON/SITIO5/z5(idz),vz5(idz),cvz5(idz),n5

C     Data for top-bridge site
      INTEGER n6
      DOUBLE PRECISION z6,vz6,cvz6
      COMMON/SITIO6/z6(idz),vz6(idz),cvz6(idz),n6

C     ----------------------------------------------------------------


C                                       Contribution from 1D potential

      IF(za.GE.z1m(1).AND.za.LE.z1p(n1p)) THEN
	IF(za.GE.0.D0) THEN
          CALL DPLCOFT(n1p,z1p,za,c,dc,ktop)
      vs(1)= c(1)*cvz1p(ktop-1) + c(2)*cvz1p(ktop)
     1     + c(3)*vz1p(ktop-1)  + c(4)*vz1p(ktop)
      dvsdz(1)= dc(1)*cvz1p(ktop-1) + dc(2)*cvz1p(ktop)
     1        + dc(3)*vz1p(ktop-1)  + dc(4)*vz1p(ktop)
        ELSE
          CALL DPLCOFT(n1m,z1m,za,c,dc,ktop)
      vs(1)= c(1)*cvz1m(ktop-1) + c(2)*cvz1m(ktop)
     1     + c(3)*vz1m(ktop-1)  + c(4)*vz1m(ktop)
      dvsdz(1)= dc(1)*cvz1m(ktop-1) + dc(2)*cvz1m(ktop)
     1        + dc(3)*vz1m(ktop-1)  + dc(4)*vz1m(ktop)
        END IF        
      ELSE
        fin=.true.
        RETURN
      END IF

      IF(za.GE.z2(1).AND.za.LE.z2(n2)) THEN
        CALL DPLCOFT(n2,z2,za,c,dc,k)
      ELSE
        fin=.true.
        RETURN
      END IF

C     ----------------------------------------------------------------
C                                   Interpolation over Z for each site

      vs(2)=c(1)*cvz2(k-1)+c(2)*cvz2(k)+c(3)*vz2(k-1)+c(4)*vz2(k)
      vs(3)=c(1)*cvz3(k-1)+c(2)*cvz3(k)+c(3)*vz3(k-1)+c(4)*vz3(k)
      vs(4)=c(1)*cvz4(k-1)+c(2)*cvz4(k)+c(3)*vz4(k-1)+c(4)*vz4(k)
      vs(5)=c(1)*cvz5(k-1)+c(2)*cvz5(k)+c(3)*vz5(k-1)+c(4)*vz5(k)
      vs(6)=c(1)*cvz6(k-1)+c(2)*cvz6(k)+c(3)*vz6(k-1)+c(4)*vz6(k)
      dvsdz(2)=dc(1)*cvz2(k-1)+dc(2)*cvz2(k)
     1        +dc(3)*vz2(k-1)+dc(4)*vz2(k)
      dvsdz(3)=dc(1)*cvz3(k-1)+dc(2)*cvz3(k)
     1        +dc(3)*vz3(k-1)+dc(4)*vz3(k)
      dvsdz(4)=dc(1)*cvz4(k-1)+dc(2)*cvz4(k)
     1        +dc(3)*vz4(k-1)+dc(4)*vz4(k)
      dvsdz(5)=dc(1)*cvz5(k-1)+dc(2)*cvz5(k)
     1        +dc(3)*vz5(k-1)+dc(4)*vz5(k)
      dvsdz(6)=dc(1)*cvz6(k-1)+dc(2)*cvz6(k)
     1        +dc(3)*vz6(k-1)+dc(4)*vz6(k)
C---------------------------------------------------------------------

C                      Calculate basis functions for Fourier expansion

      CALL FOUR3II(xa,ya,fxy,dfdx,dfdy)
      
C                         Calculation of potential and derivatives
      v=0.D0
      dvdx=0.D0
      dvdy=0.D0
      dvdz=0.D0

      DO j=1,idf3
        a(j)=0.D0
        b(j)=0.D0
        DO ksit=1,idf3
          a(j)=a(j)+fouri3D(ksit,j)*vs(ksit)
          b(j)=b(j)+fouri3D(ksit,j)*dvsdz(ksit)
        END DO
        v=v+a(j)*fxy(j)
        dvdx=dvdx+a(j)*dfdx(j)
        dvdy=dvdy+a(j)*dfdy(j)
        dvdz=dvdz+b(j)*fxy(j)
      END DO

      CALL ORIGEN100(xa,ya,delta,x0c,y0c)
      CALL MUR3D2(x0c,y0c,za,vrep,dvrepdx,dvrepdy,dvrepdz)

      v=(v+vrep)
      dvdx=(dvdx+dvrepdx)
      dvdy=(dvdy+dvrepdy)
      dvdz=(dvdz+dvrepdz)

      END
C************************ F O U R 3 I I ******************************

C     Calculates the basis functions of the Fourier expansion and 
C     their derivatives with respect to X and Y. 3D potential.

      SUBROUTINE FOUR3II(xa,ya,fxy,dfdx,dfdy)
      IMPLICIT NONE

      INTEGER idf3
      PARAMETER (idf3=6)
C     IDF3: number of reference site used for the interpolation.

      DOUBLE PRECISION pi
      PARAMETER (pi=3.141592653589793D0)

      DOUBLE PRECISION xa,ya,fac,fbigx,fbigy,s3a,
     1       s3b,s5a,s5b,s5c,s5d,s6a,s6b 
      DOUBLE PRECISION fxy(idf3),dfdx(idf3),dfdy(idf3)

      INTEGER idelta
      DOUBLE PRECISION delta
      COMMON/SURF/delta,idelta

      fac=2.D0*pi/delta
      fbigx=fac*xa
      fbigy=fac*ya
C                                                      Basis functions
C                  and derivatives with respect to crystal coordinates

      fxy(1)=1.D0
      dfdx(1)=0.D0
      dfdy(1)=0.D0

      fxy(2)=DCOS(fbigx)+DCOS(fbigy)
      dfdx(2)=-fac*DSIN(fbigx)      
      dfdy(2)=-fac*DSIN(fbigy)      

      fxy(3)=DCOS(fbigx+fbigy)+DCOS(fbigx-fbigy)
      s3a=DSIN(fbigx+fbigy)
      s3b=DSIN(fbigx-fbigy)
      dfdx(3)=-fac*(s3a+s3b)      
      dfdy(3)=-fac*(s3a-s3b)      

      fxy(4)=DCOS(2.D0*fbigx)+DCOS(2.D0*fbigy)
      dfdx(4)=-2.D0*fac*DSIN(2.D0*fbigx)      
      dfdy(4)=-2.D0*fac*DSIN(2.D0*fbigy)      

      fxy(5)=DCOS(2.D0*fbigx+fbigy)+DCOS(fbigx+2.D0*fbigy)
     #       + DCOS(2.D0*fbigx-fbigy)+DCOS(fbigx-2.D0*fbigy)
      s5a=DSIN(2.D0*fbigx+fbigy)
      s5b=DSIN(fbigx+2.D0*fbigy)
      s5c=DSIN(2.D0*fbigx-fbigy)
      s5d=DSIN(fbigx-2.D0*fbigy) 
      dfdx(5)=-fac*(2.D0*s5a+s5b+2.D0*s5c+s5d)      
      dfdy(5)=-fac*(s5a+2.D0*s5b-s5c-2.D0*s5d)

      fxy(6)=DCOS(2.D0*(fbigx+fbigy))+DCOS(2.D0*(fbigx-fbigy))
      s6a=2.D0*DSIN(2.D0*(fbigx+fbigy))
      s6b=2.D0*DSIN(2.D0*(fbigx-fbigy))
      dfdx(6)=-fac*(s6a+s6b)      
      dfdy(6)=-fac*(s6a-s6b)

      END
