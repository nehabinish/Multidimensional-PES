C     Calculation of potential at a given point

      IMPLICIT NONE

      DOUBLE PRECISION delta
      COMMON/SURF/delta
      
      LOGICAL fin
      INTEGER I,Ntot
      DOUBLE PRECISION xa,ya,za,v,dvdx,dvdy,dvdz,Zmin,Zmax
      
      CALL LECTATOM

      WRITE(*,*) 'xa,ya (units of delta)?'
      READ(*,*) xa,ya
      WRITE(*,*) 'Ntot,Zmin,Zmax (angstrom)'
      READ(*,*) Ntot,Zmin,Zmax
      
      xa=xa*delta
      ya=ya*delta
      
      OPEN(UNIT=1, FILE='SEPN.dat', STATUS='UNKNOWN')

      DO I=1,Ntot
        za=Zmin+DFLOAT(I-1)*((Zmax-Zmin)/(DFLOAT(Ntot)))
                   
      CALL AT3D2(xa,ya,za,v,dvdx,dvdy,dvdz,fin)
      
      WRITE(1,*) za,v
c      WRITE(*,*) 'Potential: ',v,' (eV)'
c      WRITE(*,*) 'Derivatives: , X:',dvdx,' Y: ',dvdy, ' Z :',dvdz,
c     1            ' (eV/Ang)'
      ENDDO
      END
