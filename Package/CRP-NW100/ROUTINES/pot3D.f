      BLOCK DATA PES
      CHARACTER*50 pesname
      COMMON/PESDEF/pesname
      DATA pesname/'N/W(100)..........................................'/
      END
C************************ P O T 3 D ******************************
      SUBROUTINE POT3D(z,v,dvdxa,dvdya,dvdza,fin)
      IMPLICIT NONE
C
      DOUBLE PRECISION evua,angua
      PARAMETER(evua=27.2D0,angua=0.529177249D0)

      LOGICAL fin 
      DOUBLE PRECISION z(3),v,dvdxa,dvdya,dvdza,xa,ya,za

      xa=z(2)*angua
      ya=z(3)*angua
      za=z(1)*angua
             
      fin=.false.
      CALL AT3D2(xa,ya,za,v,dvdxa,dvdya,dvdza,fin)
      
      V=V/evua
      DVDXA=DVDXA*angua/evua 
      DVDYA=DVDYA*angua/evua 
      DVDZA=DVDZA*angua/evua 

      END
