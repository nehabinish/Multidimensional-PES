C     Calculate the 6d potential for a given value of coordinates

      IMPLICIT NONE
      
      DOUBLE PRECISION pi
      PARAMETER(pi=3.141592653589793D0)

      INTEGER idelta
      DOUBLE PRECISION delta
      COMMON/SURF/delta,idelta

      LOGICAL switch
      DOUBLE PRECISION x,y,z,r,teta,fi,xcc,ycc,zcc,xa,xb,ya,yb,za,zb,
     1                 potv,dvdxa,dvdya,dvdza,dvdxb,dvdyb,dvdzb
      CALL LECTATOM
      CALL LECTMOL45
      CALL LECTMOL5D
      
      WRITE(*,*) 'X,Y (units of delta),Z,r (Ang.),teta,fi (deg)?'      
      READ(*,*) X,Y,Z,r,teta,fi

      WRITE(*,*) 'z:',Z,' (Ang.)'
      WRITE(*,*) 'x:',X,' y:',Y,' (*delta)' 
      WRITE(*,*) 'r:',r,' (Ang.), theta:',teta,', phi:',fi,' (deg)'
 
      xcc=x*delta
      ycc=y*delta
      zcc=z
      teta=teta*pi/180.D0
      fi=fi*pi/180.D0
      xa=xcc-r*dsin(teta)*dcos(fi)/2.D0
      xb=xcc+r*dsin(teta)*dcos(fi)/2.D0
      ya=ycc-r*dsin(teta)*dsin(fi)/2.D0
      yb=ycc+r*dsin(teta)*dsin(fi)/2.D0
      za=zcc-r*dcos(teta)/2.D0
      zb=zcc+r*dcos(teta)/2.D0

      switch=.FALSE.
      CALL INTMOL45(xa,ya,za,xb,yb,zb,xcc,ycc,zcc,r,teta,fi,
     1          potv,dvdxa,dvdya,dvdza,dvdxb,dvdyb,dvdzb,switch)

      IF(switch) THEN
        WRITE(*,*) 'Values outside allowed range'
      ELSE	
        WRITE(*,*) 'Potential: ',potv,' (eV)'
        WRITE(*,*) 'Derivatives:'
        WRITE(*,*) '/xa,ya,za',dvdxa,dvdya,dvdza,' (eV/Ang)'
        WRITE(*,*) '/xb,yb,zb',dvdxb,dvdyb,dvdzb,' (eV/Ang)'
      END IF
      WRITE(*,*) 'teta',teta,'fi',fi,'delta',delta  
      END
