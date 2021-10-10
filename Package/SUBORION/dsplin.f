C**************************** D S P L I N ***************************
C        Interpolation package. Given a set of n couples x(i),y(i),
C     DSPLIN defines a cubic spline function F(x) such that:
C        1)- F(x(i)) = y(i)
C        2)- F(x), F'(x), F''(x) are continuous in the interval
C     [x(1),x(n)].
C        Required subroutine: TRIDIA, BLKSPL.
C
C********************************************************************
C
      SUBROUTINE DSPLIN(n,x,y,cm,cm1,ic1,cmn,icn,alpha,beta,b)
C
C        Author: A. Salin   Version 3.2  20/1/86 - 10/6/97     
C
C     ****SUBROUTINE DSPLIN: calculates the vector F''(x(i)) which
C     defines the cubic spline function F(x).
C          N= number of pivots x(i).
C          X= array of abcissae. Should be stored by increasing or    
C     decreasing order.
C          Y= array of values of Y(i).
C          CM= array of F''(x(i)) with dimension equal to that of X
C     and Y.      
C          IC1,CM1,ICN,CMN: define the condition at x(1) and x(n)
C     If IC1=0, function is the same in first two intervals.
C     If IC1=1, second derivative at x(1) given in CM1
C     If IC1=2, first derivative at x(1) given in CM1
C     IF IC1>2, zero second derivative at x(1)
C          Similar definitions for ICN and CMN around x(n)
C          ALPHA, BETA, B are working arrays of dimension N at least.
C     For cyclic spline use the program TRISPL.
C
C     ****DPLCOF: the interpolated function at R is:
C       F=C(1)*CM(K-1)+C(2)*CM(K)+C(3)*Y(K-1)+C(4)*Y(K)
C       where C is obtained from:
C       CALL DPLCOF(n,x,r,c,k)
C       The vector C should be of dimension 4.
C
C     ****Functions DPL,DPLP,DPLP2: calculate respectively the function
C     F(x), its first or second derivative: 
C             N= number of pivots
C             X,Y,M= same as X,Y,CM in DSPLIN.
C             T= value for which the function (or its derivative)
C     must be determined.
C     WARNING: when T is outside the interval [X(1),X(N)], F(X) 
C     *******  (resp. F', F") is determined by linear extrapolation
C              using the value of F (resp. F', F") for the first two  
C              (last two) points in the interval.
C
C     ****Subroutine DINITI et DINTSP: calculate the integral of the
C     function F(x) from X(1) to T. The subroutine DINITI should first
C     be called before the first CALL DINTSP concerning a given
C     function F(x). DINITI calculates the array (CI) of values of the
C     integral of F(x) from X(1) to all X(i).
C
C        Arguments of DINITI:
C             N,X,Y,CM: as in subroutine DSPLIN.
C             CI: array of dimension at least N.  
C        Arguments of DINTSP:
C             N,X,Y,CM,CI: as in DINITI
C             T: value of the upper limit of the integration (X(1).LT.T.
C     LE.X(N)).
C             B: value of the integral.
C       The calculations in DINTSP are nearly as rapid as for the deter
C     mination of the function F(x) by DPL.
C
C     Parameters of COMMON/KODSPL/:
C     -----------------------------------------------------
C       -KOD= after execution of one subroutine of the package, KOD
C     takes the value:
C              *0: no error
C              *-1: values of X(i) are not stored in correct order.
C              *1: T outside the interval [X(1),X(N)]. Linear extrapo-
C     lation done using X(1) and X(2).
C              *2: T outside the interval [X(1),X(N)]. Linear extrapo-
C     lation done using X(n-1) and X(n).
C       -K= after the execution of DPL, DPLP or DPLP2, K is such that
C     T lies in the interval [X(k-1),X(k)].
C       -IW= error messages are printed on unit IW unless IW.LE.0
C     Default: 6.
C       -ILIGNE: when T is outside the interval [X(1),X(n)], DPL, DPLP
C     and DPLP2 print a message if ILIGNE.GT.0. The initial value of
C     ILIGNE (defined by DATA) is 100. For every extrapolation, the
C     value is decreased by 1.
C
C     ***************************************************************
C
      IMPLICIT REAL*8 (a-h,o-z)
      PARAMETER (six=6.D0,douze=12.D0,zero=0.D0)
C
      DIMENSION x(n),y(n),cm(n)
      DIMENSION alpha(n),beta(n),b(n)
C     
      COMMON/KODSPL/kod,k,iw,iligne
C
   91 FORMAT('Error in the data for DSPLIN',/,'Check abcissae:',
     1       1PD15.8,1X,' and ',D15.8,/,'Program stopped')
C
C                                           Contour condition at X(1)
      IF(ic1.EQ.0) THEN
        s1=zero
        fac1=(x(2)-x(1))*(1.D0+(x(2)-x(1))/(x(3)-x(2)))/six
        fab1=-(x(2)-x(1))**2/(six*(x(3)-x(2)))
      ELSE IF(ic1.EQ.1) THEN
        s1=-cm1*(x(2)-x(1))/six
        fac1=zero
        fab1=zero
      ELSE IF(ic1.EQ.2) THEN
        s1=0.5D0*(cm1-(y(2)-y(1))/(x(2)-x(1)))
        fac1=-(x(2)-x(1))/douze
        fab1=zero
      ELSE
        WRITE(iw,*) 'DSPLIN: unknown condition - use default'
        s1=zero
        fac1=zero
        fab1=zero
      END IF
C                                           Contour condition at X(n)
      IF(icn.EQ.0) THEN
        sn=zero
        facn=(x(n)-x(n-1))*(1.D0+(x(n)-x(n-1))/(x(n-1)-x(n-2)))/six
        fabn=-(x(n)-x(n-1))**2/(six*(x(n-1)-x(n-2)))
      ELSE IF(icn.EQ.1) THEN
        sn=-cmn*(x(n)-x(n-1))/six
        facn=zero
        fabn=zero
      ELSE IF(icn.EQ.2) THEN
        sn=-0.5D0*(cmn-(y(n)-y(n-1))/(x(n)-x(n-1)))
        facn=-(x(n)-x(n-1))/douze
        fabn=zero
      ELSE
        WRITE(iw,*) 'DSPLIN: unknown condition - use default'
        sn=zero
        facn=zero
        fabn=zero
      END IF
C
      kod=0
      c=x(2)-x(1)
      e=(y(2)-y(1))/c
      DO i=3,n
        i2=i-2
        a=c
        c=x(i)-x(i-1)
        IF(a*c.LE.zero) THEN   
          kod=-1
          IF(iw.GT.0) WRITE(iw,91) x(i-1),x(i)
          STOP
          END IF
        alpha(i2)=(a+c)/3.D0
        beta(i2)=c/six
        cm(i2)=beta(i2)
        d=e
        e=(y(i)-y(i-1))/c
        b(i2)=e-d
      END DO
C
      b(1)=b(1)+s1
      b(n-2)=b(n-2)+sn
      alpha(1)=alpha(1)+fac1
      alpha(n-2)=alpha(n-2)+facn
      cm(n-3)=cm(n-3)+fabn
      beta(1)=beta(1)+fab1
C                                            Solve tridiagonal system
      CALL TRIDIA(alpha,cm,beta,b,cm(2),n-2)
C     
      IF(ic1.EQ.0) THEN
        cm(1)=cm(2)*(1.D0+(x(2)-x(1))/(x(3)-x(2)))
     1        -cm(3)*(x(2)-x(1))/(x(3)-x(2))
      ELSE IF(ic1.EQ.1) THEN
        cm(1)=cm1
      ELSE IF(ic1.EQ.2) THEN
        cm(1)=-six*s1/(x(2)-x(1))-cm(2)/2.D0
      ELSE
        cm(1)=zero
      END IF
C      
      IF(icn.EQ.0) THEN
        cm(n)=cm(n-1)*(1.D0+(x(n)-x(n-1))/(x(n-1)-x(n-2)))
     1        -cm(n-2)*(x(n)-x(n-1))/(x(n-1)-x(n-2))
      ELSE IF(icn.EQ.1) THEN
        cm(n)=cmn
      ELSE IF(icn.EQ.2) THEN
        cm(n)=-six*sn/(x(n)-x(n-1))-cm(n-1)/2.D0
      ELSE
        cm(n)=zero
      END IF
C
      k=2
      RETURN
      END
