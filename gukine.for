      SUBROUTINE GUKINE
C----------------------------------------------------------------
C       Initialization of initial kinematics
C      Definition of beam particle. 
C      11-27-01 Dubna experiment
C----------------------------------------------------------------
        
      DIMENSION VERT(3),PLAB(3)
      integer NUBUF
      real urn(10),Bphi,BeamR,grn(10),ubuf(10)
      real Esample
      integer*4 a(3),b(3)
      INCLUDE 'udata.inc'
      INCLUDE 'gcflag.inc'
      INCLUDE 'gcbank.inc'
      SAVE b

C-->     Parameters of the initial beam

      PARAMETER (D2R    = 0.0174533)

      CALL RANMAR(urn,5)
      CALL RNORML(grn,5)

      call itime(a)
      IF (a(2).NE.b(2)) THEN
        PRINT *,' The time is: ',(a(i),i=1,3)
        b(1)=a(1)
        b(2)=a(2)
        b(3)=a(3)
      ENDIF

       Bphi = 2.*3.14159*urn(1)
       BeamR= ABS(dRad(1)+dRad(2)*grn(1))

       VERT(1) = BeamR * Cos(Bphi)
       VERT(2) = BeamR * Sin(Bphi)

	IF (CHAR(TState).EQ.'G') THEN
C-->	This is for Chubarian chamber ONLY!!
             VERT(3) = -32.42
        ELSEIF (CHAR(TState).EQ.'S') THEN
             VERT(3) = -18.
        ELSEIF (CHAR(TState).EQ.'L') THEN
             VERT(3) = 0.
        ENDIF

       CALL GSVERT(VERT,0,0,0,0,NVERT)

C       JPA = LQ(JPART-132)
C       AMASS = Q(JPA+7) 
       JPA = LQ(JPART-9)
       AMASS = Q(JPA+7) 
       
       BENER = BeamE + dEner*grn(2)
       
       UBUF(1) = Esample()

       IF (BENER.LE.0.) THEN
          PRINT *,'Energy of the incident particle below zero!!'
          PAUSE 
       ENDIF 

       BPTOT  = SQRT(BENER**2.+2.*AMASS*BENER)
C       BPTOT = BENER
C       print *,AMASS
       
       BTheta = D2R*ABS(theav+dtheta*grn(3))

        IF (CHAR(TState).EQ.'G') THEN
	   PLAB(1) = BPTOT*SIN(BTheta)*COS(Bphi) 
           PLAB(2) = BPTOT*SIN(BTheta)*SIN(Bphi)
           PLAB(3) = BPTOT*COS(BTheta)
	ELSEIF (CHAR(TState).EQ.'S') THEN   
           PLAB(1) = BPTOT*SIN(BTheta)*COS(Bphi) 
           PLAB(2) = BPTOT*SIN(BTheta)*SIN(Bphi)
           PLAB(3) = BPTOT*COS(BTheta)
	ELSEIF (CHAR(TState).EQ.'L') THEN
	   PLAB(1) = BPTOT*SIN(BTheta)*COS(Bphi) 
           PLAB(2) = BPTOT*SIN(BTheta)*SIN(Bphi)
           PLAB(3) = BPTOT*COS(BTheta)
	ENDIF
	
       NUBUF = 1	
C       CALL GSKINE(PLAB,132,1,UBUF,NUBUF,NT)
       CALL GSKINE(PLAB,9,1,UBUF,NUBUF,NT)
       
      IF (IDEBUG.EQ.1) THEN
C        CALL GDEBUG
        CALL GPVERT(NVERT)
        CALL GPKINE(NT)
      ENDIF
      
      RETURN
      
      END

	REAL FUNCTION Esample()
	  
	  include 'udata.inc'
	  real urn(10),grn(10)
	   
	   CALL RANMAR(urn,5)
           CALL RNORML(grn,5)
	  
	  IF (int(enkey(1)).eq.1) THEN
	     print *,'Do not know what to do. No excitation function'
	     stop
	  ELSE
	     Esample = enkey(2)+grn(3)*enkey(3)/2.34
	  ENDIF
	  
	  RETURN
	
	END 
