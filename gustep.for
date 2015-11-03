      SUBROUTINE GUSTEP
C----------------------------------------------------------------
C-     User command what to do next.
C----------------------------------------------------------------

       INCLUDE 'gckmax.inc'
       INCLUDE 'gcking.inc'
       INCLUDE 'gckine.inc'
       INCLUDE 'gcbank.inc'
       INCLUDE 'gctmed.inc'
       INCLUDE 'gctrak.inc' 
       INCLUDE 'gcvolu.inc'
       INCLUDE 'gcflag.inc'
       INCLUDE 'gcsets.inc'
       INCLUDE 'gctime.inc'
       INCLUDE 'udata.inc'
       
C-->    N-Body Monte Carlo common blocks
       real mss,cmang,WTMAX
       LOGICAL Flag
	COMMON /GENIN/ NP,TECM,MSS(18),KGENEV
        COMMON /GENOUT/ PCM(5,18),WT
	COMMON /DDD/ Flag

       INTEGER KGENEV
       INTEGER REC1,REC2,REC3,REC4,iedis,nubuf
       LOGICAL anYACK,enYACK
       INTEGER*4 enm,denm,voln1,voln2,PrNm1,PrNm2
       REAL urn(100),grn(100),M1,M2,M3,M4,M5,M6,weight
       REAL S(4),A(4),X(4),HITS(3)
       real ubuf(3)
       REAL UETOT,MPim,PPim,DeltaM,MTar,M8He,Mneut
    

C-------------------------------------------------------------

C       voln1 = 4HGTVO
C       voln2 = 4HGMVO

        voln1 = 4HTARG
	voln2 = 4HDUMP

       IF (IPART.EQ.9) THEN
          IF (voln1.EQ.NAMES(NLEVEL)) THEN

                   NP = 2
                   REC3 = 0
                   REC4 = 0
                   M5 = 0.
                   M6 = 0.
		   KGENEV = 1

       GO TO (10,20,30)RType

  10   CONTINUE
C-->    For Elastic Scattering       
	REC1 = 14
	REC2 = 132
        GO TO 100
  
  20   CONTINUE
C-->  Population of the narrow 9He state at 1.13 MeV 
        REC1 = 8 
	REC2 = 134
	GO TO 100
	
  30   CONTINUE
C-->  Population of the three particle continuum Pi+ + n + 8He
	NP   = 3
	REC1 = 8 
	REC2 = 13 
	REC3 = 133
	GO TO 100
	
  100   CONTINUE
                   JPA = LQ(JPART-9)
                   M1 = Q(JPA + 7)
                   JPA = LQ(JPART-64)
                   M2 = Q(JPA + 7)
                   JPA = LQ(JPART-REC1)
                   M3 = Q(JPA + 7)
                   JPA = LQ(JPART-REC2)
                   M4 = Q(JPA +7)
	           IF (REC3.NE.0) THEN
                        JPA = LQ(JPART-REC3)
			M5 = Q(JPA +7)
                   ENDIF
	           IF (REC4.NE.0) THEN
                        JPA = LQ(JPART-REC4)
			M6 = Q(JPA +7)
                   ENDIF

             JK = LQ(JKINE-ITRA)
	     JKU = LQ(JK-1)
C	     IF (ABS(Q(JKU+1)*(M1+M2)/M2*0.001-GEKIN).LT.DESTEP) THEN
             
	     IF (ABS(Q(JKU+1)*0.001-GEKIN).LT.DESTEP) THEN
C	          print *,Q(JKU+1),DESTEP,GEKIN
C	          print *,'I am there',IEVENT 
                  CALL RANMAR(urn,100)
                  CALL RNORML(grn,100)
	     
C-->   General part for any binary reaction (CERN Monte Carlo Phase Space)
C-->    Output in common block GENOUT
                  TECM = SQRT(M1*M1 + M2*M2 + 2*M2*GETOT)
                   IF (TECM.gt.(M3+M4+M5+M6+0.0001)) THEN
                     MSS(1) = M3
                     MSS(2) = M4
                     MSS(3) = M5
                     MSS(4) = M6

C Generating recoil particles with GENBOD according to restrictions
C given by angular distribution and Fermi energy function. 
                     anYACK = .TRUE.
                     enYACK = .TRUE.
                     in = 0
	             WTMAX = 0.
		     DO WHILE (anYACK.OR.enYACK)
                       in = in +1
                       CALL GENBOD
C			print *,'Call GENBOD'
C                       IF (WT.GT.WTMAX) WTMAX=WT
                       IF (urn(10+in).LE.WT) enYACK = .FALSE.
C Angular distribution is checking for particle number two and
C only if two particle kinematics is used.
                       IF (NP.EQ.2) THEN
                         cmang = SQRT( (PCM(1,2)**2. + 
     &                           PCM(2,2)**2.)/
     &                           (PCM(1,2)**2. + PCM(2,2)**2.
     &                           + PCM(3,2)**2.))
                         cmang = ABS(ASIN(cmang))
                         IF (PCM(3,2).LT.0.) cmang =
     &                                       3.14159 - cmang
                         weight = gtwght(cmang)
                         IF (urn(50+in).LE.weight) anYACK = .FALSE.
                       ELSE
                         anYACK = .FALSE.
                       ENDIF
                      IF (IDEBUG.EQ.1) THEN
                         print *,'Weight from GENBOD ',WT
C                         print *,'Weight from gtwght ',weight
C                         print *,'Current CM angle',cmang
                         print *,'Event number ',IEVENT
C                         pause
                      ENDIF
                     IF (in.EQ.40) in = 0
C                     CALL HFILL(10,cmang,0.,1.)
                     ENDDO
		     
		     nubuf = 3
		     ubuf(1) = VECT(1)
		     ubuf(2) = VECT(2)
		     ubuf(3) = VECT(3)
		     IADR = 1
		     CALL GSKINU(ITRA,nubuf,ubuf,IADR)
                     anincm = cmang
C-------------------------------------------------------------------

                     S(1) = 0.-VECT(4)*VECT(7)
                     S(2) = 0.-VECT(5)*VECT(7)
                     S(3) = 0.-VECT(6)*VECT(7)
C                     S(4) = VECT(7)*VECT(7)/2./(M3+M4) + M3 + M4
C	              S(4) = SQRT(M1**2 + M2**2 + 2*M2*GETOT) 
	             S(4) = GETOT + M2
		     EALOR = 0
                       DO i = 1,NP
                          DO k = 1,4
                            A(k) = PCM(k,i)
                          ENDDO

                          CALL LOREN4(S,A,X)

                          DO k = 1,4
                            GKIN(k,i) = X(k)
                          ENDDO
 
                          IF (i.eq.1) GKIN(5,i) = REC1
                          IF (i.eq.2) GKIN(5,i) = REC2
                          IF (i.eq.3) GKIN(5,i) = REC3
                          IF (i.eq.4) GKIN(5,i) = REC4
                          TOFD(i) = 0
                          
			  DO k = 1,3
                            GPOS(k,i) = VECT(k)
                          ENDDO

C			  EALOR = EALOR + X(4)

		       ENDDO

C			print *,'Total Energy after Lorenz travsform ',EALOR

                         NGKINE = NP
                         
			 CALL GSKING(0)
			
C			print *,' ' 
C			  DO i = 1,NGKINE 
C	  		    print *,GKIN(1,i),GKIN(2,i),GKIN(3,i),GKIN(4,i),GKIN(5,i)
C                          ENDDO
C			print *,' ' 

                         ISTOP = 1

	           ELSE 
                         ISTOP = 1
		   ENDIF
                   
              ENDIF
            
           ENDIF

       ELSE 
	  CALL GSKING(0)
       ENDIF

       IF (IDEBUG.EQ.1) THEN
           CALL GDEBUG
C	   CALL GPJXYZ
C           CALL GPCXYZ
C	   pause
       ENDIF
       

       IF (IPART.EQ.8) THEN
         denm = 4HMDMS 
         IF (denm.EQ.NAMES(NLEVEL)) THEN
          print *,'This is a hit!!!'
          CALL GPCXYZ
	  print *,VECT(4)*VECT(7), VECT(5)*VECT(7), VECT(6)*VECT(7)
          print *,'Event number ',IEVENT
             JK = LQ(JKINE-ITRA)
	     JKU = LQ(JK-1)
             JPA = LQ(JPART-9)
             MPim = Q(JPA + 7)
             JPA = LQ(JPART-64)
             MTar = Q(JPA + 7) 
C	     UETOT = Q(JKU+1)*0.001 + MPim 
C    178.7 MeV Pi- hardcoded
	     UETOT = 178.7*0.001 + MPim 
             PPim = SQRT(UETOT*UETOT-MPim*MPim)
             DeltaM = SQRT((UETOT + Mtar - GETOT)**2.-
     &                 (VECT(4)*VECT(7))**2.-
     &                 (VECT(5)*VECT(7))**2.-
     &                 (PPim-VECT(6)*VECT(7))**2.)
             JPA = LQ(JPART-133)
             M8He = Q(JPA + 7) 
             JPA = LQ(JPART-13)
             Mneut = Q(JPA + 7) 
             DeltaM = DeltaM - M8He - Mneut 
   	     print *,'Missing mass ',DeltaM
          CALL HFILL(10,DeltaM*1000.,0.,1.0)
	  Flag=.TRUE.
          ISTOP = 1
         ENDIF
       ENDIF
       
      IF ((IPART.EQ.134).AND.(Flag)) THEN
	CALL GPCXYZ
	print *,VECT(4)*VECT(7), VECT(5)*VECT(7), VECT(6)*VECT(7)
        print *,'Event number ',IEVENT
        Flag = .FALSE.
        ISTOP = 1 
      ENDIF

      RETURN
      END

      
        REAL FUNCTION RNDM(dummy)
        real dummy,urn(5)
      
	CALL RANMAR(urn,3)
        RNDM = urn(2)
      
	RETURN
        END

