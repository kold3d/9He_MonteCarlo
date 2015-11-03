      SUBROUTINE GUOUT
C----------------------------------------------------------------
C-     Makes an output
C-      03-21-00
C----------------------------------------------------------------   

	PARAMETER (NVDIM=1)
        PARAMETER (NHDIM=3)
        PARAMETER (NHMAX=1000)
	real E1,X1,Y1,RAD1,E2,X2,Y2,RAD2
	INCLUDE 'gcsets.inc'
	include 'gcbank.inc'
	include 'gckine.inc'
	include 'gcflag.inc'
	integer SETNM,DETNM

	COMMON /event/ E1,X1,Y1,RAD1,E2,X2,Y2,RAD2,anincm

	DIMENSION NUMVS(NVDIM),NMBV(NVDIM,NHMAX),HITS(NHDIM,NHMAX)
        DIMENSION ITRA(NHMAX)

        SETNM = 4HPSIS
        DETNM = 4HS2E1
        
	CALL GFHITS(SETNM,DETNM,NVDIM,NHDIM,NHMAX,0,0,ITRA,NMBV,
     &              HITS,NHITS)

	IF (NHITS.EQ.(NHMAX+1)) THEN
	  print *,'Number of hits in dE detector is > NHMAX'
          print *,'Change NHMAX parameter in guout.for.'
        ENDIF

        Etot = 0. 
	E1 = 0.
	X1 = 0.
	Y1 = 0.	
C	print *,'Number of hits in proton detector',NHITS
	  DO i = 1,NHITS
	    X1 = X1 + HITS(1,i)
	    Y1 = Y1 + HITS(2,i)
	    E1 = E1 + HITS(3,i)
          ENDDO
	E1 = E1*1000.
	X1 = X1/NHITS
	Y1 = Y1/NHITS
	RAD1 = SQRT(X1*X1 + Y1*Y1)
	
C        print *,dE,'=dE'
	
        SETNM = 4HHSIS
        DETNM = 4HS2E2
        
	CALL GFHITS(SETNM,DETNM,NVDIM,NHDIM,NHMAX,0,0,ITRA,NMBV,
     &              HITS,NHITS)

	IF (NHITS.EQ.(NHMAX+1)) THEN
	  print *,'Number of hits in E detector is > NHMAX'
          print *,'Change NHMAX parameter in guout.for.'
        ENDIF

	E2 = 0.
	X2 = 0.
	Y2 = 0.
C	print *,'Number of hits in HI detector',NHITS
	  DO i = 1,NHITS
	    X2 = X2 + HITS(1,i)
	    Y2 = Y2 + HITS(2,i)
	    E2 = E2 + HITS(3,i)
          ENDDO
	E2 = E2*1000.
	X2 = X2/NHITS
	Y2 = Y2/NHITS
	RAD2 = SQRT(X2*X2 + Y2*Y2)
	
C	print *,E,'=E'
     
C	Etot = E1+E2

        IF ((E1.gt.0.1).AND.(E2.gt.0.1)) THEN
C        	print *,Etot,E1,X1,Y1,RAD1
	  CALL HFNT(1)
C	  CALL HFILL(10,E,dE,1.)
C          CALL HFILL(11,dE,0.,1.)
C          CALL HFILL(12,E,0.,1.)
C          CALL HFILL(13,Etot,0.,1.)
        ENDIF
	

	JK = LQ(JKINE-1)
	JKU = LQ(JK-1)
	
	Esam = Q(JKU+1)
	Xhit = Q(JKU+2)
	Yhit = Q(JKU+3)
	
C	CALL HFILL(10,Esam,0.,1.)
C	CALL HFILL(11,Xhit,Yhit,1.) 
	
	E1  = 0.
	E2  = 0.
	X1  = 0.
	Y1  = 0.
	X2  = 0.
	Y2  = 0.
	
	anincm = 0.
	Esam = 0. 	
  
      END
