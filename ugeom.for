      SUBROUTINE UGEOM
C----------------------------------------------------------
C-       Geometry of the experiment and definition
C-       of the materials. 
C-         G.V.Rogachev 11-25-01
C-     Dubna experiment. Three different targets can be used,
C-     that is liquid, gas and solid.
C-     Modified on April 19, 2013 to include Fletchers chamber	
C----------------------------------------------------------

      real Ach2,Zch2,Wch2,pladen,D2R
      real Ach4,Zch4,Wch4,metden
      real Aaram(5),Zaram(5),Waram(5),araden
      real AHavar(5),ZHavar(5),WHavar(5),havden
      real Ah2,Zh2,Wh2
      real zdet,xdet,zdet1,xdet1,DZOFF
      real fieldm,tmaxfd,stemax,deemax,epsil,stmin,ubuf
      real PAR(3),fact,orig,TRPAR(5)
C-->   Lab. Angle of detector in degrees
      PARAMETER (D2R = 0.0174533)

C-->   Distance to the detector in cm

      DIMENSION Ach2(2),Zch2(2),Wch2(2)
      DIMENSION Ach4(2),Zch4(2),Wch4(2)
      DIMENSION Ah2(2),Zh2(2),Wh2(2)
      DIMENSION ORIG(3),FACT(3),NBITSH(3),NBITSV(3)

      CHARACTER*20 NAMATE(10)
      CHARACTER*4 CHNAMH(3),CHNMSV(4)

      INCLUDE 'udata.inc'
      include 'gcflag.inc'

      DATA Ach2/12.0,1.0/
      DATA Zch2/6.0,1.0/
      DATA Wch2/1.0,2.0/

      DATA Ach4/12.0,1.0/
      DATA Zch4/6.0,1.0/
      DATA Wch4/1.0,4.0/

      DATA Ah2/1.0,1.0/
      DATA Zh2/1.0,1.0/
      DATA Wh2/1.0,1.0/

      DATA Aaram/1.0,12.0,14.0,16.0,35.0/
      DATA Zaram/1.0,6.0,7.0,8.0,17.0/
      DATA Waram/31.0,50.0,7.0,7.0,6.0/

      DATA AHavar/1.0,12.0,14.0,16.0,35.0/
      DATA ZHavar/1.0,6.0,7.0,8.0,17.0/
      DATA WHavar/31.0,50.0,7.0,7.0,6.0/

      DATA ORIG /10.,10.,0./
      DATA FACT /1000000.,1000000.,1000000./
      DATA NBITSH /32, 32, 32/
      DATA NBITSV /  6,  6, 6/


      NAMATE(1) = 'Polethelyne CH2'
      NAMATE(2) = 'Silicon'
      NAMATE(3) = 'Aluminium'
      NAMATE(4) = 'Iron'
      NAMATE(5) = 'Vacuum'
      NAMATE(6) = 'Methane gas CH4'
      NAMATE(7) = 'Liquid Hydrogen H2'
      NAMATE(8) = 'Air'
      NAMATE(9) = 'ARAMICA'
      NAMATE(10) = 'Beryllium'
      
C-->   Density of low density polythelyne
      pladen = 0.917
      metden = 7.1499E-4 * GPre / 1000.
      hydden = LDEN
      araden = 1.4

C-->   Stores the standard GEANT materials in the data structure JMATE
      CALL GMATE

      CALL GSMIXT(21,NAMATE(1),Ach2,Zch2,pladen,-2,Wch2)

      CALL GSMATE(22,NAMATE(2),28.0,14.,2.31,10.4,1.,0,0)

      CALL GSMIXT(26,NAMATE(6),Ach4,Zch4,metden,-2,Wch4)

      CALL GSMIXT(27,NAMATE(7),Ah2,Zh2,hydden,-2,Wh2)
      
      CALL GSMIXT(29,NAMATE(9),Aaram,Zaram,araden,-5,Waram)
 
C
C     Defines tracking media parameters
C

      IFIELD=0
      FIELDM=0.
      TMAXFD=0.
      STEMAX=-1.
      DEEMAX=-1.
      EPSIL=1.E-5
      STMIN=1.E-5
       
      
C    For target materials maximum fractional energy loss
C    has to be very small
 
      DEEMAX = 5.E-3
              
      CALL GSTMED(21,NAMATE(1),21,0,IFIELD,FIELDM,TMAXFD,
     +            STEMAX,DEEMAX,EPSIL,STMIN,0.,0)
     
      CALL GSTMED(26,NAMATE(8),15,0,IFIELD,FIELDM,TMAXFD,
     +            STEMAX,DEEMAX,EPSIL,STMIN,0.,0)
 
      CALL GSTMED(27,NAMATE(6),26,0,IFIELD,FIELDM,TMAXFD,
     +            STEMAX,DEEMAX,EPSIL,STMIN,0.,0)

      CALL GSTMED(5,NAMATE(10),5,0,IFIELD,FIELDM,TMAXFD,
     +            STEMAX,DEEMAX,EPSIL,STMIN,0.,0)

      DEEMAX = -1.
       
      CALL GSTMED(22,NAMATE(2),22,1,IFIELD,FIELDM,TMAXFD,
     +            STEMAX,DEEMAX,EPSIL,STMIN,0.,0)

      CALL GSTMED(23,NAMATE(3),9,0,IFIELD,FIELDM,TMAXFD,
     +            STEMAX,DEEMAX,EPSIL,STMIN,0.,0)

      CALL GSTMED(24,NAMATE(4),10,0,IFIELD,FIELDM,TMAXFD,
     +            STEMAX,DEEMAX,EPSIL,STMIN,0.,0)

      CALL GSTMED(25,NAMATE(5),16,0,IFIELD,FIELDM,TMAXFD,
     +            STEMAX,DEEMAX,EPSIL,STMIN,0.,0)

      CALL GSTMED(29,NAMATE(9),29,0,IFIELD,FIELDM,TMAXFD,
     +            STEMAX,DEEMAX,EPSIL,STMIN,0.,0)


C
C           Defines volumes
C
c
c          Definition of the Mather volume "MVOL"
C          In all cases the mother volume is very big air box. 1x1x1 m 
c
      PAR(1)= 50.
      PAR(2)= 50.
      PAR(3)= 50.
      CALL GSVOLU('MVOL','BOX ',25,PAR,3,IVOLU)

        IF (IVOLU.LE.0) THEN
           PRINT *,'Error condition has happened!'
           PRINT *,'Mother volume has not been created.'
        ENDIF

      IF (CHAR(TState).EQ.'S') THEN
        IF (CHAR(ChName).EQ.'N') THEN
C-->   Scattering Chamber in ND.
      PAR(1) = 0.
      PAR(2) = 20.5
      PAR(3) = 16.
      CALL GSVOLU('CHSH','TUBE',23,PAR,3,IVOLU) 
      PAR(1) = 0.
      PAR(2) = 20.
      PAR(3) = 15.
      CALL GSVOLU('CHVO','TUBE',25,PAR,3,IVOLU)

C-->   Target      
      PAR(1)=0.
      PAR(2)=0.8

C-->   Thickness in cm of low desity polyethylene target 9 mg/cm^2 
C      PAR(3)=0.00981/2.
      PAR(3)=0.01150/2.

      CALL GSVOLU('TARG','TUBE',21,PAR,3,IVOLU)

C-->   Target frame.
      PAR(1)=2.
      PAR(2)=2.
      PAR(3)=0.05
      CALL GSVOLU('TFRM','BOX ',23,PAR,3,IVOLU)

C-->   Hole in target frame for target
      PAR(1)=0.
      PAR(2)=0.8
      PAR(3)=0.05
      CALL GSVOLU('HOLE','TUBE',25,PAR,3,IVOLU)

C-->   dE-E telescope 
      PAR(1)=0.0
      PAR(2)=1.0
      PAR(3)=0.001
      CALL GSVOLU('dEdt','TUBE',22,PAR,3,IVOLU)
    
      PAR(1)=0.0
      PAR(2)=1.
      PAR(3)=0.05
      CALL GSVOLU('Edt ','TUBE',22,PAR,3,IVOLU)

C-->  Calimator plate on the face of the detector
      PAR(1) = 0.00
      PAR(2) = 1.5 
      PAR(3) = 0.05
      CALL GSVOLU('CALM','TUBE',24,PAR,3,IVOLU)


C-->  Hole in the Calimator on detector
C   This settings are for 3H* run.
C      TRPAR(1)=0.05
C      TRPAR(2)=0.05
C      TRPAR(3)=0.61
C      TRPAR(4)=0.84
C      TRPAR(5)=0.19
C      CALL GSVOLU('HOCL','TRD2',25,TRPAR,5,IVOLU)

       PAR(1) = 0.
       PAR(2) = RDET
       PAR(3) = 0.05
       CALL GSVOLU('HOCL','TUBE',25,PAR,3,IVOLU)

C-->   Positioning of the objects!
      CALL GSROTM(1,90.,0.,0.,-90.,90.,90.)
      CALL GSPOS('CHSH',1,'MVOL',0.,0.,0.,1,'MANY')
      CALL GSPOS('CHVO',1,'CHSH',0.,0.,0.,0,'MANY')

      CALL GSPOS('TFRM',1,'MVOL', 0. , 0. , 0. ,0,'ONLY')
      CALL GSPOS('HOLE',1,'TFRM', 0. , 0. , 0. ,0,'ONLY')
      CALL GSPOS('TARG',1,'HOLE', 0. , 0. , 0. ,0,'ONLY')
      
      
C      zdet = Apde*cos(D2R*Aang)
      zdet = Apde
      xdet = Apde*sin(D2R*Aang)

C      zdet1 = (Apde+0.3)*cos(D2R*Aang)
      zdet1 = Apde + 0.3
C      xdet1 = (Apde+0.3)*sin(D2R*Aang)
      xdet1 = xdet
      
C      zdet2 = (Apde+0.6)*cos(D2R*Aang)
      zdet2 = Apde+0.6
C      xdet2 = (Apde+0.6)*sin(D2R*Aang)
      xdet2 = xdet

C      CALL GSROTM (1,(90.+ Aang),0.,90.,90.,Aang,0.)
      CALL GSPOS('CALM',1,'MVOL', xdet , 0. , zdet ,0,'ONLY')
      CALL GSPOS('dEdt',1,'MVOL', xdet1, 0. ,zdet1 ,0,'ONLY')
      CALL GSPOS('Edt ',1,'MVOL', xdet2, 0. ,zdet2 ,0,'ONLY')
      CALL GSPOS('HOCL',1,'MVOL', xdet, 0. ,zdet,0,'ONLY')
      
C   Positioning of the calimator for 3H experiment
C      CALL GSROTM(1,180.,0.,90.,90.,90.,0.)
C      CALL GSPOS('HOCL',1,'MVOL', xdet, 0. ,zdet,1,'ONLY')

C-->   Definition of the hit
           
	  CHNMSV(1) = 'MVOL'
C	  CHNMSV(2) = 'MCHA'
C	  CHNMSV(3) = 'GMVO'
	  CHNMSV(2) = 'dEdt'

          CALL GSDET ('TLSC','dEdt',2,CHNMSV,NBITSV,1,1000,
     &                1000,iset,idet)
          CHNAMH(1) = 'dE  '
          CALL GSDETH ('TLSC','dEdt',2,CHNAMH,NBITSH,ORIG,FACT)

          CHNMSV(2) = 'Edt '
          CALL GSDET ('TLSC','Edt ',2,CHNMSV,NBITSV,1,1000,
     &                1000,iset,idet)
          CHNAMH(1) = 'E   '
          CALL GSDETH ('TLSC','Edt ',2,CHNAMH,NBITSH,ORIG,FACT)

C Neil Fletchers's chamber. Coincidence experiment

          ELSEIF (CHAR(ChName).EQ.'F') THEN

	          PAR(1) = 0.
              PAR(2) = 30.
              PAR(3) = 10.
              CALL GSVOLU('CHSH','TUBE',23,PAR,3,IVOLU)
              CALL GSROTM(1,90.,0.,0.,-90.,90.,90.)
              CALL GSPOS('CHSH',1,'MVOL',0.,0.,0.,1,'MANY')     
             
              PAR(1) = 0.
              PAR(2) = 29.5
              PAR(3) = 9.5
              CALL GSVOLU('CHVO','TUBE',25,PAR,3,IVOLU)
              CALL GSPOS('CHVO',1,'CHSH',0.,0.,0.,0,'MANY')
              
C-->   Target      
	      PAR(1)=0.
              PAR(2)=0.5
C-->   1/2 Thickness in cm of low desity polyethylene target
              PAR(3)= 1.5/pladen*0.001/2.
              CALL GSVOLU('TARG','TUBE',21,PAR,3,IVOLU)

C-->   Target frame.
              PAR(1)=2.
              PAR(2)=2.
              PAR(3)=0.05
              CALL GSVOLU('TFRM','BOX ',23,PAR,3,IVOLU)

C-->   Hole in target frame for target
              PAR(1)=0.
              PAR(2)=0.5
              PAR(3)=0.051
              
              CALL GSVOLU('HOLE','TUBE',25,PAR,3,IVOLU)      
              
              CALL GSPOS('TFRM',1,'MVOL', 0. , 0. , 0. ,0,'ONLY')
              CALL GSPOS('HOLE',1,'TFRM', 0. , 0. , 0. ,0,'ONLY')
              CALL GSPOS('TARG',1,'HOLE', 0. , 0. , 0. ,0,'ONLY')
              
C-->    S2 Detector for protons
	      PAR(1) = 1.1
	      PAR(2) = 3.5
	      PAR(3) = 0.05
	      CALL GSVOLU('S2E1','TUBE',22,PAR,3,IVOLU)
C	      CALL GSDVN('S2E1','E1',16,2)
	      
	      
C-->    S2 Detector for heavy recoil
	      PAR(1) = 1.1
	      PAR(2) = 3.5
	      PAR(3) = 0.05
	      CALL GSVOLU('S2E2','TUBE',22,PAR,3,IVOLU)
C	      CALL GSDVN('S2E2','E2',16,2)
              
	      
	      CALL GSPOS('S2E1',1,'MVOL', 0., 0., 6.,0,'ONLY')
	      
	      CALL GSPOS('S2E2',1,'MVOL', 0., 0., 24.5,0,'ONLY')
	      
C-->   Definition of the hit
           
	      CHNMSV(1) = 'MVOL'
C	  CHNMSV(2) = 'MCHA'
C	  CHNMSV(3) = 'GMVO'
	      CHNMSV(2) = 'S2E1'
	      
              CALL GSDET ('PSIS','S2E1',2,CHNMSV,NBITSV,
     &                 1,1000,1000,iset,idet)
     
              CHNAMH(1) = 'X '
              CHNAMH(2) = 'Y '
              CHNAMH(3) = 'E '
              print *,'PSIS',iset,'S2E1',idet
              
              CALL GSDETH ('PSIS','S2E1',3,CHNAMH,NBITSH,ORIG,FACT)

              CHNMSV(2) = 'S2E2'
              
              CALL GSDET ('HSIS','S2E2',2,CHNMSV,NBITSV,
     &                 1,1000,1000,iset,idet)
              
              print *,'HSIS',iset,'S2E2',idet
     
              CALL GSDETH ('HSIS','S2E2',3,CHNAMH,NBITSH,ORIG,FACT)
              
              CALL GGCLOS

C-->    MDM Chamber

	   ELSEIF (CHAR(ChName).EQ.'M') THEN

C-->   Target
        PAR(1)=2.0
        PAR(2)=2.0
C-->   1/2 Thickness in cm of beryllium target
        PAR(3)= 0.2
        CALL GSVOLU('TARG','BOX ',5,PAR,3,IVOLU)

        CALL GSPOS('TARG',1,'MVOL', 0. , 0. , 0. ,0,'ONLY')

C-->    Effective MDM bite
        PAR(1) = 1.0
        PAR(2) = 1.0
        PAR(3) = 0.1 
        CALL GSVOLU('MDMS','BOX ',22,PAR,3,IVOLU)

        CALL GSPOS('MDMS',1,'MVOL',5.3,0.0,20.0,0,'ONLY')

        CALL GGCLOS

        ENDIF

	ELSEIF (CHAR(TState).EQ.'G') THEN
	  IF (CHAR(ChName).EQ.'C') THEN
C-->   This peace of code is for Chubarian's chamber.
C  Offset for chubarian's chamber
             DZOFF = 51.435

             TRPAR(1) = 9.2
             TRPAR(2) = 24.1
             TRPAR(3) = 14.6
             TRPAR(4) = 26.06
	     CALL GSVOLU('MCHA','TRD1',23,TRPAR,4,IVOLU)
             PAR(1) = 0.
             PAR(2) = 5.08
             PAR(3) = 3.175
             CALL GSVOLU('TCHA','TUBE',23,PAR,3,IVOLU)
             TRPAR(1) = 8.7
             TRPAR(2) = 23.6
             TRPAR(3) = 14.1
             TRPAR(4) = 25.56
	     CALL GSVOLU('GMVO','TRD1',27,TRPAR,4,IVOLU)
             PAR(1) = 0.
             PAR(2) = 4.75
             PAR(3) = 3.175
             CALL GSVOLU('GTVO','TUBE',27,PAR,3,IVOLU)
             PAR(1) = 0.
             PAR(2) = 1.75
             PAR(3) = 0.25
             CALL GSVOLU('HOLE','TUBE',27,PAR,3,IVOLU)
             PAR(1) = 0.
             PAR(2) = 1.75
             PAR(3) = 0.0001
             CALL GSVOLU('WIND','TUBE',29,PAR,3,IVOLU)
             PAR(1) = 0.
             PAR(2) = 5.00
             PAR(3) = 20.
             CALL GSVOLU('VACU','TUBE',25,PAR,3,IVOLU)

C-->   dE-E telescope (diam 3 cm); dE 100 mkm, E 3 mm
             PAR(1)=0.0
             PAR(2)=1.5
             PAR(3)=0.005
             CALL GSVOLU('dEdt','TUBE',22,PAR,3,IVOLU)
    
             PAR(1)=0.0
             PAR(2)=1.5
             PAR(3)=0.15
             CALL GSVOLU('Edt ','TUBE',22,PAR,3,IVOLU)

C-->   Colimator with diameter RDET.
             PAR(1)=RDET
             PAR(2)=2.0
             PAR(3)=0.2
             CALL GSVOLU('CALM','TUBE',23,PAR,3,IVOLU)

             CALL GSPOS('MCHA',1,'MVOL',0.,0.,0.,0,'MANY')
             CALL GSPOS('TCHA',1,'MVOL',0.,0.,-29.235,0,'MANY')
             CALL GSPOS('GMVO',1,'MCHA',0.,0.,0.,0,'MANY')
             CALL GSPOS('GTVO',1,'MVOL',0.,0.,-28.735,0,'ONLY')
             CALL GSPOS('HOLE',1,'TCHA',0.,0.,-2.925,0,'MANY')
             CALL GSPOS('WIND',1,'HOLE',0.,0.,-0.24,0,'ONLY')
             CALL GSPOS('VACU',1,'MVOL',0.,0.,-52.41,0,'MANY')

             zdet = Apde*cos(D2R*Aang) - DZOFF
             xdet = Apde*sin(D2R*Aang)

             zdet1 = (Apde+0.5)*cos(D2R*Aang) - DZOFF
             xdet1 = (Apde+0.5)*sin(D2R*Aang)
      
             zdet2 = (Apde+1.0)*cos(D2R*Aang) - DZOFF
             xdet2 = (Apde+1.0)*sin(D2R*Aang)

             CALL GSROTM (1,(90.+ Aang),0.,90.,90.,Aang,0.)
             CALL GSPOS('CALM',1,'GMVO',xdet,0.,zdet,1,'MANY')
             CALL GSPOS('dEdt',1,'GMVO',xdet1,0.,zdet1,1,'MANY')
             CALL GSPOS('Edt ',1,'GMVO',xdet2,0.,zdet2,1,'MANY')
C-->   Definition of the hit
           
	     CHNMSV(1) = 'MVOL'
	     CHNMSV(2) = 'MCHA'
	     CHNMSV(3) = 'GMVO'
	     CHNMSV(4) = 'dEdt'

             CALL GSDET ('TLSC','dEdt',1,CHNMSV,NBITSV,1,1000,
     &                1000,iset,idet)
             CHNAMH(1) = 'dE  '
             CALL GSDETH ('TLSC','dEdt',1,CHNAMH,NBITSH,ORIG,FACT)

             CHNMSV(4) = 'Edt '
             CALL GSDET ('TLSC','Edt ',1,CHNMSV,NBITSV,1,1000,
     &                 1000,iset,idet)
             CHNAMH(1) = 'E   '
             CALL GSDETH ('TLSC','Edt ',1,CHNAMH,NBITSH,ORIG,FACT) 

           ELSEIF (CHAR(ChName).EQ.'K') THEN
C -->	This peace of code is for Kirby Kempers chamber
             DZOFF = 51.435

             PAR(1) = 0.0
             PAR(2) = 25.0 
             PAR(3) = 30.0
	     CALL GSVOLU('MCHA','TUBE',23,PAR,4,IVOLU)
             PAR(1) = 0.
             PAR(2) = 5.08
             PAR(3) = 3.175
             CALL GSVOLU('TCHA','TUBE',23,PAR,3,IVOLU)
             PAR(1) = 0.0 
             PAR(2) = 24.5 
             PAR(3) = 29.0
	     CALL GSVOLU('GMVO','TUBE',27,PAR,4,IVOLU)
             PAR(1) = 0.
             PAR(2) = 4.75
             PAR(3) = 3.175
             CALL GSVOLU('GTVO','TUBE',27,PAR,3,IVOLU)
             PAR(1) = 0.
             PAR(2) = 0.35
             PAR(3) = 0.5 
             CALL GSVOLU('HOLE','TUBE',27,PAR,3,IVOLU)
             PAR(1) = 0.
             PAR(2) = 1.75
             PAR(3) = 0.0001
             CALL GSVOLU('WIND','TUBE',29,PAR,3,IVOLU)
             PAR(1) = 0.
             PAR(2) = 5.00
             PAR(3) = 20.
             CALL GSVOLU('VACU','TUBE',25,PAR,3,IVOLU)

C-->   dE-E telescope (diam 3 cm); dE 100 mkm, E 3 mm
             PAR(1)=0.0
             PAR(2)=1.5
             PAR(3)=0.005
             CALL GSVOLU('dEdt','TUBE',22,PAR,3,IVOLU)
    
             PAR(1)=0.0
             PAR(2)=1.5
             PAR(3)=0.15
             CALL GSVOLU('Edt ','TUBE',22,PAR,3,IVOLU)

C-->   Colimator with diameter RDET.
             PAR(1)=RDET
             PAR(2)=2.0
             PAR(3)=0.2
             CALL GSVOLU('CALM','TUBE',23,PAR,3,IVOLU)

             CALL GSPOS('MCHA',1,'MVOL',0.,0.,0.,0,'MANY')
             CALL GSPOS('TCHA',1,'MVOL',0.,0.,-29.235,0,'MANY')
             CALL GSPOS('GMVO',1,'MCHA',0.,0.,0.,0,'MANY')
             CALL GSPOS('GTVO',1,'MVOL',0.,0.,-28.735,0,'ONLY')
             CALL GSPOS('HOLE',1,'TCHA',0.,0.,-2.925,0,'MANY')
             CALL GSPOS('WIND',1,'HOLE',0.,0.,-0.24,0,'ONLY')
             CALL GSPOS('VACU',1,'MVOL',0.,0.,-52.41,0,'MANY')

             zdet = Apde*cos(D2R*Aang) - DZOFF
             xdet = Apde*sin(D2R*Aang)

             zdet1 = (Apde+0.5)*cos(D2R*Aang) - DZOFF
             xdet1 = (Apde+0.5)*sin(D2R*Aang)
      
             zdet2 = (Apde+1.0)*cos(D2R*Aang) - DZOFF
             xdet2 = (Apde+1.0)*sin(D2R*Aang)

             CALL GSROTM (1,(90.+ Aang),0.,90.,90.,Aang,0.)
             CALL GSPOS('CALM',1,'GMVO',xdet,0.,zdet,1,'MANY')
             CALL GSPOS('dEdt',1,'GMVO',xdet1,0.,zdet1,1,'MANY')
             CALL GSPOS('Edt ',1,'GMVO',xdet2,0.,zdet2,1,'MANY')
C-->   Definition of the hit
           
	     CHNMSV(1) = 'MVOL'
	     CHNMSV(2) = 'MCHA'
	     CHNMSV(3) = 'GMVO'
	     CHNMSV(4) = 'dEdt'

             CALL GSDET ('TLSC','dEdt',1,CHNMSV,NBITSV,1,1000,
     &                1000,iset,idet)
             CHNAMH(1) = 'dE  '
             CALL GSDETH ('TLSC','dEdt',1,CHNAMH,NBITSH,ORIG,FACT)

             CHNMSV(4) = 'Edt '
             CALL GSDET ('TLSC','Edt ',1,CHNMSV,NBITSV,1,1000,
     &                 1000,iset,idet)
             CHNAMH(1) = 'E   '
             CALL GSDETH ('TLSC','Edt ',1,CHNAMH,NBITSH,ORIG,FACT) 


		print *,'OK! I am here'
C		pause
	   ENDIF
        ELSEIF (CHAR(TState).EQ.'L') THEN
C-->	  Liquid target option is not implemented yet.
	     print *,'Geometry is not defined for Liquid target'
	     stop
	     ELSE
      ENDIF

      CALL GGCLOS

C-->  DRAWING of the experimental set-up.
       IF (IDEMIN.GE.1) THEN
         CALL GDRAW ('MVOL',50.,150.,0.,10.,10.,0.16,0.16)

         CALL GDAXIS (0.,0.,0.,25.)
         CALL GDSCAL (2.,1.)
         CALL IUWK (0,1)
       ENDIF

      RETURN

      END
