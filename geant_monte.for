C    *  PROGRAM FOR GEANT Simulation of Experiment
C    *       with very thick target and inverse kinematics.
C    *
C    *              G.V. Rogachev
C    *
      PROGRAMM Gold
      PARAMETER (NGBANK=50000,NHBOOK=20000)
      COMMON/GCBANK/Q(NGBANK)
      COMMON/PAWC/H(NHBOOK)
C-->     Initialises HBOOK and GEANT memory
      CALL GZEBRA(NGBANK)
      CALL HLIMIT(-NHBOOK)
C-->     GEANT Initialisation
      CALL UGINIT
C-->     Start events processing
      CALL GRUN
C-->     End of run
      CALL UGLAST
      END
c=========================================================================

      SUBROUTINE UGINIT
      COMMON/GCFLAG/IDEBUG,IDEMIN,IDEMAX,ITEST,IDRUN,IDEVT,IEORUN
     +        ,IEOTRI,IEVENT,ISWIT(10),IFINIT(20),NEVENT,NRNDM(2)
      COMMON/GCFLAX/BATCH, NOLOG
      LOGICAL BATCH, NOLOG

       PARAMETER (MAXMEC=30)
      COMMON/GCTRAK/VECT(7),GETOT,GEKIN,VOUT(7),NMEC,LMEC(MAXMEC)
     + ,NAMEC(MAXMEC),NSTEP ,MAXNST,DESTEP,DESTEL,SAFETY,SLENG
     + ,STEP  ,SNEXT ,SFIELD,TOFG  ,GEKRAT,UPWGHT,IGNEXT,INWVOL
     + ,ISTOP ,IGAUTO,IEKBIN, ILOSL, IMULL,INGOTO,NLDOWN,NLEVIN
     + ,NLVSAV,ISTORY
C
C      INCLUDE 'GCTRAK.INC'
      COMMON/GCKINE/IKINE,PKINE(10),ITRA,ISTAK,IVERT,IPART,ITRTYP
     +      ,NAPART(5),AMASS,CHARGE,TLIFE,VERT(3),PVERT(4),IPAOLD
C
C      INCLUDE 'GCKINE.INC'
      COMMON/UEDEP/ENE(15)
      COMMON/UTHE/THET(2)
      COMMON/UPRS/PARM(3)
      CHARACTER*30 FNAME
      COMMON/HFNAM/FNAME
      real*8 seed
      DIMENSION UHLIM(10)
      DATA UHLIM /10*0.05/
      DATA PARM /0.03,1000.,0.001/
      DATA THET /0.,0.5/
C-->     Initialise GEANT
      CALL GINIT
      IGAUTO=0
      seed= 445273452.9
      call rdmin(seed)
      PL2T=0.5
      FNAME='moi.rzd'
      CALL FFKEY('CONC',IGAUTO,1,'INTEGER')
      CALL FFKEY('ULIM',UHLIM,10,'REAL')
      CALL FFKEY('THET',THET,2,'REAL')
      CALL FFKEY('PARM',PARM,3,'REAL')
      CALL FFKEY('UPL2',PL2T,1,'REAL')
      CALL FFKEY('FNAM',FNAME,30,'MIXED')
      CALL GFFGO
      CALL GZINIT
      CALL GPART
      CALL GMATE
      CALL UGEOM
      CALL GPHYSI
      CALL GPRINT('MATE',0)
      CALL GPRINT('TMED',0)
      CALL GPRINT('VOLU',0)
C         HBOOK (ID,CHTITL,NX,XMI,XMA,ZMX)
C         Action: Books a one-dimentional histogram.
C Input parameters:
c ID    - histogram identifier, integer, non zero;
c CHTITL- histogram title (charact up to 80);
c NX    - number of channels;
c XMI   - lower edge og first channel;
c XMA   - upper edge og channel;
c VMX   - upper limit of single channel content.
c         (VMX=0 means 1 word per channel).
      CALL HBOOK1(1,'1. CENTR SINGL   ',100,0.0, UHLIM(1), 0)
      CALL HBOOK1(2,'2. ESOURCE       ',100,0.0, UHLIM(1), 0)
      CALL HBOOK1(3,'3. CENTR 1MEV-THR',100,0.0, UHLIM(1), 0)
      CALL HBOOK1(4,'4. CENTR 5MEV-THR ',100,0.0, UHLIM(1), 0)
      CALL HBOOK1(11,'11 SUM 100KEV-THR',100,0.0, UHLIM(1), 0)
      CALL HBOOK1(12,'12 SUM 200KEV-THR',100,0.0, UHLIM(1), 0)
      CALL HBOOK1(13,'13 SUM 1MEV-THR  ',100,0.0, UHLIM(1), 0)
      CALL HBOOK1(14,'14 SUM 5MEV-THR  ',100,0.0, UHLIM(1), 0)
      CALL HBOOK1(15,'15 SUM SINGLE    ',100,0.0, UHLIM(1), 0)
C
C      CALL HBOOK1(20,'10.CsJ(Tl) Perif',100,0.0, UHLIM(1), 0)
C      CALL HBOOK1(21,'11.BGO+ 7CsJ(Tl)',100,0.0, UHLIM(1), 0)
C      CALL HBOOK1(22,'12.CSJ+BGO IN TH',100,0.0, UHLIM(1), 0)
c
c    HBOOK2 (ID,CHTITL,NX,XMI,XMA,NY,YMI,YMA,VMX)
C    Action: Book a two-dimentional histogram.
C Input parameters:
c ID    - histogram identifier, integer, non zero;
c CHTITL- histogram title (charact up to 80);
c NX    - number of X channels;
c XMI   - lower edge og first X channel;
c XMA   - upper edge og X channel;
c NY    - number of Y channels;
c YMI   - lower edge og first Y channel;
c YMA   - upper edge og Y channel;
c VMX   - upper limit of single channel content.
c         (VMX=0 means 1 word per channel).
C      CALL HBOOK2(201,'E  BGO per Etotal  ',100,0.,50.,100,0.,50.,0)
      DO I=1,15
      ENE(I)=0
      ENDDO
C    TEST SOURS
      IF(ISWIT(5).EQ.1) THEN
      WRITE(*,*) 'EMISSION  ALONG Z AXIS   '
      ENDIF
      IF(ISWIT(5).EQ.2) THEN
      WRITE(*,*) 'EMISSION IZOTROPIC IN 4Pi '
      ENDIF
      IF(ISWIT(5).EQ.3) THEN
      WRITE(*,*) 'EMISSION IZOTROPIC IN CONUS   '
      WRITE(*,*) '  THETA_MIN=',THET(1),'  THETA_MAX=',THET(2)
      ENDIF
      IF(ISWIT(6).EQ.1) THEN
      WRITE(*,*) '  ENERGY IS DELTA FUNCTION  '
      WRITE(*,*) ' ENERGY =',PKINE(1),' [GeV]'
      ENDIF
      IF (ISWIT(6).EQ.2) THEN
      WRITE(*,*) 'ENERGY IS IZOTROPIC IN RANGE [E_min,E_max] '
      WRITE(*,*) '  Emin=',PKINE(1),'  Emax=',PARM(2)
      ENDIF
      END
C
c========================================================================
C
      SUBROUTINE UGEOM
      COMMON/UPL2T/PL2T
      REAL*4 FIELDM,TMAXFD,DMAXMS,DEEMAX,EPSIL,STMIN
      DIMENSION Ach2(2),Zch2(2),Wch2(2)
      DIMENSION Abgo(3),Zbgo(3),Wbgo(3)
      DIMENSION Acsj(2),Zcsj(2),Wcsj(2)
      REAL*4 PAR(10)
      CHARACTER*20 NAME_ZONE(5)
C   * Definition of the detector  materials
C   * A A mass numbers
C   * Z Z charge numbers
C   * W Weingt numbers on the molecula
c
c       plastic compound parameters
c
      DATA Ach2/12.01,1.01/
      DATA Zch2/6.,1./
      DATA Wch2/1.,1./
c
c       BGO compound parameters
c
      DATA Abgo/208.98,72.59,15.999/
      DATA Zbgo/83.,32.,8./
      DATA Wbgo/4.,3.,12./
c
c       CsJ(tl) compound parameters
c
      DATA Acsj/132.91,126.9/
      DATA Zcsj/55.,53./
      DATA Wcsj/1.,1./

      NAME_ZONE(1)='AIR$     '
      NAME_ZONE(2)='CH2$     '
      NAME_ZONE(3)='CSJ$     '
      NAME_ZONE(4)='BGO$     '
      NAME_ZONE(5)='LEAD$    '
      CALL GSMATE(1,NAME_ZONE(1),14.61,7.3,0.001205,30423.,6750.,0,0)
      CALL GSMIXT(2,NAME_ZONE(2),ACH2,ZCH2,1.032,-2,WCH2)
      CALL GSMIXT(3,NAME_ZONE(3),ACSJ,ZCSJ,4.530,-2,WCSJ)
      CALL GSMIXT(4,NAME_ZONE(4),ABGO,ZBGO,7.130,-3,WBGO)
      CALL GSMATE(5,NAME_ZONE(5),207.19,82.,11.35 ,0.56,18.5,0,0)
C
C     Defines trecking media parameterd
C
      ISVOL=1
      IFIELD=0
      FIELDM=0.
      TMAXFD=20.
      STEMAX=1000.
      DEEMAX=0.025
      EPSIL=0.001
      STMIN=0.001
      DO IM=1,5
      CALL GSTMED(IM,NAME_ZONE(IM),IM,ISVOL,IFIELD,FIELDM,TMAXFD,
     + STEMAX,DEEMAX,EPSIL,STMIN,0,0)
      ENDDO
C      WRITE(*,'('' TMAXFD='',G10.4,'' STEMAX='',G12.4,'' DEEMAX='',G12.4,
C     +          '' STMIN='',G12.4)')TMAXFD,STEMAX,DEEMAX,STMIN
C
C           Defines volumes
C
c
c          Definition of the Mather volume "ZONE"
c
      PAR(1)=200.
      PAR(2)=200.
      PAR(3)=200.
      CALL GSVOLU('ZONE','BOX ',1,PAR,3,IVOLU)
c
c          Definition of the BGO detector
c
C      PAR(1)=0.
C      PAR(2)=2.49
C      PAR(3)=1.24
C      CALL GSVOLU('BGO1','TUBE',4,PAR,3,IVOLU)
C
C           Definition og the CsJ(Tl) moduls
C
      PAR(1)=0.
      PAR(2)=360.
      PAR(3)=6.
      PAR(4)=2.
      PAR(5)=-7.
      PAR(6)=0.0
      PAR(7)=4.33
      PAR(8)=7.
      PAR(9)=0.0
      PAR(10)=4.33
      CALL GSVOLU('CS_1','PGON',3,PAR,10,IVOLU)
      CALL GSVOLU('CS_2','PGON',3,PAR,10,IVOLU)
      CALL GSVOLU('CS_3','PGON',3,PAR,10,IVOLU)
      CALL GSVOLU('CS_4','PGON',3,PAR,10,IVOLU)
      CALL GSVOLU('CS_5','PGON',3,PAR,10,IVOLU)
      CALL GSVOLU('CS_6','PGON',3,PAR,10,IVOLU)
      CALL GSVOLU('CS_7','PGON',3,PAR,10,IVOLU)
c
C                      kollimator Pb_1
      PAR(1)=5.0
      PAR(2)=15.0
      PAR(3)=2.5
      CALL GSVOLU('PB_1','TUBE',5,PAR,3,IVOLU)
c
c                      plastic  CH_1
      PAR(1)=5.0
      PAR(2)=5.0
      PAR(3)=5.0
      CALL GSVOLU('CH_1','BOX ',2,PAR,3,IVOLU)

C
      CALL GSPOS( 'CS_1',1,'ZONE', 0.    , 8.6603 ,57.1  ,0,'ONLY')
      CALL GSPOS( 'CS_2',1,'ZONE', 7.2501, 4.3302 ,57.1  ,0,'ONLY')
      CALL GSPOS( 'CS_3',1,'ZONE', 7.2501,-4.3302 ,57.1  ,0,'ONLY')
      CALL GSPOS( 'CS_4',1,'ZONE', 0.    ,-8.6603 ,57.1  ,0,'ONLY')
      CALL GSPOS( 'CS_5',1,'ZONE',-7.2501,-4.3302 ,57.1  ,0,'ONLY')
      CALL GSPOS( 'CS_6',1,'ZONE',-7.2501, 4.3302 ,57.1  ,0,'ONLY')
      CALL GSPOS( 'CS_7',1,'ZONE', 0.    , 0.     ,57.1  ,0,'ONLY')
      CALL GSPOS( 'PB_1',1,'ZONE', 0.    , 0.     ,45.0  ,0,'ONLY')
      CALL GSPOS( 'CH_1',1,'ZONE', 0.    , 0.     ,36.0  ,0,'ONLY')
      CALL GGCLOS
      END


      SUBROUTINE GUKINE
      COMMON/GCFLAG/IDEBUG,IDEMIN,IDEMAX,ITEST,IDRUN,IDEVT,IEORUN
     +        ,IEOTRI,IEVENT,ISWIT(10),IFINIT(20),NEVENT,NRNDM(2)
      COMMON/GCFLAX/BATCH, NOLOG
      LOGICAL BATCH, NOLOG
C
      COMMON/GCKINE/IKINE,PKINE(10),ITRA,ISTAK,IVERT,IPART,ITRTYP
     +      ,NAPART(5),AMASS,CHARGE,TLIFE,VERT(3),PVERT(4),IPAOLD
C
C      INCLUDE 'GCKINE.INC'
C      INCLUDE 'GCFLAG.INC'
      COMMON/UTHE/THET(2)
      COMMON/UPRS/PARM(3)
      COMMON/UMOI/EMOI,DEMOI,EGSS
      DIMENSION VERTEX(3),PLAB(3)
      CHARACTER*20 NP
      SAVE VERTEX,PLAB
c      REAL*4 UBUF(1)
      REAL EMOI, DEMOI, EGSS
      DATA VERTEX/3*0./
      DATA PLAB  /3*0./
c asaasaasa
c      WRITE(*,*) 'IKINE=', IKINE
c      WRITE(*,*) ' NP=', NP
c      WRITE(*,*) ' ITR=',ITR
c      WRITE(*,*) ' AMASS=', AMASS
c      WRITE(*,*) ' CHARGE=',CHARGE
c      WRITE(*,*) ' TLIFE=', TLIFE
c asaasaasa  
c      CALL GFPART(IKINE,NP,ITR,AMASS,CHARGE,TLIFE,0,0)
c      DEMOI = 5.0
C
C   DEFINITION GEOMETRY OF SOURS PARTICLE EMISSION
C
C                              MOIK(1).EQ.1 = ALONG Z AXIS
      IF(ISWIT(5).EQ.1) THEN
      FI=0.0
      TH=0.0
      ENDIF
C                              MOIK(1).EQ.2 = IZOTROPIC IN 4Pi
      IF(ISWIT(5).EQ.2) THEN
      FI=6.283185*RNDM(1.)
      TH=1.57075*RNDM(1.)
      ENDIF
C                            MOIK(1).EQ.3 = IZOTROPIC IN CONUS
C                      THET(1)=MIN,  THET(2) = MAX (IN RADIAN)
      IF(ISWIT(5).EQ.3) THEN
      FI=6.283185*RNDM(1.)
      TH=RNDM(1.)*(THET(1)-THET(2))+THET(2)
c      WRITE(*,*) 'th=',TH
      ENDIF
C
C              DEFINITION ENERGY OF SOURS PARTICLE EMISSION
C
C                            MOIK(2).EQ.1 = DELTA FUNCTION
C     AMASS=MASSA OF PARTICLE [GeV]   PKINE(1)=E_KIN [GeV]
      IF(ISWIT(6).EQ.1) THEN
      P=(PKINE(1)**2 - AMASS**2)**0.5
      EMOI= PKINE(1)*1000.0
      ENDIF
C
C                 MOIK(2).EQ.2 = IZOTROPIC IN RANGE [E_min,E_max]
C
C      AMASS=MASSA OF PARTICLE [GeV]
C      PKINE(1)=E_KIN_MIN [GeV]
C      PARM(1)=E_KIN_MAX [GeV]
C      EMOI=RAND*(MAX-MIN)+MIN
C
      IF (ISWIT(6).EQ.2) THEN
      EMOI=RNDM(1.)*(PARM(2)-PKINE(1))+PKINE(1)
      P=(EMOI**2 - AMASS**2)**0.5
      ENDIF
C
C                 MOIK(2).EQ.3 = AS A FUNCTION C*E* exp(-E/T)
C
C      AMASS=MASSA OF PARTICLE [GeV]
C      PKINE(1)=E_KIN_MIN [GeV]
C      PARM(1)=E_KIN_MAX [GeV]
C      PARM(2)=C (constanta)
C      PARM(3)=T (Temperature) [GeV]
C      EMOI=RAND*
C
      IF (ISWIT(6).EQ.3) THEN
      W    =10000.
      WF   =1000.
      Emin =1.0
      Emax =PKINE(1)*1000.
      C    =50000.
      C1    =1000.
      C2    =4000.
      T1    =0.8
      T2    =1.0
      Eg1   =9.0
      Eg2   =12.0
      G     =1.3
      UPPER=10000.
      DO WHILE (W.GE.WF)
      W  =RNDM(1.)*UPPER
      EMOI=RNDM(1.)*(Emax-Emin)+Emin
      e=EMOI
      F1=50000*e*e/T1/T1*exp(-e/T1)
      L1=C1*e*G/((e*e-Eg1*Eg1)**2 +G*G*Eg2*Eg2)
      L2=C2*e*G/((e*e-Eg2*Eg2)**2 +G*G*Eg2*Eg2)
      F2=2000*e*e/T2/T2*exp(-e/T2)*(L1+L2)
      WF=F1+F2
      ENDDO
      P=((EMOI*0.001)**2 - AMASS**2)**0.5
      ENDIF
      PMOI   =P*1000.
      PLAB(1)=SIN(TH)*COS(FI)*P
      PLAB(2)=SIN(TH)*SIN(FI)*P
      PLAB(3)=COS(TH)*P
C ASA
      VERTEX(1)=0.+RNDM(1.)
      VERTEX(2)=0.+RNDM(1.)
      VERTEX(3)=0.
      CALL GSVERT(VERTEX,0,0,0,0,NVERT)
      CALL GSKINE(PLAB,IKINE,NVERT,0,0,NT)
      IF (IDEBUG.EQ.1.AND.ISWIT(1).EQ.1) THEN
      CALL GPRINT('VERT',0)
      CALL GPRINT('PKINE',0)
      ENDIF
c       WRITE(*,*) 'EMIN=',EMIN,'  EMAX=',EMAX,'  EMOI=',EMOI
      EGSS =GSS(EMOI,DEMOI)
      END

      SUBROUTINE GUTREV
      COMMON/GCFLAG/IDEBUG,IDEMIN,IDEMAX,ITEST,IDRUN,IDEVT,IEORUN
     +        ,IEOTRI,IEVENT,ISWIT(10),IFINIT(20),NEVENT,NRNDM(2)
      COMMON/GCFLAX/BATCH, NOLOG
      LOGICAL BATCH, NOLOG
C
C      INCLUDE 'GCFLAG.INC'
      CALL GTREVE
      IF(IDEBUG.EQ.1.AND.ISWIT(2).EQ.1) THEN
      WRITE(*,*)'IN GUTREV'
      ENDIF
      END

      SUBROUTINE GUSTEP
      COMMON/GCFLAG/IDEBUG,IDEMIN,IDEMAX,ITEST,IDRUN,IDEVT,IEORUN
     +        ,IEOTRI,IEVENT,ISWIT(10),IFINIT(20),NEVENT,NRNDM(2)
      COMMON/GCFLAX/BATCH, NOLOG
C
      PARAMETER (MAXMEC=30)
      COMMON/GCTRAK/VECT(7),GETOT,GEKIN,VOUT(7),NMEC,LMEC(MAXMEC)
     + ,NAMEC(MAXMEC),NSTEP ,MAXNST,DESTEP,DESTEL,SAFETY,SLENG
     + ,STEP  ,SNEXT ,SFIELD,TOFG  ,GEKRAT,UPWGHT,IGNEXT,INWVOL
     + ,ISTOP ,IGAUTO,IEKBIN, ILOSL, IMULL,INGOTO,NLDOWN,NLEVIN
     + ,NLVSAV,ISTORY
C
      COMMON/GCTMED/NUMED,NATMED(5),ISVOL,IFIELD,FIELDM,TMAXFD,STEMAX
     +      ,DEEMAX,EPSIL,STMIN,CFIELD,PREC,IUPD,ISTPAR,NUMOLD
C
      COMMON/GCKINE/IKINE,PKINE(10),ITRA,ISTAK,IVERT,IPART,ITRTYP
     +      ,NAPART(5),AMASS,CHARGE,TLIFE,VERT(3),PVERT(4),IPAOLD
C
      INTEGER MXGKIN
      PARAMETER (MXGKIN=100)
      COMMON/GCKING/KCASE,NGKINE,GKIN(5,MXGKIN),
     +                           TOFD(MXGKIN),IFLGK(MXGKIN)
C
      COMMON/GCVOLU/NLEVEL,NAMES(15),NUMBER(15),
     +LVOLUM(15),LINDEX(15),INFROM,NLEVMX,NLDEV(15),LINMX(15),
     +GTRAN(3,15),GRMAT(10,15),GONLY(15),GLX(3)
C
C      INCLUDE 'GCTRAK.INC'
C      INCLUDE 'GCTMED.INC'
C      INCLUDE 'GCKINE.INC'
C      INCLUDE 'GCFLAG.INC'
C      INCLUDE 'GCKING.INC'
C      INCLUDE 'GCVOLU.INC'
      COMMON/UEDEP/ENE(15)
      CHARACTER CHNAM(15)*4
      INTEGER       KCASE,NGKINE ,IFLGK
      REAL          GKIN,TOFD
      LOGICAL BATCH, NOLOG
      EQUIVALENCE (CHNAM,NAMES)
      IF(IDEBUG.EQ.1.AND.ISWIT(3).EQ.1) THEN
      WRITE(*,*) 'I AM IN GUSTEP ',NUMED,DESTEP
      WRITE(*,*)NLEVEL,NAMES(NLEVEL),' ',CHNAM(NLEVEL),' ',
     + NUMBER(NLEVEL),' ',LVOLUM(NLEVEL),' ',LINDEX(NLEVEL),' ',
     + INFROM,' ',NLEVMX,' ',NLDEV(NLEVEL),' ',LINMX(NLEVEL),
     + GTRAN(3,NLEVEL),' ',GRMAT(10,NLEVEL),' ',GONLY(NLEVEL)
      ENDIF
      IF(DESTEP.GT.0.)THEN
      ENE(LVOLUM(NLEVEL))=ENE(LVOLUM(NLEVEL))+1000*DESTEP
      ENDIF
      IF(NGKINE.GT.0) THEN
      DO I=1,NGKINE
      CALL GSKING(I)
      ENDDO
      ENDIF
      END
C
C
      SUBROUTINE GUOUT
      COMMON/GCFLAG/IDEBUG,IDEMIN,IDEMAX,ITEST,IDRUN,IDEVT,IEORUN
     +        ,IEOTRI,IEVENT,ISWIT(10),IFINIT(20),NEVENT,NRNDM(2)
      COMMON/GCFLAX/BATCH, NOLOG
      LOGICAL BATCH, NOLOG
C
C      INCLUDE 'GCFLAG.INC'
      COMMON/UEDEP/ENE(15)
      COMMON/UMOI/EMOI,DEMOI,EGSS
c      WRITE(*,*)'event :',IEVENT,(ENE(J),J=1,10)
       REAL EIII,E_FWHM
      EIII    =0.
      E_FWHM  =0.
      ERES5MEV=0.
      ERES1MEV=0.
      ERES100 =0.
      ERES200 =0.
      ETOT    =0.
      EC5MEV  =0.
      EC1MEV  =0.
      ERING   =0.
      ECSJ    =0.
      ENELEAD =ENE(9)
      EPLAST  =ENE(10)
      EIII    =ENE(8)
       IF(EIII.LE.0.01) THEN
           E_FWHM=0.01
       ELSE
           E_FWHM = EIII*0.03/(((EIII+0.001)*0.001)**(0.25))
       ENDIF
       ECENTR = GSS(EIII,E_FWHM)
C
      DO I=2,8
       IF(ENE(I).LE.0.01) THEN
           E_FWHM=0.01
       ELSE
           E_FWHM = ENE(I)*0.03/(((ENE(I)+0.001)*0.001)**(0.25))
       ENDIF
       ETOT   = ETOT+GSS(ENE(I),E_FWHM)
C       ETOT    =ETOT+ENE(I)
      ENDDO
C
      ERING   = ETOT-ECENTR
      IF(ERING.LE.5.0) THEN
       ERES5MEV  = ETOT
       EC5MEV    = ECENTR
      ENDIF
C
      IF(ERING.LE.1.0) THEN
       ERES1MEV  = ETOT
       EC1MEV    = ECENTR
      ENDIF
C

      IF(ERING.LE.0.10) THEN
      ERES100  = ETOT
      ENDIF
C
C
      IF(ERING.LE.0.20) THEN
      ERES200  = ETOT
      ENDIF
C
c
c          HFILL (ID,X,Y,WEIGHT)
c          Action: Fills a 1 or 2 -dimentional histograms.
c   ID    -  histogram identifier
c   X     -  value of the abscissa
c   Y     -  value of the ordinate
c   WEIGHT-  event weight (positive or negative)
      CALL HFILL(1,ECENTR ,0.,1.)
      CALL HFILL(2,EMOI   ,0.,1.)
      CALL HFILL(3,EC1MEV ,0.,1.)
      CALL HFILL(4,EC5MEV ,0.,1.)
C
      CALL HFILL(11,ERES100  ,0.,1.)
      CALL HFILL(12,ERES200  ,0.,1.)
      CALL HFILL(13,ERES1MEV ,0.,1.)
      CALL HFILL(14,ERES5MEV ,0.,1.)
      CALL HFILL(15,ETOT     ,0.,1.)
C
c      CALL HFILL(17,EPLAST  ,0.,1.)
c      CALL HFILL(18,EBGO    ,0.,1.)
c      CALL HFILL(19,ECENTR  ,0.,1.)
c      CALL HFILL(20,ERING   ,0.,1.)
c      CALL HFILL(21,ETOT    ,0.,1.)
c      CALL HFILL(22,EOUTIC  ,0.,1.)
c      CALL HFILL(201,ETOT   ,EBGO,1.)
      DO I=1,15
      ENE(I)=0
      ENDDO
      END
C
      SUBROUTINE UGLAST
      COMMON/GCKINE/IKINE,PKINE(10),ITRA,ISTAK,IVERT,IPART,ITRTYP
     +      ,NAPART(5),AMASS,CHARGE,TLIFE,VERT(3),PVERT(4),IPAOLD
C
      COMMON/GCFLAG/IDEBUG,IDEMIN,IDEMAX,ITEST,IDRUN,IDEVT,IEORUN
     +        ,IEOTRI,IEVENT,ISWIT(10),IFINIT(20),NEVENT,NRNDM(2)
      COMMON/GCFLAX/BATCH, NOLOG
      LOGICAL BATCH, NOLOG
C
C      INCLUDE 'GCKINE.INC'
C      INCLUDE 'GCFLAG.INC'
      COMMON/UKEY/KEYS(2)
      COMMON/UTHE/THET(2)
      COMMON/UPRS/PARM(3)
      CHARACTER*30 FNAME
      COMMON/HFNAM/FNAME
      CALL GLAST
      CALL HRPUT(0,FNAME,'N')
      IF(ISWIT(4).EQ.1)THEN
      CALL HIDOPT(0,'BLAC')
      CALL HISTDO
      ENDIF
      END
C
C======================================================================
C
      real function rndm(dum)
c----------------------------------------------------------------------
c     real function ranf()
c     uniform random number generator from cern library
c-----------------------------------------------------------------------
      double precision    dranf,    g900gt,   g900st
      double precision    ds(2),    dm(2),    dseed
      double precision    dx24,     dx48
      double precision    dl,       dc,       du,       dr
      logical             single
      data      ds     /  1665 1885.d0, 286 8876.d0  /
      data      dm     /  1518 4245.d0, 265 1554.d0  /
      data      dx24   /  1677 7216.d0  /
      data      dx48   /  281 4749 7671 0656.d0  /
      single  =  .true.
      goto 10
      entry dranf()
      single  =  .false.
  10  dl  =  ds(1) * dm(1)
      dc  =  dint(dl/dx24)
      dl  =  dl - dc*dx24
      du  =  ds(1)*dm(2) + ds(2)*dm(1) + dc
      ds(2)  =  du - dint(du/dx24)*dx24
      ds(1)  =  dl
      dr     =  (ds(2)*dx24 + ds(1)) / dx48
      if(single)  then
         rndm  =  sngl(dr)
      else
         dranf  =  dr
      endif
      return
      entry g900gt()
      g900gt  =  ds(2)*dx24 + ds(1)
      return
      entry g900st(dseed)
      ds(2)  =  dint(dseed/dx24)
      ds(1)  =  dseed - ds(2)*dx24
      g900st =  ds(1)
      return
      end
c-----------------------------------------------------------------------
      subroutine rdmout(seed)
c     subroutine ranfgt(seed)
      double precision    seed,     g900gt,   g900st,   dummy
      seed  =  g900gt()
      return
      entry rdmin(seed)
      dummy  =  g900st(seed)
      return
      end
c======================================================================
      FUNCTION GASDEV(DUMMY)
C==================================================================
C	Returns a normally distributed deviate with zero mean and unit
C     variance using RNDM(DUMMY) as the source of uniform deviates .
C------------------------------------------------------------------
C
      DATA ISET/0/
      IF (ISET.EQ.0) THEN
1     V1=2.*RNDM(DUMMY)-1.
      V2=2.*RNDM(DUMMY)-1.
      R=V1*V1+V2*V2
      IF(R.GE.1.)GOTO 1
      FAC=SQRT(-2.*LOG(R)/R)
      GSET=V1*FAC
      GASDEV=V2*FAC
      ISET=1
      ELSE
      GASDEV=GSET
      ISET=0
      ENDIF
      RETURN
      END
C==================================================================
******************************************************************************
* REAL FUNCTION GSS(XM,DX)
* Version 1.0
*
* Returns a random number distributed according to a Gaussian around XM with
* width DX.
******************************************************************************
	REAL FUNCTION GSS(XM,DX)
	REAL DUM,XM,DX

	GSS=XM+DX*RG32(DUM)
	RETURN

	END
**
