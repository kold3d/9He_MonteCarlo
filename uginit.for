       SUBROUTINE UGINIT
C-----------------------------------------------------------------
C-     Initialization of GEANT Monte Carlo simulation process
C-     Author: Grigory Rogachev (21-11-01)
C-     Inverse geometry and very thick target experiment.
C-----------------------------------------------------------------

      real tlife,UB(2),BRATIO(6),MASSIN,AMU
      real RENE,RWID
      real enkey,bang(2),bene(2)
      integer MODE(6),Lanfl,Lenfl,R1,R2
      integer FIDX
      character flnm*20
      PARAMETER (AMU=931.494)
      
      INCLUDE 'udata.inc'
      INCLUDE 'gcflag.inc'
      INCLUDE 'gcbank.inc'

C-->     Initialise GEANT
      CALL GINIT

C-->     Read data records      
      CALL FFKEY('RTYP',RType,1,'INTE')
      CALL FFKEY('ANDI',Lanfl,1,'INTE')
      CALL FFKEY('ENDI',enkey,3,'REAL')
      CALL FFKEY('TSTA',TState,1,'INTE')
      CALL FFKEY('CHNM',ChName,1,'INTE')
      CALL FFKEY('GPRE',Gpre,1,'REAL')
      CALL FFKEY('LDEN',LDEN,1,'REAL')
      CALL FFKEY('APDE',Apde,1,'REAL')
      CALL FFKEY('AANG',Aang,1,'REAL')
      CALL FFKEY('BANG',bang,2,'REAL')
      CALL FFKEY('BENE',bene,2,'REAL')
      CALL FFKEY('BRAD',dRad,2,'REAL')
      CALL FFKEY('RDET',RDET,1,'REAL')
      CALL FFKEY('FIDX',FIDX,1,'INTE')
      CALL GFFGO

      dtheta = bang(2)
      theav  = bang(1)
      BeamE  = bene(1)*1.E-3
      dEner  = bene(2)*1.E-3

C-->    Read the angular distribution data file 'ang.dat' if Lanfl = 1
C-->    or making angular distribution uniform if Lanfl = 0.
C-->    Array andis(i,k) prepered here.
      IF (Lanfl.eq.1) THEN 
         OPEN (3,FILE='ang.dat',STATUS='OLD')
         sigmax = 0.
         DO i = 1,1000
	   READ (3,*,END=14) andis(i,1),andis(i,2)
           andis(i,1) = andis(i,1)*0.0174533
           IF (sigmax.LT.andis(i,2)) sigmax=andis(i,2)
         ENDDO
  14     CONTINUE  
      ELSE 
         DO i = 1,181
           andis(i,1) = real(i-1)*0.0174533
           andis(i,2) = 1.
           sigmax = 1.
         ENDDO
      ENDIF
         nang = i - 1

C-->  Read the excitation function data file 'energy.dat' if Lenfl = 1
C-->  or making excitation function gauss shaped with parameters given
C-->  by RENE and RWID (FWHM) if Lenfl = 0
C      Lenfl = int(enkey(1)) 
C      RENE  = enkey(2)*1.E-3
C      RWID  = enkey(3)*1.E-3
C      IF (Lenfl.eq.1) THEN 
C         OPEN (4,FILE='energy.dat',STATUS='OLD')
C         rsmax = 0.
C         DO i = 1,30000
C	   READ (4,*,END=15) endis(i,1),endis(i,2)
C           IF (rsmax.LT.endis(i,2)) rsmax=endis(i,2)
C         ENDDO
C  15     CONTINUE  
C      ELSE 
C         DO i = 1,30000
C           endis(i,1) = real(i)*1.E-6
C           endis(i,2) = exp(-(endis(i,1)-RENE)**2./(2.*(RWID/2.34)**2.))
C           rsmax = 1.
C         ENDDO
C      ENDIF
 
C         nen = i - 1


C-->     Initialise data structure
      CALL GZINIT

C-->     Initialise graphics
      IF (IDEMIN.GE.1) THEN 
        CALL HPLINT(2)
        CALL GDINIT
      ENDIF
     
      CALL GPART
      CALL GPIONS

C-->  Direct corrections to masses of standard GEANT particles
C	JPA = LQ(JPART-46)
C        Q(JPA+7) = 2.80892
C	JPA = LQ(JPART-47)
C        Q(JPA+7) = 3.72738
C        JPA = LQ(JPART-61)
C        Q(JPA+7) = 5.60152
C-->

C-->     Definition of particle unstable He5 and its mode of decay
C-->     UB(1) gives width of the GS in keV; UB(2) is mass in GeV
C      tlife = 1.E-21 
C      UB(1) = 0.6E-3
C      UB(2) = 4.6678396
C      CALL GSPART(123,'He5 ',4,UB(2),2.,tlife,UB,2)
      BRATIO(1) = 100.
C      MODE(1) = 1347
       DO i = 2,6
        BRATIO(i) = 0.
C        MODE(i)   = 0
       ENDDO
C      CALL GSDK (123,BRATIO,MODE)

C--> Definition of 6He in GS
      CALL GSPART(124,'He6 ',8,5.605538,2.,0.8067,0.,0)

C--> Definition of 6He in the first excited state, 2+
C      tlife = 2.E-20
C      UB(1) = 0.113E-3
C      UB(2) = 5.607335
C      CALL GSPART(125,'He6*',4,UB(2),2.,tlife,UB,2)
C      MODE(1) = 131347
C      CALL GSDK (125,BRATIO,MODE)

C--> Definition of 6Li in the first excited state, 3+
C      tlife = 1.0E-19
C      UB(1) = 2.4E-5
C      UB(2) = 5.603706
C      CALL GSPART(126,'Li6*',4,UB(2),3.,tlife,UB,2)
C      MODE(1) = 4547
C      CALL GSDK (126,BRATIO,MODE)

C--> Definition of 8He in GS.
C--> The number for incoming beam particle should be 132.

       MASSIN =	(AMU*8. + 31.609)/1000.
       CALL GSPART(133,'8He ',8,MASSIN,2.,0.119,0.,0)
C      CALL GSPART(132,'14O ',8,13.049007,8.,70.6,0.,0)
C	CALL GSPART(132,'8B',8,7.474411,5.,0.77,0.,0)
C	CALL GSPART(132,'11C',8,10.254084,6.,1223.4,0.,0)
C       MASSIN =	(AMU*7.+15.7689)/1000.
C       CALL GSPART(132,'7Be',8,6.5362,4.,10000.,0.,0)
        MASSIN =	(AMU*18.-0.781)/1000.
        CALL GSPART(132,'18O',8,MASSIN,8.,10000.,0.,0)

C        MASSIN = (AMU*18.+5.317)/1000.
C        CALL GSPART(133,'18Ne',8,MASSIN,10.,1.67,0.,0)

        MASSIN = (AMU*9.+40.80)/1000.
        CALL GSPART(134,'9He',8,MASSIN,2.,1000.,0.,0)


C--> Definition of 7Be in the first excited state, 1/2- at 0.429 MeV
C       tlife = 10000
C       UB(1) = 1.0E-10
C       UB(2) = MASSIN+0.429/1000.
C       CALL GSPART(133,'7Be*',4,UB(2),4.,tlife,UB,2)
C       MODE(1) = 13201
C       BRATIO(1) = 100.  
C       CALL GSDK(133,BRATIO,MODE)
       
      CALL UGEOM
      CALL GPHYSI

C-->  Initialisation of Random Numbers Generator
      CALL RMARIN(NRNDM(1),0,0)        

C-->   Output definition
       IF (RType.GE.10) THEN
        R1 = int(RType/10.)
        R2 = RType - 10
        flnm = 'events'//char(48+R1)//char(48+R2)//'_'//char(48+Lanfl)//
     &         '.hbook'
       ELSE
C        flnm = 'events'//char(48+RType)//'_'//char(48+Lanfl)//'.hbook'
         flnm = 'events'//char(48+FIDX)//'.hbook'
       ENDIF

C      CALL HROPEN(1,'LUN1',flnm,'N',1024,ISTAT)
      CALL HROPEN(1,'LUN1','histo.hbook','N',1024,ISTAT)
      CALL HBNT(1,'Event file','D')
      CALL HBNAME(1,'event',E1,'E1:r,X1:r,Y1:r,RAD1:r,E2:r,X2:r,Y2:r,
     &             RAD2:r,anincm:r')
      CALL HCDIR('//LUN1',' ')
C      CALL HBOOK1(10,'CMS angle of reaction',360,0.,3.14159,0.)
      CALL HBOOK1(10,'Spectrum of pi+',104,-10.,16.,0.)
C      CALL HBOOK2(11,'X-Y target beam hit',1000,1000,-2.,2.,-2.,2.,0.)
      
      RETURN
      END
