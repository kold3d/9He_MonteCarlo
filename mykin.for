      SUBROUTINE K2(DN,ENERGY,Q,ENEXCI,AP,EOUT1,EOUT2,CMANG1,CMANG2)
    
C=====================================================
C
C	Simple program for non-relativistic
C             two-particle kinematic.
C
C=====================================================

      IMPLICIT NONE

      INTEGER I
      DOUBLE PRECISION DN(4),ENERGY,Q,ENEXCI,EDLAB,AP,DM(4)
      DOUBLE PRECISION ETOT,A,B,C,D,Qtot,EOUT1,PHIMAX
      DOUBLE PRECISION EOUT2,APR,CMANG1,CMANG2

      DO I=1,4
         DM(I)=DN(I)*1000./931.49432
      ENDDO
      
      Qtot=Q-ENEXCI
      ETOT=ENERGY+Q
      
      APR=AP*0.01745329 
      
      A=DM(1)*DM(4)*(ENERGY/ETOT)/(DM(1)+DM(2))/(DM(3)+DM(4))
      B=DM(1)*DM(3)*(ENERGY/ETOT)/(DM(1)+DM(2))/(DM(3)+DM(4))
      C=DM(2)*DM(3)/(DM(1)+DM(2))/(DM(3)+DM(4))*
     *                              (1+DM(1)*Qtot/DM(2)/ETOT)
      D=DM(2)*DM(4)/(DM(1)+DM(2))/(DM(3)+DM(4))*
     *                              (1+DM(1)*Qtot/DM(2)/ETOT)
      
      
      IF (B.LE.D) THEN 
          EOUT1=ETOT*B*(DCOS(APR)+SQRT(D/B-(DSIN(APR))**2))**2
	  EOUT2=0.0
      ELSE
          PHIMAX=DASIN(DSQRT(D/B))
C          PRINT *,'Limit angle = ',PHIMAX
          IF (PHIMAX.GE.APR) THEN
             EOUT1=ETOT*B*(DCOS(APR)+SQRT(D/B-(DSIN(APR))**2))**2
	     EOUT2=ETOT*B*(DCOS(APR)-SQRT(D/B-(DSIN(APR))**2))**2
	  ELSE
	     PRINT *,'Attention! Limit angle was achieved.'
	     EOUT1=0.0
	     EOUT2=0.0
	  ENDIF 
	           
      ENDIF

      CMANG1=DASIN(SQRT(EOUT1/ETOT/D)*DSIN(APR))
      CMANG2=DASIN(SQRT(EOUT2/ETOT/D)*DSIN(APR))

      RETURN
      END	  
      
      
