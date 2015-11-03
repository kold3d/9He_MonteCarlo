      real E1,X1,Y1,RAD1,E2,X2,Y2,RAD2,anincm,andis,endis,sigmax,rsmax
      real Gpre,LDEN,Aang,Apde,RDET
      real dtheta,theav,dRad,BeamE,dEner
      integer nang,nen,TState,ChName
      integer RType 
      COMMON /event/ E1,X1,Y1,RAD1,E2,X2,Y2,RAD2,anincm
      COMMON /URT/ RType
      COMMON /AD/ andis(1000,2),sigmax,nang
      COMMON /EN/ enkey(3)                   
      COMMON /EPAR/ ChName,TState,Gpre,LDEN,Apde,Aang,dtheta,
     &              theav,dRad(2),BeamE,dEner,RDET
