	REAL FUNCTION gtwght(cmang)
        real andis,sigmax,tval
        integer nang
        LOGICAL there
	COMMON /AD/ andis(1000,2),sigmax,nang

         DO i = 1,nang
    	   there = 
     &     ((cmang.GE.andis(i,1)).AND.(cmang.LT.andis((i+1),1)))
           IF (there) THEN 
             tval = andis(i,2) + (andis((i+1),2)-andis(i,2))*
     &       (cmang-andis(i,1))/
     &       (andis((i+1),1)-andis(i,1))
             GO TO 10
           ENDIF
         ENDDO

  10     CONTINUE
         gtwght =  tval/sigmax

        RETURN
	END         
