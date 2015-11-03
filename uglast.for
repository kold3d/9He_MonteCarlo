      SUBROUTINE UGLAST

      CALL GLAST
 
C      CALL HCDIR('//LUN1',' ')
C      CALL HRPUT(0,FNAME,'N')

      CALL HROUT(0,ICYCLE,' ')
      CALL HREND('LUN1')
      CLOSE (UNIT=1)
C       pause 'Look at the picture!' 
      RETURN

      END




