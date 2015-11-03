      SUBROUTINE GUTREV
C----------------------------------------------------------------
C-
C-     Starting the tracking of the vertex
C-      Created by G.V.Rogachev 03-27-00
C-
C----------------------------------------------------------------
      INCLUDE 'gcflag.inc'
      INCLUDE 'gcbank.inc'
      integer pnum(7)
      DATA pnum /123,125,126,127,128,129,130/
      real er,fwhm,sigma,grn(10),newmas
      
	CALL RNORML(grn,10)

C-->  Each time new event starts new mass of unstable particles
C     is defined according
C     to avarage mass and width given in uginit.for

      DO i = 1,7 
       JPA = LQ(JPART-pnum(i))
       er = Q(JPA+11)
       fwhm = Q(JPA+10)
       sigma = fwhm/2.3
       newmas = er + sigma*grn(i)
       Q(JPA+7) = newmas
      ENDDO

      CALL GTREVE
      IF(IDEBUG.EQ.1.AND.ISWIT(2).EQ.1) THEN
         WRITE(*,*)'IN GUTREV'
      ENDIF
      END
