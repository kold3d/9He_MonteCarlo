      DOUBLE PRECISION FUNCTION HQUADF(V)
      DOUBLE PRECISION V( 1)
      INTEGER NPAR, NDIM, IMQFUN, I, J
      DOUBLE PRECISION HQDJ, VV, VCONST
      DOUBLE PRECISION SIGVMI( 1), SIGVT( 1)
      DOUBLE PRECISION SIGV(   6, 1)
      DOUBLE PRECISION SIGDEL(   6)
      DOUBLE PRECISION SIGA(   6)
      DATA NPAR, NDIM, IMQFUN /    6,    1,    2/
      DATA VCONST / 0.2000000029802    /
      DATA SIGVMI /              0.    /
      DATA SIGVT /  50.00000000000    /
      DATA SIGV / 0.1960000097752    
     +,              0.    
     +,  1.000000000000    
     +, 0.2160000056028    
     +, 0.4640000164509    
     +, 0.2880000174046    
     +/
      DATA SIGDEL / 0.9600000455976E-02
     +, 0.1200000042445E-05
     +, 0.1200000042445E-05
     +, 0.1920000091195E-01
     +, 0.3840000182390E-01
     +, 0.7680000364780E-01
     +/
      DATA SIGA / -776.0514094320    
     +,  500.5256587613    
     +, -65.21071932027    
     +,  65.11645664725    
     +, -145.3985894830    
     +,  5.715318610455    
     +/
      HQUADF = 0.
      DO 20 J = 1, NPAR
         HQDJ = 0.
         DO 10 I = 1, NDIM
            VV = (V (I) - SIGVMI (I)) / SIGVT (I)
            HQDJ = HQDJ + (VV - SIGV (J, I)) ** 2
   10    CONTINUE
         HQDJ = HQDJ + SIGDEL (J) ** 2
         HQDJ = SQRT (HQDJ)
         HQUADF = HQUADF + SIGA (J) * HQDJ
   20 CONTINUE
      IF (IMQFUN .EQ. 2) HQUADF = VCONST * EXP (HQUADF)
      END
