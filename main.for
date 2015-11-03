      PROGRAM main
C    *  PROGRAM FOR GEANT Simulation of Experiment
C    *    in Dubna, 8He+p. Possibility to use liquid hydrogen
C    *    target as well as conventional gas target is included.
C    *    Interface changed to make modification of experimental
C    *    conditions more clear.
C    *    This new version was started on 21-11-01.
C    *    Author:  Grigory Rogachev
C    *    Present Address: University of Notre Dame, IN, USA.
C    *
       INTEGER NHBOOK,PAW
      PARAMETER (NGBANK=500000,NHBOOK=50000)
      COMMON /GCBANK/ Q(NGBANK)
      COMMON /PAWC/ PAW(NHBOOK)

C-->     Will terminate the job if runs more then 24 hours.
      CALL TIMEST(11186400.)

C-->     Initialises HBOOK and GEANT memory
      CALL GZEBRA(NGBANK)

      CALL HLIMIT(-NHBOOK)

C-->     Initialise the HIGZ and HPLOT Graphycal packages

C      CALL HPLINT(1)

C-->   User initialization program.

      CALL UGINIT

C-->     Start events processing

      CALL GRUN

C-->     End of run

      CALL UGLAST

      END


