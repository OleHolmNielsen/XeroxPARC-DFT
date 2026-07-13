      SUBROUTINE EFERMI (DNEL,NSPIN,NBANDS,DEL,NSPPTS,NDIM8,NDIM10,
     +                   WEIGHT,OCC,EF,EIGVAL,SORT)
C
C-----------------------------------------------------------------------
C
C     WRITTEN BY RICHARD NEEDS ON 9TH DECEMBER 1983
C     GIVEN THE EIGENVALUES IN EIGVAL AND THE WEIGHTS OF THE
C     K-POINTS IN WEIGHT THIS SUBROUTINE CALCULATES THE FERMI LEVEL
C     EF AND THE OCCUPANCY OF THE STATES OCC.
C
C     METHOD: C-L FU AND K-M HO, PHYS. REV. B 28, 5480 (1983)
C             GAUSSIAN SMEARING OF EIGENVALUES WHEN COMPUTING OCCUPATION
C             NOTE: FOR SUM OF BANDS WE DO NOT SMEAR EIGENVALUES,
C                   AS WAS DONE BY FU AND HO. THE DIFFERENCE IS EASILY
C                   CALCULATED.
C
C     DNEL .... NUMBER OF ELECTRONS PER UNIT CELL (possibly fractional)
C     NSPIN ... 1 FOR NON-POLARIZED, 2 FOR SPIN-POLARIZED
C     NBANDS .. NUMBER OF BANDS FOR EACH K-POINT
C     DEL ..... WIDTH OF GAUSSIAN SMEARING FUNCTION
C     NSPPTS .. NUMBER OF K-POINTS
C     NDIM8 ... MAXIMUM NUMBER OF BANDS AT A K-POINT
C     WEIGHT .. THE WEIGHT OF EACH K-POINT
C     OCC ..... THE OCCUPANCY OF EACH STATE
C     EF ...... THE FERMI ENERGY
C     SORT .... THE EIGENVALUES ARE WRITTEN INTO SORT WHICH IS
C               THEN SORTED INTO ASCENDING NUMERICAL VALUE, FROM
C               WHICH BOUNDS ON EF CAN EASILY BE OBTAINED
C     EIGVAL .. IF LOWDIN CONTAINS 2ND ORDER EIGENVALUES IF NO
C               LOWDIN CONTAINS 1ST ORDER EIGENVALUES
C     NINC .... NO OF EIGENVALUE DIFFERENCES NEAR FERMI LEVEL FROM
C               WHICH THE AVERAGE BAND SEPARATION IS CALCULATED
C               NINC MUST BE ODD
C-----------------------------------------------------------------------
      DOUBLE PRECISION OCC(NDIM8,NSPIN,NSPPTS), WEIGHT(NSPPTS),
     +  EIGVAL(NDIM8,NSPIN,NSPPTS), DNEL
      DIMENSION SORT(NDIM8*NSPIN*NSPPTS)
      INTEGER NBANDS(NSPIN,NSPPTS)
      DOUBLE PRECISION EF,E1,E2,EUP,DEL,X,Z,Z1,Z2,DIV,DSOR,
     +                 TEST,WOCC,OCCMAX,ELOW
      DOUBLE PRECISION DERFC
      EXTERNAL DERFC
C.....FILES
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
C.....PHYSICAL AND MATHEMATICAL CONSTANTS
      DOUBLE PRECISION ABOHR,RYEV,RYDERG,PI,SPI
      COMMON /CONST/   ABOHR,RYEV,RYDERG,PI,SPI
C.....WARNINGS
      COMMON /WARN/ IWARN
      PARAMETER ( NINC = 9, NCYCLE = 20, NDIV = 9 )
C--------------------------------------------------------------------
C
        DSOR   = 1.0D0 / DBLE(NDIV)
        Z      = DNEL
        NEL    = DNEL + 1.0D-8
        IF (NSPIN .LT. 1 .OR. NSPIN .GT. 2) THEN
          WRITE (IOUT,*) 'EFERMI: NSPIN = ',NSPIN,' is illegal'
          CALL EXIT
          ENDIF
        IF (DEL .LT. 1.0D-6 .OR. DEL .GT. 50.0D0) THEN
          WRITE (IOUT,*) 'EFERMI: DEL = ', DEL, ' is nonsense'
          DEL = 0.01D0
          WRITE (IOUT,*) 'EFERMI: DEL is set = ', DEL
          ENDIF
C
C               COPY EIGVAL INTO SORT
        NEIG = 0
        DO 10 ISPPT = 1, NSPPTS
        DO 10 ISPIN = 1, NSPIN
          NBDS = NBANDS(ISPIN,ISPPT)
          DO 10  J = 1, NBDS
            NEIG = NEIG + 1
            SORT(NEIG) = EIGVAL(J,ISPIN,ISPPT)
10          CONTINUE
C
C               SORT THE ARRAY INTO ASCENDING ORDER OF EIGENVALUE
        CALL VSRTA (SORT,NEIG)
C
C               FIND AN UPPER BOUND E2 AND A LOWER BOUND E1 ON THE
C               FERMI ENERGY
        INT  = NEL * NSPPTS
        IOCC = 3 - NSPIN
        I1   = INT / IOCC + MOD(INT,IOCC) - 1
        I1   = MAX0(I1, 1)
        I2   = I1 + 2
        I2   = MIN0(I2, NEIG)
C
30      E1 =  DBLE(SORT(I1))
C
        Z1 = 0.0 D0
        DO 40 ISPPT = 1, NSPPTS
        WOCC = WEIGHT(ISPPT) / DBLE(NSPIN)
        DO 40 ISPIN = 1, NSPIN
          NBDS = NBANDS(ISPIN,ISPPT)
          DO 40 J = 1, NBDS
            X = (E1 - DBLE(EIGVAL(J,ISPIN,ISPPT)))/DEL
            Z1 = Z1 + WOCC * ( 2.0D0 - DERFC(X) )
40          CONTINUE
C
        IF (Z1 .GT. Z-1.0D-6) THEN
          I1 = I1 - 1
          IF (I1 .GT. 0) THEN
            GOTO 30
          ELSE
            I1 = 1
            E1 = DBLE(SORT(I1)) - DEL
            WRITE (IUNIT8,*) 'EFERMI *** lower bound on EF set to',E1
            GOTO 70
            ENDIF
          ENDIF
C
70      E2 = DBLE(SORT(I2))
C
        Z2 = 0.0 D0
        DO 60 ISPPT = 1, NSPPTS
        WOCC = WEIGHT(ISPPT) / DBLE(NSPIN)
        DO 60 ISPIN = 1, NSPIN
          NBDS = NBANDS(ISPIN,ISPPT)
          DO 80 J = 1,NBDS
            X = (E2 - DBLE(EIGVAL(J,ISPIN,ISPPT)))/DEL
            Z2 = Z2 + WOCC * ( 2.0D0 - DERFC(X) )
80          CONTINUE
60        CONTINUE
C
        IF (Z2 .LT. Z + 1.D-6) THEN
          I2 = I2 + 1
          IF (I2 .LE. NEIG) THEN
            GOTO 70
          ELSE
            I2 = NEIG
            E2 = DBLE(SORT(I2)) + DEL
            WRITE (IUNIT7,1100) E2
1100        FORMAT('0EFERMI *** WARNING ***'/
     +        ' Fermi level upper bound is > largest eigenvalue, E2=',
     +        F12.4,' EV')
            IWARN = 1
            ENDIF
          ENDIF
C
C               FIND FERMI ENERGY ENERGY BETWEEN BOUNDS E1 AND E2
        DO 90 ILOOP = 1,NCYCLE
          DIV = (E2-E1)*DSOR
          DO 100 I = 1,NDIV
            Z2 = 0.0 D0
            EUP = E1 + DBLE(I)*DIV
            DO 120 ISPPT = 1, NSPPTS
            WOCC = WEIGHT(ISPPT) / DBLE(NSPIN)
            DO 120 ISPIN = 1, NSPIN
              NBDS = NBANDS(ISPIN,ISPPT)
              DO 120 J = 1,NBDS
                X = (EUP - DBLE(EIGVAL(J,ISPIN,ISPPT)))/DEL
                Z2 = Z2 + WOCC * ( 2.0D0 - DERFC(X) )
120           CONTINUE
C           WHY WOULD YOU WANT THIS TEST? :
C           IF (Z1 .GT. Z) GOTO 1000
            IF (Z2 .GT. Z) GOTO 130
C               THIS TEST IS NEEDED TO CATCH SEMICONDUCTORS
            IF (Z2 .EQ. Z) GOTO 150
            Z1 = Z2
100         CONTINUE
          GOTO 1010
130       IF ( Z2-Z1 .LT. 1.0D-10) GOTO 160
          E1 = EUP - DIV
90        E2 = EUP
        WRITE(IOUT,*) ' * WARNING * fermi energy may not be accurate'
        WRITE(IOUT,140) NCYCLE
140     FORMAT(' After',I6,' cycles, required convergence not obtained')
        GOTO 160
C
C               WRITE OUT FERMI ENERGY
150     EF = EUP
        GOTO 170
160     EF = EUP - 0.5D0*DIV
170     WRITE (IOUT,180) EF
180     FORMAT (' Fermi energy = ',F16.8,' eV')
C
C               FORM OCCUPATIONS OCC(NBDS,NSPPTS,NSPIN)
C
        DO 190 ISPPT = 1, NSPPTS
        DO 190 ISPIN = 1, NSPIN
          NBDS = NBANDS(ISPIN,ISPPT)
          DO 190 J = 1, NBDS
            X = ( EF - DBLE(EIGVAL(J,ISPIN,ISPPT)) ) / DEL
            OCC(J,ISPIN,ISPPT) = (2.0D0 - DERFC(X)) / DBLE(NSPIN)
190       CONTINUE
C
C               TEST WHETHER OCCUPANCY ADDS UP TO Z
        TEST = 0.0D0
        DO 210 ISPPT = 1, NSPPTS
        DO 210 ISPIN = 1, NSPIN
          NBDS = NBANDS(ISPIN,ISPPT)
          DO 215 J = 1,NBDS
            TEST = TEST + WEIGHT(ISPPT) * OCC(J,ISPIN,ISPPT)
215         CONTINUE
210       CONTINUE
        IF ( DABS(TEST-Z) .GT. 1.0D-9) THEN
          WRITE(IOUT,*) ' *** WARNING ***'
          WRITE(IOUT,220) TEST,DNEL
220       FORMAT(' Sum of occupancies =',F20.12 ,' but NEL =',F20.12)
        ELSE
          WRITE(IUNIT8,230) TEST
230       FORMAT(' Sum of occupancies = ',F20.9)
          ENDIF
C
C               CALCULATE AVERAGE BAND SEPARATION NEAR EF
        DO 240 I = 1, NEIG
          IF (EF .LT. SORT(I) ) GOTO 250
240       CONTINUE
        I1 = MAX0( 1, I - 1 - NINC/2)
        I2 = MIN0( NEIG,  I + NINC/2)
250     SMEAR = (SORT(I2) - SORT(I1))/DBLE(I2 - I1 + 1)
        WRITE (IUNIT8,260) SMEAR
260     FORMAT( ' Average band separation near E-Fermi =',F8.5,' eV')
C
C       Test whether the material is a semiconductor
C
        IF ( MOD( NEL, 2) .EQ. 1) GOTO 350
        IF ( ABS(DNEL - DBLE(NEL)) .GT. 1.0D-8) GOTO 350
        INEL = NEL/2
        ELOW = 1.0D10
        DO 310 ISPPT = 1, NSPPTS
        DO 310 ISPIN = 1, NSPIN
          ELOW = DMIN1( ELOW, EIGVAL(INEL+1,ISPIN,ISPPT))
310       CONTINUE
        DO 320 ISPPT = 1,NSPPTS
        DO 320 ISPIN = 1,NSPIN
          IF (ELOW .LT. EIGVAL(INEL,ISPIN,ISPPT)) GOTO 350
320       CONTINUE
        WRITE (IOUT,*) 'This material may be a semiconductor.'
C
C       OCCUPATION OF A STATE: 1 OR 2, DEPENDING ON SPIN POLARIZATION
        OCCMAX = DBLE ( 3 - NSPIN )
        DO 330 ISPPT = 1, NSPPTS
        DO 330 ISPIN = 1, NSPIN
          DO 330 J = 1, INEL
            IF (ABS(OCC(J,ISPIN,ISPPT) - OCCMAX) .GT. 1.0E-8) GOTO 340
330         CONTINUE
        RETURN
C
340     WRITE (IOUT,345)
345     FORMAT(' However, gaussians at E-Fermi overlap,',
     +    ' so it may possibly be a metal as well')
        RETURN
C
350     WRITE (IOUT,*) 'This material is a metal.'
        RETURN
C
C----------------------------------------------------------------------
C               ERROR MESSAGES
C----------------------------------------------------------------------
1000  WRITE(IOUT,*) ' ***** ERROR *****'
      WRITE(IOUT,*) ' Fermi energy less than lower search bound set'
      EF = 1.0E20
      IWARN = 2
      RETURN
C
1010  WRITE(IOUT,*) ' ***** ERROR *****'
      WRITE(IOUT,*) ' Fermi energy greater than upper search bound set'
      EF = 1.0E20
      IWARN = 2
      RETURN
      END
