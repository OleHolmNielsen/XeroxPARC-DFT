      SUBROUTINE EXCORR (RHO,IEXCOR,NSPIN,VOLUM,ULA,EXFACT,XCE,XCMU)
C
C     Exchange-correlation functional used in routine "EXCH4"
C
C     Input:
C     RHO ...... charge density
C     IEXCOR is passed from K207:
C     0 means no XC at all
C     1 means Slater X-alpha
C     2 means Wigner interpolation
C     3 means Ceperley-Alder
C     NSPIN .... 1 for non-polarized, 2 for spin-polarized
C                (1=spin-up, 2=spin-down)
C     EXFACT ... The alpha-factor in Slater X-alpha (1 for Wigner etc.)
C
      REAL RHO(NSPIN), XCE(NSPIN), XCMU(NSPIN)
C
C     Associate LDF XC-functional names
      CHARACTER*20 LDF(0:3), POLAR(2)
      DOUBLE PRECISION BETA,ETA,XIL,PI43,ATRD,FTRD,TFTM
      LOGICAL INITLZ, SPNPOL
      REAL EX, MUX, EC, MUC, EXP, MUXP, ECP, MUCP
C
C.....FILES
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
C.....PHYSICAL AND MATHEMATICAL CONSTANTS
      DOUBLE PRECISION ABOHR,RYEV,RYDERG,PI,SPI
      COMMON /CONST/   ABOHR,RYEV,RYDERG,PI,SPI
      COMMON /WARN/ IWARN
C
      SAVE VOLAU,ATRD,FTRD,PI43,TFTM,TTRD,XFAC,SPNPOL,INITLZ
      SAVE A,B,C,D,G,B1,B2
      SAVE AP,BP,CP,DP,GP,B1P,B2P
      DATA INITLZ /.FALSE./
      DATA POLAR /'(non-polarized)', '(spin-polarized)'/
C
C-----------------------------------------------------------------------
C
      IF (.NOT. INITLZ) THEN
C
C       FIRST TIME AROUND - TELL WHICH FUNCTIONAL
C
        CALL LDFTYP (LDF)
        WRITE (IOUT,150) LDF(IEXCOR), POLAR(NSPIN)
150     FORMAT (' EXCORR - exchange-correlation is ',2A)
        VOLAU = VOLUM * ( ULA / ABOHR )**3
        ATRD  = 1.0D0 / 3.0D0
        FTRD  = 4.0D0 / 3.0D0
        PI43  = PI * FTRD
        TFTM  = 2.0D0 ** FTRD - 2.0D0
        TTRD  = 2.0D0 ** ATRD
        XFAC  = (2.0D0/PI) * ( 9.0D0*PI/4.0D0 )**ATRD
C
        IF (IEXCOR .EQ. 0) THEN
C         No XC
C
        ELSE IF (IEXCOR .EQ. 1) THEN
C         X-ALPHA
C         SCALE THE EXCHANGE TERM BY ALPHA AND 3/2:
          XFAC = EXFACT * XFAC * 3.0D0 / 2.0D0
C
        ELSE IF (IEXCOR .EQ. 2) THEN
C         WIGNER CORRELATION
          IF (NSPIN .NE. 1) THEN
            WRITE (IOUT,105)
105         FORMAT ('0EXCORR *** WARNING ***'/
     +        ' Wigner correlation is nonpolarized')
            IWARN = 1
            GOTO 400
            ENDIF
C
        ELSE IF (IEXCOR .EQ. 3) THEN
C         PARAMETERS FOR UNPOLARIZED GAS
          A =  0.0311
          B = -0.0480
          C =  0.0020
          D = -0.0116
          G = -0.1423
          B1 = 1.0529
          B2 = 0.3334
C         PARAMETERS FOR FULLY POLARIZED GAS
          AP =  0.01555
          BP = -0.0269
          CP =  0.0007
          DP = -0.0048
          GP = -0.0843
          B1P = 1.3981
          B2P = 0.2611
C
        ELSE
          WRITE (IOUT,*) 'EXCORR: illegal IEXCOR=',IEXCOR
          CALL EXIT
          ENDIF
C
        IF (NSPIN .EQ. 1) THEN
          SPNPOL = .FALSE.
        ELSE
          SPNPOL = .TRUE.
          ENDIF
        INITLZ = .TRUE.
C
        ENDIF
C
C-----------------------------------------------------------------------
C
      IF (IEXCOR .LE. 0) THEN
        DO 170 ISPIN = 1, NSPIN
          XCE (ISPIN) = 0.0
          XCMU(ISPIN) = 0.0
170       CONTINUE
        RETURN
        ENDIF
C
      IF (SPNPOL) THEN
        RHOAV = RHO(1) + RHO(2)
      ELSE
        RHOAV = RHO(1)
        ENDIF
C
C     GUARD AGAINST SLIGHTLY NEGATIVE DENSITIES
      IF (RHOAV .LT. -0.01) THEN
        WRITE (IOUT,*) 'EXCORR *** ERROR *** RHOAV =',RHOAV
        CALL EXIT
      ELSE IF (RHOAV .LE. 1.0E-9) THEN
        RHOAV = 1.0E-9
        ENDIF
C
      IF (SPNPOL) THEN
        Z     = (RHO(1) - RHO(2)) / RHOAV
        FZ    = ( (1.0+Z)**FTRD + (1.0-Z)**FTRD - 2.0 ) / TFTM
        FZP   = ( (1.0+Z)**ATRD - (1.0-Z)**ATRD ) / TFTM * FTRD
        ENDIF
C
C-----------------------------------------------------------------------
C
C     THE EXCHANGE POTENTIAL
C
      RS   = ( VOLAU / (RHOAV * PI43) )**ATRD
C     EXCHANGE-CORRELATION POTENTIAL:
      MUX  = - XFAC / RS
C     EXCHANGE-CORRELATION ENERGY:
      EX   = 0.75 * MUX
C     RELATIVISTIC CORRECTION FACTOR (MACDONALD AND VOSKO)
      BETA = 0.0140 / RS
      ETA  = DSQRT( 1.0D0 + BETA*BETA )
      XIL  = DLOG( BETA + ETA )
      MUX  = MUX * ( - 0.5D0 + 1.5D0 * XIL / (BETA*ETA) )
      EX   = EX  * (   1.0D0 - 1.5D0 * ( ETA/BETA - XIL/BETA**2 )**2 )
      IF (SPNPOL) THEN
        EXP  = TTRD * EX
        MUXP = TTRD * MUX
        ENDIF
C
C-----------------------------------------------------------------------
C
C     SLATER X-ALPHA DENSITY FUNCTIONAL
C
      IF (IEXCOR .EQ. 1) THEN
C
      EC  = 0.0
      MUC = 0.0
C
C-----------------------------------------------------------------------
C
C     WIGNER INTERPOLATION FORMULA
C
      ELSE IF (IEXCOR .EQ. 2) THEN
C
        EC  = - 0.88 / (RS + 7.8 )
        MUC = EC * (1.0 + RS/3.0/(RS + 7.8) )
C
C-----------------------------------------------------------------------
C
C     CEPERLEY-ALDER (SEE BACHELET ET AL., PHYS. REV. B 26, 4199 (1982).
C
      ELSE IF (IEXCOR .EQ. 3) THEN
C
      IF (RS .GE. 1.0) THEN
        SQRTRS = SQRT(RS)
        DENOM = 1.0 + B1*SQRTRS + B2*RS
        EC  = G / DENOM
        MUC = EC * (1.0 + 7.0/6.0*B1*SQRTRS + 4.0/3.0*B2*RS) / DENOM
      ELSE
        RSLOG = ALOG(RS)
        EC  = B + A*RSLOG + D*RS + C*RS*RSLOG
        MUC = (B - A/3.0) + A*RSLOG + (D+D-C)/3.0*RS +
     +        2.0/3.0*C*RS*RSLOG
        ENDIF
C     CONVERSION TO RYDBERGS:
      EC  = 2.0 * EC
      MUC = 2.0 * MUC
C
      IF (SPNPOL) THEN
C       THE FULLY POLARIZED DATA
        IF (RS .GE. 1.0) THEN
          DENOM = 1.0 + B1P*SQRTRS + B2P*RS
          ECP  = GP / DENOM
          MUCP = ECP * (1.0 + 7.0/6.0*B1P*SQRTRS + 4.0/3.0*B2P*RS) /
     +           DENOM
        ELSE
          ECP  = BP + AP*RSLOG + DP*RS + CP*RS*RSLOG
          MUCP = (BP - AP/3.0) + AP*RSLOG + (DP+DP-CP)/3.0*RS +
     +          2.0/3.0*CP*RS*RSLOG
          ENDIF
C       CONVERSION TO RYDBERGS:
        ECP  = 2.0 * ECP
        MUCP = 2.0 * MUCP
        ENDIF
C
C-----------------------------------------------------------------------
C
      ELSE
C
        WRITE (IOUT,140) IEXCOR
140     FORMAT('0SUBROUTINE EXCORR *** FATAL ERROR ***'/
     +    ' IEXCOR = ',I4,' IS ILLEGAL')
        CALL EXIT
        ENDIF
C
      IF (SPNPOL) THEN
        DXC = FZ * (EXP - EX) + FZ * (ECP - EC)
        EX  = EX  + FZ * (EXP - EX)
        EC  = EC  + FZ * (ECP - EC)
        MUX = MUX + FZ * (MUXP - MUX)
        MUC = MUC + FZ * (MUCP - MUC)
        DO 200 ISPIN = 1, NSPIN
C         ISPIN=1 IS UP, ISPIN=2 IS DOWN
          SIGN = FLOAT(3 - 2*ISPIN)
          XCE(ISPIN)  = EX  + EC
          XCMU(ISPIN) = MUX + MUC + DXC * FZP * (SIGN - Z)
200       CONTINUE
      ELSE
C       NON-POLARIZED
        XCE(1)  = EX  + EC
        XCMU(1) = MUX + MUC
        ENDIF
C
      RETURN
C
C-----------------------------------------------------------------------
C
C     ERROR MESSAGES
C
400   WRITE (IOUT,410)
410   FORMAT ('0EXCORR *** FATAL ERROR ***')
      WRITE (IOUT,420) NSPIN,IEXCOR
420   FORMAT ('0SPIN-POLARIZED (NSPIN=',I2,') NOT ALLOWED FOR IEXCOR =',
     +  I3)
      CALL EXIT
      RETURN
      END
