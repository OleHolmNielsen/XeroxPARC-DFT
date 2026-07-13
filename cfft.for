      SUBROUTINE CFFT (M,A,NDIMA,INV,S,MXNI,IFSET,IFERR)
C
C     Complex fast fourier transform (discrete)
C
      DIMENSION M(3), A(NDIMA), INV(MXNI), S(MXNI)
C     'A' is in reality a complex array of dimension (N1,N2,N3)
      INTEGER IFSET, IFERR
C
      DIMENSION N(3), NP(3), W(2), W2(2),W3(2)
      EQUIVALENCE (N1,N(1)), (N2,N(2)), (N3,N(3))
      SAVE MT, NT
C
      DATA PI /3.1415 92653 58979 32384 62643/
C-----------------------------------------------------------------------
C     Initialization:
      IF (IABS(IFSET) .GT. 1) GOTO 100
C
      MT = MAX0(M(1),M(2),M(3)) - 2
      IF (MT .LT. 1) THEN
        IFERR = 1
        RETURN
        ENDIF
      MT = MAX0(2,MT)
      IFERR = 0
      NT = 2**MT
      NTV2 = NT / 2
      PFNT2 = PI / FLOAT(2 * NT)
      DO 740 L = 1,NT
        S(L) = SIN(FLOAT(L) * PFNT2)
740     CONTINUE
      MTLEXP = NTV2
      LM1EXP = 1
      INV(1) = 0
      DO 760 L = 1,MT
        INV(LM1EXP + 1) = MTLEXP
CDIR$   IVDEP
        DO 750 J = 2,LM1EXP
          INV(J + LM1EXP) = INV(J) + MTLEXP
750       CONTINUE
        MTLEXP = MTLEXP / 2
        LM1EXP = LM1EXP * 2
760     CONTINUE
      IF (IFSET .EQ. 0) RETURN
C-----------------------------------------------------------------------
100   MTT = MAX0(M(1),M(2),M(3)) - 2
      IF (MTT .LT. 1 .OR. MTT .GT. MT) THEN
        IFERR = 1
        RETURN
        ENDIF
      MSUM = M(1) + M(2) + M(3)
      IFERR = 0
      ROOT2 = DSQRT(2.0D0)
      M1 = M(1)
      M2 = M(2)
      M3 = M(3)
      N1 = 2**M1
      N2 = 2**M2
      N3 = 2**M3
      IF (IFSET .GE. 0) THEN
        NX = N1 * N2 * N3
        FN = NX
        DO 160 I = 1,NX
          A(2*I - 1) =   A(2*I - 1) / FN
          A(2*I)     = - A(2*I)     / FN
160       CONTINUE
        ENDIF
      NP(1) = N1 * 2
      NP(2) = NP(1) * N2
      NP(3) = NP(2) * N3
      DO 430 ID = 1,3
        IL = NP(3) - NP(ID)
        IL1 = IL + 1
        MI = M(ID)
        IF (MI .LE. 0) GOTO 430
        IDIF = NP(ID)
        KBIT = NP(ID)
        MEV = 2 * (MI / 2)
        IF (MI .LE. MEV) GOTO 220
        KBIT = KBIT / 2
        KL = KBIT - 2
        DO 200 I = 1,IL1,IDIF
          KLAST = KL + I
CDIR$     NEXTSCALAR
          DO 200 K = I,KLAST,2
            KD = K + KBIT
            T = A(KD)
            A(KD) = A(K) - T
            A(K) = A(K) + T
            T = A(KD + 1)
            A(KD + 1) = A(K + 1) - T
            A(K + 1) = A(K + 1) + T
200         CONTINUE
        IF (MI .LE. 1) GOTO 430
        LFIRST = 3
        JLAST = 1
        GO TO 230
C
220     LFIRST = 2
        JLAST = 0
230     DO 420 L = LFIRST,MI,2
            JJDIF = KBIT
            KBIT = KBIT / 4
            KL = KBIT - 2
            DO 240 I = 1,IL1,IDIF
              KLAST = I + KL
CDIR$         NEXTSCALAR
              DO 240 K = I,KLAST,2
                K1 = K + KBIT
                K2 = K1 + KBIT
                K3 = K2 + KBIT
                T = A(K2)
                A(K2) = A(K) - T
                A(K) = A(K) + T
                T = A(K2 + 1)
                A(K2 + 1) = A(K + 1) - T
                A(K + 1) = A(K + 1) + T
                T = A(K3)
                A(K3) = A(K1) - T
                A(K1) = A(K1) + T
                T = A(K3 + 1)
                A(K3 + 1) = A(K1 + 1) - T
                A(K1 + 1) = A(K1 + 1) + T
                T = A(K1)
                A(K1) = A(K) - T
                A(K) = A(K) + T
                T = A(K1 + 1)
                A(K1 + 1) = A(K + 1) - T
                A(K + 1) = A(K + 1) + T
                R =  - A(K3 + 1)
                T = A(K3)
                A(K3) = A(K2) - R
                A(K2) = A(K2) + R
                A(K3 + 1) = A(K2 + 1) - T
                A(K2 + 1) = A(K2 + 1) + T
240             CONTINUE
              IF (JLAST .LE. 0) GOTO 410
              JJ = JJDIF + 1
              ILAST = IL + JJ
              DO 260 I = JJ,ILAST,IDIF
                KLAST = KL + I
CDIR$           NEXTSCALAR
                DO 260 K = I,KLAST,2
                  K1 = K + KBIT
                  K2 = K1 + KBIT
                  K3 = K2 + KBIT
                  R =  - A(K2 + 1)
                  T = A(K2)
                  A(K2) = A(K) - R
                  A(K) = A(K) + R
                  A(K2 + 1) = A(K + 1) - T
                  A(K + 1) = A(K + 1) + T
                  AWR = A(K1) - A(K1 + 1)
                  AWI = A(K1 + 1) + A(K1)
                  R =  - A(K3) - A(K3 + 1)
                  T = A(K3) - A(K3 + 1)
                  A(K3) = (AWR - R) / ROOT2
                  A(K3 + 1) = (AWI - T) / ROOT2
                  A(K1) = (AWR + R) / ROOT2
                  A(K1 + 1) = (AWI + T) / ROOT2
                  T = A(K1)
                  A(K1) = A(K) - T
                  A(K) = A(K) + T
                  T = A(K1 + 1)
                  A(K1 + 1) = A(K + 1) - T
                  A(K + 1) = A(K + 1) + T
                  R =  - A(K3 + 1)
                  T = A(K3)
                  A(K3) = A(K2) - R
                  A(K2) = A(K2) + R
                  A(K3 + 1) = A(K2 + 1) - T
                  A(K2 + 1) = A(K2 + 1) + T
260               CONTINUE
            IF (JLAST .LE. 1) GOTO 410
            JJ = JJ + JJDIF
            DO 400 J = 2,JLAST
              I = INV(J + 1)
              IC = NT - I
              W(1) = S(IC)
              W(2) = S(I)
              I2 = 2 * I
              I2C = NT - I2
              IF (I2C) 300,290,280
280           W2(1) = S(I2C)
              W2(2) = S(I2)
              GO TO 310
290           W2(1) = 0.
              W2(2) = 1.
              GO TO 310
300           I2CC = I2C + NT
              I2C =  - I2C
              W2(1) =  - S(I2C)
              W2(2) = S(I2CC)
310           I3 = I + I2
              I3C = NT - I3
              IF (I3C) 340,330,320
320           W3(1) = S(I3C)
              W3(2) = S(I3)
              GO TO 380
330           W3(1) = 0.
              W3(2) = 1.
              GO TO 380
340           I3CC = I3C + NT
              IF (I3CC) 370,360,350
350           I3C =  - I3C
              W3(1) =  - S(I3C)
              W3(2) = S(I3CC)
              GO TO 380
360           W3(1) =  - 1.
              W3(2) = 0.
              GO TO 380
370           I3CCC = NT + I3CC
              I3CC =  - I3CC
              W3(1) =  - S(I3CCC)
              W3(2) =  - S(I3CC)
380           ILAST = IL + JJ
              DO 390 I = JJ,ILAST,IDIF
                KLAST = KL + I
CDIR$           NEXTSCALAR
                DO 390 K = I,KLAST,2
                  K1 = K + KBIT
                  K2 = K1 + KBIT
                  K3 = K2 + KBIT
                  R = A(K2) * W2(1) - A(K2 + 1) * W2(2)
                  T = A(K2) * W2(2) + A(K2 + 1) * W2(1)
                  A(K2) = A(K) - R
                  A(K) = A(K) + R
                  A(K2 + 1) = A(K + 1) - T
                  A(K + 1) = A(K + 1) + T
                  R = A(K3) * W3(1) - A(K3 + 1) * W3(2)
                  T = A(K3) * W3(2) + A(K3 + 1) * W3(1)
                  AWR = A(K1) * W(1) - A(K1 + 1) * W(2)
                  AWI = A(K1) * W(2) + A(K1 + 1) * W(1)
                  A(K3) = AWR - R
                  A(K3 + 1) = AWI - T
                  A(K1) = AWR + R
                  A(K1 + 1) = AWI + T
                  T = A(K1)
                  A(K1) = A(K) - T
                  A(K) = A(K) + T
                  T = A(K1 + 1)
                  A(K1 + 1) = A(K + 1) - T
                  A(K + 1) = A(K + 1) + T
                  R =  - A(K3 + 1)
                  T = A(K3)
                  A(K3) = A(K2) - R
                  A(K2) = A(K2) + R
                  A(K3 + 1) = A(K2 + 1) - T
390               A(K2 + 1) = A(K2 + 1) + T
400           JJ = JJDIF + JJ
410         JLAST = 4 * JLAST + 3
420         CONTINUE
430     CONTINUE
C
      NTSQ = NT * NT
      M3MT = M3 - MT
      IF (M3MT .LT. 0) GOTO 450
      IGO3 = 1
      N3VNT = N3 / NT
      MINN3 = NT
      GO TO 460
450   IGO3 = 2
      N3VNT = 1
      NTVN3 = NT / N3
      MINN3 = N3
460   JJD3 = NTSQ / N3
      M2MT = M2 - MT
      IF (M2MT .LT. 0) GOTO 480
      IGO2 = 1
      N2VNT = N2 / NT
      MINN2 = NT
      GO TO 490
480   IGO2 = 2
      N2VNT = 1
      NTVN2 = NT / N2
      MINN2 = N2
490   JJD2 = NTSQ / N2
      M1MT = M1 - MT
      IF (M1MT .LT. 0) GOTO 510
      IGO1 = 1
      N1VNT = N1 / NT
      MINN1 = NT
      GO TO 520
510   IGO1 = 2
      N1VNT = 1
      NTVN1 = NT / N1
      MINN1 = N1
520   JJD1 = NTSQ / N1
      JJ3 = 1
      J = 1
C
      DO 660 JPP3 = 1,N3VNT
        IPP3 = INV(JJ3)
        DO 650 JP3 = 1,MINN3
          IF (IGO3 .EQ. 2) GOTO 540
C         GO TO (530,540),IGO3
530       IP3 = INV(JP3) * N3VNT
          GO TO 550
540       IP3 = INV(JP3) / NTVN3
550       I3 = (IPP3 + IP3) * N2
          JJ2 = 1
          DO 650 JPP2 = 1,N2VNT
            IPP2 = INV(JJ2) + I3
            DO 640 JP2 = 1,MINN2
              IF (IGO2 .EQ. 2) GOTO 570
C             GO TO (560,570),IGO2
560           IP2 = INV(JP2) * N2VNT
              GO TO 580
570           IP2 = INV(JP2) / NTVN2
580           I2 = (IPP2 + IP2) * N1
              JJ1 = 1
              DO 640 JPP1 = 1,N1VNT
                IPP1 = INV(JJ1) + I2
                DO 630 JP1 = 1,MINN1
                  IF (IGO1 .EQ. 2) GOTO 600
C                 GO TO (590,600),IGO1
590               IP1 = INV(JP1) * N1VNT
                  GO TO 610
600               IP1 = INV(JP1) / NTVN1
610               I = 2 * (IPP1 + IP1) + 1
                  IF (J .GE. I) GOTO 630
                  T = A(I)
                  A(I) = A(J)
                  A(J) = T
                  T = A(I + 1)
                  A(I + 1) = A(J + 1)
                  A(J + 1) = T
630               J = J + 2
640             JJ1 = JJ1 + JJD1
650         JJ2 = JJ2 + JJD2
660     JJ3 = JJ3 + JJD3
C
      IF (IFSET .LE. 0)                   RETURN
670   DO 680 I = 1,NX
        A(2 * I) =  - A(2 * I)
680     CONTINUE
      RETURN
      END
