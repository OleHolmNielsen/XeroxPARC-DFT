      SUBROUTINE VNLKK (LSTRES,GK,I1,I2,IROW,
     +  POTPAR,POTTYP,VNLIJ,NDIM5,NDIM20,NLMAX,NPARAM,
     +  STRFAC,IPLACE,NTYPMX,NDIM1,VOLUM,ULA,LREAL,
     +  COSINE,FP,FPP,
     +  FP1,FP2,XX,
     +  IGK1,V0,VX,
     +  VY,VXX,VYY,
     +  VXY,UPRIM,STRESS,VK1K2)
C
C     CALCULATE THE  N O N - L O C A L  CONTRIBUTION TO THE
C     MATRIX ELEMENT V(K,K') .
C
C     FOR NUMERICALLY GENERATED POTENTIALS.
C     (OPTIMIZED FOR VECTORIZATION)
C
C     WRITTEN    23-MAY-1984 BY OLE HOLM NIELSEN (NORDITA, COPENHAGEN)
C     REWRITTEN   6-JUL-1985 FOR CRAY VECTORIZATION
C
C     INPUT:
C     LSTRES ........ LOGICAL SWITCH FOR CALCULATING STRESS TENSOR
C     GK ............ THE LENGTHS AND VECTORS FOR(K+G)
C     I1 ............ FIRST ELEMENT OF GK CONSIDERED IN THIS CALL
C     I2 ............ LAST  ELEMENT OF GK CONSIDERED IN THIS CALL
C     IROW .......... ROW OF HAMILTONIAN IN THIS CALL
C     POTPAR ........ POTENTIAL PARAMETERS - FOR DEFINITIONS
C     SEE THE WAY POTPAR IS CONSTRUCTED IN THE
C     SETTING-UP ROUTINE K213 (OR A SIMILAR ONE)
C     VNLIJ ......... NONLOCAL CORRECTION POTENTIAL ON A UNIFORM
C                     2-DIM MESH IN (K,K')-SPACE.
C                     VNLIJ(K,K') = 4*PI * F(K,K').
C                     DIMENSION (0:NDIM20) FOR INTERPOLATION
C
C                     ALL POTENTIALS MUST BE IN RYDBERGS
C
C     POTTYP ........ POTENTIAL TYPE: MUST BE 20 FOR THIS ROUTINE
C     NPARAM ........ NUMBER OF PARAMETERS
C     STRFAC ........ ATOMIC STRUCTURE FACTORS
C     IPLACE ........ LOCATION OF (G-G') IN LIST OF STRFAC,
C     AS DETERMINED BY A PREVIOUS LOOKUP
C     ULA ........... LATTICE CONSTANT
C     COSINE-VXY..... WORKSPACES
C     OUTPUT:
C     UPRIM ......... MATRIX ELEMENTS FOR TYPES OF ATOMS (FOR FORCES)
C     STRESS ........ STRESS MATRIX ELEMENT (UPPER TRIANGLE)
C     VK1K2 ......... V(K+G,K+G')
C
C     METHOD: QUADRATIC INTERPOLATION OF (K1,K2) IN A SQUARE MESH
C
      PARAMETER ( LMAX = 2 )
      PARAMETER ( NTIPMX = 10 )
C     LMAX MUST NOT BE INCREASED OVER 2 WITHOUT REPROGRAMMING
C
C.....FILES
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
C.....PHYSICAL AND MATHEMATICAL CONSTANTS
      DOUBLE PRECISION ABOHR,RYEV,RYDERG,PI,SPI
      COMMON /CONST/   ABOHR,RYEV,RYDERG,PI,SPI
C
C.....WARNINGS
      COMMON /WARN/ IWARN
C
      DOUBLE PRECISION POTPAR(NTYPMX,NPARAM)
      REAL VNLIJ(0:NDIM20,0:NDIM20,NLMAX,NTYPMX)
      REAL GK(4,NDIM5),UPRIM(NTYPMX,NDIM5)
      INTEGER IPLACE(NDIM5)
      INTEGER POTTYP(NTYPMX)
      COMPLEX STRFAC(NTYPMX,NDIM1)
      REAL ULA,VOLUM
      DIMENSION STRESS(6,NDIM5),VK1K2(NDIM5)
      DIMENSION COSINE(NDIM5),FP(NDIM5),FPP(NDIM5),FP1(NDIM5),
     +FP2(NDIM5),XX(NDIM5),IGK1(NDIM5),V0(NDIM5),VX(NDIM5),
     +VY(NDIM5),VXX(NDIM5),VYY(NDIM5),VXY(NDIM5)
      DIMENSION LO(1+LMAX,NTIPMX)
      DIMENSION GCONV(NTIPMX),NQNL(NTIPMX)
      LOGICAL INIT,LSTRES,LREAL
C     KRONECKER DELTA
      DIMENSION DIJ(6)
C     INDEXING UPPER TRIANGLE OF STRESS TENSOR:
      INTEGER III(6),JJJ(6)
C
      SAVE LO, GCONV, NQNL, III, JJJ, DIJ, VOL, INIT
C----------------------------------------------------------------------
C     FOR COMPLEX HAMILTONIANS
C     COMPLEX STRESS,VK1K2
C     COMPLEX CMX
C     CMX(X1,X2) = CMPLX(X1,X2)
C----------------------------------------------------------------------
C     FOR REAL HAMILTONIANS
      REAL STRESS,VK1K2
      REAL CMX
      CMX(X1,X2) = X1
C----------------------------------------------------------------------
      DATA III /1,1,1, 2,2, 3/, JJJ /1,2,3, 2,3, 3/
C     KRONECKER DELTA:
      DATA DIJ /1.0,0.0,0.0, 1.0,0.0, 1.0/
      DATA INIT /.FALSE./
C----------------------------------------------------------------------
C     FOR COMPLEX HAMILTONIANS
C     IF (LREAL) WRITE(IOUT,*) '***VNLKK*** LREAL = ',LREAL
C     FOR REAL HAMILTONIANS
      IF (.NOT. LREAL) WRITE(IOUT,*) '***VNLKK*** LREAL = ',LREAL
C----------------------------------------------------------------------
C
      IF (INIT) GOTO 140
C
C     VOLUME IN ATOMIC UNITS
      VOL    = DBLE(VOLUM)*(DBLE(ULA)/ABOHR)**3
C     CONVERSION OF G INTO A.U.
      CONV   = 2.0D0*PI*ABOHR/DBLE(ULA)
      IF (NTYPMX .GT. NTIPMX) THEN
        WRITE (IOUT,110) NTYPMX,NTIPMX
110     FORMAT(' VNLKK - NTYPMX=',I5,' .GT. INTERNAL DIMENSION NTIPMX='
     +    ,I5,' - FIX IT')
        CALL EXIT
        ENDIF
C
C     CHECK POTENTIAL PARAMETERS
C
      DO 120 ITYPE = 1,NTYPMX
        IF (POTTYP(ITYPE) .EQ.  0) GOTO 120
C       INCORRECT POTENTIAL TYPE :
        IF (POTTYP(ITYPE) .NE. 20) THEN
          WRITE (IOUT,100) POTTYP
100       FORMAT('0SUBROUTINE VNLKK *** WARNING ***'/
     +      ' POTENTIAL TYPES OTHER THAN "20" ENCOUNTERED; POTTYP IS'/
     +      1X,16I5/)
          IWARN = 2
          GOTO 120
          ENDIF
C
C       ANGULAR MOMENTUM
        DO 115 NL = 1,NLMAX
          LO(NL,ITYPE) = NINT(POTPAR(ITYPE,20+NL))
115       CONTINUE
        NQNL(ITYPE) = NINT(POTPAR(ITYPE,8))
        DELQNL      =      POTPAR(ITYPE,9)
        GCONV(ITYPE) = CONV/DELQNL
120     CONTINUE
C
      INIT = .TRUE.
C
C-----------------------------------------------------------------------
C
C     THE ANGLE BETWEEN K AND K'
140   DO 150 M = I1,I2
        COSINE(M) = ( GK(1,M)*GK(1,IROW) +
     +                GK(2,M)*GK(2,IROW) +
     +                GK(3,M)*GK(3,IROW) ) /
     +              ( GK(4,M)*GK(4,IROW) + 1.0E-30 )
        VK1K2(M) = 0.0
        FP(M) = 0.0
150     CONTINUE
C
      IF (LSTRES) THEN
        DO 160 M = I1,I2
          FPP(M) = 0.0
          FP1(M) = 0.0
          FP2(M) = 0.0
160       CONTINUE
        DO 165 M = I1,I2
          STRESS(1,M) = 0.0
          STRESS(2,M) = 0.0
          STRESS(3,M) = 0.0
          STRESS(4,M) = 0.0
          STRESS(5,M) = 0.0
          STRESS(6,M) = 0.0
165       CONTINUE
        ENDIF
C
      DO 400 ITYPE = 1,NTYPMX
C
C     NO ATOM OF THIS TYPE
      IF (POTTYP(ITYPE) .NE. 20) GOTO 400
C
C     XX AND YY ARE THE G1,2-VECTOR LENGTHS IN A.U.,DIVIDED BY DELQNL
      DO 170 M = I1,I2
        XX(M)   = GK(4,M)*GCONV(ITYPE)
        IGK1(M) = MIN0(NINT(XX(M)) + 1,NQNL(ITYPE) - 1)
        XX(M)   = XX(M) - IGK1(M) + 1
170     CONTINUE
C
C     TRY TO DO: SKIP TERM IF IGK1,2 >=NQNL
      YY   = GK(4,IROW)*GCONV(ITYPE)
      IGK2 = MIN0( NINT(YY) + 1,NQNL(ITYPE) - 1)
      YY   = YY - IGK2 + 1
C
      DO 300 NL = 1,NLMAX
C
      IF (LO(NL,ITYPE) .LT. 0) GOTO 300
C
C     DO QUADRATIC INTERPOLATION IN ARRAY 'VNLIJ'
C     THE FOUR EDGES AND THE CENTRE OF THE SQUARE ARE INTERPOLATED
C     EXACTLY, AND THE FOUR CORNERS ARE FITTED.
C     THE ANALYTICAL FORMULA IS PROGRAMMED BELOW:
C
      DO 180 M = I1,I2
        V0(M)  = VNLIJ(IGK1(M)  ,IGK2  ,NL,ITYPE)
        VX(M)  = VNLIJ(IGK1(M)+1,IGK2  ,NL,ITYPE) -
     +           VNLIJ(IGK1(M)-1,IGK2  ,NL,ITYPE)
        VY(M)  = VNLIJ(IGK1(M)  ,IGK2+1,NL,ITYPE) -
     +           VNLIJ(IGK1(M)  ,IGK2-1,NL,ITYPE)
        VXX(M) = VNLIJ(IGK1(M)+1,IGK2  ,NL,ITYPE) +
     +           VNLIJ(IGK1(M)-1,IGK2  ,NL,ITYPE)
     +      -2.0*VNLIJ(IGK1(M)  ,IGK2  ,NL,ITYPE)
        VYY(M) = VNLIJ(IGK1(M)  ,IGK2+1,NL,ITYPE) +
     +           VNLIJ(IGK1(M)  ,IGK2-1,NL,ITYPE)
     +      -2.0*VNLIJ(IGK1(M)  ,IGK2  ,NL,ITYPE)
        VXY(M) = 0.25 * (
     +           VNLIJ(IGK1(M)+1,IGK2+1,NL,ITYPE) -
     +           VNLIJ(IGK1(M)-1,IGK2+1,NL,ITYPE) +
     +           VNLIJ(IGK1(M)-1,IGK2-1,NL,ITYPE) -
     +           VNLIJ(IGK1(M)+1,IGK2-1,NL,ITYPE) )
180     CONTINUE
C
C     SUM UP CONTRIBUTION:
C
      IF (LSTRES) THEN
C
      IF (LO(NL,ITYPE) .EQ. 0) THEN
C       L = 0
        DO 190 M = I1,I2
          FP(M)  = FP(M) + V0(M) + XX(M)*YY * VXY(M) +
     +      0.5 * (XX(M) * VX(M) + YY * VY(M) +
     +      XX(M)**2 * VXX(M) + YY**2 * VYY(M))
          FPP(M) = FPP(M) + 0
          FP1(M) = FP1(M) + 0.5*VX(M) + XX(M) * VXX(M) + YY * VXY(M)
          FP2(M) = FP2(M) + 0.5*VY(M) + YY * VYY(M) + XX(M) * VXY(M)
190       CONTINUE
      ELSE IF (LO(NL,ITYPE) .EQ. 1) THEN
C       L = 1
        DO 200 M = I1,I2
          FK1K2 = V0(M) + XX(M)*YY * VXY(M) +
     +      0.5 * (XX(M) * VX(M) + YY * VY(M) +
     +      XX(M)**2 * VXX(M) + YY**2 * VYY(M))
          FP(M)  = FP(M) + COSINE(M) * FK1K2
          FPP(M) = FPP(M) + FK1K2
          FP1(M) = FP1(M) +
     +      COSINE(M) * ( 0.5*VX(M) + XX(M) * VXX(M) + YY * VXY(M) )
          FP2(M) = FP2(M) +
     +      COSINE(M) * ( 0.5*VY(M) + YY * VYY(M) + XX(M) * VXY(M) )
200       CONTINUE
      ELSE IF (LO(NL,ITYPE) .EQ. 2) THEN
C       L = 2
        DO 210 M = I1,I2
          FK1K2 = V0(M) + XX(M)*YY * VXY(M) +
     +      0.5 * (XX(M) * VX(M) + YY * VY(M) +
     +      XX(M)**2 * VXX(M) + YY**2 * VYY(M))
          PL2    = 1.5 * COSINE(M)**2 - 0.5
          FP(M)  = FP(M)+ FK1K2 * PL2
          FPP(M) = FPP(M) + FK1K2 * 3.0 * COSINE(M)
          FP1(M) = FP1(M) +
     +      PL2 * ( 0.5 * VX(M) + XX(M) * VXX(M) + YY * VXY(M))
          FP2(M) = FP2(M) +
     +      PL2 * ( 0.5 * VY(M) + YY * VYY(M) +XX(M) * VXY(M) )
210       CONTINUE
        ENDIF
C
      ELSE
C
        IF(LO(NL,ITYPE) .EQ. 0) THEN
C         L = 0
          DO 220 M = I1,I2
            FP(M) = FP(M) + V0(M) + XX(M)*YY * VXY(M) +
     +        0.5 * (XX(M) * VX(M) + YY * VY(M) +
     +        XX(M)**2 * VXX(M) + YY**2 * VYY(M))
220         CONTINUE
        ELSE IF (LO(NL,ITYPE) .EQ. 1) THEN
C         L = 1
          DO 230 M = I1,I2
            FP(M) = FP(M) + COSINE(M) * ( V0(M) + XX(M)*YY * VXY(M) +
     +        0.5 * (XX(M) * VX(M) + YY * VY(M) +
     +        XX(M)**2 * VXX(M) + YY**2 * VYY(M)))
230         CONTINUE
        ELSE IF (LO(NL,ITYPE) .EQ. 2) THEN
C         L = 2
          DO 240 M = I1,I2
            FP(M)  = FP(M) + (1.5*COSINE(M)**2 - 0.5) *
     +        (V0(M) + XX(M)*YY * VXY(M) +
     +        0.5 * (XX(M) * VX(M) + YY * VY(M) +
     +        XX(M)**2 * VXX(M) + YY**2 * VYY(M)))
240         CONTINUE
          ENDIF
C
        ENDIF
C
300   CONTINUE
C
C     NEXT ATOMIC TYPE - SUM UP CONTRIBUTION
C
C     GATHER STRUCTURE FACTORS
      DO 305 M = I1,I2
        VX(M) = REAL ( STRFAC(ITYPE,IPLACE(M)) )
        VY(M) = AIMAG( STRFAC(ITYPE,IPLACE(M)) )
305     CONTINUE
      DO 310 M = I1,I2
        UPRIM(ITYPE,M) = FP(M) / VOL
        VK1K2(M) = VK1K2(M)+ FP(M) * CMX(VX(M),VY(M))
310     CONTINUE
C
      IF (LSTRES) THEN
C       MULTIPLY DERIVATIVES BY -LENGTH(K+G)
        DO 320 M = I1,I2
          FP1(M) = - FP1(M) * GK(4,  M ) * GCONV(ITYPE)
          FP2(M) = - FP2(M) * GK(4,IROW) * GCONV(ITYPE)
320       CONTINUE
        DO 350 K = 1,6
          I = III(K)
          J = JJJ(K)
          GK2IJ = GK(I,IROW) * GK(J,IROW) / ( GK(4,IROW)**2 + 1.E-30)
          DO 350 M = I1,I2
            FPPCOS      = FPP(M) * COSINE(M)
            STRESS(K,M) = STRESS(K,M) + CMX(VX(M),VY(M)) * (
     +        (FP1(M) + FPPCOS) *
     +          GK(I,M) * GK(J,M) / (GK(4,M)**2 + 1.E-30)
     +      + (FP2(M) + FPPCOS) * GK2IJ
     +      - FPP(M) *
     +         (GK(I,M) * GK(J,IROW) + GK(J,M)*GK(I,IROW)) /
     +         (GK(4,M) * GK(4,IROW) + 1.E-30)
     +      - FP(M) * DIJ(K))
350         CONTINUE
          DO 360 M = I1,I2
            FPP(M) = 0.0
            FP1(M) = 0.0
            FP2(M) = 0.0
360         CONTINUE
          ENDIF
C
        DO 370 M = I1,I2
          FP(M) = 0.0
370       CONTINUE
C
400     CONTINUE
C
      DO 410 M = I1,I2
        VK1K2(M) = VK1K2(M) / VOL
410     CONTINUE
      IF (LSTRES) THEN
        DO 420 M = I1,I2
          STRESS(1,M) = STRESS(1,M) / VOL
          STRESS(2,M) = STRESS(2,M) / VOL
          STRESS(3,M) = STRESS(3,M) / VOL
          STRESS(4,M) = STRESS(4,M) / VOL
          STRESS(5,M) = STRESS(5,M) / VOL
          STRESS(6,M) = STRESS(6,M) / VOL
420      CONTINUE
       ENDIF
C
      RETURN
      END
