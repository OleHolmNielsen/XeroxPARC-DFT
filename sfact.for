      SUBROUTINE SFACT (STRFAC,SATGR,SATGI,TY,XKAPA,DX,IGLIST,
     +  B1,B2,B3,NDIM1,NDIM9,NTYPMX,NAT,NGTOT,NG1,NG2,NG3,LREAL)
C
C     ATOMIC STRUCTURE FACTORS
C
      LOGICAL LREAL
      COMPLEX STRFAC(NTYPMX,NDIM1)
      REAL SATGR(NDIM9,NDIM1),SATGI(NDIM9,NDIM1)
      REAL XKAPA(3,NDIM9),B1(3),B2(3),B3(3)
      DOUBLE PRECISION DX(3,NDIM9),ARGMT,SATR,SATI,D1,D2,D3
      INTEGER IGLIST(3,NDIM1),TY(NDIM9)
C
C.....FILES
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
      DOUBLE PRECISION ABOHR,RYEV,RYDERG,PI,SPI
      COMMON /CONST/   ABOHR,RYEV,RYDERG,PI,SPI
C---------------------------------------------------------------------
C
      CALL USAGE (NP0,T0)
C
      DO 5 KAPA = 1,NTYPMX
      DO 5 IG   = 1,NGTOT
5       STRFAC(KAPA,IG) = (0.0,0.0)
C
C     Note: Factor of - 2 * pi
C
      DO 10 KAPA = 1,NAT
        DX(1,KAPA) = - 2.0D0*PI*
     +               (DPROD(XKAPA(1,KAPA),B1(1)) +
     +                DPROD(XKAPA(2,KAPA),B1(2)) +
     +                DPROD(XKAPA(3,KAPA),B1(3)))
        DX(2,KAPA) = - 2.0D0*PI*
     +               (DPROD(XKAPA(1,KAPA),B2(1)) +
     +                DPROD(XKAPA(2,KAPA),B2(2)) +
     +                DPROD(XKAPA(3,KAPA),B2(3)))
        DX(3,KAPA) = - 2.0D0*PI*
     +               (DPROD(XKAPA(1,KAPA),B3(1)) +
     +                DPROD(XKAPA(2,KAPA),B3(2)) +
     +                DPROD(XKAPA(3,KAPA),B3(3)))
10      CONTINUE
C
C     FIND NUMBER OF TYPES
      NTYPES = 1
      DO 20 KAPA = 1,NAT
20      NTYPES = MAX0(NTYPES,TY(KAPA))
C
C     LOOP OVER TYPES
      DO 70 J = 1,NTYPES
C     HOW MANY ATOMS OF TYPE J
      NTYPJ = 0
      DO 30 KAPA = 1,NAT
30      IF (TY(KAPA) .EQ. J) NTYPJ = NTYPJ+1
C     LOOP OVER NTYPJ ATOMS OF TYPE J
      DO 70 N = 1,NTYPJ
C       FIND NEXT ONE OF THIS TYPE
        NEXT = 0
        DO 40 KAPA = 1,NAT
          IF (TY(KAPA) .EQ. J) THEN
            NEXT = NEXT + 1
            IF (NEXT .EQ. N) THEN
              K = KAPA
              GOTO 50
              ENDIF
            ENDIF
40        CONTINUE
50      CONTINUE
C
        IF (LREAL) THEN
          DO 60 IG = 1,NGTOT
            D1 = DBLE(IGLIST(1,IG))
            D2 = DBLE(IGLIST(2,IG))
            D3 = DBLE(IGLIST(3,IG))
            ARGMT = D1*DX(1,K) + D2*DX(2,K) + D3*DX(3,K)
            SATR = DCOS(ARGMT)
            SATI = DSIN(ARGMT)
            STRFAC(J,IG) = STRFAC(J,IG) + CMPLX(SNGL(SATR),0.0)
            SATGI(K,IG) = SATI
60          CONTINUE
        ELSE
          DO 65 IG = 1,NGTOT
            D1 = DBLE(IGLIST(1,IG))
            D2 = DBLE(IGLIST(2,IG))
            D3 = DBLE(IGLIST(3,IG))
            ARGMT = D1*DX(1,K) + D2*DX(2,K) + D3*DX(3,K)
            SATR = DCOS(ARGMT)
            SATI = DSIN(ARGMT)
            STRFAC(J,IG) = STRFAC(J,IG) + CMPLX(SNGL(SATR),SNGL(SATI))
            SATGR(K,IG) = SATR
            SATGI(K,IG) = SATI
65          CONTINUE
          ENDIF
70      CONTINUE
C
      CALL USAGE (NP1,T1)
      WRITE (IUNIT8,99) T1 - T0, NP1 - NP0
99    FORMAT(' SFACT  - structure factors ',T30,F10.3,' seconds',
     +I10,' page faults ')
C
      RETURN
      END
