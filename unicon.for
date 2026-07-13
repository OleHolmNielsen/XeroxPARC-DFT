      SUBROUTINE UNICON (ENONLO,SWF2,FORCE,B1,B2,B3,ULA,NDIM9,NAT,
     +LSTRES,LFORCE)
C
C     CONVERSION FROM RYDBERGS TO EV AND TRANSFORMATION OF FORCES
C     TO CARTESIAN BASIS
C
      DIMENSION B1(3),B2(3),B3(3)
      DOUBLE PRECISION ENONLO,SWF2(6),FORCE(3,NDIM9),GVEC(3),FACTOR
      LOGICAL LSTRES,LFORCE
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
C.....PHYSICAL AND MATHEMATICAL CONSTANTS
      DOUBLE PRECISION ABOHR,RYEV,RYDERG,PI,SPI
      COMMON /CONST/   ABOHR,RYEV,RYDERG,PI,SPI
C----------------------------------------------------------------
C
      ENONLO = ENONLO*RYEV
C
      IF (LSTRES) THEN
        DO 10 I = 1,6
          SWF2(I) = SWF2(I)*RYEV
10        CONTINUE
        ENDIF
C
      IF (LFORCE) THEN
        FACTOR = 2.0D0*PI/DBLE(ULA)
        DO 20 KAPA = 1,NAT
          GVEC(1) = FORCE(1,KAPA)*FACTOR
          GVEC(2) = FORCE(2,KAPA)*FACTOR
          GVEC(3) = FORCE(3,KAPA)*FACTOR
C         TRANSFORMATION INTO CARTESIAN BASIS:
          FORCE(1,KAPA) = GVEC(1)*DBLE(B1(1)) +
     +                    GVEC(2)*DBLE(B2(1)) +
     +                    GVEC(3)*DBLE(B3(1))
          FORCE(2,KAPA) = GVEC(1)*DBLE(B1(2)) +
     +                    GVEC(2)*DBLE(B2(2)) +
     +                    GVEC(3)*DBLE(B3(2))
          FORCE(3,KAPA) = GVEC(1)*DBLE(B1(3)) +
     +                    GVEC(2)*DBLE(B2(3)) +
     +                    GVEC(3)*DBLE(B3(3))
20        CONTINUE
        ENDIF
C
      RETURN
      END
