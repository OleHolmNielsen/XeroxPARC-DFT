      SUBROUTINE KINETI (GK,NANB,KE,CEV,SWF2,WOCC,NBDS,
     +  CONV1,NDIM5,NDIM8,LREAL)
C
C     THE KINETIC ENERGY CONTRIBUTIONS TO TOTAL ENERGY AND STRESS
C
C-----------------------------------------------------------------------
      LOGICAL LREAL
      DIMENSION GK(4,NDIM5),CEV(NDIM5,NDIM8)
      DOUBLE PRECISION CONV1,KE,SUM,SWF2(6),WOCC(NDIM8)
C
C.....FILES
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
C.....PHYSICAL AND MATHEMATICAL CONSTANTS
      DOUBLE PRECISION ABOHR,RYEV,RYDERG,PI,SPI
      COMMON /CONST/   ABOHR,RYEV,RYDERG,PI,SPI
C-----------------------------------------------------------------------
C     FOR COMPLEX HAMILTONIANS
C     COMPLEX CEV,Z
C     DOUBLE PRECISION DC2ABS
C     DC2ABS(Z) = DPROD(REAL(Z),REAL(Z)) + DPROD(AIMAG(Z),AIMAG(Z))
C     IF (LREAL) WRITE(IOUT,*) '***KINETI*** LREAL =',LREAL
C-----------------------------------------------------------------------
C     FOR REAL HAMILTONIANS
      REAL CEV
      DOUBLE PRECISION DC2ABS
      DC2ABS(X) = DPROD(X,X)
      IF (.NOT. LREAL) WRITE(IOUT,*) '***KINETI*** LREAL =',LREAL
C-----------------------------------------------------------------------
C
      SUM = 0.0D0
      DO 10 K = 1,NBDS
      DO 10 M = 1,NANB
10      SUM = SUM + WOCC(K)*DPROD(GK(4,M),GK(4,M))*DC2ABS(CEV(M,K))
      KE = SUM*CONV1*RYEV
C
      INDEX = 0
      DO 30 I = 1,3
      DO 30 J = I,3
        SUM = 0.0D0
        DO 20 K = 1,NBDS
        DO 20 M = 1,NANB
20        SUM = SUM + WOCC(K)*DPROD(GK(I,M),GK(J,M))*DC2ABS(CEV(M,K))
        INDEX = INDEX + 1
30      SWF2(INDEX) = 2.0D0*SUM*CONV1*RYEV
C
      RETURN
      END
