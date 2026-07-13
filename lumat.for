      SUBROUTINE LUMAT (HAMR,HAMI,CWF1,DIAG,EBAR,
     +NDIM5,NDIM11,NA,NB,IFILEH,LREAL)
C---------------------------------------------------------------------
C
C     ADDS THE -C(D-E)C* TERM TO THE HAMILTONIAN
C     TO FORM THE LOWDIN U MATRIX.
C
C---------------------------------------------------------------------
C
      REAL HAMR(NDIM11),HAMI(NDIM11),DIAG(NDIM5)
      REAL REALP,AIMGG,RFAC,AFAC
      DIMENSION CWF1(NDIM5)
      LOGICAL LREAL
C.....FILES
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
C---------------------------------------------------------------
C     FOR COMPLEX HAMILTONIANS
C     COMPLEX CWF1,Z
C     REALP(Z) = REAL(Z)
C     AIMGG(Z) = AIMAG(Z)
C     LOC(I,J) = (I*I-I)/2 + J
C     IF (LREAL) WRITE(IOUT,*) '***LUMAT*** LREAL =',LREAL
C---------------------------------------------------------------
C     FOR REAL HAMILTONIANS
      REAL CWF1
      REALP(X) = X
      AIMGG(X) = 0.0*X
      LOC(I,J) = (I*I-I)/2 + J
      IF (.NOT. LREAL) WRITE(IOUT,*) '***LUMAT*** LREAL=',LREAL
C---------------------------------------------------------------
C
      CALL USAGE (N0,T0)
C
      CALL REWFIL (IFILEH)
      DO 30 K = NB,1,-1
      IF (      LREAL) CALL RDVEC (IFILEH,CWF1,NA)
      IF (.NOT. LREAL) CALL RDVEC (IFILEH,CWF1,2*NA)
      DENOM = EBAR-DIAG(K+NA)
      IF (DENOM .GT. -1.0) THEN
        WRITE(IOUT,*) '*** WARNING *** SUBROUTINE LUMAT DENOM =',DENOM
        ENDIF
      DO 30 I = NA,1,-1
        RFAC = REALP(CWF1(I))/DENOM
        AFAC = AIMGG(CWF1(I))/DENOM
        IF (LREAL) THEN
          LOCI1 = LOC(I,1)
          CALL SAXPY(I,RFAC,CWF1,1,HAMR(LOCI1),1)
        ELSE
          LOCI0 = LOC(I,0)
          DO 20 J = 1,I
            INDEX = LOCI0+J
            HAMR(INDEX) = HAMR(INDEX) + REALP(CWF1(J))*RFAC
     +                                + AIMGG(CWF1(J))*AFAC
            HAMI(INDEX) = HAMI(INDEX) - REALP(CWF1(J))*AFAC
     +                                + AIMGG(CWF1(J))*RFAC
20          CONTINUE
          ENDIF
30      CONTINUE
C
      CALL USAGE (N1,T1)
      WRITE(IUNIT8,40) T1-T0, N1-N0
40    FORMAT(' LUMAT  - U MATRIX',T30,F10.3,' SECONDS',
     +I10,' PAGE FAULTS')
C
      RETURN
      END
