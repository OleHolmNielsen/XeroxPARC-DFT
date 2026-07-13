      SUBROUTINE PUTNL (CWF1,VNLL,HAMR,HAMI,DIAG,IFILEH,IFILEV,
     +LRWH,WRITEV,NDIM5,NDIM11,NROW,I1,I2,M,NA,NB,IWRCUA,LREAL)
C
C     WRITE THE NROW MATRIX ELEMENTS FROM I1 TO I2
C     OF THE M TH ROW OF THE HAMILTONIAN TO FILE
C
      LOGICAL WRITEV,LREAL,LRWH
      REAL HAMR(NDIM11),HAMI(NDIM11),DIAG(NDIM5)
      DIMENSION CWF1(NDIM5),VNLL(NDIM5)
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
C--------------------------------------------------------------------
C     FOR COMPLEX HAMILTONIANS
C     COMPLEX CWF1,VNLL,Z
C     REAL REALP,AIMGG
C     REALP(Z) = REAL(Z)
C     AIMGG(Z) = AIMAG(Z)
C     LOC(I,J) = (I*I-I)/2 + J
C     IF (LREAL) WRITE(IOUT,*) '***PUTNL*** LREAL =',LREAL
C---------------------------------------------------------------------
C     FOR REAL HAMILTONIANS
      REAL CWF1,VNLL
      REALP(X) = X
      AIMGG(X) = 0.0*X
      LOC(I,J) = (I*I-I)/2 + J
      IF (.NOT. LREAL) WRITE(IOUT,*) '***PUTNL*** LREAL =',LREAL
C--------------------------------------------------------------------
C
      IF (LREAL) THEN
C
        IF (M .LE. NA) THEN
          LROW = LOC(M,0)
          DO 10 N = I1,I2
            LOCNM = LROW+N
            HAMR(LOCNM) = REALP(CWF1(N))
10          CONTINUE
          IF (IWRCUA .EQ. 3) CALL WRVEC (IFILEH,CWF1,NROW)
        ELSE IF (M .GT. NA) THEN
          IF (IWRCUA .GT. 0) CALL WRVEC (IFILEH,CWF1,NROW)
          ENDIF
C
        IF (M .GT. NA .AND. M .EQ. I2) DIAG(I2) = REALP(CWF1(I2))
C
        IF (WRITEV) CALL WRVEC (IFILEV,VNLL,NROW)
C
      ELSE
C
        IF (M .LE. NA) THEN
          LROW = LOC(M,0)
          DO 20 N = I1,I2
            LOCNM = LROW+N
            HAMR(LOCNM) = REALP(CWF1(N))
            HAMI(LOCNM) = AIMGG(CWF1(N))
20          CONTINUE
          LOCMM = LOC(M,M)
          HAMI(LOCMM) = 0.0
          IF (IWRCUA .EQ. 3) CALL WRVEC (IFILEH,CWF1,2*NROW)
        ELSE IF (M .GT. NA) THEN
          IF (IWRCUA .GT. 0) CALL WRVEC (IFILEH,CWF1,2*NROW)
          ENDIF
C
        IF (M .GT. NA .AND. M .EQ. I2) DIAG(I2) = REALP(CWF1(I2))
C
        IF (WRITEV) CALL WRVEC (IFILEV,VNLL,2*NROW)
C
        ENDIF
C
      RETURN
      END
