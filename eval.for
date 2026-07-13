      SUBROUTINE EVAL (CEV,EE1,EE2,WK,CWF1,RWK,IFILEH,
     +NDIM5,NDIM8,NA,NB,NBDS,IEVAL,IWRCUA,HRSUB,HISUB,
     +SRSUB,SISUB,VRSUB,VISUB,ALFR,ALFI,BETA,IWK,LREAL)
C
C     EVALUATE THE EIGENVALUES AFTER A LOWDIN PERTURBATION
C     CALCULATION OF THE EIGENVECTORS STORED IN CEV
C     ON ENTRY THE EIGENVECTORS M U S T BE NORMALIZED.
C
      INTEGER IWK(NDIM8)
      REAL EE1(NDIM8),EE2(NDIM8),RWK(NDIM8)
      DIMENSION CEV(NDIM5,NDIM8),WK(NDIM5,NDIM8),CWF1(NDIM5)
      DOUBLE PRECISION HRSUB(NDIM8,NDIM8),HISUB(NDIM8,NDIM8)
      DOUBLE PRECISION SRSUB(NDIM8,NDIM8),SISUB(NDIM8,NDIM8)
      DOUBLE PRECISION VRSUB(NDIM8,NDIM8),VISUB(NDIM8,NDIM8)
      DOUBLE PRECISION ALFR(NDIM8),ALFI(NDIM8),BETA(NDIM8),EPS
      LOGICAL LREAL,LMATV
      PARAMETER (LMATV = .TRUE.)
      PARAMETER (EPS = -1.0D0)
C
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
      REAL SDOT,SNRM2,SCNRM2
      COMPLEX CDOTC,DUM
      DOUBLE PRECISION DSDOT
C----------------------------------------------------------------------
C     FOR REAL HAMILTONIANS
      REAL CEV,WK,CWF1
      REALP(X) = X
      AIMGG(X) = 0.0*X
      ZERO = 0.0
      IF (.NOT. LREAL) WRITE(IOUT,*) '***EVAL*** LREAL =',LREAL
C----------------------------------------------------------------------
C     FOR COMPLEX HAMILTONIANS
C     COMPLEX CEV,WK,CWF1,Z,ZERO
C     REAL REALP,AIMGG
C     REALP(Z) = REAL(Z)
C     AIMGG(Z) = AIMAG(Z)
C     ZERO = CMPLX(0.0,0.0)
C     IF (LREAL) WRITE(IOUT,*) '***EVAL*** LREAL =',LREAL
C----------------------------------------------------------------------
C
      CALL USAGE (N0,T0)
C
      IF (IWRCUA .EQ. 3) N=1
      IF (IWRCUA .EQ. 2) N=NA+1
      IF (IWRCUA .LE. 1 .OR. IWRCUA .GT. 3) THEN
        WRITE(IOUT,*) '**ERROR** IN EVAL IWRCUA =',IWRCUA
        CALL EXIT
        ENDIF
C
      NANB = NA+NB
      DO 10 L = 1,NBDS
      DO 10 I = 1,NANB
10      WK(I,L) = ZERO
C
      CALL REWFIL (IFILEH)
      DO 20 I = NANB,N,-1
      IF (      LREAL) CALL RDVEC (IFILEH,CWF1,I)
      IF (.NOT. LREAL) CALL RDVEC (IFILEH,CWF1,2*I)
      IF (.NOT. LREAL) CALL CONVEC (I,CWF1)
      DO 20 L = 1,NBDS
        IF (LREAL) THEN
          WK(I,L) = WK(I,L) + SDOT (I,CWF1,1,CEV(1,L),1)
          CALL SAXPY(I-1,CEV(I,L),CWF1,1,WK(1,L),1)
        ELSE
          WK(I,L) = WK(I,L) + CDOTC(I,CWF1,1,CEV(1,L),1)
          CALL CAXPY(I-1,CEV(I,L),CWF1,1,WK(1,L),1)
          ENDIF
20      CONTINUE
C
      IF (IEVAL .EQ. 2 .AND. LREAL) THEN
C
        DO 30 L = 1,NBDS
        DO 30 M = 1,L
          HRSUB(M,L) = DSDOT(NANB,CEV(1,M),1, WK(1,L),1)
30        SRSUB(M,L) = DSDOT(NANB,CEV(1,M),1,CEV(1,L),1)
C
        IF (IWRCUA .EQ. 2) THEN
          DO 40 L = 1,NBDS
40          HRSUB(L,L) = HRSUB(L,L) + DPROD(RWK(L),RWK(L))*DBLE(EE1(L))
          ENDIF
C
        DO 50 L = 1,NBDS
        DO 50 M = 1,L-1
          HRSUB(L,M) = HRSUB(M,L)
50        SRSUB(L,M) = SRSUB(M,L)
C
        IFAIL = 0
        CALL F02BJF (NBDS,HRSUB,NDIM8,SRSUB,NDIM8,EPS,ALFR,ALFI,
     +  BETA,LMATV,VRSUB,NDIM8,IWK,IFAIL)
C
        IF (IFAIL .NE. 0) THEN
          WRITE(IOUT,*) '***ERROR*** F02BJF IFAIL = ',IFAIL
          CALL EXIT
          ENDIF
C
        DO 60 L = 1,NBDS
60        EE2(L) = SNGL(ALFR(L)/BETA(L))
C
        DO 80 L = 1,NBDS
          IF (ABS(ALFI(L)/ALFR(L)) .GT. 1.0E-4) THEN
          WRITE(IOUT,*) 'IMAGINARY PART OF EIGENVALUE TOO LARGE'
          WRITE(IOUT,*) 'REAL PART ',ALFR(L),' IMAG PART ',ALFI(L)
          ENDIF
80        CONTINUE
C
        DO 90 L = 1,NBDS
        DO 90 I = 1,NANB
90        WK(I,L) = ZERO
C
        DO 100 L = 1,NBDS
        DO 100 M = 1,NBDS
          VR = VRSUB(M,L)
100       CALL SAXPY(NANB,VR,CEV(1,M),1,WK(1,L),1)
C
        DO 110 L = 1,NBDS
          SNORM = SNRM2(NANB,WK(1,L),1)
          SNORM = 1.0/SNORM
          DO 110 I = 1,NANB
110         CEV(I,L) = WK(I,L)*SNORM
C
      ELSE IF (IEVAL .EQ. 2) THEN
C
        DO 120 L = 1,NBDS
        DO 120 M = 1,L
          HRSUB(M,L) = 0.0D0
          HISUB(M,L) = 0.0D0
          SRSUB(M,L) = 0.0D0
120       SISUB(M,L) = 0.0D0
C
        IF (IWRCUA .EQ. 2) THEN
          DO 130 L = 1,NBDS
130         HRSUB(L,L) = DPROD(RWK(L),RWK(L))*DBLE(EE1(L))
          ENDIF
C
        DO 140 L = 1,NBDS
        DO 140 M = 1,L
        DO 140 I = 1,NANB
          RCIM = REALP(CEV(I,M))
          ACIM = AIMGG(CEV(I,M))
          RCIL = REALP(CEV(I,L))
          ACIL = AIMGG(CEV(I,L))
          RWIL = REALP( WK(I,L))
          AWIL = AIMGG( WK(I,L))
          HRSUB(M,L) = HRSUB(M,L) + DPROD(RCIM,RWIL) + DPROD(ACIM,AWIL)
          HISUB(M,L) = HISUB(M,L) + DPROD(RCIM,AWIL) - DPROD(ACIM,RWIL)
          SRSUB(M,L) = SRSUB(M,L) + DPROD(RCIM,RCIL) + DPROD(ACIM,ACIL)
          SISUB(M,L) = SISUB(M,L) + DPROD(RCIM,ACIL) - DPROD(ACIM,RCIL)
140       CONTINUE
C
      DO 150 L = 1,NBDS
      DO 150 M = 1,L-1
        HRSUB(L,M) =  HRSUB(M,L)
        HISUB(L,M) = -HISUB(M,L)
        SRSUB(L,M) =  SRSUB(M,L)
150     SISUB(L,M) = -SISUB(M,L)
C
        IFAIL = 0
        CALL F02GJF (NBDS,HRSUB,NDIM8,HISUB,NDIM8,SRSUB,NDIM8,
     +  SISUB,NDIM8,EPS,ALFR,ALFI,BETA,LMATV,VRSUB,NDIM8,
     +  VISUB,NDIM8,IWK,IFAIL)
C
        IF (IFAIL .NE. 0) THEN
          WRITE(IOUT,*) '***ERROR*** F02GJF IFAIL = ',IFAIL
          CALL EXIT
          ENDIF
C
        DO 160 L = 1,NBDS
          IF (ABS(ALFI(L)/ALFR(L)) .GT. 1.0E-4) THEN
          WRITE(IOUT,*) 'IMAGINARY PART OF EIGENVALUE TOO LARGE'
          WRITE(IOUT,*) 'REAL PART ',ALFR(L),' IMAG PART ',ALFI(L)
          ENDIF
160       CONTINUE
C
        DO 170 L = 1,NBDS
170       EE2(L) = SNGL(ALFR(L)/BETA(L))
C
        DO 180 L = 1,NBDS
        DO 180 I = 1,NANB
180       WK(I,L) = ZERO
C
        DO 190 L = 1,NBDS
        DO 190 M = 1,NBDS
          VR = SNGL(VRSUB(M,L))
          VI = SNGL(VISUB(M,L))
          DUM = CMPLX(VR,VI)
190       CALL CAXPY(NANB,DUM,CEV(1,M),1,WK(1,L),1)
C
        DO 200 L = 1,NBDS
          SNORM = SCNRM2(NANB,WK(1,L),1)
          SNORM = 1.0/SNORM
          DO 200 I = 1,NANB
200         CEV(I,L) = WK(I,L)*SNORM
C
        ENDIF
C
      IF (IEVAL .EQ. 1) THEN
C
        IF (IWRCUA .EQ. 2) THEN
          DO 210 L = 1,NBDS
210         EE2(L) = RWK(L)*RWK(L)*EE1(L)
        ELSE IF (IWRCUA .EQ. 3) THEN
          DO 220 L = 1,NBDS
220         EE2(L) = 0.0
          ENDIF
C
        DO 230 L = 1,NBDS
          IF (LREAL) THEN
            EE2(L) = EE2(L) + SDOT(NANB,CEV(1,L),1,WK(1,L),1)
          ELSE
            DUM = CDOTC(NANB,CEV(1,L),1,WK(1,L),1)
            EE2(L) = EE2(L) + REAL(DUM)
            ENDIF
230       CONTINUE
C
        ENDIF
C
      CALL USAGE (N1,T1)
      WRITE(IUNIT8,99) T1-T0, N1-N0
99    FORMAT (' EVAL   - EIGENVALUES ',T30,F10.3,' SECONDS',
     +I10,' PAGE FAULTS')
C
      RETURN
      END
     
