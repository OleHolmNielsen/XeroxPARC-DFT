      SUBROUTINE WRVEC (IFILE,VEC,N)
      DIMENSION VEC(N)
C
C     This file keeps unit 15 (IFILEH) in memory (CRAY-2, normally).
C     It is intended that a preprocessor automatically replaces NDIM3,4
C     and sets NCMPLX=1 for real matrices, 2 for complex ones.
C
      PARAMETER (IFILEH = 15)
      PARAMETER (NDIM3 = 100, NDIM4 = 500)
      PARAMETER (NDIM5 = NDIM3+NDIM4)
      PARAMETER (MEM   = (NDIM5*NDIM5+NDIM5)*1/2+NDIM5 )
      COMMON // BUFFER(MEM)
      COMMON /INCORE/ IPT, MAXPT
C.....COMMON BLOCK FOR FILES
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
      SAVE /INCORE/
      DATA IFIRST /0/
C
      IF (IFILE .EQ. IFILEH) THEN
        IF (IFIRST .EQ. 0) THEN
          IPT = 0
          MAXPT = 0
          IFIRST = 1
          ENDIF
        IF (IPT + N + 1 .GT. MEM) THEN
C         Error: exceeded allocated memory
          WRITE (IOUT,*) 'WRVEC: exceeded memory size =', MEM
          CALL EXIT
          ENDIF
C
        IPT = IPT + 1
        BUFFER(IPT) = FLOAT(N)
        DO 100 I = 1, N
          BUFFER(IPT+I) = VEC(I)
100       CONTINUE
        IPT = IPT + N
        MAXPT = MAX( MAXPT, IPT )
C
      ELSE
        WRITE (IFILE) VEC
        ENDIF
      RETURN
      END
      SUBROUTINE RDVEC (IFILE,VEC,N)
      DIMENSION VEC(N)
C
      PARAMETER (IFILEH = 15)
      PARAMETER (NDIM3 = 100, NDIM4 = 500)
      PARAMETER (NDIM5 = NDIM3+NDIM4)
      PARAMETER (MEM   = (NDIM5*NDIM5+NDIM5)*1/2+NDIM5 )
      COMMON // BUFFER(MEM)
      COMMON /INCORE/ IPT, MAXPT
C.....COMMON BLOCK FOR FILES
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
      SAVE /INCORE/
      DATA IFIRST /0/
C
      IF (IFILE .EQ. IFILEH) THEN
        IF (IFIRST .EQ. 0) THEN
          IPT = 0
          IFIRST = 1
          ENDIF
        IF (IPT + N + 1 .GT. MAXPT) THEN
C         Error: exceeded allocated memory
          WRITE (IOUT,*) 'RDVEC: exceeded memory written =', MAXPT
          CALL EXIT
          ENDIF
C
        IPT = IPT + 1
        NBUF = NINT(BUFFER(IPT))
        IF (ABS(BUFFER(IPT) - FLOAT(NBUF)) .GT. 1.0E-6) THEN
          WRITE (IOUT,*) 'RDVEC: Something wrong with NBUF=',BUFFER(IPT)
          CALL EXIT
          ENDIF
        IF (N .GT. NBUF) THEN
          WRITE (IOUT,*) 'RDVEC: N=',N,'>NBUF=',NBUF
          CALL EXIT
          ENDIF
        DO 100 I = 1, N
          VEC(I) = BUFFER(IPT+I)
100       CONTINUE
        IPT = IPT + NBUF
      ELSE
C
        READ (IFILE) VEC
        ENDIF
      RETURN
      END
      SUBROUTINE REWFIL (IFILE)
      PARAMETER (IFILEH = 15)
      COMMON /INCORE/ IPT, MAXPT
C
      IF (IFILE .EQ. IFILEH) THEN
        IPT = 0
      ELSE
        REWIND (UNIT = IFILE)
        ENDIF
      RETURN
      END
      SUBROUTINE CONVEC (N,CVEC)
C
C     CONJUGATE THE FIRST N COMPONENTS OF THE VECTOR CVEC
C
      COMPLEX CVEC(N)
C
      IF (N .LE. 0) RETURN
      M = MOD(N,5)
      IF (M .EQ. 0) GOTO 20
      DO 10 I = 1,M
10      CVEC(I) = CONJG(CVEC(I))
      IF (N .LT. 5) RETURN
20    DO 30 I = M+1,N,5
        CVEC(I)   = CONJG(CVEC(I))
        CVEC(I+1) = CONJG(CVEC(I+1))
        CVEC(I+2) = CONJG(CVEC(I+2))
        CVEC(I+3) = CONJG(CVEC(I+3))
        CVEC(I+4) = CONJG(CVEC(I+4))
30      CONTINUE
C
      RETURN
      END
      SUBROUTINE DZERO (N,DARRAY)
C
C     ZERO THE FIRST N ELEMENTS OF THE DOUBLE PRECISION ARRAY DARRAY
C
      DOUBLE PRECISION DARRAY(N)
C
      DO 10 I = 1,N
        DARRAY(I) = 0.0D0
10      CONTINUE
C
      RETURN
      END
      SUBROUTINE CZERO (N,CARRAY)
C
C     ZERO THE FIRST N ELEMENTS OF THE COMPLEX ARRAY CARRAY
C
      COMPLEX CARRAY(N)
C
      DO 10 I = 1,N
        CARRAY(I) = (0.0,0.0)
10      CONTINUE
C
      RETURN
      END
      SUBROUTINE SZERO (N,SARRAY)
C
C     ZERO THE FIRST N ELEMENTS OF THE SCALAR ARRAY SARRAY
C
      REAL SARRAY(N)
C
      DO 10 I = 1,N
        SARRAY(I) = 0.0
10      CONTINUE
C
      RETURN
      END
