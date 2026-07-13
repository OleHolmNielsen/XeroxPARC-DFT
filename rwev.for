      SUBROUTINE RWEV (IUNIT,IREC,IRWEV,EE1,EE2,NBDS,NANB,
     +                 NDIM5,NDIM8,CEV,LREAL)
C
C     READING/WRITING EIGENVECTORS AND EIGENVALUES FROM/TO FILE
C
C     IUNIT ...... UNIT NUMBER FOR READING AND WRITING
C     IREC ....... RECORD NUMBER INDICATING READ/WRITE BEGINS
C                  THE SUBSEQUENT NDIM8 RECORDS WILL BE READ
C     IRWEV ...... -1 READ FROM FILE.  +1 WRITE TO FILE
C     EE1 ........ FIRST ORDER EIGENVALUES
C     EE2 ........ HIGHER ORDER EIGENVALUES
C     NBDS ....... NUMBER OF BANDS
C     NANB ....... NUMBER OF WAVES IN HAMILTONIAN
C     NDIM5 ...... MAXIMUM DIMENSION OF THE EIGENVECTORS
C     NDIM8 ...... MAXIMUM NUMBER OF BANDS
C     CEV ........ COMPLEX OR REAL ARRAY CONTAINING THE EIGENVECTORS
C
      REAL EE1(NDIM8),EE2(NDIM8)
C.....FILES
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
C.....MACHINE TYPE
      COMMON /MACH1/ MACTYP
      COMMON /MACH2/ MACHIN(0:20)
      CHARACTER*7 MACHIN
C
      LOGICAL INIT,LREAL
C----------------------------------------------------------------------
C     FOR COMPLEX HAMILTONIANS
C     COMPLEX CEV(NDIM5,NDIM8)
C     SAVE INIT, MAXRL
C     PARAMETER (ICRL = 2)
C     IF (LREAL) WRITE(IOUT,*) '***RWEV*** LREAL =',LREAL
C----------------------------------------------------------------------
C     FOR REAL HAMILTONIANS
      REAL CEV(NDIM5,NDIM8)
      SAVE INIT, MAXRL
      PARAMETER (ICRL = 1)
      IF (.NOT. LREAL) WRITE(IOUT,*) '***RWEV*** LREAL =',LREAL
C----------------------------------------------------------------------
      DATA INIT /.FALSE./
C
      IF (.NOT. INIT) THEN
C                MAXIMUM LENGTH OF A RECORD
        MAXRL = 3 + ICRL * NDIM5
C
        IF (MACTYP .EQ. 1) THEN
C           IBM: RECORD-LENGTH = 4 * NUMBER OF WORDS
          MAXRL = MAXRL * 4
C         WRITE (IUNIT7,*) 'RWEV   - IBM - MAXRL*4 =',MAXRL
        ELSE IF (MACTYP .EQ. 2 .OR. MACTYP .EQ. 7) THEN
C           CRAY: RECORD-LENGTH = 8 * NUMBER OF WORDS
          MAXRL = MAXRL * 8
C         WRITE (IUNIT7,*) 'RWEV   - CRAY - MAXRL*8 =',MAXRL
        ELSE IF (MACTYP .EQ. 6) THEN
C           FPS: RECORD-LENGTH = 8 * NUMBER OF WORDS
          MAXRL = MAXRL * 8
C         WRITE (IUNIT7,*) 'RWEV   - FPS - MAXRL*8 =',MAXRL
        ELSE IF (MACTYP .EQ. 10) THEN
C           UNIX: RECORD-LENGTH = 4 * NUMBER OF WORDS
          MAXRL = MAXRL * 4
C         WRITE (IUNIT7,*) 'RWEV   - UNIX - MAXRL*4 =',MAXRL
        ELSE IF (MACTYP .EQ. 11) THEN
C           Apollo: RECORD-LENGTH = 4 * NUMBER OF WORDS
          MAXRL = MAXRL * 4
          WRITE (IUNIT7,*) 'RWEV   - Apollo - MAXRL*4 =',MAXRL
          WRITE (IUNIT7,*) 'RWEV   - Apollo - IUNIT =', IUNIT
          ENDIF
C
C     Temporary solution to IBM problem: Do not name file.
        IF (MACTYP .EQ. 1) THEN
          OPEN (UNIT=IUNIT , RECL=MAXRL , ACCESS='DIRECT' ,
     +          FORM='UNFORMATTED' , STATUS='UNKNOWN')
        ELSE
          OPEN (UNIT=IUNIT , RECL=MAXRL , ACCESS='DIRECT' ,
     +          FORM='UNFORMATTED' , FILE='EVFILE' , STATUS='UNKNOWN')
          ENDIF
C
        INIT=.TRUE.
        ENDIF
C
      IF (IRWEV .EQ. 1) THEN
C                WRITE EIGENVECTORS
C
        DO 10 I = 1,NDIM8
          IREC  = IREC + 1
          WRITE (IUNIT,REC=IREC) NBDS,EE1(I),EE2(I),
     +                          (CEV(J,I),J=1,NANB)
10        CONTINUE
C
C
      ELSE IF (IRWEV .EQ. -1) THEN
C                READ EIGENVECTORS
C
          DO 20 I = 1,NDIM8
            IREC  = IREC + 1
            READ (IUNIT,REC=IREC) NBDS,EE1(I),EE2(I),
     +                          (CEV(J,I),J=1,NANB)
20          CONTINUE
C
      ELSE
C            ERROR
C
          WRITE (IOUT,30) IRWEV
30        FORMAT(' SUBROUTINE RWEV *** ERROR *** IRWEV =',I10)
          CALL EXIT
C
          ENDIF
C
      RETURN
      END
