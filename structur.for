      SUBROUTINE STRUCT(TEXT,A01,A02,A03,A1,A2,A3,NAT,
     +  TY,XKAPA,NATOM,ULA,STRAIN,NDIM9,NTYPMX)
C
C     Read lattice structure from file
C
C     FILE FORMAT:
C     Text header
C     Number of atoms
C     Translation vectors A1,2,3
C     Atomic number, (x,y,z)-position (all atoms)
C     lattice constant
C     Strain (e1-e6)
C     Atomic displacements (ux,uy,uz) (all atoms)
C     (The last two informations may be omitted entirely)
C     (The set of displacements may be repeated for additional u's)
C
      INTEGER TY(NDIM9), NATOM(NDIM9)
      CHARACTER*80 TEXT
      REAL A01(3),A02(3),A03(3),A1(3),A2(3),A3(3),
     +  XKAPA(3,NDIM9),STRAIN(6)
      REAL DISPLC(3)
      CHARACTER*2 NAME
C
C.....FILES
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
C.....PHYSICAL AND MATHEMATICAL CONSTANTS
      DOUBLE PRECISION ABOHR,RYEV,RYDERG,PI,SPI
      COMMON /CONST/   ABOHR,RYEV,RYDERG,PI,SPI
C
C-----------------------------------------------------------------------
      DO 100 I = 1,6
        STRAIN(I) = 0.0
100     CONTINUE
      WRITE (IOUT,*) 'Reading crystal structure from unit ', ISTRUC
C     Read heading card (any comment):
      READ (ISTRUC,'(A)',END=520) TEXT
      WRITE (IOUT,*) 'Structure file input:'
      WRITE (IOUT,'(1X,A)') TEXT
C     Read number of atoms per unit cell:
      READ (ISTRUC,*,END=520) NAT
      IF (NAT .GT. NDIM9 .OR. NAT .LE. 0) THEN
        WRITE (IOUT,*) 'Structur: ERROR *** illegal NAT, NDIM9 =',
     +    NAT, NDIM9
        CALL EXIT
        ENDIF
C     Read the elementary translations (unstrained crystal):
      READ  (ISTRUC,*,END=520) A01, A02, A03
      WRITE (IUNIT7,150) A01,A02,A03
150   FORMAT ('0a1',3F10.5/' a2',3F10.5/' a3',3F10.5)
      DO 170 I = 1,3
        A1(I) = A01(I)
        A2(I) = A02(I)
        A3(I) = A03(I)
170     CONTINUE
C
      WRITE (IUNIT7,160)
160   FORMAT (/'0    K   Type',14X,'X(K)',T45,'Atom no.    Name')
      ITYPE = 0
      DO 180 I = 1,NAT
C       Read atomic number and cartesian coordinates of the basis atoms
        READ (ISTRUC,*,END=520) NATOM(I),(XKAPA(J,I),J=1,3)
C       Assign an atomic type (for internal purposes)
        DO 176 J = 1, (I - 1)
          IF (NATOM(J) .EQ. NATOM(I)) THEN
C           Type was located previously:
            TY(I) = TY(J)
            GOTO 178
            ENDIF
176       CONTINUE
C       This is a new type:
        ITYPE = ITYPE + 1
        IF (ITYPE .GT. NTYPMX) GOTO 500
        TY(I) = ITYPE
C
178     CALL ELTNAM (NAME,NATOM(I))
        WRITE (IUNIT7,179) I,TY(I),(XKAPA(J,I),J=1,3),NATOM(I),NAME
179     FORMAT (1X,2I5,3F10.5,I10,5X,A)
180     CONTINUE
C
C     Lattice constant
      READ (ISTRUC,*,END=520) ULA
      WRITE (IUNIT7,190) ULA, ULA/ABOHR
190   FORMAT (/' Unit of length, ULA =',F10.5,' Angstroms =',F10.5,
     +        ' a.u.')
      IF (ULA .LT. 0.5) THEN
        WRITE (IOUT,*) 'STRUCT *** ERROR *** ULA =',ULA
        CALL EXIT
      ELSE IF (ULA .LT. 1.0 .OR. ULA .GT. 20.0) THEN
        WRITE (IOUT,*) '*** WARNING *** is ULA ok ?'
        ENDIF
C
C     Strain tensor (if any)
      READ (ISTRUC,*,END=400,ERR=400) STRAIN
      TOTSTR = 0.0
      DO 200 I = 1,6
        TOTSTR = TOTSTR + ABS(STRAIN(I))
200     CONTINUE
      IF (TOTSTR .LT. 1.0E-20) GOTO 300
      WRITE (IUNIT7,210) STRAIN
210   FORMAT('0Strain tensor elements (indexed 1-6)'/1X,6F13.5)
C     Straining the crystal:
      CALL STRNIG(A1,STRAIN)
      CALL STRNIG(A2,STRAIN)
      CALL STRNIG(A3,STRAIN)
      DO 220 I = 1,NAT
        CALL STRNIG(XKAPA(1,I),STRAIN)
220     CONTINUE
      WRITE (IUNIT7,*) ' Strained translation vectors:'
      WRITE (IUNIT7,150) A1,A2,A3
C
C     Displacing atoms in (strained) unit cell:
300   TOTDIS = 0.0
      DO 310 I = 1,NAT
        READ (ISTRUC,*,END=400,ERR=400) DISPLC
        DIS = 0.0
        DO 305 J = 1,3
          IF (ABS(DISPLC(J)) .GT. 1.0E-20) THEN
            XKAPA(J,I) = XKAPA(J,I) + DISPLC(J)
            DIS = DIS + ABS(DISPLC(J))
            ENDIF
305       CONTINUE
          IF (DIS .GT. 1.0E-20) WRITE (IUNIT7,307) I,DISPLC
307       FORMAT(' Atom ',I3,' DX=',3G15.5)
          TOTDIS = TOTDIS + DIS
310       CONTINUE
C     IF (TOTDIS .GT. 1.0E-20)
C    +    WRITE(IUNIT7,*) 'Atoms have been displaced'
C     If no strain nor displacement, continue (maybe file is ended now)
      IF (TOTSTR+TOTDIS .LT. 1.0E-20) GOTO 300
C
      WRITE (IUNIT7,160)
      DO 320 I = 1,NAT
        CALL ELTNAM (NAME,NATOM(I))
        WRITE (IUNIT7,179) I,TY(I),(XKAPA(J,I),J=1,3),NATOM(I),NAME
320     CONTINUE
C
C     Continue reading file for possible additional displacements:
      GOTO 300
C
400   RETURN
C-----------------------------------------------------------------------
500   WRITE (IOUT,510) NTYPMX,(NATOM(I),I=1,NAT)
510   FORMAT('0Subroutine STRUCT *** FATAL ERROR***'/
     +' Number of atomic types exceeds dimension (NTYPMX=',I4,')'/
     +' The array NATOM is:'/9(1X,10I7/) )
      CALL EXIT
      RETURN
C
520   WRITE (IOUT,530)
530   FORMAT('0Structure file *** premature end-of-file *** ABORT ***')
      CALL EXIT
      RETURN
      END
      SUBROUTINE STRNIG(VEC,STRAIN)
C
C     Strain VEC according to STRAIN
C
      REAL VEC(3),STRAIN(6)
C     Indexing strain tensor (see Nye, Physical properties of cryst,p134
      INTEGER INDEX(3,3)
      REAL FACTOR(3,3),TEMP(3)
      DATA ((INDEX(I,J),I=1,3),J=1,3)  /1,  6, 5,  6,2,  4,  5, 4,3/
      DATA ((FACTOR(I,J),I=1,3),J=1,3) /1.,.5,.5, .5,1.,.5, .5,.5,1./
C-----------------------------------------------------------------------
      DO 100 I  = 1,3
        TEMP(I) = VEC(I)
100     CONTINUE
      DO 110   I = 1,3
        DO 110 J = 1,3
          VEC(I) = VEC(I) + FACTOR(I,J)*STRAIN(INDEX(I,J))*TEMP(J)
110       CONTINUE
      RETURN
      END
