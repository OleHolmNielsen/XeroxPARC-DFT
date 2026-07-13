      SUBROUTINE CHKINV(LREAL,TY,XKAPA,NAT)
C
C     Check whether the lattice contains inversion around the origin.
C     If so, the Hamiltonian may be taken as real.
C     The routine prints a warning if this is inconsistent with LREAL
C
C.....FILES
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
C
      INTEGER TY(NAT)
      REAL XKAPA(3,NAT)
      LOGICAL LREAL, LINV
      PARAMETER (ACCUR = 1.0E-6)
C-----------------------------------------------------------------------
C
      LINV = .TRUE.
      DO 200 KAPA = 1, NAT
C       Locate X in the XKAPA
        DO 110 K = 1, NAT
          IF (TY(K) .NE. TY(KAPA)) GOTO 110
          DIFF = 0.0
          DO 100 I = 1, 3
            DIFF = DIFF + ABS(XKAPA(I,K) + XKAPA(I,KAPA))
100         CONTINUE
C         OK, this atom has an inverse:
          IF (DIFF .LT. ACCUR) GOTO 200
110       CONTINUE
C       The atom KAPA has no inverse
        LINV = .FALSE.
        IF (LREAL) THEN
          WRITE (IOUT,120) (XKAPA(I,KAPA), I = 1, 3)
120       FORMAT (1X,'CHKINV *** ',
     +    'No inverse atom for XKAPA = ', 3F12.6)
          ENDIF
200     CONTINUE
C
      IF (LREAL .AND. (.NOT. LINV)) THEN
        WRITE (IOUT,220)
220     FORMAT (1X,'CHKINV *** ',
     +  'ERROR: Real Hamiltonian but no inversion around origin')
        CALL EXIT
C
      ELSE IF ((.NOT. LREAL) .AND. LINV) THEN
        WRITE (IOUT,230)
230     FORMAT(' CHKINV *** WARNING: Complex Hamiltonian not necessary'/
     +  ' when the lattice has inversion around the origin')
        IWARN = 1
        ENDIF
      RETURN
      END
