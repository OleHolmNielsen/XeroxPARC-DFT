      SUBROUTINE POTGET (
     +  NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8,NDIM9,NDIM10,
     +  NDIM12,NDIM20,NLMAX,NPARAM,NTYPMX,NSPIN,NRUN,RUNPAR,
     +  ULA,VOLUM,A1,A2,A3,B1,B2,B3,NG1,NG2,NG3,EPSILO,
     +  STRAIN,XKAPA,G2TAB,TY,NTAB,IGLIST,NATOM,TEXT,
     +  POTPAR,ZKAPA,POTTYP,POTXC,POTFIL,VNLIJ,ZTOT,NEL,G2MAX,
     +  VNLL,POTWRK,VIONIC,VGTOT,ROTOT,NMIX,HKMIX,XKMIX,FKMIX,
     +  NAT,NGTOT,MINDIM,IEXCOR,EXFACT,NONLOC,LCHARG)
C
C     Read potentials from file & check their consistency
C
C.....Lattice structure
      REAL XKAPA(3,NDIM9),A1(3),A2(3),A3(3),B1(3),B2(3),B3(3),
     +     G2TAB(NDIM2),STRAIN(6)
      INTEGER TY(NDIM9),NTAB(NDIM2), IGLIST(3,NDIM1)
C.....Atomic potentials
      DOUBLE PRECISION POTPAR(NTYPMX,NPARAM)
      REAL ZKAPA(NDIM9)
      INTEGER POTTYP(NTYPMX), POTXC(NTYPMX)
      CHARACTER*80 POTFIL(NTYPMX)
      LOGICAL NONLOC, LCHARG
      COMPLEX POTWRK(*)
C.....Potentials and densities
      COMPLEX VIONIC(NDIM1,NSPIN), VGTOT(NDIM1,NSPIN),
     +        ROTOT(NDIM1,NSPIN)
      DIMENSION VNLIJ(0:NDIM20,0:NDIM20,NLMAX,NTYPMX,NSPIN)
      DOUBLE PRECISION VALPHA(2)
      INTEGER NATOM(NDIM9)
      CHARACTER*80 TEXT
C.....Running parameters
      REAL RUNPAR(NRUN)
C.....Arrays used for mixing in MIXVG
      COMPLEX HKMIX(0:NDIM12,0:NDIM12,NSPIN),
     +  XKMIX(NDIM12,NSPIN), FKMIX(NDIM12,NSPIN)
C.....Common block for keeping track of possible warnings
      COMMON /WARN/ IWARN
C
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
      DOUBLE PRECISION ABOHR,RYEV,RYDERG,PI,SPI
      COMMON /CONST/   ABOHR,RYEV,RYDERG,PI,SPI
C
C-----------------------------------------------------------------------
C
C     Reading the potential file:
C
      DO 103 ITYPE = 1, NTYPMX
        POTTYP(ITYPE) = 0
        POTXC(ITYPE)  = 0
103     CONTINUE
C     LCHARG checks whether charge density was read in:
      LCHARG = .FALSE.
      IRECTP = 11
C
C     Read the file
C
105   CALL POTRED (NPRGR,IVERSN,IDAY,MONTH,IYEAR,TEXT,
     +  IRECTP,NAT,NGTOT,MINDIM,NSPIN,
     +  ULA,A1,A2,A3,B1,B2,B3,STRAIN,TY,XKAPA,NATOM,
     +  G2MAX,NG1,NG2,NG3,EPSILO,NTAB,G2TAB,
     +  IEXCOR,EXFACT,NEL,POTTYP,POTPAR,POTFIL,
     +  IGLIST,VIONIC,VGTOT,ROTOT,RUNPAR,NMIX,HKMIX,XKMIX,FKMIX,POTWRK,
     +  NDIM1,NDIM2,NDIM9,NDIM12,NTYPMX,NPARAM,NRUN,
     +  IN213)
C
      IF (IRECTP .EQ. 13 .AND. NPRGR .EQ. 213) THEN
C       For an initial guess of screening, the density of overlapping
C       free atoms may be used.  The choice is made later.
        LCHARG = .TRUE.
        WRITE (IUNIT7,107) NPRGR
107     FORMAT (' Program K',I3,' has generated overlapping atom ',
     +          'charge density as initial guess')
        ENDIF
C
      IF (IRECTP .GT. 0) GOTO 105
C
      WRITE (ISUMRY,'(1X,A)') TEXT
C     CLOSE THE FILE:
      CLOSE (UNIT = IN213 , STATUS = 'KEEP')
C
C-----------------------------------------------------------------------
C
C     VOLUME OF THE UNIT CELL
C
      VOLUM = A1(1)*A2(2)*A3(3) + A2(1)*A3(2)*A1(3) +
     +        A3(1)*A1(2)*A2(3) - A1(3)*A2(2)*A3(1) -
     +        A2(3)*A3(2)*A1(1) - A3(3)*A1(2)*A2(1)
      VOLUM = ABS(VOLUM)
C
C     Calculate the G=0 pseudopotential contribution to total energy
C     (ALPHAZ will also calculate the ionic charges from parameters)
C
      CALL ALPHAZ (POTTYP,TY,NAT,POTPAR,NTYPMX,NDIM9,NPARAM,ULA,VOLUM,
     +             ZKAPA,VALPHA)
C
C-----------------------------------------------------------------------
C
C     Determine basic properties of potentials
C
      NONLOC = .FALSE.
C
      DO 150 ITYPE = 1, NTYPMX
C
      IF (POTTYP(ITYPE) .EQ. 0) GOTO 150
      KAPA = KFIND (ITYPE,TY,NAT)
      IF (KAPA .EQ. 0) GOTO 150
C     Atomic number:
      NNATOM = NATOM(KAPA)
C
C     Printing a summary of potentials:
C
      WRITE (IUNIT8,*) 'The pseudopotential parameters are:'
C
C     Find the number of non-zero parameters:
      DO 135 NP = NPARAM, 1, - 1
        IF (ABS(POTPAR(ITYPE,NP)) .GT. 1.0E-20) GOTO 137
135     CONTINUE
      NP = 1
137   WRITE (IUNIT8,140) ITYPE, NNATOM, POTTYP(ITYPE),
     +  (POTPAR(ITYPE,J),J=1,NP)
140   FORMAT (' Atom type ',I3,' with atomic number ',I4,
     +  ' is described by potential type ',I3/100(1X,6G13.7/) / )
C
        IF (POTTYP(ITYPE) .EQ. 10) THEN
C         Jellium potential does not care which correlation
          POTXC(ITYPE) = 0
C
        ELSE IF (POTTYP(ITYPE) .EQ. 11) THEN
C         Coulomb potential does not care which correlation
          POTXC(ITYPE) = 0
C
        ELSE IF (POTTYP(ITYPE) .EQ. 13) THEN
C         Slater X-alpha
          POTXC(ITYPE) = 1
C
        ELSE IF (POTTYP(ITYPE) .EQ. 15) THEN
C         Slater X-alpha
          POTXC(ITYPE) = 1
C
        ELSE IF (POTTYP(ITYPE) .EQ. 16) THEN
          NONLOC = .TRUE.
C         Wigner interpolation:
          POTXC(ITYPE) = 2
C
        ELSE IF (POTTYP(ITYPE) .EQ. 17) THEN
          NONLOC = .TRUE.
C         Ceperley-Alder
          POTXC(ITYPE) = 3
C
        ELSE IF (POTTYP(ITYPE) .EQ. 20) THEN
          NONLOC = .TRUE.
C
C         Read in the nonlocal matrix elements from file
C
          CALL VNLRED (ITYPE,POTFIL,POTPAR,VNLIJ,VOLUM,ULA,
     +                 NTYPMX,NPARAM,NDIM20,NLMAX,NSPIN)
          POTXC(ITYPE) = NINT( POTPAR(ITYPE,10) )
C
        ELSE
C
          WRITE (IOUT,*) 'POTGET: Illegal potential type: ',
     +      POTTYP(ITYPE)
          CALL EXIT
          ENDIF
C
150     CONTINUE
C
C-----------------------------------------------------------------------
C
C    Check consistency of correlation types (IEXCOR given from K213)
C
      DO 200 ITYPE = 1, NTYPMX
        IF (POTTYP(ITYPE) .EQ. 0) GOTO 200
C       Some potentials (e.g. Coulomb) don't care
        IF (POTXC(ITYPE) .EQ. 0) GOTO 200
        IF (IEXCOR .NE. POTXC(ITYPE)) THEN
          WRITE (IOUT,210) (I,POTTYP(I),POTXC(I), I = 1, NTYPMX)
210       FORMAT('0POTGET *** ERROR *** Correlation types disagree:'/
     +      100(' ITYPE=',I2,' POTTYP=',I3,' XC=',I2/) )
          CALL EXIT
          ENDIF
200     CONTINUE
C
C-----------------------------------------------------------------------
C
      WRITE(IUNIT8,220)
220   FORMAT(' From the parameter arrays',
     +  ' the following ionic charges were deduced:'/
     +  '  Type     Charge')
      DO 250 ITYPE = 1, NTYPMX
        KAPA = KFIND (ITYPE,TY,NAT)
        IF (KAPA .GT. 0) WRITE (IUNIT8,240) ITYPE,ZKAPA(KAPA)
240     FORMAT(1X,I5,F15.8)
250     CONTINUE
C
      ZTOT = 0.0
      DO 260 KAPA = 1, NAT
        ZTOT = ZTOT + ZKAPA(KAPA)
260     CONTINUE
      WRITE (IOUT,*) 'The total ionic charge is = ', ZTOT
      IF ( ABS( FLOAT(NEL) - ZTOT ) .GT. 1.0E-8) THEN
        WRITE (IOUT,*) '*** WARNING *** number of electrons =', NEL
        IWARN = 1
        ENDIF
C
      RETURN
      END
      INTEGER FUNCTION KFIND (ITYPE,TY,NAT)
      INTEGER TY(NAT)
C
C     Find a KFIND=KAPA for which TY(KAPA)=ITYPE
C
      KFIND = 0
      DO 100 KAPA = 1, NAT
        IF (TY(KAPA) .EQ. ITYPE) THEN
          KFIND = KAPA
          RETURN
          ENDIF
100     CONTINUE
      RETURN
      END
