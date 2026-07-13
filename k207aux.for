      SUBROUTINE QUESTN (
     +  NDIM1,NDIM2,NDIM3,NDIM4,NDIM5,NDIM6,NDIM7,NDIM8,NDIM9,
     +  NDIM10,NDIM13,NSPIN,
     +  IFIRST,LDF,IEXCOR,EXFACT,DNEL,MINBDS,NBDS,DEL,NEXTRA,ISNK,
     +  E1,E2,EBAR,CONV1,G2MAX,GMTRIC,
     +  LSEMI,LRWH,LRWV,LEXIST,LITER,LOWDIN,LUMTX,
     +  IPREV,NITMX,RES,FAC,IWRCUA,IEVAL,
     +  NUMKPT,WVLIST,WVK,ULA,B1,B2,B3,NG1,NG2,NG3,IGLIST,LIST,NGTOT,
     +  IWF5,LISTAB,GK,MAXAB,ISTRAI,
     +  NDSPL1,NDSPL2,IVGDIS,ISWCH,ISELFC,
     +  JINIT1,JINIT2,JINIT3,JSTEP1,JSTEP2,JSTEP3,N1,N2,N3,
     +  INIXC4,NTAB,ROTOT,NEWTOT,MUXCG,RSPACE,CFTWRK,VOLUM,SUMMU,SUMEXC,
     +  LSTRES,NLSTRS,LFORCE,NLFORC,NAT,INVERS,IFORC,NOFORC,NONLOC)
C
C     Ask user for values of controlling parameters
C
C.....LATTICE STRUCTURE
      REAL B1(3),B2(3),B3(3), GMTRIC(3,3)
      INTEGER NTAB(NDIM2)
      INTEGER IGLIST(3,NDIM1)
      INTEGER LIST(-NG1:NG1,-NG2:NG2,-NG3:NG3)
C.....SPECIAL POINTS
      REAL WVLIST(3,NDIM10)
C.....FORCE ARRAYS
      INTEGER IFORC(NDIM9),INVERS(NDIM9)
      LOGICAL NONLOC,LFORCE,NLFORC,LSTRES,NLSTRS,NOFORC(NDIM9)
C.....ITERATIVE DIAGONALIZATION AND LOWDIN PERTURBATION
      LOGICAL LRWH,LRWV,LEXIST,LITER,LOWDIN,LUMTX
C.....WORKSPACES
      COMPLEX RSPACE(NDIM6,NSPIN)
      REAL CFTWRK(4*NDIM7)
      INTEGER IWF5(14*NDIM5)
C.....POTENTIALS AND DENSITIES
      COMPLEX ROTOT(NDIM1,NSPIN), MUXCG(NDIM1,NSPIN)
      DOUBLE PRECISION SUMMU,SUMEXC
C.....A AND B WAVES
      DIMENSION LISTAB(NDIM5),GK(4,NDIM5)
C.....SCRATCH
      REAL WVK(3)
      LOGICAL LSEMI
      DOUBLE PRECISION DEL, CONV1, DNEL
C.....COMMON BLOCK FOR KEEPING TRACK OF POSSIBLE WARNINGS
      COMMON /WARN/ IWARN
C.....TEXT
      CHARACTER*20 LDF(0:3)
      CHARACTER*25 ASTRIX,YESNO
C
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
      DOUBLE PRECISION ABOHR,RYEV,RYDERG,PI,SPI
      COMMON /CONST/   ABOHR,RYEV,RYDERG,PI,SPI
C
C     FORMAT FOR WRITING A LINE OF - OR * ON OUTPUT
      DATA ASTRIX /'(1X,79(''*''))'/
      DATA YESNO /'('' (+1=yes; -1=no)'')'/
C
C-----------------------------------------------------------------------
C
      IF (IFIRST .GT. 0) THEN
        WRITE (IOUT,*)
     +    'Do you wish to change the parameters determining the run ?'
        WRITE (IOUT,YESNO)
        READ (INPUT,*) IYES
        WRITE (IUNIT7,*) IYES
        IF (IYES .EQ. -1) GOTO 1000
        ENDIF
C
      WRITE (IUNIT7,ASTRIX)
      WRITE (IOUT,320)
320   FORMAT(' Which of the following quantities are ok (answer: +1)'/
     +  ' and which are to be changed (answer -1) ? :'/)
      IF (IFIRST .GT. 0) GOTO 520
C
      CALL LDFTYP (LDF)
      WRITE (IOUT,325) LDF(IEXCOR)
325   FORMAT(1X,A,' local density functional - OK ?')
      READ (INPUT,*) IOK
      IF (IOK .LE. 0 .OR. IEXCOR .LE. 0) THEN
329     WRITE (IUNIT7,330) (I, LDF(I), I = 1, 3)
330     FORMAT(' Choose density-functional:'/
     +  10(I4,': ',A/))
        READ (INPUT,*) IEXCOR
        WRITE(IUNIT7,*) IEXCOR,' LDF TYPE'
        IF (IEXCOR .LT. 1 .OR. IEXCOR .GT. 3) GOTO 329
        ENDIF
C
      IF (IEXCOR .EQ. 1) THEN
C       X-ALPHA
        WRITE (IOUT,340) EXFACT
340     FORMAT (' EXFACT (ALPHA) = ',F10.5,' - O.K. ?')
        READ (INPUT,*) IOK
        IF (IOK .LE. 0) THEN
          WRITE (IOUT,*) 'Enter a new exchange factor ?'
          READ (INPUT,*) EXFACT
          WRITE (IOUT,*) EXFACT,' exchange-correlation factor'
          ENDIF
      ELSE
        EXFACT = 1.0
        ENDIF
C
C     Occupation of bands
C
      WRITE (IOUT,450)
450   FORMAT (' Do you wish to assume that the material is'/
     +  ' a semiconductor? yes=1, no=-1' )
      READ (INPUT,*) IYESNO
      IF (IYESNO .GT. 0) THEN
        LSEMI = .TRUE.
        WRITE (IUNIT7,*) 'Semiconductor material'
      ELSE
        LSEMI = .FALSE.
        WRITE (IUNIT7,*) 'Metallic material'
        ENDIF
C
      WRITE (IOUT,460) DNEL
460   FORMAT ('0Number of electrons = ',F10.6,' - O.K.?')
      READ (INPUT,*) IYES
      WRITE (IUNIT7,*) IYES,' answered'
      IF (IYES .LE. 0) THEN
465     WRITE (IOUT,470)
470     FORMAT('0How many electrons do you want to occupy the bands ?')
        READ (INPUT,*) DNEL
        IF (DNEL .LT. 1.0D0 .OR. DNEL .GT. DBLE(2*NDIM8)) THEN
          WRITE(IOUT,*) 'Unreasonable NEL =',DNEL,' (NDIM8=',NDIM8,')'
          GOTO 465
          ENDIF
        WRITE (IOUT,*) DNEL,' electrons in the occupied bands'
        ENDIF
C
C     Minimum number of bands necessary to accomodate the electrons
C     = number of electrons / 2; add a 1 in the case of a metal
      NEL = DNEL + 1.0D-8
      MINBDS = (NEL + 1) / 2
      IF (MINBDS .GT. NDIM8) THEN
        WRITE (IOUT,*) 'MINBDS = ',MINBDS,
     +    ' but dimension NDIM8 = ',NDIM8
C       Error - too many electrons in the bands
        WRITE (IOUT,2320) NEL,NEL/2,NDIM8
2320    FORMAT('0***FATAL ERROR*** NEL = ',I4/'0however, NEL/2 = ',
     +    I4,' (number of bands) > NDIM8 = ',I4)
        CALL EXIT
        ENDIF
C
C     Semiconductor:
C
      IF (LSEMI) THEN
C
        DNEL = DBLE(NEL)
        IF (MOD(NEL,2) .EQ. 1) THEN
          IF (NSPIN .EQ. 1) THEN
            WRITE (IOUT,*)
     +        'WARNING *** this cannot be semiconductor when NEL=',NEL
            IWARN = 2
          ELSE
            WRITE (IOUT,*) 'Spin-symmetry must break since NEL=',NEL
            ENDIF
          ENDIF
C
C       No. of eigenvalues
        NBDS = MINBDS
        WRITE (IOUT,490) NBDS
490     FORMAT (1X,I3,' Eigenvalues requested O.K. ?')
        READ (INPUT,*) IOK
        WRITE (IUNIT7,*) IOK,' ANSWERED'
        IF (IOK .LE. 0) THEN
495       WRITE (IOUT,*) 'Enter number of last band to be calculated'
          READ (INPUT,*) NBDS
          IF (NBDS .LT. MINBDS) THEN
            WRITE (IOUT,500) NBDS,MINBDS
500         FORMAT(' NBDS=',I4,' is < MINBDS=',I4,'; NOT ADVISABLE',
     +      ' TRY AGAIN ...')
            GOTO 495
            ENDIF
          IF (NBDS .GT. NDIM8) THEN
            WRITE (IOUT,502) NBDS,NDIM8
502         FORMAT (' NBDS=',I4,' IS > NDIM8=',I4,'; TRY AGAIN')
            GOTO 495
            ENDIF
          WRITE (IOUT,*) 'Evaluating bands up to ',NBDS
        ENDIF
C
      ELSE
C
C       A metal
C
        LSEMI = .FALSE.
        WRITE(IOUT,*) 'Enter the smearing parameter DEL in eV'
        READ (INPUT,*) DEL
        WRITE(IUNIT7,510) DEL
510     FORMAT (' Smearing parameter = ',F8.3,' eV')
        WRITE(IOUT,*) 'Enter number of bands more than MINBDS required'
        WRITE(IOUT,*) 'enter with minus sign to keep same at each cycle'
        READ (INPUT,*) NEXTRA
        ISNK = 0
        IF (NEXTRA .LT. 0) ISNK = 1
        NEXTRA = IABS(NEXTRA)
        NBDS = MINBDS + NEXTRA
        WRITE(IUNIT7,*) NEXTRA,' extra bands required'
C
        ENDIF
C------------
520   IF (IFIRST .GT. 0) THEN
        WRITE (IOUT,530)
     +    E1,E2,E1*CONV1,E2*CONV1,E1*CONV1*RYEV,E2*CONV1*RYEV
530     FORMAT(' Cut-off energies on plane waves :'/
     +    '      E1',13X,'E2',9X,'IN UNITS:'/
     +    1X,2G15.7,'  (2*PI/ULA)**2)'/1X,2G15.7,'  RYD'/
     +    1X,2G15.7,'  EV')
        WRITE (IOUT,*) ' E1,E2 - O.K. ?'
        READ (INPUT,*) IOK
        IF (IOK .GE. 0) GOTO 600
        ENDIF
540   WRITE (IOUT,550)
550   FORMAT(' Enter cut-off energies for plane waves: E1 (A waves)'
     +  ,' E2 (B waves)'/
     +' a n d  define the energy units by 1=(2PI/ULA)**2,  2=Rydberg, 3=
     +eV.'/' for anisotropic cutoff, enter units with a minus sign.')
      READ (INPUT,*) E1,E2,IUNIT
C
      IF (ABS(IUNIT) .EQ. 3) THEN
        E1 = E1/RYEV/CONV1
        E2 = E2/RYEV/CONV1
      ELSE IF (ABS(IUNIT) .EQ. 2) THEN
        E1 = E1/CONV1
        E2 = E2/CONV1
      ELSE IF (ABS(IUNIT) .NE. 1) THEN
        WRITE (IOUT,*) 'INVALID UNITS:',IUNIT
        GOTO 540
        ENDIF
      IF (E1. LT. 0.1 .OR. E2 .LT. 0.1 .OR. E2 .LT. (E1-1.0E-6)) THEN
        WRITE (IOUT,*) 'E1 AND/OR E2 UNREASONABLE:',E1,E2
        GOTO 540
        ENDIF
C
      WRITE (IOUT,530)
     +  E1,E2,E1*CONV1,E2*CONV1,E1*CONV1*RYEV,E2*CONV1*RYEV
      WRITE (ISUMRY,580) E1*CONV1,E2*CONV1
580   FORMAT(' Plane wave cutoffs: E1,E2 =',2F12.3,' RYD')
C
      IF (4.0*E2 .GT. G2MAX) THEN
        WRITE (IOUT,590) G2MAX,E2
590     FORMAT(' *** WARNING *** G2MAX = ',F10.3,' of ionic potential'/
     +  ' is too small for treating all plane waves when E2 = ',F10.3)
        WRITE (IOUT,*) 'This is potentially dangerous, so we quit'
        CALL EXIT
        IWARN = 2
        ENDIF
C
      IF (E2 .GT. E1+1.0E-4) THEN
591     WRITE(IOUT,*) 'Lowdin perturbation (1) or iterative ',
     +                'diagonalization (2) ?'
        READ (INPUT,*) ILOWIT
        IF (ILOWIT .EQ. 1) THEN
          LOWDIN = .TRUE.
          WRITE(IOUT,*) 'Lowdin perturbation will be used'
          WRITE(IOUT,*) 'Do you require the U matrix ?'
          READ(INPUT,*) IYESNO
          WRITE(IOUT,*) IYESNO,' ANSWERED'
          IF (IYESNO .EQ. 1) LUMTX = .TRUE.
          WRITE(IOUT,*) 'Enter switch for eigenvalues'
          WRITE(IOUT,*) '0=1st order  1=variational  2=subspace'
          READ(INPUT,*) IEVAL
          WRITE(IOUT,*) IEVAL,' ANSWERED'
          WRITE(IOUT,*) 'Enter EBAR in ryd (for U matrix)'
          READ(INPUT,*) EBAR
          WRITE(IOUT,*) 'EBAR = ',EBAR,' rydbergs'
        ELSE IF (ILOWIT .EQ. 2) THEN
          IF (NDIM13 .LT. NDIM3) THEN
            WRITE(IOUT,*) '***ERROR*** NDIM13 IS TOO SMALL'
            WRITE(IOUT,*) 'NDIM13= ',NDIM13,' BUT NDIM3= ',NDIM3
            WRITE(IOUT,*) 'check whether this is sufficient'
            ENDIF
          LITER = .TRUE.
          WRITE(IOUT,*) 'Iterative diagonalization will be used'
          WRITE(IOUT,*) 'Enter max number of iterations'
          READ(INPUT,*) NITMX
          WRITE(IOUT,560) NITMX
560       FORMAT(' Max no of iterations = ',I4)
          WRITE(IOUT,*) 'Enter accuracy RES of eigenvalues in Ryd'
          WRITE(IOUT,*) '(Enter RES<0 for automatic adjustment ',
     +      'starting with abs(RES))'
          READ(INPUT,*) RES
          IF (RES .LE. 0.0) THEN
            IF (RES .EQ. 0.0) RES = - 1.0E-4
            WRITE(IOUT,*) 'Automatic accuracy adjustment chosen'
            ENDIF
          WRITE(IOUT,570) ABS(RES)
570       FORMAT(' Accuracy of eigenvalues (RES) = ',E16.6,' RYD')
571       WRITE(IOUT,*) 'Enter exclusion zone parameter FAC'
          READ(INPUT,*) FAC
          IF (FAC .LE. 0.0 .OR. FAC .GT. 100.0) THEN
            WRITE(IOUT,*) 'FAC should be about 0.1-2; try again'
            GOTO 571
            ENDIF
          WRITE(IOUT,573) FAC
573       FORMAT(' Exclusion zone parameter FAC = ',F10.6)
          WRITE(IOUT,*) 'Enter cycle on which start to read ',
     +      'eigenvectors and eigenvalues from file'
          READ (INPUT,*) IPREV
          WRITE(IOUT,575) IPREV
575       FORMAT(' Eigen data from previous cycle used on cycle ',I4)
        ELSE
          WRITE(IOUT,*) '***ERROR*** bad choice - please repeat'
          GOTO 591
          ENDIF
        ENDIF
C
C................
C     USE OF FILES FOR HAMILTONIAN AND NON LOCAL MATRIX ELEMENTS
      IF (LOWDIN .OR. LITER) LRWH = .TRUE.
      IF (NONLOC) THEN
        WRITE(IOUT,*) 'Do you want a file for non-local ',
     +                'matrix elements ?'
        READ(INPUT,*) IYESNO
        WRITE(IOUT,*) IYESNO,' ANSWERED'
        IF (IYESNO .EQ. 1) LRWV = .TRUE.
        WRITE(IOUT,*) 'Do you want to read from an existing file ?'
        READ(INPUT,*) IYESNO
        WRITE(IOUT,*) IYESNO,' ANSWERED'
        IF (IYESNO .EQ. 1) LEXIST = .TRUE.
        ENDIF
C................
      IWRCUA = 0
      IF (LOWDIN) THEN
        IF (LUMTX) THEN
          IF (IEVAL .GT. 0) THEN
            IWRCUA = 3
          ELSE
            IWRCUA = 1
            ENDIF
        ELSE
          IF (IEVAL .GT. 0) THEN
            IWRCUA = 2
          ELSE
            IWRCUA = 1
            ENDIF
          ENDIF
      ELSE IF (LITER) THEN
        IWRCUA = 2
        ENDIF
C
C................
C     PLANE WAVE CUTOFF: ISOTROPIC OR ANISOTROPIC ?
      IF (IUNIT .GT. 0) THEN
C       ISOTROPIC CUTOFF (SPHERE IN THE UNSTRAINED CRYSTAL)
        DO 595 I = 1,3
        DO 595 J = 1,3
          GMTRIC(I,J) = 0.0
          IF (I .EQ. J) GMTRIC(I,J) = 1.0
595       CONTINUE
        IF (ISTRAI .GT. 0) WRITE (IOUT,*)
     +    '*** isotropic plane wave cutoff applied ***'
      ELSE IF (ISTRAI .GT. 0) THEN
        WRITE (IUNIT7,*) 'anisotropic plane wave cutoff applied'
        ENDIF
C
C.......................................................................
      WRITE (IOUT,*) 'Do you wish to test dimension of hamiltonian ?'
      WRITE (IOUT,YESNO)
      READ (INPUT,*) IYESNO
      IF (IYESNO .NE. 1) GOTO 600
C
C     Find which plane waves are to be included for this E1,E2
C
      IF (NUMKPT .LE. 0) THEN
C
C       TEST NUMBER OF PLANE WAVES FOR K=0 (IF NO SPECIAL POINTS)
C
        DO 597 I = 1,3
          WVK(I) = 0.0
597       CONTINUE
        WRITE (IUNIT7,*) 'The case of k=(0,0,0):'
        CALL GKCUT (WVK,E1,E2,B1,B2,B3,NG1,NG2,NG3,
     +    NDIM1,NDIM3,NDIM4,NDIM5,IWF5,IWF5(NDIM5+1),
     +    ISTRAI,GMTRIC,IGLIST,NGTOT,LISTAB,GK,NAMAX,NBMAX)
C
      ELSE
C
C       TEST NUMBER OF PLANE WAVES FOR SPECIAL POINTS SET
C
        NAMAX = - 1
        NBMAX = - 1
        DO 598 ISPPT = 1,NUMKPT
          CALL GKCUT (WVLIST(1,ISPPT),E1,E2,B1,B2,B3,NG1,NG2,NG3,
     +      NDIM1,NDIM3,NDIM4,NDIM5,IWF5,IWF5(NDIM5+1),
     +      ISTRAI,GMTRIC,IGLIST,NGTOT,LISTAB,GK,NA,NB)
          NAMAX = MAX(NAMAX,NA)
          NBMAX = MAX(NBMAX,NB)
598       CONTINUE
        WRITE (IUNIT7,*) 'Going through the special points set,'
C
        ENDIF
C
      WRITE (IUNIT7,599) NAMAX,NBMAX
      WRITE (ISUMRY,599) NAMAX,NBMAX
599   FORMAT(1X,I4,
     +  ' A and ',I4,' B waves are required')
      IF (NAMAX .LE. 0 .OR. NBMAX .LT. 0 .OR.
     +    NAMAX .GT. NDIM3 .OR. NBMAX .GT. NDIM4) THEN
        WRITE (IOUT,*) 'These numbers are unreasonable'
        GOTO 540
        ENDIF
      WRITE (IOUT,*) ' OK ? (YES=1, NO=-1)'
      READ (INPUT,*) IYES
      IF (IYES .LE. 0) GOTO 540
      MAXAB = MAX(NAMAX,NBMAX)
C
C...............
600   NDSPL1 = 1
      NDSPL2 = NBDS
      WRITE (IOUT,640) NDSPL1,NDSPL2
640   FORMAT (' Band structure display - NDSPL1,NDSPL2 = ',2I4,
     +  ' - O.K. ?')
      READ (INPUT,*) IOK
      IF (IOK .LE. 0) THEN
        WRITE (IOUT,650)
650     FORMAT (' Which eigenenergies do you wish to be displayed?'/
     +  ' enter the lowest and the highest order number (=from-to)'/
     +  ' (enter 0,0 for no display)')
        READ (INPUT,*) NDSPL1,NDSPL2
        IF (NDSPL1 .LT. 0 .OR. NDSPL2 .LT. 0) THEN
          WRITE (IOUT,*) 'Illegal NDSPL1,2'
          GOTO 600
          ENDIF
        NDSPL1 = MIN(NDSPL1,NDIM8)
        NDSPL2 = MIN(NDSPL2,NDIM8)
        WRITE (IOUT,*) NDSPL1,NDSPL2,' energy bands displayed'
        ENDIF
C...............
660   WRITE (IOUT,670) IVGDIS
670   FORMAT (' V(G) display: ',I4,' shells shown - O.K. ?')
      READ (INPUT,*) IOK
      IF (IOK .LE. 0) THEN
        WRITE (IOUT,680)
680     FORMAT (' Old and new potentials V(G) as well as RO(G) will be',
     +  ' displayed'/' for a fixed number of shells'/
     +  ' enter the number of shells (0 for no display)')
        READ (INPUT,*) IVGDIS
        IF (IVGDIS .LT. 0 .OR. IVGDIS .GT. NDIM2) THEN
          WRITE (IOUT,*) 'Illegal IVGDIS'
          GOTO 660
          ENDIF
        WRITE (IOUT,*) IVGDIS,' potential display'
        ENDIF
C...............
      WRITE (IOUT,700) ISWCH
700   FORMAT (' Value of ISWCH =',I3,' - OK ?')
      READ (INPUT,*) IOK
      IF (IOK .LE. 0) THEN
        WRITE (IOUT,710)
710     FORMAT(' Values of ISWCH:'/
     +' 1:  full calculation with possibly energy forces and stress',/,
     +' 2:  total energy only from bands, no forces or stress',/,
     +' 3:  eigenvalues only',/,
     +' 4:  eigenvalues and eigenvectors only',/,
     +' ENTER ISWCH')
        READ (INPUT,*) ISWCH
        WRITE (IOUT,*) ISWCH,' value of ISWCH'
        WRITE (ISUMRY,*) ISWCH,' value of ISWCH'
        ENDIF
C...............
C     THE FOLLOWING QUESTIONING IS SUPERFLUOUS IF WE ONLY WANT BAND E(K)
      IF (NUMKPT .LE. 0  .OR.  ISELFC .LE. 0) GO TO 1000
C...............
      WRITE (IUNIT7,740) JINIT1,JINIT2,JINIT3,JSTEP1,JSTEP2,JSTEP3
740   FORMAT (' RO(R) display:',10X,'JINIT1  JINIT2  JINIT3  ',
     +  'JSTEP1  JSTEP2  JSTEP3'/23X,6I8)
      WRITE (IOUT,*) 'JINIT, JSTEP - O.K.?'
      READ (INPUT,*) IOK
      IF (IOK .LE. 0) THEN
        WRITE (IOUT,750)
750     FORMAT (' For display of RO in direct space: enter 3 integers',
     +  ' defining the starting'/
     +  ' point on the real-space mesh (as defined e.g. in EXCH4)'/
     +  ' and enter 3 integers defining the increments',
     +  ' (-1,0,0,0,0,0 for no display)')
        READ (INPUT,*) JINIT1,JINIT2,JINIT3,JSTEP1,JSTEP2,JSTEP3
        WRITE (IUNIT7,740) JINIT1,JINIT2,JINIT3,JSTEP1,JSTEP2,JSTEP3
        ENDIF
C...............
760   WRITE (IOUT,770) N1,N2,N3
770   FORMAT (' Fourier transform: N1,N2,N3 = ',3I4,' - O.K. ?')
      READ (INPUT,*) IOK
      IF (IOK .LE. 0) THEN
        WRITE (IOUT,780)
780     FORMAT (' Enter 3 integers defining the mesh in real space',
     +    ' for the double fourier'/
     +    ' transform (for definition, see subroutine EXCH4)')
        READ (INPUT,*) N1,N2,N3
        WRITE (IOUT,*) N1,N2,N3,' Fourier transform parameters'
C       Having changed N1,2,3, we have to reinitialize CFFT
        INIXC4 = 0
        ENDIF
      IF (N1*N2*N3 .GT. NDIM6  .OR.  MAX(N1,N2,N3)/4 .GT. NDIM7) THEN
C       ERROR - N1,N2,N3 IS TOO BIG
        WRITE (IOUT,800) NDIM6,NDIM7
800     FORMAT(' ERROR - N1,N2,N3 exceed dimensions NDIM6,NDIM7 =',
     +    2I10/' - try again')
        GOTO 760
        ENDIF
C...............
      IF (INIXC4 .NE. 1) THEN
        NGMAX = NGTOT
        CALL EXCH4 (INIXC4,NDIM1,NDIM2,NDIM6,NDIM7,NG1,NG2,NG3,NGMAX,
     +    IGLIST,NGTOT,NTAB,ROTOT,N1,N2,N3,NEWTOT,MUXCG,RSPACE,
     +    CFTWRK,VOLUM,ULA,EXFACT,IEXCOR,NSPIN,SUMMU,SUMEXC)
        INIXC4 = 1
        ENDIF
C-----------------------------------------------------------------------
C
C     CALCULATE FORCES AND STRESS ?
C
      LSTRES = .FALSE.
      NLSTRS = .FALSE.
      LFORCE = .FALSE.
      NLFORC = .FALSE.
      IF (ISWCH .GT. 2) GOTO 1000
      WRITE (IOUT,825)
825   FORMAT (' Calculate stress tensors ? (yes=1, no=-1)')
      READ (INPUT,*) IYES
      WRITE (IOUT,*) IYES,' stress tensor calculation switch'
      IF (IYES .NE. -1) THEN
        LSTRES = .TRUE.
        IF (NONLOC) NLSTRS = .TRUE.
        ENDIF
C
      WRITE (IOUT,830)
830   FORMAT(' Calculate forces ? (1=yes, -1=no, 2=yes,+save time)')
      READ (INPUT,*) IYES
      WRITE (IOUT,*) IYES,' force calculation switch'
      IF (IYES .EQ. 1 .OR. IYES .EQ. 2) THEN
        LFORCE = .TRUE.
        IF (NONLOC) NLFORC = .TRUE.
        ENDIF
C
      IF (ISWCH .NE. 1) THEN
        WRITE (IOUT,832) ISWCH
832     FORMAT(' WARNING: Since ISWCH = ',I2/
     +    ' neither force nor stress will be calculated'/)
        LSTRES = .FALSE.
        NLSTRS = .FALSE.
        LFORCE = .FALSE.
        NLFORC = .FALSE.
        ENDIF
C...............
      DO 840 KAPA = 1,NAT
        INVERS(KAPA) = 0
        IFORC(KAPA)  = 1
        NOFORC(KAPA) = .FALSE.
840     CONTINUE
      IF (.NOT.NONLOC .OR. IYES .NE. 2) GOTO 1000
C     ASK THE USER FOR TIME SAVING REDUCTIONS IN FORCE CALCULATIONS
845   WRITE (IOUT,850)
850   FORMAT(' Enter the numbers of the atoms whose forces will be',
     +  ' calculated'/' one by one, terminate by a zero')
      DO 860 KAPA = 1,NAT
        INVERS(KAPA) = 0
        IFORC(KAPA) = 0
        NOFORC(KAPA) = .TRUE.
860     CONTINUE
      DO 870 KAPA = 1,NAT
        READ (INPUT,*) IFOR
        IF (IFOR.LE.0 .OR. IFOR.GT.NAT) GOTO 880
        IFORC(IFOR) = IFOR
        NOFORC(IFOR) = .FALSE.
870     CONTINUE
880   WRITE(IOUT,*)
     +  'for the given atoms, indicate other atoms related by inversion'
      WRITE(IOUT,*) 'Enter 0 if no atom, -1 to skip questions'
      DO 890 KAPA = 1,NAT
        IF (IFORC(KAPA) .EQ. 0) GOTO 890
        WRITE(IOUT,*) 'Atom number',KAPA
        READ (INPUT,*) INV
        IF (INV .LT. 0) GOTO 900
        IF (INV .GT. 0) INVERS(KAPA) = INV
890     CONTINUE
900   WRITE (IOUT,910) IFOR
C
910   FORMAT(' Summary: atoms included / atoms related by inversion'/
     +  1X,25I3)
      WRITE (IUNIT7,920) INVERS
920   FORMAT(1X,25I3)
      WRITE(IOUT,*) 'Satisfied (yes=1, no=-1)'
      READ (INPUT,*) IYES
      IF (IYES .EQ. -1) GOTO 845
C
1000  RETURN
      END
      SUBROUTINE TPRINT (IUNIT, TEXT, T, NPFS)
C
C     Print CPU time and page faults (if any) on unit IUNIT
C
      CHARACTER*(*) TEXT
      PARAMETER (NTEXT = 54)
      CHARACTER FORM*80
C.....FILES
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
C-----------------------------------------------------------------------
C
      IF (NPFS .GT. 0) THEN
        LEN = MIN( NTEXT, LENSTR(TEXT))
        WRITE (FORM,110) NTEXT
        WRITE (IUNIT,FORM) TEXT(1:LEN), T, NPFS
      ELSE
        NXTRA = 13
        LEN = MIN( NTEXT + NXTRA, LENSTR(TEXT))
        WRITE (FORM,110) NTEXT + NXTRA
        WRITE (IUNIT,FORM) TEXT(1:LEN), T
        ENDIF
110   FORMAT ( '(1X,A,T',I2.2,',F10.3,'' Secs'',:,I9,'' PFs'')' )
      RETURN
      END
      SUBROUTINE EDISPL (WVK,NDSPL1,NDSPL2,LITER,LOWDIN,NA,NB,IEVAL,
     +                   EE1,EE2,OCC,IOCC,ISPIN,NSPIN)
      REAL WVK(3),EE1(NDSPL2),EE2(NDSPL2)
      DOUBLE PRECISION OCC(NDSPL2)
      LOGICAL LITER,LOWDIN
C
C     BLOCK FOR DISPLAY OF SELECTED ENERGIES IN EE1 AND EE2
C
C.....FILES
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
      CHARACTER*12 SPIN(2)
      DATA SPIN /' (spin up)', ' (spin down)'/
C-----------------------------------------------------------------------
      IF (NDSPL1 .LE. 0 .OR. NDSPL2 .LE. 0) RETURN
      IF (NSPIN .EQ. 2) THEN
        WRITE (IOUT,1090) WVK, SPIN(ISPIN)
      ELSE
        WRITE (IOUT,1090) WVK, ' '
        ENDIF
1090  FORMAT ('0Wavevector: WVK =',3F14.6,A)
C
      IF (LITER) THEN
        WRITE (IUNIT7,1110) NA,NB
1110    FORMAT (1X,I4,' A waves ,',I4,' B waves ')
        IDISPL = 1
      ELSE IF (LOWDIN) THEN
        WRITE (IUNIT7,1110) NA,NB
        IDISPL = 1
        IF (IEVAL .LT. 2) IDISPL = 2
      ELSE
        WRITE (IUNIT7,*) NA,' A waves '
        IDISPL = 1
        ENDIF
C
      IF (IOCC .LT. 0) THEN
        WRITE (IUNIT7,100) ' '
      ELSE
        WRITE (IUNIT7,100) 'occupation'
        ENDIF
100   FORMAT (' Band no.   band energies in eV',T45,A)
      IF (IDISPL .EQ. 2) WRITE (IUNIT7,110)
110   FORMAT(13X,'First order   higher order ')
C
      DO 200 J = NDSPL1, NDSPL2
C
        IF (IDISPL .EQ. 1) THEN
C
          IF (IOCC .GT. 0) THEN
            WRITE (IUNIT7,120) J,EE2(J),OCC(J)
120         FORMAT (1X,I5,5X, F14.6,T45,F10.6)
          ELSE
            WRITE (IUNIT7,130) J,EE2(J)
130         FORMAT (1X,I5,5X, F14.6)
            ENDIF
C
        ELSE
C
          IF (IOCC .GT. 0) THEN
            WRITE (IUNIT7,140) J,EE1(J),EE2(J),OCC(J)
140         FORMAT (1X,I5,5X,2F14.6,T45,F10.6)
          ELSE
            WRITE (IUNIT7,150) J,EE1(J),EE2(J)
150         FORMAT (1X,I5,5X,2F14.6)
            ENDIF
C
          ENDIF
C
200     CONTINUE
C
      RETURN
      END
      SUBROUTINE GSPMTR(STRAIN,STRN,GMTRIC,ISTRAI)
      REAL STRAIN(6),STRN(3,3),GMTRIC(3,3)
C
C     CONSTRUCT THE G-SPACE METRIC TENSOR:
C
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
C-----------------------------------------------------------------------
      STRN(1,1) = 1.0 - STRAIN(1)
      STRN(2,2) = 1.0 - STRAIN(2)
      STRN(3,3) = 1.0 - STRAIN(3)
      STRN(2,3) = - 0.5*STRAIN(4)
      STRN(3,2) = STRN(2,3)
      STRN(1,3) = - 0.5*STRAIN(5)
      STRN(3,1) = STRN(1,3)
      STRN(1,2) = - 0.5*STRAIN(6)
      STRN(2,1) = STRN(1,2)
C     WRITE (IUNIT8,102) STRN
C102  FORMAT(' 1 - STRAIN TENSOR:'/3(5X,3F15.5/) )
C     INVERT
      CALL SYMINV(STRN,3,3,3,IFAIL)
      IF (IFAIL .NE. 0) GOTO 2300
      ISTRAI = 0
      SUM = - 3.0
      DO 110 I = 1,3
        DO 110 J = 1,3
          SUM = SUM + ABS(STRN(I,J))
          GMTRIC(I,J) = 0.0
          DO 110 K = 1,3
            GMTRIC(I,J) = GMTRIC(I,J) + STRN(K,I)*STRN(K,J)
110         CONTINUE
      IF (ABS(SUM) .GT. 1.0E-10) ISTRAI = 1
C     WRITE (IUNIT8,112) GMTRIC
C112  FORMAT(' G-SPACE METRIC TENSOR IS:'/3(5X,3F15.5/))
C
C     NOW DEFINE STRN AS THE TENSOR: (1+STRAIN)**(-1), FOR USE
C     IN STRAINING THE K-VECTOR.
C
      STRN(1,1) = 1.0 + STRAIN(1)
      STRN(2,2) = 1.0 + STRAIN(2)
      STRN(3,3) = 1.0 + STRAIN(3)
      STRN(2,3) = 0.5*STRAIN(4)
      STRN(3,2) = STRN(2,3)
      STRN(1,3) = 0.5*STRAIN(5)
      STRN(3,1) = STRN(1,3)
      STRN(1,2) = 0.5*STRAIN(6)
      STRN(2,1) = STRN(1,2)
C     WRITE (IUNIT8,103) STRN
C103  FORMAT(' 1 + STRAIN TENSOR:'/3(5X,3F15.5/) )
C     INVERT
      CALL SYMINV(STRN,3,3,3,IFAIL)
      IF (IFAIL .NE. 0) GOTO 2400
      RETURN
C
2300  WRITE (IOUT,*) 'COULD NOT INVERT (1 - STRAIN) MATRIX'
      CALL EXIT
      RETURN
2400  WRITE (IOUT,*) 'COULD NOT INVERT (1 + STRAIN) MATRIX'
      CALL EXIT
      RETURN
      END
      SUBROUTINE SYMINV(A,L,M,N,IFAIL)
C
C     CERNLIB F104
C     INVERT A REAL SYMMETRIC MATRIX
C     A ..... MATRIX (CONTAINS A-1 ON OUTPUT)
C     L,M ... DIMENSION OF A
C     N ..... NUMBER OF ROWS TO BE INVERTED
C     IFAIL . 0 IF OK
C             1 IF FAILURE
C
      INTEGER R
      DIMENSION A(L,M),P(100),Q(100),R(100)
C-----------------------------------------------------------------------
      IFAIL=0
      EPSILN=1.E-6
      IF(L.LT.N) GO TO 95
      IF(M.LT.N) GO TO 96
C
C                   CONSTRUCT TRUTH TABLE
C
      DO 10 I=1,N
   10 R(I)=1
C
C                   BEGIN PROGRAMME
C
      DO 65 I=1,N
      K=0
C
C                   SEARCH FOR PIVOT
C
      BIG=0.
      DO 37 J=1,N
      TEST=ABS(A(J,J))
      IF(TEST-BIG)37,37,31
   31 IF(R(J))100,37,32
   32 BIG=TEST
      K=J
   37 CONTINUE
C
C                   TEST FOR ZERO MATRIX
C
      IF(K)100,100,38
C
C                   TEST FOR LINEARITY
C
   38 IF(I.EQ.1) PIVOT1=A(K,K)
      IF(ABS(A(K,K)/PIVOT1)-EPSILN) 100,39,39
C
C                   PREPARATION FOR ELIMINATION STEP1
C
   39 R(K)=0
      Q(K)=1./A(K,K)
      P(K)=1.
      A(K,K)=0.0
      KP1=K+1
      KM1=K-1
      IF(KM1)100,50,40
   40 DO 49 J=1,KM1
      P(J)=A(J,K)
      Q(J)=A(J,K)*Q(K)
      IF(R(J))100,49,42
   42 Q(J)=-Q(J)
   49 A(J,K)=0.
   50 IF(K-N)51,60,100
   51 DO 59 J=KP1,N
      P(J)=A(K,J)
      Q(J)=-A(K,J)*Q(K)
      IF(R(J))100,52,59
   52 P(J)=-P(J)
   59 A(K,J)=0.0
C
C                   ELIMINATION PROPER
C
   60 DO 65 J=1,N
      DO 65 K=J,N
   65 A(J,K)=A(J,K)+P(J)*Q(K)
C
C                   ELEMENTS OF LEFT DIAGONAL
C
      DO 70 J=2,N
CDIR$ IVDEP
      DO 70 K=1,J-1
   70 A(J,K)=A(K,J)
      RETURN
C
C                   FAILURE RETURN
C
   95 PRINT 150,L,N
  150 FORMAT(/' L =',I5,' N =',I5,' L SHOULD BE LARGER OR EQUAL TO N')
      GO TO 100
   96 PRINT 151,M,N
  151 FORMAT(/' M =',I5,' N =',I5,' M SHOULD BE LARGER OR EQUAL TO N')
  100 IFAIL=1
      RETURN
      END
