      SUBROUTINE K213
     +  (NDIM1,NDIM2,NDIM9,NDIM12,NSPIN,NQMAX,
     +  NPARAM,NTYPMX,NRUN,RUNPAR,HKMIX,XKMIX,FKMIX,
     +  XKAPA,F0,G2TAB,TY,NTAB,IGLIST,G2LIST,NATOM,
     +  POTPAR,ZKAPA,POTTYP,POTXC,POTFIL,POTWRK,STRFAC,SATGR,SATGI,
     +  VIONIC,VGTOT,ROTOT,  DX,IWORK,ICOPY)
C
C PROGRAM K213
C                   WRITTEN ON JULY 24TH, 1981 - FROM K97
C                   MODIFIED APRIL 1982 FOR NON-LOCAL POTENTIALS
C                   BY OLE HOLM NIELSEN
C
C PREPARATION OF INPUT FILE WITH PSEUDOPOTENTIAL FOR
C RUNNING K95 (OR SIMILAR PSEUDOPOTENTIAL PROGRAMS):
C GENERATION OF RECIPROCAL LATTICE VECTORS, ESTABLISHING
C A SYNOPTIC TABLE OF RECIPROCAL SPACE, CALCULATION OF
C ATOMIC PSEUDOPOTENTIALS,
C COMBINING TOTAL IONIC PSEUDOPOTENTIAL FROM THE INDIVIDUAL
C ATOMIC ONES AND DEFINING THE INITIAL PSEUDOPOTENTIAL FOR THE
C FIRST RUN OF THE SELFCONSISTENT LOOP.
C
C THE INITIAL PSEUDOPOTENTIAL (FOR THE SELFCONSISTENT LOOP)
C IS CHOSEN TO BE THE IONIC ONE,  SCREENED
C BY THE DIELECTRIC FUNCTION OF FREE ELECTRON GAS
C
C POTENTIALS PRESENTLY AVAILABLE IN THIS PROGRAM:
C - APPELBAUM-HAMANN (SUBROUTINE POT13)
C - BERKELEY IONIC (SUBROUTINE POT15)
C - HAMANN-SCHLUTER-CHIANG NON-LOCAL (SUBROUTINE POT16)
C - BELL LABS 1982 "PERIODIC TABLE" NON-LOCAL (SUBROUTINE POT17)
C
C SUBROUTINES NEEDED: RECLAT,GCODE3,GCODE4,ALPHAZ,EPS1,
C                     POT13,POT15,POT16,PARAM(+...).
C
C     INCLUDE 'DIMENS.CMN/LIST'
C     INCLUDE '[OHN.NEW]FILES.CMN/LIST'
C     INCLUDE '[OHN.NEW]CONST.CMN/LIST'
C
      PARAMETER ( NPRGR  = 213 , IVERSN = 1 )
      REAL RUNPAR(NRUN)
      COMPLEX HKMIX(0:NDIM12,0:NDIM12,NSPIN),
     +        XKMIX(NDIM12,NSPIN), FKMIX(NDIM12,NSPIN)
      DOUBLE PRECISION POTPAR(NTYPMX,NPARAM)
C.....WORKSPACES
      DOUBLE PRECISION POTWRK(0:NQMAX,*)
      DIMENSION IWORK(NDIM1),ICOPY(NDIM1)
      COMPLEX SGVG
      DOUBLE PRECISION DG2, VG(2), VGP(2), RHOG(2),
     +  DVOLUM, VALPHA(2), DB1(3),DB2(3),DB3(3),DA1(3),DA2(3),DA3(3)
      DIMENSION A1(3), A2(3), A3(3), B1(3), B2(3), B3(3),
     +  GVEC(3),DELTA(3,3),A01(3),A02(3),A03(3),STRAIN(6)
      CHARACTER*80 TEXT
C     ADJUSTABLE DIMENSIONS:
C.....LATTICE STRUCTURE
      INTEGER TY(NDIM9), NATOM(NDIM9)
      DIMENSION R(49,3,3), IB(48), VTRANS(3,48)
      INTEGER F0(48,NDIM9)
      DIMENSION XKAPA(3,NDIM9),ZKAPA(NDIM9)
      DIMENSION IGLIST(3,NDIM1), G2LIST(NDIM1),
     +  NTAB(NDIM2), G2TAB(NDIM2)
      DOUBLE PRECISION DX(3,NDIM9)
      COMPLEX STRFAC(NTYPMX,NDIM1)
      REAL SATGR(NDIM9,NDIM1), SATGI(NDIM9,NDIM1)
C.....POTENTIALS AND DENSITIES
      INTEGER POTTYP(NTYPMX), POTXC(NTYPMX)
      CHARACTER*80 POTFIL(NTYPMX)
      COMPLEX VIONIC(NDIM1,NSPIN), VGTOT(NDIM1,NSPIN),
     +        ROTOT(NDIM1,NSPIN)
C.....MACHINE TYPE
      COMMON /MACH1/ MACTYP
      COMMON /MACH2/ MACHIN(0:20)
      CHARACTER*7 MACHIN
      CHARACTER*20 LDF(0:3)
C
C.....FILES
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
C.....PHYSICAL AND MATHEMATICAL CONSTANTS
      DOUBLE PRECISION ABOHR,RYEV,RYDERG,PI,SPI
      COMMON /CONST/   ABOHR,RYEV,RYDERG,PI,SPI
C-----------------------------------------------------------------------
C
      CALL DAY (IYEAR,MONTH,IDAY)
      WRITE (IOUT,100) MACHIN(MACTYP)
100   FORMAT (1X,80('*')/
     +  ' K213 - potential setup - ',A,' computer')
      CALL DAYPRT (IOUT,IYEAR,MONTH,IDAY)
      CALL USAGE (NPAGE0,T0)
C
C     INPUT LATTICE STRUCTURE:
      CALL STRUCT (TEXT,A01,A02,A03,A1,A2,A3,NAT,
     +             TY,XKAPA,NATOM,ULA,STRAIN,NDIM9,NTYPMX)
C-----------------------------------------------------------------------
C
C     A POSSIBLE INTRODUCTION OF GHOST ATOMIC SITES REPEATING THE
C     EXISTING ONES:
C
C     NAT IS THE TRUE NUMBER OF ATOMS IN THE CELL, NAT1 INCL. GHOSTS.
      NAT1 = NAT
220   WRITE (IUNIT7,230)
230   FORMAT ('0*** For the virtual crystal approximation only: ***'/
     +' do you wish to introduce some ghost atomic sites which would'/
     +' duplicate the real ones ? (+1=yes; -1=no)')
      READ (INPUT,*) IGHOST
      IF (IGHOST .NE. 1) GOTO 410
      WRITE (IUNIT7,250)
250   FORMAT ('0Which  of the following sites should be duplicated ?'/
     +' (-1=keep single, as it is; +1=duplicate; ',
     +'+N=duplicate N-times.)')
      DO 320 KAPA = 1,NAT
        WRITE (IUNIT7,260) KAPA
260     FORMAT (' KAPA =',I3,' ?')
        READ (INPUT,*) IDUPL
        IF (IDUPL .LE. 0) GOTO 320
        DO 310 I = 1,IDUPL
          NAT1 = NAT1 + 1
C         CHECK THAT NAT1.LE.NDIM9:
          IF (NAT1 .GT. NDIM9) THEN
            WRITE (IOUT,*) 'Too many atoms - start once more'
            GO TO 220
            ENDIF
          TY(NAT1) = TY(KAPA)
          DO 300 IALF = 1,3
            XKAPA(IALF,NAT1) = XKAPA(IALF,KAPA)
300         CONTINUE
            WRITE(IUNIT7,305) NAT1,TY(NAT1),(XKAPA(IALF,NAT1),IALF=1,3)
305         FORMAT(1X,2I5,3F10.5)
310       CONTINUE
320     CONTINUE
C     NAT1 MEANS THE TOTAL NUMBER OF THE EFFECTIVE ATOMIC SITES
      WRITE (IUNIT7,330)
330   FORMAT ('0You may want to modify the definition of atom-type',
     +'s. checkout:'/' (answer: 0=O.K., no modification; INE0=new type-'
     +,'number)')
      DO 360 KAPA = 1,NAT1
        WRITE (IUNIT7,340) KAPA,TY(KAPA)
340     FORMAT (' KAPA =',I3,': TYPE =',I3,' O.K.?')
        READ (INPUT,*) INE0
        IF (INE0 .NE. 0) TY(KAPA) = INE0
360     CONTINUE
      WRITE (IUNIT7,370)
370   FORMAT (/25X,'Summary:')
      WRITE (IUNIT7,375)
375   FORMAT(/12X,'KAPPA TYPE',12X,'X(KAPPA)')
      DO 390 KAPA = 1,NAT1
        WRITE (IUNIT7,305) KAPA,TY(KAPA),(XKAPA(J,KAPA),J=1,3)
        IF (KAPA .EQ. NAT) WRITE (IUNIT7,380)
380     FORMAT (13X,40('-'))
390     CONTINUE
      WRITE (IUNIT7,400)
400   FORMAT ('0Do you want to start the ghost-assignement',
     +' once more or rather to continue'/
     +' program?  (+1=repeat; -1=continue)')
      READ (INPUT,*) IYES
      IF (IYES .EQ. 1) GOTO 220
C-----------------------------------------------------------------------
C
C     CONVERSION OF ALL ESSENTIAL DATA TO DOUBLE PRECISION:
C
410   DO 450 IALF = 1,3
        DA1(IALF) = DBLE(A1(IALF))
        DA2(IALF) = DBLE(A2(IALF))
        DA3(IALF) = DBLE(A3(IALF))
450     CONTINUE
C VOLUME OF THE UNIT CELL
      DVOLUM = DA1(1)*DA2(2)*DA3(3) + DA2(1)*DA3(2)*DA1(3) +
     +         DA3(1)*DA1(2)*DA2(3) - DA1(3)*DA2(2)*DA3(1) -
     +         DA2(3)*DA3(2)*DA1(1) - DA3(3)*DA1(2)*DA2(1)
C RECIPROCAL LATTICE
      DB1(1) = (DA2(2)*DA3(3) - DA2(3)*DA3(2))/DVOLUM
      DB1(2) = (DA2(3)*DA3(1) - DA2(1)*DA3(3))/DVOLUM
      DB1(3) = (DA2(1)*DA3(2) - DA2(2)*DA3(1))/DVOLUM
      DB2(1) = (DA3(2)*DA1(3) - DA3(3)*DA1(2))/DVOLUM
      DB2(2) = (DA3(3)*DA1(1) - DA3(1)*DA1(3))/DVOLUM
      DB2(3) = (DA3(1)*DA1(2) - DA3(2)*DA1(1))/DVOLUM
      DB3(1) = (DA1(2)*DA2(3) - DA1(3)*DA2(2))/DVOLUM
      DB3(2) = (DA1(3)*DA2(1) - DA1(1)*DA2(3))/DVOLUM
      DB3(3) = (DA1(1)*DA2(2) - DA1(2)*DA2(1))/DVOLUM
      DVOLUM = DABS(DVOLUM)
C FOR COMPATIBILITY WITH OTHER (OLD) PARTS OF THE PROGRAM
C IT MAY BE USEFUL TO HAVE B1,2,3, ALSO IN SINGLE PRECISION:
      VOLUM = SNGL(DVOLUM)
      DO 460 IALF = 1,3
        B1(IALF) = SNGL(DB1(IALF))
        B2(IALF) = SNGL(DB2(IALF))
        B3(IALF) = SNGL(DB3(IALF))
460     CONTINUE
      WRITE (IUNIT7,470) B1,B2,B3
470   FORMAT ('0B1',3F10.5,' in units 2PI/ULA'/' B2',3F10.5,
     +' (i.e. A(I)*B(J)=1*DELTA(I,J) )'/' B3',3F10.5)
C     CHECKS:
      WRITE(IUNIT7,*) 'Unit cell volume is ',VOLUM,' times ULA**3'
      DELTA(1,1) = A1(1)*B1(1) + A1(2)*B1(2) + A1(3)*B1(3)
      DELTA(1,2) = A1(1)*B2(1) + A1(2)*B2(2) + A1(3)*B2(3)
      DELTA(1,3) = A1(1)*B3(1) + A1(2)*B3(2) + A1(3)*B3(3)
      DELTA(2,1) = A2(1)*B1(1) + A2(2)*B1(2) + A2(3)*B1(3)
      DELTA(2,2) = A2(1)*B2(1) + A2(2)*B2(2) + A2(3)*B2(3)
      DELTA(2,3) = A2(1)*B3(1) + A2(2)*B3(2) + A2(3)*B3(3)
      DELTA(3,1) = A3(1)*B1(1) + A3(2)*B1(2) + A3(3)*B1(3)
      DELTA(3,2) = A3(1)*B2(1) + A3(2)*B2(2) + A3(3)*B2(3)
      DELTA(3,3) = A3(1)*B3(1) + A3(2)*B3(2) + A3(3)*B3(3)
      CHECK = 0.0
      DO 472 I = 1, 3
      DO 472 J = 1, 3
        IF (I .EQ. J) THEN
          CHECK = CHECK + ABS(DELTA(I,J) - 1.0)
        ELSE
          CHECK = CHECK + ABS(DELTA(I,J))
          ENDIF
472     CONTINUE
        IF (CHECK .GT. 1.0E-9) WRITE (IOUT,474) DELTA
474   FORMAT ('0K213 *** WARNING *** A1,2,3 and B1,2,3 not orthonormal'/
     +        ' the product matrix is:'/3(1X,3E15.8/) )
C
C-----------------------------------------------------------------------
C
C     DETERMINE THE OPERATIONS OF THE CRYSTAL
C
      CALL GROUP1 (IUNIT7,A1,A2,A3,NAT,TY,NDIM9,XKAPA,B1,B2,B3,
     +  IHG,IHC,ISY,LI,NC,IB,VTRANS,F0,R)
C
C-----------------------------------------------------------------------
C
C     GENERATION OF THE RECIPROCAL LATTICE
C
500   WRITE (IUNIT7,510)
510   FORMAT ('0Generation of the reciprocal lattice vectors'/
     +' enter desired plane-wave cutoff in rydbergs'/
     +' and NG1,NG2,NG3, EPSILON (enter 0 for default choice)')
      READ (INPUT,*)   G2RYD,NG1,NG2,NG3,EPSILO
      WRITE (IUNIT7,*) G2RYD,NG1,NG2,NG3,EPSILO
C
C     G-VECTOR CONVERSION FACTOR FROM UNITS (2*PI/ULA)**2 TO RYDBERGS:
      G2MAX = G2RYD / (2.0D0*PI*ABOHR/DBLE(ULA))**2
C
C     CONVENTION: G2MAX IS THE G**2 OF THE PLANE-WAVES INCLUDED (IN 2PI/
C     IN THE HAMILTONIAN.  FOR MATRIX ELEMENTS AND CHARGE DENSITY
C     ONE MUST CONSIDER PLANE-WAVES UP TO 4.? * G2MAX (?=1-5)
C
      IF (G2MAX .LE. 0.0) THEN
C
C       DEDUCE G2MAX FROM DIMENSION NDIM1
C       (NG1,2,3 ARE DEDUCED BY SUBROUTINE RECLAT)
        VOLG   = B1(1)*B2(2)*B3(3) + B2(1)*B3(2)*B1(3) +
     +           B3(1)*B1(2)*B2(3) - B1(3)*B2(2)*B3(1) -
     +           B2(3)*B3(2)*B1(1) - B3(3)*B1(2)*B2(1)
        VOLG   = ABS(VOLG)
        G2MAX  = ( FLOAT(NDIM1) * VOLG * 3.0/(4.0*PI) )**(2./3.)
C       TRUNCATION (0.9 IS A FUDGE FACTOR)
        G2MAX  = FLOAT( IFIX( G2MAX * 0.9 ))
        WRITE (IUNIT7,520) G2MAX * (2.0D0*PI*ABOHR/DBLE(ULA))**2, G2MAX
520     FORMAT (' The plane-wave cutoff deduced is: ',F6.1,' Ryd',
     +    ' (G2MAX =',F10.2,')' )
        ENDIF
C
C     GENERATE RECIPROCAL LATTICE
C
      CALL RECLAT (B1,B2,B3,A1,A2,A3,NG1,NG2,NG3,G2MAX,
     +  NDIM1,NDIM2,EPSILO,IGLIST,G2LIST,NGTOT,NTAB,G2TAB,IWORK,ICOPY)
C
C     SORT IGLIST AND G2TAB INTO STARS OF EQUIVALENT G-VECTORS
C
      CALL GSHELL (IGLIST,NGTOT,NG1,NG2,NG3,A1,A2,A3,B1,B2,B3,
     +  G2TAB,NTAB,NDIM1,NDIM2, R,IB,NC,LI)
C
      DO 555 NSTARS = NDIM2, 1, -1
        IF (NTAB(NSTARS) .GT. 0) GOTO 557
555     CONTINUE
557   WRITE (IUNIT7,560) NGTOT, NSTARS - 1
560   FORMAT(I6,' Reciprocal lattice points generated in',I5,' stars')
C     CHECK ON ERROR SIGNALS
      IF (NGTOT .LE. 0) THEN
        WRITE (IOUT,*) '*** WARNING *** reciprocal lattice generation',
     +  ' ERRORS - TRY AGAIN...'
        GOTO 500
        ENDIF
C
      WRITE (IUNIT8,570)
570   FORMAT ('0Synopsis of reciprocal lattice vectors'/
     +15X,' first    squared'/' shell     occurrence    length',
     +'(in (2PI/ULA)**2)  IGLIST-value'/)
      DO 580 I = 1,NDIM2
        IF (G2TAB(I) .LT. 0.0) GOTO 600
        WRITE (IUNIT8,590) I-1,NTAB(I),G2TAB(I),
     +    (IGLIST(J,NTAB(I)), J = 1, 3)
590     FORMAT (1X,I5,I11,F15.5,15X,3I5)
580     CONTINUE
C     A MORE DETAILED OUTPRINT, IF DESIRED:
600   WRITE (IUNIT7,1100)
      READ (INPUT,*) I1,I2
C     A FOOLPROOF ASSIGNEMENT BETTER THAN ALL CHECKS:
      N1 = MIN0(IABS(I1),IABS(I2))
      N2 = MAX0(IABS(I1),IABS(I2))
      IF (N1*N2 .LE. 0) GOTO 650
      WRITE (IUNIT7,620)
620   FORMAT ('    I    IG1 IG2 IG3',15X,'G(I)',13X,
     +'G(I)**2')
      DO 640 I = N1,N2
        CALL GCODE (IGLIST(1,I),NG1,NG2,NG3,B1,B2,B3,
     +    IG1,IG2,IG3,GVEC)
        G2 = GVEC(1)**2 + GVEC(2)**2 + GVEC(3)**2
        WRITE (IUNIT7,630) I,IG1,IG2,IG3,GVEC,G2
630     FORMAT (1X,I4,3I4,4F10.5)
640     CONTINUE
      GO TO 600
C
650   WRITE (IUNIT7,660)
660   FORMAT ('0Repeat generation of reciprocal lattice vectors or',
     +' continue program ?'/' (+1=repeat; -1=continue)')
      READ (INPUT,*) IYES
      IF (IYES .EQ. 1) GOTO 500
C
C-----------------------------------------------------------------------
C
C     ASSIGNMENT  OF POTENTIALS TO ATOM TYPES:
C
670   WRITE (IUNIT7,680)
680   FORMAT('0Assignment of potential-type to individual atomic ',
     +'species'/' potentials presently available:'/
     +'0code      subroutine')
      WRITE (IUNIT7,690)
690   FORMAT (
     +'   10',8X,'POT10   (jellium)'/
     +'   11',8X,'POT11   (pure coulomb)'/
     +'   13',8X,'POT13   (appelbaum-hamann)'/
     +'   15',8X,'POT15   (berkeley ionic)'/
     +'   16',8X,'POT16   (hamann-schluter-chiang non-local)'/
     +'   17',8X,'POT17   (bell labs 1982-"periodic table")'/
     +'   20',8X,'POT20   (numerically generated)'/
     +1X)
C     COMMENT RELEVANT ONLY IN THE CASE OF GHOST ATOMS:
      IF (NAT .NE. NAT1) WRITE (IUNIT7,710)
710   FORMAT ('0When dealing with duplicated atoms remember that potenti
     +als will only be'/
     +' multiplied by structure factors and then simply added. no provis
     +ion'/
     +' was made for weighting. therefore, do not forget to modify'/
     +' potential-parameters, so as to obtain the combination desired.')
C
C     ASSIGN POTENTIAL PARAMETERS
C
      DO 740 ITYPE = 1,NTYPMX
        POTFIL(ITYPE) = '  '
        POTTYP(ITYPE) = 0
        POTXC (ITYPE) = 0
        DO 740 J = 1,NPARAM
          POTPAR(ITYPE,J) = 0.0D0
740       CONTINUE
C
      DO 870 KAPA = 1,NAT1
        ITYPE = TY(KAPA)
        IF (POTTYP(ITYPE) .NE. 0) GOTO 870
        CALL PARAM (ITYPE,NATOM(KAPA),POTPAR,POTTYP,POTXC,POTFIL,
     +              ULA,NSPIN,NTYPMX,NPARAM)
870     CONTINUE
C
C     CHECK THAT CORRELATION TYPES AGREE BETWEEN THE TYPES OF ATOMS
C
      IEXCOR = 0
      DO 880 ITYPE = 1, NTYPMX
C       THE POTENTIAL MAY NOT CARE WHICH CORRELATION (E.G. COULOMB)
        IF (POTXC(ITYPE) .EQ. 0) GOTO 880
        IF (IEXCOR .EQ. 0) THEN
          IEXCOR = POTXC(ITYPE)
        ELSE
          IF (IEXCOR .NE. POTXC(ITYPE)) THEN
            WRITE (IOUT,875) POTXC
875         FORMAT ('0K213 *** ERROR ***'/
     +        'atomic correlations disagree:',10I3)
            CALL EXIT
            ENDIF
          ENDIF
880     CONTINUE
C
      CALL LDFTYP(LDF)
      IF (IEXCOR .LE. 0) THEN
        IEXCOR = 3
        WRITE (IOUT,*) 'Default chosen:', LDF(IEXCOR)
        ENDIF
      WRITE (IOUT,885) LDF(IEXCOR)
885   FORMAT ('0The exchange-correlation formula is ',A/)
C
C-----------------------------------------------------------------------
C
910   WRITE (IUNIT7,920)
920   FORMAT ('0For linear screening of the ionic potential:'/
     +' enter the exchange factor which then will be used'/
     +' in screening (1. = slater ; 2/3 = kohn-sham)')
      READ (INPUT,*) EXFACT
      WRITE (IUNIT7,*) 'EXFACT =',EXFACT
      IF (EXFACT .LT. 0.01 .OR. EXFACT .GT. 2.0) THEN
        WRITE (IOUT,*) ' Unrealistic exfact, please repeat'
        GOTO 910
        ENDIF
C
C     DETERMINE THE NUMBER OF ELECTRONS IN THE UNIT CELL
C     (CONVENIENTLY CALCULATED BY ALPHAZ)
      CALL ALPHAZ (POTTYP,TY,NAT,POTPAR,NTYPMX,NDIM9,NPARAM,ULA,VOLUM,
     +             ZKAPA,VALPHA)
C
      Z = 0.0
      DO 931 KAPA = 1,NAT
        Z = Z + ZKAPA(KAPA)
931     CONTINUE
      NEL    = NINT(Z)
      NSCREN = NEL
      WRITE(IUNIT7,932) Z,NEL
932   FORMAT('0From the parameters, the total ionic charge = ',F12.6/
     +' the number of electrons per unit cell is set to ',I4/
     +/' for screening the ionic potentials, do you wish to use'/
     +' a different number of electrons (1=yes,-1=no) ?')
      READ (INPUT,*) IYES
      IF (IYES .EQ. 1) THEN
        WRITE (IUNIT7,*) 'Enter the number of screening electrons ?'
        READ (INPUT,*) NSCREN
        WRITE (IUNIT7,*) NSCREN,' electrons for screening'
        ENDIF
      CALL EPS1 (NSCREN,VOLUM,ULA,10.0,EXFACT,FERMIK,EPSQ)
      WVKF = 2.0*PI*FERMIK*ABOHR/ULA
      WRITE (IUNIT7,960) VOLUM,FERMIK,WVKF
960   FORMAT('0Fermi wave-vector for this number of (free) electrons'/
     +' in the volume = ',F15.6,' is'/
     +20X,'k(Fermi) = ',G15.6,' in units (2PI/ULA)'/
     +29X,'= ',G15.6,' in units ABOHR**(-1)')
C
C-----------------------------------------------------------------------
C
C     STRUCTURE FACTORS
C
      CALL SFACT (STRFAC,SATGR,SATGI,TY,XKAPA,DX,IGLIST,
     +  B1,B2,B3,NDIM1,NDIM9,NTYPMX,NAT,NGTOT,NG1,NG2,NG3, .FALSE.)
C
C-----------------------------------------------------------------------
C
C     DETERMINATION OF THE TOTAL IONIC AND SCREENED POTENTIALS.
C     IN BOTH CASES, V(G=0) WILL BE SET EQUAL TO ZERO.
C
      DO 970 ISPIN = 1, NSPIN
        ROTOT( 1,ISPIN) = FLOAT(NEL) / FLOAT(NSPIN)
        VIONIC(1,ISPIN) = 0.0
        VGTOT( 1,ISPIN) = 0.0
        DO 970 IG = 2, NDIM1
          VIONIC(IG,ISPIN) = 0.0
          VGTOT(IG,ISPIN)  = 0.0
          ROTOT(IG,ISPIN)  = 0.0
970       CONTINUE
C
C       WE SHALL MAKE USE OF THE FACT THAT THE ATOMIC (LOCAL)
C       POTENTIALS ARE SPHERICALLY SYMMETRIC AND THAT
C       THEY DEPEND ONLY ON G**2.
C       THE SAME IS TRUE FOR DIELECTRIC CONSTANT AS WELL.
C
      DO 1000 ITYPE = 1,NTYPMX
        IF (POTTYP(ITYPE) .LE. 0) GOTO 1000
        DO 1035 M = 2, NDIM2
          G2 = G2TAB(M)
          IF (G2 .LT. 0.0) GOTO 1000
          NMIN = NTAB(M)
          NMAX = NTAB(M + 1) - 1
C
          INIT = 1
          IF (M .GT. 2) INIT = 0
          L = - 1
          DG2 = DBLE( G2 )
C
          CALL POTVG (POTTYP,POTPAR,POTFIL,NTYPMX,NPARAM,ITYPE,NSPIN,
     +      POTWRK,NQMAX,VOLUM,ULA,INIT,L,DG2,VG,VGP,RHOG)
C
          CALL EPS1 (NSCREN,VOLUM,ULA,G2,EXFACT,FERMIK,EPSQ)
C
          DO 1020 ISPIN = 1, NSPIN
          DO 1020 IG = NMIN, NMAX
            SGVG = STRFAC(ITYPE,IG) * SNGL(VG(ISPIN))
C             IONIC PSEUDOPOTENTIAL
            VIONIC(IG,ISPIN) = VIONIC(IG,ISPIN) + SGVG
C             DO THE LINEAR SCREENING WITH A DIELECTRIC CONSTANT:
            VGTOT(IG,ISPIN)  = VGTOT(IG,ISPIN)  + SGVG / EPSQ
C             ADD UP THE SOLID CHARGE DENSITY
            ROTOT(IG,ISPIN)  = ROTOT(IG,ISPIN) +
     +        STRFAC(ITYPE,IG) * SNGL(RHOG(ISPIN))
1020        CONTINUE
1035    CONTINUE
1000  CONTINUE
C
C-----------------------------------------------------------------------
C
C     PRINTOUT OF SELECTED RESULTS
C
1040  WRITE (IUNIT7,1050)
1050  FORMAT ('0Total  l o c a l  potential for selected G-vectors',
     +' (in rydberg per unit cell) :'/
     +'0   I G(I)**2',10X,'V(G(I)) ionic',T54,'V(G(I)) screened'/
     +T54,'charge ROTOT(G(I))')
      NPRINT = MIN(8,NDIM2)
      DO 1080 M = 1, NPRINT
        IF (G2TAB(M) .LT. 0.0) GOTO 1090
        I = NTAB(M)
        WRITE (IUNIT7,1060) I,G2TAB(M),
     +    VIONIC(I,1), VGTOT(I,1), ROTOT(I,1)
1060    FORMAT (1X,I4,F6.2,4G16.8/43X,2G16.8)
        IF (NSPIN .EQ. 2) WRITE (IUNIT7,1070)
     +    VIONIC(I,2), VGTOT(I,2), ROTOT(I,2)
1070    FORMAT (' Spin up:  ',4G16.8/43X,2G16.8)
1080    CONTINUE
C
C     A MORE DETAILED OUTPRINT, IF DESIRED:
C
1090  WRITE (IUNIT7,1100)
1100  FORMAT ('0Do you wish to see any part of the above table ',
     +  'in some more detail ?'/
     +  ' (if yes, enter i1,i2 (=from,to); if no, enter 0,0.)')
      READ (INPUT,*) I1,I2
      WRITE (IUNIT7,*) I1,I2
C     A FOOLPROOF ASSIGNEMENT BETTER THAN ALL CHECKS:
      N1 = MAX0(1,    MIN0(IABS(I1),IABS(I2)))
      N2 = MIN0(NGTOT,MAX0(IABS(I1),IABS(I2)))
      IF (N1*N2 .LE. 0) GOTO 1150
      WRITE (IUNIT7,1120)
1120  FORMAT ('0   I  IG1 IG2 IG3',13X,'G(I)',11X,'G(I)**2   ',
     +  'V(G(I)) ionic (re,im)'/56X,'V(G(I)) screened')
      DO 1140 I = N1,N2
C       DECODING THE I-TH G-VECTOR:
        CALL GCODE (IGLIST(1,I),NG1,NG2,NG3,B1,B2,B3,
     +    IG1,IG2,IG3,GVEC)
        G2 = GVEC(1)**2 + GVEC(2)**2 + GVEC(3)**2
        WRITE (IUNIT7,1130)
     +  I,IG1,IG2,IG3,(GVEC(J),J=1,3),G2, VIONIC(I,1),VGTOT(I,1)
1130    FORMAT (1X,4I4,4F9.3,2G12.4/53X,2G12.4)
        IF (NSPIN .EQ. 2) WRITE (IUNIT7,1135) VIONIC(I,2),VGTOT(I,2)
1135    FORMAT ( 2 (53X,2G12.4/) )
1140    CONTINUE
      GO TO 1090
C-----------------------------------------------------------------------
C
C A CLEAN-UP OF POTENTIALS:
C
1150  WRITE (IUNIT7,1160)
1160  FORMAT ('0Do you wish a clean-up in the above list of resulting',
     +'potentials ?'/
     +' if yes, enter a small positive epsilon and the program will '/
     +' cancel all components of potentials (re and im separately),'/
     +' which are (in absolute value) smaller than epsilon.'/
     +' if not, enter a negative epsilon')
      READ (INPUT,*) EPSIL1
      WRITE (IOUT,*) 'Epsilon for clean-up of potential =',EPSIL1
      IF (EPSIL1 .LT. 0.0) GOTO 1280
C
      DO 1260 ISPIN = 1, NSPIN
      DO 1260 I = 1, NGTOT
        IF( ABS(REAL(VIONIC(I,ISPIN))) .LE. EPSIL1)
     +    VIONIC(I,ISPIN) = CMPLX(0.0,AIMAG(VIONIC(I,ISPIN)))
        IF( ABS(AIMAG(VIONIC(I,ISPIN))) .LE. EPSIL1)
     +    VIONIC(I,ISPIN) = CMPLX(REAL(VIONIC(I,ISPIN)),0.0)
        IF( ABS(REAL(VGTOT(I,ISPIN))) .LE. EPSIL1)
     +    VGTOT(I,ISPIN)  = CMPLX(0.0,AIMAG(VGTOT(I,ISPIN)))
        IF( ABS(AIMAG(VGTOT(I,ISPIN))) .LE. EPSIL1)
     +    VGTOT(I,ISPIN)  = CMPLX(REAL(VGTOT(I,ISPIN)),0.0)
        IF( ABS(REAL(ROTOT(I,ISPIN))) .LE. EPSIL1)
     +    ROTOT(I,ISPIN)  = CMPLX(0.0,AIMAG(ROTOT(I,ISPIN)))
        IF( ABS(AIMAG(ROTOT(I,ISPIN))) .LE. EPSIL1)
     +    ROTOT(I,ISPIN)  = CMPLX(REAL(ROTOT(I,ISPIN)),0.0)
1260    CONTINUE
      WRITE (IUNIT7,1270)
1270  FORMAT ('0Do you wish to repeat the display of total potential',
     +' for selected G-vectors ?'/' (+1=yes; -1=no)')
      READ (INPUT,*) IYES
      IF (IYES .EQ. 1) GOTO 1040
C-----------------------------------------------------------------------
C
C     SAVING THE RESULTS ON FILE
C
1280  WRITE (IOUT,1290)
1290  FORMAT ('0Do you wish to save the results on file ?'/
     +' (+1=yes; 0=no, stop the program;'/
     +'  -1=no, repeat the assignement of potentials)')
      READ (INPUT,*) IYES
      IF (IYES .EQ. 0) THEN
        GOTO 1360
      ELSE IF (IYES .EQ. -1) THEN
        GOTO 670
        ENDIF
C     HEADING OF THE FILE: CREATED BY PROGRAM K213 AND DATE
      CALL DAY (IYEAR,MONTH,IDAY)
C     HOW LONG IS THE RECORD IN TABLES NTAB,G2TAB,
C     I.E. WHAT IS THE MINIMUM DIMENSION REALLY NEEDED
      DO 1340 I = 1,NDIM2
        MINDIM = I
        IF (G2TAB(I) .LT. 0.0) GOTO 1350
1340    CONTINUE
      MINDIM = NDIM2
1350  NAT    = NAT1
      IRECTP = 11
C
      CALL POTWRI (NPRGR,IVERSN,IDAY,MONTH,IYEAR,TEXT,
     +  IRECTP,NAT,NGTOT,MINDIM,NSPIN,
     +  ULA,A1,A2,A3,B1,B2,B3,STRAIN,TY,XKAPA,NATOM,
     +  G2MAX,NG1,NG2,NG3,EPSILO,NTAB,G2TAB,
     +  IEXCOR,EXFACT,NEL,POTTYP,POTPAR,POTFIL,
     +  IGLIST,VIONIC,VGTOT,ROTOT,RUNPAR,NMIX,HKMIX,XKMIX,FKMIX,CWORK,
     +  NDIM1,NDIM2,NDIM9,NDIM12,NTYPMX,NPARAM,NRUN,
     +  ISTORE)
C
C     STORE CHARGE DENSITY
      IRECTP = 13
C
      CALL POTWRI (NPRGR,IVERSN,IDAY,MONTH,IYEAR,TEXT,
     +  IRECTP,NAT,NGTOT,MINDIM,NSPIN,
     +  ULA,A1,A2,A3,B1,B2,B3,STRAIN,TY,XKAPA,NATOM,
     +  G2MAX,NG1,NG2,NG3,EPSILO,NTAB,G2TAB,
     +  IEXCOR,EXFACT,NEL,POTTYP,POTPAR,POTFIL,
     +  IGLIST,VIONIC,VGTOT,ROTOT,RUNPAR,NMIX,HKMIX,XKMIX,FKMIX,CWORK,
     +  NDIM1,NDIM2,NDIM9,NDIM12,NTYPMX,NPARAM,NRUN,
     +  ISTORE)
C
1360  CALL USAGE (NPAGE1,T1)
      WRITE (IOUT,1370) T1 - T0, NPAGE1 - NPAGE0
1370  FORMAT(' End K213 - ',F12.3,' CPU seconds - ',I10,' page faults')
      RETURN
      END
