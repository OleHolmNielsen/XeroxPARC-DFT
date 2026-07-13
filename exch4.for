      SUBROUTINE EXCH4 (INIT,NDIM1,NDIM2,NDIM6,NDIM7,NG1,NG2,NG3,NGMAX,
     +  IGLIST,NGTOT,NTAB,ROTOT,N1,N2,N3,NEWTOT,GSPACE,RSPACE,
     +  WORK,VOLUM,ULA,EXFACT,IEXCOR,NSPIN,SUMMU,SUMEXC)
C
C                      written on january 3rd,  1980 - from exch2.
C                      modified 9-dec-81 by OLE HOLM NIELSEN
C                      to clean up the programming.
C                      slightly modified 4-may-82 by OHN
C                      to accomodate general local density-functionals.
C
C      Exchange and correlation: evaluation of Fourier series
C      of a functional FXC(RO) using double Fourier transform
C      Subroutines needed: GCODE3, CFFT,EXCORR
C
C Input data
C      INIT ... switch deciding whether we are performing
C               initialization run or normal run of the subroutine.
C
C               INIT=0 means: initialization
C
C               INIT > 0 means: normal run: 
C               Transform ROG into real space, apply XC functional,
C               transform back to reciprocal space.
C               The initialization run sets up the sine-tables
C               for the CFFT subroutine and performs number of checks.
C               At every modification of N1,N2,N3, the subroutine
C               has to be re-initialized.
C               The normal run is described below in 'numerical method'.
C
C               Special applications:
C               INIT=-1 means: perform only the first part of the normal run,
C               i.e. fourier transform to direct space, as well as the 
C               roundings and checks on RHO(r) which are mentioned below in
C               'numerical method'. Charge density on the real space mesh can 
C               then be extracted from the array RSPACE.
C
C               INIT <= -2 means: continue the normal run which
C               was started with INIT=-1 and interrupted.
C               Two runs with respectively INIT=-1 and
C               INIT <= -2 are equivalent to one run with INIT > 0
C               abs(INIT) = 2 gives mu-XC
C               abs(INIT) = 3 gives  E-XC
C               abs(INIT) = 4 gives  E-XC minus- mu-XC
C               abs(INIT) = 5 gives simply the real-space value.
C
C      NDIM1 ... variable dimension of the arrays IGLIST, ROTOT, GSPACE
C      NDIM2 ... variable dimension of the array NTAB
C      NDIM6 ... variable dimension of the array RSPACE
C               It has to be NDIM6 >= (N1*N2*N3)
C               but it is not required that NDIM6 >=NDIM1.
C      NDIM7 ... variable dimension of the array WORK
C               It has to be NDIM7 > MAX(N1,N2,N3)*12+300+2*MAX(N2,N3)
C      NG1 ... NTAB   .. See the subroutine PSPT4 or RLV4.
C      ROTOT .... Charge-density in the fourier space, i.e.
C               the elements RHO(G) corresponding to G-vectors
C               listed in IGLIST.
C               This fourier series has to have the properties
C               of charge density: its direct-space values
C               have to be real and non-negative.
C               If they are not, the subroutine itself will call
C               them to order: after performing the first
C               transform (from reciprocal to direct space),
C               the possible imaginary parts will be automatically
C               cancelled and possible negative values replaced
C               by zeroes. If these rectifications are greater
C               than what one can expect to be a consequence
C               of roundoff errors, a warning is issued.
C               (a small parameter eps3 controlling this is chosen inside 
C               the subroutine in relation to ROTOT(1).
C               Presently, eps3 = 1.E-3 * CABS(ROTOT(1)) . )
C      N1,N2,N3 ... constants defining the mesh in real and
C               reciprocal space and controlling the double fourier
C               transform (see CFFT).
C               N1,N2,N3 are required to be integer powers of 2,
C               and no one of them should be lt.4.
C               our subroutine places the origin and defines the
C               mesh density as follows:
C               R=J1*A1/N1+J2*A2/N2+J3*A3/N3 , J(I)=0,1,...N(I)-1
C               G=-0.5*(N1*B1+N2*B2+N3*B3)+K1*B1+K2*B2+K3*B3 ,
C                                       K(I)=0,1,....N(I)-1
C               F(R)=X(J1,J2,J3)*(-1)**(J1+J2+J3)
C               F(G)=A(K1,K2,K3) - When A(K1,K2,K3) is rewritten
C                  into F(G). (care is taken to get F(G) which
C                  is defined only on complete shells of
C                  G-vectors.) when, on the contrary, F(G) is
C                  written into A(K1,K2,K3), there is
C               A(K1,K2,K3)=F(G)+sum of all F("equivalent" G),
C                  where "equivalent" G is one which differs from
C                  G by N1*B1 or N2*B2 or N3*B3 or their
C                  combinations (periodicity on the supercell) .
C               Then, a and x are related as in the description
C               of CFFT and the double fourier transform is
C               performed in N1*N2*N3 points.
C
C Output data
C      NEWTOT ... number of G-vectors in IGLIST  which are
C               situated entirely inside the (extended)
C               N1 X N2 X N3 parallelopiped and which still complete a sphere.
C               There will be NEWTOT < NGTOT only if N1,N2,N3
C               are relatively small.
C      GSPACE ... Fourier components of the functional FXC of charge
C               density. elements of the array correspond to
C               G-vectors listed in IGLIST. Only the first NEWTOT
C               elements of GSPACE are meaningful - which
C               correspond to complete shells of G-vectors.
C               The rest is filled by zeros.
C      SUMMU,SUMEXC ... integral in real space of RHO*EXCORR(RHO)
C
C Working fields
C      RSPACE,WORK  ... workspaces
C               used in the subroutine CFFT  (A,S,INV in the
C               description of CFFT).
C
C Numerical method:
C               the input array ROTOT is fourier transformed into
C               the real space and stored in RSPACE. eventual
C               imaginary parts are supposed to result from round-
C               -off errors and cancelled. (if they are >EPS3,
C               a warning is issued.) the real parts are checked
C               to be nonnegative and if they are found negative,
C               they are also cancelled. (if they are <(-EPS3),
C               a warning is issued.) the function FXC is
C               then applied to each mesh point, the    array is
C               transformed back to the fourier space and stored
C               in GSPACE.
C
C Variable dimension:
      DIMENSION IGLIST(3,NDIM1), NTAB(NDIM2),
     +          WORK(4*NDIM7)
      COMPLEX ROTOT(NDIM1,NSPIN), GSPACE(NDIM1,NSPIN),
     +        RSPACE(NDIM6,NSPIN)
      DOUBLE PRECISION SUMMU, SUMEXC
      REAL RRHO(2), EXC(2), MUXC(2)
      INTEGER MMM(3),KLIST(12)
      SAVE MMM, MXNI
CDIR$ INT24 IG1,IG2,IG3, N1,N2,N3, K1,K2,K3, J1,J2,J3, INDEX,ISPIN,IG
C
C.....FILES
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
      COMMON /WARN/ IWARN
C
C Statement function for evaluation of a location
C of (J1,J2,J3) in the one-dimensional array RSPACE.
C J1,J2,J3 are those indices which run from 0 to N-1:
      LOC(J1,J2,J3) = 1 + J1 + J2*N1 + J3*N1*N2
C Note that this mapping of 3 indices onto one will be
C Singlevalued only if the J1,J2,J3 stay in the range (0,N-1).
C
C     True mathematical modulus-function (positive even if I<0)
      MODN(I,N) = MOD( N + MOD(I,N), N)
      DATA KLIST/0,0,0,1,0,0,0,1,0,0,0,1/
C-----------------------------------------------------------------------
C
      IF (INIT .NE. 0) GOTO 300
C-----------------------------------------------------------------------
C
C     Initialization
C
C     Determine M1,M2,M3: M1 = LOG(BASE=2)(N1),  etc. (see CFFT)
      CONVLG = 1.0/ALOG(2.0)
      M1 = NINT(ALOG(FLOAT(N1))*CONVLG)
      M2 = NINT(ALOG(FLOAT(N2))*CONVLG)
      M3 = NINT(ALOG(FLOAT(N3))*CONVLG)
C   Check that the M's are in the allowed limits:
      MXMI   = MAX0(M1,M2,M3)
      IF (MXMI .LT. 3) GOTO 900
C     Check the dimensions
      IF (NDIM6 .LT. N1*N2*N3) GOTO 950
      MXNI = MAX0(N1,N2,N3)
      IF (4*NDIM7 .LT. MXNI) GOTO 980
      MMM(1) = M1
      MMM(2) = M2
      MMM(3) = M3
C   Set up the sine tables:
      IFSET = 0
      CALL CFFT (MMM,RSPACE,N1*N2*N3 * 2,
     +  WORK(1),WORK(1+MXNI),MXNI,IFSET,IFERR)
      IF (IFERR .NE. 0) GOTO 1000
C   Find the greatest sphere which fits entirely inside the
C   N1 X N2 X N3 parallelopiped:
      NEWTOT = NGMAX
C     If IOVF < 0, no diagnostics will be printed by EXCH4
C     (they really just mean that sphere and box are different)
      IOVF = - 1
      DO 120 IG = 2, NGMAX
        K1 = IGLIST(1,IG) + N1/2
        K2 = IGLIST(2,IG) + N2/2
        K3 = IGLIST(3,IG) + N3/2
        IF (K1.LT.0 .OR. K1.GT.N1 .OR. K2.LT.0 .OR. K2.GT.N2 .OR.
     +      K3.LT.0 .OR. K3.GT.N3) GO TO 130
C       We are thus allowing for the "extended" parallelopiped.
120     CONTINUE
C     No overflow:
      IF (IOVF .EQ. 0) THEN
        GOTO 200
      ELSE
        GOTO 260
        ENDIF
C
C     The overflow case:
C     Locating the "IG" corresponding to the first overflowing G(I):
130   DO 140 J = 1 , (NDIM2 - 1)
        IIITAB = (IG - NTAB(J)) * (IG - NTAB(J + 1))
        IF (IIITAB .GT. 0) GOTO 140
C       The case IG=NTAB(J+1):
        IF (IIITAB .EQ. 0) NEWTOT = NTAB(J + 1) - 1
C       The case when IG is between NTAB(J) and NTAB(J+1):
C       (IG=NGMAX excluded: this would follow another routing)
        IF (IIITAB .LT. 0) NEWTOT = NTAB(J) - 1
        GOTO 150
C
140     CONTINUE
C
150   IF (IOVF .LT. 0) GOTO 260
C
      WRITE (IUNIT7,160) N1,N2,N3
160   FORMAT (' There are G-vectors in IGLIST which are too long for'/
     +  ' the ',I4,' X ',I4,' X ',I4,
     +  ' parallelopiped centered at the origin'/
     +  ' of reciprocal space. e.g. the shortest of them:')
      WRITE (IUNIT7,170) IG,IGLIST(1,IG),IGLIST(2,IG),IGLIST(3,IG)
170   FORMAT ('      I    I1    I2    I3',/,1X,4I6)
      WRITE (IUNIT7,180) NEWTOT
180   FORMAT (' The corresponding elements of RO(G) will be'/
     +' folded into parallelopiped (no truncation takes place),'/
     +' but the array gspace will be truncated after newtot =',I5,
     +' elements,'/
     +'  in order to be given only on complete shells of G-vectors')
C
C     IWARN = 1
C
C     Check whether all the points of the N1 X N2 X N3 parallelopiped
C     are among the G-vectors of IGLIST. It is enough to test 4 examples
C
200   IOVF = 0
      DO 240 I = 1,4
        K1 = (N1 - 1)*KLIST(3*I - 2)
        K2 = (N2 - 1)*KLIST(3*I - 1)
        K3 = (N3 - 1)*KLIST(3*I)
        I1 = K1 - N1/2
        I2 = K2 - N2/2
        I3 = K3 - N3/2
        DO 210 JJ = NGTOT, 2, -1
          IF (IGLIST(1,JJ) .NE. I1) GOTO 210
          IF (IGLIST(2,JJ) .NE. I2) GOTO 210
          IF (IGLIST(3,JJ) .EQ. I3) GOTO 240
210       CONTINUE
        IF (IOVF .LE. 0) WRITE (IUNIT8,220)
220     FORMAT (' Subroutine EXCH4 info: fourier transform mesh ',
     +  '(N1 X N2 X N3 parallelopiped)'/
     +  ' is larger than the sphere of G-vectors in IGLIST.',
     +  ' a few examples:'/
     +  '    K1   K2   K3    I1   I2   I3')
        IOVF = 1
        WRITE (IUNIT8,230) K1,K2,K3,I1,I2,I3
230     FORMAT (1X,6I5)
240     CONTINUE
C
      IF (IOVF .NE. 0) WRITE (IUNIT8,250)
250   FORMAT (' The values of A(K1,K2,K3) at these points will ',
     +'be neglected'/
     +' when rewriting A(K1,K2,K3) into GSPACE(G=I1*B1+I2*B2+I3*B3).')
C
C     End of the initialization run of the subroutine
C
260   IF (INIT .EQ. 0) RETURN
C-----------------------------------------------------------------------
C
C     The proper run of the subroutine
C
C     Check that N1,N2,N3 were not changed after initialization
300   IF (2**MMM(1) .NE. N1) GOTO 1020
      IF (2**MMM(2) .NE. N2) GOTO 1020
      IF (2**MMM(3) .NE. N3) GOTO 1020
C
C     Assume that we just did a call with INIT=-1, and now we are going
C     to continue the Fourier transform to calculate the exchange-
C     correlation terms:
C     Count negative or imaginary charge points:
      IF (INIT .LE. -2) GOTO 340
C-----------------------------------------------------------------------
C
C     INIT = -1 or positive
C
      DO 330 ISPIN = 1, NSPIN
C
C     Writing the coefficients RO(G) into the array RSPACE
C     in the order required by CFFT
      DO 310 I = 1, NDIM6
        RSPACE(I,ISPIN) = 0.0
310     CONTINUE
      DO 320 IG = 1, NEWTOT
C       Reducing all the indices into the limits 0, N(I)-1:
C       (MODN: statement function for "true" MOD(I,N) )
        K1 = MODN(IGLIST(1,IG) + N1/2, N1)
        K2 = MODN(IGLIST(2,IG) + N2/2, N2)
        K3 = MODN(IGLIST(3,IG) + N3/2, N3)
        INDEX         = LOC(K1,K2,K3)
        RSPACE(INDEX,ISPIN) = RSPACE(INDEX,ISPIN) + ROTOT(IG,ISPIN)
320     CONTINUE
C
C     Inverse Fourier transform in order to get the real space values:
C
      IFSET = - 2
      CALL CFFT (MMM,RSPACE(1,ISPIN),N1*N2*N3 * 2,
     +  WORK(1),WORK(1+MXNI),MXNI,IFSET,IFERR)
      IF (IFERR .NE. 0) GOTO 1040
C
330   CONTINUE
C
C     Now, RSPACE(LOC(J1,J2,J3))*(-1)**(J1+J2+J3) are the values
C     of the function on the real-space mesh
C
      IF (INIT .EQ. -1) RETURN
C-----------------------------------------------------------------------
C
C     Density functional:
C
340   NUMNEG = 0
      WORST  = 0.0
      SUMMU  = 0.0D0
      SUMEXC = 0.0D0
      SUMRO  = 0.0
      ROMAX  = - 1.0E36
      ROMIN  =   1.0E36
C     Small parameter used for fighting round-off errors.
      EPS3   = 0.001 * ABS(REAL(ROTOT(1,1)))
C
      DO 390 J1 = 0, (N1 - 1)
      DO 390 J2 = 0, (N2 - 1)
      DO 390 J3 = 0, (N3 - 1)
        INDEX = LOC(J1,J2,J3)
C       EVEN OR ODD?
        IF (MODN(J1 + J2 + J3 , 2 ) .EQ. 0) THEN
          SIGNDX =   1.0
        ELSE
          SIGNDX = - 1.0
          ENDIF
C
C       Block using specific properties of charge density - beginning
C
        DO 370 ISPIN = 1, NSPIN
C       Is the density real ?
        DENSI = ABS( AIMAG( RSPACE(INDEX,ISPIN)))
        IF (DENSI .GT. EPS3) THEN
C         Warning with fix:
          NUMNEG = NUMNEG + 1
          IF (DENSI .GT. WORST) WORST = DENSI
C         Only 3 print-outs
          IF (NUMNEG .LE. 3)
     +      WRITE (IUNIT7,360) J1,J2,J3,RSPACE(INDEX,ISPIN)
360       FORMAT(' EXCH4 *** WARNING(+FIX) *** RSPACE(',
     +      3I4,') = ',2E15.5)
          ENDIF
        RRHO(ISPIN)         = REAL(RSPACE(INDEX,ISPIN))
        RSPACE(INDEX,ISPIN) = RRHO(ISPIN)
C       Is the density positive?
        RRHO(ISPIN)   = RRHO(ISPIN) * SIGNDX
        IF (RRHO(ISPIN) .LT. 0.0) THEN
          IF (RRHO(ISPIN) .LT. - EPS3) THEN
C           Warning with fix:
            NUMNEG = NUMNEG + 1
            IF (ABS(RRHO(ISPIN)) .GT. WORST) WORST = ABS(RRHO(ISPIN))
            IF (NUMNEG .LE. 3) WRITE (IUNIT7,360) J1,J2,J3,RRHO(ISPIN)
            ENDIF
          RSPACE(INDEX,ISPIN) = (0.0,0.0)
          RRHO(ISPIN) = 0.0
          ENDIF
370     CONTINUE
C
C-----------------------------------------------------------------------
C
C     Here the local-density-functional is evaluated:
C     EXCORR calculates E-XC and MU-XC
C
      CALL EXCORR (RRHO,IEXCOR,NSPIN,VOLUM,ULA,EXFACT,EXC,MUXC)
C
      RHO = 0.0
      DO 385 ISPIN = 1, NSPIN
        IF (IABS(INIT) .EQ. 2) THEN
C         Mu-XC
          RSPACE(INDEX,ISPIN) = SIGNDX * MUXC(ISPIN)
        ELSE IF (IABS(INIT) .EQ. 3) THEN
C         E-XC
          RSPACE(INDEX,ISPIN) = SIGNDX * EXC(ISPIN)
        ELSE IF (IABS(INIT) .EQ. 4) THEN
C         E-XC - Mu-XC
          RSPACE(INDEX,ISPIN) = SIGNDX * (EXC(ISPIN) - MUXC(ISPIN))
        ELSE IF (IABS(INIT) .EQ. 5) THEN
C         The original real-space value times the (-1)**(J1+J2+J3)
          RSPACE(INDEX,ISPIN) = SIGNDX * RSPACE(INDEX,ISPIN)
        ELSE
C         Error:
          RSPACE(INDEX,ISPIN) = 0.0
          ENDIF
C       Real-space integral of XC potential and energy:
        SUMMU  = SUMMU  + DPROD(RRHO(ISPIN), MUXC(ISPIN))
        SUMEXC = SUMEXC + DPROD(RRHO(ISPIN), EXC(ISPIN))
        RHO    = RHO + RRHO(ISPIN)
385   CONTINUE
      SUMRO = SUMRO  + RHO
      ROMAX = AMAX1(ROMAX, RHO)
      ROMIN = AMIN1(ROMIN, RHO)
C
390   CONTINUE
C
C-----------------------------------------------------------------------
C
      SUMMU  = SUMMU  / DBLE(N1*N2*N3)
      SUMEXC = SUMEXC / DBLE(N1*N2*N3)
      SUMRO  = SUMRO  / DBLE(N1*N2*N3)
      WRITE (IUNIT7,395) SUMRO, ROMAX, ROMIN
395   FORMAT (' EXCH4: Average charge density = ',G20.12/
     +        '        Max,min = ',2G15.7)
C
C     Summary of errors in charge density:
      IF (NUMNEG .GT. 0) WRITE(IUNIT7,400) NUMNEG,WORST
400   FORMAT(' EXCH4 *** summary of strongly imaginary or negative densi
     +ties:'/1X,I6,' fixes made, worst case was = ',G10.4)
C
C     Now, the array is ready for fourier transform back
C     to the reciprocal space
C
      DO 430 ISPIN = 1, NSPIN
C
      IFSET = + 2
      CALL CFFT (MMM,RSPACE(1,ISPIN),N1*N2*N3 * 2,
     +  WORK(1),WORK(1+MXNI),MXNI,IFSET,IFERR)
C
      IF (IFERR .NE. 0) GOTO 1040
C   Now, we still have to sort the elements of RSPACE into the array GSPACE
C   we only need to consider NEWTOT elements in GSPACE - but note that
C   we extended our N1 X N2 X N3 parallelopiped by allowing
C   the values K(I)=N(I), so as to have it symmetric. this brings
C   an additional complication: we may be obliged to share a given
C   A(K1,K2,K3)-value between two elements of GSPACE(I).
      DO 410 IG = 1, NEWTOT
        GSPACE(IG,ISPIN) = 0.0
C       Reducing all the indices into the limits 0, N(I)-1:
        K1 = MODN(IGLIST(1,IG) + N1/2,N1)
        K2 = MODN(IGLIST(2,IG) + N2/2,N2)
        K3 = MODN(IGLIST(3,IG) + N3/2,N3)
        IF (K1.EQ.0  .OR. K2.EQ.0  .OR. K3.EQ.0) THEN
C         This is a point on the surface of the extended parallelopiped.
C         Let us examine -G(I), a point of the opposite surface
C         Reducing all the indices into the limits 0, N(I)-1:
          K1P = MODN( -IGLIST(1,IG) + N1/2,N1)
          K2P = MODN( -IGLIST(2,IG) + N2/2,N2)
          K3P = MODN( -IGLIST(3,IG) + N3/2,N3)
          INDEX = LOC(K1P,K2P,K3P)
C         Are both points related by periodicity of supercell
C         or are they completely different?
          IF (K1 .EQ. K1P .AND. K2 .EQ. K2P .AND. K3 .EQ. K3P) THEN
            GSPACE(IG,ISPIN) = 0.5 * RSPACE(INDEX,ISPIN)
          ELSE
            GSPACE(IG,ISPIN) =       RSPACE(INDEX,ISPIN)
            ENDIF
C
        ELSE
C         Here,  we have the simpler case, when no sharing is necessary
          INDEX = LOC(K1,K2,K3)
          GSPACE(IG,ISPIN) = RSPACE(INDEX,ISPIN)
          ENDIF
C
410     CONTINUE
C
C     Filling the rest by zeros:
      IF (NEWTOT .LT. NDIM1) THEN
        DO 420 IG = (NEWTOT + 1), NDIM1
          GSPACE(IG,ISPIN) = 0.0
420       CONTINUE
        ENDIF
C
430   CONTINUE
C
      RETURN
C
C There exists a better way how to sort RSPACE into GSPACE,
C without truncation at NEWTOT (or, at least, with a less
C severe truncation): by using the symmetry properties of RO(G).
C      Once we rewrite RSPACE(K1,K2,K3) into GSPACE(G) for one particular
C G, we know immediately all other GSPACE(G) of the same shell:
C RO(R G) = RO(G)*EXP(-I*R G*V(R)). This offers a possibility
C to complete those shells which only overlap (but do not
C coincide) with the parallelopiped.
C      There is, of course a (programming-) complication. all
C terms have to be attributed suitable weighting factors
C so as to produce again the original values of RSPACE, when all
C RO(G) are "folded" back to parallelopiped. a real difficulty
C is that this "folding" can mix various shells.
C      Hopefuly, we will be able to     program this in future.
C Another possible improvement: the truncated GSPACE(I) for
C I>NEWTOT could be approximated by some assymptotic formula.
C (then, however, the truncated series would not longer be
C the best short fourier series approximation to the real-space
C charge density.)
C-----------------------------------------------------------------------
C
C     Error messages
C
900   WRITE (IOUT,910)
910   FORMAT ('0Subroutine EXCH4 - abnormal end')
      WRITE (IOUT,920) N1,N2,N3,M1,M2,M3
920   FORMAT ('    N1   N2   N3   M1   M2   M3'/1X,6I5)
      WRITE (IOUT,930)
930   FORMAT (' The greatest M(I) is smaller than 3, which is'/
     +' unacceptable for the subroutine CFFT')
      CALL EXIT
      RETURN
C
950   WRITE (IOUT,910)
      WRITE (IOUT,960) NDIM6,N1,N2,N3
960   FORMAT (' NDIM6   N1   N2   N3'/1X,4I5)
      WRITE (IOUT,970)
970   FORMAT (' Dimension NDIM6 is not sufficient',
     +' to accommodate N1*N2*N3 complex numbers'/
     +' increase NDIM6 or  decrease N1,N2,N3')
      CALL EXIT
      RETURN
C
980   WRITE (IOUT,910)
      WRITE (IOUT,990) NDIM7,N1,N2,N3
990   FORMAT (' Dimension NDIM7 =',I5,' is not sufficient'/
     +' for running the subroutine CFFT with'/
     +' N1=',I5,' N2=',I5,' N3=',I5)
      CALL EXIT
      RETURN
C
1000  WRITE (IOUT,910)
      WRITE (IOUT,1010) IFERR
1010  FORMAT (' The initialization run of subroutine CFFT ended'/
     +  ' with error signal =',I3)
      CALL EXIT
      RETURN
C
1020  WRITE (IOUT,910)
      WRITE (IOUT,*) 'Contradiction between N and M:'
      WRITE (IOUT,920) N1,N2,N3,M1,M2,M3
      WRITE (IOUT,1030)
1030  FORMAT (' It is very likely that you changed N1,N2,N3'/
     +  ' and forgot to re-initialize the subroutine')
      CALL EXIT
      RETURN
C
1040  WRITE (IOUT,910)
      WRITE (IOUT,1050) IFSET,IFERR
1050  FORMAT (' The run of the subroutine CFFT ended with'/
     +  ' error signal: IFSET = ',I3,' IFERR = ',I3)
      CALL EXIT
      RETURN
      END
