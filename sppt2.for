      SUBROUTINE SPPT2 (IQ1,IQ2,IQ3,WVK0,NDIM10,A1,A2,A3,B1,B2,B3,
     +  IHC,INV,NC,IB,R,NTOT,WVKL,LWGHT,LROT,NCBRAV,IBRAV)
C                     Written on September 12-20th, 1979 by K.K.
C                     Modified 26-may-82 by Ole Holm Nielsen
C Generation of special points for an arbitrary lattice,
C following the method Monkhorst,Pack, Phys. Rev. B13 (1976) 5188
C modified by MacDonald, Phys. Rev. B18 (1978) 5897
C The subroutine is written assuming that the points are generated
C in the reciprocal space. if, however, the B1,B2,B3 are replaced
C by A1,A2,A3, then special points in the direct space can be
C produced, as well. (no multiplication by 2pi is then necessary.)
C In the case of nonsymmorphic groups, the application in the
C direct space would probably require a certain caution.
C Subroutines needed: BZDEFI,BZRDUC,INBZ,MESH
C In the cases where the point group of the crystal does not
C contain inversion. the latter may be added if we wish
C (see comment to the switch INV).
C Reduction to the 1st Brillouin zone is done by adding G-vectors
C to find the shortest wave-vector.
C
C The rotations of the Bravais lattice are applied to the
C Monkhorst/Pack mesh in order to find all k-points that are related
C by symmetry. (Ole Holm Nielsen)
C
C Input data:
C      IQ1,IQ2,IQ3 .. parameter Q of Monkhorst and Pack,
C               generalized and different for the 3 directions B1,B2,B3
C      WVK0 ... the 'arbitrary' shift of the whole mesh, denoted k0
C               in MacDonald. WVK0 = 0 corresponds to the original
C               scheme of Monkhorst and Pack.
C               Units: 2pi/(units of length  used in A1, A2, A3),
C               i.e. the same  units as the generated special points.
C      NDIM10 .. Variable dimension of the (output) arrays WVKL,
C               LWGHT,LROT, i.e. space reserved for the special
C               points and accessories.
C               It has to be NDIM10 > NTOT (total number of special
C               points), but the subroutine does not check on this.
C      The following input data may be obtained from the sbrt. GROUP1:
C      B1,B2,B3 .. Reciprocal lattice vectors, not multiplied by
C               any 2pi (in units reciprocal to those of A1,A2,A3)
C      IHC .... code distinguishing between hexagonal and cubic groups
C               IHC=0 stands for hexagonal groups
C               IHC=1 stands for cubic groups
C      INV .... Code indicating whether we wish  to add the inversion
C               to the point group of the crystal or not (in the
C               case that the point group does not contain any).
C               INV=0 means: do not add inversion
C               INV.NE.0 means: add the inversion
C               INV.NE.0 should be the standard choice when SPPT2
C               is used in recciprocal space - in order to make
C               use of the hermiticity of hamiltonian.
C               When used in direct space, the right choice of INV
C               will depend on the nature of the physical problem.
C               in the cases where the inversion is added by the
C               switch INV, the list IB will not be modified but in
C               the output list LROT some of the operations will
C               appear with negative sign; this means that they have
C               to be applied multiplied by inversion.
C      NC ..... Total number of elements in the point group of the crystal
C      IB ..... List of the rotations constituting the point group
C               of the crystal. the numbering is that defined in
C               Worlton and Warren, i.e. the one materialized in the
C               array R (see below)
C               Only the first NC elements of the array IB are meaningful
C      R ...... List of the 3 x 3 rotation matrices
C               (XYZ representation of the O(H) or D(6)H groups)
C               All 48 or 24 matrices are listed.
C      NCBRAV . Total number of elements in RBRAV
C      IBRAV .. List of NCBRAV operations of the bravais lattice
C Output data:
C      NTOT ... Total number of special points
C               If ntot appears negative, this is an error signal
C               which means that the dimension NDIM10 was chosen
C               too small so that the arrays WVKL etc. cannot
C               accomodate all the generated special points.
C               in this case the arrays will be filled up to NDIM10
C               and further generation of new points will be interrupted.
C      WVKL ... List of special points .
C               Cartesian coordinates and not multiplied by 2*pi.
C               Only the first NTOT vectors are meaningful
C               although no 2 points from the list are equivalent
C               by symmetry, this subroutine still has a kind of
C               'beauty defect': the points finally
C               selected are not necessarily situated in a
C               'compact' irreducible brill.zone; they might lie in
C               different irreducible parts of the B.Z. - but they
C               do represent an irreducible set for integration
C               over the entire B.Z.
C     LWGHT ... The list of weights of the corresponding points.
C               These weights are not normalized (just integers)
C      LROT ... For each special point the 'unfolding rotations'
C               are listed. if e.g. the weight of the i-th special
C               point is LWGHT(I), then the rotations with numbers
C               LROT(J,I), J=1,2,...,LWGHT(I) will 'spread' this
C               single point from the irreducible part of B.Z. into
C               several points in an elementary unit cell
C               (parallelopiped) of the reciprocal space.
C               some operation numbers in the list LROT may appear
C               negative, this means that the corresponding rotation
C               has to be applied with inversion (the latter having
C               been artificially added as symmetry operation in
C               case INV.NE.0).no other effort was taken,to renumber
C               the rotations with minus sign or to extend the
C               list of the point-group operations in the list NB.
C         *         *         *         *         *
      DIMENSION WVK0(3),A1(3),A2(3),A3(3),B1(3),B2(3),B3(3),
     +  IB(48),IBRAV(48),R(49,3,3)
      DIMENSION WVKL(3,NDIM10),LWGHT(NDIM10),LROT(48,NDIM10)
      PARAMETER (NRSDIR = 100)
      DIMENSION WVK(3),WVA(3),RSDIR(4,NRSDIR),PROJA(3),PROJB(3)
C     At most NMESH points in the K-point mesh
      PARAMETER (NMESH = 40000)
      DIMENSION INCLUD(NMESH)
      INTEGER YES,NO
      PARAMETER ( YES = 1, NO = 0 )
C
C.....FILES
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
C
      DATA IPLACE /-2/
C-----------------------------------------------------------------------
      NTOT = 0
C
      DO 100 I=1,NDIM10
        LROT(1,I)=1
        DO 100 J=2,48
          LROT(J,I)=0
100       CONTINUE
      DO 110 I = 1,NMESH
        INCLUD(I) = NO
110     CONTINUE
C
C     Define the 1st Brillouin zone
C
      CALL BZDEFI (B1,B2,B3,RSDIR,NRSDIR,NPLANE)
C
C-----------------------------------------------------------------------
C Generation of the mesh (they are not multiplied by 2*pi)
C by the Monkhorst/Pack algorithm, supplemented by all rotations
C
C     Initialize the list of vectors
      IPLACE = - 2
      CALL MESH(WVA,IPLACE)
      IMESH = 0
      DO 330 I1=1,IQ1
      DO 330 I2=1,IQ2
      DO 330 I3=1,IQ3
        UR1=FLOAT(1 + IQ1 - 2*I1)/FLOAT(2*IQ1)
        UR2=FLOAT(1 + IQ2 - 2*I2)/FLOAT(2*IQ2)
        UR3=FLOAT(1 + IQ3 - 2*I3)/FLOAT(2*IQ3)
        DO 170 I=1,3
          WVK(I) = UR1*B1(I) + UR2*B2(I) + UR3*B3(I) + WVK0(I)
170       CONTINUE
C       WRITE(IOUT,175) WVK
C175    FORMAT(' Basic mesh -',3F10.4)
C       Reduce WVK to the 1st Brillouin zone
        CALL BZRDUC(WVK,A1,A2,A3,B1,B2,B3,RSDIR,NRSDIR,NPLANE)
C       Apply all the Bravais lattice operations to WVK
        DO 180 IOP = 1,NCBRAV
          DO 190 I=1,3
            WVA(I) = 0.0
            DO 190 J = 1,3
              WVA(I) = WVA(I) + R(IBRAV(IOP),I,J)*WVK(J)
190           CONTINUE
C
C       Check that WVA is inside the 1 BZ.
        IF (INBZ(WVA,RSDIR,NRSDIR,NPLANE) .EQ. NO) GOTO 450
C
C       Place WVA in list
        IPLACE = 0
        CALL MESH(WVA,IPLACE)
C       If WVA was new (and therefore inserted), IPLACE is the number.
        IF (IPLACE .GT. 0) IMESH = IPLACE
        IF (IPLACE .GT. NMESH) GOTO 470
C
180     CONTINUE
330     CONTINUE
      WRITE (IUNIT7,340) IMESH
340   FORMAT('0The wavevector mesh contains ',I5,' points')
C     WRITE(IOUT,*) 'THE POINTS ARE:'
C     DO 341 I = 1,IMESH
C     CALL MESH (WVA,I)
C     WRITE (IOUT,*) WVA
C341  CONTINUE
C-----------------------------------------------------------------------
C
C     Now figure out if any special point difference (k - k') is an
C     integral multiple of a reciprocal-space vector
C
      IREMOV = 0
      DO 290 I = 1 , (IMESH - 1)
        IPLACE = I
        CALL MESH (WVA,IPLACE)
C       Project WVA onto B1,2,3:
        PROJA(1) = 0.0
        PROJA(2) = 0.0
        PROJA(3) = 0.0
        DO 200 K = 1,3
          PROJA(1) = PROJA(1) + WVA(K)*A1(K)
          PROJA(2) = PROJA(2) + WVA(K)*A2(K)
          PROJA(3) = PROJA(3) + WVA(K)*A3(K)
200       CONTINUE
C       Now loop over all the rest of the mesh points
        DO 280 J = (I + 1), IMESH
          JPLACE = J
          CALL MESH (WVK,JPLACE)
C         Project WVK onto B1,2,3:
          PROJB(1) = 0.0
          PROJB(2) = 0.0
          PROJB(3) = 0.0
          DO 210 K = 1,3
            PROJB(1) = PROJB(1) + WVK(K)*A1(K)
            PROJB(2) = PROJB(2) + WVK(K)*A2(K)
            PROJB(3) = PROJB(3) + WVK(K)*A3(K)
210         CONTINUE
C
C         Check (PROJA - PROJB): is it integral ?
          DO 220 K = 1,3
            DIFF = PROJA(K) - PROJB(K)
            IF ( ABS ( FLOAT(NINT(DIFF)) - DIFF ) .GT. 1.0E-6) GOTO 280
220         CONTINUE
C           DIFF is integral: remove WVK from mesh:
            CALL REMOVE (WVK, JPLACE )
C           If WVK actually removed, increment IREMOV
            IF (JPLACE .GT. 0) IREMOV = IREMOV + 1
280       CONTINUE
290     CONTINUE
C
C     WRITE (IOUT,*) ' NEW MESH IS:'
C     DO 250 I = 1,IMESH
C     IPLACE = I
C     CALL MESH (WVK,IPLACE)
C       WRITE (IOUT,*) WVK,' AT NO.',IPLACE
C       CALL GARBAG (WVK,IPLACE)
C       IF (IPLACE .GT. 0)
C    +    WRITE (IOUT,*) 'GARBAGE...................',WVK
C250   CONTINUE
      IF (IREMOV .GT. 0)
     +  WRITE (IUNIT7,260) IREMOV
260     FORMAT(' Some of these mesh points are related by lattice transl
     +ation vectors'/1X,I6,' of the mesh points removed.'/)
C-----------------------------------------------------------------------
C
C     In the mesh of wavevectors, now search for equivalent points:
C     the inversion (time reversal !) may be used.
C
      DO 350 IWVK = 1,IMESH
        IF (INCLUD(IWVK) .EQ. YES) GOTO 350
C       IWVK has not been encountered previously: new special point,
C       (only if WVK is not a garbage vector, however.)
        INCLUD(IWVK) = YES
        IPLACE = IWVK
        CALL MESH (WVK,IPLACE)
C
C       Find out whether WVK is in the garbage list
        CALL GARBAG (WVK,IGARBG)
        IF (IGARBG .GT. 0) GOTO 350
C
        NTOT = NTOT + 1
        IF (NTOT .GT. NDIM10) THEN
          WRITE (IOUT,*) 'SPPT2: More than max. no. of k-points ',
     +      NDIM10
          CALL EXIT
          ENDIF
C
        DO 360 I = 1,3
          WVKL(I,NTOT) = WVK(I)
360       CONTINUE
        LWGHT(NTOT) = 1
C       Find all the equivalent points
        DO 370 N = 1,NC
C         Rotate:
          DO 361 I = 1,3
            WVA(I) = 0.0
            DO 361 J = 1,3
              WVA(I) = WVA(I) + R(IB(N),I,J)*WVK(J)
361           CONTINUE
        IBSIGN = + 1
363     IPLACE = -1
        CALL MESH(WVA,IPLACE)
C
C       Find out whether WVA is in the garbage list
        CALL GARBAG (WVA,IGARBG)
        IF (IGARBG .GT. 0) GOTO 370
C
C       Was WVA encountered before ?
        IF (INCLUD(IPLACE) .EQ. YES) GOTO 364
        LWGHT(NTOT) = LWGHT(NTOT) + 1
        LROT(LWGHT(NTOT),NTOT) = IB(N)*IBSIGN
        INCLUD(IPLACE) = YES
364     IF (IBSIGN .EQ. -1 .OR. INV .EQ. 0) GOTO 370
C       The case where we also apply the inversion to WVA
C       Repeat the search, but for -WVA
        IBSIGN = - 1
        DO 365 I = 1,3
          WVA(I) = - WVA(I)
365       CONTINUE
        GOTO 363
C
370     CONTINUE
350   CONTINUE
C Total number of special points: NTOT
C before using the list WVKL as wave vectors, they have to be
C multiplied by 2*pi
C the list of weights LWGHT is not normalized
      RETURN
C-----------------------------------------------------------------------
C
C     Error messages
C
450   WRITE (IOUT,410)
410   FORMAT('0SUBROUTINE SPPT2 *** FATAL ERROR ***'/)
      WRITE (IOUT,460) WVA,WVK,IBRAV(IOP)
460   FORMAT('0The vector     ',3F10.4/' generated from ',3F10.4,
     +' in the basic mesh'/
     +' by rotation no. ',I3,' is outside the 1BZ')
      RETURN
C
470   WRITE (IOUT,410)
      WRITE (IOUT,*) 'Mesh size exceeds NMESH=',NMESH
C
      END
      SUBROUTINE MESH (WVK,IPLACE)
C
C     Mesh maintains a list of vectors for placement and/or lookup
C
C     Additional entry points: REMOVE ..... remove vector from list
C                              GARBAG ..... was vector removed ?
C
C     WVK ....... Vector
C     IPLACE .... ON INPUT:  -2 means: initialize  the list (and return)
C                            -1 means: find WVK in the list
C                             0 means: add  WVK to the list
C                            >0 means: return WVK no. IPLACE
C                 On output: the position assigned to WVK
C                            (=0 if WVK is not in the list)
C
      REAL WVK(3)
      PARAMETER ( NHASH = 2000, NMESH = 40000 , NWORDS = 4 )
      PARAMETER ( NLIST = NWORDS * NMESH)
      DIMENSION RLIST(NHASH+NLIST), LIST(NHASH+NLIST)
C
C *** WARNING ***
C
C     The equivalence of LIST and RLIST is required, but it assumes
C     that REALs and INTEGERs correspond to the same number of bits
C     in the computer
      EQUIVALENCE (RLIST(1),LIST(1))
      PARAMETER ( NIL = 0 )
      SAVE ISTORG, IGARB0, IGARBG, RLIST, LIST
C
C.....FILES
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
C-----------------------------------------------------------------------
      IF (IPLACE .GT. -2) GOTO 110
C     Initialization
      DO 100 I = 1,NHASH+NLIST
        LIST(I) = NIL
100     CONTINUE
      ISTORG = NHASH + 1
C     IGARB0 points to a linked list of removed WVK's (the garbage).
      IGARB0 = 0
      RETURN
C
110   IF (IPLACE .GT. 0) GOTO 170
C
C.....THE PARTICULAR HASH FUNCTION USED IN THIS CASE:
C
      RHASH = 0.789*WVK(1) + 0.6810*WVK(2) + 0.5811*WVK(3) + 1.0E-6
      IHASH = IFIX(ABS(RHASH) * FLOAT(NHASH))
      IHASH = MOD(IHASH,NHASH) + 1
C
C     Search for WVK in linked list
      IPOINT = LIST(IHASH)
      DO 120 I = 1,100
C       List exhausted
        IF (IPOINT .EQ. NIL) GOTO 130
C       Compare WVK with this element
        DO 112 J = 1,3
          IF(ABS(WVK(J) - RLIST(IPOINT+J)) .GT. 1.0E-5) GOTO 115
112       CONTINUE
C       WVK located
        GOTO 160
C       Next element of LIST
115     IHASH = IPOINT
        IPOINT = LIST(IHASH)
120     CONTINUE
C     List too long
      WRITE (IOUT,125)
125   FORMAT('0Subroutine MESH *** FATAL ERROR ***',
     +' linked list too long ***'/' choose a better hash-function')
      CALL EXIT
C
C.....Wvk was not found
C
130   IF (IPLACE .EQ. 0) GOTO 140
C     IPLACE=-1: search for WVK unsuccessful
      IPLACE = 0
      RETURN
C     IPLACE=0: add WVK to the list
140   LIST(IHASH) = ISTORG
      LIST(ISTORG) = NIL
      DO 150 I = 1,3
        RLIST(ISTORG+I) = WVK(I)
150     CONTINUE
      ISTORG = ISTORG + NWORDS
      IF (ISTORG .GT. NHASH + NLIST) THEN
        WRITE (IOUT,*) 'SUBROUTINE MESH *** FATAL ERROR ***'
        WRITE (IOUT,155) ISTORG,WVK
155     FORMAT('0ISTORG=',I10,' exceeds dimensions'/
     +         ' WVK = ',3F10.5)
        CALL EXIT
        ENDIF
      IPLACE = (ISTORG - (NHASH + 1) )/NWORDS
      RETURN
C
C.....WVK was found
C
160   IF (IPLACE .EQ. 0) RETURN
C     IPLACE=-1
      IPLACE = (IPOINT - (NHASH + 1))/NWORDS + 1
      RETURN
C
C.....Return a wavevector
C
170   IPOINT = (IPLACE-1)*NWORDS + (NHASH + 1)
      IF (IPOINT .GE. ISTORG) GOTO 190
      DO 180 I = 1,3
        WVK(I) = RLIST(IPOINT + I)
180     CONTINUE
      RETURN
C     ERROR - BEYOND LIST
190   WRITE (IOUT,200) IPLACE
200   FORMAT('0Subroutine MESH *** WARNING ***'/
     +' IPLACE = ',I5,' is beyond the lists - WVK set to 1.0E38'/)
      DO 210 I = 1,3
        WVK(I) = 1.0E38
210     CONTINUE
      RETURN
C-----------------------------------------------------------------------
C
C     ENTRY POINT FOR REMOVING A WAVEVECTOR
C
C     INPUT:   WVK
C     OUTPUT:  IPLACE ..... 1 IF WVK WAS REMOVED
C                           0 IF WVK WAS NOT REMOVED (WVK NOT IN THE
C                             LINKED LISTS)
C
      ENTRY REMOVE (WVK,IPLACE)
C
C.....The particular hash function used in this case:
      RHASH = 0.789*WVK(1) + 0.6810*WVK(2) + 0.5811*WVK(3) + 1.0E-6
      IHASH = IFIX(ABS(RHASH) * FLOAT(NHASH))
      IHASH = MOD(IHASH,NHASH) + 1
C
C     Search for WVK in linked list
C
      IPOINT = LIST(IHASH)
C
      DO 220 I = 1,100
C
C       List exhausted
        IF (IPOINT .EQ. NIL) THEN
C
C         WVK was not found in the mesh:
C
          IPLACE = 0
          RETURN
          ENDIF
C
C       Compare WVK with this element
C
        DO 212 J = 1,3
          IF(ABS(WVK(J) - RLIST(IPOINT+J)) .GT. 1.0E-5) GOTO 215
212       CONTINUE
C
C       WVK located, now remove it from the list:
C
        LIST(IHASH) = LIST(IPOINT)
C       LIST(IHASH) now points to the next element in the LIST,
C       and the present WVK has become garbage.
C
C       Add WVK to the list of garbage:
C
        IF (IGARB0 .EQ. 0) THEN
C       Start up the garbage list:
          IGARB0 = IPOINT
        ELSE
          LIST(IGARBG) = IPOINT
          ENDIF
        IGARBG       = IPOINT
        LIST(IGARBG) = NIL
        IPLACE = 1
C
        RETURN
C
C       Next element of list
215     IHASH = IPOINT
        IPOINT = LIST(IHASH)
220     CONTINUE
C
C     List too long
C
      WRITE (IOUT,125)
      CALL EXIT
      RETURN
C-----------------------------------------------------------------------
C
C     Entry point for checking if a wavevector is in the garbage list
C
C     INPUT:    WVK
C     OUTPUT:   IPLACE  ..... I > 0 IS THE PLACE IN THE GARBAGE LIST
C                             0 IF WVK NOT AMONG THE GARBAGE
C
      ENTRY GARBAG (WVK,IPLACE)
C
C     Search for WVK in linked list
C
C     Point to the garbage list
      IPOINT = IGARB0
      DO 320 I = 1 , NLIST/NWORDS
C
C       List exhausted
        IF (IPOINT .EQ. NIL) THEN
C
C         WVK was not found in the mesh:
C
          IPLACE = 0
          RETURN
          ENDIF
C
C       Compare WVK with this element
C
        DO 312 J = 1,3
          IF(ABS(WVK(J) - RLIST(IPOINT+J)) .GT. 1.0E-5) GOTO 315
312       CONTINUE
C
C       WVK was located in the garbage list
C
        IPLACE = I
        RETURN
C
C       Next element of list
315     IHASH = IPOINT
        IPOINT = LIST(IHASH)
320     CONTINUE
C
C     List too long
C
      WRITE (IOUT,*) ' ***SUBROUTINE GARBAG - garbage list too long***'
      CALL EXIT
      RETURN
C-----------------------------------------------------------------------
      END
      SUBROUTINE BZRDUC(WVK,A1,A2,A3,B1,B2,B3,RSDIR,NRSDIR,NPLANE)
C
C     Reduce WVK to lie entirely within the 1st Brillouin zone
C     by adding B-vectors
C
      REAL WVK(3),A1(3),A2(3),A3(3),B1(3),B2(3),B3(3)
      REAL WVA(3),WB(3),RSDIR(4,NRSDIR)
C     Look around +/- "NZONES" to locate vector
      PARAMETER ( NZONES=4, NNN=2*NZONES+1, NN=NZONES+1 )
      INTEGER YES,NO
      PARAMETER ( YES = 1, NO = 0 )
C
C.....FILES
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
C-----------------------------------------------------------------------
C     WVK already inside 1BZ
      IF (INBZ(WVK,RSDIR,NRSDIR,NPLANE) .EQ. YES) RETURN
C
C     Express WVK in the basis of B1,2,3.
C     This permits an estimate of how far WVK is from the 1BZ.
      WB(1) = WVK(1)*A1(1) + WVK(2)*A1(2) + WVK(3)*A1(3)
      WB(2) = WVK(1)*A2(1) + WVK(2)*A2(2) + WVK(3)*A2(3)
      WB(3) = WVK(1)*A3(1) + WVK(2)*A3(2) + WVK(3)*A3(3)
      NN1 = NINT(WB(1))
      NN2 = NINT(WB(2))
      NN3 = NINT(WB(3))
C
C     Look around the estimated vector for the one truly inside the 1BZ
      DO 200 N1 = 1,NNN
        I1 = NN - N1 - NN1
        DO 200 N2 = 1,NNN
          I2 = NN - N2 - NN2
          DO 200 N3 = 1,NNN
            I3 = NN - N3 - NN3
            DO 190 I = 1,3
                WVA(I) = WVK(I) + FLOAT(I1)*B1(I) + FLOAT(I2)*B2(I) +
     +            FLOAT(I3)*B3(I)
190             CONTINUE
            IF (INBZ(WVA,RSDIR,NRSDIR,NPLANE) .EQ. YES) GOTO 210
200         CONTINUE
C-----------------------------------------------------------------------
C     Fatal error
      WRITE (IOUT,205) WVK
205   FORMAT('0Subroutine BZRDUC *** FATAL ERROR ***'/
     +' wavevector ',3F10.4,' could not be reduced to the 1BZ')
      CALL EXIT
C
C-----------------------------------------------------------------------
C     The reduced vector
210   DO 220 I = 1,3
        WVK(I) = WVA(I)
220     CONTINUE
      RETURN
      END
      FUNCTION INBZ(WVK,RSDIR,NRSDIR,NPLANE)
      REAL WVK(3),RSDIR(4,NRSDIR)
C
C     Is WVK in the 1st Brillouin zone ?
C     Check whether WVK lies inside all the planes that define the 1BZ.
C
      INTEGER YES,NO
      PARAMETER ( YES = 1, NO = 0 )
      PARAMETER ( EPS = 1.0E-6 )
C
C.....FILES
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
C-----------------------------------------------------------------------
      IF (NPLANE .LT. 3) THEN
        WRITE (IOUT,*) 'INBZ *** ERROR *** too few planes: ', NPLANE
        CALL EXIT
        ENDIF
C
      INBZ = NO
      DO 100 N = 1,NPLANE
        PROJCT = (RSDIR(1,N)*WVK(1) + RSDIR(2,N)*WVK(2) +
     +            RSDIR(3,N)*WVK(3) ) / RSDIR(4,N)
C       WVK is outside the BZ
        IF (ABS(PROJCT) .GT. 0.5 + EPS) GOTO 200
100     CONTINUE
      INBZ = YES
200   RETURN
      END
      SUBROUTINE BZDEFI(B1,B2,B3,RSDIR,NRSDIR,NPLANE)
C
C     Find the vectors whose halves define the 1st Brillouin zone
C
C     On output, NPLANE tells how many elements of RSDIR contain
C                normal vectors defining the planes.
C     Method: starting with the parallelopiped spanned by B1,2,3
C     around the origin, vectors inside a sufficiently large sphere
C     are tested to see whether the planes at 1/2*B will
C     further confine the 1BZ.
C     The resulting vectors are not cleaned to avoid redundant planes.
C
      REAL B1(3),B2(3),B3(3),RSDIR(4,NRSDIR)
      REAL BVEC(3)
      PARAMETER ( EPS = 1.0E-6 )
C
C.....FILES
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
      DATA INITLZ /0/
C-----------------------------------------------------------------------
      IF (INITLZ .NE. 0) RETURN
C     Once initialized, we do not repeat the calculation
      INITLZ = 1
      B1LEN = B1(1)**2 + B1(2)**2 + B1(3)**2
      B2LEN = B2(1)**2 + B2(2)**2 + B2(3)**2
      B3LEN = B3(1)**2 + B3(2)**2 + B3(3)**2
C     Lattice containing entirely the Brillouin zone
      BMAX = B1LEN + B2LEN + B3LEN
      NB1 = IFIX( SQRT(BMAX/B1LEN) + 1.0E-6) + 1
      NB2 = IFIX( SQRT(BMAX/B2LEN) + 1.0E-6) + 1
      NB3 = IFIX( SQRT(BMAX/B3LEN) + 1.0E-6) + 1
C     PRINT *,'NB1,2,3 = ',NB1,NB2,NB3
      DO 100 I = 1,NRSDIR
      DO 100 J = 1,4
        RSDIR(J,I) = 0.0
100     CONTINUE
C     1BZ is certainly confined inside the 1/2(B1,B2,B3) parallelopiped
      DO 110 I = 1,3
        RSDIR(I,1) = B1(I)
        RSDIR(I,2) = B2(I)
        RSDIR(I,3) = B3(I)
110     CONTINUE
      RSDIR(4,1) = B1LEN
      RSDIR(4,2) = B2LEN
      RSDIR(4,3) = B3LEN
C     Starting confinement: 3 planes
      NPLANE = 3
C
      DO 150 I1 = NB1, - NB1, -1
      DO 150 I2 = NB2, - NB2, -1
      DO 150 I3 = NB3, - NB3, -1
        DO 120 I = 1,3
          BVEC(I) = FLOAT(I1)*B1(I) + FLOAT(I2)*B2(I) +
     +              FLOAT(I3)*B3(I)
120       CONTINUE
C       Does the plane of 1/2*BVEC narrow down the 1BZ ?
        DO 130 N = 1,NPLANE
          PROJCT = 0.5*(RSDIR(1,N)*BVEC(1) + RSDIR(2,N)*BVEC(2) +
     +                  RSDIR(3,N)*BVEC(3) ) / RSDIR(4,N)
C         1/2*BVEC is outside the BZ - skip this direction
C         The 1.0E-6 takes care of single points touching the BZ, and
C         of the -(plane)
          IF (ABS(PROJCT) .GT. 0.5 - EPS) GOTO 150
130       CONTINUE
        BLEN = BVEC(1)**2 + BVEC(2)**2 + BVEC(3)**2
C       The zero vector:
        IF (BLEN .LT. EPS) GOTO 150
C       1/2*BVEC further confines the 1BZ - include into RSDIR
        NPLANE = NPLANE + 1
C       PRINT *,NPLANE,' PLANE INCLUDED, I1,2,3 = ',I1,I2,I3
        IF (NPLANE .GT. NRSDIR) GOTO 470
        DO 140 I = 1,3
          RSDIR(I,NPLANE) = BVEC(I)
140       CONTINUE
        RSDIR(4,NPLANE) = BLEN
150     CONTINUE
C
C     PRINT INFORMATION
      WRITE (IUNIT7,160) 2*NPLANE,((RSDIR(I,N),I=1,3),N=1,NPLANE)
160   FORMAT('0The 1st Brillouin zone is confined by (at most)',
     +I3,' planes'/
     +' as defined by the +/- halves of the vectors:'/100(1X,3F10.4/) )
      RETURN
C-----------------------------------------------------------------------
C     ERROR MESSAGES
470   WRITE (IOUT,410)
410   FORMAT('0Subroutine BZDEFI *** FATAL ERROR ***')
      WRITE (IOUT,480) NRSDIR
480   FORMAT(' too many planes, NRSDIR = ',I5)
      CALL EXIT
C
      RETURN
      END
