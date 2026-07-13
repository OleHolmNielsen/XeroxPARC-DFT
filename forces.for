      SUBROUTINE FORLOC (TY,POTTYP,NG1,NG2,NG3,B1,B2,B3,
     +  POTPAR,POTFIL,VOLUM,ULA,NAT,NTAB,G2TAB,IGLIST,SATGR,SATGI,
     +  NDIM1,NDIM2,NDIM9,NTYPMX,NPARAM,LREAL,
     +  NSPIN,POTWRK,NQMAX,
     +  FORCWK,ROTOT,FORCE)
C
C     Local potential contribution to Hellmann-Feynman forces
C
C     Written by O. H. Nielsen (1982)
C     Subroutines needed: POTVG
C
      LOGICAL LREAL
      DOUBLE PRECISION G2,VG(2),VGP(2),RHOG(2),FORCWK(3,NDIM9),
     +  FORCE(3,NDIM9)
      INTEGER TY(NDIM9),POTTYP(NTYPMX)
      DOUBLE PRECISION POTPAR(NTYPMX,NPARAM)
      DOUBLE PRECISION POTWRK(*)
      CHARACTER*80 POTFIL(NTYPMX)
      DIMENSION B1(3),B2(3),B3(3),
     +  NTAB(NDIM2),G2TAB(NDIM2),IGLIST(3,NDIM1),
     +  SATGR(NDIM9,NDIM1),SATGI(NDIM9,NDIM1)
      COMPLEX ROTOT(NDIM1,NSPIN)
      DOUBLE PRECISION VSR,ELOCAL
C
C.....FILES
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
C.....PHYSICAL AND MATHEMATICAL CONSTANTS
      DOUBLE PRECISION ABOHR,RYEV,RYDERG,PI,SPI
      COMMON /CONST/   ABOHR,RYEV,RYDERG,PI,SPI
C
C-----------------------------------------------------------------------
C
      DO 100 I = 1,NDIM9
        FORCWK(1,I) = 0.0D0
        FORCWK(2,I) = 0.0D0
100     FORCWK(3,I) = 0.0D0
      ELOCAL = 0.0D0
C
C-----------------------------------------------------------------------
C
      DO 300 ITYPE = 1,NTYPMX
C
      IF (POTTYP(ITYPE) .LE. 0) GOTO 300
C
      DO 250 M = 2,NDIM2
C
        G2LGTH = G2TAB(M)
        IF (G2LGTH .LT. 0.0) GOTO 300
        NMIN = NTAB(M)
        NMAX = NTAB(M + 1) - 1
C
C       Decision about the potential type:
C
        INIT = 1
        IF (M .GT. 2) INIT = 0
        L = -1
        G2 = DBLE(G2LGTH)
C
        CALL POTVG (POTTYP,POTPAR,POTFIL,NTYPMX,NPARAM,ITYPE,
     +    NSPIN,POTWRK,NQMAX,VOLUM,ULA,INIT,L,G2,VG,VGP,RHOG)
C
C       Contributions to force from the M-th shell of G-vectors:
C
        DO 220 ISPIN = 1, NSPIN
        DO 220 I = NMIN, NMAX
C
          DO 200 KAPA = 1,NAT
            IF (TY(KAPA) .EQ. ITYPE) THEN
C             IMAG( V(G) * S(G) * RHO(G)* )
              IF (LREAL) THEN
                VSR = VG(ISPIN) *
     +            DPROD( SATGI(KAPA,I) , REAL (ROTOT(I,ISPIN)))
              ELSE
                VSR = VG(ISPIN) * (
     +            DPROD( SATGI(KAPA,I) , REAL (ROTOT(I,ISPIN))) -
     +            DPROD( SATGR(KAPA,I) , AIMAG(ROTOT(I,ISPIN))) )
                ENDIF
c             Would be useful, if SATGR were calculated for LREAL=.TRUE.
C             Check on the local-potential energy:
C              ELOCAL = ELOCAL + VG(ISPIN) * (
C     +          DPROD( SATGR(KAPA,I) , REAL (ROTOT(I,ISPIN))) +
C     +          DPROD( SATGI(KAPA,I) , AIMAG(ROTOT(I,ISPIN))) )
              FORCWK(1,KAPA) = FORCWK(1,KAPA) + DBLE(IGLIST(1,I))*VSR
              FORCWK(2,KAPA) = FORCWK(2,KAPA) + DBLE(IGLIST(2,I))*VSR
              FORCWK(3,KAPA) = FORCWK(3,KAPA) + DBLE(IGLIST(3,I))*VSR
              ENDIF
200         CONTINUE
C
220       CONTINUE
C
250     CONTINUE
C
300     CONTINUE
C
C-----------------------------------------------------------------------
C
C     Conversion to cartesian basis.
C     Multiply by - 2*pi/ula (note that we want the -gradient).
C     Unit conversion taken care of in main program.
C
      DO 410 KAPA = 1,NAT
      DO 410 J = 1,3
        FORCE(J,KAPA) = - 2.0D0 * PI / DBLE(ULA) *
     +    ( FORCWK(1,KAPA) * DBLE(B1(J)) +
     +      FORCWK(2,KAPA) * DBLE(B2(J)) +
     +      FORCWK(3,KAPA) * DBLE(B3(J)) )
410     CONTINUE
C
C      WRITE (IUNIT8,420) ELOCAL*RYEV
C420   FORMAT(' Subroutine FORLOC: energy term (D) =',F15.5,' eV')
C
      RETURN
      END
      SUBROUTINE FORSYM (TY,XKAPA,A1,A2,A3,NAT,ISY,NC,IB,V,F0,R,
     +  FORCWK,NDIM9,FORCE)
C
C     Apply point group operations to the forces calculated by
C     summing over special points.
C     Written 12-may-82 by Ole Holm Nielsen
C     Input: see subroutine ROSYM2.
C     FORCE ..... the forces to be symmetrized. on return force
C                 contains the correct forces
C
      INTEGER TY(NDIM9),IB(48),F0(48,NDIM9)
      REAL XKAPA(3,NDIM9),A1(3),A2(3),A3(3),V(3,48),R(49,3,3)
      DOUBLE PRECISION FORCE(3,NDIM9), FORCWK(3,NDIM9)
C
C.....FILES
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
C
C.....WARNINGS
      COMMON /WARN/ IWARN
C-----------------------------------------------------------------------
      IZERO = 0
C     Keep force in working array
      DO 100 KAPA = 1,NAT
      DO 100 J = 1,3
        IF (DABS(FORCE(J,KAPA)) .LT. 1.0E-10) IZERO = IZERO + 1
C       Divide by the number of group operations
        FORCWK(J,KAPA) = FORCE(J,KAPA) / DBLE(NC)
        FORCE(J,KAPA) = 0.0D0
100     CONTINUE
C
      IF (IZERO .GT. 0) THEN
        WRITE(IUNIT7,110) IZERO
110     FORMAT ('0Subroutine FORSYM - warning - ',I3,' forces are zero')
        IWARN = 1
        ENDIF
C     Loop over atoms
      DO 190 KAPA = 1,NAT
C       Sum over rotations
        DO 190 IROT = 1,NC
C         Note that from K290, the arrays F0 and V are arranged
C         differently from the rotation matrices.
          KAPAP = F0(IROT,KAPA)
C         Rotate the force vector by R**(-1)
          IC = IB(IROT)
          DO 190 I = 1,3
          DO 190 J = 1,3
            FORCE(I,KAPA) = FORCE(I,KAPA) + R(IC,J,I)*FORCWK(J,KAPAP)
190         CONTINUE
      RETURN
      END
