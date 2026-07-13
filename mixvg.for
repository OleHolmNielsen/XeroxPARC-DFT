      SUBROUTINE MIXVG (IDMAT,UPDATE,VCL,ALPHAS,ALPHAB,GSQUBS,VIN,VOUT,
     +  HK,CHI,INV,XK,FK,NDIM2,NDIM5,NDIM8,NDIM10,NDIM12,NDIM13,
     +  ISPIN,NSPIN,
     +  IGLIST,LIST,NGTOT,NTAB,G2TAB,NG1,NG2,NG3,A1,A2,A3,
     +  B1,B2,B3,VOLUM,ULA,NEL,LI,NC,IB,R,RB,RTABLE,RLIST,VTRANS,
     +  YK,SK,SKHK,SMHY,HKOLD,
     +  LISTAB,CEV,LWGHT,MINBDS,EE1,EE2,EF,IEVFIL,
     +  NUMKPT,LSEMI,LREAL,NHK,NPWCHI,IDIFF,IPG,INDEX,IFILEH)
C
C     MIXING THE OLD AND NEW POTENTIALS
C
C     WRITTEN BY OLE HOLM NIELSEN (NORDITA, COPENHAGEN) ON 22-MAR-1984
C
C     REFERENCES:
C          J.E.DENNIS AND J.J.MORE, SIAM REVIEW 19, 46 (1977)
C          C.G.BROYDEN, MATH. COMPUT. 19, 577 (1965)
C
C     UPDATE  (INTEGER) DETERMINES THE TYPE OF UPDATING
C     1 ..... USE PREVIOUS DIELECTRIC MATRIX (NO UPDATING)
C     2 ..... RECALCULATE DIELECTRIC MATRIX
C     3 ..... BROYDEN SCHEME, DENNIS AND MORE EQ.(4.8)
C     4 ..... SYMMETRIC SCHEME, DENNIS AND MORE EQ.(7.5)
C     5 ..... BFGS SCHEME, DENNIS AND MORE EQ.(7.25)
C
C     IDMAT . (INTEGER) CONTROLS WHERE DIELECTRIC MATRIX COMES FROM
C     1 ..... SIMPLE MIXING MATRIX
C     2 ..... FREE ELECTRON LINDHARD MATRIX
C     3 ..... FULL LINDHARD MATRIX
C     4 ..... USE PREVIOUS DIELECTRIC MATRIX
C
C     VCL ........ CLEAN-UP: V(G) < VCL ARE ZEROED
C     ALPHAS ..... MIXING COEFF. ALPHA FOR G**2 < GSQUBS
C     ALPHAB ..... MIXING COEFF. ALPHA FOR G**2 > GSQUBS
C     GSQUBS ..... G**2 SEPARATING THE MIXING COEFFICIENTS
C            ..... OR CUTOFF FOR PLANE WAVE BASIS IN LFERPA
C     VIN ........ INPUT AND NEW MIXED POTENTIAL
C     VOUT ....... OUTPUT POTENTIAL FROM SELFCONSISTENT LOOP
C     HK,XK,FK ... MIXING ARRAYS
C     IGLIST ..... LIST OF G-VECTORS
C           NGTOT,NG1,NG2,NG3: AUXILIARY TO IGLIST
C     NTAB ....... INDEX TABLE FOR THE G-VECTOR STARS
C                  N.B.:  M U S T  BE STARS, NOT SHELLS !!
C     G2TAB ...... G**2 OF G-VECTORS IN NTAB
C     A1,A2,A3 ... LATTICE TRANSLATION VECTORS
C     B1,B2,B3 ... BASIS VECTORS OF RECIPROCAL SPACE
C     VOLUM ...... UNIT CELL VOLUME
C     ULA ........ LATTICE CONSTANT
C     NEL ........ NUMBER OF ELECTRONS IN THE UNIT CELL
C     R .......... ROTATION MATRICES OF THE CRYSTAL
C     IHG,IHC,ISY,LI,NC,IB ... AUXILIARY GROUP SYMMETRY INFORMATION.
C     VTRANS ..... NON-SYMMORPHIC TRANSLATION VECTORS
C
      PARAMETER (EXFACT = 0.8)
      DOUBLE PRECISION EF
      LOGICAL LREAL,LSEMI
      COMPLEX VIN(NGTOT),VOUT(NGTOT),HKOLD(0:NDIM12,0:NDIM12)
      COMPLEX HK(0:NDIM12,0:NDIM12),XK(NDIM12),FK(NDIM12)
      INTEGER IGLIST(3,NGTOT),LIST(-NG1:NG1,-NG2:NG2,-NG3:NG3)
      COMPLEX CHI(0:NPWCHI,0:NPWCHI)
      COMPLEX INV(0:NPWCHI)
      INTEGER NTAB(NDIM2),IB(48),UPDATE,RTABLE(48,48)
      INTEGER IDIFF(3,NPWCHI),IPG(NPWCHI,NPWCHI),INDEX(NPWCHI)
      REAL A1(3),A2(3),A3(3),B1(3),B2(3),B3(3),R(49,3,3),RB(48,3,3)
      REAL EE1(NDIM13),EE2(NDIM8),G2TAB(NDIM2),VTRANS(3,48)
      COMPLEX CEV(NDIM5,NDIM8)
      INTEGER LISTAB(NDIM5),LWGHT(NDIM10),RLIST((NDIM12+1)*48)
      COMPLEX YK(NDIM12),SK(NDIM12),SKHK(NDIM12),SMHY(NDIM12)
C
C.....FILES
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
C.....PHYSICAL AND MATHEMATICAL CONSTANTS
      DOUBLE PRECISION ABOHR,RYEV,RYDERG,PI,SPI
      COMMON /CONST/   ABOHR,RYEV,RYDERG,PI,SPI
C
      SAVE IPAD,INIT,IJACOB,IHIGH
      DATA INIT,IJACOB /2*0/
C-----------------------------------------------------------------------
C
C     INITIALIZATION
C
      WRITE(IOUT,*) ' IDMAT = ',IDMAT,' UPDATE = ',UPDATE
      IF (INIT .EQ. 0) THEN
C
C       CHECK INPUT DATA
        IF (UPDATE .GT. 5 .OR. UPDATE .LE. 0) THEN
          WRITE(IOUT,*) 'UPDATE = ',UPDATE,' IS ILLEGAL'
          WRITE(IOUT,*) 'No mixing done'
          RETURN
          ENDIF
        IF (IDMAT .GT. 4 .OR. IDMAT .LE. 0) THEN
          WRITE(IOUT,*) 'IDMAT = ',IDMAT,' IS ILLEGAL'
          WRITE(IOUT,*) 'No mixing done'
          RETURN
          ENDIF
        IF (VCL .GT. 1.0E-4) THEN
          WRITE (IOUT,*) 'MIXVG *** WARNING *** VCLEAN TOO BIG ',VCL
          VCL = 0.0
          ENDIF
C
C       IF NHK IS GREATER THAN ZERO THEN HK XK FK WERE READ FROM FILE
C       ONLY IF THEY ARE NON ZERO CAN A JACOBIAN UPDATE METHOD BE USED
C       ON THE FIRST CYCLE. IJACOB FLAGS THIS.
        NHKOLD = NHK
        IF (NHKOLD .GT. 0) THEN
          HKSUM = 0.0
          XKSUM = 0.0
          FKSUM = 0.0
          DO 10 J = 1,NHKOLD
          DO 10 I = 1,NHKOLD
            HKSUM = HKSUM + ABS(HK(I,J))
10          CONTINUE
          DO 20 I = 1,NHKOLD
            XKSUM = XKSUM + ABS(XK(I))
            FKSUM = FKSUM + ABS(FK(I))
20          CONTINUE
          IF (HKSUM*XKSUM*FKSUM .GT. 1.0E-8) IJACOB = 1
          ENDIF
C
C       SET NHK AND CHECK THE NUMBER OF G-VECTORS IN EACH STAR.
        NHK = NDIM12
        DO 30 I = 1,NDIM12
          JMIN = NTAB(I+1)
          JMAX = NTAB(I+2) - 1
C         MAKE SURE THAT WE DO NOT EXCEED THE STARS CONTAINED IN NTAB
          IF (JMAX .LT. 0) THEN
            NHK = I - 1
            GOTO 40
            ENDIF
C         NO STAR CAN HAVE MORE THAN 48 G-VECTORS
          IF (JMAX - JMIN + 1 .GT. 48) THEN
            WRITE (IOUT,*) 'MIXVG *** STAR NO.', I+1,' HAS',
     +      JMAX - JMIN + 1,' > 48 VECTORS'
            WRITE (IOUT,*) 'FREE ELECTRON LINDHARD MIXING WILL BE USED'
            NHK = 0
            UPDATE = 1
            IDMAT = 2
            ENDIF
C
30        CONTINUE
C
40      CONTINUE
C
C       IF NHKOLD WAS LESS THAN NHK AND IDMAT = 4 THEN WE HAVE TO
C       PAD OUT HK TO SIZE NHK. THIS IS SIGNALLED WITH IPAD.
C       IF PADDING IS NECESSARY THEN JACOBIAN UPDATING CANNOT BE
C       DONE ON THE FIRST CYCLE.
        IPAD = 0
        IF (IDMAT .EQ. 4) IPAD = NHK - NHKOLD
        IF (IPAD .GT. 0) IJACOB = 0
C
C       SET SWITCH FOR MIXING STARS OUTSIDE OF HK
        IF (IDMAT .EQ. 1) THEN
          IHIGH = 1
        ELSE
          IHIGH = 2
          ENDIF
C
        WRITE (IUNIT8,*) 'MIXVG - JACOBIAN HAS DIMENSION = ',NHK
        INIT = 1
        ENDIF
C
C----------------------------------------------------------------------
C
C     OBTAIN NEW INPUT HK MATRIX
C
      IF (IDMAT .EQ. 1) THEN
C
C       SIMPLE MIXING FUNCTION
        WRITE(IOUT,*) ' SIMPLE MIXING FUNCTION  '
C
        DO 50 J = 1,NDIM12
        DO 50 I = 1,NDIM12
          HK(I,J) = (0.0,0.0)
50        CONTINUE
        DO 60 I = 1,NDIM12
          Q2 = G2TAB(I+1)
          IF (Q2 .LE. GSQUBS) HK(I,I) = - ALPHAS
          IF (Q2 .GT. GSQUBS) HK(I,I) = - ALPHAB
60        CONTINUE
C
        IJACOB = 0
C
      ELSE IF (IDMAT .EQ. 2) THEN
C
C       FREE ELECTRON LINDHARD FUNCTION
        WRITE(IOUT,*) ' FREE ELECTRON LINDHARD FUNCTION '
C
        DO 70 J = 1,NDIM12
        DO 70 I = 1,NDIM12
          HK(I,J) = (0.0,0.0)
70        CONTINUE
        DO 80 I = 1,NDIM12
          Q2 = G2TAB(I+1)
          CALL EPS1(NEL,VOLUM,ULA,Q2,EXFACT,FERMIK,EPSQ)
          HK(I,I) = - 1.0/EPSQ
80        CONTINUE
C
        IJACOB = 0
C
      ELSE IF (IDMAT .EQ. 3) THEN
C
C       FULL LINDHARD FUNCTION
        WRITE(IOUT,*) ' FULL LINDHARD FUNCTION '
C
        CALL USAGE(N5,T5)
        CALL LFERPA(HK,CHI,INV,NDIM12,NHK,IGLIST,NGTOT,NG1,NG2,NG3,
     +   A1,A2,A3,B1,B2,B3,R,RB,RTABLE,VTRANS,IB,NC,RLIST,LIST,LISTAB,
     +   EE1,EE2,NDIM13,IEVFIL,NUMKPT,LREAL,NDIM5,NDIM8,CEV,LWGHT,
     +   NDIM10,MINBDS,EF,LSEMI,G2TAB,NTAB,NDIM2,ULA,VOLUM,ISPIN,NSPIN,
     +   IDIFF,IPG,INDEX,NPWCHI,IFILEH)
        CALL USAGE(N6,T6)
        WRITE(IOUT,*) ' TIME SPENT IN LFERPA = ', T6 - T5
C
        IJACOB = 0
C
      ELSE IF (IDMAT .EQ. 4) THEN
C
C       USE PREVIOUS DIELECTRIC FUNCTION
        WRITE(IOUT,*) ' USE PREVIOUS DIELECTRIC FUNCTION '
C
        IF (NHK .LE. 0 .OR. NHK .GT. NDIM12) THEN
          WRITE(IOUT,*) 'NHK = ',NHK,' AND NDIM12 = ',NDIM12
          WRITE(IOUT,*) 'FREE ELECTRON LINDHARD MIXING WILL BE USED'
          NHK = 0
          UPDATE = 1
          IJACOB = 0
          ENDIF
C
C       PAD HK TO LENGTH NHK IF REQUIRED
        IF (IPAD .GT. 0) THEN
          WRITE(IOUT,*) ' IPAD = ',IPAD
          DO 105 J = NHKOLD+1,NHK
          DO 105 I = NHKOLD+1,NHK
            HK(I,J) = (0.0,0.0)
            IF (I .EQ. J) THEN
              Q2 = G2TAB(I+1)
              CALL EPS1(NEL,VOLUM,ULA,Q2,EXFACT,FERMIK,EPSQ)
              HK(I,I) = - 1.0/EPSQ
              ENDIF
105         CONTINUE
          IPAD = 0
          ENDIF
C
        ENDIF
C
C     UNLESS RECALCULATE HK ON NEXT CYCLE SET IDMAT = 4
      IF (UPDATE .NE. 2) IDMAT = 4
C
C----------------------------------------------------------------------
C
C     UPDATE JACOBIAN IF REQUIRED
C
      IF (UPDATE .EQ. 3 .OR. UPDATE .EQ. 4 .OR. UPDATE .EQ. 5) THEN
        IF (IJACOB .EQ. 1) THEN
          CALL JACOB (UPDATE,VCL,VIN,VOUT,HK,XK,FK,NDIM2,NDIM12,
     +    NGTOT,NTAB,YK,SK,SKHK,SMHY,HKOLD,NHK)
        ELSE IF (IJACOB .EQ. 0) THEN
          IJACOB = 1
          ENDIF
        ENDIF
C
C     UPDATE THE ARRAYS XK AND FK. IF JACOB WAS CALLED THEN THIS
C     WAS ALREADY DONE, BUT FOR SIMPLICITY WE DO IT AGAIN ANYWAY.
C
      DO 110 I = 1,NHK
        JMIN = NTAB(I+1)
        CALL VCLEAN (VIN(JMIN),VCL)
        CALL VCLEAN (VOUT(JMIN),VCL)
        XK(I) = VIN(JMIN)
        FK(I) = VOUT(JMIN) - VIN(JMIN)
110     CONTINUE
C
C----------------------------------------------------------------------
C
C     UPDATE THE POTENTIALS USING HK
C
      CALL UPDATV (VCL,VIN,VOUT,HK,FK,NDIM2,NDIM12,
     +IGLIST,LIST,NGTOT,NTAB,NG1,NG2,NG3,A1,A2,A3,B1,B2,B3,
     +LI,NC,IB,R,VTRANS,NHK)
C
C----------------------------------------------------------------------
C
C     MIX THE REMAINING STARS BY THE DIAGONAL DIELECTRIC METHOD
C
      DO 140 I = NHK+2,NDIM2
        Q2 = G2TAB(I)
        IF (Q2 .LE. 0.0) GOTO 150
        JMIN = NTAB(I)
        JMAX = NTAB(I+1) - 1
        IF (IHIGH .EQ. 1) THEN
          ALPHA = ALPHAB
          IF (Q2 .LE. GSQUBS) ALPHA = ALPHAS
        ELSE IF (IHIGH .EQ. 2) THEN
          CALL EPS1 (NEL,VOLUM,ULA,Q2,EXFACT,FERMIK,EPSQ)
          ALPHA = 1.0/EPSQ
          ENDIF
C
        DO 130 J = JMIN,JMAX
          VIN(J) =  ALPHA * VOUT(J) + (1.0 - ALPHA) * VIN(J)
          CALL VCLEAN (VIN(J),VCL)
130       CONTINUE
140     CONTINUE
C
150   RETURN
      END
      SUBROUTINE VCLEAN (VIN,VCL)
C
C     CLEAN UP VARIABLE VIN
C
      COMPLEX VIN
C
      IF (VCL .LE. 0.0 .OR. VCL .GT. 1.0E-4) RETURN
      VR = REAL(VIN)
      VI = AIMAG(VIN)
      IF (ABS(VR) .LT. VCL) VR = 0.0
      IF (ABS(VI) .LT. VCL) VI = 0.0
      VIN = CMPLX(VR,VI)
      RETURN
      END
      SUBROUTINE JACOB (UPDATE,VCL,VIN,VOUT,HK,XK,FK,NDIM2,NDIM12,
     +NGTOT,NTAB,YK,SK,SKHK,SMHY,HKOLD,NHK)
C
C     UPDATING THE JACOBIAN FOR POTENTIAL MIXING
C     SEE SUBROUTINE MIXVG FOR DETAILS
C
      COMPLEX VIN(NGTOT),VOUT(NGTOT)
      COMPLEX HK(0:NDIM12,0:NDIM12),XK(NDIM12),FK(NDIM12)
      COMPLEX YK(NDIM12),SK(NDIM12),SKHK(NDIM12),SMHY(NDIM12)
      COMPLEX HKOLD(0:NDIM12,0:NDIM12),SHY,SMHYY,SY
      INTEGER UPDATE,NTAB(NDIM2)
C
C.....FILES
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
C-----------------------------------------------------------------------
C
C     WRITE (IUNIT8,*) 'MIXVG - XK, FK, YK, SK (REAL,IMAG) ='
      SY = (0.0,0.0)
      DO 105 I = 1,NHK
        J = NTAB(I+1)
        CALL VCLEAN (VIN(J),VCL)
        CALL VCLEAN (VOUT(J),VCL)
        YK(I) = (VOUT(J) - VIN(J)) - FK(I)
        SK(I) = VIN(J) - XK(I)
        SY = SY + CONJG(SK(I)) * YK(I)
C       UPDATE XK AND FK
        XK(I) = VIN(J)
        FK(I) = VOUT(J) - VIN(J)
C       WRITE (IUNIT8,90) XK(I),FK(I),YK(I),SK(I)
C90     FORMAT(1X,4(2F12.5,5X))
C       FOR THE BFGS SCHEME, WE NEED TO STORE HK-OLD
        DO 100 J = 1,NHK
          HKOLD(I,J) = HK(I,J)
100       CONTINUE
105     CONTINUE
C
C     CONSTRUCT THE NEW INVERSE JACOBIAN
C
      SHY = (0.0,0.0)
      SMHYY = (0.0,0.0)
      DO 120 I = 1,NHK
        SKHK(I) = (0.0,0.0)
        SMHY(I) = SK(I)
        DO 110 J = 1,NHK
          SKHK(I) = SKHK(I) + CONJG(SK(J)) * HK(J,I)
          SMHY(I) = SMHY(I) - HK(I,J) * YK(J)
110       CONTINUE
        SHY   = SHY   +       SKHK(I)  * YK(I)
        SMHYY = SMHYY + CONJG(SMHY(I)) * YK(I)
120     CONTINUE
C     WRITE (IUNIT8,*) 'MIXVG: SHY SMHYY = ',SHY,SMHYY
C
      DO 140 I = 1,NHK
      DO 140 J = 1,NHK
C
        IF (UPDATE .EQ. 3) THEN
C
C         BROYDEN UPDATING: DENNIS AND MORE, EQ.(4.8)
          HK(I,J) = HK(I,J) + SMHY(I) * SKHK(J) / SHY
C
        ELSE IF (UPDATE .EQ. 4) THEN
C
C         SYMMETRIC UPDATING: DENNIS AND MORE, EQ.(7.5)
          HK(I,J) = HK(I,J) + SMHY(I) * CONJG(SMHY(J)) / SMHYY
C
        ELSE IF (UPDATE .EQ. 5) THEN
C
C         THE BFGS SCHEME, DENNIS AND MORE EQ.(7.25)
          HK(I,J) = SK(I) * CONJG(SK(J)) / SY
          DO 130 K = 1,NHK
          DO 130 L = 1,NHK
            DIK = 0.0
            IF (I .EQ. K) DIK = 1.0
            DLJ = 0.0
            IF (L .EQ. J) DLJ = 1.0
            HK(I,J) = HK(I,J) + (DIK - SK(I)*CONJG(YK(K))/SY) *
     +                HKOLD(K,L) * (DLJ - YK(L)*CONJG(SK(J))/SY)
130         CONTINUE
C
          ENDIF
C
140     CONTINUE
C     WRITE (IUNIT8,*) 'THE UPDATED INVERSE-JACOBIAN (10 BY 10) IS'
C     WRITE (IUNIT8,170) ((HK(I,J),J=1,10),I=1,10)
C170  FORMAT (5(1X,2F11.5,2X))
C
      RETURN
      END
      SUBROUTINE UPDATV (VCL,VIN,VOUT,HK,FK,NDIM2,NDIM12,
     +IGLIST,LIST,NGTOT,NTAB,NG1,NG2,NG3,A1,A2,A3,B1,B2,B3,
     +LI,NC,IB,R,VTRANS,NHK)
C
C     UPDATE THE POTENTIAL USING HK
C     SEE SUBROUTINE MIXVG FOR DETAILS
C
      COMPLEX VIN(NGTOT),VOUT(NGTOT),VNEW,VOLD,VCORR
      COMPLEX HK(0:NDIM12,0:NDIM12),FK(NDIM12)
      INTEGER IGLIST(3,NGTOT),LIST(-NG1:NG1,-NG2:NG2,-NG3:NG3)
      INTEGER NTAB(NDIM2),IB(48),IDONE(48),IGL(3)
      REAL A1(3),A2(3),A3(3),B1(3),B2(3),B3(3),R(49,3,3)
      REAL VTRANS(3,48),GVEC(3),GROT(3)
C
C.....FILES
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
C.....PHYSICAL AND MATHEMATICAL CONSTANTS
      DOUBLE PRECISION ABOHR,RYEV,RYDERG,PI,SPI
      COMMON /CONST/   ABOHR,RYEV,RYDERG,PI,SPI
C-----------------------------------------------------------------------
C
      NMIN = MIN(10,NHK)
      DO 20 I = 1,NMIN
        WRITE(IUNIT8,10) (HK(I,J),J=1,NMIN)
 10     FORMAT(1X,10(2F6.3,1X))
 20   CONTINUE
      WRITE (IUNIT8,*) 'THE UPDATED VNEW AND VCORR'
C
      DO 280 I = 1,NHK
C
        JMIN = NTAB(I+1)
        JMAX = NTAB(I+2) - 1
        VOLD = VIN(JMIN)
        CALL VCLEAN (VOLD,VCL)
        VNEW = VOLD
        DO 180 J = 1,NHK
          VNEW = VNEW - HK(I,J) * FK(J)
180       CONTINUE
        CALL VCLEAN (VNEW,VCL)
C       FIND THE CORRECTION FACTOR: VNEW/VOLD
        IF (ABS(VOLD) .LT. 1.0E-14) THEN
          WRITE(IUNIT8,185) JMIN
185       FORMAT(1X,I6,' TOO SMALL FOR TEST')
        ELSE
          VCORR = VNEW/VOLD
          CALL VCLEAN (VCORR,VCL)
          WRITE (IUNIT8,190) JMIN,VNEW,VCORR
190       FORMAT (1X,I4,2X,2G15.6,5X,2G15.6)
          IF (ABS(VCORR) .GT. 10.0)
     +      WRITE (IOUT,200) VCORR,I+1,VOLD,VNEW
200         FORMAT (' MIXVG - VCORR UNREALISTIC = ',2G12.5,/,
     +      ' STAR NO.',I5,' VOLD,VNEW =',2G12.5,2X,2G12.5)
          ENDIF
C
C       GIVEN THE UPDATED V-IN, WE HAVE TO APPLY THE UPDATE TO THE
C       REST OF THE STAR OF G-VECTORS
C       HERE WE RELY ON THE FACT THAT THE V-IN FOR A GIVEN G-STAR
C       ARE RELATED TO EACH OTHER SIMPLY BY A PHASE FACTOR.
C       FOR LATTICES WITHOUT INVERSION SYMMETRY, WE HAVE TO TREAT
C       (-G) SPECIALLY.
C       NOTE THAT THIS REQUIRES THAT NTAB CONTAINS G-STARS RATHER
C       THAN SHELLS OF EQUAL LENGTH (TAKEN CARE OF BY SBRT. GSHELL)
C
        DO 220 JJ = 1,48
          IDONE(JJ) = 0
220       CONTINUE
C
        INVADD = 0
        IF (LI .EQ. 0) INVADD = 1
        DO 250 INVERS = 0,INVADD
C
        IF (INVERS .EQ. 0) THEN
          IPLACE = JMIN
        ELSE
C         USE V(-G) = V(G)* FOR G=IGLIST(JMIN), AND APPLY ROTATIONS
          IGL(1) = -IGLIST(1,JMIN)
          IGL(2) = -IGLIST(2,JMIN)
          IGL(3) = -IGLIST(3,JMIN)
          CALL LOOKUP (LIST,NGTOT,IGL,IPLACE,NG1,NG2,NG3,1,1)
C         IF (IPLACE .LT. JMIN .OR. IPLACE .GT. JMAX) THEN TROUBLE
          VNEW = CONJG(VNEW)
          ENDIF
C
        CALL GCODE(IGLIST(1,IPLACE),NG1,NG2,NG3,B1,B2,B3,I1,I2,I3,GVEC)
C
C       LOOP OVER ROTATIONS
C
        DO 250 N = 1,NC
          IC = IB(N)
C         ROTATE G-VECTOR
          DO 230 II = 1,3
            GROT(II) = 0.0
            DO 230 JJ = 1,3
              GROT(II) = GROT(II) + R(IC,II,JJ) * GVEC(JJ)
230           CONTINUE
C         EXPRESS GROT IN BASIS B1,B2,B3
          IGL(1) = NINT( GROT(1)*A1(1) + GROT(2)*A1(2) + GROT(3)*A1(3) )
          IGL(2) = NINT( GROT(1)*A2(1) + GROT(2)*A2(2) + GROT(3)*A2(3) )
          IGL(3) = NINT( GROT(1)*A3(1) + GROT(2)*A3(2) + GROT(3)*A3(3) )
          CALL LOOKUP (LIST,NGTOT,IGL,IPLACE,NG1,NG2,NG3,1,1)
          IF (IPLACE .LT. JMIN .OR. IPLACE .GT. JMAX) THEN
            WRITE (IOUT,240) GROT,GVEC,IC
240         FORMAT('MIXVG - THE ROTATED VECTOR ',3F12.4/
     +        ' GENERATED FROM ',3F12.4,' BY ROT. NO. ',I5/
     +        ' IS OUTSIDE STAR')
            CALL EXIT
            ENDIF
C
          JJ = IPLACE - JMIN + 1
          IF (IDONE(JJ) .NE. 0) GOTO 250
          IDONE(JJ) = 1
C
C         PHASE FACTOR: EXP(-I*R*G*VTRANS(R))
          FAC = - 2.0*PI*( IGL(1)*VTRANS(1,N) + IGL(2)*VTRANS(2,N)
     +                   + IGL(3)*VTRANS(3,N) )
          VIN(IPLACE) = VNEW * CMPLX(COS(FAC),SIN(FAC))
          CALL VCLEAN (VIN(IPLACE),VCL)
C
250       CONTINUE
C
C       WRITE (IOUT,*) 'MIXVG - STAR NO.',I+1,' - V-MIX AND V-OLD'
        DO 270 J = JMIN,JMAX
          JJ = J - JMIN + 1
          IF (IDONE(JJ) .NE. 1) WRITE (IOUT,*) 'G=',J,' NOT UPDATED'
C         WRITE (IOUT,260) J,VIN(J),VOUT(J)
C260      FORMAT(1X,I5,2F11.5,5X,2F11.5)
270       CONTINUE
C
280   CONTINUE
C
      RETURN
      END
      SUBROUTINE LFERPA(HK,CHI,INV,NDIM12,NHK,IGLIST,NGTOT,NG1,NG2,NG3,
     +A1,A2,A3,B1,B2,B3,R,RB,RTABLE,VT,IB,NC,RLIST,LIST,LISTAB,
     +EE1,EE2,NDIM13,IEVFIL,NUMKPT,LREAL,NDIM5,NDIM8,CEV,LWGHT,
     +NDIM10,MINBDS,EF,LSEMI,G2TAB,NTAB,NDIM2,ULA,VOLUM,ISPIN,NSPIN,
     +IDIFF,IPG,INDEX,NPWCHI,IFILEH)
C
C  Calculates the dielectric response matrix in the Random Phase
C  Approximation, including Local Field Effects (off-diagonal elements).
C
      COMPLEX HK(0:NDIM12,0:NDIM12)
      COMPLEX CHI(0:NPWCHI,0:NPWCHI)
      COMPLEX INV(0:NPWCHI)
      INTEGER IGLIST(3,NGTOT),LIST(-NG1:NG1,-NG2:NG2,-NG3:NG3)
      INTEGER IDIFF(3,NPWCHI),IPG(NPWCHI,NPWCHI)
      INTEGER RTABLE(48,48),RLIST((NDIM12+1)*48),INVER(48)
      INTEGER NTAB(NDIM2)
      INTEGER LISTAB(NDIM5),INDEX(NPWCHI),LWGHT(NDIM10)
      INTEGER IB(48)
      REAL A1(3),A2(3),A3(3),RB(48,3,3),R(49,3,3)
      REAL B1(3),B2(3),B3(3),SUM(3,3)
      REAL EE1(NDIM13),EE2(NDIM8)
      REAL G2TAB(NDIM2),VT(3,48)
      DOUBLE PRECISION EF,FRED
      LOGICAL LREAL,LSEMI,LSTAT
      COMPLEX CEV(NDIM5,NDIM8)
C
C  For a calculation of the static dielectric constant, set LSTAT to
C  TRUE.
C
      PARAMETER (LSTAT=.FALSE.)
C
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
      DOUBLE PRECISION ABOHR,RYEV,RYDERG,PI,SPI
      COMMON /CONST/   ABOHR,RYEV,RYDERG,PI,SPI
C
C-----------------------------------------------------------------------
C
C  Express the rotation matrices in the basis B1,B2,B3, and reorder in
C  the order of IB.
C
      DO 40 N = 1,NC
        IC = IB(N)
        DO 10 I = 1,3
        DO 10 J = 1,3
          SUM(I,J) = 0.0
10        CONTINUE
        DO 20 I = 1,3
        DO 20 J = 1,3
          SUM(1,1) = SUM(1,1) + A1(I)*R(IC,I,J)*B1(J)
          SUM(2,1) = SUM(2,1) + A2(I)*R(IC,I,J)*B1(J)
          SUM(3,1) = SUM(3,1) + A3(I)*R(IC,I,J)*B1(J)
          SUM(1,2) = SUM(1,2) + A1(I)*R(IC,I,J)*B2(J)
          SUM(2,2) = SUM(2,2) + A2(I)*R(IC,I,J)*B2(J)
          SUM(3,2) = SUM(3,2) + A3(I)*R(IC,I,J)*B2(J)
          SUM(1,3) = SUM(1,3) + A1(I)*R(IC,I,J)*B3(J)
          SUM(2,3) = SUM(2,3) + A2(I)*R(IC,I,J)*B3(J)
          SUM(3,3) = SUM(3,3) + A3(I)*R(IC,I,J)*B3(J)
20        CONTINUE
        DO 30 I = 1,3
        DO 30 J = 1,3
          RB(N,I,J) = SUM(I,J)
30        CONTINUE
40      CONTINUE
C
      NVECHK = NTAB(NHK + 2) - 2
C  Since the actual size of CHI is not determined at the start of RUN207,
C  we use an estimate NDIM12*NGROUP which we now check is sufficient.
C
      IF (NVECHK .GT. NPWCHI) THEN
        WRITE(IOUT,*) ' PROBLEMS! SIZE OF CHI MATRIX IS ',NVECHK
        WRITE(IOUT,*) ' BUT NPWCHI IS ONLY ',NPWCHI
        CALL EXIT
        ENDIF
C
C  Initialize arrays.
C
      DO 50 J = 0,NVECHK
      DO 50 I = 0,NVECHK
50      CHI(I,J) = (0.0,0.0)
C
      DO 60 J = 0,NHK
      DO 60 I = 0,NHK
        HK(I,J) = 0.0
60      CONTINUE
C
C Set up the multiplication table for the full rotation group of the crystal.
C
      CALL MLTPLY(RTABLE,RB,INVER,NC)
C
C Compile a list such that any G-vector can be constructed by taking the
C rotations listed for it in RLIST and operating with any of them on the
C first G-vector of the star listed in IGLIST.
C
      CALL ROTATE(RLIST,NDIM12,(NHK+1),NC,IGLIST,NGTOT,NTAB,NDIM2,RB)
C
      FRED = 4*ABOHR*RYEV/(VOLUM*PI*ULA)
C
C  The top NVECHK x NVECHK left hand corner of CHI is calculated using
C the Adler-Wiser formula
C
      CALL CALCHI(CHI,NDIM12,NHK,NVECHK,IGLIST,LIST,NGTOT,CEV,
     +  ISPIN, NSPIN,
     +  LISTAB,INDEX,NDIM5,NDIM8,NTAB,NDIM2,RB,VT,IB,INVER,
     +  NPWCHI,MINBDS,EF,LSEMI,NUMKPT,LREAL,EE1,EE2,NDIM13,
     +  IDIFF,IPG,RTABLE,RLIST,NC,LWGHT,NDIM10,
     +  NG1,NG2,NG3,A1,A2,A3,IEVFIL,ULA,LSTAT)
C
C     WRITE(IOUT,*) ' CHI : '
C     NMIN = MIN(15,NVECHK)
C     DO 70 I = 0,NMIN
C       WRITE(IOUT,65) (CHI(I,II),II=0,3)
C       WRITE(IOUT,65) (CHI(I,II),II=4,7)
C       WRITE(IOUT,65) (CHI(I,II),II=8,11)
C       WRITE(IOUT,65) (CHI(I,II),II=12,NMIN)
C     WRITE(IOUT,*)
C65     FORMAT(1X,4(2F9.6,2X))
C70     CONTINUE
C
C  We use a symmetrized form of Eps to calculate the static dielectric
C  constant, defined as
C       E(G,Gi) = d(G-Gi) + X(G,Gi)*(scaling factor)/|G||Gi|.
C
      DO 90 I  = 0,NHK
      DO 90 II = 0,NHK
C
        JMIN = NTAB(II+1)
        JMAX = NTAB(II+2) - 1
C
        DO 90 J = JMIN,JMAX
          DENOM = SQRT(G2TAB(I+1))*SQRT(G2TAB(II+1))
          IF (I  .EQ. 0) DENOM = SQRT(G2TAB(II+1))
          IF (II .EQ. 0) DENOM = SQRT(G2TAB(I+1))
          IF ((I .EQ. 0) .AND. (II .EQ. 0)) DENOM = 1.0
          CHI((NTAB(I+1)-1),(J-1)) = CHI((NTAB(I+1)-1),(J-1))*
     +      SNGL(FRED)/DENOM
 90       CONTINUE
C
C  Whole CHI is constructed from certain elements using symmetry.
C
       CALL EXPAND(CHI,NPWCHI,NDIM12,NHK,NVECHK,RLIST,NC,RTABLE,
     + NTAB,NDIM2,VT,IB,IGLIST,NGTOT)
C
       DO 110 I = 0,NVECHK
         CHI(I,I) = CHI(I,I) + 1
 110     CONTINUE
C
C      WRITE(IOUT,*) ' DIELECTRIC RESPONSE MATRIX :'
C      NMIN = MIN(15,NVECHK)
C      DO 100 I = 0,NMIN
C        WRITE(IOUT,105) (CHI(I,II),II=0,3)
C        WRITE(IOUT,105) (CHI(I,II),II=4,7)
C        WRITE(IOUT,105) (CHI(I,II),II=8,11)
C        WRITE(IOUT,105) (CHI(I,II),II=12,NMIN)
C        WRITE(IOUT,*)
C105     FORMAT(1X,4(2F9.6,2X))
C100     CONTINUE
C
       WRITE(IOUT,*) ' CHI(0,0) = ',CHI(0,0)
C
       ISTAT = 1
C
       IF (LSTAT) THEN
C
         REWIND(UNIT=IFILEH)
         WRITE(IFILEH) CHI
C
         ISTAT = 0
         CALL INVERT(CHI(0,0),NPWCHI,NDIM12,NVECHK,INV(0),IDIFF(1,1),
     +     ISTAT,HK(0,0),IGLIST,NGTOT,NTAB,NDIM2,NHK,VT,IB,RLIST,NC)
         ISTAT = 1
C
         REWIND(UNIT=IFILEH)
         READ(IFILEH) CHI
C
        ENDIF
C
C  Converting to the non-symmetrized form :
C
      DO 140 I = 1,NVECHK
        CHI(0,I) = (0.0,0.0)
        CHI(I,0) = (0.0,0.0)
 140  CONTINUE
C
      DO 150 I  = 1,NHK
      DO 150 II = 1,NHK
C
        J1MIN = NTAB(I+1)
        J1MAX = NTAB(I+2) - 1
C
        J2MIN = NTAB(II+1)
        J2MAX = NTAB(II+2) - 1
C
        DO 150 J  = J1MIN,J1MAX
        DO 150 JJ = J2MIN,J2MAX
C
        CHI((J-1),(JJ-1)) = CHI((J-1),(JJ-1))*
     +    SQRT(G2TAB(II+1))/SQRT(G2TAB(I+1))
C
 150    CONTINUE
C
      WRITE(IOUT,*) ' CHI : '
      DO 1 I = 0,4
       WRITE(IOUT,2) (CHI(I,J),J=0,4)
 2     FORMAT(1X,5(2F7.4,2X))
 1    CONTINUE
      CALL INVERT(CHI(0,0),NPWCHI,NDIM12,NVECHK,INV(0),IDIFF(1,1),
     +    ISTAT,HK(0,0),IGLIST,NGTOT,NTAB,NDIM2,NHK,VT,IB,RLIST,NC)
C
      WRITE(IOUT,*) ' HK : '
      NMIN = MIN(8,NHK)
      DO 180 I = 0,NMIN
        WRITE (IOUT,170) (HK(I,J),J=0,NMIN)
 170    FORMAT (1X,8(2F7.4,2X))
 180    CONTINUE
C
      RETURN
      END
      SUBROUTINE MLTPLY (RTABLE,RB,INVER,NC)
C
C Sets up the multiplication table of a set of NC matrices.
C
      INTEGER RTABLE(48,48),INVER(48)
      REAL RB(48,3,3),MULT(3,3)
C
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
C
C
      DO 100 I = 1,NC
      DO 100 J = 1,NC
C
      DO 10 K = 1,3
      DO 10 L = 1,3
        MULT(K,L) = RB(I,K,1)*RB(J,1,L) + RB(I,K,2)*RB(J,2,L)
     +            + RB(I,K,3)*RB(J,3,L)
 10     CONTINUE
C
C
      DO 30 K = 1,NC
        A = ABS(MULT(1,1) - RB(K,1,1)) + ABS(MULT(1,2) - RB(K,1,2))
     +    + ABS(MULT(1,3) - RB(K,1,3)) + ABS(MULT(2,1) - RB(K,2,1))
     +    + ABS(MULT(2,2) - RB(K,2,2)) + ABS(MULT(2,3) - RB(K,2,3))
     +    + ABS(MULT(3,1) - RB(K,3,1)) + ABS(MULT(3,2) - RB(K,3,2))
     +    + ABS(MULT(3,3) - RB(K,3,3))
        IF (A .LE. 1.0E-6) GOTO 40
 30    CONTINUE
C
      WRITE(IOUT,*) 'HELP!  I = ',I,' J = ',J
      CALL EXIT
C
 40   RTABLE(I,J)  = K
      IF (K .EQ. 1) INVER(I) = J
C
 100  CONTINUE
C
      RETURN
      END
      SUBROUTINE ROTATE(RLIST,NDIM12,NHK,NC,IGLIST,NGTOT,NTAB,NDIM2,RB)
C
C  For a fixed group of NC symmetry operations (RB), and a set of
C vectors arranged into stars, determines which operations connect the
C first listed vector in any star to all others in it.
C
      INTEGER NTAB(NDIM2),IGLIST(3,NGTOT),IRG(3)
      REAL RB(48,3,3)
      INTEGER RLIST((NDIM12+1)*48)
C
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
C
C
      DO 5 I = 1,((NDIM12+1)*48)
5       RLIST(I) = 0
C
      DO 100 I = 1,NHK
        JMIN = NTAB(I)
        JMAX = NTAB(I+1) - 1
C
        LSTAR = NTAB(I+1) - NTAB(I)
        LEN = NC/LSTAR
C
        DO 70 NROT=1,NC
C
        DO 20 K = 1,3
          IRG(K) = 0.0
          DO 20 L = 1,3
 20         IRG(K) = IRG(K) + RB(NROT,K,L)*IGLIST(L,JMIN)
C
        DO 40 J = JMIN,JMAX
          A = ABS(IRG(1) - IGLIST(1,J)) + ABS(IRG(2) - IGLIST(2,J))
     +      + ABS(IRG(3) - IGLIST(3,J))
          IF (A .LT. 1.0E-6) GOTO 50
 40       CONTINUE
C
        WRITE(IOUT,*) ' HELP!  NROT = ',NROT,' SHELL ',I
        CALL EXIT
 50     KMIN = NC*(I-1) + LEN*(J - JMIN) + 1
        KMAX = KMIN + LEN - 1
C
        DO 60 K = KMIN,KMAX
         IF (RLIST(K) .EQ. 0) THEN
           RLIST(K) = NROT
           GOTO 70
           ENDIF
 60       CONTINUE
C
 70     CONTINUE
C
C       DO 100 J = JMIN,JMAX
C         WRITE(IOUT,*) ' FOR SHELL ',I,(J-JMIN + 1),' TH VECTOR '
C         JSTART=NC*(I-1)+LEN*(J-JMIN)
C         WRITE(IOUT,*) (RLIST(JSTART+K),K=1,LEN)
 100      CONTINUE
C
      RETURN
      END
      SUBROUTINE CALCHI (CHI,NDIM12,NHK,NVECHK,IGLIST,LIST,NGTOT,CEV,
     + ISPIN,NSPIN,
     + LISTAB,INDEX,NDIM5,NDIM8,NTAB,NDIM2,RB,VT,IB,INVER,
     + NPWCHI,MINBDS,EF,LSEMI,NUMKPT,LREAL,EE1,EE2,NDIM13,
     + IDIFF,IPG,RTABLE,RLIST,NC,LWGHT,NDIM10,
     + NG1,NG2,NG3,A1,A2,A3,IEVFIL,ULA,LSTAT)
C
C  Calculates certain elements of the polarizability matrix, X(G,G'), in
C  the Random Phase Approximation. (Calculates at one G-vector in each
C  star, for all G' vectors; all other elements can be related by
C  symmetry.)
C
      INTEGER IGLIST(3,NGTOT),LIST(-NG1:NG1,-NG2:NG2,-NG3:NG3)
      INTEGER NTAB(NDIM2),LWGHT(NDIM10),LROT(48),INVER(48)
      INTEGER RTABLE(48,48),LISTAB(NDIM5),INDEX(NPWCHI)
      INTEGER RLIST((NDIM12+1)*48),IB(48)
      INTEGER IDIFF(3,NPWCHI),IPG(NPWCHI,NPWCHI)
      REAL A1(3),A2(3),A3(3),RB(48,3,3)
      REAL VT(3,NC),WVK(3),EE1(NDIM13),EE2(NDIM8)
      COMPLEX CHI(0:NPWCHI,0:NPWCHI),PHASE
      DOUBLE PRECISION EF
C
      LOGICAL LREAL,LSEMI,LSTAT
C
      DOUBLE PRECISION ABOHR,RYEV,RYDERG,PI,SPI
      COMMON /CONST/   ABOHR,RYEV,RYDERG,PI,SPI
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
      COMMON /MACH1/ MACTYP
C
C---------------------------------------------------------------------
C For complex Hamiltonians
C     COMPLEX CEV(NDIM5,NDIM8),C1,CENTRE,Z,CONJUG,ZERO
C     CONJUG(Z) = CONJG(Z)
C     ZERO = (0.0,0.0)
C     IF (LREAL) WRITE(IOUT,*) '****WARNING**CALCHI**LREAL.TRUE****'
C---------------------------------------------------------------------
C For real Hamiltonians
      REAL CEV(NDIM5,NDIM8),C1,CENTRE
      CONJUG(X) = X
      ZERO = 0.0
      IF (.NOT. LREAL) WRITE(IOUT,*) '**WARNING**CALCHI**LREAL.FALSE*'
C---------------------------------------------------------------------
C
      IWGHT = 0
      DO 5 I = 1,NUMKPT
  5    IWGHT = IWGHT + LWGHT(I)
C
      REWIND (UNIT=IKPTS)
C
      DO 7 I = 1,NPWCHI
      DO 7 J = 1,NPWCHI
 7     IPG(I,J) = 0
C
C  Loop over special points.
C
      DO 1000 ISPPT = 1,NUMKPT
C
C     In reading we are only concerned with the case = ISPIN
C     (Omit reading GK)
      IF (NSPIN .EQ. 1) THEN
          READ (IKPTS) ICHECK,ISPCHK,WVK,NBDS,NA,NB,LISTAB
      ELSE
C       Spin-polarized case: neglect half the file.
        IF (ISPIN .EQ. 1) THEN
          READ (IKPTS) ICHECK,ISPCHK,WVK,NBDS,NA,NB,LISTAB
          READ (IKPTS)
        ELSE
          READ (IKPTS)
          READ (IKPTS) ICHECK,ISPCHK,WVK,NBDS,NA,NB,LISTAB
          ENDIF
        ENDIF
      IF (ICHECK .NE. ISPPT .OR. ISPCHK .NE. ISPIN) THEN
        WRITE (IOUT,*) 'CALCHI *** ERROR IN K-POINT FILE',
     +    ISPPT,ICHECK, '    ', ISPIN, ISPCHK
        CALL EXIT
        ENDIF
C     Read eigenvectors and eigenvalues from file
      IREC = ((ISPPT-1)*NSPIN + (ISPIN - 1)) * NDIM8
      CALL RWEV (IEVFIL,IREC,-1,EE1,EE2,NBDS,NA+NB,
     +      NDIM5,NDIM8,CEV,LREAL)
C
C  Determine limits of sums over valence and conduction bands.
C
      IF (LSEMI) THEN
            NV = MINBDS
      ELSE
            DO 8 I = 1,NBDS
              IF (EE2(I) .GT. EF) GOTO 9
 8          CONTINUE
            I = NBDS + 1
 9          NV = I - 1
      ENDIF
C
      NCOND = NBDS - NV
      IF (NCOND .LE. 0) WRITE(IOUT,*) '**WARNING**NCOND = ',NCOND
C
      IF (ISPPT .EQ. 1) WRITE(IOUT,*) ' CALCULATING WITH ',NV,
     +   ' VALENCE, AND ',NCOND,' CONDUCTION BANDS, AT ',NUMKPT,
     +   ' SPECIAL K-POINTS.  PLANE-WAVE BASIS SET RESTRICTED TO ',
     +   NHK+1, ' STARS. '
C
C
C Read unfolding rotations for ith k-point.
C
      READ(IN290,*) (LROT(J),J=1,48)
C
C Choose subset of vectors from which to calculate approximate matrix
C elements.  NO is the total number of plane waves included in the set.
C
      NO = MIN((NA+NB),NPWCHI)
C
      CALL CUTOFF(NO,LISTAB,NDIM5,NTAB((NHK+2)),INDEX,NPWCHI)
C
      IF (ISPPT .EQ. 1) WRITE(IOUT,*) NO,' PLANE WAVES USED.'
C
      IF (NO .EQ. NPWCHI) WRITE(IOUT,*) ' WORKING ARRAY SIZE ',
     +    'LIMITED BY NPWCHI ( = ',NPWCHI,' )'
C
C Calculate all possible differences between vectors of subset.
C
      IF (ISPPT .EQ. 1) WRITE(IOUT,*) ' CALCULATING A CHI MATRIX OF ',
     +     'SIZE ', NVECHK+1,' EQUIVALENT TO A MIXING MATRIX ',
     +    'OF SIZE ', NHK+1
C
      DO 25 I = 1,NO
      DO 20 J = 1,NO
      DO 20 K = 1,3
 20    IDIFF(K,J) = IGLIST(K,INDEX(I)) - IGLIST(K,INDEX(J))
C
C Look up the vector differences to see if present in IGLIST.
C
      CALL LOOKUP (LIST,NGTOT,IDIFF,IPG(1,I),NG1,NG2,NG3,NO,1)
C
 25   CONTINUE
C
      IF (LSTAT) THEN
C
C Calculate head element.
C
      CALL HEAD(CHI(0,0),IGLIST,INDEX,NGTOT,NO,CEV,NDIM5,NDIM8,
     +  IWGHT,LREAL,NV,NBDS,WVK,EE2,RB,INVER,LWGHT(ISPPT),LROT,
     +  A1,A2,A3,ULA)
C
      ENDIF
C
C Take each G-vector defined as the difference of two subset vectors in
C turn.
C
      DO 100 J = 1,NO
C
      DO 100 I = 1,NO
C
C If this defined G-vector is not present in IGLIST, ignore.
C
      IF (IPG(I,J) .EQ. 0) GOTO 100
C
C Calculate head and wing elements separately.
C
      IF ((I .EQ. 1) .AND. (J .EQ. 1)) GOTO 30
C
      IF (IPG(I,J) .EQ. 1) GOTO 100
C
C Determine the star of G, and hence its representative vector.
C
 30    CALL STAR(N,IPG(I,J),NTAB,NDIM2)
C
       IF (N .GT. (NHK+1)) GOTO 100
C
       LEN1 = NC/(NTAB(N + 1) - NTAB(N))
       ISTART = NC*(N - 1) + LEN1*(IPG(I,J) - NTAB(N))
C
C Determine which rotations will take G into the representative vector of
C its star.
C
       DO 95 IUNFLD = 1,LWGHT(ISPPT)
C
        DO 50 II = 1,LEN1
        IF(RLIST(ISTART+II) .EQ. LROT(IUNFLD)) THEN
                                    NR = II
                                    GOTO 60
        ENDIF
C
 50     CONTINUE
C
       GOTO 95
C
C Take each Gi-vector, defined as G is, in turn.
C
 60     CONTINUE
        DO 90 L = 1,NO
        DO 90 K = 1,NO
         IF (IPG(K,L) .EQ. 0) GOTO 90
         IF (IPG(K,L) .EQ. 1) GOTO 90
         IF (IPG(K,L) .GT. (NVECHK+1)) GOTO 90
C
        CALL STAR(NI,IPG(K,L),NTAB,NDIM2)
C
C  Calculating only the upper triangle of the susceptibility matrix.
        IF (NI .LT. N) GOTO 90
C
        CENTRE = ZERO
        IF ((I .EQ. 1) .AND. (J .EQ. 1)) THEN
          IF (LSTAT) CALL WING(CENTRE,IGLIST,INDEX,NGTOT,NO,CEV,
     +    NDIM5,NDIM8,WVK,RB,INVER(LROT(IUNFLD)),A1,A2,A3,
     +    ULA,NV,NBDS,LREAL,K,L,EE2)
C
        ELSE
C
          DO 70 IV = 1,NV
          E = EE2(IV)
          C1 = CEV(K,IV)*CONJUG(CEV(I,IV))
          DO 70 IC = (NV+1),NBDS
C
           CENTRE = CENTRE+C1*CEV(J,IC)*CONJUG(CEV(L,IC))/(EE2(IC)-E)
C
 70       CONTINUE
C
        ENDIF
C
        IF (ABS(CENTRE) .LT. 1.0E-15) GOTO 90
C
        LEN2 = NC/(NTAB(NI + 1) - NTAB(NI))
C
        IBEGIN = NC*(NI - 1) + LEN2*(IPG(K,L) - NTAB(NI))
C
C For each Gi, look for occurence of the compound rotations hence
C performing inverse rotation.
C
C
        DO 80 MR = 1,NC
C
         IF (RLIST(IBEGIN+1) .EQ. RTABLE(RLIST(ISTART+NR),MR)) THEN
           DO 75 IR = 1,NC
            IF (RLIST((NC*(NI-1)) + IR) .EQ. MR) GOTO 77
 75        CONTINUE
C
           WRITE(IOUT,*) ' HELP! FOR MR = ',MR
           CALL EXIT
C
 77        IPL = (IR-1)/LEN2 + NTAB(NI) - 1
C
           ARG = 0.0
           DO 78 IN = 1,3
 78         ARG = ARG+2*PI*(IGLIST(IN,IPG(I,J)) - IGLIST(IN,IPG(K,L)))
     +               *VT(IN,IB(RLIST(ISTART+NR)))
C
           PHASE = CMPLX(COS(ARG),SIN(ARG))
            CHI((NTAB(N)-1),IPL) = CHI((N-1),IPL) + CENTRE*PHASE/IWGHT
C           IF ((N .EQ. 1) .AND. (IPL .EQ. (NTAB(NI)-1)))
C    +        CHI((NTAB(NI)-1),(N-1)) = CONJG(CHI((N-1),IPL))
            GOTO 90
         ENDIF
C
 80   CONTINUE
C
      WRITE(IOUT,*) 'PROBLEMS.  I J K L = ',I,J,K,L
      CALL EXIT
 90   CONTINUE
C
 95   CONTINUE
C
100   CONTINUE
C
 1000 CONTINUE
C
C Return pointer to beginning of record on unit IN290.
C
      DO 2000 ISPPT = 1,4*NUMKPT
 2000  BACKSPACE(UNIT=IN290)
C
      RETURN
      END
      SUBROUTINE HEAD(CHI,IGLIST,INDEX,NGTOT,NO,CEV,NDIM5,NDIM8,
     +  IWGHT,LREAL,NV,NBDS,WVK,E,RB,INVER,NUNFLD,LROT,
     +  A1,A2,A3,ULA)
C
C  Calculates the element G=G'=0 of the polarizability matrix.
C
      COMPLEX CHI
      INTEGER V,C,LROT(48),INVER(48)
      INTEGER IGLIST(3,NGTOT),INDEX(NO)
      REAL A1(3),A2(3),A3(3),IWVK(3)
      REAL E(NDIM8),WVK(3),RB(48,3,3)
C
      LOGICAL LREAL
C
      DOUBLE PRECISION FRED
      DOUBLE PRECISION ABOHR,RYEV,RYDERG,PI,SPI
      COMMON /CONST/   ABOHR,RYEV,RYDERG,PI,SPI
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
      COMMON /MACH1/ MACTYP
C
C---------------------------------------------------------------------
C  For complex Hamiltonians
C     COMPLEX CEV(NDIM5,NDIM8),M,P(3),PROT(3),MSQR,SUM,Z,CONJUG,ZERO
C     CONJUG(Z) = CONJG(Z)
C     ZERO = (0.0,0.0)
C     IF (LREAL) WRITE(IOUT,*) '***HEAD*** LREAL =',LREAL
C----------------------------------------------------------------------
C  For real Hamiltonians
      REAL CEV(NDIM5,NDIM8),M,P(3),PROT(3),MSQR,SUM
      CONJUG(X) = X
      ZERO = 0.0
      IF (.NOT. LREAL) WRITE(IOUT,*) '***HEAD*** LREAL =',LREAL
C
C---------------------------------------------------------------------
C
      SUM = ZERO
C
      FRED = RYEV*ABOHR*4.0D0*PI/DBLE(ULA)
      FRED = FRED**2
C
      IWVK(1)= WVK(1)*A1(1) + WVK(2)*A1(2) + WVK(3)*A1(3)
      IWVK(2)= WVK(1)*A2(1) + WVK(2)*A2(2) + WVK(3)*A2(3)
      IWVK(3)= WVK(1)*A3(1) + WVK(2)*A3(2) + WVK(3)*A3(3)
C
C  Loop over valence and conduction bands
C
        DO 100 V = 1,NV
         NFRED = NV + 1
         DO 100 C = NFRED,NBDS
C        M = ZERO
         DO 10 I = 1,3
 10       P(I) = ZERO
C
         DO 50 I = 1,NO
           M = CONJUG(CEV(I,C)) * CEV(I,V)
           DO 50 J = 1,3
            R = FLOAT(IGLIST(J,INDEX(I)))
            P(J) = P(J) + M * (R + IWVK(J))
 50      CONTINUE
C
C Loop over unfolding rotations
C <v,Rot(k)|p|c,Rot(k)> = Rot[<v,k|p|c,k>]
C
        DO 100 IUNFLD = 1,NUNFLD
        NROT = INVER(LROT(IUNFLD))
         DO 60 I = 1,3
           PROT(I) = ZERO
           DO 60 J = 1,3
 60          PROT(I) = PROT(I) + RB(NROT,I,J)*P(J)
C
C  In the B-basis, to calculate in the limit q tends to zero in the (-1,1,1)
C  direction.
C
         MSQR = - PROT(1) + PROT(2) + PROT(3)
         MSQR  = CONJUG(MSQR)* MSQR
C
         ENERGY = E (C) - E (V)
C
         MSQR = MSQR/(ENERGY**3)
         SUM = SUM + MSQR
C
 100    CONTINUE
C
C
      CHI = CHI + SUM*SNGL(FRED)/IWGHT
C     WRITE(IOUT,*) ' CHI = ',CHI
C     CALL EXIT
      RETURN
      END
      SUBROUTINE WING(SUM,IGLIST,INDEX,NGTOT,NO,CEV,NDIM5,NDIM8,WVK,
     +RB,IUNFLD,A1,A2,A3,ULA,NV,NBDS,LREAL,K,L,E)
C
C  Calculates the elements of the polarizability matrix for which G=0,
C  but not G'.
C
      INTEGER V,C
      REAL E(NDIM8),IWVK(3)
      INTEGER IGLIST(3,NGTOT),INDEX(NO)
      LOGICAL LREAL
      REAL A1(3),A2(3),A3(3),WVK(3),RB(48,3,3)
      DOUBLE PRECISION FRED
C
      DOUBLE PRECISION ABOHR,RYEV,RYDERG,PI,SPI
      COMMON /CONST/   ABOHR,RYEV,RYDERG,PI,SPI
C
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
C---------------------------------------------------------------------
C For complex Hamiltonians
C     COMPLEX CEV(NDIM5,NDIM8),M,MI,P(3),PROT(3),SUM,Z,CONJUG,ZERO
C     CONJUG(Z)=CONJG(Z)
C     ZERO = (0.0,0.0)
C     IF (LREAL) WRITE(IOUT,*) ' WING***LREAL = ',LREAL
C------------------------------------------------------------------------
C For real Hamiltonians
      REAL CEV(NDIM5,NDIM8),M,MI,P(3),PROT(3),SUM
      CONJUG(X) = X
      ZERO = 0.0
      IF (.NOT. LREAL) WRITE(IOUT,*) ' WING***LREAL = ',LREAL
C------------------------------------------------------------------------
C
      SUM = ZERO
C
      FRED = RYEV*ABOHR*4*PI/ULA
CC    WRITE(IOUT,*)'*WING*FRED=',FRED
C
C
      IWVK(1)= WVK(1)*A1(1) + WVK(2)*A1(2) + WVK(3)*A1(3)
      IWVK(2)= WVK(1)*A2(1) + WVK(2)*A2(2) + WVK(3)*A2(3)
      IWVK(3)= WVK(1)*A3(1) + WVK(2)*A3(2) + WVK(3)*A3(3)
C
      DO 100 V = 1,NV
      DO 100 C = (NV+1),NBDS
C
C       M = ZERO
        DO 10 I = 1,3
 10       P(I) = ZERO
C
        DO 50 I = 1,NO
          M = CONJUG(CEV(I,C))*CEV(I,V)
          DO 50 J = 1,3
            R = FLOAT(IGLIST(J,INDEX(I)))
            P(J) = P(J) + M * (R + IWVK(J))
C
 50     CONTINUE
C
C       <v,Rot(k)|p|c,Rot(k)> = Rot[<v,k|p|c,k>]
        DO 60 I = 1,3
          PROT(I) = ZERO
          DO 60 J = 1,3
 60         PROT(I) = PROT(I) + RB(IUNFLD,I,J) * P(J)
        M = -PROT(1) + PROT(2) + PROT(3)
        MI = CONJUG(CEV(L,C))*CEV(K,V)
        ENERGY = E(C) - E(V)
        SUM = SUM - M*MI/(ENERGY**2)
C
 100  CONTINUE
C
      SUM = SUM*SNGL(FRED)
C     WRITE(IOUT,*) ' SUM = ',SUM
C
      RETURN
      END
      SUBROUTINE WHICHG(IPL,MR,RLIST,NC,NDIM12,LEN,IBEGIN)
C
C  Using the list compiled in ROTATE, this routine calculates the
C  position of the vector in IGLIST which is obtained by operating with
C  symmetry operation MR on the first vector of a star listed there.
C
      INTEGER RLIST((NDIM12+1)*48)
C
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
C
C
      DO 10 I = 1,NC
       IF (RLIST(IBEGIN + I) .EQ. MR) GOTO 20
 10   CONTINUE
C
      WRITE(IOUT,*) ' HELP! FOR MR = ',MR
      CALL EXIT
C
 20   IPL = (I-1)/LEN
C
      RETURN
      END
      SUBROUTINE STAR(M,I,NTAB,NDIM2)
C
C  Given the position of a G-vector in IGLIST, returns the number of the
C  star it belongs to.
C
      INTEGER NTAB(NDIM2)
C
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
C
      DO 10 J = 1,NDIM2
        IF (I .LT. NTAB(J)) GOTO 20
10      CONTINUE
C
      J = NDIM2
20    IF (NTAB(J) .LE. I .OR. NTAB(J) .LE. 0) THEN
        WRITE(IOUT,*) '***ERROR*** IN SUBROUTINE STAR'
        CALL EXIT
        ENDIF
C
      M = J - 1
C
      RETURN
      END
      SUBROUTINE EXPAND(CHI,NPWCHI,NDIM12,NHK,NVECHK,RLIST,NC,RTABLE,
     +NTAB,NDIM2,VT,IB,IGLIST,NGTOT)
C
C  From the elements of the polarizablity matrix which we have calculated
C  (see CALCHI), we can construct the rest using symmetry.
C
      COMPLEX CHI(0:NPWCHI,0:NPWCHI)
      INTEGER RLIST((NDIM12+1)*48),RTABLE(48,48)
      INTEGER NTAB(NDIM2),IB(48)
      REAL VT(3,48)
      INTEGER IGLIST(3,NGTOT)
      COMPLEX PHASE
C
      DOUBLE PRECISION ABOHR,RYEV,RYDERG,PI,SPI
      COMMON /CONST/   ABOHR,RYEV,RYDERG,PI,SPI
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
C
C
      DO 100 I = 1,NHK
C
      LSTAR1 = NTAB(I+2) - NTAB(I+1)
      LEN1 = NC/LSTAR1
      ISTRT1 = NC*I
      IPI = NTAB(I+1) - 1
C
       DO 100 J = I,NHK
C
       LSTAR2 = NTAB(J+2) - NTAB(J+1)
       LEN2 = NC/LSTAR2
       ITSRT2 = NC*J
C
        DO 100 II = 1,LSTAR1
C
        NRI = RLIST(ISTRT1 + LEN1*(II-1) + 1)
        IPL1 = NTAB(I+1) + II - 2
C
          ARG1 = 2 * PI * (IGLIST(1,(IPL1+1))*VT(1,IB(NRI))
     +                   + IGLIST(2,(IPL1+1))*VT(2,IB(NRI))
     +                   + IGLIST(3,(IPL1+1))*VT(3,IB(NRI)))
C
         DO 100 JJ = 1,LSTAR2
C
         NRJ = RLIST(ITSRT2 + LEN2*(JJ-1) + 1)
         IPJ = NTAB(J+1) + JJ - 2
C
         NR = RTABLE(NRI,NRJ)
C
         CALL WHICHG(IPL2,NR,RLIST,NC,NDIM12,LEN2,ITSRT2)
C
         IPL2 = IPL2 + NTAB(J+1) - 1
C
            ARG2 = 2 * PI * (IGLIST(1,(IPL2+1))*VT(1,IB(NRI))
     +                     + IGLIST(2,(IPL2+1))*VT(2,IB(NRI))
     +                     + IGLIST(3,(IPL2+1))*VT(3,IB(NRI)))
            PHASE = CMPLX(COS(ARG2-ARG1),SIN(ARG2-ARG1))
C
        CHI(IPL1,IPL2) = CHI(IPI,IPJ)*PHASE
C
 100   CONTINUE
C
      DO 200 I = 1,NVECHK
      DO 200 J = 0,(I-1)
 200    CHI(I,J) = CONJG(CHI(J,I))
C
      RETURN
      END
      SUBROUTINE CUTOFF(NO,LISTAB,NDIM5,NCUT,INDEX,NPWCHI)
C
C  We restrict the no. of stars of G-vectors included in the expansion
C  of any eigenvector in plane waves to be less than NHK, and
C  compile a list, INDEX, of those G-vectors satisfying this restriction
C  that belong to LISTAB.
C
      INTEGER LISTAB(NDIM5),INDEX(NPWCHI)
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
C
      M = NO
      K = 1
C
      DO 50 I = 1,M
C
        IF (LISTAB(I) .LT. NCUT)THEN
                                    INDEX(K) = LISTAB(I)
                                    K = K + 1
        ENDIF
C
C
 50   CONTINUE
C
      NO = K - 1
C
      RETURN
      END
      SUBROUTINE INVERT(CHI,NPWCHI,NDIM12,NVECHK,INV,IPIVOT,
     +  ISTAT,HK,IGLIST,NGTOT,NTAB,NDIM2,NHK,VT,IB,RLIST,NC)
C
C  Calculates the inverse of matrix CHI.
C
      COMPLEX CHI(0:NPWCHI,0:NPWCHI)
      COMPLEX HK(0:NDIM12,0:NDIM12)
      COMPLEX INV(0:NPWCHI)
      INTEGER IPIVOT(0:NPWCHI)
      INTEGER IGLIST(3,NGTOT),NTAB(NDIM2),IB(48),RLIST((NDIM12+1)*48)
      REAL VT(3,48)
C
      DOUBLE PRECISION ABOHR,RYEV,RYDERG,PI,SPI
      COMMON /CONST/   ABOHR,RYEV,RYDERG,PI,SPI
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
C
      CALL CGELIM (CHI,IPIVOT,(NPWCHI+1),(NVECHK+1),
     +   1.0E-12)
C
      IF (ISTAT .EQ. 0) THEN
              INV(0) = (1.0,0.0)
              DO 5 K = 1,NPWCHI
 5              INV(K) = (0.0,0.0)
              CALL CSUBST(CHI,IPIVOT,INV,(NPWCHI+1),
     +          (NVECHK+1))
              WRITE(IOUT,*) 'Static dielectric constant = ', 1.0/INV(0)
C
      ELSE
             HK(0,0) = 1.0/CHI(0,0)
             DO 40 I = 1,NHK
             J = NTAB(I+1) - 1
             DO 10 K = 0,NPWCHI
 10            INV(K) = (0.0,0.0)
             INV(J) = (1.0,0.0)
             CALL CSUBST(CHI,IPIVOT,INV,(NPWCHI+1),
     +         (NVECHK+1))
C
             DO 40 K = 1,NHK
             LMIN = NTAB(K+1)
             LMAX = NTAB(K+2) - 1
             LEN = NC/(LMAX - LMIN + 1)
             ISTART = NC * K
             DO 40 L = LMIN,LMAX
               NR = ISTART + (L - LMIN) * LEN + 1
               ARG = -2.0 * PI * (IGLIST(1,L)*VT(1,IB(RLIST(NR)))
     +                          + IGLIST(2,L)*VT(2,IB(RLIST(NR)))
     +                          + IGLIST(3,L)*VT(3,IB(RLIST(NR))))
               PHASE = CMPLX(COS(ARG),SIN(ARG))
               HK(I,K) = HK(I,K) - INV(L - 1)*PHASE
 40          CONTINUE
C
      ENDIF
C
      RETURN
      END
      SUBROUTINE CGELIM(AR,NT,NP,N,EMACH)
C
      COMPLEX AR,YRR,DUM
      DIMENSION AR(NP,N),NT(N)
      IF(N.LT.2)GOTO 15
      DO 12 II=2,N
      I=II-1
      YRR=AR(I,I)
      IN=I
      DO 4 J=II,N
      IF(CABS(YRR)-CABS(AR(J,I)))3,4,4
3      YRR=AR(J,I)
      IN=J
4      CONTINUE
      NT(I)=IN
      IF(IN-I)5,7,5
5      DO 6 J=I,N
      DUM=AR(I,J)
      AR(I,J)=AR(IN,J)
6      AR(IN,J)=DUM
7      IF(CABS(YRR)-EMACH)1,1,8
1      AR(I,I)=EMACH*EMACH
      GOTO 12
8      DO 11 J=II,N
      IF(CABS(AR(J,I))-EMACH)11,11,9
9      AR(J,I)=AR(J,I)/YRR
      DO 10 K=II,N
10      AR(J,K)=AR(J,K)-AR(I,K)*AR(J,I)
11      CONTINUE
12      CONTINUE
15      IF(CABS(AR(N,N))-EMACH)13,13,14
13      AR(N,N)=EMACH*EMACH
14      CONTINUE
      RETURN
      END
      SUBROUTINE CSUBS(AR,NT,XR,NP,N)
C
      COMPLEX AR,XR,DUM
      DIMENSION AR(NP,N),XR(N),NT(N)
      IF(N.LT.2)GOTO 18
      DO 20 II=2,N
      I=II-1
      IF(NT(I)-I)16,17,16
16      IN=NT(I)
      DUM=XR(IN)
      XR(IN)=XR(I)
      XR(I)=DUM
17      DO 19 J=II,N
      XR(J)=XR(J)-AR(J,I)*XR(I)
19      CONTINUE
20      CONTINUE
18      DO 25 II=1,N
      I=N-II+1
      IJ=I+1
      IF(I-N)21,25,21
21      DO 22 J=IJ,N
22      XR(I)=XR(I)-AR(I,J)*XR(J)
25      XR(I)=XR(I)/AR(I,I)
      RETURN
      END
      SUBROUTINE CSUBST(AR,NT,XR,NP,N)
      DIMENSION AR(NP,N),XR(N),NT(N)
      COMPLEX AR,XR,DUM
      XR(1)=XR(1)/AR(1,1)
      DO 2 I=2,N
      II=I-1
      DO 1 J=1,II
1     XR(I)=XR(I)-AR(J,I)*XR(J)
2     XR(I)=XR(I)/AR(I,I)
      DO 5 II=2,N
      I=N-II+1
      IJ=I+1
      DO 3 J=IJ,N
3     XR(I)=XR(I)-AR(J,I)*XR(J)
      IF(NT(I)-I)4,5,4
4     IN=NT(I)
      DUM=XR(I)
      XR(I)=XR(IN)
      XR(IN)=DUM
5     CONTINUE
      RETURN
      END
