      SUBROUTINE GROUP1 (IOUT,A1,A2,A3,NAT,TY,NDIM9,X,B1,B2,B3,IHG,IHC, 
     +ISY,LI,NC,IB,V,F0,R)  
C               WRITTEN ON SEPTEMBER 10TH - FROM THE ACMI COMPLEX   
C (WORLTON AND WARREN, COMPUT.PHYS.COMMUN. 8,71-84 (1974))  
C (AND 3,88-117 (1972)) 
C BASIC CRYSTALLOGRAPHIC INFORMATION ABOUT A GIVEN CRYSTAL STRUCTURE
C SUBROUTINES NEEDED: PGL1,ATFTM1,ROT1,RLV3 
C INPUT DATA:   
C      IOUT ... NUMBER OF THE OUTPUT UNIT FOR ON-LINE PRINTING  
C               OF VARIOUS MESSAGES 
C      A1,A2,A3 .. ELEMENTARY TRANSLATIONS OF THE LATTICE, IN SOME  
C               UNIT OF LENGTH  
C      NAT .... NUMBER OF ATOMS IN THE UNIT CELL
C      TY ..... INTEGERS DISTINGUISHING BETWEEN THE ATOMS OF
C               DIFFERENT TYPE. TY(I) IS THE TYPE OF THE I-TH ATOM  
C               OF THE BASIS
C      X ...... CARTESIAN COORDINATES OF THE NAT ATOMS OF THE BASIS 
C OUTPUT DATA   
C      B1,B2,B3 .. RECIPROCAL LATTICE VECTORS, NOT MULTIPLIED BY
C               ANY 2PI, IN UNITS RECIPROCAL TO THOSE OF A1,A2,A3   
C      IHG .... POINT GROUP OF THE PRIMITIVE LATTICE, HOLOEDRAL 
C               GROUP NUMBER:   
C               IHG=1 STANDS FOR TRICLINIC SYSTEM   
C               IHG=2 STANDS FOR MONOCLINIC SYSTEM  
C               IHG=3 STANDS FOR ORTHORHOMBIC SYSTEM
C               IHG=4 STANDS FOR TETRAGONAL SYSTEM  
C               IHG=5 STANDS FOR CUBIC SYSTEM   
C               IHG=6 STANDS FOR TRIGONAL SYSTEM
C               IHG=7 STANDS FOR HEXAGONAL SYSTEM   
C      IHC .... CODE DISTINGUISHING BETWEEN HEXAGONAL AND CUBIC 
C               GROUPS  
C               IHC=0 STANDS FOR HEXAGONAL GROUPS   
C               IHC=1 STANDS FOR CUBIC GROUPS   
C      ISY .... CODE INDICATING WHETHER THE SPACE GROUP IS  
C               SYMMORPHIC OR NONSYMMORPHIC 
C               ISY=0 MEANS NONSYMMORPHIC GROUP 
C               ISY=1 MEANS SYMMORPHIC GROUP
C               THE GROUP IS CONSIDERED SYMMORPHIC IF FOR EACH  
C               OPERATION OF THE POINT GROUP THE SUM OF THE 3   
C               COMPONENTS OF ABS(V(N)) ( NONPRIMITIVE TRANSLATION, 
C               SEE BELOW) IS LT. 0.0005
C      LI ..... CODE INDICATING WHETHER THE POINT GROUP 
C               OF THE CRYSTAL CONTAINS INVERSION OR NOT
C               (OPERATIONS 13 OR 25 IN RESPECTIVELY HEXAGONAL  
C               OR CUBIC GROUPS).   
C               LI=0 MEANS: DOES NOT CONTAIN INVERSION  
C               LI.GT.0 MEANS: THERE IS INVERSION IN THE POINT  
C                              GROUP OF THE CRYSTAL 
C      NC ..... TOTAL NUMBER OF ELEMENTS IN THE POINT GROUP OF THE  
C               CRYSTAL 
C      IB ..... LIST OF THE ROTATIONS CONSTITUTING THE POINT GROUP  
C               OF THE CRYSTAL. THE NUMBERING IS THAT DEFINED IN
C               WORLTON AND WARREN, I.E. THE ONE MATERIALIZED IN THE
C               ARRAY R (SEE BELOW) 
C               ONLY THE FIRST NC ELEMENTS OF THE ARAY IB ARE   
C               MEANINGFUL  
C      V ...... NONPRIMITIVE TRANSLATIONS (IN THE CASE OF NONSYMMOR-
C               PHIC GROUPS). V(I,N) IS THE I-TH COMPONENT  
C               OF THE TRANSLATION CONNECTED WITH THE N-TH ELEMENT  
C               OF THE POINT GROUP (I.E. WITH THE ROTATION  
C               NUMBER IB(N) ). 
C               ATTENTION: V(I) ARE NOT CARTESIAN COMPONENTS,   
C               THEY REFER TO THE SYSTEM A1,A2,A3.  
C      F0 ..... THE FUNCTION DEFINED IN MARADUDIN, IPATOVA BY   
C               EQ. (3.2.12): ATOM TRANSFORMATION TABLE.
C               THE ELEMENT F0(N,KAPA) MEANS THAT THE N-TH  
C               OPERATION OF THE SPACE GROUP (I.E. OPERATION NUMBER 
C               IB(N), TOGETHER WITH AN EVENTUAL NONPRIMITIVE   
C               TRANSLATION  V(N)) TRANSFERS THE ATOM KAPA INTO THE 
C               ATOM F0(N,KAPA).
C      R ...... LIST OF THE 3 X 3 ROTATION MATRICES 
C               (XYZ REPRESENTATION OF THE O(H) OR D(6)H GROUPS)
C               ALL 48 OR 24 MATRICES ARE LISTED.   
C PRINTED OUTPUT
C      PROGRAM PRINTS THE TYPE OF THE LATTICE (IHG, IN WORDS),  
C      LISTS THE OPERATIONS OF THE  POINT GROUP OF THE  
C      CRYSTAL, INDICATES WHETHER THE SPACE GROUP IS SYMMORPHIC OR  
C      NONSYMMORPHIC AND WHETHER THE POINT GROUP OF THE CRYSTAL 
C      CONTAINS INVERSION.  
C   
      INTEGER TY(NDIM9),F0(48,NDIM9)
      DIMENSION A1(3),A2(3),A3(3),X(3,NDIM9),B1(3),B2(3),B3(3), 
     +IB(48),V(3,48),R(49,3,3)  
      DIMENSION A(3,3), AI(3,3) 
C-----------------------------------------------------------------------
      IF (IOUT .GT. 0) WRITE (IOUT,*) 'Subroutine GROUP1 - beginning'
      IF (NAT .GT. NDIM9) THEN  
        WRITE (IOUT,120) NAT,NDIM9  
120     FORMAT ('0*** NAT =',I5,' is greater than NDIM9 =',I5)  
        CALL EXIT   
        ENDIF   
C   
      DO 140 I = 1,3
        A(I,1) = A1(I)  
        A(I,2) = A2(I)  
        A(I,3) = A3(I)  
140     CONTINUE
C     A(I,J) IS THE I-TH CARTESIAN COMPONENT OF THE J-TH PRIMITIVE  
C     TRANSLATION VECTOR OF THE DIRECT LATTICE  
C     TY(I) IS AN INTEGER DISTINGUISHING ATOMS OF DIFFERENT TYPE,I.E.,  
C     DIFFERENT ATOMIC SPECIES  
C     X(J,I) IS THE J-TH CARTESIAN COMPONENT OF THE POSITION VECTOR FOR 
C     THE I-TH ATOM IN THE UNIT CELL.   
C   
C     DETERMINE PRIMITIVE LATTICE VECTORS FOR THE RECIPROCAL LATTICE.   
C   
      DET = A(1,1)*A(2,2)*A(3,3) + A(2,1)*A(1,3)*A(3,2) +   
     +      A(3,1)*A(1,2)*A(2,3) - A(1,1)*A(2,3)*A(3,2) -   
     +      A(2,1)*A(1,2)*A(3,3) - A(3,1)*A(1,3)*A(2,2) 
      DET = 1.0/DET 
      DO 150 I = 1,3
        IL = 1  
        IU = 3  
        IF (I .EQ. 1) IL = 2
        IF (I .EQ. 3) IU = 2
        DO 150 J = 1,3  
          JL = 1
          JU = 3
          IF (J .EQ. 1) JL = 2  
          IF (J .EQ. 3) JU = 2  
          AI(J,I) = (-1.0)**(I+J) * DET *   
     +    ( A(IL,JL) * A(IU,JU) - A(IL,JU) * A(IU,JL) ) 
150       CONTINUE  
      DO 160 I = 1,3
        B1(I)  = AI(1,I)
        B2(I)  = AI(2,I)
        B3(I)  = AI(3,I)
160     CONTINUE
      CALL PGL1 (A,AI,IHC,NC,IB,IHG,R)  
      CALL ATFTM1 (IOUT,R,V,X,F0,IB,TY,NDIM9,IHG,NAT,NC,AI,LI,ISY)  
      IF (LI .GT. 0 .AND. IOUT .GT. 0) WRITE (IOUT,*) 
     +  'The point group of the crystal contains the inversion'   
      IF (IOUT .GT. 0) WRITE (IOUT,*) 'Subroutine GROUP1 - end'  
      RETURN
      END   
      SUBROUTINE PGL1 (A,AI,IHC,NC,IB,IHG,R)
C                WRITTEN ON SEPTEMBER 11TH, 1979 - FROM ACMI COMPLEX
C AUXILIARY SUBROUTINE TO GROUP1
C     SUBROUTINE PGL DETERMINES THE POINT GROUP OF THE LATTICE AND THE  
C     CRYSTAL SYSTEM.   
C SUBROUTINES NEEDED: ROT1, RLV3
C A ..... DIRECT LATTICE VECTORS
C AI .... RECIPROCAL LATTICE VECTORS
      DIMENSION R(49,3,3), IB(48), A(3,3), AI(3,3)  
      DIMENSION VR(3), XA(3)
      DOUBLE PRECISION DXA  
C   
C.....FILES 
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
C.....PHYSICAL AND MATHEMATICAL CONSTANTS   
      DOUBLE PRECISION ABOHR,RYEV,RYDERG,PI,SPI 
      COMMON /CONST/   ABOHR,RYEV,RYDERG,PI,SPI 
      PARAMETER ( EPS = 1.0E-5 )
C   
C-----------------------------------------------------------------------
      IHC = 0   
C     IHC IS 0 FOR HEXAGONAL GROUPS AND 1 FOR CUBIC GROUPS. 
      NR = 24   
100   NC = 0
      CALL ROT1 (IHC,R) 
      DO 140 N = 1,NR   
        IB(N) = 0   
        TR = 0.0
C       ROTATE THE A1,2,3 VECTORS BY ROTATION NO. N 
        DO 130 K = 1,3  
          DO 115 I = 1,3
            DXA = 0.0D0 
            DO 110 J = 1,3  
              DXA = DXA + DPROD( R(N,I,J) , A(J,K) )
110           CONTINUE  
            XA(I) = DXA 
115           CONTINUE  
          CALL RLV3 (AI,XA,VR,LX)   
          DO 120 I = 1,3
            TR = TR + ABS(VR(I))
120         CONTINUE
C         IF VR.NE.0, THEN XA CANNOT BE A MULTIPLE OF A LATTICE VECTOR  
          IF (TR .GT. EPS) GOTO 140 
130       CONTINUE  
C   
        NC = NC + 1 
        IB(NC) = N  
140     CONTINUE
C   
C     IHG STANDS FOR HOLOHEDRAL GROUP NUMBER.   
      IF (IHC .EQ. 0)  GO TO 200
C     CUBIC GROUP:  
      IF (NC  .LT. 4)  IHG = 1  
      IF (NC  .EQ. 4)  IHG = 2  
      IF (NC  .GT. 4)  IHG = 3  
      IF (NC  .EQ. 16) IHG = 4  
      IF (NC  .GT. 16) IHG = 5  
      RETURN
C   
C     HEXAGONAL GROUP:  
200   IF (NC  .EQ. 12) IHG = 6  
      IF (NC  .GT. 12) IHG = 7  
      IF (NC  .GE. 12) RETURN   
C     TOO FEW OPERATIONS, TRY CUBIC GROUP:  
      NR  = 48  
      IHC = 1   
      GOTO 100  
      END   
      SUBROUTINE RLV3 (AI,XB,VR,IL) 
      DIMENSION AI(3,3),XB(3),VR(3) 
C     WRITTEN ON SEPTEMBER 11TH, 1979 - FROM ACMI COMPLEX   
C     OPTIMIZED 24-MAY-1986 BY OLE H. NIELSEN   
C     AUXILIARY SUBROUTINE TO GROUP1
C     SUBROUTINE RLV REMOVES A DIRECT LATTICE VECTOR FROM XB LEAVING THE
C     REMAINDER IN VR.  IF A NONZERO LATTICE VECTOR WAS REMOVED, IL IS  
C     MADE NONZERO.  VR STANDS FOR V-REFERENCE. 
C     AI(I,J) ARE THE RECIPROCAL LATTICE VECTORS, B(I) = AI(I,J),J=1,2,3
C     VR IS NOT GIVEN IN CARTESIAN COORDINATES BUT  
C     IN THE SYSTEM A1,A2,A3.     K.K., 23.10.1979  
CDIR$ INT24 I   
      PARAMETER ( EPS = 1.0E-5 )
C-----------------------------------------------------------------------
      IF (ABS(XB(1))+ABS(XB(2))+ABS(XB(3)) .LE. EPS) THEN   
        VR(1) = 0.0 
        VR(2) = 0.0 
        VR(3) = 0.0 
      ELSE  
        IL = 0  
CDIR$   NOVECTOR
        DO 120 I = 1,3  
          DVR = DPROD(AI(I,1),XB(1)) + DPROD(AI(I,2),XB(2)) +   
     +          DPROD(AI(I,3),XB(3))
          IL  = IL + ABS(DVR)   
          VR(I) = NINT(DVR) - DVR   
120       CONTINUE  
          ENDIF 
      RETURN
      END   
      SUBROUTINE ATFTM1(IOUT,R,V,X,F0,IB,TY,NDIM9,IHG,NAT,NC,AI,LI,ISY) 
C           WRITTEN ON SEPTEMBER 11TH, 1979 - FROM ACMI COMPLEX 
C     AUXILIARY SUBROUTINE TO GROUP1
C     SUBROUTINE ATFTMT DETERMINES THE POINT GROUP OF THE CRYSTAL,  
C     THE ATOM TRANSFORMATION TABLE,F0, THE FRACTIONAL TRANSLATIONS,V,  
C     ASSOCIATED WITH EACH ROTATION.
C SUBROUTINES NEEDED: RLV3  
      DIMENSION  IIS(48), RX(3,50), VR(3), VT(3), XB(3) 
      DIMENSION R(49,3,3), V(3,48), X(3,NDIM9),IB(48),AI(3,3)   
      INTEGER F0(48,NDIM9),TY(NDIM9)
      PARAMETER ( EPS = 1.0E-5 )
      PARAMETER ( EPSIL = 0.0005 )  
C THE ORIGINAL VALUE EPSIL=0.0001 WAS MODIFIED  
C BY K.K. , SEPTEMBER 1979  
C   
      CHARACTER*12 ICST(7)  
      DATA ICST /'TRICLINIC','MONOCLINIC','ORTHORHOMBIC','TETRAGONAL',  
     +'CUBIC','TRIGONAL','HEXAGONAL'/   
C-----------------------------------------------------------------------
      IF (NDIM9 .LE. 50) GOTO 90
      WRITE (IOUT,80) NDIM9 
80    FORMAT('0Subroutine ATFTM1 *** FATAL ERROR *** NDIM9=',I4,'>50')  
      CALL EXIT 
C   
90    NCA = 0   
      DO 100 N = 1,48   
        IIS(N) = 0  
100     CONTINUE
C   
      DO 210 N = 1,NC   
C   
        L = IB(N)   
        IIS(L) = 1  
        DO 120 K = 1,NAT
          DO 120 I = 1,3
            RX(I,K) = 0.0   
            DO 120 J = 1,3  
              RX(I,K)=RX(I,K)+R(L,I,J)*X(J,K)   
120           CONTINUE  
        DO 180 K1 = 1,NAT   
        DO 180 K2 = 1,NAT   
          IF (TY(K1) .NE. TY(K2)) GO TO 180 
          DO 130 I=1,3  
            XB(I)=RX(I,K1)-X(I,K2)  
130         CONTINUE
          IL = 0
          CALL RLV3 (AI,XB,VR,IL)   
C     SUBROUTINE RLV REMOVES A DIRECT LATTICE VECTOR FROM XB LEAVING THE
C     REMAINDER IN VR.  IF A NONZERO LATTICE VECTOR WAS REMOVED, IL IS  
C     MADE NONZERO.  VR STANDS FOR V-REFERENCE. 
C VR IS NOT GIVEN IN CARTESIAN COORDINATES BUT  
C IN THE SYSTEM A1,A2,A3.     K.K., 23.10.1979  
          KS = 0
          DO 170 K3 = 1,NAT 
          DO 160 K4 = 1,NAT 
            IF (TY(K3) .NE. TY(K4)) GO TO 160   
CDIR$       NOVECTOR
            DO 140 I = 1,3  
              XB(I) = RX(I,K3) - X(I,K4)
140           CONTINUE  
            CALL RLV3 (AI,XB,VT,IL) 
C           VT STANDS FOR V-TEST
            DIF = 0.0   
CDIR$       NOVECTOR
            DO 150 I = 1,3  
              DA = VR(I) - VT(I)
              DIF = DIF + ABS(DA - NINT(DA))
150           CONTINUE  
            IF (DIF .GT. EPS) GO TO 160 
            F0(N,K3)=K4 
C     F0 IS THE FUNCTION DEFINED IN MARADUDIN AND VOSK0 BY EQ.(2.35).   
C     IT DEFINES THE ATOM TRANSFORMATION TABLE  
            KS = KS + K4
            IF (KS .EQ. NAT*(NAT+1)/2 ) GO TO 190   
            GO TO 170   
160         CONTINUE
          GO TO 180 
170       CONTINUE  
180     CONTINUE
      IIS(L) = 0
      GO TO 210 
190   NCA = NCA + 1 
      DO 200 I = 1,3
        V(I,NCA)=VR(I)  
200     CONTINUE
C     V(I,N) IS THE I-TH COMPONENT OF THE FRACTIONAL
C     TRANSLATION ASSOCIATED WITH THE ROTATION N.   
C ATTENTION: V(I) ARE NOT CARTESIAN COMPONENTS, THEY ARE
C GIVEN IN THE SYSTEM A1,A2,A3.     K.K., 23. 10. 1979  
C   
210     CONTINUE
C   
      I  = 0
      NI = 13   
      IF (IHG .LT. 6) NI = 25   
      LI = 0
      DO 230 N = 1,NC   
        L = IB(N)   
        IF (IIS(L) .EQ. 0) GO TO 230
        I = I + 1   
        IB(I) = IB(N)   
        IF (IB(I) .EQ. NI) LI = I   
        DO 220 K = 1,NAT
          F0(I,K) = F0(N,K) 
220       CONTINUE  
230     CONTINUE
C   
      NC = I
      IF ( (IHG .EQ. 7 .AND. NC .EQ. 24) .OR.   
     +     (IHG .EQ. 5 .AND. NC .EQ. 48) ) THEN 
        IF (IOUT .GT. 0) WRITE (IOUT,260) ICST(IHG)  
260     FORMAT (' The point group of the crystal is the full ', 
     +          A,' group') 
      ELSE  
        IF (IOUT .GT. 0) WRITE (IOUT,240) ICST(IHG),(IB(I),I=1,NC)   
240     FORMAT (' The crystal system is ',A,' with operations:'/
     +          2(3X,24I3/))
        ENDIF   
C   
      VS = 0.0  
      DO 280 N = 1,NC   
      DO 280 I = 1,3
        VS = VS + ABS(V(I,N))   
280     CONTINUE
C   
      IF (VS .LT. EPS*FLOAT(NC)) THEN   
        IF (IOUT .GT. 0)
     +  WRITE (IOUT,*) 'The space group of the crystal is symmorphic' 
        ISY = 1 
      ELSE  
        IF (IOUT .GT. 0)
     +  WRITE (IOUT,310)
310     FORMAT (' The space group is non-symmorphic,'/  
     +       ' or else a non standard origin of coordinates was used.') 
        ISY = 0 
        ENDIF   
      RETURN
      END   
      SUBROUTINE ROT1 (IHC,R)   
C                                     WRITTEN ON FEBRUARY 17TH, 1976
C GENERATION OF THE X,Y,Z-TRANSFORMATION MATRICES 3X3 FOR HEXAGONAL AND 
C CUBIC GROUPS  
C SUBROUTINES NEEDED -- NONE
C THIS IS IDENTICAL WITH THE SUBROUTINE ROT OF WORLTON-WARREN (IN THE AC
C -COMPLEX), ONLY THE WAY OF TRANSFERRING THE DATA WAS CHANGED  
C INPUT DATA
C      IHC...SWITCH DETERMINING IF WE DESIRE THE HEXAGONAL GROUP (IHC=0)
C            OR THE CUBIC GROUP (IHC=1) 
C OUTPUT DATA   
C      R...THE 3X3 MATRICES OF THE DESIRED COORDINATE REPRESENTATION
C          THEIR NUMBERING CORRESPONDS TO THE SYMMETRY ELEMENTS AS LISTE
C          IN WORLTON-WARREN
C          FOR IHC=0 THE FIRST 24 MATRICES OF THE ARRAY R REPRESENT 
C                                          THE FULL HEXAGONAL GROUP D(6H
C          FOR IHC=1 THE FIRST 48 MATRICES OF THE ARRAY R REPRESENT 
C                                          THE FULL CUBIC GROUP O(H)
C   
      DIMENSION R(49,3,3)   
C-----------------------------------------------------------------------
      DO 100 N=1,49 
      DO 100 I=1,3  
      DO 100 J=1,3  
        R(N,I,J)=0.0
100     CONTINUE
      IF (IHC .GT. 0) GOTO 160  
C   
C     DEFINE THE GENERATORS FOR THE ROTATION MATRICES--HEXAGONAL GROUP  
C   
      F = 0.5D0*DSQRT(3.0D0)
      R(2,1,1) = 0.5
      R(2,1,2) = - F
      R(2,2,1) = F  
      R(2,2,2) = 0.5
      R(7,1,1) = - 0.5  
      R(7,1,2) = - F
      R(7,2,1) = - F
      R(7,2,2) = 0.5
      DO 120 N      = 1,6   
        R(N,3,3)    = 1.0   
        R(N+18,3,3) = 1.0   
        R(N+6,3,3)  = - 1.0 
        R(N+12,3,3) = - 1.0 
120     CONTINUE
C   
C     GENERATE THE REST OF THE ROTATION MATRICES
C   
      DO 130 I = 1,2
        R(1,I,I) = 1.0  
        DO 130 J = 1,2  
          R(6,I,J) = R(2,J,I)   
          DO 130 K = 1,2
            R(3,I,J)  = R(3,I,J) +  R(2,I,K)*R(2,K,J)   
            R(8,I,J)  = R(8,I,J) +  R(2,I,K)*R(7,K,J)   
            R(12,I,J) = R(12,I,J) + R(7,I,K)*R(2,K,J)   
130         CONTINUE
      DO 140 I = 1,2
        DO 140 J = 1,2  
          R(5,I,J) = R(3,J,I)   
          DO 140 K = 1,2
            R(4,I,J)  = R(4,I,J)  + R(2,I,K)*R(3,K,J)   
            R(9,I,J)  = R(9,I,J)  + R(2,I,K)*R(8,K,J)   
            R(10,I,J) = R(10,I,J) + R(12,I,K)*R(3,K,J)  
            R(11,I,J) = R(11,I,J) + R(12,I,K)*R(2,K,J)  
140         CONTINUE
      DO 150 N = 1,12   
        NV = N + 12 
        DO 150 I = 1,2  
          DO 150 J = 1,2
            R(NV,I,J) = - R(N,I,J)  
150         CONTINUE
C   
      RETURN
C   
C     DEFINE THE GENERATORS FOR THE ROTATION MATRICES--CUBIC GROUP  
C   
160   R(9,1,3)  = 1.0   
      R(9,2,1)  = 1.0   
      R(9,3,2)  = 1.0   
      R(19,1,1) = 1.0   
      R(19,2,3) = - 1.0 
      R(19,3,2) = 1.0   
      DO 170 I = 1,3
        R(1,I,I) = 1.0  
        DO 170 J = 1,3  
          R(20,I,J) = R(19,J,I) 
          R(5,I,J)  = R(9,J,I)  
          DO 170 K  = 1,3   
            R(2,I,J)  = R(2,I,J)  + R(19,I,K)*R(19,K,J) 
            R(16,I,J) = R(16,I,J) + R(9,I,K)*R(19,K,J)  
            R(23,I,J) = R(23,I,J) + R(19,I,K)*R(9,K,J)  
170         CONTINUE
      DO 180 I = 1,3
      DO 180 J = 1,3
      DO 180 K = 1,3
        R(6,I,J)  = R(6,I,J)  + R(2,I,K)*R(5,K,J)   
        R(7,I,J)  = R(7,I,J)  + R(16,I,K)*R(23,K,J) 
        R(8,I,J)  = R(8,I,J)  + R(5,I,K)*R(2,K,J)   
        R(10,I,J) = R(10,I,J) + R(2,I,K)*R(9,K,J)   
        R(11,I,J) = R(11,I,J) + R(9,I,K)*R(2,K,J)   
        R(12,I,J) = R(12,I,J) + R(23,I,K)*R(16,K,J) 
        R(14,I,J) = R(14,I,J) + R(16,I,K)*R(2,K,J)  
        R(15,I,J) = R(15,I,J) + R(2,I,K)*R(16,K,J)  
        R(22,I,J) = R(22,I,J) + R(23,I,K)*R(2,K,J)  
        R(24,I,J) = R(24,I,J) + R(2,I,K)*R(23,K,J)  
180     CONTINUE
      DO 190 I=1,3  
      DO 190 J=1,3  
      DO 190 K=1,3  
        R(3,I,J)  = R(3,I,J)  + R(5,I,K)*R(12,K,J)  
        R(4,I,J)  = R(4,I,J)  + R(5,I,K)*R(10,K,J)  
        R(13,I,J) = R(13,I,J) + R(23,I,K)*R(11,K,J) 
        R(17,I,J) = R(17,I,J) + R(16,I,K)*R(12,K,J) 
        R(18,I,J) = R(18,I,J) + R(16,I,K)*R(10,K,J) 
        R(21,I,J) = R(21,I,J) + R(12,I,K)*R(15,K,J) 
190     CONTINUE
      DO 200 N = 1,24   
        NV = N + 24 
        DO 200 I = 1,3  
        DO 200 J = 1,3  
          R(NV,I,J) = - R(N,I,J)
200       CONTINUE  
      RETURN
      END   
