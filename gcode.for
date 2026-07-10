      SUBROUTINE GCODE (IG,NG1,NG2,NG3,B1,B2,B3,I1,I2,I3,G)
C
C INPUT DATA:
C      IG ........ G-VECTOR
C      NG1,2,3 ... LIMITS OF G SPACE GRID
C      B1,2,3 .... THE BASIS VECTORS
C OUTPUT DATA:
C      I1,2,3 .... COMPONENTS OF G IN THE SYSTEM B1,B2,B3.
C      G ......... CARTESIAN COMPONENTS OF THE G-VECTOR
C
      DIMENSION B1(3),B2(3),B3(3),G(3),IG(3)
      COMMON /FILES/INPUT,IOUT,IN290,IN213,ISTORE,IUNIT7,IUNIT8,ISTRUC,
     +               IVNLKK,ISUMRY,IKPTS
C-----------------------------------------------------------------------
C
      I1 = IG(1)
      I2 = IG(2)
      I3 = IG(3)
      IF (IABS(I1) .GT. NG1) GOTO 20
      IF (IABS(I2) .GT. NG2) GOTO 20
      IF (IABS(I3) .GT. NG3) GOTO 20
C
      DO 10 J = 1,3
        G(J) = FLOAT(I1)*B1(J) + FLOAT(I2)*B2(J) + FLOAT(I3)*B3(J)
10      CONTINUE
C
      RETURN
C
20    WRITE(IOUT,*) '***ERROR*** GCODE'
      WRITE(IOUT,*) 'VECTOR OUT OF GRID'
      CALL EXIT
C
      RETURN
      END
