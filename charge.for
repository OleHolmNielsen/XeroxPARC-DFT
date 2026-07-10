      SUBROUTINE CHARGE(NAME, Z)
      CHARACTER*2 NAME
      DOUBLE PRECISION Z
C
C     FUNCTION DETERMINES THE NUCLEAR CHARGE OF AN ELEMENT
C
      PARAMETER ( NELEM = 54 )
      CHARACTER*2 ELEMNT, PERTAB(NELEM)
      INTEGER IC(2)
      SAVE PERTAB, IOUT
      DATA IOUT /6/
C     THE PERIODIC TABLE
      DATA PERTAB /
     +  'H ','HE',
     +  'LI','BE','B ','C ','N ','O ','F ','NE',
     +  'NA','MG','AL','SI','P ','S ','CL','AR',
     +  'K ','CA',
     +       'SC','TI','V ','CR','MN','FE','CO','NI','CU','ZN',
     +            'GA','GE','AS','SE','BR','KR',
     +  'RB','SR',
     +       'Y ','ZR','NB','MO','TC','RU','RH','PD','AG','CD',
     +            'IN','SN','SB','TE','I ','XE'/
C
C
C     CONVERT THE NAME TO UPPER-CASE, AND POSSIBLY LEFT-JUSTIFY
C
C     CODE 97-122: LOWER CASE
C     CODE 65-90:  UPPER CASE
C     CODE 32:     BLANK
C
      DO 100 I = 1,2
C     GET THE ASCII VALUE
      IC(I) = ICHAR( NAME(I:I) )
      IF (IC(I) .GE. 97 .AND. IC(I) .LE. 122) THEN
C       CONVERT TO UPPER CASE
        IC(I) = IC(I) - 32
      ELSE IF (IC(I) .GE. 65 .AND. IC(I) .LE. 90) THEN
C       UPPER-CASE - DO NOTHING
      ELSE IF (IC(I) .EQ. 32) THEN
C       'SPACE' - DO NOTHING
      ELSE IF (IC(I) .EQ. 0) THEN
C       'NUL' - REPLACE BY SPACE
        IC(I) = 32
      ELSE
        WRITE (IOUT,*) 'Unrecognized element name:',NAME
        CALL EXIT
        ENDIF
100   CONTINUE
C
C     LEFT JUSTIFY
      IF (IC(1) .EQ. 32) THEN
        IC(1) = IC(2)
        IC(2) = 32
        ENDIF
C     THE STANDARD NAME OF THE ELEMENT:
      ELEMNT = CHAR(IC(1))//CHAR(IC(2))
C
C     FIND THE ELEMENT IN THE PERIODIC TABLE
C
      DO 150 I = 1, NELEM
        IF (ELEMNT .EQ. PERTAB(I)) THEN
          Z = I
          RETURN
          ENDIF
150     CONTINUE
      WRITE (IOUT,160) NAME,ELEMNT,IC
160   FORMAT (' Could not locate name in list of elements:'/
     + ' Name=',A,' converted to=',A,' ASCII codes=',2I3)
      CALL EXIT
      RETURN
C
C     ENTRY FOR ELEMENT NAMES
C
      ENTRY ELTNAM (NAME,NATOM)
C
      IF (NATOM .LT. 1 .OR. NATOM .GT. NELEM) THEN
        WRITE (IOUT,*)
     +   'Function ELTNAM: I do not know the name of element ',NATOM
        NAME = '  '
      ELSE
        NAME = PERTAB(NATOM)
        CALL LOWER(NAME(2:2))
        ENDIF
      RETURN
      END
      SUBROUTINE LOWER (CH)
C     Convert upper-case letter to lower-case
      CHARACTER CH
      ICH = ICHAR(CH)
      IF (ICH .GE. 65 .AND. ICH .LE. 90) ICH = ICH + 32
      CH = CHAR(ICH)
      RETURN
      END
