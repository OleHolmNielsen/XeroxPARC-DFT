      SUBROUTINE DAY(IYEAR,IMONTH,IDAY)
C
C     RETURN YEAR,MONTH,DAY AND MACHINE AS INTEGERS (2 DIGITS)
C
      CHARACTER*8 IC
C     For IDATE (VAX and UNIX):
      INTEGER IARRAY(3)
      EQUIVALENCE (IARRAY,IA)
C     Use machine names in all programs:
      COMMON /MACH1/ MACTYP
      COMMON /MACH2/ MACHIN(0:20)
      CHARACTER*7  MACHIN
C
C-----------------------------------------------------------------------
C
      IF (MACTYP .LE. 0 .OR. MACTYP .GT. 20) CALL MACH
C
      IF (MACTYP .EQ. 2 .OR. MACTYP .EQ. 7) THEN
C
C       CRAY: DATE RETURNS IC='MM/DD/YY' (TYPE INTEGER).
        CALL DATE(IC)
        READ (IC,'(I2,1X,I2,1X,I2)') IMONTH,IDAY,IYEAR
C
      ELSE IF (MACTYP .EQ. 3) THEN
C
C       VAX:
        CALL IDATE (IMONTH,IDAY,IYEAR)
C
      ELSE IF (MACTYP .EQ. 10) THEN
C
C       Sun UNIX:
        CALL IDATE (IA,IDUMMY,IDUMMY)
        IDAY   = IARRAY(1)
        IMONTH = IARRAY(2)
        IYEAR  = IARRAY(3) - 1900
C
      ELSE IF (MACTYP .EQ. 11) THEN
C
C       Apollo UNIX: Front-end to cal_$decode_local_time:
        CALL IDATE (IMONTH,IDAY,IYEAR)
C
        ENDIF
C
      RETURN
      END
      SUBROUTINE DAYPRT (IOUT,IYEAR,IMONTH,IDAY)
C
C     PRINT TODAY'S DATE
C
      CHARACTER*9 MONTHS(12)
      CHARACTER*2 DAYS(31)
      DATA MONTHS /'January', 'February', 'March', 'April',
     +  'May', 'June', 'July', 'August', 'September', 'October',
     +  'November', 'December'/
      DATA DAYS /'st', 'nd', 'rd', 17*'th', 'st', 'nd', 'rd',
     +           7*'th', 'st'/
C
      if (iyear .lt. 80 .or. iyear .gt. 99
     +  .or. imonth .le. 0 .or. imonth .gt. 12
     +  .or. iday .le. 0 .or. iday .gt. 31) then
        write (iout,*) 'DAYPRT: bad arguments=',IYEAR,IMONTH,IDAY
        return
        endif
      MM = LENSTR(MONTHS(IMONTH))
      if (mm .le. 0 .or. mm .gt. 20) then
        write (iout,*) 'DAY: bad MM=', MM
        return
        endif
C
      WRITE (IOUT,100) IDAY, DAYS(IDAY), 
     +  MONTHS(IMONTH)(1:MM), IYEAR
100   FORMAT (T46,'Date: ',I2,A,' of ',A,', 19',I2.2)
      RETURN
      END
c     subroutine IDATE (IMONTH,IDAY,IYEAR)
c     RETURN
c     END
