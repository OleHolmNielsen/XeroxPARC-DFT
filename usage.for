C      FUNCTION SECOND (DUMMY)
C      CALL USAGE (N, SECOND)
C      RETURN
C      END
      SUBROUTINE USAGE (NPAGES,CPUTIME)
C
C     Number of page faults (virtual machines, only) and CPU-time so far
C
C
C     Definition of machine types:
C     0 ..... Unknown
C     1 ..... IBM
C     2 ..... CRAY-1 or X-MP
C     3 ..... VAX
C     4 ..... CDC
C     5 ..... UNIVAC = Sperry = Unisys
C     6 ..... VAX + FPS attached processor
C     7 ..... CRAY-2
C     8 ..... Fujitsu VP-100/200/400 (= Amdahl 1XXX = Siemens)
C     9 ..... VAX + MAP attached processor
C     10 .... Generic UNIX (tested on SunOS3.4 and 4.0)
C     11 .... Apollo UNIX
C     12 .... Linux x86
C
C     Use machine names in all programs:
      COMMON /MACH1/ MACTYP
      COMMON /MACH2/ MACHIN(0:20)
      CHARACTER*7  MACHIN
      REAL A(2)
C-----------------------------------------------------------------------
C
C     Just to make sure that MACHIN is properly set:
      CALL DAY (IDAY,IMONTH,IYEAR)
C
      NPAGES = 0
      CPUTIME = 0.0
#ifdef x86_64
#else

C
C-----------------------------------------------------------------------
C
      IF (MACTYP .EQ. 1) THEN
C
C       IBM: Do nothing (don't know proper call)
C
C-----------------------------------------------------------------------
C
C       CRAY 1, X-MP or 2: (CFT manual p. 7-27)
C
      ELSE IF (MACTYP .EQ. 2 .OR. MACTYP .EQ. 7) THEN
C
        CPUTIME = SECOND(DUMMY)
C
C-----------------------------------------------------------------------
C
C     VAX:
C
      ELSE IF (MACTYP .EQ. 3) THEN
C
        CALL USGVAX(NPAGES,CPUTIME)
C
C-----------------------------------------------------------------------
C
C     Generic UNIX:
C
      ELSE IF (MACTYP .EQ. 10) THEN
C
        CALL UNIXUS (NPAGES,CPUTIME)
C       CALL ETIME(A)
C       CPUTIME = A(1) + A(2)
C
C-----------------------------------------------------------------------
C
C     Generic UNIX:
C
      ELSE IF (MACTYP .EQ. 11) THEN
C
        CALL APOLLO (NPAGES,CPUTIME)
C
C-----------------------------------------------------------------------
C
      ELSE
C
C       Do nothing
C
        ENDIF
C
C-----------------------------------------------------------------------
C
#endif
      RETURN
      END

#ifdef vax
C
C     Uncomment on VAX machines:
C     include 'usagevax.for/list'
#endif

#ifdef apollo
C
C     On Apollo UNIX machines, replace the following dummy code with the contents
C     of the file usapollo.for - NB: Include does NOT work for entire subroutines
C     in Apollo's fortran (how stupid can you become ....?)
C
      subroutine apollo (npages, cpu)
c
c     CPU tine and page faults on Apollo
c
      real*8 sect
      integer*2 clock(3)
      npages = 0
      return
      end
c     subroutine idate(imonth, iday, iyear)
c
c     Date subroutine for Apollo
c
c     return
c     end
#endif

#ifdef aix_xlf
      subroutine idate(imonth, iday, iyear)
 
c     Date subroutine for AIX XLF 3.X
 
        type iar
                sequence
                integer*4 iday
                integer*4 imonth
                integer*4 iyear
        end type
        type (iar) idate_struct
        call idate_(idate_struct)
        imonth = idate_struct % imonth
        iday   = idate_struct % iday
        iyear  = idate_struct % iyear
        
      return
      end
#endif
