      SUBROUTINE MACH
C
C     HERE IS DETERMINED THE TYPE OF COMPUTER WE ARE WORKING ON:
C
C     DEFINITION OF MACHINE TYPES:
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
C     11 .... Apollo Domain/OS "UNIX"
C     12 .... IBM AIX XL Fortran
C
C     Use machine names in all programs:
      COMMON /MACH1/ MACTYP
      COMMON /MACH2/ MACHIN(0:20)
      CHARACTER*7  MACHIN
C
C-----------------------------------------------------------------------
C
      MACTYP = 10
      MACHIN(0) = 'Unknown'
      MACHIN(1) = 'IBM'
      MACHIN(2) = 'CRAY-1'
      MACHIN(3) = 'VAX'
      MACHIN(4) = 'CDC'
      MACHIN(5) = 'UNIVAC'
      MACHIN(6) = 'FPS+VAX'
      MACHIN(7) = 'CRAY-2'
      MACHIN(8) = 'Fujitsu'
      MACHIN(9) = 'MAP+VAX'
      MACHIN(10) = 'SunUNIX'
      MACHIN(11) = 'Apollo'
      RETURN
      END
