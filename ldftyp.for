      SUBROUTINE LDFTYP (LDF)
C
C     Associate LDF XC-functional names
C
      CHARACTER*20 LDF(0:3)
C
      LDF(0) = 'Unknown'
      LDF(1) = 'Slater X-alpha' 
      LDF(2) = 'Wigner interpolation'
      LDF(3) = 'Ceperley-Alder'
      RETURN
      END
