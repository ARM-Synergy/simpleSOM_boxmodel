
      FUNCTION ARR(AA,BB)
      
      USE mod_MAIN, ONLY: R8,te
      
      IMPLICIT NONE
      
      REAL(r8) :: ARR
      REAL(r8) :: AA,BB
      
      ARR = AA*EXP(BB/te)
      
      RETURN
      END FUNCTION ARR
