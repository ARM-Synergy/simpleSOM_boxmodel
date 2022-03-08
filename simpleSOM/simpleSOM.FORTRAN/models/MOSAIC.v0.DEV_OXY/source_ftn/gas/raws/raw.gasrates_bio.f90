

      SUBROUTINE GASRATES_BIO(S)

      USE mod_MAIN
      USE mod_GAS

      IMPLICIT NONE

      REAL(R8) :: S(nTOT_MAX)

!FLAG1

      
      r_BIO(14) = rk_BIO(30)*S(iOH)*S(iISOP)
      
      RETURN
      END SUBROUTINE GASRATES_BIO



