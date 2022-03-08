

      SUBROUTINE GASRATES_BIO(S)

      USE mod_MAIN
      USE mod_GAS

      IMPLICIT NONE

      REAL(R8) :: S(nTOT_MAX)

!FLAG1
      r_BIO(1) = rk_BIO(17)*S(iOH)*S(iCN3)
      r_BIO(2) = rk_BIO(18)*S(iOH)*S(iCN2)
      r_BIO(3) = rk_BIO(19)*S(iOH)*S(iCN1)
      r_BIO(4) = rk_BIO(20)*S(iOH)*S(iC0)
      r_BIO(5) = rk_BIO(21)*S(iOH)*S(iC1)
      r_BIO(6) = rk_BIO(22)*S(iOH)*S(iC2)
      r_BIO(7) = rk_BIO(23)*S(iOH)*S(iC3)
      r_BIO(8) = rk_BIO(24)*S(iOH)*S(iC4)
      r_BIO(9) = rk_BIO(25)*S(iOH)*S(iC5)
      r_BIO(10) = rk_BIO(26)*S(iOH)*S(iC6)
      r_BIO(11) = rk_BIO(27)*S(iOH)*S(iC7)
      r_BIO(12) = rk_BIO(28)*S(iOH)*S(iC8)
      r_BIO(13) = rk_BIO(29)*S(iOH)*S(iC9)

      
      r_BIO(14) = rk_BIO(30)*S(iOH)*S(iISOP)
      
      RETURN
      END SUBROUTINE GASRATES_BIO



