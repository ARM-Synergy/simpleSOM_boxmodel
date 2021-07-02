      subroutine GasRates_Bio(s)
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      real(r8) :: s(ntot_max)

!FLAG1
      r_bio(1) = rk_bio(17)*s(iOH)*s(iCn3)
      r_bio(2) = rk_bio(18)*s(iOH)*s(iCn2)
      r_bio(3) = rk_bio(19)*s(iOH)*s(iCn1)
      r_bio(4) = rk_bio(20)*s(iOH)*s(iC0)
      r_bio(5) = rk_bio(21)*s(iOH)*s(iC1)
      r_bio(6) = rk_bio(22)*s(iOH)*s(iC2)
      r_bio(7) = rk_bio(23)*s(iOH)*s(iC3)
      r_bio(8) = rk_bio(24)*s(iOH)*s(iC4)
      r_bio(9) = rk_bio(25)*s(iOH)*s(iC5)
      r_bio(10) = rk_bio(26)*s(iOH)*s(iC6)
      r_bio(11) = rk_bio(27)*s(iOH)*s(iC7)
      r_bio(12) = rk_bio(28)*s(iOH)*s(iC8)
      r_bio(13) = rk_bio(29)*s(iOH)*s(iC9)
      r_bio(14) = rk_bio(30)*s(iOH)*s(iisop)
!FLAG2
      
      return
      end subroutine GasRates_Bio



