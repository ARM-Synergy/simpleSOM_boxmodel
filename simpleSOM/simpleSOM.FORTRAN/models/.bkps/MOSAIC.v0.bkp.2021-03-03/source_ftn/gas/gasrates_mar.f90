      subroutine GasRates_Mar(s)
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      real(r8) :: s(ntot_max)

      r_mar(1)    = rk_mar(1)   *s(idms)*s(ioh)
      r_mar(2)    = rk_mar(2)   *s(idms)*s(ino3)
      r_mar(3)    = rk_mar(3)   *s(idms)*s(io3p)
      r_mar(4)    = rk_mar(4)   *s(idms)*s(ioh)
      r_mar(5)    = rk_mar(5)   *s(ich3sch2oo)*s(ino)
      r_mar(6)    = rk_mar(6)   *s(ich3sch2oo)*s(ich3o2)
      r_mar(7)    = rk_mar(7)   *s(idmso)*s(ioh)
      r_mar(8)    = rk_mar(8)   *s(idmso2)*s(ioh)
      r_mar(9)    = rk_mar(9)   *s(ich3so2ch2oo)*s(ino)
      r_mar(10)   = rk_mar(10)  *s(ich3so2ch2oo)*s(ich3o2)
      r_mar(11)   = rk_mar(11)  *s(ich3so2h)*s(iho2)
      r_mar(12)   = rk_mar(12)  *s(ich3so2h)*s(ino3)
      r_mar(13)   = rk_mar(13)  *s(ich3so2h)*s(ich3o2)
      r_mar(14)   = rk_mar(14)  *s(ich3so2h)*s(ioh)
      r_mar(15)   = rk_mar(15)  *s(ich3so2h)*s(ich3so3)
      r_mar(16)   = rk_mar(16)  *s(ich3so2)
      r_mar(17)   = rk_mar(17)  *s(ich3so2)*s(ino2)
      r_mar(18)   = rk_mar(18)  *s(ich3so2)*s(io3)
      r_mar(19)   = rk_mar(19)  *s(ich3so2)*s(iho2)
      r_mar(20)   = rk_mar(20)  *s(ich3so2)*s(ich3o2)
      r_mar(21)   = rk_mar(21)  *s(ich3so2)*s(ioh)
      r_mar(22)   = rk_mar(22)  *s(ich3so2)*o2
      r_mar(23)   = rk_mar(23)  *s(ich3so2oo)
      r_mar(24)   = rk_mar(24)  *s(ich3so2oo)*s(ino)
      r_mar(25)   = rk_mar(25)  *s(ich3so2oo)*s(ich3o2)
      r_mar(26)   = rk_mar(26)  *s(ich3so3)
      r_mar(27)   = rk_mar(27)  *s(ich3so3)*s(ino2)
      r_mar(28)   = rk_mar(28)  *s(ich3so3)*s(ino)
      r_mar(29)   = rk_mar(29)  *s(ich3so3)*s(iho2)
      r_mar(30)   = rk_mar(30)  *s(ich3so3)*s(ihcho)
!
      r_mar(31)   = rk_mar(31)	*s(ich3sch2oo)*s(ich3so2)
      r_mar(32)   = rk_mar(32)  *s(ich3sch2oo)*s(ich3sch2oo)
!
      r_mar(33)   = rk_mar(33)  *s(iso2)
      r_mar(34)   = rk_mar(34)  *s(idmso)
      r_mar(35)   = rk_mar(35)  *s(idmso2)
!
      return
      end subroutine GasRates_Mar

