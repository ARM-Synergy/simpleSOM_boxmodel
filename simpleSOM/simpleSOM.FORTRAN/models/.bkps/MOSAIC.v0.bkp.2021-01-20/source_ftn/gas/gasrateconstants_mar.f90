!**************************************************************************
! subroutine GasRateConstants_Mar: generates thermal rate coefficients
!                   for the selected mechanism
! nomenclature:
! rk_mar    = reaction rate constants for marine mechanism (molec-cc-s)
! te        = ambient atmospheric temperature (K)
!
! author: Rahul A. Zaveri
! date  : february 1996
!
!-------------------------------------------------------------------------
      subroutine GasRateConstants_Mar
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      real(r8) :: rk_tot, rk_tot_num, rk_tot_den, Babs, Badd

! abstraction reaction
! Hynes et al. (1986)
      rk_tot_num  =       te * exp(-234./te) +   &
                    8.46e-10 * exp(7230./te) +   &
                    2.68e-10 * exp(7810./te)
      rk_tot_den  = 1.04e+11 * te + 88.1 * exp(7460./te)
      rk_tot	  = rk_tot_num/rk_tot_den
!
      rk_mar(1)   = 9.60e-12 * exp(-234./te) ! ch3sch3 + oh --> ch3sch2
      Babs        = rk_mar(1)/rk_tot
      Badd	  = 1. - Babs
      rk_mar(2)   = 1.40e-13 * exp(500./te)  ! ch3sch3 + no3 -->
      rk_mar(3)   = 1.26e-11 * exp(409./te)  ! ch3sch3 + o3p -->
!
! addition reaction
      rk_mar(4)   = Badd*rk_tot		     ! ch3sch3 + oh --> ch3s(oh)ch3
      rk_mar(5)   = 8.0e-12
      rk_mar(6)   = 1.8e-13
      rk_mar(7)   = 5.8e-11
      rk_mar(8)   = 1.0e-14
      rk_mar(9)   = 5.0e-12
      rk_mar(10)  = 1.8e-13
      rk_mar(11)  = 1.0e-15
      rk_mar(12)  = 1.0e-13
      rk_mar(13)  = 1.0e-15
      rk_mar(14)  = 1.6e-11
      rk_mar(15)  = 1.0e-13
      rk_mar(16)  = 2.5e+13 * exp(-8686./te) ! ch3so2 --> so2 + ch3o2
      rk_mar(17)  = 1.0e-14
      rk_mar(18)  = 5.0e-15
      rk_mar(19)  = 2.5e-13
      rk_mar(20)  = 2.5e-13
      rk_mar(21)  = 5.0e-11
      rk_mar(22)  = 2.6e-18
      rk_mar(23)  = 3.3
      rk_mar(24)  = 1.0e-11
      rk_mar(25)  = 5.5e-12
      rk_mar(26)  = 2.0e+17 * exp(-12626./te) ! ch3so3 --> h2so4 + ch3o2
      rk_mar(27)  = 3.0e-15
      rk_mar(28)  = 3.0e-15
      rk_mar(29)  = 5.0e-11
      rk_mar(30)  = 1.6e-15
!
      rk_mar(31)  = 2.5e-13	! ch3sch2oo + ch3so2 --> ch3so3 + ch3so2
      rk_mar(32)  = 8.6e-14	! 2ch3sch2oo --> .15mtf + 1.85ch3so2
!
! dry deposition
      rk_mar(33)  = 0.0 ! 2.0e-5	! 1/s
      rk_mar(34)  = 0.0 ! 2.0e-5
      rk_mar(35)  = 0.0 ! 2.0e-5

      return
      end subroutine GasRateConstants_Mar
