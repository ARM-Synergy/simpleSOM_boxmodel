!**************************************************************************
! subroutine GasRateConstants_Urb: generates thermal rate coefficients
!                   for the selected mechanism
! nomenclature:
! rk_urb    = reaction rate constants for hc1 mechanism    (molec-cc-s)
! te        = ambient atmospheric temperature (K)
!
! author: Rahul A. Zaveri
! date  : february 1996
!
!-------------------------------------------------------------------------
      subroutine GasRateConstants_Urb
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      real(r8) :: rk0, rnn, rki, rmm
      real(r8) :: ARR, Troe

      rk_urb(1) = 8.1e-13
      rk_urb(2) = rk_photo(jphoto_aone)
      rk_urb(3) = te**2*ARR(5.3d-18, -230.d0)
      rk_urb(4) = rk_photo(jphoto_mgly)
      rk_urb(5) = 1.7e-11
      rk_urb(6) = ARR(1.4d-12, -1900.d0)
      rk_urb(7) = ARR(1.2d-14, -2630.d0)

      rk0 = 1.0e-28
      rnn = 0.8
      rki = 8.8e-12
      rmm = 0.0
      rk_urb(8) = Troe(cair_mlc,te,rk0,rnn,rki,rmm)

      rk_urb(9) = ARR(4.2d-15, -1800.d0)
      rk_urb(10) = ARR(8.9d-16, -392.d0)
      rk_urb(11) = ARR(5.8d-12, 478.d0)
      rk_urb(12) = ARR(2.9d-11, 255.d0)
      rk_urb(13) = ARR(3.1d-13, -1010.d0)
      rk_urb(14) = 2.5e-12
      rk_urb(15) = ARR(2.1d-12, 322.d0)
      rk_urb(16) = ARR(1.7d-11, 116.d0)
      rk_urb(17) = 8.1e-12
      rk_urb(18) = 4.1e-11
      rk_urb(19) = 2.2e-11
      rk_urb(20) = 1.4e-11
      rk_urb(21) = 3.0e-11
      rk_urb(22) = rk_photo(jphoto_open)
      rk_urb(23) = ARR(5.4d-17, -500.d0)
      rk_urb(24) = rk_photo(jphoto_rooh)
      rk_urb(25) = ARR(3.8d-12, 200.d0)
      rk_urb(26) = ARR(1.6d-11, -540.d0)
      rk_urb(27) = rk_photo(jphoto_onit)
      rk_urb(28) = 4.0e-12
      rk_urb(29) = 4.0e-12
      rk_urb(30) = 4.0e-12
      rk_urb(31) = 4.0e-12
      rk_urb(32) = 2.5e-12
      rk_urb(33) = 1.2e-12
      rk_urb(34) = 4.0e-12
      rk_urb(35) = 2.5e-12
      rk_urb(36) = ARR(1.7d-13, 1300.d0)
      rk_urb(37) = ARR(1.2d-13, 1300.d0)
      rk_urb(38) = ARR(1.7d-13, 1300.d0)
      rk_urb(39) = ARR(1.7d-13, 1300.d0)
      rk_urb(40) = rk_param(jro2)
      rk_urb(41) = rk_param(jano2)
      rk_urb(42) = rk_param(jnap)
      rk_urb(43) = rk_param(jxo2)
      rk_urb(44) = 1.0e-11

      return
      end subroutine GasRateConstants_Urb
