!**************************************************************************
! subroutine GasRateConstants_Com: generates thermal rate coefficients
!                   for the selected mechanism
! nomenclature:
! rk_com    = reaction rate constants for common mechanism (molec-cc-s)
! te        = ambient atmospheric temperature (K)
! iregime = selected mechanism for the current chemical regime (1-6)
!
! author: Rahul A. Zaveri
! date  : february 1996
!
!-------------------------------------------------------------------------
      subroutine GasRateConstants_Com
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      real(r8) :: rk0, rnn, rki, rmm
      real(r8) :: rko, rk2, rk3
      real(r8) :: RK_OH_HNO3, RK_2HO2, RK_2HO2_H2O, RK_CO_OH
      real(r8) :: ARR, Troe

      rk_com(1) = rk_photo(jphoto_no2)
      rk_com(2) = rk_photo(jphoto_no3)
      rk_com(3) = rk_photo(jphoto_hono)
      rk_com(4) = rk_photo(jphoto_hno3)
      rk_com(5) = rk_photo(jphoto_hno4)
      rk_com(6) = rk_photo(jphoto_n2o5)
      rk_com(7) = rk_photo(jphoto_o3a)
      rk_com(8) = rk_photo(jphoto_o3b)
      rk_com(9) = rk_photo(jphoto_h2o2)
      rk_com(10) = ARR(3.2d-11, 70.d0)
      rk_com(11) = ARR(1.8d-11, 110.d0)
      rk_com(12) = 2.2e-10
      rk_com(13) = cair_mlc*6.e-34*(te/300.)**(-2.3)
      rk_com(14) = ARR(8.0d-12, -2060.d0)
      rk_com(15) = ARR(6.5d-12, -120.d0)

      rk0 = 9.0e-32
      rnn = 2.0
      rki = 2.2e-11
      rmm = 0.0
      rk_com(16) = Troe(cair_mlc,te,rk0,rnn,rki,rmm)

      rk0 = 9.0e-32
      rnn = 1.5
      rki = 3.0e-11
      rmm = 0.0
      rk_com(17) = Troe(cair_mlc,te,rk0,rnn,rki,rmm)

      rk_com(18) = ARR(2.0d-12, -1400.d0)
      rk_com(19) = ARR(1.2d-13, -2450.d0)
      rk_com(20) = ARR(1.6d-12, -940.d0)
      rk_com(21) = ARR(1.1d-14, -500.d0)
      rk_com(22) = ARR(5.5d-12, -2000.d0)

      rk0 = 7.0e-31
      rnn = 2.6
      rki = 3.6e-11
      rmm = 0.1
      rk_com(23) = Troe(cair_mlc,te,rk0,rnn,rki,rmm)

      rk0 = 2.5e-30
      rnn = 4.4
      rki = 1.6e-11
      rmm = 1.7
      rk_com(24) = Troe(cair_mlc,te,rk0,rnn,rki,rmm)
      rk_com(25) = 2.2e-11
      rk_com(26) = ARR(1.8d-11, -390.d0)

             rko = 7.2e-15 * exp(785./te)
             rk2 = 4.1e-16 * exp(1440./te)
             rk3 = 1.9e-33 * exp(725./te)*cair_mlc
      RK_OH_HNO3 = rko + rk3/(1.+rk3/rk2)
      rk_com(27) = RK_OH_HNO3
      rk_com(28) = ARR(1.3d-12, 380.d0)
      rk_com(29) = ARR(4.8d-11, 250.d0)
      rk_com(30) = ARR(2.9d-12, -160.d0)

      RK_2HO2    = 2.3e-13 * exp(600./te) + 	   &  ! ho2 + ho2 --> h2o2
                   1.7e-33 * exp(1000./te)*cair_mlc
      rk_com(31) = RK_2HO2

      RK_2HO2_H2O= RK_2HO2*1.4e-21*exp(2200./te)! ho2 + ho2 + h2o --> h2o2
      rk_com(32) = RK_2HO2_H2O

      rk_com(33) = ARR(3.5d-12, 250.d0)

      rk0 = 1.8e-31
      rnn = 3.2
      rki = 4.7e-12
      rmm = 1.4
      rk_com(34) = Troe(cair_mlc,te,rk0,rnn,rki,rmm)

      rk_com(35) = 5.0e-16
      rk_com(36) = rk_com(34)*ARR(4.8d26, -10900.d0)
      rk_com(37) = ARR(1.5d-11, 170.d0)
      rk_com(38) = ARR(4.5d-14, -1260.d0)

      rk0 = 2.2e-30
      rnn = 3.9
      rki = 1.5e-12
      rmm = 0.7
      rk_com(39) = Troe(cair_mlc,te,rk0,rnn,rki,rmm)

      rk_com(40) = ARR(8.5d-13, -2450.d0)
      rk_com(41) = 3.5e-12
      rk_com(42) = 0.0  ! 2.0e-21   ! N2O5 + H2O --> 2HNO3 (homogeneous pathway)
      rk_com(43) = rk_com(39)*ARR(3.7d26, -11000.d0)

      RK_CO_OH   = 1.5e-13 * (1.+8.18e-23*te*cair_mlc) ! co + oh --> ho2
      rk_com(44) = RK_CO_OH

      rk0 = 3.0e-31
      rnn = 3.3
      rki = 1.5e-12
      rmm = 0.0
      rk_com(45) = Troe(cair_mlc,te,rk0,rnn,rki,rmm)

      rk_com(46) = te**.667*ARR(2.8d-14, -1575.d0)
      rk_com(47) = te**2*ARR(1.5d-17, -492.d0)
      rk_com(48) = ARR(6.7d-12, -600.d0)
      rk_com(49) = rk_photo(jphoto_hchoa)
      rk_com(50) = rk_photo(jphoto_hchob)
      rk_com(51) = 1.0e-11
      rk_com(52) = ARR(3.4d-13, -1900.d0)
      rk_com(53) = rk_photo(jphoto_ch3ooh)
      rk_com(54) = rk_photo(jphoto_ethooh)
      rk_com(55) = ARR(3.8d-12, 200.d0)
      rk_com(56) = ARR(3.8d-12, 200.d0)
      rk_com(57) = ARR(3.0d-12, 280.d0)
      rk_com(58) = ARR(2.6d-12, 365.d0)
      rk_com(59) = 1.1e-12
      rk_com(60) = 2.5e-12
      rk_com(61) = ARR(3.8d-13, 800.d0)
      rk_com(62) = ARR(7.5d-13, 700.d0)
      rk_com(63) = rk_param(jch3o2)
      rk_com(64) = rk_param(jethp)
      rk_com(65) = ARR(7.0d-12, -235.d0)
      rk_com(66) = rk_photo(jphoto_ald2)
      rk_com(67) = ARR(5.6d-12, 270.d0)
      rk_com(68) = ARR(1.4d-12, -1900.d0)

      rk0 = 9.7e-29
      rnn = 5.6
      rki = 9.3e-12
      rmm = 1.5
      rk_com(69) = Troe(cair_mlc,te,rk0,rnn,rki,rmm)

      rk_com(70) = rk_com(69)*ARR(1.1d28, -14000.d0)
      rk_com(71) = ARR(5.3d-12, 360.d0)
      rk_com(72) = 4.0e-12
      rk_com(73) = ARR(4.5d-13, 1000.d0)
      rk_com(74) = rk_param(jc2o3)
!
! Heterogeneous chemistry
!      rk_com(64) = rk_het(1)	! O3 -->
!      rk_com(65) = rk_het(2)	! HO2 --> 0.5H2O2
!      rk_com(66) = rk_het(3)	! NO2 --> 0.5HONO + 0.5HNO3
!      rk_com(67) = rk_het(4)	! N2O5 --> 2HNO3
!      rk_com(68) = rk_het(5)	! HNO3 --> NO2
!      rk_com(69) = rk_het(6)	! HNO3 --> NO
!      rk_com(70) = rk_het(7)	! NO3 --> NO + O2
!
      return
      end subroutine GasRateConstants_Com
