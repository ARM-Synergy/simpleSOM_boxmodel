      subroutine GasRates_Com(s)
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      real(r8) :: s(ntot_max)

      r_com(1) = rk_com(1)*s(ino2)
      r_com(2) = rk_com(2)*s(ino3)
      r_com(3) = rk_com(3)*s(ihono)
      r_com(4) = rk_com(4)*s(ihno3)
      r_com(5) = rk_com(5)*s(ihno4)
      r_com(6) = rk_com(6)*s(in2o5)
      r_com(7) = rk_com(7)*s(io3)
      r_com(8) = rk_com(8)*s(io3)
      r_com(9) = rk_com(9)*s(ih2o2)
      r_com(10) = rk_com(10)*s(io1d)*.21*cair_mlc
      r_com(11) = rk_com(11)*s(io1d)*.79*cair_mlc
      r_com(12) = rk_com(12)*s(io1d)*H2O
      r_com(13) = rk_com(13)*s(io3p)*.21*cair_mlc
      r_com(14) = rk_com(14)*s(io3p)*s(io3)
      r_com(15) = rk_com(15)*s(io3p)*s(ino2)
      r_com(16) = rk_com(16)*s(io3p)*s(ino2)
      r_com(17) = rk_com(17)*s(io3p)*s(ino)
      r_com(18) = rk_com(18)*s(io3)*s(ino)
      r_com(19) = rk_com(19)*s(io3)*s(ino2)
      r_com(20) = rk_com(20)*s(io3)*s(ioh)
      r_com(21) = rk_com(21)*s(io3)*s(iho2)
      r_com(22) = rk_com(22)*s(ioh)*H2
      r_com(23) = rk_com(23)*s(ioh)*s(ino)
      r_com(24) = rk_com(24)*s(ioh)*s(ino2)
      r_com(25) = rk_com(25)*s(ioh)*s(ino3)
      r_com(26) = rk_com(26)*s(ioh)*s(ihono)
      r_com(27) = rk_com(27)*s(ioh)*s(ihno3)
      r_com(28) = rk_com(28)*s(ioh)*s(ihno4)
      r_com(29) = rk_com(29)*s(ioh)*s(iho2)
      r_com(30) = rk_com(30)*s(ioh)*s(ih2o2)
      r_com(31) = rk_com(31)*s(iho2)*s(iho2)
      r_com(32) = rk_com(32)*s(iho2)*s(iho2)*H2O
      r_com(33) = rk_com(33)*s(iho2)*s(ino)
      r_com(34) = rk_com(34)*s(iho2)*s(ino2)
      r_com(35) = rk_com(35)*s(iho2)*s(ino2)
      r_com(36) = rk_com(36)*s(ihno4)
      r_com(37) = rk_com(37)*s(ino3)*s(ino)
      r_com(38) = rk_com(38)*s(ino3)*s(ino2)
      r_com(39) = rk_com(39)*s(ino3)*s(ino2)
      r_com(40) = rk_com(40)*s(ino3)*s(ino3)
      r_com(41) = rk_com(41)*s(ino3)*s(iho2)
      r_com(42) = rk_com(42)*s(in2o5)*H2O
      r_com(43) = rk_com(43)*s(in2o5)
      r_com(44) = rk_com(44)*s(ico)*s(ioh)
      r_com(45) = rk_com(45)*s(iso2)*s(ioh)
      r_com(46) = rk_com(46)*s(ich4)*s(ioh)
      r_com(47) = rk_com(47)*s(ic2h6)*s(ioh)
      r_com(48) = rk_com(48)*s(ich3oh)*s(ioh)
      r_com(49) = rk_com(49)*s(ihcho)
      r_com(50) = rk_com(50)*s(ihcho)
      r_com(51) = rk_com(51)*s(ihcho)*s(ioh)
      r_com(52) = rk_com(52)*s(ihcho)*s(ino3)
      r_com(53) = rk_com(53)*s(ich3ooh)
      r_com(54) = rk_com(54)*s(iethooh)
      r_com(55) = rk_com(55)*s(ich3ooh)*s(ioh)
      r_com(56) = rk_com(56)*s(iethooh)*s(ioh)
      r_com(57) = rk_com(57)*s(ich3o2)*s(ino)
      r_com(58) = rk_com(58)*s(iethp)*s(ino)
      r_com(59) = rk_com(59)*s(ich3o2)*s(ino3)
      r_com(60) = rk_com(60)*s(iethp)*s(ino3)
      r_com(61) = rk_com(61)*s(ich3o2)*s(iho2)
      r_com(62) = rk_com(62)*s(iethp)*s(iho2)
      r_com(63) = rk_com(63)*s(ich3o2)
      r_com(64) = rk_com(64)*s(iethp)
      r_com(65) = rk_com(65)*s(ianol)*s(ioh)
      r_com(66) = rk_com(66)*s(iald2)
      r_com(67) = rk_com(67)*s(iald2)*s(ioh)
      r_com(68) = rk_com(68)*s(iald2)*s(ino3)
      r_com(69) = rk_com(69)*s(ic2o3)*s(ino2)
      r_com(70) = rk_com(70)*s(ipan)
      r_com(71) = rk_com(71)*s(ic2o3)*s(ino)
      r_com(72) = rk_com(72)*s(ic2o3)*s(ino3)
      r_com(73) = rk_com(73)*s(ic2o3)*s(iho2)
      r_com(74) = rk_com(74)*s(ic2o3)
!
! heterogeneous chemistry
!      r_com(71) = rk_com(71)*s(io3)	! O3 -->
!      r_com(72) = rk_com(72)*s(iho2)	! HO2 --> 0.5H2O2
!      r_com(73) = rk_com(73)*s(ino2)	! NO2 --> 0.5HONO + 0.5HNO3
!      r_com(74) = rk_com(74)*s(in2o5)	! N2O5 --> 2HNO3
!      r_com(75) = rk_com(75)*s(ihno3)	! HNO3 --> NO2
!      r_com(76) = rk_com(76)*s(ihno3)	! HNO3 --> NO
!      r_com(77) = rk_com(77)*s(ino3)	! NO3 --> NO

      return
      end subroutine GasRates_Com
