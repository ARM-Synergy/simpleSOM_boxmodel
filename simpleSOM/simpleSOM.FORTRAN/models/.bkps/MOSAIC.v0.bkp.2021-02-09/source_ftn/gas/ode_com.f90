!***********************************************************************
      subroutine ode_com
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      p_com(ih2so4)= r_com(45)
      d_com(ih2so4)= 0.0

      p_com(ihno3)= r_com(24)+.3*r_com(41)+r_com(42)+r_com(42)   &
               +r_com(52)+r_com(68)
      d_com(ihno3)= r_com(4)+r_com(27)

      p_com(ihcl)= 0.0
      d_com(ihcl)= 0.0

      p_com(inh3)= 0.0
      d_com(inh3)= 0.0

      p_com(ino)= r_com(1)+0.11*r_com(2)+r_com(3)+r_com(15)+r_com(38)
      d_com(ino)= r_com(17)+r_com(18)+r_com(23)+r_com(33)+r_com(37)   &
             +r_com(57)+r_com(58)+r_com(71)

      p_com(ino2)= 0.89*r_com(2)+r_com(4)+r_com(5)   &
              +r_com(6)+r_com(17)+r_com(18)    &
              +r_com(25)+r_com(26)+r_com(28)   &
              +r_com(33)+r_com(36)+r_com(37)   &
              +r_com(37)+r_com(38)+r_com(40)+r_com(40)+.7*r_com(41)   &
              +r_com(43)+r_com(57)+r_com(58)+r_com(59)   &
              +r_com(60)+r_com(70)   &
              +r_com(71)+r_com(72)
      d_com(ino2)= r_com(1)+r_com(15)+r_com(16)+r_com(19)   &
              +r_com(24)+r_com(34)   &
              +r_com(35)+r_com(38)+r_com(39)+r_com(69)

      p_com(ino3)= r_com(6)+r_com(16)+r_com(19)+r_com(27)+r_com(43)
      d_com(ino3)= r_com(2)+r_com(25)+r_com(37)+r_com(38)   &
              +r_com(39)+r_com(40)   &
              +r_com(40)+r_com(41)+r_com(52)+r_com(59)   &
              +r_com(60)+r_com(68)   &
              +r_com(72)

      p_com(in2o5)= r_com(39)
      d_com(in2o5)= r_com(6)+r_com(42)+r_com(43)

      p_com(ihono)= r_com(23)+r_com(35)
      d_com(ihono)= r_com(3)+r_com(26)

      p_com(ihno4)= r_com(34)
      d_com(ihno4)= r_com(5)+r_com(28)+r_com(36)

      p_com(io3)= r_com(13)+.4*r_com(73)
      d_com(io3)= r_com(7)+r_com(8)+r_com(14)+r_com(18)   &
             +r_com(19)+r_com(20)+r_com(21)

      p_com(io1d)= r_com(8)
      d_com(io1d)= r_com(10)+r_com(11)+r_com(12)

      p_com(io3p)= r_com(1)+0.89*r_com(2)+r_com(7)   &
             +r_com(10)+r_com(11)
      d_com(io3p)= r_com(13)+r_com(14)+r_com(15)+r_com(16)+r_com(17)

      ! edit wkc, uncomment! remove!
!      p_com(ioh)= r_com(3)+r_com(4)+2*r_com(9)   &
!             +2*r_com(12)+r_com(21)+r_com(33)    &
!             +.7*r_com(41)+r_com(53)+r_com(54)   &
!             +.3*r_com(55)+.5*r_com(56)
!      d_com(ioh)= r_com(20)+r_com(22)+r_com(23)  &
!             +r_com(24)+r_com(25)+r_com(26)      &
!             +r_com(27)+r_com(28)+r_com(29)+r_com(30)   &
!             +r_com(44)+r_com(45)   &
!             +r_com(46)+r_com(47)+r_com(48)   &
!             +r_com(51)+r_com(55)+r_com(56)   &
!             +r_com(65)+r_com(67)

      p_com(iho2)= r_com(5)+r_com(20)+r_com(22)   &
              +r_com(25)+r_com(30)+r_com(36)   &
              +r_com(44)+r_com(45)+r_com(48)   &
              +2*r_com(49)+r_com(51)           &
              +r_com(52)+r_com(53)+r_com(54)   &
              +r_com(57)+r_com(58)+r_com(59)   &
              +r_com(60)+.32*r_com(63)+.6*r_com(64)   &
              +r_com(65)+r_com(66)
      d_com(iho2)= r_com(21)+r_com(29)+r_com(31)   &
              +r_com(31)+r_com(32)+r_com(32)   &
              +r_com(33)+r_com(34)+r_com(35)   &
              +r_com(41)+r_com(61)+r_com(62)   &
              +r_com(73)

      p_com(ih2o2)= r_com(31)+r_com(32)
      d_com(ih2o2)= r_com(9)+r_com(30)

      p_com(ico)= r_com(49)+r_com(50)+r_com(51)   &
              +r_com(52)+r_com(66)
      d_com(ico)= r_com(44)

      p_com(iso2)= 0.0
      d_com(iso2)= r_com(45)

      p_com(ich4)= 0.0
      d_com(ich4)= r_com(46)

      p_com(ic2h6)= .2*r_com(64)
      d_com(ic2h6)= r_com(47)

      p_com(ich3o2)= r_com(46)+.7*r_com(55)+r_com(66)   &
                +r_com(71)+r_com(72)   &
                +r_com(74)
      d_com(ich3o2)= r_com(57)+r_com(59)+r_com(61)+r_com(63)

      p_com(iethp)= r_com(47)+.5*r_com(56)
      d_com(iethp)= r_com(58)+r_com(60)+r_com(62)+r_com(64)

      p_com(ihcho)= r_com(48)+r_com(53)+.3*r_com(55)   &
               +r_com(57)+r_com(59)   &
               +.66*r_com(63)
      d_com(ihcho)= r_com(49)+r_com(50)+r_com(51)+r_com(52)

      p_com(ich3oh)= .34*r_com(63)
      d_com(ich3oh)= r_com(48)

      p_com(ianol)= 0.0
      d_com(ianol)= r_com(65)

      p_com(ich3ooh)= r_com(61)
      d_com(ich3ooh)= r_com(53)+r_com(55)

      p_com(iethooh)= r_com(62)
      d_com(iethooh)= r_com(54)+r_com(56)

      p_com(iald2)= r_com(54)+.5*r_com(56)+r_com(58)   &
               +r_com(60)+.8*r_com(64)   &
               +r_com(65)
      d_com(iald2)= r_com(66)+r_com(67)+r_com(68)

      p_com(ihcooh)= 0.0
      d_com(ihcooh)= 0.0

      p_com(ircooh)= .4*r_com(73)
      d_com(ircooh)= 0.0

      p_com(ic2o3)= r_com(67)+r_com(68)+r_com(70)
      d_com(ic2o3)= r_com(69)+r_com(71)+r_com(72)   &
              +r_com(73)+r_com(74)

      p_com(ipan)= r_com(69)
      d_com(ipan)= r_com(70)

      return
      end subroutine ode_com
