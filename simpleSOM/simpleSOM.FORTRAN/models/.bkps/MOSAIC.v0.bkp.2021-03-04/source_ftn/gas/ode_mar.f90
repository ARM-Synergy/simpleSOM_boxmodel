!***********************************************************************
      subroutine ode_mar
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      real(r8) :: a, b
!
!
      a = 5.e+5/(5.e+5 + o2*3.e-12)
      b = 1.5e+7/(1.5e+7 + o2*1.2e-12)
!
!
      p_mar(ino)= r_mar(17)
      d_mar(ino)= r_mar(5)+r_mar(9)+r_mar(24)+r_mar(28)
!
      p_mar(ino2)= r_mar(5)+r_mar(9)+r_mar(24)
      d_mar(ino2)= r_mar(17)+r_mar(27)
!
      p_mar(ino3)=   0.0
      d_mar(ino3)= r_mar(2)+r_mar(12)
!
      p_mar(ihono)= r_mar(28)
      d_mar(ihono)=   0.0
!
      p_mar(ihno3)= r_mar(2)+r_mar(12)+r_mar(27)
      d_mar(ihno3)=   0.0
!
      p_mar(io3)=   0.0
      d_mar(io3)= r_mar(18)
!
      p_mar(io3p)=   0.0
      d_mar(io3p)= r_mar(3)
!
      p_mar(ioh)= r_mar(19)
      d_mar(ioh)= r_mar(1)+r_mar(4)+r_mar(7)+r_mar(8)+   &
                  r_mar(14)+r_mar(21)
!
      p_mar(iho2)= (1.-a)*r_mar(4)+r_mar(6)+(1.-b)*r_mar(7)+   &
                   r_mar(10)+r_mar(20)+r_mar(25)+r_mar(30)
      d_mar(iho2)= r_mar(11)+r_mar(19)+r_mar(29)
!
      p_mar(ih2o2)= r_mar(11)
      d_mar(ih2o2)=   0.0
!
      p_mar(iso2)= r_mar(16)
      d_mar(iso2)=   0.0 + r_mar(33)
!
      p_mar(ih2so4)= r_mar(26)
      d_mar(ih2so4)=   0.0
!
      p_mar(ich3o2)= r_mar(3)+a*r_mar(4)+b*r_mar(7)+r_mar(16)+   &
                     r_mar(26)
      d_mar(ich3o2)= r_mar(6)+r_mar(10)+r_mar(13)+r_mar(20)+   &
                     r_mar(25)
!
      p_mar(ich3ooh)= r_mar(13)
      d_mar(ich3ooh)=   0.0
!
      p_mar(ihcho)= r_mar(5)+r_mar(6)+r_mar(6)+r_mar(9)+   &
                    r_mar(10)+r_mar(10)+r_mar(20)+r_mar(25)
      d_mar(ihcho)= r_mar(30)
!
      p_mar(idms)= 0.0
      d_mar(idms)= r_mar(1)+r_mar(2)+r_mar(3)+r_mar(4)
!
      p_mar(imsa)= r_mar(15)+r_mar(21)+r_mar(27)+r_mar(28)+   &
                   r_mar(29)+r_mar(30)
      d_mar(imsa)=   0.0
!
      p_mar(idmso)= (1.-a)*r_mar(4)
      d_mar(idmso)= r_mar(7) + r_mar(34)
!
      p_mar(idmso2)= (1.-b)*r_mar(7)
      d_mar(idmso2)= r_mar(8) + r_mar(35)
!
      p_mar(ich3so2h)= b*r_mar(7)
      d_mar(ich3so2h)= r_mar(11)+r_mar(12)+r_mar(13)+r_mar(14)+   &
                       r_mar(15)
!
      p_mar(ich3sch2oo)= r_mar(1)+r_mar(2)
      d_mar(ich3sch2oo)= r_mar(5)+r_mar(6)+r_mar(31)+2.*r_mar(32)
!
      p_mar(ich3so2)= r_mar(3)+a*r_mar(4)+r_mar(5)+r_mar(6)+   &
                      r_mar(9)+r_mar(10)+r_mar(11)+r_mar(12)+   &
                      r_mar(13)+r_mar(14)+r_mar(15)+r_mar(23)+   &
                      r_mar(31)+1.85*r_mar(32)
      d_mar(ich3so2)= r_mar(16)+r_mar(17)+r_mar(18)+r_mar(19)+   &
                      r_mar(20)+r_mar(21)+r_mar(22)+r_mar(31)
!
      p_mar(ich3so3)= r_mar(17)+r_mar(18)+r_mar(19)+r_mar(20)+   &
                      r_mar(24)+r_mar(25)+r_mar(31)
      d_mar(ich3so3)= r_mar(15)+r_mar(26)+r_mar(27)+r_mar(28)+   &
                      r_mar(29)+r_mar(30)
!
      p_mar(ich3so2oo)= r_mar(22)
      d_mar(ich3so2oo)= r_mar(23)+r_mar(24)+r_mar(25)
!
      p_mar(ich3so2ch2oo)= r_mar(8)
      d_mar(ich3so2ch2oo)= r_mar(9)+r_mar(10)
!
      p_mar(isulfhox)= 0.15*r_mar(32)
      d_mar(isulfhox)= 0.0
!

!
      return
      end subroutine ode_mar
!***********************************************************************
