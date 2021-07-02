!***********************************************************************
      subroutine ode_urb
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      real(r8) :: yield_factor

!
      p_urb(ihno3)= r_urb(6)+r_urb(19)
      d_urb(ihno3)= 0.0
!
      p_urb(ino)= 0.0
      d_urb(ino)= r_urb(17)+r_urb(28)+r_urb(29)+r_urb(30)+r_urb(31)
!
      p_urb(ino2)= .95*r_urb(17)+r_urb(27)+.84*r_urb(28)+r_urb(29)   &
              +1.5*r_urb(30)+r_urb(31)+r_urb(32)+r_urb(33)   &
              +1.5*r_urb(34)+r_urb(35)+.5*r_urb(42)
      d_urb(ino2)= r_urb(20)
!
      p_urb(ino3)= 0.0
      d_urb(ino3)= r_urb(6)+r_urb(13)+r_urb(14)+r_urb(19)   &
              +r_urb(32)+r_urb(33)+r_urb(34)+r_urb(35)
!
      p_urb(io3)= 0.0
      d_urb(io3)= r_urb(7)+r_urb(9)+r_urb(10)+r_urb(23)
!
    ! edit wkc, remove!
!      p_urb(ioh)= .12*r_urb(7)+.33*r_urb(9)+.6*r_urb(10)   &
!             +.08*r_urb(23)+r_urb(24)+.23*r_urb(25)
!      d_urb(ioh)= r_urb(1)+r_urb(3)+r_urb(5)   &
!             +r_urb(8)+r_urb(11)+r_urb(12)+r_urb(15)   &
!             +r_urb(16)+r_urb(18)+r_urb(21)+r_urb(25)+r_urb(26)
!
      p_urb(iho2)= r_urb(4)+.22*r_urb(7)+r_urb(8)   &
              +.26*r_urb(9)+.22*r_urb(10)   &
              +r_urb(11)+r_urb(12)+.2*r_urb(15)   &
              +.55*r_urb(16)+.95*r_urb(17)   &
              +.6*r_urb(18)+2*r_urb(21)+r_urb(22)+.76*r_urb(23)   &
              +.9*r_urb(24)+.9*r_urb(27)+.76*r_urb(28)   &
              +.5*r_urb(30)+.9*r_urb(32)+.5*r_urb(34)   &
              +.54*r_urb(40)
      d_urb(iho2)= r_urb(36)+r_urb(37)+r_urb(38)+r_urb(39)
!
      p_urb(ico)= r_urb(4)+r_urb(6)+.24*r_urb(7)   &
             +.31*r_urb(9)+.3*r_urb(10)   &
             +2*r_urb(21)+r_urb(22)+.69*r_urb(23)
      d_urb(ico)= 0.0
!
      p_urb(ich3o2)= r_urb(2)+.07*r_urb(9)+.1*r_urb(10)
      d_urb(ich3o2)= 0.0
!
      p_urb(iethp)= .06*r_urb(9)+.05*r_urb(10)+.1*r_urb(24)   &
               +.1*r_urb(27)   &
               +.08*r_urb(28)+.1*r_urb(32)+.06*r_urb(40)
      d_urb(iethp)= 0.0
!
      p_urb(ihcho)= r_urb(7)+1.56*r_urb(8)+.57*r_urb(9)   &
               +r_urb(11)+r_urb(21)   &
               +.7*r_urb(23)+r_urb(29)+.5*r_urb(30)   &
               +r_urb(33)+.5*r_urb(34)   &
               +.7*r_urb(41)+.5*r_urb(42)
      d_urb(ihcho)= 0.0
!
      p_urb(ich3oh)= .03*r_urb(9)+.04*r_urb(10)
      d_urb(ich3oh)= 0.0
!
      p_urb(iald2)= .22*r_urb(8)+.47*r_urb(9)+1.03*r_urb(10)   &
               +r_urb(11)   &
               +1.77*r_urb(12)+.03*r_urb(23)+.3*r_urb(24)   &
               +.04*r_urb(25)   &
               +.3*r_urb(27)+.25*r_urb(28)+.5*r_urb(30)   &
               +.3*r_urb(32)   &
               +.5*r_urb(34)+.21*r_urb(40)+.5*r_urb(42)
      d_urb(iald2)= 0.0
!
      p_urb(ihcooh)= .52*r_urb(7)+.22*r_urb(9)
      d_urb(ihcooh)= 0.0
!
      p_urb(ircooh)= .09*r_urb(9)+.16*r_urb(10)
      d_urb(ircooh)= 0.0
!
      p_urb(ic2o3)= r_urb(2)+r_urb(4)+r_urb(5)+r_urb(6)   &
               +.13*r_urb(9)+.19*r_urb(10)   &
               +r_urb(21)+r_urb(22)+.62*r_urb(23)   &
               +r_urb(29)+r_urb(33)   &
               +.7*r_urb(41)
      d_urb(ic2o3)= 0.0
!
      p_urb(ipar)= 1.1*r_urb(16)
      d_urb(ipar)= r_urb(1)+r_urb(44)
!
      p_urb(iaone)= .07*r_urb(10)+.23*r_urb(12)   &
               +.74*r_urb(24)+.74*r_urb(27)   &
               +.62*r_urb(28)+.74*r_urb(32)   &
               +.57*r_urb(40)+.15*r_urb(41)
      d_urb(iaone)= r_urb(2)+r_urb(3)
!
      p_urb(imgly)= .04*r_urb(9)+.07*r_urb(10)   &
               +.8*r_urb(16)+.2*r_urb(23)   &
               +.19*r_urb(25)+.15*r_urb(41)
      d_urb(imgly)= r_urb(4)+r_urb(5)+r_urb(6)
!
      p_urb(ieth)= 0.0
      d_urb(ieth)= r_urb(7)+r_urb(8)
!
      p_urb(iolet)= 0.0
      d_urb(iolet)= r_urb(9)+r_urb(11)+r_urb(13)
!
      p_urb(iolei)= 0.0
      d_urb(iolei)= r_urb(10)+r_urb(12)+r_urb(14)
!
      p_urb(itol)= 0.0
      d_urb(itol)= r_urb(15)
!
      p_urb(ixyl)= 0.0
      d_urb(ixyl)= r_urb(16)
!
      p_urb(icres)= .12*r_urb(15)+.05*r_urb(16)
      d_urb(icres)= r_urb(18)+r_urb(19)
!
      p_urb(ito2)= .8*r_urb(15)+.45*r_urb(16)
      d_urb(ito2)= r_urb(17)
!
      p_urb(icro)= .4*r_urb(18)+r_urb(19)
      d_urb(icro)= r_urb(20)
!
      p_urb(iopen)= .95*r_urb(17)+.3*r_urb(18)
      d_urb(iopen)= r_urb(21)+r_urb(22)+r_urb(23)
!
      p_urb(ionit)= .05*r_urb(17)+r_urb(20)   &
               +.16*r_urb(28)+.5*r_urb(30)   &
               +.5*r_urb(34)+r_urb(38)+.5*r_urb(42)
      d_urb(ionit)= r_urb(26)+r_urb(27)
!
      p_urb(irooh)= r_urb(36)+r_urb(37)
      d_urb(irooh)= r_urb(24)+r_urb(25)
!
      p_urb(iro2)= r_urb(1)+.03*r_urb(9)+.09*r_urb(10)   &
              +.77*r_urb(25)
      d_urb(iro2)= r_urb(28)+r_urb(32)+r_urb(36)+r_urb(40)
!
      p_urb(iano2)= r_urb(3)+.11*r_urb(10)
      d_urb(iano2)= r_urb(29)+r_urb(33)+r_urb(37)+r_urb(41)
!
      p_urb(inap)= r_urb(13)+r_urb(14)+r_urb(26)
      d_urb(inap)= r_urb(30)+r_urb(34)+r_urb(38)+r_urb(42)
!
      p_urb(ixo2)= r_urb(5)+r_urb(8)+r_urb(11)+r_urb(12)   &
              +.08*r_urb(15)   &
              +.5*r_urb(16)+.6*r_urb(18)+r_urb(21)+.03*r_urb(23)   &
              +.4*r_urb(24)+.41*r_urb(27)+.34*r_urb(28)   &
              +.4*r_urb(32)+.24*r_urb(40)
      d_urb(ixo2)= r_urb(31)+r_urb(35)+r_urb(39)+r_urb(43)
!
      p_urb(ixpar)= 1.06*r_urb(9)+2.26*r_urb(10)   &
               +r_urb(11)+2.23*r_urb(12)   &
               +1.98*r_urb(24)+.42*r_urb(25)+1.98*r_urb(27)   &
               +1.68*r_urb(28)+r_urb(30)+1.98*r_urb(32)+r_urb(34)   &
               +1.25*r_urb(40)+r_urb(42)
      d_urb(ixpar)= r_urb(44)
!
      return
      end subroutine ode_urb




