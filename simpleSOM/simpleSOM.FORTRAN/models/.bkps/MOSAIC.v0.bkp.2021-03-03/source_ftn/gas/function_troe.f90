      function Troe(cair_mlc,te,rk0,rnn,rki,rmm)

      use module_data_mosaic_kind, only:  r8
      implicit none

      real(r8) :: Troe
      real(r8) :: cair_mlc, te, rk0, rnn, rki, rmm
      real(r8) :: expo

      rk0 = rk0*cair_mlc*(te/300.)**(-rnn)
      rki = rki*(te/300.)**(-rmm)
      expo= 1./(1. + (LOG10(rk0/rki))**2)
      troe  = (rk0*rki/(rk0+rki))*.6**expo

      return
      end function Troe
