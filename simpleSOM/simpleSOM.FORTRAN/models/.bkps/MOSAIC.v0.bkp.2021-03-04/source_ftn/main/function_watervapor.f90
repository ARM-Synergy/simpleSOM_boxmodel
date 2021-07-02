!**********************************************************************
! function WaterVapor
! purpose: calculates concentration of h2o using the method given in
!          Seinfeld's book, pg. 181
!----------------------------------------------------------------------

      function WaterVapor(rh, cair_mlc, te, pr_atm)

      use module_data_mosaic_kind, only:  r8
      implicit none

      real(r8) :: WaterVapor
      real(r8) :: rh, cair_mlc, te, pr_atm
      real(r8) :: t_steam, pr_std, a, arg, pr_h2o


      t_steam = 373.15			! steam temperature  [K]
      pr_std   = 1.0			! standard pressure  [atm]

      a      = 1.0 - t_steam/te
      arg    = (((-.1299*a -.6445)*a -1.976)*a +13.3185)*a
      pr_h2o = pr_std*exp(arg)  			! [atm]
      WaterVapor = RH*(pr_h2o/pr_atm)*cair_mlc/100.	! [molec/cc]

      return
      end function WaterVapor

!**********************************************************************
