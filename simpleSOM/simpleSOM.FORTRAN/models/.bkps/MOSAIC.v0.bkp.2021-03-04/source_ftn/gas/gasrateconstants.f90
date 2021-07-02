!**************************************************************************
! subroutine GasRateConstants: generates thermal rate coefficients
!                   for the selected mechanism
! nomenclature:
! rk_com    = reaction rate constants for common mechanism (molec-cc-s)
! rk_urb    = reaction rate constants for hc1 mechanism    (molec-cc-s)
! rk_bio    = reaction rate constants for hc2 mechanism    (molec-cc-s)
! rk_mar    = reaction rate constants for marine mechanism (molec-cc-s)
! te        = ambient atmospheric temperature (K)
! iregime = selected mechanism for the current chemical regime (1-6)
!
! author: Rahul A. Zaveri
! date  : february 1996
!
!-------------------------------------------------------------------------
      subroutine GasRateConstants
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none



      if(msolar.eq.1)then
       call SolarZenithAngle		! calculates cos_sza
       call PhotoConstants_Solar	! natural diurnal variation
      elseif(msolar.eq.2)then
       call PhotoConstants_Fixed	! artificial as in a smog chamber
      endif

      call gasrateconstants_het

      goto (1,2,3,4,5,6), iregime

1     call gasrateconstants_com
      return

2     call gasrateconstants_com
      call gasrateconstants_urb
      return

3     call gasrateconstants_com
      call gasrateconstants_urb
      call gasrateconstants_bio
      return

4     call gasrateconstants_com
      call gasrateconstants_mar
      return

5     call gasrateconstants_com
      call gasrateconstants_urb
      call gasrateconstants_mar
      return

6     call gasrateconstants_com
      call gasrateconstants_urb
      call gasrateconstants_bio
      call gasrateconstants_mar
      return
!
      end subroutine GasRateConstants
