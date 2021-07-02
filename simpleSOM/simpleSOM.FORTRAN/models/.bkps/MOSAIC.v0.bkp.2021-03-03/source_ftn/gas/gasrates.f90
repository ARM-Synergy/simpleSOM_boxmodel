!**************************************************************************
! subroutine gasrate: calculates reaction rates for the selected mechanism
!
! nomenclature:
! r_com, r_urb, r_bio, r_mar = reaction rates (molec/cc/sec)
! rk_com,rk_urb,rk_bio,rk_mar= rate constants in appropriate units
! s                          = species concentrations (molec/cc)
! o2                         = oxygen concentration   (molec/cc)
! cair_mlc (used for M)      = air concentration      (molec/cc)
! h2o                        = water vapor            (molec/cc)
!
! author: Rahul A. Zaveri
! date  : february 1996
!
!--------------------------------------------------------------------------
      subroutine GasRates(s)
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      real(r8) :: s(ntot_max)


      call gasrates_het(s)

      goto (1,2,3,4,5,6), iregime

1     call gasrates_com(s)
      return

2     call gasrates_com(s)
      call gasrates_urb(s)
      return

3     call gasrates_com(s)
      call gasrates_urb(s)
      call gasrates_bio(s)
      return

4     call gasrates_com(s)
      call gasrates_mar(s)
      return

5     call gasrates_com(s)
      call gasrates_urb(s)
      call gasrates_mar(s)
      return

6     call gasrates_com(s)
      call gasrates_urb(s)
      call gasrates_bio(s)
      call gasrates_mar(s)
      return

      end subroutine GasRates
