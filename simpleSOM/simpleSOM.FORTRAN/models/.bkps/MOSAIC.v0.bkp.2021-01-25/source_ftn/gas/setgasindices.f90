

!***********************************************************************
! subroutine SetGasIndices: sets up gas-phase species indices for
! the selected mechanism.
!
! input: iregime    = 1     : com
!                   = 2     : com + urb
!                   = 3     : com + urb + bio
!                   = 4     : com + mar
!                   = 5     : com + urb + mar
!                   = 6     : com + urb + bio + mar
!
! author: Rahul A. Zaveri
! date  : february 1996
!-------------------------------------------------------------------------
      subroutine SetGasIndices
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      integer ilast
!
      ilast = 0

      goto (1,2,3,4,5,6), iregime

1     call setgas_com(ilast)
      return
!
!
2     call setgas_com(ilast)
      call setgas_urb(ilast)
      return
!
!
3     call setgas_com(ilast)
      call setgas_urb(ilast)
      call setgas_bio(ilast)
      return
!
!
4     call setgas_com(ilast)
      call setgas_mar(ilast)
      return
!
!
5     call setgas_com(ilast)
      call setgas_urb(ilast)
      call setgas_mar(ilast)
      return
!
!
6     call setgas_com(ilast)
      call setgas_urb(ilast)
      call setgas_bio(ilast)
      call setgas_mar(ilast)
      return
!
      end subroutine SetGasIndices
