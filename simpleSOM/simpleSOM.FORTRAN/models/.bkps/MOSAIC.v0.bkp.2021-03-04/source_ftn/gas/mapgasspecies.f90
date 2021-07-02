

!***********************************************************************
! subroutine MapGasSpecies: maps cnn to and fro stot for the selected
!                           gas-phase mechanism.
!
! nomenclature:
! cnn       = full species concentration array.
! stot      = subset of cnn. species concentration array to be supplied to
!             lsodes. length of stot depends on the selected mechanism
! iregime   = selected chemical regime (1-6)
! imap      = 0 : map cnn to stot
!           = 1 : map stot to cnn
!
! author: Rahul A. Zaveri
! date  : february 1996
!
!-------------------------------------------------------------------------
      subroutine MapGasSpecies(stot,imap)
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      integer imap
      real(r8) :: stot(ntot_max)

      goto (1,2,3,4,5,6), iregime
!
!
1     call mapgas_com(stot,imap)
      return
!
!
2     call mapgas_com(stot,imap)
      call mapgas_urb(stot,imap)
      return
!
!
3     call mapgas_com(stot,imap)
      call mapgas_urb(stot,imap)
      call mapgas_bio(stot,imap)
      return
!
!
4     call mapgas_com(stot,imap)
      call mapgas_mar(stot,imap)
      return
!
!
5     call mapgas_com(stot,imap)
      call mapgas_urb(stot,imap)
      call mapgas_mar(stot,imap)
      return
!
!
6     call mapgas_com(stot,imap)
      call mapgas_urb(stot,imap)
      call mapgas_bio(stot,imap)
      call mapgas_mar(stot,imap)
      return
!
      end subroutine MapGasSpecies
