
!***********************************************************************
! subroutine SetGas_Urb: sets up gas-phase species indices for
! the selected mechanism.
!
! author: Rahul A. Zaveri
! date  : february 1996
!-------------------------------------------------------------------------
      subroutine SetGas_Urb(ilast)
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      integer ilast

      ipar	= ilast + 1
      iaone	= ilast + 2
      imgly	= ilast + 3
      ieth	= ilast + 4
      iolet	= ilast + 5
      iolei	= ilast + 6
      itol	= ilast + 7
      ixyl	= ilast + 8
      icres	= ilast + 9
      ito2	= ilast + 10
      icro	= ilast + 11
      iopen	= ilast + 12
      ionit	= ilast + 13
      irooh	= ilast + 14
      iro2	= ilast + 15
      iano2	= ilast + 16
      inap	= ilast + 17
      ixo2	= ilast + 18
      ixpar	= ilast + 19

      ilast	= ixpar

	  ! edit wkc
	  ! print *, 'ilast after setgas_urb = ', ilast
	  
      return
      end subroutine SetGas_Urb
