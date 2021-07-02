
!***********************************************************************
! subroutine SetGas_Bio: sets up gas-phase species indices for
! the selected mechanism.
!
! author: Rahul A. Zaveri
! date  : february 1996
!-------------------------------------------------------------------------
      subroutine SetGas_Bio(ilast)
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      integer ilast

      iisop	    = ilast + 1
      iisoprd	= ilast + 2
      iisopp	= ilast + 3
      iisopn	= ilast + 4
      iisopo2	= ilast + 5

	  icn3	 = ilast + 6
      icn2	 = ilast + 7
      icn1   = ilast + 8
      ic0    = ilast + 9
      ic1 	 = ilast + 10
      ic2	 = ilast + 11
      ic3	 = ilast + 12
      ic4	 = ilast + 13
      ic5    = ilast + 14
      ic6    = ilast + 15
      ic7	 = ilast + 16
      ic8    = ilast + 17
      ic9    = ilast + 18

      ilast	= ic9

	  ! edit wkc
	  ! print *, 'ilast after setgas_bio = ', ilast

      return
      end subroutine SetGas_Bio

