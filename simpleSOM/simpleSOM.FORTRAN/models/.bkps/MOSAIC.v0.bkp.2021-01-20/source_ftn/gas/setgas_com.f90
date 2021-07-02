
!***********************************************************************
! subroutine SetGas_Com: sets up gas-phase species indices for
! the selected mechanism.
!
! author: Rahul A. Zaveri
! date  : february 1996
!-------------------------------------------------------------------------
      subroutine SetGas_Com(ilast)
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      integer ilast

      ih2so4 = 1
      ihno3	 = 2
      ihcl	 = 3
      inh3	 = 4
      ino	 = 5
      ino2	 = 6
      ino3	 = 7
      in2o5	 = 8
      ihono	 = 9
      ihno4	 = 10
      io3	 = 11
      io1d	 = 12
      io3p	 = 13
      ioh	 = 14
      iho2	 = 15
      ih2o2	 = 16
      ico	 = 17
      iso2	 = 18
      ich4	 = 19
      ic2h6	 = 20
      ich3o2 = 21
      iethp	 = 22
      ihcho	 = 23
      ich3oh = 24
      ianol	 = 25
      ich3ooh= 26
      iethooh= 27
      iald2	 = 28
      ihcooh = 29
      ircooh = 30
      ic2o3	 = 31
      ipan	 = 32

      ilast	= ipan
	
	! edit wkc
	! print *, 'ilast after setgas_com = ', ilast	
		
      return
      end subroutine SetGas_Com

