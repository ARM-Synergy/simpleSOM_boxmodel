
!***********************************************************************
! subroutine MapGas_Bio: maps cnn to and from stot for the biogenic
!                        gas-phase mechanism.
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
      subroutine MapGas_Bio(stot,imap)
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      integer imap
      real(r8) :: stot(ntot_max)

      emit(iisop)	=emission(kisop)
      emit(iisoprd)	=emission(kisoprd)
      emit(iisopp)	=emission(kisopp)
      emit(iisopn)	=emission(kisopn)
      emit(iisopo2)	=emission(kisopo2)
      emit(iCn3)  = emission(kCn3)
      emit(iCn2)  = emission(kCn2)
      emit(iCn1)  = emission(kCn1)
      emit(iC0)  = emission(kC0)
      emit(iC1)  = emission(kC1)
      emit(iC2)  = emission(kC2)
      emit(iC3)  = emission(kC3)
      emit(iC4)  = emission(kC4)
      emit(iC5)  = emission(kC5)
      emit(iC6)  = emission(kC6)
      emit(iC7)  = emission(kC7)
      emit(iC8)  = emission(kC8)
      emit(iC9)  = emission(kC9)

      if(imap.eq.0)then    ! map cnn into stot
      stot(iisop)	=cnn(kisop)
      stot(iisoprd)	=cnn(kisoprd)
      stot(iisopp)	=cnn(kisopp)
      stot(iisopn)	=cnn(kisopn)
      stot(iisopo2)	=cnn(kisopo2)
	  stot(iCn3)    = cnn(kCn3)
      stot(iCn2)    = cnn(kCn2)
      stot(iCn1)    = cnn(kCn1)
      stot(iC0)     = cnn(kC0)
      stot(iC1)     = cnn(kC1)
      stot(iC2)     = cnn(kC2)
      stot(iC3)     = cnn(kC3)
      stot(iC4)     = cnn(kC4)
      stot(iC5)     = cnn(kC5)
      stot(iC6)     = cnn(kC6)
      stot(iC7)     = cnn(kC7)
      stot(iC8)     = cnn(kC8)
      stot(iC9)     = cnn(kC9)
!
      else                 ! map stot back into cnn
      cnn(kisop)	=stot(iisop)
      cnn(kisoprd)	=stot(iisoprd)
      cnn(kisopp)	=stot(iisopp)
      cnn(kisopn)	=stot(iisopn)
      cnn(kisopo2)	=stot(iisopo2)
	  cnn(kCn3)  = stot(iCn3)
      cnn(kCn2)  = stot(iCn2)
      cnn(kCn1)  = stot(iCn1)
      cnn(kC0)  = stot(iC0)
      cnn(kC1)  = stot(iC1)
      cnn(kC2)  = stot(iC2)
      cnn(kC3)  = stot(iC3)
      cnn(kC4)  = stot(iC4)
      cnn(kC5)  = stot(iC5)
      cnn(kC6)  = stot(iC6)
      cnn(kC7)  = stot(iC7)
      cnn(kC8)  = stot(iC8)
      cnn(kC9)  = stot(iC9)


      endif

      return
      end subroutine MapGas_Bio
