
!***********************************************************************
! subroutine MapGas_Mar: maps cnn to and fro stot for the marine
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
      subroutine MapGas_Mar(stot,imap)
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      integer imap
      real(r8) :: stot(ntot_max)

      emit(idms)	=emission(kdms)
      emit(imsa)	=emission(kmsa)
      emit(idmso)	=emission(kdmso)
      emit(idmso2)	=emission(kdmso2)
      emit(ich3so2h)	=emission(kch3so2h)
      emit(ich3sch2oo)	=emission(kch3sch2oo)
      emit(ich3so2)	=emission(kch3so2)
      emit(ich3so3)	=emission(kch3so3)
      emit(ich3so2oo)	=emission(kch3so2oo)
      emit(ich3so2ch2oo)=emission(kch3so2ch2oo)
      emit(isulfhox)	=emission(ksulfhox)

      if(imap.eq.0)then    ! map cnn into stot
      stot(idms)	=cnn(kdms)
      stot(imsa)	=cnn(kmsa)
      stot(idmso)	=cnn(kdmso)
      stot(idmso2)	=cnn(kdmso2)
      stot(ich3so2h)	=cnn(kch3so2h)
      stot(ich3sch2oo)	=cnn(kch3sch2oo)
      stot(ich3so2)	=cnn(kch3so2)
      stot(ich3so3)	=cnn(kch3so3)
      stot(ich3so2oo)	=cnn(kch3so2oo)
      stot(ich3so2ch2oo)=cnn(kch3so2ch2oo)
      stot(isulfhox)	=cnn(ksulfhox)

      else                 ! map stot back into cnn
      cnn(kdms)		=stot(idms)
      cnn(kmsa)		=stot(imsa)
      cnn(kdmso)	=stot(idmso)
      cnn(kdmso2)	=stot(idmso2)
      cnn(kch3so2h)	=stot(ich3so2h)
      cnn(kch3sch2oo)	=stot(ich3sch2oo)
      cnn(kch3so2)	=stot(ich3so2)
      cnn(kch3so3)	=stot(ich3so3)
      cnn(kch3so2oo)	=stot(ich3so2oo)
      cnn(kch3so2ch2oo)	=stot(ich3so2ch2oo)
      cnn(ksulfhox)	=stot(isulfhox)
      endif

      return
      end subroutine MapGas_Mar
