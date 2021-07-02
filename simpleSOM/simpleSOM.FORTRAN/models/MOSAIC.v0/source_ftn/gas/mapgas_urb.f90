
!***********************************************************************
! subroutine MapGas_Urb: maps cnn to and fro stot for the urban
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
      subroutine MapGas_Urb(stot,imap)
  
      USE mod_MAIN
      USE mod_GAS

      IMPLICIT NONE

      integer imap
      real(r8) :: stot(ntot_max)

      emit(ipar)	=emission(kpar)
      emit(iaone)	=emission(kaone)
      emit(imgly)	=emission(kmgly)
      emit(ieth)	=emission(keth)
      emit(iolet)	=emission(kolet)
      emit(iolei)	=emission(kolei)
      emit(itol)	=emission(ktol)
      emit(ixyl)	=emission(kxyl)
      emit(icres)	=emission(kcres)
      emit(ito2)	=emission(kto2)
      emit(icro)	=emission(kcro)
      emit(iopen)	=emission(kopen)
      emit(ionit)	=emission(konit)
      emit(irooh)	=emission(krooh)
      emit(iro2)	=emission(kro2)
      emit(iano2)	=emission(kano2)
      emit(inap)	=emission(knap)
      emit(ixo2)	=emission(kxo2)
      emit(ixpar)	=emission(kxpar)


      if(imap.eq.0)then    ! map cnn into stot
      stot(ipar)	=cnn(kpar)
      stot(iaone)	=cnn(kaone)
      stot(imgly)	=cnn(kmgly)
      stot(ieth)	=cnn(keth)
      stot(iolet)	=cnn(kolet)
      stot(iolei)	=cnn(kolei)
      stot(itol)	=cnn(ktol)
      stot(ixyl)	=cnn(kxyl)
      stot(icres)	=cnn(kcres)
      stot(ito2)	=cnn(kto2)
      stot(icro)	=cnn(kcro)
      stot(iopen)	=cnn(kopen)
      stot(ionit)	=cnn(konit)
      stot(irooh)	=cnn(krooh)
      stot(iro2)	=cnn(kro2)
      stot(iano2)	=cnn(kano2)
      stot(inap)	=cnn(knap)
      stot(ixo2)	=cnn(kxo2)
      stot(ixpar)	=cnn(kxpar)

!
      else                 ! map stot back into cnn
      cnn(kpar)		=stot(ipar)
      cnn(iaone)	=stot(kaone)
      cnn(imgly)	=stot(kmgly)
      cnn(ieth)		=stot(keth)
      cnn(iolet)	=stot(kolet)
      cnn(iolei)	=stot(kolei)
      cnn(itol)		=stot(ktol)
      cnn(ixyl)		=stot(kxyl)
      cnn(icres)	=stot(kcres)
      cnn(ito2)		=stot(kto2)
      cnn(icro)		=stot(kcro)
      cnn(iopen)	=stot(kopen)
      cnn(ionit)	=stot(konit)
      cnn(irooh)	=stot(krooh)
      cnn(iro2)		=stot(kro2)
      cnn(iano2)	=stot(kano2)
      cnn(inap)		=stot(knap)
      cnn(ixo2)		=stot(kxo2)
      cnn(ixpar)	=stot(kxpar)
      endif

      return
      end subroutine MapGas_Urb
