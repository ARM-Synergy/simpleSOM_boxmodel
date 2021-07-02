
!***********************************************************************
! subroutine MapGas_Com: maps cnn to and fro stot for the common
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
      subroutine MapGas_Com(stot,imap)

      USE mod_MAIN
      USE mod_GAS

      IMPLICIT NONE

      integer imap
      real(r8) :: stot(ntot_max)

      emit(ih2so4) =emission(kh2so4)
      emit(ihno3)	 =emission(khno3)
      emit(ihcl)	 =emission(khcl)
      emit(inh3)	 =emission(knh3)
      emit(ino)	 =emission(kno)
      emit(ino2)	 =emission(kno2)
      emit(ino3)	 =emission(kno3)
      emit(in2o5)	 =emission(kn2o5)
      emit(ihono)	 =emission(khono)
      emit(ihno4)	 =emission(khno4)
      emit(io3)	 =emission(ko3)
      emit(io1d)	 =emission(ko1d)
      emit(io3p)	 =emission(ko3p)
      emit(ioh)	 =emission(koh)
      emit(iho2)	 =emission(kho2)
      emit(ih2o2)	 =emission(kh2o2)
      emit(ico)	 =emission(kco)
      emit(iso2)	 =emission(kso2)
      emit(ich4)	 =emission(kch4)
      emit(ic2h6)	 =emission(kc2h6)
      emit(ich3o2) =emission(kch3o2)
      emit(iethp)	 =emission(kethp)
      emit(ihcho)	 =emission(khcho)
      emit(ich3oh) =emission(kch3oh)
      emit(ianol)	 =emission(kanol)
      emit(ich3ooh)=emission(kch3ooh)
      emit(iethooh)=emission(kethooh)
      emit(iald2)	 =emission(kald2)
      emit(ihcooh) =emission(khcooh)
      emit(ircooh) =emission(krcooh)
      emit(ic2o3)	 =emission(kc2o3)
      emit(ipan)	 =emission(kpan)

      if(imap.eq.0)then    ! map cnn into stot
      stot(ih2so4)   =cnn(kh2so4)
      stot(ihno3)	 =cnn(khno3)
      stot(ihcl)	 =cnn(khcl)
      stot(inh3)	 =cnn(knh3)
      stot(ino)	     =cnn(kno)
      stot(ino2)	 =cnn(kno2)
      stot(ino3)	 =cnn(kno3)
      stot(in2o5)	 =cnn(kn2o5)
      stot(ihono)	 =cnn(khono)
      stot(ihno4)	 =cnn(khno4)
      stot(io3)	     =cnn(ko3)
      stot(io1d)	 =cnn(ko1d)
      stot(io3p)	 =cnn(ko3p)
      stot(ioh)	     =cnn(koh)
      stot(iho2)	 =cnn(kho2)
      stot(ih2o2)	 =cnn(kh2o2)
      stot(ico)	     =cnn(kco)
      stot(iso2)	 =cnn(kso2)
      stot(ich4)	 =cnn(kch4)
      stot(ic2h6)	 =cnn(kc2h6)
      stot(ich3o2)   =cnn(kch3o2)
      stot(iethp)	 =cnn(kethp)
      stot(ihcho)	 =cnn(khcho)
      stot(ich3oh)   =cnn(kch3oh)
      stot(ianol)	 =cnn(kanol)
      stot(ich3ooh)  =cnn(kch3ooh)
      stot(iethooh)  =cnn(kethooh)
      stot(iald2)	 =cnn(kald2)
      stot(ihcooh)   =cnn(khcooh)
      stot(ircooh)   =cnn(krcooh)
      stot(ic2o3)	 =cnn(kc2o3)
      stot(ipan)	 =cnn(kpan)
!
      else                 ! map stot back into cnn
      cnn(kh2so4)	=stot(ih2so4)
      cnn(khno3)	=stot(ihno3)
      cnn(khcl)	    =stot(ihcl)
      cnn(knh3)	    =stot(inh3)
      cnn(kno)	    =stot(ino)
      cnn(kno2)	    =stot(ino2)
      cnn(kno3)	    =stot(ino3)
      cnn(kn2o5)	=stot(in2o5)
      cnn(khono)	=stot(ihono)
      cnn(khno4)	=stot(ihno4)
      cnn(ko3)	    =stot(io3)
      cnn(ko1d)	    =stot(io1d)
      cnn(ko3p)	    =stot(io3p)
      cnn(koh)	    =stot(ioh)
      cnn(kho2)	    =stot(iho2)
      cnn(kh2o2)	=stot(ih2o2)
      cnn(kco)	    =stot(ico)
      cnn(kso2)	    =stot(iso2)
      cnn(kch4)	    =stot(ich4)
      cnn(kc2h6)	=stot(ic2h6)
      cnn(kch3o2)	=stot(ich3o2)
      cnn(kethp)	=stot(iethp)
      cnn(khcho)	=stot(ihcho)
      cnn(kch3oh)	=stot(ich3oh)
      cnn(kanol)	=stot(ianol)
      cnn(kch3ooh)  =stot(ich3ooh)
      cnn(kethooh)  =stot(iethooh)
      cnn(kald2)	=stot(iald2)
      cnn(khcooh)	=stot(ihcooh)
      cnn(krcooh)	=stot(ircooh)
      cnn(kc2o3)	=stot(ic2o3)
      cnn(kpan)	    =stot(ipan)

      endif

      return
      end subroutine MapGas_Com
