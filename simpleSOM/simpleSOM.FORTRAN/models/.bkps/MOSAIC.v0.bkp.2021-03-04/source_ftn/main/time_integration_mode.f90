      subroutine time_integration_mode
      use module_data_mosaic_main
      use module_data_mosaic_gas
      use module_data_mosaic_aero
      use module_data_mosaic_cloud

      use module_pmcmos_subaa, only:  pmcmos_update_met, pmcmos_emit_dilu
      use module_pmcmos_aer, only:  pmcmos_print

      implicit none

! local variables
      integer mfreq_aeroptic

      REAL(R8) timer

      it = 0			! initialize time step counter
      if(mgas == mYES)then
        call print_gas		! print inital gas concentrations
      endif

!------------------------------------------------------------------
!
! main "time" loop begins...

      mfreq_aeroptic = max( 1, nint(dt_aeroptic_min/dt_min) )
      
      if ((maeroptic > 0) .and. (maer == myes)) then
        call aerosol_optical
      endif

!      call pmcmos_print( 0 )
      
      
      do 100 it = 1, nstep

         timer = time_sec

         IF (MOD(timer,600.).EQ.0) THEN
            WRITE(*,'(A,F10.1)') 'time = ', timer
         END IF

      call UpdateTime

      if(msolar.eq.1)then
       call SolarZenithAngle
      endif
      
      
      if (ipmcmos <=  0) then
         call UpdateMetFields
         call UpdateEmissions
      else
         call pmcmos_update_met( 1 )
         call pmcmos_emit_dilu

!        write(*,'(1p,a,3e14.6)') 'cos_sza,tcur,tmid', cos_sza, tcur_sec, tmid_sec
!        write(*,'(a,1p,5e13.5)') 'co,so2,ch4,c2h6', cnn(17:20)*ppb/cair_mlc
      endif

      
      call IntegrateChemistry(timer)
      

!      if ((ipmcmos >  0) .and. &
!          (mod(it,iprint) == 0)) then
!         write(lun_sect_190,'(/a,5(1x,f14.6))') 'time(h)', time_sec/3600.0, &
!            te, pr_atm*1013.25, rh, pblh
!         write(lun_sect_190,'(1p,5e14.6)') cnn(1:ngas_max)*ppb/cair_mlc
!      endif

!      call DoMassBalance

      if ((maeroptic > 0) .and. (maer == myes)) then
        if (mod(it,mfreq_aeroptic) .eq. 0) call aerosol_optical
      endif


      if(mgas == mYES)then
        	  if(mod(it,iprint) .eq. 0)call print_gas
      endif

!      if(maer == mYES)then		! UNCOMMENT THIS LINE
!        call print_aer(1)		! UNCOMMENT THIS LINE
!      endif				! UNCOMMENT THIS LINE

!      call pmcmos_print( it )

100   continue	! time loop



!------------------------------------------------------------------


      return
      end subroutine time_integration_mode




