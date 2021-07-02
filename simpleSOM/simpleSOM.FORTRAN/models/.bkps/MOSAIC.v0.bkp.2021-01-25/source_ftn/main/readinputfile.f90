
      subroutine ReadInputFile
  
      use module_data_mosaic_main
      use module_data_mosaic_gas
      use module_data_mosaic_aero
      use module_data_mosaic_asect
      use module_data_mosaic_cloud

      use module_data_mosaic_pmcmos, only:  &
         temp_profile_fname, pblh_profile_fname,  &
         aer_init_fname, aer_back_fname, aer_emit_fname,  &
         gas_init_fname, gas_back_fname, gas_emit_fname,  &
         msolar_pmcmos_dtshift
      use module_pmcmos_init, only:  &
         pmcmos_init

      implicit none

      integer i, k, ibin, idum, input_version, noffset
      character*40 dword
      real(r8) :: tmpa, tmpb

! input_version = 1 is pre-april-2008
! input_version = 2 is does input of new variables for sectional aerosol
      read(lun_inp,*)dword

      if (dword == 'INPUT_VERSION') then
         read(lun_inp,*)input_version
         read(lun_inp,*)dword
      else
         input_version = 1
      end if

      

      
      if ( (input_version /=    1) .and. &
           (input_version /=    2) .and. &
           (input_version /= 1002) .and. &
           (input_version /= 2002) .and. &
           (input_version /= 2003) .and. &
           (input_version /= 2004) ) then
         write(*,'(2a)') &
            '*** readinputfile fatal error - ', &
            'input_version must be 1/2, 1002, or 2002/3'
         stop
      end if

      read(lun_inp,*)tbeg_mo, tbeg_dd, tbeg_hh, tbeg_mm, tbeg_ss ! UTC
      read(lun_inp,*)trun_dd, trun_hh, trun_mm, trun_ss ! run time
      read(lun_inp,*)dt_min		! transport time-step [min]
      read(lun_inp,*)dt_aeroptic_min	! time-step for aerosol optics calcs [min]
      read(lun_inp,*)rlon, rlat		! lon and lat [deg]
      read(lun_inp,*)zalt_m		! altitude MSL [m]
      read(lun_inp,*)RH			! relative humidity [%]
      read(lun_inp,*)te			! temperature [K]
      read(lun_inp,*)pr_atm		! pressure [atm]

      
      if (pr_atm < 0) then		! if < 0, convert from [Pa] to [atm]
         pr_atm = -pr_atm/press0_pa
      end if

      ntype_aer = 1
      ntype_md1_aer = 0
      ntype_md2_aer = 0
      msectional_flag2 = 0
      method_bcfrac = 1
      method_kappa = 12

      if (input_version == 1) then
         read(lun_inp,*)nbin_a		! number of aerosol bins
         nsize_aer(:) = nbin_a
         msize_framework = munstructured
         maersize_init_flag1 = 1
         mhyst_method = 1
         mcoag_flag1 = 0
         ifreq_coag = 1
         mmovesect_flag1 = 0
         mnewnuc_flag1 = 0
         msectional_flag1 = 0
      else

         
         if (input_version < 2000) then
            read(lun_inp,*) ntype_aer
         else
            read(lun_inp,*) ntype_md1_aer, ntype_md2_aer
            if ((ntype_md1_aer < 1) .or. (ntype_md1_aer > maxd_atype_md1)) then
               write(*,'(2a,2i10)') &
                  '*** readinputfile fatal error - ', &
                  'bad ntype_md1_aer', ntype_md1_aer, maxd_atype_md1
               stop
            else if ((ntype_md2_aer < 1) .or. (ntype_md2_aer > maxd_atype_md2)) then
               write(*,'(2a,2i10)') &
                  '*** readinputfile fatal error - ', &
                  'bad ntype_md2_aer', ntype_md2_aer, maxd_atype_md2
               stop
            end if
            ntype_aer = ntype_md1_aer*ntype_md2_aer

            ! msectional_flag2 = 1 indicates multi-dimensional sectional
            if (ntype_aer > 1) msectional_flag2 = 1
         end if

         if (input_version >= 2003) then
            read(lun_inp,*) method_bcfrac, method_kappa
         end if

         read(lun_inp,*) nsize_aer(1)
         nsize_aer(2:) = nsize_aer(1)
         
         
         nbin_a = ntype_aer*nsize_aer(1)

         read(lun_inp,*) msize_framework
         read(lun_inp,*) maersize_init_flag1
         read(lun_inp,*) mhyst_method
         
         if (input_version < 2004) then
            read(lun_inp,*) mcoag_flag1
            ifreq_coag = 1
         else
            read(lun_inp,*) mcoag_flag1, ifreq_coag
            ifreq_coag = max( 1, ifreq_coag )
         end if
         
         read(lun_inp,*) mmovesect_flag1
         read(lun_inp,*) mnewnuc_flag1
         read(lun_inp,*) msectional_flag1
      end if

      if (nbin_a > nbin_a_max) then
         write(*,'(2a,2i10)') &
            '*** readinputfile fatal error - ', &
            'nbin_a > nbin_a_max', nbin_a, nbin_a_max
         stop
      else if (ntype_aer > maxd_atype) then
         write(*,'(2a,2i10)') &
            '*** readinputfile fatal error - ', &
            'ntype_aer > maxd_atype', ntype_aer, maxd_atype
         stop
      else if (maxval(nsize_aer(:)) > maxd_asize) then
         write(*,'(2a,2i10)') &
            '*** readinputfile fatal error - ', &
            'maxval(nsize_aer(:)) > maxd_asize', maxval(nsize_aer(:)), maxd_asize
         stop
      end if

      if (method_bcfrac /= 1) then
         write(*,'(2a,2i10)') &
            '*** readinputfile fatal error - ', &
            'bad method_bcfrac', method_bcfrac
         stop
      else if ((method_kappa /= 11) .and. (method_kappa /= 12)) then
         write(*,'(2a,2i10)') &
            '*** readinputfile fatal error - ', &
            'bad method_kappa', method_kappa
         stop
      end if

      naerbin_used = nbin_a
      ncldbin_used = 0
      
      if (ncldbin_used > 0) then
         ntot_used = ngas_max + naer_tot*naerbin      &
                              + ncld_tot*ncldbin_used
      else
         ntot_used = ngas_max + naer_tot*naerbin_used
      end if

      
      read(lun_inp,*) iprint		! freq of output
      
      
      if (input_version >= 2003) then
         read(lun_inp,*) iwrite_gas, iwrite_aer_bin, &
                         iwrite_aer_dist, iwrite_aer_species
      else if (nbin_a <= 1001) then
         iwrite_gas = 1 ; iwrite_aer_bin = 1
         iwrite_aer_dist = 1 ; iwrite_aer_species = 1
      else
         iwrite_gas = 0 ; iwrite_aer_bin = 0
         iwrite_aer_dist = 0 ; iwrite_aer_species = 0
      end if

      read(lun_inp,*)mmode		! flag
      read(lun_inp,*)mgas		! flag
      read(lun_inp,*)maer		! flag
      read(lun_inp,*)mcld		! flag
      read(lun_inp,*)maeroptic		! flag
      read(lun_inp,*)mshellcore		! flag
      read(lun_inp,*)msolar		! flag
      read(lun_inp,*)mphoto		! flag
      read(lun_inp,*)mGAS_AER_XFER	! flag
      read(lun_inp,*)mDYNAMIC_SOLVER 	! flag
      read(lun_inp,*)alpha_ASTEM	! tolerance for tau
      read(lun_inp,*)rtol_eqb_ASTEM	! relative eqb tolerance
      read(lun_inp,*)ptol_mol_ASTEM	! percent mol tolerance
!	read(lun_inp,*)rk_dil		! first order dilution rate or wall loss [1/s]

! read gas info
      read(lun_inp,*)dword ! 1st heading line
      read(lun_inp,*)dword ! 2nd heading line
      read(lun_inp,*)dword ! 3rd heading line = "GAS"

!      WRITE(*,'(A,1X,(I2))') 'nGAS_MAX = ',nGAS_MAX
!      STOP

      
      do i=1, ngas_max
! read index, species name, initial conc, and emissions
        read(lun_inp,*)k, species(k), cnn(k), emission(k)
        if (k.NE.i) then
           write(*,'(2a,2(1x,i8))') &
              '*** readinputfile fatal error', &
              ' - gas indices i & k differ', i, k
           stop
        end if
        cnn(k) = max( cnn(k), 0.0_r8 )
      enddo

! read pmcmos stuff
      if (input_version > 1000) then
         read(lun_inp,*) ipmcmos
      else
         ipmcmos = 0
      end if
      
      
      if ((input_version >= 2000) .and. (ipmcmos <= 0)) then
         write(*,'(2a)') &
            '*** readinputfile fatal error - ', &
            'input_version >= 2000 requires ipmcmos > 0'
         stop
      end if

      if (ipmcmos > 0) then
         read(lun_inp,'(a)') temp_profile_fname
         read(lun_inp,'(a)') pblh_profile_fname
         read(lun_inp,'(a)') gas_init_fname
         read(lun_inp,'(a)') gas_back_fname
         read(lun_inp,'(a)') gas_emit_fname
         read(lun_inp,'(a)') aer_init_fname
         read(lun_inp,'(a)') aer_back_fname
         read(lun_inp,'(a)') aer_emit_fname
         call pmcmos_init
         emission(:) = 0.0   ! use the pmcmos emissions
      else
         te_old = te
         pr_atm_old = pr_atm
         rh_old = rh
         pblh = 1.0e3   ! default pblh (m)
         pblh_old = pblh
      end if
      istate_pblh = 0

      msolar_pmcmos_dtshift = 0
      if (msolar >= 1000) then
         if (ipmcmos > 0) msolar_pmcmos_dtshift = 1
         msolar = mod( msolar, 1000 )
      end if

! read aerosol info
      dword = ' '
      do while (dword == ' ')
         read(lun_inp,'(a)')dword ! 1st heading line = "AEROSOL"
      end do


      dlo_aersize_init = 0.0 ; dhi_aersize_init = 0.0
      
      if (input_version > 1) then
         
         read(lun_inp,*) dlo_aersize_init, dhi_aersize_init
         
         if (input_version >= 2000) then
             read(lun_inp,*) method_atype_md1_init
             if (method_atype_md1_init <= 1) then
                read(lun_inp,*) xcut_atype_md1(0:ntype_md1_aer)
             else
                read(lun_inp,*) xcutlo_atype_md1_init, xcuthi_atype_md1_init
                tmpa = xcutlo_atype_md1_init + &
                      (xcuthi_atype_md1_init-xcutlo_atype_md1_init)/ntype_md1_aer
                tmpb = max( 0.1_r8, xcutlo_atype_md1_init+0.1 )
                if ( (tmpa < 0.0) .or. (xcuthi_atype_md1_init < tmpb) ) then
                   write(*,'(2a,1p,2e14.6)') &
                      '*** readinputfile fatal error - ', &
                      'bad xcutlo/hi_atype_md1_init', &
                      xcutlo_atype_md1_init, xcuthi_atype_md1_init
                   stop
                end if
             end if

             read(lun_inp,*) method_atype_md2_init
             if (method_atype_md2_init <= 1) then
                read(lun_inp,*) xcut_atype_md2(0:ntype_md2_aer)
             else
                read(lun_inp,*) xcutlo_atype_md2_init, xcuthi_atype_md2_init
                if ( (xcutlo_atype_md2_init < 1.0e-7) .or. &
                     (xcuthi_atype_md2_init < xcutlo_atype_md2_init*1.1) ) then
                   write(*,'(2a,1p,2e14.6)') &
                      '*** readinputfile fatal error - ', &
                      'bad xcutlo/hi_atype_md2_init', &
                      xcutlo_atype_md2_init, xcuthi_atype_md2_init
                   stop
                end if
             end if
         end if
      end if
      

      if (ipmcmos <= 0) then
! when ipmcmos > 0, this "manual" initialization of the aerosol
! is not necessary
         read(lun_inp,*)dword ! 2nd heading line
         
         read(lun_inp,*)dword ! 3rd heading line
         
         do ibin = 1, nbin_a
           noffset = ngas_max + naer_tot*(ibin - 1)

            
            
           read(lun_inp,*)idum, (cnn(k+noffset), k = 1, naer_tot)

           
           do k = 1, naer_tot
              cnn(k+noffset) = max( cnn(k+noffset), 0.0_r8 )

              !PRINT*, k,cnn(k+noffset)
              !READ(*,*)
              
           enddo
          
         enddo
      end if ! (ipmcmos <= 0)


! read cloud file
!      read(lin,*)dword ! CLOUD1
!      do i=1, ncld_tot
!      read(lin,*)k, species(k), cnn(k), tmpa
!      enddo


      write(6,*)'Finished reading all inputs...'
      return
      end subroutine ReadInputFile
