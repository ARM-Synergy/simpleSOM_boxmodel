! zz01aerchemistry.f (mosaic.25.0)
!********************************************************************************************
!   code history
!   8/28/2013 raz - added code to calculate reactive kl
!   8/27/2013 raz - bug fix in MESA_estimate_eleliquid(ibin,XT)
!   12/24/2012 raz - moved wall loss subroutine to aerchemistry.f90.
!   6/25/2008  raz - updated nh4no3 and nh4cl condensation algorithm
!   01-may-07 raz - updated CRH and hysteresis treatment for cano3 and cacl2 salts
!   09-jan-07 raz - major clean up of variables and subroutines
!   25-sep-06 raz - added kelvin effect treatment for condensing species
!   22-sep-06 raz - changed "min" to "max" in ratio_AN and ratio_AC definitions
!   21-jul-06 raz - revised and debugged kelvin effect algorithm
!   17-jun-06 raz - added MSA chemistry in particle phase
!   06-jan-05 raz - implemented revised ASTEM algorithm
!   08-oct-05 raz - debugged
!   21-sep-05 raz - revised adaptive time stepping scheme in MESA.
!   28-apr-05 raz - reversed calls to form_cacl2 and form_nacl
!                   fixed caco3 error in subr. electrolytes_to_ions
!                   renamed dens_aer to dens_aer_mac; mw_aer to mw_aer_mac
!   27-apr-05 raz - updated dry_mass calculation approach in MESA_convergence
!   22-apr-05 raz - fixed CaSO4 mass balance problem and updated algorithm to
!                   calculate phi_volatile for nh3, hno3, and hcl.
!   20-apr-05 raz - updated ASCEEM
!   19-apr-05 raz - updated the algorithm to constrain the nh4 concentration
!                   during simultaneous nh3, hno3, and hcl integration such
!                   that it does not exceed the max possible value for a given bin
!   14-apr-05 raz - fixed ASTEM_flux_wet_case3 and ASTEM_flux_dry_case3c
!   11-apr-05 raz - added SOA based on SORGAM mechanism
!   11-jan-05 raz - major updates to many subroutines
!   18-nov-04 rce - make sure that acos argument is between +/-1.0
!   28-jan-04 rce - added subr aerchem_boxtest_output;
!	eliminated some unnecessary "include v33com-"
!   01-dec-03 rce - added "implicit none" to many routines;
!	eliminated some unnecessary "include v33com-"
!   05-oct-03 raz - added hysteresis treatment
!   02-sep-03 raz - implemented ASTEM
!   10-jul-03 raz - changed ix to ixd in interp. subrs fast*_up and fast*_lo
!   08-jul-03 raz - implemented ASTEM (adaptive step time-split
!                   explicit euler method)
!   26-jun-03 raz - updated almost all the subrs. this version contains
!       options for rigorous and fast solvers (including lsode solver)
!
!   07-oct-02 raz - made zx and zm integers in activity coeff subs.
!   16-sep-02 raz - updated many subrs to treat calcium salts
!   19-aug-02 raz - inlcude v33com9a in subr aerosolmtc
!   14-aug-02 rce - "(msectional.eq.0)" changed to "(msectional.le.0)"
!   07-aug-02 rce - this is rahul's latest version from freshair
!	AFTER adding "real mean_molecular_speed" wherever it is used
!   01-apr-02 raz - made final tests and gave the code to jerome
!
!   04--14-dec-01 rce - several minor changes during initial testing/debug
!	in 3d los angeles simulation
!	(see earlier versions for details about these changes)
!-----------------------------------------------------------------------
!23456789012345678901234567890123456789012345678901234567890123456789012

!***********************************************************************
! MOSAIC (Model for Simulating Aerosol Interactions and Chemistry)
!
! author: Rahul A. Zaveri
! update: dec 2004
!-----------------------------------------------------------------------




      subroutine mosaic_box_aerchemistry( dtchem )
  
      use module_data_mosaic_main
      use module_data_mosaic_aero
      use module_sect_iface, only:  sectional_interface_1

      implicit none

!   subr arguments
      real(r8) :: dtchem


      it_mosaic = it
	iprint_mosaic = iprint
      RH_pc = RH					! RH(%)
      aH2O = 0.01*RH_pc				! aH2O (aerosol water activity)
      P_atm = pr_atm				! P(atm)
      T_K = te					! T(K)
      cair_mol_m3 = cair_molm3		! air conc in mol/m3
      cair_mol_cc = cair_mol_m3*1.e-6	! air conc in mol/cc

      !PRINT*, 'HERE3'
      !PRINT*, Dp_dry_a
      !READ(*,*)

      
      call load_mosaic_parameters		! sets up indices and other stuff once per simulation
      
      call update_thermodynamic_constants	! update T and RH dependent constants

      !PRINT*, 'HERE3-B'
      !PRINT*, Dp_dry_a
      !READ(*,*)
      
      do irepeat_mosaic = 1, 1

         call initialize_mosaic_variables

         !PRINT*, 'HERE1'
         !PRINT*, Dp_dry_a
         !READ(*,*)

         call map_mosaic_species_BOX(0)

         !PRINT*, 'HERE2'
         !PRINT*, Dp_dry_a
         !READ(*,*)
        
         call overall_massbal_in

         !PRINT*, 'HERE3-C'
         !PRINT*, Dp_dry_a
         !READ(*,*)

         
         !PRINT*, 'HERE3'
         !PRINT*, Dp_dry_a
         !READ(*,*)


         !!!!!!!!! MAKR HERE
         AER(23:24,3,:) = 0.0
        
         call MOSAIC_dynamic_solver( dtchem )

         !PRINT*, 'HERE4'
         !PRINT*, Dp_dry_a
         !READ(*,*)
        
         call overall_massbal_out(0)

         !PRINT*, 'HERE5'
         !PRINT*, Dp_dry_a
         !READ(*,*)
        
      enddo


      
      call map_mosaic_species_BOX(1)
      
      return
      end subroutine mosaic_box_aerchemistry











!***********************************************************************
! interface to dynamic gas-particle exchange solver
!
! author: Rahul A. Zaveri
! update: jan 2005
!-----------------------------------------------------------------------
      subroutine MOSAIC_dynamic_solver( dtchem )
        
      use module_data_mosaic_main, only: time_hrs
      use module_data_mosaic_aero
      use module_data_mosaic_asect

      implicit none

! subr arguments
      real(r8) :: dtchem
! local variables
      integer ibin, isize, itype, iv
      real(r8) :: XT
      real(r8) :: aerosol_water_up				! mosaic func

      do 510 ibin = 1, nbin_a

        
        call check_aerosol_mass(ibin)

        jaerosolstate_bgn(ibin) = jaerosolstate(ibin)
        if(jaerosolstate(ibin) .eq. no_aerosol)goto 500

        call conform_aerosol_number(ibin)   		! adjusts number conc so that it conforms with bin mass and diameter

        call conform_electrolytes(jtotal,ibin,XT) 	! conforms aer(jtotal) to a valid aerosol

        call check_aerosol_mass(ibin) 			! check mass again after conform_electrolytes
        jaerosolstate_bgn(ibin) = jaerosolstate(ibin)
        if(jaerosolstate(ibin) .eq. no_aerosol)goto 500	! ignore this bin

! *** should the "call conform_aerosol_number" go here instead of above ?

        ! when mhyst_method = mhyst_uporlo_waterhyst,
        ! initialize water_a_hyst at first time step using the user-input jhyst_leg
        if ( (mhyst_method == mhyst_uporlo_waterhyst) .and.   &
             (it_mosaic == 1) ) then
           if(jhyst_leg(ibin) == jhyst_lo)then
              water_a_hyst(ibin) = 0.0
           else
              water_a_up(ibin)   = aerosol_water_up(ibin)	! at 60% RH
              water_a_hyst(ibin) = water_a_up(ibin)
           endif
        end if

500     if (irepeat_mosaic == 1) then
           mass_dry_a_bgn(ibin) = mass_dry_a(ibin)
           if ( (jaerosolstate(ibin) .eq. no_aerosol) .or.   &
                (min(mass_dry_a(ibin),vol_dry_a(ibin)) .le. 1.0e-35) ) then
              call calc_aerosol_dry_density( ibin )
              dens_dry_a_bgn(ibin) = dens_dry_a(ibin)
           else
              dens_dry_a_bgn(ibin) = mass_dry_a(ibin)/vol_dry_a(ibin)
           end if
           dens_dry_a_bgn(ibin) = max( density_min_allow, &
                                  min( density_max_allow, dens_dry_a_bgn(ibin) ) )
        end if

        if (jaerosolstate(ibin) .eq. no_aerosol) then
           if (msize_framework == msectional) then
              isize = isize_of_ibin(ibin)
              itype = itype_of_ibin(ibin)
              Dp_dry_a(ibin) = dcen_sect(isize,itype)
              Dp_wet_a(ibin) = Dp_dry_a(ibin)
           end if
        end if

510   continue

      !PRINT*, 'HERE3-D'
      !PRINT*, Dp_dry_a
      !READ(*,*)

        
!cc        call save_pregrow_props	!3D

!cc	call specialoutaa( iclm_aer, jclm_aer, kclm_aer, 77, ! 3D
!cc     &		'after_conform' )
!
!-------------------------------------
! do dynamic gas-aerosol mass transfer for dtchem [s]
        
      if(mGAS_AER_XFER .eq. mON)then

         !        call wall_loss(dtchem)	! REMOVE WALL LOSS CALCULATION FROM HERE

        if(mDYNAMIC_SOLVER .eq. mASTEM)then
           call ASTEM(dtchem)
        elseif(mDYNAMIC_SOLVER .eq. mLSODE)then
          call MOSAIC_LSODE(dtchem)
        endif

      endif

      !PRINT*, 'HERE3-E'
      !PRINT*, Dp_dry_a
      !READ(*,*)
      
!-------------------------------------

! grows or shrinks size depending on mass increase or decrease

      !PRINT*, 'HERE-bCONFORM'
      !PRINT*, Dp_dry_a
      !READ(*,*)
      
      do ibin = 1, nbin_a
        if(jaerosolstate(ibin) .ne. no_aerosol)then
          call UPDATE_BINDIAM(ibin)	! BOX
        endif
      enddo

      !PRINT*, 'HERE-aCONFORM'
      !PRINT*, Dp_dry_a
      !READ(*,*)
      
      do 650 ibin = 1, nbin_a
        if(jaerosolstate(ibin).eq.no_aerosol) goto 600

        if (mhyst_method == mhyst_uporlo_jhyst) then
           if(jhyst_leg(ibin) == jhyst_lo)then
              water_a_hyst(ibin) = 0.0
           else
              water_a_up(ibin)   = aerosol_water_up(ibin)	! at 60% RH
              water_a_hyst(ibin) = water_a_up(ibin)
           endif
        elseif (mhyst_method == mhyst_uporlo_waterhyst) then
           water_a_up(ibin)   = aerosol_water_up(ibin)	! at 60% RH
           if (water_a_hyst(ibin) <= 0.5*water_a_up(ibin)) then
              jhyst_leg(ibin) = jhyst_lo
              water_a_hyst(ibin) = 0.0
           else
              jhyst_leg(ibin) = jhyst_up
              water_a_hyst(ibin) = water_a_up(ibin)
           endif
        else
           write(*,*) '*** MOSAIC_dynamic_solver - bad mhyst_method'
           stop
        endif

        call calc_dry_n_wet_aerosol_props(ibin)		! compute final mass and density

600     continue
        if ( (jaerosolstate(ibin) .eq. no_aerosol) .or.   &
             (min(mass_dry_a(ibin),vol_dry_a(ibin)) .le. 1.0e-35) ) then
           call calc_aerosol_dry_density( ibin )
        end if
        dens_dry_a(ibin) = max( density_min_allow, &
                           min( density_max_allow, dens_dry_a(ibin) ) )

650   continue

      if (irepeat_mosaic == 1) then

	  if(mod(it_mosaic,iprint_mosaic) .eq. 0)then
          call print_aer(1)	! UNCOMMENT THIS LINE
	  endif

      end if
      
      return
      end subroutine MOSAIC_dynamic_solver






!***********************************************************************
! initializes all the MOSAIC variables to zero or their default values.
!
! author: Rahul A. Zaveri
! update: jun 2003
!-----------------------------------------------------------------------
      subroutine initialize_mosaic_variables
      use module_data_mosaic_aero

      implicit none

! local variables
      integer iaer, ibin, iv, ja, jc, je



      do iv = 1, ngas_ioa
          gas(iv)        = 0.0
      enddo

! initialize to zero
      do ibin = 1, nbin_a

        num_a(ibin)          = 0.0
        mass_dry_a(ibin)     = 0.0
        mass_soluble_a(ibin) = 0.0
        dens_dry_a(ibin)     =-1.0

        do iaer = 1, naer
          aer(iaer,jtotal,ibin)  = 0.0
          aer(iaer,jsolid,ibin)  = 0.0
          aer(iaer,jliquid,ibin) = 0.0
        enddo

        do je = 1, nelectrolyte
          electrolyte(je,jtotal,ibin)  = 0.0
          electrolyte(je,jsolid,ibin)  = 0.0
          electrolyte(je,jliquid,ibin) = 0.0
          activity(je,ibin)            = 0.0
          gam(je,ibin)                 = 0.0
        enddo

          gam_ratio(ibin)   = 0.0

        do iv = 1, ngas_ioa
          flux_s(iv,ibin)   = 0.0
          flux_l(iv,ibin)   = 0.0
          kg(iv,ibin)       = 0.0
          phi_volatile_s(iv,ibin) = 0.0
          phi_volatile_l(iv,ibin) = 0.0
          df_gas(iv,ibin)   = 0.0
          volatile_s(iv,ibin) = 0.0
        enddo


        jaerosolstate(ibin) = -1	! initialize to default value
        jphase(ibin) = 0

        do jc = 1, ncation
          mc(jc,ibin) = 0.0
        enddo

        do ja = 1, nanion
          ma(ja,ibin) = 0.0
        enddo

      enddo	! ibin



      return
      end subroutine initialize_mosaic_variables






!***********************************************************************
! maps rsub(k,l,m) to and from MOSAIC arrays: gas and aer
!
! author: Rahul A. Zaveri
! update: nov 2001
!-----------------------------------------------------------------------
      subroutine map_mosaic_species_BOX(imap)
      use module_data_mosaic_main
      use module_data_mosaic_gas
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer imap
! local variables
      integer ibin, iaer, noffset
      real(r8) :: conv1, conv2, cnva1, cnva2




! define conversion factors
      conv1a = cair_mol_m3*1.e9		! converts q/mol(air) to nq/m^3 (q = mol or g)
      conv1b = 1./conv1a		! converts nq/m^3 to q/mol(air)
      conv2a = cair_mol_m3*18.*1.e-3	! converts mol(h2o)/mol(air) to kg(h2o)/m^3(air)
      conv2b = 1./conv2a		! converts kg(h2o)/m^3(air) to mol(h2o)/mol(air)


! BOX
      conv1 = 1.e15/avogad	! converts (molec/cc) to (nmol/m^3)
      conv2 = 1./conv1		! converts (nmol/m^3) to (molec/cc)

      cnva1 = 1.e3		! converts umol/m^3  to nmol/m^3 (or ug/m^3 to ng/m^3)
      cnva2 = 1./cnva1		! converts nmol/m^3) to umol/m^3 (or ng/m^3 to ug/m^3)


      if(imap.eq.0)then    ! map cnn (molec/cc) into stot (nmol/m^3)
! gas
        gas(ih2so4_g) = cnn(kh2so4)*conv1	! nmol/m^3
        gas(ihno3_g)  = cnn(khno3)*conv1
        gas(ihcl_g)   = cnn(khcl)*conv1
        gas(inh3_g)   = cnn(knh3)*conv1
        gas(imsa_g)   = cnn(kmsa)*conv1
        gas(icn3_g)  = cnn(kcn3)*conv1
        gas(icn2_g)  = cnn(kcn2)*conv1
        gas(icn1_g)  = cnn(kcn1)*conv1
        gas(ic0_g)  = cnn(kc0)*conv1
        gas(ic1_g)  = cnn(kc1)*conv1
        gas(ic2_g)  = cnn(kc2)*conv1
        gas(ic3_g)  = cnn(kc3)*conv1
        gas(ic4_g)  = cnn(kc4)*conv1
        gas(ic5_g)  = cnn(kc5)*conv1
        gas(ic6_g)  = cnn(kc6)*conv1
        gas(ic7_g)  = cnn(kc7)*conv1
        gas(ic8_g)  = cnn(kc8)*conv1
        gas(ic9_g)  = cnn(kc9)*conv1

! aerosol
        do ibin = 1, nbin_a
!!! MARK
          noffset = ngas_max + naer_tot*(ibin - 1)
          num_a(ibin)      = cnn(noffset + knum_a)! #/cc
          
          Dp_dry_a(ibin)   = cnn(noffset + kdpdry_a) * 1.e-4		! cm (dry diameter)
          
          Dpgn_a(ibin)     = Dp_dry_a(ibin)
          volume_bin(ibin) = piover6*(Dp_dry_a(ibin)**3)		! cm^3
          sigmag_a(ibin)   = cnn(noffset + ksigmag_a)
          jhyst_leg(ibin)  = INT(cnn(noffset + kjhyst_a)+0.01)
          water_a(ibin)    = cnn(noffset + kwater_a)			! kg/m^3

          do iaer = 1, naer
            aer(iaer,jtotal,ibin) = cnn(noffset+kwater_a+iaer)*cnva1	! nmol/m^3
          enddo

        enddo

      else                 ! map stot (nmol/m^3) back into cnn (umol/m^3)
! gas
        cnn(kh2so4)	= gas(ih2so4_g)*conv2
        cnn(khno3)	= gas(ihno3_g)*conv2
        cnn(khcl)		= gas(ihcl_g)*conv2
        cnn(knh3)		= gas(inh3_g)*conv2
        cnn(kmsa)		= gas(imsa_g)*conv2
        cnn(kcn3)	= gas(icn3_g)*conv2
        cnn(kcn2)	= gas(icn2_g)*conv2
        cnn(kcn1)	= gas(icn1_g)*conv2
        cnn(kc0)	= gas(ic0_g)*conv2
        cnn(kc1)	= gas(ic1_g)*conv2
        cnn(kc2)	= gas(ic2_g)*conv2
        cnn(kc3)	= gas(ic3_g)*conv2
        cnn(kc4)	= gas(ic4_g)*conv2
        cnn(kc5)	= gas(ic5_g)*conv2
        cnn(kc6)	= gas(ic6_g)*conv2
        cnn(kc7)	= gas(ic7_g)*conv2
        cnn(kc8)    = gas(ic8_g)*conv2
        cnn(kc9)    = gas(ic9_g)*conv2

! aerosol
        do ibin = 1, nbin_a

          noffset = ngas_max + naer_tot*(ibin - 1)
          cnn(noffset + knum_a)   = num_a(ibin)
          cnn(noffset + kdpdry_a) = Dp_dry_a(ibin) * 1.e4	! dry diameter (micron)
          cnn(noffset + kjhyst_a) = float(jhyst_leg(ibin))
          cnn(noffset + kwater_a) = water_a(ibin)		! kg/m^3

          do iaer = 1, naer
            cnn(noffset+kwater_a+iaer) = aer(iaer,jtotal,ibin)*cnva2	! molec/cc
          enddo

        enddo

      endif

      return
      end subroutine map_mosaic_species_BOX








      subroutine overall_massbal_in
      use module_data_mosaic_aero

      implicit none

      integer ibin

      tot_so4_in = gas(ih2so4_g)
      tot_no3_in = gas(ihno3_g)
      tot_cl_in  = gas(ihcl_g)
      tot_nh4_in = gas(inh3_g)
      tot_na_in  = 0.0
      tot_ca_in  = 0.0


      do ibin = 1, nbin_a
        tot_so4_in = tot_so4_in + aer(iso4_a,jtotal,ibin)
	tot_no3_in = tot_no3_in + aer(ino3_a,jtotal,ibin)
        tot_cl_in  = tot_cl_in  + aer(icl_a, jtotal,ibin)
        tot_nh4_in = tot_nh4_in + aer(inh4_a,jtotal,ibin)
        tot_na_in  = tot_na_in  + aer(ina_a,jtotal,ibin)
        tot_ca_in  = tot_ca_in  + aer(ica_a,jtotal,ibin)
      enddo


        total_species(inh3_g) = tot_nh4_in
        total_species(ihno3_g)= tot_no3_in
        total_species(ihcl_g) = tot_cl_in


      return
      end subroutine overall_massbal_in



      subroutine overall_massbal_out(mbin)
!      include 'v33com'
!      include 'v33com3'
!      include 'v33com9a'
!      include 'v33com9b'
      use module_data_mosaic_aero

      implicit none

! subr. agrument
      integer mbin
! local variables
      integer ibin



      tot_so4_out = gas(ih2so4_g)
      tot_no3_out = gas(ihno3_g)
      tot_cl_out  = gas(ihcl_g)
      tot_nh4_out = gas(inh3_g)
      tot_na_out  = 0.0
      tot_ca_out  = 0.0

      do ibin = 1, nbin_a
        tot_so4_out = tot_so4_out + aer(iso4_a,jtotal,ibin)
        tot_no3_out = tot_no3_out + aer(ino3_a,jtotal,ibin)
        tot_cl_out  = tot_cl_out  + aer(icl_a,jtotal,ibin)
        tot_nh4_out = tot_nh4_out + aer(inh4_a,jtotal,ibin)
        tot_na_out  = tot_na_out  + aer(ina_a,jtotal,ibin)
        tot_ca_out  = tot_ca_out  + aer(ica_a,jtotal,ibin)
      enddo

      diff_so4 = tot_so4_out - tot_so4_in
      diff_no3 = tot_no3_out - tot_no3_in
      diff_cl  = tot_cl_out  - tot_cl_in
      diff_nh4 = tot_nh4_out - tot_nh4_in
      diff_na  = tot_na_out  - tot_na_in
      diff_ca  = tot_ca_out  - tot_ca_in


      reldiff_so4 = 0.0
      if(tot_so4_in .gt. 1.e-25 .or. tot_so4_out .gt. 1.e-25)then
	reldiff_so4 = diff_so4/max(tot_so4_in, tot_so4_out)
      endif

      reldiff_no3 = 0.0
      if(tot_no3_in .gt. 1.e-25 .or. tot_no3_out .gt. 1.e-25)then
	reldiff_no3 = diff_no3/max(tot_no3_in, tot_no3_out)
      endif

      reldiff_cl = 0.0
      if(tot_cl_in .gt. 1.e-25 .or. tot_cl_out .gt. 1.e-25)then
	reldiff_cl = diff_cl/max(tot_cl_in, tot_cl_out)
      endif

      reldiff_nh4 = 0.0
      if(tot_nh4_in .gt. 1.e-25 .or. tot_nh4_out .gt. 1.e-25)then
	reldiff_nh4 = diff_nh4/max(tot_nh4_in, tot_nh4_out)
      endif

      reldiff_na = 0.0
      if(tot_na_in .gt. 1.e-25 .or. tot_na_out .gt. 1.e-25)then
	reldiff_na = diff_na/max(tot_na_in, tot_na_out)
      endif

      reldiff_ca = 0.0
      if(tot_ca_in .gt. 1.e-25 .or. tot_ca_out .gt. 1.e-25)then
	reldiff_ca = diff_ca/max(tot_ca_in, tot_ca_out)
      endif



      if( abs(reldiff_so4) .gt. 1.e-4 .or.   &
          abs(reldiff_no3) .gt. 1.e-4 .or.   &
          abs(reldiff_cl)  .gt. 1.e-4 .or.   &
          abs(reldiff_nh4) .gt. 1.e-4 .or.   &
          abs(reldiff_na)  .gt. 1.e-4 .or.   &
          abs(reldiff_ca)  .gt. 1.e-4)then


        if(iprint_input .eq. mYES)then
          write(6,*)'*** mbin = ', mbin, '  isteps = ', isteps_ASTEM
          write(6,*)'reldiff_so4 = ', reldiff_so4
          write(6,*)'reldiff_no3 = ', reldiff_no3
          write(6,*)'reldiff_cl  = ', reldiff_cl
          write(6,*)'reldiff_nh4 = ', reldiff_nh4
          write(6,*)'reldiff_na  = ', reldiff_na
          write(6,*)'reldiff_ca  = ', reldiff_ca
!          call print_input
          iprint_input = mNO
        endif

!         stop

      endif


      return
      end subroutine overall_massbal_out

















!***********************************************************************
! checks if aerosol mass is too low to be of any significance
! and determine jaerosolstate
!
! author: Rahul A. Zaveri
! update: jan 2005
!-----------------------------------------------------------------------
      subroutine check_aerosol_mass(ibin)
      USE module_data_mosaic_main, ONLY: nCOMP
      use module_data_mosaic_aero

      implicit none

      INTEGER :: i,j,k
      
      REAL(R8),PARAMETER :: NA = 6.022e23      
      
! subr arguments
      integer ibin
! local variables
      integer iaer
      real(r8) :: drymass, aer_H

      mass_dry_a(ibin) = 0.0

      aer_H = (2.*aer(iso4_a,jtotal,ibin) +   &
                  aer(ino3_a,jtotal,ibin) +   &
                  aer(icl_a,jtotal,ibin)  +   &
                  aer(imsa_a,jtotal,ibin) +   &
               2.*aer(ico3_a,jtotal,ibin))-   &
              (2.*aer(ica_a,jtotal,ibin)  +   &
                  aer(ina_a,jtotal,ibin)  +   &
                  aer(inh4_a,jtotal,ibin))
      aer_H = max(aer_H, 0.0d0)

      
      
      !PRINT*, 'aer_H = ',aer_H

      do iaer = 1, naer
        mass_dry_a(ibin) = mass_dry_a(ibin) +   &
                           aer(iaer,jtotal,ibin)*mw_aer_mac(iaer)	! ng/m^3(air)
      enddo

      
      DO k = 1,nBINS
         DO j = 1,nCOMP
            mass_dry_a(ibin) = mass_dry_a(ibin) +   &
            jk_OLGMOLE_SOM(j,k)/(1e-9*NA*1e-6)*MW_AER_MAC(icn3_a+j-1)
         END DO
      END DO

     
      mass_dry_a(ibin) = mass_dry_a(ibin) + aer_H

      drymass = mass_dry_a(ibin)			! ng/m^3(air)
      mass_dry_a(ibin) = mass_dry_a(ibin)*1.e-15	! g/cc(air)

      if(drymass .le. mass_cutoff)then			! bin mass is too small
        jaerosolstate(ibin) = no_aerosol
        jphase(ibin) = 0
        if(drymass .eq. 0.)num_a(ibin) = 0.0
      endif

      return
      end subroutine check_aerosol_mass





!***********************************************************************
! checks and conforms number according to the mass and bin size range
!
! author: Rahul A. Zaveri
! update: jan 2005
!-----------------------------------------------------------------------
      
      subroutine conform_aerosol_number(ibin)	! TOUCH
!      include 'v33com9a'
      use module_data_mosaic_aero
      use module_data_mosaic_asect

      implicit none

! subr arguments
      integer ibin
! local variables
      integer :: iaer, isize, itype
      real(r8) :: num_at_dlo, num_at_dhi, numold
      real(r8) :: aer_H


! when msize_framework = munstructured or mmodal,
!    calculate number from volume concentration and mean dry diameter
!       only when num_a(ibin) <= 0.0
!    this should only happen at the very start of the simulation
! when msize_framework = msectional,
!    check that mean dry diameter falls within the section/bin limits,
!    and adjust number is this is not true
      if (msize_framework /= msectional) then
         if (num_a(ibin) > 0.0) return
      end if

      vol_dry_a(ibin)  = 0.0		! initialize to 0.0

      if(jaerosolstate(ibin) .eq. no_aerosol) return

! calculate dry volume concentration
      aer_H = (2.*aer(iso4_a,jtotal,ibin) +   &
                  aer(ino3_a,jtotal,ibin) +   &
                  aer(icl_a,jtotal,ibin)  +   &
                  aer(imsa_a,jtotal,ibin) +   &
               2.*aer(ico3_a,jtotal,ibin))-   &
              (2.*aer(ica_a,jtotal,ibin)  +   &
                  aer(ina_a,jtotal,ibin)  +   &
                  aer(inh4_a,jtotal,ibin))

      aer_H = max(aer_H, 0.0d0)


      
      do iaer = 1, naer
         
        vol_dry_a(ibin) = vol_dry_a(ibin) +   &
             aer(iaer,jtotal,ibin)*mw_aer_mac(iaer)/dens_aer_mac(iaer)  ! ncc/m^3(air)
        
      enddo
     
     

     
      vol_dry_a(ibin) = vol_dry_a(ibin) + aer_H !!!!!!!!!!!!!!!!!!!!!!!!!!!!! MARKED

      vol_dry_a(ibin) = vol_dry_a(ibin)*1.e-15			   ! cc(aer)/cc(air)
      

      !PRINT*, Dp_dry_a(1:3)*1e7
      !PRINT*, dlo_sect(1:3,1)*1e7
      !PRINT*, dhi_sect(1:3,1)*1e7
      !STOP
      
      
      if (msize_framework /= msectional) then
! unstructured or modal - set (initialize) number
         num_a(ibin) = vol_dry_a(ibin)/volume_bin(ibin)	! #/cc(air)

      else
! sectional
         if (num_a(ibin) <= 0.0) then
! in this case, num_a has probably not yet been initialized, so do it
            num_a(ibin) = vol_dry_a(ibin)/volume_bin(ibin)	! #/cc(air)
         else
! in this case, check that bin mean size is within bounds
            isize = isize_of_ibin( ibin )
            itype = itype_of_ibin( ibin )
            num_at_dlo = vol_dry_a(ibin)/volumlo_sect(isize,itype)
            num_at_dhi = vol_dry_a(ibin)/volumhi_sect(isize,itype)
            numold = num_a(ibin)
            num_a(ibin) = min( num_a(ibin), num_at_dlo )
            num_a(ibin) = max( num_a(ibin), num_at_dhi )
         end if
      end if


      return
      end subroutine conform_aerosol_number





!***********************************************************************
! calculates dry density
!
! author: Rahul A. Zaveri
! update: apr 2010
!-----------------------------------------------------------------------
      subroutine calc_aerosol_dry_density( ibin )
!      include 'v33com9a'
      use module_data_mosaic_aero
      use module_data_mosaic_asect

      implicit none

! subr arguments
      integer ibin
! local variables
      integer :: iaer
      real(r8) :: aer_H
      real(r8) :: tmpa, tmp_volu, tmp_mass


! calculate dry volume concentration
      aer_H = ( 2.*max( 0.0_r8, aer(iso4_a,jtotal,ibin) ) +   &
                   max( 0.0_r8, aer(ino3_a,jtotal,ibin) ) +   &
                   max( 0.0_r8, aer(icl_a,jtotal,ibin) )  +   &
                   max( 0.0_r8, aer(imsa_a,jtotal,ibin) ) +   &
                2.*max( 0.0_r8, aer(ico3_a,jtotal,ibin) ) )   &
            - ( 2.*max( 0.0_r8, aer(ica_a,jtotal,ibin) )  +   &
                   max( 0.0_r8, aer(ina_a,jtotal,ibin) )  +   &
                   max( 0.0_r8, aer(inh4_a,jtotal,ibin) ) )
      aer_H = max( aer_H, 0.0_r8 )

      tmp_mass = aer_H
      tmp_volu = aer_H   ! assume density=1.0 for H+

      do iaer = 1, naer
        tmpa = max( 0.0_r8, aer(iaer,jtotal,ibin) ) * mw_aer_mac(iaer)
        tmp_mass = tmp_mass + tmpa                     !  ng/m^3(air)
        tmp_volu = tmp_volu + tmpa/dens_aer_mac(iaer)  ! ncc/m^3(air)
      enddo

! the 1.0e-20 ng/m3 cutoff here is equivalent to the
!     1.0e-35 g/cm3 cutoff used in mosaic_dynamic_solver
      if (min(tmp_mass,tmp_volu) >= 1.0e-20) then
         dens_dry_a(ibin) = tmp_mass/tmp_volu   ! g/cc
      else
         dens_dry_a(ibin) = 1.0
      end if

      return
      end subroutine calc_aerosol_dry_density











!***********************************************************************
! updates/conforms size (diameter) according to the mass and number
!
! author: Rahul A. Zaveri
! update: oct 2005
!-----------------------------------------------------------------------
      SUBROUTINE UPDATE_BINDIAM(iBIN)
        
      USE module_data_mosaic_main, ONLY:  piover6, third,nCOMP
      USE module_data_mosaic_aero

      IMPLICIT NONE
      
      INTEGER :: i,j,k
      
      REAL(R8),PARAMETER :: NA = 6.022e23      

! subr arguments
      integer ibin
! local variables
      integer iaer
      
      real(r8) :: aer_H


      vol_dry_a(ibin)  = 0.0		! initialize to 0.0

      if(jaerosolstate(ibin) .eq. no_aerosol) return

      aer_H = (2.*aer(iso4_a,jtotal,ibin) +   &
                  aer(ino3_a,jtotal,ibin) +   &
                  aer(icl_a,jtotal,ibin)  +   &
                  aer(imsa_a,jtotal,ibin) +   &
               2.*aer(ico3_a,jtotal,ibin))-   &
              (2.*aer(ica_a,jtotal,ibin)  +   &
                  aer(ina_a,jtotal,ibin)  +   &
                  aer(inh4_a,jtotal,ibin))
      
      aer_H = max(aer_H, 0.0d0)

      
      
      do iaer = 1, naer
        vol_dry_a(ibin) = vol_dry_a(ibin) +   &
        aer(iaer,jtotal,ibin)*mw_aer_mac(iaer)/dens_aer_mac(iaer)	! ng/m^3(air)
      enddo
      
      vol_dry_a(ibin) = vol_dry_a(ibin) + aer_H



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CHARLES HE
      !OLGMASS = 0.d0

      DO k = 1,nBINS
         DO j = 1,nCOMP
            vol_dry_a(ibin) = vol_dry_a(ibin) +   &
            jk_OLGMOLE_SOM(j,k)/(1e-9*NA*1e-6)*MW_AER_MAC(icn3_a+j-1)/dens_aer_mac(icn3_a+j-1)
         END DO
      END DO
      
      
      
      vol_dry_a(ibin) = vol_dry_a(ibin)*1.e-15				! cc(aer)/cc(air)


! update size
!
! Box-model only
      Dp_dry_a(ibin) = (vol_dry_a(ibin)/(piover6*num_a(ibin)))**third
      
      Dpgn_a(ibin) = Dp_dry_a(ibin)

      RETURN
      END SUBROUTINE UPDATE_BINDIAM

















!***********************************************************************
! determines phase state of an aerosol bin. includes kelvin effect.
!
! author: Rahul A. Zaveri
! update: jan 2005
!-----------------------------------------------------------------------
      subroutine aerosol_phase_state(ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      integer js, je, iaer, iv, iter_kelvin
      real(r8) :: aH2O_a_new, rel_err
      real(r8) :: aerosol_water_up, bin_molality		! mosaic func
      real(r8) :: kelvin_toler, term
      real(r8) :: aer_H


      aH2O = RH_pc*0.01
      aH2O_a(ibin) = aH2O
      kelvin(ibin) = 1.0
      do iv = 1, ngas_volatile
        kel(iv,ibin) = 1.0
      enddo

      if(RH_pc .le. 99)then
        kelvin_toler = 1.e-4
      else
        kelvin_toler = 1.e-8
      endif

! calculate dry mass and dry volume of a bin
      mass_dry_a(ibin) = 0.0		! initialize to 0.0
      vol_dry_a(ibin)  = 0.0		! initialize to 0.0

      aer_H = (2.*aer(iso4_a,jtotal,ibin) +   &
                  aer(ino3_a,jtotal,ibin) +   &
                  aer(icl_a,jtotal,ibin)  +   &
                  aer(imsa_a,jtotal,ibin) +   &
               2.*aer(ico3_a,jtotal,ibin))-   &
              (2.*aer(ica_a,jtotal,ibin)  +   &
                  aer(ina_a,jtotal,ibin)  +   &
                  aer(inh4_a,jtotal,ibin))
      aer_H = max(aer_H, 0.0d0)

      do iaer = 1, naer
        mass_dry_a(ibin) = mass_dry_a(ibin) +   &
                           aer(iaer,jtotal,ibin)*mw_aer_mac(iaer)	! ng/m^3(air)
        vol_dry_a(ibin)  = vol_dry_a(ibin) +   &
        aer(iaer,jtotal,ibin)*mw_aer_mac(iaer)/dens_aer_mac(iaer)  	! ncc/m^3(air)
      enddo
      mass_dry_a(ibin) = mass_dry_a(ibin) + aer_H
      vol_dry_a(ibin) = vol_dry_a(ibin) + aer_H

      mass_dry_a(ibin) = mass_dry_a(ibin)*1.e-15			! g/cc(air)
      vol_dry_a(ibin)  = vol_dry_a(ibin)*1.e-15				! cc(aer)/cc(air) or m^3/m^3(air)

! wet mass and wet volume
      mass_wet_a(ibin) = mass_dry_a(ibin) + water_a(ibin)*1.e-3		! g/cc(air)
      vol_wet_a(ibin)  = vol_dry_a(ibin) + water_a(ibin)*1.e-3		! cc(aer)/cc(air) or m^3/m^3(air)


      water_a_up(ibin) = aerosol_water_up(ibin)	! for hysteresis curve determination

      iter_kelvin = 0

10    iter_kelvin = iter_kelvin + 1
!!      do je = 1, nelectrolyte
!!        molality0(je) = bin_molality(je,ibin)	! compute aH2O dependent binary molalities  EFFI
!!      enddo

      call MESA(ibin)
      if(jaerosolstate(ibin) .eq. all_solid)then
        return
      endif
! new wet mass and wet volume
      mass_wet_a(ibin) = mass_dry_a(ibin) + water_a(ibin)*1.e-3		! g/cc(air)
      vol_wet_a(ibin)  = vol_dry_a(ibin) + water_a(ibin)*1.e-3		! cc(aer)/cc(air) or m^3/m^3(air)

      call calculate_kelvin(ibin)
!      kelvin(ibin) = 1.0

      aH2O_a_new = RH_pc*0.01/kelvin(ibin)

      rel_err = abs( (aH2O_a_new - aH2O_a(ibin))/aH2O_a(ibin))

      if(rel_err .gt. kelvin_toler)then
        aH2O_a(ibin) = aH2O_a_new
!        aH2O = aH2O_a_new
!        call MTEM_compute_log_gamZ	! recompute activity coeffs (for surface tension and solid-liquid equilibria)
        goto 10
      endif


      if(jaerosolstate(ibin) .eq. all_liquid)jhyst_leg(ibin) = jhyst_up

! now compute kelvin effect terms for condensing species (nh3, hno3, and hcl)
      do iv = 1,  ngas_volatile
        term = 4.*sigma_soln(ibin)*partial_molar_vol(iv)/   &
                       (8.3144e7*T_K*DpmV(ibin))
        kel(iv,ibin) = 1. + term*(1. + 0.5*term*(1. + term/3.))
      enddo

! edit wkc: the kelvin effect for organics should be calculated based on
! particle size only, and not by the terms above.
! Instead of changing that section, we'll try overwriting it
    do ibin = 1, nbin_a
      do iv = icn3_g, ngas_volatile
        if(Dp_dry_a(ibin) .eq. 0) then
            kel(iv,ibin) = 1
        else
            kel(iv,ibin) = 10.0**(3.75d-7 / Dp_dry_a(ibin))
        endif
      enddo
    enddo
      return
      end subroutine aerosol_phase_state






!***********************************************************************
! computes kelvin effect term (kelvin => 1.0)
!
! author: Rahul A. Zaveri
! update: jan 2005
!-----------------------------------------------------------------------
      subroutine calculate_kelvin(ibin)
      use module_data_mosaic_main, only:  pi
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      integer je
      real(r8) :: term, sum_dum


      volume_a(ibin) = vol_wet_a(ibin) 				! [cc/cc(air)]
      DpmV(ibin)=(6.*volume_a(ibin)/(num_a(ibin)*pi))**(1./3.)	! [cm]


! Li and Lu (2001) surface tension model:
!      sum_dum = 0.0
!      do je = 1, nelectrolyte
!        sum_dum = sum_dum + G_MX(je)*
!     &                      alog(1./(1.+K_MX(je)*activity(je,ibin)))
!      enddo
!      sigma_soln(ibin) = sigma_water + 8.3144e7*T_K*sum_dum


! simpler correlation for solution surface tension:
      sigma_soln(ibin) = sigma_water + 49.0*(1. - aH2O_a(ibin)) 	! [dyn/cm]



      term = 4.*18.*sigma_soln(ibin)/(8.3144e7*T_K*DpmV(ibin))		! [-]
!      kelvin(ibin) = exp(term)
      kelvin(ibin) = 1. + term*(1. + 0.5*term*(1. + term/3.))


      return
      end subroutine calculate_kelvin





















!***********************************************************************
! MESA: Multicomponent Equilibrium Solver for Aerosols.
! Computes equilibrum solid and liquid phases by integrating
! pseudo-transient dissolution and precipitation reactions
!
! author: Rahul A. Zaveri
! update: jan 2005
!-----------------------------------------------------------------------
      subroutine MESA(ibin)	! TOUCH
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin

! local variables
      integer idissolved, j_index, jdum, js, je
      real(r8) :: CRH, solids, sum_soluble, sum_insoluble, XT
      real(r8) :: aerosol_water				! mosaic func
      real(r8) :: drh_mutual				! mosaic func
      real(r8) :: H_ion, sum_dum


!! EFFI
!! calculate percent composition
      sum_dum = 0.0
      do je = 1, nelectrolyte
        sum_dum = sum_dum + electrolyte(je,jtotal,ibin)
      enddo

      if(sum_dum .eq. 0.)sum_dum = 1.0

      do je = 1, nelectrolyte
        epercent(je,jtotal,ibin) = 100.*electrolyte(je,jtotal,ibin)/sum_dum
      enddo


      call calculate_XT(ibin,jtotal,XT)

      CRH = 0.35

! step 1: check if aH2O is below CRH (crystallization or efflorescence point)
      if( (aH2O_a(ibin).lt.CRH)     .and. &
          (XT.gt.1.0 .or. XT.lt.0.) .and. &
          (epercent(jcano3,jtotal,ibin) .le. ptol_mol_astem) .and. &
          (epercent(jcacl2,jtotal,ibin) .le. ptol_mol_astem) )then
        jaerosolstate(ibin) = all_solid
        jphase(ibin)    = jsolid
        jhyst_leg(ibin) = jhyst_lo
        call adjust_solid_aerosol(ibin)
        return
      endif


! step 2: check for supersaturation/metastable state
      jdum = 0
      if (mhyst_method == mhyst_uporlo_jhyst) then         ! BOX method/logic
         if (jhyst_leg(ibin) == jhyst_up) jdum = 1
      elseif (mhyst_method == mhyst_uporlo_waterhyst) then ! 3-D method/logic
         if (water_a_hyst(ibin) > 0.5*water_a_up(ibin)) jdum = 1
      else
         write(*,*) '*** MESA - bad mhyst_method'
         stop
      endif
      if (jdum == 1) then
        call do_full_deliquescence(ibin)

!        call ions_to_electrolytes(jliquid,ibin,XT) ! for Li and Lu surface tension
!        call compute_activities(ibin)		   ! for Li and Lu surface tension

        sum_soluble = 0.0
        do js = 1, nsoluble
          sum_soluble = sum_soluble + electrolyte(js,jtotal,ibin)
        enddo

        solids = electrolyte(jcaso4,jtotal,ibin) +   &
                 electrolyte(jcaco3,jtotal,ibin) +   &
                 aer(ioin_a ,jtotal,ibin)


        if(sum_soluble .lt. 1.e-15 .and. solids .gt. 0.0)then

          jaerosolstate(ibin) = all_solid ! no soluble material present
          jphase(ibin) = jsolid
          call adjust_solid_aerosol(ibin)

! new wet mass and wet volume
          mass_wet_a(ibin) = mass_dry_a(ibin) + water_a(ibin)*1.e-3	! g/cc(air)
          vol_wet_a(ibin)  = vol_dry_a(ibin) + water_a(ibin)*1.e-3	! cc(aer)/cc(air) or m^3/m^3(air)
          growth_factor(ibin) = mass_wet_a(ibin)/mass_dry_a(ibin)	! mass growth factor

          return

        elseif(sum_soluble .gt. 0.0 .and. solids .eq. 0.0)then

          jaerosolstate(ibin) = all_liquid
          jhyst_leg(ibin) = jhyst_up
          jphase(ibin) = jliquid
          water_a(ibin) = aerosol_water(jtotal,ibin)

          if(water_a(ibin) .lt. 0.0)then
            jaerosolstate(ibin) = all_solid ! no soluble material present
            jphase(ibin)    = jsolid
            jhyst_leg(ibin) = jhyst_lo
            call adjust_solid_aerosol(ibin)
          else
            call adjust_liquid_aerosol(ibin)
            call compute_activities(ibin)
          endif

! new wet mass and wet volume
          mass_wet_a(ibin) = mass_dry_a(ibin) + water_a(ibin)*1.e-3	! g/cc(air)
          vol_wet_a(ibin)  = vol_dry_a(ibin) + water_a(ibin)*1.e-3	! cc(aer)/cc(air) or m^3/m^3(air)
          growth_factor(ibin) = mass_wet_a(ibin)/mass_dry_a(ibin)	! mass growth factor

          return

        endif

      endif




! step 3: diagnose MDRH
      if(XT .lt. 1. .and. XT .gt. 0. )goto 10	! excess sulfate domain - no MDRH exists

      jdum = 0
      do js = 1, nsalt
        jsalt_present(js) = 0			! default value - salt absent

        if(epercent(js,jtotal,ibin) .gt. ptol_mol_astem)then
          jsalt_present(js) = 1			! salt present
          jdum = jdum + jsalt_index(js)
        endif
      enddo

      if(jdum .eq. 0)then
        jaerosolstate(ibin) = all_solid ! no significant soluble material present
        jphase(ibin) = jsolid
        call adjust_solid_aerosol(ibin)
        return
      endif

      if(XT .ge. 2.0 .or. XT .lt. 0.0)then
        j_index = jsulf_poor(jdum)
      else
        j_index = jsulf_rich(jdum)
      endif

      MDRH(ibin) = MDRH_T(j_index)

      if(aH2O_a(ibin)*100. .lt. MDRH(ibin)) then
        jaerosolstate(ibin) = all_solid
        jphase(ibin) = jsolid
        jhyst_leg(ibin) = jhyst_lo
        call adjust_solid_aerosol(ibin)
        return
      endif


! step 4: none of the above means it must be sub-saturated or mixed-phase
10    call do_full_deliquescence(ibin)
      call MESA_PTC(ibin)	! determines jaerosolstate(ibin)



      return
      end subroutine MESA








!***********************************************************************
! this subroutine completely deliquesces an aerosol and partitions
! all the soluble electrolytes into the liquid phase and insoluble
! ones into the solid phase. It also calculates the corresponding
! aer(js,jliquid,ibin) and aer(js,jsolid,ibin) generic species
! concentrations
!
! author: Rahul A. Zaveri
! update: jan 2005
!-----------------------------------------------------------------------
      subroutine do_full_deliquescence(ibin)	! TOUCH
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      integer js




! partition all electrolytes into liquid phase
      do js = 1, nelectrolyte
       electrolyte(js,jsolid,ibin)  = 0.0
       electrolyte(js,jliquid,ibin) = electrolyte(js,jtotal,ibin)
      enddo
!
! except these electrolytes, which always remain in the solid phase
      electrolyte(jcaco3,jsolid,ibin) = electrolyte(jcaco3,jtotal,ibin)
      electrolyte(jcaso4,jsolid,ibin) = electrolyte(jcaso4,jtotal,ibin)
      electrolyte(jcaco3,jliquid,ibin)= 0.0
      electrolyte(jcaso4,jliquid,ibin)= 0.0


! partition all the generic aer species into solid and liquid phases
! solid phase
      aer(iso4_a,jsolid,ibin) = electrolyte(jcaso4,jsolid,ibin)
      aer(ino3_a,jsolid,ibin) = 0.0
      aer(icl_a, jsolid,ibin) = 0.0
      aer(inh4_a,jsolid,ibin) = 0.0
      aer(ioc_a, jsolid,ibin) = aer(ioc_a,jtotal,ibin)
      aer(imsa_a,jsolid,ibin) = 0.0
      aer(ico3_a,jsolid,ibin) = aer(ico3_a,jtotal,ibin)
      aer(ina_a, jsolid,ibin) = 0.0
      aer(ica_a, jsolid,ibin) = electrolyte(jcaco3,jsolid,ibin) +   &
                                electrolyte(jcaso4,jsolid,ibin)
      aer(ibc_a, jsolid,ibin) = aer(ibc_a,jtotal,ibin)
      aer(ioin_a,jsolid,ibin) = aer(ioin_a,jtotal,ibin)
      aer(icn3_a,jsolid,ibin)= aer(icn3_a,jtotal,ibin)
      aer(icn2_a,jsolid,ibin)= aer(icn2_a,jtotal,ibin)
      aer(icn1_a,jsolid,ibin)= aer(icn1_a,jtotal,ibin)
      aer(ic0_a,jsolid,ibin)= aer(ic0_a,jtotal,ibin)
      aer(ic1_a,jsolid,ibin)= aer(ic1_a,jtotal,ibin)
      aer(ic2_a,jsolid,ibin)= aer(ic2_a,jtotal,ibin)
      aer(ic3_a,jsolid,ibin)= aer(ic3_a,jtotal,ibin)
      aer(ic4_a,jsolid,ibin)= aer(ic4_a,jtotal,ibin)
      aer(ic5_a,jsolid,ibin)= aer(ic5_a,jtotal,ibin)
      aer(ic6_a,jsolid,ibin)= aer(ic6_a,jtotal,ibin)
	  aer(ic7_a,jsolid,ibin)= aer(ic7_a,jtotal,ibin)


! liquid-phase
      aer(iso4_a,jliquid,ibin) = aer(iso4_a,jtotal,ibin) -   &
                                 electrolyte(jcaso4,jsolid,ibin)
      aer(ino3_a,jliquid,ibin) = aer(ino3_a,jtotal,ibin)
      aer(icl_a, jliquid,ibin) = aer(icl_a,jtotal,ibin)
      aer(inh4_a,jliquid,ibin) = aer(inh4_a,jtotal,ibin)
      aer(ioc_a, jliquid,ibin) = 0.0
      aer(imsa_a,jliquid,ibin) = aer(imsa_a,jtotal,ibin)
      aer(ico3_a,jliquid,ibin) = 0.0
      aer(ina_a, jliquid,ibin) = aer(ina_a,jtotal,ibin)
      aer(ica_a, jliquid,ibin) = electrolyte(jcano3,jtotal,ibin) +   &
                                 electrolyte(jcacl2,jtotal,ibin)
      aer(ibc_a, jliquid,ibin) = 0.0
      aer(ioin_a,jliquid,ibin) = 0.0
      aer(icn3_a,jliquid,ibin)= 0.0
      aer(icn2_a,jliquid,ibin)= 0.0
      aer(ic1_a,jliquid,ibin)= 0.0
      aer(ic2_a,jliquid,ibin)= 0.0
      aer(ic3_a,jliquid,ibin)= 0.0
      aer(ic4_a,jliquid,ibin)= 0.0
      aer(ic5_a,jliquid,ibin)= 0.0
      aer(ic6_a,jliquid,ibin)= 0.0
      aer(ic7_a,jliquid,ibin)= 0.0

      return
      end subroutine do_full_deliquescence






















!***********************************************************************
! MESA: Multicomponent Equilibrium Solver for Aerosol-phase
! computes equilibrum solid and liquid phases by integrating
! pseudo-transient dissolution and precipitation reactions
!
! author: Rahul A. Zaveri
! update: jan 2005
! Reference: Zaveri R.A., R.C. Easter, and L.K. Peters, JGR, 2005b
!-----------------------------------------------------------------------
      subroutine MESA_PTC(ibin)		! TOUCH
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      integer iaer, iconverge, iconverge_flux, iconverge_mass,   &
           idissolved, itdum, js, je, jp

      real(r8) :: tau_p(nsalt), tau_d(nsalt)
      real(r8) :: frac_solid, sumflux, hsalt_min, alpha, XT, dumdum,   &
           H_ion
      real(r8) :: phi_prod, alpha_fac, sum_dum
      real(r8) :: aer_H
! function
      real(r8) :: aerosol_water



! initialize
      itdum = 0		! initialize time
      hsalt_max = 1.e25



      do js = 1, nsalt
        hsalt(js)     = 0.0
        sat_ratio(js) = 0.0
        phi_salt(js)  = 0.0
        flux_sl(js)   = 0.0
      enddo



!! EFFI calculate percent composition
      sum_dum = 0.0
      do je = 1, nelectrolyte
        sum_dum = sum_dum + electrolyte(je,jtotal,ibin)
      enddo

      if(sum_dum .eq. 0.)sum_dum = 1.0

      do je = 1, nelectrolyte
        epercent(je,jtotal,ibin) = 100.*electrolyte(je,jtotal,ibin)/sum_dum
      enddo
!! EFFI



      do js = 1, nsalt
        jsalt_present(js) = 0			! default value - salt absent
        if(epercent(js,jtotal,ibin) .gt. 1.0)then
          jsalt_present(js) = 1			! salt present
        endif
      enddo


      mass_dry_a(ibin) = 0.0

      aer_H = (2.*aer(iso4_a,jtotal,ibin) +   &
                  aer(ino3_a,jtotal,ibin) +   &
                  aer(icl_a,jtotal,ibin)  +   &
                  aer(imsa_a,jtotal,ibin) +   &
               2.*aer(ico3_a,jtotal,ibin))-   &
              (2.*aer(ica_a,jtotal,ibin)  +   &
                  aer(ina_a,jtotal,ibin)  +   &
                  aer(inh4_a,jtotal,ibin))
      aer_H = max(aer_H, 0.0d0)

      do iaer = 1, naer
       mass_dry_a(ibin) = mass_dry_a(ibin) +   &
          aer(iaer,jtotal,ibin)*mw_aer_mac(iaer) 	! [ng/m^3(air)]
        vol_dry_a(ibin)  = vol_dry_a(ibin) +   &
        aer(iaer,jtotal,ibin)*mw_aer_mac(iaer)/dens_aer_mac(iaer)  	! ncc/m^3(air)
      enddo
      mass_dry_a(ibin) = mass_dry_a(ibin) + aer_H
      vol_dry_a(ibin) = vol_dry_a(ibin) + aer_H

      mass_dry_a(ibin) = mass_dry_a(ibin)*1.e-15			! [g/cc(air)]
      vol_dry_a(ibin) = vol_dry_a(ibin)*1.e-15				! [cc(aer)/cc(air)]

      mass_dry_salt(ibin) = 0.0		! soluble salts only
      do je = 1, nsalt
        mass_dry_salt(ibin) = mass_dry_salt(ibin) +   &
              electrolyte(je,jtotal,ibin)*mw_electrolyte(je)*1.e-15	! g/cc(air)
      enddo

      jMESA_call = jMESA_call + 1

!----begin pseudo time continuation loop-------------------------------

      do 500 itdum = 1, Nmax_MESA


! compute new salt fluxes
      call MESA_flux_salt(ibin)


! check convergence
      call MESA_convergence_criterion(ibin,   &
                                      iconverge_mass,   &
                                      iconverge_flux,   &
                                      idissolved)

      if(iconverge_mass .eq. mYES)then
        iter_MESA(ibin) = iter_MESA(ibin) + itdum
        niter_MESA = niter_MESA + float(itdum)
        niter_MESA_max = max(niter_MESA_max, itdum)
        jaerosolstate(ibin) = all_solid
        call adjust_solid_aerosol(ibin)
        jhyst_leg(ibin) = jhyst_lo
        growth_factor(ibin) = 1.0
        return
      elseif(iconverge_flux .eq. mYES)then
        iter_MESA(ibin) = iter_MESA(ibin) + itdum
        niter_MESA = niter_MESA + float(itdum)
        niter_MESA_max = max(niter_MESA_max, itdum)
        jaerosolstate(ibin) = mixed
        vol_wet_a(ibin)  = vol_dry_a(ibin) + water_a(ibin)*1.e-3		! cc(aer)/cc(air) or m^3/m^3(air)
        growth_factor(ibin) = mass_wet_a(ibin)/mass_dry_a(ibin)		! mass growth factor

        if(idissolved .eq. myes)then
          jaerosolstate(ibin) = all_liquid
!          jhyst_leg(ibin) = jhyst_up  ! ! do this later (to avoid tripping kelvin iterations)
        else
          jaerosolstate(ibin) = mixed
          jhyst_leg(ibin) = jhyst_lo
        endif

! calculate epercent(jsolid) composition in mixed-phase aerosol EFFI
!!        sum_dum = 0.0
!!        jp = jsolid
!!        do je = 1, nelectrolyte
!!          electrolyte(je,jp,ibin) = max(0.d0,electrolyte(je,jp,ibin)) ! remove -ve
!!          sum_dum = sum_dum + electrolyte(je,jp,ibin)
!!        enddo
!!        electrolyte_sum(jp,ibin) = sum_dum
!!        if(sum_dum .eq. 0.)sum_dum = 1.0
!!        do je = 1, nelectrolyte
!!          epercent(je,jp,ibin) = 100.*electrolyte(je,jp,ibin)/sum_dum
!!        enddo

        return
      endif


! calculate hsalt(js)	! time step
      hsalt_min = 1.e25
      do js = 1, nsalt

        phi_prod = phi_salt(js) * phi_salt_old(js)

        if(itdum .gt. 1 .and. phi_prod .gt. 0.0)then
          phi_bar(js) = (abs(phi_salt(js))-abs(phi_salt_old(js)))/   &
                                    alpha_salt(js)
        else
          phi_bar(js) = 0.0			! oscillating, or phi_salt and/or phi_salt_old may be zero
        endif

        if(phi_bar(js) .lt. 0.0)then		! good. phi getting lower. maybe able to take bigger alphas
          phi_bar(js) = max(phi_bar(js), -10.0d0)
          alpha_fac = 3.0*exp(phi_bar(js))
          alpha_salt(js) = min(alpha_fac*abs(phi_salt(js)), 0.9d0)
        elseif(phi_bar(js) .gt. 0.0)then	! bad - phi is getting bigger. so be conservative with alpha
           alpha_salt(js) = min(abs(phi_salt(js)), 0.5d0)
        else					! very bad - phi is oscillating. be very conservative
           alpha_salt(js) = min(abs(phi_salt(js))/3.0d0, 0.5d0)
        endif

!        alpha_salt(js) = max(alpha_salt(js), 0.01)

        phi_salt_old(js) = phi_salt(js)		! update old array


        if(flux_sl(js) .gt. 0.)then

          tau_p(js) = eleliquid(js)/flux_sl(js)	! precipitation time scale
          if(tau_p(js) .eq. 0.0)then
            hsalt(js) = 1.e25
            flux_sl(js) = 0.0
            phi_salt(js)= 0.0
          else
            hsalt(js) = alpha_salt(js)*tau_p(js)
          endif

        elseif(flux_sl(js) .lt. 0.)then

          tau_p(js) = -eleliquid(js)/flux_sl(js)	! precipitation time scale
          tau_d(js) = -electrolyte(js,jsolid,ibin)/flux_sl(js) ! dissolution time scale
          if(tau_p(js) .eq. 0.0)then
            hsalt(js) = alpha_salt(js)*tau_d(js)
          else
            hsalt(js) = alpha_salt(js)*min(tau_p(js),tau_d(js))
          endif

        else

          hsalt(js) = 1.e25

        endif

          hsalt_min = min(hsalt(js), hsalt_min)

      enddo

!---------------------------------

! integrate electrolyte(solid)
      do js = 1, nsalt
        electrolyte(js,jsolid,ibin) = (   &
                         (electrolyte(js,jsolid,ibin))  +   &
                         (hsalt(js)) * (flux_sl(js)) )
      enddo


! compute aer(solid) from electrolyte(solid)
      call electrolytes_to_ions(jsolid,ibin)


! compute new electrolyte(liquid) from mass balance
      do iaer = 1, naer
        aer(iaer,jliquid,ibin) = ( (aer(iaer,jtotal,ibin)) -   &
                                   (aer(iaer,jsolid,ibin)) )
      enddo

!---------------------------------



500   continue	! end time continuation loop
!--------------------------------------------------------------------
      jMESA_fail = jMESA_fail + 1
      iter_MESA(ibin) = iter_MESA(ibin) + itdum
      niter_MESA = niter_MESA + float(itdum)
      jaerosolstate(ibin) = mixed
      jhyst_leg(ibin) = jhyst_lo
      mass_wet_a(ibin) = mass_dry_a(ibin) + water_a(ibin)*1.e-3	! g/cc(air)
      vol_wet_a(ibin)  = vol_dry_a(ibin) + water_a(ibin)*1.e-3	! cc(aer)/cc(air) or m^3/m^3(air)
      growth_factor(ibin) = mass_wet_a(ibin)/mass_dry_a(ibin)	! mass growth factor

      return
      end subroutine MESA_PTC
















!***********************************************************************
! part of MESA: calculates solid-liquid fluxes of soluble salts
!
! author: Rahul A. Zaveri
! update: jan 2005
!-----------------------------------------------------------------------
      subroutine MESA_flux_salt(ibin)	! TOUCH
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      integer js, je
      real(r8) :: XT, calcium, sum_salt, sum_dum


! compute activities and water content
      call ions_to_electrolytes(jliquid,ibin,XT)
      call compute_activities(ibin)
      activity(jna3hso4,ibin)   = 0.0

      if(water_a(ibin) .le. 0.0)then
        do js = 1, nsalt
         flux_sl(js) = 0.0
        enddo
        return
      endif


      call MESA_estimate_eleliquid(ibin,XT)

      calcium = aer(ica_a,jliquid,ibin)



!! EFFI calculate percent composition
      sum_dum = 0.0
      do je = 1, nelectrolyte
        sum_dum = sum_dum + electrolyte(je,jliquid,ibin)
      enddo

      if(sum_dum .eq. 0.)sum_dum = 1.0

      do je = 1, nelectrolyte
        epercent(je,jliquid,ibin) = 100.*electrolyte(je,jliquid,ibin)/sum_dum
      enddo
!! EFFI



! calculate % electrolyte composition in the solid and liquid phases
      sum_salt = 0.0
      do js = 1, nsalt
        sum_salt = sum_salt + electrolyte(js,jsolid,ibin)
      enddo

      if(sum_salt .eq. 0.0)sum_salt = 1.0
      do js = 1, nsalt
        frac_salt_solid(js) = electrolyte(js,jsolid,ibin)/sum_salt
        frac_salt_liq(js)   = epercent(js,jliquid,ibin)/100.
      enddo



! compute salt fluxes
      do js = 1, nsalt		! soluble solid salts

! compute new saturation ratio
        sat_ratio(js) = activity(js,ibin)/Keq_sl(js)
! compute relative driving force
        phi_salt(js)  = (sat_ratio(js) - 1.0)/max(sat_ratio(js),1.0d0)

! check if too little solid-phase salt is trying to dissolve
        if(sat_ratio(js)       .lt. 1.00 .and.   &
           frac_salt_solid(js) .lt. 0.01 .and.   &
           frac_salt_solid(js) .gt. 0.0)then
          call MESA_dissolve_small_salt(ibin,js)
          call MESA_estimate_eleliquid(ibin,XT)
          sat_ratio(js) = activity(js,ibin)/Keq_sl(js)
        endif

! compute flux
        flux_sl(js) = sat_ratio(js) - 1.0

! apply Heaviside function
        if( (sat_ratio(js)               .lt. 1.0 .and.   &
             electrolyte(js,jsolid,ibin) .eq. 0.0) .or.   &
            (calcium .gt. 0.0 .and. frac_salt_liq(js).lt.0.01).or.   &
            (calcium .gt. 0.0 .and. jsalt_present(js).eq.0) )then
          flux_sl(js) = 0.0
          phi_salt(js)= 0.0
        endif

      enddo


! force cacl2 and cano3 fluxes to zero
      sat_ratio(jcano3) = 1.0
      phi_salt(jcano3)  = 0.0
      flux_sl(jcano3)   = 0.0

      sat_ratio(jcacl2) = 1.0
      phi_salt(jcacl2)  = 0.0
      flux_sl(jcacl2)   = 0.0


      return
      end subroutine MESA_flux_salt












!***********************************************************************
! part of MESA: calculates liquid electrolytes from ions
!
! notes:
!  - this subroutine is to be used for liquid-phase or total-phase only
!  - this sub transfers caso4 and caco3 from liquid to solid phase
!
! author: Rahul A. Zaveri
! update: jan 2005
!-----------------------------------------------------------------------
      subroutine MESA_estimate_eleliquid(ibin,XT)	! TOUCH
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin, jp
      real(r8) :: XT
! local variables
      integer iaer, je, jc, ja, icase
      real(r8) :: store(naer), sum_dum, sum_naza, sum_nczc, sum_na_nh4,   &
           f_nh4, f_na, xh, xb, xl, xs, XT_d, XNa_d, XNH4_d,   &
           xdum, dum, cat_net
      real(r8) :: nc(ncation), na(nanion)
      real(r8) :: dum_ca, dum_no3, dum_cl, cano3, cacl2



! remove negative concentrations, if any
      do iaer =  1, naer
      aer(iaer,jliquid,ibin) = max(0.0d0, aer(iaer,jliquid,ibin))
      enddo


! calculate sulfate ratio
      call calculate_XT(ibin,jliquid,XT)

      if(XT .ge. 2.0 .or. XT.lt.0.)then
       icase = 1	! near neutral (acidity is caused by HCl and/or HNO3)
      else
       icase = 2	! acidic (acidity is caused by excess SO4)
      endif


! initialize to zero
      do je = 1, nelectrolyte
        eleliquid(je) = 0.0
      enddo

!
!---------------------------------------------------------
! initialize moles of ions depending on the sulfate domain

      jp = jliquid

      if(icase.eq.1)then ! XT >= 2 : SULFATE POOR DOMAIN

        dum_ca  = aer(ica_a,jp,ibin)
        dum_no3 = aer(ino3_a,jp,ibin)
        dum_cl  = aer(icl_a,jp,ibin)

        cano3   = min(dum_ca, 0.5*dum_no3)
        dum_ca  = max(0.d0, dum_ca - cano3)
        dum_no3 = max(0.d0, dum_no3 - 2.*cano3)

        cacl2   = min(dum_ca, 0.5*dum_cl)
        dum_ca  = max(0.d0, dum_ca - cacl2)
        dum_cl  = max(0.d0, dum_cl - 2.*cacl2)

        na(ja_hso4)= 0.0
        na(ja_so4) = aer(iso4_a,jp,ibin)
        na(ja_no3) = aer(ino3_a,jp,ibin)
        na(ja_cl)  = aer(icl_a, jp,ibin)
        na(ja_msa) = aer(imsa_a,jp,ibin)

        nc(jc_ca)  = aer(ica_a, jp,ibin)
        nc(jc_na)  = aer(ina_a, jp,ibin)
        nc(jc_nh4) = aer(inh4_a,jp,ibin)

        cat_net = (   &
                 (2.*na(ja_so4)+na(ja_no3)+na(ja_cl)+na(ja_msa)) -   &
                 (2.*nc(jc_ca) +nc(jc_nh4)+nc(jc_na)) )			! bug fix 8/27/2013

        if(cat_net .lt. 0.0)then

          nc(jc_h) = 0.0

        else  ! cat_net must be 0.0 or positive

          nc(jc_h) = cat_net

        endif


! now compute equivalent fractions
      sum_naza = 0.0
      do ja = 1, nanion
        sum_naza = sum_naza + na(ja)*za(ja)
      enddo

      sum_nczc = 0.0
      do jc = 1, ncation
        sum_nczc = sum_nczc + nc(jc)*zc(jc)
      enddo


      !!!! MARK HERE
      if(sum_naza .eq. 0. .or. sum_nczc .eq. 0.)then
        !write(6,*)'ionic concentrations are zero'
        !write(6,*)'sum_naza = ', sum_naza
        !write(6,*)'sum_nczc = ', sum_nczc
        return
      endif

      do ja = 1, nanion
        xeq_a(ja) = na(ja)*za(ja)/sum_naza
      enddo

      do jc = 1, ncation
        xeq_c(jc) = nc(jc)*zc(jc)/sum_nczc
      enddo

      na_Ma(ja_so4) = na(ja_so4) *MW_a(ja_so4)
      na_Ma(ja_no3) = na(ja_no3) *MW_a(ja_no3)
      na_Ma(ja_cl)  = na(ja_cl)  *MW_a(ja_cl)
      na_Ma(ja_hso4)= na(ja_hso4)*MW_a(ja_hso4)
      na_Ma(ja_msa) = na(ja_msa) *MW_a(ja_msa)

      nc_Mc(jc_ca)  = nc(jc_ca) *MW_c(jc_ca)
      nc_Mc(jc_na)  = nc(jc_na) *MW_c(jc_na)
      nc_Mc(jc_nh4) = nc(jc_nh4)*MW_c(jc_nh4)
      nc_Mc(jc_h)   = nc(jc_h)  *MW_c(jc_h)


! now compute electrolyte moles
      eleliquid(jna2so4) = (xeq_c(jc_na) *na_Ma(ja_so4) +   &
                            xeq_a(ja_so4)*nc_Mc(jc_na))/   &
                             mw_electrolyte(jna2so4)

      eleliquid(jnahso4) = (xeq_c(jc_na) *na_Ma(ja_hso4) +   &
                            xeq_a(ja_hso4)*nc_Mc(jc_na))/   &
                             mw_electrolyte(jnahso4)

      eleliquid(jnamsa)  = (xeq_c(jc_na) *na_Ma(ja_msa) +   &
                            xeq_a(ja_msa)*nc_Mc(jc_na))/   &
                             mw_electrolyte(jnamsa)

      eleliquid(jnano3)  = (xeq_c(jc_na) *na_Ma(ja_no3) +   &
                            xeq_a(ja_no3)*nc_Mc(jc_na))/   &
                             mw_electrolyte(jnano3)

      eleliquid(jnacl)   = (xeq_c(jc_na) *na_Ma(ja_cl) +   &
                            xeq_a(ja_cl) *nc_Mc(jc_na))/   &
                             mw_electrolyte(jnacl)

      eleliquid(jnh4so4) = (xeq_c(jc_nh4)*na_Ma(ja_so4) +   &
                            xeq_a(ja_so4)*nc_Mc(jc_nh4))/   &
                             mw_electrolyte(jnh4so4)

      eleliquid(jnh4hso4)= (xeq_c(jc_nh4)*na_Ma(ja_hso4) +   &
                            xeq_a(ja_hso4)*nc_Mc(jc_nh4))/   &
                             mw_electrolyte(jnh4hso4)

      eleliquid(jnh4msa) = (xeq_c(jc_nh4) *na_Ma(ja_msa) +   &
                            xeq_a(ja_msa)*nc_Mc(jc_nh4))/   &
                             mw_electrolyte(jnh4msa)

      eleliquid(jnh4no3) = (xeq_c(jc_nh4)*na_Ma(ja_no3) +   &
                            xeq_a(ja_no3)*nc_Mc(jc_nh4))/   &
                             mw_electrolyte(jnh4no3)

      eleliquid(jnh4cl)  = (xeq_c(jc_nh4)*na_Ma(ja_cl) +   &
                            xeq_a(ja_cl) *nc_Mc(jc_nh4))/   &
                             mw_electrolyte(jnh4cl)

      eleliquid(jcamsa2) = (xeq_c(jc_ca) *na_Ma(ja_msa) +   &
                            xeq_a(ja_msa)*nc_Mc(jc_ca))/   &
                             mw_electrolyte(jcamsa2)

      eleliquid(jcano3)  = (xeq_c(jc_ca) *na_Ma(ja_no3) +   &
                            xeq_a(ja_no3)*nc_Mc(jc_ca))/   &
                             mw_electrolyte(jcano3)

      eleliquid(jcacl2)  = (xeq_c(jc_ca) *na_Ma(ja_cl) +   &
                            xeq_a(ja_cl) *nc_Mc(jc_ca))/   &
                             mw_electrolyte(jcacl2)

      eleliquid(jh2so4)  = (xeq_c(jc_h)   *na_Ma(ja_hso4) +   &
                            xeq_a(ja_hso4)*nc_Mc(jc_h))/   &
                             mw_electrolyte(jh2so4)

      eleliquid(jhno3)   = (xeq_c(jc_h)  *na_Ma(ja_no3) +   &
                            xeq_a(ja_no3)*nc_Mc(jc_h))/   &
                             mw_electrolyte(jhno3)

      eleliquid(jhcl)    = (xeq_c(jc_h) *na_Ma(ja_cl) +   &
                            xeq_a(ja_cl)*nc_Mc(jc_h))/   &
                             mw_electrolyte(jhcl)

      eleliquid(jmsa)    = (xeq_c(jc_h)  *na_Ma(ja_msa) +   &
                            xeq_a(ja_msa)*nc_Mc(jc_h))/   &
                             mw_electrolyte(jmsa)

!--------------------------------------------------------------------

      elseif(icase.eq.2)then ! XT < 2 : SULFATE RICH DOMAIN

        jp = jliquid

        store(iso4_a) = aer(iso4_a,jp,ibin)
        store(imsa_a) = aer(imsa_a,jp,ibin)
        store(inh4_a) = aer(inh4_a,jp,ibin)
        store(ina_a)  = aer(ina_a, jp,ibin)
        store(ica_a)  = aer(ica_a, jp,ibin)

        call form_camsa2(store,jp,ibin)

        sum_na_nh4 = store(ina_a) + store(inh4_a)
        if(sum_na_nh4 .gt. 0.0)then
          f_nh4 = store(inh4_a)/sum_na_nh4
          f_na  = store(ina_a)/sum_na_nh4
        else
          f_nh4 = 0.0
          f_na  = 0.0
        endif

! first form msa electrolytes
        if(sum_na_nh4 .gt. store(imsa_a))then
          eleliquid(jnh4msa) = f_nh4*store(imsa_a)
          eleliquid(jnamsa)  = f_na *store(imsa_a)
          store(inh4_a)= store(inh4_a)-eleliquid(jnh4msa) ! remaining nh4
          store(ina_a) = store(ina_a) -eleliquid(jnamsa)  ! remaining na
        else
          eleliquid(jnh4msa) = store(inh4_a)
          eleliquid(jnamsa)  = store(ina_a)
          eleliquid(jmsa)    = store(imsa_a) - sum_na_nh4
          store(inh4_a)= 0.0  ! remaining nh4
          store(ina_a) = 0.0  ! remaining na
        endif

        if(store(iso4_a).eq.0.0)goto 10

        XT_d  = XT
        XNa_d = 1. + 0.5*store(ina_a)/store(iso4_a)
        xdum = store(iso4_a) - store(inh4_a)

        dum = ( (2.*store(iso4_a)) -   &
                (store(ina_a)) )
        if(store(inh4_a) .gt. 0.0 .and. dum .gt. 0.0)then
          XNH4_d = 2.*store(inh4_a)/   &
                  (2.*store(iso4_a) - store(ina_a))
        else
          XNH4_d = 0.0
        endif


        IF(store(inh4_a) .gt. 0.0)THEN


        if(XT_d .ge. XNa_d)then
          eleliquid(jna2so4) = 0.5*store(ina_a)

          if(XNH4_d .ge. 5./3.)then
            eleliquid(jnh4so4) = 1.5*store(ina_a)   &
                               - 3.*xdum - store(inh4_a)
            eleliquid(jlvcite) = 2.*xdum + store(inh4_a)   &
                               - store(ina_a)
          elseif(XNH4_d .ge. 1.5)then
            eleliquid(jnh4so4) = store(inh4_a)/5.
            eleliquid(jlvcite) = store(inh4_a)/5.
          elseif(XNH4_d .ge. 1.0)then
            eleliquid(jnh4so4) = store(inh4_a)/6.
            eleliquid(jlvcite) = store(inh4_a)/6.
            eleliquid(jnh4hso4)= store(inh4_a)/6.
          endif

        elseif(XT_d .gt. 1.0)then
          eleliquid(jnh4so4)  = store(inh4_a)/6.
          eleliquid(jlvcite)  = store(inh4_a)/6.
          eleliquid(jnh4hso4) = store(inh4_a)/6.
          eleliquid(jna2so4)  = store(ina_a)/3.
          eleliquid(jnahso4)  = store(ina_a)/3.
        elseif(XT_d .le. 1.0)then
          eleliquid(jna2so4)  = store(ina_a)/4.
          eleliquid(jnahso4)  = store(ina_a)/2.
          eleliquid(jlvcite)  = store(inh4_a)/6.
          eleliquid(jnh4hso4) = store(inh4_a)/2.
        endif

        ELSE

        if(XT_d .gt. 1.0)then
          eleliquid(jna2so4) = store(ina_a) - store(iso4_a)
          eleliquid(jnahso4) = 2.*store(iso4_a) -   &
                                  store(ina_a)
        else
          eleliquid(jna2so4) = store(ina_a)/4.
          eleliquid(jnahso4) = store(ina_a)/2.
        endif


        ENDIF



      endif
!---------------------------------------------------------


10    return
      end subroutine MESA_estimate_eleliquid










!***********************************************************************
! part of MESA: completely dissolves small amounts of soluble salts
!
! author: Rahul A. Zaveri
! update: jan 2005
!-----------------------------------------------------------------------
      subroutine MESA_dissolve_small_salt(ibin,js)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin, js, jp

      jp = jsolid


      if(js .eq. jnh4so4)then
        aer(inh4_a,jliquid,ibin) = aer(inh4_a,jliquid,ibin) +   &
                           2.*electrolyte(js,jsolid,ibin)
        aer(iso4_a,jliquid,ibin) = aer(iso4_a,jliquid,ibin) +   &
                              electrolyte(js,jsolid,ibin)

        electrolyte(js,jsolid,ibin) = 0.0

        aer(inh4_a,jp,ibin) = electrolyte(jnh4no3,jp,ibin) +   &
                             electrolyte(jnh4cl,jp,ibin)  +   &
                          2.*electrolyte(jnh4so4,jp,ibin) +   &
                          3.*electrolyte(jlvcite,jp,ibin) +   &
                             electrolyte(jnh4hso4,jp,ibin)+   &
                             electrolyte(jnh4msa,jp,ibin)

        aer(iso4_a,jp,ibin) = electrolyte(jcaso4,jp,ibin)  +   &
                             electrolyte(jna2so4,jp,ibin) +   &
                          2.*electrolyte(jna3hso4,jp,ibin)+   &
                             electrolyte(jnahso4,jp,ibin) +   &
                             electrolyte(jnh4so4,jp,ibin) +   &
                          2.*electrolyte(jlvcite,jp,ibin) +   &
                             electrolyte(jnh4hso4,jp,ibin)+   &
                             electrolyte(jh2so4,jp,ibin)
        return
      endif


      if(js .eq. jlvcite)then
        aer(inh4_a,jliquid,ibin) = aer(inh4_a,jliquid,ibin) +   &
                           3.*electrolyte(js,jsolid,ibin)
        aer(iso4_a,jliquid,ibin) = aer(iso4_a,jliquid,ibin) +   &
                           2.*electrolyte(js,jsolid,ibin)

        electrolyte(js,jsolid,ibin) = 0.0

        aer(inh4_a,jp,ibin) = electrolyte(jnh4no3,jp,ibin) +   &
                             electrolyte(jnh4cl,jp,ibin)  +   &
                          2.*electrolyte(jnh4so4,jp,ibin) +   &
                          3.*electrolyte(jlvcite,jp,ibin) +   &
                             electrolyte(jnh4hso4,jp,ibin)+   &
                             electrolyte(jnh4msa,jp,ibin)

        aer(iso4_a,jp,ibin) = electrolyte(jcaso4,jp,ibin)  +   &
                             electrolyte(jna2so4,jp,ibin) +   &
                          2.*electrolyte(jna3hso4,jp,ibin)+   &
                             electrolyte(jnahso4,jp,ibin) +   &
                             electrolyte(jnh4so4,jp,ibin) +   &
                          2.*electrolyte(jlvcite,jp,ibin) +   &
                             electrolyte(jnh4hso4,jp,ibin)+   &
                             electrolyte(jh2so4,jp,ibin)
        return
      endif


      if(js .eq. jnh4hso4)then
        aer(inh4_a,jliquid,ibin) = aer(inh4_a,jliquid,ibin) +   &
                              electrolyte(js,jsolid,ibin)
        aer(iso4_a,jliquid,ibin) = aer(iso4_a,jliquid,ibin) +   &
                             electrolyte(js,jsolid,ibin)

        electrolyte(js,jsolid,ibin) = 0.0

        aer(inh4_a,jp,ibin) = electrolyte(jnh4no3,jp,ibin) +   &
                             electrolyte(jnh4cl,jp,ibin)  +   &
                          2.*electrolyte(jnh4so4,jp,ibin) +   &
                          3.*electrolyte(jlvcite,jp,ibin) +   &
                             electrolyte(jnh4hso4,jp,ibin)+   &
                             electrolyte(jnh4msa,jp,ibin)

        aer(iso4_a,jp,ibin) = electrolyte(jcaso4,jp,ibin)  +   &
                             electrolyte(jna2so4,jp,ibin) +   &
                          2.*electrolyte(jna3hso4,jp,ibin)+   &
                             electrolyte(jnahso4,jp,ibin) +   &
                             electrolyte(jnh4so4,jp,ibin) +   &
                          2.*electrolyte(jlvcite,jp,ibin) +   &
                             electrolyte(jnh4hso4,jp,ibin)+   &
                             electrolyte(jh2so4,jp,ibin)
        return
      endif


      if(js .eq. jna2so4)then
        aer(ina_a,jliquid,ibin)  = aer(ina_a,jliquid,ibin) +   &
                           2.*electrolyte(js,jsolid,ibin)
        aer(iso4_a,jliquid,ibin) = aer(iso4_a,jliquid,ibin) +   &
                              electrolyte(js,jsolid,ibin)

        electrolyte(js,jsolid,ibin) = 0.0

        aer(ina_a,jp,ibin)  = electrolyte(jnano3,jp,ibin)  +   &
                             electrolyte(jnacl,jp,ibin)   +   &
                          2.*electrolyte(jna2so4,jp,ibin) +   &
                          3.*electrolyte(jna3hso4,jp,ibin)+   &
                             electrolyte(jnahso4,jp,ibin) +   &
                             electrolyte(jnamsa,jp,ibin)

        aer(iso4_a,jp,ibin) = electrolyte(jcaso4,jp,ibin)  +   &
                             electrolyte(jna2so4,jp,ibin) +   &
                          2.*electrolyte(jna3hso4,jp,ibin)+   &
                             electrolyte(jnahso4,jp,ibin) +   &
                             electrolyte(jnh4so4,jp,ibin) +   &
                          2.*electrolyte(jlvcite,jp,ibin) +   &
                             electrolyte(jnh4hso4,jp,ibin)+   &
                             electrolyte(jh2so4,jp,ibin)
        return
      endif


      if(js .eq. jna3hso4)then
        aer(ina_a,jliquid,ibin)  = aer(ina_a,jliquid,ibin) +   &
                           3.*electrolyte(js,jsolid,ibin)
        aer(iso4_a,jliquid,ibin) = aer(iso4_a,jliquid,ibin) +   &
                           2.*electrolyte(js,jsolid,ibin)

        electrolyte(js,jsolid,ibin) = 0.0

        aer(ina_a,jp,ibin)  = electrolyte(jnano3,jp,ibin)  +   &
                             electrolyte(jnacl,jp,ibin)   +   &
                          2.*electrolyte(jna2so4,jp,ibin) +   &
                          3.*electrolyte(jna3hso4,jp,ibin)+   &
                             electrolyte(jnahso4,jp,ibin) +   &
                             electrolyte(jnamsa,jp,ibin)

        aer(iso4_a,jp,ibin) = electrolyte(jcaso4,jp,ibin)  +   &
                             electrolyte(jna2so4,jp,ibin) +   &
                          2.*electrolyte(jna3hso4,jp,ibin)+   &
                             electrolyte(jnahso4,jp,ibin) +   &
                             electrolyte(jnh4so4,jp,ibin) +   &
                          2.*electrolyte(jlvcite,jp,ibin) +   &
                             electrolyte(jnh4hso4,jp,ibin)+   &
                             electrolyte(jh2so4,jp,ibin)
        return
      endif


      if(js .eq. jnahso4)then
        aer(ina_a,jliquid,ibin)  = aer(ina_a,jliquid,ibin) +   &
                              electrolyte(js,jsolid,ibin)
        aer(iso4_a,jliquid,ibin) = aer(iso4_a,jliquid,ibin) +   &
                              electrolyte(js,jsolid,ibin)

        electrolyte(js,jsolid,ibin) = 0.0

        aer(ina_a,jp,ibin)  = electrolyte(jnano3,jp,ibin)  +   &
                             electrolyte(jnacl,jp,ibin)   +   &
                          2.*electrolyte(jna2so4,jp,ibin) +   &
                          3.*electrolyte(jna3hso4,jp,ibin)+   &
                             electrolyte(jnahso4,jp,ibin) +   &
                             electrolyte(jnamsa,jp,ibin)

        aer(iso4_a,jp,ibin) = electrolyte(jcaso4,jp,ibin)  +   &
                             electrolyte(jna2so4,jp,ibin) +   &
                          2.*electrolyte(jna3hso4,jp,ibin)+   &
                             electrolyte(jnahso4,jp,ibin) +   &
                             electrolyte(jnh4so4,jp,ibin) +   &
                          2.*electrolyte(jlvcite,jp,ibin) +   &
                             electrolyte(jnh4hso4,jp,ibin)+   &
                             electrolyte(jh2so4,jp,ibin)
        return
      endif


      if(js .eq. jnh4no3)then
        aer(inh4_a,jliquid,ibin) = aer(inh4_a,jliquid,ibin) +   &
                              electrolyte(js,jsolid,ibin)
        aer(ino3_a,jliquid,ibin) = aer(ino3_a,jliquid,ibin) +   &
                              electrolyte(js,jsolid,ibin)

        electrolyte(js,jsolid,ibin) = 0.0

        aer(inh4_a,jp,ibin) = electrolyte(jnh4no3,jp,ibin) +   &
                             electrolyte(jnh4cl,jp,ibin)  +   &
                          2.*electrolyte(jnh4so4,jp,ibin) +   &
                          3.*electrolyte(jlvcite,jp,ibin) +   &
                             electrolyte(jnh4hso4,jp,ibin)+   &
                             electrolyte(jnh4msa,jp,ibin)

        aer(ino3_a,jp,ibin) = electrolyte(jnano3,jp,ibin)  +   &
                          2.*electrolyte(jcano3,jp,ibin)  +   &
                             electrolyte(jnh4no3,jp,ibin) +   &
                             electrolyte(jhno3,jp,ibin)
        return
      endif


      if(js .eq. jnh4cl)then
        aer(inh4_a,jliquid,ibin) = aer(inh4_a,jliquid,ibin) +   &
                              electrolyte(js,jsolid,ibin)
        aer(icl_a,jliquid,ibin)  = aer(icl_a,jliquid,ibin) +   &
                              electrolyte(js,jsolid,ibin)

        electrolyte(js,jsolid,ibin) = 0.0

        aer(inh4_a,jp,ibin) = electrolyte(jnh4no3,jp,ibin) +   &
                             electrolyte(jnh4cl,jp,ibin)  +   &
                          2.*electrolyte(jnh4so4,jp,ibin) +   &
                          3.*electrolyte(jlvcite,jp,ibin) +   &
                             electrolyte(jnh4hso4,jp,ibin)+   &
                             electrolyte(jnh4msa,jp,ibin)

        aer(icl_a,jp,ibin)  = electrolyte(jnacl,jp,ibin)   +   &
                          2.*electrolyte(jcacl2,jp,ibin)  +   &
                             electrolyte(jnh4cl,jp,ibin)  +   &
                             electrolyte(jhcl,jp,ibin)
        return
      endif


      if(js .eq. jnano3)then
        aer(ina_a,jliquid,ibin)  = aer(ina_a,jliquid,ibin) +   &
                              electrolyte(js,jsolid,ibin)
        aer(ino3_a,jliquid,ibin) = aer(ino3_a,jliquid,ibin) +   &
                              electrolyte(js,jsolid,ibin)

        electrolyte(js,jsolid,ibin) = 0.0

        aer(ina_a,jp,ibin)  = electrolyte(jnano3,jp,ibin)  +   &
                             electrolyte(jnacl,jp,ibin)   +   &
                          2.*electrolyte(jna2so4,jp,ibin) +   &
                          3.*electrolyte(jna3hso4,jp,ibin)+   &
                             electrolyte(jnahso4,jp,ibin) +   &
                             electrolyte(jnamsa,jp,ibin)

        aer(ino3_a,jp,ibin) = electrolyte(jnano3,jp,ibin)  +   &
                          2.*electrolyte(jcano3,jp,ibin)  +   &
                             electrolyte(jnh4no3,jp,ibin) +   &
                             electrolyte(jhno3,jp,ibin)
        return
      endif


      if(js .eq. jnacl)then
        aer(ina_a,jliquid,ibin)  = aer(ina_a,jliquid,ibin) +   &
                              electrolyte(js,jsolid,ibin)
        aer(icl_a,jliquid,ibin)  = aer(icl_a,jliquid,ibin) +   &
                              electrolyte(js,jsolid,ibin)

        electrolyte(js,jsolid,ibin) = 0.0

        aer(ina_a,jp,ibin)  = electrolyte(jnano3,jp,ibin)  +   &
                             electrolyte(jnacl,jp,ibin)   +   &
                          2.*electrolyte(jna2so4,jp,ibin) +   &
                          3.*electrolyte(jna3hso4,jp,ibin)+   &
                             electrolyte(jnahso4,jp,ibin) +   &
                             electrolyte(jnamsa,jp,ibin)

        aer(icl_a,jp,ibin)  = electrolyte(jnacl,jp,ibin)   +   &
                          2.*electrolyte(jcacl2,jp,ibin)  +   &
                             electrolyte(jnh4cl,jp,ibin)  +   &
                             electrolyte(jhcl,jp,ibin)
        return
      endif


      if(js .eq. jcano3)then
        aer(ica_a,jliquid,ibin)  = aer(ica_a,jliquid,ibin) +   &
                              electrolyte(js,jsolid,ibin)
        aer(ino3_a,jliquid,ibin) = aer(ino3_a,jliquid,ibin) +   &
                            2.*electrolyte(js,jsolid,ibin)

        electrolyte(js,jsolid,ibin) = 0.0

        aer(ica_a,jp,ibin)  = electrolyte(jcaso4,jp,ibin)  +   &
                             electrolyte(jcano3,jp,ibin)  +   &
                             electrolyte(jcacl2,jp,ibin)  +   &
                             electrolyte(jcaco3,jp,ibin)  +   &
                             electrolyte(jcamsa2,jp,ibin)

        aer(ino3_a,jp,ibin) = electrolyte(jnano3,jp,ibin)  +   &
                          2.*electrolyte(jcano3,jp,ibin)  +   &
                             electrolyte(jnh4no3,jp,ibin) +   &
                             electrolyte(jhno3,jp,ibin)
        return
      endif


      if(js .eq. jcacl2)then
        aer(ica_a,jliquid,ibin) = aer(ica_a,jliquid,ibin) +   &
                              electrolyte(js,jsolid,ibin)
        aer(icl_a,jliquid,ibin) = aer(icl_a,jliquid,ibin) +   &
                            2.*electrolyte(js,jsolid,ibin)

        electrolyte(js,jsolid,ibin) = 0.0

        aer(ica_a,jp,ibin)  = electrolyte(jcaso4,jp,ibin)  +   &
                             electrolyte(jcano3,jp,ibin)  +   &
                             electrolyte(jcacl2,jp,ibin)  +   &
                             electrolyte(jcaco3,jp,ibin)  +   &
                             electrolyte(jcamsa2,jp,ibin)

        aer(icl_a,jp,ibin)  = electrolyte(jnacl,jp,ibin)   +   &
                           2.*electrolyte(jcacl2,jp,ibin)  +   &
                             electrolyte(jnh4cl,jp,ibin)  +   &
                             electrolyte(jhcl,jp,ibin)
        return
      endif



      return
      end subroutine MESA_dissolve_small_salt






!***********************************************************************
! part of MESA: checks MESA convergence
!
! author: Rahul A. Zaveri
! update: jan 2005
!-----------------------------------------------------------------------
      subroutine MESA_convergence_criterion(ibin,	   &  ! TOUCH
                                       iconverge_mass,   &
                                       iconverge_flux,   &
                                       idissolved)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin, iconverge_mass, iconverge_flux, idissolved
! local variables
      integer je, js, iaer
      real(r8) :: mass_solid, mass_solid_salt,frac_solid, XT, H_ion,   &
           crustal_solids, sumflux


      idissolved = mno		! default = not completely dissolved

! check mass convergence
      iconverge_mass = mNO	! default value = no convergence

!      call electrolytes_to_ions(jsolid,ibin)
!      mass_solid = 0.0
!      do iaer = 1, naer
!        mass_solid = mass_solid +
!     &               aer(iaer,jsolid,ibin)*mw_aer_mac(iaer)*1.e-15	! g/cc(air)
!      enddo

      mass_solid_salt = 0.0
      do je = 1, nsalt
        mass_solid_salt = mass_solid_salt +   &
              electrolyte(je,jsolid,ibin)*mw_electrolyte(je)*1.e-15	! g/cc(air)
      enddo



!      frac_solid = mass_solid/mass_dry_a(ibin)

      frac_solid = mass_solid_salt/mass_dry_salt(ibin)

      if(frac_solid .ge. 0.98)then
        iconverge_mass = mYES
        return
      endif



! check relative driving force convergence
      iconverge_flux = mYES
      do js = 1, nsalt
        if(abs(phi_salt(js)).gt. rtol_mesa)then
          iconverge_flux = mNO
          return
        endif
      enddo



! check if all the fluxes are zero

      sumflux = 0.0
      do js = 1, nsalt
        sumflux = sumflux + abs(flux_sl(js))
      enddo

      crustal_solids = electrolyte(jcaco3,jsolid,ibin) +   &
                       electrolyte(jcaso4,jsolid,ibin) +   &
                       aer(ioin_a,jsolid,ibin)

      if(sumflux .eq. 0.0 .and. crustal_solids .eq. 0.0)then
        idissolved = myes
      endif



      return
      end subroutine MESA_convergence_criterion








!***********************************************************************
! called when aerosol bin is completely solid.
!
! author: Rahul A. Zaveri
! update: jan 2005
!-----------------------------------------------------------------------
      subroutine adjust_solid_aerosol(ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      integer iaer, je


      jphase(ibin)    = jsolid
      jhyst_leg(ibin) = jhyst_lo	! lower curve
      water_a(ibin)   = 0.0

! transfer aer(jtotal) to aer(jsolid)
      do iaer = 1, naer
        aer(iaer, jsolid, ibin) = aer(iaer,jtotal,ibin)
        aer(iaer, jliquid,ibin) = 0.0
      enddo

! transfer electrolyte(jtotal) to electrolyte(jsolid)
      do je = 1, nelectrolyte
        electrolyte(je,jliquid,ibin) = 0.0
        epercent(je,jliquid,ibin)    = 0.0
        electrolyte(je,jsolid,ibin)  = electrolyte(je,jtotal,ibin)
        epercent(je,jsolid,ibin)     = epercent(je,jtotal,ibin)
      enddo

! update aer(jtotal) that may have been affected above
      aer(inh4_a,jtotal,ibin) = aer(inh4_a,jsolid,ibin)
      aer(ino3_a,jtotal,ibin) = aer(ino3_a,jsolid,ibin)
      aer(icl_a,jtotal,ibin)  = aer(icl_a,jsolid,ibin)


      return
      end subroutine adjust_solid_aerosol









!***********************************************************************
! called when aerosol bin is completely liquid.
!
! author: Rahul A. Zaveri
! update: jan 2005
!-----------------------------------------------------------------------
      subroutine adjust_liquid_aerosol(ibin)	! TOUCH
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      integer je




      jphase(ibin)    = jliquid
      jhyst_leg(ibin) = jhyst_up	! upper curve

! partition all electrolytes into liquid phase
      do je = 1, nelectrolyte
        electrolyte(je,jsolid,ibin)  = 0.0
        epercent(je,jsolid,ibin)     = 0.0
        electrolyte(je,jliquid,ibin) = electrolyte(je,jtotal,ibin)
        epercent(je,jliquid,ibin)    = epercent(je,jtotal,ibin)
      enddo
! except these electrolytes, which always remain in the solid phase
      electrolyte(jcaco3,jsolid,ibin) = electrolyte(jcaco3,jtotal,ibin)
      electrolyte(jcaso4,jsolid,ibin) = electrolyte(jcaso4,jtotal,ibin)
      epercent(jcaco3,jsolid,ibin)    = epercent(jcaco3,jtotal,ibin)
      epercent(jcaso4,jsolid,ibin)    = epercent(jcaso4,jtotal,ibin)
      electrolyte(jcaco3,jliquid,ibin)= 0.0
      electrolyte(jcaso4,jliquid,ibin)= 0.0
      epercent(jcaco3,jliquid,ibin)   = 0.0
      epercent(jcaso4,jliquid,ibin)   = 0.0


! partition all the aer species into
! solid phase
      aer(iso4_a,jsolid,ibin) = electrolyte(jcaso4,jsolid,ibin)
      aer(ino3_a,jsolid,ibin) = 0.0
      aer(icl_a,jsolid,ibin)  = 0.0
      aer(inh4_a,jsolid,ibin) = 0.0
      aer(ioc_a,jsolid,ibin)  = aer(ioc_a,jtotal,ibin)
      aer(imsa_a,jsolid,ibin) = 0.0
      aer(ico3_a,jsolid,ibin) = aer(ico3_a,jtotal,ibin)
      aer(ina_a,jsolid,ibin)  = 0.0
      aer(ica_a,jsolid,ibin)  = electrolyte(jcaco3,jsolid,ibin) +   &
                                electrolyte(jcaso4,jsolid,ibin)
      aer(ibc_a,jsolid,ibin)  = aer(ibc_a,jtotal,ibin)
      aer(ioin_a,jsolid,ibin) = aer(ioin_a,jtotal,ibin)
      aer(icn3_a,jsolid,ibin) = aer(icn3_a,jtotal,ibin)
      aer(icn2_a,jsolid,ibin) = aer(icn2_a,jtotal,ibin)
      aer(icn1_a,jsolid,ibin) = aer(icn1_a,jtotal,ibin)
      aer(ic0_a,jsolid,ibin)  = aer(ic0_a,jtotal,ibin)
      aer(ic1_a,jsolid,ibin)= aer(ic1_a,jtotal,ibin)
      aer(ic2_a,jsolid,ibin)= aer(ic2_a,jtotal,ibin)
      aer(ic3_a,jsolid,ibin)= aer(ic3_a,jtotal,ibin)
      aer(ic4_a,jsolid,ibin)= aer(ic4_a,jtotal,ibin)
      aer(ic5_a,jsolid,ibin)= aer(ic5_a,jtotal,ibin)
      aer(ic6_a,jsolid,ibin)= aer(ic6_a,jtotal,ibin)
      aer(ic7_a,jsolid,ibin)= aer(ic7_a,jtotal,ibin)
      aer(ic8_a,jsolid,ibin)= aer(ic8_a,jtotal,ibin)
      aer(ic9_a,jsolid,ibin)= aer(ic9_a,jtotal,ibin)

! liquid-phase
      aer(iso4_a,jliquid,ibin) = aer(iso4_a,jtotal,ibin) -   &
                                 aer(iso4_a,jsolid,ibin)
      aer(iso4_a,jliquid,ibin) = max(0.d0, aer(iso4_a,jliquid,ibin))
      aer(ino3_a,jliquid,ibin) = aer(ino3_a,jtotal,ibin)
      aer(icl_a,jliquid,ibin)  = aer(icl_a,jtotal,ibin)
      aer(inh4_a,jliquid,ibin) = aer(inh4_a,jtotal,ibin)
      aer(ioc_a,jliquid,ibin)  = 0.0
      aer(imsa_a,jliquid,ibin) = aer(imsa_a,jtotal,ibin)
      aer(ico3_a,jliquid,ibin) = 0.0
      aer(ina_a,jliquid,ibin)  = aer(ina_a,jtotal,ibin)
      aer(ica_a,jliquid,ibin)  = aer(ica_a,jtotal,ibin) -   &
                                 aer(ica_a,jsolid,ibin)
      aer(ica_a,jliquid,ibin)  = max(0.d0, aer(ica_a,jliquid,ibin))
      aer(ibc_a,jliquid,ibin)  = 0.0
      aer(ioin_a,jliquid,ibin) = 0.0
      aer(icn3_a,jliquid,ibin)= 0.0
      aer(icn2_a,jliquid,ibin)= 0.0
      aer(icn1_a,jliquid,ibin)= 0.0
	  aer(ic0_a,jliquid,ibin)= 0.0
      aer(ic1_a,jliquid,ibin)= 0.0
      aer(ic2_a,jliquid,ibin)= 0.0
      aer(ic3_a,jliquid,ibin)= 0.0
      aer(ic4_a,jliquid,ibin)= 0.0
      aer(ic5_a,jliquid,ibin)= 0.0
      aer(ic6_a,jliquid,ibin)= 0.0
      aer(ic7_a,jliquid,ibin)= 0.0
      aer(ic8_a,jliquid,ibin)= 0.0
      aer(ic9_a,jliquid,ibin)= 0.0

      return
      end subroutine adjust_liquid_aerosol







! end of MESA package
!=======================================================================







!***********************************************************************
! computes activities
!
! author: Rahul A. Zaveri
! update: jan 2007
!-----------------------------------------------------------------------
      subroutine compute_activities(ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      integer jp, jA
      real(r8) :: XT, xmol(Nelectrolyte), sum_elec, dumK, c_bal, a_c
      real(r8) :: quad, aq, bq, cq, xq, dum
      real(r8) :: aerosol_water	! mosaic function


      water_a(ibin) = aerosol_water(jliquid,ibin)	! Kg/m^3(air)
      if(water_a(ibin) .eq. 0.0)return


      call calculate_XT(ibin,jliquid,XT)


      if(XT.gt.2.0 .or. XT.lt.0.)then
! SULFATE POOR: fully dissociated electrolytes


! anion molalities (mol/kg water)
      ma(ja_so4,ibin)  = 1.e-9*aer(iso4_a,jliquid,ibin)/water_a(ibin)
      ma(ja_hso4,ibin) = 0.0
      ma(ja_no3,ibin)  = 1.e-9*aer(ino3_a,jliquid,ibin)/water_a(ibin)
      ma(ja_cl,ibin)   = 1.e-9*aer(icl_a, jliquid,ibin)/water_a(ibin)
      ma(ja_msa,ibin)  = 1.e-9*aer(imsa_a,jliquid,ibin)/water_a(ibin)

! cation molalities (mol/kg water)
      mc(jc_ca,ibin)   = 1.e-9*aer(ica_a, jliquid,ibin)/water_a(ibin)
      mc(jc_nh4,ibin)  = 1.e-9*aer(inh4_a,jliquid,ibin)/water_a(ibin)
      mc(jc_na,ibin)   = 1.e-9*aer(ina_a, jliquid,ibin)/water_a(ibin)
      a_c              = (   &
                         (2.*ma(ja_so4,ibin)+   &
                             ma(ja_no3,ibin)+   &
                             ma(ja_cl,ibin) +   &
                                 ma(ja_msa,ibin)) -   &
                         (2.*mc(jc_ca,ibin) +   &
                             mc(jc_nh4,ibin)+   &
                             mc(jc_na,ibin)) )
      mc(jc_h,ibin) = 0.5*( (a_c) +   &
                            (sqrt(a_c**2 + 4.*Keq_ll(3))) )

      if(mc(jc_h,ibin) .eq. 0.0)then
        mc(jc_h,ibin) = 1.e-10
      endif


      jp = jliquid


      sum_elec = 2.*electrolyte(jnh4no3,jp,ibin) +   &
                 2.*electrolyte(jnh4cl,jp,ibin)  +   &
                 3.*electrolyte(jnh4so4,jp,ibin) +   &
                 3.*electrolyte(jna2so4,jp,ibin) +   &
                 2.*electrolyte(jnano3,jp,ibin)  +   &
                 2.*electrolyte(jnacl,jp,ibin)   +   &
                 3.*electrolyte(jcano3,jp,ibin)  +   &
                 3.*electrolyte(jcacl2,jp,ibin)  +   &
                 2.*electrolyte(jhno3,jp,ibin)   +   &
                 2.*electrolyte(jhcl,jp,ibin)

      if(sum_elec .eq. 0.0)then
        do jA = 1, nelectrolyte
          gam(jA,ibin) = 1.0
        enddo
        goto 10
      endif


! ionic mole fractions
      xmol(jnh4no3) = 2.*electrolyte(jnh4no3,jp,ibin)/sum_elec
      xmol(jnh4cl)  = 2.*electrolyte(jnh4cl,jp,ibin) /sum_elec
      xmol(jnh4so4) = 3.*electrolyte(jnh4so4,jp,ibin)/sum_elec
      xmol(jna2so4) = 3.*electrolyte(jna2so4,jp,ibin)/sum_elec
      xmol(jnano3)  = 2.*electrolyte(jnano3,jp,ibin) /sum_elec
      xmol(jnacl)   = 2.*electrolyte(jnacl,jp,ibin)  /sum_elec
      xmol(jcano3)  = 3.*electrolyte(jcano3,jp,ibin) /sum_elec
      xmol(jcacl2)  = 3.*electrolyte(jcacl2,jp,ibin) /sum_elec
      xmol(jhno3)   = 2.*electrolyte(jhno3,jp,ibin)  /sum_elec
      xmol(jhcl)    = 2.*electrolyte(jhcl,jp,ibin)   /sum_elec


      jA = jnh4so4
      if(xmol(jA).gt.0.0)then
      log_gam(jA) = xmol(jnh4no3)*log_gamZ(jA,jnh4no3) +   &
                    xmol(jnh4cl) *log_gamZ(jA,jnh4cl)  +   &
                    xmol(jnh4so4)*log_gamZ(jA,jnh4so4) +   &
                    xmol(jna2so4)*log_gamZ(jA,jna2so4) +   &
                    xmol(jnano3) *log_gamZ(jA,jnano3)  +   &
                    xmol(jnacl)  *log_gamZ(jA,jnacl)   +   &
                    xmol(jcano3) *log_gamZ(jA,jcano3)  +   &
                    xmol(jcacl2) *log_gamZ(jA,jcacl2)  +   &
                    xmol(jhno3)  *log_gamZ(jA,jhno3)   +   &
                    xmol(jhcl)   *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)
      activity(jnh4so4,ibin) = mc(jc_nh4,ibin)**2 * ma(ja_so4,ibin) *   &
                               gam(jnh4so4,ibin)**3
      endif



      jA = jnh4no3
      if(xmol(jA).gt.0.0)then
      log_gam(jA) = xmol(jnh4no3)*log_gamZ(jA,jnh4no3) +   &
                    xmol(jnh4cl) *log_gamZ(jA,jnh4cl)  +   &
                    xmol(jnh4so4)*log_gamZ(jA,jnh4so4) +   &
                    xmol(jna2so4)*log_gamZ(jA,jna2so4) +   &
                    xmol(jnano3) *log_gamZ(jA,jnano3)  +   &
                    xmol(jnacl)  *log_gamZ(jA,jnacl)   +   &
                    xmol(jcano3) *log_gamZ(jA,jcano3)  +   &
                    xmol(jcacl2) *log_gamZ(jA,jcacl2)  +   &
                    xmol(jhno3)  *log_gamZ(jA,jhno3)   +   &
                    xmol(jhcl)   *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)
      activity(jnh4no3,ibin) = mc(jc_nh4,ibin) * ma(ja_no3,ibin) *   &
                               gam(jnh4no3,ibin)**2
      endif


      jA = jnh4cl
      if(xmol(jA).gt.0.0)then
      log_gam(jA) = xmol(jnh4no3)*log_gamZ(jA,jnh4no3) +   &
                    xmol(jnh4cl) *log_gamZ(jA,jnh4cl)  +   &
                    xmol(jnh4so4)*log_gamZ(jA,jnh4so4) +   &
                    xmol(jna2so4)*log_gamZ(jA,jna2so4) +   &
                    xmol(jnano3) *log_gamZ(jA,jnano3)  +   &
                    xmol(jnacl)  *log_gamZ(jA,jnacl)   +   &
                    xmol(jcano3) *log_gamZ(jA,jcano3)  +   &
                    xmol(jcacl2) *log_gamZ(jA,jcacl2)  +   &
                    xmol(jhno3)  *log_gamZ(jA,jhno3)   +   &
                    xmol(jhcl)   *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)
      activity(jnh4cl,ibin)  = mc(jc_nh4,ibin) * ma(ja_cl,ibin) *   &
                               gam(jnh4cl,ibin)**2
      endif


      jA = jna2so4
      if(xmol(jA).gt.0.0)then
      log_gam(jA) = xmol(jnh4no3)*log_gamZ(jA,jnh4no3) +   &
                    xmol(jnh4cl) *log_gamZ(jA,jnh4cl)  +   &
                    xmol(jnh4so4)*log_gamZ(jA,jnh4so4) +   &
                    xmol(jna2so4)*log_gamZ(jA,jna2so4) +   &
                    xmol(jnano3) *log_gamZ(jA,jnano3)  +   &
                    xmol(jnacl)  *log_gamZ(jA,jnacl)   +   &
                    xmol(jcano3) *log_gamZ(jA,jcano3)  +   &
                    xmol(jcacl2) *log_gamZ(jA,jcacl2)  +   &
                    xmol(jhno3)  *log_gamZ(jA,jhno3)   +   &
                    xmol(jhcl)   *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)
      activity(jna2so4,ibin) = mc(jc_na,ibin)**2 * ma(ja_so4,ibin) *   &
                               gam(jna2so4,ibin)**3
      endif


      jA = jnano3
      if(xmol(jA).gt.0.0)then
      log_gam(jA) = xmol(jnh4no3)*log_gamZ(jA,jnh4no3) +   &
                    xmol(jnh4cl) *log_gamZ(jA,jnh4cl)  +   &
                    xmol(jnh4so4)*log_gamZ(jA,jnh4so4) +   &
                    xmol(jna2so4)*log_gamZ(jA,jna2so4) +   &
                    xmol(jnano3) *log_gamZ(jA,jnano3)  +   &
                    xmol(jnacl)  *log_gamZ(jA,jnacl)   +   &
                    xmol(jcano3) *log_gamZ(jA,jcano3)  +   &
                    xmol(jcacl2) *log_gamZ(jA,jcacl2)  +   &
                    xmol(jhno3)  *log_gamZ(jA,jhno3)   +   &
                    xmol(jhcl)   *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)
      activity(jnano3,ibin)  = mc(jc_na,ibin) * ma(ja_no3,ibin) *   &
                               gam(jnano3,ibin)**2
      endif



      jA = jnacl
      if(xmol(jA).gt.0.0)then
      log_gam(jA) = xmol(jnh4no3)*log_gamZ(jA,jnh4no3) +   &
                    xmol(jnh4cl) *log_gamZ(jA,jnh4cl)  +   &
                    xmol(jnh4so4)*log_gamZ(jA,jnh4so4) +   &
                    xmol(jna2so4)*log_gamZ(jA,jna2so4) +   &
                    xmol(jnano3) *log_gamZ(jA,jnano3)  +   &
                    xmol(jnacl)  *log_gamZ(jA,jnacl)   +   &
                    xmol(jcano3) *log_gamZ(jA,jcano3)  +   &
                    xmol(jcacl2) *log_gamZ(jA,jcacl2)  +   &
                    xmol(jhno3)  *log_gamZ(jA,jhno3)   +   &
                    xmol(jhcl)   *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)
      activity(jnacl,ibin)   = mc(jc_na,ibin) * ma(ja_cl,ibin) *   &
                               gam(jnacl,ibin)**2
      endif



!c      jA = jcano3
!c      if(xmol(jA).gt.0.0)then
!c      gam(jA,ibin) = 1.0
!c      activity(jcano3,ibin)  = 1.0
!c      endif



!c      jA = jcacl2
!c      if(xmol(jA).gt.0.0)then
!c      gam(jA,ibin) = 1.0
!c      activity(jcacl2,ibin)  = 1.0
!c      endif

      jA = jcano3
      if(xmol(jA).gt.0.0)then
      log_gam(jA) = xmol(jnh4no3)*log_gamZ(jA,jnh4no3) +   &
                    xmol(jnh4cl) *log_gamZ(jA,jnh4cl)  +   &
                    xmol(jnh4so4)*log_gamZ(jA,jnh4so4) +   &
                    xmol(jna2so4)*log_gamZ(jA,jna2so4) +   &
                    xmol(jnano3) *log_gamZ(jA,jnano3)  +   &
                    xmol(jnacl)  *log_gamZ(jA,jnacl)   +   &
                    xmol(jcano3) *log_gamZ(jA,jcano3)  +   &
                    xmol(jcacl2) *log_gamZ(jA,jcacl2)  +   &
                    xmol(jhno3)  *log_gamZ(jA,jhno3)   +   &
                    xmol(jhcl)   *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)
      activity(jcano3,ibin)  = mc(jc_ca,ibin) * ma(ja_no3,ibin)**2 *   &
                               gam(jcano3,ibin)**3
      endif



      jA = jcacl2
      if(xmol(jA).gt.0.0)then
      log_gam(jA) = xmol(jnh4no3)*log_gamZ(jA,jnh4no3) +   &
                    xmol(jnh4cl) *log_gamZ(jA,jnh4cl)  +   &
                    xmol(jnh4so4)*log_gamZ(jA,jnh4so4) +   &
                    xmol(jna2so4)*log_gamZ(jA,jna2so4) +   &
                    xmol(jnano3) *log_gamZ(jA,jnano3)  +   &
                    xmol(jnacl)  *log_gamZ(jA,jnacl)   +   &
                    xmol(jcano3) *log_gamZ(jA,jcano3)  +   &
                    xmol(jcacl2) *log_gamZ(jA,jcacl2)  +   &
                    xmol(jhno3)  *log_gamZ(jA,jhno3)   +   &
                    xmol(jhcl)   *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)
      activity(jcacl2,ibin)  = mc(jc_ca,ibin) * ma(ja_cl,ibin)**2 *   &
                               gam(jcacl2,ibin)**3
      endif


      jA = jhno3
      log_gam(jA) = xmol(jnh4no3)*log_gamZ(jA,jnh4no3) +   &
                    xmol(jnh4cl) *log_gamZ(jA,jnh4cl)  +   &
                    xmol(jnh4so4)*log_gamZ(jA,jnh4so4) +   &
                    xmol(jna2so4)*log_gamZ(jA,jna2so4) +   &
                    xmol(jnano3) *log_gamZ(jA,jnano3)  +   &
                    xmol(jnacl)  *log_gamZ(jA,jnacl)   +   &
                    xmol(jcano3) *log_gamZ(jA,jcano3)  +   &
                    xmol(jcacl2) *log_gamZ(jA,jcacl2)  +   &
                    xmol(jhno3)  *log_gamZ(jA,jhno3)   +   &
                    xmol(jhcl)   *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)
      activity(jhno3,ibin)   = mc(jc_h,ibin) * ma(ja_no3,ibin) *   &
                               gam(jhno3,ibin)**2


      jA = jhcl
      log_gam(jA) = xmol(jnh4no3)*log_gamZ(jA,jnh4no3) +   &
                    xmol(jnh4cl) *log_gamZ(jA,jnh4cl)  +   &
                    xmol(jnh4so4)*log_gamZ(jA,jnh4so4) +   &
                    xmol(jna2so4)*log_gamZ(jA,jna2so4) +   &
                    xmol(jnano3) *log_gamZ(jA,jnano3)  +   &
                    xmol(jnacl)  *log_gamZ(jA,jnacl)   +   &
                    xmol(jcano3) *log_gamZ(jA,jcano3)  +   &
                    xmol(jcacl2) *log_gamZ(jA,jcacl2)  +   &
                    xmol(jhno3)  *log_gamZ(jA,jhno3)   +   &
                    xmol(jhcl)   *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)
      activity(jhcl,ibin)    = mc(jc_h,ibin) * ma(ja_cl,ibin) *   &
                               gam(jhcl,ibin)**2

!----
10    gam(jlvcite,ibin) = 1.0

      gam(jnh4hso4,ibin)= 1.0

      gam(jnh4msa,ibin) = 1.0

      gam(jna3hso4,ibin) = 1.0

      gam(jnahso4,ibin) = 1.0

      gam(jnamsa,ibin)  = 1.0

      gam(jcamsa2,ibin) = 1.0

      activity(jlvcite,ibin) = 0.0

      activity(jnh4hso4,ibin)= 0.0

      activity(jnh4msa,ibin) = mc(jc_nh4,ibin) * ma(ja_msa,ibin) *   &
                               gam(jnh4msa,ibin)**2

      activity(jna3hso4,ibin)= 0.0

      activity(jnahso4,ibin) = 0.0

      activity(jnamsa,ibin) = mc(jc_na,ibin) * ma(ja_msa,ibin) *   &
                               gam(jnamsa,ibin)**2

      activity(jcamsa2,ibin) = mc(jc_ca,ibin) * ma(ja_msa,ibin)**2 *   &
                               gam(jcamsa2,ibin)**3

      gam_ratio(ibin) = gam(jnh4no3,ibin)**2/gam(jhno3,ibin)**2


      else
!  SULFATE-RICH: solve for SO4= and HSO4- ions

      jp = jliquid

      sum_elec = 3.*electrolyte(jh2so4,jp,ibin)    +   &
                 2.*electrolyte(jnh4hso4,jp,ibin)  +   &
                 5.*electrolyte(jlvcite,jp,ibin)   +   &
                 3.*electrolyte(jnh4so4,jp,ibin)   +   &
                 2.*electrolyte(jnahso4,jp,ibin)   +   &
                 5.*electrolyte(jna3hso4,jp,ibin)  +   &
                 3.*electrolyte(jna2so4,jp,ibin)   +   &
                 2.*electrolyte(jhno3,jp,ibin)     +   &
                 2.*electrolyte(jhcl,jp,ibin)


      if(sum_elec .eq. 0.0)then
        do jA = 1, nelectrolyte
          gam(jA,ibin) = 1.0
        enddo
        goto 20
      endif


      xmol(jh2so4)  = 3.*electrolyte(jh2so4,jp,ibin)/sum_elec
      xmol(jnh4hso4)= 2.*electrolyte(jnh4hso4,jp,ibin)/sum_elec
      xmol(jlvcite) = 5.*electrolyte(jlvcite,jp,ibin)/sum_elec
      xmol(jnh4so4) = 3.*electrolyte(jnh4so4,jp,ibin)/sum_elec
      xmol(jnahso4) = 2.*electrolyte(jnahso4,jp,ibin)/sum_elec
      xmol(jna3hso4)= 5.*electrolyte(jna3hso4,jp,ibin)/sum_elec
      xmol(jna2so4) = 3.*electrolyte(jna2so4,jp,ibin)/sum_elec
      xmol(jhno3)   = 2.*electrolyte(jhno3,jp,ibin)/sum_elec
      xmol(jhcl)    = 2.*electrolyte(jhcl,jp,ibin)/sum_elec


! 2H.SO4
      jA = jh2so4
      log_gam(jA) = xmol(jh2so4)  *log_gamZ(jA,jh2so4)  +   &
                    xmol(jnh4hso4)*log_gamZ(jA,jnh4hso4)+   &
                    xmol(jlvcite) *log_gamZ(jA,jlvcite) +   &
                    xmol(jnh4so4) *log_gamZ(jA,jnh4so4) +   &
                    xmol(jnahso4) *log_gamZ(jA,jnahso4) +   &
                    xmol(jna3hso4)*log_gamZ(jA,jna3hso4)+   &
                    xmol(jna2so4) *log_gamZ(jA,jna2so4) +   &
                    xmol(jhno3)   *log_gamZ(jA,jhno3)   +   &
                    xmol(jhcl)    *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)


! H.HSO4
      jA = jhhso4
      log_gam(jA) = xmol(jh2so4)  *log_gamZ(jA,jh2so4)  +   &
                    xmol(jnh4hso4)*log_gamZ(jA,jnh4hso4)+   &
                    xmol(jlvcite) *log_gamZ(jA,jlvcite) +   &
                    xmol(jnh4so4) *log_gamZ(jA,jnh4so4) +   &
                    xmol(jnahso4) *log_gamZ(jA,jnahso4) +   &
                    xmol(jna3hso4)*log_gamZ(jA,jna3hso4)+   &
                    xmol(jna2so4) *log_gamZ(jA,jna2so4) +   &
                    xmol(jhno3)   *log_gamZ(jA,jhno3)   +   &
                    xmol(jhcl)    *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)


! NH4HSO4
      jA = jnh4hso4
      log_gam(jA) = xmol(jh2so4)  *log_gamZ(jA,jh2so4)  +   &
                    xmol(jnh4hso4)*log_gamZ(jA,jnh4hso4)+   &
                    xmol(jlvcite) *log_gamZ(jA,jlvcite) +   &
                    xmol(jnh4so4) *log_gamZ(jA,jnh4so4) +   &
                    xmol(jnahso4) *log_gamZ(jA,jnahso4) +   &
                    xmol(jna3hso4)*log_gamZ(jA,jna3hso4)+   &
                    xmol(jna2so4) *log_gamZ(jA,jna2so4) +   &
                    xmol(jhno3)   *log_gamZ(jA,jhno3)   +   &
                    xmol(jhcl)    *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)


! LETOVICITE
      jA = jlvcite
      log_gam(jA) = xmol(jh2so4)  *log_gamZ(jA,jh2so4)  +   &
                    xmol(jnh4hso4)*log_gamZ(jA,jnh4hso4)+   &
                    xmol(jlvcite) *log_gamZ(jA,jlvcite) +   &
                    xmol(jnh4so4) *log_gamZ(jA,jnh4so4) +   &
                    xmol(jnahso4) *log_gamZ(jA,jnahso4) +   &
                    xmol(jna3hso4)*log_gamZ(jA,jna3hso4)+   &
                    xmol(jna2so4) *log_gamZ(jA,jna2so4) +   &
                    xmol(jhno3)   *log_gamZ(jA,jhno3)   +   &
                    xmol(jhcl)    *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)


! (NH4)2SO4
      jA = jnh4so4
      log_gam(jA) = xmol(jh2so4)  *log_gamZ(jA,jh2so4)  +   &
                    xmol(jnh4hso4)*log_gamZ(jA,jnh4hso4)+   &
                    xmol(jlvcite) *log_gamZ(jA,jlvcite) +   &
                    xmol(jnh4so4) *log_gamZ(jA,jnh4so4) +   &
                    xmol(jnahso4) *log_gamZ(jA,jnahso4) +   &
                    xmol(jna3hso4)*log_gamZ(jA,jna3hso4)+   &
                    xmol(jna2so4) *log_gamZ(jA,jna2so4) +   &
                    xmol(jhno3)   *log_gamZ(jA,jhno3)   +   &
                    xmol(jhcl)    *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)


! NaHSO4
      jA = jnahso4
      log_gam(jA) = xmol(jh2so4)  *log_gamZ(jA,jh2so4)  +   &
                    xmol(jnh4hso4)*log_gamZ(jA,jnh4hso4)+   &
                    xmol(jlvcite) *log_gamZ(jA,jlvcite) +   &
                    xmol(jnh4so4) *log_gamZ(jA,jnh4so4) +   &
                    xmol(jnahso4) *log_gamZ(jA,jnahso4) +   &
                    xmol(jna3hso4)*log_gamZ(jA,jna3hso4)+   &
                    xmol(jna2so4) *log_gamZ(jA,jna2so4) +   &
                    xmol(jhno3)   *log_gamZ(jA,jhno3)   +   &
                    xmol(jhcl)    *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)


! Na3H(SO4)2
      jA = jna3hso4
!      log_gam(jA) = xmol(jh2so4)  *log_gamZ(jA,jh2so4)  +
!     &              xmol(jnh4hso4)*log_gamZ(jA,jnh4hso4)+
!     &              xmol(jlvcite) *log_gamZ(jA,jlvcite) +
!     &              xmol(jnh4so4) *log_gamZ(jA,jnh4so4) +
!     &              xmol(jnahso4) *log_gamZ(jA,jnahso4) +
!     &              xmol(jna3hso4)*log_gamZ(jA,jna3hso4)+
!     &              xmol(jna2so4) *log_gamZ(jA,jna2so4) +
!     &              xmol(jhno3)   *log_gamZ(jA,jhno3)   +
!     &              xmol(jhcl)    *log_gamZ(jA,jhcl)
!      gam(jA,ibin) = 10.**log_gam(jA)
      gam(jA,ibin) = 1.0


! Na2SO4
      jA = jna2so4
      log_gam(jA) = xmol(jh2so4)  *log_gamZ(jA,jh2so4)  +   &
                    xmol(jnh4hso4)*log_gamZ(jA,jnh4hso4)+   &
                    xmol(jlvcite) *log_gamZ(jA,jlvcite) +   &
                    xmol(jnh4so4) *log_gamZ(jA,jnh4so4) +   &
                    xmol(jnahso4) *log_gamZ(jA,jnahso4) +   &
                    xmol(jna3hso4)*log_gamZ(jA,jna3hso4)+   &
                    xmol(jna2so4) *log_gamZ(jA,jna2so4) +   &
                    xmol(jhno3)   *log_gamZ(jA,jhno3)   +   &
                    xmol(jhcl)    *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)


! HNO3
      jA = jhno3
      log_gam(jA) = xmol(jh2so4)  *log_gamZ(jA,jh2so4)  +   &
                    xmol(jnh4hso4)*log_gamZ(jA,jnh4hso4)+   &
                    xmol(jlvcite) *log_gamZ(jA,jlvcite) +   &
                    xmol(jnh4so4) *log_gamZ(jA,jnh4so4) +   &
                    xmol(jnahso4) *log_gamZ(jA,jnahso4) +   &
                    xmol(jna3hso4)*log_gamZ(jA,jna3hso4)+   &
                    xmol(jna2so4) *log_gamZ(jA,jna2so4) +   &
                    xmol(jhno3)   *log_gamZ(jA,jhno3)   +   &
                    xmol(jhcl)    *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)


! HCl
      jA = jhcl
      log_gam(jA) = xmol(jh2so4)  *log_gamZ(jA,jh2so4)  +   &
                    xmol(jnh4hso4)*log_gamZ(jA,jnh4hso4)+   &
                    xmol(jlvcite) *log_gamZ(jA,jlvcite) +   &
                    xmol(jnh4so4) *log_gamZ(jA,jnh4so4) +   &
                    xmol(jnahso4) *log_gamZ(jA,jnahso4) +   &
                    xmol(jna3hso4)*log_gamZ(jA,jna3hso4)+   &
                    xmol(jna2so4) *log_gamZ(jA,jna2so4) +   &
                    xmol(jhno3)   *log_gamZ(jA,jhno3)   +   &
                    xmol(jhcl)    *log_gamZ(jA,jhcl)
      gam(jA,ibin) = 10.**log_gam(jA)


20    gam(jnh4no3,ibin) = 1.0
      gam(jnh4cl,ibin)  = 1.0
      gam(jnano3,ibin)  = 1.0
      gam(jnacl,ibin)   = 1.0
      gam(jcano3,ibin)  = 1.0
      gam(jcacl2,ibin)  = 1.0

      gam(jnh4msa,ibin) = 1.0
      gam(jnamsa,ibin)  = 1.0
      gam(jcamsa2,ibin) = 1.0



! compute equilibrium pH
! cation molalities (mol/kg water)
      mc(jc_ca,ibin)   = 1.e-9*aer(ica_a,jliquid,ibin)/water_a(ibin)
      mc(jc_nh4,ibin)  = 1.e-9*aer(inh4_a,jliquid,ibin)/water_a(ibin)
      mc(jc_na,ibin)   = 1.e-9*aer(ina_a, jliquid,ibin)/water_a(ibin)

! anion molalities (mol/kg water)
      mSULF            = 1.e-9*aer(iso4_a,jliquid,ibin)/water_a(ibin)
      ma(ja_hso4,ibin) = 0.0
      ma(ja_so4,ibin)  = 0.0
      ma(ja_no3,ibin)  = 1.e-9*aer(ino3_a,jliquid,ibin)/water_a(ibin)
      ma(ja_cl,ibin)   = 1.e-9*aer(icl_a, jliquid,ibin)/water_a(ibin)
      ma(ja_msa,ibin)  = 1.e-9*aer(imsa_a,jliquid,ibin)/water_a(ibin)

      gam_ratio(ibin)  = gam(jnh4hso4,ibin)**2/gam(jhhso4,ibin)**2
      dumK = Keq_ll(1)*gam(jhhso4,ibin)**2/gam(jh2so4,ibin)**3

      c_bal =  mc(jc_nh4,ibin) + mc(jc_na,ibin) + 2.*mc(jc_ca,ibin)   &
         - ma(ja_no3,ibin) - ma(ja_cl,ibin) - mSULF - ma(ja_msa,ibin)

      aq = 1.0
      bq = dumK + c_bal
      cq = dumK*(c_bal - mSULF)


!--quadratic solution
        if(bq .ne. 0.0)then
        xq = 4.*(1./bq)*(cq/bq)
        else
        xq = 1.e+6
        endif

        if(abs(xq) .lt. 1.e-6)then
          dum = xq*(0.5 + xq*(0.125 + xq*0.0625))
          quad = (-0.5*bq/aq)*dum
          if(quad .lt. 0.)then
            quad = -bq/aq - quad
          endif
        else
          quad = 0.5*(-bq+sqrt(bq*bq - 4.*cq))
        endif
!--end of quadratic solution

      mc(jc_h,ibin) = max(quad, 1.d-7)
      ma(ja_so4,ibin) = mSULF*dumK/(mc(jc_h,ibin) + dumK)
      ma(ja_hso4,ibin)= mSULF - ma(ja_so4,ibin)

      activity(jcamsa2,ibin) = mc(jc_ca,ibin) * ma(ja_msa,ibin)**2 *   &
                               gam(jcamsa2,ibin)**3

      activity(jnh4so4,ibin) = mc(jc_nh4,ibin)**2 * ma(ja_so4,ibin) *   &
                               gam(jnh4so4,ibin)**3

      activity(jlvcite,ibin) = mc(jc_nh4,ibin)**3 * ma(ja_hso4,ibin) *   &
                               ma(ja_so4,ibin) * gam(jlvcite,ibin)**5

      activity(jnh4hso4,ibin)= mc(jc_nh4,ibin) * ma(ja_hso4,ibin) *   &
                               gam(jnh4hso4,ibin)**2

      activity(jnh4msa,ibin) = mc(jc_nh4,ibin) * ma(ja_msa,ibin) *   &
                               gam(jnh4msa,ibin)**2

      activity(jna2so4,ibin) = mc(jc_na,ibin)**2 * ma(ja_so4,ibin) *   &
                               gam(jna2so4,ibin)**3

      activity(jnahso4,ibin) = mc(jc_na,ibin) * ma(ja_hso4,ibin) *   &
                               gam(jnahso4,ibin)**2

      activity(jnamsa,ibin)  = mc(jc_na,ibin) * ma(ja_msa,ibin) *   &
                               gam(jnamsa,ibin)**2

!      activity(jna3hso4,ibin)= mc(jc_na,ibin)**3 * ma(ja_hso4,ibin) *
!     &                         ma(ja_so4,ibin) * gam(jna3hso4,ibin)**5

      activity(jna3hso4,ibin)= 0.0

      activity(jhno3,ibin)   = mc(jc_h,ibin) * ma(ja_no3,ibin) *   &
                               gam(jhno3,ibin)**2

      activity(jhcl,ibin)    = mc(jc_h,ibin) * ma(ja_cl,ibin) *   &
                               gam(jhcl,ibin)**2

      activity(jmsa,ibin)    = mc(jc_h,ibin) * ma(ja_msa,ibin) *   &
                               gam(jmsa,ibin)**2


! sulfate-poor species
      activity(jnh4no3,ibin) = 0.0

      activity(jnh4cl,ibin)  = 0.0

      activity(jnano3,ibin)  = 0.0

      activity(jnacl,ibin)   = 0.0

      activity(jcano3,ibin)  = 0.0

      activity(jcacl2,ibin)  = 0.0


      endif




      return
      end subroutine compute_activities












!***********************************************************************
! computes MTEM ternary parameters only once per transport time-step
! for a given aH2O (= RH)
!
! author: Rahul A. Zaveri
! update: jan 2005
! reference: Zaveri, R.A., R.C. Easter, and A.S. Wexler,
! A new method for multicomponent activity coefficients of electrolytes
! in aqueous atmospheric aerosols, J. Geophys. Res., 2005.
!-----------------------------------------------------------------------
      subroutine MTEM_compute_log_gamZ
      use module_data_mosaic_aero

      implicit none

! local variables
      integer jA
! functions
      real(r8) :: fnlog_gamZ, bin_molality


! sulfate-poor species
      jA = jhno3
      log_gamZ(jA,jnh4so4) = fnlog_gamZ(jA,jnh4so4)
      log_gamZ(jA,jnh4no3) = fnlog_gamZ(jA,jnh4no3)
      log_gamZ(jA,jnh4cl)  = fnlog_gamZ(jA,jnh4cl)
      log_gamZ(jA,jna2so4) = fnlog_gamZ(jA,jna2so4)
      log_gamZ(jA,jnano3)  = fnlog_gamZ(jA,jnano3)
      log_gamZ(jA,jnacl)   = fnlog_gamZ(jA,jnacl)
      log_gamZ(jA,jcano3)  = fnlog_gamZ(jA,jcano3)
      log_gamZ(jA,jcacl2)  = fnlog_gamZ(jA,jcacl2)
      log_gamZ(jA,jhno3)   = fnlog_gamZ(jA,jhno3)
      log_gamZ(jA,jhcl)    = fnlog_gamZ(jA,jhcl)
      log_gamZ(jA,jh2so4)  = fnlog_gamZ(jA,jh2so4)
      log_gamZ(jA,jnh4hso4)= fnlog_gamZ(jA,jnh4hso4)
      log_gamZ(jA,jlvcite) = fnlog_gamZ(jA,jlvcite)
      log_gamZ(jA,jnahso4) = fnlog_gamZ(jA,jnahso4)
      log_gamZ(jA,jna3hso4)= fnlog_gamZ(jA,jna3hso4)


      jA = jhcl
      log_gamZ(jA,jnh4so4) = fnlog_gamZ(jA,jnh4so4)
      log_gamZ(jA,jnh4no3) = fnlog_gamZ(jA,jnh4no3)
      log_gamZ(jA,jnh4cl)  = fnlog_gamZ(jA,jnh4cl)
      log_gamZ(jA,jna2so4) = fnlog_gamZ(jA,jna2so4)
      log_gamZ(jA,jnano3)  = fnlog_gamZ(jA,jnano3)
      log_gamZ(jA,jnacl)   = fnlog_gamZ(jA,jnacl)
      log_gamZ(jA,jcano3)  = fnlog_gamZ(jA,jcano3)
      log_gamZ(jA,jcacl2)  = fnlog_gamZ(jA,jcacl2)
      log_gamZ(jA,jhno3)   = fnlog_gamZ(jA,jhno3)
      log_gamZ(jA,jhcl)    = fnlog_gamZ(jA,jhcl)
      log_gamZ(jA,jh2so4)  = fnlog_gamZ(jA,jh2so4)
      log_gamZ(jA,jnh4hso4)= fnlog_gamZ(jA,jnh4hso4)
      log_gamZ(jA,jlvcite) = fnlog_gamZ(jA,jlvcite)
      log_gamZ(jA,jnahso4) = fnlog_gamZ(jA,jnahso4)
      log_gamZ(jA,jna3hso4)= fnlog_gamZ(jA,jna3hso4)


      jA = jnh4so4
      log_gamZ(jA,jnh4so4) = fnlog_gamZ(jA,jnh4so4)
      log_gamZ(jA,jnh4no3) = fnlog_gamZ(jA,jnh4no3)
      log_gamZ(jA,jnh4cl)  = fnlog_gamZ(jA,jnh4cl)
      log_gamZ(jA,jna2so4) = fnlog_gamZ(jA,jna2so4)
      log_gamZ(jA,jnano3)  = fnlog_gamZ(jA,jnano3)
      log_gamZ(jA,jnacl)   = fnlog_gamZ(jA,jnacl)
      log_gamZ(jA,jcano3)  = fnlog_gamZ(jA,jcano3)
      log_gamZ(jA,jcacl2)  = fnlog_gamZ(jA,jcacl2)
      log_gamZ(jA,jhno3)   = fnlog_gamZ(jA,jhno3)
      log_gamZ(jA,jhcl)    = fnlog_gamZ(jA,jhcl)
      log_gamZ(jA,jh2so4)  = fnlog_gamZ(jA,jh2so4)
      log_gamZ(jA,jnh4hso4)= fnlog_gamZ(jA,jnh4hso4)
      log_gamZ(jA,jlvcite) = fnlog_gamZ(jA,jlvcite)
      log_gamZ(jA,jnahso4) = fnlog_gamZ(jA,jnahso4)
      log_gamZ(jA,jna3hso4)= fnlog_gamZ(jA,jna3hso4)


      jA = jnh4no3
      log_gamZ(jA,jnh4so4) = fnlog_gamZ(jA,jnh4so4)
      log_gamZ(jA,jnh4no3) = fnlog_gamZ(jA,jnh4no3)
      log_gamZ(jA,jnh4cl)  = fnlog_gamZ(jA,jnh4cl)
      log_gamZ(jA,jna2so4) = fnlog_gamZ(jA,jna2so4)
      log_gamZ(jA,jnano3)  = fnlog_gamZ(jA,jnano3)
      log_gamZ(jA,jnacl)   = fnlog_gamZ(jA,jnacl)
      log_gamZ(jA,jcano3)  = fnlog_gamZ(jA,jcano3)
      log_gamZ(jA,jcacl2)  = fnlog_gamZ(jA,jcacl2)
      log_gamZ(jA,jhno3)   = fnlog_gamZ(jA,jhno3)
      log_gamZ(jA,jhcl)    = fnlog_gamZ(jA,jhcl)


      jA = jnh4cl
      log_gamZ(jA,jnh4so4) = fnlog_gamZ(jA,jnh4so4)
      log_gamZ(jA,jnh4no3) = fnlog_gamZ(jA,jnh4no3)
      log_gamZ(jA,jnh4cl)  = fnlog_gamZ(jA,jnh4cl)
      log_gamZ(jA,jna2so4) = fnlog_gamZ(jA,jna2so4)
      log_gamZ(jA,jnano3)  = fnlog_gamZ(jA,jnano3)
      log_gamZ(jA,jnacl)   = fnlog_gamZ(jA,jnacl)
      log_gamZ(jA,jcano3)  = fnlog_gamZ(jA,jcano3)
      log_gamZ(jA,jcacl2)  = fnlog_gamZ(jA,jcacl2)
      log_gamZ(jA,jhno3)   = fnlog_gamZ(jA,jhno3)
      log_gamZ(jA,jhcl)    = fnlog_gamZ(jA,jhcl)


      jA = jna2so4
      log_gamZ(jA,jnh4so4) = fnlog_gamZ(jA,jnh4so4)
      log_gamZ(jA,jnh4no3) = fnlog_gamZ(jA,jnh4no3)
      log_gamZ(jA,jnh4cl)  = fnlog_gamZ(jA,jnh4cl)
      log_gamZ(jA,jna2so4) = fnlog_gamZ(jA,jna2so4)
      log_gamZ(jA,jnano3)  = fnlog_gamZ(jA,jnano3)
      log_gamZ(jA,jnacl)   = fnlog_gamZ(jA,jnacl)
      log_gamZ(jA,jcano3)  = fnlog_gamZ(jA,jcano3)
      log_gamZ(jA,jcacl2)  = fnlog_gamZ(jA,jcacl2)
      log_gamZ(jA,jhno3)   = fnlog_gamZ(jA,jhno3)
      log_gamZ(jA,jhcl)    = fnlog_gamZ(jA,jhcl)
      log_gamZ(jA,jh2so4)  = fnlog_gamZ(jA,jh2so4)
      log_gamZ(jA,jnh4hso4)= fnlog_gamZ(jA,jnh4hso4)
      log_gamZ(jA,jlvcite) = fnlog_gamZ(jA,jlvcite)
      log_gamZ(jA,jnahso4) = fnlog_gamZ(jA,jnahso4)
      log_gamZ(jA,jna3hso4)= fnlog_gamZ(jA,jna3hso4)


      jA = jnano3
      log_gamZ(jA,jnh4so4) = fnlog_gamZ(jA,jnh4so4)
      log_gamZ(jA,jnh4no3) = fnlog_gamZ(jA,jnh4no3)
      log_gamZ(jA,jnh4cl)  = fnlog_gamZ(jA,jnh4cl)
      log_gamZ(jA,jna2so4) = fnlog_gamZ(jA,jna2so4)
      log_gamZ(jA,jnano3)  = fnlog_gamZ(jA,jnano3)
      log_gamZ(jA,jnacl)   = fnlog_gamZ(jA,jnacl)
      log_gamZ(jA,jcano3)  = fnlog_gamZ(jA,jcano3)
      log_gamZ(jA,jcacl2)  = fnlog_gamZ(jA,jcacl2)
      log_gamZ(jA,jhno3)   = fnlog_gamZ(jA,jhno3)
      log_gamZ(jA,jhcl)    = fnlog_gamZ(jA,jhcl)


      jA = jnacl
      log_gamZ(jA,jnh4so4) = fnlog_gamZ(jA,jnh4so4)
      log_gamZ(jA,jnh4no3) = fnlog_gamZ(jA,jnh4no3)
      log_gamZ(jA,jnh4cl)  = fnlog_gamZ(jA,jnh4cl)
      log_gamZ(jA,jna2so4) = fnlog_gamZ(jA,jna2so4)
      log_gamZ(jA,jnano3)  = fnlog_gamZ(jA,jnano3)
      log_gamZ(jA,jnacl)   = fnlog_gamZ(jA,jnacl)
      log_gamZ(jA,jcano3)  = fnlog_gamZ(jA,jcano3)
      log_gamZ(jA,jcacl2)  = fnlog_gamZ(jA,jcacl2)
      log_gamZ(jA,jhno3)   = fnlog_gamZ(jA,jhno3)
      log_gamZ(jA,jhcl)    = fnlog_gamZ(jA,jhcl)


      jA = jcano3
      log_gamZ(jA,jnh4so4) = fnlog_gamZ(jA,jnh4so4)
      log_gamZ(jA,jnh4no3) = fnlog_gamZ(jA,jnh4no3)
      log_gamZ(jA,jnh4cl)  = fnlog_gamZ(jA,jnh4cl)
      log_gamZ(jA,jna2so4) = fnlog_gamZ(jA,jna2so4)
      log_gamZ(jA,jnano3)  = fnlog_gamZ(jA,jnano3)
      log_gamZ(jA,jnacl)   = fnlog_gamZ(jA,jnacl)
      log_gamZ(jA,jcano3)  = fnlog_gamZ(jA,jcano3)
      log_gamZ(jA,jcacl2)  = fnlog_gamZ(jA,jcacl2)
      log_gamZ(jA,jhno3)   = fnlog_gamZ(jA,jhno3)
      log_gamZ(jA,jhcl)    = fnlog_gamZ(jA,jhcl)


      jA = jcacl2
      log_gamZ(jA,jnh4so4) = fnlog_gamZ(jA,jnh4so4)
      log_gamZ(jA,jnh4no3) = fnlog_gamZ(jA,jnh4no3)
      log_gamZ(jA,jnh4cl)  = fnlog_gamZ(jA,jnh4cl)
      log_gamZ(jA,jna2so4) = fnlog_gamZ(jA,jna2so4)
      log_gamZ(jA,jnano3)  = fnlog_gamZ(jA,jnano3)
      log_gamZ(jA,jnacl)   = fnlog_gamZ(jA,jnacl)
      log_gamZ(jA,jcano3)  = fnlog_gamZ(jA,jcano3)
      log_gamZ(jA,jcacl2)  = fnlog_gamZ(jA,jcacl2)
      log_gamZ(jA,jhno3)   = fnlog_gamZ(jA,jhno3)
      log_gamZ(jA,jhcl)    = fnlog_gamZ(jA,jhcl)


! sulfate-rich species
      jA = jh2so4
      log_gamZ(jA,jh2so4)  = fnlog_gamZ(jA,jh2so4)
      log_gamZ(jA,jnh4hso4)= fnlog_gamZ(jA,jnh4hso4)
      log_gamZ(jA,jlvcite) = fnlog_gamZ(jA,jlvcite)
      log_gamZ(jA,jnh4so4) = fnlog_gamZ(jA,jnh4so4)
      log_gamZ(jA,jnahso4) = fnlog_gamZ(jA,jnahso4)
      log_gamZ(jA,jna3hso4)= fnlog_gamZ(jA,jna3hso4)
      log_gamZ(jA,jna2so4) = fnlog_gamZ(jA,jna2so4)
      log_gamZ(jA,jhno3)   = fnlog_gamZ(jA,jhno3)
      log_gamZ(jA,jhcl)    = fnlog_gamZ(jA,jhcl)


      jA = jhhso4
      log_gamZ(jA,jh2so4)  = fnlog_gamZ(jA,jh2so4)
      log_gamZ(jA,jnh4hso4)= fnlog_gamZ(jA,jnh4hso4)
      log_gamZ(jA,jlvcite) = fnlog_gamZ(jA,jlvcite)
      log_gamZ(jA,jnh4so4) = fnlog_gamZ(jA,jnh4so4)
      log_gamZ(jA,jnahso4) = fnlog_gamZ(jA,jnahso4)
      log_gamZ(jA,jna3hso4)= fnlog_gamZ(jA,jna3hso4)
      log_gamZ(jA,jna2so4) = fnlog_gamZ(jA,jna2so4)
      log_gamZ(jA,jhno3)   = fnlog_gamZ(jA,jhno3)
      log_gamZ(jA,jhcl)    = fnlog_gamZ(jA,jhcl)


      jA = jnh4hso4
      log_gamZ(jA,jh2so4)  = fnlog_gamZ(jA,jh2so4)
      log_gamZ(jA,jnh4hso4)= fnlog_gamZ(jA,jnh4hso4)
      log_gamZ(jA,jlvcite) = fnlog_gamZ(jA,jlvcite)
      log_gamZ(jA,jnh4so4) = fnlog_gamZ(jA,jnh4so4)
      log_gamZ(jA,jnahso4) = fnlog_gamZ(jA,jnahso4)
      log_gamZ(jA,jna3hso4)= fnlog_gamZ(jA,jna3hso4)
      log_gamZ(jA,jna2so4) = fnlog_gamZ(jA,jna2so4)
      log_gamZ(jA,jhno3)   = fnlog_gamZ(jA,jhno3)
      log_gamZ(jA,jhcl)    = fnlog_gamZ(jA,jhcl)


      jA = jlvcite
      log_gamZ(jA,jh2so4)  = fnlog_gamZ(jA,jh2so4)
      log_gamZ(jA,jnh4hso4)= fnlog_gamZ(jA,jnh4hso4)
      log_gamZ(jA,jlvcite) = fnlog_gamZ(jA,jlvcite)
      log_gamZ(jA,jnh4so4) = fnlog_gamZ(jA,jnh4so4)
      log_gamZ(jA,jnahso4) = fnlog_gamZ(jA,jnahso4)
      log_gamZ(jA,jna3hso4)= fnlog_gamZ(jA,jna3hso4)
      log_gamZ(jA,jna2so4) = fnlog_gamZ(jA,jna2so4)
      log_gamZ(jA,jhno3)   = fnlog_gamZ(jA,jhno3)
      log_gamZ(jA,jhcl)    = fnlog_gamZ(jA,jhcl)


      jA = jnahso4
      log_gamZ(jA,jh2so4)  = fnlog_gamZ(jA,jh2so4)
      log_gamZ(jA,jnh4hso4)= fnlog_gamZ(jA,jnh4hso4)
      log_gamZ(jA,jlvcite) = fnlog_gamZ(jA,jlvcite)
      log_gamZ(jA,jnh4so4) = fnlog_gamZ(jA,jnh4so4)
      log_gamZ(jA,jnahso4) = fnlog_gamZ(jA,jnahso4)
      log_gamZ(jA,jna3hso4)= fnlog_gamZ(jA,jna3hso4)
      log_gamZ(jA,jna2so4) = fnlog_gamZ(jA,jna2so4)
      log_gamZ(jA,jhno3)   = fnlog_gamZ(jA,jhno3)
      log_gamZ(jA,jhcl)    = fnlog_gamZ(jA,jhcl)


      jA = jna3hso4
      log_gamZ(jA,jh2so4)  = fnlog_gamZ(jA,jh2so4)
      log_gamZ(jA,jnh4hso4)= fnlog_gamZ(jA,jnh4hso4)
      log_gamZ(jA,jlvcite) = fnlog_gamZ(jA,jlvcite)
      log_gamZ(jA,jnh4so4) = fnlog_gamZ(jA,jnh4so4)
      log_gamZ(jA,jnahso4) = fnlog_gamZ(jA,jnahso4)
      log_gamZ(jA,jna3hso4)= fnlog_gamZ(jA,jna3hso4)
      log_gamZ(jA,jna2so4) = fnlog_gamZ(jA,jna2so4)
      log_gamZ(jA,jhno3)   = fnlog_gamZ(jA,jhno3)
      log_gamZ(jA,jhcl)    = fnlog_gamZ(jA,jhcl)

      return
      end subroutine MTEM_compute_log_gamZ




























!***********************************************************************
! computes sulfate ratio
!
! author: Rahul A. Zaveri
! update: dec 1999
!-----------------------------------------------------------------------
      subroutine calculate_XT(ibin,jp,XT)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin, jp
      real(r8) :: XT


      if( (aer(iso4_a,jp,ibin)+aer(imsa_a,jp,ibin)) .gt.0.0)then
        XT   = ( aer(inh4_a,jp,ibin) +   &
                 aer(ina_a,jp,ibin)  +   &
              2.*aer(ica_a,jp,ibin) )/   &
               (aer(iso4_a,jp,ibin)+0.5*aer(imsa_a,jp,ibin))
      else
        XT   = -1.0
      endif


      return
      end subroutine calculate_XT





!***********************************************************************
! computes ions from electrolytes
!
! author: Rahul A. Zaveri
! update: jan 2005
!-----------------------------------------------------------------------
      subroutine electrolytes_to_ions(jp,ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer jp, ibin
! local variables
      real(r8) :: sum_dum


      aer(iso4_a,jp,ibin) = electrolyte(jcaso4,jp,ibin)  +   &
                            electrolyte(jna2so4,jp,ibin) +   &
                         2.*electrolyte(jna3hso4,jp,ibin)+   &
                            electrolyte(jnahso4,jp,ibin) +   &
                            electrolyte(jnh4so4,jp,ibin) +   &
                         2.*electrolyte(jlvcite,jp,ibin) +   &
                            electrolyte(jnh4hso4,jp,ibin)+   &
                            electrolyte(jh2so4,jp,ibin)

      aer(ino3_a,jp,ibin) = electrolyte(jnano3,jp,ibin)  +   &
                         2.*electrolyte(jcano3,jp,ibin)  +   &
                            electrolyte(jnh4no3,jp,ibin) +   &
                            electrolyte(jhno3,jp,ibin)

      aer(icl_a,jp,ibin)  = electrolyte(jnacl,jp,ibin)   +   &
                         2.*electrolyte(jcacl2,jp,ibin)  +   &
                            electrolyte(jnh4cl,jp,ibin)  +   &
                            electrolyte(jhcl,jp,ibin)

      aer(imsa_a,jp,ibin) = electrolyte(jnh4msa,jp,ibin) +   &
                            electrolyte(jnamsa,jp,ibin)  +   &
                         2.*electrolyte(jcamsa2,jp,ibin) +   &
                            electrolyte(jmsa,jp,ibin)

      aer(ico3_a,jp,ibin) = electrolyte(jcaco3,jp,ibin)

      aer(ica_a,jp,ibin)  = electrolyte(jcaso4,jp,ibin)  +   &
                            electrolyte(jcano3,jp,ibin)  +   &
                            electrolyte(jcacl2,jp,ibin)  +   &
                            electrolyte(jcaco3,jp,ibin)  +   &
                            electrolyte(jcamsa2,jp,ibin)

      aer(ina_a,jp,ibin)  = electrolyte(jnano3,jp,ibin)  +   &
                            electrolyte(jnacl,jp,ibin)   +   &
                         2.*electrolyte(jna2so4,jp,ibin) +   &
                         3.*electrolyte(jna3hso4,jp,ibin)+   &
                            electrolyte(jnahso4,jp,ibin) +   &
                            electrolyte(jnamsa,jp,ibin)

      aer(inh4_a,jp,ibin) = electrolyte(jnh4no3,jp,ibin) +   &
                            electrolyte(jnh4cl,jp,ibin)  +   &
                         2.*electrolyte(jnh4so4,jp,ibin) +   &
                         3.*electrolyte(jlvcite,jp,ibin) +   &
                            electrolyte(jnh4hso4,jp,ibin)+   &
                            electrolyte(jnh4msa,jp,ibin)


      return
      end subroutine electrolytes_to_ions










!***********************************************************************
! combinatorial method for computing electrolytes from ions
!
! notes:
!  - to be used for liquid-phase or total-phase only
!  - transfers caso4 and caco3 from liquid to solid phase
!
! author: Rahul A. Zaveri (based on code provided by A.S. Wexler)
! update: apr 2005
!-----------------------------------------------------------------------
      subroutine ions_to_electrolytes(jp,ibin,XT)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin, jp
      real(r8) :: XT
! local variables
      integer iaer, je, jc, ja, icase
      real(r8) :: store(naer), sum_dum, sum_naza, sum_nczc, sum_na_nh4,   &
           f_nh4, f_na, xh, xb, xl, xs, cat_net, rem_nh4, rem_na
      real(r8) :: nc(ncation), na(nanion)




      if(jp .ne. jliquid)then
        write(6,*)' jp must be jliquid'
        write(6,*)' in ions_to_electrolytes sub'
        write(6,*)' wrong jp = ', jp
        stop
      endif

! remove negative concentrations, if any
!      do iaer = 1, naer
!        aer(iaer,jp,ibin) = max(0.0d0, aer(iaer,jp,ibin))	! EFFI
!      enddo


! first transfer caso4 from liquid to solid phase (caco3 should not be present here)
      store(ica_a)  = aer(ica_a, jp,ibin)
      store(iso4_a) = aer(iso4_a,jp,ibin)

      call form_caso4(store,jp,ibin)

      if(jp .eq. jliquid)then ! transfer caso4 from liquid to solid phase
        aer(ica_a,jliquid,ibin) = aer(ica_a,jliquid,ibin) -   &
                                  electrolyte(jcaso4,jliquid,ibin)

        aer(iso4_a,jliquid,ibin)= aer(iso4_a,jliquid,ibin)-   &
                                  electrolyte(jcaso4,jliquid,ibin)

        aer(ica_a,jsolid,ibin)  = aer(ica_a,jsolid,ibin) +   &
                                  electrolyte(jcaso4,jliquid,ibin)

        aer(iso4_a,jsolid,ibin) = aer(iso4_a,jsolid,ibin) +   &
                                  electrolyte(jcaso4,jliquid,ibin)

        electrolyte(jcaso4,jsolid,ibin)=electrolyte(jcaso4,jsolid,ibin)   &
                                       +electrolyte(jcaso4,jliquid,ibin)
        electrolyte(jcaso4,jliquid,ibin)= 0.0
      endif


! calculate sulfate ratio
!      call calculate_XT(ibin,jp,XT)		! EFFI

      if( (aer(iso4_a,jp,ibin)+aer(imsa_a,jp,ibin)) .gt.0.0)then
        XT   = ( aer(inh4_a,jp,ibin) +   &
                 aer(ina_a,jp,ibin)  +   &
              2.*aer(ica_a,jp,ibin) )/   &
               (aer(iso4_a,jp,ibin)+0.5*aer(imsa_a,jp,ibin))
      else
        XT   = -1.0
      endif




      if(XT .ge. 1.9999 .or. XT.lt.0.)then
       icase = 1	! near neutral (acidity is caused by HCl and/or HNO3)
      else
       icase = 2	! acidic (acidity is caused by excess SO4)
      endif


! initialize to zero
      do je = 1, nelectrolyte
        electrolyte(je,jp,ibin) = 0.0
      enddo

!
!---------------------------------------------------------
! initialize moles of ions depending on the sulfate domain

      if(icase.eq.1)then ! XT >= 2 : SULFATE POOR DOMAIN

        na(ja_hso4)= 0.0
        na(ja_so4) = aer(iso4_a,jp,ibin)
        na(ja_no3) = aer(ino3_a,jp,ibin)
        na(ja_cl)  = aer(icl_a, jp,ibin)
        na(ja_msa) = aer(imsa_a,jp,ibin)

        nc(jc_ca)  = aer(ica_a, jp,ibin)
        nc(jc_na)  = aer(ina_a, jp,ibin)
        nc(jc_nh4) = aer(inh4_a,jp,ibin)

        cat_net = (   &
            (2.*na(ja_so4)+na(ja_no3)+na(ja_cl)+na(ja_msa)) -   &
            (2.*nc(jc_ca) +nc(jc_nh4)+nc(jc_na)) )

        if(cat_net .lt. 0.0)then

          nc(jc_h) = 0.0

        else  ! cat_net must be 0.0 or positive

          nc(jc_h) = cat_net

        endif


! now compute equivalent fractions
      sum_naza = 0.0
      do ja = 1, nanion
        sum_naza = sum_naza + na(ja)*za(ja)
      enddo

      sum_nczc = 0.0
      do jc = 1, ncation
        sum_nczc = sum_nczc + nc(jc)*zc(jc)
      enddo

      !!!!! MARK HERE
      if(sum_naza .eq. 0. .or. sum_nczc .eq. 0.)then
        !write(6,*)'ionic concentrations are zero'
        !write(6,*)'sum_naza = ', sum_naza
        !write(6,*)'sum_nczc = ', sum_nczc
        return
      endif

      do ja = 1, nanion
        xeq_a(ja) = na(ja)*za(ja)/sum_naza
      enddo

      do jc = 1, ncation
        xeq_c(jc) = nc(jc)*zc(jc)/sum_nczc
      enddo

      na_Ma(ja_so4) = na(ja_so4) *MW_a(ja_so4)
      na_Ma(ja_no3) = na(ja_no3) *MW_a(ja_no3)
      na_Ma(ja_cl)  = na(ja_cl)  *MW_a(ja_cl)
      na_Ma(ja_msa) = na(ja_msa) *MW_a(ja_msa)
      na_Ma(ja_hso4)= na(ja_hso4)*MW_a(ja_hso4)

      nc_Mc(jc_ca)  = nc(jc_ca) *MW_c(jc_ca)
      nc_Mc(jc_na)  = nc(jc_na) *MW_c(jc_na)
      nc_Mc(jc_nh4) = nc(jc_nh4)*MW_c(jc_nh4)
      nc_Mc(jc_h)   = nc(jc_h)  *MW_c(jc_h)


! now compute electrolyte moles
      if(xeq_c(jc_na) .gt. 0. .and. xeq_a(ja_so4) .gt. 0.)then
        electrolyte(jna2so4,jp,ibin) = (xeq_c(jc_na) *na_Ma(ja_so4) +   &
                                        xeq_a(ja_so4)*nc_Mc(jc_na))/   &
                                         mw_electrolyte(jna2so4)
      endif

      electrolyte(jnahso4,jp,ibin) = 0.0

      if(xeq_c(jc_na) .gt. 0. .and. xeq_a(ja_msa) .gt. 0.)then
        electrolyte(jnamsa,jp,ibin)  = (xeq_c(jc_na) *na_Ma(ja_msa) +   &
                                        xeq_a(ja_msa)*nc_Mc(jc_na))/   &
                                         mw_electrolyte(jnamsa)
      endif

      if(xeq_c(jc_na) .gt. 0. .and. xeq_a(ja_no3) .gt. 0.)then
        electrolyte(jnano3,jp,ibin)  = (xeq_c(jc_na) *na_Ma(ja_no3) +   &
                                        xeq_a(ja_no3)*nc_Mc(jc_na))/   &
                                         mw_electrolyte(jnano3)
      endif

      if(xeq_c(jc_na) .gt. 0. .and. xeq_a(ja_cl) .gt. 0.)then
        electrolyte(jnacl,jp,ibin)   = (xeq_c(jc_na) *na_Ma(ja_cl) +   &
                                        xeq_a(ja_cl) *nc_Mc(jc_na))/   &
                                         mw_electrolyte(jnacl)
      endif

      if(xeq_c(jc_nh4) .gt. 0. .and. xeq_a(ja_so4) .gt. 0.)then
        electrolyte(jnh4so4,jp,ibin) = (xeq_c(jc_nh4)*na_Ma(ja_so4) +   &
                                        xeq_a(ja_so4)*nc_Mc(jc_nh4))/   &
                                         mw_electrolyte(jnh4so4)
      endif

      electrolyte(jnh4hso4,jp,ibin)= 0.0

      if(xeq_c(jc_nh4) .gt. 0. .and. xeq_a(ja_msa) .gt. 0.)then
        electrolyte(jnh4msa,jp,ibin) = (xeq_c(jc_nh4)*na_Ma(ja_msa) +   &
                                        xeq_a(ja_msa)*nc_Mc(jc_nh4))/   &
                                         mw_electrolyte(jnh4msa)
      endif

      if(xeq_c(jc_nh4) .gt. 0. .and. xeq_a(ja_no3) .gt. 0.)then
        electrolyte(jnh4no3,jp,ibin) = (xeq_c(jc_nh4)*na_Ma(ja_no3) +   &
                                        xeq_a(ja_no3)*nc_Mc(jc_nh4))/   &
                                         mw_electrolyte(jnh4no3)
      endif

      if(xeq_c(jc_nh4) .gt. 0. .and. xeq_a(ja_cl) .gt. 0.)then
        electrolyte(jnh4cl,jp,ibin)  = (xeq_c(jc_nh4)*na_Ma(ja_cl) +   &
                                        xeq_a(ja_cl) *nc_Mc(jc_nh4))/   &
                                         mw_electrolyte(jnh4cl)
      endif

      if(xeq_c(jc_ca) .gt. 0. .and. xeq_a(ja_no3) .gt. 0.0)then
        electrolyte(jcano3, jp,ibin) = (xeq_c(jc_ca) *na_Ma(ja_no3) +   &
                                        xeq_a(ja_no3)*nc_Mc(jc_ca))/   &
                                         mw_electrolyte(jcano3)
      endif

      if(xeq_c(jc_ca) .gt. 0. .and. xeq_a(ja_cl) .gt. 0.)then
        electrolyte(jcacl2,jp,ibin)  = (xeq_c(jc_ca) *na_Ma(ja_cl) +   &
                                        xeq_a(ja_cl) *nc_Mc(jc_ca))/   &
                                         mw_electrolyte(jcacl2)
      endif

      if(xeq_c(jc_ca) .gt. 0. .and. xeq_a(ja_msa) .gt. 0.)then
        electrolyte(jcamsa2,jp,ibin) = (xeq_c(jc_ca) *na_Ma(ja_msa) +   &
                                        xeq_a(ja_msa) *nc_Mc(jc_ca))/   &
                                         mw_electrolyte(jcamsa2)
      endif

      electrolyte(jh2so4, jp,ibin) = 0.0

      if(xeq_c(jc_h) .gt. 0. .and. xeq_a(ja_no3) .gt. 0.)then
      electrolyte(jhno3,jp,ibin)     = (xeq_c(jc_h)  *na_Ma(ja_no3) +   &
                                      xeq_a(ja_no3)*nc_Mc(jc_h))/   &
                                       mw_electrolyte(jhno3)
      endif

      if(xeq_c(jc_h) .gt. 0. .and. xeq_a(ja_cl) .gt. 0.)then
        electrolyte(jhcl,jp,ibin)    = (xeq_c(jc_h) *na_Ma(ja_cl) +   &
                                        xeq_a(ja_cl)*nc_Mc(jc_h))/   &
                                         mw_electrolyte(jhcl)
      endif

      if(xeq_c(jc_h) .gt. 0. .and. xeq_a(ja_msa) .gt. 0.)then
        electrolyte(jmsa,jp,ibin)    = (xeq_c(jc_h) *na_Ma(ja_msa) +   &
                                        xeq_a(ja_msa)*nc_Mc(jc_h))/   &
                                         mw_electrolyte(jmsa)
      endif

!--------------------------------------------------------------------

      elseif(icase.eq.2)then ! XT < 2 : SULFATE RICH DOMAIN

        store(imsa_a) = aer(imsa_a,jp,ibin)
        store(ica_a)  = aer(ica_a, jp,ibin)

        call form_camsa2(store,jp,ibin)

        sum_na_nh4 = aer(ina_a,jp,ibin) + aer(inh4_a,jp,ibin)

        if(sum_na_nh4 .gt. 0.0)then
          f_na  = aer(ina_a,jp,ibin)/sum_na_nh4
          f_nh4 = aer(inh4_a,jp,ibin)/sum_na_nh4
        else
          f_na  = 0.0
          f_nh4 = 0.0
        endif

! first form msa electrolytes
        if(sum_na_nh4 .gt. store(imsa_a))then
          electrolyte(jnamsa,jp,ibin)  = f_na *store(imsa_a)
          electrolyte(jnh4msa,jp,ibin) = f_nh4*store(imsa_a)
          rem_na = aer(ina_a,jp,ibin) - electrolyte(jnamsa,jp,ibin)  ! remaining na
          rem_nh4= aer(inh4_a,jp,ibin)- electrolyte(jnh4msa,jp,ibin) ! remaining nh4
        else
          electrolyte(jnamsa,jp,ibin)  = aer(ina_a,jp,ibin)
          electrolyte(jnh4msa,jp,ibin) = aer(inh4_a,jp,ibin)
          electrolyte(jmsa,jp,ibin)    = store(imsa_a) - sum_na_nh4
          rem_nh4 = 0.0  ! remaining nh4
          rem_na  = 0.0  ! remaining na
        endif


! recompute XT
        if(aer(iso4_a,jp,ibin).gt.0.0)then
          XT = (rem_nh4 + rem_na)/aer(iso4_a,jp,ibin)
        else
          goto 10
        endif

        if(XT .le. 1.0)then	! h2so4 + bisulfate
          xh = (1.0 - XT)
          xb = XT
          electrolyte(jh2so4,jp,ibin)   = xh*aer(iso4_a,jp,ibin)
          electrolyte(jnh4hso4,jp,ibin) = xb*f_nh4*aer(iso4_a,jp,ibin)
          electrolyte(jnahso4,jp,ibin)  = xb*f_na *aer(iso4_a,jp,ibin)
        elseif(XT .le. 1.5)then	! bisulfate + letovicite
          xb = 3.0 - 2.0*XT
          xl = XT - 1.0
          electrolyte(jnh4hso4,jp,ibin) = xb*f_nh4*aer(iso4_a,jp,ibin)
          electrolyte(jnahso4,jp,ibin)  = xb*f_na *aer(iso4_a,jp,ibin)
          electrolyte(jlvcite,jp,ibin)  = xl*f_nh4*aer(iso4_a,jp,ibin)
          electrolyte(jna3hso4,jp,ibin) = xl*f_na *aer(iso4_a,jp,ibin)
        else			! letovicite + sulfate
          xl = 2.0 - XT
          xs = 2.0*XT - 3.0
          electrolyte(jlvcite,jp,ibin)  = xl*f_nh4*aer(iso4_a,jp,ibin)
          electrolyte(jna3hso4,jp,ibin) = xl*f_na *aer(iso4_a,jp,ibin)
          electrolyte(jnh4so4,jp,ibin)  = xs*f_nh4*aer(iso4_a,jp,ibin)
          electrolyte(jna2so4,jp,ibin)  = xs*f_na *aer(iso4_a,jp,ibin)
        endif

        electrolyte(jhno3,jp,ibin) = aer(ino3_a,jp,ibin)
        electrolyte(jhcl,jp,ibin)  = aer(icl_a,jp,ibin)

      endif
!---------------------------------------------------------
!
! calculate % composition  EFFI
10    sum_dum = 0.0
!!      do je = 1, nelectrolyte
!!        sum_dum = sum_dum + electrolyte(je,jp,ibin)
!!      enddo
!!
!!      if(sum_dum .eq. 0.)sum_dum = 1.0
!!      electrolyte_sum(jp,ibin) = sum_dum
!!
!!      do je = 1, nelectrolyte
!!        epercent(je,jp,ibin) = 100.*electrolyte(je,jp,ibin)/sum_dum
!!      enddo
!!

      return
      end subroutine ions_to_electrolytes



























!***********************************************************************
! conforms aerosol generic species to a valid electrolyte composition
!
! author: Rahul A. Zaveri
! update: june 2000
!-----------------------------------------------------------------------
      subroutine conform_electrolytes(jp,ibin,XT)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin, jp
      real(r8) :: XT
! local variables
      integer i, iXT_case, je
      real(r8) :: sum_dum, XNa_prime, XNH4_prime, XT_prime
      real(r8) :: store(naer)

! remove negative concentrations, if any
!      do i=1,naer
!      aer(i,jp,ibin) = max(0.0d0, aer(i,jp,ibin))	! EFFI
!      enddo


!      call calculate_XT(ibin,jp,XT)	! EFFI

      if( (aer(iso4_a,jp,ibin)+aer(imsa_a,jp,ibin)) .gt.0.0)then
        XT   = ( aer(inh4_a,jp,ibin) +   &
                 aer(ina_a,jp,ibin)  +   &
              2.*aer(ica_a,jp,ibin) )/   &
               (aer(iso4_a,jp,ibin)+0.5*aer(imsa_a,jp,ibin))
      else
        XT   = -1.0
      endif


      if(XT .ge. 1.9999 .or. XT.lt.0.)then
       iXT_case = 1	! near neutral (acidity is caused by HCl and/or HNO3)
      else
       iXT_case = 2	! acidic (acidity is caused by excess SO4)
      endif

! initialize
!
! put total aer(*) into store(*)
      store(iso4_a) = aer(iso4_a,jp,ibin)
      store(ino3_a) = aer(ino3_a,jp,ibin)
      store(icl_a)  = aer(icl_a, jp,ibin)
      store(imsa_a) = aer(imsa_a,jp,ibin)
      store(ico3_a) = aer(ico3_a,jp,ibin)
      store(inh4_a) = aer(inh4_a,jp,ibin)
      store(ina_a)  = aer(ina_a, jp,ibin)
      store(ica_a)  = aer(ica_a, jp,ibin)

      do je=1,nelectrolyte
      electrolyte(je,jp,ibin) = 0.0
      enddo

!
!---------------------------------------------------------
!
      if(iXT_case.eq.1)then

! XT >= 2   : sulfate deficient

        call form_caso4(store,jp,ibin)
        call form_camsa2(store,jp,ibin)
        call form_na2so4(store,jp,ibin)
        call form_namsa(store,jp,ibin)
        call form_cano3(store,jp,ibin)
        call form_nano3(store,jp,ibin)
        call form_nacl(store,jp,ibin)
        call form_cacl2(store,jp,ibin)
        call form_caco3(store,jp,ibin)
        call form_nh4so4(store,jp,ibin)
        call form_nh4msa(store,jp,ibin)
        call form_nh4no3(store,jp,ibin)
        call form_nh4cl(store,jp,ibin)
        call form_msa(store,jp,ibin)
        call degas_hno3(store,jp,ibin)
        call degas_hcl(store,jp,ibin)
        call degas_nh3(store,jp,ibin)

      elseif(iXT_case.eq.2)then

! XT < 2   : sulfate enough or sulfate excess

        call form_caso4(store,jp,ibin)
        call form_camsa2(store,jp,ibin)
        call form_namsa(store,jp,ibin)
        call form_nh4msa(store,jp,ibin)
        call form_msa(store,jp,ibin)

        if(store(iso4_a).eq.0.0)goto 10


        XT_prime =(store(ina_a)+store(inh4_a))/   &
                        store(iso4_a)
        XNa_prime=0.5*store(ina_a)/store(iso4_a) + 1.

        if(XT_prime.ge.XNa_prime)then
          call form_na2so4(store,jp,ibin)
          XNH4_prime = 0.0
          if(store(iso4_a).gt.1.e-15)then
            XNH4_prime = store(inh4_a)/store(iso4_a)
          endif

          if(XNH4_prime .ge. 1.5)then
            call form_nh4so4_lvcite(store,jp,ibin)
          else
            call form_lvcite_nh4hso4(store,jp,ibin)
          endif

        elseif(XT_prime.ge.1.)then
          call form_nh4hso4(store,jp,ibin)
          call form_na2so4_nahso4(store,jp,ibin)
        elseif(XT_prime.lt.1.)then
          call form_nahso4(store,jp,ibin)
          call form_nh4hso4(store,jp,ibin)
          call form_h2so4(store,jp,ibin)
        endif

10      call degas_hno3(store,jp,ibin)
        call degas_hcl(store,jp,ibin)
        call degas_nh3(store,jp,ibin)

      endif ! case 1, 2


! re-calculate ions to eliminate round-off errors
      call electrolytes_to_ions(jp, ibin)
!---------------------------------------------------------
!
! calculate % composition  EFFI
!!      sum_dum = 0.0
!!      do je = 1, nelectrolyte
!!        electrolyte(je,jp,ibin) = max(0.d0,electrolyte(je,jp,ibin)) ! remove -ve
!!        sum_dum = sum_dum + electrolyte(je,jp,ibin)
!!      enddo
!!
!!      if(sum_dum .eq. 0.)sum_dum = 1.0
!!      electrolyte_sum(jp,ibin) = sum_dum
!!
!!      do je = 1, nelectrolyte
!!        epercent(je,jp,ibin) = 100.*electrolyte(je,jp,ibin)/sum_dum
!!      enddo
!!
!!
      return
      end subroutine conform_electrolytes











!***********************************************************************
! forms electrolytes from ions
!
! author: Rahul A. Zaveri
! update: june 2000
!-----------------------------------------------------------------------
      subroutine form_electrolytes(jp,ibin,XT)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin, jp
      real(r8) :: XT
! local variables
      integer i, iXT_case, j, je
      real(r8) :: sum_dum, XNa_prime, XNH4_prime, XT_prime
      real(r8) :: store(naer)

! remove negative concentrations, if any
!      do i=1,naer
!        aer(i,jp,ibin) = max(0.0d0, aer(i,jp,ibin))	! EFFI
!      enddo


!      call calculate_XT(ibin,jp,XT)	! EFFI

      if( (aer(iso4_a,jp,ibin)+aer(imsa_a,jp,ibin)) .gt.0.0)then
        XT   = ( aer(inh4_a,jp,ibin) +   &
                 aer(ina_a,jp,ibin)  +   &
              2.*aer(ica_a,jp,ibin) )/   &
               (aer(iso4_a,jp,ibin)+0.5*aer(imsa_a,jp,ibin))
      else
        XT   = -1.0
      endif




      if(XT .ge. 1.9999 .or. XT.lt.0.)then
       iXT_case = 1	! near neutral (acidity is caused by HCl and/or HNO3)
      else
       iXT_case = 2	! acidic (acidity is caused by excess SO4)
      endif

! initialize
!
! put total aer(*) into store(*)
      store(iso4_a) = aer(iso4_a,jp,ibin)
      store(ino3_a) = aer(ino3_a,jp,ibin)
      store(icl_a)  = aer(icl_a, jp,ibin)
      store(imsa_a) = aer(imsa_a,jp,ibin)
      store(ico3_a) = aer(ico3_a,jp,ibin)
      store(inh4_a) = aer(inh4_a,jp,ibin)
      store(ina_a)  = aer(ina_a, jp,ibin)
      store(ica_a)  = aer(ica_a, jp,ibin)

      do j=1,nelectrolyte
        electrolyte(j,jp,ibin) = 0.0
      enddo

!
!---------------------------------------------------------
!
      if(iXT_case.eq.1)then

! XT >= 2   : sulfate deficient
        call form_caso4(store,jp,ibin)
        call form_camsa2(store,jp,ibin)
        call form_na2so4(store,jp,ibin)
        call form_namsa(store,jp,ibin)
        call form_cano3(store,jp,ibin)
        call form_nano3(store,jp,ibin)
        call form_nacl(store,jp,ibin)
        call form_cacl2(store,jp,ibin)
        call form_caco3(store,jp,ibin)
        call form_nh4so4(store,jp,ibin)
        call form_nh4msa(store,jp,ibin)
        call form_nh4no3(store,jp,ibin)
        call form_nh4cl(store,jp,ibin)
        call form_msa(store,jp,ibin)

        if(jp .eq. jsolid)then
          call degas_hno3(store,jp,ibin)
          call degas_hcl(store,jp,ibin)
          call degas_nh3(store,jp,ibin)
        else
          call form_hno3(store,jp,ibin)
          call form_hcl(store,jp,ibin)
          call degas_nh3(store,jp,ibin)
        endif



      elseif(iXT_case.eq.2)then

! XT < 2   : sulfate enough or sulfate excess

        call form_caso4(store,jp,ibin)
        call form_camsa2(store,jp,ibin)
        call form_namsa(store,jp,ibin)
        call form_nh4msa(store,jp,ibin)
        call form_msa(store,jp,ibin)

        if(store(iso4_a).eq.0.0)goto 10


        XT_prime =(store(ina_a)+store(inh4_a))/   &
                        store(iso4_a)
        XNa_prime=0.5*store(ina_a)/store(iso4_a) + 1.

        if(XT_prime.ge.XNa_prime)then
          call form_na2so4(store,jp,ibin)
          XNH4_prime = 0.0
          if(store(iso4_a).gt.1.e-15)then
            XNH4_prime = store(inh4_a)/store(iso4_a)
          endif

          if(XNH4_prime .ge. 1.5)then
            call form_nh4so4_lvcite(store,jp,ibin)
          else
            call form_lvcite_nh4hso4(store,jp,ibin)
          endif

        elseif(XT_prime.ge.1.)then
          call form_nh4hso4(store,jp,ibin)
          call form_na2so4_nahso4(store,jp,ibin)
        elseif(XT_prime.lt.1.)then
          call form_nahso4(store,jp,ibin)
          call form_nh4hso4(store,jp,ibin)
          call form_h2so4(store,jp,ibin)
        endif

10      if(jp .eq. jsolid)then
          call degas_hno3(store,jp,ibin)
          call degas_hcl(store,jp,ibin)
          call degas_nh3(store,jp,ibin)
        else
          call form_hno3(store,jp,ibin)
          call form_hcl(store,jp,ibin)
          call degas_nh3(store,jp,ibin)
        endif

      endif ! case 1, 2


! re-calculate ions to eliminate round-off errors
      call electrolytes_to_ions(jp, ibin)
!---------------------------------------------------------
!
! calculate % composition EFFI
!!      sum_dum = 0.0
!!      do je = 1, nelectrolyte
!!        electrolyte(je,jp,ibin) = max(0.d0,electrolyte(je,jp,ibin)) ! remove -ve  EFFI
!!        sum_dum = sum_dum + electrolyte(je,jp,ibin)
!!      enddo
!!
!!      if(sum_dum .eq. 0.)sum_dum = 1.0
!!      electrolyte_sum(jp,ibin) = sum_dum
!!
!!      do je = 1, nelectrolyte
!!        epercent(je,jp,ibin) = 100.*electrolyte(je,jp,ibin)/sum_dum
!!      enddo


      return
      end subroutine form_electrolytes














!***********************************************************************
! electrolyte formation subroutines
!
! author: Rahul A. Zaveri
! update: june 2000
!-----------------------------------------------------------------------
      subroutine form_caso4(store,jp,ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer jp, ibin
      real(r8) :: store(naer)

      electrolyte(jcaso4,jp,ibin) = min(store(ica_a),store(iso4_a))
      store(ica_a)  = ( (store(ica_a)) -   &
                        (electrolyte(jcaso4,jp,ibin)) )
      store(iso4_a) = ( (store(iso4_a)) -   &
                        (electrolyte(jcaso4,jp,ibin)) )
      store(ica_a)  = max(0.d0, store(ica_a))
      store(iso4_a) = max(0.d0, store(iso4_a))

      return
      end subroutine form_caso4



      subroutine form_camsa2(store,jp,ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer jp, ibin
      real(r8) :: store(naer)

      electrolyte(jcamsa2,jp,ibin) = min(store(ica_a),0.5*store(imsa_a))
      store(ica_a)  = ( (store(ica_a)) -   &
                        (electrolyte(jcamsa2,jp,ibin)) )
      store(imsa_a) = ( (store(imsa_a)) -   &
                        (2.*electrolyte(jcamsa2,jp,ibin)) )
      store(ica_a)  = max(0.d0, store(ica_a))
      store(imsa_a) = max(0.d0, store(imsa_a))

      return
      end subroutine form_camsa2



      subroutine form_cano3(store,jp,ibin)	! Ca(NO3)2
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer jp, ibin
      real(r8) :: store(naer)

      electrolyte(jcano3,jp,ibin) = min(store(ica_a),0.5*store(ino3_a))

      store(ica_a)  = ( (store(ica_a)) -   &
                        (electrolyte(jcano3,jp,ibin)) )
      store(ino3_a) = ( (store(ino3_a)) -   &
                        (2.*electrolyte(jcano3,jp,ibin)) )
      store(ica_a)  = max(0.d0, store(ica_a))
      store(ino3_a) = max(0.d0, store(ino3_a))

      return
      end subroutine form_cano3


      subroutine form_cacl2(store,jp,ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer jp, ibin
      real(r8) :: store(naer)

      electrolyte(jcacl2,jp,ibin) = min(store(ica_a),0.5*store(icl_a))

      store(ica_a)  = ( (store(ica_a)) -   &
                        (electrolyte(jcacl2,jp,ibin)) )
      store(icl_a)  = ( (store(icl_a)) -   &
                        (2.*electrolyte(jcacl2,jp,ibin)) )
      store(ica_a)  = max(0.d0, store(ica_a))
      store(icl_a)  = max(0.d0, store(icl_a))

      return
      end subroutine form_cacl2


      subroutine form_caco3(store,jp,ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer jp, ibin
      real(r8) :: store(naer)

      if(jp.eq.jtotal .or. jp.eq.jsolid)then
      electrolyte(jcaco3,jp,ibin) = store(ica_a)

      aer(ico3_a,jp,ibin)= electrolyte(jcaco3,jp,ibin)	! force co3 = caco3

      store(ica_a) = 0.0
      store(ico3_a)= 0.0
      endif

      return
      end subroutine form_caco3


      subroutine form_na2so4(store,jp,ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer jp, ibin
      real(r8) :: store(naer)

      electrolyte(jna2so4,jp,ibin) = min(.5*store(ina_a),   &
                                            store(iso4_a))
      store(ina_a) =( (store(ina_a)) -   &
                      (2.*electrolyte(jna2so4,jp,ibin)) )
      store(iso4_a)=( (store(iso4_a)) -   &
                      (electrolyte(jna2so4,jp,ibin)) )
      store(ina_a) =max(0.d0, store(ina_a))
      store(iso4_a)=max(0.d0, store(iso4_a))

      return
      end subroutine form_na2so4



      subroutine form_nahso4(store,jp,ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer jp, ibin
      real(r8) :: store(naer)

      electrolyte(jnahso4,jp,ibin) = min(store(ina_a),   &
                                         store(iso4_a))
      store(ina_a)  = ( (store(ina_a)) -   &
                        (electrolyte(jnahso4,jp,ibin)) )
      store(iso4_a) = ( (store(iso4_a)) -   &
                        (electrolyte(jnahso4,jp,ibin)) )
      store(ina_a)  = max(0.d0, store(ina_a))
      store(iso4_a) = max(0.d0, store(iso4_a))

      return
      end subroutine form_nahso4



      subroutine form_namsa(store,jp,ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer jp, ibin
      real(r8) :: store(naer)

      electrolyte(jnamsa,jp,ibin) = min(store(ina_a),   &
                                        store(imsa_a))
      store(ina_a)  = ( (store(ina_a)) -   &
                        (electrolyte(jnamsa,jp,ibin)) )
      store(imsa_a) = ( (store(imsa_a)) -   &
                        (electrolyte(jnamsa,jp,ibin)) )
      store(ina_a)  = max(0.d0, store(ina_a))
      store(imsa_a) = max(0.d0, store(imsa_a))

      return
      end subroutine form_namsa



      subroutine form_nano3(store,jp,ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer jp, ibin
      real(r8) :: store(naer)

      electrolyte(jnano3,jp,ibin)=min(store(ina_a),store(ino3_a))
      store(ina_a)  = ( (store(ina_a)) -   &
                        (electrolyte(jnano3,jp,ibin)) )
      store(ino3_a) = ( (store(ino3_a)) -   &
                        (electrolyte(jnano3,jp,ibin)) )
      store(ina_a)  = max(0.d0, store(ina_a))
      store(ino3_a) = max(0.d0, store(ino3_a))

      return
      end subroutine form_nano3



      subroutine form_nacl(store,jp,ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer jp, ibin
      real(r8) :: store(naer)

      electrolyte(jnacl,jp,ibin) = store(ina_a)

      store(ina_a) = 0.0
      store(icl_a) = ( (store(icl_a)) -   &
                       (electrolyte(jnacl,jp,ibin)) )

      if(store(icl_a) .lt. 0.)then 				! cl deficit in aerosol. take some from gas
        aer(icl_a,jp,ibin)= aer(icl_a,jp,ibin)- store(icl_a)	! update aer(icl_a)

        if(jp .ne. jtotal)then
          aer(icl_a,jtotal,ibin)= aer(icl_a,jliquid,ibin)+ &	! update for jtotal
                                  aer(icl_a,jsolid,ibin)
        endif

        gas(ihcl_g) = gas(ihcl_g) + store(icl_a)			! update gas(ihcl_g)

        if(gas(ihcl_g) .lt. 0.0)then
          total_species(ihcl_g) = total_species(ihcl_g) - gas(ihcl_g)	! update total_species
          tot_cl_in = tot_cl_in - gas(ihcl_g)				! update tot_cl_in
        endif

        gas(ihcl_g) = max(0.d0, gas(ihcl_g))				! restrict gas(ihcl_g) to >= 0.
        store(icl_a) = 0.        					! force store(icl_a) to 0.

      endif

      store(icl_a) = max(0.d0, store(icl_a))

      return
      end subroutine form_nacl



      subroutine form_nh4so4(store,jp,ibin)	! (nh4)2so4
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer jp, ibin
      real(r8) :: store(naer)

      electrolyte(jnh4so4,jp,ibin)= min(.5*store(inh4_a),   &
                                           store(iso4_a))
      store(inh4_a)= ( (store(inh4_a)) -   &
                       (2.*electrolyte(jnh4so4,jp,ibin)) )
      store(iso4_a)= ( (store(iso4_a)) -   &
                       (electrolyte(jnh4so4,jp,ibin)) )
      store(inh4_a) = max(0.d0, store(inh4_a))
      store(iso4_a) = max(0.d0, store(iso4_a))

      return
      end subroutine form_nh4so4



      subroutine form_nh4hso4(store,jp,ibin)	! nh4hso4
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer jp, ibin
      real(r8) :: store(naer)

      electrolyte(jnh4hso4,jp,ibin) = min(store(inh4_a),   &
                                          store(iso4_a))
      store(inh4_a)= ( (store(inh4_a)) -   &
                       (electrolyte(jnh4hso4,jp,ibin)) )
      store(iso4_a)= ( (store(iso4_a)) -   &
                       (electrolyte(jnh4hso4,jp,ibin)) )
      store(inh4_a) = max(0.d0, store(inh4_a))
      store(iso4_a) = max(0.d0, store(iso4_a))

      return
      end subroutine form_nh4hso4



      subroutine form_nh4msa(store,jp,ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer jp, ibin
      real(r8) :: store(naer)

      electrolyte(jnh4msa,jp,ibin) = min(store(inh4_a),   &
                                         store(imsa_a))
      store(inh4_a) = ( (store(inh4_a)) -   &
                        (electrolyte(jnh4msa,jp,ibin)) )
      store(imsa_a) = ( (store(imsa_a)) -   &
                        (electrolyte(jnh4msa,jp,ibin)) )
      store(inh4_a) = max(0.d0, store(inh4_a))
      store(imsa_a) = max(0.d0, store(imsa_a))

      return
      end subroutine form_nh4msa



      subroutine form_nh4cl(store,jp,ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer jp, ibin
      real(r8) :: store(naer)

      electrolyte(jnh4cl,jp,ibin) = min(store(inh4_a),   &
                                        store(icl_a))
      store(inh4_a) = ( (store(inh4_a)) -   &
                        (electrolyte(jnh4cl,jp,ibin)) )
      store(icl_a)  = ( (store(icl_a)) -   &
                        (electrolyte(jnh4cl,jp,ibin)) )
      store(inh4_a) = max(0.d0, store(inh4_a))
      store(icl_a)  = max(0.d0, store(icl_a))

      return
      end subroutine form_nh4cl



      subroutine form_nh4no3(store,jp,ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer jp, ibin
      real(r8) :: store(naer)

      electrolyte(jnh4no3,jp,ibin) = min(store(inh4_a),   &
                                         store(ino3_a))
      store(inh4_a) = ( (store(inh4_a)) -   &
                        (electrolyte(jnh4no3,jp,ibin)) )
      store(ino3_a) = ( (store(ino3_a)) -   &
                        (electrolyte(jnh4no3,jp,ibin)) )
      store(inh4_a) = max(0.d0, store(inh4_a))
      store(ino3_a) = max(0.d0, store(ino3_a))

      return
      end subroutine form_nh4no3



      subroutine form_nh4so4_lvcite(store,jp,ibin) ! (nh4)2so4 + (nh4)3h(so4)2
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer jp, ibin
      real(r8) :: store(naer)

      electrolyte(jnh4so4,jp,ibin)= ( (2.*store(inh4_a)) -   &
                                      (3.*store(iso4_a)) )
      electrolyte(jlvcite,jp,ibin)= ( (2.*store(iso4_a)) -   &
                                      (store(inh4_a)) )
      electrolyte(jnh4so4,jp,ibin)= max(0.d0,   &
                                    electrolyte(jnh4so4,jp,ibin))
      electrolyte(jlvcite,jp,ibin)= max(0.d0,   &
                                    electrolyte(jlvcite,jp,ibin))
      store(inh4_a) = 0.
      store(iso4_a) = 0.

      return
      end subroutine form_nh4so4_lvcite



      subroutine form_lvcite_nh4hso4(store,jp,ibin) ! (nh4)3h(so4)2 + nh4hso4
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer jp, ibin
      real(r8) :: store(naer)

      electrolyte(jlvcite,jp,ibin) = ( (store(inh4_a)) -   &
                                       (store(iso4_a)) )
      electrolyte(jnh4hso4,jp,ibin)= ( (3.*store(iso4_a)) -   &
                                       (2.*store(inh4_a)) )
      electrolyte(jlvcite,jp,ibin) = max(0.d0,   &
                                      electrolyte(jlvcite,jp,ibin))
      electrolyte(jnh4hso4,jp,ibin)= max(0.d0,   &
                                      electrolyte(jnh4hso4,jp,ibin))
      store(inh4_a) = 0.
      store(iso4_a) = 0.

      return
      end subroutine form_lvcite_nh4hso4



      subroutine form_na2so4_nahso4(store,jp,ibin) ! na2so4 + nahso4
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer jp, ibin
      real(r8) :: store(naer)

      electrolyte(jna2so4,jp,ibin)= ( (store(ina_a)) -   &
                                      (store(iso4_a)) )
      electrolyte(jnahso4,jp,ibin)= ( (2.*store(iso4_a))-   &
                                      (store(ina_a)) )
      electrolyte(jna2so4,jp,ibin)= max(0.d0,   &
                                    electrolyte(jna2so4,jp,ibin))
      electrolyte(jnahso4,jp,ibin)= max(0.d0,   &
                                    electrolyte(jnahso4,jp,ibin))
      store(ina_a)  = 0.
      store(iso4_a) = 0.

!	write(6,*)'na2so4 + nahso4'

      return
      end subroutine form_na2so4_nahso4




      subroutine form_h2so4(store,jp,ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer jp, ibin
      real(r8) :: store(naer)

      electrolyte(jh2so4,jp,ibin) = max(0.0d0, store(iso4_a))
      store(iso4_a) = 0.0

      return
      end subroutine form_h2so4




      subroutine form_msa(store,jp,ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer jp, ibin
      real(r8) :: store(naer)

      electrolyte(jmsa,jp,ibin) = max(0.0d0, store(imsa_a))
      store(imsa_a) = 0.0

      return
      end subroutine form_msa



      subroutine form_hno3(store,jp,ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer jp, ibin
      real(r8) :: store(naer)

      electrolyte(jhno3,jp,ibin) = max(0.0d0, store(ino3_a))
      store(ino3_a) = 0.0

      return
      end subroutine form_hno3




      subroutine form_hcl(store,jp,ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer jp, ibin
      real(r8) :: store(naer)

      electrolyte(jhcl,jp,ibin) = max(0.0d0, store(icl_a))
      store(icl_a) = 0.0

      return
      end subroutine form_hcl




      subroutine degas_hno3(store,jp,ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer jp, ibin
      real(r8) :: store(naer)

      store(ino3_a) = max(0.0d0, store(ino3_a))
      gas(ihno3_g) = gas(ihno3_g) + store(ino3_a)
      aer(ino3_a,jp,ibin) = ( (aer(ino3_a,jp,ibin)) -   &
                              (store(ino3_a)) )
      aer(ino3_a,jp,ibin) = max(0.0d0,aer(ino3_a,jp,ibin))

! also do it for jtotal
      if(jp .ne. jtotal)then
        aer(ino3_a,jtotal,ibin) = aer(ino3_a,jsolid, ibin) +   &
                                  aer(ino3_a,jliquid,ibin)
      endif

      electrolyte(jhno3,jp,ibin) = 0.0
      store(ino3_a) = 0.0

      return
      end subroutine degas_hno3



      subroutine degas_hcl(store,jp,ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer jp, ibin
      real(r8) :: store(naer)

      store(icl_a) = max(0.0d0, store(icl_a))
      gas(ihcl_g) = gas(ihcl_g) + store(icl_a)
      aer(icl_a,jp,ibin) = ( (aer(icl_a,jp,ibin)) -   &
                             (store(icl_a)) )
      aer(icl_a,jp,ibin) = max(0.0d0,aer(icl_a,jp,ibin))

! also do it for jtotal
      if(jp .ne. jtotal)then
        aer(icl_a,jtotal,ibin) = aer(icl_a,jsolid, ibin) +   &
                                 aer(icl_a,jliquid,ibin)
      endif

      electrolyte(jhcl,jp,ibin) = 0.0
      store(icl_a) = 0.0

      return
      end subroutine degas_hcl



      subroutine degas_nh3(store,jp,ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer jp, ibin
      real(r8) :: store(naer)

      store(inh4_a) = max(0.0d0, store(inh4_a))
      gas(inh3_g) = gas(inh3_g) + store(inh4_a)
      aer(inh4_a,jp,ibin) = ( (aer(inh4_a,jp,ibin)) -   &
                              (store(inh4_a)) )
      aer(inh4_a,jp,ibin) = max(0.0d0,aer(inh4_a,jp,ibin))

! also do it for jtotal
      if(jp .ne. jtotal)then
        aer(inh4_a,jtotal,ibin)= aer(inh4_a,jsolid, ibin) +   &
                                 aer(inh4_a,jliquid,ibin)
      endif

      store(inh4_a) = 0.0

      return
      end subroutine degas_nh3









      subroutine degas_acids(jp,ibin,XT)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer jp, ibin
      real(r8) :: XT
! local variables
      real(r8) :: ehno3, ehcl



      if(jp .ne. jliquid)then
        write(6,*)'Error in degas_acids'
        write(6,*)'wrong jp'
      endif

      ehno3 = electrolyte(jhno3,jp,ibin)
      ehcl  = electrolyte(jhcl,jp,ibin)

! add to gas
      gas(ihno3_g) = gas(ihno3_g) + ehno3
      gas(ihcl_g)  = gas(ihcl_g)  + ehcl

! remove from aer
      aer(ino3_a,jp,ibin) = aer(ino3_a,jp,ibin) - ehno3
      aer(icl_a, jp,ibin) = aer(icl_a, jp,ibin) - ehcl

! update jtotal
      aer(ino3_a,jtotal,ibin) = aer(ino3_a,jliquid,ibin) +   &
                                aer(ino3_a,jsolid, ibin)

      aer(icl_a,jtotal,ibin)  = aer(icl_a,jliquid,ibin) +   &
                                aer(icl_a,jsolid, ibin)

      electrolyte(jhno3,jp,ibin) = 0.0
      electrolyte(jhcl,jp,ibin)  = 0.0

      return
      end subroutine degas_acids














!***********************************************************************
! subroutines to evaporate solid volatile species
!
! author: Rahul A. Zaveri
! update: sep 2004
!-----------------------------------------------------------------------
!
! nh4no3 (solid)
      subroutine degas_solid_nh4no3(ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      integer jp
      real(r8) :: a, b, c, xgas, XT
      real(r8) :: quadratic					! mosaic func


      jp = jsolid

      a = 1.0
      b = gas(inh3_g) + gas(ihno3_g)
      c = gas(inh3_g)*gas(ihno3_g) - Keq_sg(1)
      xgas = quadratic(a,b,c)

      if(xgas .ge. electrolyte(jnh4no3,jp,ibin))then ! degas all nh4no3

          gas(inh3_g) = gas(inh3_g)  + electrolyte(jnh4no3,jp,ibin)
          gas(ihno3_g)= gas(ihno3_g) + electrolyte(jnh4no3,jp,ibin)
          aer(inh4_a,jp,ibin) = aer(inh4_a,jp,ibin) -   &
                                electrolyte(jnh4no3,jp,ibin)
          aer(ino3_a,jp,ibin) = aer(ino3_a,jp,ibin) -   &
                                electrolyte(jnh4no3,jp,ibin)

      else	! degas only xgas amount of nh4no3

          gas(inh3_g) = gas(inh3_g)  + xgas
          gas(ihno3_g)= gas(ihno3_g) + xgas
          aer(inh4_a,jp,ibin) = aer(inh4_a,jp,ibin) - xgas
          aer(ino3_a,jp,ibin) = aer(ino3_a,jp,ibin) - xgas
      endif


! update jtotal
      aer(inh4_a,jtotal,ibin)  = aer(inh4_a,jsolid,ibin) +   &
                                 aer(inh4_a,jliquid,ibin)
      aer(ino3_a,jtotal,ibin)  = aer(ino3_a,jsolid,ibin) +   &
                                 aer(ino3_a,jliquid,ibin)

      return
      end subroutine degas_solid_nh4no3









! nh4cl (solid)
      subroutine degas_solid_nh4cl(ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      integer jp
      real(r8) :: a, b, c, xgas, XT
      real(r8) :: quadratic					! mosaic func


      jp = jsolid

      a = 1.0
      b = gas(inh3_g) + gas(ihcl_g)
      c = gas(inh3_g)*gas(ihcl_g) - Keq_sg(2)
      xgas = quadratic(a,b,c)

      if(xgas .ge. electrolyte(jnh4cl,jp,ibin))then ! degas all nh4cl

          gas(inh3_g) = gas(inh3_g) + electrolyte(jnh4cl,jp,ibin)
          gas(ihcl_g) = gas(ihcl_g) + electrolyte(jnh4cl,jp,ibin)
          aer(inh4_a,jp,ibin) = aer(inh4_a,jp,ibin) -   &
                                electrolyte(jnh4cl,jp,ibin)
          aer(icl_a,jp,ibin)  = aer(icl_a,jp,ibin) -   &
                                electrolyte(jnh4cl,jp,ibin)

      else	! degas only xgas amount of nh4cl

          gas(inh3_g) = gas(inh3_g) + xgas
          gas(ihcl_g) = gas(ihcl_g) + xgas
          aer(inh4_a,jp,ibin) = aer(inh4_a,jp,ibin) - xgas
          aer(icl_a,jp,ibin)  = aer(icl_a,jp,ibin)  - xgas

      endif


! update jtotal
      aer(inh4_a,jtotal,ibin)  = aer(inh4_a,jsolid,ibin) +   &
                                 aer(inh4_a,jliquid,ibin)
      aer(icl_a,jtotal,ibin)   = aer(icl_a,jsolid,ibin)  +   &
                                 aer(icl_a,jliquid,ibin)

      return
      end subroutine degas_solid_nh4cl











!***********************************************************************
! subroutines to absorb and degas small amounts of volatile species
!
! author: Rahul A. Zaveri
! update: jun 2002
!-----------------------------------------------------------------------
!
! nh4no3 (liquid)
      subroutine absorb_tiny_nh4no3(ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      integer je
      real(r8) :: small_aer, small_gas, small_amt


!! EFFI
      electrolyte_sum(jtotal,ibin) = 0.0
      do je = 1, nelectrolyte
        electrolyte_sum(jtotal,ibin) = electrolyte_sum(jtotal,ibin) + &
                                       electrolyte(je,jtotal,ibin)
      enddo
!! EFFI


      small_gas = 0.01 * min(delta_nh3_max(ibin),delta_hno3_max(ibin))
      small_aer = 0.01 * electrolyte_sum(jtotal,ibin)
      if(small_aer .eq. 0.0)small_aer = small_gas

      small_amt = min(small_gas, small_aer)

      aer(inh4_a,jliquid,ibin) = aer(inh4_a,jliquid,ibin) + small_amt
      aer(ino3_a,jliquid,ibin) = aer(ino3_a,jliquid,ibin) + small_amt

! update jtotal
      aer(inh4_a,jtotal,ibin)  = aer(inh4_a,jsolid,ibin) +   &
                                 aer(inh4_a,jliquid,ibin)
      aer(ino3_a,jtotal,ibin)  = aer(ino3_a,jsolid,ibin) +   &
                                 aer(ino3_a,jliquid,ibin)

! update gas
      gas(inh3_g)  = ((gas(inh3_g)) - (small_amt))
      gas(ihno3_g) = ((gas(ihno3_g)) - (small_amt))

      return
      end subroutine absorb_tiny_nh4no3






!--------------------------------------------------------------------
! nh4cl (liquid)
      subroutine absorb_tiny_nh4cl(ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      integer je
      real(r8) :: small_aer, small_gas, small_amt


!! EFFI
      electrolyte_sum(jtotal,ibin) = 0.0
      do je = 1, nelectrolyte
        electrolyte_sum(jtotal,ibin) = electrolyte_sum(jtotal,ibin) + &
                                       electrolyte(je,jtotal,ibin)
      enddo
!! EFFI



      small_gas = 0.01 * min(delta_nh3_max(ibin), delta_hcl_max(ibin))
      small_aer = 0.01 * electrolyte_sum(jtotal,ibin)
      if(small_aer .eq. 0.0)small_aer = small_gas

      small_amt = min(small_gas, small_aer)

      aer(inh4_a,jliquid,ibin) = aer(inh4_a,jliquid,ibin) + small_amt
      aer(icl_a,jliquid,ibin)  = aer(icl_a,jliquid,ibin)  + small_amt

! update jtotal
      aer(inh4_a,jtotal,ibin)  = aer(inh4_a,jsolid,ibin) +   &
                                 aer(inh4_a,jliquid,ibin)
      aer(icl_a,jtotal,ibin)   = aer(icl_a,jsolid,ibin)  +   &
                                 aer(icl_a,jliquid,ibin)

! update gas
      gas(inh3_g) = ((gas(inh3_g)) - (small_amt))
      gas(ihcl_g) = ((gas(ihcl_g)) - (small_amt))

      return
      end subroutine absorb_tiny_nh4cl













!--------------------------------------------------------------
! nh4no3 (liquid)
      subroutine degas_tiny_nh4no3(ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      real(r8) :: small_amt

      small_amt = 0.01 * electrolyte(jnh4no3,jliquid,ibin)

      aer(inh4_a,jliquid,ibin) = ((aer(inh4_a,jliquid,ibin))-   &
                                  (small_amt))
      aer(ino3_a,jliquid,ibin) = ((aer(ino3_a,jliquid,ibin))-   &
                                  (small_amt))

! update jtotal
      aer(inh4_a,jtotal,ibin)  = aer(inh4_a,jsolid,ibin) +   &
                                 aer(inh4_a,jliquid,ibin)
      aer(ino3_a,jtotal,ibin)  = aer(ino3_a,jsolid,ibin) +   &
                                 aer(ino3_a,jliquid,ibin)

! update gas
      gas(inh3_g)  = gas(inh3_g)  + small_amt
      gas(ihno3_g) = gas(ihno3_g) + small_amt

      return
      end subroutine degas_tiny_nh4no3




!--------------------------------------------------------------------
! nh4cl (liquid)
      subroutine degas_tiny_nh4cl(ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      real(r8) :: small_amt


      small_amt = 0.01 * electrolyte(jnh4cl,jliquid,ibin)

      aer(inh4_a,jliquid,ibin) = ((aer(inh4_a,jliquid,ibin))-   &
                                  (small_amt))
      aer(icl_a,jliquid,ibin)  = ((aer(icl_a,jliquid,ibin))-   &
                                  (small_amt))

! update jtotal
      aer(inh4_a,jtotal,ibin)  = aer(inh4_a,jsolid,ibin) +   &
                                 aer(inh4_a,jliquid,ibin)
      aer(icl_a,jtotal,ibin)   = aer(icl_a,jsolid,ibin)  +   &
                                 aer(icl_a,jliquid,ibin)

! update gas
      gas(inh3_g) = gas(inh3_g) + small_amt
      gas(ihcl_g) = gas(ihcl_g) + small_amt

      return
      end subroutine degas_tiny_nh4cl







!--------------------------------------------------------------
! hcl (liquid)
      subroutine absorb_tiny_hcl(ibin)	! and degas tiny hno3
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      real(r8) :: small_aer, small_amt, small_gas

      small_gas = 0.01 * delta_hcl_max(ibin)
      small_aer = 0.01 * aer(ino3_a,jliquid,ibin)

      small_amt = min(small_gas, small_aer)

! absorb tiny hcl
      aer(icl_a,jliquid,ibin)= aer(icl_a,jliquid,ibin) + small_amt
      aer(icl_a,jtotal,ibin) = aer(icl_a,jsolid,ibin) +   &
                               aer(icl_a,jliquid,ibin)
      gas(ihcl_g) = ((gas(ihcl_g))-(small_amt))

! degas tiny hno3
      aer(ino3_a,jliquid,ibin) = ((aer(ino3_a,jliquid,ibin))-   &
                                  (small_amt))
      aer(ino3_a,jtotal,ibin)  = aer(ino3_a,jsolid,ibin) +   &
                                 aer(ino3_a,jliquid,ibin)

! update gas
      gas(ihno3_g) = gas(ihno3_g) + small_amt

      return
      end subroutine absorb_tiny_hcl



!--------------------------------------------------------------------
! hno3 (liquid)
      subroutine absorb_tiny_hno3(ibin)	! and degas tiny hcl
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      real(r8) :: small_aer, small_amt, small_gas

      small_gas = 0.01 * delta_hno3_max(ibin)
      small_aer = 0.01 * aer(icl_a,jliquid,ibin)

      small_amt = min(small_gas, small_aer)

! absorb tiny hno3
      aer(ino3_a,jliquid,ibin) = aer(ino3_a,jliquid,ibin) + small_amt
      aer(ino3_a,jtotal,ibin)  = aer(ino3_a,jsolid,ibin) +   &
                                 aer(ino3_a,jliquid,ibin)
      gas(ihno3_g) = ((gas(ihno3_g))-(small_amt))

! degas tiny hcl
      aer(icl_a,jliquid,ibin)  = ((aer(icl_a,jliquid,ibin))-   &
                                  (small_amt))
      aer(icl_a,jtotal,ibin)   = aer(icl_a,jsolid,ibin) +   &
                                 aer(icl_a,jliquid,ibin)

! update gas
      gas(ihcl_g) = gas(ihcl_g) + small_amt

      return
      end subroutine absorb_tiny_hno3









!***********************************************************************
! subroutines to equilibrate volatile acids
!
! author: Rahul A. Zaveri
! update: may 2002
!-----------------------------------------------------------------------
      subroutine equilibrate_acids(ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables



      if(gas(ihcl_g)*gas(ihno3_g) .gt. 0.)then
        call equilibrate_hcl_and_hno3(ibin)
      elseif(gas(ihcl_g) .gt. 0.)then
        call equilibrate_hcl(ibin)
      elseif(gas(ihno3_g) .gt. 0.)then
        call equilibrate_hno3(ibin)
      endif


      return
      end subroutine equilibrate_acids








! only hcl
      subroutine equilibrate_hcl(ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      real(r8) :: a, aerH, aerHSO4, aerSO4, b, c, dum, Kdash_hcl, mH, Tcl,   &
        W, XT, Z
      real(r8) :: quadratic					! mosaic func

      aerSO4 = ma(ja_so4,ibin)*water_a(ibin)*1.e+9
      aerHSO4= ma(ja_hso4,ibin)*water_a(ibin)*1.e+9

      Tcl = aer(icl_a,jliquid,ibin) + gas(ihcl_g)	! nmol/m^3(air)
      Kdash_hcl = Keq_gl(4)*1.e+18/gam(jhcl,ibin)**2	! (nmol^2/kg^2)/(nmol/m^3(air))
      Z = (   aer(ina_a, jliquid,ibin) + 		   &  ! nmol/m^3(air)
              aer(inh4_a,jliquid,ibin) +   &
           2.*aer(ica_a, jliquid,ibin) ) -   &
          (2.*aerSO4  +   &
              aerHSO4 +   &
              aer(ino3_a,jliquid,ibin) )


      W     = water_a(ibin)				! kg/m^3(air)

      Kdash_hcl = Keq_gl(4)*1.e+18/gam(jhcl,ibin)**2	! (nmol^2/kg^2)/(nmol/m^3(air))
      a = 1.0
      b = ((Kdash_hcl*W) + (Z/W))*1.e-9
      c = Kdash_hcl*(Z - Tcl)*1.e-18


      dum = ((b*b)-(4.*a*c))
      if (dum .lt. 0.) return		! no real root


      if(c .lt. 0.)then
        mH = quadratic(a,b,c)	! mol/kg(water)
        aerH = mH*W*1.e+9
        aer(icl_a,jliquid,ibin) = ((aerH) + (Z))
      else
        mH = sqrt(Keq_ll(3))
      endif

      call form_electrolytes(jliquid,ibin,XT)

! update gas phase concentration
      gas(ihcl_g) = ( (Tcl)  - (aer(icl_a,jliquid,ibin))  )


! update the following molalities
      ma(ja_so4,ibin)  = 1.e-9*aerSO4/water_a(ibin)
      ma(ja_hso4,ibin) = 1.e-9*aerHSO4/water_a(ibin)
      ma(ja_no3,ibin)  = 1.e-9*aer(ino3_a,jliquid,ibin)/water_a(ibin)
      ma(ja_cl,ibin)   = 1.e-9*aer(icl_a, jliquid,ibin)/water_a(ibin)

      mc(jc_h,ibin)    = mH
      mc(jc_ca,ibin)   = 1.e-9*aer(ica_a, jliquid,ibin)/water_a(ibin)
      mc(jc_nh4,ibin)  = 1.e-9*aer(inh4_a,jliquid,ibin)/water_a(ibin)
      mc(jc_na,ibin)   = 1.e-9*aer(ina_a, jliquid,ibin)/water_a(ibin)


! update the following activities
      activity(jhcl,ibin)    = mc(jc_h,ibin)  *ma(ja_cl,ibin)  *   &
                               gam(jhcl,ibin)**2

      activity(jhno3,ibin)   = mc(jc_h,ibin)  *ma(ja_no3,ibin) *   &
                               gam(jhno3,ibin)**2

      activity(jnh4cl,ibin)  = mc(jc_nh4,ibin)*ma(ja_cl,ibin) *   &
                               gam(jnh4cl,ibin)**2


! also update xyz(jtotal)
      aer(icl_a,jtotal,ibin) = aer(icl_a,jliquid,ibin) +   &
                               aer(icl_a,jsolid,ibin)

      electrolyte(jhcl,jtotal,ibin) = electrolyte(jhcl,jliquid,ibin)

      return
      end subroutine equilibrate_hcl




! only hno3
      subroutine equilibrate_hno3(ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      real(r8) :: a, aerH, aerHSO4, aerSO4, b, c, dum, Kdash_hno3, mH,   &
        Tno3, W, XT, Z
      real(r8) :: quadratic					! mosaic func

      aerSO4 = ma(ja_so4,ibin)*water_a(ibin)*1.e+9
      aerHSO4= ma(ja_hso4,ibin)*water_a(ibin)*1.e+9

      Tno3 = aer(ino3_a,jliquid,ibin) + gas(ihno3_g)	! nmol/m^3(air)
      Kdash_hno3 = Keq_gl(3)*1.e+18/gam(jhno3,ibin)**2	! (nmol^2/kg^2)/(nmol/m^3(air))
      Z = (   aer(ina_a, jliquid,ibin) + 		   &  ! nmol/m^3(air)
              aer(inh4_a,jliquid,ibin) +   &
           2.*aer(ica_a, jliquid,ibin) ) -   &
          (2.*aerSO4  +   &
              aerHSO4 +   &
              aer(icl_a,jliquid,ibin) )


      W     = water_a(ibin)				! kg/m^3(air)

      Kdash_hno3 = Keq_gl(3)*1.e+18/gam(jhno3,ibin)**2	! (nmol^2/kg^2)/(nmol/m^3(air))
      a = 1.0
      b = ((Kdash_hno3*W) + (Z/W))*1.e-9
      c = Kdash_hno3*(Z - Tno3)*1.e-18

      dum = ((b*b)-(4.*a*c))
      if (dum .lt. 0.) return		! no real root



      if(c .lt. 0.)then
        mH = quadratic(a,b,c)	! mol/kg(water)
        aerH = mH*W*1.e+9
        aer(ino3_a,jliquid,ibin) = ((aerH) + (Z))
      else
        mH = sqrt(Keq_ll(3))
      endif

      call form_electrolytes(jliquid,ibin,XT)

! update gas phase concentration
      gas(ihno3_g)= ( (Tno3) - (aer(ino3_a,jliquid,ibin)) )


! update the following molalities
      ma(ja_so4,ibin)  = 1.e-9*aerSO4/water_a(ibin)
      ma(ja_hso4,ibin) = 1.e-9*aerHSO4/water_a(ibin)
      ma(ja_no3,ibin)  = 1.e-9*aer(ino3_a,jliquid,ibin)/water_a(ibin)
      ma(ja_cl,ibin)   = 1.e-9*aer(icl_a, jliquid,ibin)/water_a(ibin)

      mc(jc_h,ibin)    = mH
      mc(jc_ca,ibin)   = 1.e-9*aer(ica_a, jliquid,ibin)/water_a(ibin)
      mc(jc_nh4,ibin)  = 1.e-9*aer(inh4_a,jliquid,ibin)/water_a(ibin)
      mc(jc_na,ibin)   = 1.e-9*aer(ina_a, jliquid,ibin)/water_a(ibin)


! update the following activities
      activity(jhcl,ibin)    = mc(jc_h,ibin)  *ma(ja_cl,ibin)  *   &
                               gam(jhcl,ibin)**2

      activity(jhno3,ibin)   = mc(jc_h,ibin)  *ma(ja_no3,ibin) *   &
                               gam(jhno3,ibin)**2

      activity(jnh4no3,ibin) = mc(jc_nh4,ibin)*ma(ja_no3,ibin) *   &
                               gam(jnh4no3,ibin)**2


! also update xyz(jtotal)
      aer(ino3_a,jtotal,ibin) = aer(ino3_a,jliquid,ibin) +   &
                                aer(ino3_a,jsolid,ibin)

      electrolyte(jhno3,jtotal,ibin) = electrolyte(jhno3,jliquid,ibin)

      return
      end subroutine equilibrate_hno3










! both hcl and hno3
      subroutine equilibrate_hcl_and_hno3(ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      real(r8) :: aerH, aerHSO4, aerSO4, Kdash_hcl, Kdash_hno3,   &
        mH, p, q, r, Tcl, Tno3, W, XT, Z
      real(r8) :: cubic					! mosaic func


      aerSO4 = ma(ja_so4,ibin)*water_a(ibin)*1.e+9
      aerHSO4= ma(ja_hso4,ibin)*water_a(ibin)*1.e+9

      Tcl  = aer(icl_a,jliquid,ibin)  + gas(ihcl_g)	! nmol/m^3(air)
      Tno3 = aer(ino3_a,jliquid,ibin) + gas(ihno3_g)	! nmol/m^3(air)

      Kdash_hcl  = Keq_gl(4)*1.e+18/gam(jhcl,ibin)**2	! (nmol^2/kg^2)/(nmol/m^3(air))
      Kdash_hno3 = Keq_gl(3)*1.e+18/gam(jhno3,ibin)**2	! (nmol^2/kg^2)/(nmol/m^3(air))

      Z = (   aer(ina_a, jliquid,ibin) + 		   &  ! nmol/m^3(air)
              aer(inh4_a,jliquid,ibin) +   &
           2.*aer(ica_a, jliquid,ibin) ) -   &
          (2.*aerSO4 + aerHSO4 )


      W = water_a(ibin)

      Kdash_hcl  = Keq_gl(4)*1.e+18/gam(jhcl,ibin)**2	! (nmol^2/kg^2)/(nmol/m^3(air))
      Kdash_hno3 = Keq_gl(3)*1.e+18/gam(jhno3,ibin)**2	! (nmol^2/kg^2)/(nmol/m^3(air))

      p = (Z/W + W*(Kdash_hcl + Kdash_hno3))*1.e-9

      q = 1.e-18*Kdash_hcl*Kdash_hno3*W**2  +   &
          1.e-18*Z*(Kdash_hcl + Kdash_hno3) -   &
          1.e-18*Kdash_hcl*Tcl -   &
          1.e-18*Kdash_hno3*Tno3

      r = 1.e-18*Kdash_hcl*Kdash_hno3*W*(Z - Tcl - Tno3)*1.e-9

      mH = cubic(p,q,r)

      if(mH .gt. 0.0)then
        aerH = mH*W*1.e+9
        aer(ino3_a,jliquid,ibin) = Kdash_hno3*W*W*Tno3/   &
                                  (aerH + Kdash_hno3*W*W)
        aer(icl_a, jliquid,ibin) = Kdash_hcl*W*W*Tcl/   &
                                  (aerH + Kdash_hcl*W*W)
      else
        mH = sqrt(Keq_ll(3))
      endif

      call form_electrolytes(jliquid,ibin,XT)

! update gas phase concentration
      gas(ihno3_g)= ( (Tno3) - (aer(ino3_a,jliquid,ibin)) )
      gas(ihcl_g) = ( (Tcl)  - (aer(icl_a,jliquid,ibin))  )


! update the following molalities
      ma(ja_so4,ibin)  = 1.e-9*aerSO4/water_a(ibin)
      ma(ja_hso4,ibin) = 1.e-9*aerHSO4/water_a(ibin)
      ma(ja_no3,ibin)  = 1.e-9*aer(ino3_a,jliquid,ibin)/water_a(ibin)
      ma(ja_cl,ibin)   = 1.e-9*aer(icl_a, jliquid,ibin)/water_a(ibin)

      mc(jc_h,ibin)    = mH
      mc(jc_ca,ibin)   = 1.e-9*aer(ica_a, jliquid,ibin)/water_a(ibin)
      mc(jc_nh4,ibin)  = 1.e-9*aer(inh4_a,jliquid,ibin)/water_a(ibin)
      mc(jc_na,ibin)   = 1.e-9*aer(ina_a, jliquid,ibin)/water_a(ibin)


! update the following activities
      activity(jhcl,ibin)    = mc(jc_h,ibin)*ma(ja_cl,ibin)   *   &
                               gam(jhcl,ibin)**2

      activity(jhno3,ibin)   = mc(jc_h,ibin)*ma(ja_no3,ibin)  *   &
                               gam(jhno3,ibin)**2

      activity(jnh4no3,ibin) = mc(jc_nh4,ibin)*ma(ja_no3,ibin)*   &
                               gam(jnh4no3,ibin)**2

      activity(jnh4cl,ibin)  = mc(jc_nh4,ibin)*ma(ja_cl,ibin) *   &
                               gam(jnh4cl,ibin)**2


! also update xyz(jtotal)
      aer(icl_a,jtotal,ibin)  = aer(icl_a,jliquid,ibin) +   &
                                aer(icl_a,jsolid,ibin)

      aer(ino3_a,jtotal,ibin) = aer(ino3_a,jliquid,ibin) +   &
                                aer(ino3_a,jsolid,ibin)

      electrolyte(jhno3,jtotal,ibin) = electrolyte(jhno3,jliquid,ibin)
      electrolyte(jhcl, jtotal,ibin) = electrolyte(jhcl, jliquid,ibin)

      return
      end subroutine equilibrate_hcl_and_hno3













!***********************************************************************
! called only once per entire simulation to load gas and aerosol
! indices, parameters, physico-chemical constants, polynomial coeffs, etc.
!
! author: Rahul A. Zaveri
! update: jan 2005
!-----------------------------------------------------------------------
      subroutine load_mosaic_parameters
!      include 'v33com2'
      use module_data_mosaic_main, only:  ipmcmos
      use module_data_mosaic_aero

      implicit none

! local variables
      integer iaer, je, ja, j_index, ibin
      logical first
      save first
      data first/.true./




      if(first)then
        first=.false.

!----------------------------------------------------------------
! control settings
! *** do not change mSIZE_FRAMEWORK here ***
!      mSIZE_FRAMEWORK = mSECTIONAL	! mMODAL or mSECTIONAL
!      mDYNAMIC_SOLVER = mASTEM	! mASTEM, mLSODES
!      mGAS_AER_XFER   = mON		! mON, mOFF

! ASTEM parameters
      nmax_ASTEM      = 301		! max number of time steps in ASTEM
!      alpha_ASTEM     = 1.0		! choose a value between 0.01 and 1.0
!      rtol_eqb_ASTEM = 0.01		! equilibrium tolerance in ASTEM
!      ptol_mol_ASTEM = 0.01		! mol percent tolerance in ASTEM

! MESA parameters
      Nmax_MESA       = 80		! max number of iterations in MESA_PTC
      rtol_mesa       = 0.01		! MESA equilibrium tolerance
!----------------------------------------------------------------
!

! edit here, change aro, alk, etc. species to simpleSOM species

! set gas and aerosol indices
!
! gas (local): [ngas_volatile] = [ngas_ioa] + [ngas_soa]
      ih2so4_g	= 1	! ioa (inorganic aerosol)
      ihno3_g	= 2	! ioa
      ihcl_g	= 3	! ioa
      inh3_g	= 4	! ioa
      imsa_g	= 5	! ioa
      icn3_g	= 6	! soa (secondary organic aerosol)
      icn2_g	= 7	! soa
      icn1_g	= 8	! soa (secondary organic aerosol)
      ic0_g	= 9	! soa
      ic1_g	= 10	! soa
      ic2_g	= 11	! soa
      ic3_g	= 12	! soa
      ic4_g	= 13	! soa
      ic5_g	= 14	! soa
      ic6_g	= 15	! soa
      ic7_g	= 16	! soa
      ic8_g = 17    ! soa
      ic9_g = 18    ! soa
!      ico2_g	= 18	! currently not used
!
! aer (local): used for total species - [naer]
      iso4_a	=  1	! <-> ih2so4_g
      ino3_a	=  2	! <-> ihno3_g
      icl_a	=  3	! <-> ihcl_g
      inh4_a	=  4	! <-> inh3_g
      imsa_a	=  5	! <-> imsa_g
      icn3_a	=  6	! <-> icn3_g
      icn2_a	=  7	! <-> icn2_g
      icn1_a	=  8	! <-> icn1_g
      ic0_a	=  9	! <-> ic0_g
      ic1_a	= 10	! <-> ic1_g
      ic2_a	= 11	! <-> ic2_g
      ic3_a	= 12	! <-> ic3_g
      ic4_a	= 13	! <-> ic4_g
      ic5_a	= 14	! <-> ic5_g
      ic6_a	= 15	! <-> ic6_g
      ic7_a	= 16	! <-> ic7_g
      ic8_a = 17    ! <-> ic8_g
      ic9_a = 18    ! <-> ic9_g
      ico3_a	= 19	! <-> ico2_g
      ina_a	= 20
      ica_a	= 21
      ioin_a = 22
      ioc_a	= 23
      ibc_a	= 24


! electrolyte indices (used for water content calculations) - [naercomp]
! these indices are order sensitive
      jnh4so4	=  1	! soluble
      jlvcite	=  2	! soluble
      jnh4hso4	=  3	! soluble
      jnh4msa	=  4	! soluble: new
      jnh4no3	=  5	! soluble
      jnh4cl	=  6	! soluble
      jna2so4	=  7	! soluble
      jna3hso4	=  8	! soluble
      jnahso4	=  9	! soluble
      jnamsa	= 10	! soluble: new
      jnano3	= 11	! soluble
      jnacl	= 12	! soluble
      jcano3	= 13	! soluble
      jcacl2	= 14	! soluble
      jcamsa2	= 15	! soluble		[nsalt]
      jh2so4	= 16	! soluble
      jmsa	= 17	! soluble
      jhno3	= 18	! soluble
      jhcl	= 19	! soluble
      jhhso4	= 20	! soluble		[nsoluble]
      jcaso4	= 21	! insoluble
      jcaco3	= 22	! insoluble		[nelectrolyte]
      joc	= 23	! insoluble - part of naercomp
      jbc	= 24	! insoluble - part of naercomp
      join	= 25	! insoluble - part of naercomp
      jcn3	= 26	! insoluble - part of naercomp
      jcn2	= 27	! insoluble - part of naercomp
      jcn1	= 28	! insoluble - part of naercomp
      jc0	= 29	! insoluble - part of naercomp
      jc1	= 30	! insoluble - part of naercomp
      jc2	= 31	! insoluble - part of naercomp
      jc3 	= 32	! insoluble - part of naercomp
      jc4	= 33	! insoluble - part of naercomp
      jc5	= 34	! insoluble - part of naercomp
      jc6	= 35	! insoluble - part of naercomp
      jc7	= 36	! insoluble - part of naercomp
      jc8   = 37    ! insoluble - part of naercomp
      jc9   = 38    ! insoluble - part of naercomp
      jh2o	= 39	! water - part of naercomp


! local aerosol ions
! cations - [ncation]
      jc_h	=  1
      jc_nh4	=  2
      jc_na	=  3
      jc_ca	=  4
!
! anions - [nanion]
      ja_hso4	=  1
      ja_so4  	=  2
      ja_no3  	=  3
      ja_cl   	=  4
      ja_msa	=  5
!     ja_co3	=  6

!--------------------------------------------------------------------
! phase state names
      phasestate(no_aerosol) = "NOAERO"
      phasestate(all_solid)  = "SOLID "
      phasestate(all_liquid) = "LIQUID"
      phasestate(mixed)      = "MIXED "

! names of aer species
      aer_name(iso4_a) = "SO4"
      aer_name(ino3_a) = "NO3"
      aer_name(icl_a)  = "Cl "
      aer_name(inh4_a) = "NH4"
      aer_name(ioc_a)  = "OC "
      aer_name(imsa_a) = "MSA"
      aer_name(ico3_a) = "CO3"
      aer_name(ina_a)  = "Na "
      aer_name(ica_a)  = "Ca "
      aer_name(ibc_a)  = "BC "
      aer_name(ioin_a) = "OIN"
	  aer_name(icn3_a) = "cn3"
	  aer_name(icn2_a) = "cn2"
	  aer_name(icn1_a) = "cn1"
	  aer_name(ic0_a)  = "c0 "
	  aer_name(ic1_a)  = "c1 "
	  aer_name(ic2_a)  = "c2 "
	  aer_name(ic3_a)  = "c3 "
	  aer_name(ic4_a)  = "c4 "
	  aer_name(ic5_a)  = "c5 "
	  aer_name(ic6_a)  = "c6 "
	  aer_name(ic7_a)  = "c7 "
	  aer_name(ic8_a)  = "c8 "
	  aer_name(ic9_a)  = "c9 "

! names of gas species
      gas_name(ih2so4_g) = "H2SO4"
      gas_name(ihno3_g)  = "HNO3 "
      gas_name(ihcl_g)   = "HCl  "
      gas_name(inh3_g)   = "NH3  "
      gas_name(imsa_g)   = "MSA  "
	  gas_name(icn3_g)   = "cn3  "
	  gas_name(icn2_g)   = "cn2  "
	  gas_name(icn1_g)   = "cn1  "
	  gas_name(ic0_g)    = "c0   "
	  gas_name(ic1_g)    = "c1   "
	  gas_name(ic2_g)    = "c2   "
	  gas_name(ic3_g)    = "c3   "
	  gas_name(ic4_g)    = "c4   "
	  gas_name(ic5_g)    = "c5   "
	  gas_name(ic6_g)    = "c6   "
	  gas_name(ic7_g)    = "c7   "
	  gas_name(ic8_g)    = "c8   "
	  gas_name(ic9_g)    = "c9   "

! names of electrolytes
      ename(jnh4so4) = "AmSO4"
      ename(jlvcite) = "(NH4)3H(SO4)2"
      ename(jnh4hso4)= "NH4HSO4"
      ename(jnh4msa) = "CH3SO3NH4"
      ename(jnh4no3) = "NH4NO3"
      ename(jnh4cl)  = "NH4Cl"
      ename(jnacl)   = "NaCl"
      ename(jnano3)  = "NaNO3"
      ename(jna2so4) = "Na2SO4"
      ename(jna3hso4)= "Na3H(SO4)2"
      ename(jnamsa)  = "CH3SO3Na"
      ename(jnahso4) = "NaHSO4"
      ename(jcaso4)  = "CaSO4"
      ename(jcamsa2) = "(CH3SO3)2Ca"
      ename(jcano3)  = "Ca(NO3)2"
      ename(jcacl2)  = "CaCl2"
      ename(jcaco3)  = "CaCO3"
      ename(jh2so4)  = "H2SO4"
      ename(jhhso4)  = "HHSO4"
      ename(jhno3)   = "HNO3"
      ename(jhcl)    = "HCl"
      ename(jmsa)    = "CH3SO3H"

! molecular weights of electrolytes
      mw_electrolyte(jnh4so4) = 132.0
      mw_electrolyte(jlvcite) = 247.0
      mw_electrolyte(jnh4hso4)= 115.0
      mw_electrolyte(jnh4msa) = 113.0
      mw_electrolyte(jnh4no3) = 80.0
      mw_electrolyte(jnh4cl)  = 53.5
      mw_electrolyte(jnacl)   = 58.5
      mw_electrolyte(jnano3)  = 85.0
      mw_electrolyte(jna2so4) = 142.0
      mw_electrolyte(jna3hso4)= 262.0
      mw_electrolyte(jnahso4) = 120.0
      mw_electrolyte(jnamsa)  = 118.0
      mw_electrolyte(jcaso4)  = 136.0
      mw_electrolyte(jcamsa2) = 230.0
      mw_electrolyte(jcano3)  = 164.0
      mw_electrolyte(jcacl2)  = 111.0
      mw_electrolyte(jcaco3)  = 100.0
      mw_electrolyte(jh2so4)  = 98.0
      mw_electrolyte(jhno3)   = 63.0
      mw_electrolyte(jhcl)    = 36.5
      mw_electrolyte(jmsa)    = 96.0
      mw_electrolyte(joc)     = 1.0
      mw_electrolyte(jbc)     = 1.0
      mw_electrolyte(join)    = 1.0
	  mw_electrolyte(jcn3)    = 238.0! wkc
	  mw_electrolyte(jcn2)    = 224.0! wkc
	  mw_electrolyte(jcn1)    = 210.0! wkc
	  mw_electrolyte(jc0)     = 196.0! wkc
	  mw_electrolyte(jc1)     = 182.0! wkc
	  mw_electrolyte(jc2)     = 168.0! wkc
	  mw_electrolyte(jc3)     = 154.0! wkc
	  mw_electrolyte(jc4)     = 140.0! wkc
	  mw_electrolyte(jc5)     = 126.0! wkc
	  mw_electrolyte(jc6)     = 112.0! wkc
	  mw_electrolyte(jc7)     = 98.0! wkc
	  mw_electrolyte(jc8)     = 84.0! wkc
	  mw_electrolyte(jc9)     = 70.0! wkc
      mw_electrolyte(jh2o)    = 18.0


! molecular weights of ions [g/mol]
      MW_c(jc_h)  =  1.0
      MW_c(jc_nh4)= 18.0
      MW_c(jc_na) = 23.0
      MW_c(jc_ca) = 40.0

      MW_a(ja_so4) = 96.0
      MW_a(ja_hso4)= 97.0
      MW_a(ja_no3) = 62.0
      MW_a(ja_cl)  = 35.5
      MW_a(ja_msa) = 95.0


! magnitude of the charges on ions
      zc(jc_h)   = 1
      zc(jc_nh4) = 1
      zc(jc_na)  = 1
      zc(jc_ca)  = 2

      za(ja_hso4)= 1
      za(ja_so4) = 2
      za(ja_no3) = 1
      za(ja_cl)  = 1
      za(ja_msa) = 1


! densities of pure electrolytes in g/cc
      dens_electrolyte(jnh4so4)  = 1.8
      dens_electrolyte(jlvcite)  = 1.8
      dens_electrolyte(jnh4hso4) = 1.8
      dens_electrolyte(jnh4msa)  = 1.8 ! assumed same as nh4hso4
      dens_electrolyte(jnh4no3)  = 1.8
      dens_electrolyte(jnh4cl)   = 1.8
      dens_electrolyte(jnacl)    = 2.2
      dens_electrolyte(jnano3)   = 2.2
      dens_electrolyte(jna2so4)  = 2.2
      dens_electrolyte(jna3hso4) = 2.2
      dens_electrolyte(jnahso4)  = 2.2
      dens_electrolyte(jnamsa)   = 2.2 ! assumed same as nahso4
      dens_electrolyte(jcaso4)   = 2.6
      dens_electrolyte(jcamsa2)  = 2.6	! assumed same as caso4
      dens_electrolyte(jcano3)   = 2.6
      dens_electrolyte(jcacl2)   = 2.6
      dens_electrolyte(jcaco3)   = 2.6
      dens_electrolyte(jh2so4)   = 1.8
      dens_electrolyte(jhhso4)   = 1.8
      dens_electrolyte(jhno3)    = 1.8
      dens_electrolyte(jhcl)     = 1.8
      dens_electrolyte(jmsa)     = 1.8 ! assumed same as h2so4

!      do je = 1, nelectrolyte
!        dens_electrolyte(je) = 1.6
!      enddo

! densities of compounds in g/cc
      dens_comp_a(jnh4so4)  = 1.8
      dens_comp_a(jlvcite)  = 1.8
      dens_comp_a(jnh4hso4) = 1.8
      dens_comp_a(jnh4msa)  = 1.8	! assumed same as nh4hso4
      dens_comp_a(jnh4no3)  = 1.7
      dens_comp_a(jnh4cl)   = 1.5
      dens_comp_a(jnacl)    = 2.2
      dens_comp_a(jnano3)   = 2.2
      dens_comp_a(jna2so4)  = 2.2
      dens_comp_a(jna3hso4) = 2.2
      dens_comp_a(jnahso4)  = 2.2
      dens_comp_a(jnamsa)   = 2.2	! assumed same as nahso4
      dens_comp_a(jcaso4)   = 2.6
      dens_comp_a(jcamsa2)  = 2.6	! assumed same as caso4
      dens_comp_a(jcano3)   = 2.6
      dens_comp_a(jcacl2)   = 2.6
      dens_comp_a(jcaco3)   = 2.6
      dens_comp_a(jh2so4)   = 1.8
      dens_comp_a(jhhso4)   = 1.8
      dens_comp_a(jhno3)    = 1.8
      dens_comp_a(jhcl)     = 1.8
      dens_comp_a(jmsa)     = 1.8	! assumed same as h2so4
      dens_comp_a(joc)      = 1.4 ! 1.0
      dens_comp_a(jbc)      = 1.8
      dens_comp_a(join)     =1.0! wkc remove!= 2.6
      dens_comp_a(jcn3)    = 1.0! wkc remove!= 1.4
      dens_comp_a(jcn2)    = 1.0! wkc remove!= 1.4
      dens_comp_a(jcn1)    = 1.0! wkc remove!= 1.4
      dens_comp_a(jc0)    = 1.0! wkc remove!= 1.4
      dens_comp_a(jc1)    = 1.0! wkc remove!= 1.4
      dens_comp_a(jc2)    = 1.0! wkc remove!= 1.4
      dens_comp_a(jc3)    = 1.0! wkc remove!= 1.4 ! 1.0
      dens_comp_a(jc4)    = 1.0! wkc remove!= 1.4 ! 1.0
      dens_comp_a(jc5)    = 1.0! wkc remove!= 1.4
      dens_comp_a(jc6)    = 1.0! wkc remove!= 1.4
      dens_comp_a(jc7)    = 1.0! wkc remove!= 1.4
      dens_comp_a(jc8)    = 1.0! wkc remove!= 1.4
      dens_comp_a(jc9)    = 1.0! wkc remove!= 1.4
      dens_comp_a(jh2o)     = 1.0

!      do je = 1, naercomp
!        dens_comp_a(je) = 1.6
!      enddo

!      do je = jcn3, jc7
!        dens_comp_a(je) = 1.0
!      enddo

! mass accommodation coefficients
      accom(ih2so4_g)  = 0.1
      accom(ihno3_g)   = 0.1
      accom(ihcl_g)    = 0.1
      accom(inh3_g)    = 0.1
      accom(imsa_g)    = 0.1
      accom(icn3_g)   = 1.0
      accom(icn2_g)   = 1.0
      accom(icn1_g)   = 1.0
      accom(ic0_g)   = 1.0
      accom(ic1_g)   = 1.0
      accom(ic2_g)   = 1.0
      accom(ic3_g)   = 1.0
      accom(ic4_g)   = 1.0
      accom(ic5_g)   = 1.0
      accom(ic6_g)   = 1.0
      accom(ic7_g)   = 1.0
      accom(ic8_g)   = 1.0
      accom(ic9_g)   = 1.0

! molecular weights of generic gas species
      mw_gas_mac(ih2so4_g)= 98.0
      mw_gas_mac(ihno3_g) = 63.0
      mw_gas_mac(ihcl_g)  = 36.5
      mw_gas_mac(inh3_g) = 17.0
      mw_gas_mac(imsa_g) = 96.0	! CH3SO3
      mw_gas_mac(icn3_g) = 238.0
      mw_gas_mac(icn2_g) = 224.0
      mw_gas_mac(icn1_g) = 210.0
      mw_gas_mac(ic0_g) = 196.0
      mw_gas_mac(ic1_g) = 182.0
      mw_gas_mac(ic2_g) = 168.0
      mw_gas_mac(ic3_g) = 154.0
      mw_gas_mac(ic4_g) = 140.0
      mw_gas_mac(ic5_g) = 126.0
      mw_gas_mac(ic6_g) = 112.0
      mw_gas_mac(ic7_g) = 98.0
      mw_gas_mac(ic8_g) = 84.0
      mw_gas_mac(ic9_g) = 70.0

! molecular weights of generic aerosol species
      mw_aer_mac(iso4_a) = 96.0
      mw_aer_mac(ino3_a) = 62.0
      mw_aer_mac(icl_a)  = 35.5
      mw_aer_mac(imsa_a) = 95.0	! CH3SO3
      mw_aer_mac(ico3_a) = 60.0
      mw_aer_mac(inh4_a) = 18.0
      mw_aer_mac(ina_a)  = 23.0
      mw_aer_mac(ica_a)  = 40.0
      mw_aer_mac(ioin_a) = 1.0		! not used
      mw_aer_mac(ibc_a)  = 1.0		! not used
      mw_aer_mac(ioc_a)  = 1.0	! 200 assumed for primary organics
      mw_aer_mac(icn3_a) = 238.0
      mw_aer_mac(icn2_a) = 224.0
      mw_aer_mac(icn1_a) = 210.0
      mw_aer_mac(ic0_a) = 196.0
      mw_aer_mac(ic1_a) = 182.0
      mw_aer_mac(ic2_a) = 168.0
      mw_aer_mac(ic3_a) = 154.0
      mw_aer_mac(ic4_a) = 140.0
      mw_aer_mac(ic5_a) = 126.0
      mw_aer_mac(ic6_a) = 112.0
      mw_aer_mac(ic7_a) = 98.0
      mw_aer_mac(ic8_a) = 84.0
      mw_aer_mac(ic9_a) = 70.0

! molecular weights of compounds
      mw_comp_a(jnh4so4) = 132.0
      mw_comp_a(jlvcite) = 247.0
      mw_comp_a(jnh4hso4)= 115.0
      mw_comp_a(jnh4msa) = 113.0
      mw_comp_a(jnh4no3) = 80.0
      mw_comp_a(jnh4cl)  = 53.5
      mw_comp_a(jnacl)   = 58.5
      mw_comp_a(jnano3)  = 85.0
      mw_comp_a(jna2so4) = 142.0
      mw_comp_a(jna3hso4)= 262.0
      mw_comp_a(jnahso4) = 120.0
      mw_comp_a(jnamsa)  = 118.0
      mw_comp_a(jcaso4)  = 136.0
      mw_comp_a(jcamsa2) = 230.0
      mw_comp_a(jcano3)  = 164.0
      mw_comp_a(jcacl2)  = 111.0
      mw_comp_a(jcaco3)  = 100.0
      mw_comp_a(jh2so4)  = 98.0
      mw_comp_a(jhhso4)  = 98.0
      mw_comp_a(jhno3)   = 63.0
      mw_comp_a(jhcl)    = 36.5
      mw_comp_a(jmsa)    = 96.0
      mw_comp_a(joc)	 = 1.0
      mw_comp_a(jbc)	 = 1.0
      mw_comp_a(join)    = 1.0
      mw_comp_a(jcn3)	 = 238.0 !150.0 wkc
      mw_comp_a(jcn2)	 = 224.0 !150.0 wkc
      mw_comp_a(jcn1)	 = 210.0 !150.0 wkc
      mw_comp_a(jc0)	 = 196.0 !150.0 wkc
      mw_comp_a(jc1)	 = 182.0 !140.0 wkc
      mw_comp_a(jc2)	 = 168.0 !140.0 wkc
      mw_comp_a(jc3)	 = 154.0 !152.0 wkc
      mw_comp_a(jc4)	 = 140.0 !168.0 wkc
      mw_comp_a(jc5)	 = 126.0 !136.0 wkc
      mw_comp_a(jc6)	 = 112.0 !272.0 wkc
      mw_comp_a(jc7)	 = 98.0  !136.0 wkc
      mw_comp_a(jc8)	 = 84.0
      mw_comp_a(jc9)	 = 70.0
      mw_comp_a(jh2o)    = 18.0

! densities of generic aerosol species
      dens_aer_mac(iso4_a) = 1.8	! used
      dens_aer_mac(ino3_a) = 1.8	! used
      dens_aer_mac(icl_a)  = 2.2	! used
      dens_aer_mac(imsa_a) = 1.8	! used
      dens_aer_mac(ico3_a) = 2.6	! used
      dens_aer_mac(inh4_a) = 1.8	! used
      dens_aer_mac(ina_a)  = 2.2	! used
      dens_aer_mac(ica_a)  = 2.6	! used
      dens_aer_mac(ioin_a) = 2.6	! used
      dens_aer_mac(ioc_a)  = 1.0 !1.4 ! 1.0	! used
      dens_aer_mac(ibc_a)  = 1.8	! used
      dens_aer_mac(icn3_a) = 1.0! wkc remove!= 1.4
      dens_aer_mac(icn2_a) = 1.0! wkc remove!= 1.4
      dens_aer_mac(icn1_a) = 1.0! wkc remove!= 1.4
      dens_aer_mac(ic0_a) = 1.0! wkc remove!= 1.4
      dens_aer_mac(ic1_a) = 1.0! wkc remove!= 1.4
      dens_aer_mac(ic2_a) = 1.0! wkc remove!= 1.4
      dens_aer_mac(ic3_a) = 1.0! wkc remove!= 1.4 ! 1.0
      dens_aer_mac(ic4_a) = 1.0! wkc remove!= 1.4 ! 1.0
      dens_aer_mac(ic5_a) = 1.0! wkc remove!= 1.4
      dens_aer_mac(ic6_a) = 1.0! wkc remove!= 1.4
      dens_aer_mac(ic7_a) = 1.0! wkc remove!= 1.4
      dens_aer_mac(ic8_a) = 1.0! wkc remove!= 1.4
      dens_aer_mac(ic9_a) = 1.0! wkc remove!= 1.4

!      do iaer = icn3_a, ilc7_a  ! naer
!        dens_aer_mac(iaer) = 1.0
!      enddo

      if (ipmcmos > 0) then
	  ! comment wkc 20181108
	  ! be wary of this commented out section of partmc-mosaic densities.
	  ! I will not be updating this part with new species
! use partmc-mosaic densities
!         dens_aer_mac(1:naer) = (/ &
!            1.80, 1.80, 2.20, 1.80, 1.80, 1.40, 1.40, 1.40, 1.40, 1.40, 1.40, &
!            1.40, 1.40, 1.40, 1.40, 1.40, 1.40, 2.60, 2.20, 2.60, 2.60, 1.40, 1.80 /)
!!           so4   no3   cl    nh4   msa   cn3  cn2  cn1  c0  c1  c2
!!           c3  c4  c5  c6  c7  co3   na    ca    oin   oc    bc

      dens_aer_mac(iso4_a) = 1.8	! used
      dens_aer_mac(ino3_a) = 1.8	! used
      dens_aer_mac(icl_a)  = 2.2	! used
      dens_aer_mac(imsa_a) = 1.8	! used
      dens_aer_mac(ico3_a) = 2.6	! used
      dens_aer_mac(inh4_a) = 1.8	! used
      dens_aer_mac(ina_a)  = 2.2	! used
      dens_aer_mac(ica_a)  = 2.6	! used
      dens_aer_mac(ioin_a) = 2.6	! used
      dens_aer_mac(ioc_a)  = 1.4 ! 1.0	! used
      dens_aer_mac(ibc_a)  = 1.8	! used
      dens_aer_mac(icn3_a) = 1.0! wkc remove!= 1.4
      dens_aer_mac(icn2_a) = 1.0! wkc remove!= 1.4
      dens_aer_mac(icn1_a) = 1.0! wkc remove!= 1.4
      dens_aer_mac(ic0_a) = 1.0! wkc remove!= 1.4
      dens_aer_mac(ic1_a) = 1.0! wkc remove!= 1.4
      dens_aer_mac(ic2_a) = 1.0! wkc remove!= 1.4
      dens_aer_mac(ic3_a) = 1.0! wkc remove!= 1.4 ! 1.0
      dens_aer_mac(ic4_a) = 1.0! wkc remove!= 1.4 ! 1.0
      dens_aer_mac(ic5_a) = 1.0! wkc remove!= 1.4
      dens_aer_mac(ic6_a) = 1.0! wkc remove!= 1.4
      dens_aer_mac(ic7_a) = 1.0! wkc remove!= 1.4
      dens_aer_mac(ic8_a) = 1.0! wkc remove!= 1.4
      dens_aer_mac(ic9_a) = 1.0! wkc remove!= 1.4

      end if

! kappa values (hygroscopicities) of generic aerosol species
!
! for calculation of ccn properties, kappa of electrolytes
!    should be used
! the multi-dimensional sectional code needs a "fixed" kappa
!    for each generic aerosol species, just as the older
!    1d sectional code needs a "fixed" dry density
      kappa_aer_mac(iso4_a)  = 0.65
      kappa_aer_mac(ino3_a)  = 0.65
      kappa_aer_mac(imsa_a)  = 0.65
      kappa_aer_mac(inh4_a)  = 0.65
      kappa_aer_mac(icl_a)   = 0.65
      kappa_aer_mac(ina_a)   = 0.65
      kappa_aer_mac(ico3_a)  = 0.001  ! ??
      kappa_aer_mac(ica_a)   = 0.001  ! ??
      kappa_aer_mac(ioin_a)  = 0.001
      kappa_aer_mac(ioc_a)   = 0.001
      kappa_aer_mac(ibc_a)   = 0.001
      kappa_aer_mac(icn3_a) = 1.0!0.1
      kappa_aer_mac(icn2_a) = 1.0!0.1
      kappa_aer_mac(icn1_a) = 1.0!0.1
      kappa_aer_mac(ic0_a) = 1.0!0.1
      kappa_aer_mac(ic1_a) = 1.0!0.1
      kappa_aer_mac(ic2_a) = 1.0!0.1
      kappa_aer_mac(ic3_a) = 1.0!0.1
      kappa_aer_mac(ic4_a) = 1.0!0.1
      kappa_aer_mac(ic5_a) = 1.0!0.1
      kappa_aer_mac(ic6_a) = 1.0!0.1
      kappa_aer_mac(ic7_a) = 1.0!0.1
      kappa_aer_mac(ic8_a) = 1.0!0.1
      kappa_aer_mac(ic9_a) = 1.0!0.1

! partial molar volumes of condensing species
      partial_molar_vol(ih2so4_g) = 51.83
      partial_molar_vol(ihno3_g)  = 31.45
      partial_molar_vol(ihcl_g)   = 20.96
      partial_molar_vol(inh3_g)   = 24.03
      partial_molar_vol(imsa_g)   = 53.33
      partial_molar_vol(icn3_g)  = 238.0 ! wkc remove!= 238.0/1.4
      partial_molar_vol(icn2_g)  = 224.0 ! wkc remove!= 224.0/1.4
      partial_molar_vol(icn1_g)  = 210.0 ! wkc remove!= 210.0/1.4
      partial_molar_vol(ic0_g)  = 196.0 ! wkc remove!= 196.0/1.4
      partial_molar_vol(ic1_g)  = 182.0 ! wkc remove!= 1.0! wkc remove!= 182.0/1.4
      partial_molar_vol(ic2_g)  = 168.0 ! wkc remove!= 168.0/1.4
      partial_molar_vol(ic3_g)  = 154.0 ! wkc remove!= 154.0/1.4 !
      partial_molar_vol(ic4_g)  = 140.0 ! wkc remove!= 140.0/1.4 !
      partial_molar_vol(ic5_g)  = 126.0 ! wkc remove!= 126.0/1.4 !
      partial_molar_vol(ic6_g)  = 112.0 ! wkc remove!= 112.0/1.4 !
      partial_molar_vol(ic7_g)  = 98.0 ! wkc remove!= 98.0/1.4 !
      partial_molar_vol(ic8_g)  = 84.0 ! wkc remove!= 84.0/1.4 !
      partial_molar_vol(ic9_g)  = 70.0 ! wkc remove!= 70.0/1.4 !

! refractive index
      ref_index_a(jnh4so4) = cmplx(1.52,0.)
      ref_index_a(jlvcite) = cmplx(1.50,0.)
      ref_index_a(jnh4hso4)= cmplx(1.47,0.)
      ref_index_a(jnh4msa) = cmplx(1.50,0.)	! assumed
      ref_index_a(jnh4no3) = cmplx(1.50,0.)
      ref_index_a(jnh4cl)  = cmplx(1.50,0.)
      ref_index_a(jnacl)   = cmplx(1.45,0.)
      ref_index_a(jnano3)  = cmplx(1.50,0.)
      ref_index_a(jna2so4) = cmplx(1.50,0.)
      ref_index_a(jna3hso4)= cmplx(1.50,0.)
      ref_index_a(jnahso4) = cmplx(1.50,0.)
      ref_index_a(jnamsa)  = cmplx(1.50,0.)	! assumed
      ref_index_a(jcaso4)  = cmplx(1.56,0.006)
      ref_index_a(jcamsa2) = cmplx(1.56,0.006)	! assumed
      ref_index_a(jcano3)  = cmplx(1.56,0.006)
      ref_index_a(jcacl2)  = cmplx(1.52,0.006)
      ref_index_a(jcaco3)  = cmplx(1.68,0.006)
      ref_index_a(jh2so4)  = cmplx(1.43,0.)
      ref_index_a(jhhso4)  = cmplx(1.43,0.)
      ref_index_a(jhno3)   = cmplx(1.50,0.)
      ref_index_a(jhcl)    = cmplx(1.50,0.)
      ref_index_a(jmsa)    = cmplx(1.43,0.)	! assumed
      ref_index_a(joc)	   = cmplx(1.45,0.)
      ref_index_a(jbc)	   = cmplx(1.82,0.74)
      ref_index_a(join)    = cmplx(1.55,0.006)
      ref_index_a(jcn3)   = cmplx(1.45,0.)
      ref_index_a(jcn2)   = cmplx(1.45,0.)
      ref_index_a(jcn1)   = cmplx(1.45,0.)
      ref_index_a(jc0)   = cmplx(1.45,0.)
      ref_index_a(jc1)   = cmplx(1.45,0.)
      ref_index_a(jc2)   = cmplx(1.45,0.)
      ref_index_a(jc3)   = cmplx(1.45,0.)
      ref_index_a(jc4)   = cmplx(1.45,0.)
      ref_index_a(jc5)   = cmplx(1.45,0.)
      ref_index_a(jc6)   = cmplx(1.45,0.)
      ref_index_a(jc7)   = cmplx(1.45,0.)
      ref_index_a(jc8)   = cmplx(1.45,0.)
      ref_index_a(jc9)   = cmplx(1.45,0.)
      ref_index_a(jh2o)    = cmplx(1.33,0.)

! jsalt_index
      jsalt_index(jnh4so4) = 5		! AS
      jsalt_index(jlvcite) = 2		! LV
      jsalt_index(jnh4hso4)= 1		! AB
      jsalt_index(jnh4no3) = 2		! AN
      jsalt_index(jnh4cl)  = 1		! AC
      jsalt_index(jna2so4) = 60		! SS
      jsalt_index(jnahso4) = 10		! SB
      jsalt_index(jnano3)  = 40		! SN
      jsalt_index(jnacl)   = 10		! SC
      jsalt_index(jcano3)  = 120	! CN
      jsalt_index(jcacl2)  = 80		! CC
      jsalt_index(jnh4msa) = 0		! AM	zero for now
      jsalt_index(jnamsa)  = 0		! SM	zero for now
      jsalt_index(jcamsa2) = 0		! CM	zero for now

! Aerosol Indices
!  AC = 1, AN = 2, AS = 5, SC = 10, SN = 40, SS = 60, CC = 80, CN = 120,
!  AB = 1, LV = 2, SB = 10
!
! SULFATE-POOR DOMAIN
      jsulf_poor(1)   = 	1	! 	AC
      jsulf_poor(2)   = 	2	! 	AN
      jsulf_poor(5)   = 	3	! 	AS
      jsulf_poor(10)  = 	4	! 	SC
      jsulf_poor(40)  = 	5	! 	SN
      jsulf_poor(60)  = 	6	! 	SS
      jsulf_poor(80)  = 	7	! 	CC
      jsulf_poor(120) = 	8	! 	CN
      jsulf_poor(3)   = 	9	! 	AN + AC
      jsulf_poor(6)   = 	10	! 	AS + AC
      jsulf_poor(7)   = 	11	! 	AS + AN
      jsulf_poor(8)   =  	12	! 	AS + AN + AC
      jsulf_poor(11)  = 	13	! 	SC + AC
      jsulf_poor(41)  = 	14	! 	SN + AC
      jsulf_poor(42)  = 	15	! 	SN + AN
      jsulf_poor(43)  = 	16	! 	SN + AN + AC
      jsulf_poor(50)  = 	17	! 	SN + SC
      jsulf_poor(51)  = 	18	! 	SN + SC + AC
      jsulf_poor(61)  = 	19	! 	SS + AC
      jsulf_poor(62)  = 	20	! 	SS + AN
      jsulf_poor(63)  = 	21	! 	SS + AN + AC
      jsulf_poor(65)  = 	22	! 	SS + AS
      jsulf_poor(66)  = 	23	! 	SS + AS + AC
      jsulf_poor(67)  = 	24	! 	SS + AS + AN
      jsulf_poor(68)  = 	25	! 	SS + AS + AN + AC
      jsulf_poor(70)  = 	26	! 	SS + SC
      jsulf_poor(71)  = 	27	! 	SS + SC + AC
      jsulf_poor(100) = 	28	! 	SS + SN
      jsulf_poor(101) = 	29	! 	SS + SN + AC
      jsulf_poor(102) = 	30	! 	SS + SN + AN
      jsulf_poor(103) = 	31	! 	SS + SN + AN + AC
      jsulf_poor(110) = 	32	! 	SS + SN + SC
      jsulf_poor(111) = 	33	! 	SS + SN + SC + AC
      jsulf_poor(81)  = 	34	! 	CC + AC
      jsulf_poor(90)  = 	35	! 	CC + SC
      jsulf_poor(91)  = 	36	! 	CC + SC + AC
      jsulf_poor(121) = 	37	! 	CN + AC
      jsulf_poor(122) = 	38	! 	CN + AN
      jsulf_poor(123) = 	39	! 	CN + AN + AC
      jsulf_poor(130) = 	40	! 	CN + SC
      jsulf_poor(131) = 	41	! 	CN + SC + AC
      jsulf_poor(160) = 	42	! 	CN + SN
      jsulf_poor(161) = 	43	! 	CN + SN + AC
      jsulf_poor(162) = 	44	! 	CN + SN + AN
      jsulf_poor(163) = 	45	! 	CN + SN + AN + AC
      jsulf_poor(170) = 	46	! 	CN + SN + SC
      jsulf_poor(171) = 	47	! 	CN + SN + SC + AC
      jsulf_poor(200) = 	48	! 	CN + CC
      jsulf_poor(201) = 	49	! 	CN + CC + AC
      jsulf_poor(210) = 	50	! 	CN + CC + SC
      jsulf_poor(211) = 	51	! 	CN + CC + SC + AC
!
! SULFATE-RICH DOMAIN
      jsulf_rich(1)   = 	52	! 	AB
      jsulf_rich(2)   = 	53	! 	LV
      jsulf_rich(10)  = 	54	! 	SB
      jsulf_rich(3)   = 	55	! 	AB + LV
      jsulf_rich(7)   = 	56	! 	AS + LV
      jsulf_rich(70)  = 	57	! 	SS + SB
      jsulf_rich(62)  = 	58	! 	SS + LV
      jsulf_rich(67)  = 	59	! 	SS + AS + LV
      jsulf_rich(61)  = 	60	! 	SS + AB
      jsulf_rich(63)  = 	61	! 	SS + LV + AB
      jsulf_rich(11)  = 	62	! 	SB + AB
      jsulf_rich(71)  = 	63	! 	SS + SB + AB
      jsulf_rich(5)   = 	3	!	AS
      jsulf_rich(60)  = 	6	! 	SS
      jsulf_rich(65)  = 	22	! 	SS + AS



!
! polynomial coefficients for binary molality (used in ZSR equation)
!
!
! a_zsr for aw < 0.97
!
! (NH4)2SO4
      je = jnh4so4
      a_zsr(1,je)  =  1.30894
      a_zsr(2,je)  = -7.09922
      a_zsr(3,je)  =  20.62831
      a_zsr(4,je)  = -32.19965
      a_zsr(5,je)  =  25.17026
      a_zsr(6,je)  = -7.81632
      aw_min(je)   = 0.1
!
! (NH4)3H(SO4)2
      je = jlvcite
      a_zsr(1,je)  =  1.10725
      a_zsr(2,je)  = -5.17978
      a_zsr(3,je)  =  12.29534
      a_zsr(4,je)  = -16.32545
      a_zsr(5,je)  =  11.29274
      a_zsr(6,je)  = -3.19164
      aw_min(je)   = 0.1
!
! NH4HSO4
      je = jnh4hso4
      a_zsr(1,je)  =  1.15510
      a_zsr(2,je)  = -3.20815
      a_zsr(3,je)  =  2.71141
      a_zsr(4,je)  =  2.01155
      a_zsr(5,je)  = -4.71014
      a_zsr(6,je)  =  2.04616
      aw_min(je)   = 0.1
!
! NH4MSA (assumed same as NH4HSO4)
      je = jnh4msa
      a_zsr(1,je)  =  1.15510
      a_zsr(2,je)  = -3.20815
      a_zsr(3,je)  =  2.71141
      a_zsr(4,je)  =  2.01155
      a_zsr(5,je)  = -4.71014
      a_zsr(6,je)  =  2.04616
      aw_min(je)   = 0.1
!
! NH4NO3
      je = jnh4no3
      a_zsr(1,je)  =  0.43507
      a_zsr(2,je)  =  6.38220
      a_zsr(3,je)  = -30.19797
      a_zsr(4,je)  =  53.36470
      a_zsr(5,je)  = -43.44203
      a_zsr(6,je)  =  13.46158
      aw_min(je)   = 0.1
!
! NH4Cl: revised on Nov 13, 2003. based on Chan and Ha (1999) JGR.
      je = jnh4cl
      a_zsr(1,je)  =  0.45309
      a_zsr(2,je)  =  2.65606
      a_zsr(3,je)  = -14.7730
      a_zsr(4,je)  =  26.2936
      a_zsr(5,je)  = -20.5735
      a_zsr(6,je)  =  5.94255
      aw_min(je)   = 0.1
!
! NaCl
      je = jnacl
      a_zsr(1,je)  =  0.42922
      a_zsr(2,je)  = -1.17718
      a_zsr(3,je)  =  2.80208
      a_zsr(4,je)  = -4.51097
      a_zsr(5,je)  =  3.76963
      a_zsr(6,je)  = -1.31359
      aw_min(je)   = 0.1
!
! NaNO3
      je = jnano3
      a_zsr(1,je)  =  1.34966
      a_zsr(2,je)  = -5.20116
      a_zsr(3,je)  =  11.49011
      a_zsr(4,je)  = -14.41380
      a_zsr(5,je)  =  9.07037
      a_zsr(6,je)  = -2.29769
      aw_min(je)   = 0.1
!
! Na2SO4
      je = jna2so4
      a_zsr(1,je)  =  0.39888
      a_zsr(2,je)  = -1.27150
      a_zsr(3,je)  =  3.42792
      a_zsr(4,je)  = -5.92632
      a_zsr(5,je)  =  5.33351
      a_zsr(6,je)  = -1.96541
      aw_min(je)   = 0.1
!
! Na3H(SO4)2  added on 1/14/2004
      je = jna3hso4
      a_zsr(1,je)  =  0.31480
      a_zsr(2,je)  = -1.01087
      a_zsr(3,je)  =  2.44029
      a_zsr(4,je)  = -3.66095
      a_zsr(5,je)  =  2.77632
      a_zsr(6,je)  = -0.86058
      aw_min(je)   = 0.1
!
! NaHSO4
      je = jnahso4
      a_zsr(1,je)  =  0.62764
      a_zsr(2,je)  = -1.63520
      a_zsr(3,je)  =  4.62531
      a_zsr(4,je)  = -10.06925
      a_zsr(5,je)  =  10.33547
      a_zsr(6,je)  = -3.88729
      aw_min(je)   = 0.1
!
! NaMSA (assumed same as NaHSO4)
      je = jnamsa
      a_zsr(1,je)  =  0.62764
      a_zsr(2,je)  = -1.63520
      a_zsr(3,je)  =  4.62531
      a_zsr(4,je)  = -10.06925
      a_zsr(5,je)  =  10.33547
      a_zsr(6,je)  = -3.88729
      aw_min(je)   = 0.1
!
! Ca(NO3)2
      je = jcano3
      a_zsr(1,je)  =  0.38895
      a_zsr(2,je)  = -1.16013
      a_zsr(3,je)  =  2.16819
      a_zsr(4,je)  = -2.23079
      a_zsr(5,je)  =  1.00268
      a_zsr(6,je)  = -0.16923
      aw_min(je)   = 0.1
!
! CaCl2: Kim and Seinfeld
      je = jcacl2
      a_zsr(1,je)  =  0.29891
      a_zsr(2,je)  = -1.31104
      a_zsr(3,je)  =  3.68759
      a_zsr(4,je)  = -5.81708
      a_zsr(5,je)  =  4.67520
      a_zsr(6,je)  = -1.53223
      aw_min(je)   = 0.1
!
! H2SO4
      je = jh2so4
      a_zsr(1,je) =  0.32751
      a_zsr(2,je) = -1.00692
      a_zsr(3,je) =  2.59750
      a_zsr(4,je) = -4.40014
      a_zsr(5,je) =  3.88212
      a_zsr(6,je) = -1.39916
      aw_min(je)  = 0.1
!
! MSA (assumed same as H2SO4)
      je = jmsa
      a_zsr(1,je) =  0.32751
      a_zsr(2,je) = -1.00692
      a_zsr(3,je) =  2.59750
      a_zsr(4,je) = -4.40014
      a_zsr(5,je) =  3.88212
      a_zsr(6,je) = -1.39916
      aw_min(je)  = 0.1
!
! HHSO4
      je = jhhso4
      a_zsr(1,je) =  0.32751
      a_zsr(2,je) = -1.00692
      a_zsr(3,je) =  2.59750
      a_zsr(4,je) = -4.40014
      a_zsr(5,je) =  3.88212
      a_zsr(6,je) = -1.39916
      aw_min(je)  = 1.0
!
! HNO3
      je = jhno3
      a_zsr(1,je) =  0.75876
      a_zsr(2,je) = -3.31529
      a_zsr(3,je) =  9.26392
      a_zsr(4,je) = -14.89799
      a_zsr(5,je) =  12.08781
      a_zsr(6,je) = -3.89958
      aw_min(je)  = 0.1
!
! HCl
      je = jhcl
      a_zsr(1,je) =  0.31133
      a_zsr(2,je) = -0.79688
      a_zsr(3,je) =  1.93995
      a_zsr(4,je) = -3.31582
      a_zsr(5,je) =  2.93513
      a_zsr(6,je) = -1.07268
      aw_min(je)  = 0.1
!
! CaSO4
      je = jcaso4
      a_zsr(1,je)  =  0.0
      a_zsr(2,je)  =  0.0
      a_zsr(3,je)  =  0.0
      a_zsr(4,je)  =  0.0
      a_zsr(5,je)  =  0.0
      a_zsr(6,je)  =  0.0
      aw_min(je)   = 1.0
!
! Ca(MSA)2 (assumed same as Ca(NO3)2)
      je = jcamsa2
      a_zsr(1,je)  =  0.38895
      a_zsr(2,je)  = -1.16013
      a_zsr(3,je)  =  2.16819
      a_zsr(4,je)  = -2.23079
      a_zsr(5,je)  =  1.00268
      a_zsr(6,je)  = -0.16923
      aw_min(je)   = 0.1
!
! CaCO3
      je = jcaco3
      a_zsr(1,je)  =  0.0
      a_zsr(2,je)  =  0.0
      a_zsr(3,je)  =  0.0
      a_zsr(4,je)  =  0.0
      a_zsr(5,je)  =  0.0
      a_zsr(6,je)  =  0.0
      aw_min(je)   = 1.0



!-------------------------------------------
! b_zsr for aw => 0.97 to 0.99999
!
! (NH4)2SO4
      b_zsr(jnh4so4)  = 28.0811
!
! (NH4)3H(SO4)2
      b_zsr(jlvcite)  = 14.7178
!
! NH4HSO4
      b_zsr(jnh4hso4) = 29.4779
!
! NH4MSA
      b_zsr(jnh4msa)  = 29.4779 ! assumed same as NH4HSO4
!
! NH4NO3
      b_zsr(jnh4no3)  = 33.4049
!
! NH4Cl
      b_zsr(jnh4cl)   = 30.8888
!
! NaCl
      b_zsr(jnacl)    = 29.8375
!
! NaNO3
      b_zsr(jnano3)   = 32.2756
!
! Na2SO4
      b_zsr(jna2so4)  = 27.6889
!
! Na3H(SO4)2
      b_zsr(jna3hso4) = 14.2184
!
! NaHSO4
      b_zsr(jnahso4)  = 28.3367
!
! NaMSA
      b_zsr(jnamsa)   = 28.3367 ! assumed same as NaHSO4
!
! Ca(NO3)2
      b_zsr(jcano3)   = 18.3661
!
! CaCl2
      b_zsr(jcacl2)   = 20.8792
!
! H2SO4
      b_zsr(jh2so4)   = 26.7347
!
! HHSO4
      b_zsr(jhhso4)   = 26.7347
!
! HNO3
      b_zsr(jhno3)    = 28.8257
!
! HCl
      b_zsr(jhcl)     = 27.7108
!
! MSA
      b_zsr(jmsa)     = 26.7347 ! assumed same as H2SO4
!
! CaSO4
      b_zsr(jcaso4)   = 0.0
!
! Ca(MSA)2
      b_zsr(jcamsa2)  = 18.3661 ! assumed same as Ca(NO3)2
!
! CaCO3
      b_zsr(jcaco3)   = 0.0









!-------------------------------------------
! Li and Lu (2001) Surface tension model
! G_MX [mol/cm^2]; K_MX [-]
!
! (NH4)2SO4
      G_MX(jnh4so4)  = -8.79e-7*1.e-4
      K_MX(jnh4so4)  =  3.84e+1
!
! (NH4)3H(SO4)2
      G_MX(jlvcite)  = -8.79e-7*1.e-4	! assumed same as (NH4)2SO4
      K_MX(jlvcite)  =  3.84e+1		! assumed same as (NH4)2SO4
!
! NH4HSO4
      G_MX(jnh4hso4) = -8.79e-7*1.e-4	! assumed same as (NH4)2SO4
      K_MX(jnh4hso4) =  3.84e+1		! assumed same as (NH4)2SO4
!
! NH4MSA
      G_MX(jnh4msa)  = -8.79e-7*1.e-4	! assumed same as (NH4)2SO4
      K_MX(jnh4msa)  =  3.84e+1		! assumed same as (NH4)2SO4
!
! NH4NO3
      G_MX(jnh4no3)  = -3.08e-6*1.e-4
      K_MX(jnh4no3)  =  4.89e-1
!
! NH4Cl
      G_MX(jnh4cl)   = -1.01e-6*1.e-4
      K_MX(jnh4cl)   =  1.3
!
! NaCl
      G_MX(jnacl)    = -1.05e-6*1.e-4
      K_MX(jnacl)    =  1.2
!
! NaNO3
      G_MX(jnano3)   = -1.66e-6*1.e-4
      K_MX(jnano3)   =  1.25
!
! Na2SO4
      G_MX(jna2so4)  = -8.37e-7*1.e-4
      K_MX(jna2so4)  =  7.57e+1
!
! Na3H(SO4)2
      G_MX(jna3hso4) = -8.37e-7*1.e-4	! assumed same as Na2SO4
      K_MX(jna3hso4) =  7.57e+1		! assumed same as Na2SO4
!
! NaHSO4
      G_MX(jnahso4)  = -8.37e-7*1.e-4	! assumed same as Na2SO4
      K_MX(jnahso4)  =  7.57e+1		! assumed same as Na2SO4
!
! NaMSA
      G_MX(jnamsa)   = -8.37e-7*1.e-4
      K_MX(jnamsa)   =  7.57e+1
!
! Ca(NO3)2
      G_MX(jcano3)   = -4.88e-7*1.e-4	! assumed same as CaCl2
      K_MX(jcano3)   =  1.50e+1		! assumed same as CaCl2
!
! CaCl2
      G_MX(jcacl2)   = -4.88e-7*1.e-4
      K_MX(jcacl2)   =  1.50e+1
!
! H2SO4
      G_MX(jh2so4)   = -6.75e-8*1.e-4
      K_MX(jh2so4)   =  1.65e+3
!
! HHSO4
      G_MX(jh2so4)   = -6.75e-8*1.e-4	! assumed same as H2SO4
      K_MX(jh2so4)   =  1.65e+3		! assumed same as H2SO4
!
! HNO3
      G_MX(jhno3)    =  8.05e-7*1.e-4
      K_MX(jhno3)    =  1.06e-1
!
! HCl
      G_MX(jhcl)     =  4.12e-7*1.e-4
      K_MX(jhcl)     =  4.68e-3
!
! MSA
      G_MX(jmsa)     =  8.05e-7*1.e-4	! assumed same as HNO3
      K_MX(jmsa)     =  1.06e-1		! assumed same as HNO3
!
! CaSO4
      G_MX(jmsa)     =  0.0*1.e-4	! assumed
      K_MX(jmsa)     =  0.0		! assumed
!
! Ca(MSA)2
      G_MX(jcamsa2)  =  0.0*1.e-4	! assumed
      K_MX(jcamsa2)  =  0.0		! assumed
!
! CaCO3
      G_MX(jcaco3)   =  0.0*1.e-4	! assumed
      K_MX(jcaco3)   =  0.0		! assumed







!----------------------------------------------------------------
! parameters for MTEM mixing rule (Zaveri, Easter, and Wexler, 2005)
! log_gamZ(jA,jE)   A in E
!----------------------------------------------------------------
!
! (NH4)2SO4 in E
      jA = jnh4so4

! in (NH4)2SO4
      jE = jnh4so4
      b_mtem(1,jA,jE) = -2.94685
      b_mtem(2,jA,jE) = 17.3328
      b_mtem(3,jA,jE) = -64.8441
      b_mtem(4,jA,jE) = 122.7070
      b_mtem(5,jA,jE) = -114.4373
      b_mtem(6,jA,jE) = 41.6811

! in NH4NO3
      jE = jnh4no3
      b_mtem(1,jA,jE) = -2.7503
      b_mtem(2,jA,jE) = 4.3806
      b_mtem(3,jA,jE) = -1.1110
      b_mtem(4,jA,jE) = -1.7005
      b_mtem(5,jA,jE) = -4.4207
      b_mtem(6,jA,jE) = 5.1990

! in NH4Cl (revised on 11/15/2003)
      jE = jnh4cl
      b_mtem(1,jA,jE) = -2.06952
      b_mtem(2,jA,jE) = 7.1240
      b_mtem(3,jA,jE) = -24.4274
      b_mtem(4,jA,jE) = 51.1458
      b_mtem(5,jA,jE) = -54.2056
      b_mtem(6,jA,jE) = 22.0606

! in Na2SO4
      jE = jna2so4
      b_mtem(1,jA,jE) = -2.17361
      b_mtem(2,jA,jE) = 15.9919
      b_mtem(3,jA,jE) = -69.0952
      b_mtem(4,jA,jE) = 139.8860
      b_mtem(5,jA,jE) = -134.9890
      b_mtem(6,jA,jE) = 49.8877

! in NaNO3
      jE = jnano3
      b_mtem(1,jA,jE) = -4.4370
      b_mtem(2,jA,jE) = 24.0243
      b_mtem(3,jA,jE) = -76.2437
      b_mtem(4,jA,jE) = 128.6660
      b_mtem(5,jA,jE) = -110.0900
      b_mtem(6,jA,jE) = 37.7414

! in NaCl
      jE = jnacl
      b_mtem(1,jA,jE) = -1.5394
      b_mtem(2,jA,jE) = 5.8671
      b_mtem(3,jA,jE) = -22.7726
      b_mtem(4,jA,jE) = 47.0547
      b_mtem(5,jA,jE) = -47.8266
      b_mtem(6,jA,jE) = 18.8489

! in HNO3
      jE = jhno3
      b_mtem(1,jA,jE) = -0.35750
      b_mtem(2,jA,jE) = -3.82466
      b_mtem(3,jA,jE) = 4.55462
      b_mtem(4,jA,jE) = 5.05402
      b_mtem(5,jA,jE) = -14.7476
      b_mtem(6,jA,jE) = 8.8009

! in HCl
      jE = jhcl
      b_mtem(1,jA,jE) = -2.15146
      b_mtem(2,jA,jE) = 5.50205
      b_mtem(3,jA,jE) = -19.1476
      b_mtem(4,jA,jE) = 39.1880
      b_mtem(5,jA,jE) = -39.9460
      b_mtem(6,jA,jE) = 16.0700

! in H2SO4
      jE = jh2so4
      b_mtem(1,jA,jE) = -2.52604
      b_mtem(2,jA,jE) = 9.76022
      b_mtem(3,jA,jE) = -35.2540
      b_mtem(4,jA,jE) = 71.2981
      b_mtem(5,jA,jE) = -71.8207
      b_mtem(6,jA,jE) = 28.0758

!
! in NH4HSO4
      jE = jnh4hso4
      b_mtem(1,jA,jE) = -4.13219
      b_mtem(2,jA,jE) = 13.8863
      b_mtem(3,jA,jE) = -34.5387
      b_mtem(4,jA,jE) = 56.5012
      b_mtem(5,jA,jE) = -51.8702
      b_mtem(6,jA,jE) = 19.6232

!
! in (NH4)3H(SO4)2
      jE = jlvcite
      b_mtem(1,jA,jE) = -2.53482
      b_mtem(2,jA,jE) = 12.3333
      b_mtem(3,jA,jE) = -46.1020
      b_mtem(4,jA,jE) = 90.4775
      b_mtem(5,jA,jE) = -88.1254
      b_mtem(6,jA,jE) = 33.4715

!
! in NaHSO4
      jE = jnahso4
      b_mtem(1,jA,jE) = -3.23425
      b_mtem(2,jA,jE) = 18.7842
      b_mtem(3,jA,jE) = -78.7807
      b_mtem(4,jA,jE) = 161.517
      b_mtem(5,jA,jE) = -154.940
      b_mtem(6,jA,jE) = 56.2252

!
! in Na3H(SO4)2
      jE = jna3hso4
      b_mtem(1,jA,jE) = -1.25316
      b_mtem(2,jA,jE) = 7.40960
      b_mtem(3,jA,jE) = -34.8929
      b_mtem(4,jA,jE) = 72.8853
      b_mtem(5,jA,jE) = -72.4503
      b_mtem(6,jA,jE) = 27.7706


!-----------------
! NH4NO3 in E
      jA = jnh4no3

! in (NH4)2SO4
      jE = jnh4so4
      b_mtem(1,jA,jE) = -3.5201
      b_mtem(2,jA,jE) = 21.6584
      b_mtem(3,jA,jE) = -72.1499
      b_mtem(4,jA,jE) = 126.7000
      b_mtem(5,jA,jE) = -111.4550
      b_mtem(6,jA,jE) = 38.5677

! in NH4NO3
      jE = jnh4no3
      b_mtem(1,jA,jE) = -2.2630
      b_mtem(2,jA,jE) = -0.1518
      b_mtem(3,jA,jE) = 17.0898
      b_mtem(4,jA,jE) = -36.7832
      b_mtem(5,jA,jE) = 29.8407
      b_mtem(6,jA,jE) = -7.9314

! in NH4Cl (revised on 11/15/2003)
      jE = jnh4cl
      b_mtem(1,jA,jE) = -1.3851
      b_mtem(2,jA,jE) = -0.4462
      b_mtem(3,jA,jE) = 8.4567
      b_mtem(4,jA,jE) = -11.5988
      b_mtem(5,jA,jE) = 2.9802
      b_mtem(6,jA,jE) = 1.8132

! in Na2SO4
      jE = jna2so4
      b_mtem(1,jA,jE) = -1.7602
      b_mtem(2,jA,jE) = 10.4044
      b_mtem(3,jA,jE) = -35.5894
      b_mtem(4,jA,jE) = 64.3584
      b_mtem(5,jA,jE) = -57.8931
      b_mtem(6,jA,jE) = 20.2141

! in NaNO3
      jE = jnano3
      b_mtem(1,jA,jE) = -3.24346
      b_mtem(2,jA,jE) = 16.2794
      b_mtem(3,jA,jE) = -48.7601
      b_mtem(4,jA,jE) = 79.2246
      b_mtem(5,jA,jE) = -65.8169
      b_mtem(6,jA,jE) = 22.1500

! in NaCl
      jE = jnacl
      b_mtem(1,jA,jE) = -1.75658
      b_mtem(2,jA,jE) = 7.71384
      b_mtem(3,jA,jE) = -22.7984
      b_mtem(4,jA,jE) = 39.1532
      b_mtem(5,jA,jE) = -34.6165
      b_mtem(6,jA,jE) = 12.1283

! in Ca(NO3)2
      jE = jcano3
      b_mtem(1,jA,jE) = -0.97178
      b_mtem(2,jA,jE) = 6.61964
      b_mtem(3,jA,jE) = -26.2353
      b_mtem(4,jA,jE) = 50.5259
      b_mtem(5,jA,jE) = -47.6586
      b_mtem(6,jA,jE) = 17.5074

! in CaCl2 added on 12/22/2003
      jE = jcacl2
      b_mtem(1,jA,jE) = -0.41515
      b_mtem(2,jA,jE) = 6.44101
      b_mtem(3,jA,jE) = -26.4473
      b_mtem(4,jA,jE) = 49.0718
      b_mtem(5,jA,jE) = -44.2631
      b_mtem(6,jA,jE) = 15.3771

! in HNO3
      jE = jhno3
      b_mtem(1,jA,jE) = -1.20644
      b_mtem(2,jA,jE) = 5.70117
      b_mtem(3,jA,jE) = -18.2783
      b_mtem(4,jA,jE) = 31.7199
      b_mtem(5,jA,jE) = -27.8703
      b_mtem(6,jA,jE) = 9.7299

! in HCl
      jE = jhcl
      b_mtem(1,jA,jE) = -0.680862
      b_mtem(2,jA,jE) = 3.59456
      b_mtem(3,jA,jE) = -10.7969
      b_mtem(4,jA,jE) = 17.8434
      b_mtem(5,jA,jE) = -15.3165
      b_mtem(6,jA,jE) = 5.17123


!----------
! NH4Cl in E
      jA = jnh4cl

! in (NH4)2SO4
      jE = jnh4so4
      b_mtem(1,jA,jE) = -2.8850
      b_mtem(2,jA,jE) = 20.6970
      b_mtem(3,jA,jE) = -70.6810
      b_mtem(4,jA,jE) = 124.3690
      b_mtem(5,jA,jE) = -109.2880
      b_mtem(6,jA,jE) = 37.5831

! in NH4NO3
      jE = jnh4no3
      b_mtem(1,jA,jE) = -1.9386
      b_mtem(2,jA,jE) = 1.3238
      b_mtem(3,jA,jE) = 11.8500
      b_mtem(4,jA,jE) = -28.1168
      b_mtem(5,jA,jE) = 21.8543
      b_mtem(6,jA,jE) = -5.1671

! in NH4Cl (revised on 11/15/2003)
      jE = jnh4cl
      b_mtem(1,jA,jE) = -0.9559
      b_mtem(2,jA,jE) = 0.8121
      b_mtem(3,jA,jE) = 4.3644
      b_mtem(4,jA,jE) = -8.9258
      b_mtem(5,jA,jE) = 4.2362
      b_mtem(6,jA,jE) = 0.2891

! in Na2SO4
      jE = jna2so4
      b_mtem(1,jA,jE) = 0.0377
      b_mtem(2,jA,jE) = 6.0752
      b_mtem(3,jA,jE) = -30.8641
      b_mtem(4,jA,jE) = 63.3095
      b_mtem(5,jA,jE) = -61.0070
      b_mtem(6,jA,jE) = 22.1734

! in NaNO3
      jE = jnano3
      b_mtem(1,jA,jE) = -1.8336
      b_mtem(2,jA,jE) = 12.8160
      b_mtem(3,jA,jE) = -42.3388
      b_mtem(4,jA,jE) = 71.1816
      b_mtem(5,jA,jE) = -60.5708
      b_mtem(6,jA,jE) = 20.5853

! in NaCl
      jE = jnacl
      b_mtem(1,jA,jE) = -0.1429
      b_mtem(2,jA,jE) = 2.3561
      b_mtem(3,jA,jE) = -10.4425
      b_mtem(4,jA,jE) = 20.8951
      b_mtem(5,jA,jE) = -20.7739
      b_mtem(6,jA,jE) = 7.9355

! in Ca(NO3)2
      jE = jcano3
      b_mtem(1,jA,jE) = 0.76235
      b_mtem(2,jA,jE) = 3.08323
      b_mtem(3,jA,jE) = -23.6772
      b_mtem(4,jA,jE) = 53.7415
      b_mtem(5,jA,jE) = -55.4043
      b_mtem(6,jA,jE) = 21.2944

! in CaCl2 (revised on 11/27/2003)
      jE = jcacl2
      b_mtem(1,jA,jE) = 1.13864
      b_mtem(2,jA,jE) = -0.340539
      b_mtem(3,jA,jE) = -8.67025
      b_mtem(4,jA,jE) = 22.8008
      b_mtem(5,jA,jE) = -24.5181
      b_mtem(6,jA,jE) = 9.3663

! in HNO3
      jE = jhno3
      b_mtem(1,jA,jE) = 2.42532
      b_mtem(2,jA,jE) = -14.1755
      b_mtem(3,jA,jE) = 38.804
      b_mtem(4,jA,jE) = -58.2437
      b_mtem(5,jA,jE) = 43.5431
      b_mtem(6,jA,jE) = -12.5824

! in HCl
      jE = jhcl
      b_mtem(1,jA,jE) = 0.330337
      b_mtem(2,jA,jE) = 0.0778934
      b_mtem(3,jA,jE) = -2.30492
      b_mtem(4,jA,jE) = 4.73003
      b_mtem(5,jA,jE) = -4.80849
      b_mtem(6,jA,jE) = 1.78866


!----------
! Na2SO4 in E
      jA = jna2so4

! in (NH4)2SO4
      jE = jnh4so4
      b_mtem(1,jA,jE) = -2.6982
      b_mtem(2,jA,jE) = 22.9875
      b_mtem(3,jA,jE) = -98.9840
      b_mtem(4,jA,jE) = 198.0180
      b_mtem(5,jA,jE) = -188.7270
      b_mtem(6,jA,jE) = 69.0548

! in NH4NO3
      jE = jnh4no3
      b_mtem(1,jA,jE) = -2.4844
      b_mtem(2,jA,jE) = 6.5420
      b_mtem(3,jA,jE) = -9.8998
      b_mtem(4,jA,jE) = 11.3884
      b_mtem(5,jA,jE) = -13.6842
      b_mtem(6,jA,jE) = 7.7411

! in NH4Cl (revised on 11/15/2003)
      jE = jnh4cl
      b_mtem(1,jA,jE) = -1.3325
      b_mtem(2,jA,jE) = 13.0406
      b_mtem(3,jA,jE) = -56.1935
      b_mtem(4,jA,jE) = 107.1170
      b_mtem(5,jA,jE) = -97.3721
      b_mtem(6,jA,jE) = 34.3763

! in Na2SO4
      jE = jna2so4
      b_mtem(1,jA,jE) = -1.2832
      b_mtem(2,jA,jE) = 12.8526
      b_mtem(3,jA,jE) = -62.2087
      b_mtem(4,jA,jE) = 130.3876
      b_mtem(5,jA,jE) = -128.2627
      b_mtem(6,jA,jE) = 48.0340

! in NaNO3
      jE = jnano3
      b_mtem(1,jA,jE) = -3.5384
      b_mtem(2,jA,jE) = 21.3758
      b_mtem(3,jA,jE) = -70.7638
      b_mtem(4,jA,jE) = 121.1580
      b_mtem(5,jA,jE) = -104.6230
      b_mtem(6,jA,jE) = 36.0557

! in NaCl
      jE = jnacl
      b_mtem(1,jA,jE) = 0.2175
      b_mtem(2,jA,jE) = -0.5648
      b_mtem(3,jA,jE) = -8.0288
      b_mtem(4,jA,jE) = 25.9734
      b_mtem(5,jA,jE) = -32.3577
      b_mtem(6,jA,jE) = 14.3924

! in HNO3
      jE = jhno3
      b_mtem(1,jA,jE) = -0.309617
      b_mtem(2,jA,jE) = -1.82899
      b_mtem(3,jA,jE) = -1.5505
      b_mtem(4,jA,jE) = 13.3847
      b_mtem(5,jA,jE) = -20.1284
      b_mtem(6,jA,jE) = 9.93163

! in HCl
      jE = jhcl
      b_mtem(1,jA,jE) = -0.259455
      b_mtem(2,jA,jE) = -0.819366
      b_mtem(3,jA,jE) = -4.28964
      b_mtem(4,jA,jE) = 16.4305
      b_mtem(5,jA,jE) = -21.8546
      b_mtem(6,jA,jE) = 10.3044

! in H2SO4
      jE = jh2so4
      b_mtem(1,jA,jE) = -1.84257
      b_mtem(2,jA,jE) = 7.85788
      b_mtem(3,jA,jE) = -29.9275
      b_mtem(4,jA,jE) = 61.7515
      b_mtem(5,jA,jE) = -63.2308
      b_mtem(6,jA,jE) = 24.9542

! in NH4HSO4
      jE = jnh4hso4
      b_mtem(1,jA,jE) = -1.05891
      b_mtem(2,jA,jE) = 2.84831
      b_mtem(3,jA,jE) = -21.1827
      b_mtem(4,jA,jE) = 57.5175
      b_mtem(5,jA,jE) = -64.8120
      b_mtem(6,jA,jE) = 26.1986

! in (NH4)3H(SO4)2
      jE = jlvcite
      b_mtem(1,jA,jE) = -1.16584
      b_mtem(2,jA,jE) = 8.50075
      b_mtem(3,jA,jE) = -44.3420
      b_mtem(4,jA,jE) = 97.3974
      b_mtem(5,jA,jE) = -98.4549
      b_mtem(6,jA,jE) = 37.6104

! in NaHSO4
      jE = jnahso4
      b_mtem(1,jA,jE) = -1.95805
      b_mtem(2,jA,jE) = 6.62417
      b_mtem(3,jA,jE) = -31.8072
      b_mtem(4,jA,jE) = 77.8603
      b_mtem(5,jA,jE) = -84.6458
      b_mtem(6,jA,jE) = 33.4963

! in Na3H(SO4)2
      jE = jna3hso4
      b_mtem(1,jA,jE) = -0.36045
      b_mtem(2,jA,jE) = 3.55223
      b_mtem(3,jA,jE) = -24.0327
      b_mtem(4,jA,jE) = 54.4879
      b_mtem(5,jA,jE) = -56.6531
      b_mtem(6,jA,jE) = 22.4956


!----------
! NaNO3 in E
      jA = jnano3

! in (NH4)2SO4
      jE = jnh4so4
      b_mtem(1,jA,jE) = -2.5888
      b_mtem(2,jA,jE) = 17.6192
      b_mtem(3,jA,jE) = -63.2183
      b_mtem(4,jA,jE) = 115.3520
      b_mtem(5,jA,jE) = -104.0860
      b_mtem(6,jA,jE) = 36.7390

! in NH4NO3
      jE = jnh4no3
      b_mtem(1,jA,jE) = -2.0669
      b_mtem(2,jA,jE) = 1.4792
      b_mtem(3,jA,jE) = 10.5261
      b_mtem(4,jA,jE) = -27.0987
      b_mtem(5,jA,jE) = 23.0591
      b_mtem(6,jA,jE) = -6.0938

! in NH4Cl (revised on 11/15/2003)
      jE = jnh4cl
      b_mtem(1,jA,jE) = -0.8325
      b_mtem(2,jA,jE) = 3.9933
      b_mtem(3,jA,jE) = -15.3789
      b_mtem(4,jA,jE) = 30.4050
      b_mtem(5,jA,jE) = -29.4204
      b_mtem(6,jA,jE) = 11.0597

! in Na2SO4
      jE = jna2so4
      b_mtem(1,jA,jE) = -1.1233
      b_mtem(2,jA,jE) = 8.3998
      b_mtem(3,jA,jE) = -31.9002
      b_mtem(4,jA,jE) = 60.1450
      b_mtem(5,jA,jE) = -55.5503
      b_mtem(6,jA,jE) = 19.7757

! in NaNO3
      jE = jnano3
      b_mtem(1,jA,jE) = -2.5386
      b_mtem(2,jA,jE) = 13.9039
      b_mtem(3,jA,jE) = -42.8467
      b_mtem(4,jA,jE) = 69.7442
      b_mtem(5,jA,jE) = -57.8988
      b_mtem(6,jA,jE) = 19.4635

! in NaCl
      jE = jnacl
      b_mtem(1,jA,jE) = -0.4351
      b_mtem(2,jA,jE) = 2.8311
      b_mtem(3,jA,jE) = -11.4485
      b_mtem(4,jA,jE) = 22.7201
      b_mtem(5,jA,jE) = -22.4228
      b_mtem(6,jA,jE) = 8.5792

! in Ca(NO3)2
      jE = jcano3
      b_mtem(1,jA,jE) = -0.72060
      b_mtem(2,jA,jE) = 5.64915
      b_mtem(3,jA,jE) = -23.5020
      b_mtem(4,jA,jE) = 46.0078
      b_mtem(5,jA,jE) = -43.8075
      b_mtem(6,jA,jE) = 16.1652

! in CaCl2
      jE = jcacl2
      b_mtem(1,jA,jE) = 0.003928
      b_mtem(2,jA,jE) = 3.54724
      b_mtem(3,jA,jE) = -18.6057
      b_mtem(4,jA,jE) = 38.1445
      b_mtem(5,jA,jE) = -36.7745
      b_mtem(6,jA,jE) = 13.4529

! in HNO3
      jE = jhno3
      b_mtem(1,jA,jE) = -1.1712
      b_mtem(2,jA,jE) = 7.20907
      b_mtem(3,jA,jE) = -22.9215
      b_mtem(4,jA,jE) = 38.1257
      b_mtem(5,jA,jE) = -32.0759
      b_mtem(6,jA,jE) = 10.6443

! in HCl
      jE = jhcl
      b_mtem(1,jA,jE) = 0.738022
      b_mtem(2,jA,jE) = -1.14313
      b_mtem(3,jA,jE) = 0.32251
      b_mtem(4,jA,jE) = 0.838679
      b_mtem(5,jA,jE) = -1.81747
      b_mtem(6,jA,jE) = 0.873986


!----------
! NaCl in E
      jA = jnacl

! in (NH4)2SO4
      jE = jnh4so4
      b_mtem(1,jA,jE) = -1.9525
      b_mtem(2,jA,jE) = 16.6433
      b_mtem(3,jA,jE) = -61.7090
      b_mtem(4,jA,jE) = 112.9910
      b_mtem(5,jA,jE) = -101.9370
      b_mtem(6,jA,jE) = 35.7760

! in NH4NO3
      jE = jnh4no3
      b_mtem(1,jA,jE) = -1.7525
      b_mtem(2,jA,jE) = 3.0713
      b_mtem(3,jA,jE) = 4.8063
      b_mtem(4,jA,jE) = -17.5334
      b_mtem(5,jA,jE) = 14.2872
      b_mtem(6,jA,jE) = -3.0690

! in NH4Cl (revised on 11/15/2003)
      jE = jnh4cl
      b_mtem(1,jA,jE) = -0.4021
      b_mtem(2,jA,jE) = 5.2399
      b_mtem(3,jA,jE) = -19.4278
      b_mtem(4,jA,jE) = 33.0027
      b_mtem(5,jA,jE) = -28.1020
      b_mtem(6,jA,jE) = 9.5159

! in Na2SO4
      jE = jna2so4
      b_mtem(1,jA,jE) = 0.6692
      b_mtem(2,jA,jE) = 4.1207
      b_mtem(3,jA,jE) = -27.3314
      b_mtem(4,jA,jE) = 59.3112
      b_mtem(5,jA,jE) = -58.7998
      b_mtem(6,jA,jE) = 21.7674

! in NaNO3
      jE = jnano3
      b_mtem(1,jA,jE) = -1.17444
      b_mtem(2,jA,jE) = 10.9927
      b_mtem(3,jA,jE) = -38.9013
      b_mtem(4,jA,jE) = 66.8521
      b_mtem(5,jA,jE) = -57.6564
      b_mtem(6,jA,jE) = 19.7296

! in NaCl
      jE = jnacl
      b_mtem(1,jA,jE) = 1.17679
      b_mtem(2,jA,jE) = -2.5061
      b_mtem(3,jA,jE) = 0.8508
      b_mtem(4,jA,jE) = 4.4802
      b_mtem(5,jA,jE) = -8.4945
      b_mtem(6,jA,jE) = 4.3182

! in Ca(NO3)2
      jE = jcano3
      b_mtem(1,jA,jE) = 1.01450
      b_mtem(2,jA,jE) = 2.10260
      b_mtem(3,jA,jE) = -20.9036
      b_mtem(4,jA,jE) = 49.1481
      b_mtem(5,jA,jE) = -51.4867
      b_mtem(6,jA,jE) = 19.9301

! in CaCl2 (PSC92: revised on 11/27/2003)
      jE = jcacl2
      b_mtem(1,jA,jE) = 1.55463
      b_mtem(2,jA,jE) = -3.20122
      b_mtem(3,jA,jE) = -0.957075
      b_mtem(4,jA,jE) = 12.103
      b_mtem(5,jA,jE) = -17.221
      b_mtem(6,jA,jE) = 7.50264

! in HNO3
      jE = jhno3
      b_mtem(1,jA,jE) = 2.46187
      b_mtem(2,jA,jE) = -12.6845
      b_mtem(3,jA,jE) = 34.2383
      b_mtem(4,jA,jE) = -51.9992
      b_mtem(5,jA,jE) = 39.4934
      b_mtem(6,jA,jE) = -11.7247

! in HCl
      jE = jhcl
      b_mtem(1,jA,jE) = 1.74915
      b_mtem(2,jA,jE) = -4.65768
      b_mtem(3,jA,jE) = 8.80287
      b_mtem(4,jA,jE) = -12.2503
      b_mtem(5,jA,jE) = 8.668751
      b_mtem(6,jA,jE) = -2.50158


!----------
! Ca(NO3)2 in E
      jA = jcano3

! in NH4NO3
      jE = jnh4no3
      b_mtem(1,jA,jE) = -1.86260
      b_mtem(2,jA,jE) = 11.6178
      b_mtem(3,jA,jE) = -30.9069
      b_mtem(4,jA,jE) = 41.7578
      b_mtem(5,jA,jE) = -33.7338
      b_mtem(6,jA,jE) = 12.7541

! in NH4Cl (revised on 11/15/2003)
      jE = jnh4cl
      b_mtem(1,jA,jE) = -1.1798
      b_mtem(2,jA,jE) = 25.9608
      b_mtem(3,jA,jE) = -98.9373
      b_mtem(4,jA,jE) = 160.2300
      b_mtem(5,jA,jE) = -125.9540
      b_mtem(6,jA,jE) = 39.5130

! in NaNO3
      jE = jnano3
      b_mtem(1,jA,jE) = -1.44384
      b_mtem(2,jA,jE) = 13.6044
      b_mtem(3,jA,jE) = -54.4300
      b_mtem(4,jA,jE) = 100.582
      b_mtem(5,jA,jE) = -91.2364
      b_mtem(6,jA,jE) = 32.5970

! in NaCl
      jE = jnacl
      b_mtem(1,jA,jE) = -0.099114
      b_mtem(2,jA,jE) = 2.84091
      b_mtem(3,jA,jE) = -16.9229
      b_mtem(4,jA,jE) = 37.4839
      b_mtem(5,jA,jE) = -39.5132
      b_mtem(6,jA,jE) = 15.8564

! in Ca(NO3)2
      jE = jcano3
      b_mtem(1,jA,jE) = 0.055116
      b_mtem(2,jA,jE) = 4.58610
      b_mtem(3,jA,jE) = -27.6629
      b_mtem(4,jA,jE) = 60.8288
      b_mtem(5,jA,jE) = -61.4988
      b_mtem(6,jA,jE) = 23.3136

! in CaCl2 (PSC92: revised on 11/27/2003)
      jE = jcacl2
      b_mtem(1,jA,jE) = 1.57155
      b_mtem(2,jA,jE) = -3.18486
      b_mtem(3,jA,jE) = -3.35758
      b_mtem(4,jA,jE) = 18.7501
      b_mtem(5,jA,jE) = -24.5604
      b_mtem(6,jA,jE) = 10.3798

! in HNO3
      jE = jhno3
      b_mtem(1,jA,jE) = 1.04446
      b_mtem(2,jA,jE) = -3.19066
      b_mtem(3,jA,jE) = 2.44714
      b_mtem(4,jA,jE) = 2.07218
      b_mtem(5,jA,jE) = -6.43949
      b_mtem(6,jA,jE) = 3.66471

! in HCl
      jE = jhcl
      b_mtem(1,jA,jE) = 1.05723
      b_mtem(2,jA,jE) = -1.46826
      b_mtem(3,jA,jE) = -1.0713
      b_mtem(4,jA,jE) = 4.64439
      b_mtem(5,jA,jE) = -6.32402
      b_mtem(6,jA,jE) = 2.78202


!----------
! CaCl2 in E
      jA = jcacl2

! in NH4NO3 (PSC92: revised on 12/22/2003)
      jE = jnh4no3
      b_mtem(1,jA,jE) = -1.43626
      b_mtem(2,jA,jE) = 13.6598
      b_mtem(3,jA,jE) = -38.2068
      b_mtem(4,jA,jE) = 53.9057
      b_mtem(5,jA,jE) = -44.9018
      b_mtem(6,jA,jE) = 16.6120

! in NH4Cl (PSC92: revised on 11/27/2003)
      jE = jnh4cl
      b_mtem(1,jA,jE) = -0.603965
      b_mtem(2,jA,jE) = 27.6027
      b_mtem(3,jA,jE) = -104.258
      b_mtem(4,jA,jE) = 163.553
      b_mtem(5,jA,jE) = -124.076
      b_mtem(6,jA,jE) = 37.4153

! in NaNO3 (PSC92: revised on 12/22/2003)
      jE = jnano3
      b_mtem(1,jA,jE) = 0.44648
      b_mtem(2,jA,jE) = 8.8850
      b_mtem(3,jA,jE) = -45.5232
      b_mtem(4,jA,jE) = 89.3263
      b_mtem(5,jA,jE) = -83.8604
      b_mtem(6,jA,jE) = 30.4069

! in NaCl (PSC92: revised on 11/27/2003)
      jE = jnacl
      b_mtem(1,jA,jE) = 1.61927
      b_mtem(2,jA,jE) = 0.247547
      b_mtem(3,jA,jE) = -18.1252
      b_mtem(4,jA,jE) = 45.2479
      b_mtem(5,jA,jE) = -48.6072
      b_mtem(6,jA,jE) = 19.2784

! in Ca(NO3)2 (PSC92: revised on 11/27/2003)
      jE = jcano3
      b_mtem(1,jA,jE) = 2.36667
      b_mtem(2,jA,jE) = -0.123309
      b_mtem(3,jA,jE) = -24.2723
      b_mtem(4,jA,jE) = 65.1486
      b_mtem(5,jA,jE) = -71.8504
      b_mtem(6,jA,jE) = 28.3696

! in CaCl2 (PSC92: revised on 11/27/2003)
      jE = jcacl2
      b_mtem(1,jA,jE) = 3.64023
      b_mtem(2,jA,jE) = -12.1926
      b_mtem(3,jA,jE) = 20.2028
      b_mtem(4,jA,jE) = -16.0056
      b_mtem(5,jA,jE) = 1.52355
      b_mtem(6,jA,jE) = 2.44709

! in HNO3
      jE = jhno3
      b_mtem(1,jA,jE) = 5.88794
      b_mtem(2,jA,jE) = -29.7083
      b_mtem(3,jA,jE) = 78.6309
      b_mtem(4,jA,jE) = -118.037
      b_mtem(5,jA,jE) = 88.932
      b_mtem(6,jA,jE) = -26.1407

! in HCl
      jE = jhcl
      b_mtem(1,jA,jE) = 2.40628
      b_mtem(2,jA,jE) = -6.16566
      b_mtem(3,jA,jE) = 10.2851
      b_mtem(4,jA,jE) = -12.9035
      b_mtem(5,jA,jE) = 7.7441
      b_mtem(6,jA,jE) = -1.74821


!----------
! HNO3 in E
      jA = jhno3

! in (NH4)2SO4
      jE = jnh4so4
      b_mtem(1,jA,jE) = -3.57598
      b_mtem(2,jA,jE) = 21.5469
      b_mtem(3,jA,jE) = -77.4111
      b_mtem(4,jA,jE) = 144.136
      b_mtem(5,jA,jE) = -132.849
      b_mtem(6,jA,jE) = 47.9412

! in NH4NO3
      jE = jnh4no3
      b_mtem(1,jA,jE) = -2.00209
      b_mtem(2,jA,jE) = -3.48399
      b_mtem(3,jA,jE) = 34.9906
      b_mtem(4,jA,jE) = -68.6653
      b_mtem(5,jA,jE) = 54.0992
      b_mtem(6,jA,jE) = -15.1343

! in NH4Cl revised on 12/22/2003
      jE = jnh4cl
      b_mtem(1,jA,jE) = -0.63790
      b_mtem(2,jA,jE) = -1.67730
      b_mtem(3,jA,jE) = 10.1727
      b_mtem(4,jA,jE) = -14.9097
      b_mtem(5,jA,jE) = 7.67410
      b_mtem(6,jA,jE) = -0.79586

! in NaCl
      jE = jnacl
      b_mtem(1,jA,jE) = 1.3446
      b_mtem(2,jA,jE) = -2.5578
      b_mtem(3,jA,jE) = 1.3464
      b_mtem(4,jA,jE) = 2.90537
      b_mtem(5,jA,jE) = -6.53014
      b_mtem(6,jA,jE) = 3.31339

! in NaNO3
      jE = jnano3
      b_mtem(1,jA,jE) = -0.546636
      b_mtem(2,jA,jE) = 10.3127
      b_mtem(3,jA,jE) = -39.9603
      b_mtem(4,jA,jE) = 71.4609
      b_mtem(5,jA,jE) = -63.4958
      b_mtem(6,jA,jE) = 22.0679

! in Na2SO4
      jE = jna2so4
      b_mtem(1,jA,jE) = 1.35059
      b_mtem(2,jA,jE) = 4.34557
      b_mtem(3,jA,jE) = -35.8425
      b_mtem(4,jA,jE) = 80.9868
      b_mtem(5,jA,jE) = -81.6544
      b_mtem(6,jA,jE) = 30.4841

! in Ca(NO3)2
      jE = jcano3
      b_mtem(1,jA,jE) = 0.869414
      b_mtem(2,jA,jE) = 2.98486
      b_mtem(3,jA,jE) = -22.255
      b_mtem(4,jA,jE) = 50.1863
      b_mtem(5,jA,jE) = -51.214
      b_mtem(6,jA,jE) = 19.2235

! in CaCl2 (KM) revised on 12/22/2003
      jE = jcacl2
      b_mtem(1,jA,jE) = 1.42800
      b_mtem(2,jA,jE) = -1.78959
      b_mtem(3,jA,jE) = -2.49075
      b_mtem(4,jA,jE) = 10.1877
      b_mtem(5,jA,jE) = -12.1948
      b_mtem(6,jA,jE) = 4.64475

! in HNO3 (added on 12/06/2004)
      jE = jhno3
      b_mtem(1,jA,jE) = 0.22035
      b_mtem(2,jA,jE) = 2.94973
      b_mtem(3,jA,jE) = -12.1469
      b_mtem(4,jA,jE) = 20.4905
      b_mtem(5,jA,jE) = -17.3966
      b_mtem(6,jA,jE) = 5.70779

! in HCl (added on 12/06/2004)
      jE = jhcl
      b_mtem(1,jA,jE) = 1.55503
      b_mtem(2,jA,jE) = -3.61226
      b_mtem(3,jA,jE) = 6.28265
      b_mtem(4,jA,jE) = -8.69575
      b_mtem(5,jA,jE) = 6.09372
      b_mtem(6,jA,jE) = -1.80898

! in H2SO4
      jE = jh2so4
      b_mtem(1,jA,jE) = 1.10783
      b_mtem(2,jA,jE) = -1.3363
      b_mtem(3,jA,jE) = -1.83525
      b_mtem(4,jA,jE) = 7.47373
      b_mtem(5,jA,jE) = -9.72954
      b_mtem(6,jA,jE) = 4.12248

! in NH4HSO4
      jE = jnh4hso4
      b_mtem(1,jA,jE) = -0.851026
      b_mtem(2,jA,jE) = 12.2515
      b_mtem(3,jA,jE) = -49.788
      b_mtem(4,jA,jE) = 91.6215
      b_mtem(5,jA,jE) = -81.4877
      b_mtem(6,jA,jE) = 28.0002

! in (NH4)3H(SO4)2
      jE = jlvcite
      b_mtem(1,jA,jE) = -3.09464
      b_mtem(2,jA,jE) = 14.9303
      b_mtem(3,jA,jE) = -43.0454
      b_mtem(4,jA,jE) = 72.6695
      b_mtem(5,jA,jE) = -65.2140
      b_mtem(6,jA,jE) = 23.4814

! in NaHSO4
      jE = jnahso4
      b_mtem(1,jA,jE) = 1.22973
      b_mtem(2,jA,jE) = 2.82702
      b_mtem(3,jA,jE) = -17.5869
      b_mtem(4,jA,jE) = 28.9564
      b_mtem(5,jA,jE) = -23.5814
      b_mtem(6,jA,jE) = 7.91153

! in Na3H(SO4)2
      jE = jna3hso4
      b_mtem(1,jA,jE) = 1.64773
      b_mtem(2,jA,jE) = 0.94188
      b_mtem(3,jA,jE) = -19.1242
      b_mtem(4,jA,jE) = 46.9887
      b_mtem(5,jA,jE) = -50.9494
      b_mtem(6,jA,jE) = 20.2169


!----------
! HCl in E
      jA = jhcl

! in (NH4)2SO4
      jE = jnh4so4
      b_mtem(1,jA,jE) = -2.93783
      b_mtem(2,jA,jE) = 20.5546
      b_mtem(3,jA,jE) = -75.8548
      b_mtem(4,jA,jE) = 141.729
      b_mtem(5,jA,jE) = -130.697
      b_mtem(6,jA,jE) = 46.9905

! in NH4NO3
      jE = jnh4no3
      b_mtem(1,jA,jE) = -1.69063
      b_mtem(2,jA,jE) = -1.85303
      b_mtem(3,jA,jE) = 29.0927
      b_mtem(4,jA,jE) = -58.7401
      b_mtem(5,jA,jE) = 44.999
      b_mtem(6,jA,jE) = -11.9988

! in NH4Cl (revised on 11/15/2003)
      jE = jnh4cl
      b_mtem(1,jA,jE) = -0.2073
      b_mtem(2,jA,jE) = -0.4322
      b_mtem(3,jA,jE) = 6.1271
      b_mtem(4,jA,jE) = -12.3146
      b_mtem(5,jA,jE) = 8.9919
      b_mtem(6,jA,jE) = -2.3388

! in NaCl
      jE = jnacl
      b_mtem(1,jA,jE) = 2.95913
      b_mtem(2,jA,jE) = -7.92254
      b_mtem(3,jA,jE) = 13.736
      b_mtem(4,jA,jE) = -15.433
      b_mtem(5,jA,jE) = 7.40386
      b_mtem(6,jA,jE) = -0.918641

! in NaNO3
      jE = jnano3
      b_mtem(1,jA,jE) = 0.893272
      b_mtem(2,jA,jE) = 6.53768
      b_mtem(3,jA,jE) = -32.3458
      b_mtem(4,jA,jE) = 61.2834
      b_mtem(5,jA,jE) = -56.4446
      b_mtem(6,jA,jE) = 19.9202

! in Na2SO4
      jE = jna2so4
      b_mtem(1,jA,jE) = 3.14484
      b_mtem(2,jA,jE) = 0.077019
      b_mtem(3,jA,jE) = -31.4199
      b_mtem(4,jA,jE) = 80.5865
      b_mtem(5,jA,jE) = -85.392
      b_mtem(6,jA,jE) = 32.6644

! in Ca(NO3)2
      jE = jcano3
      b_mtem(1,jA,jE) = 2.60432
      b_mtem(2,jA,jE) = -0.55909
      b_mtem(3,jA,jE) = -19.6671
      b_mtem(4,jA,jE) = 53.3446
      b_mtem(5,jA,jE) = -58.9076
      b_mtem(6,jA,jE) = 22.9927

! in CaCl2 (KM) revised on 3/13/2003 and again on 11/27/2003
      jE = jcacl2
      b_mtem(1,jA,jE) = 2.98036
      b_mtem(2,jA,jE) = -8.55365
      b_mtem(3,jA,jE) = 15.2108
      b_mtem(4,jA,jE) = -15.9359
      b_mtem(5,jA,jE) = 7.41772
      b_mtem(6,jA,jE) = -1.32143

! in HNO3 (added on 12/06/2004)
      jE = jhno3
      b_mtem(1,jA,jE) = 3.8533
      b_mtem(2,jA,jE) = -16.9427
      b_mtem(3,jA,jE) = 45.0056
      b_mtem(4,jA,jE) = -69.6145
      b_mtem(5,jA,jE) = 54.1491
      b_mtem(6,jA,jE) = -16.6513

! in HCl (added on 12/06/2004)
      jE = jhcl
      b_mtem(1,jA,jE) = 2.56665
      b_mtem(2,jA,jE) = -7.13585
      b_mtem(3,jA,jE) = 14.8103
      b_mtem(4,jA,jE) = -21.8881
      b_mtem(5,jA,jE) = 16.6808
      b_mtem(6,jA,jE) = -5.22091

! in H2SO4
      jE = jh2so4
      b_mtem(1,jA,jE) = 2.50179
      b_mtem(2,jA,jE) = -6.69364
      b_mtem(3,jA,jE) = 11.6551
      b_mtem(4,jA,jE) = -13.6897
      b_mtem(5,jA,jE) = 7.36796
      b_mtem(6,jA,jE) = -1.33245

! in NH4HSO4
      jE = jnh4hso4
      b_mtem(1,jA,jE) = 0.149955
      b_mtem(2,jA,jE) = 11.8213
      b_mtem(3,jA,jE) = -53.9164
      b_mtem(4,jA,jE) = 101.574
      b_mtem(5,jA,jE) = -91.4123
      b_mtem(6,jA,jE) = 31.5487

! in (NH4)3H(SO4)2
      jE = jlvcite
      b_mtem(1,jA,jE) = -2.36927
      b_mtem(2,jA,jE) = 14.8359
      b_mtem(3,jA,jE) = -44.3443
      b_mtem(4,jA,jE) = 73.6229
      b_mtem(5,jA,jE) = -65.3366
      b_mtem(6,jA,jE) = 23.3250

! in NaHSO4
      jE = jnahso4
      b_mtem(1,jA,jE) = 2.72993
      b_mtem(2,jA,jE) = -0.23406
      b_mtem(3,jA,jE) = -10.4103
      b_mtem(4,jA,jE) = 13.1586
      b_mtem(5,jA,jE) = -7.79925
      b_mtem(6,jA,jE) = 2.30843

! in Na3H(SO4)2
      jE = jna3hso4
      b_mtem(1,jA,jE) = 3.51258
      b_mtem(2,jA,jE) = -3.95107
      b_mtem(3,jA,jE) = -11.0175
      b_mtem(4,jA,jE) = 38.8617
      b_mtem(5,jA,jE) = -48.1575
      b_mtem(6,jA,jE) = 20.4717


!----------
! 2H.SO4 in E
      jA = jh2so4

! in H2SO4
      jE = jh2so4
      b_mtem(1,jA,jE) = 0.76734
      b_mtem(2,jA,jE) = -1.12263
      b_mtem(3,jA,jE) = -9.08728
      b_mtem(4,jA,jE) = 30.3836
      b_mtem(5,jA,jE) = -38.4133
      b_mtem(6,jA,jE) = 17.0106

! in NH4HSO4
      jE = jnh4hso4
      b_mtem(1,jA,jE) = -2.03879
      b_mtem(2,jA,jE) = 15.7033
      b_mtem(3,jA,jE) = -58.7363
      b_mtem(4,jA,jE) = 109.242
      b_mtem(5,jA,jE) = -102.237
      b_mtem(6,jA,jE) = 37.5350

! in (NH4)3H(SO4)2
      jE = jlvcite
      b_mtem(1,jA,jE) = -3.10228
      b_mtem(2,jA,jE) = 16.6920
      b_mtem(3,jA,jE) = -59.1522
      b_mtem(4,jA,jE) = 113.487
      b_mtem(5,jA,jE) = -110.890
      b_mtem(6,jA,jE) = 42.4578

! in (NH4)2SO4
      jE = jnh4so4
      b_mtem(1,jA,jE) = -3.43885
      b_mtem(2,jA,jE) = 21.0372
      b_mtem(3,jA,jE) = -84.7026
      b_mtem(4,jA,jE) = 165.324
      b_mtem(5,jA,jE) = -156.101
      b_mtem(6,jA,jE) = 57.3101

! in NaHSO4
      jE = jnahso4
      b_mtem(1,jA,jE) = 0.33164
      b_mtem(2,jA,jE) = 6.55864
      b_mtem(3,jA,jE) = -33.5876
      b_mtem(4,jA,jE) = 65.1798
      b_mtem(5,jA,jE) = -63.2046
      b_mtem(6,jA,jE) = 24.1783

! in Na3H(SO4)2
      jE = jna3hso4
      b_mtem(1,jA,jE) = 3.06830
      b_mtem(2,jA,jE) = -3.18408
      b_mtem(3,jA,jE) = -19.6332
      b_mtem(4,jA,jE) = 61.3657
      b_mtem(5,jA,jE) = -73.4438
      b_mtem(6,jA,jE) = 31.2334

! in Na2SO4
      jE = jna2so4
      b_mtem(1,jA,jE) = 2.58649
      b_mtem(2,jA,jE) = 0.87921
      b_mtem(3,jA,jE) = -39.3023
      b_mtem(4,jA,jE) = 101.603
      b_mtem(5,jA,jE) = -109.469
      b_mtem(6,jA,jE) = 43.0188

! in HNO3
      jE = jhno3
      b_mtem(1,jA,jE) = 1.54587
      b_mtem(2,jA,jE) = -7.50976
      b_mtem(3,jA,jE) = 12.8237
      b_mtem(4,jA,jE) = -10.1452
      b_mtem(5,jA,jE) = -0.541956
      b_mtem(6,jA,jE) = 3.34536

! in HCl
      jE = jhcl
      b_mtem(1,jA,jE) = 0.829757
      b_mtem(2,jA,jE) = -4.11316
      b_mtem(3,jA,jE) = 3.67111
      b_mtem(4,jA,jE) = 3.6833
      b_mtem(5,jA,jE) = -11.2711
      b_mtem(6,jA,jE) = 6.71421


!----------
! H.HSO4 in E
      jA = jhhso4

! in H2SO4
      jE = jh2so4
      b_mtem(1,jA,jE) = 2.63953
      b_mtem(2,jA,jE) = -6.01532
      b_mtem(3,jA,jE) = 10.0204
      b_mtem(4,jA,jE) = -12.4840
      b_mtem(5,jA,jE) = 7.78853
      b_mtem(6,jA,jE) = -2.12638

! in NH4HSO4
      jE = jnh4hso4
      b_mtem(1,jA,jE) = -0.77412
      b_mtem(2,jA,jE) = 14.1656
      b_mtem(3,jA,jE) = -53.4087
      b_mtem(4,jA,jE) = 93.2013
      b_mtem(5,jA,jE) = -80.5723
      b_mtem(6,jA,jE) = 27.1577

! in (NH4)3H(SO4)2
      jE = jlvcite
      b_mtem(1,jA,jE) = -2.98882
      b_mtem(2,jA,jE) = 14.4436
      b_mtem(3,jA,jE) = -40.1774
      b_mtem(4,jA,jE) = 67.5937
      b_mtem(5,jA,jE) = -61.5040
      b_mtem(6,jA,jE) = 22.3695

! in (NH4)2SO4
      jE = jnh4so4
      b_mtem(1,jA,jE) = -1.15502
      b_mtem(2,jA,jE) = 8.12309
      b_mtem(3,jA,jE) = -38.4726
      b_mtem(4,jA,jE) = 80.8861
      b_mtem(5,jA,jE) = -80.1644
      b_mtem(6,jA,jE) = 30.4717

! in NaHSO4
      jE = jnahso4
      b_mtem(1,jA,jE) = 1.99641
      b_mtem(2,jA,jE) = -2.96061
      b_mtem(3,jA,jE) = 5.54778
      b_mtem(4,jA,jE) = -14.5488
      b_mtem(5,jA,jE) = 14.8492
      b_mtem(6,jA,jE) = -5.1389

! in Na3H(SO4)2
      jE = jna3hso4
      b_mtem(1,jA,jE) = 2.23816
      b_mtem(2,jA,jE) = -3.20847
      b_mtem(3,jA,jE) = -4.82853
      b_mtem(4,jA,jE) = 20.9192
      b_mtem(5,jA,jE) = -27.2819
      b_mtem(6,jA,jE) = 11.8655

! in Na2SO4
      jE = jna2so4
      b_mtem(1,jA,jE) = 2.56907
      b_mtem(2,jA,jE) = 1.13444
      b_mtem(3,jA,jE) = -34.6853
      b_mtem(4,jA,jE) = 87.9775
      b_mtem(5,jA,jE) = -93.2330
      b_mtem(6,jA,jE) = 35.9260

! in HNO3
      jE = jhno3
      b_mtem(1,jA,jE) = 2.00024
      b_mtem(2,jA,jE) = -4.80868
      b_mtem(3,jA,jE) = 8.29222
      b_mtem(4,jA,jE) = -11.0849
      b_mtem(5,jA,jE) = 7.51262
      b_mtem(6,jA,jE) = -2.07654

! in HCl
      jE = jhcl
      b_mtem(1,jA,jE) = 2.8009
      b_mtem(2,jA,jE) = -6.98416
      b_mtem(3,jA,jE) = 14.3146
      b_mtem(4,jA,jE) = -22.0068
      b_mtem(5,jA,jE) = 17.5557
      b_mtem(6,jA,jE) = -5.84917


!----------
! NH4HSO4 in E
      jA = jnh4hso4

! in H2SO4
      jE = jh2so4
      b_mtem(1,jA,jE) = 0.169160
      b_mtem(2,jA,jE) = 2.15094
      b_mtem(3,jA,jE) = -9.62904
      b_mtem(4,jA,jE) = 18.2631
      b_mtem(5,jA,jE) = -17.3333
      b_mtem(6,jA,jE) = 6.19835

! in NH4HSO4
      jE = jnh4hso4
      b_mtem(1,jA,jE) = -2.34457
      b_mtem(2,jA,jE) = 12.8035
      b_mtem(3,jA,jE) = -35.2513
      b_mtem(4,jA,jE) = 53.6153
      b_mtem(5,jA,jE) = -42.7655
      b_mtem(6,jA,jE) = 13.7129

! in (NH4)3H(SO4)2
      jE = jlvcite
      b_mtem(1,jA,jE) = -2.56109
      b_mtem(2,jA,jE) = 11.1414
      b_mtem(3,jA,jE) = -30.2361
      b_mtem(4,jA,jE) = 50.0320
      b_mtem(5,jA,jE) = -44.1586
      b_mtem(6,jA,jE) = 15.5393

! in (NH4)2SO4
      jE = jnh4so4
      b_mtem(1,jA,jE) = -0.97315
      b_mtem(2,jA,jE) = 7.06295
      b_mtem(3,jA,jE) = -29.3032
      b_mtem(4,jA,jE) = 57.6101
      b_mtem(5,jA,jE) = -54.9020
      b_mtem(6,jA,jE) = 20.2222

! in NaHSO4
      jE = jnahso4
      b_mtem(1,jA,jE) = -0.44450
      b_mtem(2,jA,jE) = 3.33451
      b_mtem(3,jA,jE) = -15.2791
      b_mtem(4,jA,jE) = 30.1413
      b_mtem(5,jA,jE) = -26.7710
      b_mtem(6,jA,jE) = 8.78462

! in Na3H(SO4)2
      jE = jna3hso4
      b_mtem(1,jA,jE) = -0.99780
      b_mtem(2,jA,jE) = 4.69200
      b_mtem(3,jA,jE) = -16.1219
      b_mtem(4,jA,jE) = 29.3100
      b_mtem(5,jA,jE) = -26.3383
      b_mtem(6,jA,jE) = 9.20695

! in Na2SO4
      jE = jna2so4
      b_mtem(1,jA,jE) = -0.52694
      b_mtem(2,jA,jE) = 7.02684
      b_mtem(3,jA,jE) = -33.7508
      b_mtem(4,jA,jE) = 70.0565
      b_mtem(5,jA,jE) = -68.3226
      b_mtem(6,jA,jE) = 25.2692

! in HNO3
      jE = jhno3
      b_mtem(1,jA,jE) = 0.572926
      b_mtem(2,jA,jE) = -2.04791
      b_mtem(3,jA,jE) = 2.1134
      b_mtem(4,jA,jE) = 0.246654
      b_mtem(5,jA,jE) = -3.06019
      b_mtem(6,jA,jE) = 1.98126

! in HCl
      jE = jhcl
      b_mtem(1,jA,jE) = 0.56514
      b_mtem(2,jA,jE) = 0.22287
      b_mtem(3,jA,jE) = -2.76973
      b_mtem(4,jA,jE) = 4.54444
      b_mtem(5,jA,jE) = -3.86549
      b_mtem(6,jA,jE) = 1.13441


!----------
! (NH4)3H(SO4)2 in E
      jA = jlvcite

! in H2SO4
      jE = jh2so4
      b_mtem(1,jA,jE) = -1.44811
      b_mtem(2,jA,jE) = 6.71815
      b_mtem(3,jA,jE) = -25.0141
      b_mtem(4,jA,jE) = 50.1109
      b_mtem(5,jA,jE) = -50.0561
      b_mtem(6,jA,jE) = 19.3370

! in NH4HSO4
      jE = jnh4hso4
      b_mtem(1,jA,jE) = -3.41707
      b_mtem(2,jA,jE) = 13.4496
      b_mtem(3,jA,jE) = -34.8018
      b_mtem(4,jA,jE) = 55.2987
      b_mtem(5,jA,jE) = -48.1839
      b_mtem(6,jA,jE) = 17.2444

! in (NH4)3H(SO4)2
      jE = jlvcite
      b_mtem(1,jA,jE) = -2.54479
      b_mtem(2,jA,jE) = 11.8501
      b_mtem(3,jA,jE) = -39.7286
      b_mtem(4,jA,jE) = 74.2479
      b_mtem(5,jA,jE) = -70.4934
      b_mtem(6,jA,jE) = 26.2836

! in (NH4)2SO4
      jE = jnh4so4
      b_mtem(1,jA,jE) = -2.30561
      b_mtem(2,jA,jE) = 14.5806
      b_mtem(3,jA,jE) = -55.1238
      b_mtem(4,jA,jE) = 103.451
      b_mtem(5,jA,jE) = -95.2571
      b_mtem(6,jA,jE) = 34.2218

! in NaHSO4
      jE = jnahso4
      b_mtem(1,jA,jE) = -2.20809
      b_mtem(2,jA,jE) = 13.6391
      b_mtem(3,jA,jE) = -57.8246
      b_mtem(4,jA,jE) = 117.907
      b_mtem(5,jA,jE) = -112.154
      b_mtem(6,jA,jE) = 40.3058

! in Na3H(SO4)2
      jE = jna3hso4
      b_mtem(1,jA,jE) = -1.15099
      b_mtem(2,jA,jE) = 6.32269
      b_mtem(3,jA,jE) = -27.3860
      b_mtem(4,jA,jE) = 55.4592
      b_mtem(5,jA,jE) = -54.0100
      b_mtem(6,jA,jE) = 20.3469

! in Na2SO4
      jE = jna2so4
      b_mtem(1,jA,jE) = -1.15678
      b_mtem(2,jA,jE) = 8.28718
      b_mtem(3,jA,jE) = -37.3231
      b_mtem(4,jA,jE) = 76.6124
      b_mtem(5,jA,jE) = -74.9307
      b_mtem(6,jA,jE) = 28.0559

! in HNO3
      jE = jhno3
      b_mtem(1,jA,jE) = 0.01502
      b_mtem(2,jA,jE) = -3.1197
      b_mtem(3,jA,jE) = 3.61104
      b_mtem(4,jA,jE) = 3.05196
      b_mtem(5,jA,jE) = -9.98957
      b_mtem(6,jA,jE) = 6.04155

! in HCl
      jE = jhcl
      b_mtem(1,jA,jE) = -1.06477
      b_mtem(2,jA,jE) = 3.38801
      b_mtem(3,jA,jE) = -12.5784
      b_mtem(4,jA,jE) = 25.2823
      b_mtem(5,jA,jE) = -25.4611
      b_mtem(6,jA,jE) = 10.0754


!----------
! NaHSO4 in E
      jA = jnahso4

! in H2SO4
      jE = jh2so4
      b_mtem(1,jA,jE) = 0.68259
      b_mtem(2,jA,jE) = 0.71468
      b_mtem(3,jA,jE) = -5.59003
      b_mtem(4,jA,jE) = 11.0089
      b_mtem(5,jA,jE) = -10.7983
      b_mtem(6,jA,jE) = 3.82335

! in NH4HSO4
      jE = jnh4hso4
      b_mtem(1,jA,jE) = -0.03956
      b_mtem(2,jA,jE) = 4.52828
      b_mtem(3,jA,jE) = -25.2557
      b_mtem(4,jA,jE) = 54.4225
      b_mtem(5,jA,jE) = -52.5105
      b_mtem(6,jA,jE) = 18.6562

! in (NH4)3H(SO4)2
      jE = jlvcite
      b_mtem(1,jA,jE) = -1.53503
      b_mtem(2,jA,jE) = 8.27608
      b_mtem(3,jA,jE) = -28.9539
      b_mtem(4,jA,jE) = 55.2876
      b_mtem(5,jA,jE) = -51.9563
      b_mtem(6,jA,jE) = 18.6576

! in (NH4)2SO4
      jE = jnh4so4
      b_mtem(1,jA,jE) = -0.38793
      b_mtem(2,jA,jE) = 7.14680
      b_mtem(3,jA,jE) = -38.7201
      b_mtem(4,jA,jE) = 84.3965
      b_mtem(5,jA,jE) = -84.7453
      b_mtem(6,jA,jE) = 32.1283

! in NaHSO4
      jE = jnahso4
      b_mtem(1,jA,jE) = -0.41982
      b_mtem(2,jA,jE) = 4.26491
      b_mtem(3,jA,jE) = -20.2351
      b_mtem(4,jA,jE) = 42.6764
      b_mtem(5,jA,jE) = -40.7503
      b_mtem(6,jA,jE) = 14.2868

! in Na3H(SO4)2
      jE = jna3hso4
      b_mtem(1,jA,jE) = -0.32912
      b_mtem(2,jA,jE) = 1.80808
      b_mtem(3,jA,jE) = -8.01286
      b_mtem(4,jA,jE) = 15.5791
      b_mtem(5,jA,jE) = -14.5494
      b_mtem(6,jA,jE) = 5.27052

! in Na2SO4
      jE = jna2so4
      b_mtem(1,jA,jE) = 0.10271
      b_mtem(2,jA,jE) = 5.09559
      b_mtem(3,jA,jE) = -30.3295
      b_mtem(4,jA,jE) = 66.2975
      b_mtem(5,jA,jE) = -66.3458
      b_mtem(6,jA,jE) = 24.9443

! in HNO3
      jE = jhno3
      b_mtem(1,jA,jE) = 0.608309
      b_mtem(2,jA,jE) = -0.541905
      b_mtem(3,jA,jE) = -2.52084
      b_mtem(4,jA,jE) = 6.63297
      b_mtem(5,jA,jE) = -7.24599
      b_mtem(6,jA,jE) = 2.88811

! in HCl
      jE = jhcl
      b_mtem(1,jA,jE) = 1.98399
      b_mtem(2,jA,jE) = -4.51562
      b_mtem(3,jA,jE) = 8.36059
      b_mtem(4,jA,jE) = -12.4948
      b_mtem(5,jA,jE) = 9.67514
      b_mtem(6,jA,jE) = -3.18004


!----------
! Na3H(SO4)2 in E
      jA = jna3hso4

! in H2SO4
      jE = jh2so4
      b_mtem(1,jA,jE) = -0.83214
      b_mtem(2,jA,jE) = 4.99572
      b_mtem(3,jA,jE) = -20.1697
      b_mtem(4,jA,jE) = 41.4066
      b_mtem(5,jA,jE) = -42.2119
      b_mtem(6,jA,jE) = 16.4855

! in NH4HSO4
      jE = jnh4hso4
      b_mtem(1,jA,jE) = -0.65139
      b_mtem(2,jA,jE) = 3.52300
      b_mtem(3,jA,jE) = -22.8220
      b_mtem(4,jA,jE) = 56.2956
      b_mtem(5,jA,jE) = -59.9028
      b_mtem(6,jA,jE) = 23.1844

! in (NH4)3H(SO4)2
      jE = jlvcite
      b_mtem(1,jA,jE) = -1.31331
      b_mtem(2,jA,jE) = 8.40835
      b_mtem(3,jA,jE) = -38.1757
      b_mtem(4,jA,jE) = 80.5312
      b_mtem(5,jA,jE) = -79.8346
      b_mtem(6,jA,jE) = 30.0219

! in (NH4)2SO4
      jE = jnh4so4
      b_mtem(1,jA,jE) = -1.03054
      b_mtem(2,jA,jE) = 8.08155
      b_mtem(3,jA,jE) = -38.1046
      b_mtem(4,jA,jE) = 78.7168
      b_mtem(5,jA,jE) = -77.2263
      b_mtem(6,jA,jE) = 29.1521

! in NaHSO4
      jE = jnahso4
      b_mtem(1,jA,jE) = -1.90695
      b_mtem(2,jA,jE) = 11.6241
      b_mtem(3,jA,jE) = -50.3175
      b_mtem(4,jA,jE) = 105.884
      b_mtem(5,jA,jE) = -103.258
      b_mtem(6,jA,jE) = 37.6588

! in Na3H(SO4)2
      jE = jna3hso4
      b_mtem(1,jA,jE) = -0.34780
      b_mtem(2,jA,jE) = 2.85363
      b_mtem(3,jA,jE) = -17.6224
      b_mtem(4,jA,jE) = 38.9220
      b_mtem(5,jA,jE) = -39.8106
      b_mtem(6,jA,jE) = 15.6055

! in Na2SO4
      jE = jna2so4
      b_mtem(1,jA,jE) = -0.75230
      b_mtem(2,jA,jE) = 10.0140
      b_mtem(3,jA,jE) = -50.5677
      b_mtem(4,jA,jE) = 106.941
      b_mtem(5,jA,jE) = -105.534
      b_mtem(6,jA,jE) = 39.5196

! in HNO3
      jE = jhno3
      b_mtem(1,jA,jE) = 0.057456
      b_mtem(2,jA,jE) = -1.31264
      b_mtem(3,jA,jE) = -1.94662
      b_mtem(4,jA,jE) = 10.7024
      b_mtem(5,jA,jE) = -14.9946
      b_mtem(6,jA,jE) = 7.12161

! in HCl
      jE = jhcl
      b_mtem(1,jA,jE) = 0.637894
      b_mtem(2,jA,jE) = -2.29719
      b_mtem(3,jA,jE) = 0.765361
      b_mtem(4,jA,jE) = 4.8748
      b_mtem(5,jA,jE) = -9.25978
      b_mtem(6,jA,jE) = 4.91773
!
!
!
!----------------------------------------------------------
! Coefficients for %MDRH(T) = d1 + d2*T + d3*T^2 + d4*T^3    (T in Kelvin)
! valid Temperature Range: 240 - 320 K
!----------------------------------------------------------
!
! SULFATE-POOR SYSTEMS
! AC
      j_index = 1
      d_mdrh(j_index,1) = -58.00268351
      d_mdrh(j_index,2) = 2.031077573
      d_mdrh(j_index,3) = -0.008281218
      d_mdrh(j_index,4) = 1.00447E-05

! AN
      j_index = 2
      d_mdrh(j_index,1) = 1039.137773
      d_mdrh(j_index,2) = -11.47847095
      d_mdrh(j_index,3) = 0.047702786
      d_mdrh(j_index,4) = -6.77675E-05

! AS
      j_index = 3
      d_mdrh(j_index,1) = 115.8366357
      d_mdrh(j_index,2) = 0.491881663
      d_mdrh(j_index,3) = -0.00422807
      d_mdrh(j_index,4) = 7.29274E-06

! SC
      j_index = 4
      d_mdrh(j_index,1) = 253.2424151
      d_mdrh(j_index,2) = -1.429957864
      d_mdrh(j_index,3) = 0.003727554
      d_mdrh(j_index,4) = -3.13037E-06

! SN
      j_index = 5
      d_mdrh(j_index,1) = -372.4306506
      d_mdrh(j_index,2) = 5.3955633
      d_mdrh(j_index,3) = -0.019804438
      d_mdrh(j_index,4) = 2.25662E-05

! SS
      j_index = 6
      d_mdrh(j_index,1) = 286.1271416
      d_mdrh(j_index,2) = -1.670787758
      d_mdrh(j_index,3) = 0.004431373
      d_mdrh(j_index,4) = -3.57757E-06

! CC
      j_index = 7
      d_mdrh(j_index,1) = -1124.07059
      d_mdrh(j_index,2) = 14.26364209
      d_mdrh(j_index,3) = -0.054816822
      d_mdrh(j_index,4) = 6.70107E-05

! CN
      j_index = 8
      d_mdrh(j_index,1) = 1855.413934
      d_mdrh(j_index,2) = -20.29219473
      d_mdrh(j_index,3) = 0.07807482
      d_mdrh(j_index,4) = -1.017887858e-4

! AN + AC
      j_index = 9
      d_mdrh(j_index,1) = 1761.176886
      d_mdrh(j_index,2) = -19.29811062
      d_mdrh(j_index,3) = 0.075676987
      d_mdrh(j_index,4) = -1.0116959e-4

! AS + AC
      j_index = 10
      d_mdrh(j_index,1) = 122.1074303
      d_mdrh(j_index,2) = 0.429692122
      d_mdrh(j_index,3) = -0.003928277
      d_mdrh(j_index,4) = 6.43275E-06

! AS + AN
      j_index = 11
      d_mdrh(j_index,1) = 2424.634678
      d_mdrh(j_index,2) = -26.54031307
      d_mdrh(j_index,3) = 0.101625387
      d_mdrh(j_index,4) = -1.31544547798e-4

! AS + AN + AC
      j_index = 12
      d_mdrh(j_index,1) = 2912.082599
      d_mdrh(j_index,2) = -31.8894185
      d_mdrh(j_index,3) = 0.121185849
      d_mdrh(j_index,4) = -1.556534623e-4

! SC + AC
      j_index = 13
      d_mdrh(j_index,1) = 172.2596493
      d_mdrh(j_index,2) = -0.511006195
      d_mdrh(j_index,3) = 4.27244597e-4
      d_mdrh(j_index,4) = 4.12797E-07

! SN + AC
      j_index = 14
      d_mdrh(j_index,1) = 1596.184935
      d_mdrh(j_index,2) = -16.37945565
      d_mdrh(j_index,3) = 0.060281218
      d_mdrh(j_index,4) = -7.6161E-05

! SN + AN
      j_index = 15
      d_mdrh(j_index,1) = 1916.072988
      d_mdrh(j_index,2) = -20.85594868
      d_mdrh(j_index,3) = 0.081140141
      d_mdrh(j_index,4) = -1.07954274796e-4

! SN + AN + AC
      j_index = 16
      d_mdrh(j_index,1) = 1467.165935
      d_mdrh(j_index,2) = -16.01166196
      d_mdrh(j_index,3) = 0.063505582
      d_mdrh(j_index,4) = -8.66722E-05

! SN + SC
      j_index = 17
      d_mdrh(j_index,1) = 158.447059
      d_mdrh(j_index,2) = -0.628167358
      d_mdrh(j_index,3) = 0.002014448
      d_mdrh(j_index,4) = -3.13037E-06

! SN + SC + AC
      j_index = 18
      d_mdrh(j_index,1) = 1115.892468
      d_mdrh(j_index,2) = -11.76936534
      d_mdrh(j_index,3) = 0.045577399
      d_mdrh(j_index,4) = -6.05779E-05

! SS + AC
      j_index = 19
      d_mdrh(j_index,1) = 269.5432407
      d_mdrh(j_index,2) = -1.319963885
      d_mdrh(j_index,3) = 0.002592363
      d_mdrh(j_index,4) = -1.44479E-06

! SS + AN
      j_index = 20
      d_mdrh(j_index,1) = 2841.334784
      d_mdrh(j_index,2) = -31.1889487
      d_mdrh(j_index,3) = 0.118809274
      d_mdrh(j_index,4) = -1.53007e-4

! SS + AN + AC
      j_index = 21
      d_mdrh(j_index,1) = 2199.36914
      d_mdrh(j_index,2) = -24.11926569
      d_mdrh(j_index,3) = 0.092932361
      d_mdrh(j_index,4) = -1.21774e-4

! SS + AS
      j_index = 22
      d_mdrh(j_index,1) = 395.0051604
      d_mdrh(j_index,2) = -2.521101657
      d_mdrh(j_index,3) = 0.006139319
      d_mdrh(j_index,4) = -4.43756E-06

! SS + AS + AC
      j_index = 23
      d_mdrh(j_index,1) = 386.5150675
      d_mdrh(j_index,2) = -2.4632138
      d_mdrh(j_index,3) = 0.006139319
      d_mdrh(j_index,4) = -4.98796E-06

! SS + AS + AN
      j_index = 24
      d_mdrh(j_index,1) = 3101.538491
      d_mdrh(j_index,2) = -34.19978105
      d_mdrh(j_index,3) = 0.130118605
      d_mdrh(j_index,4) = -1.66873e-4

! SS + AS + AN + AC
      j_index = 25
      d_mdrh(j_index,1) = 2307.579403
      d_mdrh(j_index,2) = -25.43136774
      d_mdrh(j_index,3) = 0.098064728
      d_mdrh(j_index,4) = -1.28301e-4

! SS + SC
      j_index = 26
      d_mdrh(j_index,1) = 291.8309602
      d_mdrh(j_index,2) = -1.828912974
      d_mdrh(j_index,3) = 0.005053148
      d_mdrh(j_index,4) = -4.57516E-06

! SS + SC + AC
      j_index = 27
      d_mdrh(j_index,1) = 188.3914345
      d_mdrh(j_index,2) = -0.631345031
      d_mdrh(j_index,3) = 0.000622807
      d_mdrh(j_index,4) = 4.47196E-07

! SS + SN
      j_index = 28
      d_mdrh(j_index,1) = -167.1252839
      d_mdrh(j_index,2) = 2.969828002
      d_mdrh(j_index,3) = -0.010637255
      d_mdrh(j_index,4) = 1.13175E-05

! SS + SN + AC
      j_index = 29
      d_mdrh(j_index,1) = 1516.782768
      d_mdrh(j_index,2) = -15.7922661
      d_mdrh(j_index,3) = 0.058942209
      d_mdrh(j_index,4) = -7.5301E-05

! SS + SN + AN
      j_index = 30
      d_mdrh(j_index,1) = 1739.963163
      d_mdrh(j_index,2) = -19.06576022
      d_mdrh(j_index,3) = 0.07454963
      d_mdrh(j_index,4) = -9.94302E-05

! SS + SN + AN + AC
      j_index = 31
      d_mdrh(j_index,1) = 2152.104877
      d_mdrh(j_index,2) = -23.74998008
      d_mdrh(j_index,3) = 0.092256654
      d_mdrh(j_index,4) = -1.21953e-4

! SS + SN + SC
      j_index = 32
      d_mdrh(j_index,1) = 221.9976265
      d_mdrh(j_index,2) = -1.311331272
      d_mdrh(j_index,3) = 0.004406089
      d_mdrh(j_index,4) = -5.88235E-06

! SS + SN + SC + AC
      j_index = 33
      d_mdrh(j_index,1) = 1205.645615
      d_mdrh(j_index,2) = -12.71353459
      d_mdrh(j_index,3) = 0.048803922
      d_mdrh(j_index,4) = -6.41899E-05

! CC + AC
      j_index = 34
      d_mdrh(j_index,1) = 506.6737879
      d_mdrh(j_index,2) = -3.723520818
      d_mdrh(j_index,3) = 0.010814242
      d_mdrh(j_index,4) = -1.21087E-05

! CC + SC
      j_index = 35
      d_mdrh(j_index,1) = -1123.523841
      d_mdrh(j_index,2) = 14.08345977
      d_mdrh(j_index,3) = -0.053687823
      d_mdrh(j_index,4) = 6.52219E-05

! CC + SC + AC
      j_index = 36
      d_mdrh(j_index,1) = -1159.98607
      d_mdrh(j_index,2) = 14.44309169
      d_mdrh(j_index,3) = -0.054841073
      d_mdrh(j_index,4) = 6.64259E-05

! CN + AC
      j_index = 37
      d_mdrh(j_index,1) = 756.0747916
      d_mdrh(j_index,2) = -8.546826257
      d_mdrh(j_index,3) = 0.035798677
      d_mdrh(j_index,4) = -5.06629E-05

! CN + AN
      j_index = 38
      d_mdrh(j_index,1) = 338.668191
      d_mdrh(j_index,2) = -2.971223403
      d_mdrh(j_index,3) = 0.012294866
      d_mdrh(j_index,4) = -1.87558E-05

! CN + AN + AC
      j_index = 39
      d_mdrh(j_index,1) = -53.18033508
      d_mdrh(j_index,2) = 0.663911748
      d_mdrh(j_index,3) = 9.16326e-4
      d_mdrh(j_index,4) = -6.70354E-06

! CN + SC
      j_index = 40
      d_mdrh(j_index,1) = 3623.831129
      d_mdrh(j_index,2) = -39.27226457
      d_mdrh(j_index,3) = 0.144559515
      d_mdrh(j_index,4) = -1.78159e-4

! CN + SC + AC
      j_index = 41
      d_mdrh(j_index,1) = 3436.656743
      d_mdrh(j_index,2) = -37.16192684
      d_mdrh(j_index,3) = 0.136641377
      d_mdrh(j_index,4) = -1.68262e-4

! CN + SN
      j_index = 42
      d_mdrh(j_index,1) = 768.608476
      d_mdrh(j_index,2) = -8.051517149
      d_mdrh(j_index,3) = 0.032342332
      d_mdrh(j_index,4) = -4.52224E-05

! CN + SN + AC
      j_index = 43
      d_mdrh(j_index,1) = 33.58027951
      d_mdrh(j_index,2) = -0.308772182
      d_mdrh(j_index,3) = 0.004713639
      d_mdrh(j_index,4) = -1.19658E-05

! CN + SN + AN
      j_index = 44
      d_mdrh(j_index,1) = 57.80183041
      d_mdrh(j_index,2) = 0.215264604
      d_mdrh(j_index,3) = 4.11406e-4
      d_mdrh(j_index,4) = -4.30702E-06

! CN + SN + AN + AC
      j_index = 45
      d_mdrh(j_index,1) = -234.368984
      d_mdrh(j_index,2) = 2.721045204
      d_mdrh(j_index,3) = -0.006688341
      d_mdrh(j_index,4) = 2.31729E-06

! CN + SN + SC
      j_index = 46
      d_mdrh(j_index,1) = 3879.080557
      d_mdrh(j_index,2) = -42.13562874
      d_mdrh(j_index,3) = 0.155235005
      d_mdrh(j_index,4) = -1.91387e-4

! CN + SN + SC + AC
      j_index = 47
      d_mdrh(j_index,1) = 3600.576985
      d_mdrh(j_index,2) = -39.0283489
      d_mdrh(j_index,3) = 0.143710316
      d_mdrh(j_index,4) = -1.77167e-4

! CN + CC
      j_index = 48
      d_mdrh(j_index,1) = -1009.729826
      d_mdrh(j_index,2) = 12.9145339
      d_mdrh(j_index,3) = -0.049811146
      d_mdrh(j_index,4) = 6.09563E-05

! CN + CC + AC
      j_index = 49
      d_mdrh(j_index,1) = -577.0919514
      d_mdrh(j_index,2) = 8.020324227
      d_mdrh(j_index,3) = -0.031469556
      d_mdrh(j_index,4) = 3.82181E-05

! CN + CC + SC
      j_index = 50
      d_mdrh(j_index,1) = -728.9983499
      d_mdrh(j_index,2) = 9.849458215
      d_mdrh(j_index,3) = -0.03879257
      d_mdrh(j_index,4) = 4.78844E-05

! CN + CC + SC + AC
      j_index = 51
      d_mdrh(j_index,1) = -803.7026845
      d_mdrh(j_index,2) = 10.61881494
      d_mdrh(j_index,3) = -0.041402993
      d_mdrh(j_index,4) = 5.08084E-05

!
! SULFATE-RICH SYSTEMS
! AB
      j_index = 52
      d_mdrh(j_index,1) = -493.6190458
      d_mdrh(j_index,2) = 6.747053851
      d_mdrh(j_index,3) = -0.026955267
      d_mdrh(j_index,4) = 3.45118E-05

! LV
      j_index = 53
      d_mdrh(j_index,1) = 53.37874093
      d_mdrh(j_index,2) = 1.01368249
      d_mdrh(j_index,3) = -0.005887513
      d_mdrh(j_index,4) = 8.94393E-06

! SB
      j_index = 54
      d_mdrh(j_index,1) = 206.619047
      d_mdrh(j_index,2) = -1.342735684
      d_mdrh(j_index,3) = 0.003197691
      d_mdrh(j_index,4) = -1.93603E-06

! AB + LV
      j_index = 55
      d_mdrh(j_index,1) = -493.6190458
      d_mdrh(j_index,2) = 6.747053851
      d_mdrh(j_index,3) = -0.026955267
      d_mdrh(j_index,4) = 3.45118E-05

! AS + LV
      j_index = 56
      d_mdrh(j_index,1) = 53.37874093
      d_mdrh(j_index,2) = 1.01368249
      d_mdrh(j_index,3) = -0.005887513
      d_mdrh(j_index,4) = 8.94393E-06

! SS + SB
      j_index = 57
      d_mdrh(j_index,1) = 206.619047
      d_mdrh(j_index,2) = -1.342735684
      d_mdrh(j_index,3) = 0.003197691
      d_mdrh(j_index,4) = -1.93603E-06

! SS + LV
      j_index = 58
      d_mdrh(j_index,1) = 41.7619047
      d_mdrh(j_index,2) = 1.303872053
      d_mdrh(j_index,3) = -0.007647908
      d_mdrh(j_index,4) = 1.17845E-05

! SS + AS + LV
      j_index = 59
      d_mdrh(j_index,1) = 41.7619047
      d_mdrh(j_index,2) = 1.303872053
      d_mdrh(j_index,3) = -0.007647908
      d_mdrh(j_index,4) = 1.17845E-05

! SS + AB
      j_index = 60
      d_mdrh(j_index,1) = -369.7142842
      d_mdrh(j_index,2) = 5.512878771
      d_mdrh(j_index,3) = -0.02301948
      d_mdrh(j_index,4) = 3.0303E-05

! SS + LV + AB
      j_index = 61
      d_mdrh(j_index,1) = -369.7142842
      d_mdrh(j_index,2) = 5.512878771
      d_mdrh(j_index,3) = -0.02301948
      d_mdrh(j_index,4) = 3.0303E-05

! SB + AB
      j_index = 62
      d_mdrh(j_index,1) = -162.8095232
      d_mdrh(j_index,2) = 2.399326592
      d_mdrh(j_index,3) = -0.009336219
      d_mdrh(j_index,4) = 1.17845E-05

! SS + SB + AB
      j_index = 63
      d_mdrh(j_index,1) = -735.4285689
      d_mdrh(j_index,2) = 8.885521857
      d_mdrh(j_index,3) = -0.033488456
      d_mdrh(j_index,4) = 4.12458E-05


      endif ! first

      return
      end subroutine load_mosaic_parameters











!***********************************************************************
! updates all temperature dependent thermodynamic parameters
!
! author: Rahul A. Zaveri
! update: jan 2005
!-----------------------------------------------------------------------
      subroutine update_thermodynamic_constants	! TOUCH
	use module_data_mosaic_main, only:  &
	    icase_soa, imodel, iraoultslaw, ireactiveuptake
      use module_data_mosaic_aero

      implicit none

! local variables
      integer iv, j_index, ibin, je
      real(r8) :: tr, rt, term
      real(r8) :: gam_nh4no3_0, gam_nh4cl_0, m_nh4no3_0, m_nh4cl_0  ! 6/25/2008
      real(r8) :: sat_factor
! function
      real(r8) :: fn_Keq, fn_Po, drh_mutual, bin_molality, molality_0, MWsoa  ! 6/25/2008


      tr = 298.15			! reference temperature
      rt = 82.056*T_K/(1.e9*1.e6)	! [m^3 atm/nmol]

! gas-liquid
      Keq_gl(1)= 1.0				        ! Kelvin Effect (default)
      Keq_gl(2)= fn_Keq(57.64d0, 13.79d0, -5.39d0,T_K)*rt     ! NH3(g)  <=> NH3(l)
      Keq_gl(3)= fn_Keq(2.63d6,  29.17d0, 16.83d0,T_K)*rt     ! HNO3(g) <=> NO3- + H+
      Keq_gl(4)= fn_Keq(2.00d6,  30.20d0, 19.91d0,T_K)*rt     ! HCl(g)  <=> Cl- + H+

! liquid-liquid
      Keq_ll(1)= fn_Keq(1.0502d-2, 8.85d0, 25.14d0,T_K)      ! HSO4- <=> SO4= + H+
      Keq_ll(2)= fn_Keq(1.805d-5, -1.50d0, 26.92d0,T_K)      ! NH3(l) + H2O = NH4+ + OH-
      Keq_ll(3)= fn_Keq(1.01d-14,-22.52d0, 26.92d0,T_K)      ! H2O(l) <=> H+ + OH-


      Kp_nh3   = Keq_ll(3)/(Keq_ll(2)*Keq_gl(2))
      Kp_nh4no3= Kp_nh3/Keq_gl(3)
      Kp_nh4cl = Kp_nh3/Keq_gl(4)


! solid-gas
      Keq_sg(1)= fn_Keq(4.72d-17,-74.38d0,6.12d0,T_K)/rt**2  ! NH4NO3<=>NH3(g)+HNO3(g)
      Keq_sg(2)= fn_Keq(8.43d-17,-71.00d0,2.40d0,T_K)/rt**2  ! NH4Cl <=>NH3(g)+HCl(g)


! solid-liquid
      Keq_sl(jnh4so4) = fn_Keq(1.040d0,-2.65d0, 38.57d0, T_K)  ! amSO4(s) = 2NH4+ + SO4=
      Keq_sl(jlvcite) = fn_Keq(11.8d0, -5.19d0, 54.40d0, T_K)  ! lvcite(s)= 3NH4+ + HSO4- + SO4=
      Keq_sl(jnh4hso4)= fn_Keq(117.0d0,-2.87d0, 15.83d0, T_K)  ! amHSO4(s)= NH4+ + HSO4-
      Keq_sl(jnh4msa) = 1.e15				       ! NH4MSA(s)= NH4+ + MSA-
      Keq_sl(jnh4no3) = fn_Keq(12.21d0,-10.4d0, 17.56d0, T_K)  ! NH4NO3(s)= NH4+ + NO3-
      Keq_sl(jnh4cl)  = fn_Keq(17.37d0,-6.03d0, 16.92d0, T_K)  ! NH4Cl(s) = NH4+ + Cl-
      Keq_sl(jna2so4) = fn_Keq(0.491d0, 0.98d0, 39.75d0, T_K)  ! Na2SO4(s)= 2Na+ + SO4=
      Keq_sl(jnahso4) = fn_Keq(313.0d0, 0.8d0,  14.79d0, T_K)  ! NaHSO4(s)= Na+ + HSO4-
      Keq_sl(jna3hso4)= 1.e15		 	               ! Na3H(SO4)2(s) = 2Na+ + HSO4- + SO4=
      Keq_sl(jnamsa)  = 1.e15				       ! NaMSA(s) = Na+ + MSA-
      Keq_sl(jnano3)  = fn_Keq(11.95d0,-8.22d0, 16.01d0, T_K)  ! NaNO3(s) = Na+ + NO3-
      Keq_sl(jnacl)   = fn_Keq(38.28d0,-1.52d0, 16.89d0, T_K)  ! NaCl(s)  = Na+ + Cl-
      Keq_sl(jcacl2)  = fn_Keq(8.0d11,  32.84d0,44.79d0, T_K)  ! CaCl2(s) = Ca++ + 2Cl-
      Keq_sl(jcano3)  = fn_Keq(4.31d5,   7.83d0,42.01d0, T_K)  ! Ca(NO3)2(s) = Ca++ + 2NO3-
      Keq_sl(jcamsa2) = 1.e15				 ! CaMSA2(s)= Ca+ + 2MSA-

! vapor pressures of VBS species, wkc added 20181022
!      Po_soa(iCn3_g) = fn_Po(1.834d-8, 156.0d0, T_K) ! [Pascal]
!      Po_soa(iCn2_g) = fn_Po(1.834d-7, 156.0d0, T_K) ! [Pascal]
!      Po_soa(iCn1_g) = fn_Po(1.834d-6, 156.0d0, T_K) ! [Pascal]
!      Po_soa(iC0_g) = fn_Po(1.834d-5, 156.0d0, T_K) ! [Pascal]
!      Po_soa(iC1_g) = fn_Po(1.834d-4, 156.0d0, T_K) ! [Pascal]
!      Po_soa(iC2_g) = fn_Po(1.834d-3, 156.0d0, T_K) ! [Pascal]
!      Po_soa(iC3_g) = fn_Po(1.834d-2, 156.0d0, T_K) ! [Pascal]
!      Po_soa(iC4_g) = fn_Po(1.834d-1, 156.0d0, T_K) ! [Pascal]
!      Po_soa(iC5_g) = fn_Po(1.834d0, 156.0d0, T_K) ! [Pascal]
!      Po_soa(iC6_g) = fn_Po(1.834d1, 156.0d0, T_K) ! [Pascal]
!      Po_soa(iC7_g) = fn_Po(1.834d2, 156.0d0, T_K) ! [Pascal]
!      Po_soa(iC8_g) = fn_Po(1.834d3, 156.0d0, T_K) ! [Pascal]
!      Po_soa(iC9_g) = fn_Po(1.834d4, 156.0d0, T_K) ! [Pascal]

!      Po_soa(iCn3_g) = fn_Po(2.595d-8, 156.0d0, T_K) ! [Pascal]
!      Po_soa(iCn2_g) = fn_Po(3.037d-7, 156.0d0, T_K) ! [Pascal]
!      Po_soa(iCn1_g) = fn_Po(3.662d-6, 156.0d0, T_K) ! [Pascal]
!      Po_soa(iC0_g) = fn_Po(1.056d-5, 156.0d0, T_K) ! [Pascal]
!      Po_soa(iC1_g) = fn_Po(1.123d-4, 156.0d0, T_K) ! [Pascal]
!      Po_soa(iC2_g) = fn_Po(1.199d-3, 156.0d0, T_K) ! [Pascal]
!      Po_soa(iC3_g) = fn_Po(1.285d-2, 156.0d0, T_K) ! [Pascal]
!      Po_soa(iC4_g) = fn_Po(1.385d-1, 156.0d0, T_K) ! [Pascal]
!      Po_soa(iC5_g) = fn_Po(1.502d0, 156.0d0, T_K) ! [Pascal]
!      Po_soa(iC6_g) = fn_Po(1.640d1, 156.0d0, T_K) ! [Pascal]
!      Po_soa(iC7_g) = fn_Po(1.806d2, 156.0d0, T_K) ! [Pascal]
!      Po_soa(iC8_g) = fn_Po(2.010d3, 156.0d0, T_K) ! [Pascal]
!      Po_soa(iC9_g) = fn_Po(2.265d4, 156.0d0, T_K) ! [Pascal]

      Po_soa(iCn3_g) = 10.0**-3 ! [ug/m^3(air)]
      Po_soa(iCn2_g) = 10.0**-2 ! [ug/m^3(air)]
      Po_soa(iCn1_g) = 10.0**-1 ! [ug/m^3(air)]
      Po_soa(iC0_g) = 10.0**0 ! [ug/m^3(air)]
      Po_soa(iC1_g) = 10.0**1 ! [ug/m^3(air)]
      Po_soa(iC2_g) = 10.0**2 ! [ug/m^3(air)]
      Po_soa(iC3_g) = 10.0**3 ! [ug/m^3(air)]
      Po_soa(iC4_g) = 10.0**4 ! [ug/m^3(air)]
      Po_soa(iC5_g) = 10.0**5 ! [ug/m^3(air)]
      Po_soa(iC6_g) = 10.0**6 ! [ug/m^3(air)]
      Po_soa(iC7_g) = 10.0**7 ! [ug/m^3(air)].
      Po_soa(iC8_g) = 10.0**8 ! [ug/m^3(air)].
      Po_soa(iC9_g) = 10.0**9 ! [ug/m^3(air)]

      sat_factor = 1.0
      do iv = icn3_g, ngas_volatile
        sat_soa(iv) = sat_factor * 1e3*Po_soa(iv)/mw_gas_mac(iv) ! [nmol/m^3(air)] Remove this and uncomment following line edit: wkc

        !sat_soa(iv) = sat_factor * 1.e9*Po_soa(iv)/(8.314*T_K)  ! [nmol/m^3(air)]
      enddo


! water surface tension
      term = (647.15 - T_K)/647.15
      sigma_water = 0.2358*term**1.256 * (1. - 0.625*term) ! surface tension of pure water in N/m

! MDRH(T)
      do j_index = 1, 63
        MDRH_T(j_index) = drh_mutual(j_index)
      enddo



! RH dependent parameters
      do ibin = 1, nbin_a
        aH2O_a(ibin) = aH2O			! initialize
      enddo

      do je = 1, nelectrolyte
        molality0(je) = bin_molality(je,1)	! compute aH2O dependent binary molalities  EFFI
      enddo

!      write(6,*)'aH2O = ', aH2O
!      write(6,*)'molality of NH4SO4 = ', molality0(jnh4so4)

      call MTEM_compute_log_gamZ		! function of aH2O and T

! 6/25/2008 - start
      gam_nh4no3_0 = 10.**log_gamZ(jnh4no3,jnh4no3)
      gam_nh4cl_0  = 10.**log_gamZ(jnh4cl,jnh4cl)

      m_nh4no3_0   = molality_0(jnh4no3,aH2O)
      m_nh4cl_0    = molality_0(jnh4cl, aH2O)

      Kp_nh4no3_0  = Kp_nh4no3*(m_nh4no3_0*gam_nh4no3_0)**2
      Kp_nh4cl_0   = Kp_nh4cl *(m_nh4cl_0 *gam_nh4cl_0 )**2
! 6/25/2008 - end


      return
      end subroutine update_thermodynamic_constants










!***********************************************************************
! computes mass transfer coefficients for each condensing species for
! all the aerosol bins
!
! author: Rahul A. Zaveri
! update: jan 2005
!-----------------------------------------------------------------------
      subroutine aerosolmtc		! TOUCH
!      include 'v33com9a'
	use module_data_mosaic_main, only:  &
	    it, icase_soa, iurbanleg, imodel, iraoultslaw, ireactiveuptake	! RAZ 12/17/2013: added it
      use module_data_mosaic_aero

      implicit none

! local variables
      integer nghq
      parameter (nghq = 2)		! gauss-hermite quadrature order
      integer ibin, iq, iv
      real(r8) :: tworootpi, root2, beta
      parameter (tworootpi = 3.5449077, root2 = 1.4142135, beta = 2.0)
      real(r8) :: cdum, Dp, Dp_avg, Fkn, Kn, lnsg, lnDpgn, lnDp, speed, sumghq
      real(r8) :: xghq(nghq), wghq(nghq)				! quadrature abscissae and weights
      real(r8) :: v_molar(ngas_volatile) 				! molar vols of inorganic species
      real(r8) :: freepath(ngas_volatile),Dg(ngas_volatile)	  	! keep local
      real(r8) :: kg_mtc, kp_mtc, cstar_by_sumA				! keep local
      real(r8) :: fuchs_sutugin						! mosaic func
      real(r8) :: gas_diffusivity					! mosaic func
      real(r8) :: mean_molecular_speed					! mosaic func
      real(r8) :: q_dum, Db						! keep local
      real(r8) :: coth
      real(r8) :: dum_fac_Dl
      real(r8) :: Dp_cutoff	! [nm]
      real(r8) :: m_c5, m_c7, kf, kc_c5, kc_c7, Vp, sum_oa_mass


      logical first
      save first
      data first/.true./


      v_molar(ih2so4_g)= 42.88
      v_molar(ihno3_g) = 24.11
      v_molar(ihcl_g)  = 21.48
      v_molar(inh3_g)  = 14.90
      v_molar(imsa_g)  = 58.00




!	do iv = 1, ngas_volatile
!	  accom(iv) = 0.01
!	enddo

! quadrature weights
      xghq(1) =  0.70710678
      xghq(2) = -0.70710678
      wghq(1) =  0.88622693
      wghq(2) =  0.88622693



! specify gas-phase and bulk-phase diffusivity and mean free path for condensing gases
! ioa
      do iv = 1, ngas_ioa
        speed  = mean_molecular_speed(T_K,mw_gas_mac(iv))	! cm/s
        Dg(iv) = gas_diffusivity(T_K,P_atm,mw_gas_mac(iv),v_molar(iv)) ! cm^2/s
        freepath(iv) = 3.*Dg(iv)/speed			! cm
      enddo

! soa
      do iv = icn3_g, ngas_volatile
        speed = mean_molecular_speed(T_K,mw_gas_mac(iv))	! cm/s
	Dg(iv) = 0.05					! gas diffusivity (cm^2/s)
	Dg(iv) = 1.38e-01 * 44.0 / mw_gas_mac(iv) ! edit wkc, may remove if desired (cm^2/s)
	freepath(iv) = 3.*Dg(iv)/speed
      enddo



! calc mass transfer coefficients for gases over various aerosol bins

      if (mSIZE_FRAMEWORK .eq. mMODAL) then

! for modal approach (DOES NOT include bulk diffusion limitation treatment)
      do 10 ibin = 1, nbin_a

        if(jaerosolstate(ibin) .eq. no_aerosol)goto 10
        call calc_dry_n_wet_aerosol_props(ibin)

        lnsg   = log(sigmag_a(ibin))
! lnDpgn = ln( wet geometric-mean Dp of number distribution )
        lnDpgn = log(Dp_wet_a(ibin)) - 1.5*lnsg*lnsg

        cdum   = tworootpi*num_a(ibin)*   &
                 exp(beta*lnDpgn + 0.5*(beta*lnsg)**2)

        do 20 iv = 1, ngas_volatile

          sumghq = 0.0
          do 30 iq = 1, nghq	! sum over gauss-hermite quadrature points
            lnDp = lnDpgn + beta*lnsg**2 + root2*lnsg*xghq(iq)
            Dp = exp(lnDp)
            Kn = 2.*freepath(iv)/Dp
            Fkn = fuchs_sutugin(Kn,accom(iv))
            sumghq = sumghq + wghq(iq)*Dp*Fkn/(Dp**beta)
30        continue

        Kg(iv,ibin) = cdum*Dg(iv)*sumghq		! 1/s

20      continue
10    continue

      elseif ((mSIZE_FRAMEWORK .eq. mSECTIONAL   ) .or. &
              (mSIZE_FRAMEWORK .eq. mUNSTRUCTURED)) then


! for sectional approach (includes bulk diffusion limitation treatment)


      Dp_cutoff = 160.0	! nm


      do 11 ibin = 1, nbin_a

        if(jaerosolstate(ibin) .eq. no_aerosol)goto 11

        call calc_dry_n_wet_aerosol_props(ibin)

        cdum  = 3.14159*(Dp_wet_a(ibin)**2)*num_a(ibin)

! specify size-dependent bulk diffusivity Db (cm2/s)

        if(Dp_wet_a(ibin)*1.e7 .lt. Dp_cutoff)then
          Db = 1.e-6	! cm2/s
        else
          Db = 1.e-15	! cm2/s
        endif
          Db = 1.e-10 ! edit wkc, remove!
        do 21 iv = 1, ngas_volatile
          Kn = 2.*freepath(iv)/Dp_wet_a(ibin)		! Knudsen number
          Fkn = fuchs_sutugin(Kn,accom(iv))		! transition regime correction factor (Fuchs and Sutugin, 1971)

          kg_mtc = (2.*Dg(iv)/Dp_wet_a(ibin))*Fkn	! gas-side mtc (cm/s)
	  Kg(iv,ibin) = cdum*kg_mtc			! overall gas-side mt rate constant (1/s), ignoring bulk diffusion limitation for inorganic species
      QQ(iv,ibin)=1.0
	  if(iv .lt. icn3_g)goto 21			! no bulk diffusion limitation for inorganic species
	  !if(gas(iv) .gt. sat_soa(iv))goto 21		! no bulk diffusion limitation for supersaturated organic species

	  ! apply bulk diffusion limitation for subsaturated organic species based on Zaveri et al. (2018) ES&T

	  kc_firstorder(iv) = 0.0				! specify first-order chemical reaction rate constant (1/s)

	  q_dum  = 0.5*Dp_wet_a(ibin)*sqrt(kc_firstorder(iv)/Db)	! dimensionless diffuso-reactive parameter
	  if(q_dum .lt. 0.1)then
	    QQ(iv,ibin) = 1.0					 ! USE THIS IN ASTEM
	  else
	    QQ(iv,ibin) = 3.0*(q_dum/tanh(q_dum) - 1.0)/q_dum**2 ! USE THIS IN ASTEM
	  endif

	  if(kc_firstorder(iv) .ge. 0.01)then	! APPROXIMATION 1

	    Kg(iv,ibin) = cdum*kg_mtc			! overall gas-side mt rate constant (1/s)

	  else ! kc_firstorder(iv) < 0.01	! APPROXIMATION 2

	    if(q_dum .lt. 0.1)then
	      kp_mtc = 10.*Db/Dp_wet_a(ibin)
	    else
	      kp_mtc = (2.0*Db/Dp_wet_a(ibin))*(q_dum/tanh(q_dum)-1.0)/(1.0-QQ(iv,ibin))
	    endif
	    cstar_by_sumA = sat_soa(iv)*mw_aer_mac(iv)/(dens_aer_mac(ioc_a)*1.e15)   	! approximation to actual C*/sumA
	    Kg(iv,ibin) = cdum*kg_mtc*kp_mtc/(kp_mtc + cstar_by_sumA*kg_mtc)		! overall gas-side mt rate constant (1/s)

          endif

21      continue ! iv loop

11    continue	 ! bin loop


      else

        write(6,*)'Error in the choice of mSIZE_FRAMEWORK'
        write(6,*)'Stopping in subr. aerosolmtc'
        stop

      endif


      return
      return
      return
      end subroutine aerosolmtc








!***********************************************************************
! calculates dry and wet aerosol properties: density, refractive indices
!
! author: Rahul A. Zaveri
! update: jan 2005
!-----------------------------------------------------------------------
      subroutine calc_dry_n_wet_aerosol_props(ibin)
!      include 'v33com9a'
      use module_data_mosaic_main, only:  &
          maeroptic, mshellcore, piover4, piover6, third
      use module_data_mosaic_aero
      use module_data_mosaic_asect, only:  &
          dcen_sect, isize_of_ibin, itype_of_ibin

      implicit none

! subr arguments
      integer ibin
! local variables
      integer isize, itype, jc, je, iaer
      real(r8) :: aer_H, duma, vol_core, vol_shell, vol_dum
      complex rixvol_tot, rixvol_core, rixvol_shell


! calculate dry mass and dry volume of a bin
      mass_dry_a(ibin) = 0.0		! initialize to 0.0
      vol_dry_a(ibin)  = 0.0		! initialize to 0.0
      area_dry_a(ibin) = 0.0		! initialize to 0.0

      if(jaerosolstate(ibin) .ne. no_aerosol)then

         aer_H = (2.*aer(iso4_a,jtotal,ibin) +   &
                     aer(ino3_a,jtotal,ibin) +   &
                     aer(icl_a,jtotal,ibin)  +   &
                     aer(imsa_a,jtotal,ibin) +   &
                  2.*aer(ico3_a,jtotal,ibin))-   &
                 (2.*aer(ica_a,jtotal,ibin)  +   &
                     aer(ina_a,jtotal,ibin)  +   &
                     aer(inh4_a,jtotal,ibin))
         aer_H = max(aer_H, 0.0d0)

         do iaer = 1, naer
            mass_dry_a(ibin) = mass_dry_a(ibin) +   &
                               aer(iaer,jtotal,ibin)*mw_aer_mac(iaer)	! ng/m^3(air)
            vol_dry_a(ibin) = vol_dry_a(ibin) +   &
            aer(iaer,jtotal,ibin)*mw_aer_mac(iaer)/dens_aer_mac(iaer)  	! ncc/m^3(air)
         enddo
         mass_dry_a(ibin) = mass_dry_a(ibin) + aer_H
         vol_dry_a(ibin) = vol_dry_a(ibin) + aer_H

         mass_dry_a(ibin) = mass_dry_a(ibin)*1.e-15			! g/cc(air)
         vol_dry_a(ibin) = vol_dry_a(ibin)*1.e-15			! cc(aer)/cc(air)

! wet mass and wet volume
         mass_wet_a(ibin) = mass_dry_a(ibin) + water_a(ibin)*1.e-3	! g/cc(air)
         vol_wet_a(ibin)  = vol_dry_a(ibin) + water_a(ibin)*1.e-3	! cc(aer)/cc(air)

! calculate mean dry and wet particle densities
         dens_dry_a(ibin) = mass_dry_a(ibin)/vol_dry_a(ibin) 		! g/cc(aerosol)
         dens_wet_a(ibin) = mass_wet_a(ibin)/vol_wet_a(ibin) 		! g/cc(aerosol)

! calculate mean dry and wet particle diameters

         !Dp_dry_a(ibin)=Dp_dry_a(ibin)!(vol_dry_a(ibin)/(piover6*num_a(ibin)))**third	! cm


         Dp_dry_a(ibin)=(vol_dry_a(ibin)/(piover6*num_a(ibin)))**third	! cm

         
         Dp_wet_a(ibin)=(vol_wet_a(ibin)/(piover6*num_a(ibin)))**third	! cm

! calculate mean dry and wet particle surface areas
         area_dry_a(ibin)= piover4*num_a(ibin)*Dp_dry_a(ibin)**2	! cm^2/cc(air)
         area_wet_a(ibin)= piover4*num_a(ibin)*Dp_wet_a(ibin)**2	! cm^2/cc(air)

! calculate volume average refractive index
!   load comp_a array with component mass concentrations

! rahul had turned this off, but it is needed
!        if(1 == 1)go to 100		! TEMP
         if (maeroptic <= 0) goto 100

         do je = 1, nelectrolyte
            comp_a(je)=electrolyte(je,jtotal,ibin)*mw_comp_a(je)*1.e-15	! g/cc(air)
         enddo
         comp_a(joc)  = aer(ioc_a,  jtotal,ibin)*mw_comp_a(joc  )*1.e-15	! g/cc(air)
         comp_a(jbc)  = aer(ibc_a,  jtotal,ibin)*mw_comp_a(jbc  )*1.e-15	! g/cc(air)
         comp_a(join) = aer(ioin_a, jtotal,ibin)*mw_comp_a(join )*1.e-15	! g/cc(air)
		 comp_a(jcn3)= aer(icn3_a,jtotal,ibin)*mw_comp_a(jcn3)*1.e-15	! g/cc(air)
		 comp_a(jcn2)= aer(icn2_a,jtotal,ibin)*mw_comp_a(jcn2)*1.e-15	! g/cc(air)
		 comp_a(jcn1)= aer(icn1_a,jtotal,ibin)*mw_comp_a(jcn1)*1.e-15	! g/cc(air)
		 comp_a(jc0)= aer(ic0_a,jtotal,ibin)*mw_comp_a(jc0)*1.e-15	! g/cc(air)
		 comp_a(jc1) = aer(ic1_a,jtotal,ibin)*mw_comp_a(jc1)*1.e-15	! g/cc(air)
		 comp_a(jc2) = aer(ic2_a,jtotal,ibin)*mw_comp_a(jc2)*1.e-15	! g/cc(air)
		 comp_a(jc3)= aer(ic3_a,jtotal,ibin)*mw_comp_a(jc3)*1.e-15	! g/cc(air)
		 comp_a(jc4)= aer(ic4_a,jtotal,ibin)*mw_comp_a(jc4)*1.e-15	! g/cc(air)
		 comp_a(jc5)= aer(ic5_a,jtotal,ibin)*mw_comp_a(jc5)*1.e-15	! g/cc(air)
		 comp_a(jc6)= aer(ic6_a,jtotal,ibin)*mw_comp_a(jc6)*1.e-15	! g/cc(air)
	 	 comp_a(jc7)= aer(ic7_a,jtotal,ibin)*mw_comp_a(jc7)*1.e-15	! g/cc(air)
	 	 comp_a(jc8)= aer(ic8_a,jtotal,ibin)*mw_comp_a(jc8)*1.e-15	! g/cc(air)
         comp_a(jc9)= aer(ic9_a,jtotal,ibin)*mw_comp_a(jc9)*1.e-15	! g/cc(air)

         comp_a(jh2o) = water_a(ibin)*1.e-3				! g/cc(air)

         rixvol_tot   = (0.0,0.0)
         do jc = 1, naercomp
            comp_a(jc) = max( 0.0d0, comp_a(jc) )
            rixvol_tot = rixvol_tot   &
                       + ref_index_a(jc)*comp_a(jc)/dens_comp_a(jc)
         enddo
         ri_avg_a(ibin) = rixvol_tot/vol_wet_a(ibin)

!
! shell/core calcs - first set values to default (corresponding to zero core)
!
         ri_shell_a(ibin) = ri_avg_a(ibin)
         ri_core_a(ibin)  = (0.0,0.0)
         Dp_core_a(ibin)  = 0.0

! sum ri*vol and vol for core species (bc and optionally oin=dust)
! currently just bc in core, but what about insoluble oin and dust species ???
         jc = jbc
         rixvol_core  = ref_index_a(jc)*comp_a(jc)/dens_comp_a(jc)
         vol_core = comp_a(jc)/dens_comp_a(jc)
         vol_core = max( 0.0d0, min( vol_core, vol_wet_a(ibin) ) )

! neglect core if (core volume) < 1.0d-9*(total volume)
!              or (core volume) < 1.0d-22 cm3 = (0.58 nm)**3
! neglect shell using similar criteria
         vol_dum = max( 1.0d-22, 1.0d-9*vol_wet_a(ibin) )
         vol_shell = vol_wet_a(ibin) - vol_core
         if (vol_core >= vol_dum) then
            if (vol_shell < vol_dum) then
               ri_shell_a(ibin)  = (0.0,0.0)
               ri_core_a(ibin) = ri_avg_a(ibin)
               Dp_core_a(ibin) = Dp_wet_a(ibin)
            else
               ri_core_a(ibin) = rixvol_core/vol_core
               Dp_core_a(ibin) = Dp_wet_a(ibin)   &
                               * (vol_core/vol_wet_a(ibin))**third

               if (vol_shell >= vol_dum) then
                  rixvol_shell = rixvol_tot - rixvol_core
                  ri_shell_a(ibin) = rixvol_shell/vol_shell
               else
                  ri_shell_a(ibin) = (0.0,0.0)
               endif
            endif
         endif

      else
! use defaults when (jaerosolstate(ibin) .eq. no_aerosol)

         dens_dry_a(ibin) = 1.0	 ! g/cc(aerosol)
         dens_wet_a(ibin) = 1.0	 ! g/cc(aerosol)
!        Dp_dry_a(ibin) = dcen_sect(ibin)	! cm
!        Dp_wet_a(ibin) = dcen_sect(ibin)	! cm
         if (msize_framework == msectional) then
            isize = isize_of_ibin(ibin)
            itype = itype_of_ibin(ibin)
            Dp_dry_a(ibin) = dcen_sect(isize,itype)
            Dp_wet_a(ibin) = Dp_dry_a(ibin)
         end if

         ri_avg_a(ibin) = (1.5,0.0)
         ri_shell_a(ibin) = (1.5,0.0)
         ri_core_a(ibin)  = (0.0,0.0)
         Dp_core_a(ibin)  = 0.0

      endif   ! if(jaerosolstate(ibin) .ne. no_aerosol)then


100   continue

      return
      end subroutine calc_dry_n_wet_aerosol_props






















!***********************************************************************
! functions used in MOSAIC
!
! author: Rahul A. Zaveri
! update: jan 2005
!-----------------------------------------------------------------------



!----------------------------------------------------------
      function fn_Keq(Keq_298, a, b, T)
      use module_data_mosaic_kind, only:  r8
      implicit none
      real(r8) :: fn_Keq
! subr. arguments
      real(r8) :: Keq_298, a, b, T
! local variables
      real(r8) :: tt


        tt = 298.15/T
        fn_Keq = Keq_298*exp(a*(tt-1.)+b*(1.+log(tt)-tt))

      return
      end function fn_Keq
!----------------------------------------------------------





!----------------------------------------------------------
      function fn_Po(Po_298, DH, T)	! TOUCH
      use module_data_mosaic_kind, only:  r8
      implicit none
      real(r8) :: fn_Po
! subr. arguments
      real(r8) :: Po_298, DH, T
! local variables

        fn_Po = Po_298*exp(-(DH/8.314e-3)*(1./T - 3.354016435e-3))

      return
      end function fn_Po
!----------------------------------------------------------





!----------------------------------------------------------
      function drh_mutual(j_index)		! TOUCH
      use module_data_mosaic_aero

      implicit none

      real(r8) :: drh_mutual
! subr. arguments
      integer j_index
! local variables
      integer j


      j = j_index

      if(j_index .eq. 7 .or. j_index .eq. 8 .or.   &
        (j_index.ge. 34 .and. j_index .le. 51))then

        drh_mutual = 10.0  ! cano3 or cacl2 containing mixtures

      else

        drh_mutual =  d_mdrh(j,1) + T_K*   &
                     (d_mdrh(j,2) + T_K*   &
                     (d_mdrh(j,3) + T_K*   &
                      d_mdrh(j,4) )) + 1.0

      endif


      return
      end function drh_mutual
!----------------------------------------------------------






!----------------------------------------------------------
! ZSR method at 60% RH
!
      function aerosol_water_up(ibin) ! kg (water)/m^3 (air)
      use module_data_mosaic_aero

      implicit none

      real(r8) :: aerosol_water_up
! subr. arguments
      integer ibin
! local variables
      integer jp, je
      real(r8) :: dum_i, dum_o, kappa_org
! function
      real(r8) :: bin_molality_60


      jp = jtotal
      dum_i = 0.0

      do je = 1, (nsalt+4)	! include hno3 and hcl in water calculation
        dum_i = dum_i + (electrolyte(je,jp,ibin)/bin_molality_60(je))*1.e-9
      enddo


      ! water associated with organics
      kappa_org = 0.1
      dum_o = (aer(ic7_a,jtotal,ibin)*mw_aer_mac(ic7_a)/dens_aer_mac(ic7_a) + &
               aer(ic3_a,jtotal,ibin)*mw_aer_mac(ic3_a)/dens_aer_mac(ic3_a) + &
               aer(ic4_a,jtotal,ibin)*mw_aer_mac(ic4_a)/dens_aer_mac(ic4_a) + &
	       aer(ic5_a,jtotal,ibin)*mw_aer_mac(ic5_a)/dens_aer_mac(ic5_a) + &
	       aer(ic6_a,jtotal,ibin)*mw_aer_mac(ic6_a)/dens_aer_mac(ic6_a) + &
	       aer(ioc_a,jtotal,ibin)*mw_aer_mac(ioc_a)/dens_aer_mac(ioc_a))*1.e-9

      dum_o = dum_o * 1.0e-3 * kappa_org * 0.6/(1.0 - 0.6)	! kg(water)/m3(air)

      aerosol_water_up = dum_i + dum_o

      return
      end function aerosol_water_up
!----------------------------------------------------------






!----------------------------------------------------------
! ZSR method
      function aerosol_water(jp,ibin) ! kg (water)/m^3 (air)
      use module_data_mosaic_aero

      implicit none

      real(r8) :: aerosol_water
! subr. arguments
      integer jp, ibin
! local variables
      integer je
      real(r8) :: dum_i, dum_o, kappa_org
! function
      real(r8) :: bin_molality



!      dum = 0.0
!      do je = 1, (nsalt+4)	! include hno3 and hcl in water calculation
!        dum = dum + electrolyte(je,jp,ibin)/molality0(je)		EFFI
!      enddo

! water associated with inorganics
      dum_i = (electrolyte(1,jp,ibin)/molality0(1) + &
               electrolyte(2,jp,ibin)/molality0(2) + &
               electrolyte(3,jp,ibin)/molality0(3) + &
               electrolyte(4,jp,ibin)/molality0(4) + &
               electrolyte(5,jp,ibin)/molality0(5) + &
               electrolyte(6,jp,ibin)/molality0(6) + &
               electrolyte(7,jp,ibin)/molality0(7) + &
               electrolyte(8,jp,ibin)/molality0(8) + &
               electrolyte(9,jp,ibin)/molality0(9) + &
               electrolyte(10,jp,ibin)/molality0(10) + &
               electrolyte(11,jp,ibin)/molality0(11) + &
               electrolyte(12,jp,ibin)/molality0(12) + &
               electrolyte(13,jp,ibin)/molality0(13) + &
               electrolyte(14,jp,ibin)/molality0(14) + &
               electrolyte(15,jp,ibin)/molality0(15) + &
               electrolyte(16,jp,ibin)/molality0(16) + &
               electrolyte(17,jp,ibin)/molality0(17) + &
               electrolyte(18,jp,ibin)/molality0(18) + &
               electrolyte(19,jp,ibin)/molality0(19))*1.e-9	! kg(water)/m3(air)

! water associated with organics
      kappa_org = 0.1
      dum_o = (aer(ic7_a,jtotal,ibin)*mw_aer_mac(ic7_a)/dens_aer_mac(ic7_a) + &
      	       aer(ic3_a,jtotal,ibin)*mw_aer_mac(ic3_a)/dens_aer_mac(ic3_a) + &
               aer(ic4_a,jtotal,ibin)*mw_aer_mac(ic4_a)/dens_aer_mac(ic4_a) + &
	       aer(ic5_a,jtotal,ibin)*mw_aer_mac(ic5_a)/dens_aer_mac(ic5_a) + &
	       aer(ic6_a,jtotal,ibin)*mw_aer_mac(ic6_a)/dens_aer_mac(ic6_a) + &
	       aer(ioc_a,jtotal,ibin)*mw_aer_mac(ioc_a)/dens_aer_mac(ioc_a))*1.e-9

      dum_o = dum_o * 1.0e-3 * kappa_org * ah2o/(1.0 - ah2o)	! kg(water)/m3(air)
    ! edit wkc: The following line must be removed! This is a testing line
    dum_i = 1.0e-20

      aerosol_water = dum_i + dum_o

!      write(6,*)dum_i, dum_o

      if(aerosol_water .le. 0.0)then
        write(6,*)'iclm  jclm  ibin  jp = ',   &
                   iclm_aer, jclm_aer, ibin, jp
        write(6,*)'aH2O, water = ', aH2O, aerosol_water
        write(6,*)'dry mass = ', mass_dry_a(ibin)
        write(6,*)'soluble mass = ', mass_soluble_a(ibin)
        write(6,*)'number = ', num_a(ibin)
        do je = 1, nsoluble
          write(6,44)ename(je), electrolyte(je,jp,ibin)
        enddo
        write(6,*)'Error in water calculation'
        write(6,*)'ibin = ', ibin
        write(6,*)'water content cannot be negative or zero'
        write(6,*)'setting jaerosolstate to all_solid'

!        call print_input

        jaerosolstate(ibin) = all_solid
        jphase(ibin)    = jsolid
        jhyst_leg(ibin) = jhyst_lo

      endif

44    format(a7, 2x, e11.3)


      return
      end function aerosol_water
!----------------------------------------------------------





!----------------------------------------------------------
      function bin_molality(je,ibin)
      use module_data_mosaic_aero

      implicit none

      real(r8) :: bin_molality
! subr. arguments
      integer je, ibin
! local variables
      real(r8) :: aw, xm


      aw = max(aH2O_a(ibin), aw_min(je))
      aw = min(aw, 0.999999d0)


      if(aw .lt. 0.97)then

        xm =     a_zsr(1,je) +   &
             aw*(a_zsr(2,je) +   &
             aw*(a_zsr(3,je) +   &
             aw*(a_zsr(4,je) +   &
             aw*(a_zsr(5,je) +   &
             aw* a_zsr(6,je) ))))

        bin_molality = 55.509*xm/(1. - xm)

      else

        bin_molality = -b_zsr(je)*log(aw)

      endif


      return
      end function bin_molality
!----------------------------------------------------------





!----------------------------------------------------------
      function bin_molality_60(je)		! TOUCH
      use module_data_mosaic_aero

      implicit none

      real(r8) :: bin_molality_60
! subr. arguments
      integer je
! local variables
      real(r8) :: aw, xm


      aw = 0.6

        xm =     a_zsr(1,je) +   &
             aw*(a_zsr(2,je) +   &
             aw*(a_zsr(3,je) +   &
             aw*(a_zsr(4,je) +   &
             aw*(a_zsr(5,je) +   &
             aw* a_zsr(6,je) ))))

      bin_molality_60 = 55.509*xm/(1. - xm)

      return
      end function bin_molality_60
!----------------------------------------------------------


!----------------------------------------------------------
! 6/25/2008 - start
      function molality_0(je,aw)
      use module_data_mosaic_aero

      implicit none

      real(r8) :: molality_0
! subr. arguments
      integer je
! local variables
      real(r8) :: aw, xm


      aw = max(aw, aw_min(je))
      aw = min(aw, 0.999999d0)


      if(aw .lt. 0.97)then

        xm =     a_zsr(1,je) +   &
             aw*(a_zsr(2,je) +   &
             aw*(a_zsr(3,je) +   &
             aw*(a_zsr(4,je) +   &
             aw*(a_zsr(5,je) +   &
             aw* a_zsr(6,je) ))))

        molality_0 = 55.509*xm/(1. - xm)

      else

        molality_0 = -b_zsr(je)*log(aw)

      endif


      return
      end function molality_0
! 6/25/2008 - end
!----------------------------------------------------------


!----------------------------------------------------------
      function fnlog_gamZ(jA,jE)	! jA in jE
      use module_data_mosaic_aero

      implicit none

      real(r8) :: fnlog_gamZ
! subr. arguments
      integer jA, jE
! local variables
      real(r8) :: aw


      aw = max(aH2O, aw_min(jE))

      fnlog_gamZ = b_mtem(1,jA,jE) + aw*   &
                  (b_mtem(2,jA,jE) + aw*   &
                  (b_mtem(3,jA,jE) + aw*   &
                  (b_mtem(4,jA,jE) + aw*   &
                  (b_mtem(5,jA,jE) + aw*   &
                   b_mtem(6,jA,jE) ))))

      return
      end function fnlog_gamZ
!----------------------------------------------------------




!----------------------------------------------------------
      function mean_molecular_speed(T, MW)	! in cm/s
      use module_data_mosaic_kind, only:  r8
      implicit none
      real(r8) :: mean_molecular_speed
! subr. arguments
      real(r8) :: T, MW	! T(K)

        mean_molecular_speed = 1.455e4 * sqrt(T/MW)

      return
      end function mean_molecular_speed
!----------------------------------------------------------




!----------------------------------------------------------
      function gas_diffusivity(T, P, MW, Vm)	! in cm^2/s
      use module_data_mosaic_kind, only:  r8
      use module_data_mosaic_main, only:  third
      implicit none
      real(r8) :: gas_diffusivity
! subr. arguments
      real(r8) :: MW, Vm, T, P	! T(K), P(atm)


      gas_diffusivity = (1.0e-3 * T**1.75 * sqrt(1./MW + 0.035))/   &
                             (P * (Vm**third + 2.7189)**2)


      return
      end function gas_diffusivity
!----------------------------------------------------------




!----------------------------------------------------------
      function fuchs_sutugin(rkn,a)
      use module_data_mosaic_kind, only:  r8
      implicit none
      real(r8) :: fuchs_sutugin
! subr. arguments
      real(r8) :: rkn, a
! local variables
      real(r8) :: rnum, denom


      rnum  = 0.75*a*(1. + rkn)
      denom = rkn**2 + rkn + 0.283*rkn*a + 0.75*a
      fuchs_sutugin = rnum/denom

      return
      end function fuchs_sutugin
!----------------------------------------------------------





!----------------------------------------------------------
! solution to x^3 + px^2 + qx + r = 0
!
      function cubic( psngl, qsngl, rsngl )
      use module_data_mosaic_kind, only:  r8
      implicit none
      real(r8) :: cubic
! subr arguments
      real(r8) :: psngl, qsngl, rsngl
! local variables
      real(r8) :: p, q, r, A, B, D, M, N, third, y
      real(r8) :: k, phi, thesign, x(3), duma
      integer icase, kk

      third = 1.d0/3.d0

      q = (qsngl)
      p = (psngl)
      r = (rsngl)

      A = (1.d0/3.d0)*((3.d0*q) - (p*p))
      B = (1.d0/27.d0)*((2.d0*p*p*p) - (9.d0*p*q) + (27.d0*r))

      D = ( ((A*A*A)/27.d0) + ((B*B)/4.d0) )

      if(D .gt. 0.)then	!	=> 1 real and 2 complex roots
        icase = 1
      elseif(D .eq. 0.)then !	=> 3 real roots, atleast 2 identical
        icase = 2
      else	! D < 0		=> 3 distinct real roots
        icase = 3
      endif


      goto (1,2,3), icase

! case 1: D > 0
1     thesign = 1.
      if(B .gt. 0.)then
        B = -B
        thesign = -1.
      endif

      M = thesign*((-B/2.d0) + (sqrt(D)))**(third)
      N = thesign*((-B/2.d0) - (sqrt(D)))**(third)

      cubic = ( (M) + (N) - (p/3.d0) )
      return

! case 2: D = 0
2     thesign = 1.
      if(B .gt. 0.)then
        B = -B
        thesign = -1.
      endif

      M = thesign*(-B/2.d0)**third
      N = M

      x(1) = ( (M) + (N) - (p/3.d0) )
      x(2) = ( (-M/2.d0) + (-N/2.d0) - (p/3.d0) )
      x(2) = ( (-M/2.d0) + (-N/2.d0) - (p/3.d0) )

      cubic = 0.
      do kk = 1, 3
        if(x(kk).gt.cubic) cubic = x(kk)
      enddo
      return

! case 3: D < 0
3     if(B.gt.0.)then
        thesign = -1.
      elseif(B.lt.0.)then
        thesign = 1.
      endif

! rce 18-nov-2004 -- make sure that acos argument is between +/-1.0
!     phi = acos(thesign*sqrt( (B*B/4.d0)/(-A*A*A/27.d0) ))	! radians
      duma = thesign*sqrt( (B*B/4.d0)/(-A*A*A/27.d0) )
      duma = min( duma, +1.0d0 )
      duma = max( duma, -1.0d0 )
      phi  = acos( duma )	! radians


      cubic = 0.
      do kk = 1, 3
        k = kk-1
        y = 2.*Sqrt(-A/3.)*cos(phi + 120.*k*0.017453293)
        x(kk) = ((y) - (p/3.d0))
        if(x(kk).gt.cubic) cubic = x(kk)
      enddo
      return

      end function cubic
!----------------------------------------------------------




!----------------------------------------------------------
      function quadratic(a,b,c)
      use module_data_mosaic_kind, only:  r8
      implicit none
      real(r8) :: quadratic
! subr. arguments
      real(r8) :: a, b, c
! local variables
      real(r8) :: x, dum, quad1, quad2


        if(b .ne. 0.0)then
        x = 4.*(a/b)*(c/b)
        else
        x = 1.e+6
        endif

        if(abs(x) .lt. 1.e-6)then
          dum = ( (0.5*x) +   &
                  (0.125*x**2) +   &
                  (0.0625*x**3) )

          quadratic = (-0.5*b/a)*dum

          if(quadratic .lt. 0.)then
            quadratic = -b/a - quadratic
          endif

        else
          quad1 = ((-b)+sqrt((b*b)-(4.*a*c)))/   &
                                     (2.*a)
          quad2 = ((-b)-sqrt((b*b)-(4.*a*c)))/   &
                                     (2.*a)

          quadratic = max(quad1, quad2)
        endif

      return
      end function quadratic
!----------------------------------------------------------



!----------------------------------------------------------
! currently not used
!
! two roots of a quadratic equation
!
      subroutine quadratix(a,b,c, qx1,qx2)
      use module_data_mosaic_kind, only:  r8
      implicit none
! subr. arguments
      real(r8) :: a, b, c, qx1, qx2
! local variables
      real(r8) :: x, dum


      if(b .ne. 0.0)then
        x = 4.*(a/b)*(c/b)
        else
        x = 1.e+6
      endif

      if(abs(x) .lt. 1.e-6)then
        dum = ( (0.5*x) +   &
                (0.125*x**2) +   &
                (0.0625*x**3) )

        qx1 = (-0.5*b/a)*dum
        qx2 = -b/a - qx1

      else

        qx1 = ((-b)+sqrt((b*b)-(4.*a*c)))/   &
                                     (2.*a)
        qx2 = ((-b)-sqrt((b*b)-(4.*a*c)))/   &
                                     (2.*a)

      endif

      return
      end subroutine quadratix


!=====================================================================

















!***********************************************************************
! computes aerosol optical properties
!
! author: Rahul A. Zaveri
! update: jan 2005
!-----------------------------------------------------------------------
      subroutine aerosol_optical_properties(iclm,jclm)
      use module_data_mosaic_main
      use module_data_mosaic_gas
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer iclm, jclm

! local variables
      integer iaer, ibin, je, k

      real(r8) :: XT
      real(r8) :: number_bin(nbin_a,kmaxd)
      real(r8) :: radius_wet(nbin_a, kmaxd)

      complex refindx(nbin_a, kmaxd)

	integer kmax 			! maximum number of levels in k direction
	parameter (kmax=100)
        integer nspint  		! Num of spectral intervals across solar spectrum
        parameter ( nspint = 4 )	! for FAST-J

	real(r8) :: wavmid, sizeaer, extaer, waer, gaer,   &
          l2, l3, l4, l5, l6, l7, tauaer1
        common /mie/  wavmid(nspint), sizeaer(nspint,kmax),   &
          extaer(nspint,kmax),   &
          waer(nspint,kmax), gaer(nspint,kmax), l2(nspint,kmax),   &
          l3(nspint,kmax), l4(nspint,kmax), l5(nspint,kmax),   &
          l6(nspint,kmax), l7(nspint,kmax), tauaer1(nspint,kmax)




      call load_mosaic_parameters

!      iclm_aer = iclm
!      jclm_aer = jclm

!      do 100 k = 1, ktot

!        cair_mol_m3 = cairclm(k)*1.e6		! cairclm(k) is in mol/cc
!        cair_mol_cc = cairclm(k)

        conv1a = cair_mol_m3*1.e9		! converts q/mol(air) to nq/m^3 (q = mol or g)
        conv1b = 1./conv1a			! converts nq/m^3 to q/mol(air)
        conv2a = cair_mol_m3*18.*1.e-3		! converts mol(h2o)/mol(air) to kg(h2o)/m^3(air)
        conv2b = 1./conv2a			! converts kg(h2o)/m^3(air) to mol(h2o)/mol(air)


! initialize to zero
        do ibin = 1, nbin_a
          do iaer = 1, naer
            aer(iaer,jtotal,ibin)  = 0.0
          enddo

          do je = 1, nelectrolyte
            electrolyte(je,jtotal,ibin)  = 0.0
          enddo

          jaerosolstate(ibin) = -1	! initialize to default value

        enddo

! map rclm(k,l) to aer(i,jtotal,ibin)
        call map_mosaic_species_BOX(0)

        do 90 ibin = 1, nbin_a

          call check_aerosol_mass(ibin)
          if(jaerosolstate(ibin) .eq. no_aerosol)goto 90	! ignore this bin
          call conform_electrolytes(jtotal,ibin,XT) 		! conforms aer(jtotal) to a valid aerosol
          call check_aerosol_mass(ibin) 			! check mass again after conform_electrolytes
          if(jaerosolstate(ibin) .eq. no_aerosol)goto 90	! ignore this bin
          call conform_aerosol_number(ibin)   			! adjusts number conc so that it conforms with bin mass and diameter
          call calc_dry_n_wet_aerosol_props(ibin)		! calc Dp_wet, ref index

!          refindx(ibin,k)    = ri_avg_a(ibin)			! vol avg ref index
!          radius_wet(ibin,k) = Dp_wet_a(ibin)/2.0		! wet radius (cm)
!          number_bin(ibin,k) = num_a(ibin)			! #/cc air

90      continue

! 100   continue	! k levels


!      call p2cputim( 14 )
!      call mieaer(nbin_a,ktot,number_bin,radius_wet,refindx)
!      call p2cputim( 15 )


      return
      end subroutine aerosol_optical_properties
















