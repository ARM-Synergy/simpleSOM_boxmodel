! 1/1/2013 raz set idiagbb_coag = 0 (turn off printing diagnostics)

	module module_sect_iface


	use module_data_mosaic_kind, only:  r8


	implicit none


	integer, parameter :: iunits_flagaa = 1
	integer, parameter :: iunits_flagbb = 1
	integer, parameter :: iunits_diagout_flagaa = 1


	contains


!-----------------------------------------------------------------------
	subroutine sectional_interface_1( dtchem, cnn_sv1 )
	use module_data_mosaic_main
	use module_data_mosaic_aero
	use module_data_mosaic_asect, only:  maxd_asize, maxd_atype, &
		nsize_aer, ntype_aer
	use module_movesect, only:  move_sections_x3
	use module_movesect_3d, only:  move_sect_3d_x1
	use module_coag,     only:  mosaic_coag_1box
	use module_coag_3d,  only:  mosaic_coag_3d_1box
	use module_newnuc,   only:  mosaic_newnuc_1box

	implicit none

!   subr arguments
	real(r8) :: dtchem
	real(r8) :: cnn_sv1(ntot_used)

!   local variables
	integer  :: idiagbb_coag, istat_coag
	integer  :: idiagbb_newnuc, istat_newnuc, itype_newnuc
	integer  :: idiagcc
	integer  :: jsv
	integer  :: method_coag, method_movesect

	real(r8) :: dens_nh4so4a_newnuc
	real(r8) :: fact_apmassmr, fact_apnumbmr,   &
	            fact_apdens, fact_apdiam
	real(r8) :: rh_box, rhoair_g_cc
	real(r8) :: tmpa

	real(r8) :: cnn_sv2(ntot_used)
	real(r8) :: rbox(ntot_used)
	real(r8) :: rbox0(ntot_used)
	real(r8) :: rbox_sv(ntot_used,4)

	real(r8) :: drydens_pregrow(maxd_asize,maxd_atype)
	real(r8) :: drydens_aftgrow(maxd_asize,maxd_atype)
	real(r8) :: drymass_pregrow(maxd_asize,maxd_atype)
	real(r8) :: drymass_aftgrow(maxd_asize,maxd_atype)
! drydens_xxxgrow = dry density (g/cm3) of bin
! drymass_xxxgrow = all-component dry-mass mixing ratio (g-AP/m3-air) of bin
!    pregrow (aftgrow) means before (after) the last call to the
!    gas-aerosol mass-transfer solver
! these arrays are used by the movesect routine to calculate transfer 
!    between bins due to particle growth/shrinkage

	real(r8) :: adrydens_box(maxd_asize,maxd_atype)
	real(r8) :: adrydpav_box(maxd_asize,maxd_atype)
	real(r8) :: adryqmas_box(maxd_asize,maxd_atype)
	real(r8) :: awetdens_box(maxd_asize,maxd_atype)
	real(r8) :: awetdpav_box(maxd_asize,maxd_atype)
! adrydens_box = current dry density (g/cm3) of bin
! adrydpav_box = current mean dry diameter (cm) of bin
! adryqmas_box = current all-component dry-mass mixing ratio (g-AP/m3-air) of bin
! awetdens_box = current wet/ambient density (g/cm3) of bin
! awetdpav_box = current mean wet/ambient diameter (cm) of bin
!
	real(r8) :: adrydens_box_sv(maxd_asize,maxd_atype,4)
	real(r8) :: adrydpav_box_sv(maxd_asize,maxd_atype,4)
	real(r8) :: adryqmas_box_sv(maxd_asize,maxd_atype,4)
	real(r8) :: awetdens_box_sv(maxd_asize,maxd_atype,4)
	real(r8) :: awetdpav_box_sv(maxd_asize,maxd_atype,4)


	cnn_sv2(1:ntot_used) = cnn(1:ntot_used)


! the movesect routine puts "initial" values into the axxxxxxx_box arrays
!    so need to have movesect turned (mmovesect_flag1 > 0) 
!    for any of the sectional routines to function properly
	if (mmovesect_flag1 <= 0) then
	    write(*,*)   &
	    '*** skipping sectional_interface_1 -- mmovesect_flag1 <= 0'
	    return
	end if


! set conversion factors
	tmpa = cair_mol_m3*mw_air   ! air_density in (g/m^3)
	if (iunits_flagaa <= 1) then
	    fact_apmassmr = 1.0e-6_r8/tmpa  ! converts aerosol mass mrs  in rbox
	                                    ! from ug-AP/m^3-air to g-AP/g-air
	    fact_apnumbmr = 1.0_r8/tmpa     ! converts aerosol numb mrs  in rbox
	                                    ! from #/m^3-air to #/g-air
	else if (iunits_flagaa == 2) then
	    fact_apmassmr = 1.0_r8/tmpa     ! converts aerosol mass mrs  in rbox
	                                    ! from g-AP/m^3-air to g-AP/g-air
	    fact_apnumbmr = 1.0_r8/tmpa     ! converts aerosol numb mrs  in rbox
	                                    ! from #/m^3-air to #/g-air
	else
	    fact_apmassmr = 1.0_r8          ! aerosol mass mrs already are g-AP/g-air
	    fact_apnumbmr = 1.0_r8          ! aerosol numb mrs already are #/g-air
	end if
	fact_apdens = 1.0_r8             ! converts mosaic aerosol densities from g/cm^3 to g/cm^3
	fact_apdiam = 1.0_r8             ! converts mosaic aerosol diameters from cm to cm

! convert densities and sizes for testing purposes
	dens_nh4so4a_newnuc = dens_aer_mac(iso4_a)
	if (iunits_flagbb == 2) call sect_iface_dens_size_testbb( 1, &
	    fact_apdens, fact_apdiam, dens_nh4so4a_newnuc )


! map gases only from cnn_sv1 to rbox0
!   cnn_sv1 (and rbox0) values are before the mosaic-astem mass-xfer calcs, 
!   and are used in the newnuc routine
!
!   first load cnn with cnn_sv1 values
	cnn(1:ntot_used) = cnn_sv1(1:ntot_used)
!   then map gases from cnn to rbox0
	call map_mosaic_cnn_rbox( 0, rbox0, &
	    drydens_aftgrow, drydens_pregrow,   &
	    drymass_aftgrow, drymass_pregrow,   &
	    adrydens_box, awetdens_box, adrydpav_box, awetdpav_box,   &
	    adryqmas_box,   &
	    fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam )
!   then restore cnn from cnn_sv2
	cnn(1:ntot_used) = cnn_sv2(1:ntot_used)


! map from cnn to rbox
	call map_mosaic_cnn_rbox( 1, rbox, &
	    drydens_aftgrow, drydens_pregrow,   &
	    drymass_aftgrow, drymass_pregrow,   &
	    adrydens_box, awetdens_box, adrydpav_box, awetdpav_box,   &
	    adryqmas_box,   &
	    fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam )
	jsv = 1
	rbox_sv(:,jsv) = rbox(:)
	adrydens_box_sv(:,:,jsv) = adrydens_box(:,:)
	awetdens_box_sv(:,:,jsv) = awetdens_box(:,:)
	adrydpav_box_sv(:,:,jsv) = adrydpav_box(:,:)
	awetdpav_box_sv(:,:,jsv) = awetdpav_box(:,:)
	adryqmas_box_sv(:,:,jsv) = adryqmas_box(:,:)


! move_sections transfers particles between sections
! following condensational growth / evaporative shrinking
!
	method_movesect = mod( max(0,mmovesect_flag1), 100 )
	if (method_movesect < 50) then
!	  call move_sections_x3( iflag, iclm, jclm, k, m, rbox,   &
!	    drydens_aftgrow, drydens_pregrow,   &
!	    drymass_aftgrow, drymass_pregrow,   &
!	    adrydens_tmp, awetdens_tmp, adrydbar_tmp, awetdbar_tmp,   &
!	    adryqmas_tmp, adryqvol_tmp )
	  call move_sections_x3( 1, 1, 1, 1, 1, rbox,   &
	    fact_apmassmr, fact_apnumbmr,   &
	    fact_apdens, fact_apdiam,   &
	    drydens_aftgrow, drydens_pregrow,   &
	    drymass_aftgrow, drymass_pregrow,   &
	    adrydens_box, awetdens_box, adrydpav_box, awetdpav_box,   &
	    adryqmas_box )
	else
	  call move_sect_3d_x1( 1, 1, 1, 1, 1, rbox,   &
	    fact_apmassmr, fact_apnumbmr,   &
	    fact_apdens, fact_apdiam,   &
	    drydens_aftgrow, drydens_pregrow,   &
	    drymass_aftgrow, drymass_pregrow,   &
	    adrydens_box, awetdens_box, adrydpav_box, awetdpav_box,   &
	    adryqmas_box )
	end if

	jsv = 2
	rbox_sv(:,jsv) = rbox(:)
	adrydens_box_sv(:,:,jsv) = adrydens_box(:,:)
	awetdens_box_sv(:,:,jsv) = awetdens_box(:,:)
	adrydpav_box_sv(:,:,jsv) = adrydpav_box(:,:)
	awetdpav_box_sv(:,:,jsv) = awetdpav_box(:,:)
	adryqmas_box_sv(:,:,jsv) = adryqmas_box(:,:)


! do new particle nucleation
!
	if (mnewnuc_flag1 > 0) then
	    idiagbb_newnuc = +200
	    itype_newnuc = 1
	    rh_box = rh*0.01
	    call mosaic_newnuc_1box(   &
		istat_newnuc, idiagbb_newnuc,   &
		it_mosaic, itype_newnuc, dtchem,   &
		te, pr_atm, cair_mol_cc, rh_box,   &
		dens_nh4so4a_newnuc, mw_aer_mac(iso4_a), mw_aer_mac(inh4_a),   &
		mw_comp_a(jh2o), mw_air,  &
		fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam,   &
		rbox0, rbox,   &
		adrydens_box, awetdens_box,   &
		adrydpav_box, awetdpav_box, adryqmas_box )
	end if

	jsv = 3
	rbox_sv(:,jsv) = rbox(:)
	adrydens_box_sv(:,:,jsv) = adrydens_box(:,:)
	awetdens_box_sv(:,:,jsv) = awetdens_box(:,:)
	adrydpav_box_sv(:,:,jsv) = adrydpav_box(:,:)
	awetdpav_box_sv(:,:,jsv) = awetdpav_box(:,:)
	adryqmas_box_sv(:,:,jsv) = adryqmas_box(:,:)


! do particle coagulation
!
	if (mcoag_flag1 > 0) then
	    idiagbb_coag = 0 ! +200
	    rhoair_g_cc = cair_mol_cc*mw_air
	    method_coag = mod( max(0,mcoag_flag1), 100 )
	    if (method_coag < 50) then
		call mosaic_coag_1box( istat_coag,   &
		    idiagbb_coag, it_mosaic,   &
		    dtchem, te, pr_atm, rhoair_g_cc, rbox,   &
		    fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam,   &
		    adrydens_box, awetdens_box,   &
		    adrydpav_box, awetdpav_box, adryqmas_box )
	    else
		call mosaic_coag_3d_1box( istat_coag,   &
		    idiagbb_coag, it_mosaic,   &
		    dtchem, te, pr_atm, rhoair_g_cc, rbox,   &
		    fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam,   &
		    adrydens_box, awetdens_box,   &
		    adrydpav_box, awetdpav_box, adryqmas_box )
	    end if
	end if

	jsv = 4
	rbox_sv(:,jsv) = rbox(:)
	adrydens_box_sv(:,:,jsv) = adrydens_box(:,:)
	awetdens_box_sv(:,:,jsv) = awetdens_box(:,:)
	adrydpav_box_sv(:,:,jsv) = adrydpav_box(:,:)
	awetdpav_box_sv(:,:,jsv) = awetdpav_box(:,:)
	adryqmas_box_sv(:,:,jsv) = adryqmas_box(:,:)


! map from rbox back to cnn
	call map_mosaic_cnn_rbox( 2, rbox, &
	    drydens_aftgrow, drydens_pregrow,   &
	    drymass_aftgrow, drymass_pregrow,   &
	    adrydens_box, awetdens_box, adrydpav_box, awetdpav_box,   &
	    adryqmas_box,   &
	    fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam )


	idiagcc = 1
	if (nsize_aer(1)*ntype_aer > 10000) idiagcc = 0
	if (nsize_aer(1)*ntype_aer >   100) idiagcc = 0
! this was needed for early testing -- should not need it for large runs
	if (idiagcc > 0) &
	  call mbox_sectional_diagnostics( cnn_sv1, cnn_sv2, rbox_sv, &
	    drydens_aftgrow, drydens_pregrow,   &
	    drymass_aftgrow, drymass_pregrow,   &
	    adrydens_box_sv, awetdens_box_sv, adrydpav_box_sv,   &
	    awetdpav_box_sv, adryqmas_box_sv,   &
	    fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam )


! unconvert densities and sizes
	if (iunits_flagbb == 2) call sect_iface_dens_size_testbb( 2, &
	    fact_apdens, fact_apdiam, dens_nh4so4a_newnuc )


	return
	end subroutine sectional_interface_1


!-----------------------------------------------------------------------
	subroutine map_mosaic_cnn_rbox( iflagaa, rbox, &
	    drydens_aftgrow, drydens_pregrow,   &
	    drymass_aftgrow, drymass_pregrow,   &
	    adrydens_box, awetdens_box, adrydpav_box,   &
	    awetdpav_box, adryqmas_box,   &
	    fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam )
!
!   maps between cnn array and rbox array
!
	use module_data_mosaic_main
	use module_data_mosaic_aero
	use module_data_mosaic_asect
	use module_movesect, only:  test_move_sections

	implicit none

! subr parameters
	integer,  intent(in)    :: iflagaa
	real(r8), intent(inout) :: rbox(ntot_used)

	real(r8), intent(inout) :: drydens_aftgrow(maxd_asize,maxd_atype)
	real(r8), intent(inout) :: drydens_pregrow(maxd_asize,maxd_atype)
	real(r8), intent(inout) :: drymass_aftgrow(maxd_asize,maxd_atype)
	real(r8), intent(inout) :: drymass_pregrow(maxd_asize,maxd_atype)

	real(r8), intent(inout) :: adrydens_box(maxd_asize,maxd_atype)
	real(r8), intent(inout) :: awetdens_box(maxd_asize,maxd_atype)
	real(r8), intent(inout) :: adrydpav_box(maxd_asize,maxd_atype)
	real(r8), intent(inout) :: awetdpav_box(maxd_asize,maxd_atype)
	real(r8), intent(inout) :: adryqmas_box(maxd_asize,maxd_atype)

	real(r8), intent(in) :: fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam


! local variables
	integer :: ibin, iphase, isize, itype, jhyst_tmp
	integer :: l, ll, lunaa, mtmp, noffset
	real(r8) :: conv_diam, conv_drym, conv_gas, conv_numb, conv_watr
	real(r8) :: tmpa, tmph, tmpj, tmpr


! on first time step, call test_move_sections
! note that test_move_sections only executes when mmovesect_flag1 = 80nn
	if ((it_mosaic == 1) .and. (iflagaa == 0)) then
	    call test_move_sections( 1, 1, 1, 1, 1 )
	end if


	if ((iflagaa < 0) .or. (iflagaa > 2)) return


	lunaa = lun_sect_172

	tmpa = cair_mol_m3*mw_air   ! air_density in (g/m^3)
	conv_gas  = 1.0_r8/(avogad*cair_mol_cc)     ! convert molecules/cm^3 to mol/mol-air
	conv_numb = 1.0e6_r8 /(fact_apnumbmr*tmpa)  ! convert #/cm^3 to 
	                                            ! (1)  #/m^3 OR (2) #/m^3 OR (3) #/g-air
	conv_drym = 1.0e-6_r8/(fact_apmassmr*tmpa)  ! convert umol-xxx/m^3 to 
	                                            ! (1) ug/m^3 OR (2) g/m^3 OR (3) g/g-air
	                                            ! (also multiply by molecular weight in g)
	conv_watr = 1.0e3_r8 /(fact_apmassmr*tmpa)  ! convert kg-water/m^3 to 
	                                            ! (1) ug/m^3 OR (2) g/m^3 OR (3) g/g-air
	conv_diam = 1.0e-4_r8/fact_apdiam           ! convert um to cm (or m)


	if (iflagaa == 2) goto 20000


!
! iflagaa = 1 --> map from cnn to rbox
! iflagaa = 0 --> map gases only from cnn to rbox
!
	do l = 1, ngas_max
	    rbox(l) = cnn(l)*conv_gas
	end do
	if (iflagaa == 0) return

	adrydens_box(:,:) = 0.0
	awetdens_box(:,:) = 0.0
	adrydpav_box(:,:) = 0.0
	awetdpav_box(:,:) = 0.0
	adryqmas_box(:,:) = 0.0

	do ibin = 1, nbin_a
	    noffset = ngas_max + (ibin-1)*naer_tot
	    rbox(knum_a   +noffset) = cnn(knum_a   +noffset)*conv_numb
	    rbox(kdpdry_a +noffset) = cnn(kdpdry_a +noffset)*conv_diam
	    rbox(ksigmag_a+noffset) = cnn(ksigmag_a+noffset)
	    if (mhyst_method == mhyst_uporlo_waterhyst) then
		! units same as for water
		rbox(kjhyst_a +noffset) = water_a_hyst(ibin)*conv_watr
	    else
		! do not apply units conversion here as this will
		! be modified (below)
		rbox(kjhyst_a +noffset) = cnn(kjhyst_a +noffset)
	    end if
	    rbox(kwater_a +noffset) = cnn(kwater_a +noffset)*conv_watr
	    do ll = 1, naer
		l = noffset + (naer_tot-naer) + ll
		rbox(l) = cnn(l)*(conv_drym*mw_aer_mac(ll))
	    end do

	    isize = isize_of_ibin(ibin)
	    itype = itype_of_ibin(ibin)
	    if (jaerosolstate_bgn(ibin) /= no_aerosol) then
		drydens_pregrow(isize,itype) = dens_dry_a_bgn(ibin)
		drymass_pregrow(isize,itype) = mass_dry_a_bgn(ibin)
	    else
		drydens_pregrow(isize,itype) = -1.0
		drymass_pregrow(isize,itype) =  0.0
	    end if
	    if (jaerosolstate(ibin) /= no_aerosol) then
		drydens_aftgrow(isize,itype) = dens_dry_a(ibin)
		drymass_aftgrow(isize,itype) = mass_dry_a(ibin)
	    else
		drydens_aftgrow(isize,itype) = -1.0
		drymass_aftgrow(isize,itype) =  0.0
	    end if

	    ! convert drymass_xxxgrow from (g/cm^3) to
	    ! (1) ug/m^3 OR (2) g/m^3 OR (3) g/g-air
	    tmpa = 1.0e6_r8/(cair_mol_m3*mw_air*fact_apmassmr)
	    drymass_pregrow(isize,itype) = drymass_pregrow(isize,itype) * tmpa
	    drymass_aftgrow(isize,itype) = drymass_aftgrow(isize,itype) * tmpa

	    ! convert drydens_xxxgrow from (g/cm^3) to (kg/m^3) when iunits_flagbb == 2
	    tmpa = 1.0_r8/fact_apdens
	    drydens_pregrow(isize,itype) = drydens_pregrow(isize,itype) * tmpa
	    drydens_aftgrow(isize,itype) = drydens_aftgrow(isize,itype) * tmpa

	    if (mhyst_method /= mhyst_uporlo_waterhyst) then
! here rbox(hyswptr_aer(:,:)) holds jhyst_leg rather than water_a_hyst
! convert it from [0/1 flag] to [(0/1 flag)*(hygro*volume mixing ratio)]
!    so it can be treated reasonably in movesect and coag routines
! note - the tmph (below) units are essentially the same as for mass species in rbox
		iphase = ai_phase
		tmph = 0.0
		do ll = 1, ncomp_aer(itype)
		    l = massptr_aer(ll,isize,itype,iphase)
		    tmph = tmph + max(0.0_r8,rbox(l))*   &
				(hygro_aer(ll,itype)/dens_aer(ll,itype))
		end do
		tmph = max( tmph*dens_water_aer, 1.0e-30_r8 )
		l = hyswptr_aer(isize,itype)
		rbox(l) = rbox(l)*tmph
	    end if

	end do

	return


20000	continue
!
! iflagaa = 2 --> map from rbox to cnn
!
	do l = 1, ngas_max
	    cnn(l) = rbox(l)/conv_gas
	end do

	do ibin = 1, nbin_a
	    isize = isize_of_ibin(ibin)
	    itype = itype_of_ibin(ibin)
	    noffset = ngas_max + (ibin-1)*naer_tot
	    cnn(knum_a   +noffset) = rbox(knum_a   +noffset)/conv_numb
!	    cnn(kdpdry_a +noffset) = rbox(kdpdry_a +noffset)/conv_diam
	    cnn(kdpdry_a +noffset) = adrydpav_box(isize,itype)/conv_diam
	    cnn(ksigmag_a+noffset) = rbox(ksigmag_a+noffset)
	    if (mhyst_method == mhyst_uporlo_waterhyst) then
		! units same as for water
		water_a_hyst(ibin) = rbox(kjhyst_a +noffset)/conv_watr
	    else
		cnn(kjhyst_a +noffset) = rbox(kjhyst_a +noffset)
	    end if
	    cnn(kwater_a +noffset) = rbox(kwater_a +noffset)/conv_watr
	    do ll = 1, naer
		l = noffset + (naer_tot-naer) + ll
		cnn(l) = rbox(l)/(conv_drym*mw_aer_mac(ll))
	    end do

	    if (mhyst_method /= mhyst_uporlo_waterhyst) then
! convert rbox(hyswptr_aer(:,:)) from [(0/1 flag)*(hygro*volume mixing ratio)] 
!    back to [0/1 flag]
! note - the tmph (below) units are essentially the same as for mass species in rbox
		iphase = ai_phase
		tmph = 0.0
		do ll = 1, ncomp_aer(itype)
		    l = massptr_aer(ll,isize,itype,iphase)
		    tmph = tmph + max(0.0_r8,rbox(l))*   &
			              (hygro_aer(ll,itype)/dens_aer(ll,itype))
		end do
		tmph = max( tmph*dens_water_aer, 1.0e-30_r8 )
		l = hyswptr_aer(isize,itype)
		tmpr = max( rbox(l), 0.0_r8 )
		if (tmpr > tmph*1.0e6) then   ! here tmpr/tmph > 1.0e6
		    tmpj = 1.0e6
		else
		    tmpj = tmpr/tmph
		end if
! *** at this point, tmpj can be anything between 0 and 1 (or even > 1)
!	 because movesect and coag can "mix" upper and lower curve
!	 particles from different bins
!     need to decide on a reasonable "threshold" that determines how a
!	 "mixed" particles is be classified
!	 (0.99 probably NOT a good threshold)
		if (tmpj > 0.5) then
		    jhyst_tmp = 1
		else
		    jhyst_tmp = 0
		end if
		cnn(kjhyst_a +noffset) = jhyst_tmp
	    end if

	end do

	return
	end subroutine map_mosaic_cnn_rbox


!-----------------------------------------------------------------------
	subroutine mbox_sectional_diagnostics( cnn_sv1, cnn_sv2, rbox_sv, &
	    drydens_aftgrow, drydens_pregrow,   &
	    drymass_aftgrow, drymass_pregrow,   &
	    adrydens_box_sv, awetdens_box_sv, adrydpav_box_sv,   &
	    awetdpav_box_sv, adryqmas_box_sv,   &
	    fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam )
!
!   writes test diagnostics for the combination of movesect, newnuc, coag
!
	use module_data_mosaic_main
	use module_data_mosaic_aero
	use module_data_mosaic_asect

	implicit none

! subr parameters
	real(r8), intent(in) :: cnn_sv1(ntot_used), cnn_sv2(ntot_used)
	real(r8), intent(in) :: rbox_sv(ntot_used,4)
	real(r8), intent(in) :: drydens_aftgrow(maxd_asize,maxd_atype)
	real(r8), intent(in) :: drydens_pregrow(maxd_asize,maxd_atype)
	real(r8), intent(in) :: drymass_aftgrow(maxd_asize,maxd_atype)
	real(r8), intent(in) :: drymass_pregrow(maxd_asize,maxd_atype)

	real(r8), intent(in) :: adrydens_box_sv(maxd_asize,maxd_atype,4)
	real(r8), intent(in) :: awetdens_box_sv(maxd_asize,maxd_atype,4)
	real(r8), intent(in) :: adrydpav_box_sv(maxd_asize,maxd_atype,4)
	real(r8), intent(in) :: awetdpav_box_sv(maxd_asize,maxd_atype,4)
	real(r8), intent(in) :: adryqmas_box_sv(maxd_asize,maxd_atype,4)

	real(r8), intent(in) :: fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam

! local variables
	integer  :: ibin, iphase, isize, itype
	integer  :: jsv,  l, ll, lunaa, mtmp, noffset
	real(r8) :: conv_watr
	real(r8) :: tmpa, tmpb, tmph, tmpveca(4), tmpvech(4)


	lunaa = lun_sect_172

	write(lunaa,'(//a,2i5)') 'mbox_sectional_diagnostics -- it', it_mosaic
	write(lunaa,'(a,1p,4e12.4)') 'te, pr, cair1,2', te, pr_atm, cair_mol_m3,  cair_mol_cc
	write(lunaa,'(a,1p,4e12.4)') 'time_sec,dt_sec', time_sec, dt_sec

	write(lunaa,'(a)')
	do l = 1, 4
		write(lunaa,'(i4,2x,a,1p,8e12.4)') &
		    l, species(l)(1:12), cnn_sv1(l), cnn_sv2(l), cnn(l),   &
		    rbox_sv(l,1:4)
	end do

	do ibin = 1, nbin_a
	    isize = isize_of_ibin(ibin)
	    itype = itype_of_ibin(ibin)

	    iphase = ai_phase
	    do jsv = 1, 4
		if (mhyst_method == mhyst_uporlo_waterhyst) then
		    if (iunits_diagout_flagaa > 0) then
			tmpa = cair_mol_m3*mw_air   ! air_density in (g/m^3)
			conv_watr = 1.0e3_r8 /(fact_apmassmr*tmpa)
			tmpvech(jsv) = conv_watr
		    else
			tmpvech(jsv) = 1.0
		    end if
		else
		    tmph = 0.0
		    do ll = 1, ncomp_aer(itype)
			l = massptr_aer(ll,isize,itype,iphase)
			tmph = tmph + max(0.0_r8,rbox_sv(l,jsv))*   &
		                  (hygro_aer(ll,itype)/dens_aer(ll,itype))
		    end do
		    tmph = max( tmph*dens_water_aer, 1.0e-30_r8 )
		    tmpvech(jsv) = tmph
		end if
	    end do

	    write(lunaa,'(a)')
	    write(lunaa,'(4x,2x,a,1p,2i12)') &
		'jaero_state ', jaerosolstate_bgn(ibin), jaerosolstate(ibin)
	    noffset = ngas_max + (ibin-1)*naer_tot
	    do ll = 1, naer_tot
!		if ((ll >= kmsa_a) .and. (ll <= klim2_a)) cycle
		l = noffset + ll
		if (ll == kjhyst_a) then
		    tmpveca(1:4) = rbox_sv(l,1:4)/tmpvech(1:4)
!		    write(lunaa,'(i4,2x,a,1p,8e12.4)') &
!		    -1, 'water_jhyst ', water_a_hyst(ibin)
		else
		    tmpveca(1:4) = rbox_sv(l,1:4)
		    if (ll == kdpdry_a) tmpveca(1:4) = tmpveca(1:4)*fact_apdiam
		    if (iunits_diagout_flagaa > 0) then
		    if (ll == knum_a) then
		        ! convert (1)  #/m^3 OR (2) #/m^3 OR (3) #/g-air
		        ! to #/mol-air
			tmpa = mw_air*fact_apnumbmr
			tmpveca(1:4) = tmpveca(1:4)*tmpa
		    else if (ll >= kwater_a) then
		        ! convert (1) ug/m^3 OR (2) g/m^3 OR (3) g/g-air
		        ! to mol/mol-air
			tmpb = mw_comp_a(jh2o)
			if (ll > kwater_a) tmpb = mw_aer_mac(ll-kwater_a)
			tmpa = (mw_air/tmpb)*fact_apmassmr
			tmpveca(1:4) = tmpveca(1:4)*tmpa
		    end if
		    end if ! (iunits_diagout_flagaa > 0)
		end if
		write(lunaa,'(i4,2x,a,1p,8e12.4)') &
		    l, species(l)(1:12), cnn_sv1(l), cnn_sv2(l), cnn(l), &
		    tmpveca(1:4)
	    end do

	    tmpa = fact_apdiam
	    write(lunaa,'(4x,2x,a,1p,36x,5e12.4)') &
		'adrydpav    ', &
		    adrydpav_box_sv(isize,itype,1:4)*tmpa
	    write(lunaa,'(4x,2x,a,1p,36x,5e12.4)') &
		'awetdpav    ', &
		    awetdpav_box_sv(isize,itype,1:4)*tmpa

	    tmpa = fact_apdens
	    write(lunaa,'(4x,2x,a,1p,2e12.4,12x,5e12.4)') &
		'adrydens    ', &
		    drydens_pregrow(isize,itype)*tmpa, drydens_aftgrow(isize,itype)*tmpa, &
		    adrydens_box_sv(isize,itype,1:4)*tmpa
	    write(lunaa,'(4x,2x,a,1p,36x,5e12.4)') &
		'awetdens    ', &
		    awetdens_box_sv(isize,itype,1:4)*tmpa

	    if (iunits_diagout_flagaa <= 0) then
		tmpa = 1.0_r8
	    else
	        ! convert (1) ug/m^3 OR (2) g/m^3 OR (3) g/g-air
	        ! to g/mol-air
		tmpa = mw_air*fact_apmassmr
	    end if
	    write(lunaa,'(4x,2x,a,1p,2e12.4,12x,6e12.4)') &
		'adryqmas    ', &
		tmpa*drymass_pregrow(isize,itype), tmpa*drymass_aftgrow(isize,itype), &
		tmpa*adryqmas_box_sv(isize,itype,1:4)
	end do

	if (it_mosaic == 1) then
	write(lunaa,'(a)')
	do ll = 1, 10
	    l = ntot_used + ll
	    write(lunaa,'(i4,2x,a,1p,2e12.4)') &
		l, species(l)(1:12)
	end do
	end if


	return
	end subroutine mbox_sectional_diagnostics


!-----------------------------------------------------------------------
	subroutine sect_iface_dens_size_testbb( iflagaa, &
	    fact_apdens, fact_apdiam, dens_nh4so4a_newnuc )
!
!   when iunits_flagbb = 2, densities and sizes in module_data_mosaic_asect
!	are temporarily changes from (g/cm^3) to (kg/m^3), and (cm) to (m)
!
	use module_data_mosaic_asect

	implicit none

! subr parameters
	integer,  intent(in)    :: iflagaa
	real(r8), intent(inout) :: fact_apdens, fact_apdiam
	real(r8), intent(inout) :: dens_nh4so4a_newnuc

! local variables
	real(r8) :: tmpa

	real(r8), save ::   &
	  dens_nh4so4a_newnuc_tmp,  &
	  dens_water_aer_tmp,  &
	  dens_aer_tmp( maxd_acomp, maxd_atype ),  &
	  dens_mastercomp_aer_tmp( maxd_acomp )
	real(r8), save ::   &
	  volumcut_sect_tmp( 0:maxd_asize, maxd_atype ),  &
	  volumcen_sect_tmp(   maxd_asize, maxd_atype ),  &
	  volumlo_sect_tmp(    maxd_asize, maxd_atype ),   &
	  volumhi_sect_tmp(    maxd_asize, maxd_atype ),   &
	  dcut_sect_tmp( 0:maxd_asize, maxd_atype ),      &
	  dcen_sect_tmp(   maxd_asize, maxd_atype ),      &
	  dlo_sect_tmp(    maxd_asize, maxd_atype ),       &
	  dhi_sect_tmp(    maxd_asize, maxd_atype )


	if (iunits_flagbb /= 2) return
	if (iflagaa == 1) goto 10000
	if (iflagaa == 2) goto 20000
	write(*,'(2a,2(1x,i10))') &
	    '*** sect_iface_dens_size_testbb fatal error - ', &
	    'bad iflagaa = ', iflagaa
	stop


10000	continue
	fact_apdens = 0.001_r8   ! converts (kg/m^3) to (g/cm^3)
	fact_apdiam = 100.0_r8   ! converts (m) to (cm)

	dens_mastercomp_aer_tmp(:) = dens_mastercomp_aer(:)
	dens_aer_tmp(:,:) = dens_aer(:,:)
	dens_water_aer_tmp = dens_water_aer
	dens_nh4so4a_newnuc_tmp = dens_nh4so4a_newnuc
	dlo_sect_tmp( :,:) = dlo_sect( :,:)
	dhi_sect_tmp( :,:) = dhi_sect( :,:)
	dcen_sect_tmp(:,:) = dcen_sect(:,:)
	dcut_sect_tmp(:,:) = dcut_sect(:,:)
	volumlo_sect_tmp( :,:) = volumlo_sect( :,:)
	volumhi_sect_tmp( :,:) = volumhi_sect( :,:)
	volumcen_sect_tmp(:,:) = volumcen_sect(:,:)
	volumcut_sect_tmp(:,:) = volumcut_sect(:,:)

	tmpa = 1.0_r8/fact_apdens
	dens_mastercomp_aer(:) = dens_mastercomp_aer(:)*tmpa
	dens_aer(:,:) = dens_aer(:,:)*tmpa
	dens_water_aer = dens_water_aer*tmpa
	dens_nh4so4a_newnuc = dens_nh4so4a_newnuc*tmpa

	tmpa = 1.0_r8/fact_apdiam
	dlo_sect( :,:) = dlo_sect( :,:)*tmpa
	dhi_sect( :,:) = dhi_sect( :,:)*tmpa
	dcen_sect(:,:) = dcen_sect(:,:)*tmpa
	dcut_sect(:,:) = dcut_sect(:,:)*tmpa

	tmpa = tmpa*tmpa*tmpa
	volumlo_sect( :,:) = volumlo_sect( :,:)*tmpa
	volumhi_sect( :,:) = volumhi_sect( :,:)*tmpa
	volumcen_sect(:,:) = volumcen_sect(:,:)*tmpa
	volumcut_sect(:,:) = volumcut_sect(:,:)*tmpa

	return


20000	continue
	fact_apdens = 1.0_r8
	fact_apdiam = 1.0_r8

	dens_mastercomp_aer(:) = dens_mastercomp_aer_tmp(:)
	dens_aer(:,:) = dens_aer_tmp(:,:)
	dens_water_aer = dens_water_aer_tmp
	dens_nh4so4a_newnuc = dens_nh4so4a_newnuc_tmp

	dlo_sect( :,:) = dlo_sect_tmp( :,:)
	dhi_sect( :,:) = dhi_sect_tmp( :,:)
	dcen_sect(:,:) = dcen_sect_tmp(:,:)
	dcut_sect(:,:) = dcut_sect_tmp(:,:)

	volumlo_sect( :,:) = volumlo_sect_tmp( :,:)
	volumhi_sect( :,:) = volumhi_sect_tmp( :,:)
	volumcen_sect(:,:) = volumcen_sect_tmp(:,:)
	volumcut_sect(:,:) = volumcut_sect_tmp(:,:)

	return


	end subroutine sect_iface_dens_size_testbb


!-----------------------------------------------------------------------


	end module module_sect_iface

