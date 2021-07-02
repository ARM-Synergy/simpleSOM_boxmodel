	module module_pmcmos_aer


	use module_data_mosaic_kind, only:  r8
	use module_data_mosaic_pmcmos


	implicit none


	contains


!-----------------------------------------------------------------------
	subroutine pmcmos_init_aerosol
!
!   allocates and loads "cnn" arrays for background and emissions
!      size distributions
!   loads the main mosaic cnn array with initial size distribution
!      info
!
	use module_data_mosaic_main, only:  &
	    cnn, knum_a, lun_sect_180, naer_tot, ngas_max, ntot_max, pi
	use module_data_mosaic_aero, only:  &
	    dens_aer_mac, ibc_a, kappa_aer_mac,  &
	    method_bcfrac, method_kappa,  &
	    msectional_flag2, mw_aer_mac, naer
	use module_data_mosaic_asect, only:  &
	    dlo_sect, dhi_sect, itype_of_itype_md1md2,  &
            nsize_aer, ntype_aer, ntype_md1_aer, ntype_md2_aer,  &
            xcut_atype_md1, xcut_atype_md2

!   subr arguments (none)

!   local variables
	integer  :: i, ibin, idiagaa, imode, isize, itype, itype_use, it1, it2
	integer  :: j, jj
	integer  :: l, luna
	integer  :: noffseta, noffsetb

	real(r8) :: cnnx(ntot_max)
	real(r8) :: densmean(0:ad_maxmode), dvolmean
	real(r8) :: fracnum, fracvol
	real(r8) :: mtot(0:ad_maxmode,2)
	real(r8) :: ntot(0:ad_maxmode,2)
	real(r8) :: tmp_dens, tmp_kappa, tmp_wbc
	real(r8) :: tmpa, tmpb, tmpc
	real(r8) :: tmp_c, tmp_m, tmp_mtot, tmp_n, tmp_ntot, tmp_v, tmp_vtot
	real(r8) :: vtot(0:ad_maxmode,2)

	type( aer_dist_t ) :: ad


	allocate( aer_back_cnn( ntot_max, aer_back_nt ) )
	allocate( aer_emit_cnn( ntot_max, aer_emit_nt ) )

jj_loop_01: &
	do jj = 1, 1 + aer_back_nt + aer_emit_nt

	luna = 6

	if (jj == 1) then
	    j = 1
	    ad = aer_init_dist
	    idiagaa = -1
	    if (j == 1) write(luna,'(/2a)') 'pmcmos_init_aerosol -- init', &
		'--------------------------------------------------'
	else if (jj <= 1+aer_back_nt) then
	    j = jj - 1
	    ad = aer_back_dist(j)
	    idiagaa = -1
	    if (j == 1) write(luna,'(/2a)') 'pmcmos_init_aerosol -- back', &
		'--------------------------------------------------'
	else
	    j = jj - 1 - aer_back_nt
	    ad = aer_emit_dist(j)
	    idiagaa = -1
	    if (j == 1) write(luna,'(/2a)') 'pmcmos_init_aerosol -- emit', &
		'--------------------------------------------------'
	end if


	cnnx(:) = 0.0
	ntot(:,:) = 0.0
	vtot(:,:) = 0.0
	mtot(:,:) = 0.0

	do imode = 1, ad%n_mode

            ! calculate number, volume, and mass for the mode
	    ! ntot(imode,1) = number conc of the mode (#/cm^3)
	    if (jj <= 1+aer_back_nt) then
		ntot(imode,1) = ad%numden(imode)*1.0e-6_r8
	    else
		! for emissions, numden is #/m^2/s.  Mult by (3600 s)/(1000 m)
		! to get #/cm^3 emitted over 1 h into a 1000 m pbl
		ntot(imode,1) = ad%numden(imode)*1.0e-6_r8*(3600.0_r8/1000.0_r8)
	    end if
	    dvolmean = ad%dgnum_cm(imode) * exp( 1.5_r8 * (ad%lnsg(imode)**2) )  * 1.0e4_r8  ! um
	    ! vtot(imode) = volume conc of the mode (um^3/cm^3)
	    vtot(imode,1) = ntot(imode,1) * (pi/6.0_r8) * (dvolmean**3)
	    tmp_mtot = 0.0
	    tmp_vtot = 0.0
!	    write(*,'(/a,2i5)') 'tmp_m,v for jj,j =', jj, j
!	    write(*,'(a,1p,4e11.3)') 'dgnum, lnsg, massfrac', ad%dgnum_cm(imode), ad%lnsg(imode)
!	    write(*,'(10f8.3)') ad%massfrac(1:naer,imode)
	    do l = 1, naer
		tmp_m = ad%massfrac(l,imode)
		tmp_v = tmp_m/dens_aer_mac(l)
!		write(*,'(a,i5,1p,2e11.3)') 'l,tmp_m,v =', l, tmp_m, tmp_v
		tmp_mtot = tmp_mtot + tmp_m
		tmp_vtot = tmp_vtot + tmp_v
	    end do
	    ! (tmp_mtot/tmp_vtot) = mean dry density (g/cm^3)
	    densmean(imode) = 1.0
	    if (tmp_vtot > 1.0e-35) densmean(imode) = (tmp_mtot/tmp_vtot)
	    ! mtot(imode) = mass conc of the mode (ug/m^3)
	    mtot(imode,1) = vtot(imode,1) * densmean(imode)

	    ntot(0,1) = ntot(0,1) + ntot(imode,1)
	    vtot(0,1) = vtot(0,1) + vtot(imode,1)
	    mtot(0,1) = mtot(0,1) + mtot(imode,1)


	    ! determine which type will receive this mode
	    if (msectional_flag2 <= 0) then
		itype_use = 1
		it1 = -1 ; it2 = -1
	    else
		! bc mass fraction
		tmp_wbc = ad%massfrac(ibc_a,imode)

		! volume weighted hygroscopicity of non-bc species
		tmpa = 0.0
		tmpb = 0.0
		do l = 1, naer
		    if (method_kappa == 12) then
			if (l == ibc_a) cycle
		    end if
		    tmpc = ad%massfrac(l,imode)/dens_aer_mac(l)
		    tmpa = tmpa + tmpc*kappa_aer_mac(l)
		    tmpb = tmpb + tmpc
		end do
		tmp_kappa = max(tmpa,0.0_r8)/max(tmpb,1.0e-35_r8)

		it1 = ntype_md1_aer
		do i = 1, ntype_md1_aer-1
		    if (tmp_wbc < xcut_atype_md1(i)) then
			! if code never gets here, then it1 = ntype_md1_aer
			it1 = i
			exit
		    end if
		end do

		it2 = ntype_md2_aer
		do i = 1, ntype_md2_aer-1
		    if (tmp_kappa < xcut_atype_md2(i)) then
			it2 = i
			exit
		    end if
		end do
		itype_use = itype_of_itype_md1md2(it1,it2)
	    end if
	    if (imode == 1) write(luna,'(a)')
	    write(luna,'(a,6i5)') 'pmcmos_init_aerosol -- jj,j,imode,itype_use,it1,it2', &
		jj, j, imode, itype_use, it1, it2


            ! now do size bin calculations
	    ibin = 0
	    noffseta = ngas_max - naer_tot
	    noffsetb = ngas_max + (naer_tot-naer) - naer_tot

	    do itype = 1, itype_use
	    do isize = 1, nsize_aer(itype)
		noffseta = noffseta + naer_tot
		noffsetb = noffsetb + naer_tot
		ibin = ibin + 1
		if (itype /= itype_use) cycle

		call modal_to_bin_aa( ad%dgnum_cm(imode), ad%lnsg(imode), &
		    dlo_sect(isize,itype), dhi_sect(isize,itype), fracnum, fracvol )
		tmp_n = ntot(imode,1)*fracnum
		tmp_v = vtot(imode,1)*fracvol
		tmp_m = mtot(imode,1)*fracvol
		ntot(imode,2) = ntot(imode,2) + tmp_n
		vtot(imode,2) = vtot(imode,2) + tmp_v
		mtot(imode,2) = mtot(imode,2) + tmp_m

		! cnn for dry-mass species are umol/m^3, so just divide by mw
		do l = 1, naer
		    tmp_c = tmp_m*ad%massfrac(l,imode)/mw_aer_mac(l)
		    cnnx(noffsetb+l) = cnnx(noffsetb+l) + tmp_c
		end do

		! cnn for number is #/cm^3
		cnnx(noffseta+knum_a) = cnnx(noffseta+knum_a) + tmp_n

	    end do ! isize
	    end do ! itype

	    ntot(0,2) = ntot(0,2) + ntot(imode,2)
	    vtot(0,2) = vtot(0,2) + vtot(imode,2)
	    mtot(0,2) = mtot(0,2) + mtot(imode,2)
	end do ! imode

! diagnostics
	luna = lun_sect_180
	luna = 6
	write(luna,'(/a,2i5)') 'pmcmos_init_aerosol -- jj,j', jj, j
	do imode = 0, ad%n_mode
	    write(luna,'(i5,1p,4(2x,2e11.3))') imode, ntot(imode,1:2), &
		vtot(imode,1:2), mtot(imode,1:2), densmean(imode)
	end do


	if (jj == 1) then
	    cnn(ngas_max+1:ntot_max) = cnnx(ngas_max+1:ntot_max)
	else if (jj <= 1+aer_back_nt) then
	    j = jj - 1
	    aer_back_cnn(:,j) = cnnx(:)
	else
	    j = jj - 1 - aer_back_nt
	    aer_emit_cnn(:,j) = cnnx(:)
	end if


	end do jj_loop_01


	return
	end subroutine pmcmos_init_aerosol


!-----------------------------------------------------------------------
	subroutine modal_to_bin_aa( dgnum, lnsg, dlo, dhi, fracnum, fracvol )
!
! dgnum = number median diameter of a log-normal size distribution
! lnsg = natural log of the geometric standard deviation
! dlo, dhi = lower and upper diameters of a size bin
! fracnum = fraction of the mode's number in the size bin
! fracvol = fraction of the mode's volume in the size bin
!
	implicit none
	real(r8), intent(in)  :: dgnum, lnsg, dlo, dhi
	real(r8), intent(out) :: fracnum, fracvol

	real(r8) :: sxroot2, tlo, thi, x0, x3, xlo, xhi

	real(r8) :: derfc

	sxroot2 = lnsg * sqrt( 2.0_r8 )

	xlo = log( dlo )
	xhi = log( dhi )
	x0 = log( dgnum )
	tlo = (xlo - x0)/sxroot2
	thi = (xhi - x0)/sxroot2
	if (tlo <= 0.0) then
	    fracnum = 0.5d0*( derfc(-thi) - derfc(-tlo) )
	else
	    fracnum = 0.5d0*( derfc(tlo) - derfc(thi) )
	end if

	x3 = x0 + 3.0d0*lnsg*lnsg
	tlo = (xlo - x3)/sxroot2
	thi = (xhi - x3)/sxroot2
	if (tlo <= 0.0) then
	    fracvol = 0.5d0*( derfc(-thi) - derfc(-tlo) )
	else
	    fracvol = 0.5d0*( derfc(tlo) - derfc(thi) )
	end if

	return
	end subroutine modal_to_bin_aa


!-----------------------------------------------------------------------
	subroutine pmcmos_print( it )
!
!   allocates and loads "cnn" arrays for background and emissions
!      size distributions
!   loads the main mosaic cnn array with initial size distribution
!      info
!
	use module_data_mosaic_main, only:  &
	    cair_mlc, cnn, inputfile, iprint, knum_a, kwater_a,  &
	    kso4_a, kno3_a, kcl_a, knh4_a, koc_a, kbc_a, &
        kcn3_a,     kcn2_a,     kcn1_a,     kc0_a,    &
        kc1_a,      kc2_a,      kc3_a,      kc4_a,    &
        kc5_a,      kc6_a,      kc7_a,                &
		kcn3,       kcn2,       kcn1,       kc0,      &
		kc1,        kc2,        kc3,        kc4,      &
		kc5,        kc6,        kc7,        kc8,      &
		kc9,    &
	    koh, ko3, kno, kno2, kno3, kn2o5, khno3, kpan, khono, khno4, konit, knap, kisopn, &
	    lun_sect_183, lun_sect_184, lun_sect_185, lun_sect_186,  &
	    lun_sect_188,  &
	    naer_tot, ngas_max, ntot_max, pi, ppb, species, time_sec

	use module_data_mosaic_aero, only:  &
	    aer_name, dens_aer_mac, jaerosolstate, &
	    method_bcfrac, method_kappa, &
	    msectional, msectional_flag2, msize_framework, &
	    mw_aer_mac, naer, nbin_a, nbin_a_max, &
	    iso4_a, ino3_a, icl_a, inh4_a, ioc_a, ibc_a, &
		icn3_a,     icn2_a,    icn1_a,     ic0_a,     &
        ic1_a,      ic2_a,      ic3_a,     ic4_a,      ic5_a,     &
        ic6_a,      ic7_a,      ic8_a,      ic9_a

	use module_data_mosaic_asect, only:  &
	    dcen_sect, dlo_sect, dhi_sect,  &
	    itype_of_ibin, itype_md1_of_itype, itype_md2_of_itype,  &
	    nsize_aer, ntype_aer, ntype_md1_aer, ntype_md2_aer,  &
	    xcut_atype_md1, xcut_atype_md2

!   subr arguments (none)
	integer  :: it

!   local variables
	integer  :: i, ibin, idiagaa, imode, isize, itype, it1, it2
	integer  :: j, jj
	integer  :: l, ltmp, lun
	integer  :: n, noffseta, noffsetb, noffsetc

	real(r8) :: cnnx(ntot_max)
	real(r8) :: dlog10dp, dvolmean
	real(r8) :: fracnum, fracvol
	real(r8) :: mtot(0:ad_maxmode,2)
	real(r8) :: ntot(0:ad_maxmode,2)
	real(r8) :: tmpa, tmpb, tmpc, tmpd, tmpe
	real(r8) :: tmp_dens
	real(r8) :: tmp_m, tmp_mtot, tmp_n, tmp_ntot, tmp_v, tmp_vtot
	real(r8) :: tmpveca(naer_tot), tmpvecb(naer_tot)
	real(r8) :: vtot(0:ad_maxmode,2)
	real(r8) :: xdcen(0:nbin_a_max)
	real(r8) :: xmass(0:nbin_a_max,0:naer+1)
	real(r8) :: xnumb(0:nbin_a_max)
	real(r8) :: xtot_vol

	character(len=80) :: txtaa

	type( aer_dist_t ) :: ad



!	if (time_sec <= 1.0e-10) then
!	if (it == 0) then
!	    write(lun_sect_183,*) 'dens_aer_mac*1000'
!	    write(lun_sect_183,'(a,i6)') (aer_name(l), &
!		nint( dens_aer_mac(l)*1000.0 ), l=1,naer)
!	end if

	if (mod(it,iprint) /= 0) return

	lun = lun_sect_183 ! output size distributed number and masses
!	  write(lun,'(/f10.1, a)') time_sec/60, 'min'
	  write(lun,'(/a,16(a,7x))') '     time      bin s   dcen   ', 'ntot', 'mdry', 'mh2o', &
	        aer_name(iso4_a)(1:4), aer_name(inh4_a)(1:4), aer_name(ioc_a)(1:4), &
	        aer_name(icn3_a)(1:4),&
			aer_name(icn2_a)(1:4),&
	        aer_name(icn1_a)(1:4),&
			aer_name(ic0_a)(1:4),&
			aer_name(ic1_a)(1:4),&
			aer_name(ic2_a)(1:4),&
			aer_name(ic3_a)(1:4),&
			aer_name(ic4_a)(1:4),&
		    aer_name(ic5_a)(1:4),&
			aer_name(ic6_a)(1:4),&
			aer_name(ic7_a)(1:4),&
			aer_name(ic8_a)(1:4),&
			aer_name(ic9_a)(1:4),&
			'SOAt'

	if(it.eq.32 .or. it.eq.64 .or. it.eq.128 .or. it.eq.256 .or. it.eq.512 .or. it.eq.1024)then
	  lun = it
	  write(lun,'(/a,16(a,7x))') '     time      bin s   dcen   ', 'ntot', 'mdry', 'mh2o', &
	        aer_name(iso4_a)(1:4), aer_name(inh4_a)(1:4), aer_name(ioc_a)(1:4), &
	        aer_name(icn3_a)(1:4),&
			aer_name(icn2_a)(1:4),&
	        aer_name(icn1_a)(1:4),&
			aer_name(ic0_a)(1:4),&
			aer_name(ic1_a)(1:4),&
			aer_name(ic2_a)(1:4),&
			aer_name(ic3_a)(1:4),&
			aer_name(ic4_a)(1:4),&
			aer_name(ic5_a)(1:4),&
			aer_name(ic6_a)(1:4),&
			aer_name(ic7_a)(1:4),&
			aer_name(ic8_a)(1:4),&
			aer_name(ic9_a)(1:4),&
			'SOAt'
	endif

	lun = lun_sect_184 ! output total number and species masses (summed over all bins)
	  if (it == 0) &
	    write(lun,'(14(a,7x))') 'time_min', 'ntot', 'vtot', 'mdry', 'mh2o', &
	        aer_name(iso4_a)(1:4), aer_name(inh4_a)(1:4), aer_name(ioc_a)(1:4), &
	        aer_name(icn3_a)(1:4),&
			aer_name(icn2_a)(1:4),&
	        aer_name(icn1_a)(1:4),&
			aer_name(ic0_a)(1:4),&
			aer_name(ic1_a)(1:4),&
			aer_name(ic2_a)(1:4),&
			aer_name(ic3_a)(1:4),&
			aer_name(ic4_a)(1:4),&
			aer_name(ic5_a)(1:4),&
			aer_name(ic6_a)(1:4),&
			aer_name(ic7_a)(1:4),&
			aer_name(ic8_a)(1:4),&
			aer_name(ic9_a)(1:4),&
			'SOAt'

!	lun = lun_sect_185 ! output total number and total aerosol mass
!	  if (it == 0) then
!	    write(lun,'(a)') 'numb (#/cm3) and mass(1:25) (ug/m3)'
!	    write(lun,'(24(a,5x))') (aer_name(i), i=1,23), '  WATER'
!        endif


!	if (msize_framework == msectional) then
!	    write(lun_sect_183,'(a)') &
!		'masses = dm/dlog10dp (ug/m3) ; number = dn/dlog10dp (#/cm3)'
!	else
!	    write(lun_sect_183,'(a)') &
!		'masses = mass in bin (ug/m3) ; number = number in bin (#/cm3)'
!	end if


	xnumb(:) = 0.0
	xmass(:,:) = 0.0

	noffseta = ngas_max
	noffsetb = ngas_max + (naer_tot - naer)
	i = 1

	do itype = 1, ntype_aer
	do isize = 1, nsize_aer(itype)
	    do l = 1, naer
            xmass(i,l) = cnn(noffsetb+l)*mw_aer_mac(l) ! ug/m^3
	    end do
	    xnumb(i) = cnn(noffseta+knum_a) ! #/cm^3
	    xmass(i,naer+1) = cnn(noffseta+kwater_a)*1.0e9 ! kg/m^3 --> ug/m^3

	    xmass(i,0) = sum( xmass(i,1:naer) )

	    xdcen(i) = dcen_sect(isize,itype)*1.0e4 ! um
	    i = i + 1
	    noffseta = noffseta + naer_tot
	    noffsetb = noffsetb + naer_tot
	end do ! isize
	end do ! itype

	do l = 0, naer+1
	    xmass(0,l) = sum( xmass(1:nbin_a,l) )
	end do
	xnumb(0) = sum( xnumb(1:nbin_a) )

! output size distributed number and species masses
	lun = lun_sect_183
	do i = 1, nbin_a
	    tmpa = 1.0
	    if (msize_framework == msectional) then
		itype = itype_of_ibin( i )
		n = nsize_aer(itype)
		if (n > 1) then
		    dlog10dp = log10( dcen_sect(n,itype) &
		                    / dcen_sect(1,itype) ) / (n-1)
!		    tmpa = 1.0/dlog10dp		! use this to report num and masses/dlog10dp
		    tmpa = 1.0			! use this to report num and masses in bin
		end if
	    end if
	    write(lun,'(f10.1,2x,i6,i2,f8.2,1p,12e11.3)') &
		time_sec/60., i, jaerosolstate(i), xdcen(i)*1000.0, &
		xnumb(i)*tmpa,    xmass(i, 0)*tmpa, &
		xmass(i,20)*tmpa,                   &
		xmass(i,iso4_a)*tmpa, xmass(i,inh4_a)*tmpa, xmass(i,ioc_a)*tmpa, &
		xmass(i,icn3_a)*tmpa, &
		xmass(i,icn2_a)*tmpa, &
		xmass(i,icn1_a)*tmpa, &
		xmass(i,ic0_a)*tmpa, &
		xmass(i,ic1_a)*tmpa, &
		xmass(i,ic2_a)*tmpa, &
		xmass(i,ic3_a)*tmpa, &
		xmass(i,ic4_a)*tmpa, &
		xmass(i,ic5_a)*tmpa, &
		xmass(i,ic6_a)*tmpa, &
		xmass(i,ic7_a)*tmpa, &
		xmass(i,ic8_a)*tmpa, &
		xmass(i,ic9_a)*tmpa, &
		sum(xmass(i,6:17))*tmpa

	   if(it.eq.32 .or. it.eq.64 .or. it.eq.128 .or. it.eq.256 .or. it.eq.512 .or. it.eq.1024)then
	     write(it,'(f10.1,2x,i6,i2,f8.2,1p,12e11.3)') &
		time_sec/60., i, jaerosolstate(i), xdcen(i)*1000.0, &
		xnumb(i)*tmpa,    xmass(i, 0)*tmpa, &
		xmass(i,20)*tmpa,                   &
		xmass(i,iso4_a)*tmpa, xmass(i,inh4_a)*tmpa, xmass(i,ioc_a)*tmpa, &
		xmass(i,icn3_a)*tmpa, &
		xmass(i,icn2_a)*tmpa, &
		xmass(i,icn1_a)*tmpa, &
		xmass(i,ic0_a)*tmpa, &
		xmass(i,ic1_a)*tmpa, &
		xmass(i,ic2_a)*tmpa, &
		xmass(i,ic3_a)*tmpa, &
		xmass(i,ic4_a)*tmpa, &
		xmass(i,ic5_a)*tmpa, &
		xmass(i,ic6_a)*tmpa, &
		xmass(i,ic7_a)*tmpa, &
        xmass(i,ic8_a)*tmpa, &
		xmass(i,ic9_a)*tmpa, &
		sum(xmass(i,6:17))*tmpa
	   endif

	end do

! output total number and species masses (summed over all bins)
	lun = lun_sect_184
      xtot_vol = 0.0
	do i = 1, nbin_a
	  xtot_vol = xtot_vol + xnumb(i)*(pi/6.0)*(xdcen(i)*1000.0)**3
	enddo

	i = 0
	write(lun,'(f10.1,1p,13e11.3)') time_sec/60.0, &
		xnumb(i),    xtot_vol, xmass(i,0), &
		xmass(i,20),                   &
		xmass(i,iso4_a), xmass(i,inh4_a), xmass(i,ioc_a), &
		xmass(i,icn3_a), &
		xmass(i,icn2_a), &
		xmass(i,icn1_a), &
		xmass(i,ic0_a), &
		xmass(i,ic1_a), &
		xmass(i,ic2_a), &
		xmass(i,ic3_a), &
		xmass(i,ic4_a), &
		xmass(i,ic5_a), &
		xmass(i,ic6_a), &
		xmass(i,ic7_a), &
        xmass(i,ic8_a), &
		xmass(i,ic9_a), &
		sum(xmass(i,6:17))

! output total number and total aerosol mass
!	lun = lun_sect_185
!	write(lun,'(1p,8e12.4)') &
!		xnumb(i),    xmass(i, 1:20)


	lun = lun_sect_188
	if (msectional_flag2 > 0) then
!	if (msectional_flag2 > -999888999) then
	    write(lun,'(a)')
	    write(lun,'(a)') 'IVERSION 0004'
	    write(lun,'(a)') '10 params = num, so4, no3, cl, nh4, oc, bc, api1, other_soa, water'
	    write(lun,'(a)') trim(inputfile)
	    write(lun,'(3i4,i10,1p,2e14.6,a)') nsize_aer(1), &
		ntype_md1_aer, ntype_md2_aer, -1, &
		(time_sec/3600.0), -1.0, &
		'     ndpbin, nbcbin, nkappabin, npart, time (h), ntot (#/cm3)'

	    write(lun,'(a)') 'Dp_cut values'
	    tmpd = 1.0e4
	    write(lun,'(1p,5e14.6)') tmpd*dlo_sect(1,1), (tmpd*dhi_sect(i,1), i=1,nsize_aer(1))

	    if (method_bcfrac == 1) then
		txtaa = ' (BC dry-mass fraction)'
	    else
		write(*,'(/a,1x,i20/)') &
		    '*** pmcmos_print - bad method_bcfrac = ', method_bcfrac
		stop
	    end if
	    write(lun,'(a,i2,a)') 'wBC_cut values     method_bcfrac=', method_bcfrac, trim(txtaa)
	    write(lun,'(1p,5e14.6)') xcut_atype_md1(0:ntype_md1_aer)

	    if (method_kappa == 11) then
		txtaa = ' (all species using ion/aer kappas)'
	    else if (method_kappa == 12) then
		txtaa = ' (non-bc species using ion/aer kappas)'
	    else
		write(*,'(/a,1x,i20/)') &
		    '*** pmcmos_print - bad method_kappa = ', method_kappa
		stop
	    end if
	    write(lun,'(a,i2,a)') 'kappa_cut values     method_kappa=', method_kappa, trim(txtaa)
	    write(lun,'(1p,5e14.6)') xcut_atype_md2(0:ntype_md2_aer)

	    noffseta = ngas_max
	    noffsetc = naer_tot - naer
	    do itype = 1, ntype_aer
		if (msectional_flag2 > 0) then
		    it1 = itype_md1_of_itype(itype)
		    it2 = itype_md2_of_itype(itype)
		else
		    it1 = 1 ; it2 = 1
		end if
		tmpd =        xcut_atype_md1(it1) - xcut_atype_md1(it1-1)
		tmpe = log10( xcut_atype_md2(it2) / xcut_atype_md2(it2-1) )
		do isize = 1, nsize_aer(itype)
		    tmpc = log10( dhi_sect(isize,itype) / dlo_sect(isize,itype) )
		    tmpa = 1.0_r8/(tmpc*tmpd*tmpe)

		    tmpveca(1:naer_tot) = max( 0.0_r8, cnn(noffseta+1:noffseta+naer_tot) )
		    do l = 1, naer
		        tmpveca(noffsetc+l) = tmpveca(noffsetc+l)*mw_aer_mac(l)   ! umol/m^3 to ug/m^3
		    end do

		    tmpvecb( 1) = tmpveca(knum_a )   ! #/cm^3
		    tmpvecb( 2) = tmpveca(kso4_a )   ! ug/m^3
		    tmpvecb( 3) = tmpveca(kno3_a )
		    tmpvecb( 4) = tmpveca(kcl_a  )
		    tmpvecb( 5) = tmpveca(knh4_a )
		    tmpvecb( 6) = tmpveca(koc_a  )
		    tmpvecb( 7) = tmpveca(kbc_a  )
		    ! soa species order is aro1, aro2, alk1, ole1, api1, api2, lim1, lim2
		    tmpvecb( 8) = sum( tmpveca(kcn3_a:kc7_a) )
		    tmpvecb(9) = tmpveca(kwater_a)*1.0e9   ! kg/m^3 to ug/m^3
		    tmpvecb(1:9) = tmpvecb(1:9)*tmpa
		    tmpb = tmp_n/tmpc
		    write(lun,'(3i4,1p,30e12.5)') isize, it1, it2, tmpvecb(1:10)

		    noffseta = noffseta + naer_tot
		end do
	    end do
	end if ! (msectional_flag2 > 0)


! some gases
! h2so4=1    hno3=2    nh3=4    no2=6    o3=11    so2=18    hcho=23
	lun = lun_sect_186
	if (it == 0) &
	write(lun,'(2x,19(a,7x))') &
	    'time_min', &
	    species(koh)(1:4), species(ko3)(1:4), &
	    species(kcn3)(1:4), &
		species(kcn2)(1:4), &
	    species(kcn1)(1:4), &
		species(kc0)(1:4),  &
	    species(kc1)(1:4),  &
		species(kc2)(1:4),  &
		species(kc3)(1:4),  &
		species(kc4)(1:4),  &
		species(kc5)(1:4),  &
		species(kc6)(1:4),  &
		species(kc7)(1:4),  &
        species(kc8)(1:4),  &
		species(kc9)(1:4),  &
	    species(kno)(1:4), species(kno2)(1:4), species(kno3)(1:4), species(kn2o5)(1:4), &
	    species(khno3)(1:4), species(konit)(1:4), species(khono)(1:4), species(khno4)(1:4), &
          'NOy'
!	write(lun,'(/f10.1,a,5x,a,i3.3)') time_sec, ' = time', 'mosaic_box fort.', lun
	tmpa = ppb/cair_mlc
	write(lun,'(f10.1, 2x, 19e11.3)') &
          time_sec/60.0, &
	    cnn(koh)*tmpa, cnn(ko3)*tmpa, &
	    cnn(kcn3)*tmpa, &
		cnn(kcn2)*tmpa, &
	    cnn(kcn1)*tmpa, &
		cnn(kc0)*tmpa,  &
	    cnn(kc1)*tmpa,  &
		cnn(kc2)*tmpa,  &
		cnn(kc3)*tmpa,  &
		cnn(kc4)*tmpa,  &
		cnn(kc5)*tmpa,  &
		cnn(kc6)*tmpa,  &
		cnn(kc7)*tmpa,  &
        cnn(kc8)*tmpa,  &
		cnn(kc9)*tmpa,  &
	    cnn(kno)*tmpa, cnn(kno2)*tmpa, cnn(kno3)*tmpa, cnn(kn2o5)*tmpa, &
	    cnn(khno3)*tmpa, cnn(konit)*tmpa, cnn(khono)*tmpa, cnn(khno4)*tmpa, &
	   (cnn(kno)+cnn(kno2)+cnn(kno3)+2.*cnn(kn2o5)+cnn(khno3)+cnn(khono)+ &
          cnn(kpan)+cnn(khno4)+cnn(konit)+cnn(knap)+cnn(kisopn))*tmpa


	return
	end subroutine pmcmos_print


!-----------------------------------------------------------------------


	end module module_pmcmos_aer

