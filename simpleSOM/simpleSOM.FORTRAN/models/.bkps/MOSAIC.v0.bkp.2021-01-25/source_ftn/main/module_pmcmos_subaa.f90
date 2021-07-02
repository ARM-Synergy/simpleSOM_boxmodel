	module module_pmcmos_subaa


	use module_data_mosaic_kind, only:  r8
	use module_data_mosaic_pmcmos


	implicit none


	contains


!-----------------------------------------------------------------------
	subroutine pmcmos_emit_dilu
!
!   applies changes to gas and aerosol concentrations due to
!      air density changes
!      emissions
!      dilution/mixing with background
!   using same approach as in partmc-mosaic
!
	use module_data_mosaic_main, only:  &
	    avogad, cair_mlc, cair_mlc_old, cair_molm3, cnn, dt_sec, &
	    istate_pblh, knum_a, kwater_a, naer_tot, ngas_max, &
	    pblh, pblh_old, ppb, te, time_hrs, time_sec
	use module_data_mosaic_aero, only:  &
	    naer
	use module_data_mosaic_asect, only:  &
	    nsize_aer, ntype_aer


!   subr arguments (none)

!   local variables
	integer :: i, ibin, isize, itype
	integer :: jt_aer_back, jt_aer_emit
	integer :: jt_gas_back, jt_gas_emit
	integer :: l, ll
	integer :: noffseta, noffsetb

	real(r8) :: ecur_tmp, esat_tmp
	real(r8) :: press_tmp, rh_tmp, te_tmp
	real(r8) :: tmpa
	real(r8) :: tmp_dilu, tmp_horzdilu, tmp_vertdilu
	real(r8) :: tmp_back1, tmp_back2
	real(r8) :: tmp_emit1, tmp_emit2, tmp_emitfac


!   do changes from air density
	tmpa = cair_pmcmos/cair_pmcmos_old
	tmpa = cair_mlc/cair_mlc_old
	do i = 1, ngas_max
	    cnn(i) = cnn(i)*tmpa
	end do

	noffseta = ngas_max
	noffsetb = ngas_max + (naer_tot - naer)
	ibin = 1

	do itype = 1, ntype_aer
	do isize = 1, nsize_aer(itype)
	    do l = -1, naer
		if (l == -1) then
		    ll = noffseta + knum_a
		else if (l == 0) then
		    ll = noffseta + kwater_a
		else
		    ll = noffsetb + l
		end if
		cnn(ll) = cnn(ll)*tmpa
	    end do

	    ibin = ibin + 1
	    noffseta = noffseta + naer_tot
	    noffsetb = noffsetb + naer_tot
	end do ! isize
	end do ! itype


!   do changes to gases from emissions and mixing
	call pmcmos_time_interp2( gas_emit_tim, gas_emit_nt, &
	    time_sec, jt_gas_emit )
	tmp_emitfac = gas_emit_fac(jt_gas_emit)
	! minor adjustment to better match partmc-mosaic
	tmp_emitfac = tmp_emitfac * (cair_molm3/cair_pmcmos)

	call pmcmos_time_interp2( gas_back_tim, gas_back_nt, &
	    time_sec, jt_gas_back )
	tmp_horzdilu = gas_back_dil(jt_gas_back)
	if ((istate_pblh == 1) .and. (pblh > pblh_old)) then
	    tmp_vertdilu = (pblh - pblh_old)/(pblh_old*dt_sec)
	else
	    tmp_vertdilu = 0.0
	end if
	tmp_dilu = tmp_horzdilu + tmp_vertdilu

	do i = 1, ngas_max
	    tmp_emit1 = gas_emit_val(jt_gas_emit,i)
	    tmp_back1 = gas_back_val(jt_gas_back,i)

	    tmp_back2 = tmp_back1*cair_mlc/ppb   ! convert from ppb to molecu/cm3
	    tmp_emit2 = (tmp_emit1/(pblh*1.0e6))*avogad   ! convert from mole/m2/s to molecu/cm3/s

	    cnn(i) = cnn(i) + tmp_emitfac*tmp_emit2*dt_sec + (tmp_back2 - cnn(i))*tmp_dilu*dt_sec
	end do


!   do changes to aerosols from emissions and mixing
	call pmcmos_time_interp2( aer_emit_tim, aer_emit_nt, &
	    time_sec, jt_aer_emit )
	tmp_emitfac = aer_emit_fac(jt_aer_emit)

	call pmcmos_time_interp2( aer_back_tim, aer_back_nt, &
	    time_sec, jt_aer_back )
	tmp_horzdilu = aer_back_dil(jt_aer_back)
	tmp_dilu = tmp_horzdilu + tmp_vertdilu

! *** need to calc current aerosol water 
!     for background and emitted particles,
!     consistent with partmc

	noffseta = ngas_max
	noffsetb = ngas_max + (naer_tot - naer)
	ibin = 1

	do itype = 1, ntype_aer
	do isize = 1, nsize_aer(itype)
	    do l = -1, naer
		if (l == -1) then
		    ll = noffseta + knum_a
		else if (l == 0) then
		    ll = noffseta + kwater_a
		else
		    ll = noffsetb + l
		end if

		tmp_emit1 = aer_emit_cnn(ll,jt_aer_emit)
		! the aer_emit_ccn values are concentrations (ug/m^3 or #/cm^3)
		!    resulting from emission over 1 h into a 1000 m pbl
		! the next line adjusts for the current pblh and 
		!    removes the (3600 s) factor
		tmp_emit2 = tmp_emit1*(1000.0_r8/pblh)/3600.0_r8

		tmp_back2 = aer_back_cnn(ll,jt_aer_back)

		cnn(ll) = cnn(ll) + tmp_emitfac*tmp_emit2*dt_sec + (tmp_back2 - cnn(ll))*tmp_dilu*dt_sec
	    end do

	    ibin = ibin + 1
	    noffseta = noffseta + naer_tot
	    noffsetb = noffsetb + naer_tot
	end do ! isize
	end do ! itype



	return
	end subroutine pmcmos_emit_dilu


!-----------------------------------------------------------------------
	subroutine pmcmos_update_met( iflagaa )
!
!   updates temperature, air density, and relative humidity
!   using same approach as in partmc-mosaic
!
	use module_data_mosaic_main, only:  &
	    avogad, cair_mlc, cair_mlc_old, cair_molm3, cair_molm3_old, &
	    istate_pblh, &
	    pblh, pblh_old, pr_atm, pr_atm_old, press0_pa, &
	    rh, rh_old, te, te_old, time_hrs, time_sec


!   subr arguments
	integer, intent(in) :: iflagaa

!   local variables
	real(r8) :: ecur_tmp, esat_tmp
	real(r8) :: press_tmp, rh_tmp, te_tmp

!   set temperature
	te_old = te
	call pmcmos_time_interp1( temp_tim, temp_val, temp_nt, &
	    time_sec, te_tmp )
	if (iflagaa <= 0) te_old = te_tmp
	te = te_tmp   ! *** eventually this should be active

!   currently pressure is constant
	if (iflagaa <= 0) then
	    pr_atm_old = pr_atm
	end if
	press_tmp = pr_atm*press0_pa   ! [Pa]

!   set air density 
!   change this to dry air density ???
	cair_mlc_old = cair_mlc
	cair_molm3_old = cair_molm3
	cair_mlc = avogad*pr_atm/(82.056*te)	! air conc [molec/cc]
	cair_molm3 = 1.e6*pr_atm/(82.056*te)	! air conc [mol/m^3]
	cair_pmcmos_old = cair_pmcmos
	cair_pmcmos = press_tmp/(univ_gas_const_pmcmos*te_tmp)
	if (iflagaa <= 0) then
	    cair_mlc_old = cair_mlc
	    cair_molm3_old = cair_molm3
	    cair_pmcmos_old = cair_pmcmos
	end if

!   set relative humidity
!
!   method 1 -- rh = qv/qvsat = (ambient mixing ratio)/(saturation ...)
!   method 2 -- rh =  e/esat  = (ambient vapor pressure)/(saturation ...)
!   method 2 appears to best match the partmc rh
!
	rh_old = rh
	esat_tmp = pmcmos_esat_liq( te_tmp )
	if (iflagaa <= 0) then
	    ecur_tmp = esat_tmp  * (rh*0.01_r8)
	    qh2o = mw_h2o_air * ecur_tmp / (press_tmp - ecur_tmp)
	    rh_tmp = rh
	else
	    ecur_tmp = qh2o * press_tmp / (mw_h2o_air + qh2o)
	    rh_tmp = 100.0_r8 * ecur_tmp / esat_tmp
	    rh_tmp = max( 0.0_r8, min( 100.0_r8, rh_tmp ) )
	end if

	if (iflagaa <= 0) rh_old = rh_tmp
	rh = rh_tmp   ! *** eventually this should be active

!   set pblh
	pblh_old = pblh
	call pmcmos_time_interp1( pblh_tim, pblh_val, pblh_nt, &
	    time_sec, pblh )
	if (iflagaa <= 0) pblh_old = pblh

!   istate_pblh starts at 0
!   when pblh is increasing and istate_pblh==0, it changes to 1
!   when pblh is decreasing and istate_pblh==1, it changes to 2
	if (iflagaa <= 0) then
	    istate_pblh = 0
	else if (pblh > pblh_old) then
	    if (istate_pblh == 0) istate_pblh = 1
	else if (pblh < pblh_old) then
	    if (istate_pblh == 1) istate_pblh = 2
	end if

!	write(*,'(/a)') 't_cur (s,h), temp, press, rh, pblh'
!	write(*,'(1p,10e13.5)') time_sec, time_hrs, te_tmp, &
!	    (pr_atm*press0_pa), rh_tmp, pblh
	write(*,'(a,f10.1)') 'time (s)', time_sec


	return
	end subroutine pmcmos_update_met


!-----------------------------------------------------------------------
	subroutine pmcmos_time_interp1( x_tim, x_val, nt, &
	    cur_tim, cur_val )
!
!   does linear interpolation of a time-varying array,
!      with no extrapolation beyond end points
!

!   subr arguments
	integer, intent(in) :: nt
	real(r8), intent(in) :: cur_tim, x_tim(nt), x_val(nt)
	real(r8), intent(out) :: cur_val

!   local variables
	integer :: j
	real(r8) :: tmpa

	if ((nt <= 1) .or. (cur_tim <= x_tim(1))) then
	    cur_val = x_val(1)
	else if (cur_tim >= x_tim(nt)) then
	    cur_val = x_val(nt)
	else
	    do j = 1, nt-1
		if (cur_tim <= x_tim(j+1)) then
		    tmpa = (cur_tim - x_tim(j))/(x_tim(j+1) - x_tim(j))
		    tmpa = max( 0.0_r8, min( 1.0_r8, tmpa ) )
		    cur_val = x_val(j)*(1.0_r8 - tmpa) + x_val(j+1)*tmpa
		    exit
		end if
	    end do
	end if

	return
	end subroutine pmcmos_time_interp1


!-----------------------------------------------------------------------
	subroutine pmcmos_time_interp2( x_tim, nt, &
	    cur_tim, cur_jt )
!
!   does "step-wise" interpolation of a time-varying array
!   *** step-wise means that the data being interpolated has
!      the form of a step function
!   *** this routine returns a pointer (index) to the 
!      interpolated value, rather than the value itself
!

!   subr arguments
	integer, intent(in) :: nt
	integer, intent(out) :: cur_jt
	real(r8), intent(in) :: cur_tim, x_tim(nt)

!   local variables
	integer :: j
	real(r8) :: tmpa

	if ((nt <= 1) .or. (cur_tim <= x_tim(1))) then
	    cur_jt = 1
	else if (cur_tim >= x_tim(nt)) then
	    cur_jt = nt
	else
	    do j = 1, nt-1
		if (cur_tim < x_tim(j+1)) then
		    cur_jt = j
		    exit
		end if
	    end do
	end if

	return
	end subroutine pmcmos_time_interp2


!----------------------------------------------------------------------
	real(r8) function pmcmos_qvsat( temp, press )
! returns saturation mixing ratio over (liquid) water (kg/kg)
	implicit none
! arguments
	real(r8), intent(in) :: temp  ! temperature (K)
	real(r8), intent(in) :: press ! pressure (Pa)  
! local
	real(r8) :: esat

	esat = pmcmos_esat_liq( temp )
	pmcmos_qvsat = mw_h2o_air * esat/(press-esat)

	end function pmcmos_qvsat


!----------------------------------------------------------------------
	real(r8) function pmcmos_esat_liq( temp )
! returns saturation vapor pressure over (liquid) water (Pa)
	implicit none
! arguments
	real(r8), intent(in) :: temp ! temperature (K)

	pmcmos_esat_liq = 611.0d0 * &
	 10d0**( 7.45d0*(temp - 273.15d0)/(temp - 38d0) )

	end function pmcmos_esat_liq


!-----------------------------------------------------------------------


	end module module_pmcmos_subaa

