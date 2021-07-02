!   code history
!   6/25/2008  raz - updated nh4no3 and nh4cl condensation algorithm
!   feb 22. new flux_mix

!***********************************************************************
! ASTEM: Adaptive Step Time-Split Euler Method
!
! author: Rahul A. Zaveri
! update: jan 2007
!-----------------------------------------------------------------------
      SUBROUTINE ASTEM(dtchem)
      
      use module_data_mosaic_main
      use module_data_mosaic_aero

      implicit none

! subr arguments
      real(r8) :: dtchem
! local variables
      integer ibin

      logical first
      save first
      data first/.true./

      
! update ASTEM call counter
      jASTEM_call  = jASTEM_call

! reset input print flag
      iprint_input = mYES

      !PRINT*, 'HERE4-A'
      !PRINT*, Dp_dry_a
      !READ(*,*)

! compute aerosol phase state before starting integration
      do ibin = 1, nbin_a
         if(jaerosolstate(ibin) .ne. no_aerosol)then

           
            
           call aerosol_phase_state(ibin)

           !PRINT*, 'HERE5-A'
           !PRINT*, Dp_dry_a
           !READ(*,*)
           
           ! MARK HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           call calc_dry_n_wet_aerosol_props(ibin)

           !PRINT*, 'HERE5-B'
           !PRINT*, Dp_dry_a
           !READ(*,*)
           
        endif
      enddo

      !PRINT*, 'HERE4-B'
      !PRINT*, Dp_dry_a
      !READ(*,*)
      

! BOX
!      if(first)then		! UNCOMMENT THIS LINE
!        first=.false.		! UNCOMMENT THIS LINE
!        call print_aer(0)	! UNCOMMENT THIS LINE
!      endif			! UNCOMMENT THIS LINE


! compute new gas-aerosol mass transfer coefficients
      call aerosolmtc

! condense h2so4, msa, and nh3 only
      call ASTEM_non_volatiles(dtchem)	! analytical solution

! condense inorganic semi-volatile gases hno3, hcl, nh3, and co2
      call ASTEM_semi_volatiles(dtchem)	! semi-implicit + explicit euler

! condense secondary organic gases (8 sorgam species)
      call ASTEM_secondary_organics(dtchem) ! semi-implicit euler

      !PRINT*, 'HERE4-C'
      !PRINT*, Dp_dry_a
      !READ(*,*)
      
      return
      end subroutine ASTEM














!***********************************************************************
! part of ASTEM: integrates semi-volatile inorganic gases
!
! author: Rahul A. Zaveri
! update: jan 2007
!-----------------------------------------------------------------------
      subroutine ASTEM_semi_volatiles(dtchem)
      use module_data_mosaic_main
      use module_data_mosaic_aero

      implicit none

! subr arguments
      real(r8) :: dtchem
! local variables
      integer ibin, iv, jp
      real(r8) :: dtmax, t_new, t_old, t_out, XT
      real(r8) :: sum1, sum2, sum3, sum4, sum4a, sum4b, h_flux_s


! initialize time
      t_old = 0.0
      t_out = dtchem

! reset ASTEM time steps and MESA iterations counters to zero
      isteps_ASTEM = 0
      do ibin = 1, nbin_a
        iter_MESA(ibin) = 0
      enddo

!--------------------------------
! overall integration loop begins over dtchem seconds

10    isteps_ASTEM = isteps_ASTEM + 1

! compute new fluxes
      phi_nh4no3_s = 0.0
      phi_nh4cl_s  = 0.0
      ieqblm_ASTEM = mYES			! reset to default

      do 501 ibin = 1, nbin_a

        idry_case3a(ibin) = mNO			! reset to default
! default fluxes and other stuff
        do iv = 1, ngas_ioa
          sfc_a(iv)                  = gas(iv)
          df_gas_s(iv,ibin)          = 0.0
          df_gas_l(iv,ibin)          = 0.0
          flux_s(iv,ibin)            = 0.0
          flux_l(iv,ibin)            = 0.0
          Heff(iv,ibin)              = 0.0
          volatile_s(iv,ibin)        = 0.0
          phi_volatile_s(iv,ibin)    = 0.0
          phi_volatile_l(iv,ibin)    = 0.0
          integrate(iv,jsolid,ibin)  = mNO	! reset to default
          integrate(iv,jliquid,ibin) = mNO	! reset to default
        enddo


        if(jaerosolstate(ibin) .eq. all_solid)then
          jphase(ibin) = jsolid
          call ASTEM_flux_dry(ibin)
        elseif(jaerosolstate(ibin) .eq. all_liquid)then
          jphase(ibin) = jliquid
          call ASTEM_flux_wet(ibin)
        elseif(jaerosolstate(ibin) .eq. mixed)then
          call ASTEM_flux_mix(ibin)	! jphase(ibin) will be determined in this subr.
        endif

501   continue

      if(ieqblm_ASTEM .eq. mYES)goto 30	! all bins have reached eqblm, so quit.

!-------------------------


! calculate maximum possible internal time-step
11    call ASTEM_calculate_dtmax(dtchem, dtmax)
      t_new = t_old + dtmax	! update time
      if(t_new .gt. t_out)then	! check if the new time step is too large
        dtmax = t_out - t_old
        t_new = t_out*1.01
      endif


!------------------------------------------
! do internal time-step (dtmax) integration

      do 20 iv = 2, 4

        sum1 = 0.0
        sum2 = 0.0
        sum3 = 0.0
        sum4 = 0.0
        sum4a= 0.0
        sum4b= 0.0

        do 21 ibin = 1, nbin_a
          if(jaerosolstate(ibin) .eq. no_aerosol)goto 21

          jp = jliquid
          sum1 = sum1 + aer(iv,jp,ibin)/   &
          (1. + dtmax*kg(iv,ibin)*Heff(iv,ibin)*integrate(iv,jp,ibin))

          sum2 = sum2 + kg(iv,ibin)*integrate(iv,jp,ibin)/   &
          (1. + dtmax*kg(iv,ibin)*Heff(iv,ibin)*integrate(iv,jp,ibin))

          jp = jsolid
          sum3 = sum3 + aer(iv,jp,ibin)

          if(flux_s(iv,ibin) .gt. 0.)then
            h_flux_s = dtmax*flux_s(iv,ibin)
            sum4a = sum4a + h_flux_s
            aer(iv,jp,ibin) = aer(iv,jp,ibin) + h_flux_s
          elseif(flux_s(iv,ibin) .lt. 0.)then
            h_flux_s = min(h_s_i_m(iv,ibin),dtmax)*flux_s(iv,ibin)
            sum4b = sum4b + h_flux_s
            aer(iv,jp,ibin) = aer(iv,jp,ibin) + h_flux_s
            aer(iv,jp,ibin) = max(aer(iv,jp,ibin), 0.0d0)
          endif

21      continue

        sum4 = sum4a + sum4b


! first update gas concentration
        gas(iv) = (total_species(iv) - (sum1 + sum3 + sum4) )/   &
                              (1. + dtmax*sum2)
        gas(iv) = max(gas(iv), 0.0d0)

!        if(gas(iv) .lt. 0.)write(6,*) gas(iv)

! now update aer concentration in the liquid phase
        do 22 ibin = 1, nbin_a

          if(integrate(iv,jliquid,ibin) .eq. mYES)then
            aer(iv,jliquid,ibin) =   &
             (aer(iv,jliquid,ibin) + dtmax*kg(iv,ibin)*gas(iv))/   &
                  (1. + dtmax*kg(iv,ibin)*Heff(iv,ibin))
          endif

22      continue


20    continue
!------------------------------------------
! sub-step integration done


!------------------------------------------
! now update aer(jtotal) and update internal phase equilibrium
! also do integration of species by mass balance if necessary
!
      do 40 ibin = 1, nbin_a
        if(jaerosolstate(ibin) .eq. no_aerosol)goto 40

        if(jphase(ibin) .eq. jsolid)then
          call form_electrolytes(jsolid,ibin,XT)  ! degas excess nh3 (if present)
        elseif(jphase(ibin) .eq. jliquid)then
          call form_electrolytes(jliquid,ibin,XT) ! degas excess nh3 (if present)
        elseif(jphase(ibin) .eq. jtotal)then
          call form_electrolytes(jsolid,ibin,XT)  ! degas excess nh3 (if present)
          call form_electrolytes(jliquid,ibin,XT) ! degas excess nh3 (if present)
        endif

!========================
! now update jtotal
        do iv = 2, ngas_ioa
          aer(iv,jtotal,ibin)=aer(iv,jsolid,ibin)+aer(iv,jliquid,ibin)
        enddo
!========================


        call form_electrolytes(jtotal,ibin,XT)	! for MDRH diagnosis



! update internal phase equilibrium
        if(jhyst_leg(ibin) .eq. jhyst_lo)then
          call ASTEM_update_phase_eqblm(ibin)
        else
          call do_full_deliquescence(ibin)		! simply do liquid <-- total
        endif


40    continue
!------------------------------------------

! update time
      t_old = t_new


      if(isteps_ASTEM .ge. nmax_ASTEM)then
        jASTEM_fail = jASTEM_fail + 1
        write(6,*)'ASTEM internal steps exceeded', nmax_ASTEM
        write(6,*)'ibin =', ibin
        if(iprint_input .eq. mYES)then
!          call print_input
          iprint_input = mNO
        endif
        goto 30
      elseif(t_new .lt. t_out)then
        goto 10
      endif


! check if end of dtchem reached
      if(t_new .lt. 0.9999*t_out) goto 10

30    cumul_steps_ASTEM = cumul_steps_ASTEM + isteps_ASTEM
      isteps_ASTEM_max = max(isteps_ASTEM_max, isteps_ASTEM )
!================================================
! end of overall integration loop over dtchem seconds


!
! call subs to calculate fluxes over mixed-phase particles to update H+ ions,
! which were wiped off during update_phase_eqblm
      do ibin = 1, nbin_a

        if(jaerosolstate(ibin) .eq. mixed)then
          if( electrolyte(jnh4no3,jsolid,ibin).gt. 0.0 .or.   &
              electrolyte(jnh4cl, jsolid,ibin).gt. 0.0 )then
            call ASTEM_flux_mix(ibin)		! jphase(ibin) will be determined in this subr.
          else
            jphase(ibin) = jliquid
            call ASTEM_flux_wet(ibin)
          endif
        endif

      enddo



      return
      end subroutine ASTEM_semi_volatiles












!***********************************************************************
! part of ASTEM: computes max time step for gas-aerosol integration
!
! author: Rahul A. Zaveri
! update: jan 2005
!-----------------------------------------------------------------------
      subroutine ASTEM_calculate_dtmax(dtchem, dtmax)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      real(r8) :: dtchem, dtmax
! local variables
      integer ibin, iv
      real(r8) :: alpha, h_gas, h_sub_max,   &
           h_gas_i(ngas_ioa), h_gas_l, h_gas_s,   &
           sum_kg_phi, sumflux_s


      h_sub_max = 60.0	! sec




! GAS-SIDE

! solid-phase
! calculate h_gas_i and h_gas_l

      h_gas_s = 2.e16

      do 5 iv = 2, ngas_ioa
        h_gas_i(iv) = 1.e16
        sumflux_s = 0.0
        do ibin = 1, nbin_a
          if(flux_s(iv,ibin) .gt. 0.0)then
            sumflux_s = sumflux_s + flux_s(iv,ibin)
          endif
        enddo

        if(sumflux_s .gt. 0.0)then
          h_gas_i(iv) = 0.1*gas(iv)/sumflux_s
          h_gas_s     = min(h_gas_s, h_gas_i(iv))
        endif

5     continue


! liquid-phase
! calculate h_gas_s and h_gas_l

      h_gas_l = 2.e16

      do 6 iv = 2, ngas_ioa
        h_gas_i(iv) = 1.e16
        sum_kg_phi = 0.0
        do ibin = 1, nbin_a
          if(integrate(iv,jliquid,ibin) .eq. mYES)then
          sum_kg_phi = sum_kg_phi +   &
                       abs(phi_volatile_l(iv,ibin))*kg(iv,ibin)
          endif
        enddo

        if(sum_kg_phi .gt. 0.0)then
          h_gas_i(iv) = alpha_astem/sum_kg_phi
          h_gas_l     = min(h_gas_l, h_gas_i(iv))
        endif

6     continue

      h_gas = min(h_gas_s, h_gas_l)
      h_gas = min(h_gas, h_sub_max)




! AEROSOL-SIDE: solid-phase

! first load volatile_solid array
      do ibin = 1, nbin_a

        volatile_s(ino3_a,ibin) = electrolyte(jnh4no3,jsolid,ibin)
        volatile_s(inh4_a,ibin) = electrolyte(jnh4cl,jsolid,ibin) +   &
                                  electrolyte(jnh4no3,jsolid,ibin)

        if(idry_case3a(ibin) .eq. mYES)then
          volatile_s(icl_a,ibin)  = aer(icl_a,jsolid,ibin)
        else
          volatile_s(icl_a,ibin)  = electrolyte(jnh4cl,jsolid,ibin)
        endif

      enddo


! next calculate weighted avg_df_gas_s
      do iv = 2, ngas_ioa

        sum_bin_s(iv) = 0.0
        sum_vdf_s(iv) = 0.0
        sum_vol_s(iv) = 0.0

        do ibin = 1, nbin_a
          if(flux_s(iv,ibin) .lt. 0.)then	! aer -> gas
            sum_bin_s(iv) = sum_bin_s(iv) + 1.0
            sum_vdf_s(iv) = sum_vdf_s(iv) +   &
                            volatile_s(iv,ibin)*df_gas_s(iv,ibin)
            sum_vol_s(iv) = sum_vol_s(iv) + volatile_s(iv,ibin)
          endif
        enddo

        if(sum_vol_s(iv) .gt. 0.0)then
          avg_df_gas_s(iv) = sum_vdf_s(iv)/sum_vol_s(iv)
        else
          avg_df_gas_s(iv) = 1.0 ! never used, but set to 1.0 just to be safe
        endif

      enddo


! calculate h_s_i_m


      do 20 ibin = 1, nbin_a

        if(jaerosolstate(ibin) .eq. no_aerosol) goto 20

        do 10 iv = 2, ngas_ioa

          if(flux_s(iv,ibin) .lt. 0.)then				! aer -> gas

            alpha = abs(avg_df_gas_s(iv))/   &
                   (volatile_s(iv,ibin)*sum_bin_s(iv))
            alpha = min(alpha, 1.0d0)

            if(idry_case3a(ibin) .eq. mYES)alpha = 1.0

            h_s_i_m(iv,ibin) =   &
                 -alpha*volatile_s(iv,ibin)/flux_s(iv,ibin)

          endif

10      continue


20    continue


      dtmax = min(dtchem, h_gas)


      if(dtmax .eq. 0.0)then
        write(6,*)' dtmax = ', dtmax
      endif

      return
      end subroutine ASTEM_calculate_dtmax















!***********************************************************************
! part of ASTEM: updates solid-liquid partitioning after each gas-aerosol
! mass transfer step
!
! author: Rahul A. Zaveri
! update: jan 2005
!-----------------------------------------------------------------------
      subroutine ASTEM_update_phase_eqblm(ibin)	! TOUCH
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      integer jdum, js, j_index, je
      real(r8) :: XT, sum_dum



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


! calculate overall sulfate ratio
      call calculate_XT(ibin,jtotal,XT)		! calc updated XT

! now diagnose MDRH
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

      if(aH2O*100. .lt. MDRH(ibin)) then
        jaerosolstate(ibin) = all_solid
        jphase(ibin) = jsolid
        call adjust_solid_aerosol(ibin)
        return
      endif


! none of the above means it must be sub-saturated or mixed-phase
10    if(jphase(ibin) .eq. jsolid)then
        call do_full_deliquescence(ibin)
        call MESA_PTC(ibin)
      else
        call MESA_PTC(ibin)
      endif



      return
      end subroutine ASTEM_update_phase_eqblm












!==================================================================
!
! LIQUID PARTICLES
!
!***********************************************************************
! part of ASTEM: computes fluxes over wet aerosols
!
! author: Rahul A. Zaveri
! update: Jan 2007
!-----------------------------------------------------------------------
      subroutine ASTEM_flux_wet(ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      integer iv, iadjust, iadjust_intermed
      real(r8) :: XT, g_nh3_hno3, g_nh3_hcl, a_nh4_no3, a_nh4_cl



      call ions_to_electrolytes(jliquid,ibin,XT)  	! for water content calculation
      call compute_activities(ibin)

      if(water_a(ibin) .eq. 0.0)then
	write(6,*)'Water is zero in liquid phase'
	write(6,*)'Stopping in ASTEM_flux_wet'
        stop
      endif

!-------------------------------------------------------------------
! CASE 1: caco3 > 0 absorb acids (and indirectly degas co2)

      if(electrolyte(jcaco3,jsolid,ibin) .gt. 0.0)then
        call ASTEM_flux_wet_case1(ibin)
        return
      endif

!-------------------------------------------------------------------
! CASE 2: Sulfate-Rich Domain

      if(XT.lt.1.9999 .and. XT.ge.0.)then
        call ASTEM_flux_wet_case2(ibin)
        return
      endif

!-------------------------------------------------------------------

      if( (gas(inh3_g)+aer(inh4_a,jliquid,ibin)) .lt. 1.e-25)goto 10  ! no ammonia in the system

!-------------------------------------------------------------------
! CASE 3: nh4no3 and/or nh4cl maybe active
! do some small adjustments (if needed) before deciding case 3

      iadjust = mNO		! default
      iadjust_intermed = mNO	! default

! nh4no3
      g_nh3_hno3 = gas(inh3_g)*gas(ihno3_g)
      a_nh4_no3  = aer(inh4_a,jliquid,ibin)*aer(ino3_a,jliquid,ibin)

      if(g_nh3_hno3 .gt. 0. .and. a_nh4_no3 .eq. 0.)then
        call absorb_tiny_nh4no3(ibin)
        iadjust = mYES
        iadjust_intermed = mYES
      endif

      if(iadjust_intermed .eq. mYES)then
        call ions_to_electrolytes(jliquid,ibin,XT)  	! update after adjustments
        iadjust_intermed = mNO	! reset
      endif

! nh4cl
      g_nh3_hcl = gas(inh3_g)*gas(ihcl_g)
      a_nh4_cl  = aer(inh4_a,jliquid,ibin)*aer(icl_a,jliquid,ibin)

      if(g_nh3_hcl .gt. 0. .and. a_nh4_cl .eq. 0.)then
        call absorb_tiny_nh4cl(ibin)
        iadjust = mYES
        iadjust_intermed = mYES
      endif

      if(iadjust_intermed .eq. mYES)then
        call ions_to_electrolytes(jliquid,ibin,XT)  	! update after adjustments
      endif

      if(iadjust .eq. mYES)then
        call compute_activities(ibin)			! update after adjustments
      endif


! all adjustments done...

!--------
      kelvin_nh4no3 = kel(inh3_g,ibin)*kel(ihno3_g,ibin)
      Keq_nh4no3 = kelvin_nh4no3*activity(jnh4no3,ibin)*Kp_nh4no3	! = [NH3]s * [HNO3]s

      kelvin_nh4cl = kel(inh3_g,ibin)*kel(ihcl_g,ibin)
      Keq_nh4cl = kelvin_nh4cl*activity(jnh4cl,ibin)*Kp_nh4cl	! = [NH3]s * [HCl]s

      call ASTEM_flux_wet_case3(ibin)

      return


!-------------------------------------------------------------------
! CASE 4: ammonia = 0. hno3 and hcl exchange may happen here
! do small adjustments (if needed) before deciding case 4

10    iadjust = mNO		! default
      iadjust_intermed = mNO	! default

! hno3
      if(gas(ihno3_g).gt.0. .and. aer(ino3_a,jliquid,ibin).eq.0. .and.   &
         aer(icl_a,jliquid,ibin) .gt. 0.0)then
        call absorb_tiny_hno3(ibin)	! and degas tiny hcl
        iadjust = mYES
        iadjust_intermed = mYES
      endif

      if(iadjust_intermed .eq. mYES)then
        call ions_to_electrolytes(jliquid,ibin,XT)  	! update after adjustments
        iadjust_intermed = mNO	! reset
      endif

! hcl
      if(gas(ihcl_g).gt.0. .and. aer(icl_a,jliquid,ibin) .eq. 0. .and.   &
         aer(ino3_a,jliquid,ibin) .gt. 0.0)then
        call absorb_tiny_hcl(ibin)	! and degas tiny hno3
        iadjust = mYES
        iadjust_intermed = mYES
      endif

      if(iadjust_intermed .eq. mYES)then
        call ions_to_electrolytes(jliquid,ibin,XT)  	! update after adjustments
      endif

      if(iadjust .eq. mYES)then
        call compute_activities(ibin)			! update after adjustments
      endif

! all adjustments done...

      call ASTEM_flux_wet_case4(ibin)


      return
      end subroutine ASTEM_flux_wet












!***********************************************************************
! part of ASTEM: subroutines for flux_wet cases
!
! author: Rahul A. Zaveri
! update: Jan 2007
!-----------------------------------------------------------------------

! CASE 1: CaCO3 > 0 absorb all acids (and indirectly degas co2)

      subroutine ASTEM_flux_wet_case1(ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      integer iv

      mc(jc_h,ibin) = sqrt(Keq_ll(3))

! same as dry case1
      if(gas(ihno3_g) .gt. 1.e-6)then
        sfc_a(ihno3_g) = 0.0
        df_gas_s(ihno3_g,ibin) = gas(ihno3_g)
        phi_volatile_s(ihno3_g,ibin) = 1.0
        flux_s(ihno3_g,ibin) = kg(ihno3_g,ibin)*df_gas_s(ihno3_g,ibin)
        integrate(ihno3_g,jsolid,ibin) = mYES
        jphase(ibin) = jsolid
        ieqblm_ASTEM = mNO
      endif

      if(gas(ihcl_g) .gt. 1.e-6)then
        sfc_a(ihcl_g)  = 0.0
        df_gas_s(ihcl_g,ibin) = gas(ihcl_g)
        phi_volatile_s(ihcl_g,ibin) = 1.0
        flux_s(ihcl_g,ibin) = kg(ihcl_g,ibin)*df_gas_s(ihcl_g,ibin)
        integrate(ihcl_g,jsolid,ibin)  = mYES
        jphase(ibin) = jsolid
        ieqblm_ASTEM = mNO
      endif

      return
      end subroutine ASTEM_flux_wet_case1



!--------------------------------------------------------------------
! CASE 2: Sulfate-Rich Domain

      subroutine ASTEM_flux_wet_case2(ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      real(r8) :: dum_hno3, dum_hcl, dum_nh3


      sfc_a(inh3_g)  = kel(inh3_g,ibin)*   &
                       gam_ratio(ibin)*mc(jc_nh4,ibin)*Keq_ll(3)/   &
                        (mc(jc_h,ibin)*Keq_ll(2)*Keq_gl(2))

      sfc_a(ihno3_g) = kel(ihno3_g,ibin)*   &
                   mc(jc_h,ibin)*ma(ja_no3,ibin)*gam(jhno3,ibin)**2/   &
                   Keq_gl(3)

      sfc_a(ihcl_g)  = kel(ihcl_g,ibin)*   &
                   mc(jc_h,ibin)*ma(ja_cl,ibin)*gam(jhcl,ibin)**2/   &
                   Keq_gl(4)

      dum_hno3 = max(sfc_a(ihno3_g), gas(ihno3_g))
      dum_hcl  = max(sfc_a(ihcl_g), gas(ihcl_g))
      dum_nh3  = max(sfc_a(inh3_g), gas(inh3_g))


! compute relative driving forces
      if(dum_hno3 .gt. 0.0)then
        df_gas_l(ihno3_g,ibin) = gas(ihno3_g) - sfc_a(ihno3_g)
        phi_volatile_l(ihno3_g,ibin)= df_gas_l(ihno3_g,ibin)/dum_hno3
      else
        phi_volatile_l(ihno3_g,ibin)= 0.0
      endif

      if(dum_hcl .gt. 0.0)then
        df_gas_l(ihcl_g,ibin)  = gas(ihcl_g)  - sfc_a(ihcl_g)
        phi_volatile_l(ihcl_g,ibin) = df_gas_l(ihcl_g,ibin)/dum_hcl
      else
        phi_volatile_l(ihcl_g,ibin) = 0.0
      endif

      if(dum_nh3 .gt. 0.0)then
        df_gas_l(inh3_g,ibin)  = gas(inh3_g)  - sfc_a(inh3_g)
        phi_volatile_l(inh3_g,ibin) = df_gas_l(inh3_g,ibin)/dum_nh3
      else
        phi_volatile_l(inh3_g,ibin) = 0.0
      endif


!      if(phi_volatile_l(ihno3_g,ibin) .le. rtol_eqb_astem .and.
!     &   phi_volatile_l(ihcl_g,ibin)  .le. rtol_eqb_astem .and.
!     &   phi_volatile_l(inh3_g,ibin)  .le. rtol_eqb_astem)then
!
!        return
!
!      endif


! compute Heff
      if(dum_hno3 .gt. 0.0)then
        Heff(ihno3_g,ibin)=   &
          kel(ihno3_g,ibin)*gam(jhno3,ibin)**2*mc(jc_h,ibin)*1.e-9/   &
                       (water_a(ibin)*Keq_gl(3))
        integrate(ihno3_g,jliquid,ibin)= mYES
        ieqblm_ASTEM = mNO
      endif

      if(dum_hcl .gt. 0.0)then
        Heff(ihcl_g,ibin)=   &
          kel(ihcl_g,ibin)*gam(jhcl,ibin)**2*mc(jc_h,ibin)*1.e-9/   &
                       (water_a(ibin)*Keq_gl(4))
        integrate(ihcl_g,jliquid,ibin) = mYES
        ieqblm_ASTEM = mNO
      endif

      if(dum_nh3 .gt. 0.0)then
        Heff(inh3_g,ibin) =   &
             kel(inh3_g,ibin)*gam_ratio(ibin)*1.e-9*Keq_ll(3)/   &
             (water_a(ibin)*mc(jc_h,ibin)*Keq_ll(2)*Keq_gl(2))
        integrate(inh3_g,jliquid,ibin) = mYES
        ieqblm_ASTEM = mNO
      endif


      return
      end subroutine ASTEM_flux_wet_case2








!---------------------------------------------------------------------
! CASE 3: nh4no3 and/or nh4cl may be active

      subroutine ASTEM_flux_wet_case3(ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      real(r8) :: a, b, c, dum_hno3, dum_hcl, dum_nh3
! function
      real(r8) :: quadratic

      a =   kg(inh3_g,ibin)
      b = - kg(inh3_g,ibin)*gas(inh3_g)   &
          + kg(ihno3_g,ibin)*gas(ihno3_g)   &
          + kg(ihcl_g,ibin)*gas(ihcl_g)
      c = -(kg(ihno3_g,ibin)*Keq_nh4no3 + kg(ihcl_g,ibin)*Keq_nh4cl)

      sfc_a(inh3_g)  = quadratic(a,b,c)
      sfc_a(ihno3_g) = Keq_nh4no3/max(sfc_a(inh3_g),1.d-20)
      sfc_a(ihcl_g)  = Keq_nh4cl/max(sfc_a(inh3_g),1.d-20)


! diagnose mH+
      if(sfc_a(ihno3_g).gt.0.0 .and. ma(ja_no3,ibin).gt.0.0)then
        mc(jc_h,ibin) = Keq_gl(3)*sfc_a(ihno3_g)/   &
        (kel(ihno3_g,ibin)*gam(jhno3,ibin)**2 * ma(ja_no3,ibin))
      elseif(sfc_a(ihcl_g).gt.0.0 .and. ma(ja_cl,ibin).gt.0.0)then
        mc(jc_h,ibin) = Keq_gl(4)*sfc_a(ihcl_g)/   &
        (kel(ihcl_g,ibin)*gam(jhcl,ibin)**2 * ma(ja_cl,ibin))
      else
        call equilibrate_acids(ibin)	! hno3 and/or hcl may be > 0 in the gas phase
        mc(jc_h,ibin)  = max(mc(jc_h,ibin), sqrt(Keq_ll(3)))

        sfc_a(inh3_g)  = kel(inh3_g,ibin)*   &
                         gam_ratio(ibin)*mc(jc_nh4,ibin)*Keq_ll(3)/   &
                        (mc(jc_h,ibin)*Keq_ll(2)*Keq_gl(2))

        sfc_a(ihno3_g) = kel(ihno3_g,ibin)*   &
                   mc(jc_h,ibin)*ma(ja_no3,ibin)*gam(jhno3,ibin)**2/   &
                   Keq_gl(3)
        sfc_a(ihcl_g)  = kel(ihcl_g,ibin)*   &
                   mc(jc_h,ibin)*ma(ja_cl,ibin)*gam(jhcl,ibin)**2/   &
                   Keq_gl(4)
      endif



      dum_hno3 = max(sfc_a(ihno3_g), gas(ihno3_g))
      dum_hcl  = max(sfc_a(ihcl_g), gas(ihcl_g))
      dum_nh3  = max(sfc_a(inh3_g), gas(inh3_g))

! compute relative driving forces
      if(dum_hno3 .gt. 0.0)then
        df_gas_l(ihno3_g,ibin) = gas(ihno3_g) - sfc_a(ihno3_g)
        phi_volatile_l(ihno3_g,ibin)= df_gas_l(ihno3_g,ibin)/dum_hno3
      else
        phi_volatile_l(ihno3_g,ibin)= 0.0
      endif

      if(dum_hcl .gt. 0.0)then
        df_gas_l(ihcl_g,ibin)  = gas(ihcl_g)  - sfc_a(ihcl_g)
        phi_volatile_l(ihcl_g,ibin) = df_gas_l(ihcl_g,ibin)/dum_hcl
      else
        phi_volatile_l(ihcl_g,ibin) = 0.0
      endif

      if(dum_nh3 .gt. 0.0)then
        df_gas_l(inh3_g,ibin)  = gas(inh3_g)  - sfc_a(inh3_g)
        phi_volatile_l(inh3_g,ibin) = df_gas_l(inh3_g,ibin)/dum_nh3
      else
        phi_volatile_l(inh3_g,ibin) = 0.0
      endif



!      if(phi_volatile_l(ihno3_g,ibin) .le. rtol_eqb_astem .and.
!     &   phi_volatile_l(ihcl_g,ibin)  .le. rtol_eqb_astem .and.
!     &   phi_volatile_l(inh3_g,ibin)  .le. rtol_eqb_astem)then
!
!        return
!
!      endif


! compute Heff
      if(dum_hno3 .gt. 0.0)then
        Heff(ihno3_g,ibin)=   &
          kel(ihno3_g,ibin)*gam(jhno3,ibin)**2*mc(jc_h,ibin)*1.e-9/   &
                       (water_a(ibin)*Keq_gl(3))
        integrate(ihno3_g,jliquid,ibin)= mYES
        ieqblm_ASTEM = mNO
      endif

      if(dum_hcl .gt. 0.0)then
        Heff(ihcl_g,ibin)=   &
          kel(ihcl_g,ibin)*gam(jhcl,ibin)**2*mc(jc_h,ibin)*1.e-9/   &
                       (water_a(ibin)*Keq_gl(4))
        integrate(ihcl_g,jliquid,ibin) = mYES
        ieqblm_ASTEM = mNO
      endif

      if(dum_nh3 .gt. 0.0)then
        Heff(inh3_g,ibin) =   &
             kel(inh3_g,ibin)*gam_ratio(ibin)*1.e-9*Keq_ll(3)/   &
             (water_a(ibin)*mc(jc_h,ibin)*Keq_ll(2)*Keq_gl(2))
        integrate(inh3_g,jliquid,ibin) = mYES
        ieqblm_ASTEM = mNO
      endif



      return
      end subroutine ASTEM_flux_wet_case3









!--------------------------------------------------------------------
! CASE 3a: only NH4NO3 (aq) active

      subroutine ASTEM_flux_wet_case3a(ibin)	! NH4NO3 (aq)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      real(r8) :: a, b, c, dum_hno3, dum_nh3
! function
      real(r8) :: quadratic


      a =   kg(inh3_g,ibin)
      b = - kg(inh3_g,ibin)*gas(inh3_g)   &
          + kg(ihno3_g,ibin)*gas(ihno3_g)
      c = -(kg(ihno3_g,ibin)*Keq_nh4no3)

      sfc_a(inh3_g)  = quadratic(a,b,c)
      sfc_a(ihno3_g) = Keq_nh4no3/sfc_a(inh3_g)


! diagnose mH+
      if(sfc_a(ihno3_g).gt.0.0 .and. ma(ja_no3,ibin).gt.0.0)then
        mc(jc_h,ibin) = Keq_gl(3)*sfc_a(ihno3_g)/   &
          (kel(ihno3_g,ibin)*gam(jhno3,ibin)**2 * ma(ja_no3,ibin))
      else
        mc(jc_h,ibin) = sqrt(Keq_ll(3))
      endif


! compute Heff
      dum_hno3 = max(sfc_a(ihno3_g), gas(ihno3_g))
      dum_nh3  = max(sfc_a(inh3_g), gas(inh3_g))

! compute relative driving forces
      if(dum_hno3 .gt. 0.0)then
        df_gas_l(ihno3_g,ibin) = gas(ihno3_g) - sfc_a(ihno3_g)
        phi_volatile_l(ihno3_g,ibin)= df_gas_l(ihno3_g,ibin)/dum_hno3
      else
        phi_volatile_l(ihno3_g,ibin)= 0.0
      endif

      if(dum_nh3 .gt. 0.0)then
        df_gas_l(inh3_g,ibin)  = gas(inh3_g)  - sfc_a(inh3_g)
        phi_volatile_l(inh3_g,ibin) = df_gas_l(inh3_g,ibin)/dum_nh3
      else
        phi_volatile_l(inh3_g,ibin) = 0.0
      endif


!      if(phi_volatile_l(ihno3_g,ibin) .le. rtol_eqb_astem .and.
!     &   phi_volatile_l(inh3_g,ibin)  .le. rtol_eqb_astem)then
!
!        return
!
!      endif


! compute Heff
      Heff(ihno3_g,ibin)=   &
        kel(ihno3_g,ibin)*gam(jhno3,ibin)**2*mc(jc_h,ibin)*1.e-9/   &
                     (water_a(ibin)*Keq_gl(3))
      integrate(ihno3_g,jliquid,ibin)= mYES


      Heff(inh3_g,ibin) =   &
           kel(inh3_g,ibin)*gam_ratio(ibin)*1.e-9*Keq_ll(3)/   &
           (water_a(ibin)*mc(jc_h,ibin)*Keq_ll(2)*Keq_gl(2))
      integrate(inh3_g,jliquid,ibin) = mYES


      ieqblm_ASTEM = mNO


      return
      end subroutine ASTEM_flux_wet_case3a









!--------------------------------------------------------------------
! CASE 3b: only NH4Cl (aq) active

      subroutine ASTEM_flux_wet_case3b(ibin)	! NH4Cl (aq)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      real(r8) :: a, b, c, dum_hcl, dum_nh3
! function
      real(r8) :: quadratic


      a =   kg(inh3_g,ibin)
      b = - kg(inh3_g,ibin)*gas(inh3_g)   &
          + kg(ihcl_g,ibin)*gas(ihcl_g)
      c = -(kg(ihcl_g,ibin)*Keq_nh4cl)

      sfc_a(inh3_g)  = quadratic(a,b,c)
      sfc_a(ihcl_g)  = Keq_nh4cl /sfc_a(inh3_g)


! diagnose mH+
      if(sfc_a(ihcl_g).gt.0.0 .and. ma(ja_cl,ibin).gt.0.0)then
        mc(jc_h,ibin) = Keq_gl(4)*sfc_a(ihcl_g)/   &
          (kel(ihcl_g,ibin)*gam(jhcl,ibin)**2 * ma(ja_cl,ibin))
      else
        mc(jc_h,ibin) = sqrt(Keq_ll(3))
      endif


! compute Heff
      dum_hcl  = max(sfc_a(ihcl_g), gas(ihcl_g))
      dum_nh3  = max(sfc_a(inh3_g), gas(inh3_g))


! compute relative driving forces
      if(dum_hcl .gt. 0.0)then
        df_gas_l(ihcl_g,ibin)  = gas(ihcl_g)  - sfc_a(ihcl_g)
        phi_volatile_l(ihcl_g,ibin) = df_gas_l(ihcl_g,ibin)/dum_hcl
      else
        phi_volatile_l(ihcl_g,ibin) = 0.0
      endif

      if(dum_nh3 .gt. 0.0)then
        df_gas_l(inh3_g,ibin)  = gas(inh3_g)  - sfc_a(inh3_g)
        phi_volatile_l(inh3_g,ibin) = df_gas_l(inh3_g,ibin)/dum_nh3
      else
        phi_volatile_l(inh3_g,ibin) = 0.0
      endif



!      if(phi_volatile_l(ihcl_g,ibin)  .le. rtol_eqb_astem .and.
!     &   phi_volatile_l(inh3_g,ibin)  .le. rtol_eqb_astem)then
!
!        return
!
!      endif



! compute Heff
      Heff(ihcl_g,ibin)=   &
          kel(ihcl_g,ibin)*gam(jhcl,ibin)**2*mc(jc_h,ibin)*1.e-9/   &
                       (water_a(ibin)*Keq_gl(4))
      integrate(ihcl_g,jliquid,ibin) = mYES


      Heff(inh3_g,ibin) =   &
             kel(inh3_g,ibin)*gam_ratio(ibin)*1.e-9*Keq_ll(3)/   &
             (water_a(ibin)*mc(jc_h,ibin)*Keq_ll(2)*Keq_gl(2))
      integrate(inh3_g,jliquid,ibin) = mYES


      ieqblm_ASTEM = mNO



      return
      end subroutine ASTEM_flux_wet_case3b









!-----------------------------------------------------------------------
! CASE 4: NH3 = 0 (in gas and aerosol). hno3 and hcl exchange may happen here

      subroutine ASTEM_flux_wet_case4(ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      real(r8) :: dum_numer, dum_denom, gas_eqb_ratio, dum_hno3, dum_hcl


      dum_numer = kel(ihno3_g,ibin)*Keq_gl(4)*ma(ja_no3,ibin)*   &
                  gam(jhno3,ibin)**2
      dum_denom = kel(ihcl_g,ibin)*Keq_gl(3)*ma(ja_cl ,ibin)*   &
                  gam(jhcl,ibin)**2


      if(dum_denom .eq. 0.0 .or. dum_numer .eq. 0.0)then
        mc(jc_h,ibin) = sqrt(Keq_ll(3))
        return
      endif

      gas_eqb_ratio = dum_numer/dum_denom	! Ce,hno3/Ce,hcl


! compute equilibrium surface concentrations
      sfc_a(ihcl_g) =   &
       ( kg(ihno3_g,ibin)*gas(ihno3_g) + kg(ihcl_g,ibin)*gas(ihcl_g) )/   &
           ( kg(ihcl_g,ibin) + gas_eqb_ratio*kg(ihno3_g,ibin) )
      sfc_a(ihno3_g)= gas_eqb_ratio*sfc_a(ihcl_g)


! diagnose mH+
      if(sfc_a(ihno3_g).gt.0.0 .and. ma(ja_no3,ibin).gt.0.0)then
        mc(jc_h,ibin) = Keq_gl(3)*sfc_a(ihno3_g)/   &
        (kel(ihno3_g,ibin)*gam(jhno3,ibin)**2 * ma(ja_no3,ibin))
      elseif(sfc_a(ihcl_g).gt.0.0 .and. ma(ja_cl,ibin).gt.0.0)then
        mc(jc_h,ibin) = Keq_gl(4)*sfc_a(ihcl_g)/   &
        (kel(ihcl_g,ibin)*gam(jhcl,ibin)**2 * ma(ja_cl,ibin))
      else
        mc(jc_h,ibin) = sqrt(Keq_ll(3))
      endif


! compute Heff
      dum_hno3 = max(sfc_a(ihno3_g), gas(ihno3_g))
      dum_hcl  = max(sfc_a(ihcl_g), gas(ihcl_g))

! compute relative driving forces
      if(dum_hno3 .gt. 0.0)then
        df_gas_l(ihno3_g,ibin) = gas(ihno3_g) - sfc_a(ihno3_g)
        phi_volatile_l(ihno3_g,ibin)= df_gas_l(ihno3_g,ibin)/dum_hno3
      else
        phi_volatile_l(ihno3_g,ibin)= 0.0
      endif

      if(dum_hcl .gt. 0.0)then
        df_gas_l(ihcl_g,ibin)  = gas(ihcl_g)  - sfc_a(ihcl_g)
        phi_volatile_l(ihcl_g,ibin)= df_gas_l(ihcl_g,ibin)/dum_hcl
      else
        phi_volatile_l(ihcl_g,ibin)= 0.0
      endif


!      if(phi_volatile_l(ihno3_g,ibin) .le. rtol_eqb_astem .and.
!     &   phi_volatile_l(ihcl_g,ibin)  .le. rtol_eqb_astem)then
!
!        return
!
!      endif



! compute Heff
      Heff(ihno3_g,ibin)=   &
          kel(ihno3_g,ibin)*gam(jhno3,ibin)**2*mc(jc_h,ibin)*1.e-9/   &
                       (water_a(ibin)*Keq_gl(3))
      integrate(ihno3_g,jliquid,ibin)= mYES


      Heff(ihcl_g,ibin)=   &
          kel(ihcl_g,ibin)*gam(jhcl,ibin)**2*mc(jc_h,ibin)*1.e-9/   &
                       (water_a(ibin)*Keq_gl(4))
      integrate(ihcl_g,jliquid,ibin) = mYES


      ieqblm_ASTEM = mNO



      return
      end subroutine ASTEM_flux_wet_case4














!===========================================================
!
! DRY PARTICLES
!
!===========================================================
!***********************************************************************
! part of ASTEM: computes gas-aerosol fluxes over dry aerosols
!
! author: Rahul A. Zaveri
! update: dec 2006
!-----------------------------------------------------------------------
      subroutine ASTEM_flux_dry(ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      integer iv
      real(r8) :: XT, prod_nh4no3, prod_nh4cl, volatile_cl




      call calculate_XT(ibin,jsolid,XT)

!-----------------------------------------------------------------
! CASE 1:  caco3 > 0 absorb all acids (and indirectly degas co2)

      if(electrolyte(jcaco3,jsolid,ibin) .gt. 0.0)then

        call ASTEM_flux_dry_case1(ibin)

        return
      endif

!-----------------------------------------------------------------
! CASE 2: Sulfate-Rich Domain

      if(XT.lt.1.9999 .and. XT.ge.0.)then	! excess sulfate (acidic)

	call ASTEM_flux_dry_case2(ibin)

        return
      endif

!-------------------------------------------------------------------
! CASE 3: hno3 and hcl exchange may happen here and nh4cl may form/evaporate

      volatile_cl  = electrolyte(jnacl,jsolid,ibin) +   &
                     electrolyte(jcacl2,jsolid,ibin)


      if(volatile_cl .gt. 0.0 .and. gas(ihno3_g).gt. 0.0 )then

        call ASTEM_flux_dry_case3a(ibin)

        Keq_nh4cl_0  = min(Kp_nh4cl_0,  Keq_sg(2))	! 6/25/2008

        prod_nh4cl = max( (gas(inh3_g)*gas(ihcl_g)-Keq_nh4cl_0), 0.0d0) +   &
                     electrolyte(jnh4cl, jsolid,ibin)	! 6/25/2008

        if(prod_nh4cl .gt. 0.0)then
          call ASTEM_flux_dry_case3b(ibin)
        endif


        return
      endif

!-----------------------------------------------------------------
! CASE 4: nh4no3 or nh4cl or both may be active

      Keq_nh4no3_0 = min(Kp_nh4no3_0, Keq_sg(1))	! 6/25/2008
      Keq_nh4cl_0  = min(Kp_nh4cl_0,  Keq_sg(2))	! 6/25/2008

      prod_nh4no3 = max( (gas(inh3_g)*gas(ihno3_g)-Keq_nh4no3_0), 0.0d0) +   &
                    electrolyte(jnh4no3,jsolid,ibin)	! 6/25/2008
      prod_nh4cl  = max( (gas(inh3_g)*gas(ihcl_g) -Keq_nh4cl_0), 0.0d0) +   &
                    electrolyte(jnh4cl, jsolid,ibin)	! 6/25/2008

      if(prod_nh4no3 .gt. 0.0 .or. prod_nh4cl .gt. 0.0)then
        call ASTEM_flux_dry_case4(ibin)
        return
      endif

!-----------------------------------------------------------------

      return
      end subroutine ASTEM_flux_dry

!----------------------------------------------------------------------













!***********************************************************************
! part of ASTEM: subroutines for flux_dry cases
!
! author: Rahul A. Zaveri
! update: dec 2006
!-----------------------------------------------------------------------

! CASE 1:  caco3 > 0 absorb all acids (and indirectly degas co2)

      subroutine ASTEM_flux_dry_case1(ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin


      if(gas(ihno3_g) .gt. 1.e-6)then
        sfc_a(ihno3_g) = 0.0
        df_gas_s(ihno3_g,ibin) = gas(ihno3_g)
        phi_volatile_s(ihno3_g,ibin) = 1.0
        flux_s(ihno3_g,ibin) = kg(ihno3_g,ibin)*df_gas_s(ihno3_g,ibin)
        integrate(ihno3_g,jsolid,ibin) = mYES
        ieqblm_ASTEM = mNO
      endif

      if(gas(ihcl_g) .gt. 1.e-6)then
        sfc_a(ihcl_g)  = 0.0
        df_gas_s(ihcl_g,ibin) = gas(ihcl_g)
        phi_volatile_s(ihcl_g,ibin) = 1.0
        flux_s(ihcl_g,ibin)  = kg(ihcl_g,ibin)*df_gas_s(ihcl_g,ibin)
        integrate(ihcl_g,jsolid,ibin)  = mYES
        ieqblm_ASTEM = mNO
      endif


      return
      end subroutine ASTEM_flux_dry_case1



!---------------------------------------------------------------------
! CASE 2: Sulfate-Rich Domain

      subroutine ASTEM_flux_dry_case2(ibin) ! TOUCH
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin


      if(gas(inh3_g).gt.1.e-6)then
        sfc_a(inh3_g) = 0.0
        df_gas_s(inh3_g,ibin) = gas(inh3_g)
        phi_volatile_s(inh3_g,ibin)  = 1.0
        flux_s(inh3_g,ibin) = kg(inh3_g,ibin)*gas(inh3_g)
        integrate(inh3_g,jsolid,ibin) = mYES
        ieqblm_ASTEM = mNO
      endif


      return
      end subroutine ASTEM_flux_dry_case2




!---------------------------------------------------------------------
! CASE 3a: degas hcl from nacl or cacl2 by flux_s balance with hno3

      subroutine ASTEM_flux_dry_case3a(ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin


      if(gas(ihno3_g) .gt. 1.e-6)then
        sfc_a(ihno3_g) = 0.0
        sfc_a(ihcl_g)  = gas(ihcl_g) + aer(icl_a,jsolid,ibin)

        df_gas_s(ihno3_g,ibin) = gas(ihno3_g)
        df_gas_s(ihcl_g,ibin)  = -aer(icl_a,jsolid,ibin)

        flux_s(ihno3_g,ibin) = kg(ihno3_g,ibin)*gas(ihno3_g)
        flux_s(ihcl_g,ibin)  = -flux_s(ihno3_g,ibin)

        phi_volatile_s(ihno3_g,ibin) = 1.0
        phi_volatile_s(ihcl_g,ibin)=df_gas_s(ihcl_g,ibin)/sfc_a(ihcl_g)

        integrate(ihno3_g,jsolid,ibin) = mYES
        integrate(ihcl_g,jsolid,ibin)  = mYES

        idry_case3a(ibin) = mYES
        ieqblm_ASTEM = mNO
      endif

      return
      end subroutine ASTEM_flux_dry_case3a




!---------------------------------------------------------------------
! CASE 3b: nh4cl may form/evaporate here

      subroutine ASTEM_flux_dry_case3b(ibin)	! TOUCH
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      integer iactive_nh4cl, js
      real(r8) :: a, b, c
      real(r8) :: sum_dum
! function
      real(r8) :: quadratic



!! EFFI calculate percent composition
      sum_dum = 0.0
      do js = 1, nsalt
        sum_dum = sum_dum + electrolyte(js,jsolid,ibin)
      enddo

      if(sum_dum .eq. 0.)sum_dum = 1.0

      epercent(jnh4cl,jsolid,ibin) = 100.*electrolyte(jnh4cl,jsolid,ibin)/sum_dum
!! EFFI



!-------------------
! set default values for flags
      iactive_nh4cl  = 1


! compute relative driving force
      phi_nh4cl_s = (gas(inh3_g)*gas(ihcl_g) - Keq_sg(2))/   &
                    max(gas(inh3_g)*gas(ihcl_g),Keq_sg(2))


!-------------------
! now determine if nh4cl is active or significant
! nh4cl
      if( abs(phi_nh4cl_s) .lt. rtol_eqb_ASTEM )then
        iactive_nh4cl = 0
      elseif(gas(inh3_g)*gas(ihcl_g) .lt. Keq_sg(2) .and.   &
             epercent(jnh4cl, jsolid,ibin) .le. ptol_mol_ASTEM)then
        iactive_nh4cl = 0
        if(epercent(jnh4cl, jsolid,ibin) .gt. 0.0)then
          call degas_solid_nh4cl(ibin)
        endif
      endif


! check the outcome
      if(iactive_nh4cl .eq. 0)return


!-----------------
! nh4cl is active


      a =   kg(inh3_g,ibin)
      b = - kg(inh3_g,ibin)*gas(inh3_g)   &
          + kg(ihcl_g,ibin)*gas(ihcl_g)
      c = -(kg(ihcl_g,ibin)*Keq_sg(2))

      sfc_a(inh3_g) = quadratic(a,b,c)
      sfc_a(ihcl_g) = Keq_sg(2)/sfc_a(inh3_g)

      df_gas_s(ihcl_g,ibin) = gas(ihcl_g) - sfc_a(ihcl_g)
      df_gas_s(inh3_g,ibin) = gas(inh3_g) - sfc_a(inh3_g)

      flux_s(inh3_g,ibin) = kg(inh3_g,ibin)*df_gas_s(inh3_g,ibin)
      flux_s(ihcl_g,ibin) = flux_s(ihcl_g,ibin) + flux_s(inh3_g,ibin)

      phi_volatile_s(inh3_g,ibin) = phi_nh4cl_s

      if(flux_s(ihcl_g,ibin) .gt. 0.0)then
        df_gas_s(ihcl_g,ibin) = flux_s(ihcl_g,ibin)/kg(ihcl_g,ibin)	! recompute df_gas
        phi_volatile_s(ihcl_g,ibin) = phi_nh4cl_s
      else
        sfc_a(ihcl_g)  = gas(ihcl_g) + aer(icl_a,jsolid,ibin)
        df_gas_s(ihcl_g,ibin) = -aer(icl_a,jsolid,ibin)
        phi_volatile_s(ihcl_g,ibin)=df_gas_s(ihcl_g,ibin)/sfc_a(ihcl_g)  ! not to be used
      endif

      integrate(inh3_g,jsolid,ibin) = mYES
      integrate(ihcl_g,jsolid,ibin) = mYES	! integrate HCl with explicit euler

      ieqblm_ASTEM = mNO

      return
      end subroutine ASTEM_flux_dry_case3b




!---------------------------------------------------------------------
! Case 4: NH4NO3 and/or NH4Cl may be active

      subroutine ASTEM_flux_dry_case4(ibin)	! TOUCH
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      integer iactive_nh4no3, iactive_nh4cl, iactive, js
      real(r8) :: a, b, c
      real(r8) :: sum_dum
! function
      real(r8) :: quadratic



!! EFFI calculate percent composition
      sum_dum = 0.0
      do js = 1, nsalt
        sum_dum = sum_dum + electrolyte(js,jsolid,ibin)
      enddo

      if(sum_dum .eq. 0.)sum_dum = 1.0

      epercent(jnh4no3,jsolid,ibin) = 100.*electrolyte(jnh4no3,jsolid,ibin)/sum_dum
      epercent(jnh4cl, jsolid,ibin) = 100.*electrolyte(jnh4cl, jsolid,ibin)/sum_dum
!! EFFI


!-------------------
! set default values for flags
      iactive_nh4no3 = 1
      iactive_nh4cl  = 2


! compute diagnostic products and ratios
      phi_nh4no3_s = (gas(inh3_g)*gas(ihno3_g) - Keq_sg(1))/   &
                     max(gas(inh3_g)*gas(ihno3_g),Keq_sg(1))
      phi_nh4cl_s  = (gas(inh3_g)*gas(ihcl_g) - Keq_sg(2))/   &
                     max(gas(inh3_g)*gas(ihcl_g),Keq_sg(2))


!-------------------
! now determine if nh4no3 and/or nh4cl are active or significant

! nh4no3
      if( abs(phi_nh4no3_s) .lt. rtol_eqb_ASTEM )then
        iactive_nh4no3 = 0
      elseif(gas(inh3_g)*gas(ihno3_g) .lt. Keq_sg(1) .and.   &
             epercent(jnh4no3,jsolid,ibin) .le. ptol_mol_ASTEM)then
        iactive_nh4no3 = 0
        if(epercent(jnh4no3,jsolid,ibin) .gt. 0.0)then
          call degas_solid_nh4no3(ibin)
        endif
      endif

! nh4cl
      if( abs(phi_nh4cl_s) .lt. rtol_eqb_ASTEM )then
        iactive_nh4cl = 0
      elseif(gas(inh3_g)*gas(ihcl_g) .lt. Keq_sg(2) .and.   &
             epercent(jnh4cl, jsolid,ibin) .le. ptol_mol_ASTEM)then
        iactive_nh4cl = 0
        if(epercent(jnh4cl, jsolid,ibin) .gt. 0.0)then
          call degas_solid_nh4cl(ibin)
        endif
      endif


      iactive = iactive_nh4no3 + iactive_nh4cl

! check the outcome
      if(iactive .eq. 0)return


      goto (1,2,3),iactive

!---------------------------------
! only nh4no3 solid is active
1     call ASTEM_flux_dry_case4a(ibin)

      return


!-----------------
! only nh4cl solid is active
2     call ASTEM_flux_dry_case4b(ibin)

      return


!-----------------
! both nh4no3 and nh4cl are active
3     call ASTEM_flux_dry_case4ab(ibin)




      return
      end subroutine ASTEM_flux_dry_case4







!---------------------------------------------------------------------
! Case 4a

      subroutine ASTEM_flux_dry_case4a(ibin) ! NH4NO3 solid
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      real(r8) :: a, b, c
! function
      real(r8) :: quadratic



      a =   kg(inh3_g,ibin)
      b = - kg(inh3_g,ibin)*gas(inh3_g)   &
          + kg(ihno3_g,ibin)*gas(ihno3_g)
      c = -(kg(ihno3_g,ibin)*Keq_nh4no3_0)	! 6/25/2008

      sfc_a(inh3_g)  = quadratic(a,b,c)
      sfc_a(ihno3_g) = Keq_nh4no3_0/sfc_a(inh3_g)	! 6/25/2008

      integrate(ihno3_g,jsolid,ibin) = mYES
      integrate(inh3_g,jsolid,ibin)  = mYES

      df_gas_s(ihno3_g,ibin)=gas(ihno3_g)-sfc_a(ihno3_g)
      df_gas_s(inh3_g,ibin) =gas(inh3_g) -sfc_a(inh3_g)

      phi_volatile_s(ihno3_g,ibin)= phi_nh4no3_s
      phi_volatile_s(inh3_g,ibin) = phi_nh4no3_s

      flux_s(ihno3_g,ibin) = kg(ihno3_g,ibin)*df_gas_s(ihno3_g,ibin)
      flux_s(inh3_g,ibin)  = flux_s(ihno3_g,ibin)

      ieqblm_ASTEM = mNO

      return
      end subroutine ASTEM_flux_dry_case4a




!----------------------------------------------------------------
! Case 4b

      subroutine ASTEM_flux_dry_case4b(ibin) ! NH4Cl solid
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      real(r8) :: a, b, c
! function
      real(r8) :: quadratic


      a =   kg(inh3_g,ibin)
      b = - kg(inh3_g,ibin)*gas(inh3_g)   &
          + kg(ihcl_g,ibin)*gas(ihcl_g)
      c = -(kg(ihcl_g,ibin)*Keq_nh4cl_0)	! 6/25/2008

      sfc_a(inh3_g) = quadratic(a,b,c)
      sfc_a(ihcl_g) = Keq_nh4cl_0 /sfc_a(inh3_g)	! 6/25/2008

      integrate(ihcl_g,jsolid,ibin) = mYES
      integrate(inh3_g,jsolid,ibin) = mYES

      df_gas_s(ihcl_g,ibin) = gas(ihcl_g)-sfc_a(ihcl_g)
      df_gas_s(inh3_g,ibin) = gas(inh3_g)-sfc_a(inh3_g)

      phi_volatile_s(ihcl_g,ibin) = phi_nh4cl_s
      phi_volatile_s(inh3_g,ibin) = phi_nh4cl_s

      flux_s(ihcl_g,ibin) = kg(ihcl_g,ibin)*df_gas_s(ihcl_g,ibin)
      flux_s(inh3_g,ibin) = flux_s(ihcl_g,ibin)

      ieqblm_ASTEM = mNO

      return
      end subroutine ASTEM_flux_dry_case4b




!-------------------------------------------------------------------
! Case 4ab

      subroutine ASTEM_flux_dry_case4ab(ibin)	! NH4NO3 + NH4Cl (solid)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      real(r8) :: a, b, c,   &
           flux_nh3_est, flux_nh3_max, ratio_flux
! function
      real(r8) :: quadratic

      call ASTEM_flux_dry_case4a(ibin)
      call ASTEM_flux_dry_case4b(ibin)


! estimate nh3 flux and adjust hno3 and/or hcl if necessary

      flux_nh3_est = flux_s(ihno3_g,ibin)+flux_s(ihcl_g,ibin)
      flux_nh3_max = kg(inh3_g,ibin)*gas(inh3_g)


      if(flux_nh3_est .le. flux_nh3_max)then

        flux_s(inh3_g,ibin) = flux_nh3_est			! all ok - no adjustments needed
        sfc_a(inh3_g)       = gas(inh3_g) -  			   &  ! recompute sfc_a(ihno3_g)
                              flux_s(inh3_g,ibin)/kg(inh3_g,ibin)
        phi_volatile_s(inh3_g,ibin) = max(abs(phi_nh4no3_s),   &
                                          abs(phi_nh4cl_s))

      else			! reduce hno3 and hcl flux_ses as necessary so that nh3 flux_s = flux_s_nh3_max

        ratio_flux          = flux_nh3_max/flux_nh3_est
        flux_s(inh3_g,ibin) = flux_nh3_max
        flux_s(ihno3_g,ibin)= flux_s(ihno3_g,ibin)*ratio_flux
        flux_s(ihcl_g,ibin) = flux_s(ihcl_g,ibin) *ratio_flux

        sfc_a(inh3_g) = 0.0
        sfc_a(ihno3_g)= gas(ihno3_g) -  			   &  ! recompute sfc_a(ihno3_g)
                        flux_s(ihno3_g,ibin)/kg(ihno3_g,ibin)
        sfc_a(ihcl_g) = gas(ihcl_g)  -  			   &  ! recompute sfc_a(ihcl_g)
                        flux_s(ihcl_g,ibin)/kg(ihcl_g,ibin)

        df_gas_s(inh3_g,ibin) =gas(inh3_g) -sfc_a(inh3_g)
        df_gas_s(ihno3_g,ibin)=gas(ihno3_g)-sfc_a(ihno3_g)
        df_gas_s(ihcl_g,ibin) =gas(ihcl_g) -sfc_a(ihcl_g)

        phi_volatile_s(inh3_g,ibin) = max(abs(phi_nh4no3_s),   &
                                          abs(phi_nh4cl_s))


      endif

      ieqblm_ASTEM = mNO

      return
      end subroutine ASTEM_flux_dry_case4ab








!=======================================================================
!
! MIXED-PHASE PARTICLES
!
!***********************************************************************
! part of ASTEM: computes gas-aerosol fluxes over mixed-phase aerosols
!
! author: Rahul A. Zaveri
! update: apr 2006
!-----------------------------------------------------------------------

      subroutine ASTEM_flux_mix(ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      integer iv, iadjust, iadjust_intermed, js
      real(r8) :: XT, g_nh3_hno3, g_nh3_hcl,   &
           a_nh4_no3, a_nh4_cl, a_no3, a_cl,   &
           prod_nh4no3, prod_nh4cl
      real(r8) :: volatile_cl, sum_dum


      call ions_to_electrolytes(jliquid,ibin,XT)  	! for water content calculation
      call compute_activities(ibin)

      if(water_a(ibin) .eq. 0.0)then
	write(6,*)'Water is zero in liquid phase'
	write(6,*)'Stopping in ASTEM_flux_wet'
        stop
      endif


!! EFFI calculate percent composition
      sum_dum = 0.0
      do js = 1, nsalt
        sum_dum = sum_dum + electrolyte(js,jsolid,ibin)
      enddo

      if(sum_dum .eq. 0.)sum_dum = 1.0

      epercent(jcaco3,jsolid,ibin) = 100.*electrolyte(jcaco3,jsolid,ibin)/sum_dum
!! EFFI


! reset
	Keq_nh4no3_0 = Keq_sg(1)	! 6/25/2008
	Keq_nh4cl_0  = Keq_sg(2)	! 6/25/2008

!-----------------------------------------------------------------
! MIXED CASE 1:  caco3 > 0 absorb all acids (and indirectly degas co2)

      if(epercent(jcaco3,jsolid,ibin) .gt. 0.0)then
        jphase(ibin) = jliquid
        call ASTEM_flux_wet_case1(ibin)
        return
      endif

!-----------------------------------------------------------------
! MIXED CASE 2: Sulfate-Rich Domain

      if(XT.lt.1.9999 .and. XT.ge.0.)then	! excess sulfate (acidic)
        jphase(ibin) = jliquid
	call ASTEM_flux_wet_case2(ibin)
        return
      endif

!-------------------------------------------------------------------
! MIXED CASE 3: hno3 and hcl exchange may happen here and nh4cl may form/evaporate

      volatile_cl  = electrolyte(jnacl,jsolid,ibin) +   &
                     electrolyte(jcacl2,jsolid,ibin)


      if(volatile_cl .gt. 0.0 .and. gas(ihno3_g).gt. 0.0 )then

        call ASTEM_flux_dry_case3a(ibin)

        prod_nh4cl = max( (gas(inh3_g)*gas(ihcl_g)-Keq_sg(2)), 0.0d0) +   &
                     electrolyte(jnh4cl, jsolid,ibin)

        if(prod_nh4cl .gt. 0.0)then
          call ASTEM_flux_dry_case3b(ibin)
        endif

        jphase(ibin) = jsolid

        return
      endif

!-------------------------------------------------------------------
! MIXED CASE 4: nh4no3 or nh4cl or both may be active

      if( electrolyte(jnh4no3,jsolid,ibin).gt.0. .and.   &
          electrolyte(jnh4cl,jsolid,ibin) .gt.0. )then
        jphase(ibin) = jsolid
        call ASTEM_flux_dry_case4(ibin)

        if(sfc_a(ihno3_g).gt.0.0 .and. ma(ja_no3,ibin).gt.0.0)then
          mc(jc_h,ibin) = Keq_gl(3)*sfc_a(ihno3_g)/   &
          (kel(ihno3_g,ibin)*gam(jhno3,ibin)**2 * ma(ja_no3,ibin))
        elseif(sfc_a(ihcl_g).gt.0.0 .and. ma(ja_cl,ibin).gt.0.0)then
          mc(jc_h,ibin) = Keq_gl(4)*sfc_a(ihcl_g)/   &
          (kel(ihcl_g,ibin)*gam(jhcl,ibin)**2 * ma(ja_cl,ibin))
        else
          mc(jc_h,ibin) = sqrt(Keq_ll(3))
        endif

        return

      elseif( electrolyte(jnh4no3,jsolid,ibin).gt.0. )then
! do small adjustments for nh4cl aq
        g_nh3_hcl= gas(inh3_g)*gas(ihcl_g)
        a_nh4_cl = aer(inh4_a,jliquid,ibin)*aer(icl_a,jliquid,ibin)

        iadjust = mNO		! initialize
        if(g_nh3_hcl .gt. 0.0 .and. a_nh4_cl .eq. 0.0)then
          call absorb_tiny_nh4cl(ibin)
          iadjust = mYES
        elseif(g_nh3_hcl .eq. 0.0 .and. a_nh4_cl .gt. 0.0)then
          call degas_tiny_nh4cl(ibin)
          iadjust = mYES
        endif

        if(iadjust .eq. mYES)then
          call ions_to_electrolytes(jliquid,ibin,XT)  	! update after adjustments
          call compute_activities(ibin)			! update after adjustments
        endif

        call ASTEM_flux_mix_case4a(ibin)	! nh4no3 solid + nh4cl aq
        jphase(ibin) = jtotal
        return

      elseif( electrolyte(jnh4cl,jsolid,ibin).gt.0.)then
! do small adjustments for nh4no3 aq
        g_nh3_hno3= gas(inh3_g)*gas(ihno3_g)
        a_nh4_no3 = aer(inh4_a,jliquid,ibin)*aer(ino3_a,jliquid,ibin)

        iadjust = mNO		! initialize
        if(g_nh3_hno3 .gt. 0.0 .and. a_nh4_no3 .eq. 0.0)then
          call absorb_tiny_nh4no3(ibin)
          iadjust = mYES
        elseif(g_nh3_hno3 .eq. 0.0 .and. a_nh4_no3 .gt. 0.0)then
          call degas_tiny_nh4no3(ibin)
          iadjust = mYES
        endif

        if(iadjust .eq. mYES)then
          call ions_to_electrolytes(jliquid,ibin,XT)  	! update after adjustments
          call compute_activities(ibin)			! update after adjustments
        endif

        kelvin_nh4no3 = kel(inh3_g,ibin)*kel(ihno3_g,ibin)
        Keq_nh4no3 = kelvin_nh4no3*activity(jnh4no3,ibin)*Kp_nh4no3	! = [NH3]s * [HNO3]s

        call ASTEM_flux_mix_case4b(ibin)	! nh4cl solid + nh4no3 aq
        jphase(ibin) = jtotal
        return
      endif


!-------------------------------------------------------------------

      if( (gas(inh3_g)+aer(inh4_a,jliquid,ibin)) .lt. 1.e-25)goto 10  ! no ammonia in the system

!-------------------------------------------------------------------
! MIXED CASE 5: liquid nh4no3 and/or nh4cl maybe active
! do some small adjustments (if needed) before deciding case 3

      iadjust = mNO		! default
      iadjust_intermed = mNO	! default

! nh4no3
      g_nh3_hno3 = gas(inh3_g)*gas(ihno3_g)
      a_nh4_no3  = aer(inh4_a,jliquid,ibin)*aer(ino3_a,jliquid,ibin)

      if(g_nh3_hno3 .gt. 0. .and. a_nh4_no3 .eq. 0.)then
        call absorb_tiny_nh4no3(ibin)
        iadjust = mYES
        iadjust_intermed = mYES
      endif

      if(iadjust_intermed .eq. mYES)then
        call ions_to_electrolytes(jliquid,ibin,XT)  	! update after adjustments
        iadjust_intermed = mNO	! reset
      endif

! nh4cl
      g_nh3_hcl = gas(inh3_g)*gas(ihcl_g)
      a_nh4_cl  = aer(inh4_a,jliquid,ibin)*aer(icl_a,jliquid,ibin)

      if(g_nh3_hcl .gt. 0. .and. a_nh4_cl .eq. 0.)then
        call absorb_tiny_nh4cl(ibin)
        iadjust = mYES
        iadjust_intermed = mYES
      endif

      if(iadjust_intermed .eq. mYES)then
        call ions_to_electrolytes(jliquid,ibin,XT)  	! update after adjustments
      endif

      if(iadjust .eq. mYES)then
        call compute_activities(ibin)			! update after adjustments
      endif


! all adjustments done...

!--------
      kelvin_nh4no3 = kel(inh3_g,ibin)*kel(ihno3_g,ibin)
      Keq_nh4no3 = kelvin_nh4no3*activity(jnh4no3,ibin)*Kp_nh4no3	! = [NH3]s * [HNO3]s

      kelvin_nh4cl = kel(inh3_g,ibin)*kel(ihcl_g,ibin)
      Keq_nh4cl = kelvin_nh4cl*activity(jnh4cl,ibin)*Kp_nh4cl	! = [NH3]s * [HCl]s

      call ASTEM_flux_wet_case3(ibin)
      jphase(ibin) = jliquid

      return


!-------------------------------------------------------------------
! MIXED CASE 6: ammonia = 0. liquid hno3 and hcl exchange may happen here
! do small adjustments (if needed) before deciding case 4

10    iadjust = mNO		! default
      iadjust_intermed = mNO	! default

! hno3
      if(gas(ihno3_g).gt.0. .and. aer(ino3_a,jliquid,ibin).eq.0. .and.   &
         aer(icl_a,jliquid,ibin) .gt. 0.0)then
        call absorb_tiny_hno3(ibin)	! and degas tiny hcl
        iadjust = mYES
        iadjust_intermed = mYES
      endif

      if(iadjust_intermed .eq. mYES)then
        call ions_to_electrolytes(jliquid,ibin,XT)  	! update after adjustments
        iadjust_intermed = mNO	! reset
      endif

! hcl
      if(gas(ihcl_g).gt.0. .and. aer(icl_a,jliquid,ibin) .eq. 0. .and.   &
         aer(ino3_a,jliquid,ibin) .gt. 0.0)then
        call absorb_tiny_hcl(ibin)			! and degas tiny hno3
        iadjust = mYES
        iadjust_intermed = mYES
      endif

      if(iadjust_intermed .eq. mYES)then
        call ions_to_electrolytes(jliquid,ibin,XT)  	! update after adjustments
      endif

      if(iadjust .eq. mYES)then
        call compute_activities(ibin)			! update after adjustments
      endif

! all adjustments done...

      call ASTEM_flux_wet_case4(ibin)
      jphase(ibin) = jliquid



      return
      end subroutine ASTEM_flux_mix

!----------------------------------------------------------------------








!------------------------------------------------------------------
! Mix Case 4a: NH4NO3 solid maybe active. NH4Cl aq maybe active

      subroutine ASTEM_flux_mix_case4a(ibin)	! TOUCH
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      integer iactive_nh4no3, iactive_nh4cl, js
      real(r8) :: sum_dum


! set default values for flags
      iactive_nh4no3 = mYES
      iactive_nh4cl  = mYES


!! EFFI calculate percent composition
      sum_dum = 0.0
      do js = 1, nelectrolyte
        sum_dum = sum_dum + electrolyte(js,jsolid,ibin)
      enddo

      if(sum_dum .eq. 0.)sum_dum = 1.0

      epercent(jnh4no3,jsolid,ibin) = 100.*electrolyte(jnh4no3,jsolid,ibin)/sum_dum
!! EFFI



! nh4no3 (solid)
      phi_nh4no3_s = (gas(inh3_g)*gas(ihno3_g) - Keq_sg(1))/   &
                     max(gas(inh3_g)*gas(ihno3_g),Keq_sg(1))

! nh4cl (liquid)
      kelvin_nh4cl = kel(inh3_g,ibin)*kel(ihcl_g,ibin)
      Keq_nh4cl = kelvin_nh4cl*activity(jnh4cl,ibin)*Kp_nh4cl	! = [NH3]s * [HCl]s


!-------------------
! now determine if nh4no3 and/or nh4cl are active or significant
! nh4no3 solid
      if( abs(phi_nh4no3_s) .le. rtol_eqb_ASTEM )then
        iactive_nh4no3 = mNO
      elseif(gas(inh3_g)*gas(ihno3_g) .lt. Keq_sg(1) .and.   &
             epercent(jnh4no3,jsolid,ibin) .le. ptol_mol_ASTEM)then
        iactive_nh4no3 = mNO
        if(epercent(jnh4no3,jsolid,ibin) .gt. 0.0)then
          call degas_solid_nh4no3(ibin)
        endif
      endif

! nh4cl aq
      if( gas(inh3_g)*gas(ihcl_g).eq.0. .or. Keq_nh4cl.eq.0. )then
        iactive_nh4cl = mNO
      endif


!---------------------------------
      if(iactive_nh4no3 .eq. mYES)then

        jphase(ibin) = jsolid
        call ASTEM_flux_dry_case4a(ibin)	! NH4NO3 (solid)

        if(sfc_a(ihno3_g).gt.0.0 .and. ma(ja_no3,ibin).gt.0.0)then
          mc(jc_h,ibin) = Keq_gl(3)*sfc_a(ihno3_g)/   &
          (kel(ihno3_g,ibin)*gam(jhno3,ibin)**2 * ma(ja_no3,ibin))
        elseif(sfc_a(ihcl_g).gt.0.0 .and. ma(ja_cl,ibin).gt.0.0)then
          mc(jc_h,ibin) = Keq_gl(4)*sfc_a(ihcl_g)/   &
          (kel(ihcl_g,ibin)*gam(jhcl,ibin)**2 * ma(ja_cl,ibin))
        else
          mc(jc_h,ibin) = sqrt(Keq_ll(3))
        endif

      endif


      if(iactive_nh4cl .eq. mYES)then

        jphase(ibin) = jliquid
        call ASTEM_flux_wet_case3b(ibin)	! NH4Cl (liquid)

        if(sfc_a(ihcl_g).gt.0.0 .and. ma(ja_cl,ibin).gt.0.0)then
          mc(jc_h,ibin) = Keq_gl(4)*sfc_a(ihcl_g)/   &
          (kel(ihcl_g,ibin)*gam(jhcl,ibin)**2 * ma(ja_cl,ibin))
        else
          mc(jc_h,ibin) = sqrt(Keq_ll(3))
        endif

      endif


      if(iactive_nh4cl .eq. mYES .and. iactive_nh4no3 .eq. mYES)then
        jphase(ibin) = jtotal
      endif



      return
      end subroutine ASTEM_flux_mix_case4a








!------------------------------------------------------------------
! Mix Case 4b: NH4Cl solid maybe active. NH4NO3 aq may or maybe active

      subroutine ASTEM_flux_mix_case4b(ibin)	! TOUCH
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      integer iactive_nh4no3, iactive_nh4cl, js
      real(r8) :: sum_dum


! set default values for flags
      iactive_nh4cl  = mYES
      iactive_nh4no3 = mYES


!! EFFI calculate percent composition
      sum_dum = 0.0
      do js = 1, nsalt
        sum_dum = sum_dum + electrolyte(js,jsolid,ibin)
      enddo

      if(sum_dum .eq. 0.)sum_dum = 1.0

      epercent(jnh4cl,jsolid,ibin) = 100.*electrolyte(jnh4cl,jsolid,ibin)/sum_dum
!! EFFI


! nh4cl (solid)
      phi_nh4cl_s  = (gas(inh3_g)*gas(ihcl_g) - Keq_sg(2))/   &
                     max(gas(inh3_g)*gas(ihcl_g),Keq_sg(2))

! nh4no3 (liquid)
      kelvin_nh4no3 = kel(inh3_g,ibin)*kel(ihno3_g,ibin)
      Keq_nh4no3 = kelvin_nh4no3*activity(jnh4no3,ibin)*Kp_nh4no3	! = [NH3]s * [HNO3]s


!-------------------
! now determine if nh4no3 and/or nh4cl are active or significant
! nh4cl (solid)
      if( abs(phi_nh4cl_s) .le. rtol_eqb_ASTEM )then
        iactive_nh4cl = mNO
      elseif(gas(inh3_g)*gas(ihcl_g) .lt. Keq_sg(2) .and.   &
             epercent(jnh4cl,jsolid,ibin) .le. ptol_mol_ASTEM)then
        iactive_nh4cl = mNO
        if(epercent(jnh4cl,jsolid,ibin) .gt. 0.0)then
          call degas_solid_nh4cl(ibin)
        endif
      endif

! nh4no3 (liquid)
      if( gas(inh3_g)*gas(ihno3_g).eq.0. .or. Keq_nh4no3.eq.0. )then
        iactive_nh4no3 = mNO
      endif


!---------------------------------
      if(iactive_nh4cl .eq. mYES)then

        jphase(ibin) = jsolid
        call ASTEM_flux_dry_case4b(ibin)	! NH4Cl (solid)

        if(sfc_a(ihcl_g).gt.0.0 .and. ma(ja_cl,ibin).gt.0.0)then
          mc(jc_h,ibin) = Keq_gl(4)*sfc_a(ihcl_g)/   &
          (kel(ihcl_g,ibin)*gam(jhcl,ibin)**2 * ma(ja_cl,ibin))
        elseif(sfc_a(ihno3_g).gt.0.0 .and. ma(ja_no3,ibin).gt.0.0)then
          mc(jc_h,ibin) = Keq_gl(3)*sfc_a(ihno3_g)/   &
          (kel(ihno3_g,ibin)*gam(jhno3,ibin)**2 * ma(ja_no3,ibin))
        else
          mc(jc_h,ibin) = sqrt(Keq_ll(3))
        endif

      endif


      if(iactive_nh4no3 .eq. mYES)then

        jphase(ibin) = jliquid
        call ASTEM_flux_wet_case3a(ibin)	! NH4NO3 (liquid)

        if(sfc_a(ihno3_g).gt.0.0 .and. ma(ja_no3,ibin).gt.0.0)then
          mc(jc_h,ibin) = Keq_gl(3)*sfc_a(ihno3_g)/   &
          (kel(ihno3_g,ibin)*gam(jhno3,ibin)**2 * ma(ja_no3,ibin))
        else
          mc(jc_h,ibin) = sqrt(Keq_ll(3))
        endif

      endif


      if(iactive_nh4cl .eq. mYES .and. iactive_nh4no3 .eq. mYES)then
        jphase(ibin) = jtotal
      endif



      return
      end subroutine ASTEM_flux_mix_case4b











!***********************************************************************
! part of ASTEM: condenses h2so4, msa, and nh3 analytically over dtchem [s]
!
! author: Rahul A. Zaveri
! update: jan 2007
!-----------------------------------------------------------------------

      subroutine ASTEM_non_volatiles(dtchem) ! TOUCH
      use module_data_mosaic_aero

      implicit none

! subr arguments
      real(r8) :: dtchem
! local variables
      integer ibin, iupdate_phase_state
      real(r8) :: decay_h2so4, decay_msa,   &
           delta_h2so4, delta_tmsa, delta_nh3, delta_hno3, delta_hcl,   &
           delta_so4(nbin_a), delta_msa(nbin_a),   &
           delta_nh4(nbin_a)
      real(r8) :: XT




      sumkg_h2so4 = 0.0
      sumkg_msa   = 0.0
      sumkg_nh3   = 0.0
      sumkg_hno3  = 0.0
      sumkg_hcl   = 0.0
      do ibin = 1, nbin_a
        sumkg_h2so4 = sumkg_h2so4 + kg(ih2so4_g,ibin)
        sumkg_msa   = sumkg_msa   + kg(imsa_g,ibin)
        sumkg_nh3   = sumkg_nh3   + kg(inh3_g,ibin)
        sumkg_hno3  = sumkg_hno3  + kg(ihno3_g,ibin)
        sumkg_hcl   = sumkg_hcl   + kg(ihcl_g,ibin)
      enddo



!--------------------------------------
! H2SO4
      if(gas(ih2so4_g) .gt. 1.e-14)then

! integrate h2so4 condensation analytically
        decay_h2so4   = exp(-sumkg_h2so4*dtchem)
        delta_h2so4   = gas(ih2so4_g)*(1.0 - decay_h2so4)
        gas(ih2so4_g) = gas(ih2so4_g)*decay_h2so4


! now distribute delta_h2so4 to each bin and conform the particle (may degas by massbal)
        do ibin = 1, nbin_a
          if(jaerosolstate(ibin) .ne. no_aerosol)then
            delta_so4(ibin) = delta_h2so4*kg(ih2so4_g,ibin)/sumkg_h2so4
            aer(iso4_a,jtotal,ibin) = aer(iso4_a,jtotal,ibin) +   &
                                      delta_so4(ibin)
          endif
        enddo

      else

        delta_h2so4 = 0.0
        do ibin = 1, nbin_a
            delta_so4(ibin) = 0.0
        enddo

      endif
! h2so4 condensation is now complete
!--------------------------------------



! MSA
      if(gas(imsa_g) .gt. 1.e-14)then

! integrate msa condensation analytically
        decay_msa   = exp(-sumkg_msa*dtchem)
        delta_tmsa  = gas(imsa_g)*(1.0 - decay_msa)
        gas(imsa_g) = gas(imsa_g)*decay_msa

! now distribute delta_msa to each bin and conform the particle (may degas by massbal)
        do ibin = 1, nbin_a
          if(jaerosolstate(ibin) .ne. no_aerosol)then
            delta_msa(ibin) = delta_tmsa*kg(imsa_g,ibin)/sumkg_msa
            aer(imsa_a,jtotal,ibin) = aer(imsa_a,jtotal,ibin) +   &
                                      delta_msa(ibin)
          endif
        enddo

      else

        delta_tmsa = 0.0
        do ibin = 1, nbin_a
            delta_msa(ibin) = 0.0
        enddo

      endif
! msa condensation is now complete
!-------------------------------------



! compute max allowable nh3, hno3, and hcl condensation
      delta_nh3 = gas(inh3_g) *(1.0 - exp(-sumkg_nh3*dtchem))
      delta_hno3= gas(ihno3_g)*(1.0 - exp(-sumkg_hno3*dtchem))
      delta_hcl = gas(ihcl_g) *(1.0 - exp(-sumkg_hcl*dtchem))

! compute max possible nh4 condensation for each bin
      do ibin = 1, nbin_a
        if(jaerosolstate(ibin) .ne. no_aerosol)then
          delta_nh3_max(ibin) = delta_nh3*kg(inh3_g,ibin)/sumkg_nh3
          delta_hno3_max(ibin)= delta_hno3*kg(ihno3_g,ibin)/sumkg_hno3
          delta_hcl_max(ibin) = delta_hcl*kg(ihcl_g,ibin)/sumkg_hcl
        endif
      enddo


      if(delta_h2so4 .eq. 0.0 .and. delta_tmsa .eq. 0.0)then
        iupdate_phase_state = mNO
        goto 100
      endif


! now condense appropriate amounts of nh3 to each bin  EFFI
      do ibin = 1, nbin_a

        if(electrolyte(jnacl,jtotal,ibin)  .eq. 0.0 .and.   &
           electrolyte(jcacl2,jtotal,ibin) .eq. 0.0 .and.   &
           electrolyte(jnano3,jtotal,ibin) .eq. 0.0 .and.   &
           electrolyte(jcano3,jtotal,ibin) .eq. 0.0 .and.   &
           electrolyte(jcaco3,jtotal,ibin) .eq. 0.0 .and.   &
           jaerosolstate(ibin) .ne. no_aerosol)then

          delta_nh4(ibin) = min( (2.*delta_so4(ibin)+delta_msa(ibin)),   &
                                 delta_nh3_max(ibin) )

          aer(inh4_a,jtotal,ibin) = aer(inh4_a,jtotal,ibin) +	   &  ! update aer-phase
                                    delta_nh4(ibin)

          gas(inh3_g) = gas(inh3_g) - delta_nh4(ibin)		! update gas-phase

        else

          delta_nh4(ibin)     = 0.0

        endif

      enddo

      iupdate_phase_state = mYES


! recompute phase equilibrium
100   if(iupdate_phase_state .eq. mYES)then
        do ibin = 1, nbin_a
          if(jaerosolstate(ibin) .ne. no_aerosol)then
            call conform_electrolytes(jtotal,ibin,XT)
            call aerosol_phase_state(ibin)
          endif
        enddo
      endif

      return
      end subroutine ASTEM_non_volatiles








!=================================================================
! SOA module




!***********************************************************************
! part of ASTEM: condenses secondary organic species over TSI time interval
! mechanism adapted from SORGAM
!
! author: Rahul A. Zaveri
! update: apr 2005
!-----------------------------------------------------------------------
      subroutine ASTEM_secondary_organics(dtchem)
      use module_data_mosaic_main
      use module_data_mosaic_aero

      implicit none

! subr arguments
      real(r8) :: dtchem
! local variables
      integer ibin, iv, jp, &
              ieqblm, nsteps_max
      integer iforce_eqblm
      parameter(nsteps_max = 400)
      real(r8) :: dtmax, t_new, t_old, t_out
      real(r8) :: sum1, sum2, sumkg, gas_save, decay_iv, delta_iv
      real(r8) :: numr, denr, Keq_dimer1, Keq_dimer2, Vp, m_molm3, d_molm3, tm_molm3, sum_oa_mass
      real(r8) :: kf, kd1, kd2, m_0, d_0, m_f, d_f
      real(r8) :: SOA_tot, NVO_tot, ratio_SOA_to_NVO
      real(r8) :: alpha


! initialize time
      t_old = 0.0
      t_out = dtchem
      isteps_SOA = 0

      dtmax = dtchem/6.0

! overall integration loop begins over dtchem seconds
10    isteps_SOA = isteps_SOA + 1

      
      do iv = icn3_g, ngas_volatile
        total_species(iv) = gas(iv)
        do ibin = 1, nbin_a
          if(t_since_start .eq. 0.0)then
            aer0(iv,ibin) = aer(iv,jtotal,ibin)
          endif

          if (jaerosolstate(ibin) .eq. no_aerosol) cycle
          total_species(iv) = total_species(iv) + aer(iv,jtotal,ibin)
        enddo
      enddo


! update time_since_start
	t_since_start = t_since_start + dtmax

! compute new fluxes
      ieqblm_soa = mYES			! reset to default


      do 501 ibin = 1, nbin_a
        if (jaerosolstate(ibin) .eq. no_aerosol) goto 501

        call ASTEM_flux_soa(ibin)

501   continue


      if(ieqblm_soa .eq. mYES)goto 30 ! all bins have reached equilibrium

!-----------------------


! calculate maximum possible internal time-step
!11    call ASTEM_dtmax_soa(dtchem, dtmax)
      t_new = t_old + dtmax	! update time
      if(t_new .gt. t_out)then	! check if the new time step is too large
        dtmax = t_out - t_old
        t_new = t_out*1.01
      endif


!	write(6,*)'dtmax = ', dtmax

!------------------------------------------
! do internal time-step (dtmax) integration

      jp = jtotal
      do 20 iv = icn3_g, ngas_volatile
      ! right now ngas_volatile is the number for the highest volatility bin.
        sum1 = 0.0
        sum2 = 0.0
        do 21 ibin = 1, nbin_a
          if(jaerosolstate(ibin) .eq. no_aerosol)goto 21

        ! wkc note: kc_firstorder is currently 0.0 (20190513)
	    numr=aer(iv,jp,ibin)
	    denr = 1. + dtmax*(integrate(iv,jp,ibin)*kg(iv,ibin)*Heff(iv,ibin)/QQ(iv,ibin) + kc_firstorder(iv))
	    sum1 = sum1 + numr/denr

	    ! wkc note: kc_firstorder is currently 0.0 (20190513)
	    numr = kg(iv,ibin)*integrate(iv,jp,ibin)
	    denr = 1. + dtmax*(integrate(iv,jp,ibin)*kg(iv,ibin)*Heff(iv,ibin)/QQ(iv,ibin) + kc_firstorder(iv))
	    sum2 = sum2 + numr/denr


21      continue ! ibin

! first update gas concentration
! wkc note: kc_firstorder is currently 0.0 (20190513)
        gas(iv) = (total_species(iv) - (1.0 + dtmax*kc_firstorder(iv))*sum1)/   &
                       (1.0 + (1.0 + dtmax*kc_firstorder(iv))*dtmax*sum2)

! now update aer concentration in the jp phase
        do 22 ibin = 1, nbin_a
          if (jaerosolstate(ibin) .eq. no_aerosol) goto 22

          if(integrate(iv,jp,ibin) .eq. mYES)then
            aer(iv,jp,ibin) =   &
             (aer(iv,jp,ibin) + dtmax*kg(iv,ibin)*integrate(iv,jp,ibin)*gas(iv))/ &
                  (1. + dtmax*(integrate(iv,jp,ibin)*kg(iv,ibin)*Heff(iv,ibin)/QQ(iv,ibin) + kc_firstorder(iv)))
          endif

22      continue ! ibin

20    continue ! iv
!
!if (ihet_gas(iOH) .eq. 1)
!! edit wkc: Add in heterogeneous reactions
!    v_molar(iOH) = .08206 * te / P_atm ! molar volume in L/mol
!    Dg(iOH) = gas_diffusivity(T_K,P_atm,mw_gas(iOH),v_molar(iOH)) ! [cm2/s], following Rahul
!    vel_gas(iOH) = 1.455e4 * sqrt(te/mw_gas(iOH)) ! avg. molec speed [cm/s]
!    freepath(iOH) = 3 * Dg(iOH) / vel_gas(iOH) ! mean free path [cm]
!    alpha = 1 ! accommodation coefficient
!
!
!    do ibin = 1, nbin_a ! loop through size bins
!        Kn = 2 * freepath(iOH) / Dp(ibin) ! Knudsen number
!        beta = Kn * (1+Kn) / (Kn**2+Kn+0.283*Kn*alpha+0.75*alpha) ! Fuchs correction factor from kinetic regime
!    ! compute het chem rate constant with OH
!        sfc_area = (3.14159*Dp(ibin)**2) * nsize_aer(ibin) ! total surface area for each size bin, not sure about nsize_aer
!
!        aer_sum_size(ibin) = 0
!        ! sum of organics at this size bin
!        do iv = icn3_a, ngas_volatile ! loop through volatiliy bins
!            aer_sum_size(ibin) = aer_sum_size(ibin) + aer(iv,jp,ibin)
!        enddo
!
!        do iv = icn3_a, ngas_volatile ! loop through volatility bins
!            k_het = (vel_gas(iOH)/4) * sfc_area * uptake_gas(iOH) * beta * cnn(iOH) ! het rate constant, [molec/cc/s]
!            ! reaction_rate(iv, jp, ibin) = k_het * aer(iv, jp, ibin) / aer_sum_size(ibin) ! reaction rate for species iv, [molec/cc/s]
!
!            if (iv .eq. icn3_a) ! lowest volatility species
!            ! This is incompatible with multicomponent
!            ! There will be multiple lowest volatility species
!                aer_next_step(iv,jp,ibin) = &
!                    (aer(iv,jp,ibin) + (dtmax*k_het*aer(iv-1,jp,ibin)*(1-p_frag(iv-1)) / aer_sum_size(ibin)) / &
!                    (1 + dtmax*k_het/aer_sum_size(ibin))
!
!            else if (iv .eq. ngas_volatile) ! highest volatility species
!            ! This is incompatible with multicomponent.
!            ! There will be multiple highest volatility species.
!                aer_next_step(iv,jp,ibin) = &
!                    aer(iv,jp,ibin) / (1 + dtmax*k_het/aer_sum_size(ibin))
!
!            else
!                aer_next_step(iv,jp,ibin) = &
!                    (aer(iv,jp,ibin) + dtmax*k_het*aer(iv-1,jp,ibin)/aer_sum_size(ibin)*(1-p_frag(iv-1))) / &
!                    (1 + dtmax*k_het/aer_sum_size(ibin))
!            endif
!
!
!        enddo
!    enddo
!endif

    ! edit wkc
    !if(t_since_start .gt. 2.0) then
    !    print *, t_since_start
    !endif

!------------------------------------------
! sub-step integration done


! update jtotal
!      do iv = icn3_g, ngas_volatile
!        aer(iv,jtotal,ibin)=aer(iv,jsolid,ibin)+aer(iv,jliquid,ibin)
!      enddo


! update time
      t_old = t_new

      if(t_new .lt. 0.9999*t_out) goto 10
!================================================
! end of integration

30    continue

      return
      end subroutine ASTEM_secondary_organics











!***********************************************************************
! part of ASTEM: computes fluxes of soa species
!
! author: Rahul A. Zaveri
! update: mar 2013 - added code for reactive mass transfer to semi-solids
!
!-----------------------------------------------------------------------
      subroutine ASTEM_flux_soa(ibin)		! TOUCH
      use module_data_mosaic_main
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      integer iv, jp, in
      real(r8) :: dum, sum_dum, sum_soa, sum_soa_only, small_oc, gamma_org
	real(r8) :: ndum, sum_n, ydum, &
			q(ngas_volatile)
	real(r8) :: coth

      REAL(R8),PARAMETER :: NA = 6.022e23
      REAL(R8),ALLOCATABLE :: k_PARMOLE_SOM(:)
 
      small_oc  = 1.e-15		! ng/m^3


! default fluxes and other stuff
      do iv = icn3_g, ngas_volatile
        sfc_org_a(iv,ibin)      = gas(iv)
        df_gas_o(iv,ibin)       = 0.0
        flux_o(iv,ibin)         = 0.0
        phi_volatile_o(iv,ibin) = 0.0
      enddo


      jp = jtotal

! compute mole fractions of soa species
      sum_soa = 0.0
      do iv = icn3_g, ngas_volatile
        sum_soa = sum_soa + aer(iv,jp,ibin)
      enddo
      sum_soa = sum_soa + aer(ioc_a,jp,ibin)/200.0		! 200 is assumed MW of primary OC

!!!!! BY CHARLES HE
      ALLOCATE(k_PARMOLE_SOM(nBINS))
      k_PARMOLE_SOM(:) = 0.0
      
      k_PARMOLE_SOM = SUM(jk_PARMOLE_SOM,1)
      
      sum_soa = sum_soa + k_PARMOLE_SOM(iBIN)/(1e-9*NA*1e-6)

      
      sum_soa_only = sum_soa

!      sum_soa = sum_soa + water_a(ibin)*1.e12/18.0		! add water moles



! check threshold concentration for SOA formation in the absence of primary OC
      if(aer(ioc_a,jp,ibin) .eq. 0.0)then
        sum_dum = 0.0
        do iv = icn3_g, ngas_volatile
          sum_dum = sum_dum + (gas(iv)+aer(iv,jp,ibin))/sat_soa(iv)
        enddo

        if(sum_dum .le. 1.0)then	! transfer all aer to gas and quit
          do iv = icn3_g, ngas_volatile
            gas(iv)         = gas(iv) + aer(iv,jp,ibin)
            aer(iv,jp,ibin) = 0.0
            integrate(iv,jp,ibin) = 0.0
          enddo
          return
        endif

        sum_soa = max(sum_soa, 1.d-10)

      endif




! compute Heff
      do iv = icn3_g, ngas_volatile

	  gamma_org = 1.0 ! activity coefficient
          Heff(iv,ibin) = kel(iv,ibin)*gamma_org*sat_soa(iv)/sum_soa

! check if reached equilibrium
          sfc_org_a(iv,ibin) = aer(iv,jp,ibin)*Heff(iv,ibin)/QQ(iv,ibin) ! nmol/m^3
          df_gas_o(iv,ibin) = gas(iv) - sfc_org_a(iv,ibin)

          dum = max(sfc_org_a(iv,ibin),gas(iv))
          if(dum .gt. 0.0)then
            phi_volatile_o(iv,ibin) = df_gas_o(iv,ibin)/dum
          else
            phi_volatile_o(iv,ibin) = 0.0
          endif

!        if(abs(phi_volatile_o(iv,ibin)) .le. rtol_eqb_ASTEM)then
!          integrate(iv,jp,ibin) = 0.0
!        else
          integrate(iv,jp,ibin) = 1.0
          ieqblm_soa = mNO
!        endif

      enddo


      return
      end subroutine ASTEM_flux_soa




      function coth(xx_arg)
      use module_data_mosaic_aero
	implicit none
	real(r8) :: xx_arg, yy, coth

	yy = xx_arg
      if(xx_arg .gt. 30)yy=30.

	if(yy .gt. 0.6)then
        coth = (exp(yy) + exp(-yy))/(exp(yy) - exp(-yy))
	else
	  coth = (1.0 + 0.5*yy*yy + (yy**4)/24. + (yy**6)/720.)/ &
               (yy + (yy**3)/6.0 + (yy**5)/120. + (yy**7)/5040.)
	endif

      return
      end function coth





!***********************************************************************
! part of ASTEM: computes fluxes of soa species
!
! author: Rahul A. Zaveri
! update: apr 2005
!-----------------------------------------------------------------------
      subroutine ASTEM_dtmax_soa(dtchem, dtmax)		! TOUCH
      use module_data_mosaic_aero

      implicit none

! subr arguments
      real(r8) :: dtchem, dtmax
! local variables
      integer ibin, iv, jp
      real(r8) :: h_gas, h_gas_i(ngas_volatile), h_sub_max,   &
           sum_kg_phi


      h_sub_max = dtchem/6.0	! sec

      jp = jtotal

! GAS-SIDE
! calculate h_gas_i and h_gas

      h_gas = 2.e16

      do 6 iv = icn3_g, ngas_volatile

        h_gas_i(iv) = 1.e16
        sum_kg_phi = 0.0

        do ibin = 1, nbin_a
          if(integrate(iv,jtotal,ibin) .eq. mYES)then
          sum_kg_phi = sum_kg_phi +   &
                       abs(phi_volatile_o(iv,ibin))*kg(iv,ibin)
          endif
        enddo

        if(sum_kg_phi .gt. 0.0)then
          h_gas_i(iv) = alpha_astem/sum_kg_phi
          h_gas       = min(h_gas, h_gas_i(iv))
        endif

6     continue


      dtmax = min(h_gas, h_sub_max)


!      if(dtmax .le. 1.0e-10)then
!        write(6,*)' SOA dtmax = ', dtmax, h_gas_i(iapi1_g)
!      endif


      return
      end subroutine ASTEM_dtmax_soa





