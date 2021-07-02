!
! double precision lsodes solver
!
      subroutine MOSAIC_LSODE(dtchem)
      use module_data_mosaic_main
      use module_data_mosaic_gas
      use module_data_mosaic_aero

      implicit none


! subr arguments
      real(r8) :: dtchem

      write(*,'(//a//)') &
      '*** error - mosaic_lsode has been deactivated ***'
      stop

      return
      end subroutine MOSAIC_LSODE











      subroutine ode_aer(ntot,tt,stot,sdot)
!
      use module_data_mosaic_main
      use module_data_mosaic_gas
      use module_data_mosaic_aero

      implicit none

      integer ntot
      real(r8) :: tt, stot(ntot_max), sdot(ntot_max)

      integer iaer, ibin, iv, istot


! compute fluxes

      call make_gas_from_stot(stot)	! distribute stot into gas and aer

! compute aerosol phase state
      do ibin = 1, nbin_a
        if(jaerosolstate(ibin) .ne. no_aerosol)then

          if(jaerosolstate(ibin).eq.all_solid)then
            call make_aer_from_stot(stot,ibin)
            call aerosol_phase_state(ibin) ! loads aer(jsolid) & aer(jliquid)
          else
           call make_aer_from_stot(stot,ibin)
          endif

        endif
      enddo


! calculate fluxes
      do 500 ibin = 1, nbin_a

        if(jaerosolstate(ibin) .eq. all_solid)then
          jphase(ibin) = jsolid
          call LSODE_flux_dry(ibin)
        elseif(jaerosolstate(ibin) .eq. all_liquid)then
          jphase(ibin) = jliquid
          call LSODE_flux_wet(ibin)
        elseif(jaerosolstate(ibin) .eq. mixed)then
          jphase(ibin) = jliquid
          call LSODE_flux_wet(ibin)
        endif

500   continue




! load sdot array
! gas
      do 600 iv = 1, ngas_ioa
        sdot(iv) = 0.0
        do ibin = 1, nbin_a
          flux(iv,ibin) = flux_s(iv,ibin)+flux_l(iv,ibin)
          sdot(iv) = sdot(iv) - flux(iv,ibin)
        enddo
600   continue

! aerosol
      do 700 ibin = 1, nbin_a

        do iaer = 1, ngas_ioa
          istot = ngas_ioa*ibin + iaer
          sdot(istot) = flux(iaer,ibin)
        enddo

700   continue

      return
      end subroutine ode_aer









      subroutine make_aer_from_stot(stot,ibin)
      use module_data_mosaic_main
      use module_data_mosaic_gas
      use module_data_mosaic_aero

      implicit none

      real(r8) :: stot(ntot_max)
      integer jp, iaer, ibin, istot

      if(jaerosolstate(ibin).eq.all_solid)then
        jp = jtotal
      elseif(jaerosolstate(ibin).eq.all_liquid)then
        jp = jliquid
      elseif(jaerosolstate(ibin).eq.mixed)then
        jp = jliquid
      else
        jp = jtotal
      endif

      do iaer = 1, ngas_ioa
        istot= ngas_ioa*ibin + iaer
        aer(iaer,jp,ibin) = stot(istot)	! nmol/m^3
      enddo

      return
      end subroutine make_aer_from_stot







      subroutine make_gas_from_stot(stot)
      use module_data_mosaic_main
      use module_data_mosaic_gas
      use module_data_mosaic_aero

      implicit none

      real(r8) :: stot(ntot_max)

!
      gas(ih2so4_g) = stot(ih2so4_g)
      gas(ihno3_g)  = stot(ihno3_g)
      gas(ihcl_g)   = stot(ihcl_g)
      gas(inh3_g)   = stot(inh3_g)
      gas(imsa_g)   = stot(imsa_g)

      return
      end subroutine make_gas_from_stot




      subroutine make_stot_from_aer(stot,ibin)
      use module_data_mosaic_main
      use module_data_mosaic_gas
      use module_data_mosaic_aero

      implicit none

      real(r8) :: stot(ntot_max)
      integer jp, iaer, istot, ibin

      if(jaerosolstate(ibin).eq.all_solid)then
        jp = jtotal
      elseif(jaerosolstate(ibin).eq.all_liquid)then
        jp = jliquid
      elseif(jaerosolstate(ibin).eq.mixed)then
        jp = jliquid
      else
        jp = jtotal
      endif

      do iaer = 1, ngas_ioa
        istot= ngas_ioa*ibin + iaer
        stot(istot) = aer(iaer,jp,ibin)	! nmol/m^3
      enddo

      return
      end subroutine make_stot_from_aer





      subroutine make_stot_from_gas(stot)
      use module_data_mosaic_main
      use module_data_mosaic_gas
      use module_data_mosaic_aero

      implicit none

      real(r8) :: stot(ntot_max)

!
      stot(ih2so4_g) = gas(ih2so4_g)
      stot(ihno3_g)  = gas(ihno3_g)
      stot(ihcl_g)   = gas(ihcl_g)
      stot(inh3_g)   = gas(inh3_g)
      stot(imsa_g)   = gas(imsa_g)

      return
      end subroutine make_stot_from_gas




!***********************************************************************
! part of ASTEM: computes fluxes over wet aerosols
!
! author: Rahul A. Zaveri
! update: may 2006
!-----------------------------------------------------------------------
      subroutine LSODE_flux_wet(ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      integer iv
      real(r8) :: XT



      call ions_to_electrolytes(jliquid,ibin,XT)  	! for water content calculation
      call compute_activities(ibin)

      if(water_a(ibin) .eq. 0.0)then
	write(6,*)'Water is zero in liquid phase'
	write(6,*)'Stopping in LSODE_flux_wet'
        stop
      endif



! default fluxes and other stuff
      do iv = 1, ngas_ioa
        sfc_a(iv)               = gas(iv)
        df_gas(iv,ibin)         = 0.0
        flux_s(iv,ibin)         = 0.0
        flux_l(iv,ibin)         = 0.0
      enddo

      if(gas(ih2so4_g) .gt. 1.e-14)then
        flux_l(ih2so4_g,ibin) = kg(ih2so4_g,ibin)*gas(ih2so4_g)
      endif

!-------------------------------------------------------------------
      if(electrolyte(jcaco3,jsolid,ibin) .gt. 0.0)then
        df_gas(inh3_g,ibin)  = 0.0
        df_gas(ihno3_g,ibin) = gas(ihno3_g)
        df_gas(ihcl_g,ibin)  = gas(ihcl_g)
      else
        sfc_a(inh3_g) = kel(inh3_g,ibin)*   &
                        gam_ratio(ibin)*mc(jc_nh4,ibin)*Keq_ll(3)/   &
                        (mc(jc_h,ibin)*Keq_ll(2)*Keq_gl(2))

        sfc_a(ihno3_g) = kel(ihno3_g,ibin)*activity(jhno3,ibin)/   &
                         Keq_gl(3)

        sfc_a(ihcl_g)  = kel(ihcl_g,ibin)*activity(jhcl,ibin)/   &
                         Keq_gl(4)

        df_gas(inh3_g,ibin) = gas(inh3_g) - sfc_a(inh3_g)
        df_gas(ihno3_g,ibin)= gas(ihno3_g)- sfc_a(ihno3_g)
        df_gas(ihcl_g,ibin) = gas(ihcl_g) - sfc_a(ihcl_g)

      endif

       flux_l(inh3_g,ibin) = kg(inh3_g,ibin)*df_gas(inh3_g,ibin)
       flux_l(ihno3_g,ibin)= kg(ihno3_g,ibin)*df_gas(ihno3_g,ibin)
       flux_l(ihcl_g,ibin) = kg(ihcl_g,ibin)*df_gas(ihcl_g,ibin)

      return
      end subroutine LSODE_flux_wet






!***********************************************************************
! part of ASTEM: computes gas-aerosol fluxes over dry aerosols
!
! author: Rahul A. Zaveri
! update: jan 2005
!-----------------------------------------------------------------------
      subroutine LSODE_flux_dry(ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      integer iv
      real(r8) :: XT, prod_nh4no3, prod_nh4cl, volatile_cl, volatile_no3




      call calculate_XT(ibin,jsolid,XT)

! default fluxes and other stuff
      do iv = 1, ngas_ioa
        sfc_a(iv)               = gas(iv)
        df_gas(iv,ibin)         = 0.0
        flux_s(iv,ibin)         = 0.0
        flux_l(iv,ibin)         = 0.0
        phi_volatile_s(iv,ibin) = 0.0
        phi_volatile_l(iv,ibin) = 0.0
      enddo

! h2so4
      if(gas(ih2so4_g) .gt. 1.e-14)then
        sfc_a(ih2so4_g)  = 0.0
        flux_s(ih2so4_g,ibin)       = kg(ih2so4_g,ibin)*gas(ih2so4_g)
        phi_volatile_s(ih2so4_g,ibin) = 1.0
      endif

!-----------------------------------------------------------------
! CASE 1: Sulfate-Rich Domain

      if(XT.lt.1.9999 .and. XT.ge.0.)then	! excess sulfate (acidic)

	call LSODE_flux_dry_case1(ibin)

        return
      endif

!-----------------------------------------------------------------
! CASE 2:  caco3 > 0 absorb all acids (and indirectly degas co2)

      if(electrolyte(jcaco3,jsolid,ibin) .gt. 0.0)then

        call LSODE_flux_dry_case2(ibin)

        return
      endif

!-------------------------------------------------------------------
! CASE 3: hno3 and hcl exchange may happen here and nh4cl may form/evaporate

      volatile_cl  = epercent(jnacl,jsolid,ibin) +   &
                     epercent(jcacl2,jsolid,ibin)


      if( volatile_cl .gt. 0.001 .and.   &
         (gas(ihno3_g).gt. 1.e-15 .or. gas(ih2so4_g).gt.1.e-15))then

        call LSODE_flux_dry_case3a(ibin)

        prod_nh4cl = max( (gas(inh3_g)*gas(ihcl_g) -Keq_sg(2)), 0.0d0) +   &
                     epercent(jnh4cl, jsolid,ibin)

        if(prod_nh4cl .gt. 0.0)then
          call LSODE_flux_dry_case3b(ibin)
        endif

        return
      endif

!-----------------------------------------------------------------
! CASE 4: nh4no3 or nh4cl or both may be active

      prod_nh4no3 = max( (gas(inh3_g)*gas(ihno3_g)-Keq_sg(1)), 0.0d0) +   &
                    epercent(jnh4no3,jsolid,ibin)
      prod_nh4cl  = max( (gas(inh3_g)*gas(ihcl_g) -Keq_sg(2)), 0.0d0) +   &
                    epercent(jnh4cl, jsolid,ibin)

      if(prod_nh4no3 .gt. 0.0 .or. prod_nh4cl .gt. 0.0)then
        call LSODE_flux_dry_case4(ibin)
        return
      endif

!-----------------------------------------------------------------
! CASE 5: condense h2so4 and degas hno3
      volatile_no3 = epercent(jnano3,jsolid,ibin) +   &
                     epercent(jcano3,jsolid,ibin)

      if(volatile_no3 .gt. 0.001 .and.   &
         gas(ih2so4_g).gt. 1.e-15 )then

        call LSODE_flux_dry_case5(ibin)

        return
      endif

!-------------------------------------------------------------------
! CASE 6: probably pure (nh4)2so4 particle.
        flux_s(ihno3_g,ibin)  = 0.0
        flux_s(ihcl_g,ibin)   = 0.0
        flux_s(inh3_g,ibin)   = min( kg(inh3_g,ibin)*gas(inh3_g),   &
                                 2.*flux_s(ih2so4_g,ibin) )
        return

      end subroutine LSODE_flux_dry

!----------------------------------------------------------------------







!***********************************************************************
! part of ASTEM: subroutines for flux_dry cases
!
! author: Rahul A. Zaveri
! update: jan 2005
!-----------------------------------------------------------------------
!
!
!
! CASE 1: Sulfate-Rich Domain
!
      subroutine LSODE_flux_dry_case1(ibin) ! TOUCH
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
      real(r8) :: beta_nh4, ppbv_nh3

      ppbv_nh3    = gas(inh3_g)/cair_mol_m3

      sfc_a(inh3_g)  = 0
      df_gas(inh3_g,ibin)   = gas(inh3_g)

      if(ppbv_nh3 .lt. 1.e-6)then
        flux_s(inh3_g,ibin) = 0.0
      else
        flux_s(inh3_g,ibin) = kg(inh3_g,ibin)*gas(inh3_g)
      endif


      return
      end subroutine LSODE_flux_dry_case1




! CASE 2:  caco3 > 0 absorb all acids (and indirectly degas co2)
!
      subroutine LSODE_flux_dry_case2(ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      real(r8) :: ppbv_hno3, ppbv_hcl


      ppbv_hno3 = gas(ihno3_g)/cair_mol_m3
      ppbv_hcl  = gas(ihcl_g)/cair_mol_m3

      if(ppbv_hno3 .gt. 1.e-6)then
        sfc_a(ihno3_g)       = 0.0
        df_gas(ihno3_g,ibin) = gas(ihno3_g)
        flux_s(ihno3_g,ibin) = kg(ihno3_g,ibin)*df_gas(ihno3_g,ibin)
      endif

      if(ppbv_hcl .gt. 1.e-6)then
        sfc_a(ihcl_g)       = 0.0
        df_gas(ihcl_g,ibin) = gas(ihcl_g)
        flux_s(ihcl_g,ibin) = kg(ihcl_g,ibin)*df_gas(ihcl_g,ibin)
      endif


      return
      end subroutine LSODE_flux_dry_case2





! CASE 3a: hno3 and hcl exchange may happen here
!
      subroutine LSODE_flux_dry_case3a(ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      real(r8) :: ppbv_hno3

! degas hcl from nacl or cacl2 by flux_s balance with 2 h2so4 and hno3

      ppbv_hno3 = gas(ihno3_g)/cair_mol_m3

      if(ppbv_hno3 .gt. 1.e-6)then
        flux_s(ihno3_g,ibin) = kg(ihno3_g,ibin)*gas(ihno3_g)
        phi_volatile_s(ihno3_g,ibin) = 1.0

        flux_s(ihcl_g,ibin) = - flux_s(ihno3_g,ibin)   &
                              - 2.*flux_s(ih2so4_g,ibin)
      endif

      if( df_gas(ihcl_g,ibin) .lt. 0. .and.   &
           epercent(jnacl,jsolid,ibin) .le. 0.001 )then
        df_gas(ihcl_g,ibin) = 0.0
        flux_s(ihcl_g,ibin) = 0.0
        flux_s(ihno3_g,ibin) = 0.0
      endif


      return
      end subroutine LSODE_flux_dry_case3a





! CASE 3b: nh4cl may form/evaporate here
!
      subroutine LSODE_flux_dry_case3b(ibin)	! TOUCH
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      integer iv, iactive_nh4cl
      real(r8) :: gnh3_hcl, pcnt_nh4cl,   &
           beta_nh4,   &
           a, b, c
! function
      real(r8) :: quadratic


!-------------------
! set default values for flags
      iactive_nh4cl  = 1

!-------------------
! compute diagnostic products and ratios
      gnh3_hcl    = gas(inh3_g)*gas(ihcl_g)

      phi_nh4cl_s = (gnh3_hcl - Keq_sg(2))/   &
                               Keq_sg(2)

      pcnt_nh4cl  = epercent(jnh4cl, jsolid,ibin)


!-------------------
! now determine if nh4cl is active or significant
! nh4cl
      if( abs(phi_nh4cl_s) .lt. rtol_eqb_astem )then
        iactive_nh4cl = 0
      elseif(gnh3_hcl.lt.Keq_sg(2) .and. pcnt_nh4cl.lt.0.001)then
        iactive_nh4cl = 0
      endif


! check the outcome
      if(iactive_nh4cl .eq. 0)then

        return
      endif



!-----------------
! nh4cl is active


      a =   kg(inh3_g,ibin)
      b = - kg(inh3_g,ibin)*gas(inh3_g)   &
          + kg(ihcl_g,ibin)*gas(ihcl_g)   &
          + 2.0*flux_s(ih2so4_g,ibin)
      c = -(kg(ihcl_g,ibin)*Keq_sg(2))

      sfc_a(inh3_g) = quadratic(a,b,c)
      sfc_a(ihcl_g) = Keq_sg(2) /sfc_a(inh3_g)

      df_gas(ihcl_g,ibin) = gas(ihcl_g) - sfc_a(ihcl_g)
      df_gas(inh3_g,ibin) = gas(inh3_g) - sfc_a(inh3_g)

      flux_s(inh3_g,ibin) = kg(inh3_g,ibin)*df_gas(inh3_g,ibin)
      flux_s(ihcl_g,ibin) = kg(ihcl_g,ibin)*df_gas(ihcl_g,ibin)

      if( df_gas(ihcl_g,ibin) .lt. 0. .and.   &
           epercent(jnh4cl,jsolid,ibin) .le. 0.001 )then
        df_gas(ihcl_g,ibin) = 0.0
        flux_s(ihcl_g,ibin) = 0.0
        flux_s(inh3_g,ibin) = min(2.0*flux_s(ih2so4_g,ibin),   &
                                  kg(inh3_g,ibin)*gas(inh3_g))
      endif


      return
      end subroutine LSODE_flux_dry_case3b









! Case 4: NH4NO3 and/or NH4Cl may be active
      subroutine LSODE_flux_dry_case4(ibin)	! TOUCH
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      integer iv, iactive_nh4no3, iactive_nh4cl, iactive
      real(r8) :: gnh3_hno3, gnh3_hcl, pcnt_nh4no3, pcnt_nh4cl,   &
           beta_nh4,   &
           a, b, c,   &
           flux_nh3_est, flux_nh3_max, ratio_flux
! function
      real(r8) :: quadratic


!-------------------
! set default values for flags
      iactive_nh4no3 = 1
      iactive_nh4cl  = 2

!-------------------
! compute diagnostic products and ratios
      gnh3_hno3   = gas(inh3_g)*gas(ihno3_g)
      gnh3_hcl    = gas(inh3_g)*gas(ihcl_g)

      phi_nh4no3_s = (gnh3_hno3 - Keq_sg(1))/   &
                                  Keq_sg(1)
      phi_nh4cl_s  = (gnh3_hcl - Keq_sg(2))/   &
                                  Keq_sg(2)

      pcnt_nh4no3 = epercent(jnh4no3,jsolid,ibin)
      pcnt_nh4cl  = epercent(jnh4cl, jsolid,ibin)


!-------------------
! now determine if nh4no3 and/or nh4cl are active or significant

! nh4no3
      if( abs(phi_nh4no3_s) .lt. rtol_eqb_astem )then
        iactive_nh4no3 = 0
      elseif(gnh3_hno3.lt.Keq_sg(1) .and. pcnt_nh4no3.lt.0.001)then
        iactive_nh4no3 = 0
      endif

! nh4cl
      if( abs(phi_nh4cl_s) .lt. rtol_eqb_astem )then
        iactive_nh4cl = 0
      elseif(gnh3_hcl.lt.Keq_sg(2) .and. pcnt_nh4cl.lt.0.001)then
        iactive_nh4cl = 0
      endif


      iactive = iactive_nh4no3 + iactive_nh4cl

! check the outcome
      if(iactive .eq. 0)then

        return
      endif

      goto (1,2,3),iactive

!---------------------------------
! only nh4no3 solid is active
1     call LSODE_flux_dry_case4a(ibin)

      return


!-----------------
! only nh4cl solid is active
2     call LSODE_flux_dry_case4b(ibin)

      return


!-----------------
! both nh4no3 and nh4cl are active
3      call LSODE_flux_dry_case4ab(ibin)

!      call LSODE_flux_dry_case4c(ibin)

      return
      end subroutine LSODE_flux_dry_case4









      subroutine LSODE_flux_dry_case4a(ibin) ! NH4NO3 solid
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
          + kg(ihno3_g,ibin)*gas(ihno3_g)   &
          + 2.0*flux_s(ih2so4_g,ibin)
      c = -(kg(ihno3_g,ibin)*Keq_sg(1))

      sfc_a(inh3_g)  = quadratic(a,b,c)
      sfc_a(ihno3_g) = Keq_sg(1)/sfc_a(inh3_g)

      df_gas(ihno3_g,ibin)=gas(ihno3_g)-sfc_a(ihno3_g)
      df_gas(inh3_g,ibin) =gas(inh3_g) -sfc_a(inh3_g)

      phi_volatile_s(ihno3_g,ibin)= phi_nh4no3_s
      phi_volatile_s(inh3_g,ibin) = phi_nh4no3_s

      flux_s(inh3_g,ibin)  = kg(inh3_g,ibin)*df_gas(inh3_g,ibin)
      flux_s(ihno3_g,ibin) = kg(ihno3_g,ibin)*df_gas(ihno3_g,ibin)

      if( df_gas(ihno3_g,ibin) .lt. 0..and.   &
           epercent(jnh4no3,jsolid,ibin) .le.0.001 )then
        df_gas(ihno3_g,ibin) = 0.0
        flux_s(ihno3_g,ibin) = 0.0
        flux_s(inh3_g,ibin) = min(2.0*flux_s(ih2so4_g,ibin),   &
                                  kg(inh3_g,ibin)*gas(inh3_g))
      endif


      return
      end subroutine LSODE_flux_dry_case4a






      subroutine LSODE_flux_dry_case4b(ibin) ! NH4Cl solid
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
          + kg(ihcl_g,ibin)*gas(ihcl_g)   &
          + 2.0*flux_s(ih2so4_g,ibin)
      c = -(kg(ihcl_g,ibin)*Keq_sg(2))

      sfc_a(inh3_g) = quadratic(a,b,c)
      sfc_a(ihcl_g) = Keq_sg(2) /sfc_a(inh3_g)

      df_gas(ihcl_g,ibin) = gas(ihcl_g)-sfc_a(ihcl_g)
      df_gas(inh3_g,ibin) = gas(inh3_g)-sfc_a(inh3_g)

      flux_s(inh3_g,ibin) = kg(inh3_g,ibin)*df_gas(inh3_g,ibin)
      flux_s(ihcl_g,ibin) = kg(ihcl_g,ibin)*df_gas(ihcl_g,ibin)

      if( df_gas(ihcl_g,ibin) .lt. 0. .and.   &
           epercent(jnh4cl,jsolid,ibin) .le. 0.001 )then
        df_gas(ihcl_g,ibin) = 0.0
        flux_s(ihcl_g,ibin) = 0.0
        flux_s(inh3_g,ibin) = min(2.0*flux_s(ih2so4_g,ibin),   &
                                  kg(inh3_g,ibin)*gas(inh3_g))
      endif


      return
      end subroutine LSODE_flux_dry_case4b














      subroutine LSODE_flux_dry_case4ab(ibin)	! NH4NO3 + NH4Cl (solid)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      real(r8) :: pcnt_nh4no3, pcnt_nh4cl,   &
           a, b, c,   &
           flux_nh3_est, flux_nh3_max, ratio_flux
! function
      real(r8) :: quadratic

      call LSODE_flux_dry_case4a(ibin)
      call LSODE_flux_dry_case4b(ibin)

      phi_volatile_s(inh3_g,ibin)  = max( abs(phi_nh4cl_s),   &
                                          abs(phi_nh4no3_s) )


! estimate nh3 flux and adjust hno3 and/or hcl if necessary

      flux_nh3_est = flux_s(ihno3_g,ibin)+flux_s(ihcl_g,ibin)+   &
                     2.*flux_s(ih2so4_g,ibin)
      flux_nh3_max = kg(inh3_g,ibin)*gas(inh3_g)


      if(flux_nh3_est .le. flux_nh3_max)then

        flux_s(inh3_g,ibin) = flux_nh3_est		! all ok - no adjustments needed

      else			! reduce hno3 and hcl flux_ses as necessary so that nh3 flux_s = flux_s_nh3_max

        flux_s(inh3_g,ibin) = flux_nh3_max
        flux_nh3_max = max(0.d0, (flux_nh3_max-2.0d0*flux_s(ih2so4_g,ibin)))
        ratio_flux          = flux_nh3_max/flux_nh3_est
        flux_s(ihno3_g,ibin)= flux_s(ihno3_g,ibin)*ratio_flux
        flux_s(ihcl_g,ibin) = flux_s(ihcl_g,ibin) *ratio_flux

        sfc_a(inh3_g)   = 0.0
        sfc_a(ihno3_g)  = gas(ihno3_g) -  			   &  ! recompute sfc_a(ihno3_g)
                          flux_s(ihno3_g,ibin)/kg(ihno3_g,ibin)
        sfc_a(ihcl_g)   = gas(ihcl_g)  -  			   &  ! recompute sfc_a(ihcl_g)
                          flux_s(ihcl_g,ibin)/kg(ihcl_g,ibin)

        df_gas(inh3_g,ibin) =gas(inh3_g) -sfc_a(inh3_g)
        df_gas(ihno3_g,ibin)=gas(ihno3_g)-sfc_a(ihno3_g)
        df_gas(ihcl_g,ibin) =gas(ihcl_g) -sfc_a(ihcl_g)

      endif

      return
      end subroutine LSODE_flux_dry_case4ab










! this is not used
      subroutine LSODE_flux_dry_case4c(ibin)	! NH4NO3 + NH4Cl (solid)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin
! local variables
      real(r8) :: pcnt_nh4no3, pcnt_nh4cl,   &
           a, b, c,   &
           flux_nh3_est, flux_nh3_max, ratio_flux
! function
      real(r8) :: quadratic



      a =   kg(inh3_g,ibin)
      b = - kg(inh3_g,ibin)*gas(inh3_g)   &
          + kg(ihno3_g,ibin)*gas(ihno3_g)   &
          + kg(ihcl_g,ibin)*gas(ihcl_g)   &
          + 2.0*flux_s(ih2so4_g,ibin)
      c = -( kg(ihno3_g,ibin)*Keq_sg(1) + kg(ihcl_g,ibin)*Keq_sg(2) )

      sfc_a(inh3_g)  = quadratic(a,b,c)
      sfc_a(ihno3_g) = Keq_sg(1)/sfc_a(inh3_g)
      sfc_a(ihcl_g)  = Keq_sg(2)/sfc_a(inh3_g)

      df_gas(ihno3_g,ibin)  = gas(ihno3_g) - sfc_a(ihno3_g)
      df_gas(ihcl_g,ibin)   = gas(ihcl_g)  - sfc_a(ihcl_g)
      df_gas(inh3_g,ibin)   = gas(inh3_g)  - sfc_a(inh3_g)

      flux_s(ihno3_g,ibin) = kg(ihno3_g,ibin)*df_gas(ihno3_g,ibin)
      flux_s(ihcl_g,ibin) = kg(ihcl_g,ibin)*df_gas(ihcl_g,ibin)


      if( df_gas(ihno3_g,ibin).lt.0..and.   &
           pcnt_nh4no3.lt.0.001  )then
        df_gas(ihno3_g,ibin) = 0.0
      endif

      if( df_gas(ihcl_g,ibin).lt.0..and.   &
           pcnt_nh4cl.lt.0.001 )then
        df_gas(ihcl_g,ibin) = 0.0
      endif



! estimate nh3 flux and adjust hno3 and/or hcl if necessary

      flux_nh3_est = flux_s(ihno3_g,ibin)+flux_s(ihcl_g,ibin)+   &
                     2.*flux_s(ih2so4_g,ibin)
      flux_nh3_max = kg(inh3_g,ibin)*gas(inh3_g)


      if(flux_nh3_est .le. flux_nh3_max .or. 1.eq.1)then

        flux_s(inh3_g,ibin) = flux_nh3_est		! all ok - no adjustments needed

      else			! reduce hno3 and hcl flux_ses as necessary so that nh3 flux_s = flux_s_nh3_max

        ratio_flux          = flux_nh3_max/flux_nh3_est
        flux_s(inh3_g,ibin) = flux_nh3_max
        flux_s(ihno3_g,ibin)= flux_s(ihno3_g,ibin)*ratio_flux
        flux_s(ihcl_g,ibin) = flux_s(ihcl_g,ibin) *ratio_flux

      endif

      return
      end subroutine LSODE_flux_dry_case4c














      subroutine LSODE_flux_dry_case5(ibin)
      use module_data_mosaic_aero

      implicit none

! subr arguments
      integer ibin


! just degas hno3 from nano3 or cano3 by flux balance with h2so4

      flux_s(ihno3_g,ibin) = -2.*flux_s(ih2so4_g,ibin)

      if(epercent(jnano3,jsolid,ibin) .lt. 0.001 .or.   &
         epercent(jcano3,jsolid,ibin) .lt. 0.001)then
        flux_s(ihno3_g,ibin) = 0.0
      endif


      return
      end subroutine LSODE_flux_dry_case5


































