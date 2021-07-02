! growth: (wet mass)/(dry mass)
      subroutine parametric_analysis_mode
      use module_data_mosaic_main
      use module_data_mosaic_gas
      use module_data_mosaic_aero
      use module_data_mosaic_cloud

      implicit none

      integer jdum, idum, ibin, iter_kelvin
      real(r8) :: solvated_water_mass, water_mass, wet_mass, xt



      write(6,*)'    '
      write(6,*)'         Parametric Analysis Begins'

      RH_pc = RH				! RH(%)
      aH2O = 0.01*RH_pc				! aH2O (aerosol water activity)
      P_atm = pr_atm				! P(atm)
      T_K = te					! T(K)
      cair_mol_m3 = cair_molm3			! air conc in mol/m3
      cair_mol_cc = cair_mol_m3*1.e-6		! air conc in mol/cc

!      call print_aer				! inital concentrations
        write(66,300)

300    format('   RH(%)     aH2O        MGF       sigma',   &
              '   iterations caso4',   &
              '      na2so4     nahso4     nh4so4     lvcite',   &
              '     nh4hso4    cano3      nano3     nh4no3',   &
              '      cacl2       nacl      nh4cl')
!------------------------------------------------------------------
!
! parameter loop begins...

      do 100 idum = 10, 100, 10
!      do 50 jdum = 1, 10

        RH_pc = float(idum) ! + float(jdum - 1)/10.
        RH = RH_pc
        aH2O = 0.01*RH_pc

        aH2O = min(aH2O, 0.99_r8)

        call load_mosaic_parameters		! sets up indices and other stuff once per simulation

        call update_thermodynamic_constants	! update temperature dependent constants

        call initialize_mosaic_variables

        call map_mosaic_species_BOX(0)

        ibin = 1
        iter_mesa(ibin) = 0

        aH2O_a(ibin) = aH2O			! initialize

!        call cputime(5)
        call MTEM_compute_log_gamZ			! function of RH
!        call cputime(3)

!        call cputime(5)
        call conform_electrolytes(jtotal,ibin,XT)
        call conform_aerosol_number(ibin)   		! adjusts number conc so that it conforms with bin mass and diameter

        call do_full_deliquescence(ibin)
        call ions_to_electrolytes(jliquid,ibin,XT) ! for Li and Lu surface tension
        call compute_activities(ibin)

!        call aerosol_phase_state(ibin)


!        solvated_water_mass=6.*electrolyte(jcacl2,jsolid,ibin)*1.e-9*18.

        solvated_water_mass = 0.0

!        water_mass = water_a(ibin)*1.e-3 + solvated_water_mass ! [g/cc(air)]
!        wet_mass = mass_dry_a(ibin) + water_mass		   ! [g/cc(air)]
!        growth_factor(ibin) = wet_mass/mass_dry_a(ibin)


!        r_na = aer(ina_a,jsolid,ibin)/aer(ina_a,jtotal,ibin)
!        r_no3= aer(ino3_a,jsolid,ibin)/aer(ino3_a,jtotal,ibin)
!        r_cl = aer(icl_a,jsolid,ibin)/aer(icl_a,jtotal,ibin)
!        r_nh4= aer(inh4_a,jsolid,ibin)/aer(inh4_a,jtotal,ibin)


!        write(66,228) RH_pc, aH2O_a(ibin), growth_factor(ibin),   &
!        sigma_soln(ibin),   &
!        iter_kelvin,   &
!        electrolyte(jcaso4, jsolid,ibin),   &
!        electrolyte(jna2so4,jsolid,ibin),   &
!        electrolyte(jnahso4,jsolid,ibin),   &
!        electrolyte(jnh4so4,jsolid,ibin),   &
!        electrolyte(jlvcite,jsolid,ibin),   &
!        electrolyte(jnh4hso4,jsolid,ibin),   &
!        electrolyte(jcano3, jsolid,ibin),   &
!        electrolyte(jnano3, jsolid,ibin),   &
!        electrolyte(jnh4no3,jsolid,ibin),   &
!        electrolyte(jcacl2, jsolid,ibin),   &
!        electrolyte(jnacl,  jsolid,ibin),   &
!        electrolyte(jnh4cl, jsolid,ibin)


        write(6,230)aH2O, gam_ratio(ibin)
230     format(f8.2, 2x, f11.4)

228     format(f7.1, 2x, f9.4,2x,f10.3, 1x, f10.4, 2x, i5, 12(2x, f9.3))
229     format(f7.3, 5(2x,f9.3), 2x, i5)

!50    continue
100   continue	! time loop

!      call cputime(-1)


!------------------------------------------------------------------

      write(6,*)'   '
      write(6,*)'         End of Parametric Analysis'
!      write(6,*)'******************************************************'

      return
      end subroutine parametric_analysis_mode




