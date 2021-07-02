! Aerosol-phase Species
! subroutine for printing output at iprint time steps
!
!--------------------------------------------------------------------
!
      subroutine print_aer(idum)
  
      use module_data_mosaic_main
      use module_data_mosaic_gas
      use module_data_mosaic_aero
      use module_data_mosaic_cloud

	use module_data_mosaic_asect, only:  &
	    dcen_sect, dlo_sect, dhi_sect,  &
	    itype_of_ibin, itype_md1_of_itype, itype_md2_of_itype,  &
	    nsize_aer, ntype_aer, ntype_md1_aer, ntype_md2_aer,  &
	    xcut_atype_md1, xcut_atype_md2

      implicit none

      integer idum, lun
      integer l, ibin, lbin, je, iv, ifix
      real(r8) :: timeofday, time_dum, cppb(ntot_max), conv1, conv2, conv3, conv4
      real(r8) :: dlogDp
      real(r8) :: SOA_tot, SOA_tot1, SOA_tot2, SOA_tot3, SOA_tot4, SOA_tot5, SOA_tot6, &
	              SOA_tot7 ,SOA_tot8, SOA_tot9, SOA_tot10, SOA_tot11, SOA_tot12, SOA_tot13, &
	              OC_tot, &
		          avol_tot, anum_tot, Dp_cutoff



      conv1 = 1.e15/avogad	 ! converts (molec/cc) to (nmol/m^3)
      conv2 = 1./conv1		 ! converts (nmol/m^3) to (molec/cc)
      conv3 = conv2/cair_mlc*ppb ! converts (nmol/m^3) to (ppbv)
!      conv4 = 1.e-3		 ! converts (nmol/m^3) to (umol/m^3)
      conv4 = 1.0		 ! converts (nmol/m^3) to (nmol/m^3)



!      dlogDp = 1.0/64	         ! 110 SMPS bins
      dlogDp = 0.055243448       ! 44 FIMS bins
!     dlogDp = 1.0		 ! monodisperse

      do l=1,ngas_max
        cppb(l) = (cnn(l)/cair_mlc)*ppb	! converting (molecules/cc) to (ppb)
      enddo

      if(idum .eq. 0)then
        time_dum = 0.0
        timeofday = time_UTC_beg
      else
        time_dum = time_hrs
        timeofday = time_UTC
      endif


      SOA_tot1  = 0.0
      SOA_tot2  = 0.0
      SOA_tot3  = 0.0
      SOA_tot4  = 0.0
      SOA_tot5  = 0.0
      SOA_tot6  = 0.0
	  SOA_tot7  = 0.0
	  SOA_tot8  = 0.0
	  SOA_tot9  = 0.0
	  SOA_tot10 = 0.0
	  SOA_tot11 = 0.0
	  SOA_tot12 = 0.0
	  SOA_tot13 = 0.0
      OC_tot   = 0.0
      avol_tot = 0.0
      anum_tot = 0.0

      
      do 30 ibin = 1, nbin_a
	    SOA_tot1  = SOA_tot1 + aer(icn3_a,jtotal,ibin)*1.e-3*mw_aer_mac(icn3_a) ! ug/m^3
            SOA_tot2  = SOA_tot2 + aer(icn2_a,jtotal,ibin)*1.e-3*mw_aer_mac(icn2_a) ! ug/m^3
            SOA_tot3  = SOA_tot3 + aer(icn1_a,jtotal,ibin)*1.e-3*mw_aer_mac(icn1_a) ! ug/m^3
            SOA_tot4  = SOA_tot4 + aer(ic0_a,jtotal,ibin)*1.e-3*mw_aer_mac(ic0_a) ! ug/m^3
            SOA_tot5  = SOA_tot5 + aer(ic1_a,jtotal,ibin)*1.e-3*mw_aer_mac(ic1_a) ! ug/m^3
            SOA_tot6  = SOA_tot6 + aer(ic2_a,jtotal,ibin)*1.e-3*mw_aer_mac(ic2_a) ! ug/m^3
	    SOA_tot7  = SOA_tot7 + aer(ic3_a,jtotal,ibin)*1.e-3*mw_aer_mac(ic3_a) ! ug/m^3
	    SOA_tot8  = SOA_tot8 + aer(ic4_a,jtotal,ibin)*1.e-3*mw_aer_mac(ic4_a) ! ug/m^3
	    SOA_tot9  = SOA_tot9 + aer(ic5_a,jtotal,ibin)*1.e-3*mw_aer_mac(ic5_a) ! ug/m^3
	    SOA_tot10 = SOA_tot10 + aer(ic6_a,jtotal,ibin)*1.e-3*mw_aer_mac(ic6_a) ! ug/m^3
	    SOA_tot11 = SOA_tot11 + aer(ic7_a,jtotal,ibin)*1.e-3*mw_aer_mac(ic7_a) ! ug/m^3
	    SOA_tot12 = SOA_tot12 + aer(ic8_a,jtotal,ibin)*1.e-3*mw_aer_mac(ic8_a) ! ug/m^3
            SOA_tot13 = SOA_tot13 + aer(ic9_a,jtotal,ibin)*1.e-3*mw_aer_mac(ic9_a) ! ug/m^3

	    OC_tot   = OC_tot   + aer(ioc_a,jtotal,ibin)  *1.e-3	   ! ug/m^3
            avol_tot = avol_tot + (num_a(ibin)*0.523598776*(Dp_dry_a(ibin)*1.e4)**3)  ! um^3/cc(air)
	    anum_tot = anum_tot + num_a(ibin)
30    continue

        SOA_tot = SOA_tot1 + SOA_tot2 + SOA_tot3 + SOA_tot4 + SOA_tot5 + SOA_tot6 &
					+ SOA_tot7 + SOA_tot8 + SOA_tot9 + SOA_tot10 + SOA_tot11 + SOA_tot12 &
					+ SOA_tot13
       
      do 20 ibin = 1, nbin_a

        if(mc(jc_h,ibin) .gt. 0.0)then
          pH(ibin) = -log10(mc(jc_h,ibin))
        else
          pH(ibin) = 0.0
        endif


! write to "soa.aer.txt" file
    if (iwrite_aer_species > 0 .and. it .eq. nstep)   &
        write(lun_species,101)    &
        time_dum, 		  &  ! hours since start of simulation
        Dp_dry_a(ibin)*1.e4*1.e3, &  ! nm (same as Dp_dry_a)
!        dcen_sect(ibin,1)*1.0e4, &
        num_a(ibin)/dlogDp,	  &  ! #/cc(air)
        (num_a(ibin)*0.523598776*(Dp_dry_a(ibin)*1.e4)**3)/dlogDp, &  ! um^3/cc(air)
        SOA_tot,                  &  ! ug/m3
        aer(iso4_a,jtotal,ibin),  &  ! nmol/m3
        aer(inh4_a,jtotal,ibin),  &  ! nmol/m3
        aer(inh4_a,jtotal,ibin),  &  ! nmol/m3
        water_a(ibin)*1.e12/18.0     ! nmol/m3

    write(lun_species,101)    &
    time_dum, 		  &  ! hours since start of simulation
        Dp_dry_a(ibin)*1.e4*1.e3, &  ! nm (same as Dp_dry_a)
    !    dcen_sect(ibin,1)*1.0e4, &
        num_a(ibin)/dlogDp,	  &  ! #/cc(air)
        (num_a(ibin)*0.523598776*(Dp_dry_a(ibin)*1.e4)**3)/dlogDp, &  ! um^3/cc(air)
        SOA_tot,                  &  ! ug/m3
        aer(iso4_a,jtotal,ibin),  &  ! nmol/m3
            aer(inh4_a,jtotal,ibin),  &  ! nmol/m3
        aer(inh4_a,jtotal,ibin),  &  ! nmol/m3
        water_a(ibin)*1.e12/18.0,  &  ! nmol/m3
        num_a(ibin), &
        cnn(kcn3), &
        cnn(kcn2), &
        cnn(kcn1), &
        cnn(kc0), &
        cnn(kc1), &
        cnn(kc2), &
        cnn(kc3), &
        cnn(kc4), &
        cnn(kc5), &
        cnn(kc6), &
        cnn(kc7), &
        cnn(kc8), &
        cnn(kc9), &
        cnn(kisop), &
        cnn(kOH), &
        aer(icn3_a, jtotal, 1)/conv1, &
        aer(icn2_a, jtotal, 1)/conv1, &
        aer(icn1_a, jtotal, 1)/conv1, &
        aer(ic0_a, jtotal, 1)/conv1, &
        aer(ic1_a, jtotal, 1)/conv1, &
        aer(ic2_a, jtotal, 1)/conv1, &
        aer(ic3_a, jtotal, 1)/conv1, &
        aer(ic4_a, jtotal, 1)/conv1, &
        aer(ic5_a, jtotal, 1)/conv1, &
        aer(ic6_a, jtotal, 1)/conv1, &
        aer(ic7_a, jtotal, 1)/conv1, &
        aer(ic8_a, jtotal, 1)/conv1, &
        aer(ic9_a, jtotal, 1)/conv1

20    continue

101   format(f8.4,2x,72(2x,e13.6))


      return
      end subroutine print_aer

