! 12/24/2012 raz - moved wall loss / dilution subroutine to this file
!
      subroutine aerchemistry(t_in, t_out)
  
      use module_data_mosaic_main
      use module_data_mosaic_gas
      use module_data_mosaic_aero
      use module_sect_iface, only:  sectional_interface_1

      implicit none

!   subr arguments
      real(r8) :: t_in, t_out
!   local variables
      integer ibin, iaer, noffset
      real(r8) :: cnn_sv1(ntot_used)
      real(r8) :: dtchem


      dtchem = t_out - t_in

! do first-order wall loss or dilution here
!      call wall_loss( dtchem )
!     call dilution( dtchem )

      cnn_sv1(1:ntot_used) = cnn(1:ntot_used)


! calculate gas-aerosol exchange over timestep dtchem
      ! (during this calculation there is no transfer of particles between bins)


      
      call mosaic_box_aerchemistry( dtchem )


! for sectional framework, calculate
!    transfer of particles betweens due to growth/shrinkage
!    new particle nucleation (optional)
!    particle coagulation (optional)
      if (msize_framework == msectional) then
         call sectional_interface_1( dtchem, cnn_sv1 )

!     else if (msize_framework == mmodal) then
! for modal framework, do similar calculations
! (not yet implemented)
      end if



! add source of particles in the smallest bin
!        do ibin = 1, 1

!          noffset = ngas_max + naer_tot*(ibin - 1)
!          cnn(noffset + knum_a)   = cnn_sv1(noffset + knum_a)
!          cnn(noffset + kdpdry_a) = cnn_sv1(noffset + kdpdry_a)	! dry diameter (micron)
!          cnn(noffset + kjhyst_a) = cnn_sv1(noffset + kjhyst_a)
!          cnn(noffset + kwater_a) = cnn_sv1(noffset + kwater_a)		! kg/m^3

!          do iaer = 1, naer
!            cnn(noffset+kwater_a+iaer) = cnn_sv1(noffset+kwater_a+iaer)	! molec/cc
!          enddo
!
!        enddo

!            write(6,*)'-------------------------------------'
!            PAUSE



      return
      end subroutine aerchemistry




!***********************************************************************
! applies first-order wall loss to number and mass
!
! author: Rahul A. Zaveri
! update: dec 2012
!-----------------------------------------------------------------------
      subroutine wall_loss(dtchem)
      use module_data_mosaic_main
      use module_data_mosaic_gas
      use module_data_mosaic_aero
      use module_data_mosaic_asect

      implicit none
! subr arguments
      real(r8) :: dtchem, kwall
! local variables
      integer iv, iaer, ibin, isize, itype
! local variables
	real(r8) :: kloss_aer, kloss_gas, num_a_bkg(nbin_a), aer_so4_bkg(nbin_a), aer_nh4_bkg(nbin_a), aer_oc_bkg(nbin_a)
	real(r8) :: dum_factor
! function
	real(r8) :: lamda_wall

      P_atm = pr_atm				! P(atm)
      T_K = te					! T(K)
      cair_mol_m3 = cair_molm3		! air conc in mol/m3
      cair_mol_cc = cair_mol_m3*1.e-6	! air conc in mol/cc

      call initialize_mosaic_variables
      call map_mosaic_species_BOX(0)



! initialize background concentrations
        num_a_bkg( 1) =      29858.00000   ! #/cc
!        num_a_bkg( 2) =      0.00000   ! #/cc
!        num_a_bkg( 3) =      0.00000   ! #/cc
!        num_a_bkg( 4) =      0.00000   ! #/cc
!        num_a_bkg( 5) =      0.00000   ! #/cc
!        num_a_bkg( 6) =      0.00003   ! #/cc
!        num_a_bkg( 7) =      0.00026   ! #/cc
!        num_a_bkg( 8) =      0.00170   ! #/cc
!        num_a_bkg( 9) =      0.00918   ! #/cc
!        num_a_bkg(10) =      0.04141   ! #/cc
!        num_a_bkg(11) =      0.15599   ! #/cc
!        num_a_bkg(12) =      0.49101   ! #/cc
!        num_a_bkg(13) =      1.29137   ! #/cc
!        num_a_bkg(14) =      2.83787   ! #/cc
!        num_a_bkg(15) =      5.43272   ! #/cc
!        num_a_bkg(16) =      7.79294   ! #/cc
!        num_a_bkg(17) =     10.48545   ! #/cc
!        num_a_bkg(18) =     11.69284   ! #/cc
!        num_a_bkg(19) =     13.53134   ! #/cc
!        num_a_bkg(20) =     15.27858   ! #/cc
!        num_a_bkg(21) =     18.18119   ! #/cc
!        num_a_bkg(22) =     22.93771   ! #/cc
!        num_a_bkg(23) =     29.85795   ! #/cc
!        num_a_bkg(24) =     35.51799   ! #/cc
!        num_a_bkg(25) =     40.55044   ! #/cc
!        num_a_bkg(26) =     43.37145   ! #/cc
!        num_a_bkg(27) =     43.55781   ! #/cc
!        num_a_bkg(28) =     41.32913   ! #/cc
!        num_a_bkg(29) =     32.73955   ! #/cc
!        num_a_bkg(30) =     28.92083   ! #/cc
!        num_a_bkg(31) =     25.36255   ! #/cc
!        num_a_bkg(32) =     21.88892   ! #/cc
!        num_a_bkg(33) =     18.53446   ! #/cc
!        num_a_bkg(34) =     15.69874   ! #/cc
!        num_a_bkg(35) =     11.48562   ! #/cc
!        num_a_bkg(36) =      8.64723   ! #/cc
!        num_a_bkg(37) =      5.72210   ! #/cc
!        num_a_bkg(38) =      3.00147   ! #/cc
        !num_a_bkg(39) =      1.22834   ! #/cc
        !num_a_bkg(40) =      0.46336   ! #/cc
        !num_a_bkg(41) =      0.17894   ! #/cc
        !num_a_bkg(42) =      0.05340   ! #/cc
        !num_a_bkg(43) =      0.01381   ! #/cc
        !num_a_bkg(44) =      0.00310   ! #/cc

        aer_so4_bkg( 1) =  0.22355E+03 !        aer_so4_bkg( 1) =  0.00000E+00   ! nmol/m3
!        aer_so4_bkg( 2) =  0.00000E+00   ! nmol/m3
!        aer_so4_bkg( 3) =  0.00000E+00   ! nmol/m3
!        aer_so4_bkg( 4) =  0.00000E+00   ! nmol/m3
!        aer_so4_bkg( 5) =  0.00000E+00   ! nmol/m3
!        aer_so4_bkg( 6) =  0.34199E-11   ! nmol/m3
!        aer_so4_bkg( 7) =  0.43411E-10   ! nmol/m3
!        aer_so4_bkg( 8) =  0.41572E-09   ! nmol/m3
!        aer_so4_bkg( 9) =  0.32880E-08   ! nmol/m3
!        aer_so4_bkg(10) =  0.21723E-07   ! nmol/m3
!        aer_so4_bkg(11) =  0.11985E-06   ! nmol/m3
!        aer_so4_bkg(12) =  0.55255E-06   ! nmol/m3
!        aer_so4_bkg(13) =  0.21284E-05   ! nmol/m3
!        aer_so4_bkg(14) =  0.68506E-05   ! nmol/m3
!        aer_so4_bkg(15) =  0.19208E-04   ! nmol/m3
!        aer_so4_bkg(16) =  0.40355E-04   ! nmol/m3
!        aer_so4_bkg(17) =  0.79527E-04   ! nmol/m3
!        aer_so4_bkg(18) =  0.12989E-03   ! nmol/m3
!        aer_so4_bkg(19) =  0.22016E-03   ! nmol/m3
!        aer_so4_bkg(20) =  0.36409E-03   ! nmol/m3
!        aer_so4_bkg(21) =  0.63456E-03   ! nmol/m3
!        aer_so4_bkg(22) =  0.11726E-02   ! nmol/m3
!        aer_so4_bkg(23) =  0.22355E-02   ! nmol/m3
!        aer_so4_bkg(24) =  0.38949E-02   ! nmol/m3
!        aer_so4_bkg(25) =  0.65128E-02   ! nmol/m3
!        aer_so4_bkg(26) =  0.10203E-01   ! nmol/m3
!        aer_so4_bkg(27) =  0.15007E-01   ! nmol/m3
!        aer_so4_bkg(28) =  0.20855E-01   ! nmol/m3
!        aer_so4_bkg(29) =  0.24197E-01   ! nmol/m3
!        aer_so4_bkg(30) =  0.31307E-01   ! nmol/m3
!        aer_so4_bkg(31) =  0.40211E-01   ! nmol/m3
!        aer_so4_bkg(32) =  0.50829E-01   ! nmol/m3
!        aer_so4_bkg(33) =  0.63037E-01   ! nmol/m3
!        aer_so4_bkg(34) =  0.78201E-01   ! nmol/m3
!        aer_so4_bkg(35) =  0.83797E-01   ! nmol/m3
!        aer_so4_bkg(36) =  0.92403E-01   ! nmol/m3
!        aer_so4_bkg(37) =  0.89556E-01   ! nmol/m3
!        aer_so4_bkg(38) =  0.68802E-01   ! nmol/m3
        !aer_so4_bkg(39) =  0.41240E-01   ! nmol/m3
        !aer_so4_bkg(40) =  0.22785E-01   ! nmol/m3
        !aer_so4_bkg(41) =  0.12887E-01   ! nmol/m3
        !aer_so4_bkg(42) =  0.56329E-02   ! nmol/m3
        !aer_so4_bkg(43) =  0.21336E-02   ! nmol/m3
        !aer_so4_bkg(44) =  0.70147E-03   ! nmol/m3

        aer_nh4_bkg( 1) =  0.44710E+03   !        aer_nh4_bkg( 1) =  0.00000E+00   ! nmol/m3
!        aer_nh4_bkg( 2) =  0.00000E+00   ! nmol/m3
!        aer_nh4_bkg( 3) =  0.00000E+00   ! nmol/m3
!        aer_nh4_bkg( 4) =  0.00000E+00   ! nmol/m3
!        aer_nh4_bkg( 5) =  0.00000E+00   ! nmol/m3
!        aer_nh4_bkg( 6) =  0.68399E-11   ! nmol/m3
!        aer_nh4_bkg( 7) =  0.86822E-10   ! nmol/m3
!        aer_nh4_bkg( 8) =  0.83145E-09   ! nmol/m3
!        aer_nh4_bkg( 9) =  0.65760E-08   ! nmol/m3
!        aer_nh4_bkg(10) =  0.43446E-07   ! nmol/m3
!        aer_nh4_bkg(11) =  0.23970E-06   ! nmol/m3
!        aer_nh4_bkg(12) =  0.11051E-05   ! nmol/m3
!        aer_nh4_bkg(13) =  0.42569E-05   ! nmol/m3
!        aer_nh4_bkg(14) =  0.13701E-04   ! nmol/m3
!        aer_nh4_bkg(15) =  0.38416E-04   ! nmol/m3
!        aer_nh4_bkg(16) =  0.80711E-04   ! nmol/m3
!        aer_nh4_bkg(17) =  0.15905E-03   ! nmol/m3
!        aer_nh4_bkg(18) =  0.25978E-03   ! nmol/m3
!        aer_nh4_bkg(19) =  0.44031E-03   ! nmol/m3
!        aer_nh4_bkg(20) =  0.72817E-03   ! nmol/m3
!        aer_nh4_bkg(21) =  0.12691E-02   ! nmol/m3
!        aer_nh4_bkg(22) =  0.23451E-02   ! nmol/m3
!        aer_nh4_bkg(23) =  0.44710E-02   ! nmol/m3
!        aer_nh4_bkg(24) =  0.77897E-02   ! nmol/m3
!        aer_nh4_bkg(25) =  0.13026E-01   ! nmol/m3
!        aer_nh4_bkg(26) =  0.20405E-01   ! nmol/m3
!        aer_nh4_bkg(27) =  0.30014E-01   ! nmol/m3
!        aer_nh4_bkg(28) =  0.41711E-01   ! nmol/m3
!        aer_nh4_bkg(29) =  0.48395E-01   ! nmol/m3
!        aer_nh4_bkg(30) =  0.62613E-01   ! nmol/m3
!        aer_nh4_bkg(31) =  0.80423E-01   ! nmol/m3
!        aer_nh4_bkg(32) =  0.10166E+00   ! nmol/m3
!        aer_nh4_bkg(33) =  0.12607E+00   ! nmol/m3
!        aer_nh4_bkg(34) =  0.15640E+00   ! nmol/m3
!        aer_nh4_bkg(35) =  0.16759E+00   ! nmol/m3
!        aer_nh4_bkg(36) =  0.18481E+00   ! nmol/m3
!        aer_nh4_bkg(37) =  0.17911E+00   ! nmol/m3
!        aer_nh4_bkg(38) =  0.13760E+00   ! nmol/m3
        !aer_nh4_bkg(39) =  0.82480E-01   ! nmol/m3
        !aer_nh4_bkg(40) =  0.45570E-01   ! nmol/m3
        !aer_nh4_bkg(41) =  0.25775E-01   ! nmol/m3
        !aer_nh4_bkg(42) =  0.11266E-01   ! nmol/m3
        !aer_nh4_bkg(43) =  0.42672E-02   ! nmol/m3
        !aer_nh4_bkg(44) =  0.14029E-02   ! nmol/m3

        aer_oc_bkg( 1) =  0.19033E+04    !        aer_oc_bkg( 1) =  0.00000E+00   ! ng/m3
!        aer_oc_bkg( 2) =  0.00000E+00   ! ng/m3
!        aer_oc_bkg( 3) =  0.00000E+00   ! ng/m3
!        aer_oc_bkg( 4) =  0.00000E+00   ! ng/m3
!        aer_oc_bkg( 5) =  0.00000E+00   ! ng/m3
!        aer_oc_bkg( 6) =  0.29117E-08   ! ng/m3
!        aer_oc_bkg( 7) =  0.36960E-07   ! ng/m3
!        aer_oc_bkg( 8) =  0.35395E-06   ! ng/m3
!        aer_oc_bkg( 9) =  0.27994E-05   ! ng/m3
!        aer_oc_bkg(10) =  0.18495E-04   ! ng/m3
!        aer_oc_bkg(11) =  0.10204E-03   ! ng/m3
!        aer_oc_bkg(12) =  0.47044E-03   ! ng/m3
!        aer_oc_bkg(13) =  0.18121E-02   ! ng/m3
!        aer_oc_bkg(14) =  0.58326E-02   ! ng/m3
!        aer_oc_bkg(15) =  0.16354E-01   ! ng/m3
!        aer_oc_bkg(16) =  0.34359E-01   ! ng/m3
!        aer_oc_bkg(17) =  0.67710E-01   ! ng/m3
!        aer_oc_bkg(18) =  0.11059E+00   ! ng/m3
!        aer_oc_bkg(19) =  0.18744E+00   ! ng/m3
!        aer_oc_bkg(20) =  0.30998E+00   ! ng/m3
!        aer_oc_bkg(21) =  0.54027E+00   ! ng/m3
!        aer_oc_bkg(22) =  0.99831E+00   ! ng/m3
!        aer_oc_bkg(23) =  0.19033E+01   ! ng/m3
!        aer_oc_bkg(24) =  0.33161E+01   ! ng/m3
!        aer_oc_bkg(25) =  0.55450E+01   ! ng/m3
!        aer_oc_bkg(26) =  0.86864E+01   ! ng/m3
!        aer_oc_bkg(27) =  0.12777E+02   ! ng/m3
!        aer_oc_bkg(28) =  0.17756E+02   ! ng/m3
!        aer_oc_bkg(29) =  0.20602E+02   ! ng/m3
!        aer_oc_bkg(30) =  0.26654E+02   ! ng/m3
!        aer_oc_bkg(31) =  0.34236E+02   ! ng/m3
!        aer_oc_bkg(32) =  0.43276E+02   ! ng/m3
!        aer_oc_bkg(33) =  0.53670E+02   ! ng/m3
!        aer_oc_bkg(34) =  0.66580E+02   ! ng/m3
!        aer_oc_bkg(35) =  0.71345E+02   ! ng/m3
!        aer_oc_bkg(36) =  0.78671E+02   ! ng/m3
!        aer_oc_bkg(37) =  0.76248E+02   ! ng/m3
!        aer_oc_bkg(38) =  0.58578E+02   ! ng/m3
        !aer_oc_bkg(39) =  0.35112E+02   ! ng/m3
        !aer_oc_bkg(40) =  0.19399E+02   ! ng/m3
        !aer_oc_bkg(41) =  0.10972E+02   ! ng/m3
        !aer_oc_bkg(42) =  0.47958E+01   ! ng/m3
        !aer_oc_bkg(43) =  0.18165E+01   ! ng/m3
        !aer_oc_bkg(44) =  0.59723E+00   ! ng/m3


      do ibin = 1, nbin_a

        if (msize_framework == msectional) then
           isize = isize_of_ibin(ibin)
           itype = itype_of_ibin(ibin)
           Dp_dry_a(ibin) = dcen_sect(isize,itype)
        end if

!        kloss_aer = 4.0e-6 ! 0.35*lamda_wall(Dp_dry_a(ibin)*1.e7)  ! 1/s  -- 9/19/2012
!	 kloss_aer = 7.0e-6 ! 12/7/2012
!	 kloss_aer = 1.e-5  ! 12/19/2012
!        kloss_aer = 0.0 ! 0.35*lamda_wall(Dp_dry_a(ibin)*1.e7)  ! 1/s  -- Alla pyrene evap
!	  write(6,*)Dp_dry_a(ibin)*1.e7, dlog(Dp_dry_a(ibin)*1.e7), kloss_aer


!        kloss_aer = 0.0 ! 1/s	! default



        kloss_aer = 0.0*1.30e-4  ! 1/s



!	write(6,*)'kloss_aer = ', kloss_aer

        do iaer = 1, naer

          if(iaer .eq. iso4_a)then
            aer(iaer,jtotal,ibin)  = aer_so4_bkg(ibin) + (aer(iaer,jtotal,ibin)-aer_so4_bkg(ibin))*exp(-kloss_aer*dtchem)
	  elseif(iaer .eq. inh4_a)then
            aer(iaer,jtotal,ibin)  = aer_nh4_bkg(ibin) + (aer(iaer,jtotal,ibin)-aer_nh4_bkg(ibin))*exp(-kloss_aer*dtchem)
	  elseif(iaer .eq. ioc_a)then
	    aer(iaer,jtotal,ibin)  = aer_oc_bkg(ibin) + (aer(iaer,jtotal,ibin)-aer_oc_bkg(ibin))*exp(-kloss_aer*dtchem)
	  else
            aer(iaer,jtotal,ibin)  = aer(iaer,jtotal,ibin)*exp(-kloss_aer*dtchem)
	  endif


        enddo

        num_a(ibin) = num_a_bkg(ibin) + (num_a(ibin) - num_a_bkg(ibin))*exp(-kloss_aer*dtchem)

      enddo


! source of small particles
	dum_factor = 0.
        num_a(1)            = 5.32479 * dum_factor			! #/cc
        aer(ioc_a,jtotal,1) = 0.46535E-02 * dum_factor	! ng/m3




	kloss_gas = 0.0*1.0* 1.e-6

	do iv = icn3_g, ngas_volatile
	  gas(iv) = gas(iv)*exp(-kloss_gas*dtchem)
	enddo


      call map_mosaic_species_BOX(1)

      return
      end subroutine wall_loss




      function lamda_wall(Dp)		! Dp in nm
      use module_data_mosaic_kind, only:  r8
	implicit none
	real(r8) :: lamda_wall
! subr argument
	real(r8) :: Dp, x0, y0, a, b

	x0 = 4.7909
	y0 = 9.175e-6
	a  = 9.0134e-6
	b = 0.2531

	lamda_wall = y0 + a*exp(-0.5*((dlog(Dp)-x0)/b)**2)

	return
	end function lamda_wall
