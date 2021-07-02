      module module_data_mosaic_asect

      use module_data_mosaic_kind, only:  r8
      use module_data_mosaic_aero, only:  nbin_a_max

      implicit none

!-----------------------------------------------------------------------
!
!   The variables in this module provide a means of organizing and accessing
!   aerosol species by their chemical component, size bin (or mode), "type", and "phase"
!
!   Their purpose is to allow flexible coding of process modules,
!   compared to "hard-coding" using specify indices.
!   Most (if not all) of these variables are usxwed in the i
!   WRF-chem MOSAIC implementation.
!
!-----------------------------------------------------------------------
!
!   maxd_atype = maximum allowable number of aerosol types
!   maxd_asize = maximum allowable number of aerosol size bins
!   maxd_acomp = maximum allowable number of chemical components
!	in each aerosol size bin
!   maxd_aphase = maximum allowable number of aerosol phases
!	(gas, cloud, ice, rain, ...)
!
!   ntype_aer = number of aerosol types
!	The aerosol type will allow treatment of an externally mixed
!	aerosol.  For a traditional internally-mixed sectional approach,
!	ntype_aer=1.  Eventually, multiple types
!	could treat fresh primary BC/OC, fresh SO4 from nucleation,
!	aged BC/OC/SO4/... mixture, soil dust, sea salt, ...
!
!   nphase_aer = number of aerosol phases
!
!   ai_phase = phase (p) index for interstitial (unactivated) aerosol particles
!   cw_phase = phase (p) index for aerosol particles in cloud water
!   ci_phase = phase (p) index for aerosol particles in cloud ice
!   rn_phase = phase (p) index for aerosol particles in rain
!   sn_phase = phase (p) index for aerosol particles in snow
!   gr_phase = phase (p) index for aerosol particles in graupel
!   [Note:  the value of "xx_phase" will be between 1 and nphase_aer
!	for phases that are active in a simulation.  The others
!	will have non-positive values.]
!
!   nsize_aer(t) = number of aerosol size bins for aerosol type t
!
!   ncomp_aer(t) = number of "regular" chemical components for aerosol type t
!   ncomp_plustracer_aer(t) = number of "regular" plus "tracer"
!	chemical components for aerosol type t
!   [Note:  only "regular" components are used for calculating
!	aerosol physical (mass, volume) and chemical properties.
!	"Tracer" components are optional, and can be used to track source
!	regions, source mechanisms, etc.]
!   [Note:  for aerosol type t, all phases have the same number of size
!	bins, and all size bins have the same number of
!	both regular and tracer components.]
!
!   ntot_mastercomp_aer = number of aerosol chemical components defined
!	in the "master component list".
!   [Note:  each aerosol type will use some but not necessarily all
!	of the components in the "master component list".]
!
!   mastercompptr_aer(c,t) = the position/index/i.d. in the
!       "master component list" for chemical component c of aerosol type t.
!	(1=sulfate, others to be defined by user.)
!
!   massptr_aer(c,s,t,p) = the position/index in the chem array for mixing-
!	ratio for chemical component c, size bin s, type t, and phase p.
!
!   lptr_so4_aer(s,t,p) = the position/index in the chem array for mixing-
!	ratio for sulfate for aerosol size bin s, type t, and phase p
!   (similar lptr's are defined for no3, cl, msa, co3,
!	nh4, na, ca, oin, oc, bc, ...)
!   [Note:  the massptr_aer allow you to loop over all species of
!	an aerosol type.  The lptr_so4_aer, etc., allow you to access
!	a specific chemical component.]
!
!   waterptr_aer(s,t) = the position/index in the chem array for mixing-
!	ratio of aerosol water content for size bin s, type t.
!	[Note:  water content is only carried for the interstitial aerosol
!	phase, so there is no p dimension.]
!
!   hyswptr_aer(s,t) = the position/index in the chem array for mixing-
!	ratio of aerosol "hysteresis water" content for size bin s, type t.
!	This is used to determine if aerosol is in the dry or wet state, when
!	the ambient RH is between the crystallization and deliquescence RH.
!	[Note:  hysteresis water content is only carried for the
!	interstitial aerosol phase, so there is no p dimension.]
!
!   numptr_aer(s,t,p) = the position/index in the chem array for mixing-
!	ratio of particle number for size bin s, type t, and phase p.
!
!   mprognum_aer(s,t,p) - if positive, number mixing-ratio for size s, type t,
!       and phase p will be prognosed.  Otherwise, it is diagnosed using
!	mass mixing-ratio add assumed/prescribed size.
!
!	mixing ratio (mol-water/mol-air) for water
!       associated with aerosol size bin s and type t
!
!   ibin_of_isize_itype(s,t) - maps from the sectional isize,itype
!	to the mosaic "1-D bin index"
!   isize_of_ibin(b) - maps from the mosaic "1-D bin index" ibin
!	to the sectional isize
!   itype_of_ibin(b) - maps from the mosaic "1-D bin index" ibin
!	to the sectional itype
!
!   itype_of_itype_md1md2(t1,t2) - maps from the "new 3d sectional"
!	itype_md1,itype_md2 to the "old sectional" itype
!   itype_md1_of_itype(t) - maps from the "old sectional" itype
!	to the "new 3d sectional" itype_md1
!   itype_md2_of_itype(t) - maps from the "old sectional" itype
!	to the "new 3d sectional" itype_md2
!
!   mastercompindx_so4_aer = the position/index in the
!       "master component list" for sulfate.
!   (similar lptr's are defined for no3, cl, msa, co3,
!	nh4, na, ca, oin, oc, bc, ...)
!   [Note:  the mastercompindx_xxx_aer are used primarily in
!	initialization routines, and generally aren't needed elsewhere.]
!
!-----------------------------------------------------------------------
!
!   dens_mastercomp_aer(mc) = dry density (g/cm^3) of component mc
!	of the master component list.
!   dens_aer(c,t) = dry density (g/cm^3) of aerosol chemical component
!	c of type t
!   [Note:  dens_aer(c,t) == dens_mastercomp_aer(mastercompptr_aer(c,t))
!	The dens_mastercomp_aer is used in some initialization routines.
!	The dens_aer is used in most other places because of convenience.]
!
!   mw_mastercomp_aer(mc) = molecular weight (g/mol) of component mc
!	of the master component list.
!   mw_aer(c,t) = molecular weight (g/mol) of aerosol chemical component
!	c of type t
!   [Note:  mw_aer(c,t) == mw_mastercomp_aer(mastercompptr_aer(c,t)) ]
!
!   name_mastercomp_aer(mc) = name of component mc of the
!	master component list (e.g., "sulfate", "nitrate", ...).
!   name_aer(c,t) = name of aerosol chemical component c of type t
!   [Note:  name_aer(c,t) == name_mastercomp_aer(mastercompptr_aer(c,t)) ]
!
!   hygro_mastercomp_aer(mc) = bulk hygroscopicity (--) at dilute conditions
!	(RH near 100%) of component mc of the master component list.
!   hygro_aer(c,t) = bulk hygroscopicity (--) at dilute conditions
!	(RH near 100%) of aerosol chemical component c of type t
!   [For definition of bulk hygroscopicity,
!	see Abdul-Razzak and Ghan, 2004, J Geophys Res, V105, p. 6837-6844.]
!   [*** this bulk hygroscopicity is equivalent to the "kappa" of
!       Peters and Kreidenweis, 2007, Atmos Chem Phys, V7, p. 1961-1971.]
!   [Note:  hygro_aer(c,t) == hygro_mastercomp_aer(mastercompptr_aer(c,t)) ]
!
!-----------------------------------------------------------------------
!
!   volumlo_sect(s,t) = 1-particle volume (cm^3) at lower boundary of section m
!   volumhi_sect(s,t) = 1-particle volume (cm^3) at upper boundary of section m
!   volumcen_sect(s,t)= 1-particle volume (cm^3) at "center" of section m
!
!   dlo_sect(s,t) = 1-particle diameter (cm) at lower boundary of section m
!   dhi_sect(s,t) = 1-particle diameter (cm) at upper boundary of section m
!   dcen_sect(s,t) = 1-particle diameter (cm) at "center" section m
!
!   [Note:  the "center" values are defined as follows:
!       volumcen_sect == 0.5*(volumlo_sect + volumhi_sect)
!                     == (pi/6) * (dcen_sect**3) ]
!
!-----------------------------------------------------------------------

!
! the sectional mosaci uses a 2d bin structure
!    dimension 1 (dry diameter) = "size"
!    dimension 2 (composition ) = "type"
!
! for the newer 3d bin structure,
!    dimension 1 (dry diameter)     = "size" still
!    dimension 2 (bc mass fraction) = "type_md1"
!    dimension 3 (hygroscopicity  ) = "type_md2"
!
!    (dimension 2) x (dimension 3) is mapped to the old 1d "type"
!
	integer, parameter :: maxd_atype_md1 = 29
	integer, parameter :: maxd_atype_md2 = 30
	integer, parameter :: maxd_atype = 1
!	integer, parameter :: maxd_atype = maxd_atype_md1*maxd_atype_md2

	integer, parameter :: maxd_asize = 1000 ! 110

	integer, parameter :: maxd_acomp = 25
	integer, parameter :: maxd_aphase = 1

	integer, parameter :: lunerr = 6
	integer, save :: lunout = 170

	integer, save :: ai_phase = -999888777
	integer, save :: cw_phase = -999888777
	integer, save :: ci_phase = -999888777
	integer, save :: rn_phase = -999888777
	integer, save :: sn_phase = -999888777
	integer, save :: gr_phase = -999888777

	integer, save :: ntype_aer = 0 ! number of types
	integer, save :: ntype_md1_aer = 0 ! number of md1 types
	integer, save :: ntype_md2_aer = 0 ! number of md2 types

	integer, save :: ntot_mastercomp_aer = 0 ! number of master components
	integer, save :: nphase_aer = 0 ! number of phases

	integer, save ::   &
	  nsize_aer( maxd_atype ),   & ! number of size bins
	  ncomp_aer( maxd_atype ),   & ! number of chemical components
	  ncomp_plustracer_aer( maxd_atype ),   &
	  mastercompptr_aer(maxd_acomp, maxd_atype), &   !  mastercomp index
	  massptr_aer( maxd_acomp, maxd_asize, maxd_atype, maxd_aphase ), &
		! index for mixing ratio
	  waterptr_aer( maxd_asize, maxd_atype ), & ! index for aerosol water
	  hyswptr_aer( maxd_asize, maxd_atype ), &
	  numptr_aer( maxd_asize, maxd_atype, maxd_aphase ), &
		! index for the number mixing ratio
	  mprognum_aer(maxd_asize,maxd_atype,maxd_aphase)

	integer, save ::   &
	  ibin_of_isize_itype( maxd_asize, maxd_atype ), &
	  isize_of_ibin( nbin_a_max ), &
	  itype_of_ibin( nbin_a_max ), &
	  itype_of_itype_md1md2(maxd_atype_md1,maxd_atype_md2), &
	  itype_md1_of_itype(maxd_atype), &
	  itype_md2_of_itype(maxd_atype)


!   these indices give the location in the "mastercomp list" of
!   the different aerosol chemical (or tracer) components
	integer, save :: mastercompindx_so4_aer = -999888777
	integer, save :: mastercompindx_no3_aer = -999888777
	integer, save :: mastercompindx_cl_aer  = -999888777
	integer, save :: mastercompindx_msa_aer = -999888777
	integer, save :: mastercompindx_co3_aer = -999888777
	integer, save :: mastercompindx_nh4_aer = -999888777
	integer, save :: mastercompindx_na_aer  = -999888777
	integer, save :: mastercompindx_ca_aer  = -999888777
	integer, save :: mastercompindx_oin_aer = -999888777
	integer, save :: mastercompindx_oc_aer  = -999888777
	integer, save :: mastercompindx_bc_aer  = -999888777


	real(r8), save ::   &
	  dens_aer( maxd_acomp, maxd_atype ),  &
	  dens_mastercomp_aer( maxd_acomp ),   &
      	  mw_mastercomp_aer( maxd_acomp ),     &
      	  mw_aer( maxd_acomp, maxd_atype ),    &
      	  hygro_mastercomp_aer( maxd_acomp ),  &
      	  hygro_aer( maxd_acomp, maxd_atype )

	real(r8), save ::   &
	  volumcut_sect( 0:maxd_asize, maxd_atype ),  &
	  volumcen_sect(   maxd_asize, maxd_atype ),  &
	  volumlo_sect(    maxd_asize, maxd_atype ),   &
	  volumhi_sect(    maxd_asize, maxd_atype ),   &
	  dcut_sect( 0:maxd_asize, maxd_atype ),      &
	  dcen_sect(   maxd_asize, maxd_atype ),      &
	  dlo_sect(    maxd_asize, maxd_atype ),       &
	  dhi_sect(    maxd_asize, maxd_atype ),       &
	  sigmag_aer(maxd_asize, maxd_atype)

! these are the cut values that separate the "md1" and "md2" types
! for example, the bc mass fraction (itype_md1) is between
!    xcut_atype_md1(0) and xcut_atype_md1(1) for itype_md1=1
!    xcut_atype_md1(2) and xcut_atype_md1(2) for itype_md1=2
	real(r8), save ::    &
	  xcut_atype_md1( 0:maxd_atype_md1 ),   &
	  xcut_atype_md2( 0:maxd_atype_md2 )

	character*10, save ::   &
      	  name_mastercomp_aer( maxd_acomp ),  &
      	  name_aer( maxd_acomp, maxd_atype )

	integer, save ::                     &
      	  lptr_so4_aer(maxd_asize, maxd_atype, maxd_aphase),      &
      	  lptr_msa_aer(maxd_asize, maxd_atype, maxd_aphase),      &
      	  lptr_no3_aer(maxd_asize, maxd_atype, maxd_aphase),      &
      	  lptr_cl_aer(maxd_asize, maxd_atype, maxd_aphase),       &
	  lptr_co3_aer(maxd_asize, maxd_atype, maxd_aphase),      &
      	  lptr_nh4_aer(maxd_asize, maxd_atype, maxd_aphase),      &
      	  lptr_na_aer(maxd_asize, maxd_atype, maxd_aphase),       &
      	  lptr_ca_aer(maxd_asize, maxd_atype, maxd_aphase),       &
      	  lptr_oin_aer(maxd_asize, maxd_atype, maxd_aphase),      &
      	  lptr_oc_aer(maxd_asize, maxd_atype, maxd_aphase),       &
      	  lptr_bc_aer(maxd_asize, maxd_atype, maxd_aphase)

!   in the mosaic box model, the molecular weight, densities,
!       and hygroscopities in module_data_mosaic_aero are the
!       correct ones to use
!   those values are copied into the mw_aer and dens_aer arrays
!   *** the individual "parameter" values below should not
!       be used and thus are left undefined
!
!   molecular weights (g/mol)
!	real(r8), parameter :: mw_so4_aer = 96.066
!	real(r8), parameter :: mw_no3_aer = 62.007
!	real(r8), parameter :: mw_cl_aer  = 35.450
!	real(r8), parameter :: mw_msa_aer = 96.109
!	real(r8), parameter :: mw_co3_aer = 60.007
!	real(r8), parameter :: mw_nh4_aer = 18.042
!	real(r8), parameter :: mw_na_aer  = 22.990
!	real(r8), parameter :: mw_ca_aer  = 40.080
!	real(r8), parameter :: mw_oin_aer = 1.0
!	real(r8), parameter :: mw_oc_aer  = 1.0
!	real(r8), parameter :: mw_bc_aer  = 1.0

!   dry densities (g/cm3)
!	real(r8), parameter :: dens_so4_aer = 1.80
!	real(r8), parameter :: dens_no3_aer = 1.80
!	real(r8), parameter :: dens_cl_aer  = 2.20
!	real(r8), parameter :: dens_msa_aer = 1.80
!	real(r8), parameter :: dens_co3_aer = 2.60
!	real(r8), parameter :: dens_nh4_aer = 1.80
!	real(r8), parameter :: dens_na_aer  = 2.20
!	real(r8), parameter :: dens_ca_aer  = 2.60
!	real(r8), parameter :: dens_oin_aer = 2.60
!	real(r8), parameter :: dens_oc_aer  = 1.00
!	real(r8), parameter :: dens_bc_aer  = 1.70

!   water molecular weights (g/mol) and density (g/cm3)
	real(r8), parameter :: mw_water_aer  = 18.0
!   29-mar-2010 - change to allow testing of non-cgs densities
!	real(r8), parameter :: dens_water_aer  = 1.0
	real(r8), save      :: dens_water_aer  = 1.0

!   hygroscopicities (dimensionless)
!	real(r8), parameter :: hygro_so4_aer = 0.5
!	real(r8), parameter :: hygro_no3_aer = 0.5
!	real(r8), parameter :: hygro_ca_aer  = 0.1
!	real(r8), parameter :: hygro_co3_aer = 0.1
!	real(r8), parameter :: hygro_nh4_aer = 0.5
!	real(r8), parameter :: hygro_msa_aer = 0.58
!	real(r8), parameter :: hygro_cl_aer  = 1.16
!	real(r8), parameter :: hygro_na_aer  = 1.16
!	real(r8), parameter :: hygro_oin_aer = 0.14
!	real(r8), parameter :: hygro_oc_aer  = 0.14
!	real(r8), parameter :: hygro_bc_aer  = 1.e-6


!-----------------------------------------------------------------------
!   following are used in movesect, newnuc, and coag routines
!   to identify bins with essentially negligible mass
!
!   if bin mass mixrat < smallmassaa (1.0e-22 g/g-air),
!   then assume no growth AND no water AND conform number so that size is within bin limits
	real(r8), parameter :: smallmassaa = 1.0e-22_r8
!   if bin mass mixrat < smallmassab (1.0e-32 g/g-air),
!   then assume default density to avoid divide by zero
	real(r8), parameter :: smallmassbb = 1.0e-32_r8
!
!   with single-particle diameter = 1 nm and mass ~1e-21 g,
!	and number = 1e-4 #/cm3 ~= 1e-1 #/g-air, the mass mixing ratio ~= 1e-22 g/g-air
!   for simulations focusing on nucleation and ultrafine particles,
!      one might want to use reduce smallmassaa


!-----------------------------------------------------------------------
!
!   following are used by coag, movesect, and newnuc routines
!   which were adapted from wrf-chem version of mosaic
!
!-----------------------------------------------------------------------
!	integer, parameter :: lunerr = 6
!	integer, save :: lunout = 6


!-----------------------------------------------------------------------
!
!   following are used by coag, movesect, and newnuc routines
!   which were adapted from wrf-chem version of mosaic,
!   BUT are subr parameters in the mosaic box code
!
!	drymass_pregrow(s,t) = dry-mass (g/mol-air) before gas-aerosol mass transfer
!	drymass_aftgrow(s,t) = dry-mass (g/mol-air) after   "     "     "      "
!	drydens_pregrow(s,t) = dry-density (g/cm3)  before  "     "     "      "
!	drydens_aftgrow(s,t) = dry-density (g/cm3)  after   "     "     "      "
!
!       aqvoldry_box(s,t)  = dry-volume mixing ratio (cm^3-aerosol/mol-air)
!       aqmassdry_box(s,t) = dry-mass mixing ratio (g-aerosol/mol-air)
!
!       adrydens_box(s,t)  = dry-density (g-aerosol/cm^3-aerosol)
!                            == amassdry_box/avoldry_box
!       awetdens_box(s,t)  = wet-density (g-aerosol/cm^3-aerosol)
!
!       admeandry_box(s,t) = current mean dry-diameter (cm) for unactivated aerosol
!       admeanwet_box(s,t) = current mean wet-diameter (cm) for unactivated aerosol
!
!-----------------------------------------------------------------------



	end module module_data_mosaic_asect
