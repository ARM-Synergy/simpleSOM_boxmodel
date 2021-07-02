      module module_data_mosaic_pmcmos

      use module_data_mosaic_kind, only:  r8
      use module_data_mosaic_main, only:  naerbin, ngas_max
      use module_data_mosaic_aero, only:  naer

      implicit none

                                                                                                                                            
      integer, parameter :: ad_maxmode = 6  ! max # of modes per distribution
      integer, parameter :: ad_fname_maxlen= 128  ! max length of an fname
      integer, parameter :: ad_other_maxlen=  64  ! max length of another

      type :: aer_dist_t
         character(len=ad_fname_maxlen) :: dist_fname
         character(len=ad_other_maxlen) :: mode_name(ad_maxmode)
         character(len=ad_other_maxlen) :: mode_type(ad_maxmode)
         character(len=ad_other_maxlen) :: frac_type(ad_maxmode)
         character(len=ad_fname_maxlen) :: frac_fname(ad_maxmode)
         integer :: n_mode
         real(r8) :: dgnum_cm(ad_maxmode) ! number geo-mean diameter (cm)
         real(r8) :: rgnum_m(ad_maxmode)  ! number geo-mean radius (m)
         real(r8) :: sg(ad_maxmode)       ! geo std dev
         real(r8) :: lnsg(ad_maxmode)     ! ln( geo std dev )
         real(r8) :: numden(ad_maxmode)   ! total number (#/m3)
         real(r8) :: massfrac(naer,ad_maxmode) ! mass frac of each species
      end type aer_dist_t


! max number of times for temp, pblh, and gas input arrays
      integer, parameter :: ntmax_pmcmos = 100
! max number of times for aer input arrays
      integer, parameter :: ntmax_aerdist = 10

! Universal gas constant (J mole^{-1} K^{-1}).
      real(r8) :: univ_gas_const_pmcmos = 8.314472d0
! (water vapor molec wght)/(air molec wght)
      real(r8), parameter :: mw_h2o_air = 18.016d0/28.966d0


      integer,  save :: lun_pmcmos_in1   ! unit number for input
      integer,  save :: lun_pmcmos_in2   ! unit number for input


! prescribed background aerosol versus time
      integer,  save :: aer_back_nt                    ! number of times
      real(r8), allocatable, save :: aer_back_tim(:)   ! time (h)
      real(r8), allocatable, save :: aer_back_dil(:)   ! aer dilution rate (1/s)
      type( aer_dist_t ), allocatable, save :: &
                                     aer_back_dist(:)  ! aer distribution (***)
      real(r8), allocatable, save :: aer_back_cnn(:,:) ! 

! prescribed aerosol emissions versus time
      integer,  save :: aer_emit_nt                    ! number of times
      real(r8), allocatable, save :: aer_emit_tim(:)   ! time (h)
      real(r8), allocatable, save :: aer_emit_fac(:)   ! adjustment factor (--_)
      type( aer_dist_t ), allocatable, save :: &
                                     aer_emit_dist(:)  ! aer distribution (***)
      real(r8), allocatable, save :: aer_emit_cnn(:,:) ! 

! initial aerosol
      type( aer_dist_t ) :: aer_init_dist                ! aer distribution (***)


! prescribed background gases versus time
      integer,  save :: gas_back_nt                         ! number of times
      real(r8), allocatable, save :: gas_back_tim(:)        ! time (h)
      real(r8), allocatable, save :: gas_back_dil(:)        ! gas dilution rate (1/s)
      real(r8), allocatable, save :: gas_back_val(:,:)      ! gas mixrat (ppb)

! prescribed gas emissions versus time
      integer,  save :: gas_emit_nt                         ! number of times
      real(r8), allocatable, save :: gas_emit_tim(:)        ! time (h)
      real(r8), allocatable, save :: gas_emit_fac(:)        ! adjustment factor (--)
      real(r8), allocatable, save :: gas_emit_val(:,:)      ! gas mixrat (ppb)

! initial gases
      real(r8), save :: gas_init_val(ngas_max)              ! gas mixrat (ppb)


! prescribed mixing height versus time
      integer,  save :: pblh_nt                  ! number of times
      real(r8), allocatable, save :: pblh_tim(:) ! time (h)
      real(r8), allocatable, save :: pblh_val(:) ! pblh (k)

! prescribed temperature versus time
      integer,  save :: temp_nt                  ! number of times
      real(r8), allocatable, save :: temp_tim(:) ! time (h)
      real(r8), allocatable, save :: temp_val(:) ! temp (k)

! constant in time water vapor mixing ratio
      real(r8), save :: qh2o ! kg/kg
! air molar densities
      real(r8), save :: cair_pmcmos, cair_pmcmos_old ! mol/m3


! names of pmcmos input files
      character(len=128), save :: aer_back_fname
      character(len=128), save :: aer_emit_fname
      character(len=128), save :: aer_init_fname
      character(len=128), save :: gas_back_fname
      character(len=128), save :: gas_emit_fname
      character(len=128), save :: gas_init_fname
      character(len=128), save :: pblh_profile_fname
      character(len=128), save :: temp_profile_fname


! flag for shifting tmid_sec for photol rates to match partmc-mosaic
      integer, save :: msolar_pmcmos_dtshift


!----------------------------------------------------------------------

      end module module_data_mosaic_pmcmos
