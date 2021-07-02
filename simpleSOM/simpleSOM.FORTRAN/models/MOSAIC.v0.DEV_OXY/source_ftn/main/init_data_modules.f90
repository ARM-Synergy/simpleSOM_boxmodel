      subroutine init_data_modules
!
!   place various constant or initial values into common
      USE mod_MAIN
      USE mod_GAS
      use module_data_mosaic_aero
      use module_data_mosaic_cloud

      implicit none


!--------------------------------------------------
!     define fundamental constants...
       avogad	= 6.02217e+23

!      pi	= 3.141592654
!      deg2rad	= 0.017453293
       pi		= 4.0_r8 * atan( 1.0_r8 )
       piover4	= pi/4.0_r8
       piover6	= pi/6.0_r8
       deg2rad	= pi/180.0_r8
       third	= 1.0_r8/3.0_r8

!---------------------------------
! define species indices
!
! species in inorganic chemistry
       kh2so4		=  1
       khno3		=  2
       khcl		=  3
       knh3		=  4
       kno		=  5
       kno2		=  6
       kno3		=  7
       kn2o5		=  8
       khono		=  9
       khno4		= 10
       ko3		= 11
       ko1d		= 12
       ko3p		= 13
       koh		= 14
       kho2		= 15
       kh2o2		= 16
       kco		= 17
       kso2		= 18
!
! species in methane, ethane, formaldehyde chemistry
       kch4		= 19
       kc2h6		= 20
       kch3o2		= 21
       kethp		= 22
       khcho		= 23
       kch3oh		= 24
       kanol		= 25
       kch3ooh		= 26
       kethooh		= 27
       kald2		= 28
       khcooh		= 29
       krcooh		= 30
       kc2o3		= 31
       kpan		= 32
!
! species in hc1 mechanism. initialize indices to zero
       kpar		= 33
       kaone		= 34
       kmgly		= 35
       keth		= 36
       kolet		= 37
       kolei		= 38
       ktol		= 39
       kxyl		= 40
       kcres		= 41
       kto2		= 42
       kcro		= 43
       kopen	 	= 44
       konit		= 45
       krooh		= 46
       kro2		= 47
       kano2		= 48
       knap		= 49
       kxo2		= 50
       kxpar		= 51
!
! species in hc2 mechanism. initialize indices to zero
       kisop		= 52
       kisoprd		= 53
       kisopp		= 54
       kisopn		= 55
       kisopo2		= 56

! species in organic chemistry; edit wkc
	kCn3     = 57    ! soa prec 1
	kCn2     = 58    ! soa prec 2
	kCn1     = 59    ! soa prec 3
	kC0      = 60    ! soa prec 4
	kC1      = 61    ! soa prec 5
	kC2      = 62    ! soa prec 6
	kC3      = 63    ! soa prec 7
	kC4      = 64    ! soa prec 8
	kC5      = 65    ! soa prec 9
	kC6      = 66    ! soa prec 10
	kC7      = 67    ! soa prec 11
	kC8      = 68    ! soa prec 12
	kC9      = 69    ! soa prec 13

! species in dms mechanism. initialize indices to zero
       kdms		= 70
       kmsa		= 71
       kdmso		= 72
       kdmso2		= 73
       kch3so2h		= 74
       kch3sch2oo   = 75
       kch3so2		= 76
       kch3so3		= 77
       kch3so2oo	= 78
       kch3so2ch2oo	= 79
       ksulfhox		= 80


! aerosol bin: (24 species=bin)
       knum_a		=  1
       kdpdry_a		=  2
       ksigmag_a	=  3
       kjhyst_a		=  4
       kwater_a		=  5
       kso4_a		=  6
       kno3_a		=  7
       kcl_a		=  8
       knh4_a		=  9
       kmsa_a		= 10
       kco3_a		= 11
       kna_a		= 12
       kca_a		= 13
       koin_a		= 14
       koc_a		= 15
       kbc_a		= 16

	  ! organic aerosol species; edit wkc
      kCn3_a      = 17
      kCn2_a      = 18
      kCn1_a      = 19
      kC0_a      = 20
      kC1_a      = 21
      kC2_a      = 22
      kC3_a      = 23
      kC4_a      = 24
      kC5_a      = 25
      kC6_a      = 26
      kC7_a      = 27
      kC8_a      = 28
      kC9_a      = 29


! cloud mode: 1-ncldmode (13=mode)
       knum_c		=  1
       kwater_c		=  2
       kso4_c		=  3
       kno3_c		=  4
       kcl_c		=  5
       knh4_c		=  6
       koc_c		=  7
       kmsa_c		=  8
       kco3_c		=  9
       kna_c		= 10
       kca_c		= 11
       kbc_c		= 12
       koin_c		= 13
!
!
! regime-dependent chemistry definitions
!
       iregime		=  3 ! edit wkc
!
!     GAS
!
       ih2so4		=  1
       ihno3		=  2
       ihcl		=  3
       inh3		=  4
       ino		=  5
       ino2		=  6
       ino3		=  7
       in2o5		=  8
       ihono		=  9
       ihno4		= 10
       io3		= 11
       io1d		= 12
       io3p		= 13
       ioh		= 14
       iho2		= 15
       ih2o2		= 16
       ico		= 17
       iso2		= 18
!
! species in methane, ethane, formaldehyde chemistry
       ich4		= 19
       ic2h6		= 20
       ich3o2		= 21
       iethp		= 22
       ihcho		= 23
       ich3oh		= 24
       ianol		= 25
       ich3ooh		= 26
       iethooh		= 27
       iald2		= 28
       ihcooh		= 29
       ircooh		= 30
       ic2o3		= 31
       ipan		= 32
!
! species in hc1 mechanism. initialize indices to zero
       ipar		= 33
       iaone		= 34
       imgly		= 35
       ieth		= 36
       iolet		= 37
       iolei		= 38
       itol		= 39
       ixyl		= 40
       icres		= 41
       ito2		= 42
       icro		= 43
       iopen		= 44
       ionit		= 45
       irooh		= 46
       iro2		= 47
       iano2		= 48
       inap		= 49
       ixo2		= 50
       ixpar		= 51
!
! species in hc2 mechanism. initialize indices to zero
       iisop		= 52
       iisoprd		= 53
       iisopp		= 54
       iisopn		= 55
       iisopo2		= 56

! gas species in organics mechanism; edit wkc
      iCn3    = 57    ! soa prec 13
      iCn2    = 58    ! soa prec 14
      iCn1    = 59    ! soa prec 15
      iC0     = 60    ! soa prec 16
      iC1     = 61    ! soa prec 17
      iC2     = 62    ! soa prec 18
      iC3     = 63    ! soa prec 19
      iC4     = 64    ! soa prec 20
      iC5     = 65    ! soa prec 21
      iC6     = 66    ! soa prec 22
      iC7     = 67    ! soa prec 23
      iC8     = 68    ! ?
      iC9     = 69    ! ?

! species in dms mechanism. initialize indices to zero
       idms		= 70
       imsa		= 71
       idmso		= 72
       idmso2		= 73
       ich3so2h		= 74
       ich3sch2oo	= 75
       ich3so2		= 76
       ich3so3		= 77
       ich3so2oo	= 78
       ich3so2ch2oo	= 79
       isulfhox		= 80




! alkylperoxy radical indices for parameterized permutation reactions
       jch3o2		=  1
       jethp		=  2
       jro2		=  3
       jc2o3		=  4
       jano2		=  5
       jnap		=  6
       jisopp		=  7
       jisopn		=  8
       jisopo2		=  9
       jxo2		= 10

! photolyzing species indices
       jphoto_no2	=  1
       jphoto_no3	=  2
       jphoto_hono	=  3
       jphoto_hno3	=  4
       jphoto_hno4	=  5
       jphoto_n2o5	=  6
       jphoto_o3a	=  7
       jphoto_o3b	=  8
       jphoto_h2o2	=  9
       jphoto_hchoa	= 10
       jphoto_hchob	= 11
       jphoto_ch3ooh	= 12
       jphoto_ethooh	= 13
       jphoto_ald2	= 14
       jphoto_aone	= 15
       jphoto_mgly	= 16
       jphoto_open	= 17
       jphoto_rooh	= 18
       jphoto_onit	= 19
       jphoto_isoprd	= 20


!
!     CLOUD
!
! cloud (local): used for total and undissociated species
       iso4_c		=  1
       ino3_c		=  2
       icl_c		=  3
       inh4_c		=  4
       ioc_c		=  5
       imsa_c		=  6
       ico2_c		=  7
       ina_c		=  8
       ica_c		=  9
       ibc_c		= 10
       ioin_c		= 11
       iso2_c		= 12
       ihono_c		= 13
       ih2o2_c		= 14
       ich3ooh_c	= 15
       ihcooh_c		= 16
       ircooh_c		= 17
       ihcho_c		= 18
       io3_c		= 19
       iho2_c		= 20	! don't actually exist in cloud (surf rxn only)
       ino2_c		= 21	! don't actually exist in cloud (surf rxn only)
       ino3r_c		= 22	! don't actually exist in cloud (surf rxn only)
       in2o5_c		= 23	! don't actually exist in cloud (surf rxn only)

! cloud ionic species
       jh_c		=  1
       jnh4_c		=  2
       jna_c		=  3
       jhso4_c		=  4
       jso4_c		=  5
       jno3_c		=  6
       jcl_c		=  7
       jno2_c		=  8
       jhso3_c		=  9
       jso3_c		= 10
       jhco3_c		= 11
       jco3_c		= 12
       jho2_c		= 13
       jhcoo_c		= 14
       jrcoo_c		= 15
       jmsa_c		= 16
       joh_c		= 17
       jch2oh2_c	= 18


      RETURN
      END SUBROUTINE
