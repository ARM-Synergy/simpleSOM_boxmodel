!     ==========================================================================
!                    THIS MODULE KEEPS GAS-PHASE INDICES AND ARRAYS
!     ==========================================================================

      MODULE mod_GAS
      
      USE mod_REALKIND, ONLY: R8
      USE mod_MAIN,     ONLY: nGAS_MAX, nGAS_COM, nGAS_URB, nGAS_BIO, nGAS_MAR, nCOMP
      
      IMPLICIT NONE
      
      INTEGER,PARAMETER :: &
                nRXN_HET = nGAS_MAX,       &
                nRXN_COM = 75,             &
                nRXN_URB = 44,             &
                nRXN_BIO = 16 + nCOMP + 1, &
                nRXN_MAR = 35

      INTEGER,PARAMETER :: &
                nREG1 = nGAS_COM,                       &
                nREG2 = nGAS_COM + nGAS_URB,            &
                nREG3 = nGAS_COM + nGAS_URB + nGAS_BIO, &
                nREG4 = nGAS_COM + nGAS_MAR,            &
                nREG5 = nGAS_COM + nGAS_URB + nGAS_MAR, &
                nREG6 = nGAS_COM + nGAS_URB + nGAS_BIO + nGAS_MAR

      real(r8), parameter ::   &
                foh  = 0.228,   &
                fo3  = 0.772,   &
                fno3 = 0.0

      integer, parameter ::   &
      		nperox = 10,	   &  ! total number of alkylperoxy radicals
      		nphoto = 20	! total number of photolyzing species
!------------------------------------------------------------------------
      INTEGER,SAVE :: iREGIME

      REAL(R8),SAVE :: VOC,VOC_MW,VOC_kOH,VOC_CSAT
      REAL(R8),SAVE :: VOC_COEFF(nCOMP)

      REAL(R8),SAVE :: AX1,BX1,AX2,BX2

      REAL(R8),SAVE :: GAS_ON_WALL(nCOMP)

      REAL(R8),SAVE :: kvap_on
      INTEGER,SAVE :: VWL
      
      real(r8), save ::   &
      		rk_com(nrxn_com),   &
      		rk_urb(nrxn_urb),   &
      		rk_bio(nrxn_bio),   &
      		rk_mar(nrxn_mar),   &
      		rk_het(nrxn_het),   &
			rk_dil,		  &
      		rk_param(nperox),   &
      		rk_photo(nphoto),   &
      		Aperox(nperox,nperox),   &
      		Bperox(nperox,nperox)

      real(r8), save ::   &
      		r_com(nrxn_com),   &
      		r_urb(nrxn_urb),   &
      		r_bio(nrxn_bio),   &
      		r_mar(nrxn_mar),   &
      		r_het(nrxn_het)

      REAL(R8),SAVE :: &
               p_com(ngas_max), d_com(ngas_max),   &
      	       p_urb(ngas_max), d_urb(ngas_max),   &
      	       p_bio(ngas_max), d_bio(ngas_max),   &
      	       p_mar(ngas_max), d_mar(ngas_max),   &
      	       p_het(ngas_max), d_het(ngas_max)
!
      real(r8), save ::	Npcasp(15), NOy_in, SO2_in, sum_sfc_area
!
!
!------------------------------------------------------------------------
      integer, save ::   &
       ih2so4,      ihno3,       ihcl,        inh3,        ino,       &
       ino2,        ino3,        in2o5,       ihono,       ihno4,     &
       io3,         io1d,        io3p,        ioh,         iho2,      &
       ih2o2,       ico,         iso2,        ich4,        ic2h6,     &
       ich3o2,      iethp,       ihcho,       ich3oh,      ianol,     &
       ich3ooh,     iethooh,     iald2,       ihcooh,      ircooh,    &
       ic2o3,       ipan,                                             &
       ipar,        iaone,       imgly,       ieth,        iolet,     &
       iolei,       itol,        ixyl,        icres,       ito2,      &
       icro,        iopen,       ionit,       irooh,       iro2,      &
       iano2,       inap,        ixo2,        ixpar,                  &
       iisop,       iisoprd,     iisopp,      iisopn,      iisopo2,   &
       idms,        imsa,        idmso,       idmso2,      ich3so2h,  &
       ich3sch2oo,  ich3so2,     ich3so3,     ich3so2ch2oo,ich3so2oo, &
       isulfhox

      integer, save ::   &
       iCn3,    iCn2,    iCn1,    iC0,        &
       iC1,     iC2,     iC3,     iC4,        &
       iC5,     iC6,     iC7,     iC8,        &
       iC9

      integer, save ::   &
       jch3o2,      jethp,       jro2,        jc2o3,       jano2,   &
       jnap,        jisopp,      jisopn,      jisopo2,     jxo2

      integer, save ::   &
       jphoto_no2,    jphoto_no3,   jphoto_hono,   jphoto_hno3,   &
       jphoto_hno4,   jphoto_n2o5,  jphoto_o3a,    jphoto_o3b,   &
       jphoto_h2o2,   jphoto_hchoa, jphoto_hchob,  jphoto_ch3ooh,   &
       jphoto_ethooh, jphoto_ald2,  jphoto_aone,   jphoto_mgly,   &
       jphoto_open,   jphoto_rooh,  jphoto_onit,   jphoto_isoprd

      real(r8), save ::   &
       mw_gas(ngas_max),   &
       uptake_gas(ngas_max),   &
       D_gas(ngas_max),   &
       vel_gas(ngas_max),   &
       k_gas(ngas_max),   &
       ihet_gas(ngas_max), &
	 delta_api1


      END MODULE mod_GAS
