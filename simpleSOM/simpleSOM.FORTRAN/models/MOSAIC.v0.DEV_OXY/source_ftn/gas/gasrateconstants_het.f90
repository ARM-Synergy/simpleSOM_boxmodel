!**************************************************************************
! subroutine GasRateConstants_Het: generates thermal rate coefficients
!                   for the selected mechanism
! nomenclature:
! rk_het    = reaction rate constants for heterogeneous rxns (1/s)
!
! author: Rahul A. Zaveri
! date  : June 2006
!
!-------------------------------------------------------------------------
      subroutine GasRateConstants_Het
  
      USE mod_MAIN
      USE mod_GAS

      IMPLICIT NONE

      logical first
      save first
      data first/.true./
      integer k, igas
      real(r8) :: Dp(15), sfc_area, sum_k_A	! Dp in cm

      data Dp /.12e-4, .14e-4, .17e-4, .20e-4, .25e-4, .31e-4,   &
      	       .40e-4, .54e-4, .76e-4, .95e-4,   &
               1.20e-4, 1.32e-4, 1.75e-4, 2.25e-4, 2.75e-4/


      if(first)then
        first=.false.

        do igas = 1, ngas_max
          mw_gas(igas)     = 1.0	! molecular weight
          uptake_gas(igas) = 1.0	! reaction probablity or uptake coefficient
          ihet_gas(igas)   = 0		! flag to turn on/off reaction
          D_gas(igas)      = 0.1	! gas-phase diffusivity in air [cm^2/s]
        enddo

        mw_gas(io3)   = 48.0
        mw_gas(in2o5) = 108.0
        mw_gas(ihno3) = 63.0
        mw_gas(ino3)  = 62.0
        mw_gas(iho2)  = 33.0
        mw_gas(ino2)  = 46.0
        mw_gas(ino)   = 30.0
        mw_gas(iro2)  = 75.0	! assumed as C3H7O2
        mw_gas(iso2)  = 64.0

        uptake_gas(io3)   = 1.e-3	! O3 -->
        uptake_gas(in2o5) = 0.1	! N2O5 --> 2HNO3
        uptake_gas(ihno3) = 0.1	! HNO3 --> NO2
        uptake_gas(ino3)  = 0.1	! NO3 --> NO + O2
        uptake_gas(iho2)  = 0.1	! HO2 --> 0.5H2O2
        uptake_gas(ino2)  = 0.1	! NO2 --> 0.5HONO + 0.5HNO3
        uptake_gas(ino)   = 0.1	! NO --> ?
        uptake_gas(iro2)  = 0.1	! RO2 -->
        uptake_gas(iso2)  = 0.1	! SO2 --> H2SO4

        ihet_gas(io3)   = 0	! O3 -->
        ihet_gas(in2o5) = 0	! N2O5 --> 2HNO3
        ihet_gas(ihno3) = 0	! HNO3 --> NO2
        ihet_gas(ino3)  = 0	! NO3 --> NO + O2
        ihet_gas(iho2)  = 0	! HO2 --> 0.5H2O2
        ihet_gas(ino2)  = 0	! NO2 --> 0.5HONO + 0.5HNO3
        ihet_gas(ino)   = 0	! NO --> ?
        ihet_gas(iro2)  = 0	! RO2 -->
        ihet_gas(iso2)  = 0	! SO2 --> H2SO4

        ! edit wkc: parameters for heterogeneous chemistry
        mw_gas(iOH) = 17.0 ! g/mol
        uptake_gas(iOH) = 0.1 ! gamma, uptake coefficient
        ihet_gas(iOH)   = 0 ! switch for OH causing hetchem

      endif


! update temperature dependent molecular speed of gases
        do igas = 1, ngas_max
          vel_gas(igas) = 1.455e4 * sqrt(te/mw_gas(igas))	! avg. molec speed [cm/s]
        enddo

! compute integral first-order heterogeneous reaction rate constants
      do igas = 1, ngas_max

        sum_k_A = 0.0
        do k = 1, 15
          sfc_area   = Npcasp(k)*0.25*3.14159*Dp(k)**2	! cm^2/cm^3(air)
          k_gas(igas) =           ihet_gas(igas)/	   &  ! mass transfer coefficient [cm/s]
          (0.5*Dp(k)/D_gas(igas)+4./(vel_gas(igas)*uptake_gas(igas)))
          sum_k_A = sum_k_A + k_gas(igas)*sfc_area	! integral first order mtc [1/s]
        enddo

        rk_het(igas) = sum_k_A		! [1/s]

      enddo

      return
      end subroutine GasRateConstants_Het


