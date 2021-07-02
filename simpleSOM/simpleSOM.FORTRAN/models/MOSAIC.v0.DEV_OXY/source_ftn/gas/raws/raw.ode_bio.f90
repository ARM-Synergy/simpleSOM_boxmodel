!     ==========================================================================
!                 THIS SUBROUTINE DEFINES ODEs FOR BIO-REGIME CHEMISTRY
!     ==========================================================================

      SUBROUTINE ODE_BIO
  
      USE mod_MAIN
      USE mod_GAS

      IMPLICIT NONE
      
      REAL(R8) :: OXY_AVE

!FLAG1
      

      p_BIO(iOH)= 0.0
      d_BIO(iOH)= 0.0

      p_BIO(iISOP) = 0.0
      d_BIO(iISOP) = r_BIO(14)

      RETURN
      END SUBROUTINE ODE_BIO






