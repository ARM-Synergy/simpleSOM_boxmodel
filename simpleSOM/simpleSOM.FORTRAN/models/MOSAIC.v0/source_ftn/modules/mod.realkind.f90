!     ==========================================================================
!                           THIS MODULE DEFINES REAL TYPES
!     ==========================================================================

      MODULE mod_REALKIND
      
      IMPLICIT NONE
      
!     REAL OF 8 AND 4 BYTES:      
      INTEGER,PARAMETER :: R8 = SELECTED_REAL_KIND(12)
      INTEGER,PARAMETER :: R4 = SELECTED_REAL_KIND(6)
      
      END MODULE mod_REALKIND
