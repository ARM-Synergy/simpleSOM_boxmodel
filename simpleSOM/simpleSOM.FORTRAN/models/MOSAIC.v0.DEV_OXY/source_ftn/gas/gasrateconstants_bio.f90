!**************************************************************************
! subroutine GasRateConstants_Bio: generates thermal rate coefficients
!                   for the selected mechanism
! nomenclature:
! rk_bio    = reaction rate constants for hc2 mechanism    (molec-cc-s)
! te        = ambient atmospheric temperature (K)
!
! author: Rahul A. Zaveri
! date  : february 1996
!
!-------------------------------------------------------------------------

      SUBROUTINE GASRATECONSTANTS_BIO
  
      USE mod_MAIN
      USE mod_GAS

      IMPLICIT NONE

      REAL(R8) :: ARR

!     ISOPRENE CHEMISTRY:
      rk_bio(1)  = ARR(2.6d-11, 409.d0)
      rk_bio(2)  = ARR(1.2d-14, -2013.d0)
      rk_bio(3)  = ARR(3.0d-12, -446.d0)
      rk_bio(4)  = rk_photo(jphoto_isoprd)
      rk_bio(5)  = 3.3e-11
      rk_bio(6)  = 7.0e-18
      rk_bio(7)  = 1.0e-15
      rk_bio(8)  = 4.0e-12
      rk_bio(9)  = 4.0e-12
      rk_bio(10) = 4.0e-12
      rk_bio(11) = ARR(1.7d-13, 1300.d0)
      rk_bio(12) = ARR(1.7d-13, 1300.d0)
      rk_bio(13) = ARR(1.7d-13, 1300.d0)
      rk_bio(14) = rk_param(jisopp)
      rk_bio(15) = rk_param(jisopn)
      rk_bio(16) = rk_param(jisopo2)
      
!FLAG1
      rk_BIO(17) = 5.53e-11
      rk_BIO(18) = 5.51e-11
      rk_BIO(19) = 5.43e-11
      rk_BIO(20) = 5.29e-11
      rk_BIO(21) = 5.08e-11
      rk_BIO(22) = 4.82e-11
      rk_BIO(23) = 4.49e-11
      rk_BIO(24) = 4.10e-11
      rk_BIO(25) = 3.65e-11
      rk_BIO(26) = 3.14e-11
      rk_BIO(27) = 2.57e-11
      rk_BIO(28) = 1.94e-11
      rk_BIO(29) = 1.24e-11
      
      
      rk_BIO(30) = 1.00e-10

      RETURN
      END SUBROUTINE GASRATECONSTANTS_BIO



