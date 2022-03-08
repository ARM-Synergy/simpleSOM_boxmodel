
!***********************************************************************
! subroutine SetGas_Mar: sets up gas-phase species indices for
! the selected mechanism.
!
! author: Rahul A. Zaveri
! date  : february 1996
!-------------------------------------------------------------------------
      subroutine SetGas_Mar(ilast)

      USE mod_MAIN
      USE mod_GAS

      IMPLICIT NONE

      integer ilast

      idms         = ilast + 1
      imsa         = ilast + 2
      idmso        = ilast + 3
      idmso2       = ilast + 4
      ich3so2h     = ilast + 5
      ich3sch2oo   = ilast + 6
      ich3so2      = ilast + 7
      ich3so3      = ilast + 8
      ich3so2oo    = ilast + 9
      ich3so2ch2oo = ilast + 10
      isulfhox     = ilast + 11

      ilast	   = isulfhox

      return
      end subroutine SetGas_Mar

