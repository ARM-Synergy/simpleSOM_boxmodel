
!**************************************************************************
! subroutine PeroxyRateConstants: calculates parameterized thermal rate
!                     constants for the alkylperoxy radical permutation
!                     reactions for the entire mechanism.
! nomenclature:
! rk_param  = parameterized reaction rate constants (1/s)
! rk_perox  = individual permutation reaction rate constants (molec-cc-s)
! te        = ambient atmospheric temperature (K)
!
! author: Rahul A. Zaveri
! date  : June 1998
!
!-------------------------------------------------------------------------
      subroutine PeroxyRateConstants
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      integer i, j
      real(r8) :: sperox(nperox), rk_perox(nperox,nperox)
      real(r8) :: ARR
!
      sperox(jch3o2)  = cnn(kch3o2)
      sperox(jethp)   = cnn(kethp)
      sperox(jro2)    = cnn(kro2)
      sperox(jc2o3)   = cnn(kc2o3)
      sperox(jano2)   = cnn(kano2)
      sperox(jnap)    = cnn(knap)
      sperox(jisopp)  = cnn(kisopp)
      sperox(jisopn)  = cnn(kisopn)
      sperox(jisopo2) = cnn(kisopo2)
      sperox(jxo2)    = cnn(kxo2)

!
! initialize to zero
      do i = 1, nperox
      rk_param(i) = 0.0
      enddo

      do i = 1, nperox
      do j = 1, nperox
      rk_perox(i,j) = ARR(Aperox(i,j),Bperox(i,j))
      rk_param(i) = rk_param(i) + rk_perox(i,j)*sperox(j)
      enddo
      enddo
!
      return
      end subroutine PeroxyRateConstants
