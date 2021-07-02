!**************************************************************************
! subroutine LoadPeroxyParameters: loads thermal rate coefficients
!                                  for peroxy-peroxy permutation reactions
!
! nomenclature:
! Aperox  = Pre-exponential factor (molec-cc-s)
! Bperox  = activation energy (-E/R)  (K)
!
! author: Rahul A. Zaveri
! date  : June 1998
!
!-------------------------------------------------------------------------
      subroutine LoadPeroxyParameters
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      integer i, j

      Aperox(jch3o2,jch3o2)   = 2.5e-13
      Aperox(jethp,jethp)     = 6.8e-14
      Aperox(jc2o3,jc2o3)     = 2.9e-12
      Aperox(jano2,jano2)     = 8.0e-12
      Aperox(jnap,jnap)       = 1.0e-12
      Aperox(jro2,jro2)       = 5.3e-16
      Aperox(jisopp,jisopp)   = 3.1e-14
      Aperox(jisopn,jisopn)   = 3.1e-14
      Aperox(jisopo2,jisopo2) = 3.1e-14
      Aperox(jxo2,jxo2)       = 3.1e-14

      Bperox(jch3o2,jch3o2)   = 190.
      Bperox(jethp,jethp)     = 0.0
      Bperox(jc2o3,jc2o3)     = 500.
      Bperox(jano2,jano2)     = 0.0
      Bperox(jnap,jnap)       = 0.0
      Bperox(jro2,jro2)       = 1980.
      Bperox(jisopp,jisopp)   = 1000.
      Bperox(jisopn,jisopn)   = 1000.
      Bperox(jisopo2,jisopo2) = 1000.
      Bperox(jxo2,jxo2)       = 1000.

      do i = 1, nperox
      do j = 1, nperox
        if(i.ne.j)then
          Aperox(i,j) = 2.0*sqrt(Aperox(i,i)*Aperox(j,j))
          Bperox(i,j) = 0.5*(Bperox(i,i) + Bperox(j,j))
        endif
      enddo
      enddo
!
! except for
      Aperox(jc2o3,jch3o2) = 1.3e-12
      Aperox(jch3o2,jc2o3) = 1.3e-12
      Bperox(jc2o3,jch3o2) = 640.
      Bperox(jch3o2,jc2o3) = 640.
!
      return
      end subroutine LoadPeroxyParameters

