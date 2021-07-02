! external dummy jacobian evaluation for LSODES (when mf=222)
      subroutine jac(ngas,tt,s,j,ian,jan,pdj)

      use module_data_mosaic_kind, only:  r8
      implicit none

      integer ngas, j, ian(1), jan(1)
      real(r8) :: tt, s(1), pdj(1)

      return
      end subroutine jac
