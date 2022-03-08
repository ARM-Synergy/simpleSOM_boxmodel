

      subroutine CldChemistry(t_in, t_out)
  
      USE mod_MAIN
      use module_data_mosaic_cloud

      IMPLICIT NONE

      real(r8) :: t_in, t_out

!      call ThermodynamicConstants

!      call SetGasAerCldIndices

!      call GasAerCldIntegrator(t_in,t_out)

      return
      end subroutine CldChemistry
