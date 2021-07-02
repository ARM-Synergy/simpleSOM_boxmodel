      subroutine CldChemistry(t_in, t_out)
      use module_data_mosaic_main
      use module_data_mosaic_cloud

      implicit none

      real(r8) :: t_in, t_out

!      call ThermodynamicConstants

!      call SetGasAerCldIndices

!      call GasAerCldIntegrator(t_in,t_out)

      return
      end subroutine CldChemistry
