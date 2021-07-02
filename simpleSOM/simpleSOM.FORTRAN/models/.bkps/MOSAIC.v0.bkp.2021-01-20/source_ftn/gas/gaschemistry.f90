



      subroutine GasChemistry(t_in, t_out)

  
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      real(r8) :: t_in, t_out
      integer ntot
      real(r8) :: stot(ntot_max)		! local species array
      real(r8) :: WaterVapor

      o2       = 0.21*cair_mlc
      h2       = 0.58*1.e-6*cair_mlc
      
      h2o = WaterVapor(RH, cair_mlc, te, pr_atm)
      
      call SelectGasRegime(ntot)	! selects iregime and calculates ntot
      
      call PeroxyRateConstants

      call GasRateConstants

      call SetGasIndices		! set gas indices for selected iregime

      call MapGasSpecies(stot,0)	! map cnn into stot for selected iregime

      call GasIntegrator(ntot,stot,t_in,t_out)

      call MapGasSpecies(stot,1)	! map stot back into cnn


      return
      end subroutine GasChemistry
