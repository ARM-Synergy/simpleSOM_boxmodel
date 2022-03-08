



      subroutine GasChemistry(t_in, t_out)
      
      USE mod_MAIN        
      USE mod_GAS

      IMPLICIT NONE

      real(r8) :: t_in, t_out
      integer ntot
      real(r8) :: stot(ntot_max)		! local species array
      real(r8) :: WaterVapor

      REAL(R8) :: OH,mdt
      
      o2       = 0.21*cair_mlc
      h2       = 0.58*1.e-6*cair_mlc
      
      h2o = WaterVapor(RH, cair_mlc, te, pr_atm)


      CNN(kOH) = AX1*EXP(-BX1*time_hrs) + AX2*EXP(-BX2*time_hrs)
      
      OH = CNN(kOH)
      
!!!      PRINT*, 'OH = ',OH

      mdt = t_out - t_in

      CALL STEP_VOCCHEM(VOC_kOH,OH,mdt)
      
      CALL STEP_WALLLOSS(kvap_on,VWL,mdt)
      
      call SelectGasRegime(ntot)	! selects iregime and calculates ntot
      
      call PeroxyRateConstants

      call GasRateConstants

      call SetGasIndices		! set gas indices for selected iregime

      call MapGasSpecies(stot,0)	! map cnn into stot for selected iregime

      call GasIntegrator(ntot,stot,t_in,t_out)
      
      
      j_GASMOLE_SOM(:) = CNN(kCn3:kCn3+nCOMP-1)

      !PRINT*, 'O:C = ',SUM(j_GASMOLE_OXY)/SUM(j_GASMOLE_SOM)/10.
      !READ(*,*)
      
      CALL STEP_OXYCHEM(t_out-t_in)


      call MapGasSpecies(stot,1)	! map stot back into cnn


      return
      end subroutine GasChemistry
