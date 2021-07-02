      subroutine IntegrateChemistry
      use module_data_mosaic_main
      use module_data_mosaic_gas
      use module_data_mosaic_aero
      use module_data_mosaic_cloud

      implicit none

      real(r8) :: t_in, t_out

      t_in = told_sec
      t_out= tcur_sec

      if(mgas.eq.1)then
        call GasChemistry(t_in, t_out)
      endif

      if(maer.eq.1)then
        call AerChemistry(t_in, t_out)
      endif

      if(mcld.eq.1)then
        call CldChemistry(t_in, t_out)
      endif
      
      
      return
      end subroutine IntegrateChemistry
