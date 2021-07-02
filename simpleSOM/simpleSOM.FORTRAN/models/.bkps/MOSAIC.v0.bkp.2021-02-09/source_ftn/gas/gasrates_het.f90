      subroutine GasRates_Het(s)
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      real(r8) :: s(ntot_max)
      integer igas


! heterogeneous chemistry
      do igas = 1, ngas_max
        r_het(igas) = rk_het(igas)*s(igas)
      enddo

      return
      end subroutine GasRates_Het
