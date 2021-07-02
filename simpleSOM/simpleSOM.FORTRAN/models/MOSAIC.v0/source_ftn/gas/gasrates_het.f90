


      subroutine GasRates_Het(s)
  
      USE mod_MAIN
      USE mod_GAS

      IMPLICIT NONE

      real(r8) :: s(ntot_max)
      integer igas


! heterogeneous chemistry
      do igas = 1, ngas_max
        r_het(igas) = rk_het(igas)*s(igas)
      enddo

      return
      end subroutine GasRates_Het
