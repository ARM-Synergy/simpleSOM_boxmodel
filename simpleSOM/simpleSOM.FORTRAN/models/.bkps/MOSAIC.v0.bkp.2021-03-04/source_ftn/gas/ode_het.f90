!***********************************************************************
      subroutine ode_het
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      integer igas


      do igas = 1, ngas_max
        p_het(igas) = 0.0
        d_het(igas) = r_het(igas)
      enddo

!      p_het(ihno3)  = 2.*r_het(in2o5)
!      p_het(ino)    = r_het(ino3)
!      p_het(ih2so4) = r_het(iso2)

      return
      end subroutine ode_het


