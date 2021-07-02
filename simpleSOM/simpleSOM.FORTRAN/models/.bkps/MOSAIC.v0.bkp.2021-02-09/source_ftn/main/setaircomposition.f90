      subroutine SetAirComposition
      use module_data_mosaic_main
      use module_data_mosaic_gas
      use module_data_mosaic_aero

      implicit none

      integer l, ibin, ib, i1, i2, iaer, kaerstart
      real(r8) :: WaterVapor

!-----------------------------------------------------------------------
      if (ipmcmos <= 0) then
! set bulk air density in (molec/cc)
         cair_mlc = avogad*pr_atm/(82.056*te)	! air conc [molec/cc]
         cair_molm3 = 1.e6*pr_atm/(82.056*te)	! air conc [mol/m^3]
         cair_mlc_old = cair_mlc
         cair_molm3_old = cair_molm3
      end if

      h2o      = WaterVapor(RH, cair_mlc, te, pr_atm)

      o2       = 0.21*cair_mlc
      h2       = 0.58*1.e-6*cair_mlc

!      write(6,60)cair_mlc, h2o
60    format(' air conc. = ', e12.4, '  h2o = ',e12.4)

!-----------------------------------------------------------------------
! conversion factor for converting [ppb] to [molec/cc]
!
      ppb = 1.e+9
!
!-------------------------------------------------------------
! converting gas-phase conc from [ppb] to [molec/cc]
      do l=1, ngas_max
        cnn(l) = cnn(l)*cair_mlc/ppb
      enddo

!
!
! convert from ppb/hr to molec/cc/s
      do l=1, ngas_max
        emission(l) = emission(l)*cair_mlc/ppb/3600.
      enddo
!
      return
      end subroutine SetAirComposition


