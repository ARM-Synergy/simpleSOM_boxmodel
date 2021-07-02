! aeroptical.f90
!-----------------------------------------------------------------------
      subroutine aerosol_optical
      use module_data_mosaic_main
      use module_data_mosaic_aero

      implicit none

!   subr arguments
!     none yet

!   local variables
      integer ibin
      real(r8), save :: 	 &
      		vlambc, 	 &  ! wavelength in microns
      		qextc, 		 &  ! extinction efficiency, dimensionless
      		qscatc,		 &  ! scattering efficiency, dimensionless
      		gscac,		 &  ! asymmetry parameter for given size, dimensionless
      		extc,		 &  ! extinction cross section, cm^2
      		scatc		    ! scattering cross section, cm^2



!
!   NOTES
!
!   mshellcore = 0 - do non-shell/core mie calc
!              = 1 - do shell/core mie calc, and shell = { bc }
!             ??? what about insoluble oin and dust species ???
!
!   core volume == 0 if either
!      ri_core_a(ibin) == (0.0, 0.0)
!      dp_core_a(ibin) == 0.0
!   shell volume == 0 if either
!      ri_shell_a(ibin) == (0.0, 0.0)
!

      write( lun_aeroptic, '(a,2f12.3)' )   &
         'in aerosol_optical - run hrs, utc hrs = ', time_hrs, time_utc

      if (it == 0) then
         write( lun_aeroptic, '(a)' )   &
            '   *** currently cannot do optical calcs when it=0'
         write( lun_aeroptic, '(a)' )
         return
      endif

!
! calc optical properties of each bin
!
      write( lun_aeroptic, '(a)' )   &
         '   calling aerosol_optical_properties'
      call aerosol_optical_properties

!
! do mie calculations 
!
      write( lun_aeroptic, '(a)' )   &
           '   doing mie calculations'
      do ibin = 1, nbin_a

! set vlambc to convention "550 nm" wavelength
         vlambc= 0.550  ! wavelength in microns

         call miedriver( dp_wet_a(ibin), dp_core_a(ibin), &
            ri_shell_a(ibin), ri_core_a(ibin), &
            vlambc, qextc, qscatc, gscac, extc, scatc )

! save cross sections & asymmetry parameter for each particle
         ext_cross(ibin) = extc
         scat_cross(ibin) = scatc
         asym_particle(ibin) = gscac

      end do

!
! output
!
      write( lun_aeroptic, '(3a)' )   &
         '   bin, num (#/cm3), dpwet, dpcore (cm)',   &
         ', ri_tot, ri_shell, ri_core',   &
         ', ext_cross, scat_cross, asymmetry_particle'
      do ibin = 1, nbin_a
         write( lun_aeroptic, '(i7,1p,e10.2,4(2x,2e10.2)2x,3e10.2)' )   &
            ibin, num_a(ibin), dp_wet_a(ibin), dp_core_a(ibin),   &
            ri_avg_a(ibin), ri_shell_a(ibin), ri_core_a(ibin),   &
            ext_cross(ibin), scat_cross(ibin), asym_particle(ibin)
      end do
      write( lun_aeroptic, '(a)' )
     
      return
      end subroutine aerosol_optical

