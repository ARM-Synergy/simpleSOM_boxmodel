!*************************************************************************
! subroutine PhotoConstants_Solar: calculates photochemical rate constants (1/s)
!
! input: cos_sza (cosine of solar zenith angle from SolarZenithAngle.f)
!        zalt_m (altitude above sea level in meters)
!
!-------------------------------------------------------------------------

      subroutine PhotoConstants_Solar
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      integer jphoto
      real(r8) :: sza_cut, cos_sza_cut, factcld
      parameter (sza_cut = 89.00)		! cutoff solar zenith angle
      parameter (cos_sza_cut = 0.017452406)	! cos of sza_cut

      factcld = 1.0

      if(cos_sza .ge. cos_sza_cut)then	! daytime

         idaytime = 1

         if(mphoto.eq.1)then
           call PhotoParam1
         elseif(mphoto.eq.2)then
           call PhotoParam2
         endif

! apply cloudiness correction factor
         do jphoto = 1, nphoto
           rk_photo(jphoto) = rk_photo(jphoto)*factcld
         enddo

      else				! nighttime

         idaytime = 0

         do jphoto = 1, nphoto
           rk_photo(jphoto) = 0.0
         enddo

      endif

      return
      end subroutine PhotoConstants_Solar
