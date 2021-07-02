!*******************************************************************
!      SUBROUTINE SOLAR - calculates cosine of zenith angle
!                         for use in photochemical rate coefficient
!                         calculations.
!
!      Nomenclature:
!
!      tmid_sec = time in seconds from Greenich Noon March 21
!
!      cos_za = cosine of zenith angle
!
!      rlon   = local longitude, W hemis. values are negative (radians)
!
!      rlat   = local latitude, S hemis. values are negative (radians)
!
!*******************************************************************

       subroutine SolarZenithAngle
      use module_data_mosaic_main

      implicit none

       real(r8) :: tlocal, tdec, sidec, codec, tloc, thou
!
       tlocal=tmid_sec
       tdec=0.4092797*sin(1.992385E-7*tlocal)
       sidec=sin(tdec)
       codec=cos(tdec)
       tloc=7.272205E-5*tlocal
       thou=cos(rlon+tloc)
       cos_sza=sin(rlat)*sidec+cos(rlat)*codec*thou

       return
       end subroutine SolarZenithAngle
