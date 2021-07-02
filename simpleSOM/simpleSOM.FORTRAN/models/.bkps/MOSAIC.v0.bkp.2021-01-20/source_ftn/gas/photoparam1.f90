!*************************************************************************
! subroutine PhotoParam1: calculates photochemical rate constants (1/s)
!
! input: cos_sza (cosine of solar zenith angle from SolarZenithAngle.f)
!        zalt_m (altitude above sea level in meters)
!
!-------------------------------------------------------------------------

      subroutine PhotoParam1
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      real(r8) :: z,   &
           alpha0,  beta0,   &
           alpha1,  beta1,   &
           alpha7,  beta7,   &
           alpha8,  beta8,   &
           alpha49, beta49,   &
           alpha50, beta50

      z = min(1.d-4*zalt_m, 1.1d0)	! zalt in meters

! no2 + hv --> no + o3p
      alpha0=2.10223e-2-1.6e-3*(z-1.15)**2.
      alpha1=-12.2e-2+alpha0**.5
      beta0=1.258-.16*(z-2.15)**2.
      beta1=-1.2+beta0**.5
      rk_photo(jphoto_no2)    = alpha1*exp(beta1/cos_sza)

! no3 + hv --> .89 no2 + .89 o3p + .11 no
      rk_photo(jphoto_no3)    = 23.8*rk_photo(jphoto_no2)

! hono + hv --> oh + no
      rk_photo(jphoto_hono)   = 0.197*rk_photo(jphoto_no2)

! hno3 + hv --> oh + no2
      rk_photo(jphoto_hno3)   = 3.3e-5*rk_photo(jphoto_no2)

! hno4 + hv --> ho2 + no2
      rk_photo(jphoto_hno4)   = 5.9e-4*rk_photo(jphoto_no2)

! n2o5 + hv --> no2 + no3
      rk_photo(jphoto_n2o5)   = 0.0*rk_photo(jphoto_no2)

! o3 + hv --> o3p
      rk_photo(jphoto_o3a)    = 0.053*rk_photo(jphoto_no2)

! o3 + hv --> o1d
      alpha0=2.69924e-7-4.0e-8*(z-.375)**2.
      alpha7=-3.25e-4+sqrt(alpha0)
      beta0=4.173-.64*(z-2.0)**2.
      beta7=-3.2+sqrt(beta0)
      rk_photo(jphoto_o3b)    = alpha7*exp(beta7/cos_sza)

! h2o2 + hv --> 2 oh
      alpha0=2.540534e-9-4.e-11*(z-0.75)**2.
      alpha8=-3.e-5+sqrt(alpha0)
      beta0=.539284-.16*(z-1.75)**2.
      beta8=-1+sqrt(beta0)
      rk_photo(jphoto_h2o2)   = alpha8*exp(beta8/cos_sza)

! hcho + hv ---> 2 ho2 + co
      alpha0=7.1747e-9-1.6e-9*(z-.75)**2.
      alpha49=-4.e-5+sqrt(alpha0)
      beta0=.7631144-.16*(z-2.)**2.
      beta49=-1.2+sqrt(beta0)
      rk_photo(jphoto_hchoa)  = alpha49*exp(beta49/cos_sza)

! hcho + hv ---> co
      alpha0=1.693813e-7-4.e-8*(z-.875)**2.
      alpha50=-2.5e-4+sqrt(alpha0)
      beta0=.7631144-.16*(z-1.875)**2.
      beta50=-1.1+sqrt(beta0)
      rk_photo(jphoto_hchob)  = alpha50*exp(beta50/cos_sza)

      rk_photo(jphoto_ch3ooh) = 0.7   *rk_photo(jphoto_h2o2)
      rk_photo(jphoto_ethooh) = 0.7   *rk_photo(jphoto_h2o2)
      rk_photo(jphoto_ald2)   = 4.6e-4*rk_photo(jphoto_no2)
      rk_photo(jphoto_aone)   = 7.8e-5*rk_photo(jphoto_no2)
      rk_photo(jphoto_mgly)   = 9.64  *rk_photo(jphoto_hchoa)
      rk_photo(jphoto_open)   = 9.04  *rk_photo(jphoto_hchoa)
      rk_photo(jphoto_rooh)   = 0.7   *rk_photo(jphoto_h2o2)
      rk_photo(jphoto_onit)   = 1.0e-4*rk_photo(jphoto_no2)
      rk_photo(jphoto_isoprd) = .025  *rk_photo(jphoto_hchob)

      return
      end subroutine PhotoParam1
