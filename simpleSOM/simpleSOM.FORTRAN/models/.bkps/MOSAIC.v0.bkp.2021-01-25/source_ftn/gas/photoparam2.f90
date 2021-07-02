!*************************************************************************
! subroutine PhotoParam2: calculates photochemical rate constants (1/s)
!
! input: cos_sza (cosine of solar zenith angle from SolarZenithAngle.f)
!        zalt_m (altitude above sea level in meters)
!
!-------------------------------------------------------------------------

      subroutine PhotoParam2
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      integer idum, jdum, alpha, k, km1, a, a1
      real(r8) :: kr(3,nphoto)	! kr(level, species #):
      real(r8) :: cz

      cz = cos_sza

      kr = 0d0

! surface photolysis rates:
! NO2
	KR(1,jphoto_no2)=   &
        -1.0184E-3+cz*1.8542E-2-cz**2*9.5368E-3+cz**3*1.8165E-3
! NO3
	KR(1,jphoto_no3)=   &
        4.3945E-3+0.556446*cz-0.71996*cz**2+0.34253*cz**3
! O3
        KR(1,jphoto_o3b)=(-3.2168E-8+cz*2.4588E-6-cz**2*2.6274E-5+   &
       	1.2005E-4*cz**3-7.6216E-5*cz**4+1.4628E-5*cz**5)
! HONO
	KR(1,jphoto_hono)=   &
        -1.7863E-4+3.2272E-3*cz-8.5989E-4*cz**2-1.8987E-4*cz**3
! HNO3
	KR(1,jphoto_hno3)=   &
        1.9592E-8-2.8147E-7*cz+1.3533E-6*cz**2-4.2010E-7*cz**3
! H2O2
	KR(1,jphoto_h2o2)=   &
        -2.1672E-7+4.0070E-6*cz+1.1189E-5*cz**2-6.4306E-6*cz**3
! HNO4
	KR(1,jphoto_hno4)=   &
        2.1392E-8-2.0854E-7*cz+1.6131E-5*cz**2-7.2297E-6*cz**3
! HCHOB
	KR(1,jphoto_hchoa)=   &
        -5.4160E-8+9.4694E-7*cz+6.4697E-5*cz**2-3.2594E-5*cz**3
! HCHOA
	KR(1,jphoto_hchob)=   &
        -2.6855E-6+4.9102E-5*cz+5.7086E-5*cz**2-4.3040E-5*cz**3


! photolysis rates at 4 km:
! NO2
  	KR(2,jphoto_no2)=   &
        -1.3136E-3+2.4948E-2*cz-1.9513E-2*cz**2+6.611E-3*cz**3
! NO3
	KR(2,jphoto_no3)=   &
        1.59E-2+0.54202*cz-0.72079*cz**2+0.34898*cz**3
! O3
        KR(2,jphoto_o3b)=   &
          (1.6295E-7+cz*4.9940E-7-cz**2*2.9055E-5+1.8187E-4*cz**3   &
                    -1.5627E-4*cz**4+4.5975E-5*cz**5)
! HONO
	KR(2,jphoto_hono)=   &
        -2.6339E-4+4.6997E-3*cz-2.9408E-3*cz**2+7.4996E-4*cz**3
! HNO3
	KR(2,jphoto_hno3)=   &
        2.2106E-8-3.4422E-7*cz+1.8449E-6*cz**2-6.7994E-7*cz**3
! H2O2
	KR(2,jphoto_h2o2)=   &
        -4.73E-7+7.4881E-6*cz+9.7183E-6*cz**2-6.4955E-6*cz**3
! HNO4
	KR(2,jphoto_hno4)=   &
        -1.0672E-7+1.1698E-6*cz+1.9044E-5*cz**2-9.4072E-6*cz**3
! HCHOB
	KR(2,jphoto_hchoa)=   &
        -7.4493E-7+8.7149E-6*cz+7.1885E-5*cz**2-3.9526E-5*cz**3
! HCHOA
	KR(2,jphoto_hchob)=   &
        -5.1681E-6+8.4398E-5*cz+2.6478E-5*cz**2-3.4452E-5*cz**3


! photolysis rates at 8 km:
! NO2
       KR(3,jphoto_no2)=   &
       -1.3748E-3+2.9757E-2*cz-2.8355E-2*cz**2+1.1168E-2*cz**3
! NO3
       KR(3,jphoto_no3)=   &
       2.80132E-2+0.51381*cz-0.68839*cz**2+0.33448*cz**3
! O3
       KR(3,jphoto_o3b)=   &
       (1.6295E-7+cz*4.9940E-7-cz**2*2.9055E-5+1.8187E-4*cz**3   &
                    -1.5627E-4*cz**4+4.5975E-5*cz**5)
! HONO
       KR(3,jphoto_hono)=   &
       -3.1944E-4+6.0983E-3*cz-5.2694E-3*cz**2+1.9111E-3*cz**3
! HNO3
       KR(3,jphoto_hno3)=   &
       1.9176E-8-3.4083E-7*cz+2.1560E-6*cz**2-8.7941E-7*cz**3
! H2O2
       KR(3,jphoto_h2o2)=   &
       -7.6642E-7+1.1717E-5*cz+5.3611E-6*cz**2-4.9358E-6*cz**3
! HNO4
       KR(3,jphoto_hno4)=   &
       -3.2131E-7+3.6898E-6*cz+1.8481E-5*cz**2-9.826E-6*cz**3
! HCHOB
       KR(3,jphoto_hchoa)=   &
       -1.7563E-6+2.0714E-5*cz+6.5668E-5*cz**2-3.9386E-5*cz**3
! HCHOA
       KR(3,jphoto_hchob)=   &
       -7.9124E-6+1.258E-4*cz-2.8767E-5*cz**2-1.0505E-5*cz**3


! force all the kr to be non-negative
	do jdum = 1, nphoto
	do idum = 1, 3
	    kr(idum,jdum) = max( 0.0d0, kr(idum,jdum) )
	end do
	end do


! above 8km, use values at 8km
	if (zalt_m .ge. 8000) then
	 alpha = 1
	 k = 3
	 km1 = 3
! linear interpolation
	else if (zalt_m .lt. 8000. .and.   &
                 zalt_m .ge. 4000.) then
	 alpha = (8000. - zalt_m)/4000.
	 k = 3
	 km1 = 2
	else if (zalt_m .lt. 4000) then
	 alpha = (4000. - zalt_m)/4000
	 k = 2
	 km1 = 1
	end if

        a = alpha
        a1= 1. - alpha

      rk_photo(jphoto_no2)   = a1*kr(k,jphoto_no2)   +   &
                               a*kr(km1,jphoto_no2)
      rk_photo(jphoto_no3)   = a1*kr(k,jphoto_no3)   +   &
                               a*kr(km1,jphoto_no3)
      rk_photo(jphoto_o3a)   = 0.053*rk_photo(jphoto_no2)

      rk_photo(jphoto_o3b)   = a1*kr(k,jphoto_o3b)   +   &
                               a*kr(km1,jphoto_o3b)
      rk_photo(jphoto_hono)  = a1*kr(k,jphoto_hono)  +   &
                               a*kr(km1,jphoto_hono)
      rk_photo(jphoto_hno3)  = a1*kr(k,jphoto_hno3)  +   &
                               a*kr(km1,jphoto_hno3)
      rk_photo(jphoto_h2o2)  = a1*kr(k,jphoto_h2o2)  +   &
                               a*kr(km1,jphoto_h2o2)
      rk_photo(jphoto_hno4)  = a1*kr(k,jphoto_hno4)  +   &
                               a*kr(km1,jphoto_hno4)
      rk_photo(jphoto_hchoa) = a1*kr(k,jphoto_hchoa) +   &
                               a*kr(km1,jphoto_hchoa)
      rk_photo(jphoto_hchob) = a1*kr(k,jphoto_hchob) +   &
                               a*kr(km1,jphoto_hchob)

      rk_photo(jphoto_n2o5)   = 0.0   *rk_photo(jphoto_no2)
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
      end subroutine PhotoParam2
