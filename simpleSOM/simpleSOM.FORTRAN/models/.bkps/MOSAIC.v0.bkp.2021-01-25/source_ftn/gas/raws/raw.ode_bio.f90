! 1/23/2013 raz: fixed bugs in d_bio(ioh), d_bio(ino3).
! 1/23/2013 raz: fixed bugs in p_bio(ionit)
! 1/23/2013 raz: added OH, HO2, HCHO, and ALD2 sources from a-pinene + O3 reaction
!***********************************************************************
      subroutine ode_bio
  
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      ! real(r8) :: yield_factor, alphaO3_1, alphaO3_2, alphaO3_3, alphaO3_4
      ! real(r8) :: alphaOH_1, alphaOH_2, alphaOH_3, alphaOH_4
      ! real(r8) :: xonit

      ! xonit = 0.5

      ! p_bio(ino)= 0.0
      ! d_bio(ino)= r_bio(8)+r_bio(9)+r_bio(10)

      ! p_bio(ino2)= .91*r_bio(8)+1.2*r_bio(9)+r_bio(10) + (1.0 - xonit)*r_bio(19)*0.0 ! 1/23/2013 raz: added xonit ! r_bio(19) temporarily disabled here
      ! d_bio(ino2)= 0.0

      ! p_bio(ino3)= 0.0
      ! d_bio(ino3)= r_bio(3)+r_bio(7) + 0.0*r_bio(19)+r_bio(22)	! 1/23/2013 raz: fixed bug (19 and 22 were missing) ! r_bio(19) temporarily disabled here

      ! p_bio(ihno3)= .07*r_bio(7)
      ! d_bio(ihno3)= 0.0

      ! p_bio(io3)= 0.0
      ! d_bio(io3)= r_bio(2)+r_bio(6)+r_bio(18)+r_bio(21)

      ! p_bio(ioh)= .27*r_bio(2)+.27*r_bio(6) &
			! + 0.85*r_bio(18) 					! 1/23/2013 raz: added 0.85*r_bio(18)
      ! d_bio(ioh)= r_bio(1)+r_bio(5)+r_bio(17)+r_bio(20)	! 1/23/2013 raz: fixed bug (17 and 20 were missing)

      ! p_bio(iho2)= .07*r_bio(2)+.33*r_bio(4)+.1*r_bio(6)   &
              ! +.93*r_bio(7)+.91*r_bio(8)+.8*r_bio(9)+r_bio(10) &
		  ! + 0.24*r_bio(18) 					! 1/23/2013 raz: added 0.24*r_bio(18)
      ! d_bio(iho2)= r_bio(11)+r_bio(12)+r_bio(13)

      ! p_bio(ih2o2)= 0.0
      ! d_bio(ih2o2)= 0.0

      ! p_bio(ico)= .07*r_bio(2)+.33*r_bio(4)+.16*r_bio(6)   &
             ! +.64*r_bio(7)+.59*r_bio(10)
      ! d_bio(ico)= 0.0

      ! p_bio(ihcho)= .6*r_bio(2)+.2*r_bio(4)+.15*r_bio(6)   &
               ! +.28*r_bio(7)+.63*r_bio(8)+.25*r_bio(10)  &
		   ! + 0.37*r_bio(18)					! 1/23/2013 raz: added 0.37*r_bio(18)
      ! d_bio(ihcho)= 0.0

      ! p_bio(iald2)= .15*r_bio(2)+.07*r_bio(4)+.02*r_bio(6)   &
               ! +.28*r_bio(7)+.8*r_bio(9)+.55*r_bio(10)+r_bio(15)   &
               ! +.5*r_bio(16)  &
		   ! + 0.0*r_bio(18)					! 1/23/2013 raz: added 0.49*r_bio(18)
      ! d_bio(iald2)= 0.0

      ! p_bio(ipar)= 1.86*r_bio(7)+.18*r_bio(8)+1.6*r_bio(9)+   &
                   ! 2*r_bio(12)+2*r_bio(15)
      ! d_bio(ipar)= 0.0

      ! p_bio(iaone)= .03*r_bio(4)+.09*r_bio(6)+.63*r_bio(10)   &
               ! +.5*r_bio(16)
      ! d_bio(iaone)= 0.0

      ! p_bio(imgly)= .85*r_bio(6)+.34*r_bio(10)
      ! d_bio(imgly)= 0.0

      ! p_bio(ionit)= .93*r_bio(7)+.09*r_bio(8)+.8*r_bio(9)+r_bio(12)   &
               ! +r_bio(15) + xonit*r_bio(19)*0.0+r_bio(22)			! 1/23/2013 raz: added 19 and 22 ! r_bio(19) temporarily disabled here
      ! d_bio(ionit)= 0.0

      ! p_bio(ircooh)= .39*r_bio(2)+.46*r_bio(6)
      ! d_bio(ircooh)= 0.0

      ! p_bio(irooh)= r_bio(11)+r_bio(13)
      ! d_bio(irooh)= 0.0

      ! p_bio(ich3o2)= .7*r_bio(4)+.05*r_bio(6)
      ! d_bio(ich3o2)= 0.0

      ! p_bio(ic2o3)= .2*r_bio(2)+.97*r_bio(4)+.5*r_bio(5)   &
               ! +.11*r_bio(6)+.07*r_bio(7)
      ! d_bio(ic2o3)= 0.0

      ! p_bio(ixo2)= .08*r_bio(1)+.2*r_bio(2)+.2*r_bio(5)+.07*r_bio(6)   &
              ! +.93*r_bio(7)
      ! d_bio(ixo2)= 0.0

      ! p_bio(iisop)= 0.0
      ! d_bio(iisop)= r_bio(1)+r_bio(2)+r_bio(3)

      ! p_bio(iisoprd)= .65*r_bio(2)+.91*r_bio(8)+.2*r_bio(9)+r_bio(14)
      ! d_bio(iisoprd)= r_bio(4)+r_bio(5)+r_bio(6)+r_bio(7)

      ! p_bio(iisopp)= r_bio(1)
      ! d_bio(iisopp)= r_bio(8)+r_bio(11)+r_bio(14)

      ! p_bio(iisopn)= r_bio(3)
      ! d_bio(iisopn)= r_bio(9)+r_bio(12)+r_bio(15)

      ! p_bio(iisopo2)= .5*r_bio(5)
      ! d_bio(iisopo2)= r_bio(10)+r_bio(13)+r_bio(16)

!      p_bio(iapi)= 0.0
!      d_bio(iapi)= r_bio(17)+r_bio(18)+r_bio(19)*0.0 ! r_bio(19) temporarily disabled here

!      p_bio(ilim)= 0.0
!      d_bio(ilim)= r_bio(20)+r_bio(21)+r_bio(22)


! SORGAM
      ! yield_factor = 1.0 ! = original SORGAM

!FLAG1
p_bio(iCn3)= 9.27d-6*r_bio(1)       &
              +4.07d-5*r_bio(1)+1.29d-4*r_bio(1)+2.50d-4*r_bio(1)             &
              +7.96d-5*r_bio(1)+3.72d-6*r_bio(1)+2.44d-4*r_bio(1)             &
              +1.67d-4*r_bio(1)+5.55d-7*r_bio(1)+1.66d-5*r_bio(2)             &
              +7.29d-5*r_bio(2)+2.31d-4*r_bio(2)+4.47d-4*r_bio(2)             &
              +1.42d-4*r_bio(2)+6.65d-6*r_bio(2)+4.37d-4*r_bio(2)             &
              +2.99d-4*r_bio(2)+9.93d-7*r_bio(2)+2.97d-5*r_bio(3)             &
              +1.30d-4*r_bio(3)+4.14d-4*r_bio(3)+8.01d-4*r_bio(3)             &
              +2.55d-4*r_bio(3)+1.19d-5*r_bio(3)+7.82d-4*r_bio(3)             &
              +1.78d-6*r_bio(3)+5.32d-5*r_bio(4)+2.34d-4*r_bio(4)             &
              +7.41d-4*r_bio(4)+1.43d-3*r_bio(4)+4.56d-4*r_bio(4)             &
              +2.13d-5*r_bio(4)+3.18d-6*r_bio(4)+9.52d-5*r_bio(5)             &
              +4.18d-4*r_bio(5)+1.33d-3*r_bio(5)+2.57d-3*r_bio(5)             &
              +8.17d-4*r_bio(5)+5.70d-6*r_bio(5)+1.70d-4*r_bio(6)             &
              +7.48d-4*r_bio(6)+2.37d-3*r_bio(6)+4.59d-3*r_bio(6)             &
              +1.02d-5*r_bio(6)+3.05d-4*r_bio(7)+1.34d-3*r_bio(7)             &
              +4.25d-3*r_bio(7)+1.83d-5*r_bio(7)+5.46d-4*r_bio(8)             &
              +2.40d-3*r_bio(8)+3.27d-5*r_bio(8)+9.77d-4*r_bio(9)             &
              +5.85d-5*r_bio(9)+1.05d-4*r_bio(10)+1.88d-4*r_bio(11)           &
              +3.36d-4*r_bio(12)+6.01d-4*r_bio(13)+6.01d-4*r_bio(14)
      d_bio(iCn3)= r_bio(1)

      p_bio(iCn2)= 2.50d-3*r_bio(1)+5.35d-4*r_bio(3)+1.40d-3*r_bio(4)             &
              +3.82d-5*r_bio(5)+1.46d-3*r_bio(6)+8.22d-3*r_bio(7)             &
              +7.61d-3*r_bio(8)+4.29d-3*r_bio(9)+1.75d-3*r_bio(10)
      d_bio(iCn2)= r_bio(2)

      p_bio(iCn1)= 2.50d-3*r_bio(1)+2.50d-3*r_bio(2)+9.58d-4*r_bio(4)             &
              +2.51d-3*r_bio(5)+6.83d-5*r_bio(6)+2.62d-3*r_bio(7)             &
              +1.47d-2*r_bio(8)+1.36d-2*r_bio(9)+7.69d-3*r_bio(10)            &
              +3.13d-3*r_bio(11)
      d_bio(iCn1)= r_bio(3)

      p_bio(iC0)= 2.50d-3*r_bio(2)+2.49d-3*r_bio(3)+1.71d-3*r_bio(5)              &
             +4.49d-3*r_bio(6)+1.22d-4*r_bio(7)+4.69d-3*r_bio(8)              &
             +2.64d-2*r_bio(9)+2.44d-2*r_bio(10)+1.38d-2*r_bio(11)            &
             +5.61d-3*r_bio(12)
      d_bio(iC0)= r_bio(4)

      p_bio(iC1)= 2.49d-3*r_bio(3)+2.49d-3*r_bio(4)+3.07d-3*r_bio(6)              &
             +8.03d-3*r_bio(7)+2.19d-4*r_bio(8)+8.39d-3*r_bio(9)              &
             +4.72d-2*r_bio(10)+4.36d-2*r_bio(11)+2.46d-2*r_bio(12)           &
             +1.00d-2*r_bio(13)+1.00d-2*r_bio(14)
      d_bio(iC1)= r_bio(5)

      p_bio(iC2)= 2.49d-3*r_bio(4)+2.48d-3*r_bio(5)+5.49d-3*r_bio(7)              &
             +1.44d-2*r_bio(8)+3.92d-4*r_bio(9)+1.50d-2*r_bio(10)             &
             +8.45d-2*r_bio(11)+7.81d-2*r_bio(12)+4.41d-2*r_bio(13)+4.41d-2*r_bio(14)
      d_bio(iC2)= r_bio(6)

      p_bio(iC3)= 2.48d-3*r_bio(5)+2.46d-3*r_bio(6)+9.83d-3*r_bio(8)              &
             +2.57d-2*r_bio(9)+7.02d-4*r_bio(10)+2.69d-2*r_bio(11)            &
             +1.51d-1*r_bio(12)+1.40d-1*r_bio(13)+1.40d-1*r_bio(14)
      d_bio(iC3)= r_bio(7)

      p_bio(iC4)= 2.46d-3*r_bio(6)+2.42d-3*r_bio(7)+1.76d-2*r_bio(9)              &
             +4.61d-2*r_bio(10)+1.26d-3*r_bio(11)+4.81d-2*r_bio(12)           &
             +2.71d-1*r_bio(13)+2.71d-1*r_bio(14)
      d_bio(iC4)= r_bio(8)

      p_bio(iC5)= 2.42d-3*r_bio(7)+2.36d-3*r_bio(8)+3.15d-2*r_bio(10)             &
             +8.25d-2*r_bio(11)+2.25d-3*r_bio(12)+8.62d-2*r_bio(13)+8.62d-2*r_bio(14)
      d_bio(iC5)= r_bio(9)

      p_bio(iC6)= 2.36d-3*r_bio(8)+2.26d-3*r_bio(9)+5.64d-2*r_bio(11)             &
             +1.48d-1*r_bio(12)+4.02d-3*r_bio(13)+4.02d-3*r_bio(14)
      d_bio(iC6)= r_bio(10)

      p_bio(iC7)= 2.26d-3*r_bio(9)+2.06d-3*r_bio(10)+1.01d-1*r_bio(12)            &
             +2.64d-1*r_bio(13)+2.64d-1*r_bio(14)
      d_bio(iC7)= r_bio(11)

      p_bio(iC8)= 2.06d-3*r_bio(10)+1.72d-3*r_bio(11)+1.81d-1*r_bio(13)+1.81d-1*r_bio(14)
      d_bio(iC8)= r_bio(12)

      p_bio(iC9)= 1.72d-3*r_bio(11)+1.10d-3*r_bio(12)!+1.10d-3*r_bio(12)
      d_bio(iC9)= r_bio(13)
!FLAG2

      p_bio(iOH)= 0.0
      d_bio(iOH)= 0.0 !r_bio(1)+r_bio(2)+r_bio(3)+r_bio(4)+r_bio(5)+r_bio(6)+r_bio(7)    &
       !+r_bio(8)+r_bio(9)+r_bio(10)+r_bio(11)

      p_bio(iisop) = 0.0
      d_bio(iisop) = r_bio(14)

      return
      end subroutine ode_bio






