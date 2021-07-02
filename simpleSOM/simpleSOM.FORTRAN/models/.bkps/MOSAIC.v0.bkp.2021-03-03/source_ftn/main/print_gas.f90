! Gas-phase species
! subroutine for printing output at iprint time steps

      subroutine print_gas
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      integer l
      real(r8) :: cppb(ntot_max)
	real(r8) :: cppb_api1_0
	save cppb_api1_0


      if (iwrite_gas <= 0) return


      do l=1,ngas_max
        cppb(l) = (cnn(l)/cair_mlc)*ppb	! converting (molecules/cc) to (ppb)
      enddo


!      if(it.eq.0)write(lun_gas,201)(species(l), l=1, ngas_max)
!201   format(' UTC(hr)    t(hr)   Temp(K)    Pr(atm)   RH(%)',   &
!             '  AIR(molec/cc) H2O(molec/cc)', 100(2x,a11))

!      if(it.eq.0)write(lun_gas,201)species(kapi1)

!	if(it .eq. 0)cppb_api1_0 = cppb(kapi1)
!	if(it .eq. 0)write(6,*)'cppb_api1_0 = ', cppb_api1_0*pr_atm*1.e5/(82.056*te)
!201   format(' t(hr)   ', 100(2x,a11))



!      write(lun_gas,202)time_UTC,time_hrs,te,pr_atm,RH,cair_mlc,h2o,   &
!                   (cppb(l),l=1,ngas_max)
!202   format(f8.4,4(2x,f8.4),100(2x,e11.5))




!      write(lun_gas,202)time_hrs, cppb(kapi1)*pr_atm*1.e5/(82.056*te)
!202   format(f8.4,2x, f12.5, 100(2x,e11.5))

!	delta_api1 = (cppb_api1_0-cppb(kapi1))*pr_atm*1.e5/(82.056*te)
!	if(it .eq. nstep)write(6,*)'cppb_api1_0 = ',cppb_api1_0*pr_atm*1.e5/(82.056*te)
!	if(it .eq. nstep)write(6,*)'cppb_api1_f = ',cppb(kapi1)*pr_atm*1.e5/(82.056*te)
!	if(it .eq. nstep)write(6,*)nstep, 'Delta api1 (ug/m^3) = ', delta_api1
!
!	if(it .eq. nstep)write(6,*)' '
!      close(20)


      return
      end subroutine print_gas

