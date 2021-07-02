
!***************************************************************************
! subroutine ode_gas
!
! purpose: computes time derivatives of species concentrations ds/dt = sdot(*).
!          calls ode_com, ode_urb, ode_bio, ode_mar depending on the
!          chemical regime (iregime)
!
! author: Rahul A. Zaveri
!
!---------------------------------------------------------------------------
!
      subroutine ode_gas(ntot,tt,s,sdot)
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      integer ntot, i
      real(r8) :: s(ntot_max),sdot(ntot_max), tt
!
      do i=1,nrxn_com
        r_com(i) = 0.
      enddo

      do i=1,nrxn_urb
        r_urb(i) = 0.
      enddo

      do i=1,nrxn_bio
        r_bio(i) = 0.
      enddo

      do i=1,nrxn_mar
        r_mar(i) = 0.
      enddo

      do i=1,nrxn_het
        r_het(i) = 0.
      enddo

      call GasRates(s)
!
      do i=1,ngas_max
        p_com(i) = 0.
        p_urb(i) = 0.
        p_bio(i) = 0.
        p_mar(i) = 0.
        p_het(i) = 0.

        d_com(i) = 0.
        d_urb(i) = 0.
        d_bio(i) = 0.
        d_mar(i) = 0.
        d_het(i) = 0.
      enddo
!
!
      goto (1,2,3,4,5,6), iregime

1     call ode_com
      call ode_het

      do i=1,nreg1
        sdot(i) =     ( (p_com(i)+p_het(i)) -   &
                        (d_com(i)+d_het(i)) )   &
                + emit(i)
      enddo

      return
!
!----------------------------------------------------------
2     call ode_com
      call ode_urb
      call ode_het

      do i=1,nreg2
        sdot(i) =     ( (p_com(i)+p_urb(i)+p_het(i)) -   &
                        (d_com(i)+d_urb(i)+d_het(i)) )   &
                + emit(i)
      enddo

      return
!
!----------------------------------------------------------
3     call ode_com
      call ode_urb
      call ode_bio
      call ode_het

      do i=1,nreg3
        sdot(i) =     ( (p_com(i)+p_urb(i)+p_bio(i)+p_het(i)) -   &
                        (d_com(i)+d_urb(i)+d_bio(i)+d_het(i)) )   &
                + emit(i)
      enddo

      return
!
!----------------------------------------------------------
4     call ode_com
      call ode_mar
      call ode_het

      do i=1,nreg4
        sdot(i) =     ( (p_com(i)+p_mar(i)+p_het(i)) -   &
                        (d_com(i)+d_mar(i)+d_het(i)) )   &
                + emit(i)
      enddo

      return
!
!----------------------------------------------------------
5     call ode_com
      call ode_urb
      call ode_mar
      call ode_het

      do i=1,nreg5
        sdot(i) =     ( (p_com(i)+p_urb(i)+p_mar(i)+p_het(i)) -   &
                        (d_com(i)+d_urb(i)+d_mar(i)+d_het(i)) )   &
                + emit(i)
      enddo

      return
!
!----------------------------------------------------------
6     call ode_com
      call ode_urb
      call ode_bio
      call ode_mar
      call ode_het

      do i=1,nreg6
        sdot(i) =     ( (p_com(i)+p_urb(i)+p_bio(i)+p_mar(i)+   &
                             p_het(i)) -   &
                        (d_com(i)+d_urb(i)+d_bio(i)+d_mar(i)+   &
                             d_het(i)) )   &
                + emit(i)
      enddo

      return

      end subroutine ode_gas
