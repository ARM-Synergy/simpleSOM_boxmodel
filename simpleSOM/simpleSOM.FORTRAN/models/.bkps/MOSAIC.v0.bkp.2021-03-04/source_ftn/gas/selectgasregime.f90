!***********************************************************************
! subroutine SelectGasRegime: selects an optimum combination of gas-phase
!                             mechanisms
!
! input : cnn       = full species concentrations array (molec/cc)
!
! output: iregime   = 1     : com
!                   = 2     : com + urb
!                   = 3     : com + urb + bio
!                   = 4     : com + mar
!                   = 5     : com + urb + mar
!                   = 6     : com + urb + bio + mar
!         ngas      = number of gas-phase species in the selected mechanism
!
! author: Rahul A. Zaveri
! date  : february 1996
!
!---------------------------------------------------------------------
      subroutine SelectGasRegime(ntot)
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      integer ntot
      integer m_com, m_urb, m_bio, m_mar, k
      real(r8) :: cutoff_conc

      cutoff_conc = 5.e+6     ! [molec/cc]

!
!---------------------------------------------
! initialize regime flags to zero...
      m_com = 1 ! 1 (always)
      m_urb = 0	! 0 or 1
      m_bio = 2	! 0 or 2 ! edit wkc
      m_mar = 0	! 0 or 3

!
! decide mechanism flags...
      do k = kpar, kxpar
      if( (cnn(k) .gt. cutoff_conc)  .or.   &
          (emission(k) .gt. 0.0)    ) m_urb=1
      enddo
!
      do k = kisop, kC9 ! edit wkc to add species
      if( (cnn(k) .gt. cutoff_conc)  .or.   &
          (emission(k) .gt. 0.0)    ) m_bio=2
      enddo
!
      do k = kdms, ksulfhox
      if( (cnn(k) .gt. cutoff_conc)  .or.   &
          (emission(k) .gt. 0.0)     ) m_mar=3
      enddo
!
      iregime = m_com + m_urb*((2-m_bio)/2) + m_bio + m_mar

!	  print *, 'm_bio = ', m_bio
!	  print *, 'iregime = ', iregime

!-------------------------------------------------------------------------
      goto (1,2,3,4,5,6), iregime

1     ntot = nreg1
      return

2     ntot = nreg2
      return

3     ntot = nreg3
      return

4     ntot = nreg4
      return

5     ntot = nreg5
      return

6     ntot = nreg6
      return

      end subroutine SelectGasRegime

