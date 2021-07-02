


      subroutine GasIntegrator(ntot,stot,t_in,t_out)
  
      USE mod_MAIN
      USE mod_GAS

      IMPLICIT NONE

      integer ntot
      real(r8) :: stot(ngas_max), t_in, t_out

! lsodes parameters
      integer itoler, itask, iopt, mf,   &
              nnz, lwm, nrdim, nidim
      parameter(itoler = 2, itask = 1, iopt = 1, mf = 222)
      parameter(nnz = ngas_max*ngas_max)
      parameter(lwm = 3*ngas_max*ngas_max + 12*ngas_max)
      parameter(nrdim = 20 + 9*ngas_max + lwm)
      parameter(nidim = 31 + ngas_max + ngas_max*ngas_max)

      real(r8) :: rwork(nrdim), atol(ngas_max)
      integer iwork(nidim), istate

!      dimension rwork(nrdim),iwork(nidim), atol(ngas_max)
      real(r8) :: rtol, t0, t1
      integer i, iflag
      external ode_gas,jac
	  
!	  print *, 'now in gasintegrator' ! edit wkc
	  
! set LSODE parameters...
! relative tolerance
      rtol=1.e-6
!
! absolute tolerances for gas-phase species
      do i=1,ngas_max
        atol(i)=1.e+1	! [molec/cc]
      enddo

      rwork = 0d0
      iwork = 0

      iwork(6) = 1000
      iwork(7) = 1
      istate   = 1
      rwork(6) = dt_sec

      iflag = 0
      if(iflag .eq. 0)then
      call xsetf(iflag)
      endif

      t0 = 0.0
      t1 = dt_sec

      call DLSODES(ode_gas,ntot,stot,t0,t1,itoler,   &
       rtol,atol,itask,istate,iopt,rwork,nrdim,iwork,nidim,jac,mf)

!
! error message from LSODES...
      if(istate.lt.0)then
        write(6,*)'WARNING from DLSODES: istate=',istate
        write(6,*)'Calling DLSODE to solve this gas integration step'
!
! reset time
        t0 = 0.0
        t1 = dt_sec
        istate = 1
!
        call MapGasSpecies(stot,0)	! map cnn into stot
        print *, 'this is timestep ', t1    ! edit wkc
        call DLSODE(ode_gas,ntot,stot,t0,t1,itoler,   &
        rtol,atol,itask,istate,iopt,rwork,nrdim,iwork,nidim,jac,22)

! error message from LSODE...
        if(istate.lt.0)then
        write(6,*)'WARNING from DLSODE: istate=',istate
        write(6,*)'Stopping in subroutine GasIntegrator'
        stop
        else
        write(6,*)'successful solution from DLSODE'
        endif
!
      else
!        write(6,*)'successful solution from DLSODES'
      endif
!
! successful solution...
      if(istate.eq.2)then
        do i=1,ntot
        stot(i)=max(0.0d0,stot(i))
        enddo
      endif

      return
      end subroutine GasIntegrator
