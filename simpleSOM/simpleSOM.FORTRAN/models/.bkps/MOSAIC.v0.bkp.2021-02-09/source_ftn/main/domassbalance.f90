      subroutine DoMassBalance
      use module_data_mosaic_main

      implicit none

      real(r8) :: totS, totN, totCl, totNH4
      real(r8) :: tNid, tSid, tClid, tNH4id
      integer i

!
!=================================================================
! Mass Balances
!=================================================================
!     initialize...
      totS = 0.0
      totN = 0.0
      totCl= 0.0
!
!__________________________________________________________________
! Sulfur Balance
!
      do i=kdms,ksulfhox
      totS = totS + cnn(i)
      enddo
      totS = totS + cnn(kso2) + cnn(kh2so4)
!     &     + cnn(ktso4_a1)+cnn(ktso4_a2)+cnn(ktso4_a3)+cnn(ktso4_a4)

!__________________________________________________________________
! Nitrogen Balance
!
      do i=kno,khno4
      totN = totN + cnn(i)
      enddo
      totN = totN + cnn(kn2o5)+cnn(kpan)+cnn(konit)+cnn(knap)
!     &     + cnn(ktno3_a1)+cnn(ktno3_a2)+cnn(ktno3_a3)+cnn(ktno3_a4)
!
!__________________________________________________________________
! Chlorine Balance
!
      totCl = cnn(khcl)
!     &      + cnn(ktcl_a1)+cnn(ktcl_a2)+cnn(ktcl_a3)+cnn(ktcl_a4)
!
!__________________________________________________________________
!
! Ammonium Balance
      totNH4= cnn(knh3)
!     &     + cnn(ktnh4_a1)+cnn(ktnh4_a2)+cnn(ktnh4_a3)+cnn(ktnh4_a4)


!
!  Initial total elements (N,S,Cl)
      if(it.eq.0) then

          if(totN.gt.1.E-30)then
          tNi = totN
          tNid= totN
          else
          tNi = totN
          tNid= 1.0
          endif

          if(totS.gt.1.E-30)then
          tSi = totS
          tSid= totS
          else
          tSi = totS
          tSid= 1.0
          endif

          if(totCl.gt.1.E-30)then
          tCli = totCl
          tClid= totCl
          else
          tCli = totCl
          tClid= 1.0
          endif


          if(totNH4.gt.1.E-30)then
          tNH4i = totNH4
          tNH4id= totNH4
          else
          tNH4i = totNH4
          tNH4id= 1.0
          endif


      end if

!
!  Calculate percent deviation in elemental mass balance
      DN = 100.*(totN-tNi)/tNid
      DS = 100.*(totS-tSi)/tSid
      DCl= 100.*(totCl-tCli)/tClid
      DNH4= 100.*(totNH4-tNH4i)/tNH4id

      return
      end subroutine DoMassBalance
