!    ======================================================================================
!               THIS SUBROUTINE STEPS FORWARD FOR OXYGEN CHEMICAL TRANSFORMATION
!    ======================================================================================

     SUBROUTINE STEP_OXYCHEM(mdt)

     USE mod_MAIN, ONLY: R8,nCOMP
     USE mod_GAS
     
     REAL(R8) :: mdt
     
     INTEGER :: i,j,k
     
!     PRINT*, p_OXY
!     PRINT*, d_OXY

!     READ(*,*)


     DO i = 1,nCOMP
        
        j_GASMOLE_OXY(i) = &
        j_GASMOLE_OXY(i) + p_OXY(i)*mdt - d_OXY(i)*mdt

     END DO
     
     RETURN
     END SUBROUTINE STEP_OXYCHEM
