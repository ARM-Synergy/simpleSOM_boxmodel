!    ======================================================================================
!                        THIS SUBROUTINE STEPS FORWARD FOR VOC OXIDATION
!    ======================================================================================

     SUBROUTINE STEP_VOCCHEM(VOC_kOH,OH,mdt)
     
     USE mod_MAIN, ONLY: R8,nCOMP,CNN,kCN3
     USE mod_GAS,  ONLY: VOC,VOC_COEFF,VOC_OXY,j_GASMOLE_OXY
     
     IMPLICIT NONE

     REAL(R8) :: VOC_kOH
     REAL(R8) :: OH
     REAL(R8) :: mdt
     REAL(R8) :: dVOC

     REAL(R8) :: dSOM(nCOMP)
     REAL(R8) :: dOXY(nCOMP)
     
     INTEGER :: i,j,k
     
!    REACTED VOC:
     dVOC = VOC_kOH*OH*mdt*VOC

!    UPDATE VOC CONC.:
     VOC = VOC - dVOC

!    UPDATE SIMPLE-SOM PRODUCTS:
     dSOM = dVOC*VOC_COEFF*1e-9*(101325./8.314/298.)*1e-6*6.022e23
     
     DO j = 1,nCOMP
        CNN(kCN3+j-1) = &
        CNN(kCN3+j-1) + dSOM(nCOMP-j+1)*0.4
     END DO
     
     DO j = 1,nCOMP-1
        CNN(kCN3+j-1) = &
        CNN(kCN3+j-1) + dSOM(nCOMP-j-1+1)*0.6
     END DO
     
!    UPDATE SIMPLE-SOM PRODUCTS:
     dOXY = dVOC*VOC_OXY*1e-9*(101325./8.314/298.)*1e-6*6.022e23
     
     DO j = 1,nCOMP
        j_GASMOLE_OXY(j) = &
        j_GASMOLE_OXY(j) + dOXY(nCOMP-j+1)*0.4
     END DO

     DO j = 1,nCOMP-1
        j_GASMOLE_OXY(j) = &
        j_GASMOLE_OXY(j) + dOXY(nCOMP-j-1+1)*0.6
     END DO

     RETURN
     END SUBROUTINE STEP_VOCCHEM
