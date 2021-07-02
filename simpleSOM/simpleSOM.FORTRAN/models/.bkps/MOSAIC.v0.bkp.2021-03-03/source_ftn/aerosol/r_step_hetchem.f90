!    ======================================================================================
!                   THIS SUBROUTINE STEPS FORWARD FOR HETEROGENEOUS CHEMISTRY
!    ======================================================================================

     SUBROUTINE STEP_HETCHEM(HETCHEM,Gamma_OH,OH,Temp,mdt)

     USE module_data_mosaic_main
     USE module_data_mosaic_aero
     
     IMPLICIT NONE

     INTEGER  :: HETCHEM
     REAL(R8) :: Gamma_OH
     REAL(R8) :: OH
     REAL(R8) :: Temp
     REAL(R8) :: mdt
     
     REAL(R8),PARAMETER :: DG_OH = 3e-5
     REAL(R8),PARAMETER :: MW_OH = 17.
     REAL(R8),PARAMETER :: RGAS  = 8.314

     REAL(R8),PARAMETER :: j_PFRAG = 0.9 !0.5
     REAL(R8),PARAMETER :: j_PFUNC = 1 - j_PFRAG
     
     REAL(R8) :: MFP_OH
     REAL(R8) :: RMS_OH
     
     REAL(R8) :: k_KNUD_OH(nBINS)
     REAL(R8) :: k_FUCH_OH(nBINS)
     REAL(R8) :: k_COEF_OH(nBINS)

     REAL(R8),ALLOCATABLE :: jk_PARMOLE_SOM_PREV(:,:)
     REAL(R8),ALLOCATABLE :: k_PARMOLE_SOM(:)
     
     INTEGER :: i,j,k
     
!    SKIPPING CRITERIA:
!    ======================================================================================
     IF (HETCHEM.EQ.0) RETURN
     
!    CALCULATE THE HET. CHEM. RATE COEFFICIENT:
!    ======================================================================================     
!    MEAN FREE PATH OF OH:
     MFP_OH = 3.*DG_OH/SQRT(8.*RGAS*Temp/PI/(MW_OH/1000.))
     
!    ROOT-MEAN-SQUARE SPEED OF OH:
     RMS_OH = SQRT(RGAS*Temp/(MW_OH/1000.)/2./PI)
     !RMS_OH = SQRT(RGAS*Temp/(MW_OH/1000.))
     
!    KNUDSEN NUMBER:
     k_KNUD_OH = 2.*MFP_OH/k_DIAM_DRY

!    FUCHS-SUTUGIN FACTOR: 
     k_FUCH_OH = (1. + k_KNUD_OH)/(1. + k_KNUD_OH + 0.75/k_KNUD_OH)

!    HET. CHEM. RATE COEFFICIENT [# cm-3 s-1]:
     k_COEF_OH = Gamma_OH*(k_NUM*PI*(k_DIAM_DRY**2.)*k_FUCH_OH)*(0.25*OH*1e6*RMS_OH)
     
!    SOLVE HET. CHEM.:
!    ======================================================================================
!    KEEP PREVIOUS VALUES: 
     ALLOCATE(jk_PARMOLE_SOM_PREV(nCOMP,nBINS))
     jk_PARMOLE_SOM_PREV(:,:) = 0.d0

     DO j = 1,nCOMP
        DO k = 1,nBINS
           jk_PARMOLE_SOM_PREV(j,k) = jk_PARMOLE_SOM(j,k)
        END DO
     END DO

!    TOTAL SOM MOLES IN EACH BIN:
     ALLOCATE(k_PARMOLE_SOM(nBINS))
     k_PARMOLE_SOM(:) = 0.d0
     
     k_PARMOLE_SOM = MAX(SUM(jk_PARMOLE_SOM,1),1.)
     
!    LOOP OVER SPECIES:       
     DO 200 j = 1,nCOMP
     DO 201 k = 1,nBINS
        IF (j.EQ.nCOMP) THEN
           jk_PARMOLE_SOM(j,k) = (jk_PARMOLE_SOM_PREV(j,k) + &
                                  jk_PARMOLE_SOM_PREV(j-1,k)*k_COEF_OH(k)*mdt/k_PARMOLE_SOM(k)*j_PFRAG)/ &
                                 (1. + k_COEF_OH(k)*mdt/k_PARMOLE_SOM(k))
        ELSE IF (j.EQ.1) THEN
           jk_PARMOLE_SOM(j,k) = (jk_PARMOLE_SOM_PREV(j,k) + &
                                  jk_PARMOLE_SOM_PREV(j,k)*k_COEF_OH(k)*mdt/k_PARMOLE_SOM(k)*j_PFUNC +   &
                                  jk_PARMOLE_SOM_PREV(j+1,k)*k_COEF_OH(k)*mdt/k_PARMOLE_SOM(k)*j_PFUNC)/ &
                                 (1. + k_COEF_OH(k)*mdt/k_PARMOLE_SOM(k))
        ELSE
           jk_PARMOLE_SOM(j,k) = (jk_PARMOLE_SOM_PREV(j,k) + &
                                  jk_PARMOLE_SOM_PREV(j+1,k)*k_COEF_OH(k)*mdt/k_PARMOLE_SOM(k)*j_PFUNC + &
                                  jk_PARMOLE_SOM_PREV(j-1,k)*k_COEF_OH(k)*mdt/k_PARMOLE_SOM(k)*j_PFRAG)/ &
                                 (1. + k_COEF_OH(k)*mdt/k_PARMOLE_SOM(k))
        END IF
201  CONTINUE
200  CONTINUE
           
     RETURN
     END SUBROUTINE STEP_HETCHEM
