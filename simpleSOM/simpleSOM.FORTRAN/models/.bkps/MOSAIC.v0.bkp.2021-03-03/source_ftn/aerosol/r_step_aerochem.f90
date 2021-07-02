!    ======================================================================================
!                THIS SUBROUTINE STEPS FORWARD FOR AEROSOL-PHASE OLIGOMERIZATION
!    ======================================================================================

     SUBROUTINE STEP_AEROCHEM(OLIG,ko_f,ko_d,mdt)
     
     USE module_data_mosaic_main
     USE module_data_mosaic_aero
     
     IMPLICIT NONE

     INTEGER  :: OLIG
     REAL(R8) :: ko_f
     REAL(R8) :: ko_d
     REAL(R8) :: mdt

     REAL(R8),ALLOCATABLE :: k_ROLIG(:)
     REAL(R8),ALLOCATABLE :: k_PARMOLE_SOM(:)

     REAL(R8),ALLOCATABLE :: jk_PARMOLE_SOM_NEXT(:,:)
     REAL(R8),ALLOCATABLE :: jk_OLGMOLE_SOM_NEXT(:,:)
     REAL(R8),ALLOCATABLE :: jk_TOTMOLE_SOM(:,:)

     INTEGER :: i,j,k
     
!    SKIPPING CRITERIA:
!    ======================================================================================     
     IF (OLIG.EQ.0) RETURN
     
!    EFFECTIVE OLIGOMERIZATION RATE:
!    ======================================================================================
     ALLOCATE(k_ROLIG(nBINS))
     
     k_ROLIG = ko_f/(k_NUM*(1e6*PI/6.0*(k_DIAM_DRY**3.0)))
     
!    STEP FOR OLIGOMERIZATION:
!    ======================================================================================
     ALLOCATE(k_PARMOLE_SOM(nBINS))
     k_PARMOLE_SOM(:) = 0.d0
     
     k_PARMOLE_SOM = SUM(jk_PARMOLE_SOM,1)
     
     ALLOCATE(jk_PARMOLE_SOM_NEXT(nCOMP,nBINS))
     ALLOCATE(jk_OLGMOLE_SOM_NEXT(nCOMP,nBINS))
     jk_PARMOLE_SOM_NEXT(:,:) = 0.d0
     jk_OLGMOLE_SOM_NEXT(:,:) = 0.d0
     
     ALLOCATE(jk_TOTMOLE_SOM(nCOMP,nBINS))
     jk_TOTMOLE_SOM(:,:) = 0.d0
     
     jk_TOTMOLE_SOM(:,:) = jk_PARMOLE_SOM(:,:) + jk_OLGMOLE_SOM(:,:)

     DO 200 j = 1,nCOMP
     DO 201 k = 1,nBINS
        !IF (k_NUM(k).LE.10) GOTO 201
        
        jk_PARMOLE_SOM_NEXT(j,k) = (jk_PARMOLE_SOM(j,k) + mdt*0.5d0*ko_d*jk_TOTMOLE_SOM(j,k))/ &
                                   (1.d0 + mdt*k_ROLIG(k)*k_PARMOLE_SOM(k) + mdt*0.5d0*ko_d)

        jk_OLGMOLE_SOM_NEXT(j,k) = (jk_TOTMOLE_SOM(j,k)*(1.d0 + mdt*k_ROLIG(k)*k_PARMOLE_SOM(k)) - jk_PARMOLE_SOM(j,k))/ &
                                   (1.d0 + mdt*k_ROLIG(k)*k_PARMOLE_SOM(k) + mdt*0.5d0*ko_d)
           
201  CONTINUE
200  CONTINUE  
     
     jk_PARMOLE_SOM(:,:) = jk_PARMOLE_SOM_NEXT(:,:)
     jk_OLGMOLE_SOM(:,:) = jk_OLGMOLE_SOM_NEXT(:,:)
     
     RETURN
     END SUBROUTINE STEP_AEROCHEM
