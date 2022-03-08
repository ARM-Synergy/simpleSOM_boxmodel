!    ======================================================================================
!                   THIS SUBROUTINE STEPS FORWARD FOR HOMOGENEOUS NUCLEATION
!    ======================================================================================

     SUBROUTINE STEP_NUCL(NUC,NUC_T0,NUC_T1,mdt)
     
     USE mod_MAIN
     USE module_data_mosaic_aero
     
     IMPLICIT NONE
     
     REAL(R8) :: NUC
     REAL(R8) :: NUC_T0
     REAL(R8) :: NUC_T1
     REAL(R8) :: mdt

     REAL(R8) :: nOFFSET
     REAL(R8) :: iBIN

     IF (time_hrs.LT.(NUC_T0/3600.)) RETURN
     IF (time_hrs.GT.(NUC_T1/3600.)) RETURN
     
     iBIN = 1
     
     noffset = ngas_max + naer_tot*(ibin - 1)
     
     cnn(noffset + knum_a) = &
     cnn(noffset + knum_a) + NUC*mdt     
                        
     cnn(nOFFSET+kwater_a+iCN3_a+1-1) = &
     cnn(nOFFSET+kwater_a+iCN3_a+1-1) + &
     NUC*mdt*(1./6.*PI*(Dp_dry_a(iBIN)*1e-2)**3.)*1000./(350.*1.0958*1e-3)*1e6*1e6
     
     
     RETURN
     END SUBROUTINE STEP_NUCL
