
      SUBROUTINE IntegrateChemistry(timer)
  
      USE mod_MAIN
      USE mod_GAS
      use module_data_mosaic_aero
      use module_data_mosaic_cloud
      
      implicit none

      REAL(R8) :: timer
      
      real(r8) :: t_in, t_out
      
      REAL(R8) :: Gamma_OH,OH,Temp,mdt
      
      INTEGER :: i,j,k
      
      REAL(R8),PARAMETER :: NA = 6.022d23
      
      INTEGER :: HETCHEM
      INTEGER :: OLIG
      
      INTEGER :: OFFSET
      
      REAL(R8) :: ko_f
      REAL(R8) :: ko_d
      REAL(R8) :: NUC
      REAL(R8) :: NUC_T0,NUC_T1

      INTEGER nOFFSET,iBIN

      REAL(R8) TOT
      
      t_in = told_sec
      t_out= tcur_sec
      
      if(mgas.eq.1)then
        call GasChemistry(t_in, t_out)
      endif


      IBIN = 1
      
      noffset = ngas_max + naer_tot*(ibin - 1)

      !PRINT*, 'HERE in INTEGRATE'
      !PRINT*, 'mass = ',cnn(nOFFSET+kwater_a+iCN3_a+1-1)
      !READ(*,*)



      
      
      if(maer.eq.1)then
        call AerChemistry(t_in, t_out)
      endif




      
      OLIG = OLIG_IN
      
      ko_f = ko_f_in
      ko_d = ko_d_in

      HETCHEM = HETCHEM_IN
      
      Gamma_OH = 1.

      OH = CNN(kOH)
      
      Temp = 298.

      mdt = t_out - t_in

      NUC = NUC_in

      NUC_T0 = NUC_T0_in
      NUC_T1 = NUC_T1_in
      
      DO j = 1,nCOMP
         DO k = 1,nBINS
            OFFSET = ngas_max + naer_tot*(k - 1)
            
            
            jk_PARMOLE_SOM(j,k) = cnn(OFFSET+kwater_a+iCN3_a+j-1)*1d3*1d-9*1d-6*NA
            
         END DO
      END DO
            
      
      CALL STEP_HETCHEM(HETCHEM,Gamma_OH,OH,Temp,mdt)
      
      CALL STEP_AEROCHEM(OLIG,ko_f,ko_d,mdt)

      
      
      DO j = 1,nCOMP
         DO k = 1,nBINS
            OFFSET = ngas_max + naer_tot*(k - 1)
                        
            cnn(OFFSET+kwater_a+iCN3_a+j-1) = jk_PARMOLE_SOM(j,k)*(1d-3*1d9*1d6)/NA 
            
         END DO
      END DO

      CALL STEP_NUCL(NUC,NUC_T0,NUC_T1,mdt)

      !TOT = 0.0

      !DO k = 1,nBINS
     
      !   noffset = ngas_max + naer_tot*(k - 1)

      !   TOT = TOT + cnn(noffset + knum_a)

      !END DO

      !PRINT*, 'HERE in INTEGRATE'
      !PRINT*, 'TOT = ',TOT
      !READ(*,*)
      
      if(mcld.eq.1)then
        call CldChemistry(t_in, t_out)
      endif
      
      
      return
      end subroutine IntegrateChemistry
