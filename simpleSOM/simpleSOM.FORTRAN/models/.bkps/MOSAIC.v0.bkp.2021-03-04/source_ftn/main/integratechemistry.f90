
      SUBROUTINE IntegrateChemistry(timer)
  
      use module_data_mosaic_main
      use module_data_mosaic_gas
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
      
      t_in = told_sec
      t_out= tcur_sec
      
      if(mgas.eq.1)then
        call GasChemistry(t_in, t_out)
      endif


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

      
      if(mcld.eq.1)then
        call CldChemistry(t_in, t_out)
      endif
      
      
      return
      end subroutine IntegrateChemistry
