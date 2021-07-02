      subroutine UpdateTime
        
      use module_data_mosaic_main
      use module_data_mosaic_pmcmos, only:  msolar_pmcmos_dtshift
      
      implicit none
      
      tsav_sec = tcur_sec
      told_sec = tcur_sec
      tcur_sec = tcur_sec + dt_sec
      tcur_min = tcur_sec/60.
      tcur_hrs = tcur_min/60.
      tmid_sec = told_sec + 0.5*dt_sec
      if (msolar_pmcmos_dtshift > 0) tmid_sec = tcur_sec + 0.5*dt_sec
      
      time_UTC = time_UTC + dt_sec/3600.
      
      if(time_UTC .ge. 24.0)then
        time_UTC = time_UTC - 24.0
      endif
      
      time_sec_old = time_sec
      time_sec = time_sec + dt_sec

      
      
      time_min = time_sec/60.
      time_hrs = time_min/60.

      if(tmid_sec .ge. tmar21_sec)then
        tmid_sec = tmid_sec - tmar21_sec		! seconds since noon, march 21
      else
        tmid_sec = tmid_sec + ((365-79)*24 - 12)*3600	! seconds since noon, march 21
      endif


      return
      end subroutine UpdateTime
