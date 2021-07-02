      subroutine SetRunParameters
      use module_data_mosaic_main

      implicit none

      real(r8) :: trunm, dtmin
      integer month, imonth, num_days
      real(r8) :: days(12)

      days(1)  = 31	! jan
      days(2)  = 28	! feb
      days(3)  = 31	! mar
      days(4)  = 30	! apr
      days(5)  = 31	! may
      days(6)  = 30	! jun
      days(7)  = 31	! jul
      days(8)  = 31	! aug
      days(9)  = 30	! sep
      days(10) = 31	! oct
      days(11) = 30	! nov
      days(12) = 31	! dec

      tmar21_sec = (79*24 + 12)*3600	! noon, mar 21

      month = tbeg_mo-1
      num_days = 0
      do imonth = 1, month
        num_days = num_days + days(imonth)
      enddo
      tbeg_sec = (num_days+tbeg_dd-1)*24*3600 +   &
                 tbeg_hh*3600 + tbeg_mm*60 + tbeg_ss	! time since beginning of year in seconds

      trun_sec = ((trun_dd*24+trun_hh)*60+trun_mm)*60+trun_ss 	! total runtime duration in seconds

      dt_sec   = 60.*dt_min				! time step [seconds]
      nstep = INT(float(trun_sec)/dt_sec + 0.01)	! number of time steps

      tcur_sec = tbeg_sec		! initialize current time to tbeg_sec
      tcur_min = tcur_sec/60.
      tcur_hrs = tcur_min/60.

      time_UTC = float(tbeg_hh)      +	   &  ! initialize time of day in hours (UTC)
                 float(tbeg_mm)/60.  +   &
                 float(tbeg_ss)/3600.
      time_UTC_beg = time_UTC

      time_sec = 0.0		! initialize time since start to zero
      time_sec_old = time_sec
      time_min = 0.0
      time_hrs = 0.0
      t_since_start = 0.0	! sec
      it       = 0		! initialize step counter to zero

!
! convert rlon and rlat to radians
      rlon = rlon*deg2rad
      rlat = rlat*deg2rad
!
      if(msolar .eq. 2)then
	write(6,*)'  '
        write(6,*)' msolar = 2: Constant photolysis rates set for the entire simulation'
        write(6,*)' Set the values in file: PhotoRateConstants.f90'
      endif
!
      return
      end subroutine SetRunParameters
