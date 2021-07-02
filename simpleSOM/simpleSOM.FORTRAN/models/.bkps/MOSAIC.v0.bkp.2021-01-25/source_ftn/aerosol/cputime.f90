	subroutine cputime( icpuflag )

	use module_data_mosaic_kind, only:  r8

	implicit none

!	include 'chemistry.com'

	integer icpuflag

	integer, parameter :: ncpugrps = 5
	real(r8), save :: tcpugrp(ncpugrps)
	character(len=8), save :: cpugrpnm(ncpugrps)
	data cpugrpnm /   &
              'PSC ', 'KM  ', 'TEM ', 'MET ', 'Other' /

	integer i
	real(r8) :: tcpu, tcpu00, tcpuold, tcputot

	call cpusecnd( tcpu )

	if (icpuflag .eq. 0) then
	    tcpu00 = tcpu
	    do 1000 i = 1, ncpugrps
		tcpugrp(i) = 0.
1000	    continue

	else if (icpuflag .gt. 0) then
	    i = min( icpuflag, ncpugrps )
	    tcpugrp(i) = tcpugrp(i) + (tcpu - tcpuold)

	else
	    tcputot = max( (tcpu - tcpu00), 1.d-10 )

	    write(6,9100) tcputot, (tcputot/60.)
	    do 2000 i = 1, ncpugrps
	         write(6,9110) i, cpugrpnm(i),   &
      			tcpugrp(i), (tcpugrp(i)/tcputot),   &
                        tcpugrp(i)/tcpugrp(3)
2000	    continue

	end if

	tcpuold = tcpu

9100	format( / ' MOSAIC cpu time statistics' /   &
           '     total time (secs, mins)', 2f12.2 /   &
           '     module        secs     secs/total  Time wrt TEM' )
9110	format( i3, 2x, a, f10.4, 2x, f11.5, 2x, f9.3 )

	return
	end subroutine cputime




	subroutine cpusecnd( tcpu )
!
!   returns elapsed cpu time in seconds
!   uses ibm function mclock which returns cpu time in hundreths of seconds
!
	use module_data_mosaic_kind, only:  r4, r8

	implicit none

!   subr. parameters
	real(r8) :: tcpu

!   local variables
	real(r4) :: dum, tcpu_2d(2)
	real(r4) :: etime
	external etime

	dum = etime(tcpu_2d)
	tcpu = tcpu_2d(1) + tcpu_2d(2)

	return
	end subroutine cpusecnd
