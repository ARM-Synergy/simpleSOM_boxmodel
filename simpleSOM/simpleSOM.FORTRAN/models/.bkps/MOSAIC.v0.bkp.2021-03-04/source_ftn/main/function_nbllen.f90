!
!   returns the position of the last non-blank character in str
!
	integer function nbllen( str )
        implicit none
	character*(*) str
        integer j

	j = len(str)

	if (j .gt. 0) then
1000	    if (str(j:j) .eq. ' ') then
		j = j - 1
		if (j .gt. 0) goto 1000
	    end if
	end if
	nbllen = max0( j, 0 )

	return
	end function nbllen
