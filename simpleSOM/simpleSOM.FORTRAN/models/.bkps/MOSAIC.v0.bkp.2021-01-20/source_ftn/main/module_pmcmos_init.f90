	module module_pmcmos_init


	use module_data_mosaic_kind, only:  r8
	use module_data_mosaic_pmcmos


	implicit none


	contains


!-----------------------------------------------------------------------
	subroutine pmcmos_init

	use module_data_mosaic_main, only:  &
	    cnn, lun_sect_190, ngas_max, pblh, pr_atm,   &
	    rh, species, te, time_sec
	use module_pmcmos_subaa, only:  pmcmos_update_met

!   subr arguments (none)

!   local variables
	integer  :: i


	call pmcmos_read_inputs

	call pmcmos_update_met( 0 )

!   initial gas mixing ratios (ppb)
	do i = 1, ngas_max
	    cnn(i) = gas_init_val(i)
	end do
!	write(*,'(a,1p,5e13.5)') 'co,so2,ch4,c2h6', cnn(17:20)
	write(lun_sect_190,'(i3,2x,a)') (i, species(i)(1:12), i=1,ngas_max)
	write(lun_sect_190,'(/a,5(1x,f14.6))') 'time(h)', time_sec/3600.0, &
		te, pr_atm*1013.25, rh, pblh
	write(lun_sect_190,'(1p,5e14.6)') cnn(1:ngas_max)

!   initial aerosol concentrations


	return
	end subroutine pmcmos_init


!-----------------------------------------------------------------------
	subroutine pmcmos_read_inputs

	use module_data_mosaic_main, only:  &
	    ngas_max, species

!   subr arguments (none)

!   local variables
	integer, parameter :: nparam_max = 2+ngas_max
	integer  :: i, idiagbb, idiagcc, iok, itmpa
	integer  :: imustfind(nparam_max), ifound(nparam_max)
	integer  :: itmpvecaa(nparam_max)
	integer  :: j
	integer  :: nparam
	integer  :: ptype_aa(nparam_max)
	real(r8) :: pval_aa(ntmax_pmcmos,nparam_max)
	character(len=20) :: pname_aa(nparam_max)
	character(len=80) :: ptxt_aa(ntmax_pmcmos,nparam_max)


	idiagbb = -1
	idiagcc =  1

!   temperature versus time
	nparam = 2
	pname_aa(1) = 'time'
	pname_aa(2) = 'temp'
	ptype_aa(:) = 2
	imustfind(1:2) = 1
	call pmcmos_readfile_aa( 0, iok, &
	    temp_profile_fname, lun_pmcmos_in1, &
	    nparam, temp_nt, ntmax_pmcmos, &
	    imustfind, ifound, pname_aa, ptype_aa, ptxt_aa, pval_aa )
	allocate( temp_tim(1:temp_nt) )
	allocate( temp_val(1:temp_nt) )
	temp_tim(1:temp_nt) = pval_aa(1:temp_nt,1) 
	temp_val(1:temp_nt) = pval_aa(1:temp_nt,2) 
	    
	if (idiagbb > 0) then
	write(*,'(/a,i5)'  ) 'temp profile', temp_nt
	write(*,'(10f10.1)') temp_tim(1:temp_nt)
	write(*,'(10f10.2)') temp_val(1:temp_nt)
	end if ! (idiagbb > 0)


!   pblh versus time
	nparam = 2
	pname_aa(1) = 'time'
	pname_aa(2) = 'height'
	ptype_aa(:) = 2
	imustfind(1:2) = 1
	call pmcmos_readfile_aa( 0, iok, &
	    pblh_profile_fname, lun_pmcmos_in1, &
	    nparam, pblh_nt, ntmax_pmcmos, &
	    imustfind, ifound, pname_aa, ptype_aa, ptxt_aa, pval_aa )
	allocate( pblh_tim(1:pblh_nt) )
	allocate( pblh_val(1:pblh_nt) )
	pblh_tim(1:pblh_nt) = pval_aa(1:pblh_nt,1) 
	pblh_val(1:pblh_nt) = pval_aa(1:pblh_nt,2) 
	    
	if (idiagbb > 0) then
	write(*,'(/a,i5)'  ) 'pblh profile', pblh_nt
	write(*,'(10f10.1)') pblh_tim(1:pblh_nt)
	write(*,'(10f10.2)') pblh_val(1:pblh_nt)
	end if ! (idiagbb > 0)


!   initial gas mixing ratios
	nparam = ngas_max
	pname_aa(1:ngas_max) = species(1:ngas_max)
	ptype_aa(:) = 2
	imustfind(1:nparam) = 0
	call pmcmos_readfile_aa( 0, iok, &
	    gas_init_fname, lun_pmcmos_in1, &
	    nparam, itmpa, ntmax_pmcmos, &
	    imustfind, ifound, pname_aa, ptype_aa, ptxt_aa, pval_aa )
	if (itmpa > 1) then
	    write(*,*) '*** pmcmos_read_inputs error for gas_init'
	    write(*,*) '    bad nt value = ', itmpa
	    stop
	end if
	gas_init_val(1:ngas_max) = pval_aa(1,1:ngas_max) 
	    
	if (idiagbb > 0) then
	write(*,'(/a,i5)'  ) 'gas_init found'
	itmpa = 0
	do i = 1, ngas_max
	    if (ifound(i) > 0) then
		itmpa = itmpa + 1
		itmpvecaa(itmpa) = i
	    end if
	end do
	write(*,'(4(2x,a,1p,e12.4))') (species(itmpvecaa(i))(1:12), &
		gas_init_val(itmpvecaa(i)), i=1,itmpa)

	write(*,'(a,i5)'  ) 'gas_init not found'
	itmpa = 0
	do i = 1, ngas_max
	    if (ifound(i) <= 0) then
		itmpa = itmpa + 1
		itmpvecaa(itmpa) = i
	    end if
	end do
	write(*,'(4(2x,a,1p,e12.4))') (pname_aa(itmpvecaa(i))(1:12), &
		gas_init_val(itmpvecaa(i)), i=1,itmpa)
	end if ! (idiagbb > 0)


!   gas background mixing ratios
	nparam = ngas_max + 2
	pname_aa(1:ngas_max) = species(1:ngas_max)
	pname_aa(nparam-1) = 'time'
	pname_aa(nparam  ) = 'rate'
	ptype_aa(:) = 2
	imustfind(1:nparam) = 0
	imustfind(nparam-1:nparam) = 1
	call pmcmos_readfile_aa( 0, iok, &
	    gas_back_fname, lun_pmcmos_in1, &
	    nparam, gas_back_nt, ntmax_pmcmos, &
	    imustfind, ifound, pname_aa, ptype_aa, ptxt_aa, pval_aa )
	allocate( gas_back_tim(1:gas_back_nt) )
	allocate( gas_back_dil(1:gas_back_nt) )
	allocate( gas_back_val(1:gas_back_nt,1:ngas_max) )
	gas_back_val(1:gas_back_nt,1:ngas_max) = pval_aa(1:gas_back_nt,1:ngas_max) 
	gas_back_tim(1:gas_back_nt) = pval_aa(1:gas_back_nt,nparam-1) 
	gas_back_dil(1:gas_back_nt) = pval_aa(1:gas_back_nt,nparam  ) 
	    
	if (idiagbb > 0) then
	write(*,'(/a,i5)' ) 'gas_back found', gas_back_nt
	write(*,'(a,1p,10e10.2)') 'time        ', gas_back_tim(1:gas_back_nt)
	write(*,'(a,1p,10e10.2)') 'dilu        ', gas_back_dil(1:gas_back_nt)
	do i = 1, ngas_max
	    if (ifound(i) <= 0) cycle
	    write(*,'(a,1p,10e10.2)') pname_aa(i)(1:12), gas_back_val(1:gas_back_nt,i)
	end do

	write(*,'(a,i5)'  ) 'gas_init not found'
	do i = 1, ngas_max
	    if (ifound(i) > 0) cycle
	    write(*,'(a,1p,10e10.2)') pname_aa(i)(1:12), gas_back_val(1:gas_back_nt,i)
	end do
	end if ! (idiagbb > 0)


!   gas emissions 
	nparam = ngas_max + 2
	pname_aa(1:ngas_max) = species(1:ngas_max)
	pname_aa(nparam-1) = 'time'
	pname_aa(nparam  ) = 'rate'
	ptype_aa(:) = 2
	imustfind(1:nparam) = 0
	imustfind(nparam-1:nparam) = 1
	call pmcmos_readfile_aa( 0, iok, &
	    gas_emit_fname, lun_pmcmos_in1, &
	    nparam, gas_emit_nt, ntmax_pmcmos, &
	    imustfind, ifound, pname_aa, ptype_aa, ptxt_aa, pval_aa )
	allocate( gas_emit_tim(1:gas_emit_nt) )
	allocate( gas_emit_fac(1:gas_emit_nt) )
	allocate( gas_emit_val(1:gas_emit_nt,1:ngas_max) )
	gas_emit_val(1:gas_emit_nt,1:ngas_max) = pval_aa(1:gas_emit_nt,1:ngas_max) 
	gas_emit_tim(1:gas_emit_nt) = pval_aa(1:gas_emit_nt,nparam-1) 
	gas_emit_fac(1:gas_emit_nt) = pval_aa(1:gas_emit_nt,nparam  ) 
	    
	if (idiagbb > 0) then
	j = min( 10, gas_emit_nt )
	write(*,'(/a,i5)' ) 'gas_emit found', gas_emit_nt
	write(*,'(a,1p,10e10.2)') 'time        ', gas_emit_tim(1:j)
	if (j < gas_emit_nt) &
	write(*,'(12x,1p,10e10.2)') gas_emit_tim(j+1:gas_emit_nt)
	write(*,'(a,1p,10e10.2)') 'factor      ', gas_emit_fac(1:j)
	if (j < gas_emit_nt) &
	write(*,'(12x,1p,10e10.2)') gas_emit_fac(j+1:gas_emit_nt)

	do i = 1, ngas_max
	    if (ifound(i) <= 0) cycle
	    write(*,'(a,1p,10e10.2)') pname_aa(i)(1:12), gas_emit_val(1:j,i)
	    if (j < gas_emit_nt) &
	    write(*,'((12x,1p,10e10.2))') gas_emit_val(j+1:gas_emit_nt,i)
	end do

	write(*,'(a,i5)'  ) 'gas_init not found'
	do i = 1, ngas_max
	    if (ifound(i) > 0) cycle
	    write(*,'(a,1p,10e10.2)') pname_aa(i)(1:12), gas_emit_val(1:j,i)
	    if (j < gas_emit_nt) &
	    write(*,'((12x,1p,10e10.2))') gas_emit_val(j+1:gas_emit_nt,i)
	end do
	end if ! (idiagbb > 0)


!   aerosol background distributions
	nparam = 3
	pname_aa(1) = 'time'
	pname_aa(2) = 'rate'
	pname_aa(3) = 'dist'
	ptype_aa(:) = 2
	ptype_aa(3) = 3
	imustfind(:) = 1
	call pmcmos_readfile_aa( 0, iok, &
	    aer_back_fname, lun_pmcmos_in1, &
	    nparam, aer_back_nt, ntmax_pmcmos, &
	    imustfind, ifound, pname_aa, ptype_aa, ptxt_aa, pval_aa )
	allocate( aer_back_tim(1:aer_back_nt) )
	allocate( aer_back_dil(1:aer_back_nt) )
	allocate( aer_back_dist(1:aer_back_nt) )
	aer_back_tim(1:aer_back_nt) = pval_aa(1:aer_back_nt,1) 
	aer_back_dil(1:aer_back_nt) = pval_aa(1:aer_back_nt,2) 
	aer_back_dist(1:aer_back_nt)%dist_fname = ptxt_aa(1:aer_back_nt,3) 
	    
	if (idiagcc > 0) then
	write(*,'(/a,i5)' ) 'aer_back found', aer_back_nt
	write(*,'(a,1p,10e10.2)') 'time        ', aer_back_tim(1:aer_back_nt)
	write(*,'(a,1p,10e10.2)') 'dilu        ', aer_back_dil(1:aer_back_nt)
!	write(*,'(a,1p,10(1x,a)') 'dist_fname  ', ( trim(aer_back_dist(i)%dist_fname), i=1,aer_back_nt )
	
	end if ! (idiagbb > 0)


!   aerosol emissions
	nparam = 3
	pname_aa(1) = 'time'
	pname_aa(2) = 'rate'
	pname_aa(3) = 'dist'
	ptype_aa(:) = 2
	ptype_aa(3) = 3
	imustfind(:) = 1
	call pmcmos_readfile_aa( 0, iok, &
	    aer_emit_fname, lun_pmcmos_in1, &
	    nparam, aer_emit_nt, ntmax_pmcmos, &
	    imustfind, ifound, pname_aa, ptype_aa, ptxt_aa, pval_aa )
	allocate( aer_emit_tim(1:aer_emit_nt) )
	allocate( aer_emit_fac(1:aer_emit_nt) )
	allocate( aer_emit_dist(1:aer_emit_nt) )
	aer_emit_tim(1:aer_emit_nt) = pval_aa(1:aer_emit_nt,1) 
	aer_emit_fac(1:aer_emit_nt) = pval_aa(1:aer_emit_nt,2) 
	aer_emit_dist(1:aer_emit_nt)%dist_fname = ptxt_aa(1:aer_emit_nt,3) 
	    
	if (idiagcc > 0) then
	write(*,'(/a,i5)' ) 'aer_emit found', aer_emit_nt
	write(*,'(a,1p,10e10.2)') 'time        ', aer_emit_tim(1:aer_emit_nt)
	write(*,'(a,1p,10e10.2)') 'rate        ', aer_emit_fac(1:aer_emit_nt)
!	write(*,'(a,1p,10(1x,a)') 'dist_fname  ', ( trim(aer_emit_dist(i)%dist_fname), i=1,aer_emit_nt )
	end if ! (idiagbb > 0)


	call pmcmos_read_aer_distribs


	return
	end subroutine pmcmos_read_inputs


!-----------------------------------------------------------------------
	subroutine pmcmos_readfile_aa( idiagaa, iok, &
	    fname, lun, &
	    np, nt, ntmax, &
	    imustfind, ifound, pname, ptype, ptxt, pval )
!
!   reads data from a temp_profile or pblh_profile or ... file
!   the data are equal length 1d arrays of values for several parameters
!

!   subr arguments (none)
	integer, intent(in)  :: idiagaa, lun, np, ntmax
	integer, intent(in)  :: imustfind(np), ptype(np)
	integer, intent(out) :: ifound(np), iok, nt
	real(r8), intent(out) :: pval(ntmax,np)
	character(len=*), intent(in) :: fname, pname(np)
	character(len=80), intent(out) :: ptxt(ntmax,np)

!   local variables
	integer, parameter :: maxch = 128, maxlnlen=2000
	integer :: i, iline, ip, itmpa
	integer :: ipostn(ntmax_pmcmos), itype(ntmax_pmcmos), ival(ntmax_pmcmos)
	integer :: j
	integer :: nchqq(ntmax_pmcmos)
	integer :: nfound, ntok

	real(r8) :: xval(ntmax_pmcmos)
	real(r8) :: tmpa

	character(len=maxlnlen) :: ast
	character(len=maxch) :: strqq(ntmax_pmcmos)
        character(len=maxch) :: txtaa, txtee

 !!!!!!!!
        iok = 0
        
	open( unit=lun, file=fname, status='old', iostat=itmpa )
	if (itmpa /= 0) then
	    txtee = 'open error'
	    goto 7200
	end if

	ifound(1:np) = 0
	pval(1:ntmax,1:np) = 0.0
	ptxt(1:ntmax,1:np) = ' '

	nfound = 0

loop_iline: &
	do iline = 1, 999888777
	if (idiagaa > 0) write(*,'(i5,2x,a)') iline, trim(fname)

	read(lun,'(a)',end=7100) ast
	if (idiagaa > 0) write(*,'(3a)') 'ast = |', trim(ast), '|'
	call parselnb( ast, ntok, itype, ival, xval,   &
		nchqq, strqq, ipostn, maxlnlen, maxch, ntmax_pmcmos )
	if (idiagaa >= 10) write(*,'((i4,3a,2x,i3,1p,e10.2))') (j, ' strqq=|', &
	    trim(strqq(j)), '|', itype(j), xval(j), j=1,ntok)

	if (ntok <= 0) cycle loop_iline
	if (strqq(1)(1:1) == '#') cycle loop_iline

	ip = -1
	do i = 1, np
	    if (strqq(1) == pname(i)) then
		ip = i
		exit
	    end if
	end do
	if (idiagaa > 0) write(*,'(i5,2x,a)') ip, ' = ip'
	if (ip <= 0) cycle loop_iline
	if (ifound(ip) > 0) cycle loop_iline

! check nt = number of times
	if (nfound == 0) then
	    nt = ntok - 1
	    if ((nt < 1) .or. (nt > ntmax_pmcmos)) then
		txtee = 'bad nt'
		goto 7200
	    end if
	else
	    if (ntok-1 /= nt) then
		txtee = 'bad nt'
		goto 7200
	    end if
	end if

! load parameter values into pval array
	do j = 1, nt
	    tmpa = -9.99e35
	    if (itype(j+1) == 1) then
		tmpa = ival(j+1)
	    else if (itype(j+1) == 2) then
		tmpa = xval(j+1)
	    else if (ptype(ip) <= 2) then
		txtee = 'bad token type'
		goto 7200
	    end if
	    if (ptype(ip) <= 2) pval(j,ip) = tmpa
	    ptxt(j,ip) = strqq(j+1)
	end do
	ifound(ip) = 1
	nfound = nfound + 1

	end do loop_iline

7100	close( unit=lun )
	if (ifound(ip) <= 0) then 
	    txtee = 'did not find anything (nfound == 0)'
	    goto 7200
	end if
	do ip = 1, np
	    if (imustfind(ip) <= 0) cycle
	    if (ifound(ip) <= 0) then 
		txtee = 'did not find ' // pname(ip)
		goto 7200
	    end if
	end do

	return


7200	write(*,'(2a)') '*** pmcmos_readfile_aa error = ', trim(txtee)
	write(*,'(2a)') 'filename = ', trim(fname)
	stop

	end subroutine pmcmos_readfile_aa


!-----------------------------------------------------------------------
	subroutine parselnb( ast, ntok, itype, ival, xval,   &
		nchqq, strqq, ipostn, maxlnlen, maxch, maxp )
!
!   parses a text string, using blanks as delims.
!
!	ast - string that is parsed
!
!	maxlnlen - length of ast
!
!	ntok - no. of tokens found.
!		on end-of-file, ntok = -1
!
!	for each token, 5 params are returned:
!
!    	itype - type code
!		1 = integer
!		2 - real
!		3 - string
!	
!	ival - integer value for a numeric token
!	xval - f.p. value for a numeric token
!
!	nchqq - no. of chars in the token
!	strqq - the token as a string
!
!	ipostn - the position of the start of string in the line
!
!    calls:
!	isdigit, issign, isperiod
!
!   sub arguments

	character*(maxlnlen) :: ast
!   ast holds input line,  bst holds token for numeric decoding

	integer :: maxlnlen, maxch, maxp, ntok
	integer :: itype(maxp), ival(maxp), nchqq(maxp), ipostn(maxp)
	real(r8) :: xval(maxp)
	character(len=maxch) :: strqq(maxp)

!   local variables
	integer :: j, ja, jz, lenj, nr
	character(len=30) :: bst
	character(len=1) :: ach1,ach2,ach3
	character(len=1) :: tabch

	tabch = char(9)
9000	format( a )


!... count non-blank, non-tab chars
!... also convert tabs to blanks
	do j = 1, maxlnlen
	    if (ichar(ast(j:j)) .le. 25) ast(j:j) = ' '
	end do

	nr = maxlnlen
1100	if (ast(nr:nr) .eq. ' ') then
	    nr = nr - 1
	    if (nr .gt. 1) goto 1100
	end if

!... locate and identify tokens
	ntok = 0
	j = 1

2100	if (j .le. nr) then

!... ignore blanks
2200	    if (ast(j:j) .eq. ' ') then
		j = j+1
		if (j .gt. nr) go to 4900
		goto 2200
	    end if

!... now count non-blanks
	    ja = j
2300	    if ((ast(j:j).ne.' ') .and. (j.le.nr)) then
		j = j+1
		goto 2300
	    end if

	    ntok = ntok + 1
	    if (ntok .gt. maxp) go to 8000
	    lenj = j-ja
	    nchqq(ntok) = lenj
	    ipostn(ntok) = ja

	    lenj = min0( maxch, j - ja)
	    jz = j - 1
	    strqq(ntok) = ' '
	    strqq(ntok)(1:lenj) = ast(ja:jz)

	    lenj = min0( 30, j-ja )
	    bst = ' '
	    bst(31-lenj : 30) = ast(ja:jz)
!	    write(*,'(a,i5,3a)') 'ntok, bst', ntok, '  >', bst, '<'

!    bypass the following
	goto 3400

!... try to identify numeric - must start as
!...	digit - integer or real
!...	sign digit - integer or real
!...	period digit - real
!...	sign period digit - real
	    ach1 = ast(ja:ja)
	    ach2 = ast(ja+1:ja+1)
	    ach3 = ast(ja+2:ja+2)
	    if ( isdigit(ach1) .or.   &
      		(issign(ach1).and.isdigit(ach2)) ) then
		go to 3400
	    else if ( (isperiod(ach1).and.isdigit(ach2)) .or.   &
      		    (issign(ach1).and.isperiod(ach2).and.   &
      		    isdigit(ach3)) ) then
		go to 3500
	    else
		go to 3600
	    end if

!... try for integer - non-integer causes read error
!...
3400	    read( bst, 9100, err=3500 ) ival(ntok)
9100	    format( i30 )
	    xval(ntok) = float( ival(ntok) )
	    itype(ntok) = 1
	    go to 3900

!... try for real - non-real causes read error
!...
3500	    read( bst, 9200, err=3600 ) xval(ntok)
9200	    format( e30.5 )
	    read( bst,    *, err=3600 ) xval(ntok)
!	    ival(ntok) = ifix( xval(ntok) )
	    ival(ntok) = 0
	    itype(ntok) = 2
	    go to 3900

!... must be symbol - only 4 chars allowed, and check for 'show' and 'all'
!...
3600	    itype(ntok) = 3
	    ival(ntok) = 0
	    xval(ntok) = 0.

3900	    goto 2100
	end if

4900	continue

8000	return
	end subroutine parselnb


!-----------------------------------------------------------------------
	logical function isdigit( ach )
	character*1 ach
	isdigit = (ach .ge. '0') .and. (ach .le. '9')
	return
	end function isdigit


!-----------------------------------------------------------------------
	logical function issign( ach )
	character*1 ach
	issign = (ach .eq. '+') .or. (ach .eq. '-')
	return
	end function issign


!-----------------------------------------------------------------------
	logical function isperiod( ach )
	character*1 ach
	isperiod = (ach .eq. '.')
	return
	end function isperiod


!-----------------------------------------------------------------------
	subroutine pmcmos_read_aer_distribs

	use module_data_mosaic_aero, only:  &
	    aer_name, naer

!   subr arguments (none)

!   local variables
	integer  :: i, idiagaa, idiagbb, idiagcc, iok, itmpa
	integer  :: j, jj
	integer  :: l

	type( aer_dist_t ) :: ad

	idiagbb = -1
	idiagcc = -1

!   aerosol initial distrib
	aer_init_dist%dist_fname = aer_init_fname
	call pmcmos_read_1_aer_distrib( 0, iok, &
	    lun_pmcmos_in1, aer_init_dist )

!   aerosol background distributions
	do j = 1, aer_back_nt
	    call pmcmos_read_1_aer_distrib( 0, iok, &
		lun_pmcmos_in1, aer_back_dist(j) )
	end do

!   aerosol emissions
	do j = 1, aer_emit_nt
	    call pmcmos_read_1_aer_distrib( 0, iok, &
		lun_pmcmos_in1, aer_emit_dist(j) )
	end do


!   ***do this to set the aer_name array
	call load_mosaic_parameters


!   read composition files and write diagnostics
	do jj = 1, 1 + aer_back_nt + aer_emit_nt
	    if (jj == 1) then
		ad = aer_init_dist
		idiagaa = -1
	    else if (jj <= 1+aer_back_nt) then
		j = jj - 1
		ad = aer_back_dist(j)
		idiagaa = -1
	    else
		j = jj - 1 - aer_back_nt
		ad = aer_emit_dist(j)
		idiagaa = -1
	    end if

	    call pmcmos_read_1_ad_comps( idiagaa, iok, &
		jj, ad )

	    write(*,'(/2a)') 'filename = ', trim(ad%dist_fname)
	    write(*,'(a,i5)') 'n_mode =', ad%n_mode
	    do i = 1, ad%n_mode
		ad%sg(i) = 10.0d0**ad%sg(i)
		ad%lnsg(i) = log( max( 1.0d0, ad%sg(i) ) )
		ad%dgnum_cm(i) = 200.0d0*ad%rgnum_m(i)
		write(*,'(4(a,3x))') trim(ad%mode_name(i)), trim(ad%mode_type(i)), &
			trim(ad%frac_type(i)), trim(ad%frac_fname(i)) 
		write(*,'(i5,1p,3e12.4)') i, ad%numden(i), ad%dgnum_cm(i), ad%sg(i)
		do l = 1, naer
		    if (ad%massfrac(l,i) > 0.0) write(*,'(f12.6,2x,a)') &
			ad%massfrac(l,i), aer_name(l)
		end do
	    end do

	    if (jj == 1) then
		aer_init_dist = ad
	    else if (jj <= 1+aer_back_nt) then
		j = jj - 1
		aer_back_dist(j) = ad
	    else
		j = jj - 1 - aer_back_nt
		aer_emit_dist(j) = ad
	    end if

	end do
	write(*,'(/)')


	return
	end subroutine pmcmos_read_aer_distribs


!-----------------------------------------------------------------------
	subroutine pmcmos_read_1_aer_distrib( idiagaa, iok, &
	    lun, ad )
!
!   reads data from a temp_profile or pblh_profile or ... file
!   the data are equal length 1d arrays of values for several parameters
!

!   subr arguments (none)
	integer, intent(in)  :: idiagaa, lun
	integer, intent(out) :: iok
	type( aer_dist_t ), intent(inout) :: ad

!   local variables
	integer, parameter :: maxch = 128, maxlnlen=2000
	integer :: i, iline, imode, imodeb, itmpa
	integer :: ipostn(ntmax_pmcmos), itype(ntmax_pmcmos), ival(ntmax_pmcmos)
	integer :: j, jp
	integer :: nchqq(ntmax_pmcmos)
	integer :: ntok

	real(r8) :: xval(ntmax_pmcmos)
	real(r8) :: tmpa

	character(len=maxlnlen) :: ast
	character(len=maxch) :: strqq(ntmax_pmcmos)
	character(len=maxch) :: str_expected(6)
	character(len=128) :: txtaa, txtee

        iok = 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
	open( unit=lun, file=ad%dist_fname, status='old', iostat=itmpa )
	if (itmpa /= 0) then
	    txtee = 'open error'
	    goto 7200
	end if

	str_expected(1) = 'frac_type   '
	str_expected(2) = 'frac        '
	str_expected(3) = 'num_den     '
	str_expected(4) = 'mode_type   '
	str_expected(5) = 'mean_radius '
	str_expected(6) = 'log_std_dev '

	ad%n_mode = 0
	imodeb = 0
	iline = 0
	if (idiagaa > 0) write(*,'(/a)') &
		'** pmcmos_read_1_aer_distrib diagnostics'

loop_imode: &
	do imode = 1, 999888777
	jp = -1

	strqq(1) = ' '
	do while (strqq(1) /= 'mode_name')
	    iline = iline + 1
	    if (idiagaa > 0) write(*,'(i5,2x,a)') iline, trim(ad%dist_fname)
	    read(lun,'(a)',end=5100) ast
	    if (idiagaa > 0) write(*,'(3a)') 'ast = |', trim(ast), '|'
	    call parselnb( ast, ntok, itype, ival, xval,   &
		nchqq, strqq, ipostn, maxlnlen, maxch, ntmax_pmcmos )
	    if (idiagaa >= 10) write(*,'((i4,3a,2x,i3,1p,e10.2))') (j, ' strqq=|', &
		trim(strqq(j)), '|', itype(j), xval(j), j=1,ntok)
	end do
	if (imode > ad_maxmode) then
	    txtee = 'too many modes'
	    goto 7200
	end if
	ad%n_mode = imode
	ad%mode_name(imode) = strqq(2)

	do jp = 1, 6
	    iline = iline + 1
	    if (idiagaa > 0) write(*,'(i5,2x,a)') iline, trim(ad%dist_fname)
	    read(lun,'(a)',end=7100) ast
	    if (idiagaa > 0) write(*,'(3a)') 'ast = |', trim(ast), '|'
	    call parselnb( ast, ntok, itype, ival, xval,   &
		nchqq, strqq, ipostn, maxlnlen, maxch, ntmax_pmcmos )
	    if (idiagaa >= 10) write(*,'((i4,3a,2x,i3,1p,e10.2))') (j, ' strqq=|', &
		trim(strqq(j)), '|', itype(j), xval(j), j=1,ntok)

	    if (strqq(1) /= str_expected(jp)) then
		txtee = 'expecting strqq(1) = ' // str_expected(jp)
	        write(*,'(3a)') ' strqq=|', trim(strqq(1)), '|'
	        write(*,'(3a)') ' expec=|', trim(str_expected(jp)), '|'
		goto 7200
	    else if (ntok < 2) then
		txtee = 'expecting at least 2 tokens' ; goto 7200
	    end if

	    tmpa = -9.99e35
	    if (itype(2) == 1) then
		tmpa = ival(2)
	    else if (itype(2) == 2) then
		tmpa = xval(2)
	    end if
	
	    if      (jp == 1) then
		ad%frac_type(imode) = strqq(2)
	    else if (jp == 2) then
		ad%frac_fname(imode) = strqq(2)
	    else if (jp == 3) then
		ad%numden(imode) = tmpa
	    else if (jp == 4) then
		ad%mode_type(imode) = strqq(2)
	    else if (jp == 5) then
		ad%rgnum_m(imode) = tmpa
	    else
		ad%sg(imode) = tmpa
	    end if

	end do ! jp
	jp = -2

	imodeb = imode

	end do loop_imode

5100	close( unit=lun )

	return


7100	txtee = 'eof before done'

7200	write(*,'(2a)') '*** pmcmos_read_1_aer_distrib error = ', trim(txtee)
	write(*,'(2a)') 'filename = ', trim(ad%dist_fname)
	write(*,'(a,4i8)') 'imode, imodeb, iline, jp =', imode, imodeb, iline, jp
	stop

	end subroutine pmcmos_read_1_aer_distrib


!-----------------------------------------------------------------------
	subroutine pmcmos_read_1_ad_comps( idiagaa, iok, &
	    jj, ad )
!
!   reads data from a temp_profile or pblh_profile or ... file
!   the data are equal length 1d arrays of values for several parameters
!
	use module_data_mosaic_aero, only:  &
	    aer_name, naer

!   subr arguments (none)
	integer, intent(in)  :: idiagaa, jj
	integer, intent(out) :: iok
	type( aer_dist_t ), intent(inout) :: ad

!   local variables
	integer, parameter :: nparam_max = 2+naer
	integer, parameter :: ntmax_adcomp = 2
	integer  :: imode, itmpa
	integer  :: imustfind(nparam_max), ifound(nparam_max)
	integer  :: jok
	integer  :: nparam
	integer  :: ptype_aa(nparam_max)
	real(r8) :: pval_aa(ntmax_adcomp,nparam_max)
	real(r8) :: tmpa
	character(len=20) :: pname_aa(nparam_max)
	character(len=80) :: ptxt_aa(ntmax_adcomp,nparam_max)

        iok = 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	do imode = 1, ad%n_mode

	nparam = naer
	pname_aa(1:naer) = aer_name(1:naer)
	ptype_aa(:) = 2
	imustfind(1:nparam) = 0
	call pmcmos_readfile_aa( idiagaa, jok, &
	    ad%frac_fname(imode), lun_pmcmos_in1, &
	    nparam, itmpa, ntmax_adcomp, &
	    imustfind, ifound, pname_aa, ptype_aa, ptxt_aa, pval_aa )
	if (itmpa > 1) then
	    write(*,*) '*** pmcmos_read_1_ad_comps error'
	    write(*,*) '    bad nt value = ', itmpa
	    write(*,*) '    jj, imode    = ', jj, imode
	    stop
	end if
	tmpa = sum( pval_aa(1,1:naer) )
	if (tmpa < 0.1) then
	    write(*,*) '*** pmcmos_read_1_ad_comps error'
	    write(*,*) '    massfrac sum < 0.1', tmpa
	    write(*,*) '    jj, imode    = ', jj, imode
	    stop
	end if
	ad%massfrac(1:naer,imode) = pval_aa(1,1:naer)/tmpa
	    
	end do ! imode


	return
	end subroutine pmcmos_read_1_ad_comps


!-----------------------------------------------------------------------


	end module module_pmcmos_init

