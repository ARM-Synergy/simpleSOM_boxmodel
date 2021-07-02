c Program zap analyzes large reaction mechanisms and creates
c ordinary differential equations for the dynamic species.
c 
c Author: Rahul A. Zaveri
c Date: Feb 18, 1998
c Email bugs, comments, and suggestions to zaveri@vtaix.cc.vt.edu
c
c June 19, 1998: fixed bug in analyze_reactions - dimension dumline(maxlines)
c
c----------------------------------------------------------------------- 
c summary of usage.
c
c zap requires two input files containing:
c 1. list of chemical species, and 
c 2. list of reaction mechanism
c
c example:
c
c file 1: species.list (user is free to choose the name)
c 
c
c file 2: mechanism.list (user is free to choose the name)
c
c
c
c
c----------------------------------------------------------------------- 
      program zap
      parameter(maxspecies=100, maxfixed=20, maxelements=10,
     &          maxrxn=200, maxreact=5, maxprod=40, maxterms=500,
     &          maxextralines=100)
c
      parameter (maxlines = 3*maxrxn    , maxlinech = 80, 
     &           maxwords = 4*maxprod+10, maxwordch = 40,
     &           linecutoff = 50)
c
c
      dimension nreacts(maxrxn),nprods(maxrxn),
     &          ryield(maxprod, maxrxn),
     &          nddlines(maxspecies),
     &          npplines(maxspecies),
     &          nelements(maxspecies),	! number of elements in ith species
     &          elemfactor(maxspecies,maxelements)
c
      character*(maxwordch) reactnt(maxreact, maxrxn), 
     &                      product(maxprod,  maxrxn), 
     &                      yield(maxprod, maxrxn),
     &                      rate(maxrxn),
     &                      species(maxspecies), 
     &                      fixed(maxfixed),
     &                      index(maxspecies),
     &                      comp(maxspecies,maxelements),
     &                      atom(maxelements),
     &                      sign(maxprod, maxrxn)
c
      character*(maxlinech) bigblank,
     &                      dd(maxspecies,maxextralines),
     &                      pp(maxspecies,maxextralines),
     &                      rconst(maxrxn),
     &                      rateexpr(maxrxn),
     &			    species_list, mechanism
c
      data bigblank(1:80)/' '/
c
c-----------------------------------------------------------------
c
c
      write(6,*)'Enter species file name'
      read(5,*)species_list
      write(6,*)'Enter mechanism file name'
      read(5,*)mechanism
c
      open(9 ,file=species_list)	! input: user-supplied
      open(10,file=mechanism)		! input: user-supplied
      open(21,file='rconst.f')		! output
      open(22,file='rate.f')		! output
      open(23,file='ode.f')		! output
c-----------------------------------------------------------------
c
      if(maxspecies.gt.maxrxn)then
        write(6,*)'Please set maxrxn at least equal to maxspecies'
        write(6,*)'in the parameter statements even if the'
        write(6,*)'actual no. of rxns is less than no. of species'
        write(6,*)'****stopping execution****'
        write(6,*)'   '
        stop
      endif
c
c------------------------------------------------------------------------
c
      call analyze_species(
     &		nspecies,	! number of species            : int : out
     &		index,		! species indices              : char: out
     & 		species,	! species names                : char: out
     &		comp,		! species composition          : char: out
     &          elemfactor,	! factor to comp	       : real: out 
     &		nelements,	! no. elements in each species : int : out 
     &		nfixed,		! number of fixed species      : int : out
     &		fixed,		! fixed species names	       : char: out
     &          natoms,		! no. of critical atoms	       : int : out
     &          atom)		! critical atoms for mass bal  : char: out
c
c------------------------------------------------------------------------
c
      call analyze_reactions(
     &		nrxn,		! number of reactions          : int : out
     &		nreacts,	! no. of reactants in rxn ir   : int : out
     &		nprods,		! no. of products in rxn ir    : int : out
     &		reactnt,	! reactant species in rxn ir   : char: out
     &          sign,		! sign of product in rxn ir    : char: out
     &		yield,		! yield of product in rxn ir   : char: out
     &		ryield,		! yield of product in rxn ir   : real: out
     &		product,	! product species in rxn ir    : char: out
     &		rate,		! rate array       	       : char: out
     & 		rconst)		! rate constant expression     : char: out
c
c------------------------------------------------------------------------
c
      call detect_phantoms(
     &		nspecies,	! number of species            : int : in
     & 		species,	! species names                : char: in
     &		nfixed,		! number of fixed species      : int : in
     &		fixed,		! fixed species names	       : char: in
     &		nrxn,		! number of reactions          : int : in
     &		nreacts,	! no. of reactants in rxn ir   : int : in
     &		nprods,		! no. of products in rxn ir    : int : in
     &		reactnt,	! reactant species in rxn ir   : char: in
     &		product)	! product species in rxn ir    : char: in
c
c------------------------------------------------------------------------
c
      call check_mass_balance(
     &		nspecies,	! number of species            : int : in
     & 		species,	! species names                : char: in
     &		nrxn,		! number of reactions          : int : in
     &		nreacts,	! no. of reactants in rxn ir   : int : in
     &		nprods,		! no. of products in rxn ir    : int : in
     &		reactnt,	! reactant species in rxn ir   : char: in
     &		product,	! product species in rxn ir    : char: in
     &          sign,		! yield sign                   : char: in
     &          ryield,         ! product yield		       : real: in
     &          nelements,	! no. of elements in a species : int : in
     &          comp,           ! elemental composition        : char: in
     &          elemfactor,     ! factor to comp               : real: in
     &          atom,		! critical atoms for mass bal  : char: in
     &          natoms)		! no. of atoms		       : int : in
c
c------------------------------------------------------------------------
c
      call create_ppdd(
     &		nspecies,	! number of species            : int : in
     &		index,		! species indices              : char: in
     & 		species,	! species names                : char: in
     &		nfixed,		! number of fixed species      : int : in
     &		fixed,		! fixed species names	       : char: in
     &		nrxn,		! number of reactions          : int : in
     &		nreacts,	! no. of reactants in rxn ir   : int : in
     &		nprods,		! no. of products in rxn ir    : int : in
     &		reactnt,	! reactant species in rxn ir   : char: in
     &          sign,		! sign of product in rxn ir    : char: in
     &		yield,		! yield of product in rxn ir   : char: in
     &		product,	! product species in rxn ir    : char: in
     &		rate,		! rate array       	       : char: in
     &          nddlines,	! number of destruction lines  : int : out
     &          npplines,	! number of production lines   : int : out
     &		rateexpr,	! rate expresion array	       : char: out
     &          pp,		! production for isp species   : char: out
     &		dd)		! destruction for isp species  : char: out
c
c------------------------------------------------------------------------
c
c create output FORTRAN codes
c
      write(6,*)'  '
      write(6,*)'Creating FORTRAN codes:'
c
c
      write(6,*)'(1) rconst.f'
      do irxn = 1, nrxn
        nchar = nbllen( rconst(irxn) )
        write(21,9001)rconst(irxn)
      enddo
c
c
      write(6,*)'(2) rate.f' 
      do irxn = 1, nrxn
        nchar = nbllen( rateexpr(irxn) )
        write(22,9001)rateexpr(irxn)
      enddo
c
c
      write(6,*)'(3) ode.f'
      do isp = 1, nspecies
        do il = 1,  npplines(isp)
          nchar = nbllen(pp(isp,il))
          write(23,9001)pp(isp,il)
        enddo

        do il = 1,  nddlines(isp)
          nchar = nbllen(dd(isp,il))
          write(23,9001)dd(isp,il)
        enddo
c
        nchar = 1
        write(23,9001)'c'
      enddo
c
      write(6,*)'  '
      write(6,*)'  '
c
c--------------------------------------------------------------------
c formatting statements
c
9000  format( a )
9001  format(A) 
      stop
      end
c
c
c
c
c==========================================================================
c
c reads and analyze user-supplied species.list
c
      subroutine analyze_species(
     &		nspecies,	! number of species            : int : out
     &		index,		! species indices              : char: out
     & 		species,	! species names                : char: out
     &		comp,		! species elemental composition: char: out
     &          elemfactor,	! factor to comp	       : real: out 
     &		nelements,	! no. elements in each species : int : out 
     &		nfixed,		! number of fixed species      : int : out
     &		fixed,		! fixed species names	       : char: out
     &          natoms,		! no. of critical atoms	       : int : out
     &          atom)		! critical atoms for mass bal  : char: out
c
c--------------------------------------------------------------------------
c
      parameter(maxspecies=100, maxfixed=20, maxelements=10,
     &          maxrxn=200, maxreact=5, maxprod=40, maxterms=500,
     &          maxextralines=100)
c
      parameter (maxlines = 3*maxrxn    , maxlinech = 80, 
     &           maxwords = 4*maxprod+10, maxwordch = 40,
     &           linecutoff = 50)
c
      dimension itype(maxwords), ival(maxwords), xval(maxwords), 
     &          nchqq(maxwords), ipostn(maxwords), nwords(maxlines)
c
      dimension nelements(maxspecies),	! no. of elements in the ith species
     &          elemfactor(maxspecies,maxelements)
c
      character*(maxwordch) strqq(maxwords),
     &                      species(maxspecies), 
     &                      index(maxspecies),
     &                      comp(maxspecies,maxelements),
     &                      fixed(maxfixed),
     &                      atom(maxelements)
c
      character*(maxlinech) dumline(maxspecies+maxfixed+10), 
     &                      dumch, bigblank, header
      data bigblank(1:80)/' '/
c-----------------------------------------------------------------------
c
      read(9,*)header, nspecies
c
c reading dynamic species...
      iline = 0
700   read(9,9000)dumch
      if(dumch(1:75).eq.bigblank(1:75) .or.
     &   dumch(1:1) .eq.'#')     goto 700 ! ignore blank/comment line
      if(dumch(1:5) .eq."*end*") goto 701 ! end of dynamic species
      iline = iline + 1
      dumline(iline) = dumch
      goto 700
c
c error check for dynamic species...
701   if(iline.eq.nspecies)then
        write(6,*)'      '
        write(6,*)'Number of dynamic species ok'
      else
        write(6,*)'      '
        write(6,*)'Error in dynamic species listing'
        write(6,*)'number of species expected by user = ', nspecies
        write(6,*)'number of species actually found   = ', iline
        write(6,*)'****stopping execution****'
        write(6,*)'   '
      stop
      endif
c
c analyze dynamic species list
      do 500 il = 1, nspecies
        call linefields(maxwords,dumline(il),nwords(il),itype,  
     &		        ival,xval,nchqq,strqq,ipostn)
c
          nelements(il) = 1 ! initialize
          index(il)     = strqq(2)
          species(il)   = strqq(3)

c set default value
          do ie = 1, maxelements
            elemfactor(il,ie) = 1.0	! default value
          enddo
c
c analyze elemental composition
      do 501 iw = 5, nwords(il)
c
        if(strqq(iw).eq.';')goto 500
        if(strqq(iw).eq.'+' .or.
     &     strqq(iw).eq.'IGNORE')goto 501
c
          if(xval(iw).gt.0.0)then
            elemfactor(il,nelements(il)) = xval(iw)
            goto 501
          else
            comp(il,nelements(il)) = strqq(iw)
            nelements(il) = nelements(il) + 1
          endif

501   continue
            nelements(il) = nelements(il) - 1

500   continue
c
c=========================================================================
c
c reading fixed species...
      read(9,*)header, nfixed
c
c reading fixed species...
      iline = 0
800   read(9,9000)dumch
      if(dumch(1:75).eq.bigblank(1:75) .or.
     &   dumch(1:1) .eq.'#')     goto 800 ! ignore blank/comment line
      if(dumch(1:5) .eq."*end*") goto 801 ! end of fixed species
      iline = iline + 1
      dumline(iline) = dumch
      goto 800

c error check for dynamic species...
801   if(iline.eq.nfixed)then
        write(6,*)'      '
        write(6,*)'Number of fixed species ok'
      else
        write(6,*)'      '
        write(6,*)'Error in fixed species listing'
        write(6,*)'number of species expected by user = ', nfixed
        write(6,*)'number of species actually found   = ', iline
        write(6,*)'****stopping execution****'
        write(6,*)'   '
      stop
      endif
c
      do 550 il = 1, nfixed
        call linefields(maxwords,dumline(il),nwords(il),itype,  
     &		        ival,xval,nchqq,strqq,ipostn)
         fixed(il) = strqq(2)
550   continue
c
c
c=========================================================================
c
c reading atoms...
      read(9,*)header, natoms

      do ia = 1, natoms
      read(9,*)idum, atom(ia)
      enddo
c
c
9000  format( a )
1000  return
      end






c-------------------------------------------------------------------------
c reads and analyzes mechanism
c
      subroutine analyze_reactions(
     &		nrxn,		! number of reactions          : int : out
     &		nreacts,	! no. of reactants in rxn ir   : int : out
     &		nprods,		! no. of products in rxn ir    : int : out
     &		reactnt,	! reactant species in rxn ir   : char: out
     &          sign,		! sign of product in rxn ir    : char: out
     &		yield,		! yield of product in rxn ir   : char: out
     &		ryield,		! yield of product in rxn ir   : real: out
     &		product,	! product species in rxn ir    : char: out
     &		rate,		! rate array       	       : char: out
     & 		rconst)		! rate constant expression     : char: out
c
c-------------------------------------------------------------------------
      parameter(maxspecies=100, maxfixed=20, maxelements=10,
     &          maxrxn=200, maxreact=5, maxprod=40, maxterms=500,
     &          maxextralines=100)
c
      parameter (maxlines = 3*maxrxn    , maxlinech = 80, 
     &           maxwords = 4*maxprod+10, maxwordch = 40,
     &           linecutoff = 50)
c
c
      dimension itype(maxwords), ival(maxwords), xval(maxwords), 
     &          nchqq(maxwords), ipostn(maxwords), nwords(maxlines)
c
      dimension valu(maxwords,maxlines),
     &          nchr(maxwords,maxlines),
     &          ityp(maxwords,maxlines),
     &          ipos(maxwords,maxlines)

c
      dimension nreacts(maxrxn),nprods(maxrxn),
     &          ryield(maxprod, maxrxn)

c
      character*(maxwordch) strqq(maxwords), 
     &                      word(maxwords,maxlines),
     &                      reactnt(maxreact, maxrxn), 
     &                      product(maxprod,  maxrxn), 
     &                      yield(maxprod, maxrxn),
     &                      rate(maxrxn),
     &                      sign(maxprod, maxrxn),
     &                      dumwd, irch
c
      character*(maxlinech) dumline(maxlines), 
     &                      dumch, bigblank, rconst(maxrxn)
c
      character*20 header
      data bigblank(1:80)/' '/
c------------------------------------------------------------------------
c
c read user-supplied mechanism.txt...
      read(10,*)header, nrxn
c  
      iline = 0
1000  read(10,9000)dumch
c
      if(dumch(1:75).eq.bigblank(1:75) .or.
     &   dumch(1:1) .eq.'#')     goto 1000	! ignore blank/comment line
      if(dumch(1:5) .eq."*end*") goto 1001	! end of file
c
      iline = iline + 1
      dumline(iline) = dumch
c      write(6,*)dumline(iline)	! Comment later
      goto 1000
c
1001  nlines = iline	! total no. of non-blank lines read
c
c
c create linewise words...
      do 500 il = 1, nlines
      call linefields(maxwords, dumline(il), nwords(il),itype,ival,  
     +		      xval, nchqq, strqq, ipostn)
c
        do iw = 1, nwords(il)
        word(iw,il) = strqq(iw)
        valu(iw,il) = xval(iw)
        nchr(iw,il) = nchqq(iw)
        ityp(iw,il) = itype(iw)
        ipos(iw,il) = ipostn(iw)
c
        enddo
c
500   continue
c
c
c------------------------------------------------------------------------
c create reactnt, sign, yield, and product arrays for each reaction...
c
c initializing nreacts and nprods
      do irxn=1,maxrxn
      nreacts(irxn) = 0
      nprods(irxn) = 1
      enddo

      ir = 0
c
      do 600 il = 1, nlines		! lines loop
c
      if(word(1,il)(1:1).eq.'{')then ! 1st word of new line indicates new rxn
        ir = ir + 1
        mp1 = 1		! reactant part
        mp2 = 0		! product part
        mp3 = 0		! rate constant part
        iwstart = 2	! 1st word is rate no., so start with 2nd word
      else
        iwstart = 1	! rxn cont'd on the next line, so start with 2nd word
      endif

      do 601 iw = iwstart, nwords(il)	! words loop
      dumwd = word(iw,il)

      if(dumwd.eq.'=')then
        mp2 = 1
        goto 601
      endif
c
      if(dumwd.eq.':')then
        mp3 = 1
        goto 601
      endif
c
      if(dumwd.eq.';')goto 601
c
c determine ipart
      ipart = 100 + 10*mp2 + mp3
      
      if(ipart.eq.100)goto 100		! part 1: reactants
      if(ipart.eq.110)goto 110		! part 2: products
      if(ipart.eq.111)goto 111		! part 3: rate constant
c
c
c---reactants-------------------------------------
100   if(dumwd.ne.'+')then
      nreacts(ir) = nreacts(ir) + 1
      reactnt(nreacts(ir),ir) = dumwd
      endif
      goto 601
c
c
c---products--------------------------------------
110   continue
c
c check for sign
      if(isign_state.eq.0)then
         if(dumwd.eq.'+' .or. dumwd.eq.'-')then
           sign(nprods(ir),ir) = dumwd(1:1)
           isign_state = 1
           goto 601
         elseif(nprods(ir) .eq. 1)then
           sign(nprods(ir),ir) = '+'
           isign_state = 1
         else
           write(6,*)'    '
           write(6,*)'Formatting error in reaction ', ir
           write(6,*)'* no space between sign and species or'
           write(6,*)'* sign missing for a product'
           write(6,*)'Hint: the error is at or near: ',dumwd
           write(6,*)'***stopping execution***'
           write(6,*)'    '
           stop
         endif
      endif
c
c check for yield
      if(iyield_state.eq.0)then
         if(valu(iw,il).gt.0.0)then
           ryield(nprods(ir),ir) = valu(iw,il)
           yield(nprods(ir),ir) = dumwd(1:nchr(iw,il))//'*'
           iyield_state= 1
           goto 601
         elseif(dumwd(1:1).eq.'#')then
           ryield(nprods(ir),ir) = 999.
           yield(nprods(ir),ir) = dumwd(2:nchr(iw,il))//'*'
           iyield_state= 1
           goto 601
         else
           yield(nprods(ir),ir) = ''
           ryield(nprods(ir),ir) = 1.0
           iyield_state= 1
         endif
      endif
c
c check for product species
      if(ityp(iw,il).eq.3  .and. 
     &   dumwd.ne.'+' .and.
     &   dumwd.ne.'-' .and.
     &   dumwd(1:1).ne.'#' )then
         product(nprods(ir),ir) = dumwd
         nprods(ir) = nprods(ir) + 1
         isign_state = 0
         iyield_state= 0
      endif
c
      goto 601
c
c
c---rate constant---------------------------------
c
111   write(irch,*)ir
      call linefields(maxwords, irch, ndum, itype, ival,  
     &		      xval, nchqq, strqq, ipostn)
      irch = strqq(1)(1:nchqq(1))
      rconst(ir) = '      rk('//irch(1:nbllen(irch))//') = '//
     &        dumline(il)(ipos(iw,il):ipos(nwords(il),il)-1)
      goto 600
c
601   continue
600   continue
c
c
      if(ir.eq.nrxn)then
        write(6,*)'      '
        write(6,*)'Number of reactions ok'
        do ir = 1, nrxn
          nprods(ir) =  nprods(ir) - 1
        enddo 
      else
        write(6,*)'Fatal error in reaction numbering'
        write(6,*)'number of rxn expected by user = ', nrxn
        write(6,*)'number of rxn actually found   = ', ir
        write(6,*)'****stopping execution****'
        write(6,*)'   '
        stop
      endif

c
c load rate array: r(1), r(2), etc.
        do ir = 1, nrxn
          write(irch,*)ir
          call linefields(maxwords, irch, ndum, itype, ival,  
     &		          xval, nchqq, strqq, ipostn)
          irch = strqq(1)(1:nchqq(1))
          rate(ir) = 'r('//irch(1:nbllen(irch))//')'
        enddo
c
c
9000  format( a )
      return
      end






c------------------------------------------------------------------------
c creates production and destruction arrays
c
      subroutine create_ppdd(
     &		nspecies,	! number of species            : int : in
     &		index,		! species indices              : char: in
     & 		species,	! species names                : char: in
     &		nfixed,		! number of fixed species      : int : in
     &		fixed,		! fixed species names	       : char: in
     &		nrxn,		! number of reactions          : int : in
     &		nreacts,	! no. of reactants in rxn ir   : int : in
     &		nprods,		! no. of products in rxn ir    : int : in
     &		reactnt,	! reactant species in rxn ir   : char: in
     &          sign,		! sign of product in rxn ir    : char: in
     &		yield,		! yield of product in rxn ir   : char: in
     &		product,	! product species in rxn ir    : char: in
     &		rate,		! rate array       	       : char: in
     &          nddlines,	! number of destruction lines  : int : out
     &          npplines,	! number of production lines   : int : out
     &		rateexpr,	! rate expresion array	       : char: out
     &          pp,		! production for isp species   : char: out
     &		dd)		! destruction for isp species  : char: out
c
c-------------------------------------------------------------------------
      parameter(maxspecies=100, maxfixed=20, maxelements=10,
     &          maxrxn=200, maxreact=5, maxprod=40, maxterms=500,
     &          maxextralines=100)
c
      parameter (maxlines = 3*maxrxn    , maxlinech = 80, 
     &           maxwords = 4*maxprod+10, maxwordch = 40,
     &           linecutoff = 50)
c
c
      dimension itype(maxwords), ival(maxwords), xval(maxwords), 
     &          nchqq(maxwords), ipostn(maxwords)
c
      dimension nreacts(maxrxn),nprods(maxrxn),
     &          nddlines(maxspecies),
     &          npplines(maxspecies),
     &          ndd(maxterms), npp(maxterms)
c
      character*(maxwordch) strqq(maxwords), 
     &                      reactnt(maxreact, maxrxn), 
     &                      product(maxprod,  maxrxn), 
     &                      yield(maxprod, maxrxn),
     &                      rate(maxrxn),
     &                      rk(maxrxn),
     &                      species(maxspecies), 
     &                      fixed(maxfixed),
     &                      index(maxspecies),
     &                      sign(maxprod, maxrxn),
     &                      irch,
     &                      dd_rate(maxspecies,maxterms),
     &                      pp_rate(maxspecies,maxterms),
     &                      pp_sign(maxspecies,maxterms),
     &                      pp_yield(maxspecies,maxterms)
c
      character*(maxlinech) bigblank,
     &                      dd(maxspecies,maxextralines),
     &                      pp(maxspecies,maxextralines),
     &                      rateexpr(maxrxn)
      data bigblank(1:80)/' '/
c
c------------------------------------------------------------------------
c
c create various production and destruction arrays
c
      do 200 isp  = 1, nspecies
      do 200 irxn = 1, nrxn

      do irct = 1, nreacts(irxn)
        if(reactnt(irct,irxn).eq.species(isp))then
          ndd(isp) = ndd(isp) + 1
          dd_rate(isp,ndd(isp)) = rate(irxn)
        endif
      enddo

      do iprd = 1, nprods(irxn)
        if(product(iprd,irxn).eq.species(isp))then
          npp(isp) = npp(isp) + 1
          pp_sign(isp ,npp(isp)) = sign(iprd, irxn)
          pp_yield(isp,npp(isp)) = yield(iprd, irxn)
          pp_rate(isp ,npp(isp)) = rate(irxn)
        endif
      enddo
          if(pp_sign(isp,1)(1:1).eq.'+')then
             pp_sign(isp,1)(1:1) = ''
          endif
c
200   continue
c
c build dd and pp character strings
      do 300 isp = 1, nspecies


c dd strings
           nlen = nbllen(index(isp))
           nddlines(isp) = 1		! no. of lines: 1st line
        if(ndd(isp).eq.0)then
           dd(isp,1)=
     &     '      d('//index(isp)(1:nlen)//')= 0.0'
        else
           dd(isp,1)=
     &     '      d('//index(isp)(1:nlen)//')= '//dd_rate(isp,1)
        endif
c
        do in = 2, ndd(isp)
         nnew = nbllen(dd_rate(isp,in))
         nddchar = nbllen( dd(isp,nddlines(isp)) )
         ntot = nddchar + nnew
 
         if(ntot.gt.linecutoff)then
           nddlines(isp)=nddlines(isp)+1
           indent = nlen + 4
           dd(isp,nddlines(isp))='     &'//bigblank(1:indent)//
     &                               '+'//dd_rate(isp,in)(1:nnew)
         else
           dd(isp,nddlines(isp))=dd(isp,nddlines(isp))(1:nddchar)//
     &                                         '+'             //
     &                                  dd_rate(isp,in)(1:nnew)
         endif
        enddo
c
c
c pp strings
           nlen = nbllen(index(isp))
           npplines(isp) = 1		! no. of lines: 1st line
        if(npp(isp).eq.0)then
           pp(isp,1)=
     &     '      p('//index(isp)(1:nlen)//')= 0.0'
        else
           pp(isp,1)=
     &     '      p('//index(isp)(1:nlen)//')= '//
c     &            pp_sign(isp,1)(1:nbllen( pp_sign(isp,1) ))//
     &            pp_yield(isp,1)(1: nbllen( pp_yield(isp,1) ))//
     &            pp_rate(isp,1)
        endif
c
        do in = 2, npp(isp)
         nnew = nbllen(pp_rate(isp,in))
         nppchar = nbllen( pp(isp,npplines(isp)) )
         ntot = nppchar + nnew
 
         if(ntot.gt.linecutoff)then
           npplines(isp)=npplines(isp)+1
           indent = nlen + 4
           pp(isp,npplines(isp))='     &'//bigblank(1:indent)//
     &                                 pp_sign(isp,in)(1:1)//
     &       pp_yield(isp,in)(1:nbllen( pp_yield(isp,in) ))//
     &                              pp_rate(isp,in)(1:nnew)
         else
           pp(isp,npplines(isp))=pp(isp,npplines(isp))(1:nppchar)//
     &                                   pp_sign(isp,in)(1:1)//
     &         pp_yield(isp,in)(1:nbllen( pp_yield(isp,in) ))//
     &                                pp_rate(isp,in)(1:nnew)
         endif
        enddo

300   continue
c
c
c create rateexpr array: r(1) = rk(1)*s(ia)*s(ib), etc.
      do 400 ir = 1, nrxn
          write(irch,*)ir
          call linefields(maxwords, irch, ndum, itype, ival,  
     &		          xval, nchqq, strqq, ipostn)
          irch = strqq(1)(1:nchqq(1))
          rk(ir) = 'rk('//irch(1:nbllen(irch))//')'

        rateexpr(ir) = '      '//
     &                   rate(ir)(1:nbllen( rate(ir) ))//
     &                   ' = '//
     &                   rk(ir)(1:nbllen( rk(ir) ))
c
      do 401 irct = 1, nreacts(ir)
c
c dynamic species
        do isp = 1, nspecies
          if(reactnt(irct,ir).eq.species(isp))then
            nrat = nbllen(rateexpr(ir))
            nind = nbllen(index(isp))
            rateexpr(ir) = rateexpr(ir)(1:nrat) // '*' //
     &                     's('//index(isp)(1:nind)//')'
          endif
        enddo
c
c fixed species
        do ifix = 1, nfixed
          if(reactnt(irct,ir).eq.fixed(ifix) .and.
     &       reactnt(irct,ir).ne.'hv')then
            nfix = nbllen(fixed(ifix))
            nrat = nbllen(rateexpr(ir))
            rateexpr(ir) = rateexpr(ir)(1:nrat) // '*' //
     &                     fixed(ifix)(1:nfix)
          endif
        enddo
401   continue
400   continue

      return
      end






c------------------------------------------------------------------------
c
      subroutine detect_phantoms(
     &		nspecies,	! number of species            : int : in
     & 		species,	! species names                : char: in
     &		nfixed,		! number of fixed species      : int : in
     &		fixed,		! fixed species names	       : char: in
     &		nrxn,		! number of reactions          : int : in
     &		nreacts,	! no. of reactants in rxn ir   : int : in
     &		nprods,		! no. of products in rxn ir    : int : in
     &		reactnt,	! reactant species in rxn ir   : char: in
     &		product)	! product species in rxn ir    : char: in
c
c------------------------------------------------------------------------
      parameter(maxspecies=100, maxfixed=20, maxelements=10,
     &          maxrxn=200, maxreact=5, maxprod=40, maxterms=500,
     &          maxextralines=100)
c
      parameter (maxlines = 3*maxrxn    , maxlinech = 80, 
     &           maxwords = 4*maxprod+10, maxwordch = 40,
     &           linecutoff = 50)
c
c
      dimension itype(maxwords), ival(maxwords), xval(maxwords), 
     &          nchqq(maxwords), ipostn(maxwords)
c
      dimension nreacts(maxrxn),nprods(maxrxn)
c
      character*(maxwordch) strqq(maxwords), 
     &                      reactnt(maxreact, maxrxn), 
     &                      product(maxprod,  maxrxn), 
     &                      species(maxspecies), 
     &                      fixed(maxfixed),
     &                      dumwd, irch
c
c-----------------------------------------------------------------
c
        write(6,*)'  '

        iphantom_tot = 0

      do 200 irxn = 1, nrxn
      do 200 irct = 1, nreacts(irxn)

        iphantom = 1

      do isp = 1, nspecies
        if(reactnt(irct,irxn).eq.species(isp))then
        iphantom = 0
        endif
      enddo

      do ifix = 1, nfixed
        if(reactnt(irct,irxn).eq.fixed(ifix))then
        iphantom = 0
        endif
      enddo


      if(iphantom.eq.1)then
        iphantom_tot = iphantom_tot + 1
        write(irch,*)irxn
        call linefields(maxwords, irch, ndum, itype, ival,  
     +		      xval, nchqq, strqq, ipostn)
        irch = strqq(1)(1:nchqq(1))
        nchar1 = nchqq(1)
        nchar2 = nbllen( reactnt(irct,irxn) )
        dumwd = 'reactant side'
        nchar3 = 13
        write(6,8001)irch, reactnt(irct,irxn), dumwd
      endif

200   continue
c
c
c
      do 300 irxn = 1, nrxn
      do 300 iprd = 1, nprods(irxn)

        iphantom = 1
      do isp = 1, nspecies
        if(product(iprd,irxn).eq.species(isp))then
        iphantom = 0
        endif
      enddo

      do ifix = 1, nfixed
        if(product(iprd,irxn).eq.fixed(ifix))then
        iphantom = 0
        endif
      enddo
c
      if(iphantom.eq.1)then
        iphantom_tot = iphantom_tot + 1
        write(irch,*)irxn
        call linefields(maxwords, irch, ndum, itype, ival,  
     +		      xval, nchqq, strqq, ipostn)
        irch = strqq(1)(1:nchqq(1))
        nchar1 = nchqq(1)
        nchar2 = nbllen( product(iprd,irxn) )
        dumwd = 'product side'
        nchar3 = 12
        write(6,8001)irch, product(iprd,irxn), dumwd
      endif

300   continue

        if(iphantom_tot.gt.0)then
      write(6,*)'Hint: there is a typo in species name or a space'
      write(6,*)'is missing between the sign and the species name'
      write(6,*)'***stopping execution***'
      write(6,*)'   '
      stop
        else
      write(6,*)'No phantom species were detected'
        endif
c
c8001  format(' Reaction ', A<nchar1>, ' contains phantom species < ',
c    &        A<nchar2> ' > on the ', A<nchar3>)

8001  format(' Reaction ', A, ' contains phantom species < ',
     &        A ' > on the ', A)

      return
      end






c------------------------------------------------------------------------
      subroutine check_mass_balance(
     &		nspecies,	! number of species            : int : in
     & 		species,	! species names                : char: in
     &		nrxn,		! number of reactions          : int : in
     &		nreacts,	! no. of reactants in rxn ir   : int : in
     &		nprods,		! no. of products in rxn ir    : int : in
     &		reactnt,	! reactant species in rxn ir   : char: in
     &		product,	! product species in rxn ir    : char: in
     &          sign,		! yield sign                   : char: in
     &          ryield,         ! product yield		       : real: in
     &          nelements,	! no. of elements in a species : int : in
     &          comp,           ! elemental composition        : char: in
     &          elemfactor,     ! factor to comp               : real: in
     &          atom,		! critical atoms for mass bal  : char: in
     &          natoms)         ! no. of atoms		       : int : in
c------------------------------------------------------------------------
      parameter(maxspecies=100, maxfixed=20, maxelements=10,
     &          maxrxn=200, maxreact=5, maxprod=40, maxterms=500,
     &          maxextralines=100)
c
      parameter (maxlines = 3*maxrxn    , maxlinech = 80, 
     &           maxwords = 4*maxprod+10, maxwordch = 40,
     &           linecutoff = 50)
c
      dimension nelements(maxspecies),
     &          elemfactor(maxspecies,maxelements),
     &          nreacts(maxrxn),nprods(maxrxn),
     &          ryield(maxprod, maxrxn)
c
      character*(maxwordch) species(maxspecies), 
     &                      comp(maxspecies,maxelements),
     &                      reactnt(maxreact, maxrxn), 
     &                      product(maxprod,  maxrxn), 
     &                      sign(maxprod, maxrxn),
     &                      atom(maxelements)
c
c------------------------------------------------------------------------
c
      iunbalance = 0 
c
      do 200 ir = 1, nrxn
      do 100 ia = 1, natoms
      elhs = 0.0
      erhs = 0.0
c
c left hand side (reactants)
        do 101 irct = 1, nreacts(ir)
        do isp= 1, nspecies
          if(reactnt(irct,ir).eq.species(isp))then
            do ie = 1, nelements(isp)
            if(comp(isp,ie).eq.atom(ia))then
            elhs = elhs + elemfactor(isp,ie)
            endif
            enddo
          endif 
        enddo
101     continue

c
c right hand side (products)
        do 102 iprd = 1, nprods(ir)
        do isp= 1, nspecies
          if(product(iprd,ir).eq.species(isp))then
            do ie = 1, nelements(isp)
            if(comp(isp,ie).eq.atom(ia))then
              if(sign(iprd,ir)(1:1).eq.'+')then
                erhs = erhs + elemfactor(isp,ie)*ryield(iprd,ir)
              elseif(sign(iprd,ir)(1:1).eq.'-')then
                erhs = erhs - elemfactor(isp,ie)*ryield(iprd,ir)
              else
                write(6,*)'Serious error in reaction ', ir
                write(6,*)'***stopping execution***'
                write(6,*)'  '
                stop
              endif
            endif
            enddo
          endif 
        enddo
102     continue


      if(abs(elhs-erhs).gt.1.e-5)then
        if(iunbalance.eq.0)then
        write(6,*)'  '
        endif
      write(6,7001)ir, atom(ia)(1:1)
      iunbalance = iunbalance + 1
      endif

100   continue
200   continue
c
7001  format(' Reaction',I6,' is unbalanced for element ',a5)
c7001  format(' Reaction',I6,' is unbalanced for element ',a<5>)
c
      if(iunbalance.gt.0)then
      write(6,*)'  '
      write(6,*)'One or more unbalanced reactions were detected'
      write(6,*)'***stopping execution***'
      write(6,*)'  '
      stop
      else
      write(6,*)'  '
      write(6,*)'Mass balance ok'
      endif
c
      return
      end






























c---------------------------------------------------------------------
c
	integer function nbllen( str )
c
c   returns the position of the last non-blank character in str
c
	character*(*) str

	j = len(str)

	if (j .gt. 0) then
1000	    if (str(j:j) .eq. ' ') then
		j = j - 1
		if (j .gt. 0) goto 1000
	    end if
	end if
	nbllen = max0( j, 0 )

	return
	end
c---------------------------------------------------------------------


c---------------------------------------------------------------------
c   subroutine linefields(ast,nfield,itype,ival,xval,nchqq,strqq,ipostn)
c---------------------------------------------------------------------
c   file linefields.for - from zvax on 7-nov-86
c   changes - maxwordch from 40 to 20
c
c   9-may-89 - cleaned up code ("end do", "do while", end-of-line comments)
c	deleted readline and bell subroutines
c
c   26-jun-90 - Sun version - deactivated conversion to uppercase
c
c------------------------
c   reads a line from unit lun and breaks it up into individual fields, 
c	using blanks as delims.
c
c	nfield - no. of fields found.  
c		on end-of-file, nfield = -1
c
c	for each field, 5 params are returned:
c
c    	itype - type code
c		1 = integer
c		2 - real
c		3 - string
c	
c	ival - integer value for a numeric field
c	xval - f.p. value for a numeric field
c
c	nchqq - no. of chars in the field
c	strqq - the field as a string
c
c	ipostn - the position of the start of string in the line 
c
c    calls:
c	isdigit, issign, isperiod
c---------------------------------------------------------------------

	subroutine linefields(maxwords, ast, nfield, itype, ival,  
     +		              xval, nchqq, strqq, ipostn )

	parameter (maxlinech = 80, maxwordch = 40)
	logical isdigit, issign, isperiod

c   ast holds input line,  bst holds field for numeric decoding
	character*(maxlinech) ast
	character*20 bst
	character*1 ach1,ach2,ach3

	dimension itype(*), ival(*), xval(*), nchqq(*), ipostn(*)
	character*(maxwordch) strqq(*)

c----------------------------------------------------------------------
	nfield = -1
	nr = maxlinech
1100	if (ast(nr:nr) .eq. ' ') then
	    nr = nr - 1
	    if (nr .gt. 1) goto 1100
	end if

c... convert tabs to spaces
	call tabfix( nr, maxlinech, ast )

c...convert to upper case (system library routine)
c	call str$upcase( ast, ast )

c... locate and identify fields
	nfield = 0
	j = 1

2100	if (j .le. nr) then

c... ignore blanks
2200	    if (ast(j:j) .eq. ' ') then
		j = j+1
		if (j .gt. 80) go to 4900
		goto 2200
	    end if

c... now count non-blanks
	    ja = j
2300	    if ((ast(j:j).ne.' ') .and. (j.le.nr)) then
		j = j+1
		goto 2300
	    end if

	    nfield = nfield + 1
	    if (nfield .gt. maxwords) go to 8000
	    len = j-ja
	    nchqq(nfield) = len
	    ipostn(nfield) = ja

	    len = min0( maxwordch, j - ja)
	    jz = j - 1
	    strqq(nfield) = ' '
	    strqq(nfield)(1:len) = ast(ja:jz)

	    len = min0( 20, j-ja )
	    bst = ' '
	    bst(21-len : 20) = ast(ja:jz)

c... try to identify numeric - must start as 
c...	digit - integer or real
c...	sign digit - integer or real
c...	period digit - real
c...	sign period digit - real
	    ach1 = ast(ja:ja)
	    ach2 = ast(ja+1:ja+1)
	    ach3 = ast(ja+2:ja+2)
	    if ( isdigit(ach1) .or. 
     +		(issign(ach1).and.isdigit(ach2)) ) then
		go to 3400
	    else if ( (isperiod(ach1).and.isdigit(ach2)) .or.
     +		    (issign(ach1).and.isperiod(ach2).and.
     +		    isdigit(ach3)) ) then
		go to 3500
	    else
		go to 3600
	    end if

c... try for integer - non-integer causes read error
c...
3400	    read( bst, 9100, err=3500 ) ival(nfield)
9100	    format( i20 )
	    xval(nfield) = float( ival(nfield) )
	    itype(nfield) = 1
	    go to 3900

c... try for real - non-real causes read error
c...
3500	    read( bst, 9200, err=3600 ) xval(nfield)
9200	    format( f20.5 )
c	    ival(nfield) = ifix( xval(nfield) )
	    ival(nfield) = 0
	    itype(nfield) = 2
	    go to 3900

c... must be symbol - only 4 chars allowed, and check for 'show' and 'all'
c...
3600	    itype(nfield) = 3
	    ival(nfield) = 0
	    xval(nfield) = 0.

3900	    goto 2100
	end if

4900	continue

8000	return
	end



	logical function isdigit( ach )
	character*1 ach
	isdigit = (ach .ge. '0') .and. (ach .le. '9')
	return
	end


	logical function issign( ach )
	character*1 ach
	issign = (ach .eq. '+') .or. (ach .eq. '-')
	return
	end


	logical function isperiod( ach )
	character*1 ach
	isperiod = (ach .eq. '.')
	return
	end


	subroutine tabfix( nr, maxlinech, ast )
c
c   changes tabs to blanks in a string, using standard vax tabbing
c
	character*(*) ast
	character*200 bst
	character*1 ch

	bst(1:nr) = ast(1:nr)

	nout = 0
	do 1900 ir = 1, nr
	    ch = bst(ir:ir)
	    if (nout .ge. maxlinech) goto 1910
	    nout = nout + 1
	    if (ch .ne. char(9)) then
		ast(nout:nout) = ch
	    else
		ntab = 9
1600		if (ntab .le. nout) then
		    ntab = ntab + 8
		    goto 1600
		end if
		do 1700 i = nout, ntab-1
		    if (i .gt. maxlinech) goto 1910
		    ast(i:i) = ' '
		    nout = i
1700		continue
	    end if
1900	continue

1910	nr = nout
	return
	end
c
