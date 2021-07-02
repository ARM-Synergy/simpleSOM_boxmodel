      subroutine init_aerosol
!
!   does initialization of variables in module_data_mosaic_asect
!
      use module_data_mosaic_main
      use module_data_mosaic_gas
      use module_data_mosaic_aero
      use module_data_mosaic_asect
      use module_pmcmos_aer, only:  pmcmos_init_aerosol
      
      implicit none
      
      integer :: ibin, iphase, isize, itype, it1, it2
      integer ::  l, ll, lunaa, noffset
      real(r8) :: tmpa, tmpb
      character(len=5) :: tmpch5

      LOGICAL COND1,COND2


      call load_mosaic_parameters

      
! check for valid msize_framework
      if ((msize_framework /= mmodal        ) .and. &
          (msize_framework /= munstructured ) .and. &
          (msize_framework /= msectional    )) then
         write(*,'(2a,2(1x,i10))') &
         '*** init_aerosol fatal error - ', &
         'bad msize_framework = ', msize_framework
         stop
      end if

! check misc flags
      if ((mhyst_method < 1) .or. (mhyst_method > 2)) then
         write(*,'(2a,2(1x,i10))') &
         '*** init_aerosol fatal error - ', &
         'bad mhyst_method = ', mhyst_method
         stop
      end if


! check for valid ntype_aer and nsize_aer(:)
      if ((ntype_aer < 1) .or. (ntype_aer > maxd_atype)) then
         write(*,'(2a,2(1x,i10))') &
         '*** init_aerosol fatal error - ', &
         'too big/small ntype_aer = ', ntype_aer
         stop
      end if

      do itype = 1, ntype_aer
         if ((nsize_aer(itype) < 1) .or. &
             (nsize_aer(itype) > maxd_asize)) then
            write(*,'(2a,2(1x,i10))') &
            '*** init_aerosol fatal error - ', &
            'too big/small nsize_aer = ', nsize_aer(itype), itype
            stop
         end if
         if (nsize_aer(itype) /= nsize_aer(1)) then
            write(*,'(2a,2(1x,i10))') &
            '*** init_aerosol fatal error - ', &
            'non-constant nsize_aer = ', nsize_aer(itype), itype
            stop
         end if
      end do

! phases -- only 1 for now
      nphase_aer = 1
      ai_phase = 1


! species pointers and properties
      if ((naer /= 24) .or. (naer_tot /= 27)) then
         write(*,'(2a,2(1x,i10))') &
         '*** init_aerosol fatal error - ', &
         'expecting naer=22 & naer_tot=27 but they = ', naer, naer_tot
         stop
      end if
      if (naer > maxd_acomp) then
         write(*,'(2a,2(1x,i10))') &
         '*** init_aerosol fatal error - ', &
         'naer > maxd_acomp = ', naer, maxd_acomp
         stop
      end if


      dens_mastercomp_aer(:) = 0.0
      mw_mastercomp_aer(:) = 0.0
      hygro_mastercomp_aer(:) = 0.0

      ntot_mastercomp_aer = naer
      do ll = 1, naer
         dens_mastercomp_aer(ll) = dens_aer_mac(ll)
         mw_mastercomp_aer(  ll) = mw_aer_mac(ll)
         name_mastercomp_aer(ll) = aer_name(ll)
      end do
      mastercompindx_so4_aer = iso4_a
      mastercompindx_no3_aer = ino3_a
      mastercompindx_cl_aer  = icl_a
      mastercompindx_msa_aer = imsa_a
      mastercompindx_co3_aer = ico3_a
      mastercompindx_nh4_aer = inh4_a
      mastercompindx_na_aer  = ina_a
      mastercompindx_ca_aer  = ica_a
      mastercompindx_oin_aer = ioin_a
      mastercompindx_oc_aer  = ioc_a
      mastercompindx_bc_aer  = ibc_a

! set hygroscopicities that are (temporarily) used
!    by sectional routines for hysteresis water
! for consistency with the mosaic water uptake,
!    only the inorganic salts/acids have non-negligible values
      hygro_mastercomp_aer(:) = 1.0e-10
      hygro_mastercomp_aer( mastercompindx_so4_aer ) = 0.5
      hygro_mastercomp_aer( mastercompindx_nh4_aer ) = 0.5
      hygro_mastercomp_aer( mastercompindx_no3_aer ) = 0.5
      hygro_mastercomp_aer( mastercompindx_msa_aer ) = 0.58
      hygro_mastercomp_aer( mastercompindx_na_aer  ) = 1.16
      hygro_mastercomp_aer( mastercompindx_cl_aer  ) = 1.16
      hygro_mastercomp_aer( mastercompindx_ca_aer  ) = 0.1
      hygro_mastercomp_aer( mastercompindx_co3_aer ) = 0.1

! 2010-02-26 - get hygro from mosaic kappa, always
! 2010-03-28 - only do this when ipmcmos > 0,
!    to allow backwards compatibility with mos27a
      if (ipmcmos > 0) then
         do ll = 1, naer
            hygro_mastercomp_aer(ll) = kappa_aer_mac(ll)
         end do
      end if


      ibin_of_isize_itype(:,:) = -999888777
      isize_of_ibin(:) = -999888777
      itype_of_ibin(:) = -999888777
      itype_of_itype_md1md2(:,:) = -999888777
      itype_md1_of_itype(:) = -999888777
      itype_md2_of_itype(:) = -999888777

      ncomp_aer(:) = 0
      ncomp_plustracer_aer(:) = 0
      numptr_aer(:,:,:) = -999888777
      waterptr_aer(:,:) = -999888777
      hyswptr_aer(:,:) = -999888777
      massptr_aer(:,:,:,:) = -999888777
      mastercompptr_aer(:,:) = -999888777

      dens_aer(:,:) = 0.0
      mw_aer(:,:) = 0.0
      hygro_aer(:,:) = 0.0

      iphase = ai_phase
      ibin = 0
      do itype = 1, ntype_aer
         ncomp_aer(itype) = naer
         ncomp_plustracer_aer(itype) = naer
         do isize = 1, nsize_aer(itype)
            ibin = ibin + 1
            noffset = ngas_max + naer_tot*(ibin - 1)
            numptr_aer(isize,itype,iphase) = noffset + knum_a
            waterptr_aer(isize,itype) = noffset + kwater_a
            hyswptr_aer(isize,itype) = noffset + kjhyst_a
            do ll = 1, ncomp_aer(itype)
               massptr_aer(ll,isize,itype,iphase) = &
                  noffset + (naer_tot-naer) + ll
               mastercompptr_aer(ll,itype) = ll
               name_aer( ll,itype) = name_mastercomp_aer(ll)
               dens_aer( ll,itype) = dens_mastercomp_aer(ll)
               mw_aer(   ll,itype) = mw_mastercomp_aer(ll)
               hygro_aer(ll,itype) = hygro_mastercomp_aer(ll)
            end do

            ibin_of_isize_itype(isize,itype) = ibin
            isize_of_ibin(ibin) = isize
            itype_of_ibin(ibin) = itype
         end do ! isize
      end do ! itype


! generate initial dry diameters automatically
      dlo_sect( :,:) = 0.0
      dhi_sect( :,:) = 0.0
      dcen_sect(:,:) = 0.0
      dcut_sect(:,:) = 0.0
      volumlo_sect( :,:) = 0.0
      volumhi_sect( :,:) = 0.0
      volumcen_sect(:,:) = 0.0
      volumcut_sect(:,:) = 0.0
      
      if ((maersize_init_flag1 == 2) .or. &
           (msize_framework == msectional)) then


         
         if ((dlo_aersize_init <= 0.0) .or. &
             (dlo_aersize_init >= dhi_aersize_init)) then
            write(*,'(2a,1p,2e11.3)') &
               '*** init_aerosol fatal error - ', &
               'bad dlo/hi_aersize_init =', &
                dlo_aersize_init, dhi_aersize_init
            stop
         end if

         !PRINT*, 'dlo = ',dlo_aersize_init
         !PRINT*, 'dhi = ',dhi_aersize_init
         !STOP
         
         
         
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         tmpa = LOG(dlo_aersize_init)
         !tmpb = (log(dhi_aersize_init) - tmpa)/nsize_aer(1)

         tmpb = LOG(dhi_aersize_init/dlo_aersize_init)/nsize_aer(1)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
         
         ibin = 0
         
         do itype = 1, ntype_aer
         do isize = 1, nsize_aer(itype)
               
! dlo_sect & dhi_sect are bin lower/upper bound dry diam (cm)
            dlo_sect( isize,itype) = exp(tmpa + (isize-1  )*tmpb)*1.0e-4_r8
            dhi_sect( isize,itype) = exp(tmpa + (isize    )*tmpb)*1.0e-4_r8
            
            dcen_sect(isize,itype) = exp(tmpa + (isize-0.5)*tmpb)*1.0e-4_r8
            dcut_sect(isize,itype) = dhi_sect(isize,itype)

            volumlo_sect( isize,itype) = (dlo_sect( isize,itype)**3)*pi/6.0_r8
            volumhi_sect( isize,itype) = (dhi_sect( isize,itype)**3)*pi/6.0_r8
            
            volumcen_sect(isize,itype) = (dcen_sect(isize,itype)**3)*pi/6.0_r8
            volumcut_sect(isize,itype) = volumhi_sect(isize,itype)

            ibin = ibin + 1
            
            if (maersize_init_flag1 == 2) then
               noffset = ngas_max + naer_tot*(ibin - 1)
! cnn(noffset + kdpdry_a) is bin initial dry diameter (micron)
               cnn(noffset + kdpdry_a) = exp(tmpa + (isize-0.5)*tmpb)
            end if
         end do ! isize

         dcut_sect(0,itype) = dlo_sect(1,itype)
         volumcut_sect(0,itype) = volumlo_sect(1,itype)
         end do ! itype

      end if


! sigmag
      do itype = 1, ntype_aer
      do isize = 1, nsize_aer(itype)
            
         ibin = ibin_of_isize_itype(isize,itype)
         noffset = ngas_max + naer_tot*(ibin - 1)
         
         sigmag_aer(isize,itype) = cnn(noffset + ksigmag_a)
         
         if (msize_framework == munstructured ) then
!           assume monodisperse for each "bin"
            sigmag_aer(isize,itype) = 1.0
         else if (msize_framework == msectional) then
!           sigmag based on width of bin - following is same as
!           log(sigmag) = [log(dhi) - log(dlo)]/sqrt(12)
            sigmag_aer(isize,itype) = &
               (dhi_sect(isize,itype)/dlo_sect(isize,itype))**0.289
         end if
         
         cnn(noffset + ksigmag_a) = sigmag_aer(isize,itype)
      end do ! isize
      end do ! itype


! number
!    initialize it to zero, which will cause conform_aerosol_number
!       to initialize it correctly on the first time step
!*** this seems to work, but code should really be changed so that
!       aerosol number is initialized BEFORE the first time step
      do itype = 1, ntype_aer
      do isize = 1, nsize_aer(itype)
         cnn(numptr_aer(isize,itype,iphase)) = 0.0
      end do ! isize
      end do ! itype


! names for cnn array
      do l = ngas_max+1, ntot_max
         write( species(l), '(a,i9.9)' ) 's', l
      end do
      
      
      do ibin = 1, nbin_a
         noffset = ngas_max + naer_tot*(ibin - 1)
         write(tmpch5,'(a,i3.3)') '_a', mod(ibin,1000)
         species(noffset+knum_a   ) = 'num'    // tmpch5
         species(noffset+kdpdry_a ) = 'dpdry'  // tmpch5
         species(noffset+ksigmag_a) = 'sigmag' // tmpch5
         species(noffset+kjhyst_a ) = 'jhyst'  // tmpch5
         species(noffset+kwater_a ) = 'water'  // tmpch5
         do l = 1, naer
             species(noffset+5+l) = trim(aer_name(l)) // tmpch5
         end do
      end do


! 3d sectional stuff
      if (msectional_flag2 > 0) then

! mappings between itype and (itype_md1,itype_md2)
         itype = 0
         do it2 = 1, ntype_md2_aer
         do it1 = 1, ntype_md1_aer
            itype = itype + 1
            itype_of_itype_md1md2(it1,it2) = itype
            itype_md1_of_itype(itype) = it1
            itype_md2_of_itype(itype) = it2
         end do
         end do

! cut values for bc mass fraction
! when method_atype_md1_init <= 1, they are read from input file
! otherwise, they have a uniform linear spacing between the
!    specified xcutlo & xcuthi values
         if (method_atype_md1_init >= 2) then
            tmpa = xcutlo_atype_md1_init
            tmpb = (xcuthi_atype_md1_init - tmpa) / ntype_md1_aer
            do it1 = 0, ntype_md1_aer
               xcut_atype_md1(it1) = tmpa + tmpb*it1
            end do
         else
            tmpb = 1.0
         end if
         if ( (tmpb <= 0.0) .or. &
              (xcut_atype_md1(1) < 0.0) .or. &
              (xcut_atype_md1(ntype_md1_aer) < 0.1) ) then
            write(*,'(2a)') &
               '*** init_aerosol fatal error - ', &
               'bad xcut_atype_md1 ='
            write(*,'(1p,7e11.3)') &
               xcut_atype_md1(0:ntype_md1_aer)
            stop
         end if

! cut values for hygroscopicity
! when method_atype_md2_init <= 1, they are read from input file
! otherwise, they have a uniform logarithmic spacing between the
!    specified xcutlo & xcuthi values
         if (method_atype_md2_init >= 2) then
            tmpa = log(xcutlo_atype_md2_init)
            tmpb = (log(xcuthi_atype_md2_init) - tmpa) / ntype_md2_aer
            do it2 = 0, ntype_md2_aer
               xcut_atype_md2(it2) = exp( tmpa + tmpb*it2 )
            end do
         else
            tmpb = 1.0
         end if
         if ( (tmpb <= 0.0) .or. &
              (xcut_atype_md2(0) <= 0.0) .or. &
              (xcut_atype_md2(1) < 0.0001) .or. &
              (xcut_atype_md2(ntype_md2_aer) < 0.1) ) then
            write(*,'(2a)') &
               '*** init_aerosol fatal error - ', &
               'bad xcut_atype_md2 ='
            write(*,'(1p,7e11.3)') &
               xcut_atype_md2(0:ntype_md2_aer)
            stop
         end if

      else
         xcut_atype_md1(0:1) = (/ 0.0_r8, 1.0_r8 /)
         xcut_atype_md2(0:1) = (/ 0.1_r8, 1.0_r8 /)

      end if ! (msectional_flag2 > 0)


! diagnostic output
      lunaa = lun_sect_171

      lunaa = 0 ! ADDED BY CHARLES

      
      if (lunaa > 0) then

      write(lunaa,'(//a,i5)') 'ntot_mastercomp_aer', ntot_mastercomp_aer
      do ll = 1, ntot_mastercomp_aer
         write(lunaa,'(a,i4,2x,a,1p,4e12.4)') 'll,name,dens,mw,hy', &
            ll, name_mastercomp_aer(ll), dens_mastercomp_aer(ll), &
            mw_mastercomp_aer(ll), hygro_mastercomp_aer(ll)
      end do

      iphase = ai_phase
      do itype = 1, ntype_aer
      write(lunaa,'(//a,2i5)') 'itype, ncomp', itype, ncomp_aer(itype)
      do ll = 1, ncomp_aer(itype)
         write(lunaa,'(a,i4,2x,a,1p,4e12.4)') 'll,name,dens,mw', &
            ll, name_aer(ll,itype), dens_aer(ll,itype), mw_aer(ll,itype)
      end do
      do isize = 1, nsize_aer(itype)
         write(lunaa,'(/a,5i5)')     'itype,isize,numptr,waterptr', &
            itype, isize, numptr_aer(isize,itype,iphase), &
            waterptr_aer(isize,itype)
         do ll = 1, ncomp_aer(itype)
            write(lunaa,'(a,3i5,i7,2x,a)') 'itype,isize,ll,massptr,name', &
               itype, isize, ll, massptr_aer(ll,isize,itype,iphase), &
               name_aer(ll,itype)
         end do
      end do
      end do

      end if   ! (lunaa > 0)


      ! pmcmos stuff
      
      if (ipmcmos > 0) call pmcmos_init_aerosol


      return
      end subroutine init_aerosol

