

      subroutine SetIOfiles

      use module_data_mosaic_main
      use module_data_mosaic_gas
      use module_data_mosaic_aero
      use module_data_mosaic_cloud
      use module_data_mosaic_pmcmos, only: &
         lun_pmcmos_in1, lun_pmcmos_in2
      use module_data_mosaic_asect, only: lunout

      implicit none

      integer ibin, lun
      integer nchar, nbllen, nbb
      character*7 bb

! INPUT FILES
      write(6,*)'Enter gas input filename. Example: case1.inp'

      READ(5,*) INPUTFILE

      lun_pmcmos_in1 = 16 ; lun_pmcmos_in2 = 17

      lun_sect_170 = 170
      lun_sect_171 = 171
      lun_sect_172 = 172
      lun_sect_180 = 180
      lun_sect_183 = 183
      lun_sect_184 = 184
      lun_sect_185 = 185
      lun_sect_186 = 186
      lun_sect_188 = 188
      lun_sect_190 = 190

      lunout = lun_sect_170

      lun_inp = 10
      
      OPEN(10,FILE='inputs/'//INPUTFILE)
      CALL ReadInputFile
      CLOSE(10)

!------------------------------------------------------------------
! OUTPUT FILES
!
        nchar = nbllen(inputfile) - 4

! GAS output file
      if(mgas .eq. mYES)then

        lun_gas = 19
        gas_output = 'outputs/'//inputfile(1:nchar)//'.gas.txt'
        open(lun_gas, file = gas_output)

      endif




! AEROSOL output files
      if(maer .eq. mYES)then	! UNCOMMENT THIS LINE


! species distribution output file
        lun_species = 22
        species_output= 'outputs/'//inputfile(1:nchar)//'.aero.txt'

        
        OPEN(lun_species, file = species_output)


! aerosol optical info output file
        lun_aeroptic = 23
        aeroptic_output=   &
                   'outputs/'//inputfile(1:nchar)//'.aeroptic.txt'
        if (maeroptic > 0) then
           open(lun_aeroptic, file = aeroptic_output)
        end if

      endif

      return
      end subroutine SetIOfiles


!---------------------------------------------------------------------
      subroutine lun_aer_open( ibin )
      use module_data_mosaic_aero
      use module_data_mosaic_main
      implicit none
      integer, intent(in) :: ibin

      if (nbin_a <= 20) return
      if ((ibin < 1) .or. (ibin > nbin_a)) then
         write(*,*) '*** lun_aer_open -- bad ibin = ', ibin
         stop
      end if

! open/write/close these files on each entry to avoid
! problems with lun_aer being too large
      if (lun_aer_status(ibin) <= 0) then
         open( unit=lun_aer(ibin), &
            file=aer_output(ibin), status='unknown' )
         lun_aer_status(ibin) = 1
      else
         open( unit=lun_aer(ibin), &
            file=aer_output(ibin), status='old', position='append' )
      endif
      return
      end subroutine lun_aer_open


!---------------------------------------------------------------------
      subroutine lun_aer_close( ibin )
      use module_data_mosaic_main
      use module_data_mosaic_aero
      implicit none
      integer, intent(in) :: ibin

      if (nbin_a <= 20) return
      if ((ibin < 1) .or. (ibin > nbin_a)) then
         write(*,*) '*** lun_aer_close -- bad ibin = ', ibin
         stop
      end if

! open/write/close these files on each entry to avoid
! problems with lun_aer being too large
      close( unit=lun_aer(ibin) )
      return
      end subroutine lun_aer_close


