!**********************************************************************************
! This computer software was prepared by Battelle Memorial Institute, hereinafter
! the Contractor, under Contract No. DE-AC05-76RL0 1830 with the Department of
! Energy (DOE). NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Copyright (c) 2018 Battelle Memorial Institute
! Written by Rahul A. Zaveri
!
!********************************************************************************************
! Model         : MOSAIC (Model for Simulating Aerosol Interactions & Chemistry)
!
!
! Sub-Modules   : CBM-Z  (gas-phase photochemistry)
!                 ASTEM  (Adaptive Step Time-split Euler Method)
!                 MESA   (Multicomponent Equilibrium Solver for Aerosols)
!                 MTEM   (Multicomponent Taylor Expansion Method for activity coefficients)
!
! Purpose       : CBM-Z simulates lower tropospheric trace gas photochemistry for
!                 background, urban, biogenic, and marine (DMS) sub-regimes.
!
!                 ASTEM solves the dynamic partitioning of semi-volatile
!                 species between gas and particle phases.
!
!                 MESA solves the multicomponent solid-liquid equilibria within
!                 the aerosol phase.
!
!                 MTEM computes the multicomponent activity coefficients of
!                 electrolytes in aqueous atmospheric aerosols.
!
! Author        : Rahul A. Zaveri, PhD
!                 Senior Research Scientist
!                 Pacific Northwest National Laboratory
!                 Atmospheric Sciences Technical Group
!                 P.O. Box 999, MSIN K9-30
!                 Richland, WA 99352
!                 Phone: (509) 372-6159, Fax: (509) 372-6168
!                 Email: Rahul.Zaveri@pnl.gov
!
! Bugs/Problems : Please report any bugs or problems to Rahul.Zaveri@pnnl.gov
!
! Terms of Use  : (1) MOSAIC and its submodules CBM-Z, ASTEM, MESA, and MTEM may not be
!                     included in any commercial package, or used for any commercial
!                     applications without prior authorization from the author.
!                 (2) The MOSAIC code may be used for educational or non-profit purposes
!                     only. Any other usage must be first approved by the author.
!                 (3) The MOSAIC code cannot be modified in any way or form or distributed
!                     without the author's prior consent.
!                 (4) No portion of the MOSAIC source code can be used in other codes
!                     without the author's prior consent.
!                 (5) The MOSAIC code is provided on an as-is basis, and the author
!                     bears no liability from its usage.
!                 (6) Publications resulting from the usage of MOSAIC must cite
!                     the references below for proper acknowledgment.
!
! References    : Zaveri R.A., R.C. Easter, J.D. Fast, and L.K. Peters, Model
!                   for simulating aerosol interactions and chemistry (MOSAIC),
!                   J. Geophys. Res., in review.
!
!                 Zaveri R.A., R.C. Easter, and L.K. Peters (2005a) A computationally
!                   efficient multicomponent equilibrium solver for aerosols (MESA),
!                   J. Geophys. Res, 110, D24203, doi:10.1029/2004JD005618.
!
!                 Zaveri R.A., R.C. Easter, and A.S. Wexler (2005b) A new method for
!                   multicomponent activity coefficients of electrolytes in aqueous
!                   atmospheric aerosols, J. Geophys. Res., 110, D02201,
!                   doi:10.1029/2004JD004681.
!
!                 Zaveri R.A. and L.K. Peters (1999) A new lumped structure photochemical
!                   mechanism for large-scale applications
!
! Support       : Funding for the development and evaluation of MOSAIC and
!                 its sub-modules was provided by:
!                 (a) the U.S. Department of Energy (DOE) under the auspices of the
!                     Atmospheric Science Program (ASP) of the Office of Biological and
!                     Environmental Research
!                 (b) the NASA Aerosol Program and NASA Earth Science Enterprise
!                 (c) the U.S. Environmental Protection Agency (EPA) Aerosol Program
!                 (d) PNNL Laboratory Directed Research and Development (LDRD) Program
!------------------------------------------------------------------------

      program main
      
      USE mod_MAIN
      USE mod_GAS

      use module_data_mosaic_aero
      use module_data_mosaic_cloud

      IMPLICIT NONE


      write(6,*)'   '
      write(6,*)'*****************************************************'
      write(6,*)'                     MOSAIC'
      write(6,*)'Model for Simulating Aerosol Interactions & Chemistry'
      write(6,*)'    Copyright (c) 2018 Battelle Memorial Institute'
      write(6,*)'  '
      write(6,*)'    Contact: Rahul A. Zaveri (rahul.zaveri@pnl.gov)'
      write(6,*)'       Pacific Northwest National Laboratory'
      write(6,*)'*****************************************************'
      write(6,*)'   '
      write(6,*)'simulation begins...'
      
      call init_data_modules			! initializes various indices
      
      call SetIOfiles				! reads inputfile
      
      call SetRunParameters
      
      call SetAirComposition
      
      call init_aerosol
      
      call LoadPeroxyParameters			! Aperox and Bperox
      
!!      call DoMassBalance			! initial elemental mass balance
      
      IF (mmode.EQ.1) THEN
         
        call time_integration_mode
         
      elseif(mmode .eq. 2)then
        call parametric_analysis_mode
      endif

      
      WRITE(6,'(/,A,/)') 'SIMULATION COMPLETE.'      
      END PROGRAM
