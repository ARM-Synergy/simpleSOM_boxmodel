#=============================================================================================
#                             THIS IS THE MAKEFILE FOR MOSAIC.EXE
#=============================================================================================

COMP1 = ifort
COMP2 = icc

FLAGS1 = -debug -O0
FLAGS2 = #-Mfree
FLAGS3 =

PATH0 = objects/
PATH1 = source_ftn/datamodules/
PATH2 = source_ftn/main/
PATH3 = source_ftn/aerosol/
PATH4 = source_ftn/cloud/
PATH5 = source_ftn/gas/
PATH6 = source_ftn/solver/
PATH7 = source_ftn/main/

DATAMODU =  \
$(PATH0)module_data_mosaic_kind.o  \
$(PATH0)module_data_mosaic_main.o  \
$(PATH0)module_data_mosaic_aero.o  \
$(PATH0)module_data_mosaic_asect.o \
$(PATH0)module_data_mosaic_gas.o   \
$(PATH0)module_data_mosaic_cloud.o \
$(PATH0)module_data_mosaic_pmcmos.o

PMCMOS =  \
$(PATH0)module_pmcmos_aer.o   \
$(PATH0)module_pmcmos_subaa.o \
$(PATH0)module_pmcmos_init.o

AEROSOL =  \
$(PATH0)module_peg_util.o  \
$(PATH0)module_coag.o  \
$(PATH0)module_coag_3d.o  \
$(PATH0)module_movesect.o  \
$(PATH0)module_movesect_3d.o  \
$(PATH0)module_newnuc.o  \
$(PATH0)module_sect_iface.o  \
$(PATH0)aerchemistry.o  \
$(PATH0)aeroptical.o  \
$(PATH0)ASTEM.o  \
$(PATH0)cputime.o  \
$(PATH0)init_aerosol.o  \
$(PATH0)miedriver.o  \
$(PATH0)mosaic_box.25.o  \
$(PATH0)mosaic_lsode.o \
$(PATH0)r_step_hetchem.o

CLOUD =  \
$(PATH0)cldchemistry.o

GAS =  \
$(PATH0)function_arr.o  \
$(PATH0)function_troe.o  \
$(PATH0)gaschemistry.o  \
$(PATH0)gasintegrator.o  \
$(PATH0)gasrateconstants_bio.o  \
$(PATH0)gasrateconstants_com.o  \
$(PATH0)gasrateconstants.o  \
$(PATH0)gasrateconstants_het.o  \
$(PATH0)gasrateconstants_mar.o  \
$(PATH0)gasrateconstants_urb.o  \
$(PATH0)gasrates_bio.o  \
$(PATH0)gasrates_com.o  \
$(PATH0)gasrates.o  \
$(PATH0)gasrates_het.o  \
$(PATH0)gasrates_mar.o  \
$(PATH0)gasrates_urb.o  \
$(PATH0)jac.o  \
$(PATH0)loadperoxyparameters.o  \
$(PATH0)mapgas_bio.o  \
$(PATH0)mapgas_com.o  \
$(PATH0)mapgas_mar.o  \
$(PATH0)mapgasspecies.o  \
$(PATH0)mapgas_urb.o  \
$(PATH0)ode_bio.o  \
$(PATH0)ode_com.o  \
$(PATH0)ode_gas.o  \
$(PATH0)ode_het.o  \
$(PATH0)ode_mar.o  \
$(PATH0)ode_urb.o  \
$(PATH0)peroxyrateconstants.o  \
$(PATH0)photoconstants_fixed.o  \
$(PATH0)photoconstants_solar.o  \
$(PATH0)photoparam1.o  \
$(PATH0)photoparam2.o  \
$(PATH0)selectgasregime.o  \
$(PATH0)setgas_bio.o  \
$(PATH0)setgas_com.o  \
$(PATH0)setgasindices.o  \
$(PATH0)setgas_mar.o  \
$(PATH0)setgas_urb.o  \
$(PATH0)solarzenithangle.o

SOLVER =  \
$(PATH0)dlsode.o

MAIN =  \
$(PATH0)domassbalance.o  \
$(PATH0)function_nbllen.o  \
$(PATH0)function_watervapor.o  \
$(PATH0)init_data_modules.o  \
$(PATH0)integratechemistry.o  \
$(PATH0)main.o  \
$(PATH0)parametric_analysis_mode.o  \
$(PATH0)print_aer.o  \
$(PATH0)print_gas.o  \
$(PATH0)readinputfile.o  \
$(PATH0)setaircomposition.o  \
$(PATH0)setiofiles.o  \
$(PATH0)setrunparameters.o  \
$(PATH0)time_integration_mode.o  \
$(PATH0)updateemissions.o  \
$(PATH0)updatemetfields.o  \
$(PATH0)updatetime.o

# EXECUTABLE:
#=============================================================================================
mosaic.exe:  $(DATAMODU)  $(PMCMOS)  $(AEROSOL)  $(CLOUD)  $(GAS)  $(SOLVER)  $(MAIN)
	$(COMP1)  $(FLAGS1)  -o $@  $(DATAMODU)  $(PMCMOS)  $(AEROSOL)  $(CLOUD)  $(GAS)  $(SOLVER)  $(MAIN)

# GROUP-1 OBJECTS:
#=============================================================================================
$(PATH0)module_data_mosaic_kind.o :: $(PATH1)module_data_mosaic_kind.f90
	$(COMP1) -c $(FLAGS1) $< -o $@ -module $(PATH0)

$(PATH0)module_data_mosaic_main.o :: $(PATH1)module_data_mosaic_main.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)module_data_mosaic_aero.o : $(PATH1)module_data_mosaic_aero.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)module_data_mosaic_asect.o : $(PATH1)module_data_mosaic_asect.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)module_data_mosaic_gas.o : $(PATH1)module_data_mosaic_gas.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)module_data_mosaic_cloud.o : $(PATH1)module_data_mosaic_cloud.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)module_data_mosaic_pmcmos.o : $(PATH1)module_data_mosaic_pmcmos.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

# GROUP-2 OBJECTS:
#=============================================================================================
$(PATH0)module_pmcmos_aer.o : $(PATH2)module_pmcmos_aer.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)module_pmcmos_subaa.o : $(PATH2)module_pmcmos_subaa.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)module_pmcmos_init.o : $(PATH2)module_pmcmos_init.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

# GROUP-3 OBJECTS:
#=============================================================================================
$(PATH0)aerchemistry.o : $(PATH3)aerchemistry.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@

$(PATH0)aeroptical.o : $(PATH3)aeroptical.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@

$(PATH0)ASTEM.o : $(PATH3)ASTEM.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@

$(PATH0)cputime.o : $(PATH3)cputime.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@

$(PATH0)init_aerosol.o : $(PATH3)init_aerosol.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@

$(PATH0)miedriver.o : $(PATH3)miedriver.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@

$(PATH0)module_coag.o : $(PATH3)module_coag.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)module_coag_3d.o : $(PATH3)module_coag_3d.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)module_movesect.o : $(PATH3)module_movesect.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)module_movesect_3d.o : $(PATH3)module_movesect_3d.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)module_newnuc.o : $(PATH3)module_newnuc.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)module_peg_util.o : $(PATH3)module_peg_util.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)module_sect_iface.o : $(PATH3)module_sect_iface.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)mosaic_box.25.o : $(PATH3)mosaic_box.25.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)mosaic_lsode.o : $(PATH3)mosaic_lsode.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)r_step_hetchem.o : $(PATH3)r_step_hetchem.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

# GROUP-4 OBJECTS:
#=============================================================================================
$(PATH0)cldchemistry.o : $(PATH4)cldchemistry.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

# GROUP-5 OBJECTS:
#=============================================================================================
$(PATH0)function_arr.o : $(PATH5)function_arr.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)function_troe.o : $(PATH5)function_troe.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)gaschemistry.o : $(PATH5)gaschemistry.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)gasintegrator.o : $(PATH5)gasintegrator.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)gasrateconstants_bio.o : $(PATH5)gasrateconstants_bio.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)gasrateconstants_com.o : $(PATH5)gasrateconstants_com.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)gasrateconstants.o : $(PATH5)gasrateconstants.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)gasrateconstants_het.o : $(PATH5)gasrateconstants_het.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)gasrateconstants_mar.o : $(PATH5)gasrateconstants_mar.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)gasrateconstants_urb.o : $(PATH5)gasrateconstants_urb.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)gasrates_bio.o : $(PATH5)gasrates_bio.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)gasrates_com.o : $(PATH5)gasrates_com.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)gasrates.o : $(PATH5)gasrates.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)gasrates_het.o : $(PATH5)gasrates_het.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)gasrates_mar.o : $(PATH5)gasrates_mar.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)gasrates_urb.o : $(PATH5)gasrates_urb.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)jac.o : $(PATH5)jac.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)loadperoxyparameters.o : $(PATH5)loadperoxyparameters.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)mapgas_bio.o : $(PATH5)mapgas_bio.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)mapgas_com.o : $(PATH5)mapgas_com.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)mapgas_mar.o : $(PATH5)mapgas_mar.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)mapgasspecies.o : $(PATH5)mapgasspecies.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)mapgas_urb.o : $(PATH5)mapgas_urb.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)ode_bio.o : $(PATH5)ode_bio.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)ode_com.o : $(PATH5)ode_com.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)ode_gas.o : $(PATH5)ode_gas.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)ode_het.o : $(PATH5)ode_het.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)ode_mar.o : $(PATH5)ode_mar.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)ode_urb.o : $(PATH5)ode_urb.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)peroxyrateconstants.o : $(PATH5)peroxyrateconstants.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)photoconstants_fixed.o : $(PATH5)photoconstants_fixed.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)photoconstants_solar.o : $(PATH5)photoconstants_solar.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)photoparam1.o : $(PATH5)photoparam1.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)photoparam2.o : $(PATH5)photoparam2.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)selectgasregime.o : $(PATH5)selectgasregime.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)setgas_bio.o : $(PATH5)setgas_bio.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)setgas_com.o : $(PATH5)setgas_com.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)setgasindices.o : $(PATH5)setgasindices.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)setgas_mar.o : $(PATH5)setgas_mar.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)setgas_urb.o : $(PATH5)setgas_urb.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)solarzenithangle.o : $(PATH5)solarzenithangle.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

# GROUP-6 OBJECTS:
#=============================================================================================
$(PATH0)dlsode.o : $(PATH6)dlsode.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

# GROUP-7 OBJECTS:
#=============================================================================================
$(PATH0)domassbalance.o : $(PATH7)domassbalance.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)function_nbllen.o : $(PATH7)function_nbllen.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)function_watervapor.o : $(PATH7)function_watervapor.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)init_data_modules.o : $(PATH7)init_data_modules.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)integratechemistry.o : $(PATH7)integratechemistry.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)main.o : $(PATH7)main.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)parametric_analysis_mode.o : $(PATH7)parametric_analysis_mode.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)print_aer.o : $(PATH7)print_aer.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)print_gas.o : $(PATH7)print_gas.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)readinputfile.o : $(PATH7)readinputfile.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)setaircomposition.o : $(PATH7)setaircomposition.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)setiofiles.o : $(PATH7)setiofiles.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)setrunparameters.o : $(PATH7)setrunparameters.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)time_integration_mode.o : $(PATH7)time_integration_mode.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)updateemissions.o : $(PATH7)updateemissions.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)updatemetfields.o : $(PATH7)updatemetfields.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

$(PATH0)updatetime.o : $(PATH7)updatetime.f90
	$(COMP1) -c $(FLAGS1) $< -I$(PATH0) -o $@ -module $(PATH0)

# CLEAN:
#=============================================================================================
clean:
	rm mosaic.exe >/dev/null
	rm objects/*.o >/dev/null
	rm objects/*.mod >/dev/null

