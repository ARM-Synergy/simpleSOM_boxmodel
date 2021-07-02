'''================================================================
                THIS FUNCTION WRITES THE INPUT FILE
================================================================'''

import numpy  as np
import pandas as pd

def write_inputs(ns):

    # RUN NAME:
    run = ns.run_name
    
    # READ TEMPLATE:
    with open('../inputs/raw.inp','r') as ff:
        lines = ff.readlines()
    
    # WRITE SIMULATION TIME:    
    lines[4] = '%s %s %s %s\n'%(tuple(ns.t_end.split(':')))

    # WRITE PRESSURE:
    lines[11] = '%.1f\n'%ns.Pres

    # WRITE TEMPERATURE:
    lines[10] = '%.1f\n'%ns.Temp

    # WRITE RELATIVE HUMIDITY:
    lines[9] = '%.1f\n'%ns.RH

    # WRITE GAS CONCENTRATIONS:
    nSPECIES = ns.df_GAS.shape[0]
    
    for i in range(nSPECIES):
        lines[39+i] = '%i\t%s\t%.6f\t%.6f\t%.6f\n'%(tuple(list(ns.df_GAS.iloc[i].values)))

    # WRITE INITIAL AEROSOL DISTRIBUTION:
    nBINS = ns.nBINS

    for i in range(nBINS):
        
        add = ['0.00']*29
        
        add[0] = '%i'%(i + 1)
        add[2] = '%.3f'%ns.DIAM[i]
        add[3] = '1.0'
        add[6] = '%.2e'%ns.INIT_MOLE[i]
        
        lines[123+i] = '%s\n'%(' '.join(add))

    # WRITE NUMBER OF BINS:
    lines[13] = '%i\n'%nBINS

    # WRITE THE SIZE RANGE:
    lines[120] = '%.3f %.3f\n'%(ns.BD_LOWER,ns.BD_UPPER)
        
    # WRITE TO THE INPUT FILE:
    with open('../inputs/%s.inp'%run,'w') as ff:
        ff.writelines(lines)
    
    # WRITE TO THE INPUT FILE:
    with open('../inputs/input.b','w') as ff:
        ff.writelines('%i\n'%ns.HET)
        ff.writelines('%i\n'%ns.OLIG)
        ff.writelines('%.2e\n'%ns.ko_f)
        ff.writelines('%.2e\n'%ns.ko_d)
        ff.writelines('%i\n'%ns.VWL)
        ff.writelines('%.2e\n'%ns.kvap_on)
        ff.writelines('%.2e\n'%ns.NUC)
        ff.writelines('%.2e\n'%ns.NUC_T0)
        ff.writelines('%.2e\n'%ns.NUC_T1)
    
    # WRITE TO THE INPUT FILE:
    with open('../inputs/in.VOC','w') as ff:

        p1 = ns.VOC['Conc0']
        p2 = ns.VOC['MW']
        p3 = ns.VOC['kOH']
        p4 = ns.VOC['Csat']
        
        ff.writelines((4*'%.4e ')%(p1,p2,p3,p4))


    # WRITE TO THE INPUT FILE:
    with open('../inputs/in.OH','w') as ff:

        p1 = ns.OH['AX1']
        p2 = ns.OH['BX1']
        p3 = ns.OH['AX2']
        p4 = ns.OH['BX2']
        
        ff.writelines((4*'%.4e ')%(p1,p2,p3,p4))
    
