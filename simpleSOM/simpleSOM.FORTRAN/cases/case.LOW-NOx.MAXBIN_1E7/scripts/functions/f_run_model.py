'''===================================================================
                  THIS FUNCTION RUNS THE MOSAIC MODEL
==================================================================='''

import os
import numpy  as np
import pandas as pd

from functions.f_write_inputs import write_inputs
from functions.f_read_aeout   import read_aeout

def run_model(ns):

    # RUN NAME:
    # ================================================================
    run = ns.run_name
    
    # WRITE THE INPUT FILE:
    # ================================================================
    write_inputs(ns)
    
    # RUN THE MODEL:
    # ================================================================
    os.system('cd ../outputs/ ; rm *')
    os.system('cd ../ ; echo %s.inp | ./mosaic.exe'%run)
    
    # READ AEROSOL OUTPUT FILE:
    # ================================================================
    df_AEOUT,df_SOA,df_SIZE,df_AENUM = read_aeout(ns)

    # READ VOC OUTPUT FILE:
    # ================================================================
    df_VOC = pd.read_csv('../outputs/out.VOC',header=None,delim_whitespace=True)

    # READ NUMBER CONC. OUTPUT FILE:
    # ================================================================
    df_NCONC = pd.read_csv('../outputs/out.NCONC',header=None,delim_whitespace=True)

    # READ O:C OUTPUT FILE:
    # ================================================================
    df_O2C = pd.read_csv('../outputs/out.O2C',header=None,delim_whitespace=True)
    
    # OUTPUT OBJECT:
    # ================================================================
    outs = {}
    
    outs['df_AEOUT'] = df_AEOUT
    outs['df_SOA']   = df_SOA
    outs['df_SIZE']  = df_SIZE
    outs['df_AENUM'] = df_AENUM
    outs['df_VOC']   = df_VOC
    outs['df_NCONC'] = df_NCONC
    outs['df_O2C']   = df_O2C
    
    return outs
