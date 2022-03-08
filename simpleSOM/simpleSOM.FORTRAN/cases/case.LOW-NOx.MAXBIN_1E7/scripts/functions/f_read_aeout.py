'''====================================================================================
                       THIS FUNCTION READS THE AEROSOL OUTPUT FILE
===================================================================================='''

import numpy             as np
import pandas            as pd
import matplotlib.pyplot as plt

def read_aeout(ns):
    
    # RUN NAME:
    # =================================================================================
    run = ns.run_name
    
    # READ MODEL OUTFILE:
    # =================================================================================
    df = pd.read_csv('../outputs/%s.aero.txt'%run,header=None,delim_whitespace=True)
    
    # DATAFRAME FOR SOA:
    # =================================================================================
    df_SOA = df[[0,5]].drop_duplicates().rename(columns={0:'time',5:'SOA'})
    
    # DATAFRAME FOR AEROSOL BIN SIZE:
    # =================================================================================
    df_TEMP = df[[0,1,2]].pivot_table(index=0,columns=1,values=2)
    df_SIZE = pd.DataFrame({'time':df_TEMP.index})
    
    for col in df_TEMP.columns:
       df_SIZE['BIN_%02i'%col] = df_TEMP[col].values
    
    # DATAFRAME FOR BIN PARTICLE NUMBER:
    # =================================================================================
    df_TEMP  = df[[0,1,3]].pivot_table(index=0,columns=1,values=3)
    df_AENUM = pd.DataFrame({'time':df_TEMP.index})
    
    for col in df_TEMP.columns:
       df_AENUM['BIN_%02i'%col] = df_TEMP[col].values

    # DATAFRAME FOR BIN SO4 MOLES:
    # =================================================================================
    df_TEMP = df[[0,1,6]].pivot_table(index=0,columns=1,values=6)
    df_SO4  = pd.DataFrame({'time':df_TEMP.index})
    
    for col in df_TEMP.columns:
       df_SO4['BIN_%02i'%col] = df_TEMP[col].values

    # DATAFRAME FOR BIN NH4 MOLES:
    # =================================================================================
    df_TEMP = df[[0,1,7]].pivot_table(index=0,columns=1,values=7)
    df_NH4  = pd.DataFrame({'time':df_TEMP.index})
    
    for col in df_TEMP.columns:
       df_NH4['BIN_%02i'%col] = df_TEMP[col].values

    # DATAFRAME FOR BIN H2O MOLES:
    # =================================================================================
    df_TEMP = df[[0,1,9]].pivot_table(index=0,columns=1,values=9)
    df_H2O  = pd.DataFrame({'time':df_TEMP.index})
    
    for col in df_TEMP.columns:
       df_H2O['BIN_%02i'%col] = df_TEMP[col].values

    # READ OLIGOMER MASS:
    # =================================================================================
    df_OLIG = pd.read_csv('../outputs/dat.olgmass',header=None,delim_whitespace=True)
    df_OLIG = df_OLIG.rename(columns={0:'time',1:'OLIG'})
    
    df_SOA['SOA']  = df_OLIG['OLIG'].values + df_SOA['SOA'].values  
    df_SOA['OLIG'] = df_OLIG['OLIG'].values

    
    return df,df_SOA,df_SIZE,df_AENUM
