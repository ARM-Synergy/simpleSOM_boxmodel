'''===================================================================================
        THIS IS THE OBJECTIVE FUNCTION FOR FITTING THE SIMPLESOM-MOSAIC MODEL
==================================================================================='''

import os
import numpy             as np
import pandas            as pd
import matplotlib.pyplot as plt

from data_mod.counter import ctr

import data_mod.mdata as mdata
import variables  as ns

from functions.f_comp_model import comp_model
from functions.f_run_model  import run_model

def fit_model(x00,dex=None,x01=None):
    
    # RE-COMBINE PARAMETERS:
    # ================================================================================
    dex1 = dex[0]
    dex2 = dex[1]
    
    xx = np.zeros(max(max(dex1),max(dex2))+1)

    for i,dex in enumerate(dex1):
        xx[dex] = x01[i]

    for i,dex in enumerate(dex2):
        xx[dex] = x00[i]
    
    # RECORD PARAMETER VALUES AT EACH EVAL.:
    # ================================================================================
    global ctr
    
    with open('outputs/diag.params','a') as ff:
        
        ff.write('%03i... \n'%(ctr))
        ff.write(('%.4f '*8)%(xx[0],np.exp(-xx[1]),xx[2]/100.,xx[3],
                              xx[4]/sum(xx[4:]),
                              xx[5]/sum(xx[4:]),
                              xx[6]/sum(xx[4:]),
                              xx[7]/sum(xx[4:]))+'\n\n'); ff.close()
    
    ctr += 1

    # COMPILE THE MODEL:
    # ================================================================================
    comp_model('MOSAIC.v0.DEV_OXY',xx)
    
    # RUN THE MODEL:
    # ================================================================================
    os.system('rm ../outputs/* 2>/dev/null')
    
    output = run_model(ns)
    df_SOA = output['df_SOA']
    df_O2C = output['df_O2C']

    t1 = df_SOA['time'].values
    y1 = df_SOA['SOA'].values

    t2 = df_O2C[0].values
    y2 = df_O2C[3].values
        
    # CALCULATE OBJECTIVE FUNCTION:
    # ================================================================================
    x_obs = mdata.x_obs
    y_obs = mdata.y_obs

    sample = np.linspace(0.,t1[-1],num=100)
    y_mod  = np.interp(sample,t1,y1)
    y_obs  = np.interp(sample,x_obs,y_obs)
    
    diff1 = abs(y_mod - y_obs)
    diff1 = diff1[~np.isnan(diff1)]
    diff1 = diff1/np.max(diff1)*100.

    # CALCULATE OBJECTIVE FUNCTION:
    # ================================================================================
    x_obs = mdata.x_obs
    y_obs = mdata.z_obs

    sample = np.linspace(0.5,t2[-1],num=100)
    y_mod  = np.interp(sample,t2,y2)
    y_obs  = np.interp(sample,x_obs,y_obs)
    
    diff2 = abs(y_mod - y_obs)
    diff2 = diff2[~np.isnan(diff2)]
    diff2 = diff2/np.max(diff2)*100.

    return diff1 + diff2
