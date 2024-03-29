'''==============================================================
            THIS IS THE WRAPPER FOR THE MOSAIC MODEL
=============================================================='''

import os
import numpy             as np
import pandas            as pd
import scipy.optimize    as sio
import matplotlib.pyplot as plt

from functions.f_feed4ward import feed4ward
from functions.f_fit_model import fit_model

# SIMPLE-SOM PARAMETERS:
# ===============================================================
#xx = np.array([3.51346,0.0106,34.00,1.6303,10.,8348.,3087.,10.])
#x0 = np.array([3.513,0.0111,34.00,1.63,10.,7040.,2600.,10.])
xx = np.array([3.513,0.0111,34,1.63,0.001,.704,.26,0.001]) #fwd
#xx = np.array([3.46173,0.00620,34.00000,1.63339,6.12317,7380.90358,2542.18336,7.81189]) #fit
# FIT CONFIGURATIONS:
# ===============================================================
# CONTROL SWITCH
# 0 FOR FORWARD
# 1 FOR FITTING:
switch = 0
fit_file = 'base_fwd.csv'
# RUN THE MODEL (FIT/FORWARD):
# ===============================================================
if switch == 1:
    
    # REFRESH RECORD FILE:
    f1 = open('outputs/diag.params','w').close()

    # INDEX OF VALUE(S) TO HOLD CONST.:
    dex1 = [2]
    dex2 = [i for i in range(len(x0)) if i not in dex1] 
    
    # GET SLICES:
    x00 = np.array(x0[dex2])
    x01 = np.array(x0[dex1])
    
    # LOWER AND UPPER BOUNDS:
    bL = [ 0.01,  0.001,   0.0, 1.000, 1e-2, 1e-2, 1e-2, 1e-2]
    bU = [10.00, 10.000, 100.0, 3.000,  1e4,  1e4,  1e4,  1e4]
    
    bL = tuple(np.array(bL)[dex2])
    bU = tuple(np.array(bU)[dex2])

    # RELATIVE DIFFERENTIAL STEP FOR FITTING:
    ds = [0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001]
    
    ds = tuple(np.array(ds)[dex2])
    
    # RUN FITTING:
    results = sio.least_squares(fit_model,x00,
                diff_step=ds,bounds=(bL,bU),loss='linear',
                ftol=1e-6,xtol=1e-6,verbose=0,max_nfev=10,
                kwargs={'dex':(dex1,dex2),'x01':x01})
    
    # PRINT FITS:
    fits = results.x

    oo = np.zeros(max(max(dex1),max(dex2))+1)

    for i,dex in enumerate(dex1):
        oo[dex] = x01[i]

    for i,dex in enumerate(dex2):
        oo[dex] = fits[i]

    with open('outputs/diag.params','a') as ff:
        ff.write('xx = np.array([%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f])'%tuple(list(oo)))
    
elif switch == 0:
    # VISUALIZE:
    feed4ward(xx, fit_file); print('FINISHED')

else:
    pass





