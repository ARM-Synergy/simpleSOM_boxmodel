'''==================================================================
     THIS FUNCTION RUNS THE PREPROCESSOR TO GENERATE MODEL FILES 
=================================================================='''

import os
import numpy  as np
import pandas as pd

from functions.f_write_mmss import write_mmss

def use_prep(path,params,ns):

    # CREATE THE SPECIES AND MECHANISM FILES:
    # ===============================================================
    write_mmss(path,params,ns)

    '''
    # GENERATE THE MODEL FILES:
    # ===============================================================
    cmd1 = "cd %s/source_ftn/gas/OrganicsMech/"%path
    cmd2 = "printf 'simplesom.ss\nsimplesom.mm' | ./zap.exe"
    
    os.system('%s;%s > /dev/null'%(cmd1,cmd2))

    # WRITE TO THE SOURCE CODE:
    # ===============================================================
    # ODE FILE:
    with open('%s/source_ftn/gas/raws/raw.ode_bio.f90'%path,'r') as f1:
        lines1 = np.array(f1.readlines())

    x1 = np.where(lines1=='!FLAG1\n')[0][0]
    x2 = np.where(lines1=='!FLAG2\n')[0][0]
    
    with open('%s/source_ftn/gas/OrganicsMech/ode.f'%path,'r') as f2:
        lines2 = np.array(f2.readlines())
    
    lines3 = np.hstack((lines1[:x1+1],lines2,lines1[x2:]))

    with open('%s/source_ftn/gas/ode_bio.f90'%path,'w') as f3:
        f3.writelines(lines3)

    # RATE FILE:
    with open('%s/source_ftn/gas/raws/raw.gasrates_bio.f90'%path,'r') as f1:
        lines1 = np.array(f1.readlines())

    x1 = np.where(lines1=='!FLAG1\n')[0][0]
    x2 = np.where(lines1=='!FLAG2\n')[0][0]

    with open('%s/source_ftn/gas/OrganicsMech/rate.f'%path,'r') as f2:
        lines2 = np.array(f2.readlines())

    lines3 = np.hstack((lines1[:x1+1],lines2,lines1[x2:]))
        
    with open('%s/source_ftn/gas/gasrates_bio.f90'%path,'w') as f3:
        f3.writelines(lines3)

    # RATE CONSTANTS FILE:
    with open('%s/source_ftn/gas/raws/raw.gasrateconstants_bio.f90'%path,'r') as f1:
        lines1 = np.array(f1.readlines())
    
    x1 = np.where(lines1=='!FLAG1\n')[0][0]
    x2 = np.where(lines1=='!FLAG2\n')[0][0]

    with open('%s/source_ftn/gas/OrganicsMech/rconst.f'%path,'r') as f2:
        lines2 = np.array(f2.readlines())

    lines3 = np.hstack((lines1[:x1+1],lines2,lines1[x2:]))

    with open('%s/source_ftn/gas/gasrateconstants_bio.f90'%path,'w') as f3:
        f3.writelines(lines3)
    '''
