'''==================================================================
                  THIS FUNCTION RUNS THE MODEL FORWARD
=================================================================='''

import numpy             as np
import pandas            as pd
import matplotlib.pyplot as plt

from functions.f_comp_model import comp_model
from functions.f_run_model  import run_model

import variables  as ns
import data_mod.mdata as mdata

def feed4ward(params, df_file):

    # COMPILE THE MODEL:
    # ===============================================================
    comp_model('MOSAIC.v0.DEV_OXY',params)

    # RUN THE MODEL AND GET OUTPUTS:
    # ===============================================================
    outs = run_model(ns)

    df_AEOUT = outs['df_AEOUT']
    df_SOA   = outs['df_SOA']
    df_SIZE  = outs['df_SIZE']
    df_AENUM = outs['df_AENUM']
    df_VOC   = outs['df_VOC']
    df_NCONC = outs['df_NCONC']
    df_O2C   = outs['df_O2C']

    # SAVE DATA:
    # ===============================================================
    with pd.ExcelWriter('outputs/model_outs.xlsx') as xls:
        df_SOA.to_excel(xls,sheet_name='df_SOA')
        xls.close()

    df = pd.DataFrame({'time': df_SOA['time'].values,
                       'SOA': df_SOA['SOA'].values,
                       'O2C': df_O2C[3].values})
    df.to_csv(df_file)
    
    # MAKE PLOTS:
    # ===============================================================
    plt.figure()
    ax = plt.gca()

    tt = df_SOA['time'].values
    yy = df_SOA['SOA'].values
    zz = df_SOA['OLIG'].values

    ttt = mdata.x_obs
    yyy = mdata.y_obs

    print('SOA = ',yy[-1])
    print('OLG = ',zz[-1])
    
    ax.plot(tt,yy,color='r',label='SOA')
    ax.plot(tt,zz,color='r',linestyle='--',label='DIMER')
    ax.plot(ttt,yyy,color='grey',marker='o',linestyle='None')

    ax.set_xlabel('Time (hours)')
    ax.set_ylabel('SOA ($\mu$g m$^{-3}$)')

    ax.set_xlim((0.,None))
    ax.set_ylim((0.,None))

    ax.legend(loc=2)
    
    plt.savefig('outputs/soa.png',bbox_inches='tight')
    
    # MAKE PLOTS:
    # ===============================================================
    plt.figure()
    ax = plt.gca()

    tt = df_SIZE['time'].values

    for col in df_SIZE.columns[1:]:
        
        yy = df_SIZE[col].values

        ax.plot(tt,yy,label=col)

    ax.set_xlabel('Time (hours)')
    ax.set_ylabel('')

    ax.set_xlim((0.,None))
    ax.set_ylim((0.,None))

    ax.legend()
    
    # MAKE PLOTS:
    # ===============================================================
    plt.figure()
    ax = plt.gca()

    time = df_SIZE['time'].values
    
    for i in range(np.shape(df_SIZE)[0]):
        
        if i in [0,600,1200,1800]:
        
            xx = df_SIZE.values[i,1:]
            yy = df_AENUM.values[i,1:]

            ax.plot(xx,yy,label='t = %.1f hrs'%time[i])

    ax.set_xscale('log')

    ax.legend(loc=1)

    ax.set_xlabel('Size (nm)')
    ax.set_ylabel('dNdLogDp (cm$^{-3}$)')

    ax.ticklabel_format(axis='y',style='sci',scilimits=(0.,0.),useMathText=True)
    
    plt.savefig('outputs/ndist.png',bbox_inches='tight')

    # MAKE PLOTS:
    # ===============================================================
    plt.figure()
    ax = plt.gca()
    
    tt = df_VOC[0]
    yy = df_VOC[1]
    
    ax.plot(tt,yy)

    ax.set_xlim((0.,None))
    ax.set_ylim((0.,None))

    ax.set_xlabel('Time (hours)')
    ax.set_ylabel('VOC (ppb)')
    plt.close()

    # MAKE PLOTS:
    # ===============================================================
    plt.figure()
    ax = plt.gca()
    
    tt = df_NCONC[0]
    yy = df_NCONC[1]

    print('NCONC = ',yy.values[-1])
    
    ax.plot(tt,yy)

    ax.set_xlim((0.,None))
    ax.set_ylim((0.,None))

    ax.set_xlabel('Time (hours)')
    ax.set_ylabel('Number Conc. (cm$^{-3}$)')
    
    ax.ticklabel_format(axis='y',style='sci',scilimits=(0.,0.),useMathText=True)
    plt.savefig('outputs/nconc.png',bbox_inches='tight')
    
    # MAKE PLOTS:
    # ===============================================================
    plt.figure()
    ax = plt.gca()
    
    tt = mdata.x_obs
    yy = mdata.z_obs

    ax.plot(tt,yy,color='k',marker='o',linestyle='None')
    
    tt = df_O2C[0].values
    yy = df_O2C[3].values

    print(yy)
    ax.plot(tt,yy)

    ax.set_xlim((0.,None))
    ax.set_ylim((0.,1.0))
    
    ax.set_xlabel('Elapsed Time (hrs)')
    ax.set_ylabel('O:C')
    
    plt.show()


    df = pd.DataFrame({'time': df_SOA['time'].values,
                       'SOA': df_SOA['SOA'].values,
                       'O2C': yy})
    df.to_csv(df_file)
