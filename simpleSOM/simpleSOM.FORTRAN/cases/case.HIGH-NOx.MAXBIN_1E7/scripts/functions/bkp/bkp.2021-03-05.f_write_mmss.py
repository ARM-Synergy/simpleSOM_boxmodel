'''=====================================================================
         THIS FUNCTION WRITES THE MECHANISM AND SPECIES FILES
====================================================================='''

import numpy  as np
import pandas as pd

def write_mmss(path,params):

    # GET PARAMETERS:
    # ==================================================================
    mfrag  = params[0]
    f_LOSS = np.exp(-params[1])
    f_HOM  = params[2]
    dLVP   = params[3]
    
    PO1 = params[4]
    PO2 = params[5]
    PO3 = params[6]
    PO4 = params[7]

    POs = [i/sum([PO1,PO2,PO3,PO4])*(1. - f_HOM) for i in [PO1,PO2,PO3,PO4]]

    PO1 = POs[0]
    PO2 = POs[1]
    PO3 = POs[2]
    PO4 = POs[3]
    
    # VOLATILITY DIST. OF ADDING 1, 2, 3 or 4 OXYGENS:
    # ==================================================================
    SIGMAX = 1.

    DIST_PO1 = np.array([np.exp(-(1.*dLVP-(i+1))**2./SIGMAX) for i in range(2)])
    DIST_PO2 = np.array([np.exp(-(2.*dLVP-(i+2))**2./SIGMAX) for i in range(3)])
    DIST_PO3 = np.array([np.exp(-(3.*dLVP-(i+3))**2./SIGMAX) for i in range(4)])
    DIST_PO4 = np.array([np.exp(-(4.*dLVP-(i+4))**2./SIGMAX) for i in range(5)])
    
    DIST_PO1 = DIST_PO1/np.sum(DIST_PO1)
    DIST_PO2 = DIST_PO2/np.sum(DIST_PO2)
    DIST_PO3 = DIST_PO3/np.sum(DIST_PO3)
    DIST_PO4 = DIST_PO4/np.sum(DIST_PO4)
    
    # GENERIC FUNCTIONALIZATION DISTRIBUTION:
    # ==================================================================
    nFRAG = 2
    nFUNC = 8

    CHEMDIST_SOM = np.zeros(nFUNC+1+nFRAG)

    CHEMDIST_SOM[nFRAG]   = 0.
    CHEMDIST_SOM[nFRAG+1] = PO1*DIST_PO1[0]
    CHEMDIST_SOM[nFRAG+2] = PO1*DIST_PO1[1] + PO2*DIST_PO2[0]
    CHEMDIST_SOM[nFRAG+3] = PO2*DIST_PO2[1] + PO3*DIST_PO3[0]
    CHEMDIST_SOM[nFRAG+4] = PO2*DIST_PO2[2] + PO3*DIST_PO3[1] + PO4*DIST_PO4[0]
    CHEMDIST_SOM[nFRAG+5] = PO3*DIST_PO3[2] + PO4*DIST_PO4[1]
    CHEMDIST_SOM[nFRAG+6] = PO3*DIST_PO3[3] + PO4*DIST_PO4[2]
    CHEMDIST_SOM[nFRAG+7] = PO4*DIST_PO4[3]
    CHEMDIST_SOM[nFRAG+8] = PO4*DIST_PO4[4]

    # FRAGMENTATION PROBABILITY:
    # ==================================================================
    maxCSAT = 9
    minCSAT = -3
    nCOMP   = maxCSAT - minCSAT + 1

    # SATURATION CONC.:
    j_CSATLOG = np.arange(maxCSAT,minCSAT-1,-1)

    # FRAG. PROB.:
    j_PFRAG = 1. - np.exp(mfrag*(j_CSATLOG - maxCSAT)/maxCSAT)

    # PROB. DISTRIBUTION FOR ALL SPECIES:
    # ==================================================================
    CHEMDIST_SOM_MATRIX = np.zeros((nFUNC+1+nFRAG,nCOMP))
    
    for j in range(nCOMP):
        
        for jj in range(nFRAG+1,nFUNC+1+nFRAG):
            CHEMDIST_SOM_MATRIX[jj,j] = CHEMDIST_SOM[jj]*(1. - j_PFRAG[j])
            
        for jj in range(nFRAG):
            CHEMDIST_SOM_MATRIX[jj,j] = j_PFRAG[j]/float(nFRAG)*(1. - f_LOSS)
    
    # GAS-PHASE CHEMISTRY KERNEL:
    # ==================================================================
    jj_CHEMKERNEL_SOM = np.zeros((nCOMP,nCOMP))

    for j in range(nCOMP):

        SLICE_SOM = CHEMDIST_SOM_MATRIX[:,j]

        for k in range(nFRAG):
            if ((j-k-1)<0):
                pass
            else:
                jj_CHEMKERNEL_SOM[j-k-1,j] = SLICE_SOM[nFRAG-k-1]

        for k in range(nFUNC):
            if ((j+k+1)<nCOMP-1):
                jj_CHEMKERNEL_SOM[j+k+1,j] = SLICE_SOM[nFRAG+k+1]
            else:
                jj_CHEMKERNEL_SOM[nCOMP-1,j] += SLICE_SOM[nFRAG+k+1]

    # ADD HOM FORMATION:
    jj_CHEMKERNEL_SOM[nCOMP-1,:] += f_HOM * (1 - j_PFRAG)
                
    
    
    
    
    OFFSET = 17
    
    # FUNCTIONALIZATION:
    # ==================================================================
    maxAddOxygen = 4

    width_nO = 1 # Currently set to unity

    frac   = np.zeros((maxAddOxygen,maxAddOxygen*2))
    prob   = np.zeros((maxAddOxygen,maxAddOxygen*2))

    # Number of fragments. This is a user input
    numFragments = 2

    for n in range(maxAddOxygen):
        for m in range(n, 2*n+2):
            frac[n][m] = np.exp(-1 * ((n+1)*dLVP - (m+1))**2 / width_nO) 
    
        prob[n] = frac[n] / frac[n].sum() * pfs[n]

    # Sum rows to get the probability for each volatility bin
    prob_prod = sum(prob)
    
    # Range of volatility space
    maxLogCStar = 9
    minLogCStar = -3

    maxCSAT = maxLogCStar
    minCSAT = minLogCStar

    nSPEC = maxLogCStar - minLogCStar + 1

    
    
    # Array of C*
    logCStar25C = np.array(range(minLogCStar, maxLogCStar+1))

    j_CSATLOG = logCStar25C 
    
    # Initialize products by bin
    prob_prod_bins = pd.DataFrame(np.zeros((maxAddOxygen*2, maxLogCStar-minLogCStar+1)))
    prob_frag_bin_temp = pd.DataFrame(np.zeros((1,logCStar25C.size)))

    # Probability of fragmentation for each volatility bin
    pFrag_bins = 1 - np.exp(mfrag * (logCStar25C-maxLogCStar)/maxLogCStar) # Checked against SOM

    # Probability of functionalization for each volatility bin
    ## One would expect that the part of functionalization that produces
    ## ELVOCs would be subtracted here then considered separately, 
    ## but this is not how it is implemented in simpleSOM.
    ## We follow simpleSOM here for consistency.
    pFunc_bins = (1 - pFrag_bins)
    # Probability of fragmentation for each volatility bin
    prob_frag_bins = np.zeros((logCStar25C.size))
    prob_zero_bins = np.zeros((1,logCStar25C.size))

    # Loop through each volatility bin to calculate the func/frag change
    # IMPORTANT NOTE: prob_frag_bin_temp does NOT sum to 1
    # frac_mass_loss assumes that at each fragmentation step, mass is 
    # completely lost from the system, and is no longer tracked.
    for i in list(range(logCStar25C.size)):
        prob_prod_bins[i] = pFunc_bins[i] * prob_prod # Checked against SOM
        prob_frag_bin_temp[i] = pFrag_bins[i] * (1-fLOSS) / numFragments
        
    # Repeat the array the same number of times as numFragments, because there are 
    # exactly that many fragments in the array
    prob_frag_bins = np.tile(prob_frag_bin_temp, (numFragments,1))

    # Combine all the reaction results together
    full_reaction_bins = pd.DataFrame(np.vstack([prob_frag_bins, prob_zero_bins, prob_prod_bins]))
    df_full_reaction_bins = pd.DataFrame(full_reaction_bins)

    # Change the column and row labels on the reaction dataframe
    df_full_reaction_bins.columns = np.array(range(minLogCStar, maxLogCStar+1))
    df_full_reaction_bins.index = np.flip(np.array(range(-maxAddOxygen*2, numFragments+1)), axis=0)

    
    # HOMS FORMATION YIELD:
    # ==================================================================
    j_fHOMS = fHOMS * pFunc_bins
    
    # OH-OXIDATION REACTION RATES:
    # ==================================================================
    kxA =  1.56341E-13
    kxB = -7.12119E-13
    kxC = -8.22305E-12
    kxD = -5.62410E-13
    kxE = -5.69473E-13
    kxF =  6.62928E-11

    j_kRXN = (kxA*dLVP + kxD)*(j_CSATLOG**2) + (kxB*dLVP + kxE)*j_CSATLOG + (kxC*dLVP + kxF)

    
    # ==================================================================
    def rxn_prod_string(reactant_slice, reactant_cstar):
        prod_string = ""
        # Loop through the reactions for 1 volatility bin
        for j in list(range(-maxAddOxygen*2, numFragments+1)):
            # This is to enforce upper and lower limits of volatility bins
            if j+reactant_cstar < minLogCStar:
                prod_bin = minLogCStar
            elif j+reactant_cstar > maxLogCStar:
                prod_bin = maxLogCStar
            else:
                prod_bin = j + reactant_cstar
            # Write each reaction to a string, only if the coefficient is > 0
            if reactant_slice[j] > 0:
                #prod_string = ( prod_string + " + " + str(sigfig.sci_notation(reactant_slice[j], 3)) 
                #                + " C" + str(prod_bin).replace("-","n"))
                
                prod_string = ( prod_string + " + " + '%.2e'%reactant_slice[j] 
                                + " C" + str(prod_bin).replace("-","n"))
            # Adds line breaks every 2 indices, 
            # because Fortran has trouble with lines going over 80 characters (or something)
                if j % 2 == 1:
                    prod_string = prod_string + '\n'
    # Return the reaction string for this CStar
        return prod_string

    # CREATE THE SPECIES AND MECHANISM LINES:
    # ==================================================================
    lines_mm = ''
    lines_ss = ''
    
    for i,CSAT in enumerate(np.arange(minCSAT,maxCSAT+1)):

        # STRINGS:
        str1 = str(i + OFFSET)
        str2 = str(CSAT).replace('-','n')
        str3 = str(minCSAT).replace('-','n')

        str4 = rxn_prod_string(df_full_reaction_bins[CSAT], CSAT).replace('e','d')
        str5 = ('%.2e'%(j_fHOMS[i])).replace('e','d')
        str6 = ('%.2e'%(j_kRXN[i])).replace('e','d')
        
        # FOR SPECIES:
        lines_ss = lines_ss + '%s iC%s C%s   = IGNORE ;\n'%(str1,str2,str2)
        
        # FOR MECHANISM:
        lines_mm = lines_mm + '{%s} OH + C%s =%s + %s C%s : %s ; \n\n'%(str1,str2,str4,str5,str3,str6)
        
        
    lines_ss = lines_ss + str(nSPEC + OFFSET) + ' iOH OH  = IGNORE ;\n'
    
    
    # WRITE THE SPECIES AND MECHANISM FILES:
    # ==================================================================
    wpath = '%s/source_ftn/gas/OrganicsMech/'%path
    
    with open('%s/simplesom.ss'%wpath,'w') as f1:
        
        f1.write("*dynamic* %i\n"%(nSPEC+1))
        f1.write(lines_ss)
        f1.write("*end*\n\n")

        f1.write("*fixed* 0\n")
        f1.write("*end*\n\n")

        f1.write("*atoms* 0\n")
        f1.write("*end*")
        
        f1.close()
    
    with open('%s/simplesom.mm'%wpath,'w') as f2:

        f2.write("*Reactions* %i\n\n"%nSPEC)
        f2.write(lines_mm)
        f2.write("*end*\n\n")

        f2.close()
