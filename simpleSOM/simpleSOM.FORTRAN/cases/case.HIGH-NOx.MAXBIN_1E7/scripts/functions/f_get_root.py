'''=========================================================================
      THIS FUNCTION FINDS THE ROOT DIRECTORY OF THE CURRENT EXPERIMENT
========================================================================='''

import os

def get_root():
    
    i = 0
    x = '..'

    while i <= 10:
        
        xx = '/'.join([x]*i)
        pp = './%s/.root'%xx
        
        path = 0

        if os.path.isfile(pp):
            path = os.popen('./%s/.root'%xx).read().strip('\n')
            break
        
        i += 1

    if path == 0:
        exit('ROOT DIRECTORY NOT FOUND!')
    else:
        return path
