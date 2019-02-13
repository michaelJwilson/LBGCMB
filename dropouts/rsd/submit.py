import os
import sys

'''
##  Euclid run. 
print('\n\nEuclid run.\n\n')

for zmin in [1.4, 1.6, 1.8]:
    ##  Must set _bz, nz in code:                                                                                                                
    ## 
    ##  _bz = euclid_bz
    ##  nz = euclid_nz

    survey = 'euclid'
    fsky   = 15000. / 41253. 

    zmax   = zmin + 0.2

    os.system('python rsd_fish.py %s %lf %lf %lf' % (survey, fsky, zmin, zmax))
'''

print('\n\nBEAST run.\n\n')                                                                                                                              

for [zmin, zmax] in [[2.0, 2.5], [2.5, 3.5], [3.5, 4.25], [4.25, 5.2]]:                                                                                  
    ##  Must set _bz, nz in code:
    ## 
    ##  _bz =  beast_bz                                                                                                                                
    ##  nz  =  beast_nz  

    survey  = 'beast'                                                                                                                                   
    fsky    =  14000. / 41253.                                                                                                                          
                                                                                                                                         
    os.system('python rsd_fish.py %s %lf %lf %lf' % (survey, fsky, zmin, zmax))

print('\n\nDone.\n\n')
