#!/usr/bin/env python3
import os 
import sys
#--------------------------------------------------------
# INPUTS
#--------------------------------------------------------
runs=['s1','c1','d1','b1','s2','c2','d2','b2']             #the model runs, there must be corresponding particle.bp file in bp/
if len(sys.argv)>1: runs=sys.argv[1:]                      #allow you to define the runs with command line arguments

#--------------------------------------------------------
# remove the sybomlic link of nc files
#--------------------------------------------------------
for brun in runs:
    run='run'+brun
    print(f'\n=================== {run} =====================')
    if os.path.exists(run):
        cmd=f'''cd {run}; rm *.nc'''
        print(cmd); os.system(cmd)
    
