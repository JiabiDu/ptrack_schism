#!/usr/bin/env python3
#sample call: ./run.multiple
import os 
import sys
#--------------------------------------------------------
# INPUTS
#--------------------------------------------------------
settlings=[0,1,2,3]
locs={'g':'galveston','s':'sanJacinto','c':'clearLake','d':'dickingson','b':'buffalo'} #short and full name of release locations
check_run=True #check if run exists

#--------------------------------------------------------
# creat run directories and submit jobs
#--------------------------------------------------------
for loc in locs:
    loc_name=locs[loc] #the full name
    if not os.path.exists(f'Inputs/particle.bp_{loc_name}'): print(f'Inputs/particle.bp_{loc_name} not eists'); continue
    for settling in settlings:
        run='run'+loc+'_settling{settling}'
        print(f'\n=================== {run} =====================')
        if os.path.exists(run) and check_run:
            x=input(f'{run} already exists; do you want to overwrite?\nType Y/y to confirm, other keys to skip this run: ')
            if x not in ['Y','y']: continue
    
        #following are a bunch of linux command to set-up a run and then submit the job
        cmd=f'''rm -rf {run}; mkdir {run}; cd {run}; ln -sf ../Inputs/particle.bp_{loc_name} particle.bp; 
             ln -sf ../base/ptrack3.WW .; ln -sf ../base/hgrid.ll .; ln -sf ../base/vgrid.in .;
             cp ../base/param.in .; sed -i 's/settling_velocity=0/settling_velocity={settling}/g' param.in; 
             cp ../base/run_ptrack .; sed -i 's/job-name=ptrack/job-name=ptrack-{run}/g' run_ptrack; qsub run_ptrack'''
        print(cmd)
        
        #run the command on your front node
        os.system(cmd)
