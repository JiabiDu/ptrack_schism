#!/usr/bin/env python3
''' This code is to compile the ptrack3.f90 and add tag to the compiled code '''
import os 
import sys
home='/home/jdu'

#compile command
cmd=f'''cd ~/ptrack_schism; ifort -mcmodel=medium -CB -assume byterecl -O2 -o ptrack3.WW ptrack3.f90 compute_zcor.f90 schism_geometry.f90 -I$NETCDF/include -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff'''
print(cmd); os.system(cmd)
if not os.path.exists(f'{home}/ptrack_schism/ptrack3.WW'): sys.exit('no success in compiling')
#add git tag
tag=os.popen('cd ~/ptrack_schism; git log').read().split('\n')[0].split()[1][:8]
host=os.popen('echo $HOST').read().split()[0].upper()
cmd=f'mv ~/ptrack_schism/ptrack3.WW ./ptrack3.WW.{host}.{tag}; ln -sf ptrack3.WW.{host}.{tag} ptrack3.WW;' # ./ptrack3.WW'
print(cmd); os.system(cmd)


