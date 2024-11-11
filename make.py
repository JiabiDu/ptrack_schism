#!/usr/bin/env python3
''' This code is to compile the ptrack3.f90 and add tag to the compiled code '''
import os 
import sys

#compile the code
bdir='~/data10/ptrack_schism'
cmd=f'''cd {bdir}; ifort -mcmodel=medium -CB -assume byterecl -O2 -o ptrack3.WW ptrack3.f90 compute_zcor.f90 schism_geometry.f90 -I$NETCDF/include -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff'''
print(cmd); os.system(cmd)
if len(os.popen(f'ls {bdir}/ptrack3.WW').read())==0: sys.exit(f'no success in compiling; no ptrack3.WW in {bdir}')

#add git tag
tag=os.popen(f'cd {bdir}; git log').read().split('\n')[0].split()[1][:8]

#get the host name
host=os.popen('echo $HOST').read().split()[0].upper()

#add git tag and host name to the ptrack3.WW
cmd=f'mv {bdir}/ptrack3.WW ./ptrack3.WW.{host}.{tag}; ln -sf ptrack3.WW.{host}.{tag} ptrack3.WW;' # ./ptrack3.WW'
print(cmd); os.system(cmd)


