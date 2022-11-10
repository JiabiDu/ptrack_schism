#!/usr/bin/env python3
''' This code is to compile the ptrack3.f90 and add tag to the compiled code '''
import os 
#compile command
cmd=f'''cd ~/ptrack_schism; ifort -mcmodel=medium -CB -assume byterecl -O2 -o ptrack3.WW ptrack3.f90 compute_zcor.f90 schism_geometry.f90 -I$NETCDF/include -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff; rm *.mod fort.11'''
print(cmd); os.system(cmd)

#add git tag
tag=os.popen('git log').read().split('\n')[0].split()[1][:8]
cmd='mv ~/ptrack_schism/ptrack3.WW ./ptrack3.WW.{}; ln -sf ptrack3.WW.{} ptrack3.WW; ./ptrack3.WW'.format(tag,tag)
print(cmd); os.system(cmd)



