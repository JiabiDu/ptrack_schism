#!/usr/bin/env python3
from pylib import *
#read data
C=ReadNC('out.nc')
plon=C.lon.val
plat=C.lat.val

#read schism hgrid
gd=read_schism_hgrid('hgrid.ll')
gd.plot_bnd()

#plot the first 10 particles
for i arange(10):
    plot(plon[i,:],plat[i,:])
show()

