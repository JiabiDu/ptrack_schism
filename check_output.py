#!/usr/bin/env python3
from pylib import read_schism_hgrid, ReadNC, plot, show, arange
#read data
C=ReadNC('out.nc')
plon=C.lon.val
plat=C.lat.val

#read schism hgrid
gd=read_schism_hgrid('hgrid.ll')
gd.plot_bnd()

#plot the first 10 particles
for i in arange(10):
    plot(plon[:,i],plat[:,i])
show()

