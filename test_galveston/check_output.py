#!/usr/bin/env python3
from pylib import read_schism_hgrid, ReadNC, plot, show, arange,xlim,ylim
#read data
C=ReadNC('out.nc')
plon=C.lon.val
plat=C.lat.val

#read schism hgrid
gd=read_schism_hgrid('/scratch/user/jdu/GoM/run19q3_20yr/hgrid.ll')
gd.plot_bnd()

#plot the first 10 particles
for i in arange(0,500,10):
    plot(plon[:,i],plat[:,i])
xlim([-95.4,-94.4])
ylim([29,30])
show()

