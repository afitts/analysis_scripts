import numpy as np
import sys 
import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.cm as cm
import matplotlib.colors as co
import pylab
import time
import pandas as pd
import pygadgetreader as pg
from mikecm import mikecm

hnum = sys.argv[1]
snum = np.int(sys.argv[2])

pathname = '/nobackup/afitts/Gadget-2.0.7/production/mfm%s_13_giz5_12_16_raw_output'%(hnum)

hdf = pd.HDFStore('/nobackup/afitts/Gadget-2.0.7/production/mfm%s_13_giz5_12_16_raw_output/analysis/dataframes/halo%s_13_giz5_12_16_snap%03d.h5'%(hnum,hnum,snum))
xcen = np.float(hdf['props']['halox'])*1000
ycen = np.float(hdf['props']['haloy'])*1000
zcen = np.float(hdf['props']['haloz'])*1000
rvir = np.float(hdf['props']['rvir'])*1000
red = np.float(hdf['props']['redshift'])
drdf = hdf['particles/dm']['r'].as_matrix()*1000
grdf = hdf['particles/gas']['r'].as_matrix()*1000
dmdf = hdf['particles/dm']['mass'].as_matrix()
gmdf = hdf['particles/gas']['mass'].as_matrix()
dpos = pg.readsnap('%s/snapdir_%03d/snapshot_%03d.0.hdf5'%(pathname,snum,snum),'pos','dm')
dm = pg.readsnap('%s/snapdir_%03d/snapshot_%03d.0.hdf5'%(pathname,snum,snum),'mass','dm')
gpos = pg.readsnap('%s/snapdir_%03d/snapshot_%03d.0.hdf5'%(pathname,snum,snum),'pos','gas')
gm = pg.readsnap('%s/snapdir_%03d/snapshot_%03d.0.hdf5'%(pathname,snum,snum),'mass','gas')

dx = dpos[:,0]-xcen
dy = dpos[:,1]-ycen
dz = dpos[:,2]-zcen
gx = gpos[:,0]-xcen
gy = gpos[:,1]-ycen
gz = gpos[:,2]-zcen
dr = np.sqrt(dx**2+dy**2+dz**2)/.71
gr = np.sqrt(gx**2+gy**2+gz**2)/.71
dsc = mikecm(fname = dpos[dr<rvir],nofile=True,pmass=dm[dr<rvir])
ddiff = np.sqrt((dpos[:,0]-dsc[0])**2+(dpos[:,1]-dsc[1])**2+(dpos[:,2]-dsc[2])**2)/.71
gdiff = np.sqrt((gpos[:,0]-dsc[0])**2+(gpos[:,1]-dsc[1])**2+(gpos[:,2]-dsc[2])**2)/.71
bfrac = sum(gm[gr<rvir*1])/(sum(dm[dr<rvir*1])+sum(gm[gr<rvir*1]))
bfrac1 = sum(gm[gdiff<rvir])/sum(dm[ddiff<rvir])
bfracdf = sum(gmdf[grdf<rvir])/sum(dmdf[drdf<rvir])
print 'OLD CENTER:',xcen,ycen,zcen
print 'NEW CENTER:',dsc
print 'bfrac from snap is %.3e for halo %s at z=%f, Mvir = %.3e, Mgas = %.3e, rvir = %.3e (comoving), %.3e (physical'%(bfrac,hnum,red,sum(dm[dr<rvir])*1e10,sum(gm[gr<rvir])*1e10,rvir,rvir/(red+1))
print 'bfrac from df is %.3e for halo %s at z=%f, Mvir = %.3e, Mgas = %.3e'%(bfracdf,hnum,red,sum(dmdf[drdf<rvir]),sum(gmdf[grdf<rvir]))
print 'bfrac from snap new cen is %.3e for halo %s at z=%f, Mvir = %.3e, Mgas = %.3e'%(bfrac1,hnum,red,sum(dm[ddiff<rvir])*1e10,sum(gm[gdiff<rvir])*1e10)

