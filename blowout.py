import pygadgetreader as pg
import pandas as pd
import numpy as np
from pyd import mikecm
import sys

i = int(sys.argv[1])
hdf = pd.HDFStore('/nobackup/afitts/Gadget-2.0.7/production/mfm1016SNe_13_giz5_12_16_raw_output/analysis/dataframes/halo1016SNe_13_giz5_12_16_snap%03d.h5'%(i))
red = np.float(hdf['props']['redshift'])
pos = pg.readsnap('/nobackupp8/afitts/Gadget-2.0.7/production/mfm1016SNe_13_giz5_12_16_raw_output/snapdir_%03d/snapshot_%03d.0.hdf5'%(i,i),'pos','gas')/1000
vel = pg.readsnap('/nobackupp8/afitts/Gadget-2.0.7/production/mfm1016SNe_13_giz5_12_16_raw_output/snapdir_%03d/snapshot_%03d.0.hdf5'%(i,i),'vel','gas')*np.sqrt(red+1)
mass = pg.readsnap('/nobackupp8/afitts/Gadget-2.0.7/production/mfm1016SNe_13_giz5_12_16_raw_output/snapdir_%03d/snapshot_%03d.0.hdf5'%(i,i),'mass','gas')

hvx = np.float(hdf['props']['halo_vx'])*1000
hvy = np.float(hdf['props']['halo_vy'])*1000
hvz = np.float(hdf['props']['halo_vz'])*1000
rvir = np.float(hdf['props']['rvir'])*1000*np.sqrt(red+1)

vx = vel[:,0]-hvx
vy = vel[:,1]-hvy
vz = vel[:,2]-hvz

dm = hdf['particles/dm']['mass'].as_matrix()
dx = hdf['particles/dm']['x'].as_matrix()
dy = hdf['particles/dm']['y'].as_matrix()
dz = hdf['particles/dm']['z'].as_matrix()
dp = np.column_stack((dx,dy,dz))
smm = hdf['particles/star']['mass'].as_matrix()
sx = hdf['particles/star']['x'].as_matrix()
sy = hdf['particles/star']['y'].as_matrix()
sz = hdf['particles/star']['z'].as_matrix()
sp = np.column_stack((sx,sy,sz))
dsp = np.vstack((dp,sp))
dsm = np.append(dm,smm)
dsc = mikecm(fname = dsp,nofile=True,pmass=dsm)
dist = np.sqrt((pos[:,0]-dsc[0])**2+(pos[:,1]-dsc[1])**2+(pos[:,2]-dsc[2])**2)/.71*1000*np.sqrt(red+1)
pos = pos/(red+1)/.71
vr = (pos[:,0]*vx+pos[:,1]*vy+pos[:,2]*vz)/dist
v = np.sqrt(vx**2+vy**2+vz**2)
ke = 0.5*mass*v**2*1e10*1.99e33*(1000*100)**2
