import numpy as np
import sys
import glob
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as co
import pylab
import pygadgetreader as pg
import scipy.integrate
import time
import pandas as pd

def HI_mass(hnum,res,snum):
	pathname = '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz5_12_16_raw_output'%(hnum,res)
	pos = pg.readsnap('%s/snapdir_%03d/snapshot_%03d.0.hdf5'%(pathname,snum,snum),'pos','gas')
	gx = pos[:,0]
	gy = pos[:,1]
	gz = pos[:,2]
	nh = pg.readsnap('%s/snapdir_%03d/snapshot_%03d.0.hdf5'%(pathname,snum,snum),'nh','gas')
	gm = pg.readsnap('%s/snapdir_%03d/snapshot_%03d.0.hdf5'%(pathname,snum,snum),'mass','gas')
	hdf = pd.HDFStore('%s/analysis/dataframes/halo%s%s_giz5_12_16_snap%03d.h5'%(pathname,hnum,res,snum))
	xcen = np.float(hdf['props']['halox'])*1000
	ycen = np.float(hdf['props']['haloy'])*1000
	zcen = np.float(hdf['props']['haloz'])*1000
	rvir =  np.float(hdf['props']['rvir'])*1000
	gr = np.sqrt((gx-xcen)**2+(gy-ycen)**2+(gz-zcen)**2)
	himass = sum(gm[gr<0.1*rvir]*nh[gr<0.1*rvir])*1e10
	hdf.close()
	print 'HI mass for halo %s is %.3e'%(hnum,himass)

if __name__ == "__main__":
        hnum =['1084']#11707','32503','12596','007','848','796','20192','20910','897','948','1016','32257','10q','10v']
        res = ['_13','_13','_13','_11','_13','_13','_13','_13','_13','_13','_13','_13','_11','_11','_13']
        ver = ['5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16']
        snap = [184,184,184,184,184,184,184,184,184,184,184,184,600,600]
        date = time.strftime("%m_%d_%Y")
        for i in np.arange(len(hnum)):
		HI_mass(hnum[i],res[i],snap[i])
