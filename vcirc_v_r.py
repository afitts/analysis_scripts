import numpy as np
import sys 
import glob
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as co
import pylab
import pygadgetreader as pg
import scipy.integrate
import time
import scipy.stats as stats
import scipy.special as sp
import pandas as pd
import scipy.optimize as opt
from matplotlib import rcParams
from matplotlib import rc

def load_vcirc(pathname,hnum,res,ver):
	G = 4.3e-6 #in kpc/M_sun (km/s)^2 
	bincount = 50
	if hnum == '10q' or hnum == '10v':
		i = 600
	else:
		i = 184
	hdf = pd.HDFStore("%s/halo%s%s_giz%s_snap%s.h5"%(pathname,hnum,res,ver,i))
	mass = hdf['particles/dm']['mass'].as_matrix()
	gmass = hdf['particles/gas']['mass'].as_matrix()
	smass = hdf['particles/star']['mass'].as_matrix()
	r = hdf['particles/dm']['r'].as_matrix()*1000
	gr = hdf['particles/gas']['r'].as_matrix()*1000
	sr = hdf['particles/star']['r'].as_matrix()*1000
	mass = np.concatenate((mass,gmass,smass))
	r = np.concatenate((r,gr,sr))
	mass,bin_edges = np.histogram(r,bins=np.logspace(-2,np.log10(6),bincount),weights=mass)
	bin_mids = 10**(np.log10(bin_edges[1:])-np.log10(bin_edges[2]/bin_edges[1])/2)
	cmass = np.cumsum(mass)
	vcirc = np.sqrt(G*cmass/bin_mids)
	hdf.close()
	return bin_edges[:bincount-1],bin_mids,bin_edges[1:],vcirc

if __name__ = "__main__":

	hnum =['11707','32503','12596','007','848','796','20192','20910','897','948','1016','32257','10q','10v','1084']
	res = ['_13','_13','_13','_11','_13','_13','_13','_13','_13','_13','_13','_13','_11','_11','_13']
	#res = ['','','','','','','','','','','','','_11','_11','']
	ver = ['5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16']

	for w in np.arange(len(hnum)):
		date = time.strftime("%m_%d_%Y")

		pathname = '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/analysis/dataframes'%(hnum[w],res[w],ver[w])

		lbin,mbin,rbin,vcirc = load_vcirc(pathname,hnum[w],res[w],ver[w])
		np.savetxt('halo%s_vcirc_v_r.txt'%(hnum[w]),np.column_stack((lbin,mbin,rbin,vcirc)),header = '(0) L side of bin	(1) Mid of bin	(2) R side of bin	(3) Vcirc')
	return True




