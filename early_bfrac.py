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
if hnum == '007' or hnum == '10q' or hnum == '10v':
	res = '_11'
else:
	res = '_13'
snums = np.linspace(5,10,6)
diff = np.zeros(len(snums)-1)
cen = np.zeros((len(snums),3))
rvir = np.zeros((len(snums)))
red = np.zeros((len(snums)))
pathname = '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz5_12_16_raw_output'%(hnum,res)
for i,snum in enumerate(snums):
	hdf = pd.HDFStore('/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz5_12_16_raw_output/analysis/dataframes/halo%s%s_giz5_12_16_snap%03d.h5'%(hnum,res,hnum,res,snum))
	cen[i,0] = np.float(hdf['props']['halox'])*1000
	cen[i,1] = np.float(hdf['props']['haloy'])*1000
	cen[i,2] = np.float(hdf['props']['haloz'])*1000
	rvir[i] = np.float(hdf['props']['rvir'])*1000
	red[i] = np.float(hdf['props']['redshift'])
	if snum != snums[0]:
		diff[i-1] = np.sqrt(sum(cen[i]-cen[i-1])**2)/.71


