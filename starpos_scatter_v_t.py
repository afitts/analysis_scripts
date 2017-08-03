# IPython log file
#Import necessary modules
#import matplotlib
#matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LogNorm
import numpy as np
import pygadgetreader as pg

#Turn on interactive plots
#plt.ion()

def loader(hnum,ver,res,snap,i):
	try:
		#Load up the halo center and virial radius
		hdf = pd.HDFStore('/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/analysis/dataframes/halo%s%s_giz%s_snap%03d.h5'%(hnum,res,ver,hnum,res,ver,snap))
		hx = np.float(hdf['props']['halox'])
		hy = np.float(hdf['props']['haloy'])
		hz = np.float(hdf['props']['haloz'])
		rvir = np.float(hdf['props']['rvir'])
		red = np.float(hdf['props']['redshift'])
		smass = hdf['particles/star']['mass'].as_matrix()
		sx = hdf['particles/star']['x'].as_matrix()
		sy = hdf['particles/star']['y'].as_matrix()
		sz = hdf['particles/star']['z'].as_matrix()
		spos = hdf['particles/star']['r'].as_matrix()
		time = np.float(hdf['props']['time'])

		#Plot it!
		fig = plt.figure()
		halopic = fig.add_subplot(111)
		halopic.scatter(sx,sy)
		halopic.add_patch(patches.Circle((hx, hy),0.1*rvir,fill=False,color='r',linewidth = 6))
		halopic.set_xlabel('x')
		halopic.set_ylabel('y')

		halopic.set_xlim(hx-.2*rvir,hx+.2*rvir)
		halopic.set_ylim(hy-.2*rvir,hy+.2*rvir)
		halopic.set_title(r'$M_\star=%.2e,\:\:\rm Time:\:%.2e\:Gyr$'%(sum(smass[spos<.1*rvir]),time))
		hdf.close()
		fig.savefig('Halo%s%s_starpos_xy/Halo%s%s_starpos_xy_%03d.png'%(hnum,res,hnum,res,i),transparent = True)
	except Exception,f:
		print f,i
		hdf.close()
	return True

hnum =['948']#'11707','32503','12596','007','848','796','20192','20910','897','948','1016','32257','10q','10v','1084']
res = ['_13','_13','_13','_11','_13','_13','_13','_13','_13','_13','_13','_13','_11','_11','_13']
ver = ['11_13','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16']
dmover = ['11_13','11_13','11_13','11_13','11_13','11_13','5_12_16','5_12_16','11_13','11_13','11_13','11_13','5_12_16','5_12_16','5_12_16']
snum = [184]#[184,184,184,184,184,184,184,184,184,184,184,184,600,600,184]


for w in np.arange(len(hnum)):
	snap = np.linspace(0,snum[w],snum[w]+1)
	for i in np.arange(len(snap)):
  		loader(hnum[w],ver[w],res[w],snap[i],i)



