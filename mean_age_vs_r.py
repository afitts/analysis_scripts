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
from mikecm import mikecm

class mean_age_vs_r(object):
	"""Plots the radial age gradient of in situ vs external stars (e.g.Middle panel of Fig 7. in El-Badry et al. 2016).
	Requires starmerger files (list of all stars and whether or not they participated in a merger, from star_and_halo_particleplot.py)"""
	def __init__(self):
		self.fig = plt.figure()
		self.sub = self.fig.add_subplot(111)
		self.sub.set_xlabel(r'$\mathrm{Radius}$ $(\mathrm{kpc})$',fontsize=20)
		self.sub.set_ylabel(r'$\mathrm{Mean}$ $\mathrm{Age}$ $(\mathrm{Gyr})$',fontsize=20)#Stellar Mass Fraction',fontsize=20)#
		self.sub.set_xlim(0,5)

	def compile_and_plot(self,pathname,hnum,res,ver,snap):
		hdf = pd.HDFStore('%sdataframes/halo%s%s_giz%s_snap%03d.h5'%(pathname,hnum,res,ver,snap))
		red = np.float(hdf['props']['redshift'])
		rvir = np.float(hdf['props']['rvir'])*1000/(red+1)
		xcen = np.float(hdf['props']['halox'])*1000
		ycen = np.float(hdf['props']['haloy'])*1000
		zcen = np.float(hdf['props']['haloz'])*1000
		binz = np.logspace(np.log(.0014),np.log(rvir/10),60,base=np.e)
		x = np.e**(np.log(binz[1:])-np.log(binz[1]/binz[0])/2)
		dlogr = np.log(binz[1]/binz[0])
		dm = hdf['particles/dm']['mass']
		dx = hdf['particles/dm']['x'].as_matrix()*1000/(red+1)
		dy = hdf['particles/dm']['y'].as_matrix()*1000/(red+1)
		dz = hdf['particles/dm']['z'].as_matrix()*1000/(red+1)
		dp = np.column_stack((dx,dy,dz))
		smm = hdf['particles/star']['mass']
		sx = hdf['particles/star']['x'].as_matrix()*1000/(red+1)
		sy = hdf['particles/star']['y'].as_matrix()*1000/(red+1)
		sz = hdf['particles/star']['z'].as_matrix()*1000/(red+1)
		sp = np.column_stack((sx,sy,sz))
		dsp = np.vstack((dp,sp))
		dsm = np.append(dm,smm)
		dsc = mikecm(fname = dsp,nofile=True,pmass=dsm)
		starfile = np.genfromtxt('starmerger_out/Halo%s_starmerger_test.out'%hnum)
		insitu = (starfile[:,1]==4)&(starfile[:,8]==1)
		insitu_y = (starfile[:,1]==4)&(starfile[:,8]==1)&(starfile[:,9]>3.09)
		insitu_o = (starfile[:,1]==4)&(starfile[:,8]==1)&(starfile[:,9]<3.09)
		exsitu = (starfile[:,1]==4)&(starfile[:,8]==0)
		exsitu_y = (starfile[:,1]==4)&(starfile[:,8]==0)&(starfile[:,9]>3.09)
		exsitu_o = (starfile[:,1]==4)&(starfile[:,8]==0)&(starfile[:,9]<3.09)
		sp_in = starfile[insitu,2:5]
		sp_ex = starfile[exsitu,2:5]

		sp_in_y = starfile[insitu_y,2:5]
		sp_ex_y = starfile[exsitu_y,2:5]
		sp_in_o = starfile[insitu_o,2:5]
		sp_ex_o = starfile[exsitu_o,2:5]
		spos_in = np.sqrt((sp_in[:,0]-dsc[0])**2+(sp_in[:,1]-dsc[1])**2+(sp_in[:,2]-dsc[2])**2)/.71/(red+1)
		spos_in_y = np.sqrt((sp_in_y[:,0]-dsc[0])**2+(sp_in_y[:,1]-dsc[1])**2+(sp_in_y[:,2]-dsc[2])**2)/.71/(red+1)
		spos_in_o = np.sqrt((sp_in_o[:,0]-dsc[0])**2+(sp_in_o[:,1]-dsc[1])**2+(sp_in_o[:,2]-dsc[2])**2)/.71/(red+1)
		spos_ex = np.sqrt((sp_ex[:,0]-dsc[0])**2+(sp_ex[:,1]-dsc[1])**2+(sp_ex[:,2]-dsc[2])**2)/.71/(red+1)
		spos_ex_y = np.sqrt((sp_ex_y[:,0]-dsc[0])**2+(sp_ex_y[:,1]-dsc[1])**2+(sp_ex_y[:,2]-dsc[2])**2)/.71/(red+1)
		spos_ex_o = np.sqrt((sp_ex_o[:,0]-dsc[0])**2+(sp_ex_o[:,1]-dsc[1])**2+(sp_ex_o[:,2]-dsc[2])**2)/.71/(red+1)
		sr_in_bin, bin_edge = np.histogram(spos_in,bins=binz, weights =13.73-starfile[insitu,11])
		sr_in_bin_o, bin_edge = np.histogram(spos_in_o,bins=binz, weights =starfile[insitu_o,12])
		sr_in_bin_num, bin_edge = np.histogram(spos_in,bins=binz)
		sr_ex_bin, bin_edge = np.histogram(spos_ex,bins=binz, weights =13.73-starfile[exsitu,11])
		sr_ex_bin_o, bin_edge = np.histogram(spos_ex_o,bins=binz, weights =starfile[exsitu_o,12])
		sr_ex_bin_num, bin_edge = np.histogram(spos_ex,bins=binz)
		sr_in_bin = sr_in_bin/sr_in_bin_num
		sr_ex_bin = sr_ex_bin/sr_ex_bin_num
		sr_in_bin_o = np.cumsum(sr_in_bin_o)/(np.cumsum(sr_in_bin_o)+np.cumsum(sr_ex_bin_o))
		#sr_in_bin_y = 1-np.cumsum(sr_in_bin_y)/(np.cumsum(sr_in_bin_y)+np.cumsum(sr_ex_bin_y))
		#print sr_in_bin_y
		sr_ex_bin_o = 1-sr_in_bin_o#np.cumsum(sr_ex_bin)/(sum(sr_in_bin)+sum(sr_ex_bin))
		#sr_ex_bin_y = 1-sr_in_bin_y
		self.sub.plot(x,sr_in_bin,'k',label='insitu')
		self.sub.plot(x,sr_ex_bin,'r',label='merger')
		#self.sub.plot(x,sr_in_bin_o,'k',label='insitu_o')
		#self.sub.plot(x,sr_ex_bin_o,'r',label='merger_o')
		hdf.close()
	
	def save(self,hnum,date):
		#self.sub.legend(loc=1)#,prop={'size':10})
		self.fig.savefig('mean_age_vs_r_halo%s_%s.pdf'%(hnum,date),transparent=True)#stellar_mass_frac_old_vs_r_all_%s.pdf'%date)
		self.fig.show()

if __name__ == "__main__":
	hnum =['11707','32503','12596','007','848','796','20192','20910','897','948','1016','32257','10q','10v']
	res = ['_13','_13','_13','_11','_13','_13','_13','_13','_13','_13','_13','_13','_11','_11','_13']
	ver = ['5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16']
	snap = [184,184,184,184,184,184,184,184,184,184,184,184,600,600]
	date = time.strftime("%m_%d_%Y")
	plot = mean_age_vs_r()
	for i in np.arange(len(hnum)):
		#plot = mean_age_vs_r()
		pathname = '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz5_12_16_raw_output/analysis/'%(hnum[i],res[i])
		plot.compile_and_plot(pathname,hnum[i],res[i],ver[i],snap[i])
	plot.save('all',date)
		


