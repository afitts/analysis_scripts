import numpy as np
import sys 
import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as co
import pylab
import pygadgetreader as pg
import scipy.integrate
import time as dtime
import scipy.stats as stats
import pandas as pd
import scipy.optimize as opt
from profiles import *
from mikecm import mikecm
from matplotlib import rcParams
import matplotlib.animation as animation
from matplotlib.ticker import FixedLocator
from matplotlib.ticker import FixedFormatter
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

rcParams['lines.linewidth'] = 4
rcParams['axes.linewidth'] = 2
rcParams['xtick.major.size'] = 10
rcParams['xtick.minor.size'] = 5
rcParams['xtick.major.width'] = 2
rcParams['xtick.minor.width'] = 2
rcParams['ytick.major.size'] = 10
rcParams['ytick.minor.size'] = 5
rcParams['ytick.major.width'] = 2
rcParams['ytick.minor.width'] = 2
rcParams['font.size'] = 20
rcParams['xtick.labelsize']= '16'
rcParams['ytick.labelsize']= '16'
rcParams['savefig.bbox'] = 'tight' 

def fit_profile(hnum,res,ver,dmo,snum):
	grain = 60
	if dmo == 0: ###For hydro runs ###
		### Load in dataframe ###
		pathname = '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/'%(hnum,res,ver)
		hdf = pd.HDFStore('%sanalysis/dataframes/halo%s%s_giz%s_snap%03d.h5'%(pathname,hnum,res,ver,snum))
		red = np.float(hdf['props']['redshift'])
		time = np.float(hdf['props']['time'])
		rhalf = np.float(hdf['props']['rhalf'])*1000/(red+1)
		rvir = np.float(hdf['props']['rvir'])*1000/(red+1)
		dm = hdf['particles/dm']['mass'].as_matrix()
		gm = hdf['particles/gas']['mass'].as_matrix()
		sm = hdf['particles/star']['mass'].as_matrix()

		### Recentering procedure ###
		dx = hdf['particles/dm']['x'].as_matrix()
		dy = hdf['particles/dm']['y'].as_matrix()
		dz = hdf['particles/dm']['z'].as_matrix()
		gx = hdf['particles/gas']['x'].as_matrix()
		gy = hdf['particles/gas']['y'].as_matrix()
		gz = hdf['particles/gas']['z'].as_matrix()
		sx = hdf['particles/star']['x'].as_matrix()
		sy = hdf['particles/star']['y'].as_matrix()
		sz = hdf['particles/star']['z'].as_matrix()
		dp = np.column_stack((dx,dy,dz))
		sp = np.column_stack((sx,sy,sz))
		dsp = np.vstack((dp,sp))
		dsm = np.append(dm,sm)
		dsc = mikecm(fname = dsp,nofile=True,pmass=dsm)
		dpos = np.sqrt((dp[:,0]-dsc[0])**2+(dp[:,1]-dsc[1])**2+(dp[:,2]-dsc[2])**2)/.71*1000/(red+1)
		gpos = np.sqrt((gp[:,0]-dsc[0])**2+(gp[:,1]-dsc[1])**2+(gp[:,2]-dsc[2])**2)/.71*1000/(red+1)
		spos = np.sqrt((sp[:,0]-dsc[0])**2+(sp[:,1]-dsc[1])**2+(sp[:,2]-dsc[2])**2)/.71*1000/(red+1)

		### Binning procedure ###
		binz = np.logspace(np.log(.0014),np.log(rvir),grain,base=np.e)
		x = np.e**(np.log(binz[1:])-np.log(binz[1]/binz[0])/2)
		dlogr = np.log(binz[1]/binz[0])
		massall, bin_edge = np.histogram(dpos,bins=binz, weights =dm)
		gmassall, bin_edge = np.histogram( gpos,bins=binz, weights =gm) 
		smassall, bin_edge = np.histogram( spos,bins=binz, weights =sm)
		totden = (massall+gmassall+smassall)/(4*3.14159*x**3)/dlogr

		### Chi-squared minimization fitting routine ###
		chimin =10000
		sig = 0.1 * totden
		cenden = totden[np.sqrt((x-.3)**2).argmin()]
		for p in np.logspace(-1,np.log10(10),100):
			(fit, cmatrix)= opt.curve_fit(psuedoiso,x,totden,p0=(cenden,p),sigma=sig)
			chisq = sum((totden-psuedoiso(x,*fit))**2/(sig)**2)
			if chisq < chimin:
				chimin = chisq
				bestfit = fit
	else: ### For DMO runs ###
		### Load in dataframe ###
		pathname = '/nobackup/afitts/Gadget-2.0.7/production/gizdmsi%s%s_raw_output/'%(hnum,res)
		hdf = pd.HDFStore('%sanalysis/dataframes/halo%s%s_giz%s_snap%03d.h5'%(pathname,hnum,res,ver,snum))
		red = np.float(hdf['props']['redshift'])
		time = np.float(hdf['props']['time'])
		c = np.float(hdf['props']['cNFW'])
		rhalf = np.float(hdf['props']['rhalf'])*1000/(red+1)
		rvir = np.float(hdf['props']['rvir'])*1000/(red+1)
		dm = hdf['particles/dm']['mass'].as_matrix()
		mvir = sum(dm)

		### Recentering procedure ###
		dx = hdf['particles/dm']['x'].as_matrix()
		dy = hdf['particles/dm']['y'].as_matrix()
		dz = hdf['particles/dm']['z'].as_matrix()
		dp = np.column_stack((dx,dy,dz))
		dsc = mikecm(fname = dp,nofile=True,pmass=dm)
		dpos = np.sqrt((dp[:,0]-dsc[0])**2+(dp[:,1]-dsc[1])**2+(dp[:,2]-dsc[2])**2)/.71*1000/(red+1)

		### Binning procedure ###
		binz = np.logspace(np.log(.0014),np.log(rvir),grain,base=np.e)
		x = np.e**(np.log(binz[1:])-np.log(binz[1]/binz[0])/2)
		dlogr = np.log(binz[1]/binz[0])
		massall, bin_edge = np.histogram(dpos,bins=binz, weights =dm)
		totden = (massall)/(4*3.14159*x**3)/dlogr

		### Chi-squared minimization fitting routine ###
		chimin =10000
		sig = 0.1 * totden
		prad = .3 #power radius (approximation for now)
		cenden0 = totden[np.sqrt((x-prad)**2).argmin()]

		param_bounds=([cenden0*1e-1,.499],[cenden0*1e1,10])
		outxlim = 50
		fit_lims = (x<outxlim) & (x>prad)
		### Burkert fit ###
		for p in np.logspace(np.log10(.5),np.log10(2),20):
			for cenden in np.logspace(np.log10(cenden0*.5),np.log10(cenden0*1.5),20):
				(fit, cmatrix)= opt.curve_fit(burkert,x[fit_lims],totden[fit_lims],p0=(cenden,p),bounds=param_bounds)
				chisq = sum((totden[fit_lims]-burkert(x[fit_lims],*fit))**2/(sig[fit_lims])**2)
				if chisq < chimin:
					chimin = chisq
					bestfit = fit
		### Psuedo iso fit ###
		chimin =10000
		for p in np.logspace(np.log10(.5),np.log10(2),20):
			for cenden in np.logspace(np.log10(cenden0*.5),np.log10(cenden0*1.5),20):
				(fit, cmatrix)= opt.curve_fit(psuedoiso,x[fit_lims],totden[fit_lims],p0=(cenden,p),bounds=param_bounds)
				chisq = sum((totden[fit_lims]-psuedoiso(x[fit_lims],*fit))**2/(sig[fit_lims])**2)
				if chisq < chimin:
					chimin = chisq
					bestfit1 = fit
		### Cored-NFW fit ###
		chimin =10000
                def CNFW(mvir,c): return lambda r,n,rc : corenfw(r,mvir,c,n,rc)
                nfwcored = CNFW(mvir,c)
		param_bounds=([0.1,0.1],[1.1,5])
                for rc in np.logspace(np.log10(0.1),np.log10(5),20):
                        for n in np.logspace(np.log10(0.1),np.log10(1.1),20):
				(bestfit2, cmatrix)= opt.curve_fit(nfwcored,x[fit_lims],np.log(totden[fit_lims]),p0 = (n,rc),bounds=param_bounds)
                                chisq = sum((totden[fit_lims]-psuedoiso(x[fit_lims],*bestfit2))**2/(sig[fit_lims])**2)
                                if chisq < chimin:
                                        chimin = chisq
                                        bestfit2 = fit
		test_plot = radpro()
		test_plot.add_dmoline(x,totden,'SIDM')
		test_plot.add_fit(x,burkert(x,*bestfit),'bur')
		test_plot.add_fit(x,psuedoiso(x,*bestfit1),'piso')
		test_plot.add_fit(x,np.exp(corenfw(x,mvir,c,*bestfit2)),'cnfw')
		date = dtime.strftime("%m_%d_%Y")
		test_plot.save(date,hnum)
	return time, bestfit


hnum = '11707'
res = '_13'
ver = '5_12_16'
dmo = 1
time, coresize = fit_profile(hnum,res,ver,dmo,89)
time = comm.gather(time, root=0)
coresize = comm.gather(coresize[1], root=0)
if rank == 0:
	print 'HI',coresize
	a = plt.plot(time,coresize)
	plt.savefig('coresize_test.pdf',transparent=True)
