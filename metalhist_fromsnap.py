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

def load_metal(pathname,i):
	sname = "%ssnapdir_%03d/snapshot_%03d"%(pathname,i,i)
	smass = pg.readsnap(sname, 'mass', 'star')*1.e10/.71
	smass = smass.astype('float64')
	print len(smass)
	sz = pg.readsnap(sname, 'zarray','star')
	sz = sz.astype('float64')
	he = sz[:,1]
	fe = sz[:,10]
	metal = sz[:,0]
	numH = smass*(1-(metal+he))/(1.6733e-24/1.99e33)
	numFe = smass*fe/(9.27e-23/1.99e33)
	meta = numFe/numH
	nofloor = meta#[meta>4.54877795e-09]
	avgnum = np.mean(nofloor)
	fullmetal = np.log10(nofloor)+4.5
	return fullmetal

def load_metal_df(pathname,hnum,i):
	hdf = pd.HDFStore('%sanalysis/dataframes/halo%s_13_giz5_12_16_snap%03d.h5'%(pathname,hnum,i))
	smass = hdf['particles/star']['mass'].as_matrix()
	he = hdf['particles/star']['metal_He'].as_matrix()
	fe = hdf['particles/star']['metal_Fe'].as_matrix()
	metal = hdf['particles/star']['metal_tot'].as_matrix()
	numH = smass*(1-(metal+he))/(1.6733e-24/1.99e33)
	numFe = smass*fe/(9.27e-23/1.99e33)
	meta = numFe/numH
	nofloor = meta#[meta>4.54877795e-09]
	avgnum = np.mean(nofloor)
	fullmetal = np.log10(nofloor)+4.5
	metal = np.log10(avgnum)+4.5
	return fullmetal

class metal_dist(object):
	""" Metallicity Distributions (Cumulative Metalicity Distribution under construction)"""
	def __init__(self):
		self.fig = plt.figure()
		self.sub = self.fig.add_subplot(111)
		self.sub.set_xlabel('$\mathrm{[Fe/H]}$')
		self.sub.set_ylabel('# $\mathrm{of}$ $\mathrm{Stars}$')
		self.sub.set_xlim(-4,0)
		self.sub.set_ylim(0,200)
		#figm1 = plt.figure()
		#cumumetal = figm1.add_subplot(111)

  #metalall, bin_edge = np.histogram(fullmetal,bins=binz)
  #cumumetal.plot((binz[1:]-(binz[1]-binz[0])/2),np.cumsum(metalall)/float(max(np.cumsum(metalall))),linewidth=4)
  #cumumetal.set_xlabel('$\mathrm{[Fe/H]}$')
  #cumumetal.set_ylabel('$\mathrm{Cumulative}$ $\mathrm{Distribution}$')
  #cumumetal.set_ylim(-.05,1.05)
  #figm1.savefig('cumu_metal_halo%s_%dbin_%s_z0.pdf'%(hnum,grain,date),Transparent=True)

	def add_dist(self,fullmetal,grain):
		binz = np.linspace(min(fullmetal),max(fullmetal),grain)
		self.sub.hist(fullmetal,binz)
	def save(self,hnum,grain,date):
		self.fig.savefig('metalhist_halo%s_%dbin_%s_z0.pdf'%(hnum,grain,date),Transparent=True)
		self.fig.show()

hnum = sys.argv[1]
snum = np.float(sys.argv[2])
date = time.strftime("%m_%d_%Y")

pathname = '/nobackup/afitts/Gadget-2.0.7/production/mfm%s_13_giz5_12_16_raw_output/'%(hnum)
grain = 100

fullmetal = load_metal_df(pathname,hnum,snum)
mhist = metal_dist()
mhist.add_dist(fullmetal,grain)
mhist.save(hnum,grain,date)
