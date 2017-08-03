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
import scipy.stats as stats
import pandas as pd
import scipy.optimize as opt
from matplotlib import rcParams
import matplotlib.animation as animation
from matplotlib.ticker import FixedLocator
from matplotlib.ticker import FixedFormatter

rcParams['lines.linewidth'] = 2
rcParams['axes.linewidth'] = 2
rcParams['xtick.major.size'] = 10
rcParams['xtick.minor.size'] = 5
rcParams['xtick.major.width'] = 2
rcParams['xtick.minor.width'] = 2
rcParams['ytick.major.size'] = 10
rcParams['ytick.minor.size'] = 5
rcParams['ytick.major.width'] = 2
rcParams['ytick.minor.width'] = 2
rcParams['font.size'] = 10
rcParams['xtick.labelsize']= '14'
rcParams['ytick.labelsize']= '14'
rcParams['savefig.bbox'] = 'tight' 

class radpro(object):

	def __init__(self, pathname, hnum, res, ver, dmo, i, grain):
		""" Loads a H5 pandas object and all relevant quantities""" 
		self.df = pd.HDFStore('%sdataframes/halo%s%s_giz%s_snap%03d.h5'%(pathname,hnum,res,ver,i))
		self.red = np.float(self.df['props']['redshift'])
		self.rvir = np.float(self.df['props']['rvir'])*1000/(self.red+1)
		self.rhalf = np.float(self.df['props']['rhalf'])*1000/(self.red+1)
		self.vmax = np.float(self.df['props']['vmax'])
		self.dmass = self.df['particles/dm']['mass'].as_matrix()
		self.dpos =  self.df['particles/dm']['r'].as_matrix()*1000/(self.red+1)
		binz = np.logspace(np.log10(.06),np.log10(self.rvir),grain)
		self.massall, bin_edge = np.histogram(self.dpos,bins=binz, weights = self.dmass) 
		self.x = 10**(np.log10(binz[1:])-np.log10(binz[1]/binz[0])/2)
		self.dlogr = np.log(binz[1]/binz[0])	
		if dmo == 0:
		  self.gmass = self.df['particles/gas']['mass'].as_matrix()
		  self.gpos = self.df['particles/gas']['r'].as_matrix()*1000/(self.red+1)
		  self.smass = self.df['particles/star']['mass'].as_matrix()
		  self.spos = self.df['particles/star']['r'].as_matrix()*1000/(self.red+1)
		  self.temp = self.df['particles/gas']['temp'].as_matrix()
		  self.den = self.df['particles/gas']['rho'].as_matrix()*1e10*1.99e33/1.67e-24/(3.806e24)**3
		  self.sfr = self.df['particles/gas']['sfr'].as_matrix()
		  he = self.df['particles/star']['metal_He'].as_matrix()
		  fe = self.df['particles/star']['metal_Fe'].as_matrix()
		  metal = self.df['particles/star']['metal_tot'].as_matrix()
		  numH = self.smass*(1-(metal+he))/(1.6733e-24/1.99e33)
		  numFe = self.smass*fe/(9.27e-23/1.99e33)
		  avgnum = np.mean(numFe/numH)
		  self.fullmetal = np.log10(numFe/numH)+4.5
		  self.Z = np.log10(avgnum)+4.5
		  if hnum == '007' or res == '5_12_16':
		    den= den*(1e3)**3
		  self.gmassall, bin_edge = np.histogram( gpos,bins=binz, weights =gmass) 
		  self.smassall, bin_edge = np.histogram( spos,bins=binz, weights =smass)
		  self.tempall, bin_edge, binnum= stats.binned_statistic(gpos, np.log10(temp), statistic='mean', bins=binz)
		  self.tbin = np.logspace(np.log10(min(temp)),np.log10(max(temp)),len(temp+1))
		  cumutempall, bin_edge = np.histogram(temp,bins=tbin, weights =np.ones(len(temp)))   
		  self.cumutempall = np.cumsum(cumutempall)
		  self.df.close()
    return massall, gmassall, smassall, x, rvir,rhalf,rmax, temp, den, sfr,red,dlogr, tempall,tbin[1:],cumutempall,cNFW,vmax,fullmetal,metal
		else:
		  count += 1
		self.df.close()
    return massall, x,cNFW,vmax	
	def plot(self)




def radpro_df(pathname,sname,hist,dmo,rhalf,rmax,i):
  global grain, clr, count
  massp = 1.67262178e-24 #in grams
  gamma = 5./3.
  kb = 1.3806e-26 #in km^2*g/s^2 
  G = 4.3e-6 #in kpc/M_sun (km/s)^2   
  switch = 0
  hdf = pd.HDFStore('%sdataframes/halo%s%s_giz%s_snap%03d.h5'%(pathname,hnum[w],res[w],ver[w],i))
  red = np.float(self.df['props']['redshift'])
  rvir = np.float(self.df['props']['rvir'])*1000/(red+1)
  rhalf = np.float(self.df['props']['rhalf'])*1000/(red+1)
  vmax = np.float(self.df['props']['vmax'])
  dmass = self.df['particles/dm']['mass'].as_matrix()
  dpos =  self.df['particles/dm']['r'].as_matrix()*1000/(red+1)
  binz = np.logspace(np.log10(.06),np.log10(rvir),grain)
  massall, bin_edge = np.histogram(dpos,bins=binz, weights =dmass) 
  x = 10**(np.log10(binz[1:])-np.log10(binz[1]/binz[0])/2)
  dlogr = np.log(binz[1]/binz[0])
  if dmo == 0:
    gmass = self.df['particles/gas']['mass'].as_matrix()
    gpos = self.df['particles/gas']['r'].as_matrix()*1000/(red+1)
    smass = self.df['particles/star']['mass'].as_matrix()
    spos = self.df['particles/star']['r'].as_matrix()*1000/(red+1)
    temp = self.df['particles/gas']['temp'].as_matrix()
    den = self.df['particles/gas']['rho'].as_matrix()*1e10*1.99e33/1.67e-24/(3.806e24)**3
    sfr = self.df['particles/gas']['sfr'].as_matrix()
    he = self.df['particles/star']['metal_He'].as_matrix()
    fe = self.df['particles/star']['metal_Fe'].as_matrix()
    metal = self.df['particles/star']['metal_tot'].as_matrix()
    numH = smass*(1-(metal+he))/(1.6733e-24/1.99e33)
    numFe = smass*fe/(9.27e-23/1.99e33)
    avgnum = np.mean(numFe/numH)
    fullmetal = np.log10(numFe/numH)+4.5
    metal = np.log10(avgnum)+4.5
    if hnum[w] == '007' or res[w] == '5_12_16':
      den= den*(1e3)**3
    gmassall, bin_edge = np.histogram( gpos,bins=binz, weights =gmass) 
    smassall, bin_edge = np.histogram( spos,bins=binz, weights =smass)
    tempall, bin_edge, binnum= stats.binned_statistic(gpos, np.log10(temp), statistic='mean', bins=binz)
    tbin = np.logspace(np.log10(min(temp)),np.log10(max(temp)),len(temp+1))
    cumutempall, bin_edge = np.histogram(temp,bins=tbin, weights =np.ones(len(temp)))   
    cumutempall = np.cumsum(cumutempall)
    hdf.close()
    return massall, gmassall, smassall, x, rvir,rhalf,rmax, temp, den, sfr,red,dlogr, tempall,tbin[1:],cumutempall,cNFW,vmax,fullmetal,metal
  else:
    count += 1
    hdf.close()
    return massall, x,cNFW,vmax
