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

def load_vel(pathname,hnum,i):

	massp = 1.67262178e-24 #in grams
	gamma = 5./3.
	kb = 1.3806e-26 #in km^2*g/s^2
	cutoff = 20 #kpc
	solartog = 1.989e33

	sname = "%ssnapdir_%03d/snapshot_%03d"%(pathname,i,i)
	gid = pg.readsnap(sname, 'pid', 'gas')
	pos = pg.readsnap(sname, 'pos', 'gas') #in kpc
	gv = pg.readsnap(sname, 'vel', 'gas') #in km/s
	gm = pg.readsnap(sname, 'mass', 'gas')*1.e10/.71
	gm = gm.astype('float64')
	zarray = pg.readsnap(sname, 'zarray', 'gas')
	he_mf = zarray[:,1] #Helium mass fraction
	y_he = he_mf/(4*(1-he_mf))
	ne = pg.readsnap(sname, 'ne', 'gas') #Electron Abundance
	mu = (1+4*y_he)/(1+y_he+ne)
	mmw = mu*massp #mean molecular weight
	u = pg.readsnap(sname,'u','gas') #specific internal energy in km^2/s^2
	temp = mmw * (gamma-1.)*u/kb #temperature of gas
	
	hdf = pd.HDFStore('/nobackup/afitts/Gadget-2.0.7/production/mfm%s_13_giz5_12_16_raw_output/analysis/dataframes/halo%s_13_giz5_12_16_snap%03d.h5'%(hnum,hnum,i))
	vx = np.float(hdf['props']['halo_vx'])*1000
	vy = np.float(hdf['props']['halo_vy'])*1000
	vz = np.float(hdf['props']['halo_vz'])*1000
	hx = np.float(hdf['props']['halox'])*1000
	hy = np.float(hdf['props']['haloy'])*1000
	hz = np.float(hdf['props']['haloz'])*1000

	gr = np.sqrt((pos[:,0]-hx)**2+(pos[:,1]-hy)**2+(pos[:,2]-hz)**2)/.71
	gv = gv[gr<cutoff] #Only focus on particles within cutoff kpc of halo center
	gm = gm[gr<cutoff]*solartog
	temp = temp[gr<cutoff]
	mmw = mmw[gr<20]

	gv[:,0] = gv[:,0]-vx
	gv[:,1] = gv[:,1]-vy
	gv[:,2] = gv[:,2]-vz
	gv = gv*1e5 #convert to cm/s
	ke = sum(.5*gm*(gv[:,0]**2+gv[:,1]**2+gv[:,2]**2))

	N = gm/mmw
	te = sum(3/2.*N*(kb*1e5**2)*temp)
	hdf.close()
	return gid,gv,ke,te

def load_vel_df(pathname,hnum,i):

	massp = 1.67262178e-24 #in grams
	gamma = 5./3.
	kb = 1.3806e-26 #in km^2*g/s^2
	cutoff = 20 #kpc
	solartog = 1.989e33

	#sname = "%ssnapdir_%03d/snapshot_%03d"%(pathname,i,i)
	#gid = pg.readsnap(sname, 'pid', 'gas')
	#pos = pg.readsnap(sname, 'pos', 'gas') #in kpc
	#gv = pg.readsnap(sname, 'vel', 'gas') #in km/s
	#gm = pg.readsnap(sname, 'mass', 'gas')*1.e10/.71
	#gm = gm.astype('float64')
	#zarray = pg.readsnap(sname, 'zarray', 'gas')
	#he_mf = zarray[:,1] #Helium mass fraction
	#y_he = he_mf/(4*(1-he_mf))
	#ne = pg.readsnap(sname, 'ne', 'gas') #Electron Abundance
	#mu = (1+4*y_he)/(1+y_he+ne)
	#mmw = mu*massp #mean molecular weight
	#u = pg.readsnap(sname,'u','gas') #specific internal energy in km^2/s^2
	#temp = mmw * (gamma-1.)*u/kb #temperature of gas
	
	hdf = pd.HDFStore('%sanalysis/dataframes/halo%s_13_giz5_12_16_snap%03d.h5'%(pathname,hnum,i))
	gvx = hdf['particles/gas']['vx'].as_matrix()#*1000
	gvy = hdf['particles/gas']['vx'].as_matrix()#*1000
	gvz = hdf['particles/gas']['vx'].as_matrix()#*1000
	gv = np.column_stack((gvx,gvy,gvz))
	gm = hdf['particles/gas']['mass'].as_matrix()
	gr = hdf['particles/gas']['r'].as_matrix()*1000
	temp = hdf['particles/gas']['temp'].as_matrix()

	gv = gv[gr<cutoff] #Only focus on particles within cutoff kpc of halo center
	gm = gm[gr<cutoff]*solartog
	temp = temp[gr<cutoff]
	#mmw = mmw[gr<20]

	gv = gv*1e5 #convert to cm/s
	ke = sum(.5*gm*(gv[:,0]**2+gv[:,1]**2+gv[:,2]**2))

	#N = gm/mmw
	#te = sum(3/2.*N*(kb*1e5**2)*temp)
	hdf.close()
	return gv,ke

class metal_dist(object):
	""" Metallicity Distributions (Cumulative Metalicity Distribution under construction)"""
	def __init__(self):
		self.fig = plt.figure()
		self.sub = self.fig.add_subplot(111)
		self.sub.set_xlabel('$\mathrm{[Fe/H]}$')
		self.sub.set_ylabel('# $\mathrm{of}$ $\mathrm{Stars}$')
		self.sub.set_xlim(-4,0)
		self.sub.set_ylim(0,350)
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
new_gid,new_gv,new_ke,new_te = load_vel(pathname,hnum,snum)
old_gid,old_gv,old_ke,old_te = load_vel(pathname,hnum,snum-1)
print '%d KE %.2e, TE, %.2e'%(snum-1,old_ke,old_te)
print '%d KE %.2e, TE, %.2e'%(snum,new_ke,new_te)
#mhist = metal_dist()
#mhist.add_dist(fullmetal,grain)
#mhist.save(hnum,grain,date)
#new_gv,new_ke = load_vel_df(pathname,hnum,snum)
#old_gv,old_ke = load_vel_df(pathname,hnum,snum-1)
#print '%d KE %.2e, %d KE, %.2e'%(snum-1,old_ke,snum,new_ke)
