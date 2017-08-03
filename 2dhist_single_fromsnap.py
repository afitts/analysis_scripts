import numpy as np
import sys
import os
import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as co
from matplotlib.colors import LogNorm
import pylab
import pygadgetreader as pg
import scipy.integrate
import time
import pandas as pd
import scipy.stats as stats
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

class hist2d_plot(object):
	def __init__(self):
		self.fig = plt.figure()
		self.sub = self.fig.add_subplot(111)
		self.sub.set_xlabel('x (kpc)',fontsize=12)
		self.sub.set_ylabel('y (kpc)',fontsize=12)
	def add_plot(self,dx,dy,sx,sy,red):#,rho,red):
		#counts, xedges, yedges = np.histogram2d(dx, dy, bins=5e2, weights = np.log10(rho),normed=True)
                counts, xedges, yedges, Image = self.sub.hist2d(dx, dy, bins=1e3,norm=LogNorm())
		#counts, xedges, yedges, Image = self.sub.hist2d(dx, dy, bins=5e3, weights = rho,norm=LogNorm(),cmap=plt.cm.Greys)
		#im = mpl.image.NonUniformImage(self.sub, interpolation='bilinear',cmap=plt.cm.Greys)
		#xcenters = xedges[:-1] + 0.5 * (xedges[1:] - xedges[:-1])
		#ycenters = yedges[:-1] + 0.5 * (yedges[1:] - yedges[:-1])
		#im.set_data(xcenters, ycenters, counts)
		#self.sub.images.append(im)
		#ax.set_xlim(xedges[0], xedges[-1])
		#ax.set_ylim(yedges[0], yedges[-1])
		#self.sub.imshow(counts, origin = "center", interpolation = "none",cmap=plt.cm.Greys, extent=[-505,505,-505,505])
                #cbar = self.fig.colorbar(im, ax = self.sub)
                #cbar.set_label(r'$\mathrm{Density}$ $\left(\frac{\mathrm{g}}{\mathrm{cm}^3}\right)$')
		self.sub.scatter(sx,sy,s=10,lw=0.0)
		#self.sub.text(0.15, 0.15,'$\mathrm{z}$ = %.3f'%red, ha='center', va='center', transform=self.sub.transAxes, fontsize = 15,color='k')
        def save(self,hnum,snum,red):
		self.sub.set_xlim(-30,30)
		self.sub.set_ylim(-30,30)
                newpath = '2dhist_%s_plots/2dhist_halo%s'%('Z14',hnum)
                if not os.path.exists(newpath):
                        os.makedirs(newpath)
                #plt.savefig('%s/2dhist_halo%s_%03d.png'%(newpath,'both796',snum))
                plt.savefig('%s/2dhist_halo%s_z%f.pdf'%(newpath,hnum,red),transparent=True)
                plt.show()



if __name__ == "__main__":
	i=rank
	hnum = '848'
	res = '_14'
	ver = '5_12_16'
	pathname = '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output'%(hnum,res,ver)
	#hist = np.genfromtxt('%s/analysis/halo%smerger_hist.txt'%(pathname,hnum))
	snum = 172#hist[i,0]
	hdf = pd.HDFStore('%s/analysis/dataframes/halo%s%s_giz%s_snap%03d.h5'%(pathname,hnum,res,ver,snum))
	xcen = np.float(hdf['props']['halox'])*1e3
	ycen = np.float(hdf['props']['haloy'])*1e3
	zcen = np.float(hdf['props']['haloz'])*1e3
	red = np.float(hdf['props']['redshift'])
        rvir = np.float(hdf['props']['rvir'])*1e3
	dx = hdf['particles/dm']['x'].as_matrix()*1e3-xcen
	dy = hdf['particles/dm']['y'].as_matrix()*1e3-ycen
	dz = hdf['particles/dm']['z'].as_matrix()*1e3-zcen
	#dr = hdf['particles/dm']['r'].as_matrix()
	#dx = dx[dr<10*rvir]
        #dy = dy[dr<10*rvir]
	#gx = hdf['particles/gas']['x'].as_matrix()*1e3-xcen
	#gy = hdf['particles/gas']['y'].as_matrix()*1e3-ycen
	#gz = hdf['particles/gas']['z'].as_matrix()*1e3-zcen
	sx = hdf['particles/star']['x'].as_matrix()*1e3-xcen
	sy = hdf['particles/star']['y'].as_matrix()*1e3-ycen
	sz = hdf['particles/star']['z'].as_matrix()*1e3-zcen
	hinv = 1/.71
	ascale = 1/(red+1)
	#rho = hdf['particles/gas']['rho'].as_matrix()*1e10*1.99e33/(3.806e21)**3*(hinv/((ascale*hinv)**3))
	test_plot = hist2d_plot()
	test_plot.add_plot(dx,dy,sx,sy,red)
	test_plot.save(hnum,snum,red)
	#switch =0
	#numfiles = 16
	#for j in np.arange(numfiles): #Bring together all the AHF files for the snapshot
	#	temph = glob.glob(pathname+'analysis/ahf_snap%03d/ahf.snap_%03d.%04d.z*.*.AHF_halos'%(snum,snum,j))
	#	temph = str(temph).strip('[]').replace("'","")
	#	h = np.genfromtxt(temph)
	#	if switch == 0 and len(h) >0:
	#		halo = h
	#		switch = 1
	#	if switch == 1:
	#		try:
	#			halo = np.vstack((halo,h))
	#		except:
	#			print "nothing there"
	#for j in np.arange(len(halo)):
	#	if halo[j,0] == 
