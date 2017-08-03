import numpy as np
#import matplotlib 
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.cm as cm
import matplotlib.colors as co
from matplotlib import rcParams
from matplotlib import rc
import pandas as pd
import scipy.stats as st
import os
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
rcParams['savefig.bbox'] = 'tight' 

class hist2d_plot(object):
	def __init__(self):
		self.fig, (self.trip1, self.trip2) = plt.subplots(2,2,sharex=True,sharey=True,figsize=plt.figaspect(.9))
		#self.fig = plt.figure()
		#self.trip2 = self.fig.add_subplot(111)
		self.trip2[0].set_xlabel('x (kpc)',fontsize=12)
		self.trip2[1].set_xlabel('x (kpc)',fontsize=12)
		self.trip1[0].set_ylabel('y (kpc)',fontsize=12)
		self.trip2[0].set_ylabel('y (kpc)',fontsize=12)
		
	def add_plots(self,hnum,res,snum,z):
		hdf = pd.HDFStore('/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz5_12_16_raw_output/analysis/dataframes/halo%s%s_giz5_12_16_snap%03d.h5'%(hnum,res,hnum,res,snum))
		try:
			xcen = hdf['props']['halox'].values[0]
			ycen = hdf['props']['haloy'].values[0]
			x = hdf['particles/dm']['x'].values-xcen
			y = hdf['particles/dm']['y'].values-ycen
			print x*1e3,y*1e3
			counts, xedges, yedges, Image = self.trip1[0].hist2d(x*1e3, y*1e3, bins=1e3, norm=LogNorm())
		except Exception,e:
			print 'No DM'
		try:
			x = hdf['particles/star']['x'].values-xcen
			y = hdf['particles/star']['y'].values-ycen
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
			sft = hdf['particles/star']['sft'].as_matrix()
			#self.trip1[1].hist2d(x*1e3, y*1e3, bins=[xedges,yedges], norm=LogNorm())
			self.trip1[1].scatter(x*1e3,y*1e3,s=40,color='b',edgecolor='k')
			#self.trip2.scatter(x[fullmetal<-3.83]*1e3,y[fullmetal<-3.83]*1e3,s=40,color='r')
			#self.trip2.scatter(x[(fullmetal<-3.83)&(sft>0.689176245211)]*1e3,y[(fullmetal<-3.83)&(sft>0.689176245211)]*1e3,s=40,color='g')
			self.trip1[1].set_xlim(min(xedges),max(xedges))
			self.trip1[1].set_ylim(min(yedges),max(yedges))#-3,3)
		except Exception,e:
			print 'No Stars'
		hdf = pd.HDFStore('/nobackup/afitts/Gadget-2.0.7/production/mfmw%s%s_giz5_12_16_raw_output/analysis/dataframes/halow%s%s_giz5_12_16_snap%03d.h5'%(hnum,res,hnum,res,snum))
		try:
			xcen = hdf['props']['halox'].values[0]
			ycen = hdf['props']['haloy'].values[0]
			x = hdf['particles/dm']['x'].values-xcen
			y = hdf['particles/dm']['y'].values-ycen
			counts, xedges, yedges, Image = self.trip2[0].hist2d(x, y, bins=1e3, norm=LogNorm())
		except Exception,e:
			print 'No DM',e
		try:
			x = hdf['particles/star']['x'].values-xcen
			y = hdf['particles/star']['y'].values-ycen
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
			sft = hdf['particles/star']['sft'].as_matrix()
			#self.trip2[1].hist2d(x, y, bins=[xedges,yedges], norm=LogNorm())
			self.trip2[1].scatter(x,y,s=40,color='b',edgecolor='k')
			#self.trip2.scatter(x[fullmetal<-3.83]*1e3,y[fullmetal<-3.83]*1e3,s=40,color='r')
			#self.trip2.scatter(x[(fullmetal<-3.83)&(sft>0.689176245211)]*1e3,y[(fullmetal<-3.83)&(sft>0.689176245211)]*1e3,s=40,color='g')
			self.trip2[1].set_xlim(min(xedges),max(xedges))
			self.trip2[1].set_ylim(min(yedges),max(yedges))#-3,3)
		except Exception,e:
			print 'No Stars'
		#self.trip2.set_title('z=%.3f'%z)
		try:
			bob = bob1
			x = hdf['particles/gas']['x'].values-xcen
			y = hdf['particles/gas']['y'].values-ycen
			self.trip3.hist2d(x*1e3, y*1e3, bins=[xedges,yedges], norm=LogNorm())
		
		except Exception,e:
			print 'No Gas'
		#plt.setp(self.trip1, aspect=0.75, adjustable='box-forced')
		#plt.setp(self.trip2, aspect=0.75, adjustable='box-forced')
		self.fig.subplots_adjust(wspace=0.00)
		self.fig.subplots_adjust(hspace=0.00)
		self.trip1[0].set_xlim(-10,10)
		self.trip1[0].set_ylim(-10,10)
		self.trip2[1].set_xticklabels(['','-5','0','5','10'])
		#plt.setp(self.trip3, aspect=0.75, adjustable='box-forced')
		hdf.close()
	def save(self,hnum,snum):
		#plt.colorbar()
		newpath = '2dhist_%s_plots/2dhist_halo%s'%('wdm',hnum)
		if not os.path.exists(newpath):
			os.makedirs(newpath)
		#plt.savefig('%s/2dhist_halo%s_%03d.png'%(newpath,'both796',snum))
		plt.savefig('%s/2dhist_halo%s_%03d.pdf'%(newpath,'both796',snum),transparent=True)
		plt.show()
		plt.close()

class josetemp_plot(object):
	def __init__(self):
		self.fig= plt.figure()
		self.sub = self.fig.add_subplot(111)
		self.sub.set_xlabel('x (kpc)',fontsize=12)
		self.sub.set_ylabel('y (kpc)',fontsize=12)
		#my_cmap=plt.get_cmap('jet')
		#self.sm = cm.ScalarMappable(cmap=my_cmap,norm=co.Normalize(vmin=4, vmax=7))#vmin=5.874853, vmax=7.11806)) 11_13 limits
		#self.sm._A = []
		#cb = self.fig.colorbar(self.sm)
		#cb.set_label(r'$\mathrm{log}$ T (K)')
	def add_plots(self,hnum,res,snum,z):
		hdf = pd.HDFStore('/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz5_12_16_raw_output/analysis/dataframes/halo%s%s_giz5_12_16_snap%03d.h5'%(hnum,res,hnum,res,snum))
		self.sub.set_title('z=%.3f'%z)
		xcen = hdf['props']['halox'].values[0]/.71
		ycen = hdf['props']['haloy'].values[0]/.71
		zcen = hdf['props']['haloz'].values[0]/.71
		x = hdf['particles/gas']['x'].values/.71-xcen
		y = hdf['particles/gas']['y'].values/.71-ycen
		z = np.absolute(hdf['particles/gas']['z'].values/.71-zcen)
		temp = hdf['particles/gas']['temp'].values
		#xx, yy = np.mgrid[-40:40:100j, -40:40:100j]
		#positions = np.vstack([xx.ravel(), yy.ravel()])
		#values = np.vstack([x*1e3, y*1e3])
		#kernel = st.gaussian_kde(values)
		#f = np.reshape(kernel(positions).T, xx.shape)

		#cax = self.trip2.imshow(np.rot90(f), cmap='jet',vmin = 4,vmax=7,
           #extent=[-40, 40, -40, 40])
		#counts, xedges, yedges, Image =self.trip1.hist2d(x[r<0.002]*1e3, y[r<0.002]*1e3, bins=1e3,weights = temp[r<0.002])
		#self.trip2.hist2d(x*1e3, y*1e3, bins=[xedges,yedges], weights=temp)
		#self.trip2.pcolor(x, y,temp[r<0.002],cmap = 'plasma')	
		n = np.histogram2d(x[z<0.001]*1e3, y[z<0.001]*1e3,bins = 1e2)[0]
		temp,xedges,yedges = np.histogram2d(x[z<0.001]*1e3, y[z<0.001]*1e3,bins = 1e2,weights=temp[z<0.001])
		cax = self.sub.imshow(np.log10(temp/n), cmap='jet',vmin = 4,vmax=6,
           extent=[-40, 40, -40, 40],
           interpolation='gaussian', origin='lower')
		cb = plt.colorbar(cax)
		cb.set_label(r'$\mathrm{log}$ T (K)')
		plt.show()
		hdf.close()
	def save(self,hnum,snum):
		plt.savefig('2dhist_josetemp_halo%s_2kpc.pdf'%(hnum))
		#plt.show()

def main(hnum,res,snum):
	if hnum == '10q' or hnum == '10v':
		a = np.genfromtxt('snapshot_scale-factors.txt')
	else:
		a = np.genfromtxt('output_times.txt')
	z = (1/a)-1
	#josetemp = josetemp_plot()
	#josetemp.add_plots(hnum,res,184,z[184])
	#josetemp.save(hnum,184)
	snum = [184]
	for j in np.arange(len(snum)):
		h2dplot = hist2d_plot()
		print hnum,snum[j],'%.3f'%z[snum[j]]
		h2dplot.add_plots(hnum,res,snum[j],z[int(snum[j])])
		h2dplot.save(hnum,snum[j])

hnum =['848']#,'w20192','w32503','w796','w848','w897','w948','w1016']#11707','32503','12596','007','848','796','20192','20910','897','948','1016','32257']#,'10q','10v','1084']#
res = ['_13','_13','_13','_13','_13','_13','_13','_13','','','','','','','']#['_11','_11']#,'_13','_11','_13','_13','_13','_13','_13','_13','_13','_13','_11','_11','_13']
snum = np.linspace(160,184,25)
if hnum[rank] == '10q' or hnum[rank] == '10v':
	snum = np.linspace(0,600,601)
#for i in np.arange(len(hnum)):
	#print hnum[i]
main(hnum[rank],res[rank],snum)
