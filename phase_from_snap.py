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
import scipy.stats as stats
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
class plot_phase(object):
	def __init__(self):
		self.fig = plt.figure()
		self.sub = self.fig.add_subplot(111)
		self.sub.set_xlabel(r'$\mathrm{log(\rho)}$ $\mathrm{(n_H/cm^3)}$',fontsize=15, labelpad=-2)
		self.sub.set_ylabel(r'$\mathrm{log(Temperature)}$ $\mathrm{(K)}$',fontsize=15, labelpad=-2)
		self.sub.set_xlim(-6,3)
		self.sub.set_ylim(1,6)
		#plt.ylim(1e-1,2e0)
		plt.legend(loc='best',prop={'size':10})

	def add_plot(self,hymass,hygmass,hysmass,temp,den,red,grain):
		hytot = hymass+hygmass+hysmass
		nbins = 100
		#temp = np.nan_to_num(temp)
		H, denedges, tempedges = np.histogram2d(np.log10(den),np.log10(temp),bins=nbins)
		H = np.rot90(H)
		H = np.flipud(H)
 
		# Mask zeros
		Hmasked = np.ma.masked_where(H==0,H)
		phase = self.sub.pcolormesh(denedges,tempedges,np.log10(Hmasked))
		cbar = self.fig.colorbar(phase, ax = self.sub)
		cbar.set_label('$\mathrm{Log(Counts)}$')
		phase.set_clim(vmin=0,vmax=5)
                #self.sub.spines['bottom'].set_linewidth(4)
                #self.sub.spines['top'].set_linewidth(4)
                #self.sub.spines['right'].set_linewidth(4)
                #self.sub.spines['left'].set_linewidth(4)
		self.sub.text(0.15, 0.3,'$\mathrm{z}$ = %.3f'%red, ha='center', va='center', transform=self.sub.transAxes, fontsize = 15)
		self.sub.text(0.25, 0.2,'$\mathrm{M_{baryon}}$ = %.3e'%hygmass, ha='center', va='center', transform=self.sub.transAxes, fontsize = 15)
		#self.sub.text(0.225, 0.1,'$\mathrm{SFR}$ = %.3e'%sum(sfr), ha='center', va='center', transform=self.sub.transAxes, fontsize = 15)

	def save(self,hnum,snum,date,ex):
		self.fig.savefig('phase_halo%s_%dx_%03d.pdf'%(hnum,ex,snum),Transparent=True)
		self.fig.show()

class temp_vs_r(object):
	def __init__(self,sm):
		self.fig = plt.figure()
		self.sub = self.fig.add_subplot(111)
		self.sub1 = self.sub.twinx()
		self.sub.set_xlabel(r'$\mathrm{r/R_{vir}} (kpc)$',fontsize=15, labelpad=-2)
		self.sub.set_ylabel(r'$\mathrm{log(Temperature)}$ $\mathrm{(K)}$',fontsize=15, labelpad=-2)
		self.sm = sm
                cb = self.fig.colorbar(self.sm)
                cb.set_label(r'$\mathrm{Redshift}$')
	def add_line(self,radbin,tempbin,gmassbin,bfracbin,red):
		self.sub.plot(radbin,tempbin)
		self.sub1.plot(radbin,gmassbin/sum(gmassbin))
		self.sub.axvline(radbin[abs(bfracbin-0.9*(1-0.83120300751)).argmin()],color = 'k',linestyle='--')
	def save(self,hnum,snum,date,ex):
		self.fig.savefig('temp_vs_r_halo%s_%dx_%03d.pdf'%(hnum,ex,snum),Transparent=True)
		self.fig.show()
class cumutempfunc(object):
        """ Cumulative temperature function"""
        def __init__(self,sm):
                self.fig = plt.figure()
                self.sub = self.fig.add_subplot(111)
                self.sm = sm
                self.sub.set_xlabel(r'$\mathrm{Temperature}$ $\mathrm{(K)}$')
                self.sub.set_ylabel(r'$\mathrm{Cumulative}$ $\mathrm{Fraction}$ $(<T)$')
                cb = self.fig.colorbar(self.sm)
                cb.set_label(r'f$_\mathrm{b,max}$')
		self.sub.set_xlim(1e1,1e6)
        def add_line(self,smass,cumutempx,cumutempy,red,hygmass,hnum,tvir,maxbfrac,name):
                self.sub.semilogx(cumutempx,cumutempy/max(cumutempy),linewidth=4,label= name,color=self.sm.to_rgba(maxbfrac))
		if hnum == '007':
                	self.sub.text(0.15, 0.3,'$\mathrm{z}$ = %.3f'%red, ha='center', va='center', transform=self.sub.transAxes, fontsize = 15)
		self.sub.axvline(tvir,color=self.sm.to_rgba(maxbfrac),linestyle = '--')
                #self.sub.text(0.25, 0.2,'$\mathrm{M_{baryon}}$ = %.3e'%hygmass, ha='center', va='center', transform=self.sub.transAxes, fontsize = 15)
        def save(self,hnum,date):
		self.sub.legend(loc=4)
                self.fig.savefig('cumutempfunc_halo%s_%s.pdf'%(hnum,date),transparent=True)
                self.fig.show()
                #plt.close()

def tvir(mvir,rvir):
    G = 4.3e-6 #in kpc/M_sun (km/s)^2
    kb = 1.3806e-26 #in km^2*g/s^2
    mu = 0.59
    mp = 1.6726219e-24
    tv = G*mvir/rvir*(0.5*mu*mp)/kb
    return tv

if __name__ == "__main__":
        hnum =['1084','007']#['32503','20910','1084','007','20192','948','11707','10q','10v']#32503','12596','007','848','796','20192','20910','897','948','1016','32257','10q','10v']
        res = ['_13','_11','_13','_11','_13','_13','_13','_11','_11','_13','_13','_13','_11','_11','_13']
        ver = ['5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16']
        snap = [10]#np.linspace(5,10,6)#[9,184,184,184,184,184,184,184,184,184,184,184,600,600]
        date = time.strftime("%m_%d_%Y")
	clr = ['b','g','r','c','m','y','k']
	name = ['m10k','m10l','m10a','m10g','m10m','m10f','m10h','m10q','m10v']
	ex = 1
	grain = 100
	massp = 1.67262178e-24 #in grams
	gamma = 5./3.
	kb = 1.3806e-26 #in km^2*g/s^2
	G = 4.3e-6 #in kpc/M_sun (km/s)^2
	my_cmap=plt.get_cmap('viridis')
	sm1 = cm.ScalarMappable(cmap=my_cmap,norm=co.Normalize(vmin=0.01, vmax=0.168))
	sm1._A = []
	i = rank
	if hnum[i] == '10q' or hnum[i] == '10v':
		snap = [5]
	for j in np.arange(len(snap)):
                phase_plot = plot_phase()
		pathname = '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output'%(hnum[i],res[i],ver[i])
		sname = '%s/snapdir_%03d/snapshot_%03d.0.hdf5'%(pathname,snap[j],snap[j])
		zarray = pg.readsnap(sname, 'zarray', 'gas')
		he_mf = zarray[:,1] #Helium mass fraction
		y_he = he_mf/(4*(1-he_mf))
		ne = pg.readsnap(sname, 'ne', 'gas') #Electron Abundance
		mu = (1+4*y_he)/(1+y_he+ne)
		mmw = mu*massp #mean molecular weight
		u = pg.readsnap(sname,'u','gas') #specific internal energy in km^2/s^2
		temp = mmw * (gamma-1.)*u/kb #temperature of gas
		hdf = pd.HDFStore('%s/analysis/dataframes/halo%s%s_giz%s_snap%03d.h5'%(pathname,hnum[i],res[i],ver[i],snap[j]))
		xcen = np.float(hdf['props']['halox'])*1000
		ycen = np.float(hdf['props']['haloy'])*1000
		zcen = np.float(hdf['props']['haloz'])*1000
		rvir = np.float(hdf['props']['rvir'])*1000
		red = np.float(hdf['props']['redshift'])
		hinv = 1/.71
		ascale = 1/(red+1)
		den = pg.readsnap(sname, 'rho','gas')*1e10*1.99e33/1.67e-24/(3.806e21)**3*(hinv/((ascale*hinv)**3))
		dpos = pg.readsnap(sname,'pos','dm')
		dm = pg.readsnap(sname,'mass','dm')*1e10
		gpos = pg.readsnap(sname,'pos','gas')
		gm = pg.readsnap(sname,'mass','gas')*1e10
		sm = 0
		dx = dpos[:,0]-xcen
		dy = dpos[:,1]-ycen
		dz = dpos[:,2]-zcen
		gx = gpos[:,0]-xcen
		gy = gpos[:,1]-ycen
		gz = gpos[:,2]-zcen
		dr = np.sqrt(dx**2+dy**2+dz**2)/.71
		gr = np.sqrt(gx**2+gy**2+gz**2)/.71
		ex = 15
		phase_plot.add_plot(sum(dm[dr<ex*rvir]),sum(gm[gr<ex*rvir]),sm,temp[gr<ex*rvir],den[gr<ex*rvir],red,grain)
		phase_plot.save(hnum[i],snap[j],date,ex)
		binz = np.logspace(np.log(.0014),np.log(rvir),48,base=np.e)
		x = np.e**(np.log(binz[1:])-np.log(binz[1]/binz[0])/2)
		tempbin, bin_edge, binnum= stats.binned_statistic(gr, np.log10(temp), statistic='mean', bins=binz)
		tbin = np.logspace(np.log10(min(temp)),np.log10(max(temp)),len(temp)/1e5)
		dmassbin, bin_edge = np.histogram(dr,bins=binz, weights =dm)			
		gmassbin, bin_edge = np.histogram(gr,bins=binz, weights =gm)
		#temp_vs_r_plot = temp_vs_r(sm1)
		#temp_vs_r_plot.add_line(x/rvir*(red+1),tempbin,gmassbin,np.cumsum(gmassbin)/sum(np.cumsum(gmassbin)+np.cumsum(dmassbin)),red)
		#temp_vs_r_plot.save(hnum[i],snap[j],date,ex)
		cumutempall, bin_edge = np.histogram(temp,bins=tbin, weights =np.ones(len(temp)))
		cumutempall = np.cumsum(cumutempall)
		tv = tvir(sum(dm[dr<rvir])+sum(gm[gr<rvir])/.71,rvir/(red+1))
		print tv
		mgas = np.genfromtxt('mvir_out/Halo%s%s_%s_mgas.out'%(hnum[i],res[i],ver[i]))
		mvir = np.genfromtxt('mvir_out/Halo%s%s_%s_mvir.out'%(hnum[i],res[i],ver[i]))
		maxbfrac = max(np.nan_to_num(mgas[:,1]/mvir[:,1]))
		print maxbfrac
		mbary = sum(gm[gr<ex*rvir])
		hdf.close()

	cumutempall = comm.gather(cumutempall, root=0)
	mbary = comm.gather(mbary, root=0)
	tv = comm.gather(tv, root=0)
	maxbfrac = comm.gather(maxbfrac, root=0)
	tbin = comm.gather(tbin[1:], root=0)
        cumutemp_plot = cumutempfunc(sm1)
	if rank == 0:
		for i in np.arange(len(hnum)):
			cumutemp_plot.add_line(sm,tbin[i],cumutempall[i],red,mbary[i],hnum[i],tv[i],maxbfrac[i],name[i])
		cumutemp_plot.save('all',date)



