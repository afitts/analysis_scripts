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
import time
import scipy.stats as stats
import pandas as pd
import scipy.optimize as opt
from mikecm import *
from matplotlib import rcParams
import matplotlib.animation as animation
from matplotlib.ticker import FixedLocator
from matplotlib.ticker import FixedFormatter
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def binner(pathname,hnum,res,ver,dmo,snap,grain=60):
	"""
	Loads particle data from dataframes and bins the particles for radial profiles
	
	INPUTS:
	pathname -- path to snaps or dataframes
	hnum -- halo number
	res -- resolution
	ver -- gizmo version 
	dmo -- 0 if hydro run, 1 if dmo run
	snap -- snap number

	OPTIONAL INPUTS:
	grain -- How many bins for the binning procedure

	OUTPUTS:
	massall -- radial binned mass
	x -- centers of radial bins
	dlogr -- delta log r, THE SPACE BETWEEEEEEN bins"""

	massp = 1.67262178e-24 #in grams
	gamma = 5./3.
	kb = 1.3806e-26 #in km^2*g/s^2 
	G = 4.3e-6 #in kpc/M_sun (km/s)^2   
	switch = 0
	
	### Load in particle data from dataframes ###
	hdf = pd.HDFStore('%sdataframes/halo%s%s_giz%s_snap%03d.h5'%(pathname,hnum,res,ver,snap))
	red = np.float(hdf['props']['redshift'])
	rhalf = np.float(hdf['props']['rhalf'])*1000/(red+1)
	rvir = np.float(hdf['props']['rvir'])*1000/(red+1)
	vmax = np.float(hdf['props']['vmax'])
	dmass = hdf['particles/dm']['mass'].as_matrix()
	dpos =  hdf['particles/dm']['r'].as_matrix()*1000/(red+1)

	### Binning procedure ####
	binz = np.logspace(np.log(.0014),np.log(rvir),12,base=np.e)
	if dmo == 0:
		binz = np.logspace(np.log(.0014),np.log(rvir),100,base=np.e)
	x = np.e**(np.log(binz[1:])-np.log(binz[1]/binz[0])/2)
	dlogr = np.log(binz[1]/binz[0])

	### Recentering procedure ###
	dx = hdf['particles/dm']['x'].as_matrix()
	dy = hdf['particles/dm']['y'].as_matrix()
	dz = hdf['particles/dm']['z'].as_matrix()
	dp = np.column_stack((dx,dy,dz))
	#dsc = mikecm(fname = dp,nofile=True,pmass=dmass)
	#dpos = np.sqrt((dp[:,0]-dsc[0])**2+(dp[:,1]-dsc[1])**2+(dp[:,2]-dsc[2])**2)/.71*1000/(red+1)
	massall, bin_edge = np.histogram(dpos,bins=binz, weights =dmass)



	sortdmp = dpos[np.argsort(dpos)]
	sortdm = dmass[np.argsort(dpos)]
  
	if dmo == 0:
		### Load starmerger files to get merger vs insitu info ###
		insitu_file = np.genfromtxt('/nobackupp8/afitts/analysis_scripts/starmerger_out/Halo%s_starmerger_test.out'%hnum)
		insitu = insitu_file[(insitu_file[:,8]==1)&(insitu_file[:,1]==4)]
		merger = insitu_file[(insitu_file[:,8]==0)&(insitu_file[:,1]==4)]
		gmass = hdf['particles/gas']['mass'].as_matrix()
		gpos = hdf['particles/gas']['r'].as_matrix()*1000/(red+1)
		smass = hdf['particles/star']['mass'].as_matrix()
		spos = hdf['particles/star']['r'].as_matrix()*1000/(red+1)
		svx = insitu[:,5]
		svy = insitu[:,6]
		svz = insitu[:,7]
		msvx = merger[:,5]
		msvy = merger[:,6]
		msvz = merger[:,7]
		xcen = np.float(hdf['props']['halox'])*1000
		ycen = np.float(hdf['props']['haloy'])*1000
		zcen = np.float(hdf['props']['haloz'])*1000
		hvx = np.float(hdf['props']['halo_vx'])
		hvy = np.float(hdf['props']['halo_vy'])
		hvz = np.float(hdf['props']['halo_vz'])
		dm = hdf['particles/dm']['mass']
		dx = hdf['particles/dm']['x'].as_matrix()
		dy = hdf['particles/dm']['y'].as_matrix()
		dz = hdf['particles/dm']['z'].as_matrix()
		ddpos = pg.readsnap('/nobackupp8/afitts/Gadget-2.0.7/production/mfm%s%s_giz5_12_16_raw_output/snapdir_%03d/snapshot_%03d.0.hdf5'%(hnum,res,snap,snap),'pos','dm')
		ddr = np.sqrt((ddpos[:,0]-xcen)**2+(ddpos[:,1]-ycen)**2+(ddpos[:,2]-zcen)**2)/.71
		dvx = hdf['particles/dm']['vx'].as_matrix()
		dvy = hdf['particles/dm']['vy'].as_matrix()
		dvz = hdf['particles/dm']['vz'].as_matrix()
		print 'dvx is ',dvx,dvy,dvz
		dx = hdf['particles/dm']['x'].as_matrix()*1000/(red+1)
		dy = hdf['particles/dm']['y'].as_matrix()*1000/(red+1)
		dz = hdf['particles/dm']['z'].as_matrix()*1000/(red+1)
		dp = np.column_stack((dx,dy,dz))
		gm = hdf['particles/gas']['mass']
		gx = hdf['particles/gas']['x'].as_matrix()
		gy = hdf['particles/gas']['y'].as_matrix()
		gz = hdf['particles/gas']['z'].as_matrix() 
		gp = np.column_stack((gx,gy,gz))
		spid = hdf['particles/star'].index.values
		smm = hdf['particles/star']['mass']
		sx = hdf['particles/star']['x'].as_matrix()*1000/(red+1)
		sy = hdf['particles/star']['y'].as_matrix()*1000/(red+1)
		sz = hdf['particles/star']['z'].as_matrix()*1000/(red+1)
		svx = hdf['particles/star']['vx'].as_matrix()
		svy = hdf['particles/star']['vy'].as_matrix()
		svz = hdf['particles/star']['vz'].as_matrix()
		he = hdf['particles/star']['metal_He'].as_matrix()
		fe = hdf['particles/star']['metal_Fe'].as_matrix()
		metal = hdf['particles/star']['metal_tot'].as_matrix()
		numH = smass*(1-(metal+he))/(1.6733e-24/1.99e33)
		numFe = smass*fe/(9.27e-23/1.99e33)
		meta = numFe/numH
		nofloor = meta#[meta>4.54877795e-09]
		avgnum = np.mean(nofloor)
		fullmetal = np.log10(nofloor)+4.5
		sp = np.column_stack((sx,sy,sz))
		dsp = np.vstack((dp,sp))
		dsm = np.append(dm,smm)
		dsc = mikecm(fname = dsp,nofile=True,pmass=dsm)
		#dpos = np.sqrt((dp[:,0]-dsc[0])**2+(dp[:,1]-dsc[1])**2+(dp[:,2]-dsc[2])**2)/.71*1000/(red+1)
		#gpos = np.sqrt((gp[:,0]-dsc[0])**2+(gp[:,1]-dsc[1])**2+(gp[:,2]-dsc[2])**2)/.71*1000/(red+1)
		#spos = np.sqrt((sp[:,0]-dsc[0])**2+(sp[:,1]-dsc[1])**2+(sp[:,2]-dsc[2])**2)/.71*1000/(red+1)

		massall, bin_edge = np.histogram(dpos,bins=binz, weights =dmass)
		gmassall, bin_edge = np.histogram( gpos,bins=binz, weights =gmass) 
		smassall, bin_edge = np.histogram( spos,bins=binz, weights =smass)
		merger_id = np.zeros(len(spid))
		#for i in np.arange(len(spid)):
		#	merger_id[i] = insitu_file[insitu_file[:,0]==spid[i],8]
		#merger_smassall,bin_edge = np.histogram( spos[merger_id==0],bins=binz, weights =smass[merger_id==0])
		#insitu_smassall,bin_edge = np.histogram( spos[merger_id==1],bins=binz, weights =smass[merger_id==1])
		#merger_smassall,bin_edge = np.histogram( spos[fullmetal<-1.7],bins=binz, weights =smass[fullmetal<-1.7])
		#insitu_smassall,bin_edge = np.histogram( spos[fullmetal>-1.7],bins=binz, weights =smass[fullmetal>-1.7])
		merger_smassall=0
		insitu_smassall=0
		sortdmp = dpos[np.argsort(dpos)]
		sortgp = gpos[np.argsort(gpos)]
		sortsp = spos[np.argsort(spos)]
		sortdm = dmass[np.argsort(dpos)]
		sortgm = gmass[np.argsort(gpos)]
		sortsm = smass[np.argsort(spos)]
		insitu_vdisp = np.zeros(len(binz)-1)
		merger_vdisp = np.zeros(len(binz)-1)
		svrdisp = np.zeros(len(binz)-1)
		svdisp = np.zeros(len(binz)-1)
		svrdisp_insitu = np.zeros(len(binz)-1)
		svrdisp_merger = np.zeros(len(binz)-1)
		dvrdisp = np.zeros(len(binz)-1)
		dvdisp = np.zeros(len(binz)-1)
		sr = spos
		svr = (sx*svx+sy*svy+sz*svz)/sr
		dr = dpos
		dvr = (dx*dvx+dy*dvy+dz*dvz)/dr



		insitu_r = np.sqrt((insitu[:,2]-dsc[0])**2+(insitu[:,3]-dsc[1])**2+(insitu[:,4]-dsc[2])**2)
		merger_r = np.sqrt((merger[:,2]-dsc[0])**2+(merger[:,3]-dsc[1])**2+(merger[:,4]-dsc[2])**2)
		sr_0 = np.sqrt((insitu_file[:,2]-dsc[0])**2+(insitu_file[:,3]-dsc[1])**2+(insitu_file[:,4]-dsc[2])**2)
		svx_0 = insitu_file[:,5]-hvx
		svy_0 = insitu_file[:,6]-hvy
		svz_0 = insitu_file[:,7]-hvz
		svr_0 = ((insitu_file[:,2]-dsc[0])*svx_0+(insitu_file[:,3]-dsc[1])*svy_0+(insitu_file[:,4]-dsc[2])*svz_0)/sr_0
		print 'sr_0 is',min(sr_0),max(sr_0)
		print 'svr_0 is',min(svr_0),max(svr_0)
		for i in np.arange(len(binz)-1):
			insitu_vdisp[i] = np.sqrt(np.std(svx_0[(sr_0>binz[i])&(sr_0<binz[i+1])&(insitu_file[:,8]==1)])**2+np.std(svy_0[(sr_0>binz[i])&(sr_0<binz[i+1])&(insitu_file[:,8]==1)])**2+np.std(svz_0[(sr_0>binz[i])&(sr_0<binz[i+1])&(insitu_file[:,8]==1)])**2)/np.sqrt(3.)
			merger_vdisp[i] = np.sqrt(np.std(svx_0[(sr_0>binz[i])&(sr_0<binz[i+1])&(insitu_file[:,8]==0)])**2+np.std(svy_0[(sr_0>binz[i])&(sr_0<binz[i+1])&(insitu_file[:,8]==0)])**2+np.std(svz_0[(sr_0>binz[i])&(sr_0<binz[i+1])&(insitu_file[:,8]==0)])**2)/np.sqrt(3.)
			svrdisp[i] = np.sqrt(np.std(svr[(sr>binz[i])&(sr<binz[i+1])])**2)
			svrdisp_insitu[i] = np.sqrt(np.std(svr_0[(sr_0>binz[i])&(sr_0<binz[i+1])&(insitu_file[:,8]==1)])**2)
			svrdisp_merger[i] = np.sqrt(np.std(svr_0[(sr_0>binz[i])&(sr_0<binz[i+1])&(insitu_file[:,8]==0)])**2)
			svdisp[i] = np.sqrt(np.std(svx[(sr>binz[i])&(sr<binz[i+1])])**2+np.std(svy[(sr>binz[i])&(sr<binz[i+1])])**2+np.std(svz[(sr>binz[i])&(sr<binz[i+1])])**2)/np.sqrt(3.)
			dvrdisp[i] = np.sqrt(np.std(dvr[(dr>binz[i])&(dr<binz[i+1])])**2)
			dvdisp[i] = np.sqrt(np.std(dvx[(dr>binz[i])&(dr<binz[i+1])])**2+np.std(dvy[(dr>binz[i])&(dr<binz[i+1])])**2+np.std(dvz[(dr>binz[i])&(dr<binz[i+1])])**2)/np.sqrt(3.)
			
		print dvdisp,dvrdisp
		hdf.close()
		return massall, gmassall, smassall, x,merger_smassall,insitu_smassall,np.nan_to_num(merger_vdisp),np.nan_to_num(insitu_vdisp),svrdisp,np.nan_to_num(svrdisp_insitu),np.nan_to_num(svrdisp_merger),svdisp,dvrdisp,dvdisp
	else:
		hdf.close()
		return massall, x,dlogr 

class radpro(object):
	"""Radial density profile """
	def __init__(self):
		self.fig = plt.figure()
		self.sub = self.fig.add_subplot(111)
		self.sub.spines['bottom'].set_linewidth(4)
		self.sub.spines['top'].set_linewidth(4)
		self.sub.spines['right'].set_linewidth(4)
		self.sub.spines['left'].set_linewidth(4)
		self.sub.tick_params('both',length=5,width=2,which='minor')
		self.sub.tick_params('both',length=10,width=2,which='major')
		self.sub.xaxis.set_tick_params(labelsize=20)
		self.sub.yaxis.set_tick_params(labelsize=20)
		self.sub.xaxis.set_label_coords(.48,-.07)
		self.sub.set_xlabel(r'$Radius$ $(kpc)$',fontsize=20, labelpad=-10)
		self.sub.set_ylabel(r'$\rho$ $(M_\odot/kpc^3$)',fontsize=20, labelpad=-5)
		self.sub.set_xlim(.04,60)
		self.sub.set_ylim(1e2,1e6)
		pylab.rcParams['xtick.major.pad']='6'
		pylab.rcParams['ytick.major.pad']='6'
  
	def add_line(self,x,totden,res):
		#if res == '':	# lvl 12
		#	clr = 'r'
		#	lbl = 'LOW '
		#else:		# lvl 13
		#	clr = 'k'
		#	lbl = 'HI '
		if res == '':	# turb diff
			clr = 'r'
			lbl = ''#'in situ'
		else:		# regular
			clr = 'k'
			lbl = 'low metal'#'merger'
		self.sub.loglog(x,totden, color = '%s'%clr,linewidth=4,label = 'Star %s'%lbl)


	def add_dmoline(self,x,totden,res):
		if res == '':	# turb diff
			clr = 'k'
			lbl = 'CDM'#'TD'
		else:		# regular
			clr = 'r'
			lbl = 'SIDM'
		barycorrect = 0.83120300751
		self.sub.loglog(x,totden, color = '%s'%clr,linestyle='',linewidth=4,label = 'DMO %s'%lbl,marker = 's')#self.sub.plot(dmox,dmomass/(4*3.14159*dmox**3)*barycorrect/dlogr, color = 'r',linestyle='',linewidth=4, label = '%sDMO'%lbl,marker='s')


	def save(self,date,hnum):
		self.sub.legend(loc=1)#,prop={'size':10})
		self.fig.savefig('radden_halo%s_metal_%s.pdf'%(hnum,date),transparent=True)
		self.fig.show()

class vdisp(object):
	"""3D velocity dispersion profile """
	def __init__(self):
		self.fig = plt.figure()
		self.sub = self.fig.add_subplot(111)
		self.sub.spines['bottom'].set_linewidth(4)
		self.sub.spines['top'].set_linewidth(4)
		self.sub.spines['right'].set_linewidth(4)
		self.sub.spines['left'].set_linewidth(4)
		self.sub.tick_params('both',length=5,width=2,which='minor')
		self.sub.tick_params('both',length=10,width=2,which='major')
		self.sub.xaxis.set_tick_params(labelsize=20)
		self.sub.yaxis.set_tick_params(labelsize=20)
		self.sub.xaxis.set_label_coords(.48,-.07)
		self.sub.set_xlabel(r'$Radius$ $(kpc)$',fontsize=20, labelpad=-10)
		self.sub.set_ylabel(r'$\sigma_{\star,3D}$',fontsize=20, labelpad=-5)#r'$\sigma_\mathrm{merger,3D}/\sigma_\mathrm{insitu,3D}$',fontsize=20, labelpad=-5)
		self.sub.set_xlim(.1,20)
		#self.sub.set_ylim(1e2,1e6)
		pylab.rcParams['xtick.major.pad']='6'
		pylab.rcParams['ytick.major.pad']='6'
		my_cmap=plt.get_cmap('plasma')
		#self.sm = cm.ScalarMappable(cmap=my_cmap,norm=co.Normalize(vmin=5.004, vmax=7.176))#vmin=5.874853, vmax=7.11806)) 11_13 limits
		self.sm = cm.ScalarMappable(cmap=my_cmap,norm=co.Normalize(vmin=0, vmax=9))
		self.sm._A = []
  
	def add_line(self,x,vdisp,vdisp_std,res,smass):
		if res == '':	# turb diff
			clr = 'k'
			lbl = 'in situ'#$\sigma_\mathrm{r}$'#
		else:		# regular
			clr = 'r'
			lbl = 'merger'#'$\sigma_\mathrm{3D,LOS}$'#
		self.sub.semilogx(x,vdisp, color = '%s'%clr,linewidth=4)#,label = 'Star %s'%lbl)self.sm.to_rgba(smass)
		#self.sub.fill_between(x,vdisp+vdisp_std,vdisp-vdisp_std,color = '%s'%clr,alpha=0.5)

	def save(self,date):
		self.sub.legend(loc=1)#,prop={'size':10})
		self.fig.savefig('vdisp_pro_3D_948_10q_%s.pdf'%(date),transparent=True)
		self.fig.show()

### VV PSEUDO CODE VV ######
#if dmo == 0:
#	dmass,gmass,smass,x = binner(pathname,hnum,res,ver,dmo,snap)
#	totmass = dmass+gmass+smass
#else:
#	totmass,x,dlogr = binner(pathname,hnum,res,ver,dmo,snap)


hnum =['948','10q']#,'12596','007','848','796','20192','20910','897','948','1016','32257','10q','10v']
nummerger = [6,6,4,3,5,5,7,9,1,6,7,5,2,0]
res = ['_13','_11','_13','_11','_13','_13','_13','_13','_13','_13','_13','_13','_11','_11','_13']
#res = ['','','','','','','','','','','','','_11','_11','']
ver = ['5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16']
snum = [184,600,184,184,184,184,184,184,184,184,184,184,600,600]
date = time.strftime("%m_%d_%Y")
vdisp_plot = vdisp()
svrdisp_insitu = np.zeros(60)
svrdisp_merger = np.zeros(60)
i = rank
#for i in np.arange(len(hnum)): for serial
#radpro_plot = radpro()
dmass,gmass,smass,x,merger_smass,insitu_smass,merger_vdisp,insitu_vdisp,svrdisp,svrdisp_insitu,svrdisp_merger,svdisp,dvrdisp,dvdisp = binner('/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz5_12_16_raw_output/analysis/'%(hnum[i],res[i]),hnum[i],res[i],'5_12_16',0,snum[i],60)
#radpro_plot.add_line(x,merger_smass,'_13')
#radpro_plot.add_line(x,insitu_smass,'')
#radpro_plot.save(date,hnum)

#vdisp_plot.add_line(x,merger_vdisp,'_13')

svrdisp_insitu = comm.gather(svrdisp_insitu,root=0)
svrdisp_merger = comm.gather(svrdisp_merger,root=0)
insitu_vdisp = comm.gather(insitu_vdisp,root=0)
merger_vdisp = comm.gather(merger_vdisp,root=0)

if rank == 0:
	svrdisp_inmean = np.mean(([i for i in svrdisp_insitu]),axis = 0)
	svrdisp_instd = np.std(([i for i in svrdisp_insitu]),axis = 0)
	svrdisp_exmean = np.mean(([i for i in svrdisp_merger]),axis = 0)
	svrdisp_exstd = np.std(([i for i in svrdisp_merger]),axis = 0)

	svdisp_inmean = np.mean(([i for i in insitu_vdisp]),axis = 0)
	svdisp_instd = np.std(([i for i in insitu_vdisp]),axis = 0)
	svdisp_exmean = np.mean(([i for i in merger_vdisp]),axis = 0)
	svdisp_exstd = np.std(([i for i in merger_vdisp]),axis = 0)

	for i,e in enumerate(insitu_vdisp):
		vdisp_plot.add_line(x,insitu_vdisp[i],svrdisp_instd,'',nummerger)
                vdisp_plot.add_line(x,merger_vdisp[i],svrdisp_exstd,'0',nummerger)
	#vdisp_plot.add_line(x,svrdisp_inmean,svrdisp_instd,'',nummerger)
	#vdisp_plot.add_line(x,svrdisp_exmean,svrdisp_exstd,'0',nummerger)
	#vdisp_plot.add_line(x,svdisp_inmean,svdisp_instd,'',nummerger)
	#vdisp_plot.add_line(x,svdisp_exmean,svdisp_exstd,'0',nummerger)
	vdisp_plot.save(date)
##x = np.genfromtxt('radden_out/Halo%s_raddenZ13_hydro.out'%hnum)
##res = ''
##radpro_plot.add_line(x[:,0],x[:,1],res)
#x = np.genfromtxt('radden_out/Halo%s_raddenZ13_dmo.out'%hnum)
#radpro_plot.add_dmoline(x[:,0],x[:,1],res)
##x = np.genfromtxt('radden_out/Halo%s_raddenZ13_hydro.out'%hnum1)
##res = '_13'
##radpro_plot.add_line(x[:,0],x[:,1],res)
##radpro_plot.save(date,hnum1)
