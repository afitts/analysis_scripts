import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as co
import time
import scipy.stats as st
import pygadgetreader as pg
from pyd import mikecm
import glob
from matplotlib import rcParams
from matplotlib import gridspec
from scipy import stats

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
rcParams['xtick.labelsize']= '20'
rcParams['ytick.labelsize']= '20'
rcParams['savefig.bbox'] = 'tight'

def behroozi_shmf(mhalo,a):
	z = 1./a-1.
	Mh0 = 10**(11.9852-2.8458*(1-a)-3.0443*np.log(a)-0.5396*z)
	Mstar0true = 10**(np.log10(Mh0)-1.1976-1.4796*(1-a)-0.7104*np.log(a)-0.0044*z)
	Mstar0obs = 10**(np.log10(Mh0)-1.3386-1.6644*(1-a)-0.7104*np.log(a)-0.0044*z)
	m = mhalo/Mh0
	alphaPrime = 2.1463-0.3725*(1-a)-0.0879*np.log(a)-0.0108*z
	betaPrime = 0.4371+0.4399*(1-a)-0.2610*z
	deltaPrime = 0.4017
	gammaPrime = 10**(-1.1953+3.1139*(1-a)-0.7971*z)
	MstarTrue = Mstar0true*np.exp(np.log(10)*gammaPrime*np.exp(-np.log10(m)**2/(2*deltaPrime**2)))/(m**(-alphaPrime)+m**(-betaPrime))
	MstarObs = Mstar0obs*np.exp(np.log(10)*gammaPrime*np.exp(-np.log10(m)**2/(2*deltaPrime**2)))/(m**(-alphaPrime)+m**(-betaPrime))
	return MstarTrue,MstarObs

def contaminator(hnum,res,snap,percent,dmo):
	""" Creates a file that includes all halos/galaxies above a certain contamination threshold"""
	global haloid,hysub,mvir,mstar,mgas,num,snum,gnum,mfrac,rvir,xcen,ycen,zcen,dmohaloid,dmosub,dmomvir,dmonum,dmomfrac,dmodist,dmorvir,dmoxcen,dmoycen,dmozcen,gmvir,gmstar,gmfrac,gsub
	if dmo ==1:
		f = open('halo%s_contam95dmo_halo_catalogue.txt'%(hnum), 'w')
		f.write('(0) HaloID	(1) Subhalo (0=no,#=yes) (2) Mvir	(3) # of particles	(4) Mass frac. of high res	(5) Dist to main progen\n')
		for i in xrange(16):
			ahf = np.genfromtxt('/nobackup/afitts/Gadget-2.0.7/production/gizdm%s%s_raw_output/analysis/ahf_snap%s/ahf.snap_%s.00%02d.z0.000.AHF_halos'%(hnum,res,snap,snap,i))
			ahf = ahf[ahf[:,37]>=percent]
			dmohaloid = np.append(dmohaloid,ahf[:,0])
			dmosub = np.append(dmosub,ahf[:,1])
			dmomvir = np.append(dmomvir,ahf[:,3]/0.71*0.83120300751)
			dmonum = np.append(dmonum,ahf[:,4])
			dmomfrac = np.append(dmomfrac,ahf[:,37])
			dmoxcen = np.append(dmoxcen,ahf[:,5])
			dmoycen = np.append(dmoycen,ahf[:,6])
			dmozcen = np.append(dmozcen,ahf[:,7])
			dmorvir = np.append(dmorvir,ahf[:,11]/.71)
	       
		dmohaloid = dmohaloid[dmomvir.argsort()[::-1]]
		dmosub = dmosub[dmomvir.argsort()[::-1]]
		dmomfrac = dmomfrac[dmomvir.argsort()[::-1]]
		dmoxcen = dmoxcen[dmomvir.argsort()[::-1]]
		dmoycen = dmoycen[dmomvir.argsort()[::-1]]
		dmozcen = dmozcen[dmomvir.argsort()[::-1]]
		dmorvir = dmorvir[dmomvir.argsort()[::-1]]
		dmomvir = dmomvir[dmomvir.argsort()[::-1]]
		xmain = dmoxcen[dmomfrac==1][0]
		ymain = dmoycen[dmomfrac==1][0]
		zmain = dmozcen[dmomfrac==1][0]
		dmodist = np.sqrt((dmoxcen-xmain)**2+(dmoycen-ymain)**2+(dmozcen-zmain)**2)
		cata = np.column_stack((dmohaloid,dmosub,dmomvir,dmonum,dmomfrac,dmodist))
		np.savetxt(f,cata)
		f.close()
		#print hnum, cata[0],cata[(cata[:,3]>0) & (cata[:,1]>0)]
	else:
		f = open('halo%s_%s_contam95_halo_catalogue.txt'%(hnum,snap),'w')#,np.int(100-percent*100)), 'w')
		f.write('(0) HaloID	(1) Subhalo (0=no,#=yes) (2) Mvir	(3) Mstar	(4) Mgas	(5)# of particles	(6)# of stars	(7)# of gas	(8) Mass frac. of high res	(9) Dist to main progen\n')
		for i in xrange(16):
			temph = glob.glob('/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz5_12_16_raw_output/analysis/ahf_snap%s/ahf.snap_%s.%04d.z*.*.AHF_halos'%(hnum,res,snap,snap,i))
			temph = str(temph).strip('[]').replace("'","")
			ahf = np.genfromtxt(temph)

			ahf = ahf[ahf[:,37]>=percent]
			haloid = np.append(haloid,ahf[:,0])
			hysub = np.append(hysub,ahf[:,1])
			mvir = np.append(mvir,ahf[:,3]/0.71)
			mstar = np.append(mstar,ahf[:,64]/0.71)
			mgas = np.append(mgas,ahf[:,44]/0.71)
			num = np.append(num,ahf[:,4])
			snum = np.append(snum,ahf[:,63])
			gnum = np.append(gnum,ahf[:,43])
			mfrac = np.append(mfrac,ahf[:,37])
			xcen = np.append(xcen,ahf[:,5])
			ycen = np.append(ycen,ahf[:,6])
			zcen = np.append(zcen,ahf[:,7])
			rvir = np.append(rvir,ahf[:,11]/.71)
	       
		haloid = haloid[mvir.argsort()[::-1]]
		hysub = hysub[mvir.argsort()[::-1]]
		mstar = mstar[mvir.argsort()[::-1]]
		mgas = mgas[mvir.argsort()[::-1]]
		num = num[mvir.argsort()[::-1]]
		snum = snum[mvir.argsort()[::-1]]
		gnum = gnum[mvir.argsort()[::-1]]
		mfrac = mfrac[mvir.argsort()[::-1]]
		xcen = xcen[mvir.argsort()[::-1]]
		ycen = ycen[mvir.argsort()[::-1]]
		zcen = zcen[mvir.argsort()[::-1]]
		rvir = rvir[mvir.argsort()[::-1]]
		mvir = mvir[mvir.argsort()[::-1]]
		xmain = xcen[mfrac==1][0]
		ymain = ycen[mfrac==1][0]
		zmain = zcen[mfrac==1][0]
		dist = np.sqrt((xcen-xmain)**2+(ycen-ymain)**2+(zcen-zmain)**2)
		gmvir = np.append(gmvir,mvir)
		gmstar = np.append(gmstar,mstar)
		gmfrac = np.append(gmfrac,mfrac)
		gsub = np.append(gsub,hysub)
		cata = np.column_stack((haloid,hysub,mvir,mstar,mgas,num,snum,gnum,mfrac,dist))
		np.savetxt(f,cata)
		f.close()
		print hnum, cata[0],cata[(cata[:,3]>0) & (cata[:,1]>0)]

def denpro(hnum,res,ver,snap,mvir,dmomvir,xcen,dmoxcen,ycen,dmoycen,zcen,dmozcen,cut,x,rho,dmorho):
	""" Finds the radial density profiles for all the halos"""
	switch = [1,0]
	for dmo in switch:
		if dmo == 1:
			pathname = '/nobackup/afitts/Gadget-2.0.7/production/gizdm%s%s_raw_output/snapdir_%s/snapshot_%s.0.hdf5'%(hnum,res,snap,snap)
			boxsize = pg.readheader(pathname, 'boxsize')
			k = 1
			if boxsize < 5000:
				k = 1000
			pos = dpos = pg.readsnap(pathname,'pos','dm')*k
			mass = dmass = pg.readsnap(pathname,'mass','dm')*1e10/.71
			print 'DMO LENGTH IS ',len(dmomvir[dmomvir>cut])
			for i in np.arange(len(dmomvir[dmomvir>cut])):
				ldmorvir = dmorvir[dmomvir>cut]
				lxcen = dmoxcen[dmomvir>cut]
				print ldmorvir[0],pos[0],lxcen[0]
				lycen = dmoycen[dmomvir>cut]
				lzcen = dmozcen[dmomvir>cut]
				binz = np.logspace(np.log(.0014),np.log(ldmorvir[i]),48,base=np.e)
				x.append(np.e**(np.log(binz[1:])-np.log(binz[1]/binz[0])/2))
				xx = np.e**(np.log(binz[1:])-np.log(binz[1]/binz[0])/2)
				diff = np.sqrt((pos[:,0]-lxcen[i])**2+(pos[:,1]-lycen[i])**2+(pos[:,2]-lzcen[i])**2)/.71
				cen_cm = mikecm(fname = dpos[diff<ldmorvir[i]],nofile=True,pmass=dmass[diff<ldmorvir[i]])
				diff = np.sqrt((pos[:,0]-cen_cm[0])**2+(pos[:,1]-cen_cm[1])**2+(pos[:,2]-cen_cm[2])**2)/.71	
				massall, bin_edge = np.histogram(diff[diff<ldmorvir[i]],bins=binz, weights = mass[diff<ldmorvir[i]])
				dlogr = np.log(binz[1]/binz[0])
				dmorho.append(massall/(4*3.14159*xx**3)/dlogr*0.83120300751)
		else:
			pathname = '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/snapdir_%s/snapshot_%s.0.hdf5'%(hnum,res,ver,snap,snap)
			print pathname
			dpos = pg.readsnap(pathname,'pos','dm')
			gpos = pg.readsnap(pathname,'pos','gas')
			spos = pg.readsnap(pathname,'pos','star')
			pos = np.concatenate((dpos,gpos,spos))
			dspos = np.concatenate((dpos,spos))
			dm = pg.readsnap(pathname,'mass','dm')*1e10/.71
			gm = pg.readsnap(pathname,'mass','gas')*1e10/.71
			sm = pg.readsnap(pathname,'mass','star')*1e10/.71
			mass = np.concatenate((dm,gm,sm))
			dsmass = np.concatenate((dm,sm))		
			print 'HYDRO LENGTH IS ',len(mvir[mvir>cut])
			args = np.zeros(len(mvir[mvir>cut]))
			args1 = np.zeros(len(mvir[mvir>cut]))
			for i in np.arange(len(mvir[mvir>cut])):
				ldmorvir = dmorvir[dmomvir>cut]
				lxcen = xcen[mvir>cut]
				lycen = ycen[mvir>cut]
				lzcen = zcen[mvir>cut]
				binz = np.logspace(np.log(.0014),np.log(ldmorvir[i]),48,base=np.e)	
				diff = np.sqrt((dspos[:,0]-lxcen[i])**2+(dspos[:,1]-lycen[i])**2+(dspos[:,2]-lzcen[i])**2)/.71
				cen_cm = mikecm(fname = dspos[diff<ldmorvir[i]],nofile=True,pmass=dsmass[diff<ldmorvir[i]])
				diff = np.sqrt((pos[:,0]-cen_cm[0])**2+(pos[:,1]-cen_cm[1])**2+(pos[:,2]-cen_cm[2])**2)/.71	
				massall, bin_edge = np.histogram(diff[diff<ldmorvir[i]],bins=binz, weights = mass[diff<ldmorvir[i]])
				dlogr = np.log(binz[1]/binz[0])
				rho.append(massall/(4*3.14159*xx**3)/dlogr)
				dist = np.sqrt((dmoxcen[dmomvir>cut][i]-lxcen)**2+(dmoycen[dmomvir>cut][i]-lycen)**2+(dmozcen[dmomvir>cut][i]-lzcen)**2)
				print 'distxyz is ',dist
				dist1 = np.sqrt((dmoxcen[dmomvir>cut][i]-lxcen)**2+(dmoycen[dmomvir>cut][i]-lycen)**2)
				print 'distxy is ',dist1

				args[i] = dist.argmin()
				args1[i] = dist1.argmin()
	print args, args1
	return x, rho, dmorho, args


class rhoratiovr(object):
	""" Rho_hydro/rho_dmo vs r/rhalf """
	def __init__(self,sm):
		self.sm = sm
		self.fig = plt.figure()
		self.sub = self.fig.add_subplot(111)
		self.sub.axhline(1,linestyle='--',color='k',linewidth=1)
		self.sub.set_xlabel(r'$\rm Radius\:(kpc)$',fontsize=16, labelpad=-10)#'r/R$_{vir}$'
		self.sub.set_ylabel(r'$\rho_{\mathrm{hydro}}/\rho_{\mathrm{dmo}}$',fontsize=16, labelpad=-5)
		self.sub.set_xlim(.06,10)
		self.sub.set_ylim(1e-1,2e0)
		self.sub.xaxis.set_label_coords(.48,-.06)
		self.sub.set_xticklabels([0.1,0.1,1,10])
		cb  = self.fig.colorbar(self.sm)
		cb.set_label(r'$\mathrm{log}$ $M_\star$',fontsize=16)

	def add_line(self,x,rho,dmorho,args):
		#self.sub.loglog(x[x<prad]/rhalf,(totmass[x<prad])/(4*3.14159*(x[x<prad]/rhalf)**3)/(dmomass[x<prad]*0.83120300751/(4*3.14159*(dmox[x<prad]/rhalf)**3)), linestyle='--',color=self.sm.to_rgba(np.log10(sum(smass))),linewidth=2)
		for i in np.arange(len(x)):
			if mstar[i] == 0:
				clr = 'k'
				wid = 1
				zo = 0
			else:
				clr = self.sm.to_rgba(np.log10(mstar[i]))
				wid = 3
				zo = 10
			self.sub.loglog(x[i],rho[int(args[i])]/dmorho[i], color=clr,linewidth=wid,zorder=zo)  


	def save(self,date,grain):
		self.fig.savefig('radden_ratio_allhalos_%s_z0.pdf'%(grain,date),transparent=True)
		self.fig.show()




x = []
rho =  []
dmorho = []
#hnum = sys.argv[1]
hnum = ['848','897','1016','796','948','007','11707','12596','32257','32503','20910','20192']#,'10q','10v','1084']
res =  ['_13','_13','_13','_13','_13','_11','_13','_13','_13','_13','_13','_13','_11','_11','_13']
snap = ['184','184','184','184','184','184','184','184','184','184','184','184','600','600','184']
ver = ['5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16']
cut = 1e8
my_cmap=plt.get_cmap('plasma')
sm = cm.ScalarMappable(cmap=my_cmap,norm=co.Normalize(vmin=5.00, vmax=7.158))
sm._A = []
rhoratio_plot = rhoratiovr(sm)
for q in np.arange(2):
	gmvir = np.array([])
	gmstar = np.array([])
	gmfrac = np.array([])
	gsub = np.array([])
	if q == 0:
		hnum = ['848','897','1016','796','948','007','11707','12596','32257','32503','20910','20192','10q','10v','1084']
		res = ['_13','_13','_13','_13','_13','_11','_13','_13','_13','_13','_13','_13','_11','_11','_13']
	else:
		hnum = ['848','897','1016','796','948','007','11707','12596','32257','32503','20910','20192']#,'10q','10v','1084']
		res =  ['','','','','','','','','','','','','','','']
	for i in xrange(len(hnum)):

	  haloid = np.array([])
	  hysub = np.array([])
	  mvir = np.array([])
	  mstar = np.array([])
	  mgas = np.array([])
	  num = np.array([])
	  snum = np.array([])
	  gnum = np.array([])
	  mfrac = np.array([])
	  xcen = np.array([])
	  ycen = np.array([])
	  zcen = np.array([])
	  rvir = np.array([])
	  dmohaloid = np.array([])
	  dmosub = np.array([])
	  dmomvir = np.array([])
	  dmonum = np.array([])
	  dmomfrac = np.array([])
	  dmodist = np.array([])
	  dmoxcen = np.array([])
	  dmoycen = np.array([])
	  dmozcen = np.array([])
	  dmorvir = np.array([])

	 # contaminator(hnum[i],res[i],snap[i],.95,1)
	  contaminator(hnum[i],res[i],snap[i],.95,0)
	  #x, rho, dmorho, args = denpro(hnum[i],res[i],ver[i],snap[i],mvir,dmomvir,xcen,dmoxcen,ycen,dmoycen,zcen,dmozcen,cut,x,rho,dmorho)
	  #rhoratio_plot.add_line(x,rho,dmorho,args)
	if q == 0:
		fig = plt.figure(figsize=(8, 6)) #, ax1 = plt.subplots()
		gs = gridspec.GridSpec(2, 1, height_ratios=[1, 3],hspace=0)
		sub = plt.subplot(gs[1])
		sub.set_xlabel(r'$\mathrm{M_{halo}}$ $(\mathrm{M}_\odot)$')
		sub.set_ylabel(r'$\mathrm{M_{\star}}$ $(\mathrm{M}_\odot)$')
		sub.set_xscale('log')
		sub.set_yscale('log')
		sub.set_xlim(1e8,2e10)
		sub.set_ylim(1e2,2e8)
		#sub.axes.get_yaxis().set_ticklabels(['',r'$10^{2}$',r'$10^{3}$',r'$10^{4}$',r'$10^{5}$',r'$10^{6}$',r'$10^{7}$',''])
		#cb = fig.colorbar(sm)
		#cb.set_label(r'High res mass frac')
		#fig2 = plt.figure()
		sub2 = plt.subplot(gs[0])
		#sub2.set_xlabel(r'$\mathrm{M_{vir}}$ $(\mathrm{M}_\odot)$')
		sub2.set_ylabel(r'$f_\mathrm{dark}$')
		sub2.set_xscale('log')
		#sub2.set_yscale('log')
		sub2.set_xlim(1e8,2e10)
		sub2.set_ylim(-0.2,1.2)
		sub2.axes.get_xaxis().set_ticklabels([])
		sub2.axes.get_yaxis().set_ticklabels(['','0','','0.4','','0.8','','1.2'])

	gggmvir = gmvir[(gmstar<120) & (gmvir>2e7)]
	gmvir = gmvir[gmstar>0]
	gmfrac = gmfrac[gmstar>0]
	gsub = gsub[gmstar>0]
	#mkr = [i for i in gsub[gsub > 0]]
	mkr = np.where(gsub>0,'^','o')

	#gmstar[gmstar == 0] = 120
	    
	gmstar = gmstar[gmstar>0]
	my_cmap=plt.get_cmap('viridis')
	sm = cm.ScalarMappable(cmap=my_cmap,norm=co.Normalize(vmin=0.95, vmax=1))
	sm._A = []
	symbols = [u'\u2193']
	#ggmvir = gmvir[gmstar>120]
	#ggmfrac = gmfrac[gmstar>120]
	#ggmkr = mkr[gmstar>120]
	#ggmstar = gmstar[gmstar>120]
	bin = np.logspace(np.log10(2e7),np.log10(2e10),20)
	x = 10**(np.log10(bin[1:])-np.log10(bin[1]/bin[0])/2)
	gnum,bin_edge = np.histogram(gmvir,bins = bin)
	gggnum,bin_edge = np.histogram(gggmvir,bins = bin)
	print gggnum,gnum
	if q == 0:
		sub2.plot(x,np.asfarray(gggnum)/np.asfarray(gnum+gggnum), color='k')
	else:
		sub2.plot(x,np.asfarray(gggnum)/np.asfarray(gnum+gggnum), color='k',linestyle='--')
	mhalo = np.logspace(np.log10(1e10),np.log10(2e10),10)
	mstar,mstarobs = behroozi_shmf(mhalo,1)
	if q == 0:
		sub.fill_between(mhalo,10**(.5+np.log10(mstar)),10**(-.5+np.log10(mstar)),color = 'lightgray')
		sub.fill_between(mhalo,10**(.2+np.log10(mstar)),10**(-.2+np.log10(mstar)),color = 'gray')
		z0slope, z0intercept, r_value, p_value, std_err = st.linregress(np.log10(mhalo),np.log10(mstar))
		mhalo = np.logspace(np.log10(2e7),np.log10(1e10),100)
		mstar,mstarobs = behroozi_shmf(mhalo,1)
		sub.loglog(mhalo,mstar,color='k',linestyle='--')
		slopehi = 2.6
		slopelow = 1.8
		z0intercepthi = (z0slope-slopehi)*11+z0intercept
		z0interceptlow = (z0slope-slopelow)*11+z0intercept
		sub.fill_between(mhalo,10**z0interceptlow*mhalo**slopelow,10**z0intercepthi*mhalo**slopehi,color = 'cyan')

	#alpha0 = st.linregress(np.log10(gmvir),np.log10(gmstar))
	#alpha1 = st.linregress(np.log10(gmvir[gmvir>1e8]),np.log10(gmstar[gmvir>1e8]))
	#print alpha0[0],alpha1[0]
	if q == 0:
		for i in xrange(len(gmvir)):
		  sub.scatter(gmvir[i],gmstar[i],s =100,marker = mkr[i],edgecolors='k',color='r')#color=sm.to_rgba(gmfrac[i]),label='%s'%hnum)
	else:
		for i in xrange(len(gmvir)):
		  sub.scatter(gmvir[i],gmstar[i],s =100,marker = mkr[i],edgecolors='k',color='b')#color=sm.to_rgba(gmfrac[i]),label='%s'%hnum)
	slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(gmvir),np.log10(gmstar))
	print "r-squared:", r_value**2
	print "slope:",slope
	#sub.loglog(gmvir,10**(slope*np.log10(gmvir)+intercept))
	if q == 0:
		slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(gmvir[gmvir>1e8]),np.log10(gmstar[gmvir>1e8]))
	else:
		gmvir = np.append(gmvir,np.ones(10)*1e10)
		gmstar = np.append(gmstar,np.ones(10)*3e6)
		slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(gmvir[gmvir>3e8]),np.log10(gmstar[gmvir>3e8]))

	print "r-squared:", r_value**2
	print "slope:",slope
	if q == 0:
		sub.loglog(gmvir[gmvir>1e8],10**(slope*np.log10(gmvir[gmvir>1e8])+intercept),color = 'r')
	else:
		sub.loglog(gmvir[gmvir>1e8],10**(slope*np.log10(gmvir[gmvir>1e8])+intercept),color='b')
	if q == 0:
		sub.scatter(gggmvir,np.zeros(len(gggmvir))+120,s =50,marker = 'v',edgecolors='k',color='k',label='%s'%hnum)
	#sub.plot(10**(np.linspace(8,np.log10(4e11))),10**(np.linspace(8,np.log10(4e11))*alpha0[0]+alpha0[1]))
	#sub.plot(10**(np.linspace(8,np.log10(4e11))),10**(np.linspace(8,np.log10(4e11))*alpha1[0]+alpha1[1]))
	#print len(gggmvir)


date = time.strftime("%m_%d_%Y")
fig.savefig('mstar_fracofdarkhalos_v_mvir_%s.pdf'%date,Transparent=True)
#fig2.savefig('FracOfDarkHalos_v_mvir_catalogue_95limit_%s.pdf'%date,Transparent=True)
plt.show()
