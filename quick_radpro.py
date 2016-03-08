import numpy as np
import sys 
import glob
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pylab
import pygadgetreader as pg
import scipy.integrate
import time
import scipy.stats as stats
import pandas as pd
import scipy.optimize as opt

def _a_dot(a, h0, om_m, om_l):
    om_k = 1.0 - om_m - om_l
    return h0 * a * np.sqrt(om_m * (a ** -3) + om_k * (a ** -2) + om_l)


def _a_dot_recip(*args):
    return 1. / _a_dot(*args)

def psuedoiso(r,cenden,rcore):
    return cenden*(1+(r/rcore)**2)**-1

def radpro(pathname,sname,hist,dmo,rhalf,rmax):
  global i, grain, clr, count,conv
  numcores = 16
  massp = 1.67262178e-24 #in grams
  gamma = 5./3.
  kb = 1.3806e-26 #in km^2*g/s^2 
  G = 4.3e-6 #in kpc/M_sun (km/s)^2   
  switch = 0
  for j in np.arange(numcores): #Bring together all the AHF files for the snapshot
    temph = glob.glob(pathname+'ahf_snap%03d/ahf.snap_%03d.%04d.z*.*.AHF_halos'%(i,i,j)) 
    temph = str(temph).strip('[]').replace("'","")
    h = np.genfromtxt(temph)
    temppr = glob.glob(pathname+'ahf_snap%03d/ahf.snap_%03d.%04d.z*.*.AHF_profiles'%(i,i,j))
    temppr = str(temppr).strip('[]').replace("'","")
    p = np.genfromtxt(temppr)
    if switch == 0:
      halo = h
      switch = 1
    if switch == 1:
      try:
        halo = np.vstack((halo,h))
      except:
        print "nothing there"
  for j in np.arange(len(halo)): #Find our halo from all the halos of the snapshot
    if halo[j,0] == hist[184-i,1]:
      mstar = halo[j,64]
      mvir = halo[j,3]
      vmax = halo[j,16]
      xcen = halo[j,5] #Grab relevant data (can always add more later)
      ycen = halo[j,6]
      zcen = halo[j,7]
      rvir = halo[j,11]/.71
      bsize = pg.readheader(sname, 'boxsize')
      if bsize > 25: #mfm007 is in kpc while everything else (even gizdm007) is in Mpc. This standardizes all that.
        kfix = 1
      else:
	kfix = 1000
      mass = pg.readsnap(sname, 'mass', 'dm')*1.e10/.71
      pos = pg.readsnap(sname, 'pos', 'dm')*kfix
      tmp = np.array([0,0])
      diff = np.sqrt((pos[:,0]-xcen)**2+(pos[:,1]-ycen)**2+(pos[:,2]-zcen)**2)/.71
      mass = mass[diff<=rvir]
      diff = diff[diff<=rvir]
      cum = np.cumsum(mass[np.argsort(diff)])

      print xcen,ycen,zcen, rvir,pos
      binz = np.logspace(np.log10(.06),np.log10(rvir),grain)
      massall, bin_edge = np.histogram( diff,bins=binz, weights =mass) 
      x = 10**(np.log10(binz[1:])-np.log10(binz[1]/binz[0])/2)
      if dmo == 0:
        rmax = halo[j,12]/.71
        gmass = pg.readsnap(sname, 'mass', 'gas')*1.e10/.71
        smass = pg.readsnap(sname, 'mass', 'star')*1.e10/.71
        gpos = pg.readsnap(sname, 'pos', 'gas')*kfix
        spos = pg.readsnap(sname, 'pos', 'star')*kfix
        sdiff = np.sqrt((spos[:,0]-xcen)**2+(spos[:,1]-ycen)**2+(spos[:,2]-zcen)**2)/.71
        gdiff = np.sqrt((gpos[:,0]-xcen)**2+(gpos[:,1]-ycen)**2+(gpos[:,2]-zcen)**2)/.71
	smass = smass[sdiff<=rvir]
	gmass = gmass[gdiff<=rvir]
	sdiff = sdiff[sdiff<=rvir]
	gdiff = gdiff[gdiff<=rvir]
	scum = np.cumsum(smass[np.argsort(sdiff)])
	gcum = np.cumsum(gmass[np.argsort(gdiff)])
	totdiff = np.append(np.append(diff,sdiff),gdiff)
	totmass = np.append(np.append(mass,smass),gmass)
	totcum = np.cumsum(totmass[np.argsort(totdiff)])
	irhalf = abs(scum-sum(smass)/2).argmin()
	rhalf = sdiff[np.argsort(sdiff)][irhalf]
	itotrhalf = abs(np.sort(totdiff)-rhalf).argmin()
	print 'rhalf is %f'%rhalf
	irmax = np.sqrt(G*totcum/np.sort(totdiff)).argmax()
	rmax = totdiff[np.argsort(totdiff)][irmax]
	print 'rmax is %f'%rmax


	#print len(sdiff), binz, len(smass)
        gmassall, bin_edge = np.histogram( gdiff,bins=binz, weights =gmass) 
        smassall, bin_edge = np.histogram( sdiff,bins=binz, weights =smass)
	#irhalf = (abs(bin_edge-max(bin_edge[np.cumsum(smassall)<=sum(smassall)/2]))).argmin()

        #print bin_edge[np.cumsum(smassall)<=sum(smassall)/2],rvir,sum(smassall)#sum(smassall), np.cumsum(smassall) #x[irhalf],np.sqrt(G*sum(massall[:irhalf+1])/x[irhalf])
        ax1.scatter(rhalf,np.sqrt(G*totcum[itotrhalf]/rhalf),marker='s',s = 80,color='%s'%clr[count],label='%s'%(hnum[w]))
        #ax1.scatter(rmax,np.sqrt(G*totcum[irmax]/rmax),marker='D',s = 80,color='%s'%clr[count],label='%s Rmax'%(hnum))
        return massall, gmassall, smassall, x, rhalf,rmax
      else:
        massall = 0.83120300751*massall
	irhalf = (abs(np.sort(diff)-rhalf)).argmin()
	irmax = (abs(np.sort(diff)-rmax)).argmin()
	print np.sqrt(G*cum[irhalf]/rhalf), np.sqrt(G*cum[irmax]/rmax)
	ax1.scatter(rhalf,np.sqrt(G*cum[irhalf]/rhalf),marker='^',s = 80,color='%s'%clr[count],label='%s DMO'%(hnum[w]))
	#ax1.scatter(rmax,np.sqrt(G*cum[irmax]/rmax),marker='v',s = 80,color='%s'%clr[count],label='%s DMO Rmax'%(hnum))
        count +=1
        return massall, x

def radpro_df(pathname,sname,hist,dmo,rhalf,rmax):
  global i, grain, clr, count
  massp = 1.67262178e-24 #in grams
  gamma = 5./3.
  kb = 1.3806e-26 #in km^2*g/s^2 
  G = 4.3e-6 #in kpc/M_sun (km/s)^2   
  switch = 0
  hdf = pd.HDFStore('%sdataframes/halo%s%s_giz%s_snap%03d.h5'%(pathname,hnum[w],res[w],ver[w],i))
  red = np.float(hdf['props']['redshift'])
  rvir = np.float(hdf['props']['rvir'])*1000/(red+1)
  rhalf = np.float(hdf['props']['rhalf'])*1000/(red+1)
  dmass = hdf['particles/dm']['mass'].as_matrix()
  dpos =  hdf['particles/dm']['r'].as_matrix()*1000/(red+1)
  binz = np.logspace(np.log10(.06),np.log10(rvir),grain)
  print binz
  massall, bin_edge = np.histogram(dpos,bins=binz, weights =dmass) 
  x = 10**(np.log10(binz[1:])-np.log10(binz[1]/binz[0])/2)
  if dmo == 0:
    gmass = hdf['particles/gas']['mass'].as_matrix()
    gpos = hdf['particles/gas']['r'].as_matrix()*1000/(red+1)
    smass = hdf['particles/star']['mass'].as_matrix()
    spos = hdf['particles/star']['r'].as_matrix()*1000/(red+1)
    temp = hdf['particles/gas']['temp'].as_matrix()
    gmassall, bin_edge = np.histogram( gpos,bins=binz, weights =gmass) 
    smassall, bin_edge = np.histogram( spos,bins=binz, weights =smass)
    print gpos
    tempall, bin_edge, binnum= stats.binned_statistic(gpos, temp, statistic='mean', bins=binz)
    return massall, gmassall, smassall, x, rhalf,rmax, tempall, red
  else:
    count += 1
    return massall, x

def plot_phase(hnum,hymass,hygmass,hysmass,hyx,temp,date, red):
  global count,grain
  hytot = hymass+hygmass+hysmass
  print temp, hygmass*1.99e33/(4*3.14159/3*(hyx*3.09e21)**3)/1.67e-24
  pylab.rcParams['xtick.major.pad']='6'
  pylab.rcParams['ytick.major.pad']='6'
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.loglog(hygmass*1.99e33/1.67e-24/(4*3.14159/3*(hyx*3.09e21)**3),temp, '%s'%clr[count],linewidth=4,label = '11_13(ff) MFM')#(4*3.14159/3*hyx**3)
  ax.spines['bottom'].set_linewidth(4)
  ax.spines['top'].set_linewidth(4)
  ax.spines['right'].set_linewidth(4)
  ax.spines['left'].set_linewidth(4)
  ax.tick_params('both',length=5,width=2,which='minor')
  ax.tick_params('both',length=10,width=2,which='major')
  ax.xaxis.set_tick_params(labelsize=20)
  ax.yaxis.set_tick_params(labelsize=20)
  plt.xlabel('Radius (kpc)',fontsize=20, labelpad=-10)
  ax.xaxis.set_label_coords(.48,-.07)
  plt.xlabel(r'$\rho (n_H$/cm$^3$)',fontsize=20, labelpad=-5)
  plt.ylabel(r'Temperature (K)',fontsize=20, labelpad=-5)
  #plt.xlim(.06,60)
  #plt.ylim(1e4,1e9)
  #plt.ylim(1e-1,2e0)
  plt.legend(loc=3,prop={'size':10})
  #fname = 'm_enc_ratio_halo%s_%dbin_12_16_z0.pdf'%(hnum,grain)
  fname = 'phase_halo%s_%dbin_%s_z%f.pdf'%(hnum,grain,date,red)
  plt.savefig(fname,transparent=True)
  plt.show()
  plt.close()
  return True

def plot_radpro(hnum,hymass,hygmass,hysmass,hyx,dmomass,dmox,date,dm):

  pylab.rcParams['xtick.major.pad']='6'
  pylab.rcParams['ytick.major.pad']='6'
  fig = plt.figure()
  ax = fig.add_subplot(111)
  hytot = hymass+hygmass+hysmass
  chimin =1000
  for p in np.logspace(-1,np.log10(50)):
    (fit, cmatrix)= opt.curve_fit(psuedoiso,hyx,hytot/(4*3.14159/3*hyx**3),p0=(hytot[0]/(4*3.14159/3*hyx[0]**3),p))
    chisq = sum((hytot/(4*3.14159/3*hyx**3)-psuedoiso(hyx,*fit))**2/(.3*hytot/(4*3.14159/3*hyx**3))**2)
    if chisq < chimin:
      chimin = chisq
      bestfit = fit

  dof = len(hytot)-len(fit)
  print fit, chisq, dof
  if dm == 0:
    ax.loglog(hyx,hytot/(4*3.14159/3*hyx**3), 'k',linewidth=4,label = '11_13(ff) MFM')#(4*3.14159/3*hyx**3)
    ax.loglog(hyx,psuedoiso(hyx,*bestfit), 'g-',linewidth=2, label='psuedo isothermal')
  #ax.loglog(x9_17,y9_17, 'b',linewidth=4,label = '9_17 MFM')
  #ax.loglog(x6_2,y6_2, 'g',linewidth=4,label = '6_2 MFM')
  #ax.loglog(x11_13,y11_13-g11_13-s11_13, 'b',linewidth=4,label = 'MFM DM')
  #ax.loglog(x11_13,g11_13, 'r',linewidth=4,label = 'MFM gas')
  #ax.loglog(x11_13,s11_13, 'y',linewidth=4,label = 'MFM star')
  #ax.loglog(dmox,dmoy, 'g--',linewidth=4, label = 'DMO (corrected)')
  else:
    ax.loglog(hyx,hytot/(4*3.14159/3*hyx**3), 'k',linewidth=4,label = '11_13(ff) MFM')
    ax.loglog(dmox,dmomass/(4*3.14159/3*dmox**3), 'r--',linewidth=4, label = 'DMO')
    ax.loglog(hyx,psuedoiso(hyx,*bestfit), 'g-',linewidth=2, label='psuedo isothermal')
  ax.spines['bottom'].set_linewidth(4)
  ax.spines['top'].set_linewidth(4)
  ax.spines['right'].set_linewidth(4)
  ax.spines['left'].set_linewidth(4)
  ax.tick_params('both',length=5,width=2,which='minor')
  ax.tick_params('both',length=10,width=2,which='major')
  ax.xaxis.set_tick_params(labelsize=20)
  ax.yaxis.set_tick_params(labelsize=20)
  plt.xlabel('Radius (kpc)',fontsize=20, labelpad=-10)
  ax.xaxis.set_label_coords(.48,-.07)
  plt.ylabel(r'$\rho (M_\odot$/kpc$^3$)',fontsize=20, labelpad=-5)
  #plt.ylabel(r'$M_\odot/M_{DM,corrected}$',fontsize=20, labelpad=-5)
  plt.title('Halo %s Radial M(<r) (%d bins) Profile'%(hnum,grain))
  plt.xlim(.06,60)
  plt.ylim(1e4,1e9)
  #plt.ylim(1e-1,2e0)
  plt.legend(loc=3,prop={'size':10})
  #fname = 'm_enc_ratio_halo%s_%dbin_12_16_z0.pdf'%(hnum,grain)
  fname = 'radden_halo%s_%dbin_%s_z0.pdf'%(hnum,grain,date)
  plt.savefig(fname,transparent=True)
  plt.show()
  plt.close()
  return True

grain = 30
i = 184
h0 = 71
om_l = 0.734
om_m = 0.266
conv = 3.085677581e+19
date = time.strftime("%m_%d_%Y")
hnum = ['007','007','2','2','897','897','1016','1016','796','796','948','948'] #dmo should be placed after hydro run for irhalf variable.
res = ['_11','_11','_13','_13','_13','_13','_13','_13','_13','_13','_13','_13']
ver = ['11_13','11_13','11_13','11_13','11_13discover','11_13','11_13discover','11_13','11_13','11_13','11_13','11_13']
dm = [0,1,0,1,0,1,0,1,0,1,0,1]
snap = [184,184,184,184,184,184,184,184,184,184,184,184]
count = 0
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
clr = ['r','b','k','g','y','m']
for w in np.arange(len(hnum)):
  if dm[w] == 0:
    pathname = '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/analysis/'%(hnum[w],res[w],ver[w])
    sname = "/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/snapdir_%03d/snapshot_%03d"%(hnum[w],res[w],ver[w],snap[w],snap[w])
    hist = np.genfromtxt(pathname+'halo%smerger_hist.txt'%hnum[w])
    hymass,hygmass,hysmass,hyx,rhalf,rmax= radpro(pathname,sname,hist,dm[w],0,0)
    hymass,hygmass2,hysmass,hyx2,rhalf,rmax, temp = radpro_df(pathname,sname,hist,dm[w],0,0)
    plot_phase(hnum[w],hymass,hygmass,hysmass,hyx,temp,date)
    i = 
    hymass,hygmass2,hysmass,hyx2,rhalf,rmax, temp = radpro_df(pathname,sname,hist,dm[w],0,0)
    plot_phase(hnum[w],hymass,hygmass,hysmass,hyx,temp,date)

  else:
    pathname = '/nobackup/afitts/Gadget-2.0.7/production/gizdm%s%s_raw_output/analysis/'%(hnum[w],res[w])
    sname = "/nobackup/afitts/Gadget-2.0.7/production/gizdm%s%s_raw_output/snapdir_%03d/snapshot_%03d"%(hnum[w],res[w],snap[w],snap[w])
    hist = np.genfromtxt(pathname+'halo%smerger_hist.txt'%hnum[w])
    dmomass,dmox = radpro(pathname,sname,hist,dm[w],rhalf,rmax)
    #plot_radpro(hnum[w],hymass,hygmass,hysmass,hyx,dmomass,dmox,date,dm[w])


ax1.set_ylabel("$v_{circ}$ (km/s)")
ax1.set_xlabel("$r_{1/2}$ (kpc)")
ax1.set_title("Vcirc v rhalf")
ax1.set_ylim(0,30)
ax1.legend(loc=4,prop={'size':8},ncol=1)
plt.show()
fig1.savefig('vcirc_vs_rhalf_%s.pdf'%date,transparent=True)


