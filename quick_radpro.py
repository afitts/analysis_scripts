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
#import seaborn as sns

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

##################################################
# Includes scripts to make mass density profiles #
# as well as phase diagrams (rho v temp)         #
##################################################


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
        #ax1.scatter(rhalf,np.sqrt(G*totcum[itotrhalf]/rhalf),marker='s',s = 80,color='%s'%clr[count],label='%s'%(hnum[w]))
        #ax1.scatter(rmax,np.sqrt(G*totcum[irmax]/rmax),marker='D',s = 80,color='%s'%clr[count],label='%s Rmax'%(hnum))
        return massall, gmassall, smassall, x, rhalf,rmax
      else:
        massall = 0.83120300751*massall
	irhalf = (abs(np.sort(diff)-rhalf)).argmin()
	irmax = (abs(np.sort(diff)-rmax)).argmin()
	#ax1.scatter(rhalf,np.sqrt(G*cum[irhalf]/rhalf),marker='^',s = 80,color='%s'%clr[count],label='%s DMO'%(hnum[w]))
	#ax1.scatter(rmax,np.sqrt(G*cum[irmax]/rmax),marker='v',s = 80,color='%s'%clr[count],label='%s DMO Rmax'%(hnum))
        count +=1
        return massall, x

def radpro_df(pathname,sname,hist,dmo,rhalf,rmax,i):
  global grain, clr, count
  massp = 1.67262178e-24 #in grams
  gamma = 5./3.
  kb = 1.3806e-26 #in km^2*g/s^2 
  G = 4.3e-6 #in kpc/M_sun (km/s)^2   
  switch = 0
  hdf = pd.HDFStore('%sdataframes/halo%s%s_giz%s_snap%03d.h5'%(pathname,hnum[w],res[w],ver[w],i))
  red = np.float(hdf['props']['redshift'])
  rvir = np.float(hdf['props']['rvir'])*1000/(red+1)
  rhalf = np.float(hdf['props']['rhalf'])*1000/(red+1)
  vmax = np.float(hdf['props']['vmax'])
  dmass = hdf['particles/dm']['mass'].as_matrix()
  dpos =  hdf['particles/dm']['r'].as_matrix()*1000/(red+1)
  binz = np.logspace(np.log10(.06),np.log10(rvir),grain)
  massall, bin_edge = np.histogram(dpos,bins=binz, weights =dmass) 
  x = 10**(np.log10(binz[1:])-np.log10(binz[1]/binz[0])/2)
  dlogr = np.log(binz[1]/binz[0])
########### Only necessary until all dataframes have been updated to contain cNFW######
#  switch = 0
#  for j in np.arange(16): #Bring together all the AHF files for the snapshot
#    temph = glob.glob(pathname+'ahf_snap%03d/ahf.snap_%03d.%04d.z*.*.AHF_halos'%(i,i,j)) 
#    temph = str(temph).strip('[]').replace("'","")
#    h = np.genfromtxt(temph)
#    if switch == 0:
#      halo = h
#      switch = 1
#    if switch == 1:
#      try:
#        halo = np.vstack((halo,h))
#      except:
#        print "nothing there"
#  for j in np.arange(len(halo)): #Find our halo from all the halos of the snapshot
#    if halo[j,0] == hist[184-i,1]:
#      cNFW = halo[j,42]
#####################################################################################
  cNFW = 0
  if dmo == 0:
    gmass = hdf['particles/gas']['mass'].as_matrix()
    gpos = hdf['particles/gas']['r'].as_matrix()*1000/(red+1)
    smass = hdf['particles/star']['mass'].as_matrix()
    spos = hdf['particles/star']['r'].as_matrix()*1000/(red+1)
    temp = hdf['particles/gas']['temp'].as_matrix()
    den = hdf['particles/gas']['rho'].as_matrix()*1e10*1.99e33/1.67e-24/(3.806e24)**3
    sfr = hdf['particles/gas']['sfr'].as_matrix()
    if hnum[w] == '125961':
     he = 0
     fe = 0
     metal = 0
     fullmetal = 0
    else:
     he = hdf['particles/star']['metal_He'].as_matrix()
     fe = hdf['particles/star']['metal_Fe'].as_matrix()
     metal = hdf['particles/star']['metal_tot'].as_matrix()
     numH = smass*(1-(metal+he))/(1.6733e-24/1.99e33)
     numFe = smass*fe/(9.27e-23/1.99e33)
     avgnum = np.mean(numFe/numH)
     fullmetal = np.log10(numFe/numH)+4.5
     metal = np.log10(avgnum)+4.5
    if hnum[w] == '007':
      den= den*(1e3)**3
    print 'mstar is %f'%(np.log10(sum(smass)))
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

def plot_phase(hnum,hymass,hygmass,hysmass,hyx,temp,den,sfr,date,red):
  global count,grain
  hytot = hymass+hygmass+hysmass
  #print temp, hygmass*1.99e33/(4*3.14159*(hyx*3.09e21)**3)/1.67e-24
  nbins = 100
  #temp = np.nan_to_num(temp)
  H, denedges, tempedges = np.histogram2d(np.log10(den),np.log10(temp),bins=nbins)
  H = np.rot90(H)
  H = np.flipud(H)
 
  # Mask zeros
  Hmasked = np.ma.masked_where(H==0,H)
  #pylab.rcParams['xtick.major.pad']='6'
  #pylab.rcParams['ytick.major.pad']='6'
  fig = plt.figure()
  ax = fig.add_subplot(111)
  #cou = cm.ScalarMappable(norm=co.Normalize(vmin=0, vmax=5e3))
  #cou._A = []
  #cbar = fig.colorbar(cou)
  phase = ax.pcolormesh(denedges,tempedges,np.log10(Hmasked))
  #cbar.set_label(r'$\mathrm{Counts}$')
  cbar = fig.colorbar(phase, ax = ax)
  cbar.set_label('$\mathrm{Log(Counts)}$')
  phase.set_clim(vmin=0,vmax=5)
  #cbar.ScalarMappable()
  #ax.loglog(hygmass*1.99e33/1.67e-24/(4*3.14159/3*(hyx*3.09e21)**3),temp, '%s'%clr[count],linewidth=4,label = '11_13(ff) MFM')#(4*3.14159/3*hyx**3)
  #ax.spines['bottom'].set_linewidth(4)
  #ax.spines['top'].set_linewidth(4)
  #ax.spines['right'].set_linewidth(4)
  #ax.spines['left'].set_linewidth(4)
  #ax.tick_params('both',length=5,width=2,which='minor')
  #ax.tick_params('both',length=10,width=2,which='major')
  #ax.xaxis.set_tick_params(labelsize=20)
  #ax.yaxis.set_tick_params(labelsize=20)
  #plt.xlabel('Radius (kpc)',fontsize=20, labelpad=-10)
  #ax.xaxis.set_label_coords(.48,-.07)
  plt.xlabel(r'$\mathrm{log(\rho)}$ $\mathrm{(n_H/cm^3)}$',fontsize=15, labelpad=-2)
  plt.ylabel(r'$\mathrm{log(Temperature)}$ $\mathrm{(K)}$',fontsize=15, labelpad=-2)
  plt.xlim(-6,2)
  plt.ylim(1,6)
  #plt.ylim(1e-1,2e0)
  plt.legend(loc='best',prop={'size':10})
  ax.text(0.15, 0.3,'$\mathrm{z}$ = %.3f'%red, ha='center', va='center', transform=ax.transAxes, fontsize = 15)
  ax.text(0.25, 0.2,'$\mathrm{M_{baryon}}$ = %.3e'%sum(hysmass+hygmass), ha='center', va='center', transform=ax.transAxes, fontsize = 15)
  ax.text(0.225, 0.1,'$\mathrm{SFR}$ = %.3e'%sum(sfr), ha='center', va='center', transform=ax.transAxes, fontsize = 15)
  #fname = 'm_enc_ratio_halo%s_%dbin_12_16_z0.pdf'%(hnum,grain)
  fname = 'phase_halo%s_%dbin_%03d.png'%(hnum,grain,i-11)
  print i
  plt.savefig(fname,transparent=True)
  #plt.show()
  plt.close()
  return True
def cmike(fname, nofile=False, centered=False, pmass=[0e0], r_stop=None, **args):
    """compute center of mass of a GADGET snapshot.  

    INPUTS:
    fname  -- name of GADGET snapshot file, or, if nofile=True, an array of
    particle positions. 

    OPTIONAL INPUTS:
    nofile=False  -- if True, then assume that 'fname' is actually an Nx3 array
    with positions. 
    centered=False  -- if True, then return Nx3 array in CM frame rather than
    the COM. 
    pmass=[0e0]  -- use a vector containing mass of each particle in computing
    the center of mass 

    **args:  any argument that can be passed to get_pos
    """
    # only do things iteratively for np >= nplim
    nplim=1000

    if nofile == True:
	pos=fname.copy()
    else:
	pos=get_pos(fname, **args)

    if r_stop is None:
        r_stop=1e-10

    ntot=pos.shape[0]
    # if length of pmass is >1, then presumably we are trying to specify the
    # particle masses explicitly:
    if len(pmass) != 1:
	if len(pmass) != ntot:
	    print('if specified, pmass must be an array with '+
                  'the same length as the array of positions')
	    return
	else:
	    # need to reshape particle mass array so that it can be multiplied
	    # by position array in a consistent manner
	    pmass=pmass.reshape(-1,1)

    # get center of mass of array, as defined by (sum(x_i))/N_i
    if len(pmass) == 1: 
	cm0=pos.astype('float64').mean(axis=0).astype('float32')
    else: 
	cm0=(pos.astype('float64')*pmass).mean(axis=0).astype('float32')/pmass.mean()

    if ntot < nplim:
        print('N_p is below limit for iterative CM determination; ' + 
              'computing CM for all particles directly.')
        print 'Center of mass is ', cm0
        cmtemp2=cm0
    else:
        # get radii in new CM frame
        radi=((pos-cm0)**2).sum(axis=1)**0.5

        # get radius of sphere holding 85% of particles, and locations of
        # particles within this sphere:
        tempi=radi.argsort()[int(ntot*0.85)]
        rmax=radi[tempi]
        locs=np.ravel((radi < rmax).nonzero())
        nel=np.size(locs)
        print 'number of elements within ' , rmax , ' kpc is ', nel

        # compute CM of this subset of particles
        if len(pmass) == 1:
            cmn=(pos[locs,:]).mean(axis=0).astype('float32')
        else:
            cmn=(pos[locs,:].astype('float64')*pmass[locs]).mean(axis=0).astype('float32')/(pmass[locs]).mean()

        print cm0; print cmn
        cmtemp1=cm0
        cmtemp2=cmn

        # once radius is below limitc, reduce search radius by 5% rather than 15%.
        rmax1=rmax/4.
        # iteratively determine COM:
        count=0
        # while nel >= 1000 and nel >= ntot/100:
        while nel >= 1000:
            if rmax >= rmax1:
                rmax=rmax*0.85
            else:
                rmax=rmax*0.95

            rad=((pos-cmtemp2)**2).sum(axis=1)**0.5
            locs=np.ravel((rad < rmax).nonzero())
            nel=np.size(locs)
            cmtemp1=cmtemp2
            if len(pmass) == 1:
                cmtemp2=pos[locs,:].astype('float64').mean(axis=0).astype('float32')
            else:
                cmtemp2=(pos[locs,:].astype('float64')*pmass[locs]).mean(axis=0).astype('float32')/pmass[locs].mean()
            count += 1
            if nel <= ntot/500:
                break
            if rmax < r_stop:
                break
        # del radi, tempi, locs, rad
        print 'number of iterations was ', count
        print 'Center of mass is ', cmtemp2, ' within ', rmax
        print 'change from previous iteration was ', cmtemp2-cmtemp1

    if centered == True:
	pos=pos-cmtemp2
	return pos
    else:
	return cmtemp2

def plot_radpro(hnum,hymass,hygmass,hysmass,hyx,dmomass,dmox,date,dm,rhalf,rvir,quench,core,cnfw,clr,dlogr,tempall,cumutempx,cumutempy,cNFW,vmax,fullmetal,metal):
  global clrcount,my_cmap,sm,ex,coreden,vcirc,mhalo,dfrhalf,dfmstar
  G = 4.3e-6 #in kpc/M_sun (km/s)^2   
  pylab.rcParams['xtick.major.pad']='6'
  pylab.rcParams['ytick.major.pad']='6'
  fig = plt.figure()
  ax = fig.add_subplot(111)
  #figm = plt.figure()
  #metalhist = figm.add_subplot(111)
  #figm1 = plt.figure()
  #cumumetal = figm1.add_subplot(111)
  hytot = hymass+hygmass+hysmass
  #print 'here',hyx, hytot/(4*3.14159/3*hyx**3)
  chimin =10000
  #for p in np.logspace(-1,np.log10(50)):
  #  (fit, cmatrix)= opt.curve_fit(psuedoiso,hyx,hytot/(4*3.14159/3*hyx**3),p0=(hytot[0]/(4*3.14159/3*hyx[0]**3),p))
  #  chisq = sum((hytot/(4*3.14159/3*hyx**3)-psuedoiso(hyx,*fit))**2/(.3*hytot/(4*3.14159/3*hyx**3))**2)
  #  if chisq < chimin:
  #    chimin = chisq
  #    bestfit = fit

  #dof = len(hytot)-len(fit)
  #print fit, chisq, dof
  if dm == 0:
    ax.loglog(hyx,hytot/(4*3.14159/3*hyx**3)/dlogr, 'k',linewidth=4,label = '11_13')#(4*3.14159/3*hyx**3)
    #ax.loglog(hyx,psuedoiso(hyx,*bestfit), 'g-',linewidth=2, label='psuedo isothermal')
  #ax.loglog(x9_17,y9_17, 'b',linewidth=4,label = '9_17 MFM')
  #ax.loglog(x6_2,y6_2, 'g',linewidth=4,label = '6_2 MFM')
  #ax.loglog(x11_13,y11_13-g11_13-s11_13, 'b',linewidth=4,label = 'MFM DM')
  #ax.loglog(x11_13,g11_13, 'r',linewidth=4,label = 'MFM gas')
  #ax.loglog(x11_13,s11_13, 'y',linewidth=4,label = 'MFM star')
  #ax.loglog(dmox,dmoy, 'g--',linewidth=4, label = 'DMO (corrected)')
  else:
   if hnum != '125961':
    #if maxstellar<sum(hysmass):  
    # maxstellar = sum(hysmass)
    #if minstellar>sum(hysmass): 
    # minstellar = sum(hysmass)
    #print np.log10(sum(hysmass))
    cumutempfunc.loglog(cumutempx,cumutempy, color=sm.to_rgba(np.log10(sum(hysmass))),linewidth=4)
    avgtempvr.loglog(hyx,10**tempall, color=sm.to_rgba(np.log10(sum(hysmass))),linewidth=4)
    #ax.axhline(y=1,color='k',linestyle='--',linewidth=1.5)
    #ax.loglog(hyx,hytot/(4*3.14159*hyx**3)/dlogr, 'k',linewidth=4,label = '11_13_15')
    if hnum == '007': 
     ress = '_11'
    else:
     ress = '_13'
     if hnum == '2':
      hnum = '848'
    pathname1 = '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz5_12_16_raw_output/analysis/'%(hnum,ress)
    hdf = pd.HDFStore('%sdataframes/halo%s%s_giz5_12_16_snap%03d.h5'%(pathname1,hnum,ress,i))
    red = np.float(hdf['props']['redshift'])
    rhalf1 = np.float(hdf['props']['rhalf'])*1000/(red+1)
    rvir = np.float(hdf['props']['rvir'])*1000/(red+1)
    dmass = hdf['particles/dm']['mass'].as_matrix()
    dpos =  hdf['particles/dm']['r'].as_matrix()*1000/(red+1)
    gmass = hdf['particles/gas']['mass'].as_matrix()
    gpos = hdf['particles/gas']['r'].as_matrix()*1000/(red+1)
    smass = hdf['particles/star']['mass'].as_matrix()
    spos = hdf['particles/star']['r'].as_matrix()*1000/(red+1)
    binz = np.logspace(np.log10(.06),np.log10(rvir),grain)
    hyx = 10**(np.log10(binz[1:])-np.log10(binz[1]/binz[0])/2)
    dlogr = np.log(binz[1]/binz[0])

    dm = hdf['particles/dm']['mass']
    dx = hdf['particles/dm']['x'].as_matrix()
    dy = hdf['particles/dm']['y'].as_matrix()
    dz = hdf['particles/dm']['z'].as_matrix()
    dp = np.column_stack((dx,dy,dz))
    gm = hdf['particles/gas']['mass']
    gx = hdf['particles/gas']['x'].as_matrix()
    gy = hdf['particles/gas']['y'].as_matrix()
    gz = hdf['particles/gas']['z'].as_matrix()
    gp = np.column_stack((gx,gy,gz))
    smm = hdf['particles/star']['mass']
    sx = hdf['particles/star']['x'].as_matrix()
    sy = hdf['particles/star']['y'].as_matrix()
    sz = hdf['particles/star']['z'].as_matrix()
    sp = np.column_stack((sx,sy,sz))
    dsp = np.vstack((dp,sp))
    dsm = np.append(dm,smm)
    dsc = cmike(fname = dsp,nofile=True,pmass=dsm)
    dpos = np.sqrt((dp[:,0]-dsc[0])**2+(dp[:,1]-dsc[1])**2+(dp[:,2]-dsc[2])**2)/.71*1000
    gpos = np.sqrt((gp[:,0]-dsc[0])**2+(gp[:,1]-dsc[1])**2+(gp[:,2]-dsc[2])**2)/.71*1000
    spos = np.sqrt((sp[:,0]-dsc[0])**2+(sp[:,1]-dsc[1])**2+(sp[:,2]-dsc[2])**2)/.71*1000

    massall, bin_edge = np.histogram(dpos,bins=binz, weights =dmass)
    gmassall, bin_edge = np.histogram(gpos,bins=binz, weights =gmass)  
    smassall, bin_edge = np.histogram( spos,bins=binz, weights =smass)
    ax.loglog(hyx,(massall+gmassall+smassall)/(4*3.14159*hyx**3)/dlogr, 'r',linewidth=4,label = '5_12_16')
    ax.text(1,1e8,r'5_12_16 M$_\star=$%.2e'%sum(smassall),fontsize=16)
    ax.axvline(x=rhalf1,color='r',linestyle='--',linewidth=1)
    ax.loglog(dmox,dmomass/(4*3.14159*dmox**3)*0.83120300751/dlogr, color = 'lightgrey',linestyle='--',linewidth=4, label = 'DMO')
    rhohyrhodmovrrhalf.loglog(hyx/rhalf1,(massall+gmassall+smassall)/(4*3.14159*hyx**3)/(dmomass*0.83120300751/(4*3.14159*dmox**3)), color=sm.to_rgba(np.log10(sum(smassall))),linewidth=4,label = r'$\rho_star/\rho_{dm}$') 
    rhohyrhodmovrrhalf.axhline(1,linestyle='--',color='k',linewidth=1)
    #binz = np.linspace(min(fullmetal),max(fullmetal),grain)
    #metalhist.hist(fullmetal,binz)
    #metalhist.set_xlabel('$\mathrm{[Fe/H]}$')
    #metalhist.set_ylabel('# $\mathrm{of}$ $\mathrm{Stars}$')
    #figm.savefig('metalhist_halo%s_%dbin_%s_z0.pdf'%(hnum,grain,date),Transparent=True)
    #metalall, bin_edge = np.histogram(fullmetal,bins=binz)
    #cumumetal.plot((binz[1:]-(binz[1]-binz[0])/2),np.cumsum(metalall)/float(max(np.cumsum(metalall))),linewidth=4)
    #cumumetal.set_xlabel('$\mathrm{[Fe/H]}$')
    #cumumetal.set_ylabel('$\mathrm{Cumulative}$ $\mathrm{Distribution}$')
    #cumumetal.set_ylim(-.05,1.05)
    #figm1.savefig('cumu_metal_halo%s_%dbin_%s_z0.pdf'%(hnum,grain,date),Transparent=True)
    if hnum == '1016':
      trip1.loglog(hyx,hytot/(4*3.14159*hyx**3)/dlogr, 'k',linewidth=2,label = '11_13(ff) MFM')
      trip1.loglog(dmox,dmomass/(4*3.14159*dmox**3)*0.83120300751/dlogr, color = 'darkgrey',linestyle='--',linewidth=2, label = 'DMO')
      trip1.axvline(x=rhalf,color='k',linestyle=':',linewidth=1)
      #trip1.text(1.5,2e8,r'M$_star=$%.2e'%sum(hysmass),fontsize=16)
    if hnum == '11707':
      trip2.loglog(hyx,hytot/(4*3.14159*hyx**3)/dlogr, 'k',linewidth=2,label = '11_13(ff) MFM')
      trip2.loglog(dmox,dmomass/(4*3.14159*dmox**3)*0.83120300751/dlogr, color = 'darkgrey',linestyle='--',linewidth=2, label = 'DMO')
      trip2.axvline(x=rhalf,color='k',linestyle=':',linewidth=1)
      #f.show()
      #trip2.text(1.5,2e8,r'M$_star=$%.2e'%sum(hysmass),fontsize=16)
    if hnum == '20910':
      trip3.loglog(hyx,hytot/(4*3.14159*hyx**3)/dlogr, 'k',linewidth=2,label = '11_13(ff) MFM')
      trip3.loglog(dmox,dmomass/(4*3.14159*dmox**3)*0.83120300751/dlogr, color = 'darkgrey',linestyle='--',linewidth=2, label = 'DMO')
      trip3.axvline(x=rhalf,color='k',linestyle=':',linewidth=1)
      #trip2.set_xlim(.1,10)#(.06,60)
      #trip2.set_ylim(1e5,1e9)
      trip1.set_xlim(.1,20)
      trip1.set_ylim(1e5,1e9)
      #labels = ['%.5f' % float(t.get_text()) for t in ax.get_xticklabels()]
      trip1.set_xticklabels([0.1,0.1,1,10])
      f.subplots_adjust(wspace=0.05)
      plt.setp(trip1, aspect=0.5, adjustable='box-forced')
      plt.setp(trip2, aspect=0.5, adjustable='box-forced')
      plt.setp(trip3, aspect=0.5, adjustable='box-forced')
      f.savefig('trip_rhoden_2_1016_11707_%s.pdf'%date,Transparent=True)
      #trip3.set_xlim(.1,10)#(.06,60)
      #trip3.set_ylim(1e5,1e9)
      #trip3.text(1.5,2e8,r'M$_star=$%.2e'%sum(hysmass),fontsize=16)
    #ax.loglog(hyx,hysmass/(4*3.14159*hyx**3)/dlogr, 'b',linewidth=4,label = '11_13(ff) Star')
    #ax.loglog(hyx,hygmass/(4*3.14159*hyx**3)/dlogr, 'g',linewidth=4,label = '11_13(ff) Gas')
    ax.text(1,2e8,r'11_13_15 M$_\star=$%.2e'%sum(hysmass),fontsize=16)
    #ax.text(1.5,1e8,r'M$_{gas}=$%.2e'%sum(hygmass),fontsize=16)
    #ax.text(1.5,0.3,r'max |$\rho_{dm}$|=%.2e'%(max(hymass/(4*3.14159/3*hyx**3))),fontsize=16)
    #ax.text(1.5,0.2,r'max |$\rho_{*}$|=%.2e'%(max(hysmass/(4*3.14159/3*hyx**3))),fontsize=16)
    #ax.loglog(hyx,psuedoiso(hyx,*bestfit), 'g-',linewidth=2, label='psuedo isothermal') #Core Fitting
    
    ax.axvline(x=rhalf,color='k',linestyle='--',linewidth=1)
  ax.spines['bottom'].set_linewidth(4)
  ax.spines['top'].set_linewidth(4)
  ax.spines['right'].set_linewidth(4)
  ax.spines['left'].set_linewidth(4)
  ax.tick_params('both',length=5,width=2,which='minor')
  ax.tick_params('both',length=10,width=2,which='major')
  ax.xaxis.set_tick_params(labelsize=20)
  ax.yaxis.set_tick_params(labelsize=20)
  ax.set_xlabel(r'$Radius$ $(kpc)$',fontsize=20, labelpad=-10)#'r/R$_{vir}$'
  ax.xaxis.set_label_coords(.48,-.07)
  plt.ylabel(r'$\rho$ $(M_\odot/kpc^3$)',fontsize=20, labelpad=-5)
  #plt.ylabel(r'$\rho_{hydro}/\rho_{dmo}$',fontsize=20, labelpad=-5)
  #plt.title('Halo %s Radial M(<r) (%d bins) Profile'%(hnum,grain))
  ax.set_xlim(.1,10)#(.06,60)
  ax.set_ylim(1e5,1e9)
  #plt.ylim(1e-1,1e1)
  ax.legend(loc=3,prop={'size':10})
  #fname = 'm_enc_ratio_halo%s_%dbin_12_16_z0.pdf'%(hnum,grain)
  fname = 'radden_halo%s_%dbin_%s_z%03d.png'%(hnum,grain,date,i-11)
  fig.savefig(fname,transparent=True)
  #plt.show()
  plt.close()
  #ax1.scatter(rhalf,fit[1],marker='s',s = 80,color='r')
  if quench == 0:
    mrkr = '*'
  else:
    mrkr = 'o'
  w = 0
  while hyx[w] <=.5:
    w+=1
  dw = 0
  while dmox[dw] <=.5:
    dw+=1
  v = 0
  while hyx[v] <=rhalf:
    v+=1
  cmass = np.cumsum(hytot)
  if hnum != '125961':
   rhalfvmstar.scatter(rhalf,sum(hysmass),marker=mrkr,s=90,color=sm.to_rgba(np.log10(sum(hysmass))),label='%s'%hnum)
   trip4.scatter(np.log10(sum(hysmass)),rhalf,s=30,color=sm.to_rgba(np.log10(sum(hysmass))),label='%s'%hnum)
   trip5.scatter(np.log10(sum(hysmass)),np.log10(sum(hytot[:v+1])/sum(hysmass[:v+1])),s=30,color=sm.to_rgba(np.log10(sum(hysmass))),label='%s'%hnum)
   trip6.scatter(np.log10(sum(hysmass)),metal,s=30,color=sm.to_rgba(np.log10(sum(hysmass))),label='%s'%hnum)
   plt.setp(trip4, aspect=(1.6/1.05/1.25), adjustable='box-forced')
   #trip5.set_yticks([25,50,75,100])
   plt.setp(trip5, aspect=1.6/1/1.25, adjustable='box-forced')
   plt.setp(trip6, aspect=(1.6/1.35/1.25), adjustable='box-forced')
   trip4.set_xticks([6,6.5,7.])
   f1.subplots_adjust(hspace=0.06)
   f1.savefig('trip_correlation_mstar_%s.pdf'%date,Transparent=True)
   print 'metal is %.2e, cNFW is %f'%(metal,cNFW)
   coreden = np.append(coreden,hytot[w]/(4*3.14159*hyx[w]**3)/dlogr/(dmomass[w]*0.83120300751/(4*3.14159*hyx[w]**3)/dlogr))
   vcirc = np.append(vcirc, np.sqrt(G*cmass[v]/hyx[v]))
   mhalo = np.append(mhalo,sum(hytot[:v+1])/sum(hysmass[:v+1]))
   dfrhalf = np.append(dfrhalf,rhalf)
   dfmstar = np.append(dfmstar,sum(hysmass))
   coredenvrhalf.scatter(rhalf,(hymass[v]/(4*3.14159*hyx[v]**3)/dlogr),marker=mrkr,s=80,color=sm.to_rgba(np.log10(sum(hysmass))),label='%s'%hnum)
   coredenvmstar.scatter(sum(hysmass),(hymass[v]/(4*3.14159*hyx[v]**3)/dlogr),marker=mrkr,s=80,color=sm.to_rgba(np.log10(sum(hysmass))),label='%s'%hnum)
   vcircvrhalf.scatter(rhalf,np.sqrt(G*cmass[v]/hyx[v]),marker = mrkr, s = 80, color=sm.to_rgba(np.log10(sum(hysmass))), label='%s'%hnum)
   mhalovrhalf.scatter(rhalf,sum(hytot[:w+1])/sum(dmomass[:dw+1]*0.83120300751),marker=mrkr,s=80,color=sm.to_rgba(np.log10(sum(hysmass))),label='%s'%hnum)
   #rho200vrhalf.scatter(sum(hysmass),(dmomass[dw]/(4*3.14159*dmox[dw]**3)/dlogr)*0.83120300751,marker = mrkr,facecolors='none',s=80,edgecolors=sm.to_rgba(np.log10(sum(hysmass))))#,label='%s dmo'%hnum)
   rho200vrhalf.scatter(sum(hysmass),(hymass[w]/(4*3.14159*hyx[w]**3)/dlogr)/((dmomass[dw]*0.83120300751/(4*3.14159*dmox[dw]**3)/dlogr)*0.83120300751),marker = 'o',s=80,color=sm.to_rgba(np.log10(sum(hysmass))),label='%s'%hnum)#(80*sum(hysmass)/5e6)
   vcircvr.plot(hyx,np.sqrt(G*cmass/hyx),linewidth = 2,color=sm.to_rgba(np.log10(sum(hysmass))),label='%s'%hnum)
   vcircvr.scatter(rhalf,np.sqrt(G*cmass[v]/hyx[v]),marker='s', s = 40, color=sm.to_rgba(np.log10(sum(hysmass))))
   mstarvr.semilogx(hyx/rvir,np.cumsum(hysmass)/sum(hysmass),linewidth = 2,color=sm.to_rgba(np.log10(sum(hysmass))),label='%s'%hnum)
   metalvmstar.scatter(np.log10(sum(hysmass)),metal,linewidth = 2,color=sm.to_rgba(np.log10(sum(hysmass))),label='%s'%hnum)
   mdynmdmovmdmo.scatter(sum(dmomass)/1e10,sum(hytot[:w+1])/sum(dmomass[:dw+1]*0.83120300751),linewidth = 2,color=sm.to_rgba(np.log10(sum(hysmass))),label='%s'%hnum)
   mdynmdmovcnfw.scatter(cNFW,sum(hytot[:w+1])/sum(dmomass[:dw+1]*0.83120300751),linewidth = 2,color=sm.to_rgba(np.log10(sum(hysmass))),label='%s'%hnum)
   mdynmdmovvmax.scatter(vmax,sum(hytot[:w+1])/sum(dmomass[:dw+1]*0.83120300751),linewidth = 2,color=sm.to_rgba(np.log10(sum(hysmass))),label='%s'%hnum)
   mdynmstarvmstar.scatter(sum(hysmass),sum(hytot[:v+1])/sum(hysmass[:v+1]),linewidth = 2,color=sm.to_rgba(np.log10(sum(hysmass))),label='%s'%hnum)
  return True

grain = 100
sp = np.linspace(11,184,174)
i=184
h0 = 71
om_l = 0.734
om_m = 0.266
conv = 3.085677581e+19
date = time.strftime("%m_%d_%Y")
hnum = ['007','007','2','2','897','897','1016','1016','796','796','948','948'] #dmo should be placed after hydro run for irhalf variable.
res = ['_11','_11','_13','_13','_13','_13','_13','_13','_13','_13','_13','_13']
ver = ['11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13']
dm = [0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1]
snap = [184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184]
cnfw = [15.39,0,15.24,0,0,0,18.815,0,13.47,0,9.4239,0,8.93069,0,11.7743,0,7.6973,0,16.0826,0,2,2,2,2,2,2,2,2]#125961 16.6265,0 
quench = [1,0,0,0,1,0,0,0,0,0,1,0,2,2,2,2,2,2,2,2]
core = [1,0,0,0,1,0,1,0,1,0,1,0,2,2,2,2,2,2,2,2]


hnum = ['32257','32257','11707','11707','32503','32503','125961','125961','12596','12596','007','007','2','2','897','897','1016','1016','796','796','948','948'] #dmo should be placed after hydro run for irhalf variable.
res = ['_13','_13','_13','_13','_13','_13','_13','_13','_13','_13','_11','_11','_13','_13','_13','_13','_13','_13','_13','_13','_13','_13']
ver = ['11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13']
dm = [0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1]
snap = [184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184]
cnfw = [15.39,0,15.24,0,0,0,18.815,0,13.47,0,9.4239,0,8.93069,0,11.7743,0,7.6973,0,16.0826,0,2,2]#125961 16.6265,0
quench = [1,0,1,0,1,0,1,0,1,0,0,0,1,0,0,0,0,0,1,0,2,2]
core = [2,2,2,2,2,2,2,2,1,0,1,0,0,0,1,0,1,0,1,0,1,0,2,2]
hnum =['11707','11707','32503','32503','12596','12596','007','007','2','2','1016','1016','796','796','20192','20192','20910','20910','897','897','948','948']#['32257','11707','32503','125961','12596','007','2','897','1016','796','948'] #dmo should be placed after hydro run for irhalf variable.
res = ['_13','_13','_13','_13','_13','_13','_11','_11','_13','_13','_13','_13','_13','_13','_13','_13','_13','_13','_13','_13','_13','_13','_13','_13','_13','_13','_13','_13','_13']
ver = ['11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','5_12_16','5_12_16','5_12_16','5_12_16','11_13','11_13','11_13','11_13']#['11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16']
dm = [0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1]
count = 0
clrcount = 0

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
fig2 = plt.figure()
concen = fig2.add_subplot(111)
fig3 = plt.figure()
formt = fig3.add_subplot(111)
fig4 = plt.figure()
smass = fig4.add_subplot(111)
fig5 = plt.figure()
maxsfr = fig5.add_subplot(111)
fig6 = plt.figure()
rhalfplot = fig6.add_subplot(111)
fig7 = plt.figure()
rhalfvmstar = fig7.add_subplot(111)
fig8 = plt.figure()
coredenvmstar = fig8.add_subplot(111)
fig9 = plt.figure()
coredenvrhalf = fig9.add_subplot(111)
fig10 = plt.figure()
vcircvrhalf = fig10.add_subplot(111)
fig11 = plt.figure()
mhalovrhalf = fig11.add_subplot(111)
fig12 = plt.figure()
rho200vrhalf = fig12.add_subplot(111)
fig13 = plt.figure()
vcircvr = fig13.add_subplot(111)
fig14 = plt.figure()
rhohyrhodmovrrhalf = fig14.add_subplot(111)
fig15 = plt.figure()
mstarvr = fig15.add_subplot(111)
fig16 = plt.figure()
metalvmstar = fig16.add_subplot(111)
fig17 = plt.figure()
avgtempvr = fig17.add_subplot(111)
fig18 = plt.figure()
cumutempfunc = fig18.add_subplot(111)
fig19 = plt.figure()
mdynmdmovmdmo = fig19.add_subplot(111)
fig20 = plt.figure()
mdynmdmovcnfw = fig20.add_subplot(111)
fig21 = plt.figure()
mdynmdmovvmax = fig21.add_subplot(111)
fig22 = plt.figure()
mdynmstarvmstar = fig22.add_subplot(111)


clr = ['brown','brown','olive','olive','darkviolet','darkviolet','lime','lime','olive','olive','gray','gray','r','r','b','b','k','k','g','g','y','y','m','m']
my_cmap=plt.get_cmap('plasma')
sm = cm.ScalarMappable(cmap=my_cmap,norm=co.Normalize(vmin=6.0969, vmax=7.176))#vmin=5.874853, vmax=7.11806)) 11_13 limits
sm._A = []
ex = 1

#Initialize triple horizontal radial density profiles
f, (trip1, trip2, trip3) = plt.subplots(1,3,sharex=True,sharey=True)
trip2.set_xlabel('$\mathrm{Radius}$ $\mathrm{(kpc)}$',fontsize=10)
trip1.set_ylabel(r'$\rho$ $(\mathrm{M_\odot/kpc^3})$',fontsize=10)

#Initialize triple vertical M* correlation plots with one colorbar
f1, (trip4, trip5, trip6) = plt.subplots(3,sharex=True)
trip6.set_xlabel('$\mathrm{log}$ $M_\star$',fontsize=12)
trip4.set_ylabel(r'$r_\mathrm{{1/2}}$ $(\mathrm{kpc})$',fontsize=12)
f1.sca(trip4)
plt.yticks(fontsize=10)
trip5.set_ylabel(r'$(M_{\mathrm{dyn}}/M_{\star})(<r_{1/2})$',fontsize=12)
trip5.yaxis.set_label_coords(-.2,.5)
trip6.set_ylabel(r'$\mathrm{[Fe/H]}$',fontsize=12)
trip6.yaxis.set_label_coords(-.2,.5)
f1.sca(trip6)
plt.yticks(fontsize=10)
plt.xticks(fontsize=10)
trip4.set_ylim(.15,1.2)
trip4.set_xlim(np.log10(5e5),np.log10(2e7))
#trip5.set_yscale('log')
trip5.set_ylim(1,2)
trip6.set_ylim(-2.4,-1.05)
f1.sca(trip5)
plt.tick_params(which='minor',labelsize=10,pad=2)
ticks = np.linspace(1,2,2)
plt.yticks(ticks,['10','100'],fontsize=10)#[r'$10^1$',r'$10^2$'],fontsize=10)
minor_ticks=[]
for j in range(2,10):
  minor_ticks.append(1+np.log10(j))
trip5.yaxis.set_minor_locator(FixedLocator(minor_ticks))
minor_labels = ['','30 ','','','60 ','','','']
trip5.yaxis.set_minor_formatter(FixedFormatter(minor_labels))
#f.subplots_adjust(right=0.8)
#cbar_ax = f1.add_axes([0.65, 0.15, 0.025, 0.7])
#cobar = f1.colorbar(sm, cax=cbar_ax)
#cobar.set_label(r'$\mathrm{log}$ $\mathrm{M_star}$')
rhohy = np.zeros(185)
rhodmo = np.zeros(185)
coreden = np.array([])
vcirc = np.array([])
mhalo = np.array([])
dfrhalf = np.array([])
dfmstar = np.array([])


for w in np.arange(len(hnum)):
  if dm[w] == 0:
    pathname = '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/analysis/'%(hnum[w],res[w],ver[w])
    sname = "/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/snapdir_%03d/snapshot_%03d"%(hnum[w],res[w],ver[w],snap[w],snap[w])
    hist = np.genfromtxt(pathname+'halo%smerger_hist.txt'%hnum[w])
    #hymass,hygmass,hysmass,hyx,rhalf,rmax= radpro(pathname,sname,hist,dm[w],0,0)
    hymass,hygmass,hysmass,hyx,rvir,rhalf,rmax, temp, den, sfr,red,incre,tempall,cumutempx,cumutempy,cNFW,vmax,fullmetal,metal = radpro_df(pathname,sname,hist,dm[w],0,0,i)

    #for q in np.arange(len(sp)): #For phase plots
    #  i = sp[q]
    #  try:
    #   pathname = '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/analysis/'%(hnum[w],res[w],ver[w])
    #   sname = "/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/snapdir_%03d/snapshot_%03d"%(hnum[w],res[w],ver[w],i,i)
    #   hymass,hygmass,hysmass,hyx,rvir,rhalf,rmax, temp ,den,sfr,red,incre,tempall,cumutempx,cumutempy,cNFW,vmax,fullmetal,metal = radpro_df(pathname,sname,hist,dm[w],0,0,i)
    #   pathname = '/nobackup/afitts/Gadget-2.0.7/production/gizdm%s%s_raw_output/analysis/'%(hnum[w],res[w])
    #   sname = "/nobackup/afitts/Gadget-2.0.7/production/gizdm%s%s_raw_output/snapdir_%03d/snapshot_%03d"%(hnum[w],res[w],i,i)
    #   print i, hnum[w]
    #   hist = np.genfromtxt(pathname+'halo%smerger_hist.txt'%hnum[w])
    #   dmomass,dmox,dmocNFW,dmovmax = radpro_df(pathname,sname,hist,1,rhalf,rmax,i)
    #   v = 0
    #   while hyx[v] <=.5:
    #     v+=1
    #   rhohy[q+11] = ((hymass[v]+hygmass[v]+hysmass[v])/(4*3.14159*hyx[v]**3)/incre)
    #   dv = 0
    #   while dmox[dv] <=.5:
    #     dv+=1
    #   rhodmo[q+11] = (dmomass[dv]*0.83120300751/(4*3.14159*dmox[dv]**3)/incre)
    #  except Exception,e:
    #   print e, 'no gas'
    #  time = np.genfromtxt('Halo897_13_mstar.out')
    #np.savetxt('Halo%s%s_rhohy_halfrhalf.out'%(hnum[w],res[w]),np.column_stack((time[:,0],rhohy)), fmt= '%.8e')
    #np.savetxt('Halo%s%s_rhodmo_halfrhalf.out'%(hnum[w],res[w]),np.column_stack((time[:,0],rhodmo)), fmt= '%.8e')
    print '%s DONE'%hnum[w]


  else:
    pathname = '/nobackup/afitts/Gadget-2.0.7/production/gizdm%s%s_raw_output/analysis/'%(hnum[w],res[w])
    sname = "/nobackup/afitts/Gadget-2.0.7/production/gizdm%s%s_raw_output/snapdir_%03d/snapshot_%03d"%(hnum[w],res[w],snap[w],snap[w])
    print pathname, hnum[w]
    hist = np.genfromtxt(pathname+'halo%smerger_hist.txt'%hnum[w])
    dmomass,dmox,dmocNFW,dmovmax = radpro_df(pathname,sname,hist,dm[w],rhalf,rmax,i)
    print 'dmocNFW is %f\ncNFW is %f\n'%(dmocNFW,cNFW)
    plot_radpro(hnum[w],hymass,hygmass,hysmass,hyx,dmomass,dmox,date,dm[w],rhalf,rvir,quench[w-1],core[w-1],cnfw[w-1],clr[w-1],incre,tempall,cumutempx,cumutempy,cNFW,dmovmax,fullmetal,metal)


ax1.set_ylabel("$r_{core}$ (kpc)")#"$v_{circ}$ (km/s)")
ax1.set_xlabel("$r_{1/2}$ (kpc)")
#ax1.set_title("Vcirc v rhalf")
#ax1.set_ylim(0,30)
#ax1.legend(loc=4,prop={'size':8},ncol=1)
#plt.show()
#fig1.savefig('vcirc_vs_rhalf_%s.pdf'%date,transparent=True)
fig1.savefig('rcore_vs_rhalf_%s.pdf'%date,transparent=True)

rhalfvmstar.set_xlabel(r'$\mathrm{r}_{1/2}$ $(kpc)$')#$\mathrm{M}_{\mathrm{baryon}}(<5 \mathrm{kpc})$ $(\mathrm{M}_\odot)$')#"$v_{circ}$ (km/s))
rhalfvmstar.set_ylabel(r'$\mathrm{M}_{\star}$ $(\mathrm{M}_\odot)$')#r_{1/2}$ $(kpc)$')
#rhalfvmstar.legend(loc=2,prop={'size':10})
rhalfvmstar.set_xscale('log')
rhalfvmstar.set_yscale('log')
rhalfvmstar.set_xlim(2e-1,2e0)
rhalfvmstar.set_ylim(6e5,2e7)
rhalfvmstar.xaxis.set_label_coords(.48,-.07)
cb7 = fig7.colorbar(sm)
cb7.set_label(r'$\mathrm{log}$ $M_\star$')
#rhalfvmstar.set_xticklabels([])
#rhalfvmstar.set_yticklabels([])
fig7.savefig('mstar_vs_rhalf_%s.pdf'%date,transparent=True)

coredenvmstar.set_ylabel(r'$\rho_{\mathrm{dm}}\mathrm{(<r_{1/2})}$ $(\mathrm{M}_\odot/\mathrm{kpc}^3)$')#"$v_{circ}$ (km/s)")
coredenvmstar.set_xlabel(r'$\mathrm{M}_{\star}$ $(\mathrm{M}_\odot)$')
#coredenvmstar.legend(loc=2,prop={'size':10})
coredenvmstar.set_xscale('log')
coredenvmstar.set_yscale('log')
coredenvrhalf.set_xlim(8e4,2e7)
coredenvrhalf.set_ylim(4e6,1e8)
coredenvrhalf.xaxis.set_label_coords(.48,-.06)
cb8 = fig8.colorbar(sm)
cb8.set_label(r'$\mathrm{logM}_star$')
fig8.savefig('rhodm_vs_mstar_loglog_%s.pdf'%date,transparent=True)

coredenvrhalf.set_ylabel(r'$\rho_{\mathrm{dm}}\mathrm{(<r_{1/2})}$ $(\mathrm{M}_\odot/\mathrm{kpc}^3)$')#"$v_{circ}$ (km/s)")
coredenvrhalf.set_xlabel(r'$\mathrm{r}_{1/2}$ $(kpc)$')
#coredenvrhalf.legend(loc=2,prop={'size':10})
coredenvrhalf.set_xscale('log')
coredenvrhalf.set_yscale('log')
coredenvrhalf.set_xlim(2e-1,2e0)
coredenvrhalf.set_ylim(4e6,1e8)
coredenvrhalf.xaxis.set_label_coords(.48,-.06)
cb9 = fig9.colorbar(sm)
cb9.set_label(r'$\mathrm{log}$ $M_\star$')
fig9.savefig('rhodm_vs_rhalf_loglog_%s.pdf'%date,transparent=True)

vcircvrhalf.set_ylabel(r'$\mathrm{v}_{\mathrm{circ,r}_{1/2}}$ $\mathrm{(km/s)}$')
vcircvrhalf.set_xlabel(r'$\mathrm{r}_{1/2}$ $\mathrm{(kpc)}$')
#vcircvrhalf.legend(loc=2,prop={'size':10})
vcircvrhalf.set_xlim(0,1.2e0)
vcircvrhalf.set_ylim(5,3e1)
cb10 = fig10.colorbar(sm)
cb10.set_label(r'$\mathrm{log}$ $M_\star$')
fig10.savefig('vcircrhalf_vs_rhalf_%s.pdf'%date,transparent=True)

mhalovrhalf.set_ylabel(r'$\mathrm{M}_{\mathrm{hydro}}(<%d kpc)/\mathrm{M}_{\mathrm{dmo,corrected}}(<%d \mathrm{kpc})$'%(ex,ex))# ($M_\odot$)')#"$v_{circ}$ (km/s))
mhalovrhalf.set_xlabel(r'$\mathrm{r}_{1/2}$ $\mathrm{(kpc)}$')#'$log(M_star)$ ($M_\odot$)')
#mhalovrhalf.legend(loc=1,prop={'size':10})
#mhalovrhalf.set_xscale('log')
#mhalovrhalf.set_yscale('log')
mhalovrhalf.set_xlim(1e-1,1.2e0)#set_xlim(1e4,1e8)
mhalovrhalf.set_ylim(4.5e-1,1.1e0)
mhalovrhalf.xaxis.set_label_coords(.48,-.06)
#mhalovrhalf.set_ylim(1e7,3e8)
cb11 = fig11.colorbar(sm)
cb11.set_label(r'$\mathrm{log}$ $M_\star$')
fig11.savefig('mhydro_mdmo_ratio_1kpc_vs_rhalf_loglog_%s.pdf'%date,transparent=True)

mdynmdmovmdmo.set_ylabel(r'$\mathrm{M}_{\mathrm{hydro}}/\mathrm{M}_{\mathrm{dmo,corrected}}(<%d \mathrm{kpc})$'%(ex))# ($M_\odot$)')#"$v_{circ}$ (km/s))
mdynmdmovmdmo.set_xlabel(r'$M_{\mathrm{dmo,corrected}}(<%d \mathrm{kpc})/$$(\mathrm{M_\odot}/10\mathrm{M_\odot})$'%ex)
#mhalovrhalf.legend(loc=1,prop={'size':10})
#mdynmdmovmdmo.set_xscale('log')
#mdynmdmovmdmo.set_yscale('log')
mdynmdmovmdmo.set_xlim(8e-1,1.6e0)#set_xlim(1e4,1e8)
mdynmdmovmdmo.set_ylim(4.5e-1,1.1e0)
mdynmdmovmdmo.xaxis.set_label_coords(.48,-.06)
#mhalovrhalf.set_ylim(1e7,3e8)
cb19 = fig19.colorbar(sm)
cb19.set_label(r'$\mathrm{log}$ $M_\star$')
fig19.savefig('mhydro_mdmo_ratio_1kpc_vs_mdmo_loglog_%s.pdf'%date,transparent=True)

mdynmdmovcnfw.set_ylabel(r'$\mathrm{M}_{\mathrm{hydro}}/\mathrm{M}_{\mathrm{dmo,corrected}}(<%d \mathrm{kpc})$'%(ex))# ($M_\odot$)')#"$v_{circ}$ (km/s))
mdynmdmovcnfw.set_xlabel(r'$\mathrm{cNFW}$')#'$log(M_star)$ ($M_\odot$)')
#mhalovrhalf.legend(loc=1,prop={'size':10})
#mdynmdmovcnfw.set_xscale('log')
#mdynmdmovcnfw.set_yscale('log')
mdynmdmovcnfw.set_xlim(6e0,2e1)#set_xlim(1e4,1e8)
mdynmdmovcnfw.set_ylim(4.5e-1,1.1e0)
mdynmdmovcnfw.xaxis.set_label_coords(.48,-.06)
#mhalovrhalf.set_ylim(1e7,3e8)
cb20 = fig20.colorbar(sm)
cb20.set_label(r'$\mathrm{log}$ $M_\star$')
fig20.savefig('mhydro_mdmo_ratio_1kpc_vs_cNFW_loglog_%s.pdf'%date,transparent=True)

mdynmdmovvmax.set_ylabel(r'$\mathrm{M}_{\mathrm{hydro}}/\mathrm{M}_{\mathrm{dmo,corrected}}(<%d \mathrm{kpc})$'%(ex))# ($M_\odot$)')#"$v_{circ}$ (km/s))
mdynmdmovvmax.set_xlabel(r'$\mathrm{v}_{\mathrm{max,dmo}}$ $\mathrm{(km/s)}$')#'$log(M_star)$ ($M_\odot$)')
#mhalovrhalf.legend(loc=1,prop={'size':10})
#mdynmdmovvmax.set_xscale('log')
#mdynmdmovvmax.set_yscale('log')
mdynmdmovvmax.set_xlim(30,50)#set_xlim(1e4,1e8)
mdynmdmovvmax.set_ylim(4.5e-1,1.1e0)
mdynmdmovvmax.xaxis.set_label_coords(.48,-.06)
#mhalovrhalf.set_ylim(1e7,3e8)
cb21 = fig21.colorbar(sm)
cb21.set_label(r'$\mathrm{log}$ $M_\star$')
fig21.savefig('mhydro_mdmo_ratio_1kpc_vs_vmaxdmo_loglog_%s.pdf'%date,transparent=True)

#rho200vrhalf.set_ylabel(r'$\rho_{dm,200pc}$')
rho200vrhalf.set_ylabel(r'$\mathrm{\rho}_{\mathrm{dm(hydro)}}/\mathrm{\rho}_{\mathrm{dmo(corrected)}}\mathrm{(500pc)}$',fontsize=20)#"$v_{circ}$ (km/s))
rho200vrhalf.set_xlabel(r'$M_{\mathrm{\star}}$ $(M_\odot)$',fontsize=20)#r_{\mathrm{1/2}}$ $\mathrm{(kpc)}$',fontsize=20)
#rho200vrhalf.legend(loc=2,prop={'size':10})
rho200vrhalf.set_xscale('log')
#rho200vrhalf.set_yscale('log')
rho200vrhalf.set_xlim (5e5,2e7) #(1e-1,1.2e0) rhalf
rho200vrhalf.set_ylim(4e-1,1.25e0)#(0,1e0) 200 pc
rho200vrhalf.xaxis.set_label_coords(.48,-.06)
cb12 = fig12.colorbar(sm)
cb12.set_label(r'$\mathrm{log}$ $M_\star$',fontsize=16)
fig12.savefig('rhohydrorhodmo_500pc_vs_mstar_loglog_%s.pdf'%date,transparent=True)

vcircvr.set_ylabel(r'$\mathrm{v}_{\mathrm{circ}}$ $\mathrm{(km/s)}$')
vcircvr.set_xlabel(r'$\mathrm{Radius}$ $\mathrm{(kpc)}$')
#vcircvr.legend(loc=4,prop={'size':10})
vcircvr.set_xlim(1e-1,1.1e0)
vcircvr.set_ylim(5,3e1)
#cb13 = fig13.colorbar(sm)
#cb13.set_label(r'$\mathrm{log}$ $M_\star$')
fig13.savefig('vcirc_profile_%s.pdf'%date,transparent=True)

mstarvr.set_ylabel(r'$\mathrm{M}_star(<\mathrm{r})/\mathrm{M}_{*,\mathrm{tot}}$')
mstarvr.set_xlabel(r'$\mathrm{r}/\mathrm{r}_{\mathrm{vir}}$')
#vcircvr.legend(loc=4,prop={'size':10})
mstarvr.set_xlim(1e-2,1e0)
mstarvr.set_ylim(1e-1,1.2e0)
mstarvr.xaxis.set_label_coords(.48,-.06)
cb15 = fig15.colorbar(sm)
cb15.set_label(r'$\mathrm{log}$ $M_\star$')
fig15.savefig('cumu_mstar_ratio_rvirscale_%s.pdf'%date,transparent=True)

rhohyrhodmovrrhalf.set_xlabel(r'$\mathrm{r}/\mathrm{r}_{1/2}$',fontsize=20, labelpad=-10)#'r/R$_{vir}$'
rhohyrhodmovrrhalf.set_ylabel(r'$\rho_{\mathrm{hydro}}/\rho_{\mathrm{dmo}}$',fontsize=20, labelpad=-5)
rhohyrhodmovrrhalf.set_xlim(.1,10)#(.06,60)
rhohyrhodmovrrhalf.set_ylim(1e-1,2e0)
rhohyrhodmovrrhalf.xaxis.set_label_coords(.48,-.06)
rhohyrhodmovrrhalf.set_xticklabels([0.1,0.1,1,10])
cb14 = fig14.colorbar(sm)
cb14.set_label(r'$\mathrm{log}$ $M_\star$',fontsize=16)
fname = 'radden_ratio_%dbin_%s_z0.pdf'%(grain,date)
fig14.savefig(fname,transparent=True)

metalvmstar.set_xlabel(r'$\mathrm{M}_\star$ $(\mathrm{M}_\odot)$')
metalvmstar.set_ylabel(r'$\mathrm{[Fe/H]}$')
metalvmstar.set_xlim(2,10)
metalvmstar.set_ylim(-3,0)
#cb16 = fig16.colorbar(sm)
#cb16.set_label(r'$\mathrm{log}$ $M_\star$')
metalvmstar.set_xticklabels([])
metalvmstar.set_yticklabels([])
fname = 'metal_vs_mstar_%s.pdf'%date
fig16.savefig(fname,transparent=True)

avgtempvr.set_xlabel(r'$\mathrm{Radius}$ $\mathrm{(kpc)}$')
avgtempvr.set_ylabel(r'$\mathrm{Mean}$ $\mathrm{Temperature}$ $\mathrm{(K)}$')
cb17 = fig17.colorbar(sm)
cb17.set_label(r'$\mathrm{log}$ $M_\star$')
fig17.savefig('avgtemp_vs_r_%s.pdf'%date,transparent=True)

cumutempfunc.set_xlabel(r'$\mathrm{Temperature}$ $\mathrm{(K)}$')
cumutempfunc.set_ylabel(r'$\mathrm{N (<T)}$')
cb18 = fig18.colorbar(sm)
cb18.set_label(r'$\mathrm{log}$ $M_\star$')
fig18.savefig('cumutempfunc_%s.pdf'%date,transparent=True)

mdynmstarvmstar.set_ylabel(r'$\mathrm{M}_{\mathrm{dyn}}/\mathrm{M}_{\mathrm{*}}\mathrm{(<r_{1/2})}$')# ($M_\odot$)')#"$v_{circ}$ (km/s))
mdynmstarvmstar.set_xlabel(r'$\mathrm{M}_{*}$ $(\mathrm{M}_\odot)$')
#mhalovrhalf.legend(loc=1,prop={'size':10})
mdynmstarvmstar.set_xscale('log')
mdynmstarvmstar.set_yscale('log')
mdynmstarvmstar.set_xlim(5e5,2e7)
#mdynmstarvmstar.set_ylim(4.5e-1,1.1e0)
mdynmstarvmstar.xaxis.set_label_coords(.48,-.06)
cb22 = fig22.colorbar(sm)
cb22.set_label(r'$\mathrm{log}$ $M_\star$')
fig22.savefig('mdynmstar_ratio_vs_mstar_loglog_%s.pdf'%date,transparent=True)

#sns.set_style("ticks")
#data = {'mstar':dfmstar,'rhalf':dfrhalf,'coreden':coreden,'mhalo':mhalo}
#df = pd.DataFrame(data,columns=['mstar','rhalf','coreden','mhalo'])
#g = sns.PairGrid(df,x_vars=['mstar','rhalf'],y_vars=['coreden','mhalo','mstar'], size=5)
#g.map(plt.scatter, s=100)
#axes = g.axes
#axes[0,1].set_xlim(.3,1.6)
#axes[2,1].set_xlabel(r'$r_{1/2}$',fontsize=14)
#axes[0,0].set_xlim(5e5,5e7)
#axes[1,0].set_ylim(1e1,1e2)
#axes[2,0].set_ylim(5e5,5e7)
#axes[2,0].set_xlabel(r'$M_\star$',fontsize=14)
#axes[0,0].set_xscale('log')
#axes[1,0].set_yscale('log')
#axes[2,0].set_yscale('log')
#axes[2,0].set_ylabel(r'$M_\star$',fontsize=14)
#axes[0,0].set_ylabel(r'$\rho_{\mathrm{hydro}}/\rho_{DMO}$',fontsize=14)#
#axes[1,0].set_ylabel(r'$M_{\mathrm{dyn}}/M_\star$ $(<r_{1/2})$',fontsize=14)
#g.savefig('scatter.pdf')
#plt.show()

