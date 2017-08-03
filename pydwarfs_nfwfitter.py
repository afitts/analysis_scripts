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
import matplotlib.animation as animation
from matplotlib.ticker import FixedLocator
from matplotlib.ticker import FixedFormatter
#import seaborn as sns

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
rcParams['xtick.labelsize']= '16'
rcParams['ytick.labelsize']= '16'
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

def nfw(r,mvir,c):
    gc = 1./(np.log(1.+c)-(c/(1.+c)))
    rhocrit = 139.94 #M_sun/kpc^3
    rvir = (3/4.*mvir*(1./(np.pi*96.45*rhocrit)))**(1/3.)
    rs = rvir/c
    #return mvir*gc*(np.log(1+r/rs)-r/rs*(1+r/rs)**-1)
    rho0 = rhocrit*96.45/3*gc*c**3.
    return np.log(rho0*(r/rs)**-1*(1+(r/rs))**-2)#np.log(rho0)-np.log(r/rs)-2*np.log(1+(r/rs))

def corenfw(r,mvir,c,n,rc):
    G = 4.3e-6 #in kpc/M_sun (km/s)^2
    gc = 1./(np.log(1.+c)-(c/(1.+c)))
    rhocrit = 139.94 #M_sun/kpc^3
    rvir = (3/4.*mvir*(1./(np.pi*96.45*rhocrit)))**(1/3.)
    rs = rvir/c
    rho0 = rhocrit*96.45/3*gc*c**3.
    nfw = rho0*(r/rs)**-1*(1+(r/rs))**-2
    mnfw = lambda r:mvir*gc*(np.log(1+r/rs)-r/rs*(1+r/rs)**-1)
    mnfw_rs = mnfw(rs)
    #rc = eta*rhalf
    f = np.tanh(r/rc)
    #tdyn = 2*np.pi*np.sqrt(rs**3/(G*mnfw_rs))
    #q = kappa_tsf/tdyn
    #n = np.tanh(q)
    print 'n is ',n
    return np.log(f**n*nfw+n*f**(n-1)*(1-f**2)/(4*np.pi*r**2*rc)*mnfw(r))

def einasto(r,mvir,c):
    rhocrit = 139.94 #M_sun/kpc^3
    alpha = 0.17
    rvir = (3/4.*mvir*(1./(np.pi*96.45*rhocrit)))**(1/3.)
    rs = rvir/c
    h = rs/((2/alpha)**(1/alpha))
    s = ((2/alpha)**(1/alpha))*r/rs
    return np.log(mvir/(4*np.pi*h**3/alpha*sp.gamma(3./alpha))*np.exp(-s**alpha))

def einasto1(r,mvir,c):
    rhocrit = 139.94 #M_sun/kpc^3
    alpha = 0.17
    rvir = (3/4.*mvir*(1./(np.pi*96.45*rhocrit)))**(1/3.)
    x = r*c/rvir
    M = mvir*(1-sp.gammaincc(3./alpha,((2*x)**alpha)/alpha))
    return np.log(M/(4*np.pi*r**3))

def einasto2(r,mvir,c):
    rhocrit = 139.94 #M_sun/kpc^3
    alpha = 0.17
    rvir = (3/4.*mvir*(1./(np.pi*96.45*rhocrit)))**(1/3.)
    x = r*c/rvir
    M = mvir*(sp.gamma(3./alpha)-sp.gamma(3./alpha)*sp.gammaincc(3./alpha,((2*x)**alpha)/alpha))/sp.gamma(3./alpha)
    return np.log(M)

def radpro_AHF(pathname,sname,hist,dmo,rhalf,rmax):
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

def radpro_df(pathname,hnum,res,ver,hist,dmo,rhalf,rmax,i,grain,dmover):
  massp = 1.67262178e-24 #in grams
  gamma = 5./3.
  kb = 1.3806e-26 #in km^2*g/s^2 
  G = 4.3e-6 #in kpc/M_sun (km/s)^2   
  switch = 0
  print '%sdataframes/halo%s%s_giz%s_snap%03d.h5'%(pathname,hnum,res,ver,i)
  hdf = pd.HDFStore('%sdataframes/halo%s%s_giz%s_snap%03d.h5'%(pathname,hnum,res,ver,i))
  red = np.float(hdf['props']['redshift'])
  rvir = np.float(hdf['props']['rvir'])*1000/(red+1)
  if dmo == 0:
    if hnum == '848':
      hdf1 = pd.HDFStore('/nobackup/afitts/Gadget-2.0.7/production/gizdm%s%s_raw_output/analysis/dataframes/halo%s%s_giz%s_snap%03d.h5'%('2',res,'2',res,dmover,i))
    else:
      hdf1 = pd.HDFStore('/nobackup/afitts/Gadget-2.0.7/production/gizdm%s%s_raw_output/analysis/dataframes/halo%s%s_giz%s_snap%03d.h5'%(hnum,res,hnum,res,dmover,i))
    rvir = np.float(hdf1['props']['rvir'])*1000/(red+1)  
  rhalf = np.float(hdf['props']['rhalf'])*1000/(red+1)
  vmax = np.float(hdf['props']['vmax'])
  dmass = hdf['particles/dm']['mass'].as_matrix()
  dpos =  hdf['particles/dm']['r'].as_matrix()*1000/(red+1)
  binz = np.logspace(np.log(.0014),np.log(rvir),48,base=np.e)
  x = np.e**(np.log(binz[1:])-np.log(binz[1]/binz[0])/2)
  print binz
  dx = hdf['particles/dm']['x'].as_matrix()
  dy = hdf['particles/dm']['y'].as_matrix()
  dz = hdf['particles/dm']['z'].as_matrix()
  dp = np.column_stack((dx,dy,dz))
  dsc = mikecm(fname = dp,nofile=True,pmass=dmass)
  dpos = np.sqrt((dp[:,0]-dsc[0])**2+(dp[:,1]-dsc[1])**2+(dp[:,2]-dsc[2])**2)/.71*1000/(red+1)
  massall, bin_edge = np.histogram(dpos,bins=binz, weights =dmass)
  dmnum, bin_edge = np.histogram(dpos,bins=binz)
  dlogr = np.log(binz[1]/binz[0])
  cNFW = np.float(hdf['props']['cNFW'])
  sortdmp = dpos[np.argsort(dpos)]
  sortdm = dmass[np.argsort(dpos)]
  
  if dmo == 0:
    gmass = hdf['particles/gas']['mass'].as_matrix()
    gpos = hdf['particles/gas']['r'].as_matrix()*1000/(red+1)
    smass = hdf['particles/star']['mass'].as_matrix()
    spos = hdf['particles/star']['r'].as_matrix()*1000/(red+1)
    temp = hdf['particles/gas']['temp'].as_matrix()
    den = hdf['particles/gas']['rho'].as_matrix()*1e10*1.99e33/1.67e-24/(3.806e24)**3
    sfr = hdf['particles/gas']['sfr'].as_matrix()
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
    dsc = mikecm(fname = dsp,nofile=True,pmass=dsm)
    dpos = np.sqrt((dp[:,0]-dsc[0])**2+(dp[:,1]-dsc[1])**2+(dp[:,2]-dsc[2])**2)/.71*1000/(red+1)
    gpos = np.sqrt((gp[:,0]-dsc[0])**2+(gp[:,1]-dsc[1])**2+(gp[:,2]-dsc[2])**2)/.71*1000/(red+1)
    spos = np.sqrt((sp[:,0]-dsc[0])**2+(sp[:,1]-dsc[1])**2+(sp[:,2]-dsc[2])**2)/.71*1000/(red+1)

    if hnum == '125961':
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
     meta = numFe/numH
     nofloor = meta[meta>4.54877795e-09]
     avgnum = np.mean(nofloor)
     fullmetal = np.log10(nofloor)+4.5
     metal = np.log10(avgnum)+4.5
    if hnum == '007':
      den= den*(1e3)**3 # this was for the old runs when 007 was the only kpc run. need to dbl check this now.
    print 'mstar is %f'%(np.log10(sum(smass)))
    massall, bin_edge = np.histogram(dpos,bins=binz, weights =dmass)
    gmassall, bin_edge = np.histogram( gpos,bins=binz, weights =gmass) 
    smassall, bin_edge = np.histogram( spos,bins=binz, weights =smass)
    dmnum, bin_edge = np.histogram(dpos,bins=binz)
    gnum, bin_edge = np.histogram( gpos,bins=binz) 
    snum, bin_edge = np.histogram( spos,bins=binz)

    tempall, bin_edge, binnum= stats.binned_statistic(gpos, np.log10(temp), statistic='mean', bins=binz)
    tbin = np.logspace(np.log10(min(temp)),np.log10(max(temp)),len(temp+1))
    cumutempall, bin_edge = np.histogram(temp,bins=tbin, weights =np.ones(len(temp)))   
    cumutempall = np.cumsum(cumutempall)

    sortdmp = dpos[np.argsort(dpos)]
    sortgp = gpos[np.argsort(gpos)]
    sortsp = spos[np.argsort(spos)]
    sortdm = dmass[np.argsort(dpos)]
    sortgm = gmass[np.argsort(gpos)]
    sortsm = smass[np.argsort(spos)]
    sortsmcumu = np.cumsum(sortsm)

    gvx = hdf['particles/gas']['vx'].as_matrix()*1000
    gvy = hdf['particles/gas']['vy'].as_matrix()*1000
    gvz = hdf['particles/gas']['vz'].as_matrix()*1000
    gv = np.sqrt(gvx**2+gvy**2+gvz**2)
    vgbin = stats.binned_statistic(gpos,gv,bins = 1000)
    vmaxx = max(np.sqrt(G*np.cumsum(sortdm)/sortdmp))
    vgmax = max(vgbin[0])
    rout = (vgbin[1][vgbin[0].argmax()+1]-vgbin[1][vgbin[0].argmax()])/2+vgbin[1][vgbin[0].argmax()]
    print vmaxx,vgmax,rout

    w = 0
    while sortsmcumu[w] <= 0.5*sum(sortsm):
      rhalf = sortsp[w]
      w+=1

    print 'rhalf is ',rhalf
    w = 0
    while sortdmp[w] <=rhalf:
      w+=1
    dmrhalf=sum(sortdm[:w])
    w = 0
    while sortgp[w] <=rhalf:
      w+=1
    gmrhalf=sum(sortgm[:w])
    w = 0
    while sortsp[w] <=rhalf:
      w+=1
    smrhalf=sum(sortsm[:w])

    w = 0
    while sortdmp[w] <=.5:
      w+=1
    dm500pc=sum(sortdm[:w])
    w = 0
    while sortgp[w] <=.5:
      w+=1
    gm500pc=sum(sortgm[:w])
    w = 0
    while sortsp[w] <=.5:
      w+=1
    sm500pc=sum(sortsm[:w])

    hdf.close()
    return massall, gmassall, smassall, x, rvir,rhalf,rmax, temp, den, sfr,red, tempall,tbin[1:],cumutempall,cNFW,vmax,fullmetal,metal,dmrhalf,gmrhalf,smrhalf,dm500pc,gm500pc,sm500pc,sum(dm),sum(gm),sum(smm),vmaxx,vgmax,rout,dmnum,gnum,snum
  else:
    w = 0
    power = 1
    while np.sqrt((power-0.6)**2) >1e-2:
      w+=1
      N = len(sortdmp[:w+1])
      q = np.sqrt((x-sortdmp[w])**2).argmin()
      rhomean = sum(sortdm[:w+1])/(4/3.*3.14159*sortdmp[w]**3)*0.83120300751
      #print N
      power = np.sqrt(96.45)/8*N/np.log(N)*(rhomean/139.94)**(-1/2.)/0.77905265114710698 #Pulled from eq 20 in Power et al 2003. This term at the end = t_circ(r_vir)/hubble time
      #print sortdmp[w],power
    prad = sortdmp[w]
    w = 0
    while sortdmp[w] <=.5:
      w+=1
    dm500pc=sum(sortdm[:w])
    hdf.close()
    return massall, x,cNFW,vmax,dlogr,dm500pc,dmnum,prad 

def mikecm(fname, nofile=False, centered=False, pmass=[0e0], r_stop=None, **args):
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

#def plot_phase(hnum,hymass,hygmass,hysmass,hyx,temp,den,sfr,date,red):
#	""" UNDER CONSTRUCTION"""
#	def __init__(self,sm):
#		self.fig = plt.figure()
#		self.sub = self.fig.add_subplot(111)
#		self.sm = sm
#  global count,grain
#  hytot = hymass+hygmass+hysmass
#  #print temp, hygmass*1.99e33/(4*3.14159*(hyx*3.09e21)**3)/1.67e-24
#  nbins = 100
#  #temp = np.nan_to_num(temp)
#  H, denedges, tempedges = np.histogram2d(np.log10(den),np.log10(temp),bins=nbins)
#  H = np.rot90(H)
#  H = np.flipud(H)
# 
#  # Mask zeros
#  Hmasked = np.ma.masked_where(H==0,H)
#  #pylab.rcParams['xtick.major.pad']='6'
#  #pylab.rcParams['ytick.major.pad']='6'
#  fig = plt.figure()
#  ax = fig.add_subplot(111)
#  #cou = cm.ScalarMappable(norm=co.Normalize(vmin=0, vmax=5e3))
#  #cou._A = []
#  #cbar = fig.colorbar(cou)
#  phase = ax.pcolormesh(denedges,tempedges,np.log10(Hmasked))
#  #cbar.set_label(r'$\mathrm{Counts}$')
#  cbar = fig.colorbar(phase, ax = ax)
#  cbar.set_label('$\mathrm{Log(Counts)}$')
#  phase.set_clim(vmin=0,vmax=5)
#  #cbar.ScalarMappable()
#  #ax.loglog(hygmass*1.99e33/1.67e-24/(4*3.14159/3*(hyx*3.09e21)**3),temp, '%s'%clr[count],linewidth=4,label = '11_13(ff) MFM')#(4*3.14159/3*hyx**3)
#  #ax.spines['bottom'].set_linewidth(4)
#  #ax.spines['top'].set_linewidth(4)
#  #ax.spines['right'].set_linewidth(4)
#  #ax.spines['left'].set_linewidth(4)
#  #ax.tick_params('both',length=5,width=2,which='minor')
#  #ax.tick_params('both',length=10,width=2,which='major')
#  #ax.xaxis.set_tick_params(labelsize=20)
#  #ax.yaxis.set_tick_params(labelsize=20)
#  #plt.xlabel('Radius (kpc)',fontsize=20, labelpad=-10)
#  #ax.xaxis.set_label_coords(.48,-.07)
#  plt.xlabel(r'$\mathrm{log(\rho)}$ $\mathrm{(n_H/cm^3)}$',fontsize=15, labelpad=-2)
#  plt.ylabel(r'$\mathrm{log(Temperature)}$ $\mathrm{(K)}$',fontsize=15, labelpad=-2)
#  plt.xlim(-6,2)
#  plt.ylim(1,6)
#  #plt.ylim(1e-1,2e0)
#  plt.legend(loc='best',prop={'size':10})
#  ax.text(0.15, 0.3,'$\mathrm{z}$ = %.3f'%red, ha='center', va='center', transform=ax.transAxes, fontsize = 15)
#  ax.text(0.25, 0.2,'$\mathrm{M_{baryon}}$ = %.3e'%sum(hysmass+hygmass), ha='center', va='center', transform=ax.transAxes, fontsize = 15)
#  ax.text(0.225, 0.1,'$\mathrm{SFR}$ = %.3e'%sum(sfr), ha='center', va='center', transform=ax.transAxes, fontsize = 15)
#  #fname = 'm_enc_ratio_halo%s_%dbin_12_16_z0.pdf'%(hnum,grain)
#  fname = 'phase_halo%s_%dbin_%03d.png'%(hnum,grain,i-11)
#  print i
#  plt.savefig(fname,transparent=True)
#  #plt.show()
#  plt.close()
#  return True

class metal_dist(object):
	""" Metallicity Distributions (Cumulative Metalicity Distribution under construction)"""
	def __init__(self,sm):
		self.fig = plt.figure()
		self.sub = self.fig.add_subplot(111)
		self.sm = sm
		self.sub.set_xlabel('$\mathrm{[Fe/H]}$')
		self.sub.set_ylabel('# $\mathrm{of}$ $\mathrm{Stars}$')
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
	def save(self,hnum,date,grain):
		self.fig.savefig('metalhist_halo%s_%dbin_%s_z0.pdf'%(hnum,grain,date),Transparent=True)
		self.fig.show()

class cumutempfunc(object):
	""" Cumulative temperature function"""
	def __init__(self,sm):
		self.fig = plt.figure()
		self.sub = self.fig.add_subplot(111)
		self.sm = sm
		self.sub.set_xlabel(r'$\mathrm{Temperature}$ $\mathrm{(K)}$')
		self.sub.set_ylabel(r'$\mathrm{N (<T)}$')
		cb = self.fig.colorbar(self.sm)
		cb.set_label(r'$\mathrm{log}$ $M_\star$')

	def add_point(self,smass,cumutempx,cumutempy):
		self.sub.loglog(cumutempx,cumutempy, color=sm.to_rgba(np.log10(sum(smass))),linewidth=4)

	def save(self,date):
		self.fig.savefig('cumutempfunc_%s.pdf'%date,transparent=True)
		self.fig.show()
		plt.close()

class avgtempvr(object):
	""" Average Temperature vs r"""
	def __init__(self,sm):
		self.fig = plt.figure()
		self.sub = self.fig.add_subplot(111)
		self.sm = sm
		self.sub.set_xlabel(r'$\mathrm{Radius}$ $\mathrm{(kpc)}$')
		self.sub.set_ylabel(r'$\mathrm{Mean}$ $\mathrm{Temperature}$ $\mathrm{(K)}$')
		cb = self.fig.colorbar(self.sm)
		cb.set_label(r'$\mathrm{log}$ $M_\star$')

	def add_point(self,x,smass,temp):
		self.sub.loglog(x,10**temp, color=self.sm.to_rgba(np.log10(sum(smass))),linewidth=4)
	def save(self,date):
		self.fig.savefig('avgtemp_vs_r_%s.pdf'%date,transparent=True)
		self.fig.show()
		plt.close()

class rhoratiovr(object):
	""" Rho_hydro/rho_dmo vs r/rhalf """
	def __init__(self,sm):
		self.sm = sm
		self.fig = plt.figure()
		self.sub = self.fig.add_subplot(111)
		self.sub.axhline(1,linestyle='--',color='k',linewidth=1)
		self.sub.set_xlabel(r'$\mathrm{r}/\mathrm{r}_{1/2}$',fontsize=20, labelpad=-10)#'r/R$_{vir}$'
		self.sub.set_ylabel(r'$\rho_{\mathrm{hydro}}/\rho_{\mathrm{dmo}}$',fontsize=20, labelpad=-5)
		self.sub.set_xlim(.1,10)#(.06,60)
		self.sub.set_ylim(1e-1,2e0)
		self.sub.xaxis.set_label_coords(.48,-.06)
		self.sub.set_xticklabels([0.1,0.1,1,10])
		cb  = self.fig.colorbar(self.sm)
		cb.set_label(r'$\mathrm{log}$ $M_\star$',fontsize=16)

	def add_line(self,x,dmox,totmass,smass,dmomass,rhalf):
		#print 'totmass is',totmass
		#print 'dmomass is',dmomass
		#print 'x is',x
		#print 'dmox is',dmox
		#print 'ratio is', (totmass)/(4*3.14159*x**3)/(dmomass*0.83120300751/(4*3.14159*dmox**3))
		self.sub.loglog(x/rhalf,(totmass)/(4*3.14159*(x/rhalf)**3)/(dmomass*0.83120300751/(4*3.14159*(dmox/rhalf)**3)), color=self.sm.to_rgba(np.log10(sum(smass))),linewidth=4) 

	def save(self,date,grain):
		self.fig.savefig('radden_ratio_%dbin_%s_z0.pdf'%(grain,date),transparent=True)
		self.fig.show()
		#plt.close()

class trip_radpro(object):
	""" 1 x 3 subplot figure of radial density profiles. Meant to display the transition of cuspy to cored. """
	def __init__(self,sm):
		self.fig = plt.figure()
		self.sub = self.fig.add_subplot(111)
		self.sm = sm
		#Initialize triple horizontal radial density profiles
		self.fig, (self.trip1, self.trip2, self.trip3) = plt.subplots(1,3,sharex=True,sharey=True)
		self.trip2.set_xlabel('$\mathrm{Radius}$ $\mathrm{(kpc)}$',fontsize=10)
		self.trip1.set_ylabel(r'$\rho$ $(\mathrm{M_\odot/kpc^3})$',fontsize=10)

	def add_line(self,x,dmox,totmass,dmomass,dlogr,rhalf,panel,smass):
		if panel == 0:
			self.trip1.loglog(x,totmass/(4*3.14159*x**3)/dlogr, 'k',linewidth=2,label = 'HYDRO')
			self.trip1.loglog(dmox,dmomass/(4*3.14159*dmox**3)*0.83120300751/dlogr, color = 'darkgrey',linestyle='--',linewidth=2, label = 'DMO')
			self.trip1.axvline(x=rhalf,color='k',linestyle=':',linewidth=1)
			self.trip1.text(0.6,2e8,r'M$_\star=$%.2e M$_\odot$'%sum(smass),fontsize=8)
		elif panel == 1:
			self.trip2.loglog(x,totmass/(4*3.14159*x**3)/dlogr, 'k',linewidth=2,label = 'HYDRO')
			self.trip2.loglog(dmox,dmomass/(4*3.14159*dmox**3)*0.83120300751/dlogr, color = 'darkgrey',linestyle='--',linewidth=2, label = 'DMO')
			self.trip2.axvline(x=rhalf,color='k',linestyle=':',linewidth=1)
			self.trip2.text(0.6,2e8,r'M$_\star=$%.2e M$_\odot$'%sum(smass),fontsize=8)
		elif panel == 2:
			self.trip3.loglog(x,totmass/(4*3.14159*x**3)/dlogr, 'k',linewidth=2,label = 'HYDRO')
			self.trip3.loglog(dmox,dmomass/(4*3.14159*dmox**3)*0.83120300751/dlogr, color = 'darkgrey',linestyle='--',linewidth=2, label = 'DMO')
			self.trip3.axvline(x=rhalf,color='k',linestyle=':',linewidth=1)
			self.trip3.text(0.6,2e8,r'M$_\star=$%.2e M$_\odot$'%sum(smass),fontsize=8)
			#trip2.set_xlim(.1,10)#(.06,60)
			#trip2.set_ylim(1e5,1e9)
			self.trip1.set_xlim(.1,20)
			self.trip1.set_ylim(1e5,1e9)
			self.trip1.set_xticklabels([0.1,0.1,1,10])
			self.fig.subplots_adjust(wspace=0.05)
			plt.setp(self.trip1, aspect=0.5, adjustable='box-forced')
			plt.setp(self.trip2, aspect=0.5, adjustable='box-forced')
			plt.setp(self.trip3, aspect=0.5, adjustable='box-forced')
	def save(self,date):
		#labels = ['%.5f' % float(t.get_text()) for t in ax.get_xticklabels()]
		self.fig.savefig('trip_rhoden_848_948_32503_%s.pdf'%date,Transparent=True)
		self.fig.show()
		plt.close()


class trip_correl(object):
	""" 3 x 1 Subplot figure containing scatter plots of:
	    Rhalf vs Mstar
	    Mdyn/Mstar vs Mstar
	    Metallicity vs Mstar """
	def __init__(self,sm):
		self.fig = plt.figure()
		self.sub = self.fig.add_subplot(111)
		self.sm = sm
		#Initialize triple vertical M* correlation plots with one colorbar
		self.fig, (self.trip1, self.trip2) = plt.subplots(2,sharex=True)
		#self.fig, (self.trip1) = plt.subplots(1)
		self.trip2.set_xlabel('$M_\star(<r_{1/2})$',fontsize=16)
		self.trip1.set_ylabel(r'$r_\mathrm{1/2}$ $(\mathrm{kpc})$',fontsize=16)
		self.trip2.set_ylabel(r'$M/L(<r_{1/2})$',fontsize=16)
		self.trip2.yaxis.set_label_coords(-.175,.5)
		#self.trip3.set_ylabel(r'$\mathrm{[Fe/H]}$',fontsize=12)
		#self.trip3.yaxis.set_label_coords(-.2,.5)
		#self.fig.sca(self.trip3)
		plt.yticks(fontsize=16)
		plt.xticks(fontsize=16)
		self.trip1_yliml = -.1#.15
		self.trip1_ylimh = 2.2#1.3
		self.trip1.set_ylim(self.trip1_yliml,self.trip1_ylimh)
		self.trip2_xliml = np.log10(7e4)
		self.trip2_xlimh = np.log10(2e8)
		self.trip2.set_xlim(self.trip2_xliml,self.trip2_xlimh)#np.log10(1.7e5),np.log10(2e7))
		#trip5.set_yscale('log')
		self.trip2_yliml = .001#1
		self.trip2_ylimh = 3.1#2.1
		self.trip2.set_ylim(self.trip2_yliml,self.trip2_ylimh)
		#self.trip3_yliml = -2.7#-2.5
		#self.trip3_ylimh = -1.3#-1.05 hi only#-1.8 median
		#self.trip3.set_ylim(self.trip3_yliml,self.trip3_ylimh)
		self.fig.sca(self.trip2)
		plt.tick_params(which='minor',labelsize=16,pad=2)
		ticks = np.linspace(0,3,4)
		plt.yticks(ticks,['1','10','100','1000'],fontsize=16)#[r'$10^1$',r'$10^2$'],fontsize=10) #[r'$10^6$',r'$10^7$',r'$10^8$',r'$10^9$'],fontsize=16)#
		minor_ticks=[]
		for j in range(2,10):
		  minor_ticks.append(0+np.log10(j))
		for j in range(2,10):
		  minor_ticks.append(1+np.log10(j))
		for j in range(2,10):
		  minor_ticks.append(2+np.log10(j))
		self.trip2.yaxis.set_minor_locator(FixedLocator(minor_ticks))
		minor_labels = ['','','','','','','','','','','','','','','','','','','','','','','','','','','',]#['','30 ','','','60 ','','','']
		self.trip2.yaxis.set_minor_formatter(FixedFormatter(minor_labels))
		ticks = np.linspace(5,8,4)
		plt.xticks(ticks,['$10^5$','$10^6$','$10^7$','$10^8$'],fontsize=16)#[r'$10^1$',r'$10^2$'],fontsize=10)
		minor_ticks=[]
		for j in range(2,10):
		  minor_ticks.append(4+np.log10(j))
		for j in range(2,10):
		  minor_ticks.append(5+np.log10(j))
		for j in range(2,10):
		  minor_ticks.append(6+np.log10(j))
		for j in range(2,10):
		  minor_ticks.append(7+np.log10(j))
		self.trip1.xaxis.set_minor_locator(FixedLocator(minor_ticks))
		minor_labels = ['','','','','','','','','','','','','','','','','','','','','','','','','','','',]#['','30 ','','','60 ','','','']
		self.trip1.xaxis.set_minor_formatter(FixedFormatter(minor_labels))
		#self.fig.subplots_adjust(right=0.8)
		#cbar_ax = self.fig.add_axes([0.65, 0.15, 0.025, 0.7])
		#cb = self.fig.colorbar(self.sm,cax=cbar_ax)
		#cb.set_label(r'$\mathrm{\rho}_{\mathrm{dm(hydro)}}/\mathrm{\rho}_{\mathrm{dmo(corrected)}}\mathrm{(500pc)}$',fontsize=16)

	def add_point(self,hnum,x,dmox,rhalf,metal,totrhalf,smrhalf,smass,res,hismass=0,hirhalf=0,totmass=0,dmomass=0):#tot500=0,dmo500=0):
		v = 0
		while x[v] <=rhalf:
			v+=1
		if res == '':
			mkr = '^'
			self.trip1.plot([np.log10(sum(smass)),np.log10(sum(hismass))],[rhalf,hirhalf],'k')
		else:
			mkr = 'o'
		q = 0
		while x[q] <=.5:
			q+=1
		w = 0
		while dmox[w] <=.5:
			w+=1
		if x[v]-.5>.5-x[v-1]:
		 v-=1
		if dmox[w]-.5>.5-dmox[w-1]:
		 w-=1
		self.trip1.scatter(np.log10(smrhalf),rhalf,marker = 's',s=60,color='b')#self.sm.to_rgba((totmass[q]/(4*3.14159*(x[q])**3))/(dmomass[w]*0.83120300751/(4*3.14159*(dmox[w])**3))),label='%s'%hnum)
		self.trip2.scatter(np.log10(smrhalf),np.log10(totrhalf/smrhalf),marker = 's',s=60,color='b')#self.sm.to_rgba((totmass[q]/(4*3.14159*(x[q])**3))/(dmomass[w]*0.83120300751/(4*3.14159*(dmox[w])**3))),label='%s'%hnum)
		print hnum
		if hnum == '848':
			isodwarf = np.genfromtxt('isoDwarfs.dat')
			isodwarf = isodwarf[:,1:16]
			satmwdwarf = np.genfromtxt('satMWDwarfs.dat')
			satmwdwarf = satmwdwarf[:,1:16]
			self.trip1.errorbar(isodwarf[:,6],4/3.*isodwarf[:,0]/1000,yerr=[4/3.*isodwarf[:,2]/1000,4/3.*isodwarf[:,1]/1000],xerr=[isodwarf[:,8],isodwarf[:,7]],color='k',label='%s'%hnum, fmt='o',markersize=5,elinewidth=2)
			self.trip1.errorbar(satmwdwarf[:,6],4/3.*satmwdwarf[:,0]/1000,yerr=[4/3.*satmwdwarf[:,2]/1000,4/3.*satmwdwarf[:,1]/1000],xerr=[satmwdwarf[:,8],satmwdwarf[:,7]],color='0.6',label='%s'%hnum, fmt='o',markersize=5,elinewidth=2)
			# vv for mdyn/mstar vv
			#yerp = np.sqrt((10**isodwarf[:,5]/10**isodwarf[:,3])**2+(10**isodwarf[:,8]/10**isodwarf[:,6])**2)*(10**isodwarf[:,3])/((10**isodwarf[:,6])/2)
			#yerm = np.sqrt((10**isodwarf[:,4]/10**isodwarf[:,3])**2+(10**isodwarf[:,7]/10**isodwarf[:,6])**2)*(10**isodwarf[:,3])/((10**isodwarf[:,6])/2)
			self.trip2.errorbar(isodwarf[:,6],np.log10(isodwarf[:,12]),yerr=[np.log10((isodwarf[:,12]+isodwarf[:,14])/isodwarf[:,12]),np.log10(isodwarf[:,12]/(isodwarf[:,12]-isodwarf[:,13]))],xerr=[isodwarf[:,8],isodwarf[:,7]],color='k',label='%s'%hnum, fmt='o',markersize=5,elinewidth=2)
			#yerp = np.sqrt((10**satmwdwarf[:,5]/10**satmwdwarf[:,3])**2+(10**satmwdwarf[:,8]/10**satmwdwarf[:,6])**2)*10**satmwdwarf[:,3]/((10**satmwdwarf[:,6])/2)
			#yerm = np.sqrt((10**satmwdwarf[:,4]/10**satmwdwarf[:,3])**2+(10**satmwdwarf[:,7]/10**satmwdwarf[:,6])**2)*10**satmwdwarf[:,3]/((10**satmwdwarf[:,6])/2)
			self.trip2.errorbar(satmwdwarf[:,6],np.log10(satmwdwarf[:,12]),yerr=[np.log10((satmwdwarf[:,12]+satmwdwarf[:,14])/satmwdwarf[:,12]),np.log10(satmwdwarf[:,12]/(satmwdwarf[:,12]-satmwdwarf[:,13]))],xerr=[satmwdwarf[:,8],satmwdwarf[:,7]],color='0.6',label='%s'%hnum, fmt='o',markersize=5,elinewidth=2)
			## vv for mdyn vv ##
			##self.trip2.errorbar(isodwarf[:,6],isodwarf[:,3],yerr=[isodwarf[:,5],isodwarf[:,4]],xerr=[isodwarf[:,8],isodwarf[:,7]],color='k',label='%s'%hnum, fmt='o',markersize=5,elinewidth=2)
			##self.trip2.errorbar(satmwdwarf[:,6],satmwdwarf[:,3],yerr=[satmwdwarf[:,5],satmwdwarf[:,4]],xerr=[satmwdwarf[:,8],satmwdwarf[:,7]],color='0.6',label='%s'%hnum, fmt='o',markersize=5,elinewidth=2)
		#print 'trip 2 is ',satmwdwarf[:,6],10**(satmwdwarf[:,3])/((10**satmwdwarf[:,6])/2)
		#self.trip3.scatter(np.log10(sum(smass)),metal,marker = '%s'%mkr,s=30,color=self.sm.to_rgba(np.log10(sum(smass))),label='%s'%hnum)
		plt.setp(self.trip1, aspect=((self.trip2_xlimh-self.trip2_xliml)/(self.trip1_ylimh-self.trip1_yliml)/1.25), adjustable='box-forced')#2.1/(self.trip1_ylimh-self.trip1_yliml)/1.25
		#trip5.set_yticks([25,50,75,100])
		#self.trip2.set_yscale('log')
		plt.setp(self.trip2, aspect=((self.trip2_xlimh-self.trip2_xliml)/(self.trip2_ylimh-self.trip2_yliml)/1.25), adjustable='box-forced')
		#plt.setp(self.trip3, aspect=(2.1/(self.trip3_ylimh-self.trip3_yliml)/1.25), adjustable='box-forced')
		#^ aspect=(len(xaxis)/len(yaxis)/'squareness')
		#self.trip2.set_xticks([5.0,5.5,6,6.5,7.,7.5,8.0])
		plt.xticks(fontsize=16)
		#self.fig.subplots_adjust(hspace=0.06)

	def save(self,date):
		self.fig.savefig('dbl_correlation_mstar_%s.pdf'%date,Transparent=True)
		self.fig.show()
		#plt.close()

class rhalfvmstar(object):
	""" Scatter plot for rhalf v M_star """
	def __init__(self,sm):
		self.fig = plt.figure()
		self.sub = self.fig.add_subplot(111)
		self.sm = sm
		self.sub.set_xlabel(r'$\mathrm{r}_{1/2}$ $(kpc)$')#$\mathrm{M}_{\mathrm{baryon}}(<5 \mathrm{kpc})$ $(\mathrm{M}_\odot)$')#"$v_{circ}$ (km/s))
		self.sub.set_ylabel(r'$\mathrm{M}_{\star}$ $(\mathrm{M}_\odot)$')#r_{1/2}$ $(kpc)$')
		#self.sub.legend(loc=2,prop={'size':10})
		self.sub.set_xscale('log')
		self.sub.set_yscale('log')
		self.sub.set_xlim(2e-1,2e0)
		self.sub.set_ylim(6e5,2e7)
		self.sub.xaxis.set_label_coords(.48,-.07)
		cb = fig.colorbar(self.sm)
		cb.set_label(r'$\mathrm{log}$ $M_\star$')
		#self.sub.set_xticklabels([])
		#self.sub.set_yticklabels([])

	def add_point(self,hnum,rhalf,smass):
		self.sub.scatter(rhalf,sum(smass),s=90,color=self.sm.to_rgba(np.log10(sum(smass))),label='%s'%hnum)

	def save(self,date):
		self.fig.savefig('mstar_vs_rhalf_%s.pdf'%date,transparent=True)
		self.fig.show()
		plt.close()

class vcircvr (object):
	""" Circular velocity plots v r. Includes points where rhalf is."""
	def __init__(self,sm):
		self.fig = plt.figure()
		self.sub = self.fig.add_subplot(111)
		self.sm = sm
		self.sub.set_ylabel(r'$\mathrm{v}_{\mathrm{circ}}$ $\mathrm{(km/s)}$',fontsize = 20)
		self.sub.set_xlabel(r'$\mathrm{Radius}$ $\mathrm{(kpc)}$',fontsize = 20)
		#self.sub.legend(loc=4,prop={'size':10})
		self.sub.set_xlim(1e-1,1.3e0)
		self.sub.set_ylim(5,3e1)
		#cb13 = fig13.colorbar(sm)
		#cb13.set_label(r'$\mathrm{log}$ $M_\star$')

	def add_line(self,hnum,x,rhalf,totmass,smass):
		G = 4.3e-6 #in kpc/M_sun (km/s)^2 
		v = 0
		while x[v] <=rhalf:
			v+=1
		if (x[v]-rhalf)>(rhalf-x[v-1]):
			v-=1
		cmass = np.cumsum(totmass)
		print 'vmax is \n',max(np.sqrt(G*cmass/x))
		self.sub.plot(x,np.sqrt(G*cmass/x),linewidth = 2,color=self.sm.to_rgba(np.log10(sum(smass))),label='%s'%hnum)
		self.sub.scatter(x[v],np.sqrt(G*cmass[v]/x[v]),marker='s', s = 40, color=self.sm.to_rgba(np.log10(sum(smass))))

	def save(self,date):
		self.fig.savefig('vcirc_profile_%s.pdf'%date,transparent=True)
		self.fig.show()
		plt.close()

class coreden (object):
	""" Scatter plot for core density vs rhalf(ver = 0) or M_star(ver = 1)"""
	def __init__(self,sm,ver):
		self.fig = plt.figure()
		self.sub = self.fig.add_subplot(111)
		self.sm = sm
		self.ver = ver
		cb = self.fig.colorbar(self.sm)
		cb.set_label(r'$\mathrm{logM}_star$')
		if ver == 0:
			self.sub.set_ylabel(r'$\rho_{\mathrm{dm}}\mathrm{(<r_{1/2})}$ $(\mathrm{M}_\odot/\mathrm{kpc}^3)$')#"$v_{circ}$ (km/s)")
			self.sub.set_xlabel(r'$\mathrm{r}_{1/2}$ $(kpc)$')
			#self.sub.legend(loc=2,prop={'size':10})
			self.sub.set_xscale('log')
			self.sub.set_yscale('log')
			self.sub.set_xlim(2e-1,2e0)
			self.sub.set_ylim(4e6,1e8)
			self.sub.xaxis.set_label_coords(.48,-.06)
		elif ver == 1:
			self.sub.set_ylabel(r'$\rho_{\mathrm{dm}}\mathrm{(<r_{1/2})}$ $(\mathrm{M}_\odot/\mathrm{kpc}^3)$')#"$v_{circ}$ (km/s)")
			self.sub.set_xlabel(r'$\mathrm{M}_{\star}$ $(\mathrm{M}_\odot)$')
			#self.sub.legend(loc=2,prop={'size':10})
			self.sub.set_xscale('log')
			self.sub.set_yscale('log')
			self.sub.set_xlim(8e4,2e7)
			self.sub.set_ylim(4e6,1e8)
			self.sub.xaxis.set_label_coords(.48,-.06)

	def add_point(self,hnum,x,rhalf,dmass,smass,dlogr):
		v = 0
		while x[v] <=rhalf:
			v+=1
		if self.ver == 0:
			self.sub.scatter(rhalf,(dmass[v]/(4*3.14159*x[v]**3)/dlogr),s=80,color=self.sm.to_rgba(np.log10(sum(smass))),label='%s'%hnum)
		elif self.ver ==1:
			self.sub.scatter(sum(smass),(dmass[v]/(4*3.14159*x[v]**3)/dlogr),s=80,color=self.sm.to_rgba(np.log10(sum(smass))),label='%s'%hnum)

	def save(self,date):
		if self.ver == 0:
			self.fig.savefig('rhodm_vs_rhalf_loglog_%s.pdf'%date,transparent=True)
		if self.ver == 1:
			self.fig.savefig('rhodm_vs_mstar_loglog_%s.pdf'%date,transparent=True)
		self.fig.show()
		plt.close()

class vcircvrhalf(object):
	"""Scatter plot for circular velocity vs rhalf """
	def __init__(self,sm):
		self.fig = plt.figure()
		self.sub = self.fig.add_subplot(111)
		self.sm = sm
		self.sub.set_ylabel(r'$\mathrm{v}_{\mathrm{circ,r}_{1/2}}$ $\mathrm{(km/s)}$')
		self.sub.set_xlabel(r'$\mathrm{r}_{1/2}$ $\mathrm{(kpc)}$')
		#self.sub.legend(loc=2,prop={'size':10})
		self.sub.set_xlim(0,1.2e0)
		self.sub.set_ylim(5,3e1)
		cb = fig.colorbar(self.sm)
		cb.set_label(r'$\mathrm{log}$ $M_\star$')

	def add_point(self,hnum,x,rhalf,totmass,smass):
		global G
		v = 0
		while x[v] <=rhalf:
			v+=1
		cmass = np.cumsum(totmass)
		self.sub.scatter(rhalf,np.sqrt(G*cmass[v]/x[v]), s = 80, color=self.sm.to_rgba(np.log10(sum(smass))), label='%s'%hnum)
	def save(self,date):
		self.fig.savefig('vcircrhalf_vs_rhalf_%s.pdf'%date,transparent=True)
		self.fig.show()
		plt.close()

class rhoratio(object):
	""" Scatter plot of Rho_hydro/Rho_dmo(corrected) vs mstar (ver = 0) vs rhalf (ver = 1) """
	def __init__(self,sm,ver):
		self.fig = plt.figure()
		self.sub = self.fig.add_subplot(111)
		self.sm = sm
		self.ver = ver
		if ver == 0:
			self.sub.set_xlabel(r'$M_{\mathrm{\star}}$ $(M_\odot)$',fontsize=20)
			self.sub.set_xscale('log')
			self.sub.set_xlim(1.7e5,2e7)
		elif ver == 1:
			self.sub.set_xlabel(r'$r_{\mathrm{1/2}}$ $\mathrm{(kpc)}$',fontsize=20)
			self.sub.set_xlim(1e-1,1.3e0)
		self.sub.set_ylabel(r'$\mathrm{\rho}_{\mathrm{dm(hydro)}}/\mathrm{\rho}_{\mathrm{dmo(corrected)}}\mathrm{(500pc)}$',fontsize=20)#"$v_{circ}$ (km/s))

		#self.sub.legend(loc=2,prop={'size':10})
		#self.sub.set_yscale('log')
		self.sub.set_ylim(.5,1.05)
		self.sub.xaxis.set_label_coords(.48,-.06)
		cb = self.fig.colorbar(self.sm)
		cb.set_label(r'$\mathrm{log}$ $M_\star$',fontsize=16)

	def add_point(self,hnum,x,dmox,totmass,dmomass,smass,dlogr,res,rhalf=0):
		if res == '':
			mkr = '^'
		else:
			mkr = 'o'
		#self.sub.scatter(sum(hysmass),(dmomass[dw]/(4*3.14159*dmox[dw]**3)/dlogr)*0.83120300751,marker = mrkr,facecolors='none',s=80,edgecolors=sm.to_rgba(np.log10(sum(hysmass))))#,label='%s dmo'%hnum)
		v = 0
		while x[v] <=.5:
			v+=1
		w = 0
		while dmox[w] <=.5:
			w+=1
		if x[v]-.5>.5-x[v-1]:
		 v-=1
		if dmox[w]-.5>.5-dmox[w-1]:
		 w-=1
		print '500 pc is ',x[v], dmox[w]
		if self.ver == 0:
			self.sub.scatter(sum(smass),(totmass[v]/(4*3.14159*(x[v])**3)/dlogr)/((dmomass[w]*0.83120300751/(4*3.14159*(dmox[w])**3)/dlogr)),marker = '%s'%mkr,s=80,color=self.sm.to_rgba(np.log10(sum(smass))),label='%s'%hnum)#(80*sum(hysmass)/5e6)

		elif self.ver == 1:
			self.sub.scatter(rhalf,(totmass[v]/(4*3.14159*(x[v])**3)/dlogr)/((dmomass[w]*0.83120300751/(4*3.14159*(dmox[w])**3)/dlogr)),marker = '%s'%mkr,s=80,color=self.sm.to_rgba(np.log10(sum(smass))),label='%s'%hnum)#(80*sum(hysmass)/5e6)

	def save(self,date):
		if self.ver == 0:
		  self.fig.savefig('rhohydrorhodmo_500pc_vs_mstar_loglog_%s.pdf'%date,transparent=True)
		elif self.ver == 1:
		  self.fig.savefig('rhohydrorhodmo_500pc_vs_rhalf_loglog_%s.pdf'%date,transparent=True)
		self.fig.show()
		#plt.close()

class cumumstar_ratio(object):
	""" Cumulative sum of Mstar/ Mstar,tot vs r/rvir ??"""
	def __init__(self,sm):
		self.fig = plt.figure()
		self.sub = self.fig.add_subplot(111)
		self.sm = sm
		self.sub.set_ylabel(r'$\mathrm{M}_star(<\mathrm{r})/\mathrm{M}_{*,\mathrm{tot}}$')
		self.sub.set_xlabel(r'$\mathrm{r}/\mathrm{r}_{\mathrm{vir}}$')
		#self.sub.legend(loc=4,prop={'size':10})
		self.sub.set_xlim(1e-2,1e0)
		self.sub.set_ylim(1e-1,1.2e0)
		self.sub.xaxis.set_label_coords(.48,-.06)
		cb = self.fig.colorbar(self.sm)
		cb.set_label(r'$\mathrm{log}$ $M_\star$')

	def add_point(self,hnum,x,rvir,smass):
		self.sub.semilogx(x/rvir,np.cumsum(smass)/sum(smass),linewidth = 2,color=self.sm.to_rgba(np.log10(sum(smass))),label='%s'%hnum)
	def save(self,date):
		self.fig.savefig('cumu_mstar_ratio_rvirscale_%s.pdf'%date,transparent=True)
		self.fig.show()
		plt.close()

class metalvmstar(object):
	""" Scatter plot of Metallicity vs Mstar """
	def __init__(self,sm):
		self.fig = plt.figure()
		self.sub = self.fig.add_subplot(111)
		self.sm = sm
		self.sub.set_xlabel(r'$\mathrm{M}_\star$ $(\mathrm{M}_\odot)$')
		self.sub.set_ylabel(r'$\mathrm{[Fe/H]}$')
		self.sub.set_xlim(2,10)
		self.sub.set_ylim(-3,0)
		#cb16 = fig16.colorbar(sm)
		#cb16.set_label(r'$\mathrm{log}$ $M_\star$')
		self.sub.set_xticklabels([])
		self.sub.set_yticklabels([])

	def add_point(self,hnum,smass,metal):
		self.sub.scatter(np.log10(sum(smass)),metal,linewidth = 2,color=self.sm.to_rgba(np.log10(sum(smass))),label='%s'%hnum)
	def save(self,date):
		self.fig.savefig('metal_vs_mstar_%s.pdf'%date,transparent=True)
		self.fig.show()
		plt.close()

class mdynmdmo(object):
	"""Scatter plot for M_dyn/M_dmo vs rhalf (ver = 0), M_dmo (ver = 1), cnfw (ver = 2) or vmax (ver = 3) """
	def __init__(self,sm,ver,ex):
		self.fig = plt.figure()
		self.sub = self.fig.add_subplot(111)
		self.sm = sm

		if ver == 0:
			self.sub.set_ylabel(r'$\mathrm{M}_{\mathrm{hydro}}/\mathrm{M}_{\mathrm{dmo,corrected}} \mathrm{(<r_{1/2})}$')# ($M_\odot$)')#"$v_{circ}$ (km/s))
			self.sub.set_xlabel(r'$\mathrm{r}_{1/2}$ $\mathrm{(kpc)}$')#'$log(M_star)$ ($M_\odot$)')
			#self.sub.legend(loc=1,prop={'size':10})
			#self.sub.set_xscale('log')
			#self.sub.set_yscale('log')
			self.sub.set_xlim(1e-1,1.2e0)#set_xlim(1e4,1e8)
			self.sub.set_ylim(4.5e-1,1.1e0)
			self.sub.xaxis.set_label_coords(.48,-.06)
			#self.sub.set_ylim(1e7,3e8)
			cb = self.fig.colorbar(self.sm)
			cb.set_label(r'$\mathrm{log}$ $M_\star$')
		if ver == 1:
			self.sub.set_ylabel(r'$\mathrm{M}_{\mathrm{hydro}}/\mathrm{M}_{\mathrm{dmo,corrected}}(<%d \mathrm{kpc})$'%(ex))# ($M_\odot$)')#"$v_{circ}$ (km/s))
			self.sub.set_xlabel(r'$M_{\mathrm{dmo,corrected}}(<%d \mathrm{kpc})/$$(\mathrm{M_\odot}/10\mathrm{M_\odot})$'%ex)
			#self.sub.legend(loc=1,prop={'size':10})
			#self.sub.set_xscale('log')
			#self.sub.set_yscale('log')
			self.sub.set_xlim(8e-1,1.6e0)#set_xlim(1e4,1e8)
			self.sub.set_ylim(4.5e-1,1.1e0)
			self.sub.xaxis.set_label_coords(.48,-.06)
			#self.sub.set_ylim(1e7,3e8)
			cb = self.fig.colorbar(self.sm)
			cb.set_label(r'$\mathrm{log}$ $M_\star$')

		if ver == 2:
			self.sub.set_ylabel(r'$\mathrm{M}_{\mathrm{hydro}}/\mathrm{M}_{\mathrm{dmo,corrected}}(<%d \mathrm{kpc})$'%(ex))# ($M_\odot$)')#"$v_{circ}$ (km/s))
			self.sub.set_xlabel(r'$\mathrm{cNFW}$')#'$log(M_star)$ ($M_\odot$)')
			#self.sub.legend(loc=1,prop={'size':10})
			#self.sub.set_xscale('log')
			#self.sub.set_yscale('log')
			self.sub.set_xlim(6e0,2e1)#set_xlim(1e4,1e8)
			self.sub.set_ylim(4.5e-1,1.1e0)
			self.sub.xaxis.set_label_coords(.48,-.06)
			#self.sub.set_ylim(1e7,3e8)
			cb = self.fig.colorbar(self.sm)
			cb.set_label(r'$\mathrm{log}$ $M_\star$')

		if ver == 3:
			self.sub.set_ylabel(r'$\mathrm{M}_{\mathrm{hydro}}/\mathrm{M}_{\mathrm{dmo,corrected}}(<%d \mathrm{kpc})$'%(ex))# ($M_\odot$)')#"$v_{circ}$ (km/s))
			self.sub.set_xlabel(r'$\mathrm{v}_{\mathrm{max,dmo}}$ $\mathrm{(km/s)}$')#'$log(M_star)$ ($M_\odot$)')
			#self.sub.legend(loc=1,prop={'size':10})
			#self.sub.set_xscale('log')
			#self.sub.set_yscale('log')
			self.sub.set_xlim(30,50)#set_xlim(1e4,1e8)
			self.sub.set_ylim(4.5e-1,1.1e0)
			self.sub.xaxis.set_label_coords(.48,-.06)
			#self.sub.set_ylim(1e7,3e8)
			cb = self.fig.colorbar(self.sm)
			cb.set_label(r'$\mathrm{log}$ $M_\star$')

	def add_point(self,hnum,x,dx,totmass,dmomass,smass,extent,rhalf=0,cNFW=0,vmax=0):
		v = 0
		while x[v] <=extent:
			v+=1
		dv = 0
		while dx[v] <=extent:
			dv+=1
		if ver == 0:
			self.sub.scatter(rhalf,sum(totmass[:v+1])/sum(dmomass[:dv+1]*0.83120300751),s=80,color=self.sm.to_rgba(np.log10(sum(smass))),label='%s'%hnum)
		elif ver == 1:
			self.sub.scatter(sum(dmomass)/1e10,sum(totmass[:v+1])/sum(dmomass[:dv+1]*0.83120300751),linewidth = 2,color=self.sm.to_rgba(np.log10(sum(smass))),label='%s'%hnum)
		elif ver == 2:
			self.sub.scatter(cNFW,sum(totmass[:v+1])/sum(dmomass[:dv+1]*0.83120300751),linewidth = 2,color=self.sm.to_rgba(np.log10(sum(smass))),label='%s'%hnum)
		elif ver == 3:
			self.sub.scatter(vmax,sum(totmass[:v+1])/sum(dmomass[:dv+1]*0.83120300751),linewidth = 2,color=self.sm.to_rgba(np.log10(sum(smass))),label='%s'%hnum)

	def save(self,date,ver):
		if ver == 0:
			self.fig.savefig('mhydro_mdmo_ratio_1kpc_vs_rhalf_loglog_%s.pdf'%date,transparent=True)
		elif ver == 1:
			self.fig.savefig('mhydro_mdmo_ratio_1kpc_vs_mdmo_loglog_%s.pdf'%date,transparent=True)
		elif ver == 2:
			self.fig.savefig('mhydro_mdmo_ratio_1kpc_vs_cNFW_loglog_%s.pdf'%date,transparent=True)
		elif ver == 3:
			self.fig.savefig('mhydro_mdmo_ratio_1kpc_vs_vmaxdmo_loglog_%s.pdf'%date,transparent=True)
		self.fig.show()
		plt.close()

class mdynmstar(object):
	"""Scatter plot for mdyn/mstar vs mstar """
	def __init__(self,sm):
		self.fig = plt.figure()
		self.sub = self.fig.add_subplot(111)
		self.sm = sm
		self.sub.set_ylabel(r'$\mathrm{M}_{\mathrm{dyn}}/\mathrm{M}_{\mathrm{*}}\mathrm{(<r_{1/2})}$')# ($M_\odot$)')#"$v_{circ}$ (km/s))
		self.sub.set_xlabel(r'$\mathrm{M}_{\star}$ $(\mathrm{M}_\odot)$')
		#self.sub.legend(loc=1,prop={'size':10})
		self.sub.set_xscale('log')
		self.sub.set_yscale('log')
		self.sub.set_xlim(5e5,2e7)
		#self.sub.set_ylim(4.5e-1,1.1e0)
		self.sub.xaxis.set_label_coords(.48,-.06)
		cb = self.fig.colorbar(self.sm)
		cb.set_label(r'$\mathrm{log}$ $M_\star$')

	def add_point(self,hnum,x,totmass,smass,extent):
		v = 0
		while x[v] <=extent:
			v+=1
		self.sub.scatter(sum(smass),sum(totmass[:v+1])/sum(smass[:v+1]),linewidth = 2,color=self.sm.to_rgba(np.log10(sum(smass))),label='%s'%hnum)

	def save(self,date):
		self.fig.savefig('mdynmstar_ratio_vs_mstar_loglog_%s.pdf'%date,transparent=True)
		self.fig.show()
		plt.close()

class mhalomstar(object):
	"""Scatter plot for mhalo vs mstar """
	def __init__(self,sm):
		self.fig = plt.figure()
		self.sub = self.fig.add_subplot(111)
		self.sm = sm
		self.sub.set_ylabel(r'$\mathrm{M}_{\mathrm{\star}}$ $(\mathrm{M}_\odot)$')# ($M_\odot$)')#"$v_{circ}$ (km/s))
		self.sub.set_xlabel(r'$\mathrm{M}_{\mathrm{halo}}$ $(\mathrm{M}_\odot)$')
		#self.sub.legend(loc=1,prop={'size':10})
		self.sub.set_xscale('log')
		self.sub.set_yscale('log')
		self.sub.set_xlim(1e1,2e7)
		self.sub.set_ylim(1e7,2e10)
		self.sub.xaxis.set_label_coords(.48,-.06)
		cb = self.fig.colorbar(self.sm)
		cb.set_label(r'$\mathrm{log}$ a')

	def add_point(self,hnum,dmass,smass,red):
		self.sub.scatter(smass,dmass,linewidth = 2,color=self.sm.to_rgba(np.log10(1/(1+red))),label='%s'%hnum)

	def save(self,date):
		self.fig.savefig('mhalomstar_ratio_vs_mstar_loglog_allred_%s.pdf'%date,transparent=True)
		self.fig.show()

class alpha_vs_mstarmhalo(object):
	"""Money plot from Di Cintio et al. 2014"""
	def __init__(self,sm):
		self.fig = plt.figure()
		self.sub = self.fig.add_subplot(111)
		self.sm = sm
		self.sub.set_ylabel(r'$\alpha(0.01<\mathrm{r/R_{vir}}<0.02)$')# ($M_\odot$)')#"$v_{circ}$ (km/s))
		self.sub.set_xlabel(r'$\mathrm{log}_{10}$ $(\mathrm{M}_\star/\mathrm{M_{halo}})$')
		#self.sub.legend(loc=1,prop={'size':10})
		#self.sub.set_xscale('log')
		self.sub.set_xlim(-4,-1)
		self.sub.set_ylim(-2,4e-1)
		self.sub.xaxis.set_label_coords(.48,-.06)
		cb = self.fig.colorbar(self.sm)
		cb.set_label(r'$\mathrm{log}$ $\mathrm{M}_\star$')

	def add_point(self,hnum,x,dmox,totmass,smass,dmomass,rvir,dlogr,cNFW):
		totmass1 = totmass
		x1 = x
		totmass = totmass[(x>0.01*rvir) & (x<0.02*rvir)]
		x = x[(x>0.01*rvir) & (x<0.02*rvir)]
		print 'x is ',x
		alpha,intercept, r_value, p_value, std_err = stats.linregress(np.log10(x),np.log10(totmass/(4*3.14159*(x)**3)/dlogr))
		print 'alpha is ',alpha
		self.sub.scatter(np.log10(sum(smass)/sum(totmass1)),alpha,linewidth = 2,color=self.sm.to_rgba(np.log10(sum(smass))),label='%s'%hnum)
		fig = plt.figure()
		bob = fig.add_subplot(111)
		bob.loglog(x1,totmass1/(4*3.14159*(x1)**3)/dlogr)
		bob.loglog(dmox,dmomass/(4*3.14159*(dmox)**3)*0.83120300751/dlogr)
		#bob.loglog(x1,(x1)**(alpha)*(10)**(intercept))
		G = 4.3e-6 #in kpc/M_sun (km/s)^2 
		print
		denc = 200/3.*cNFW**3/(np.log(1+cNFW)-(cNFW/(1+cNFW)))
		print 'denc is ',denc
		rhocrit = 3*71e-3**2/(8*3.14159*G)#orders of magnitude are weird. Figure this out
		rs = rvir/cNFW
		nfw = denc*rhocrit/((x1/rs)*(1+x1/rs)**2)*0.83120300751
		#bob.loglog(x1,nfw,'k--')
		jose = np.genfromtxt('f11den.txt')
		bob.loglog(jose[:,0],jose[:,1],'m')
		fig.savefig('jose_denpro_08_10_2016.pdf',transparent=True)
		plt.show()

	def save(self,date):
		self.fig.savefig('alpha_vs_mstarmhalo_%s.pdf'%date,transparent=True)
		self.fig.show()


class radpro(object):
	"""Radial density profile """
	def __init__(self,sm):
		self.fig = plt.figure()
		self.sub = self.fig.add_subplot(111)
		self.fig1 = plt.figure()
		self.sub1 = self.fig1.add_subplot(111)
		self.sub.set_xlabel(r'$M_{\mathrm{halo}}$ $(M_\odot)$')
		self.sub.set_ylabel(r'concentration',fontsize=16, labelpad=-5)
		self.sm = sm
		self.sub.spines['bottom'].set_linewidth(4)
		self.sub.spines['top'].set_linewidth(4)
		self.sub.spines['right'].set_linewidth(4)
		self.sub.spines['left'].set_linewidth(4)
		self.sub.tick_params('both',length=5,width=2,which='minor')
		self.sub.tick_params('both',length=10,width=2,which='major')
		self.sub.xaxis.set_tick_params(labelsize=20)
		self.sub.yaxis.set_tick_params(labelsize=20)
		self.sub.xaxis.set_label_coords(.48,-.07)
		self.sub.set_xlabel(r'$Radius$ $(kpc)$',fontsize=20, labelpad=-10)#'r/R$_{vir}$'
		self.sub.set_ylabel(r'$\rho$ $(M_\odot/kpc^3$)',fontsize=20, labelpad=-5)
		self.sub.set_xlim(.1,10)#.1,10)#(.06,60)
		self.sub.set_ylim(1e5,1e9)#1e5,1e9)
		pylab.rcParams['xtick.major.pad']='6'
		pylab.rcParams['ytick.major.pad']='6'
		cb = self.fig.colorbar(self.sm)
		cb.set_label(r'$\mathrm{log}$ $\mathrm{M}_\star$')

  
	def add_line(self,x,rhalf,totmass,dlogr,res,red,hnum,mvir,c,prad,f,subr,num,smass):
		if res == '':
			clr = 'r'
			lbl = 'LOW '
		#	np.savetxt('Halo%s_raddenZ12_hydro.out'%(hnum),np.column_stack((x,totmass/(4*3.14159*x**3)/dlogr,(dmnum+gnum+snum))))
		else:
			clr = 'k'
			lbl = 'HI '
		#	np.savetxt('Halo%s_raddenZ13_hydro.out'%(hnum),np.column_stack((x,totmass/(4*3.14159*x**3)/dlogr,(dmnum+gnum+snum))))
		#!self.sub.loglog(x,totmass/(4*3.14159*x**3)/dlogr, color = '%s'%clr,linewidth=4,label = 'HYDRO')
		if hnum =='1084':
			self.sub.loglog(x,totmass,linewidth=2,label = 'HYDRO',color='k')
		else:
			self.sub.loglog(x,totmass,linewidth=2,label = 'HYDRO',color=self.sm.to_rgba(np.log10(smass)))#, color = '%s'%clr,linewidth = 4
		#self.sub.text(1,1e8,r'%s M$_\star=$%.2e'%(ver,sum(smassall)),fontsize=16)
		#!self.sub.text(1,1e8,r'z=%.2f'%(red),fontsize=16)
		#self.sub.text(1,2e8,r'%s M$_\star=$%.2e'%(ver,sum(hysmass)),fontsize=16)
		#self.sub.axvline(x=rhalf,color='b',linestyle='--',linewidth=1)
		llim = 0.215828318088
		param_bounds=([0.0,0.0],[1,10])
		totden = totmass#!/dlogr/(4*3.14159*(x)**3)
		print 'x is ', x
		print 'totden is ', totden
		#mvir = 1.14296249e+10   
		#c = 1.83482319e+01
		print mvir,c
		#def CNFW(mvir,c): return lambda r,n,rc : corenfw(r,mvir,c,n,rc)
		#nfwcored = CNFW(mvir,c)
		#(fit, cmatrix)= opt.curve_fit(nfwcored,x[x>llim],np.log(totden[x>llim]),bounds = param_bounds,sigma = x[x>llim])
		#print 'n is ', fit[0]
		#print 'core radius is ',fit[1]
		#residuals = np.log10(totmass[x>llim]/np.exp(corenfw(x[x>llim],mvir,c, *fit)))#!/(4*3.14159*x[x>llim]**3)/dlogr-nfwcored(x[x>llim],*fit))**2)
		#RSS = sum(residuals**2)/len(residuals)
		#f.write('Halo %s n = %f Rcore = %f RSS = %f\n'%(hnum,fit[0],fit[1],RSS))#kappa*tsf
		#self.sub.loglog(x,np.exp(corenfw(x,mvir,c, *fit)),'b',linewidth=2, label='nfw_core') #Core Fitting
		#subr.semilogx(x[x>llim],residuals,marker='s')


	def add_dmoline(self,dmox,dmomass,dlogr,res,c,hnum,prad,g,smass,mvir):
		if res == '':
			clr = 'r'
			lbl = 'Z12 '
		elif res == '_14':
			clr = 'b'
			lbl = 'Z14 '
		elif res == '_13' or res == '_11':
			clr = 'k'
			lbl = 'Z13 '

		#!self.sub.plot(dmox,dmomass/(4*3.14159*dmox**3)*0.83120300751/dlogr, color = '%s'%clr,linestyle='--',linewidth=4, label = '%sDMO'%lbl)#'lightgrey'
		if hnum != '1084':
			self.sub.loglog(dmox,dmomass,linewidth=2, label = '%sDMO'%lbl,color=self.sm.to_rgba(np.log10(smass)))
		else:
			self.sub.loglog(dmox,dmomass,linewidth=2, label = '%sDMO'%lbl,color='k')
		#self.sub.fill_between(np.linspace(1e-3,prad),1e3,1e11,facecolor='k',alpha=0.35)
		#if res == '':
		#	self.sub.axvline(x=prad,color='r',linestyle='--',linewidth=1)
		#	np.savetxt('Halo%s_raddenZ12_dmo.out'%(hnum),np.column_stack((dmox,dmomass/(4*3.14159*dmox**3)*0.83120300751/dlogr,dmonum)))
		#elif res == '_14':
		#	self.sub.axvline(x=prad,color='b',linestyle='--',linewidth=1)
		#elif res == '_13' or res == '_11':
		#	self.sub.axvline(x=prad,color='k',linestyle='--',linewidth=1)
		#	np.savetxt('Halo%s_raddenZ13_dmo.out'%(hnum),np.column_stack((dmox,dmomass/(4*3.14159*dmox**3)*0.83120300751/dlogr,dmonum)))

  #===========================================================================#
  #====For fitting a radial density profile with a psuedo-isothermal curve====#
  #===========================================================================#
  #  chimin =10000
  #for p in np.logspace(-1,np.log10(50)):
  #  (fit, cmatrix)= opt.curve_fit(psuedoiso,hyx,hytot/(4*3.14159/3*hyx**3),p0=(hytot[0]/(4*3.14159/3*hyx[0]**3),p))
  #  chisq = sum((hytot/(4*3.14159/3*hyx**3)-psuedoiso(hyx,*fit))**2/(.3*hytot/(4*3.14159/3*hyx**3))**2)
  #  if chisq < chimin:
  #    chimin = chisq
  #    bestfit = fit

  #dof = len(hytot)-len(fit)
  #print fit, chisq, dof
  #ax.loglog(hyx,psuedoiso(hyx,*bestfit), 'g-',linewidth=2, label='psuedo isothermal') #Core Fitting
  #===========================================================================#
  #====For fitting a radial density profile with a NFW profile====#
  #===========================================================================#
		chimin =1e8
		llim = .25
		hlim = 100
		prad = .2
		#jose = np.genfromtxt('f11den.txt')
		param_bounds=([4e9,0.5],[1.4e10,26])#0.5,26)
		#dmonum[dmonum==0]=1
		dmoden = dmomass#!*0.83120300751/dlogr/(4*3.14159*(dmox)**3)
		#def NFW(mvir): return lambda r,fitc : nfw(r,mvir,fitc)
		#nfww = NFW(mvir)
		#fit =[0,0]
		#!print dmox[dmox>prad]
		(fit, cmatrix)= opt.curve_fit(einasto1,dmox[(dmox>prad) & (dmox<10)],np.log(dmoden[(dmox>prad) &(dmox<10)]),bounds = param_bounds)#dmox,dmomass/(4*3.14159*dmox**3),p0=(mvir,c))
		print 'fit is ',fit
		#!self.sub.loglog(dmox,np.exp(nfw(dmox, *fit)),'g-',linewidth=2, label='nfw') #Core Fitting
		self.sub.loglog(dmox,np.exp(einasto1(dmox, *fit)),'g-',linewidth=2, label=r'einasto $\alpha=0.17$') #Core Fitting
		#!residuals = np.log10(dmomass[dmox>prad]/np.exp(nfw(dmox[dmox>prad],*fit)))
		#!RSS = sum(residuals**2)/len(residuals)
		#!g.write('Halo %s mvir = %.2e cNFW = %f RSS = %f\n'%(hnum,fit[0],fit[1],RSS))
		#self.sub.loglog(jose[:,0],jose[:,1],'m')
		#if hnum == '11707':
		#	print smass,c,fit
		#	self.sub1.scatter(smass,c,c='k',label = 'AHF')
		#	self.sub1.scatter(smass,fit[0],c='r',label = 'Fit')
		#	self.sub1.legend()
		#else:
		#	self.sub1.scatter(smass,c,c='k')
		#	self.sub1.scatter(smass,fit[0],c='r')
		
		return fit

	def save(self,date,hnum,grain):
		#self.sub.legend(loc=3,prop={'size':10})
		self.fig.savefig('raddenpros_dmo_%s.pdf'%(date),transparent=True)#radden_halo%s_%dbin_%s.pdf'%(hnum,grain,date),transparent=True)
		self.fig.show()
		self.fig1.savefig('cNFWcomp_vs_mstar.pdf',transparent=True)
		self.fig1.show()		
		#plt.close()

class vmaxvgmax(object):
	""" Max velocity of gas curve vs max velocity of dm curve. Same in form as fig 4a in Papastergis 2016"""
	def __init__(self,sm):
		self.fig = plt.figure()
		self.sub = self.fig.add_subplot(111)
		self.sm = sm
		self.sub.set_xlim(10,90)
		self.sub.set_xscale('log')
		self.sub.set_xlabel(r'$\mathrm{V_{dm,max}}$ $\mathrm{km }$ $\mathrm{s^{-1})}$')
		self.sub.set_ylim(7,90)
		self.sub.set_yscale('log')
		self.sub.set_ylabel(r'$\mathrm{V_{gas,max}}$ $\mathrm{(km }$ $\mathrm{s^{-1})}$')
		self.sub.plot(np.logspace(np.log10(7),np.log10(90),20),np.logspace(np.log10(7),np.log10(90),20),'k--')
		cb = self.fig.colorbar(self.sm)
		cb.set_label(r'$\mathrm{R_{gas,max}}$')
		#plt.tick_params(which='minor',labelsize=10,pad=2)
		#ticks = np.linspace(1,2,2)
		#plt.yticks(ticks,['10','100'],fontsize=10)#[r'$10^1$',r'$10^2$'],fontsize=10)
		#minor_ticks=[]
		#for j in range(2,10):
		#  minor_ticks.append(1+np.log10(j))
		#self.trip2.yaxis.set_minor_locator(FixedLocator(minor_ticks))
		#minor_labels = ['','30 ','','','60 ','','','']
		#self.trip2.yaxis.set_minor_formatter(FixedFormatter(minor_labels))
	def add_point(self,vmax,vgmax,rout):
		self.sub.scatter(vmax,vgmax,color=self.sm.to_rgba(rout))

	def save(self,date):
		self.fig.savefig('vgmax_vs_vmax_%s.pdf'%(date),transparent=True)
		self.fig.show()

class mstarvmax(object):
	""" Stellar mass vs max velocity. Same in form as fig 9 in MBK 2012"""
	def __init__(self):
		self.fig = plt.figure()
		self.sub = self.fig.add_subplot(111)
		self.sub.set_xlim(10,200)
		self.sub.set_xscale('log')
		self.sub.set_xlabel(r'$\mathrm{V_{max}}$ $\mathrm{(km }$ $\mathrm{s^{-1})}$')
		self.sub.set_ylim(7e4,1e11)
		self.sub.set_yscale('log')
		self.sub.set_ylabel(r'$\mathrm{M}_{\mathrm{\star}}$ $(\mathrm{M}_\odot)$')
		abun = np.genfromtxt('vacc_SHAM.txt')
		self.sub.plot(abun[:,0],abun[:,1])
		#plt.tick_params(which='minor',labelsize=10,pad=2)
		#ticks = np.linspace(1,2,2)
		#plt.yticks(ticks,['10','100'],fontsize=10)#[r'$10^1$',r'$10^2$'],fontsize=10)
		#minor_ticks=[]
		#for j in range(2,10):
		#  minor_ticks.append(1+np.log10(j))
		#self.trip2.yaxis.set_minor_locator(FixedLocator(minor_ticks))
		#minor_labels = ['','30 ','','','60 ','','','']
		#self.trip2.yaxis.set_minor_formatter(FixedFormatter(minor_labels))
	def add_point(self,vmax,mstar):
		print vmax,mstar
		self.sub.scatter(vmax,mstar,c='k',s=80)

	def save(self,date):
		self.fig.savefig('mstar_vs_vmax_%s.pdf'%(date),transparent=True)
		self.fig.show()

def main():

  grain = 100
  sp = np.linspace(30,184,155)
  i=184
  bob = 0
  h0 = 71
  om_l = 0.734
  om_m = 0.266
  conv = 3.085677581e+19
  G = 4.3e-6 #in kpc/M_sun (km/s)^2 
  date = time.strftime("%m_%d_%Y")
  clr = ['olive','olive','darkviolet','darkviolet','lime','lime','olive','olive','gray','gray','r','r','b','b','k','k','g','g','y','y','m','m']
  my_cmap=plt.get_cmap('plasma')
  sm = cm.ScalarMappable(cmap=my_cmap,norm=co.Normalize(vmin=5.004, vmax=7.176))#vmin=5.874853, vmax=7.11806)) 11_13 limits
  sm._A = []
  my_cmap1=plt.get_cmap('viridis')
  sm1 = cm.ScalarMappable(cmap=my_cmap1,norm=co.Normalize(vmin=np.log10(0.007936507936507936), vmax=np.log10(1)))#vmin=5.874853, vmax=7.11806)) 11_13 limits
  sm1._A = []
  my_cmap2=plt.get_cmap('cool')
  sm2 = cm.ScalarMappable(cmap=my_cmap2,norm=co.Normalize(vmin=0.6, vmax=4))
  sm2._A = []
  my_cmap3=plt.get_cmap('viridis_r')
  sm3 = cm.ScalarMappable(cmap=my_cmap3,norm=co.Normalize(vmin=0.6, vmax=1.0))
  sm3._A = []
  hnum =['11707']#,'32503','12596','007','848','796','20192','20910','897','948','1016','32257','10q','10v','1084']
  res = ['_13','_13','_13','_11','_13','_13','_13','_13','_13','_13','_13','_13','_11','_11','_13']
  ver = ['5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16']
  dmover = ['11_13','11_13','11_13','11_13','11_13','11_13','5_12_16','5_12_16','11_13','11_13','11_13','11_13','5_12_16','5_12_16','5_12_16']
  #!rhoratio_plot = rhoratiovr(sm)
  #!trip_correl_plot= trip_correl(sm3)
  #!rhoratiovmstar_plot= rhoratio(sm,0)
  #!rhoratiovrhalf_plot= rhoratio(sm,1)
  #mhalomstar_plot = mhalomstar(sm1)
  #vcirc_plot = vcircvr(sm)
  #trip_radpro_plot = trip_radpro(sm)
  #vmaxvgmax_plot = vmaxvgmax(sm2)
  #alpha_plot = alpha_vs_mstarmhalo(sm)
  ##mstarvmax_plot = mstarvmax()
  f = open('core_fits.out', 'w')
  g = open('nfw_fits.out', 'w')
  figr = plt.figure()
  subr = figr.add_subplot(111)
  subr.set_xlabel('Radius (kpc)')
  subr.set_ylabel('Residual')
  subr.set_xlim(.04,60)#.1,10)#(.06,60)
  radpro_plot = radpro(sm)
  qsq = np.array([],dtype = [('name','|S10'),('qsq','>f4')])
  for w in np.arange(len(hnum)):
   #for i in sp:

    try:
     rhalf = rmax = 0
     if hnum[w] == '10q' or hnum[w] == '10v':
       i = 600
     if hnum[w] == '1084':
       i = 184

#     if hnum[w] == '10q' or hnum[w] == '10v':
#       bob = hnum[w]
#       hnum[w] = '32257'
#       res[w] = '_13'
#       i = 184
     pathname = '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/analysis/'%(hnum[w],res[w],ver[w])
     hdf = pd.HDFStore('%sdataframes/halo%s%s_giz%s_snap%03d.h5'%(pathname,hnum[w],res[w],ver[w],i))
     if hnum[w] == '1084':
      smass = 0
     else:
      smass = sum(hdf['particles/star']['mass'].as_matrix())
     if hnum[w] == '848':
       hnum[w] = '2'
     pathname = '/nobackup/afitts/Gadget-2.0.7/production/gizdm%s%s_raw_output/analysis/'%(hnum[w],res[w])
     #!hist = np.genfromtxt(pathname+'halo%smerger_hist.txt'%hnum[w])
     #!dmomass,dmox,dmocNFW,dmovmax,dlogr,dmo500pc,dmonum,prad = radpro_df(pathname,hnum[w],res[w],dmover[w],hist,1,rhalf,rmax,i,grain,dmover[w])
     #!dmo500pc = dmo500pc *0.83120300751
#     if bob == '10q' or bob == '10v':
#       hnum[w] = bob
#       res[w] = '_11'

     dmox = np.genfromtxt('Halo%s_raddenZ13_dmo.out'%(hnum[w])) 
     dmomass = dmox[:,1]
     dmox = dmox[:,0]
     hdf = pd.HDFStore('%sdataframes/halo%s%s_giz%s_snap%03d.h5'%(pathname,hnum[w],res[w],dmover[w],i))
     rvir = np.float(hdf['props']['rvir'])*1000. 
     binz = np.logspace(np.log(.0014),np.log(rvir),48,base=np.e)
     dlogr = np.log(binz[1]/binz[0])
     dmocNFW = np.float(hdf['props']['cNFW'])
     prad = .22  
     mvir = sum(hdf['particles/dm']['mass'].as_matrix())*0.83120300751
     print mvir
     mvir,c = radpro_plot.add_dmoline(dmox,dmomass,dlogr,res[w],dmocNFW,hnum[w],prad,g,smass,mvir)
     if hnum[w] == '2':
       hnum[w] = '848'
     pathname = '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/analysis/'%(hnum[w],res[w],ver[w])
     #!hist = np.genfromtxt(pathname+'halo%smerger_hist.txt'%hnum[w])

     #!dmass,gmass,smass,x,rvir,rhalf,rmax, temp, den, sfr,red,tempall,cumutempx,cumutempy,cNFW,vmax,fullmetal,metal,dmrhalf,gmrhalf,smrhalf,dm500pc,gm500pc,sm500pc,dmm,gmm,smm,vmaxx,vgmax,rout,dmnum,gnum,snum = radpro_df(pathname,hnum[w],res[w],ver[w],hist,0,0,0,i,grain,dmover[w])
     #!totmass = dmass+gmass+smass
     #!hismass = smass
     #!hirhalf = rhalf
     hydro = np.genfromtxt('Halo%s_raddenZ13_hydro.out'%(hnum[w]))
     totmass = hydro[:,1]
     num = hydro[:,2]
     x = hydro[:,0]
     print 1/48.*sum((np.log(totmass[x>.22])-np.log(dmomass[x>.22]))**2)
     qsq = np.append(qsq,np.array([('Halo%s'%(hnum[w]),1/48.*sum((np.log(totmass[x>.22])-np.log(dmomass[x>.22]))**2))],dtype = qsq.dtype))
     hdf = pd.HDFStore('%sdataframes/halo%s%s_giz%s_snap%03d.h5'%(pathname,hnum[w],res[w],ver[w],i))
     rvir = np.float(hdf['props']['rvir'])*1000
     rhalf = np.float(hdf['props']['rhalf'])*1000 
     red = np.float(hdf['props']['redshift'])
     #mvir = sum(hdf['particles/dm']['mass'].as_matrix())+sum(hdf['particles/gas']['mass'].as_matrix()) + sum(hdf['particles/star']['mass'].as_matrix())
     #smass = sum(hdf['particles/star']['mass'].as_matrix())
     binz = np.logspace(np.log(.0014),np.log(rvir),48,base=np.e)
     dlogr = np.log(binz[1]/binz[0])
     c = np.float(hdf['props']['cNFW'])
     #radpro_plot.add_line(x,rhalf,totmass,dlogr,res[w],red,hnum[w],mvir,c,prad,f,subr,num,smass)
     ##if hnum == '948'
     ## radpro_plot = radpro(sm)
     ## radpro_plot.add_dmoline(dmox,dmomass,dlogr,res[w])
     ## res[w] = '_12'
     #pathname = '/nobackup/afitts/Gadget-2.0.7/production/gizdm%s%s_raw_output/analysis/'%(hnum[w],res[w])
     #hist = np.genfromtxt(pathname+'halo%smerger_hist.txt'%hnum[w])
     #dmomass,dmox,dmocNFW,dmovmax,dmo500pc = radpro_df(pathname,hnum[w],res[w],dmover[w],hist,1,rhalf,rmax,i,grain, dmover[w])

     ## res[w] = '_14'
     ## pathname = '/nobackup/afitts/Gadget-2.0.7/production/gizdm%s%s_raw_output/analysis/'%(hnum[w],res[w])
     ## hist = np.genfromtxt(pathname+'halo%smerger_hist.txt'%hnum[w])
     ## dmomass,dmox,dmocNFW,dmovmax,dmo500pc = radpro_df(pathname,hnum[w],res[w],dmover[w],hist,1,rhalf,rmax,i,grain)  
     ## radpro_plot.add_dmoline(dmox,dmomass,dlogr,res[w])    
     ## radpro_plot.save(date,hnum[w],red,grain,i)
     #!radpro_plot.save(date,hnum[w],grain)
     #mhalomstar_plot.add_point(hnum,dmm,smm,red)
     #vcirc_plot.add_line(hnum[w],x,rhalf,dmass,smass)
     #vmaxvgmax_plot.add_point(vmaxx,vgmax,rout)
     #alpha_plot.add_point(hnum,x,dmox,totmass,smass,dmomass,rvir,dlogr,cNFW)
     #mstarvmax_plot.add_point(vmaxx,smm)
     #!totrhalf = dmrhalf+gmrhalf+smrhalf
     #!tot500pc = dm500pc+gm500pc+sm500pc
     #trip_correl_plot.add_point(hnum[w],x,dmox,rhalf,metal,totrhalf,smrhalf,smass,res[w],0,0,totmass,dmomass)
    except Exception,e:
     print e, 'no gas'


    #if hnum[w] == '2':
    #  trip_radpro_plot.add_line(x,dmox,totmass,dmomass,dlogr,rhalf,0,smass)
    #if hnum[w] == '948':
    #  trip_radpro_plot.add_line(x,dmox,totmass,dmomass,dlogr,rhalf,1,smass)
    #if hnum[w] == '32503':
    #  trip_radpro_plot.add_line(x,dmox,totmass,dmomass,dlogr,rhalf,2,smass)
    #metal_plot = metal_dist(sm)
    #metal_plot.add_dist(fullmetal,grain)
    #metal_plot.save(hnum[w],date,grain)


    #!trip_correl_plot.add_point(hnum[w],x,dmox,rhalf,metal,totrhalf,smrhalf,smass,res[w],0,0,totmass,dmomass)
    #!rhoratio_plot.add_line(x,dmox,totmass,smass,dmomass,rhalf)
    #rhoratio_plot.save(date,grain)
    #if hnum[w] != '10q' and hnum[w] != '10v':
    #!rhoratiovmstar_plot.add_point(hnum[w],x,dmox,totmass,dmomass,smass,dlogr,res[w])
    #!rhoratiovrhalf_plot.add_point(hnum[w],x,dmox,totmass,dmomass,smass,dlogr,res[w],rhalf)
    #print hnum[w]
    print i,hnum[w]

    ##### LOW RES #####
    #res[w] = ''
    #dmover[w] = '5_12_16'
    #try:
     #if hnum[w] == '848':
     #  hnum[w] = '2'
     #  dmover[w] = '5_12_16'
#     if hnum[w] == '10q' or hnum[w] == '10v':
#       bob = hnum[w]
#       hnum[w] = '32257'
#       res[w] = '_13'
#       i = 184
    # pathname = '/nobackup/afitts/Gadget-2.0.7/production/gizdm%s%s_raw_output/analysis/'%(hnum[w],res[w])
    # hist = np.genfromtxt(pathname+'halo%smerger_hist.txt'%hnum[w])
    # dmomass,dmox,dmocNFW,dmovmax,dlogr,dmo500pc,dmonum,prad = radpro_df(pathname,hnum[w],res[w],dmover[w],hist,1,rhalf,rmax,i,grain,dmover[w])
    # dmo500pc = dmo500pc *0.83120300751
#     if bob == '10q' or bob == '10v':
#       hnum[w] = bob
#       res[w] = '_11'    
    # mvir,c = radpro_plot.add_dmoline(dmox,dmomass,dlogr,res[w],dmocNFW,hnum[w],dmonum,prad)
    # if hnum[w] == '2':
    #   hnum[w] = '848'
    # pathname = '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/analysis/'%(hnum[w],res[w],ver[w])
    # hist = np.genfromtxt(pathname+'halo%smerger_hist.txt'%hnum[w])
    # dmass,gmass,smass,x,rvir,rhalf,rmax, temp, den, sfr,red,tempall,cumutempx,cumutempy,cNFW,vmax,fullmetal,metal,dmrhalf,gmrhalf,smrhalf,dm500pc,gm500pc,sm500pc,dmm,gmm,smm,vmaxx,vgmax,rout,dmnum,gnum,snum = radpro_df(pathname,hnum[w],res[w],ver[w],hist,0,0,0,i,grain,dmover[w])
    # totmass = dmass+gmass+smass
    # hismass = smass
    # hirhalf = rhalf
    # radpro_plot.add_line(x,rhalf,totmass,dlogr,res[w],red,hnum[w],dmnum,gnum,snum,mvir,c,prad)
    # radpro_plot.save(date,hnum[w],grain)
    #except Exception,e:
    # print e, 'low res fail'


    ##### HIGHEST RES #####
#@    res[w] = '_14'
 #@   try:
 #@    if hnum[w] == '848':
  #@     hnum[w] = '2'
#     if hnum[w] == '10q' or hnum[w] == '10v':
#       bob = hnum[w]
#       hnum[w] = '32257'
#       res[w] = '_13'
#       i = 184
 #@    pathname = '/nobackup/afitts/Gadget-2.0.7/production/gizdm%s%s_raw_output/analysis/'%(hnum[w],res[w])
     #hist = np.genfromtxt(pathname+'halo%smerger_hist.txt'%hnum[w])
     #dmomass,dmox,dmocNFW,dmovmax,dlogr,dmo500pc,dmonum,prad = radpro_df(pathname,hnum[w],res[w],dmover[w],hist,1,rhalf,rmax,i,grain,dmover[w])
  #@   dmo500pc = dmo500pc *0.83120300751
#     if bob == '10q' or bob == '10v':
#       hnum[w] = bob
#       res[w] = '_11'    
     #mvir,c = radpro_plot.add_dmoline(dmox,dmomass,dlogr,res[w],dmocNFW,hnum[w],dmonum,prad)
 #@    if hnum[w] == '2':
 #@      hnum[w] = '848'
     #pathname = '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/analysis/'%(hnum[w],res[w],ver[w])
     #hist = np.genfromtxt(pathname+'halo%smerger_hist.txt'%hnum[w])
     #dmass,gmass,smass,x,rvir,rhalf,rmax, temp, den, sfr,red,tempall,cumutempx,cumutempy,cNFW,vmax,fullmetal,metal,dmrhalf,gmrhalf,smrhalf,dm500pc,gm500pc,sm500pc,dmm,gmm,smm,vmaxx,vgmax,rout,dmnum,gnum,snum = radpro_df(pathname,hnum[w],res[w],ver[w],hist,0,0,0,i,grain)
     #totmass = dmass+gmass+smass
     #hismass = smass
     #hirhalf = rhalf
     #radpro_plot.add_line(x,rhalf,totmass,dlogr,res[w],red,hnum[w],dmnum,gnum,snum,mvir,c)
 #@   except Exception,e:
 #@    print e, 'highest res fail'
 #! trip_correl_plot.save(date)
  #trip_radpro_plot.save(date)
 #! rhoratio_plot.save(date,grain)
 #! rhoratiovmstar_plot.save(date)
 #! rhoratiovrhalf_plot.save(date)
  #vcirc_plot.save(date)
  #mhalomstar_plot.save(date)
  #vmaxvgmax_plot.save(date)
  #alpha_plot.save(date)
  #mstarvmax_plot.save(date)
  #f.close()
  #figr.savefig('residuals_n_to_1.pdf',transparent=True)
  radpro_plot.save(date,hnum[w],grain)
  np.savetxt('hydro_dmo_fits.out',qsq,fmt = ['%s','%f'])
  return True

main()




