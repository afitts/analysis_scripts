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
import seaborn as sns

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

def table_df(pathname,sname,hist,dmo,rhalf,rmax,hnum,res,ver,snap):
 if hnum != '125961':
  global i, grain, clr, count,numcores,coreden, concen,cgas_rhalf, mhalo, mhalo_rhalf, dfrhalf, dfmstar, dfmstar_rhalf
  massp = 1.67262178e-24 #in grams
  gamma = 5./3.
  kb = 1.3806e-26 #in km^2*g/s^2 
  G = 4.3e-6 #in kpc/M_sun (km/s)^2   
  switch = 0
  print '%sdataframes/halo%s%s_giz%s_snap%03d.h5'%(pathname,hnum,res,ver,snap)
  hdf = pd.HDFStore('%sdataframes/halo%s%s_giz%s_snap%03d.h5'%(pathname,hnum,res,ver,snap))
  red = np.float(hdf['props']['redshift'])
  hx = np.float(hdf['props']['halox'])
  hy = np.float(hdf['props']['haloy'])
  hz = np.float(hdf['props']['haloz'])
  rvir = np.float(hdf['props']['rvir'])*1000/(red+1)
  rhalf = np.float(hdf['props']['rhalf'])*1000/(red+1)
  dmass = hdf['particles/dm']['mass'].as_matrix()
  dmpos =  hdf['particles/dm']['r'].as_matrix()*1000/(red+1)
  vmax = hdf['props']['vmax']
  c = hdf['props']['cNFW']
  i = snap
  sortdmp = dmpos[np.argsort(dmpos)]
  sortdm = dmass[np.argsort(dmpos)]
  w = 0
  while sortdmp[w] <=.5:
    w+=1
  dm500pc=sum(sortdm[:w])
 # if hist[184-i,1] != 0:
 #   switch = 0
 #   for j in np.arange(numcores): #Bring together all the AHF files for the snapshot
 #     temph = glob.glob(pathname+'ahf_snap%03d/ahf.snap_%03d.%04d.z*.*.AHF_halos'%(i,i,j))
 #     temph = str(temph).strip('[]').replace("'","")
 #     h = np.genfromtxt(temph)

#      if switch == 0 and len(h) >0:
#        halo = h
#        switch = 1
#      if switch == 1:
#        try:
#          halo = np.vstack((halo,h))
#        except:
#	  print "nothing there"
    #for j in np.arange(len(halo)): #Find our halo from all the halos of the snapshot
    #  if halo[j,0] == hist[184-i,1]:
#	c = halo[j,42]
  if dmo == 0:
    gmass = hdf['particles/gas']['mass'].as_matrix()
    cgas = hdf['particles/gas']['mass'][(hdf['particles/gas']['temp']<1e4) & (hdf['particles/gas']['r']<(rhalf/1000))]
    cgas = cgas.as_matrix()
    gpos = hdf['particles/gas']['r'].as_matrix()*1000/(red+1)
    smass = hdf['particles/star']['mass'].as_matrix()
    exsmass = hdf['particles/star']['mass'][hdf['particles/star']['r']<(rhalf/1000)]
    exsmass = exsmass.as_matrix()
    spos = hdf['particles/star']['r'].as_matrix()*1000/(red+1)
    temp = hdf['particles/gas']['temp'].as_matrix()
    den = hdf['particles/gas']['rho'].as_matrix()*1e10*1.99e33/1.67e-24/(3.806e24)**3
    sfr = hdf['particles/gas']['sfr'].as_matrix()
    nh = pg.readsnap('/nobackupp8/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/snapdir_%03d/snapshot_%03d.0.hdf5'%(hnum,res,ver,snap,snap),'nh','gas')
    ggpos = pg.readsnap('/nobackupp8/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/snapdir_%03d/snapshot_%03d.0.hdf5'%(hnum,res,ver,snap,snap),'pos','gas')/1000.
    ggmass = pg.readsnap('/nobackupp8/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/snapdir_%03d/snapshot_%03d.0.hdf5'%(hnum,res,ver,snap,snap),'mass','gas')*1e10/.71
    r = np.sqrt((ggpos[:,0]-hx)**2+(ggpos[:,1]-hy)**2+(ggpos[:,2]-hz)**2)/.71
    nhtot = sum(nh[r<rvir/1000]*ggmass[r<rvir/1000])
    #gtot = sum(ggmass[r<rvir/1000])    
    if hnum == '007':
      den= den*(1e3)**3
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
    dmpos = np.sqrt((dp[:,0]-dsc[0])**2+(dp[:,1]-dsc[1])**2+(dp[:,2]-dsc[2])**2)/.71*1000/(red+1)
    gpos = np.sqrt((gp[:,0]-dsc[0])**2+(gp[:,1]-dsc[1])**2+(gp[:,2]-dsc[2])**2)/.71*1000/(red+1)
    spos = np.sqrt((sp[:,0]-dsc[0])**2+(sp[:,1]-dsc[1])**2+(sp[:,2]-dsc[2])**2)/.71*1000/(red+1)

    sortdmp = dmpos[np.argsort(dmpos)]
    sortgp = gpos[np.argsort(gpos)]
    sortsp = spos[np.argsort(spos)]
    sortdm = dmass[np.argsort(dmpos)]
    sortgm = gmass[np.argsort(gpos)]
    sortsm = smass[np.argsort(spos)]
    sorttemp = temp[np.argsort(gpos)]
    sortsmcumu = np.cumsum(sortsm)
    w=0
    while sortsmcumu[w] <= 0.5*sum(sortsm):
      rhalf = sortsp[w]
      w+=1

    vmax = max(np.sqrt(G*np.cumsum(sortdm)/sortdmp))
    vgmax = max(np.sqrt(G*np.cumsum(sortgm)/sortgp))
    #w = 0
    #while np.sqrt(G*np.cumsum(sortgm[w])/sortgp[w]) <=vgmax:
    #  rgout = sortgp[w]
    #  w+=1

    w = 0
    while sortdmp[w] <=rhalf:
      w+=1
    dmrhalf=sum(sortdm[:w])
    w = 0
    while sortgp[w] <=rhalf:
      w+=1
    gmrhalf=sum(sortgm[:w])
    cgasrhalf = sum(sortgm[:w][sorttemp[:w]<1e4])
    w = 0
    while sortgp[w] <=.5:
      w+=1
    gm500pc=sum(sortgm[:w])
    cgas500pc = sum(sortgm[:w][sorttemp[:w]<1e4])
    w = 0
    while sortsp[w] <=rhalf:
      w+=1
    smrhalf=sum(sortsm[:w])
    w = 0
    while sortsp[w] <=.5:
      w+=1
    sm500pc=sum(sortsm[:w])

    coreden = np.append(coreden,(dm500pc+gm500pc+sm500pc))
    concen = np.append(concen, hdf['props']['cNFW'])
    cgas_rhalf = np.append(cgas_rhalf,cgas500pc)
    mhalo = np.append(mhalo,sum(dmass)+sum(smass)+sum(gmass))
    mhalo_rhalf = np.append(mhalo_rhalf,dmrhalf+smrhalf+gmrhalf)
    dfrhalf = np.append(dfrhalf,rhalf)
    dfmstar = np.append(dfmstar,sum(smass))
    dfmstar_rhalf = np.append(dfmstar_rhalf,smrhalf)
    print '\nFor Halo %s:\nmvir is %.2e Msun\nmstar is %.2e Msun\nrhalf is %f kpc\ncNFW is %f\nvmax is %f km/s\nM_dyn/M*(<rhalf) is %f\nCold gas is %.2e\nmgas is %.2e\nrvir is %.2e\nHI tot is %.2e\nGas frac is %.2e\nneutral frac is %.2e\n'%(hnum,sum(dmass),sum(smass),rhalf,c,vmax,(smrhalf+gmrhalf+dmrhalf)/smrhalf,cgasrhalf,sum(gmass),rvir,nhtot,nhtot/(nhtot+sum(smass)),nhtot/sum(gmass))
    #print 'Center (Mpc/h):\n %f, %f, %f'%(hx,hy,hz)
    print max(nh)

    f.write('For Halo %s:\nColdgas/M* (<r_1/2) is %.2e\n'%(hnum,sum(cgas)/sum(exsmass)))#mhalo is %.2e Msun\nmstar is %.2e Msun\nCenter (Mpc/h):\n %f, %f, %f'%(hnum,sum(dmass),sum(smass),hx,hy,hz))
    return True
  elif dmo == 2 :
    gmass = hdf['particles/gas']['mass'].as_matrix()
    cgas = hdf['particles/gas']['mass'][(hdf['particles/gas']['temp']<1e4) & (hdf['particles/gas']['r']<(rhalf/1000))]
    cgas = cgas.as_matrix()
    gpos = hdf['particles/gas']['r'].as_matrix()*1000/(red+1)
    temp = hdf['particles/gas']['temp'].as_matrix()
    den = hdf['particles/gas']['rho'].as_matrix()*1e10*1.99e33/1.67e-24/(3.806e24)**3
    sfr = hdf['particles/gas']['sfr'].as_matrix()
    nh = pg.readsnap('/nobackupp8/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/snapdir_184/snapshot_184.0.hdf5'%(hnum,res,ver),'nh','gas')
    ggpos = pg.readsnap('/nobackupp8/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/snapdir_184/snapshot_184.0.hdf5'%(hnum,res,ver),'pos','gas')/1000.
    ggmass = pg.readsnap('/nobackupp8/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/snapdir_184/snapshot_184.0.hdf5'%(hnum,res,ver),'mass','gas')*1e10/.71
    r = np.sqrt((ggpos[:,0]-hx)**2+(ggpos[:,1]-hy)**2+(ggpos[:,2]-hz)**2)/.71
    nhtot = sum(nh[r<rvir/1000]*ggmass[r<rvir/1000])
    #gtot = sum(ggmass[r<rvir/1000])    
    if hnum == '007':
      den= den*(1e3)**3
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

    dsc = mikecm(fname = dp,nofile=True,pmass=dm)
    dmpos = np.sqrt((dp[:,0]-dsc[0])**2+(dp[:,1]-dsc[1])**2+(dp[:,2]-dsc[2])**2)/.71*1000/(red+1)
    gpos = np.sqrt((gp[:,0]-dsc[0])**2+(gp[:,1]-dsc[1])**2+(gp[:,2]-dsc[2])**2)/.71*1000/(red+1)


    sortdmp = dmpos[np.argsort(dmpos)]
    sortgp = gpos[np.argsort(gpos)]

    sortdm = dmass[np.argsort(dmpos)]
    sortgm = gmass[np.argsort(gpos)]

    sorttemp = temp[np.argsort(gpos)]


    vmax = max(np.sqrt(G*np.cumsum(sortdm)/sortdmp))
    vgmax = max(np.sqrt(G*np.cumsum(sortgm)/sortgp))
    #w = 0
    #while np.sqrt(G*np.cumsum(sortgm[w])/sortgp[w]) <=vgmax:
    #  rgout = sortgp[w]
    #  w+=1

    w = 0
    while sortgp[w] <=.5:
      w+=1
    gm500pc=sum(sortgm[:w])
    cgas500pc = sum(sortgm[:w][sorttemp[:w]<1e4])

    coreden = np.append(coreden,(dm500pc+gm500pc))
    concen = np.append(concen, hdf['props']['cNFW'])
    mhalo = np.append(mhalo,sum(dmass)+sum(gmass))


    print '\nFor Halo %s:\nmvir is %.2e Msun\nrhalf is %f kpc\ncNFW is %f\nvmax is %f km/s\nCold gas in 500pc is %.2e\nmgas is %.2e\nrvir is %.2e\nHI tot is %.2e\nneutral frac is %.2e\n'%(hnum,sum(dmass),rhalf,c,vmax,cgas500pc,sum(gmass),rvir,nhtot,nhtot/sum(gmass))
    #print 'Center (Mpc/h):\n %f, %f, %f'%(hx,hy,hz)
    print max(nh)

    #f.write('For Halo %s:\nColdgas/M* (<r_1/2) is %.2e\n'%(hnum,sum(cgas)/sum(exsmass)))#mhalo is %.2e Msun\nmstar is %.2e Msun\nCenter (Mpc/h):\n %f, %f, %f'%(hnum,sum(dmass),sum(smass),hx,hy,hz))
    return True
  else:
    coreden[len(coreden)-1] = coreden[len(coreden)-1]/(dm500pc*0.83120300751)
    print 'DMO cNFW is %f\nvmax is %f\nM_hydro/M_dmo-0.83120300751 is %.2e\n'%(c,vmax,(mhalo[len(mhalo)-1]/sum(dmass)))
    return True

#hnum = ['32257','32257','11707','11707','32503','32503','125961','125961','12596','12596','007','007','2','2','897','897','1016','1016','796','796','948','948'] #dmo should be placed after hydro run for irhalf variable.
res = ['_13','_13','_13','_13','_13','_13','_11','_11','_13','_13','_13','_13','_13','_13','_13','_13','_13','_13','_13','_13','_13','_13']
#ver = ['11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13']
hnum =['848','2','897','897','1016','1016','796','796','948','948','007','007','11707','11707','12596','12596','32503','32503','20910','20910','20192','20192','32257','32257','1084','1084']#,'10q','10q','10v','10v']#['32257','11707','32503','125961','12596','007','2','897','1016','796','948'] #dmo should be placed after hydro run for irhalf variable.
res = ['','','','','','','','','','','','','','','','','','','','','','','','','','','','','','']#['_13','_13','_13','_13','_13','_13','_13','_13','_13','_13','_11','_11','_13','_13','_13','_13','_13','_13','_13','_13','_13','_13','_13','_13','_13','_13','_11','_13','_11','_13']
ver = ['5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16']#['5_12_16','11_13','5_12_16','11_13','5_12_16','11_13','5_12_16','11_13','5_12_16','11_13','5_12_16','11_13','5_12_16','11_13','5_12_16','11_13','5_12_16','11_13','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','11_13','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16']
snap = [184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,184,600,600,600,600]
dm = [0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,2,1,0,1,0,1]
numcores = 16
coreden = np.array([])
concen = np.array([])
cgas_rhalf = np.array([])
mhalo = np.array([])
mhalo_rhalf = np.array([])
dfrhalf = np.array([])
dfmstar = np.array([])
dfmstar_rhalf = np.array([])
f = open('dwarf_coldgas_mstar_rhalf.txt', 'w')
for w in np.arange(len(hnum)):
  if dm[w] == 0 or dm[w] == 2:
    pathname = '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/analysis/'%(hnum[w],res[w],ver[w])
    sname = "/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/snapdir_%03d/snapshot_%03d"%(hnum[w],res[w],ver[w],snap[w],snap[w])
    hist = np.genfromtxt(pathname+'halo%smerger_hist.txt'%hnum[w])
    table_df(pathname,sname,hist,dm[w],0,0,hnum[w],res[w],ver[w],snap[w])

  else:
    pathname = '/nobackup/afitts/Gadget-2.0.7/production/gizdm%s%s_raw_output/analysis/'%(hnum[w],res[w])
    sname = "/nobackup/afitts/Gadget-2.0.7/production/gizdm%s%s_raw_output/snapdir_%03d/snapshot_%03d"%(hnum[w],res[w],snap[w],snap[w])
    hist = np.genfromtxt(pathname+'halo%smerger_hist.txt'%hnum[w])
    table_df(pathname,sname,hist,dm[w],0,0,hnum[w],res[w],ver[w],snap[w]) 



  #else:
    #pathname = '/nobackup/afitts/Gadget-2.0.7/production/gizdm%s%s_raw_output/analysis/'%(hnum[w],res[w])
    #sname = "/nobackup/afitts/Gadget-2.0.7/production/gizdm%s%s_raw_output/snapdir_%03d/snapshot_%03d"%(hnum[w],res[w],snap[w],snap[w])
    #hist = np.genfromtxt(pathname+'halo%smerger_hist.txt'%hnum[w])
    #dmomass,dmox = radpro_df(pathname,sname,hist,dm[w],rhalf,rmax)
f.close()
sns.set_style("ticks")
data = {'mstar':dfmstar,'mstar_rhalf':dfmstar_rhalf,'rhalf':dfrhalf,'coreden':coreden,'mhalo':mhalo,'mhalo_rhalf':mhalo,'cgas_500pc':cgas_rhalf,'cNFW':concen}
df = pd.DataFrame(data,columns=['mstar','mstar_rhalf','rhalf','coreden','mhalo','mhalo_rhalf','cgas_500pc','cNFW'])
#g = sns.pairplot(df,kind='reg')
g = sns.PairGrid(df,x_vars=['mstar_rhalf','rhalf','coreden','mhalo_rhalf','cgas_500pc','cNFW'],y_vars=['mstar_rhalf','rhalf','coreden','mhalo_rhalf','cgas_500pc','cNFW'], size=5)
g.map(plt.scatter)# s=100)
axes = g.axes
#axes[0,1].set_xlim(.3,1.6)
#axes[2,1].set_xlabel(r'$r_{1/2}$',fontsize=14)
#axes[0,0].set_xlim(5e5,5e7)
#axes[1,0].set_ylim(1e1,1e2)
#axes[2,0].set_ylim(5e5,5e7)
#axes[2,0].set_xlabel(r'$M_\star$',fontsize=14)
axes[0,0].set_xlim(1e5,2e7)
axes[0,0].set_ylim(1e5,2e7)
axes[4,0].set_ylim(1e6,3e7)
axes[0,4].set_xlim(1e6,3e7)
axes[0,0].set_xscale('log')
#axes[0,1].set_xscale('log')
#axes[2,0].set_xscale('log')
axes[0,2].set_xscale('log')
#axes[0,4].set_xscale('log')
#axes[0,5].set_xscale('log')
axes[0,4].set_xscale('log')
#axes[0,7].set_xscale('log')
axes[0,0].set_yscale('log')
axes[2,0].set_yscale('log')
#axes[2,0].set_yscale('log')
#axes[4,0].set_yscale('log')
axes[4,0].set_yscale('log')
#axes[5,0].set_yscale('log')
#axes[6,0].set_yscale('log')
#axes[7,0].set_yscale('log')

#axes[2,0].set_ylabel(r'$M_\star$',fontsize=14)
#axes[0,0].set_ylabel(r'$\rho_{\mathrm{hydro}}/\rho_{DMO}$',fontsize=14)
#axes[1,0].set_ylabel(r'$M_{\mathrm{dyn}}/M_\star$ $(<r_{1/2})$',fontsize=14)
g.savefig('scatter_reg.pdf')
plt.show()

