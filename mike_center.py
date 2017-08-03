import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as co
from matplotlib import rcParams
import pylab
import pygadgetreader as pg
import scipy.integrate
import time

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

def df(fname):

    hdf = pd.HDFStore(fname)
    dm = hdf['particles/dm']['mass']
    dx = hdf['particles/dm']['x']
    dy = hdf['particles/dm']['y']
    dz = hdf['particles/dm']['z']
    dr = hdf['particles/dm']['r']
    dpos = np.column_stack((dx,dy,dz))
    gm = hdf['particles/gas']['mass']
    gx = hdf['particles/gas']['x']
    gy = hdf['particles/gas']['y']
    gz = hdf['particles/gas']['z']
    gr = hdf['particles/gas']['r']
    gpos = np.column_stack((gx,gy,gz))
    sm = hdf['particles/star']['mass']
    sx = hdf['particles/star']['x']
    sy = hdf['particles/star']['y']
    sz = hdf['particles/star']['z']
    sr = hdf['particles/star']['r']
    spos = np.column_stack((sx,sy,sz))
    allmass = np.append(dm,gm)
    allmass = np.append(allmass,sm)
    allpos = np.vstack((dpos,gpos))
    allpos = np.vstack((allpos,spos))
    hx = hdf['props']['halox'].as_matrix()
    hy = hdf['props']['haloy'].as_matrix()
    hz = hdf['props']['haloz'].as_matrix()
    hcen = np.append(hx,hy)
    hcen = np.append(hcen,hz)
    rvir = np.float(hdf['props']['rvir'])*1000

    return dm,dpos,dr,gm,gpos,gr,sm,spos,sr,allmass,allpos,hcen,rvir

def dfdmo(fname):

    hdf = pd.HDFStore(fname)
    dm = hdf['particles/dm']['mass']
    dx = hdf['particles/dm']['x']
    dy = hdf['particles/dm']['y']
    dz = hdf['particles/dm']['z']
    dr = hdf['particles/dm']['r']
    dpos = np.column_stack((dx,dy,dz))
    hx = hdf['props']['halox'].as_matrix()
    hy = hdf['props']['haloy'].as_matrix()
    hz = hdf['props']['haloz'].as_matrix()
    hcen = np.append(hx,hy)
    hcen = np.append(hcen,hz)
    rvir = np.float(hdf['props']['rvir'])*1000

    return dm,dpos,dr,allmass,allpos,hcen,rvir

def cm(fname, nofile=False, centered=False, pmass=[0e0], r_stop=None, **args):
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

def radpro(pos,mass,cen,label):
  global rvir
  diff = np.sqrt((pos[:,0]-cen[0])**2+(pos[:,1]-cen[1])**2+(pos[:,2]-cen[2])**2)/.71*1000
  binz = np.logspace(np.log10(.06),np.log10(rvir),100)
  massall, bin_edge = np.histogram( diff,bins=binz,weights =mass) 
  x = 10**(np.log10(binz[1:])-np.log10(binz[1]/binz[0])/2)
  dlogr = np.log(binz[1]/binz[0])

  hdf = pd.HDFStore('/nobackup/afitts/Gadget-2.0.7/production/gizdm%s_13_raw_output/analysis/dataframes/halo%s_13_giz11_13_snap184.h5'%(hnum,hnum))
  diffdm = hdf['particles/dm']['r']*1000
  massdm = hdf['particles/dm']['mass']
  massalldm, bin_edge = np.histogram( diffdm,bins=binz,weights =massdm) 
  ax.loglog(x,massall/(4*3.14159*x**3)/dlogr,linewidth=4,label = '%s'%label)
  ax.loglog(x,0.83120300751*massalldm/(4*3.14159*x**3)/dlogr,linestyle='--',color='lightgrey',linewidth=4,label = 'DMO %s'%label)
  #diffdmo = np.sqrt((pos[:,0]-cen[0])**2+(pos[:,1]-cen[1])**2+(pos[:,2]-cen[2])**2)/.71*1000
  #binzdmo = np.logspace(np.log10(.06),np.log10(rvir),100)
  #massalldmo, bin_edge = np.histogram( diff,bins=binz, weights =mass) 
  #xdmo = 10**(np.log10(binz[1:])-np.log10(binz[1]/binz[0])/2)
  #ax.loglog(xdmo,massalldmo*0.83120300751/(4*3.14159*xdmo**3)/dlogr,'r--',linewidth=4,label = '%s DMO'%label)

  return True

def den2d (pos,cen,label,clr):
  diff = np.sqrt((pos[:,0]-cen[0,0])**2+(pos[:,1]-cen[0,1])**2+(pos[:,2]-cen[0,2])**2)/.71*1000
  H, xedges, yedges = np.histogram2d(pos[:,0],pos[:,1],bins=1e3)
  H = np.rot90(H)
  H = np.flipud(H)
 
  # Mask zeros
  Hmasked = np.ma.masked_where(H==0,H)
  #pylab.rcParams['xtick.major.pad']='6'
  #pylab.rcParams['ytick.major.pad']='6'
  ax1.pcolormesh(xedges,yedges,np.log10(Hmasked))
  ax1.set_xlim(cen[0,0]-1e-3,cen[0,0]+1e-3)
  ax1.set_ylim(cen[0,1]-1e-3,cen[0,1]+1e-3)
  for i in np.arange(len(label)):
    ax1.scatter(cen[i,0],cen[i,1],marker='*',s=500,color ='%s'%clr[i],label = '%s'%label[i])

  return True





hnum = sys.argv[1]
dm,dp,dr,gm,gp,gr,sm,sp,sr,am,ap,hc,rvir = df('/nobackup/afitts/Gadget-2.0.7/production/mfm%s_13_giz5_12_16_raw_output/analysis/dataframes/halo%s_13_giz5_12_16_snap184.h5'%(hnum,hnum))

dsp = np.vstack((dp,sp))
dsm = np.append(dm,sm)
dc = cm(fname = dp,nofile=True,pmass=dm)
sc = cm(fname = sp,nofile=True,pmass=sm)
dsc = cm(fname = dsp,nofile=True,pmass=dsm)
ac = cm(fname = ap,nofile=True,pmass=am)
f = open('old_vs_new_centers.txt', 'a')
f.write('%s\t%f\t%f\t%f\t%f\t%f\t%f\n'%(hnum,hc[0],hc[1],hc[2],dsc[0],dsc[1],dsc[2]))
f.close()
#vc = np.array((13.3221399999,11.4755,12.68309)) # rockstar center for 848
print "CENTERS FOR %s"%hnum
print '================='
print 'AHF %f, %f, %f'%(hc[0],hc[1],hc[2])
print 'DM %f, %f, %f'%(dc[0],dc[1],dc[2])
print 'STARS %f, %f, %f'%(sc[0],sc[1],sc[2])
print 'DM + STARS %f, %f, %f'%(dsc[0],dsc[1],dsc[2])
print 'ALL PARTICLES %f, %f, %f'%(ac[0],ac[1],ac[2])
clr = ['b','g','r','c','m']
fig = plt.figure()
ax = fig.add_subplot(111)
#radpro(ap,am,dc,'DM ONLY')
#radpro(ap,am,sc,'STARS ONLY')
radpro(ap,am,dsc,'DM + STARS')
radpro(ap,am,hc,'AHF (ALL PART)')
#radpro(ap,am,vc,'ROCKSTAR')
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
ax.set_xlim(.1,10)#(.06,60)
ax.set_ylim(1e5,1e9)
ax.legend(loc=3,prop={'size':10})
date = time.strftime("%m_%d_%Y")
fname = 'radpro_halo%s_centers_%s.pdf'%(hnum,date)
fig.savefig(fname,transparent=True)


fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
den2d(ap,np.array((dc,sc,dsc,hc)),['DM ONLY','STARS ONLY', 'DM + STARS', 'AHF (ALL PART)'],clr)
ax1.legend()
fname = 'den2d_halo%s_centers_%s.pdf'%(hnum,date)
fig1.savefig(fname,transparent=True)
plt.show()

#cou = cm.ScalarMappable(norm=co.Normalize(vmin=0, vmax=5e3))
#cou._A = []
#cbar = fig.colorbar(cou)

#cbar.set_label(r'$\mathrm{Counts}$')
#cbar = fig.colorbar(phase, ax = ax)
#cbar.set_label('$\mathrm{Log(Counts)}$')
#phase.set_clim(vmin=0,vmax=5)

#diffd = np.sqrt((dc[0]-hc[0])**2+(dc[1]-hc[1])**2+(dc[2]-hc[2])**2)
#diffa = np.sqrt((ac[0]-hc[0])**2+(ac[1]-hc[1])**2+(ac[2]-hc[2])**2)
#diffdm = np.sqrt((dc[0]-vc[0])**2+(dc[1]-vc[1])**2+(dc[2]-vc[2])**2)
#diffam = np.sqrt((ac[0]-vc[0])**2+(ac[1]-vc[1])**2+(ac[2]-vc[2])**2)
#print "The difference using dm only is %f pc"%(diffd*1e6) 
#print "The difference using all particles is %f pc"%(diffa*1e6)
#print "The difference using dmo between mike's and rockstar is %f pc"%(diffdm*1e6) 
#print "The difference using all particles between mike's and rockstar is %f pc"%(diffam*1e6)



