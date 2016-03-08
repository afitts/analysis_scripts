import numpy as np
import sys 
import glob
import pygadgetreader as pg
import scipy.integrate
import pandas as pd

def _a_dot(a, h0, om_m, om_l):
    om_k = 1.0 - om_m - om_l
    return h0 * a * np.sqrt(om_m * (a ** -3) + om_k * (a ** -2) + om_l)


def _a_dot_recip(*args):
    return 1. / _a_dot(*args)

def create_df(pathname,hist,hnum,dm):
  h0 = 71
  om_l = 0.734
  om_m = 0.266
  conv = 3.085677581e+19
  numcores = 16
  numsnaps = 185
  massp = 1.67262178e-24 #in grams
  gamma = 5./3.
  kb = 1.3806e-26 #in km^2*g/s^2
  for i in np.arange(len(hist)):
    print i
    if hist[184-i,1] != 0:
      switch = 0
      for j in np.arange(numcores): #Bring together all the AHF files for the snapshot
        temph = glob.glob(pathname+'analysis/ahf_snap%03d/ahf.snap_%03d.%04d.z*.*.AHF_halos'%(i,i,j))
        temph = str(temph).strip('[]').replace("'","")
        h = np.genfromtxt(temph)
        temppr = glob.glob(pathname+'analysis/ahf_snap%03d/ahf.snap_%03d.%04d.z*.*.AHF_profiles'%(i,i,j))
        temppr = str(temppr).strip('[]').replace("'","")
        p = np.genfromtxt(temppr)
        if switch == 0 and len(h) >0:
          halo = h
          switch = 1
        if switch == 1:
	  try:
            halo = np.vstack((halo,h))
          except:
	    print "nothing there"
      for j in np.arange(len(halo)): #Find our halo from all the halos of the snapshot
       try:
        if halo[j,0] == hist[184-i,1]:
          mvir = halo[j,3]/.71
	  vmax = halo[j,16]
          xcen = halo[j,5]/1000. #Grab relevant data (can always add more later)
          ycen = halo[j,6]/1000.
          zcen = halo[j,7]/1000.
          rvir = halo[j,11]/1000./.71
	  rmax = halo[j,12]/1000./.71
	  rhalf = 0
          sname = "%ssnapdir_%03d/snapshot_%03d"%(pathname,i,i)
	  bsize = pg.readheader(sname, 'boxsize')
	  if bsize > 25: #mfm007 is in kpc while everything else (even gizdm007) is in Mpc. This standardizes all that.
	    kfix = 1000
	  else:
	    kfix = 1
          red = pg.readheader(sname, 'redshift')
	  a = pg.readheader(sname, 'time')
          yr = scipy.integrate.quad(_a_dot_recip, 0,a, (h0, om_m, om_l))[0]*conv/(60*60*24*365*1e9)
	  dpid = pg.readsnap(sname, 'pid', 'dm')
          dmass = pg.readsnap(sname, 'mass', 'dm')*1.e10/.71
          dpos = pg.readsnap(sname, 'pos', 'dm')/kfix
	  dvel = pg.readsnap(sname, 'vel', 'dm')
          ddiff = np.sqrt((dpos[:,0]-xcen)**2+(dpos[:,1]-ycen)**2+(dpos[:,2]-zcen)**2)/.71
	  dpos = dpos[ddiff<rvir]
	  dvel = dvel[ddiff<rvir]
	  dmass = dmass[ddiff<rvir]
          dpid = dpid[ddiff<rvir]
	  ddiff = ddiff[ddiff<rvir]
	  rhalf = 0
	  hdf = pd.HDFStore('%sanalysis/dataframes/halo%s%s_giz%s_snap%03d.h5'%(pathname,hnum,res,ver,i))
	  hdf.put('particles/dm',pd.DataFrame(np.hstack((dpos,np.matrix(ddiff).T,np.matrix(dmass).T)), index = dpid, columns = (['x','y','z','r','mass'])))
	  if dm == 0:
            mstar = halo[j,64]/.71
            zarray = pg.readsnap(sname, 'zarray', 'gas')
            he_mf = zarray[:,1] #Helium mass fraction
            y_he = he_mf/(4*(1-he_mf))
            ne = pg.readsnap(sname, 'ne', 'gas') #Electron Abundance
            mu = (1+4*y_he)/(1+y_he+ne)
            mmw = mu*massp #mean molecular weight
            u = pg.readsnap(sname,'u','gas') #specific internal energy in km^2/s^2
            temp = mmw * (gamma-1.)*u/kb #temperature of gas
	    gpid = pg.readsnap(sname, 'pid', 'gas')
	    spid = pg.readsnap(sname, 'pid', 'star')
            gmass = pg.readsnap(sname, 'mass', 'gas')*1.e10/.71
            smass = pg.readsnap(sname, 'mass', 'star')*1.e10/.71
            gpos = pg.readsnap(sname, 'pos', 'gas')/kfix
            spos = pg.readsnap(sname, 'pos', 'star')/kfix
            gvel = pg.readsnap(sname, 'vel', 'gas')/kfix
            svel = pg.readsnap(sname, 'vel', 'star')/kfix
            sfr = pg.readsnap(sname, 'sfr', 'gas')
	    sft = pg.readsnap(sname, 'age', 'star')
            tmp = np.array([0,0])
	    gdiff = np.sqrt((gpos[:,0]-xcen)**2+(gpos[:,1]-ycen)**2+(gpos[:,2]-zcen)**2)/.71
            sdiff = np.sqrt((spos[:,0]-xcen)**2+(spos[:,1]-ycen)**2+(spos[:,2]-zcen)**2)/.71
	    gpos = gpos[gdiff<rvir]
	    gvel = gvel[gdiff<rvir]
	    gmass = gmass[gdiff<rvir]
	    temp = temp[gdiff<rvir]
	    sfr = sfr[gdiff<rvir]
            gpid = gpid[gdiff<rvir]
	    gdiff = gdiff[gdiff<rvir]
	    spos = spos[sdiff<rvir]
  	    svel = svel[sdiff<rvir]
	    smass = smass[sdiff<rvir]
	    sft = sft[sdiff<rvir]
            spid = spid[sdiff<rvir]
	    sdiff = sdiff[sdiff<rvir]
            sortsdiff = sdiff[np.argsort(sdiff)]
            cumsmass = np.cumsum(smass[np.argsort(sdiff)])
            q = 0
            print len(sortsdiff), mstar, max(cumsmass)
            while (cumsmass[q]<=mstar/2):
              rhalf = sortsdiff[q]
              q = q+1 
	    hdf.put('particles/gas',pd.DataFrame(np.hstack((gpos,np.matrix(gdiff).T,np.matrix(gmass).T,np.matrix(temp).T,np.matrix(sfr).T)), index = gpid, columns = (['x','y','z','r','mass','temp','sfr'])))
	    hdf.put('particles/star',pd.DataFrame(np.hstack((spos,np.matrix(sdiff).T,np.matrix(smass).T,np.matrix(sft).T)), index = spid, columns = (['x','y','z','r','mass','sft'])))
	  props = np.array([[xcen,ycen,zcen,rvir,rhalf,rmax,vmax,red,yr]])
	  hdf.put('props',pd.DataFrame(props,index = [halo[j,0]], columns = (['halox','haloy','haloz','rvir','rhalf','rmax','vmax','redshift','time'])))
	  hdf.close()

       except Exception,f:
          print f
	  props = np.array([[xcen,ycen,zcen,rvir,rhalf,rmax,vmax,red,yr]])
	  print halo[j,0], props
	  hdf = pd.HDFStore('%sanalysis/dataframes/halo%s%s_giz%s_snap%03d.h5'%(pathname,hnum,res,ver,i))
	  hdf.put('props',pd.DataFrame(props,index = [halo[j,0]], columns = (['halox','haloy','haloz','rvir','rhalf','rmax','vmax','redshift','time'])))
	  hdf.close()
          print "MISS"
  return True



#sys.stdout = open('df.log', 'w')
hnum = sys.argv[1]
res = sys.argv[2]
if res == '12' or res== '10':
  res = ''
elif res == '13':
  res = '_13'
elif res == '11':
  res = '_11'
ver = sys.argv[3]
dm = int(sys.argv[4])

if dm == 0:
  pathname = '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/'%(hnum,res,ver)
else:
  pathname = '/nobackup/afitts/Gadget-2.0.7/production/gizdm%s%s_raw_output/'%(hnum,res)
print pathname
hist = np.genfromtxt(pathname+'analysis/halo%smerger_hist.txt'%hnum)

create_df(pathname,hist,hnum,dm)





