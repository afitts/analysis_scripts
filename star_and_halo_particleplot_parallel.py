import numpy as np
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import asciitable
import h5py
import pylab
import glob
import pandas as pd
import types
import pygadgetreader as pg
import scipy
import scipy.integrate
import fileinput
import time
import sys
import pickle
from matplotlib.colors import LogNorm
from matplotlib import rcParams
from matplotlib.ticker import FixedLocator
from matplotlib.ticker import FixedFormatter
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

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

def save_obj(obj, name ):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)

def readahf(af):
       data=asciitable.read(af)
       return data

def get_cen_rvir(bob):
    xcen = bob[5]/h
    ycen = bob[6]/h
    zcen = bob[7]/h
    rvir = bob[11]/h
    return xcen,ycen,zcen,rvir

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

def get_starsnhalos(pathname,hnum,res,ver,snum): #Simple function to get 2D scatter plots of star particles with circles to denote rvir of halos at z=0
  global conv,h
  fname = '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/snapdir_%03d/snapshot_%03d.0.hdf5'%(hnum,res,ver,snum,snum)
  dpos = pg.readsnap(fname,'pos','dm')
  dmass = pg.readsnap(fname,'mass','dm')
  xtot = dpos[:,0]/h
  ytot = dpos[:,1]/h
  ztot = dpos[:,2]/h
  hdf = pd.HDFStore('%sanalysis/dataframes/halo%s%s_giz%s_snap%03d.h5'%(pathname,hnum,res,ver,snum))
  xcen = np.float(hdf['props']['halox'])
  ycen = np.float(hdf['props']['haloy'])
  zcen = np.float(hdf['props']['haloz'])
  rvir = np.float(hdf['props']['rvir'])*1000
  for w in np.arange(16):
    temph = glob.glob(pathname+'analysis/ahf_snap%03d/ahf.snap_%03d.%04d.z*.*.AHF_halos'%(snum,snum,w))
    temph = str(temph).strip('[]').replace("'","")
    halo = np.loadtxt(temph,usecols=(0,3,5,6,7,11,64,63)) #Halo ID, Mvir, halox, haloy, haloz, rvir, mstar, # of stars
    for k in np.arange(len(halo)):
      halodiff = np.sqrt((xcen-halo[k,2]/h)**2+(ycen-halo[k,3]/h)**2+(zcen-halo[k,4]/h)**2)
      try:
       if halodiff<100 or halodiff==0:#halo[k,5] > 5 and halo[k,7]>10 and halodiff>20 and halodiff<100 or halodiff==0:# and halodiff<rvir or halodiff==0:
        ans = np.vstack((ans,(halo[k,1],halo[k,6])))
        x = np.append(x,halo[k,2]/h)
        y = np.append(y,halo[k,3]/h)
        z = np.append(z,halo[k,4]/h)
        rad = np.append(rad,halo[k,5]/h)
        bobob += halo[k,6]
	print halo[k,0], halo[k,1],halo[k,6]
      except:
       if halodiff<100 or halodiff==0:#and halodiff<rvir or halodiff ==0:
        ans = np.array([[halo[k,1],halo[k,6]]])
        x = np.array([halo[k,2]])/h
        y = np.array([halo[k,3]])/h
        z = np.array([halo[k,4]])/h
        rad = np.array([halo[k,5]])/h
        bobob = halo[k,6]
	print halo[k,0], halo[k,1],halo[k,6]

  return xtot,ytot,ztot,x,y,z,rad,xcen,ycen,zcen 

def my_circle_scatter(axes, x_array, y_array, radius, **kwargs):
    for x, y, rad in zip(x_array, y_array,radius):
        circle = pylab.Circle((x,y), radius=rad, **kwargs)
        axes.add_patch(circle)
    return True

def plot_star_and_halo(x,y,z,halox,haloy,haloz,halorad,xcen,ycen,zcen,hnum,res,snum): #Plotting function for 2d scatter plots at z=0.
  global date
  plt.figure()
  axes=pylab.axes()
  my_circle_scatter(axes, halox-xcen, haloy-ycen, radius=halorad, fill=False,color = 'b')
  pylab.axis('scaled')
  plt.hist2d(x,y,100,norm = LogNorm(),weights = mass)
  plt.scatter(x[abs(z-zcen)<100]-xcen,y[abs(z-zcen)<100]-ycen,marker = '.',color='r')
  plt.xlim(-100,100)
  plt.ylim(-100,100)
  plt.xlabel('X Position (kpc)')
  plt.ylabel('Y Position (kpc)')
  plt.ticklabel_format(useOffset=False)
  plt.title("Halo %s XY Stellar Particle + Halo Positions at z=0"%(hnum))
  plt.savefig('2den_plots/halo%s%s/halo%s%s_xy_2dden_%d_within100kpc_%s.pdf'%(hnum,res,hnum,res,snum,date),transparent=True)
  plt.show()
  plt.close()

  plt.figure()
  axes=pylab.axes()
  my_circle_scatter(axes, haloy-ycen, haloz-zcen, radius=halorad, fill=False,color = 'b')
  pylab.axis('scaled')
  plt.scatter(y[abs(x-xcen)<100]-ycen,z[abs(x-xcen)<100]-zcen,marker = '.',color='r')
  plt.xlim(-100,100)
  plt.ylim(-100,100)
  plt.xlabel('Y Position (kpc)')
  plt.ylabel('Z Position (kpc)')
  plt.ticklabel_format(useOffset=False)
  plt.title("Halo %s YZ Stellar Particle + Halo Positions at z=0"%hnum)
  plt.savefig('2den_plots/halo%s%s/halo%s%s_yz_2dden_%d_within100kpc_%s.pdf'%(hnum,res,hnum,res,snum,date),transparent=True)
  plt.show()
  plt.close()

  plt.figure()
  axes=pylab.axes()
  my_circle_scatter(axes, halox-xcen, haloz-zcen, radius=halorad, fill=False,color = 'b')
  pylab.axis('scaled')
  plt.scatter(x[abs(y-ycen)<100]-xcen,z[abs(y-ycen)<100]-zcen,marker = '.',color='r')
  plt.xlim(-100,100)
  plt.ylim(-100,100)
  plt.xlabel('X Position (kpc)')
  plt.ylabel('Z Position (kpc)')
  plt.ticklabel_format(useOffset=False)
  plt.title("Halo %s XZ Stellar Particle + Halo Positions at z=0"%hnum)
  plt.savefig('2den_plots/halo%s%s/halo%s%s_xz_2dden_%d_within100kpc_%s.pdf'%(hnum,res,hnum,res,snum,date),transparent=True)
  plt.show()
  plt.close()

  return True

def gather_progen(pathname,snap,hnum): #Gather progenitors of halo from mtree file. Requires Mergertree and my_merger_hist.py to have been already run (for merger_hist.txt)

  mainid = np.genfromtxt(pathname+'analysis/halo%smerger_hist.txt'%hnum)
  if hnum == '10q' or hnum == '10v':
    totsnap = 600
  else:
    totsnap = 184
  try:
    f = np.genfromtxt(pathname+'analysis/merger/halo%s_snap%03d_snap%03d_mtree'%(hnum,snap,snap-1))
    mainid = mainid[totsnap-snap,1]
    maindex = np.where(f[:,0] == mainid) #Find our halo
    maindex = int(maindex[0])
    numprogen = f[maindex,2]
    haloid = f[maindex,0]
    progen = f[maindex+1:(maindex+1)+(numprogen+1),1] #Grab progenitors from table beneath Halo id
    mergerrank = f[maindex+1:(maindex+1)+(numprogen+1),0]**2/(f[maindex+1:(maindex+1)+(numprogen+1),2]*f[maindex,1])
    mergerrank_sort = np.argsort(mergerrank)
    mainmergerid = progen[mergerrank_sort[len(mergerrank_sort)-2]]
    mainprogenid = progen[0]
  except Exception,e:
    print e
    haloid = 'skip'
    mainprogenid = 0
    mainmergerid = 0
    progen = np.zeros(9)

  return haloid,mainprogenid,progen,mainmergerid

def _a_dot(a, h0, om_m, om_l):
    om_k = 1.0 - om_m - om_l
    return h0 * a * np.sqrt(om_m * (a ** -3) + om_k * (a ** -2) + om_l)


def _a_dot_recip(*args):
    return 1. / _a_dot(*args)

def stellar_age(pathname,snum,main_part): #Add stellar creation time
  h0 = 71
  om_l = 0.734
  om_m = 0.266
  conv = 3.085677581e+19
  age = pg.readsnap(pathname+'snapdir_%03d/snapshot_%03d'%(snum,snum),'age','star')
  test = pg.readheader(pathname+'snapdir_%03d/snapshot_%03d'%(snum,snum),'time')
  print age[age>test]
  for i in np.arange(len(age)):
    age[i] = scipy.integrate.quad(_a_dot_recip, 0,age[i], (h0, om_m, om_l))[0]*conv
  test = scipy.integrate.quad(_a_dot_recip, 0,test, (h0, om_m, om_l))[0]*conv
  #age = np.array([age]) #So that we can easily take the transpose of age to concatenate with mass_progen
  pid = pg.readsnap(pathname+'snapdir_%03d/snapshot_%03d'%(snum,snum),'pid','star')
  bob = age[np.in1d(pid, main_part[:,0].astype('int'))]
  main_part[np.argsort(main_part[:,0].astype('int')),11] = bob[np.argsort(pid[np.in1d(pid, main_part[:,0].astype('int'))])]/3600/24/365/1e9
  #print main_part[np.argsort(main_part[:,0].astype('int')),0],main_part[np.argsort(main_part[:,0].astype('int')),11]
  #main_part[np.argsort(main_part[np.in1d(main_part[:,0].astype('int'),bob[:,0].astype('int')),0])= bob[np.argsort(bob[np.in1d(bob[:,0].astype('int'),bob[:,0].astype('int')),0])

  return main_part

def get_projden(pathname,snum,main_part):

  global h
  mass = pg.readsnap(pathname+'snapdir_%03d/snapshot_%03d'%(snum,snum),'mass','star')/h
  pid = pg.readsnap(pathname+'snapdir_%03d/snapshot_%03d'%(snum,snum),'pid','star')
  bob = mass[np.in1d(pid, main_part[:,0].astype('int'))]
  main_part[np.argsort(main_part[np.in1d(main_part[:,0].astype('int'),pid),0]),12] = bob[np.argsort(pid[np.in1d(pid, main_part[:,0].astype('int'))])]*1e10 

  return main_part

def formation_radius(pathname,snum,i,mainprogenid,halostars): #Add formation radius for main progenitor star particles
  
  global h,conv

  temp = glob.glob(pathname+'analysis/ahf_snap%03d/ahf.snap_%03d.%04d.z*.*.AHF_halos'%(snum-1,snum-1,i))
  temp = str(temp).strip('[]').replace("'","")
  ahfdata = readahf(temp)
  for i in np.arange(len(ahfdata)):
    if ahfdata[i][0] == mainprogenid:
      ahfhalo = ahfdata[i]
  xcen,ycen,zcen,rvir = get_cen_rvir(ahfhalo)
  halostars[:,9] = np.sqrt((halostars[:,2].astype('float')/h*conv-xcen)**2+(halostars[:,3].astype('float')/h*conv-ycen)**2+(halostars[:,4].astype('float')/h*conv-zcen)**2) #Formation Radius

  return halostars

def gather_particles(main_progen,pathname,snum,hnum,fnum): #Use progenitors array to gather all their star particles in AHF particle files into the main particle array (now we have each particles' ID, position, velocity and whether they were part of the main progenitor or from a merger)

  global switch
  h0 = 71
  om_l = 0.734
  om_m = 0.266
  conv = 3.085677581e+19
  fix = 1
  mainhaloid, mainprogenid, progen = gather_progen(pathname,snum,hnum)
  my_cols = ['id', 'type', 'x', 'y', 'z', 'vx', 'vy', 'vz']
  for i in np.arange(fnum):
    temp = glob.glob(pathname+'analysis/ahf_snap%03d/ahf.snap_%03d.%04d.z*.*.AHF_particles'%(snum-1,snum-1,i))
    temp = str(temp).strip('[]').replace("'","")
    print temp
    prev_particles = pd.read_csv(temp, names=my_cols, sep='\t')
    prev_particles['main'] = pd.Series(2*np.ones(len(prev_particles)),index=prev_particles.index) # append a new column that denotes whether the particle is in the main halo (1), is from a merger(0) or is undetermined (2). 
    prev_particles['r_form'] = pd.Series(np.zeros(len(prev_particles)),index=prev_particles.index)
    prev_particles['atime'] = pd.Series(np.zeros(len(prev_particles)),index=prev_particles.index)
    prev_particles['ctime'] = pd.Series(np.zeros(len(prev_particles)),index=prev_particles.index)
    prev_particles['mass'] = pd.Series(np.zeros(len(prev_particles)),index=prev_particles.index)

    for j,e in enumerate(prev_particles['id']): # Search through particle file to find instances where the first column contains two numbers (acts as a header for each halo). 
      if type(e)==types.StringType and len(e.split())==2 and (int(e.split()[1]) in progen)==True: # Make sure the halo id is in the progenitor array
        haloid = int(e.split()[1])# e.split()[0] is the number of particles, [1] is the halo ID
	numpart = int(e.split()[0])
	halostars =np.array(prev_particles.ix[j+1:j+numpart][prev_particles.ix[j+1:j+numpart,'type']==4])
	if len(halostars>0):
	  try:
	    if haloid == mainprogenid: #Main progenitor?
	      if switch == 0:
	        halostars[:,8] = 1
		halostars = formation_radius(pathname,snum,i,mainprogenid,halostars)
	        switch = 1
	      else:
	        halostars[np.logical_not(np.in1d(halostars[:,0].astype('int'),main_progen[:,0].astype('int'))),8] = 3 #Mark particles in halostars that arent in main_progen as 3.
	        main_progen[np.logical_not(np.in1d(main_progen[:,0].astype('int'),halostars[:,0].astype('int'))),8] = 3 # Mark particles in main_progen that arent in halostars as 3. This is just in case the main progenitor did not contribute ALL of its particles to the daughter.
	        hstar_und = halostars[:,8]==2 #Undetermined particles remaining in the halostars array (need their flags from main_progen) <---this might be overkill but its good to have regardless
	        main_need = main_progen[:,8]!=3 #Cuts out all the particles that did not participate in merger
	        halostars[np.argsort(halostars[hstar_und,0].astype('int')),8] = main_progen[np.argsort(main_progen[main_need,0].astype('int')),8]
		halostars[np.argsort(halostars[hstar_und,0].astype('int')),9] = main_progen[np.argsort(main_progen[main_need,0].astype('int')),9]
		halostars[np.argsort(halostars[hstar_und,0].astype('int')),10] = main_progen[np.argsort(main_progen[main_need,0].astype('int')),10]
		halostars[np.argsort(halostars[hstar_und,0].astype('int')),11] = main_progen[np.argsort(main_progen[main_need,0].astype('int')),11]
		halostars[np.argsort(halostars[hstar_und,0].astype('int')),12] = main_progen[np.argsort(main_progen[main_need,0].astype('int')),12]
	    else:# Then Merger.
	      halostars[:,8] = 0 #Now the last column becomes a zero to denote that the particle came from a merger MAY NEED TO EDIT THIS. CLEARLY THE COLUMNS ARE BEING PUT IN OUT OF ORDER.
	      atime = pg.readheader(pathname+'snapdir_%03d/snapshot_%03d'%(snum,snum),'time')
	      atime = scipy.integrate.quad(_a_dot_recip, 0,atime, (h0, om_m, om_l))[0]*conv/3600/24/365/1e9
	      halostars[:,10] =  atime# SET EQUAL TO ACCRETION TIME. ALTER C TIME TO 11 AND MASS FOR 12!!!!!
	    progen_part = np.vstack((progen_part, halostars)) #Add the id of particle, the type (only allow stars at this point), position, velocity and main or merger flag. 
	  except Exception,f:
            print f
	    if haloid == mainprogenid: 
	      if switch == 0:
	        halostars[:,8] = 1
		halostars = formation_radius(pathname,snum,i,mainprogenid,halostars)
	        switch = 1
	      else:
	        halostars[np.logical_not(np.in1d(halostars[:,0].astype('int'),main_progen[:,0].astype('int'))),8] = 3 
	        main_progen[np.logical_not(np.in1d(main_progen[:,0].astype('int'),halostars[:,0].astype('int'))),8] = 3
	        hstar_und = halostars[:,8]==2 #Undetermined particles remaining in the halostars array (need their flags from main_progen) <---this might be overkill but its good to have regardless
	        main_need = main_progen[:,8]!=3 #Cuts out all the particles that did not participate in merger
		if len(halostars[hstar_und])>0:
	          halostars[np.argsort(halostars[hstar_und,0].astype('int')),8] = main_progen[np.argsort(main_progen[main_need,0].astype('int')),8]
	    else:
	      halostars[:,8] = 0
	    progen_part = halostars
  try:
    progen_part = progen_part[np.unique(progen_part[:,0].astype('int'),return_index=True)[1]] #Get rid of duplicates (In AHF, host halos calim ownership of sub halos' particles.
  except:
    progen_part = np.zeros((2,12))
  for i in np.arange(fnum):
    temp = glob.glob(pathname+'analysis/ahf_snap%03d/ahf.snap_%03d.%04d.z*.*.AHF_particles'%(snum,snum,i))
    temp = str(temp).strip('[]').replace("'","")
    current_particles = pd.read_csv(temp, names=my_cols, sep='\t')
    current_particles['main'] = pd.Series(2*np.ones(len(current_particles)),index=current_particles.index)
    current_particles['r_form'] = pd.Series(np.zeros(len(current_particles)),index=current_particles.index)
    current_particles['atime'] = pd.Series(np.zeros(len(current_particles)),index=current_particles.index)
    current_particles['ctime'] = pd.Series(np.zeros(len(current_particles)),index=current_particles.index)
    current_particles['mass'] = pd.Series(np.zeros(len(current_particles)),index=current_particles.index)
    for j,e in enumerate(current_particles['id']): 
       if type(e)==types.StringType and len(e.split())==2 and int(e.split()[1]) == mainhaloid:
         haloid = int(e.split()[1])# e.split()[0] is the number of particles, [1] is the halo ID
	 numpart = int(e.split()[0])
         main_part = np.array(current_particles.ix[j+1:j+numpart][current_particles.ix[j+1:j+numpart,'type']==4])
	 temp1 = glob.glob(pathname+'analysis/ahf_snap%03d/ahf.snap_%03d.%04d.z*.*.AHF_halos'%(snum,snum,i))
	 temp1 = str(temp1).strip('[]').replace("'","")
	 ahfdata = readahf(temp1)
	 for i in np.arange(len(ahfdata)):
	  if ahfdata[i][0] == mainhaloid:
	    ahfhalo = ahfdata[i]
	 xcen,ycen,zcen,rvir = get_cen_rvir(ahfhalo)
	 if len(main_part>0):
           main_part[np.logical_not(np.in1d(main_part[:,0].astype('int'),progen_part[:,0].astype('int'))),8] = 1 # If there are particles in main part that aren't in the progen_part, they must have formed in main halo and must be given that distinction
           main_part[np.logical_not(np.in1d(main_part[:,0].astype('int'),progen_part[:,0].astype('int'))),9] = np.sqrt(((main_part[np.logical_not(np.in1d(main_part[:,0].astype('int'),progen_part[:,0].astype('int'))),2]/h*fix-xcen)**2+(main_part[np.logical_not(np.in1d(main_part[:,0].astype('int'),progen_part[:,0].astype('int'))),3]/h*fix-ycen)**2+(main_part[np.logical_not(np.in1d(main_part[:,0].astype('int'),progen_part[:,0].astype('int'))),4]/h*fix-zcen)**2).astype('float')) # If they formed in the main halo, mark down their formation radius.
	   progen_part[np.logical_not(np.in1d(progen_part[:,0].astype('int'),main_part[:,0].astype('int'))),8] = 3 #This flag is set so that these particles are ignored in the next step.  These are stars that either died or didn't take part in the merger.
	   main_und = main_part[main_part[:,8]==2] #Undetermined particles remaining in the main_part array (need their flags from progen_part)
	   progen2main = progen_part[progen_part[:,8]!=3] #particles in progen_part that contribute to main.
	   if len(main_part[0]) >9 and len(progen_part[0])>9:
             main_und[np.argsort(main_und[:,0].astype('int')),9]=progen2main[np.argsort(progen2main[:,0].astype('int')),9]
	     main_part[main_part[:,8]==2,9]=main_und[:,9]
	   if len(main_part[0]) >10 and len(progen_part[0])>12:
             main_und[np.argsort(main_und[:,0].astype('int')),10]=progen2main[np.argsort(progen2main[:,0].astype('int')),10]
	     main_part[main_part[:,8]==2,10]=main_und[:,10]
             main_und[np.argsort(main_und[:,0].astype('int')),11]=progen2main[np.argsort(progen2main[:,0].astype('int')),11]
	     main_part[main_part[:,8]==2,10]=main_und[:,10]
             main_und[np.argsort(main_und[:,0].astype('int')),12]=progen2main[np.argsort(progen2main[:,0].astype('int')),12]
	     main_part[main_part[:,8]==2,10]=main_und[:,10]
           main_und[np.argsort(main_und[:,0].astype('int')),8]=progen2main[np.argsort(progen2main[:,0].astype('int')),8] #Sort particles by id (in both main and progen). Set main_part flags (main or merger) by/
	   main_part[main_part[:,8]==2,8]=main_und[:,8]										 # gathering them from progen (thanks to our previous slices, the sorted arrays should line up now.)

	 try:
	   main_part = stellar_age(pathname,snum,main_part)
	   main_part = get_projden(pathname,snum,main_part)
	 except Exception,f:
            print f
	 stars_in_halo = main_part
	 print snum, len(stars_in_halo[stars_in_halo[:,8]==0]),len(stars_in_halo[stars_in_halo[:,8]==1])
  try:
    b = stars_in_halo[stars_in_halo[:,10]!=0]
    b = b[(stars_in_halo[stars_in_halo[:,10]!=0,10])<(stars_in_halo[stars_in_halo[:,10]!=0,11]),0]
    if len(b)>0:
      print b
      return False
    else:
      return stars_in_halo
  except Exception,f:
    print f
    return np.zeros((2,9))

def plot_starmerge(dwarf,hnum,res,halox,haloy,halorad,xcen,ycen):
  global h, date
  agecut = 3.09 #2.12 (z=3) 2.58 (z=2.5) 3.09(z=2.1) 3.24 (z=2) 5.77 (z=1) 
  figprops = dict(figsize=(8., 8. / 1.618), dpi=128)                                          # Figure properties
  adjustprops = dict(left=0.1, bottom=0.1, right=0.97, top=0.93, wspace=0.001, hspace=0.001)       # Subplot properties

  fig = pylab.figure(**figprops)                                                              # New figure
  fig.subplots_adjust(**adjustprops)                                                          # Tunes the subplot layout

  ax = fig.add_subplot(2, 3, 1)
  bx = fig.add_subplot(2, 3, 2, sharex=ax, sharey=ax)
  cx = fig.add_subplot(2, 3, 3, sharex=ax, sharey=ax)
  dx = fig.add_subplot(2, 3, 4, sharex=ax, sharey=ax)
  ex = fig.add_subplot(2, 3, 5, sharex=ax, sharey=ax)
  fx = fig.add_subplot(2, 3, 6, sharex=ax, sharey=ax)

  merger = dwarf[:,8]==0
  mergery = (dwarf[:,8]==0) & (dwarf[:,11]>agecut)
  mergero = (dwarf[:,8]==0) & (dwarf[:,11]<agecut)
  main = dwarf[:,8]==1
  mainy = (dwarf[:,8]==1) & (dwarf[:,11]>agecut)
  maino = (dwarf[:,8]==1) & (dwarf[:,11]<agecut)
  ax.scatter(dwarf[main,2]/h-xcen,dwarf[main,3]/h-ycen,marker = 'o',s = 10,c='b',edgecolors='none',label = 'merger')
  ax.text(-17,17, 'All in situ', size=10)
  ax.text(-17,-17, '$M_{*,tot}$=%.3g'%sum(dwarf[main,12]),size=8)
  bx.scatter(dwarf[maino,2]/h-xcen,dwarf[maino,3]/h-ycen,marker = 'o',s = 10,c='b',edgecolors='none',label = 'merger')
  bx.text(-17,17, 'in situ, old', size=10)
  bx.text(-17,-17, '$M_{*,tot}$=%.3g'%sum(dwarf[maino,12]),size=8)
  cx.scatter(dwarf[mainy,2]/h-xcen,dwarf[mainy,3]/h-ycen,marker = 'o',s = 10,c='b',edgecolors='none',label = 'merger')
  cx.text(-17,17, 'in situ, young', size=10)
  cx.text(-17,-17, '$M_{*,tot}$=%.3g'%sum(dwarf[mainy,12]),size=8)
  dx.scatter(dwarf[merger,2]/h-xcen,dwarf[merger,3]/h-ycen,marker = 'o',s = 10,c='b',edgecolors='none',label = 'merger')
  dx.text(-17,17, 'All merger', size=10)
  dx.text(-17,-17, '$M_{*,tot}$=%.3g'%sum(dwarf[merger,12]),size=8)
  ex.scatter(dwarf[mergero,2]/h-xcen,dwarf[mergero,3]/h-ycen,marker = 'o',s = 10,c='b',edgecolors='none',label = 'merger')
  #ex.text(-17,17, 'merger, early accre.', size=10)
  ex.text(-17,17, 'merger, old (>%.2fGyr)'%agecut, size=10)
  ex.text(-17,-17, '$M_{*,tot}$=%.3g'%sum(dwarf[mergero,12]),size=8)
  fx.scatter(dwarf[mergery,2]/h-xcen,dwarf[mergery,3]/h-ycen,marker = 'o',s = 10,c='b',edgecolors='none',label = 'merger')
  #fx.text(-17,17, 'merger, late accre.', size=10)
  fx.text(-17,17, 'merger, young(<%.2fGyr)'%agecut, size=10)
  fx.text(-17,-17, '$M_{*,tot}$=%.3g'%sum(dwarf[mergery,12]),size=8)

  pylab.setp(ax.get_xticklabels(), visible=False)
  pylab.setp(bx.get_xticklabels(), visible=False)
  pylab.setp(bx.get_yticklabels(), visible=False)
  pylab.setp(cx.get_xticklabels(), visible=False)
  pylab.setp(cx.get_yticklabels(), visible=False)
  pylab.setp(ex.get_yticklabels(), visible=False)
  pylab.setp(fx.get_yticklabels(), visible=False)

  ax.set_ylabel('Y Position (kpc)           ', fontsize=14)
  ex.set_xlabel('X Position (kpc)', fontsize=14)
  bx.set_title('Halo %s 9_17 XY Star Particle Positions at z=0'%hnum)
  ax.set_xlim(-20,20)
  ax.set_ylim(-20,20)
  fig.savefig('halo%snew%s_xy_starmerger_cutatz=2p1_%s.pdf'%(hnum,res,date),transparent=True)
  plt.show()

   

  return True 

def kristy_plots(dwarf,hnum,res,ver,halox,halo,halorad,xcen,ycen,zcen,thick):
  
  global h, PI,date
  starhalf = 0.83
  dwarf[:,11]=13.424-dwarf[:,11]
  inter = (dwarf[:,11]>0.5) & (dwarf[:,11]<3) #take only intermediately aged stars
  rstar = np.sqrt((dwarf[:,2]/h-xcen)**2+(dwarf[:,3]/h-ycen)**2+(dwarf[:,4]/h-zcen)**2)<3*starhalf #cut only stars within 3 halflight radii
  print len(dwarf[:,11]),len(dwarf[(inter) & (rstar),2]),len(dwarf[(inter) & (rstar),2])
  plt.figure()
  axes=pylab.axes()
  my_circle_scatter(axes, halox, haloy, radius=3*starhalf*np.ones(len(halox)), fill=False,color = 'b')
  pylab.axis('scaled')
  plt.scatter(dwarf[(inter) & (rstar),2]/h-xcen,dwarf[(inter) & (rstar),3]/h-ycen,marker = 'o',s = 10,edgecolors='none',label = 'intermediate')
  plt.xlabel('X Position (kpc)')
  plt.ylabel('Y Position (kpc)')
  plt.xlim(-5,5)
  plt.ylim(-5,5)
  plt.ticklabel_format(useOffset=False)
  plt.title("Halo %d XY Stellar Particle + Halo Positions at z=0"%hnum)
  plt.savefig('halo%d_xy_kristy_11_12.pdf'%hnum,transparent=True)
  #plt.show()
  plt.close()

  conv = 1
  x = dwarf[:,2]/h* conv
  y = dwarf[:,3]/h* conv
  z = dwarf[:,4]/h* conv
  mass = dwarf[:,12]
  xbins = np.logspace(-1,np.log10(thick),50)
  bindiff = np.log(xbins[1]/xbins[0])
  R = np.sqrt((x-xcen)**2+(y-ycen)**2)
  inter_stars, xbin_edges= np.histogram(R[(inter) & (rstar)], bins = xbins,weights=mass[(inter) & (rstar)])
  projR = xbins/np.exp(bindiff/2)
  projdeny = inter_stars/(4*PI*np.exp(np.log(xbins[1:])-(bindiff/2)))
  starhalf = 0.83
  i = find_nearest(projR,starhalf)
  Rhalf = projR[i]
  pylab.rcParams['xtick.major.pad']='6'
  pylab.rcParams['ytick.major.pad']='6'
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.loglog(projR[1:],projdeny, 'r',linewidth=4,label = 'Intermediate Age (0.5<t<3 Gyr)')
  ax.spines['bottom'].set_linewidth(4)
  ax.spines['top'].set_linewidth(4)
  ax.spines['right'].set_linewidth(4)
  ax.spines['left'].set_linewidth(4)
  ax.tick_params('both',length=5,width=2,which='minor')
  ax.tick_params('both',length=10,width=2,which='major')
  ax.xaxis.set_tick_params(labelsize=15)
  ax.yaxis.set_tick_params(labelsize=15)
  plt.xlabel('Projected Radius (in z) (kpc)',fontsize=20, labelpad=-10)
  ax.xaxis.set_label_coords(.48,-.07)
  plt.ylabel(r'$\Sigma (M_\odot$/kpc$^2$)',fontsize=20, labelpad=-5)
  plt.xlim(.1,100)
  #plt.ylim(1e3,1e9)
  plt.legend(loc=3,prop={'size':10})
  fname = 'projden_kristy_halo%d_11_12_z0.pdf'%hnum
  plt.savefig(fname,transparent=True)
  #plt.show()
  plt.close()    

  return True

def plot_projden(dwarf,xcen,ycen,zcen,thick, hnum,res,ver):

  global h, PI,date
  x = dwarf[:,2]/h
  y = dwarf[:,3]/h
  z = dwarf[:,4]/h
  mass = dwarf[:,12]
  xbins = np.logspace(-1,np.log10(thick),50)
  bindiff = np.log(xbins[1]/xbins[0])
  R = np.sqrt(((x-xcen)**2+(y-ycen)**2).astype('float'))
  young, xbin_edges= np.histogram(R[dwarf[:,11]>6.5], bins = xbins,weights=mass[dwarf[:,11]>6.5])
  old, xbin_edges= np.histogram(R[dwarf[:,11]<6.5], bins = xbins,weights=mass[dwarf[:,11]<6.5])
  young1, xbin_edges= np.histogram(R[(dwarf[:,11]>6.5) & (dwarf[:,8]==0)], bins = xbins,weights=mass[(dwarf[:,11]>6.5) & (dwarf[:,8]==0)])
  old1, xbin_edges= np.histogram(R[(dwarf[:,11]<6.5) & (dwarf[:,8]==0)], bins = xbins,weights=mass[(dwarf[:,11]<6.5) & (dwarf[:,8]==0)])
  young2, xbin_edges= np.histogram(R[(dwarf[:,11]>6.5) & (dwarf[:,8]==1)], bins = xbins,weights=mass[(dwarf[:,11]>6.5) & (dwarf[:,8]==1)])
  old2, xbin_edges= np.histogram(R[(dwarf[:,11]<6.5) & (dwarf[:,8]==1)], bins = xbins,weights=mass[(dwarf[:,11]<6.5) & (dwarf[:,8]==1)])
  projR = xbins/np.exp(bindiff/2)
  projdeny = young/(4*PI*np.exp(np.log(xbins[1:])-(bindiff/2)))
  projdeno = old/(4*PI*np.exp(np.log(xbins[1:])-(bindiff/2)))
  projdeny1 = young1/(4*PI*np.exp(np.log(xbins[1:])-(bindiff/2)))
  projdeno1 = old1/(4*PI*np.exp(np.log(xbins[1:])-(bindiff/2)))
  projdeny2 = young2/(4*PI*np.exp(np.log(xbins[1:])-(bindiff/2)))
  projdeno2 = old2/(4*PI*np.exp(np.log(xbins[1:])-(bindiff/2)))
  starhalf = 0.83#1.07#
  i = find_nearest(projR,starhalf)
  Rhalf = projR[i]
  pylab.rcParams['xtick.major.pad']='6'
  pylab.rcParams['ytick.major.pad']='6'
  fig = plt.figure()
  print projdeny
  ax = fig.add_subplot(111)
  ax.semilogy(projR[1:],projdeny/projdeny[i], 'r',linewidth=4,label = 'young (after 6.5 Gyr)')
  ax.semilogy(projR[1:],projdeno/projdeno[i], 'k',linewidth=4, label = 'old (before 6.5 Gyr)')
  ax.semilogy(projR[1:],projdeny1/projdeny1[i], 'r--',linewidth=4,label = 'young (merger)')
  ax.semilogy(projR[1:],projdeno1/projdeno1[i], 'k--',linewidth=4, label = 'old (merger)')
  ax.semilogy(projR[1:],projdeny2/projdeny2[i], 'r-.',linewidth=4,label = 'young (main)')
  ax.semilogy(projR[1:],projdeno2/projdeno2[i], 'k-.',linewidth=4, label = 'old (main)')
  ax.spines['bottom'].set_linewidth(4)
  ax.spines['top'].set_linewidth(4)
  ax.spines['right'].set_linewidth(4)
  ax.spines['left'].set_linewidth(4)
  ax.tick_params('both',length=5,width=2,which='minor')
  ax.tick_params('both',length=10,width=2,which='major')
  ax.xaxis.set_tick_params(labelsize=15)
  ax.yaxis.set_tick_params(labelsize=15)
  plt.xlabel('Projected Radius (in z) (kpc)',fontsize=20, labelpad=-10)
  ax.xaxis.set_label_coords(.48,-.07)
  plt.ylabel(r'$\Sigma (M_\odot$/kpc$^2$)',fontsize=20, labelpad=-5)
  plt.xlim(.1,8)
  plt.ylim(1e-3,100)
  plt.legend(loc=1,prop={'size':10})
  fname = 'projden_halo%dnew_11_12_z0.pdf'%hnum
  plt.savefig(fname,transparent=True)
  plt.show()
  plt.close()
  
  return True

def meanstellarage_vs_r(dwarf,xcen,ycen,zcen,thick, hnum):

  global h, hage
  x = dwarf[:,2]/h
  y = dwarf[:,3]/h
  z = dwarf[:,4]/h
  R = np.sqrt(((x-xcen)**2+(y-ycen)**2+(z-zcen)**2).astype('float'))
  xbins = np.logspace(-1,np.log10(thick),50)
  bindiff = np.log(xbins[1]/xbins[0])
  pR = xbins/np.exp(bindiff/2)
  allstar, xbin_edges= np.histogram(R, bins = xbins)
  allstar_age, xbin_edges= np.histogram(R, bins = xbins,weights=hage-dwarf[:,11])
  allstar_age2, xbin_edges= np.histogram(R, bins = xbins,weights=(hage-dwarf[:,11])**2)
  allstar_mean = allstar_age/allstar
  allstar_err = np.sqrt(allstar_age2/allstar-allstar_mean**2)
  mainstar, xbin_edges= np.histogram(R[dwarf[:,8]==1], bins = xbins)
  mainstar_age, xbin_edges= np.histogram(R[dwarf[:,8]==1], bins = xbins,weights=hage-dwarf[dwarf[:,8]==1,11])
  mainstar_age2, xbin_edges= np.histogram(R[dwarf[:,8]==1], bins = xbins,weights=(hage-dwarf[dwarf[:,8]==1,11])**2)
  mainstar_mean = mainstar_age/mainstar
  mainstar_err = np.sqrt(mainstar_age2/mainstar-mainstar_mean**2)
  mergestar, xbin_edges= np.histogram(R[dwarf[:,8]==0], bins = xbins)
  mergestar_age, xbin_edges= np.histogram(R[dwarf[:,8]==0], bins = xbins,weights=hage-dwarf[dwarf[:,8]==0,11])
  mergestar_age2, xbin_edges= np.histogram(R[dwarf[:,8]==0], bins = xbins,weights=(hage-dwarf[dwarf[:,8]==0,11])**2)
  mergestar_mean = mergestar_age/mergestar
  mergestar_err = np.sqrt(mergestar_age2/mergestar-mergestar_mean**2)
  print mergestar_mean,mergestar_err, mergestar
  plt.figure()
  plt.errorbar(np.log10(pR[1:]), allstar_mean, yerr=allstar_err, label='all')
  plt.errorbar(np.log10(pR[1:]), mainstar_mean, yerr=mainstar_err, label='in situ')
  plt.errorbar(np.log10(pR[1:]), mergestar_mean, yerr=mergestar_err, label='merger')
  plt.xlim(-1,1.2)
  plt.xlabel('Log(Radius) (kpc)')
  plt.ylabel(r'$\langle age\rangle$ (Gyr)')
  plt.title('Mean stellar age vs Radius')
  plt.legend(loc=4)
  plt.savefig('meanstarage_halo%dnew_11_12_z0.pdf'%hnum,transparent=True)
  plt.show()
  plt.close()

  return True

def merger_tracker(pathname,hnum,res,ver,dmo,snum,time,mergers,final_tab,final_star,size,rank,allahfdata,z,mvir_mhist,mstar_mhist):

	""" Creates an entire history of all the progenitors of the main halo."""

	mainhaloid, mainprogenid, progen, mainmergerid = gather_progen(pathname,snum,hnum)
	for i in np.arange(16):
		temp = glob.glob(pathname+'analysis/ahf_snap%03d/ahf.snap_%03d.%04d.z*.*.AHF_halos'%(snum-1,snum-1,i))
		temp = str(temp).strip('[]').replace("'","")
		temp = np.genfromtxt(temp)
		if len(temp) != 0:
			try:
				ahfdata = np.vstack((ahfdata,temp))
			except:
				ahfdata = temp
	if len(ahfdata) != 0 and mainprogenid != 0:
		progen = ahfdata[np.in1d(ahfdata[:,0],progen)]
		rvir = progen[progen[:,0]==mainprogenid,11]
		print 'rvir is ',rvir
		xcen = progen[progen[:,0]==mainprogenid,5]
		ycen = progen[progen[:,0]==mainprogenid,6]
		zcen = progen[progen[:,0]==mainprogenid,7]
		dist_to_main = np.sqrt((progen[:,5]-xcen)**2+(progen[:,6]-ycen)**2+(progen[:,7]-zcen)**2)
		progen = progen[:,[0,0,3,64,3,3,64,3,37,3]]
		progen[:,0] = time
		progen[progen[:,1]==mainprogenid,4] = 1
		progen[progen[:,1]==mainmergerid,4] = 2
		mvir = progen[progen[:,1]==mainprogenid,2]
		mstar = progen[progen[:,1]==mainprogenid,3]
		progen[progen[:,1]!=mainprogenid,4] = 0
		progen[:,5] = progen[:,5]/mvir
		if mstar > 0:
			progen[:,6] = progen[:,6]/mstar
		else:
			progen[progen[:,3]>0,6] = 2 
			progen[progen[:,3]==0,6] = 0 
		progen[:,7] = dist_to_main
		progen[:,9] = 0

		if rank == 0:
			data = progen[((progen[:,2]>.01*progen[progen[:,1]==mainprogenid,2])|(progen[:,3]>0)) & ((progen[:,7]>rvir)|(progen[:,7]==0)),1]
			print 'dist to main ', progen[(progen[:,2]>.01*progen[progen[:,1]==mainprogenid,2]) & ((progen[:,7]>rvir)|(progen[:,7]==0)),7]
			outputData = np.zeros(len(data)) #Create output array of same size
			split = np.array_split(data,size,axis = 0) #Split input array by the number of available cores

			split_sizes = []

			for i in range(0,len(split),1):
				split_sizes = np.append(split_sizes, len(split[i]))
			split_sizes_input = split_sizes
			displacements_input = np.insert(np.cumsum(split_sizes),0,0)[0:-1]
			split_sizes_output = split_sizes
			displacements_output = np.insert(np.cumsum(split_sizes),0,0)[0:-1]


			print("Input data split into vectors of sizes %s" %split_sizes)
			print("Input data split with displacements of %s" %displacements_input)

		else:
			#Create variables on other cores
			split_sizes_input = None
			displacements_input = None
			split_sizes_output = None
			displacements_output = None
			split = None
			data = None
			outputData = None

		split = comm.bcast(split, root=0) #Broadcast split array to other cores
		split_sizes = comm.bcast(split_sizes_input, root = 0)
		displacements = comm.bcast(displacements_input, root = 0)
		split_sizes_output = comm.bcast(split_sizes_output, root = 0)
		displacements_output = comm.bcast(displacements_output, root = 0)		

		output_chunk = np.zeros(np.shape(split[rank])) #Create array to receive subset of data on each core, where rank specifies the core
		print("Rank %d with output_chunk shape %s" %(rank,output_chunk.shape))
		print data,split_sizes_input, displacements_input,output_chunk
		comm.Scatterv([data,split_sizes_input, displacements_input,MPI.DOUBLE],output_chunk,root=0)
		output = np.zeros(len(output_chunk))

		#comm.Scatterv([data,sendcounts,displs,MPI.DOUBLE],block,root=0)
		print 'Success'
		for q,haloid in enumerate(output_chunk): #only do merger tree/peak mass analysis for big enough halos
			print snum,' q is ',q
			peaklist,snums = peakmass(hnum,res,ver,dmo,haloid,snum-1)
			snums = snums[::-1].astype(int)
			peaklist = peaklist[::-1]
			print 'peaklist is',len(peaklist)
			w = 0
			while peaklist[w] == 0:
				w+=1
			print 'w is ',w
			output[q] = max(allahfdata[np.in1d(np.round(allahfdata[:,0],10),np.round(peaklist,10)),1])
			mvirdata = allahfdata[np.in1d(np.round(allahfdata[:,0],10),np.round(peaklist,10)),1]
			mstardata = allahfdata[np.in1d(np.round(allahfdata[:,0],10),np.round(peaklist,10)),2]
			print 'mvirdata is',len(mvirdata)
			if hnum == '10q' or hnum == '10v':
				if len(mvirdata) == 601 and mvirdata[len(mvirdata)-1]>6e9:
					peakpoint = 600
				else:
					peakpoint = len(mvirdata)-1#mvirdata.argmax() 
			else:
				if len(mvirdata) > 135 and mvirdata[len(mvirdata)-1]>6e9: #This prevents the main progenitor of 11707 or 848 from being cut off too early.
					peakpoint = len(mvirdata)-1
				else:
					peakpoint = len(mvirdata)-1#mvirdata.argmax()

			if output_chunk.shape[0] > 0 and rank != 0:
				rankrank = int('%d%d'%(rank,rank))
				rankrankrank = int('%d%d%d'%(rank,rank,rank))
				rankrankrankrank = int('%d%d%d%d'%(rank,rank,rank,rank))
				#comm.Send([z[(snum+1-len(mvirdata[:peakpoint+1])):(snum+1)],MPI.DOUBLE],dest=0, tag = rank)
				comm.Send([z[w:len(mvirdata[:peakpoint+1])+w],MPI.DOUBLE],dest=0, tag = rank)
				comm.Send([mvirdata[:peakpoint+1],MPI.DOUBLE],dest=0, tag = rankrank)
				comm.Send([mstardata[:peakpoint+1],MPI.DOUBLE],dest=0, tag = rankrankrank)
				comm.send(len(mvirdata[:peakpoint+1]),dest=0,tag = rankrankrankrank)

		comm.Gatherv(output,[outputData,split_sizes_output,displacements_output,MPI.DOUBLE], root=0) #Gather output data together
		if rank == 0:
			KEYS = final_tab.keys()
			KEYS_array = np.zeros((2,len(KEYS)))
			KEYS_star = final_star.keys()
			KEYS_star_array = np.zeros((2,len(KEYS_star)))
			for i in np.arange(len(KEYS)):
				KEYS_array[0,i] = float(KEYS[i].split(',')[0]) #Contains all the z's of the merger points so far
				KEYS_array[1,i] = float(KEYS[i].split(',')[1]) #Contains all the mvir's of the merger points so far
			for i in np.arange(len(KEYS_star)):
				KEYS_star_array[0,i] = float(KEYS_star[i].split(',')[0]) #Contains all the z's of the merger points so far
				KEYS_star_array[1,i] = float(KEYS_star[i].split(',')[1]) #Contains all the mvir's of the merger points so far

			for i in np.linspace(1,len(data)-1,len(data)-1): #Recieve merger histories from other computers
				lenlen = comm.recv(source=int('%d'%i),tag=int('%d%d%d%d'%(i,i,i,i)))
				z_for_plot = np.empty(lenlen,dtype='d')
				mvir_for_plot = np.empty(lenlen,dtype='d')
				mstar_for_plot = np.empty(lenlen,dtype='d')
				comm.Recv([z_for_plot,MPI.DOUBLE],source=int('%d'%i),tag=int('%d'%i))
				comm.Recv([mvir_for_plot,MPI.DOUBLE],source=int('%d'%i),tag=int('%d%d'%(i,i)))
				comm.Recv([mstar_for_plot,MPI.DOUBLE],source=int('%d'%i),tag=int('%d%d%d'%(i,i,i)))

				for j in np.arange(len(KEYS)):
					if np.round(KEYS_array[0,j],5) in np.round(z_for_plot,5) and np.round(KEYS_array[1,j],5) in np.round(mvir_for_plot,5): #Are any of these new merger points just continuations from an old point? 
						final_tab.pop('%s,%s'%(KEYS_array[0,j],KEYS_array[1,j])) #Remove old merger point from final tab. 
						KEYS_array[1,j] = mvir_for_plot[-1] #Replace old merger point with latest merger point.
						KEYS_array[0,j] = z_for_plot[-1]
						final_tab['%s,%s'%(KEYS_array[0,j],KEYS_array[1,j])] = np.vstack((z_for_plot,mvir_for_plot)) #Make new home for merger point and merger history
					else: #Brand new merger
						final_tab['%s,%s'%(z_for_plot[-1],mvir_for_plot[-1])] = np.vstack((z_for_plot,mvir_for_plot)) #Make new home for merger point and merger history
				if len(final_tab) == 0:
					final_tab['%s,%s'%(z_for_plot[-1],mvir_for_plot[-1])] = np.vstack((z_for_plot,mvir_for_plot)) #Make new home for merger point and merger history

				#Now the same procedure but for stars
				for j in np.arange(len(KEYS_star)):
					if np.round(KEYS_star_array[0,j],5) in np.round(z_for_plot,5) and np.round(KEYS_star_array[1,j],5) in np.round(mstar_for_plot,5) and np.round(KEYS_star_array[1,j],5) > 0: #Are any of these new merger points just continuations from an old point? 
						final_star.pop('%s,%s'%(KEYS_star_array[0,j],KEYS_star_array[1,j])) #Remove old merger point from final tab. 
						KEYS_star_array[1,j] = mstar_for_plot[-1] #Replace old merger point with latest merger point.
						KEYS_star_array[0,j] = z_for_plot[-1]
						final_star['%s,%s'%(KEYS_star_array[0,j],KEYS_star_array[1,j])] = np.vstack((z_for_plot,mstar_for_plot)) #Make new home for merger point and merger history
					else: #Brand new merger
						final_star['%s,%s'%(z_for_plot[-1],mstar_for_plot[-1])] = np.vstack((z_for_plot,mstar_for_plot)) #Make new home for merger point and merger history
				if len(final_star) == 0:
					final_star['%s,%s'%(z_for_plot[-1],mstar_for_plot[-1])] = np.vstack((z_for_plot,mstar_for_plot)) #Make new home for merger point and merger history


			z_for_plot = z[w:len(mvirdata[:peakpoint+1])+w] #Now repeat previous KEY loop except for the data on rank 0 proc.
			mvir_for_plot = mvirdata[:peakpoint+1]
			mstar_for_plot = mstardata[:peakpoint+1]
			for j in np.arange(len(KEYS)): 
				if np.round(KEYS_array[0,j],5) in np.round(z_for_plot,5) and np.round(KEYS_array[1,j],5) in np.round(mvir_for_plot,5):
					print final_tab 
					final_tab.pop('%s,%s'%(KEYS_array[0,j],KEYS_array[1,j]))  
					KEYS_array[1,j] = mvir_for_plot[-1] 
					KEYS_array[0,j] = z_for_plot[-1]
					final_tab['%s,%s'%(KEYS_array[0,j],KEYS_array[1,j])] = np.vstack((z_for_plot,mvir_for_plot)) 
				else: #Brand new merger
					final_tab['%s,%s'%(z_for_plot[-1],mvir_for_plot[-1])] = np.vstack((z_for_plot,mvir_for_plot)) #Make new home for merger point and merger history
			if len(final_tab) == 0:
				final_tab['%s,%s'%(z_for_plot[-1],mvir_for_plot[-1])] = np.vstack((z_for_plot,mvir_for_plot)) #Make new home for merger point and merger history

			for j in np.arange(len(KEYS_star)):
				if np.round(KEYS_star_array[0,j],5) in np.round(z_for_plot,5) and np.round(KEYS_star_array[1,j],5) in np.round(mstar_for_plot,5) and np.round(KEYS_star_array[1,j],5) > 0: #Are any of these new merger points just continuations from an old point? 
					final_star.pop('%s,%s'%(KEYS_star_array[0,j],KEYS_star_array[1,j])) #Remove old merger point from final tab. 
					KEYS_star_array[1,j] = mstar_for_plot[-1] #Replace old merger point with latest merger point.
					KEYS_star_array[0,j] = z_for_plot[-1]
					final_star['%s,%s'%(KEYS_star_array[0,j],KEYS_star_array[1,j])] = np.vstack((z_for_plot,mstar_for_plot)) #Make new home for merger point and merger history
				else: #Brand new merger
					final_star['%s,%s'%(z_for_plot[-1],mstar_for_plot[-1])] = np.vstack((z_for_plot,mstar_for_plot)) #Make new home for merger point and merger history
			if len(final_star) == 0:
				final_star['%s,%s'%(z_for_plot[-1],mstar_for_plot[-1])] = np.vstack((z_for_plot,mstar_for_plot)) #Make new home for merger point and merger history

			#mvir_mhist.add_line(z[w:len(mvirdata[:peakpoint+1])+w],mvirdata[:peakpoint+1])
			#mstar_mhist.add_line(z[w:len(mvirdata[:peakpoint+1])+w],mstardata[:peakpoint+1])
			progen[((progen[:,2]>.01*progen[progen[:,1]==mainprogenid,2])|(progen[:,3]>0)) & ((progen[:,7]>rvir)|(progen[:,7]==0)),9] = outputData
			try:
				mergers = np.vstack((mergers,progen))
				print "length is ",len(mergers)
			except Exception,e:
				print e
				mergers = progen
	return mergers,final_tab,final_star

class merger_hist(object):
	"""Merger history plot """
	def __init__(self):
		self.fig = plt.figure(figsize=plt.figaspect(.9))
		self.sub = self.fig.add_subplot(111)
		fs = 20
		twin = self.sub.twinx()
		twin.set_ylim(1e-5,1e1)
		twin.set_yscale('log')
		twin.set_ylabel(r'$M_{\rm \star,merger}/M_{\rm \star,main}$',fontsize=fs, labelpad=5)
		#self.sm = sm
		self.sub.spines['bottom'].set_linewidth(4)
		self.sub.spines['top'].set_linewidth(4)
		self.sub.spines['right'].set_linewidth(4)
		self.sub.spines['left'].set_linewidth(4)
		self.sub.tick_params('both',length=5,width=2,which='minor')
		self.sub.tick_params('both',length=10,width=2,which='major')
		self.sub.xaxis.set_tick_params(labelsize=fs)
		self.sub.yaxis.set_tick_params(labelsize=fs)
		self.sub.xaxis.set_label_coords(.48,-.07)
		self.sub.set_xlabel(r'Time (Gyr)',fontsize=fs, labelpad=-10)
		self.sub.set_ylabel(r'$M_{\rm vir,merger}/M_{\rm vir,main}$',fontsize=fs, labelpad=5)
		self.sub.set_xlim(0,13.8)
		self.sub.set_ylim(1e-5,1e1)
	def add_line(self,mergers):
		main = np.where(mergers[:,4]==1)[0]
		t = mergers[main,0]
		topmerge = np.zeros(len(main))
		topmerge = np.where(mergers[:,4]==2)[0]
		if len(topmerge)<len(main):
			topmerge = np.append(topmerge[0],topmerge)
		#for i in np.arange(len(main)-1):
		#	if main[i+1]-(main[i]+1)>0:
		#		topmerge[i] = mergers[mergers[main[i]+1:main[i+1],4]==2]
		#for i in np.arange(len(main)-1):
		#	if main[i+1]-(main[i]+1)>0:
		#		topmerge[i] = max(mergers[main[i]+1:main[i+1],5])
		self.sub.semilogy(t,topmerge[:,5],'k')
		self.sub.scatter(mergers[(mergers[:,3]!=0) & (mergers[:,4]!=1),0],mergers[(mergers[:,3]!=0) & (mergers[:,4]!=1),6],color='r',zorder=100)
	def save(self,hnum,date):
		#self.sub.legend(loc=3,prop={'size':10})
		self.sub.text(8,3,'Halo %s'%hnum)
		self.fig.savefig('Merger_hist_%s_%s.pdf'%(hnum,date),transparent=True)#radden_halo%s_%dbin_%s.pdf'%(hnum,grain,date),transparent=True)

class mvir_progen(object):
	""" Mvir for all progen """
	def __init__(self):
		fs = 20
		self.fig = plt.figure(figsize=plt.figaspect(.9))
		self.sub = self.fig.add_subplot(111)
		self.sub.set_xlabel(r'z',fontsize=fs, labelpad=0)
		self.sub.set_ylabel(r'$M_{\rm vir}$',fontsize=fs, labelpad=5)
		self.sub.set_xlim(1,0)
		self.sub.set_ylim(1e7,2e10)
		plt.xticks([1,0],['9','0'],fontsize=fs)
		minor_ticks=[]
		for j in np.linspace(9,2,8):
		  minor_ticks.append(0+np.log10(j))
		self.sub.xaxis.set_minor_locator(FixedLocator(minor_ticks))
		minor_labels = ['8','7','6','5','4','3','2','1']
		self.sub.xaxis.set_minor_formatter(FixedFormatter(minor_labels))
	def add_line(self,z,mvirs):
		print 'hi'
		self.sub.semilogy(np.log10(z),mvirs,'k',linewidth = 1)
		#if len(mvirs) != 185:
		#	self.sub.scatter(np.log10(z[len(z)-1]),mvirs[len(mvirs)-1],color='b',s=40)
	def save(self,hnum,date):
		self.fig.savefig('Mvir_merger_hist_%s_%s.pdf'%(hnum,date),transparent=True)
class mvirmstar_progen(object):
	""" Mvir for all progen """
	def __init__(self):
		fs = 20
		self.fig = plt.figure(figsize=plt.figaspect(.9))
		self.sub = self.fig.add_subplot(111)
		self.sub2 = self.sub.twinx()
		self.sub.set_xlabel(r'z',fontsize=fs, labelpad=0)
		self.sub.set_ylabel(r'$M_{\rm vir}$',fontsize=fs, labelpad=5)
		self.sub2.set_ylabel(r'$M_{\star}$',fontsize=fs, labelpad=5)
		self.sub.set_xlim(1,0)
		self.sub.set_ylim(1e7,2e10)
		self.sub2.set_ylim(1e3,2e7)
		plt.xticks([1,0],['9','0'],fontsize=fs)
		minor_ticks=[]
		for j in np.linspace(9,2,8):
		  minor_ticks.append(0+np.log10(j))
		self.sub.xaxis.set_minor_locator(FixedLocator(minor_ticks))
		minor_labels = ['8','7','6','5','4','3','2','1']
		self.sub.xaxis.set_minor_formatter(FixedFormatter(minor_labels))
	def add_virline(self,z,mvirs):
		self.sub.semilogy(np.log10(z),mvirs/.71,'k',linewidth = 3)
	def add_starline(self,z,mstars):
		self.sub2.semilogy(np.log10(z),mstars/.71,'k',linewidth = 2, color='darkgray')
	def add_point(self,z,mvirs,mainz,mainvir,ver):
		if ver == 0:
			ratio_index = np.where(mainz == z[-1])[0][0]
			self.sub.scatter(np.log10(z[-1]),mvirs[-1]/.71,color='b',s=80+30*np.log10(mvirs[len(mvirs)-1]/mainvir[ratio_index]))
		if ver == 1:
			ratio_index = np.where(mainz == z[-1])[0][0]
			self.sub2.scatter(np.log10(z[-1]),mvirs[-1]/.71,color='r',marker='*',s=80+30*np.log10(mvirs[len(mvirs)-1]/mainvir[ratio_index]))
	def save(self,hnum,date):
		self.fig.savefig('Mvir_merger_hist_%s_%s.pdf'%(hnum,date),transparent=True)

class mstar_progen(object):
	""" Mstar for all progen """
	def __init__(self):
		fs = 20
		self.fig = plt.figure(figsize=plt.figaspect(.9))
		self.sub = self.fig.add_subplot(111)
		self.sub.set_xlabel(r'z',fontsize=fs, labelpad=0)
		self.sub.set_ylabel(r'$M_{\star}$',fontsize=fs, labelpad=5)
		self.sub.set_xlim(1,0)
		self.sub.set_ylim(1e3,2e7)
		plt.xticks([1,0],['9','0'],fontsize=fs)
		minor_ticks=[]
		for j in np.linspace(9,2,8):
		  minor_ticks.append(0+np.log10(j))
		self.sub.xaxis.set_minor_locator(FixedLocator(minor_ticks))
		minor_labels = ['8','7','6','5','4','3','2','1']
		self.sub.xaxis.set_minor_formatter(FixedFormatter(minor_labels))
	def add_line(self,z,mstars):
		self.sub.semilogy(np.log10(z),mstars,'k',linewidth = 1)
		#if len(mstars) != 185:
		#	self.sub.scatter(np.log10(z[len(z)-1]),mstars[len(mstars)-1],color='b',s=40)
	def save(self,hnum,date):
		self.fig.savefig('Mstar_merger_hist_%s_%s.pdf'%(hnum,date),transparent=True)

class merger_histogram(object):
	""" Multi-layered histogram that depicts when mergers are occuring"""
	def ___init___(self):
		fs = 20
		self.fig = plt.figure(figsize=plt.figaspect(.9))
		self.sub = self.fig.add_subplot(111)
		self.sub.set_xlabel(r'z',fontsize=fs, labelpad=0)
		self.sub.set_ylabel(r'Number of Mergers',fontsize=fs, labelpad=5)
		self.sub.set_xlim(13.8,0)
		#self.sub.set_ylim(1e3,2e7)
		#plt.xticks([1,0],['9','0'],fontsize=fs)
		#minor_ticks=[]
		#for j in np.linspace(9,2,8):
		#  minor_ticks.append(0+np.log10(j))
		#self.sub.xaxis.set_minor_locator(FixedLocator(minor_ticks))
		#minor_labels = ['8','7','6','5','4','3','2','1']
		#self.sub.xaxis.set_minor_formatter(FixedFormatter(minor_labels))
	def add_hist(self,z,zmain,virmain,starmain):
		h0 = 71
		om_l = 0.734
		om_m = 0.266
		conv = 3.085677581e+19
		above50 = []
		above10 = []
		for i in np.arange(len(z)):
			if mvir[i] > 1./50*virmain[np.in1d(zmain,z[i])]:
				above50.append(z[i])
			if mvir[i] > 1./10*virmain[np.in1d(zmain,z[i])]:
				above10.append(z[i])
		time = np.zeros(len(z))
		for i in np.arange(len(z)):
			time[i] = scipy.integrate.quad(_a_dot_recip, 0,1/z[i], (h0, om_m, om_l))[0]*conv/3600/24/365/1e9
		hist,bin_edges = np.histogram(time,bins=[3,6,9,12,15])
		self.sub.bar(bin_edges,hist,color = 'blue',edgecolor = 'black')
		time = np.zeros(len(above50))
		for i in np.arange(len(above50)):
			time[i] = scipy.integrate.quad(_a_dot_recip, 0,1/above50[i], (h0, om_m, om_l))[0]*conv/3600/24/365/1e9
		hist,bin_edges = np.histogram(time,bins=[3,6,9,12,15])
		self.sub.bar(bin_edges,hist,color = 'blue',edgecolor = 'black',hatch = '//')
		time = np.zeros(len(above10))
		for i in np.arange(len(above10)):
			time[i] = scipy.integrate.quad(_a_dot_recip, 0,1/above10[i], (h0, om_m, om_l))[0]*conv/3600/24/365/1e9
		hist,bin_edges = np.histogram(time,bins=[3,6,9,12,15])
		self.sub.bar(bin_edges,hist,color = 'blue',edgecolor = 'black',hatch = '\\')
	def save(self,date):
		self.fig.savefig('mergers_histogram_%s.pdf'%(date),tranparent=True)

class mergers_vs_mstar(object):
	""" Merger rate (ver 0) or In-situ/merger SF (ver 1) vs. Mstar"""
	def ___init___(self,ver):
		fs = 20
		self.fig = plt.figure(figsize=plt.figaspect(.9))
		self.sub = self.fig.add_subplot(111)
		self.sub.set_xlabel(r'Time (Gyr)',fontsize=fs, labelpad=0)
		#self.sub.set_xlim(0,0)
		if ver == 0:
			self.sub.set_ylabel(r'$M_{\star}$',fontsize=fs, labelpad=5)
			self.sub.set_ylim(1e3,2e7)
		if ver == 1:
			self.sub.set_ylabel(r'$M_{\star}$',fontsize=fs, labelpad=5)
			self.sub.set_ylim(1e3,2e7)
	def add_point(self,mstar,ver):
		print 'filler'

	def save(self,date):
		self.fig.savefig('mergers_histogram_%s.pdf'%(date),tranparent=True)

def peakmass(hnum,res,ver,dmo,haloid,snapf):

	#Import some files
	if dmo == 1:
		path_pref = "/nobackup/afitts/Gadget-2.0.7/production/gizdm%s%s_raw_output/analysis/merger"%(hnum,res)
	else:
		path_pref = "/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/analysis/merger"%(hnum,res,ver)
	numoffiles = 16
	setname = "MWa"
	setpref = "halo%s"%hnum
	mtree_pref= "halo%s_snap"%hnum #halo12596_13_snap184_snap183_mtree_idx
	totsnap = snapf


	pathname = path_pref+"/"
	candidate_file = haloid
	filein_pref = setpref #wdm_s7_bin_200.0000.z0.000.AHF_halos
	fileout_oref = setpref
	#candidate_file = open(pathname+"wdm_s7_bin_200_LGsubs_MWb_v2.txt")
	#fileout_pref = pathname + 'density_profiles_L7_snap200_'


	candidates = []
	cand_id_list = []
	candidates.append(0)#not halo id, but number sequentially from 0
	cand_id_list.append(haloid)#halo id

	#need to do this for host first? - identify host via r_dist

	cand_float = np.array(cand_id_list,dtype=np.float64)
	cand_id = cand_float.astype(np.uint64)
	#cand_int = np.uint64(cand_float)
	#cand_id = np.array(cand_int,dtype=np.uint64)
	num_orig_cand = len(cand_id)
	if hnum == '10q' or hnum == '10v':
		all_ids = np.zeros((num_orig_cand,601)) #includes host
	else:
		all_ids = np.zeros((num_orig_cand,185)) #includes host
	all_ids[:,0] = cand_id 
	print cand_id
	#snapf = 184
	#snapi = 183
	snapi = snapf-1
	i = 0
	snoops = np.linspace(snapf,0,snapf+1)
	terp = []
	terp.append(str(cand_id_list).strip("[']"))
	while snapi >= 1:
		if snapi >= 100:
		        snapi_str = str(snapi)
		elif (snapi < 100) & (snapi >= 10):
		        snapi_str = '0'+str(snapi)
		else:
		        snapi_str = '00'+str(snapi)
		if snapf >= 100:
		        snapf_str = str(snapf)
		elif (snapf < 100) & (snapf >= 10):
		        snapf_str = '0'+str(snapf)
		else:
		        snapf_str = '00'+str(snapf)
		mtree_file = open(pathname+mtree_pref+snapf_str+"_snap"+snapi_str+"_mtree_idx")
		print snapf,snapi
		mtree = np.loadtxt(mtree_file,dtype=np.uint64)
		#print np.shape(mtree)
		consider = all_ids[:,i]
		#print consider
		nonwhere_ind = np.nonzero(consider)
		#print nonwhere_ind
		#print nonwhere_ind[0][3]
		k = np.count_nonzero(consider)
		print k
		i += 1
		if k > 0:
			for a in np.arange(k):
				halo = consider[a]
		                #print halo
				try:
				  parent_ind = np.nonzero(halo == mtree[:,0])
				except:
				  break;
		                #print parent_ind[0]
		                parent = mtree[parent_ind[0],1]
		                #print parent
		                if len(parent) > 0:
		                        all_ids[nonwhere_ind[0][a],i] = parent
					terp.append(str(parent).strip('[]'))
			snapi -= 1
			snapf -= 1
			finalsnap = snapf
		else:
			finalsnap = snapf
			snapi = 0
		try:
			print parent
		except Exception,E:
			print E
	while len(terp)<totsnap+1:
	  terp.append(0)
	return np.array(terp).astype(np.float),snoops
	#bf = open(pathname+setpref+'merger_hist.txt','w')
	#for i in np.arange(len(terp)):
	#  bf.write('%d %s\n'%((snoops[i],terp[i]))) 
	#bf.close()


PI = 3.14159265359
h = 0.71
hnum =['948']#,'20910','32257','32503']#'948','20910','20192','007','1016']#'848','32503','32257','10q','10v','11707','12596','796','897']#'11707','32503','12596','007','848','796','20192','20910','897','948']#
res = ['_13','','','','_13','_13','_13','_13','_13','_13','_13','_11','_11','_13']
ver = ['5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16']
dmo = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
snum = [185,185,185,185,185,185,185,185,185,185,185,601,601,185]
date = time.strftime("%m_%d_%Y")
conv = 1000
snap = 184
hage = 13.424
#for j in np.linspace(100,184,85):
#  pathname =  '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/'%(hnum,res[0],ver[0])
#  starx,stary,starz,halox,haloy,haloz,halorad,xcen,ycen,zcen = get_starsnhalos(pathname,hnum,res[0],ver[0],j)
#  plot_star_and_halo(starx,stary,starz,halox,haloy,haloz,halorad,xcen,ycen,zcen,hnum,res[0])

switch = 0 #switch that determines if a main progenitor with stars has been found. When it is first found, all star particles will be denoted as from the main progenitor. From then on switch =1 and the complied list of particles will be carried through each snapshot.
main_progen = np.zeros(1)
mergers = []
fnum = 16
#mhist_plot= merger_hist()
mvir_mhist = mvir_progen()
mstar_mhist = mstar_progen()


def parallel_main(hnum,res,ver,snum):
	mergers = []
	final_tab = {}
	final_star = {}
	pathname =  '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/'%(hnum,res,ver)
	mainid_file = np.genfromtxt(pathname+'analysis/halo%smerger_hist.txt'%hnum)
	if hnum == '10v' or hnum == '10q':
		time_file = np.genfromtxt('/nobackup/afitts/analysis_scripts/Halo10q_11_mstar.out')
	else:
		time_file = np.genfromtxt('/nobackup/afitts/analysis_scripts/Halo848_13_mstar.out')
	allahfdata = pd.read_csv(pathname+'analysis/all_ahf_halo%s.out'%hnum,sep=' ',header=None)
	allahfdata = allahfdata.values.astype('float')
	if hnum == '10q' or hnum == '10v':
		a = np.genfromtxt('/nobackup/afitts/analysis_scripts/snapshot_scale-factors.txt')
	else:
		a = np.genfromtxt('/nobackup/afitts/analysis_scripts/output_times.txt')
	z = 1./a
	for i,e in enumerate(np.arange(snum)):
		print e
		if e != 0:
			mainid = mainid_file[snum-e,1]
		else:
			mainid = mainid_file[snum-e-1,1]
		if mainid != 0:
			#main_progen = gather_particles(main_progen,pathname,e,hnum,fnum)
			current_time = time_file[i,0]
			print current_time
			#hdf.close()
			mergers, final_tab, final_star = merger_tracker(pathname,hnum,res,ver,dmo,e,current_time,mergers,final_tab,final_star,size,rank,allahfdata,z,mvir_mhist,mstar_mhist)
	if rank == 0:
		np.savetxt('Halo%s_mergers1.out'%(hnum),mergers,header='(0)Time (Gyr)\t(1)Halo ID\t(2)M_vir\t(3)M_star\t(4)Main Progenitor (0=no,1=yes,2=Main Merger)\t(5)M/M_vir,mainprogen\t(6)M/M_star,mainprogen\t(7)Distance to Main Progenitor (kpc)\t(8) Frac of Mass in high res\t(9) M_peak')
		save_obj(final_tab,'Mvir%s_hist_allprogen_at'%(hnum))
		save_obj(final_star,'Mstar%s_hist_allprogen_at'%(hnum))
	#mergers = np.array(mergers)
	#print mergers
	return True
def plotter(hnum):
	for i in np.arange(len(hnum)):
		mergers = np.genfromtxt('/nobackup/afitts/analysis_scripts/Halo%s_mergers.out'%hnum[i])
		mhist_plot = merger_hist()
		mhist_plot.add_line(mergers)
		mhist_plot.save(hnum[i],date)
	return True
def mhist_plots(hnum):
	for i in np.arange(len(hnum)):
		print hnum[i]
		mvirmstar_mhist = mvirmstar_progen()
		raw = load_obj('Mvir%s_hist_allprogen_at'%(hnum[i]))
		bob = np.zeros((2,len(raw)))
		for j in np.arange(len(raw)):
			bob[0,j]= float(raw.keys()[j].split(',')[0])
			bob[1,j]= float(raw.keys()[j].split(',')[1])
		maxpoint = bob[1,:].argmax() # Find the main progenitor by searching for the largest endpoint
		raw_star = load_obj('Mstar%s_hist_allprogen_at'%(hnum[i]))
		bob_star= np.zeros((2,len(raw_star)))
		for j in np.arange(len(raw_star)):
			bob_star[0,j]= float(raw_star.keys()[j].split(',')[0])
			bob_star[1,j]= float(raw_star.keys()[j].split(',')[1])
		maxpoint_star = bob_star[1,:].argmax() # Find the main progenitor by searching for the largest endpoint
		####MHIST PLOTS####
		for k in np.arange(len(bob[0])):
			mvirmstar_mhist.add_virline(raw['%s,%s'%(bob[0,k],bob[1,k])][0],raw['%s,%s'%(bob[0,k],bob[1,k])][1])
			#mvirmstar_mhist.add_point(raw['%s,%s'%(bob[0,k],bob[1,k])][0],raw['%s,%s'%(bob[0,k],bob[1,k])][1],raw['%s,%s'%(bob[0,maxpoint],bob[1,maxpoint])][0],raw['%s,%s'%(bob[0,maxpoint],bob[1,maxpoint])][1],0)
		for k in np.arange(len(bob_star[0])):
			mvirmstar_mhist.add_point(raw_star['%s,%s'%(bob_star[0,k],bob_star[1,k])][0],raw_star['%s,%s'%(bob_star[0,k],bob_star[1,k])][1],raw_star['%s,%s'%(bob_star[0,maxpoint_star],bob_star[1,maxpoint_star])][0],raw_star['%s,%s'%(bob_star[0,maxpoint_star],bob_star[1,maxpoint_star])][1],1)
		mvirmstar_mhist.add_virline(raw['%s,%s'%(bob[0,maxpoint],bob[1,maxpoint])][0],raw['%s,%s'%(bob[0,maxpoint],bob[1,maxpoint])][1])
		mvirmstar_mhist.add_starline(raw_star['%s,%s'%(bob_star[0,maxpoint_star],bob_star[1,maxpoint_star])][0],raw_star['%s,%s'%(bob_star[0,maxpoint_star],bob_star[1,maxpoint_star])][1])		
		#mstar_mhist.add_line(z[w:len(mvirdata[:peakpoint+1])+w],mstardata[:peakpoint+1])
		mvirmstar_mhist.save(hnum[i],date)
		

def starmerge_plots(hnum,res,ver):
	for j in np.arange(len(hnum)):
		main_progen = np.genfromtxt('Halo848_starmerger_test_11_18_2016.out')
		pathname =  '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/'%(hnum[j],res[j],ver[j])
		hdf = pd.HDFStore('%sanalysis/dataframes/halo%s%s_giz%s_snap%03d.h5'%(pathname,hnum[j],res[j],ver[j],184))
		xcen = np.float(hdf['props']['halox'])*1e3
		ycen = np.float(hdf['props']['haloy'])*1e3
		zcen = np.float(hdf['props']['haloz'])*1e3
		halorad = np.float(hdf['props']['rvir'])*1000
		plot_starmerge(main_progen,hnum[j],res[j],xcen,ycen,halorad,xcen,ycen)

plt.ion()
j = int(sys.argv[1])
print j

parallel_main(hnum[j],res[j],ver[j],snum[j])
#mvir_mhist.save(hnum[j],date)
#mstar_mhist.save(hnum[j],date)
mhist_plots(hnum)
#starmerge_plots(hnum,res,ver)
#plotter(hnum)
#mhist_plot.add_line(mergers)
#mhist_plot.save(date)
#main_progen[:,0] = main_progen[:,0].astype('int')
#main_progen[:,2] = main_progen[:,2]*conv
#main_progen[:,3] = main_progen[:,3]*conv
#main_progen[:,4] = main_progen[:,4]*conv
#np.savetxt('Halo%d_starmerger_test_%s.out'%(hnum,date),main_progen)

#headers = 'id type x y z vx vy vz main r_form(kpc) atime(Gyr) ctime(Gyr) mass (M_sun)'.split()
#for line in fileinput.input(['Halo%d_starmerger_hist_11_3.out'%hnum], inplace=True):
#    if fileinput.isfirstline():
#        print '\t'.join(headers)
#main_progen = np.genfromtxt('Halo7new_starmerger_hist_11_12.out')
#plot_starmerge(main_progen,hnum,halox,haloy,halorad,xcen,ycen)
#plot_projden(main_progen,xcen,ycen,zcen,rvir, hnum)
#kristy_plots(main_progen,hnum,halox,haloy,halorad,xcen,ycen,zcen,rvir)
#meanstellarage_vs_r(main_progen,xcen,ycen,zcen,rvir, hnum)
