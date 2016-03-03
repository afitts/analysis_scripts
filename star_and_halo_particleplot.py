import numpy as np
import matplotlib 
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

def get_starsnhalos(pathname,hnum,res,ver,snap_ind,xcen,ycen,zcen,rvir): #Simple function to get 2D scatter plots of star particles with circles to denote rvir of halos at z=0
  global conv
  for i in np.arange(8):
    fname = '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/snapdir_184/snapshot_184.%d.hdf5'%(hnum,res,ver,i)
    a = h5py.File(fname,'r')
    if i == 0:
      pos = a['PartType4']['Coordinates']
      masstot = a['PartType4']['Masses']
    else:
      pos = np.vstack((pos,a['PartType4']['Coordinates']))
      masstot = np.concatenate((masstot,a['PartType4']['Masses']),0) 
  xtot = pos[:,0]/h* conv
  ytot = pos[:,1]/h* conv
  ztot = pos[:,2]/h* conv
  for w in np.arange(16):
    temph = glob.glob(pathname+'analysis/ahf_snap%03d/ahf.snap_%03d.%04d.z*.*.AHF_halos'%(snap_ind,snap_ind,w))
    temph = str(temph).strip('[]').replace("'","")
    halo = np.loadtxt(temph,usecols=(0,3,5,6,7,11,64,63))
    for k in np.arange(len(halo)):
      halodiff = np.sqrt((xcen-halo[k,2]/h)**2+(ycen-halo[k,3]/h)**2+(zcen-halo[k,4]/h)**2)
      try:
       if halo[k,5] > 5 and halo[k,7]>10 and halodiff>20 or halodiff==0:# and halodiff<rvir or halodiff==0:
        ans = np.vstack((ans,(halo[k,1],halo[k,6])))
        x = np.append(x,halo[k,2]/h)
        y = np.append(y,halo[k,3]/h)
        z = np.append(z,halo[k,4]/h)
        rad = np.append(rad,halo[k,5]/h)
        bobob += halo[k,6]
      except:
       if halo[k,5] > 5 and halo[k,7]>10 and halodiff>20 or halodiff==0:#and halodiff<rvir or halodiff ==0:
        ans = np.array([[halo[k,1],halo[k,6]]])
        x = np.array([halo[k,2]])/h
        y = np.array([halo[k,3]])/h
        z = np.array([halo[k,4]])/h
        rad = np.array([halo[k,5]])/h
        bobob = halo[k,6]

  return xtot,ytot,ztot,x,y,z,rad 

def my_circle_scatter(axes, x_array, y_array, radius, **kwargs):
    for x, y, rad in zip(x_array, y_array,radius):
        circle = pylab.Circle((x,y), radius=rad, **kwargs)
        axes.add_patch(circle)
    return True

def plot_star_and_halo(starx,stary,starz,halox,haloy,haloz,halorad,xcen,ycen,zcen,hnum,res): #Plotting function for 2d scatter plots at z=0.
  global date
  plt.figure()
  axes=pylab.axes()
  my_circle_scatter(axes, halox-xcen, haloy-ycen, radius=halorad, fill=False,color = 'b')
  pylab.axis('scaled')
  plt.scatter(starx[abs(starz-zcen)<100]-xcen,stary[abs(starz-zcen)<100]-ycen,marker = '.',color='r')
  plt.xlim(-100,100)
  plt.ylim(-100,100)
  plt.xlabel('X Position (kpc)')
  plt.ylabel('Y Position (kpc)')
  plt.ticklabel_format(useOffset=False)
  plt.title("Halo %s XY Stellar Particle + Halo Positions at z=0"%(hnum))
  plt.savefig('halo%s%s_xy_starnhalo_z=0_satellite20_%s.pdf'%(hnum,res,date),transparent=True)
  plt.show()
  plt.close()

  plt.figure()
  axes=pylab.axes()
  my_circle_scatter(axes, haloy-ycen, haloz-zcen, radius=halorad, fill=False,color = 'b')
  pylab.axis('scaled')
  plt.scatter(stary[abs(starx-xcen)<100]-ycen,starz[abs(starx-xcen)<100]-zcen,marker = '.',color='r')
  plt.xlim(-100,100)
  plt.ylim(-100,100)
  plt.xlabel('Y Position (kpc)')
  plt.ylabel('Z Position (kpc)')
  plt.ticklabel_format(useOffset=False)
  plt.title("Halo %s YZ Stellar Particle + Halo Positions at z=0"%hnum)
  plt.savefig('halo%s%s_yz_starnhalo_z=0_satellite20_%s.pdf'%(hnum,res,date),transparent=True)
  plt.show()
  plt.close()

  plt.figure()
  axes=pylab.axes()
  my_circle_scatter(axes, halox-xcen, haloz-zcen, radius=halorad, fill=False,color = 'b')
  pylab.axis('scaled')
  plt.scatter(starx[abs(stary-ycen)<100]-xcen,starz[abs(stary-ycen)<100]-zcen,marker = '.',color='r')
  plt.xlim(-100,100)
  plt.ylim(-100,100)
  plt.xlabel('X Position (kpc)')
  plt.ylabel('Z Position (kpc)')
  plt.ticklabel_format(useOffset=False)
  plt.title("Halo %s XZ Stellar Particle + Halo Positions at z=0"%hnum)
  plt.savefig('halo%s%s_xz_starnhalo_z=0_satellite20_%s.pdf'%(hnum,res,date),transparent=True)
  plt.show()
  plt.close()

  return True

def gather_progen(pathname,snap,hnum): #Gather progenitors of halo from mtree file. Requires Mergertree and my_merger_hist.py to have been already run (for merger_hist.txt)

  mainid = np.genfromtxt(pathname+'analysis/halo%smerger_hist.txt'%hnum)
  f = np.genfromtxt(pathname+'analysis/merger/halo%s_snap%03d_snap%03d_mtree'%(hnum,snap,snap-1))
  mainid = mainid[184-snum,1]
  try:
    maindex = np.where(f[:,0] == mainid) #Find our halo
    maindex = int(maindex[0])
    numprogen = f[maindex,2]
    haloid = f[maindex,0]
    progen = f[maindex+1:(maindex+1)+(numprogen+1),1] #Grab progenitors from table beneath Halo id
    mainprogenid = progen[0]
  except:
    haloid = 'skip'
    mainprogenid = 0
    progen = np.zeros(9)

  return haloid,mainprogenid,progen

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

PI = 3.14159265359
h = 0.71
hnum = '2'
res = '_13'
ver = '11_13'
snum = 10
date = '12_18'
conv = 1000
snap = 184
hage = 13.424
pathname =  '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/'%(hnum,res,ver)
ahf_fname = '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/analysis/ahf_snap184/ahf.snap_184.00%02d.z0.000.AHF_halos'%(hnum,res,ver,snum)
ahfdata = readahf(ahf_fname)
ahfdata = ahfdata[0]
xcen,ycen,zcen,rvir = get_cen_rvir(ahfdata)
starx,stary,starz,halox,haloy,haloz,halorad = get_starsnhalos(pathname,hnum,res,ver,snap,xcen,ycen,zcen,rvir)
plot_star_and_halo(starx,stary,starz,halox,haloy,haloz,halorad,xcen,ycen,zcen,hnum,res)

switch = 0 #switch that determines if a main progenitor with stars has been found. When it is first found, all star particles will be denoted as from the main progenitor. From then on switch =1 and the complied list of particles will be carried through each snapshot.
main_progen = np.zeros(1)
fnum = 16

#for i,e in enumerate(np.array([182,183])):#np.arange(snap+1)):
#  if hnum == 7:
#    mainid = np.genfromtxt(pathname+'analysis/halo%03dmerger_hist.txt'%hnum)
#  else:
#    mainid = np.genfromtxt(pathname+'analysis/halo%dmerger_hist.txt'%hnum)
#  mainid = mainid[184-e,1]
#  if mainid != 0:
#    main_progen = gather_particles(main_progen,pathname,e,hnum,fnum)
#main_progen[:,0] = main_progen[:,0].astype('int')
#main_progen[:,2] = main_progen[:,2]*conv
#main_progen[:,3] = main_progen[:,3]*conv
#main_progen[:,4] = main_progen[:,4]*conv
#np.savetxt('Halo%dnew_starmerger_test_11_16.out'%hnum,main_progen)

#headers = 'id type x y z vx vy vz main r_form(kpc) atime(Gyr) ctime(Gyr) mass (M_sun)'.split()
#for line in fileinput.input(['Halo%d_starmerger_hist_11_3.out'%hnum], inplace=True):
#    if fileinput.isfirstline():
#        print '\t'.join(headers)
#main_progen = np.genfromtxt('Halo7new_starmerger_hist_11_12.out')
#plot_starmerge(main_progen,hnum,halox,haloy,halorad,xcen,ycen)
#plot_projden(main_progen,xcen,ycen,zcen,rvir, hnum)
#kristy_plots(main_progen,hnum,halox,haloy,halorad,xcen,ycen,zcen,rvir)
#meanstellarage_vs_r(main_progen,xcen,ycen,zcen,rvir, hnum)
