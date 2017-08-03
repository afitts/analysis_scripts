import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as co
from matplotlib import rcParams
from matplotlib import rc
import pylab
import time
import scipy.integrate
from collections import Counter
#import astropy as ap
#from astropy.cosmology import FlatLambdaCDM
#import astropy.units as u

rcParams['lines.linewidth'] = 2
rcParams['axes.linewidth'] = 3
rcParams['xtick.major.size'] = 10
rcParams['xtick.minor.size'] = 5
rcParams['xtick.major.width'] = 2
rcParams['xtick.minor.width'] = 2
rcParams['ytick.major.size'] = 10
rcParams['ytick.minor.size'] = 5
rcParams['ytick.major.width'] = 2
rcParams['ytick.minor.width'] = 2
rcParams['font.size'] = 18
rcParams['xtick.labelsize']= '14'
rcParams['ytick.labelsize']= '14'
rcParams['savefig.bbox'] = 'tight' 
rcParams['axes.labelsize'] = 16
rcParams['xtick.labelsize']= 16
rcParams['ytick.labelsize']= 16

#rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#rc('text', usetex=True)

def scinot(x,pos=None):
    if x == 0:
        s = '0'
    else:
        xp = int(np.floor(np.log10(np.abs(x))))
        mn = x/10.**xp
        # Here we truncate to 2 significant digits -- may not be enough 
        # in all cases
        s = '$'+str('%.1f'%mn) +'\\times 10^{'+str(xp)+'}$'
    return s
#ax.yaxis.set_major_formatter(tick.FuncFormatter(scinot)) how to use

def _a_dot(a, h0, om_m, om_l):
    om_k = 1.0 - om_m - om_l
    return h0 * a * np.sqrt(om_m * (a ** -3) + om_k * (a ** -2) + om_l)


def _a_dot_recip(*args):
    return 1. / _a_dot(*args)

def arch_SFH(pathname,hist,hnum,snum,res,ver,date,extent,colors):
  global j
  hnum1 = hnum
  time = np.zeros(len(hist))
  rvir = np.zeros(len(hist))
  rhalf = np.zeros(len(hist))
  red = np.zeros(len(hist))
  hdf = pd.HDFStore('%sanalysis/dataframes/halo%s%s_giz%s_snap%03d.h5'%(pathname,hnum,res,ver,snum))
  red = np.float(hdf['props']['redshift'])
  rhalf = np.float(hdf['props']['rhalf'])
  rvir = np.float(hdf['props']['rvir'])
  if extent == 0 :
    extent = rhalf
  elif extent == 1 :
    extent = 3*rhalf
  elif extent == 2 :
    extent = rvir
  sm = hdf['particles/star']['mass']
  dm_mass = sum(hdf['particles/dm']['mass'][hdf['particles/dm']['r']<extent])
  cg_mass = sum(hdf['particles/gas']['mass'][hdf['particles/gas']['r']<extent][hdf['particles/gas']['temp']<2e4])
  g_mass = sum(hdf['particles/gas']['mass'][hdf['particles/gas']['r']<extent])
  ct = hdf['particles/star']['sft'][hdf['particles/star']['r']<extent].as_matrix()
  #print ct, extent,len(ct)
  S_mas = hdf['particles/star']['mass'][hdf['particles/star']['r']<extent].as_matrix()
  s_id = hdf['particles/star'].index[hdf['particles/star']['r']<extent]

  YEAR = 60*60*24*365
  h0 = 71
  om_l = 0.734
  om_m = 0.266
  conv = 3.085677581e+19
  #cosmo = FlatLambdaCDM(H0=71,Om0=0.266)
  current_time = np.float(hdf['props']['time'])*YEAR*1e9 #4.331972446834427e+17
  hdf.close()
#begin dm mass within rhalf calc

  print "dm_rhalf is %.3e"%dm_mass

#begin sfr calculations
  C_age = np.zeros(len(ct))
  for i in np.arange(len(ct)):
    C_age[i] = scipy.integrate.quad(_a_dot_recip, 0, ct[i], (h0, om_m, om_l))[0]*conv
  #ct = (1/np.array(ct))-1
  #lookback = cosmo.lookback_time(ct)
  #print lookback
  #C_age = current_time-lookback*YEAR*1e9
  bin_count = 137#186

  # Find the oldest stars in units of code time.
  tmin = min(C_age)
  # Multiply the end to prevent numerical issues.
  time_bins = np.linspace(tmin*1.01,current_time, num = bin_count+1)#np.genfromtxt('Halo897_13_mstar.out')#
  #time_bins = time_bins[:,0]*1e9*YEAR
  # Figure out which bins the stars go into.
  inds = np.digitize(C_age, time_bins) - 1
  # Sum up the stars created in each time bin.
  mass_bins = np.zeros(bin_count + 1, dtype ='float64')
  for i,e in enumerate(inds):
    mass_bins[e] += S_mas[i]
  # Calculate the cumulative mass sum over time by forward adding.
  cum_mass_bins = mass_bins.copy()
  for index in xrange(bin_count):
    cum_mass_bins[index+1] += cum_mass_bins[index]
  # We will want the time taken between bins.
  time_bins_dt = time_bins[1:] - time_bins[:-1]

  tc = 1#ds["Time"] IS THIS CORRECT FOR ALL SNAPS OR JUST 184?
  time = []
  lookback_time = []
  redshift = []
  Msol_yr = np.zeros(len(time_bins)-1)
  Msol = []
  Msol_cumulative = []
# Use the center of the time_bin, not the left edge.
  for i, t in enumerate((time_bins[1:] + time_bins[:-1])/2.):
    time = np.append(time,t * tc / YEAR)
    lookback_time = np.append(lookback_time,(current_time - t * tc)/YEAR)
    redshift = np.append(redshift,1/tc-1)
    Msol_yr[i] = mass_bins[i] /(time_bins_dt[i] * tc /YEAR)
    Msol = np.append(Msol,mass_bins[i])
    Msol_cumulative = np.append(Msol_cumulative,cum_mass_bins[i])
    if (Msol_yr[i]==0):
      Msol_yr[i] = 1e-7

  w, h = plt.figaspect(.9)
  fig = plt.figure(figsize=(w,h))
  plt.ylim(1e-5,1e-2)
  #plt.xlim(0,0.8)
  plt.xlabel("Time (Gyr)")
  plt.ylabel(r"Star Formation Rate (M$_\odot$/yr)")
  plt.semilogy(time/1e9,Msol_yr, label =r'Y%s 100 Myr'%(res))
  #plt.legend()
  plt.savefig("Halo%s%s_giz%s_SFR_rhalf_z0_bin100Myr_%s.pdf"%(hnum1,res,ver,date), transparent=True)
  #plt.show()
  #plt.close()
  if hnum != '125961':
   if ver=='11_13':
    ax.plot(time/1e9,Msol_cumulative,'%s'%colors, label = 'Halo %s %s'%(hnum, ver))
    ax1.plot(time/1e9,Msol_cumulative/max(Msol_cumulative),color=smm.to_rgba(np.log10(max(Msol_cumulative))),linewidth=4,label = 'Halo %s %s'%(hnum, res))
    if hnum == '1016':    
      ax2.semilogy(time/1e9,Msol_yr, color='k',label ='Halo %s %s'%(hnum, ver))
    elif hnum == '2' or hnum == '848':
      ax3.semilogy(time/1e9,Msol_yr, color='k',label ='Halo %s %s'%(hnum, ver))  
   else:
    ax.plot(time/1e9,Msol_cumulative,'%s'%colors, label = 'Halo %s %s'%(hnum, ver))
    ax1.plot(time/1e9,Msol_cumulative/max(Msol_cumulative),color=smm.to_rgba(np.log10(max(Msol_cumulative))),linewidth=4,label = 'Halo %s %s'%(hnum, res))
    if hnum == '1016':    
      ax2.semilogy(time/1e9,Msol_yr, color='r',label ='Halo %s %s'%(hnum, ver))
    elif  hnum == '2'or hnum == '848':
      ax3.semilogy(time/1e9,Msol_yr, color='r',label ='Halo %s %s'%(hnum, ver))      
  plt.figure()
  plt.plot(time/1e9,Msol_cumulative)
  plt.xlabel("Time (Gyr)")
  plt.ylabel(r"Cumulative SFR ($10^6$ M$_\odot$)")
  plt.legend(loc=4)
  #print time_bins/1e9/YEAR, Msol_cumulative/max(Msol_cumulative)
  np.savetxt('Halo%s_13_sfr_rhalf.out'%hnum1,np.column_stack((time_bins[1:]/1e9/YEAR,Msol_yr)), fmt= '%.8e')  
  np.savetxt('Halo%s_13_cumusfr_rhalf.out'%hnum1,np.column_stack((time_bins[1:]/1e9/YEAR,Msol_cumulative)), fmt= '%.8e')
  plt.savefig("Halo%s_arch_SFH_rhalf_z0_bin100Myr_%s.pdf"%(hnum,date), transparent=True)
  #plt.show()
  plt.close()
  j+=1
  return True

class merger_SFH(object):
	def __init__(self):
		w, h = plt.figaspect(.9)
		self.fig = plt.figure(figsize=(w,h))
		self.sub = self.fig.add_subplot(111)
		self.sub.set_ylim(0,1.1)
		#plt.xlim(0,0.8)
		self.sub.set_xlabel("Time (Gyr)")
		self.sub.set_ylabel('Cumulative Stellar Mass Fraction')
	def add_line(self,hnum):
		hnum = '1016'
		stars = np.genfromtxt('starmerger_out/Halo%s_starmerger_test.out'%hnum)
		insitu = stars[stars[:,8]==1]
		merger = stars[stars[:,8]==0]
		time = np.linspace(min(stars[:,8]),13.75,137) 
		insitu_sfh,bin_edges = np.histogram(insitu[:,11], bins = time, weights=insitu[:,12])
		merger_sfh,bin_edges = np.histogram(merger[:,11], bins = time, weights=merger[:,12])
		insitu_sfh = np.cumsum(insitu_sfh)/sum(insitu_sfh)
		merger_sfh = np.cumsum(merger_sfh)/sum(merger_sfh)
		self.sub.plot(time[1:],insitu_sfh,color = 'k',linewidth = 3)
		self.sub.plot(time[1:],merger_sfh,color = 'r',linewidth = 3)
	def save(self,hnum,date):
		self.fig.savefig("Halo%s_arch_SFH_merger_v_insitu_z0_bin100Myr_%s.pdf"%(hnum,date), transparent=True)
class single_merger_SFH(object):
	def __init__(self):
		w, h = plt.figaspect(.9)
		self.fig = plt.figure(figsize=(w,h))
		self.sub = self.fig.add_subplot(111)
		#self.sub.set_ylim(0,1.1)
		#plt.xlim(0,0.8)
		self.sub.set_xlabel("Time (Gyr)")
		self.sub.set_ylabel('Cumulative Stellar Mass')
	def add_line(self,hnum):
		#hnum = '1016'
		stars = np.genfromtxt('starmerger_out/Halo%s_starmerger_test.out'%hnum)
		insitu = stars[stars[:,8]==1]
		merger = stars[stars[:,8]==0]
		time = np.linspace(min(stars[:,8]),13.75,137) 
		insitu_sfh,bin_edges = np.histogram(insitu[:,11], bins = time, weights=insitu[:,12])
		insitu_sfh = np.cumsum(insitu_sfh)
		self.sub.semilogy(time[1:],insitu_sfh,color = 'k',linewidth = 3)

		mylist = np.sort(merger[:,10])
		binx = np.sort([k for k,v in Counter(mylist).items() if v>1])
		binx = binx+0.001
		for i in np.arange(len(binx)):
			merger_sfh,bin_edges = np.histogram(merger[np.digitize(merger[:,10],binx)==i,11], bins = time, weights=merger[np.digitize(merger[:,10],binx)==i,12])
			merger_sfh = np.cumsum(merger_sfh)
			self.sub.semilogy(time[1:],merger_sfh,color = 'r',linewidth = 3)

	def save(self,hnum,date):
		self.fig.savefig("Halo%s_arch_SFH_merger_v_insitu_z0_bin100Myr_%s.pdf"%(hnum,date), transparent=True)
class merger_SFH(object):
	def __init__(self):
		w, h = plt.figaspect(.9)
		self.fig = plt.figure(figsize=(w,h))
		self.sub = self.fig.add_subplot(111)
		#self.sub.set_ylim(0,1.1)
		self.sub.set_ylim(1e4,5e6)
		#plt.xlim(0,0.8)
		self.sub.set_xlabel("Time (Gyr)")
		self.sub.set_ylabel('Cumulative Stellar Mass')
	def add_line(self,hnum):
		if hnum == '10166':
			stars = np.genfromtxt('Halo1016rerun_13_cumusfr_rhalf.out')
			clr = 'r'
			hnum = '1016NEW'
		elif hnum == '1016og':
			stars = np.genfromtxt('Halo1016_13_cumusfr_rhalf.out')
			clr = 'k'
			hnum = '1016'
		elif hnum == '848td':
			stars = np.genfromtxt('Halo848td_13_cumusfr_rhalf.out')
			clr = 'r'
		elif hnum == '848':
			stars = np.genfromtxt('Halo848_13_cumusfr_rhalf.out')
			clr = 'k'
		elif hnum == '11707td':
			stars = np.genfromtxt('Halo11707td_13_cumusfr_rhalf.out')
			clr = 'r'
		elif hnum == '11707':
			stars = np.genfromtxt('Halo11707_13_cumusfr_rhalf.out')
		elif hnum == '12596':
			stars = np.genfromtxt('Halo12596_13_cumusfr_rhalf.out')
			clr = 'r'
			hnum = '12596NEW'
		elif hnum == '12596OG':
			stars = np.genfromtxt('Halo12596OG_13_cumusfr_rhalf.out')
			clr = 'k'
			hnum = '12596'
		else:
			stars = np.genfromtxt('Halo%s_13_cumusfr_rhalf.out'%hnum)
		time = stars[:,0]
		self.sub.semilogy(time,stars[:,1],linewidth = 3,label=hnum)
	def save(self,hnum,date):
		self.sub.legend(loc = 4)
		self.fig.savefig("Halo1016%s_arch_SFH_z0_bin100Myr_%s.pdf"%(hnum,date), transparent=True)

class single_merger_SFH1(object):
	def __init__(self):
		w, h = plt.figaspect(.9)
		self.fig = plt.figure(figsize=(w,h))
		self.sub = self.fig.add_subplot(111)
		#self.sub.set_ylim(0,1.1)
		#plt.xlim(0,0.8)
		self.sub.set_xlabel("Time (Gyr)")
		self.sub.set_ylabel('Cumulative Stellar Mass')
	def add_line(self,hnum):
		if hnum == '1016':
			stars = np.genfromtxt('Halo1016rerun_13_cumusfr_rhalf.out')
		elif hnum == '1016og':
			stars = np.genfromtxt('Halo1016_13_cumusfr_rhalf.out')
		elif hnum == '848td':
			stars = np.genfromtxt('Halo848td_13_cumusfr_rhalf.out')
		time = stars[:,0]
		self.sub.semilogy(time,stars[:,1],linewidth = 3)

	def save(self,hnum,date):
		self.fig.savefig("Halo%s_arch_SFH_z0_bin100Myr_%s.pdf"%(hnum,date), transparent=True)

date = time.strftime("%m_%d_%Y")
my_cmap=plt.get_cmap('plasma')
smm = cm.ScalarMappable(cmap=my_cmap,norm=co.Normalize(vmin=5.6739, vmax=7.176))#vmin=5.874853, vmax=7.11806))
smm._A = []
fig = plt.figure()
ax = fig.add_subplot(111)
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
cb1 = fig1.colorbar(smm)
cb1.set_label(r'$\mathrm{log M}_*$')
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
colors = ['r','r','k','k','b','b','g','g','m','m','y','y','olive','darkviolet']#['olive','darkviolet','lime','b','gray','r','b','k','g','y','m']

j = 0
hnum = '2'
res = '_13'
ver = '5_12_16'
hnum = ['1016','1016REDO','1016SNe']#'11707','32503','12596','007','848','796','20192','20910','897','948','1016']#,'32257','10q','10v'] #dmo should be placed after hydro run for irhalf variable.
res = ['_13','_13','_13','_11','_13','_13','_13','_13','_13','_13','_13','_13','_11','_11']
ver = ['5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16']#['11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13']
#hnum = ['2','948','796','897','1016']# No 007 for now as the dataframe from 184 apparently has 4 star particles under the same index for some reason
#res = ['_13','_13','_13','_13','_13']
#ver = ['11_13','11_13','11_13','11_13','11_13']
snum = 184
extent = 2
SFH_plot = merger_SFH()
for w in np.arange(len(hnum)):
	#if hnum[w] == '1016':
	pathname = '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/'%(hnum[w],res[w],ver[w])
	hist = np.genfromtxt(pathname+'analysis/halo%smerger_hist.txt'%hnum[w])
	arch_SFH(pathname,hist,hnum[w],snum,res[w],ver[w],date,extent,colors[j])
	single_SFH_plot = single_merger_SFH1()
	#single_SFH_plot.add_line(hnum[w])
	#single_SFH_plot.save(hnum[w],date)
	#if w == 0:
	#	hnum[w] = '12596OG'
	SFH_plot.add_line(hnum[w])
#  pathname = '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/'%(hnum[w],res[w],ver[w])

#  hist = np.genfromtxt(pathname+'analysis/halo%smerger_hist.txt'%hnum[w])
#  if hnum == '10q' or hnum == '10v':
#   snum = 600
#  arch_SFH(pathname,hist,hnum[w],snum,res[w],ver[w],date,extent,colors[j])
SFH_plot.save('all',date)
ax.set_ylabel(r'$\mathrm{M}_\star$ $\mathrm{M}_\odot$')
ax.set_xlabel('Time (Gyr)')
ax.legend(loc=2,prop={'size':8})
#fig.savefig('mstar_v_t_%s.pdf'%(date))
ax1.set_ylabel('Cumulative Stellar Mass Fraction')
ax1.set_xlabel('Time (Gyr)')
ax1.set_ylim(-.05,1.05)
#ax1.legend(loc=4,prop={'size':10})
fig1.savefig('arch_cumu_sfr_%s.pdf'%(date))
ax2.set_ylim(1e-5,1e-2)
ax2.set_xlabel("Time (Gyr)")
ax2.set_ylabel(r"Star Formation Rate (M$_\odot$/yr)")
ax2.legend(loc=1,prop={'size':10})
#fig2.savefig('arch_sfr_rhalf_halo1016comp_%s.pdf'%(date))
ax3.set_ylim(1e-5,1e-2)
ax3.set_xlabel("Time (Gyr)")
ax3.set_ylabel(r"Star Formation Rate (M$_\odot$/yr)")
ax3.legend(loc=1,prop={'size':10})
#fig3.savefig('arch_sfr_rhalf_%s.pdf'%(date))
#plt.show()

