import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pylab
import time
import scipy.integrate
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
rcParams['figure.figsize'] = [12,9]

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
  ct = hdf['particles/star']['sft'][hdf['particles/star']['r']<extent]
  S_mas = hdf['particles/star']['mass'][hdf['particles/star']['r']<extent]
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
      C_age[i] = scipy.integrate.quad(_a_dot_recip, 0, ct[ct.index[i]], (h0, om_m, om_l))[0]*conv
  #ct = (1/np.array(ct))-1
  #lookback = cosmo.lookback_time(ct)
  #print lookback
  #C_age = current_time-lookback*YEAR*1e9
  bin_count = 137

  # Find the oldest stars in units of code time.
  tmin = min(C_age)
  # Multiply the end to prevent numerical issues.
  time_bins = np.linspace(tmin*1.01,current_time, num = bin_count+1)
  # Figure out which bins the stars go into.
  inds = np.digitize(C_age, time_bins) - 1
  # Sum up the stars created in each time bin.
  mass_bins = np.zeros(bin_count + 1, dtype ='float64')
  for i,e in enumerate(inds):
    mass_bins[e] += S_mas[S_mas.index[i]]
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

  plt.figure()
  plt.ylim(1e-5,1e-2)
  #plt.xlim(0,0.8)
  plt.xlabel("Time (Gyr)")
  plt.ylabel(r"Star Formation Rate (M$_\odot$/yr)")
  plt.semilogy(time/1e9,Msol_yr, label =r'Y%s 100 Myr'%(res))
  plt.legend()
  plt.savefig("Halo%s%s_giz%s_SFR_rhalf_z%s_%s.pdf"%(hnum,res,ver,red,date), transparent=True)
  #plt.show()
  #plt.close()
  ax.plot(time/1e9,Msol_cumulative,'%s'%colors, label = 'Halo %s %s'%(hnum, ver))
  ax1.plot(time/1e9,Msol_cumulative/max(Msol_cumulative),'%s'%colors,label = 'Halo %s %s'%(hnum, ver))
  ax2.semilogy(time/1e9,Msol_yr, '%s'%colors,label ='Halo %s %s'%(hnum, ver))
  plt.figure()
  plt.plot(time/1e9,Msol_cumulative)
  plt.xlabel("Time (Gyr)")
  plt.ylabel(r"Cumulative SFR ($10^6$ M$_\odot$)")
  plt.legend(loc=4)
  plt.savefig("Halo%s%s_giz%s_cum_SFR_rhalf_z%f_%s.pdf"%(hnum,res,ver,red,date), transparent=True)
  plt.show()
  #plt.close()
  j+=1
  return True

date = time.strftime("%m_%d_%Y")
fig = plt.figure()
ax = fig.add_subplot(111)
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
colors = ['r','b','g','c','m','y','y']
j = 0
hnum = '2'
res = '_13'
ver = '11_13'
snum = 184
extent = 1
pathname = '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/'%(hnum,res,ver)

hist = np.genfromtxt(pathname+'analysis/halo%smerger_hist.txt'%hnum)

arch_SFH(pathname,hist,hnum,snum,res,ver,date,extent,colors[j])

