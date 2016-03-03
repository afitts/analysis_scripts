import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pylab
import time

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

def gas_shell_plot(pathname,hist,hnum,res,ver,date):
  sixr = np.zeros(len(hist))
  sevenr = np.zeros(len(hist))
  eightr = np.zeros(len(hist))
  time = np.zeros(len(hist))
  rvir = np.zeros(len(hist))
  rhalf = np.zeros(len(hist))
  red = np.zeros(len(hist))
  for i in np.arange(len(hist)):
    print i
    try:
      hdf = pd.HDFStore('%sanalysis/dataframes/halo%s%s_giz%s_snap%03d.h5'%(pathname,hnum,res,ver,i))
      red[i] = np.float(hdf['props']['redshift'])
      gcum = np.cumsum(hdf['particles/gas'].sort('r')['mass'])
      isix = abs(gcum - 1e6).argmin()
      iseven = abs(gcum - 1e7).argmin()
      ieight = abs(gcum - 1e8).argmin()
      sixr[i] =  np.float(hdf['particles/gas']['r'][isix])*1000/(red[i]+1)
      if (max(gcum) <= 1e6):
	sixr[i] = np.nan
      if (max(gcum) >= 1e7):
        sevenr[i] =  np.float(hdf['particles/gas']['r'][iseven])*1000/(red[i]+1)
      else:
	sevenr[i] = np.nan
      if (max(gcum) >= 1e8):
        eightr[i] =  np.float(hdf['particles/gas']['r'][ieight])*1000/(red[i]+1)
      else:
	eightr[i] = np.nan
      time[i] = np.float(hdf['props']['time'])
      rvir[i] = np.float(hdf['props']['rvir'])*1000/(red[i]+1)
      rhalf[i] = np.float(hdf['props']['rhalf'])*1000/(red[i]+1)
      hdf.close()
    except Exception,f:
      print f
      hdf.close()
  plt.figure()
  plt.semilogy(time,sixr, 'k', label = r'$10^6 M_\odot$')
  plt.semilogy(time,sevenr, 'b', label = r'$10^7 M_\odot$')
  plt.semilogy(time,eightr, 'r', label = r'$10^8 M_\odot$')
  plt.semilogy(time,rvir, 'k--', label = r'$R_{vir}$')
  plt.semilogy(time,rhalf,'g--', label = r'$R_{1/2}$',linewidth = 1)
  plt.xlabel('Time (Gyr)')
  plt.ylabel('Radial Extenet (kpc)')
  plt.xlim(min(time[time != 0]),max(time))
  plt.ylim(.1,60)
  plt.legend(loc = 2,prop={'size':14})
  plt.savefig('%sanalysis/gas_shells_semilogy_halo%s%s_ver%s_%s.pdf'%(pathname,hnum,res,ver,date),transparent = True)
  plt.show()
  return True

date = time.strftime("%m_%d_%Y")
hnum = '897'
res = '_13'
ver = '11_13'
pathname = '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/'%(hnum,res,ver)

hist = np.genfromtxt(pathname+'analysis/halo%smerger_hist.txt'%hnum)

gas_shell_plot(pathname,hist,hnum,res,ver,date)

