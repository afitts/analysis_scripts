import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import FixedLocator
from matplotlib.ticker import FixedFormatter
from matplotlib import gridspec
import pylab
import time
import sys

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
#rcParams['figure.figsize'] = [12,9]

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
  mstarhalf = np.zeros(len(hist))
  for i in np.arange(len(hist)):
    print i
    try:
      hdf = pd.HDFStore('%sanalysis/dataframes/halo%s%s_giz%s_snap%03d.h5'%(pathname,hnum,res,ver,i))
      red[i] = np.float(hdf['props']['redshift'])
  #    gcum = np.cumsum(hdf['particles/gas'].sort('r')['mass'])
  #    isix = abs(gcum - 1e6).argmin()
  #    iseven = abs(gcum - 1e7).argmin()
  #    ieight = abs(gcum - 1e8).argmin()
  #    sixr[i] =  np.float(hdf['particles/gas']['r'][isix])*1000/(red[i]+1)
  #    if (max(gcum) <= 1e6):
	#sixr[i] = np.nan
      #if (max(gcum) >= 1e7):
      #  sevenr[i] =  np.float(hdf['particles/gas']['r'][iseven])*1000/(red[i]+1)
      #else:
#	sevenr[i] = np.nan
#      if (max(gcum) >= 1e8):
#        eightr[i] =  np.float(hdf['particles/gas']['r'][ieight])*1000/(red[i]+1)
#      else:
#	eightr[i] = np.nan
      time[i] = np.float(hdf['props']['time'])
      rvir[i] = np.float(hdf['props']['rvir'])*1000/(red[i]+1)
      rhalf[i] = np.float(hdf['props']['rhalf'])*1000/(red[i]+1)
      mu = 0.59 # mean molecular weight for a fully ionized gas w/ primordial composition
      mp = 1.67e-24 #mass of proton
      G = 6.67e-8 #Grav
      k = 1.38e-16 #Boltzmann
      #vvir = np.sqrt(G*sum(hdf['particles/dm']['mass']*1.99e33)/(rvir[i]*3.09e21))
      #mstarhalf[i] = rhalf[i]#0.5 * mu * mp *vvir**2/k <--virial temp
      hdf.close()
    except Exception,f:
      print f
      hdf.close()
  bob = np.loadtxt('halo12596_og_extent.txt')
  sixr = bob[:,0]
  sevenr = bob[:,1]
  eightr = bob[:,2]
  fig = plt.figure(figsize=(8, 6)) #, ax1 = plt.subplots()
  gs = gridspec.GridSpec(2, 1, height_ratios=[1, 3],hspace=0)
  ax2 = plt.subplot(gs[0])
  ax1 = plt.subplot(gs[1])
  ax1.semilogy(time,sixr, 'k', label = r'$10^6 M_\odot$',zorder=100)
  vers = ['a']#,'b','c','d','e','f','g','h','i','j']
  for k in np.arange(len(vers)):
    redo = np.genfromtxt('halo12596%s_extended_extent.txt'%vers[k])
    ax1.semilogy(time[121:161],redo[:,0], 'k',linewidth=1,zorder=100)
    #ax1.semilogy(time[121:129],redo[:,1], '0.6',linewidth=1,zorder=100)
    #ax1.semilogy(time[121:129],redo[:,2], 'r',linewidth=1,zorder=100)
  redo = np.genfromtxt('halo12596_redo_extent.txt')
  ax1.semilogy(time[120:150],redo[:,0], 'k',linewidth=1,zorder=100)
  #ax1.semilogy(time,sevenr, '0.6', label = r'$10^7 M_\odot$',zorder=100)
  #ax1.semilogy(time[120:150],redo[:,1], '0.6',linewidth=1,zorder=100)
  #eightr = np.genfromtxt('halo12596_1e8_extent.txt')
  #ax1.semilogy(time,eightr, 'r', label = r'$10^8 M_\odot$',zorder=100)
  #ax1.semilogy(time[120:150],redo[:,2], 'r',linewidth=1,zorder=100)
  ax1.semilogy(time,rvir, 'k--', label = r'$R_{vir}$',zorder=100)
  ax1.semilogy(time,rhalf,'m--', label = r'$R_{1/2}$',linewidth = 2,zorder=100)
  np.savetxt('1e6_gas_shell.txt',np.column_stack((time,sixr)))
  plt.xlabel('Time (Gyr)')
  plt.ylabel('Radial Extenet (kpc)')
  plt.xlim(8,max(time))#min(time[time != 0]),max(time))
  plt.ylim(.1,2e2)
  l = plt.legend(loc = 2,prop={'size':10},frameon=False)
  l.set_zorder(20)
  #ax2 = ax1.twinx()
  ax2.axes.get_xaxis().set_ticklabels([])
  ax2.set_ylim(-5.0,np.log10(2e-2))
  ax2.set_xlim(8,max(time))
  ticks = np.linspace(-5.,-2.,4)
  plt.axes(ax2)
  plt.yticks(ticks,[r'$10^{-5}$',r'$10^{-4}$',r'$10^{-3}$',r'$10^{-2}$'])#[r'$10^1$',r'$10^2$'],fontsize=10) #[r'$10^6$',r'$10^7$',r'$10^8$',r'$10^9$'],fontsize=16)#
  minor_ticks=[]
  for j in range(2,10):
    minor_ticks.append(-5+np.log10(j))
  for j in range(2,10):
    minor_ticks.append(-4+np.log10(j))
  for j in range(2,10):
    minor_ticks.append(-3+np.log10(j))
  ax2.yaxis.set_minor_locator(FixedLocator(minor_ticks))
  ax2.set_ylabel(r"$\dot{M}_\star\:\:\:$(M$_\odot$/yr)")
  sfrbin = 5
  a = np.genfromtxt('sfr%s0_out/Halo%s%s_sfr_rvir.out'%(sfrbin,hnum,res))
  sout = np.genfromtxt('mvir_out/Halo%s%s_5_12_16_mstar.out'%(hnum,res))
  ax2.plot(a[:,0],np.log10(a[:,1]), color='k',linewidth=2,zorder=1)


  #ax1.set_zorder(ax2.get_zorder()+1) # put ax in front of ax2 
  #ax1.patch.set_visible(False) # hide the 'canvas' 

  #ax2.plot(time,mstarhalf, 'm--',marker ='o', label = r'$T(<r_{vir})$',linewidth=4)
  #ax2.set_ylim(0,4e4)
  plt.savefig('gas_shells_semilogy_1e6only_halo%s%s_ver%s_%s.pdf'%(hnum,res,ver,date),transparent = True)
  plt.show()
  return True

w = int(sys.argv[1])
date = time.strftime("%m_%d_%Y")
hnum =['12596OG']#11707','32503','12596','007','848','796','20192','20910','897','948','1016','32257','10q','10v','1084']
res = ['_13','_13','_13','_11','_13','_13','_13','_13','_13','_13','_13','_13','_11','_11','_13']
#res = ['','','','','','','','','','','','','_11','_11','']
ver = ['5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16']


pathname = '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/'%(hnum[w],res[w],ver[w])
hnum[w] = '12596'
hist = np.genfromtxt(pathname+'analysis/halo%smerger_hist.txt'%hnum[w])

gas_shell_plot(pathname,hist,hnum[w],res[w],ver[w],date)

