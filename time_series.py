import numpy as np
import sys 
import glob
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pylab
import time
import pandas as pd

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

littleh = 0.71
omegam = .266
rhocrit = 277.5

def add_vt(hnum,res,ver,extent,dmo):
  global count
  pathname = '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/'%(hnum,res,ver)
  hist = np.genfromtxt(pathname+'analysis/halo%smerger_hist.txt'%hnum)
  dm_mass = np.zeros(len(hist))
  g_mass = np.zeros(len(hist))
  s_mass = np.zeros(len(hist))
  time = np.zeros(len(hist))
  vmax = np.zeros(len(hist))
  GYR = 60*60*24*365*1e9
  for i in np.arange(len(hist)):
    try:
      hdf = pd.HDFStore('%sanalysis/dataframes/halo%s%s_giz%s_snap%03d.h5'%(pathname,hnum,res,ver,i))
      red = float(hdf['props']['redshift'])
      rhalf = float(hdf['props']['rhalf'])/(red+1)
      rvir = float(hdf['props']['rvir'])/(red+1)
      if extent == 0 :
        ex = rhalf
      elif extent == 1 :
        ex = 3*rhalf
      elif extent == 2 :
        ex = rvir
      print ex
      dm_mass[i] = sum(hdf['particles/dm']['mass'][hdf['particles/dm']['r']<ex])
      time[i] = hdf['props']['time']
      vmax[i] = hdf['props']['vmax'] 
      if dmo == 0:
        g_mass[i] = sum(hdf['particles/gas']['mass'][hdf['particles/gas']['r']<ex])
        s_mass[i] = sum(hdf['particles/star']['mass'][hdf['particles/star']['r']<ex])
    except Exception,f:
      f = str(f)
      f.replace("'","")
      print f
      if f != "'No object named props in the file'" and f != "'No object named particles/dm in the file'":
        dm_mass[i] = sum(hdf['particles/dm']['mass'][hdf['particles/dm']['r']<ex])
        time[i] = hdf['props']['time']
        vmax[i] = hdf['props']['vmax']
  print s_mass
  if dmo == 0:
    mvirvt.loglog(time,1.42*dm_mass, '%s'%clr[count],linewidth=2,label='%s tot'%(hnum))
    #ax.loglog(time/GYR,1.42*(prop[:,1]-prop[:,2]-prop[:,11]), '%s'%clr[count+2],linewidth=2,label='%s dm'%(hnum))
    #ax.loglog(prop[:,0]/GYR,1.42*(prop[:,11]), '%s--'%clr[count+2],linewidth=2,label='%s gas'%(hnum))
    #ax.loglog(prop[:,0]/GYR,1.42*(prop[:,2]), '%s-.'%clr[count+2],linewidth=2,label='%s star'%(hnum))
    #ax1.plot(prop[:,0]/GYR,1.42*prop[:,2]/max(prop[:,2])/1.42, '%s'%clr[count],label='%s'%(hnum))
    mstarvt.plot(time,1.42*s_mass, '%s'%clr[count],label='%s'%(hnum))
    vmaxvt.plot(time,vmax, '%s'%clr[count],label='%s'%(hnum))
    bfracvt.semilogy(time,(g_mass+s_mass)/dm_mass, '%s'%clr[count],label='%s'%(hnum))
    bfracvt.axhline(y=1-0.83120300751)
    #ax4.plot(time/GYR,prop[:,3]/dm_mass[:,2], '%s'%clr[count],label='%s'%(hnum))
  else:
    mvirvt.loglog(time,1.42*0.83120300751*dm_mass, '%s--'%clr[count],linewidth=2,label='%s dmo (corrected)'%(hnum))
    #ax.loglog(time/GYR,1.42*prop[:,1], '%s:'%clr[count-2],linewidth=2,label='%s dmo'%(hnum))
    vmaxvt.plot(time,vmax, '%s--'%clr[count],label='%s'%(hnum))
  count +=1
  return True

fig = plt.figure(figsize=(10*1.5,10))
fig1 = plt.figure()
fig2 = plt.figure()
fig3 = plt.figure()
fig4 = plt.figure()
mvirvt = fig.add_subplot(111)
mstarvt = fig1.add_subplot(111)
vmaxvt = fig2.add_subplot(111)
bfracvt = fig3.add_subplot(111)
ax4 = fig4.add_subplot(111)

mvirvt.set_ylabel("$M_{tot}$($M_\odot$)")
mvirvt.set_xlabel("Time (Gyr)")
mvirvt.set_xlim(5e-1,13.7)
mvirvt.set_ylim(1e7,2e10)

mstarvt.set_ylabel("$M_*$($M_\odot$)")
mstarvt.set_xlabel("Time (Gyr)")

vmaxvt.set_ylabel("$V_{max}$($km/s$)")
vmaxvt.set_xlabel("Time (Gyr)")

bfracvt.set_ylabel("$f_{bar}$")
bfracvt.set_xlabel("Time (Gyr)")
bfracvt.set_ylim(1e-3,.5)

ax4.set_ylabel("$(Hydro V_{max})/(DM V_{max})$")
ax4.set_xlabel("Time (Gyr)")
ax4.set_title("Halos Vmax v t")
clr = ['r','b','k','g','y','m','c','r','b','k','g','y','m']

extent = 2
hnum = '007'#['007','2','948','796','897','1016','007','2','948','796','897','1016']
res = '_11'#['_11','_13','_13','_13','_13','_13','_11','_13','_13','_13','_13','_13']
ver = '11_13'#['11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13']
dm = 0#[0,0,0,0,0,0,1,1,1,1,1,]
count = 0
#for w in np.arange(len(hnum)):
  #add_vt(hnum[w],res[w],ver[w],extent,dm[w])
  
add_vt(hnum,res,ver,extent,dm)
#hnum = '1084'
#res = '_13'
#ver = '11_13'
#fname = '%s/mfm%s%s_giz%s_raw_output/analysis/Halo%s%s%s_props.out'%(pathname,hnum,res,ver,hnum,res,ver)
#fname1 = '%s/gizdm%s%s_raw_output/analysis/Halo_dm%s_%s_props.out'%(pathname,hnum,res,hnum,res)
#add_vt(hnum,res,fname,fname1,0)

mvirvt.legend(loc=4,prop={'size':10})
mstarvt.legend(loc=2,prop={'size':10})
vmaxvt.legend(loc=4, prop={'size':10})
bfracvt.legend(loc=1)
ax4.legend(loc=4)
date = time.strftime("%m_%d_%Y")
fig.savefig('mvir_v_t_%s.pdf'%date,transparent=True)
fig1.savefig('mstar_v_t_linear_%s.pdf'%date,transparent=True)
fig2.savefig('vmax_v_t_%s.pdf'%date,transparent=True)
fig3.savefig('mbarymvir_ratio_all_%s.pdf'%date,transparent=True)
fig4.savefig('vcirc_ratio_all_1_11.pdf',transparent=True)
plt.show('all')
plt.close('all')

