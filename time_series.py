import numpy as np
import sys 
import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.cm as cm
import matplotlib.colors as co
import pylab
import time
import pandas as pd
from pyd import mikecm
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

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
rcParams['font.size'] = 16
rcParams['savefig.bbox'] = 'tight'


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

def add_vt(hnum,res,ver,extent,dmo,i):
  global count,sm,date
  if dmo == 0:
    pathname = '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/'%(hnum,res,ver)
  else:
    pathname = '/nobackup/afitts/Gadget-2.0.7/production/gizdm%s%s_raw_output/'%(hnum,res)
    pathnameh = '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/'%(hnum,res,ver)
  hist = np.genfromtxt(pathname+'analysis/halo%smerger_hist.txt'%hnum)
  #dm_mass = np.zeros(len(hist))
  #g_mass = np.zeros(len(hist))
  #cgas = np.zeros(len(hist))
  #sfr = np.zeros(len(hist))
  #s_mass = np.zeros(len(hist))
  #time = np.zeros(len(hist))
  #vmax = np.zeros(len(hist))
  GYR = 60*60*24*365*1e9
  try:
    hdf = pd.HDFStore('%sanalysis/dataframes/halo%s%s_giz%s_snap%03d.h5'%(pathname,hnum,res,ver,i))
    print '%sanalysis/dataframes/halo%s%s_giz%s_snap%03d.h5'%(pathname,hnum,res,ver,i)
    red = float(hdf['props']['redshift'])
    rvir = float(hdf['props']['rvir'])
    #if dmo == 1 and extent != 2:
    #  hdf1 = pd.HDFStore('%sanalysis/dataframes/halo%s%s_giz%s_snap%03d.h5'%(pathnameh,hnum,'_11',ver,i))
    #  rhalf = float(hdf1['props']['rhalf'])/(red+1)
    #  hdf1.close()
    #else:
    #rhalf = float(hdf['props']['rhalf'])/(red+1)
    dmass = hdf['particles/dm']['mass'].as_matrix()
    dmm = hdf['particles/dm']['mass']
    dx = hdf['particles/dm']['x'].as_matrix()
    dy = hdf['particles/dm']['y'].as_matrix()
    dz = hdf['particles/dm']['z'].as_matrix()
    dp = np.column_stack((dx,dy,dz))
    if dm == 1:
      dsc = mikecm(fname = dp,nofile=True,pmass=dmm)
      dpos = np.sqrt((dp[:,0]-dsc[0])**2+(dp[:,1]-dsc[1])**2+(dp[:,2]-dsc[2])**2)/.71/(red+1) 
    else:     
     try:
      gmass = hdf['particles/gas']['mass'].as_matrix()
      gpos = hdf['particles/gas']['r'].as_matrix()
      smass = hdf['particles/star']['mass'].as_matrix()
      spos = hdf['particles/star']['r'].as_matrix()
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
      dsm = np.append(dmm,smm)
      dsc = mikecm(fname = dsp,nofile=True,pmass=dsm)
      dpos = np.sqrt((dp[:,0]-dsc[0])**2+(dp[:,1]-dsc[1])**2+(dp[:,2]-dsc[2])**2)/.71
      gpos = np.sqrt((gp[:,0]-dsc[0])**2+(gp[:,1]-dsc[1])**2+(gp[:,2]-dsc[2])**2)/.71
      spos = np.sqrt((sp[:,0]-dsc[0])**2+(sp[:,1]-dsc[1])**2+(sp[:,2]-dsc[2])**2)/.71
     except:
      gmass = hdf['particles/gas']['mass'].as_matrix()
      gpos = hdf['particles/gas']['r'].as_matrix()
      gm = hdf['particles/gas']['mass']
      gx = hdf['particles/gas']['x'].as_matrix()
      gy = hdf['particles/gas']['y'].as_matrix()
      gz = hdf['particles/gas']['z'].as_matrix() 
      gp = np.column_stack((gx,gy,gz))

      dsc = mikecm(fname = dp,nofile=True,pmass=dmass)
      dpos = np.sqrt((dp[:,0]-dsc[0])**2+(dp[:,1]-dsc[1])**2+(dp[:,2]-dsc[2])**2)/.71
      gpos = np.sqrt((gp[:,0]-dsc[0])**2+(gp[:,1]-dsc[1])**2+(gp[:,2]-dsc[2])**2)/.71
      spos = np.zeros(10)
      smass = np.zeros(10)
    if extent == 0 :
      ex = rhalf
    elif extent == 1 :
      ex = 3*rhalf
    elif extent == 2 :
      ex = rvir
    elif extent == 3 :
      ex = .0005
    print ex
    dm_mass = sum(dmass[dpos<=rvir])
    time = hdf['props']['time']
    vmax = hdf['props']['vmax'] 
    if dmo == 0:
      g_mass = sum(gmass[gpos<=ex])
      cgas = sum(gmass[(gpos<=ex) & (hdf['particles/gas']['temp']<1e4)])
      #sfr[i] = sum(hdf['particles/gas']['sfr'][hdf['particles/gas']['r']/(red+1)<=ex])
      sfr=0
      #s_mass[i]=0
      s_mass = sum(smass[spos<=ex])
    hdf.close()
    return time,dm_mass,s_mass,g_mass,vmax,cgas,rvir
  except Exception,f:
    f = str(f)
    f.replace("'","")
    print f,i
    if f != "'No object named props in the file'" and f != "'No object named particles/dm in the file'" and f != "Reindexing only valid with uniquely valued Index objects":
      dm_mass = sum(hdf['particles/dm']['mass'][hdf['particles/dm']['r']/(red+1)<=rvir])
      time = hdf['props']['time']
      vmax = hdf['props']['vmax']
      g_mass = 0
      s_mass = 0
      cgas = 0
      hdf.close()
      return time,dm_mass,s_mass,g_mass,vmax,cgas,rvir
    else:
      hdf.close()
      return 0,0,0,0,0,0,0

def graph_and_compile(dmo,time,dm_mass,s_mass,g_mass,vmax,cgas,rvir):
  if dmo == 0:
   if hnum != '125961':
    ##sfr[sfr<1e-7] = 1e-7
    ##mvirvt.semilogy(time,dm_mass,linewidth=2,color=sm.to_rgba(np.log10(s_mass[184])))
    #mvirvt1 = mvirvt.twinx()
    #mvirvt1.loglog(time,1.42*s_mass, '%s'%clr[count+1],linewidth=2,label=r'%s M$_*$'%(hnum))
    #mvirvt1.set_xlim(5e-1,13.7)
    #mvirvt1.legend(loc=4,prop={'size':20})
    #ax.loglog(time/GYR,1.42*(prop[:,1]-prop[:,2]-prop[:,11]), '%s'%clr[count+2],linewidth=2,label='%s dm'%(hnum))
    #ax.loglog(prop[:,0]/GYR,1.42*(prop[:,11]), '%s--'%clr[count+2],linewidth=2,label='%s gas'%(hnum))
    #ax.loglog(prop[:,0]/GYR,1.42*(prop[:,2]), '%s-.'%clr[count+2],linewidth=2,label='%s star'%(hnum))
    #ax1.plot(prop[:,0]/GYR,1.42*prop[:,2]/max(prop[:,2])/1.42, '%s'%clr[count],label='%s'%(hnum))
    ##mstarvt.semilogy(time,s_mass, color=sm.to_rgba(np.log10(s_mass[184])))
    ##vmaxvt.plot(time,vmax, color=sm.to_rgba(np.log10(s_mass[184])),label='%s'%(hnum))
    ##bfracvt.semilogy(time,(g_mass+s_mass)/dm_mass, color=sm.to_rgba(np.log10(s_mass[184])))
    ##bfracvt.axhline(y=1-0.83120300751)
    #sfrvt.semilogy(time,sfr/ex**2, color=sm.to_rgba(np.log10(s_mass[184])),label='%s'%(hnum))
    #plt.figure()
    #plt.semilogy(time,sfr/ex**2, color=sm.to_rgba(np.log10(s_mass[184])),label='%s'%(hnum))
    #plt.ylabel(r"$\mathrm{SFR}/\mathrm{r}_{1/2}^2$ $(\mathrm{M}_\odot/\mathrm{yr})$")
    #plt.xlabel("$\mathrm{Time (Gyr)}$")
    #plt.ylim(1e1,1e4)
    #plt.savefig('sfr_v_t_halo%s_%s.pdf'%(hnum,date),transparent=True)
    if extent == 2:
      np.savetxt('Halo%s%s_5_12_16_mvir.out'%(hnum,res),np.column_stack((time,dm_mass)), fmt= '%.8e')
      np.savetxt('Halo%s%s_5_12_16_mstar.out'%(hnum,res),np.column_stack((time,s_mass)), fmt= '%.8e')
      np.savetxt('Halo%s%s_5_12_16_mgas.out'%(hnum,res),np.column_stack((time,g_mass)), fmt= '%.8e')
      np.savetxt('Halo%s%s_5_12_16_vmax.out'%(hnum,res),np.column_stack((time,vmax)), fmt= '%.8e')
      np.savetxt('Halo%s%s_5_12_16_coldgas.out'%(hnum,res),np.column_stack((time,cgas)), fmt= '%.8e')
      np.savetxt('Halo%s%s_5_12_16_rvir.out'%(hnum,res),np.column_stack((time,rvir)), fmt= '%.8e')
    elif extent == 3:
      np.savetxt('Halo%s%s_500pc_mvir.out'%(hnum,res),np.column_stack((time,dm_mass)), fmt= '%.8e')
      np.savetxt('Halo%s%s_500pc_mstar.out'%(hnum,res),np.column_stack((time,s_mass)), fmt= '%.8e')
      np.savetxt('Halo%s%s_500pc_mgas.out'%(hnum,res),np.column_stack((time,g_mass)), fmt= '%.8e')
      np.savetxt('Halo%s%s_500pc_vmax.out'%(hnum,res),np.column_stack((time,vmax)), fmt= '%.8e')
      np.savetxt('Halo%s%s_500pc_coldgas.out'%(hnum,res),np.column_stack((time,cgas)), fmt= '%.8e')
  else:
    #mvirvt.loglog(time,0.83120300751*dm_mass,linestyle='--',linewidth=2,label='%s dmo (corrected)'%(hnum))
    if extent == 2:
      np.savetxt('Halo%s%sdmo_central_mvir.out'%(hnum,res),np.column_stack((time,dm_mass)), fmt= '%.8e')
      np.savetxt('Halo%s%sdmo_central_vmax.out'%(hnum,res),np.column_stack((time,vmax)), fmt= '%.8e')
    if extent == 3:
      np.savetxt('Halo%s%sdmo_500pc_mvir.out'%(hnum,res),np.column_stack((time,dm_mass)), fmt= '%.8e')
      np.savetxt('Halo%s%sdmo_500pc_vmax.out'%(hnum,res),np.column_stack((time,vmax)), fmt= '%.8e')
    #ax.loglog(time/GYR,1.42*prop[:,1], '%s:'%clr[count-2],linewidth=2,label='%s dmo'%(hnum))
    #vmaxvt.plot(time,vmax, '%s'%clr[count],linestyle='--',label='%s'%(hnum))
  #count +=1
  return True

fig = plt.figure()
mvirvt = fig.add_subplot(111)
fig1 = plt.figure()
mstarvt = fig1.add_subplot(111)
fig2 = plt.figure()
vmaxvt = fig2.add_subplot(111)
fig3 = plt.figure()
bfracvt = fig3.add_subplot(111)
fig5 = plt.figure()
sfrvt = fig5.add_subplot(111)
fig6 = plt.figure()
rhorhovmstar = fig6.add_subplot(111)
fig7 = plt.figure()
rhorhovarch = fig7.add_subplot(111)

mvirvt.set_ylabel("$\mathrm{M}_{\mathrm{tot}}$($\mathrm{M}_\odot$)")
mvirvt.set_xlabel("$\mathrm{Time (Gyr)}$")
mvirvt.set_xlim(5e-1,13.7)
mvirvt.set_ylim(1e7,2e10)#(1e6,2e9) for extent =0

mstarvt.set_ylabel("$\mathrm{M_{baryon}}$($\mathrm{M}_\odot$)")
mstarvt.set_xlabel("$\mathrm{Time (Gyr)}$")
mstarvt.set_ylim(1e4,5e7)

vmaxvt.set_ylabel("$\mathrm{V}_{\mathrm{max}}$($\mathrm{km/s}$)")
vmaxvt.set_xlabel("$\mathrm{Time (Gyr)}$")

bfracvt.set_ylabel("$\mathrm{f}_{\mathrm{bar}}$")
bfracvt.set_xlabel("$\mathrm{Time (Gyr)}$")
bfracvt.set_ylim(1e-3,5e-1)

sfrvt.set_ylabel(r"$\mathrm{SFR}/\mathrm{r}_{1/2}^2$ $(\mathrm{M}_\odot/\mathrm{yr})$")
sfrvt.set_xlabel("$\mathrm{Time (Gyr)}$")
sfrvt.set_ylim(1e-1,1e2)
clr = ['olive','darkviolet','lime','olive','gray','r','b','k','g','y','m']

extent = 2
hnum = '10v'#'11707','32503','12596','007','848','796','897','948','1016','20192','20910','32257','10q','10v']#'11707','11707','32503','32503','12596','12596','007','007','848','848','897','897','1016','1016','796','796','948','948','20192','20192','20910','20910','32257','32257','10q','10q','10v','10v'] #dmo should be placed after hydro run for irhalf variable.
res = '_11'#,'_13','_13','_11','_13','_13','_13','_13','_13','_13','_13','_13','_11','_11','_13']#['_13','_13','_13','_13','_13','_13','_11','_11','_13','_13','_13','_13','_13','_13','_13','_13','_13','_13','_13','_13','_13','_13','_13','_13','_11','_11','_11','_11']
ver = '5_12_16'#,'5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16']#'11_13','11_13','5_12_16','11_13','5_12_16','11_13','5_12_16','11_13','5_12_16','5_12_16','5_12_16','11_13','5_12_16','11_13','5_12_16','11_13','5_12_16','11_13','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16']#['11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13']
dm = [0,0,0,0,0,0,0,0,0,0,0,0,0,0]#[0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1]#[0,0,0,0,0,0,0,0,0,0,0,0]

count = 0
my_cmap=plt.get_cmap('plasma')
sm = cm.ScalarMappable(cmap=my_cmap,norm=co.Normalize(vmin=5.874853, vmax=7.11806))
sm._A = []
date = time.strftime("%m_%d_%Y")

time,dm_mass,s_mass,g_mass,vmax,cgas,rvir=add_vt(hnum,res,ver,extent,0,rank)
time = comm.gather(time, root=0)
dm_mass = comm.gather(dm_mass, root=0)
s_mass = comm.gather(s_mass, root=0)
g_mass = comm.gather(g_mass, root=0)
vmax = comm.gather(vmax, root=0)
cgas = comm.gather(cgas, root=0)
rvir = comm.gather(rvir, root=0)
if rank == 0:
	graph_and_compile(0,time,dm_mass,s_mass,g_mass,vmax,cgas,rvir)
	cb = fig.colorbar(sm)
	cb.set_label(r'$\mathrm{log M}_*$')
	cb1 = fig1.colorbar(sm)
	cb1.set_label(r'$\mathrm{log M}_*$')
	cb2 = fig2.colorbar(sm)
	cb2.set_label(r'$\mathrm{log M}_*$')
	cb3 = fig3.colorbar(sm)
	cb3.set_label(r'$\mathrm{log M}_*$')
	cb5 = fig5.colorbar(sm)
	cb5.set_label(r'$\mathrm{log M}_*$')

	#fig.savefig('mvir_v_t_all_%s.pdf'%(date),transparent=True)
	#fig1.savefig('mstar_v_t_all_%s.pdf'%(date),transparent=True)
	#fig2.savefig('vmax_v_t_%s.pdf'%date,transparent=True)
	#fig3.savefig('baryfrac_jose_%s.pdf'%(date),transparent=True)
	#fig5.savefig('sfr_v_t_%s.pdf'%date,transparent=True)
	#plt.show('all')
	#plt.close('all')
  
#add_vt(hnum,res,ver,extent,dm)
#hnum = '12596'
#res = '_13'
#ver = '11_13'
#fname = '%s/mfm%s%s_giz%s_raw_output/analysis/Halo%s%s%s_props.out'%(pathname,hnum,res,ver,hnum,res,ver)
#fname1 = '%s/gizdm%s%s_raw_output/analysis/Halo_dm%s_%s_props.out'%(pathname,hnum,res,hnum,res)
#add_vt(hnum,res,ver,extent,0)

#mvirvt.legend(loc=4,prop={'size':10})
#mstarvt.legend(loc=4,prop={'size':10})
#vmaxvt.legend(loc=4, prop={'size':10})
#bfracvt.legend(loc='best')
#sfrvt.legend(loc=4)


