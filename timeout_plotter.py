import numpy as np
import sys 
import glob
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as co
import pylab
import pygadgetreader as pg
import scipy.integrate
import time
import scipy.stats as stats
import pandas as pd
import scipy.optimize as opt
from matplotlib import rcParams
from matplotlib.ticker import FixedLocator
from matplotlib.ticker import FixedFormatter

rcParams['lines.linewidth'] = 1.5
rcParams['axes.linewidth'] = 3
rcParams['xtick.major.size'] = 10
rcParams['xtick.minor.size'] = 5
rcParams['xtick.major.width'] = 2
rcParams['xtick.minor.width'] = 2
rcParams['ytick.major.size'] = 10
rcParams['ytick.minor.size'] = 5
rcParams['ytick.major.width'] = 2
rcParams['ytick.minor.width'] = 2
rcParams['font.size'] = 20
rcParams['savefig.bbox'] = 'tight'

hnum = ['11707','32503','12596','007','848','796','20192','20910','897','948','1016','32257','10q','10v','1084'] #dmo should be placed after hydro run for irhalf variable.
res = ['_13','_13','_13','_11','_13','_13','_13','_13','_13','_13','_13','_13','_11','_11','_13','_13']
clr = ['olive','darkviolet','lime','olive','gray','r','b','k','g','y','m','k','k','k']
fig = plt.figure()
fig1 = plt.figure()
fig2 = plt.figure()
fig3 = plt.figure()
fig4 = plt.figure()
fig5 = plt.figure()
#xpos = [6e9,8e9,1e10,1.2e10,1.4e10]
#labels=[r'$6\times10^9$',r'$8\times10^9$',r'$1\times10^{10}$',r'$1.2\times10^{10}$',r'$1.4\times10^{10}$']
#plt.xticks(xpos,labels,fontsize=16)
fig6 = plt.figure()
fig7 = plt.figure()
fig8 = plt.figure()
mvirvt = fig.add_subplot(111)
mstarvt = fig1.add_subplot(111)
mstarmhalo = fig2.add_subplot(111)
rhohyvt = fig3.add_subplot(111)
vmaxvt = fig4.add_subplot(111)
mhalovmstar = fig5.add_subplot(111)
corevt = fig6.add_subplot(111)
bfracvt = fig7.add_subplot(111)
mvirvz = fig8.add_subplot(111)
#mvirvz1 = fig9.add_subplot(111)
my_cmap=plt.get_cmap('plasma')
sm = cm.ScalarMappable(cmap=my_cmap,norm=co.Normalize(vmin=5.0025, vmax=7.176))
sm._A = []
my_cmap1=plt.get_cmap('viridis')
sm1 = cm.ScalarMappable(cmap=my_cmap1,norm=co.Normalize(vmin=np.log10(0.00793650793651), vmax=np.log10(1)))
sm1._A = []
sm2 = cm.ScalarMappable(cmap=my_cmap,norm=co.Normalize(vmin=0, vmax=1))
sm2._A = []
for i in np.arange(len(hnum)):
 if hnum[i] != '125961' and hnum[i] !='209102':
  out = np.genfromtxt('mvir_out/Halo%s%s_5_12_16_mvir.out'%(hnum[i],res[i]))
  sout = np.genfromtxt('mvir_out/Halo%s%s_5_12_16_mstar.out'%(hnum[i],res[i])) 
  ssout = np.genfromtxt('mvir_out/Halo%s%s_5_12_16_mstar.out'%(hnum[i],res[i])) 
  gout = np.genfromtxt('mvir_out/Halo%s%s_5_12_16_mgas.out'%(hnum[i],res[i]))
  vout = np.genfromtxt('mvir_out/Halo%s%s_5_12_16_vmax.out'%(hnum[i],res[i]))
  #dmoout = np.genfromtxt('Halo%s%sdmo_5_12_16_mvir.out'%(hnum[i],res[i]))
  #archout = np.genfromtxt('Halo%s%s_cumusfr_rhalf.out'%(hnum[i],res[i]))
  #rhohy = np.genfromtxt('Halo%s%s_rhohy_halfrhalf.out'%(hnum[i],res[i]))
  #rhodmo = np.genfromtxt('Halo%s%s_rhodmo_halfrhalf.out'%(hnum[i],res[i]))
  a = np.genfromtxt('output_times.txt')
  w = 0
  while ssout[w,1]<.5*ssout[184,1]:
    w+=1
  print hnum[i],a[w]
  v = 0
  #while ssout[v,1]<.9*ssout[184,1]:
  #  v+=1
  if hnum[i] == '1084':
   mvirvt.semilogy(out[:,0],out[:,1],color='k',linewidth=2,label=r'%s'%(hnum[i]))
   mvirvz.loglog(1/a,out[:,1],color='k',linewidth=2)
   vmaxvt.plot(vout[:,0],vout[:,1], color='k',label='%s'%(hnum))
  else:
   if hnum[i] == '10q' or hnum[i] == '10v':
     a = np.genfromtxt('snapshot_scale-factors.txt')
   mvirvt.semilogy(out[:,0],out[:,1],color=sm.to_rgba(np.log10(ssout[len(ssout)-1,1])),linewidth=2,label=r'%s'%(hnum[i]),alpha=1)
   mvirvz.loglog(1./a,out[:,1],color=sm.to_rgba(np.log10(ssout[len(ssout)-1,1])),linewidth=2,alpha=1)
   #mvirvz1.semilogx(a,out[:,1]/ok,color=sm.to_rgba(np.log10(ssout[len(ssout)-1,1])),linewidth=2,label=r'%s'%(hnum[i]),alpha=1)
   vmaxvt.plot(vout[:,0],vout[:,1], color=sm.to_rgba(np.log10(ssout[len(ssout)-1,1])),label='%s'%(hnum))
  mstarvt.semilogy(sout[:,0],sout[:,1],color=sm.to_rgba(np.log10(ssout[len(ssout)-1,1])),linewidth=2,label=r'%s'%(hnum[i]))
  print out[len(ssout)-1,1],sout[len(ssout)-1,1]
  mhalovmstar.scatter(out[len(ssout)-1,1],sout[len(ssout)-1,1],color='k',s=80,label='%s'%(hnum))
  mstarmhalo.loglog(out[:,0],sout[:,1]/out[:,1],color=sm.to_rgba(np.log10(ssout[len(ssout)-1,1])),linewidth=2,label=r'%s'%(hnum[i]))
  #print hnum[i], (out[184,1]+sout[184,1]+gout[184,1])/(dmoout[184,1]*0.83120300751)
  #corevt.scatter(a[w],(out[184,1])/(dmoout[184,1]*0.83120300751),s=100,color=sm.to_rgba(np.log10(ssout[184,1])),label='%s'%(hnum))
  bfracvt.plot(out[:,0],(gout[:,1]+sout[:,1])/out[:,1], color=sm.to_rgba(np.log10(sout[len(ssout)-1,1])))
  bfracvt.axhline(y=1-0.83120300751)
  corevt.scatter(out[184,0]-out[w,0],out[184,0]-out[v,0],s=100,color=sm.to_rgba(np.log10(ssout[184,1])),label='%s'%(hnum))
  #rhohyvt.semilogy(rhohy[:,0],rhohy[:,1]/rhodmo[:,1],color=sm.to_rgba(np.log10(sout[184,1])),linewidth=1,label=r'%s'%(hnum[i]))
 #if hnum[i] == '1084':
 # a = np.genfromtxt('output_times.txt')
 # out = np.genfromtxt('mvir_out/Halo%s%s_5_12_16_mvir.out'%(hnum[i],res[i]))
 # mvirvt.semilogy(out[:,0],out[:,1],color='darkgrey',linewidth=2,label=r'%s'%(hnum[i]))
 # vout = np.genfromtxt('mvir_out/Halo%s%s_5_12_16_vmax.out'%(hnum[i],res[i]))
 # vmaxvt.plot(vout[:,0],vout[:,1], color='darkgrey')
 if hnum[i] == '209102':
  a = np.genfromtxt('output_times.txt')
  out = np.genfromtxt('mvir_out/Halo%s%s_5_12_16_mvir.out'%(hnum[i],res[i]))
  mvirvt.semilogy(out[:,0],out[:,1],color='k',linewidth=2,label=r'%s'%(hnum[i]))
  vout = np.genfromtxt('mvir_out/Halo%s%s_5_12_16_vmax.out'%(hnum[i],res[i]))
  vmaxvt.plot(vout[:,0],vout[:,1], color='k')
 if hnum[i] == '1016':
  oka = np.load('okamoto_dat.npy')
  mvirvz.loglog(1+oka[0],oka[3],color = '0.5',linewidth = 3,linestyle='--',alpha=1,label='Fakhouri 2010')
  mvirvz.loglog(1+oka[0],oka[2],color = 'k',linewidth = 3,linestyle='--',alpha =1,label='Okamoto 2008')
  mvirvz.legend(loc=4,prop={'size':15})
#mstarvt.legend(loc=4,prop={'size':9},ncol=2)

#mvirvt.set_xlim(1/0.0909090909,1)

mvirvt.set_ylim(1e7,2e10)#1e6,2e9)#
mvirvt.set_xlabel('$\mathrm{Time}$ $\mathrm{(Gyr)}$')#expansion factor (a)')#'1+z'
mvirvt.xaxis.set_label_coords(.48,-.06)

cb = fig.colorbar(sm)
cb.set_label(r'$\mathrm{log M}_\star(z=0)$')
mvirvt.set_ylabel('$\mathrm{M_{vir}}$ $\mathrm{(M_\odot)}$')
#mstarvt.set_xlim(0,13.7)
mstarvt.set_ylim(2e3,2e7)
mstarvt.set_xlabel('$\mathrm{Time}$ $\mathrm{(Gyr)}$')
mstarvt.set_ylabel('$M_{\star}$ $(M_\odot)$')#/M_{\star,arch}}$
cb1 = fig1.colorbar(sm)
cb1.set_label(r'$\mathrm{log M}_\star(z=0)$')
mstarmhalo.set_xlabel('expansion factor (a)')#\mathrm{Time}$ $\mathrm{(Gyr)}$')
mstarmhalo.set_ylabel('$\mathrm{M_\star/M_{vir}}$')
mstarmhalo.set_xlim(0.0909090909,1)
cb2 = fig2.colorbar(sm)
cb2.set_label(r'$\mathrm{log M}_\star(z=0)$')
rhohyvt.set_xlabel('$\mathrm{Time}$ $\mathrm{(Gyr)}$')
rhohyvt.set_ylabel(r'$\rho_\mathrm{hydro}/\rho_{\mathrm{dmo}}$ $(1/2*r_\mathrm{half})$')
cb3 = fig3.colorbar(sm)
cb3.set_label(r'$\mathrm{log M}_\star(z=0)$')
vmaxvt.set_ylabel("$\mathrm{V}_{\mathrm{max}}$($\mathrm{km/s}$)")
vmaxvt.set_xlabel(r'$\mathrm{Time}$ $\mathrm{(Gyr)}$')#'expansion factor (a)')#
#vmaxvt.set_xlim(0.0909090909,1)
cb4 = fig4.colorbar(sm)
cb4.set_label(r'$\mathrm{log M}_\star(z=0)$')
mhalovmstar.set_xlabel('$\mathrm{M_{vir}}$ $\mathrm{(M_\odot)}$')
mhalovmstar.set_ylabel('$M_{\star}$ $(M_\odot)$')#/M_{\star,arch}}$
mhalovmstar.set_xscale('log')
mhalovmstar.set_yscale('log')
mhalovmstar.set_ylim(7e4,2e7)
mhalovmstar.set_xlim(5e9,1.5e10)
cb5 = fig5.colorbar(sm1)
cb5.set_label(r'log a')
corevt.set_yscale('log')
corevt.set_ylabel(r'$\mathrm{\rho}_{\mathrm{dm(hydro)}}/\mathrm{\rho}_{\mathrm{dmo(corrected)}}\mathrm{(500pc)}$',fontsize=20)#$\tau_q$ $\mathrm{(Gyr)}$',fontsize=14)#
corevt.set_xlabel(r'$a_{50\%}$',fontsize=20)#'$\tau_{50}$ $\mathrm{(Gyr)}$',fontsize=14)#
corevt.set_xscale('log')
mvirvz.set_ylabel('$\mathrm{M_{vir}}$ $\mathrm{(M_\odot)}$')
mvirvz.set_xlabel('$(1+z)$')
mvirvz.set_xlim(11,1)
mvirvz.set_ylim(1e7,2e10)
#plt.tick_params(which='minor',labelsize=16,pad=5)
#ticks = np.linspace(5,10,2)
#plt.xticks(ticks,['5','10'],fontsize=16)#[r'$10^1$',r'$10^2$'],fontsize=10)
#minor_ticks=[]
#for j in np.linspace(0,14,15):
#  minor_ticks.append(j)
#corevt.xaxis.set_minor_locator(FixedLocator(minor_ticks))
#minor_labels = ['',' ','','',' ','','','']
#corevt.xaxis.set_minor_formatter(FixedFormatter(minor_labels))
#ticks = np.linspace(5,10,2)
#plt.yticks(ticks,['5','10'],fontsize=16)#[r'$10^1$',r'$10^2$'],fontsize=10)
#minor_ticks=[]
#for j in np.linspace(0,14,15):
#  minor_ticks.append(j)
#corevt.yaxis.set_minor_locator(FixedLocator(minor_ticks))
#minor_labels = ['',' ','','',' ','','','']
#corevt.yaxis.set_minor_formatter(FixedFormatter(minor_labels))
#ticks = np.linspace(.1,1,2)
#plt.xticks(ticks,['0.1','1'],fontsize=16)#[r'$10^1$',r'$10^2$'],fontsize=10)
#minor_ticks=[]
#for j in np.linspace(0.2,0.9,8):
#  minor_ticks.append(j)
#corevt.xaxis.set_minor_locator(FixedLocator(minor_ticks))
#minor_labels = ['','0.3 ','','','0.6 ','','','']
#corevt.xaxis.set_minor_formatter(FixedFormatter(minor_labels))
#corevt.set_ylim(1e-1,1.1)
#corevt.set_xlim(1.5e-1,1)
#cb6 = fig6.colorbar(sm)
#cb6.set_label(r'$\mathrm{log M}_\star$ at z=0')#% of Stellar Mass Formed')
bfracvt.set_ylabel("$\mathrm{f}_{\mathrm{bar}}$")
bfracvt.set_xlabel("$\mathrm{Time (Gyr)}$")
bfracvt.set_ylim(1e-3,5e-1)
cb7 = fig7.colorbar(sm)
cb7.set_label(r'$\mathrm{log M}_\star$(z=0)')
date = time.strftime("%m_%d_%Y")
fig.savefig('mvir_v_t_%s.pdf'%(date),transparent=True)
fig1.savefig('mstar%s_v_t_%s.pdf'%(hnum[0],date),transparent=True)
fig2.savefig('mstarmhalo%s_v_a_all_%s.pdf'%(hnum[0],date),transparent=True)
fig3.savefig('rhohyrhodmo%s_v_t_%s.pdf'%(hnum[0],date),transparent=True)
fig4.savefig('vmax_v_t_%s.pdf'%(date),transparent=True)
fig5.savefig('mhalo_v_mstar%s_%s.pdf'%(hnum[0],date),transparent=True)
fig6.savefig('core_v_loga_%s.pdf'%(date),transparent=True)
fig7.savefig('baryfrac_jose_%s.pdf'%(date),transparent=True)
fig8.savefig('mvir_v_z_%s.pdf'%(date),transparent=True)
plt.show()
