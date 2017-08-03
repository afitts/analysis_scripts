import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator
from matplotlib.ticker import FixedFormatter
from matplotlib import rcParams
rcParams['savefig.bbox'] = 'tight' 
fs=24
plt.ion()
w, h = plt.figaspect(.9)
fig = plt.figure(figsize=(w,h))
sub = fig.add_subplot(111)
sub.set_xlabel(r'# of cores',fontsize=fs, labelpad=5)
sub.set_ylabel(r'Node hrs',fontsize=fs, labelpad=-5)
sub.spines['bottom'].set_linewidth(4)
sub.spines['top'].set_linewidth(4)
sub.spines['right'].set_linewidth(4)
sub.spines['left'].set_linewidth(4)
sub.tick_params('both',length=7.5,width=2,which='minor')
sub.tick_params('both',length=15,width=2,which='major')
sub.xaxis.set_tick_params(labelsize=fs)
sub.yaxis.set_tick_params(labelsize=fs)

ticks = np.linspace(1,3,3)
plt.xticks(ticks,['','100',''])
minor_ticks=[]
for j in range(2,10):
  minor_ticks.append(1+np.log10(j))
for j in range(2,10):
  minor_ticks.append(2+np.log10(j))
for j in range(2,10):
  minor_ticks.append(3+np.log10(j))
sub.xaxis.set_minor_locator(FixedLocator(minor_ticks))
minor_labels = ['','','','','60','','','','','','','500','','','','','','','','','','','','','','','',]#['','30 ','','','60 ','','','']
sub.xaxis.set_minor_formatter(FixedFormatter(minor_labels))
sub.set_xlim(np.log10(50),np.log10(700))
sub.set_ylim(4*1e1,4*2.5e1)
sub.tick_params(which='both',axis='x',labelsize=24)
#ticks = np.linspace(1e4,1e5,2)
#plt.yticks(ticks,['10,000','100,000'])
sub.text(np.log10(70),4*13.5,'Cold Dark Matter',color='r',fontsize=fs)
sub.text(np.log10(120),4*23,'Self-Interacting',color='b',fontsize=fs)
sub.text(np.log10(120),4*22,'Dark Matter',color='b',fontsize=fs)
ms = 12.5
nodes = np.array([64,128,256,512])
cdmSU = np.array([29.8243,15.2581*2,7.94*4,4.15*8])/2
sidmSU = np.array([29.8243*1.2,15.2581*2*1.3,7.64*4*1.25,4.15*8*1.3])/2
sub.plot(np.log10(nodes),4*cdmSU,color = 'r',marker='o',linewidth=2,markersize=ms)
sub.plot(np.log10(nodes),4*sidmSU,color = 'b',marker='o',linewidth=2,markersize=ms)
plt.savefig('XSEDE_strong_scaling_4_15_2017.pdf',transparent=True)

w, h = plt.figaspect(.9)
fig = plt.figure(figsize=(w,h))
sub = fig.add_subplot(111)

sub.set_xlabel(r'# of particles',fontsize=fs, labelpad=5)
sub.set_ylabel(r'Node hrs',fontsize=fs, labelpad=-5)
sub.spines['bottom'].set_linewidth(4)
sub.spines['top'].set_linewidth(4)
sub.spines['right'].set_linewidth(4)
sub.spines['left'].set_linewidth(4)
sub.tick_params('both',length=7.5,width=2,which='minor')
sub.tick_params('both',length=15,width=2,which='major')
sub.xaxis.set_tick_params(labelsize=24)
sub.yaxis.set_tick_params(labelsize=24)
sub.set_xlim(3e6,1.5e9)
sub.set_ylim(1e1,1.5e4)
#sub.text(np.log10(70),2.5e4,'Cold Dark Matter',color='r',fontsize=fs)
#sub.text(np.log10(120),5e4,'Self-Interacting Dark Matter',color='b',fontsize=fs)

parts = np.array([6.7e6,5.2e7,4.1e8])
cdmSU = np.array([420*1.1,2820,5.15e4/1.47])/16/2
cdmSU_stam1 = np.array([650*1.2,5700,9e4/1.47])/16/2
sidmSU = np.array([420*1.2*1.2,2820*1.3,5.15e4*1.35/1.47])/16/2
sub.text(1.5e8,2.74*2.5e2,'Projection',color='k',fontsize=fs-4)
sub.text(1.05e8,2.74*1.7e2,'for ultra-high',color='k',fontsize=fs-4)
sub.text(8e7,2.74*1.2e2,'resolution runs',color='k',fontsize=fs-4)
sub.loglog(parts,4*cdmSU_stam1,color = 'g',marker='*',linewidth=2,markersize=ms,label = 'CDM (Stampede)')
sub.loglog(parts[:2],4*sidmSU[:2],color = 'b',marker='o',linewidth=2,markersize=ms,label = 'SIDM (Stampede2)')
sub.loglog(parts[1:3],4*sidmSU[1:3],color = 'b',marker='o',linewidth=2,markersize=ms,linestyle='--')
sub.loglog(parts[:2],4*cdmSU[:2],color = 'r',marker='o',linewidth=2,markersize=ms,label = 'CDM (Stampede2)')
sub.loglog(parts[1:3],4*cdmSU[1:3],color = 'r',marker='o',linewidth=2,markersize=ms,linestyle='--')
sub.tick_params(which='both',axis='x',pad=5)
leg = sub.legend(loc=4,frameon=False,prop={'size':15})
lcolors = ['g','b','r']
for color,text,item in zip(lcolors,leg.get_texts(),leg.legendHandles):
    text.set_color(color)
    text.set_weight('heavy')
    item.set_visible(False)
plt.savefig('XSEDE_weak_scaling_4_15_2017.pdf',transparent=True)
