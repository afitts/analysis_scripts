import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['lines.linewidth'] = 2
rcParams['axes.linewidth'] = 3
rcParams['xtick.major.size'] = 10
rcParams['xtick.minor.size'] = 5
rcParams['xtick.major.width'] = 1.5
rcParams['xtick.minor.width'] = 1.5
rcParams['ytick.major.size'] = 10
rcParams['ytick.minor.size'] = 5
rcParams['ytick.major.width'] = 1.5
rcParams['ytick.minor.width'] = 1.5
rcParams['font.size'] = 16
rcParams['figure.figsize'] = [9,3.8]

hx1 = np.array([11468800,91750400])
hx2 = np.array([91750400,734003200])
dx = np.array([6553600,52428800,419430400])

dy = np.array([620,5480,84480])
hy1 = np.array([5300,69000])
hy2 = np.array([69000,1.2e6])

ddy = np.array([16.76,18.73,36])
hhy1 = np.array([82.81,134.77])
hhy2 = np.array([134.77,292.97])


fig, (ax1, ax2) = plt.subplots(1,2, sharex=False, sharey=False)
fig.tight_layout()
fig.subplots_adjust(wspace=0.25)
ax1.loglog(hx1,hy1, color='b', marker = '*',markersize = 15,linewidth=2)
ax1.loglog(hx2,hy2, linestyle='--',color='b', marker = '*',markersize = 15,linewidth=2)
ax1.loglog(dx,dy, color='r', marker = 'o',markersize = 8,linewidth=2)
ax1.text(5.5e6,1e5,'full gravity + ',fontsize=15, color='b')
ax1.text(5.5e6,5e4,'hydro + stars',fontsize=15, color='b')
ax1.text(4e7,8e5,'projection for ultra-',fontsize=9, color='k')
ax1.text(4e7,5e5,'high resolution run',fontsize=9, color='k')
ax1.text(4.5e7,1.243e3,'dark matter only',fontsize=15, color='r')
ax1.set_xlim(4e6,1.5e9)
ax1.set_ylim(3e2,3e6)
ax1.set_xlabel('paritcle number',labelpad=-2, fontsize=14)
ax1.set_ylabel('CPU time to z=0 [hr]',labelpad=-2, fontsize=14)

ax2.loglog(hx1,hhy1, color='b', marker = '*',markersize = 15,linewidth=2)
ax2.loglog(hx2,hhy2, linestyle='--',color='b', marker = '*',markersize = 15,linewidth=2)
ax2.loglog(dx,ddy, color='r', marker = 'o',markersize = 8,linewidth=2)
ax2.text(7.5e6,174,'full gravity + ',fontsize=15, color='b')
ax2.text(7.5e6,127,'hydro + stars',fontsize=15, color='b')
ax2.text(6.72e7,287,'projection for ultra-',fontsize=9, color='k')
ax2.text(6.72e7,240,'high resolution run',fontsize=9, color='k')
ax2.text(2.7e7,12,'dark matter only',fontsize=15, color='r')
ax2.set_xlim(5e6,1e9)
ax2.set_ylim(8e0,4e2)
ax2.set_xlabel('paritcle number',labelpad=-2, fontsize=14)
ax2.set_ylabel('wall time to z=0 [hr]', fontsize=14)
plt.savefig('weakcputime.pdf',pad_inches=0)
plt.show()

