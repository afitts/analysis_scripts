import numpy as np
import matplotlib.pyplot as plt
import pylab 
import scipy.optimize as opt
import scipy.special as sp
def nfw(r,mvir,c):
    gc = 1./(np.log(1.+c)-(c/(1.+c)))
    rhocrit = 139.94 #M_sun/kpc^3
    rvir = (3/4.*mvir*(1./(np.pi*96.45*rhocrit)))**(1/3.)
    rs = rvir/c
    #return mvir*gc*(np.log(1+r/rs)-r/rs*(1+r/rs)**-1)
    rho0 = rhocrit*96.45/3*gc*c**3.
    return np.log(rho0*(r/rs)**-1*(1+(r/rs))**-2)#np.log(rho0)-np.log(r/rs)-2*np.log(1+(r/rs))

def corenfw(r,mvir,c,n,rc):
    G = 4.3e-6 #in kpc/M_sun (km/s)^2
    gc = 1./(np.log(1.+c)-(c/(1.+c)))
    rhocrit = 139.94 #M_sun/kpc^3
    rvir = (3/4.*mvir*(1./(np.pi*96.45*rhocrit)))**(1/3.)
    rs = rvir/c
    rho0 = rhocrit*96.45/3*gc*c**3.
    nfw = rho0*(r/rs)**-1*(1+(r/rs))**-2
    mnfw = lambda r:mvir*gc*(np.log(1+r/rs)-r/rs*(1+r/rs)**-1)
    mnfw_rs = mnfw(rs)
    #rc = eta*rhalf
    f = np.tanh(r/rc)
    #tdyn = 2*np.pi*np.sqrt(rs**3/(G*mnfw_rs))
    #q = kappa_tsf/tdyn
    #n = np.tanh(q)
    return np.log(f**n*nfw+n*f**(n-1)*(1-f**2)/(4*np.pi*r**2*rc)*mnfw(r))

def psuedoiso(r,cenden,rcore):
    return cenden*(1+(r/rcore)**2)**-1

def burkert(r,cenden,rcore):
	return cenden/((1+(r/rcore))*(1+(r/rcore)**2))

def fit_einasto_mass(x, mvir, c0):
    rhocrit = 139.94 #M_sun/kpc^3
    my_alpha = 0.17
    rvir = (3/4.*mvir*(1./(np.pi*96.45*rhocrit)))**(1/3.)
    return np.log(mvir/(4*np.pi*rvir**3)*c0**3/gfac(c0, my_alpha)) - 2/my_alpha*((x/(rvir/c0))**(my_alpha)-1)

def gfac(xx, ind, po=3e0):
    return (1e0/ind*np.exp(2e0/ind)*(2e0/ind)**(-po/ind)*sp.gammainc(po/ind, 2e0/ind*xx**ind)*sp.gamma(po/ind))

class radpro(object):
        """Radial density profile """
        def __init__(self):
                self.fig = plt.figure()
                self.sub = self.fig.add_subplot(111)
                self.sub.spines['bottom'].set_linewidth(4)
                self.sub.spines['top'].set_linewidth(4)
                self.sub.spines['right'].set_linewidth(4)
                self.sub.spines['left'].set_linewidth(4)
                self.sub.tick_params('both',length=5,width=2,which='minor')
                self.sub.tick_params('both',length=10,width=2,which='major')
                self.sub.xaxis.set_tick_params(labelsize=20)
                self.sub.yaxis.set_tick_params(labelsize=20)
                self.sub.xaxis.set_label_coords(.48,-.07)
                self.sub.set_xlabel(r'$Radius$ $(kpc)$',fontsize=20, labelpad=-10)
                self.sub.set_ylabel(r'$\rho$ $(M_\odot/kpc^3$)',fontsize=20, labelpad=-5)
                self.sub.set_xlim(.06,60)
                self.sub.set_ylim(1e4,1e9)
                pylab.rcParams['xtick.major.pad']='6'
                pylab.rcParams['ytick.major.pad']='6'

        def add_line(self,x,totden,res):
                #if res == '':  # lvl 12
                #       clr = 'r'
                #       lbl = 'LOW '
                #else:          # lvl 13
                #       clr = 'k'
                #       lbl = 'HI '
                if res == '_14':   # turb diff
                        clr = 'k'
                        lbl = 'Z14'#'in situ'
                else:           # regular
                        clr = 'm'
                        lbl = 'Z13'#'merger'
                self.sub.loglog(x,totden, color = '%s'%clr,linewidth=4,label = 'HYDRO %s'%lbl)
	def fit_einasto(self):
		dmo = np.genfromtxt('radden_out/Halo848_raddenZ14_dmo.out')
		prad = .095
		param_bounds=([4e9,0.5],[1.4e10,26])
		(fit, cmatrix)= opt.curve_fit(fit_einasto_mass,dmo[:,0][(dmo[:,0]>prad) & (dmo[:,0]<10)],np.log(dmo[:,1][(dmo[:,0]>prad) &(dmo[:,0]<10)]),bounds = param_bounds)
		self.sub.loglog(dmo[:,0],np.exp(fit_einasto_mass(dmo[:,0],*fit)),color = 'darkgrey',linestyle='--',linewidth=2,label = 'Einasto',zorder=100)

        def add_dmoline(self,x,totden,res):
                if res == '_14':   # turb diff
                        clr = 'k'
                        lbl = 'Z14'#'TD'
			dmo = np.genfromtxt('radden_out/Halo848_raddenZ14_dmo.out')
			x = dmo[:,0]
			totden = dmo[:,1]
                else:           # regular
                        clr = 'm'
                        lbl = 'Z13'
                        dmo = np.genfromtxt('radden_out/Halo848_raddenZ13_dmo.out')
                        x = dmo[:,0]
                        totden = dmo[:,1]
                barycorrect = 0.83120300751
                self.sub.loglog(x,totden, color = '%s'%clr,linestyle='--',linewidth=4,label = 'DMO %s'%lbl)#self.sub.plot(dmox,dmomass/(4*3.14159*dmox**3)*barycorrect/dlogr, color = 'r',linestyle='',linewidth=4, label = '%sDMO'%lbl,marker='s')
	def add_fit(self,x,fit,ver):
		if ver =='bur': #burkert
			clr = 'b'
			lbl = 'Burkert'
                elif ver =='piso': #psuedoiso
                        clr = 'y'
                        lbl = 'Psuedo Iso'
                elif ver =='cnfw': #cored NFW
                        clr = 'g'
                        lbl = 'Cored NFW'
                self.sub.loglog(x,fit, color = '%s'%clr,linewidth=4,label = '%s fit'%lbl)

        def save(self,date,hnum):
                self.sub.legend(loc=1,prop={'size':10})
                self.fig.savefig('radden_halo%s_metal_%s.pdf'%(hnum,date),transparent=True)
                self.fig.show()

class vcircvr (object):
        """ Circular velocity plots v r. Includes points where rhalf is."""
        def __init__(self):
                self.fig = plt.figure()
                self.sub = self.fig.add_subplot(111)
                #self.sm = sm
                self.sub.set_ylabel(r'$V_{\mathrm{circ}}$ $\mathrm{[km\:s^{-1}]}$',fontsize = 20)
                self.sub.set_xlabel(r'$r$ $\mathrm{[kpc]}$',fontsize = 20)
                self.sub.set_xlim(0,5e0)
		self.sub.set_ylim(0,3.5e1)
                #cb13 = fig13.colorbar(sm)
                #cb13.set_label(r'$\mathrm{log}$ $M_\star$')

        def add_line(self,x,totmass,res):
                G = 4.3e-6 #in kpc/M_sun (km/s)^2
                if res == '_14':  
                        clr = 'k'
                        lbl = 'Z14'
                elif res == '_13':           
                        clr = 'm'
                        lbl = 'Z13'
                elif res == '':           
                        clr = 'c'
                        lbl = 'Z12'
                #v = 0
                #while x[v] <=rhalf:
                #        v+=1
                #if (x[v]-rhalf)>(rhalf-x[v-1]):
                #        v-=1
                cmass = np.cumsum(totmass)
                print 'vmax is \n',max(np.sqrt(G*cmass/x))
                self.sub.plot(x,np.sqrt(G*cmass/x),linewidth = 2,color='%s'%clr,label='%s'%lbl)
                #self.sub.scatter(x[v],np.sqrt(G*cmass[v]/x[v]),marker='s', s = 40, color=self.sm.to_rgba(np.log10(sum(smass))))

        def save(self,date):
                self.sub.legend(loc=2,frameon=False,prop={'size':10})
                self.fig.savefig('vcirc_profile_halo848_hydro_%s.pdf'%date,transparent=True)
                self.fig.show()
                plt.close()
