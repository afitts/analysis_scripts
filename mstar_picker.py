import numpy as np
import sys 
import glob
#import matplotlib
#matplotlib.use('Agg')
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
import matplotlib.animation as animation
from matplotlib.ticker import FixedLocator
from matplotlib.ticker import FixedFormatter

def _a_dot(a, h0, om_m, om_l):
    om_k = 1.0 - om_m - om_l
    return h0 * a * np.sqrt(om_m * (a ** -3) + om_k * (a ** -2) + om_l)


def _a_dot_recip(*args):
    return 1. / _a_dot(*args)

hnum =['848','1016','32503']#['11707','32503','12596','007','848','796','20192','20910','897','948','1016','32257','10q','10v','1084']
res = ['_13','_13','_13','_11','_13','_13','_13','_13','_13','_13','_13','_13','_11','_11','_13']
ver = ['5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16']
snap = [184,184,184,184,184,184,184,184,184,184,184,184,600,600,184]

for i in np.arange(len(hnum)):
	print 'HALO ',hnum[i]
	for k in np.arange(185):
		a = pd.HDFStore('/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/analysis/dataframes/halo%s%s_giz%s_snap%03d.h5'%(hnum[i],res[i],ver[i],hnum[i],res[i],ver[i],k))
		try:
			pid = a['particles/star'].index.values
			mass = a['particles/star']['mass']
			x = a['particles/star']['x'].as_matrix()*1000/.71
			y = a['particles/star']['y'].as_matrix()*1000/.71
			z = a['particles/star']['z'].as_matrix()*1000/.71
			vx = a['particles/star']['vx'].as_matrix()*1000
			vy = a['particles/star']['vy'].as_matrix()*1000
			vz = a['particles/star']['vz'].as_matrix()*1000
			he = a['particles/star']['metal_He'].as_matrix()
			fe = a['particles/star']['metal_Fe'].as_matrix()
			metal = a['particles/star']['metal_tot'].as_matrix()
			numH = mass*(1-(metal+he))/(1.6733e-24/1.99e33)
			numFe = mass*fe/(9.27e-23/1.99e33)
			meta = numFe/numH
			nofloor = meta[meta>4.54877795e-09]
			avgnum = np.mean(nofloor)
			fullmetal = np.log10(meta)+4.5
			metal = np.log10(avgnum)+4.5
			sft = a['particles/star']['sft'].as_matrix()
			YEAR = 60*60*24*365.
			h0 = 71
			om_l = 0.734
			om_m = 0.266
			conv = 3.085677581e+19
			for j in np.arange(len(sft)):
				sft[j] = 13.736583-scipy.integrate.quad(_a_dot_recip, 0, sft[j], (h0, om_m, om_l))[0]*conv/1e9/YEAR
			time = np.zeros(len(sft))
			time += np.float(a['props']['time'])
			np.savetxt('Halo%s_%03d_stars.out'%(hnum[i],k),np.column_stack((pid,mass,x,y,z,vx,vy,vz,sft,fullmetal,time)),header = '(0) ID	(1) MASS (M_sun)	(2) X (kpc)	(3) Y	(4) Z	(5) Vx (km/s)	(6) Vy	(7) Vz	(8) Age (Gyr)	(9) [Fe/H]	(10) Time of snapshot (Gyr)')
			print k
		except Exception,e:
			print e,'No stars'
			try:
				time = np.float(a['props']['time'])
				np.savetxt('Halo%s_%03d_stars.out'%(hnum[i],k),np.column_stack((0,0,0,0,0,0,0,0,0,0,time)),header = '(0) ID	(1) MASS (M_sun)	(2) X (kpc)	(3) Y	(4) Z	(5) Vx (km/s)	(6) Vy	(7) Vz	(8) Age (Gyr)	(9) [Fe/H]	(10) Time of snapshot (Gyr)')
			except:
				print 'No time'
		a.close()
