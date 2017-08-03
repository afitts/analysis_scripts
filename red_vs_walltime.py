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

clr = ['olive','darkviolet','lime','olive','gray','r','b','k','g','y','m']
hnum = ['32257','11707','32503','125961','12596','007','2','897','1016','796','948'] #dmo should be placed after hydro run for irhalf variable.
res = ['_13','_13','_13','_13','_13','_11','_13','_13','_13','_13','_13','_13','_13','_13','_13','_13']
ver = ['11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13','11_13']
dm = [0,0,0,0,0,0,0,0,0,0,0,0]
for w in np.arange(len(hnum)):
 if hnum[w] == '32257' or hnum[w] == '11707' or hnum[w] == '32503' or hnum[w] == '12596':
  if dm[w] == 0:
    pathname = '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz%s_raw_output/'%(hnum[w],res[w],ver[w])
  else:
    pathname = '/nobackup/afitts/Gadget-2.0.7/production/gizdm%s%s_raw_output/'%(hnum[w],res[w])
  print hnum[w]
  f = '%scpu.txt'%pathname
  a = np.zeros(1)
  wt = np.zeros(1)
  for line in open(f):
    if "Step" in line:
      a = np.append(a,float(line.split()[3].strip(',')))
    if "total" in line:
      #print line,float(line.split()[1])
      wt = np.append(wt,float(line.split()[1]))
  np.savetxt('Halo%s_scale_vs_totwalltime.out'%hnum[w],np.column_stack((a,wt)))
  plt.semilogy(a,wt)
plt.xlabel('Scale Factor a')
plt.ylabel('walltime (s)')
date = time.strftime("%m_%d_%Y")
plt.savefig('scale_vs_walltime_%s.pdf'%date,transparent=True)
plt.show()
