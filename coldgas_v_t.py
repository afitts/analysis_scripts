import pygadgetreader as pg
import numpy as np
import pandas as pd
import sys
from pyd import mikecm
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

hnum = sys.argv[1]
extent = np.linspace(0,184,185)
i = extent[rank]
my_N = 1
N = my_N * comm.size

#if rank == 0:
#    A = np.arange(N, dtype=np.float64)
#else:
#    A = np.empty(N, dtype=np.float64)
#
#my_A = np.empty(my_N, dtype=np.float64)
print i
try:
	hdf = pd.HDFStore('/nobackup/afitts/Gadget-2.0.7/production/mfm1016og_13_giz5_12_16_raw_output/analysis/dataframes/halo1016_13_giz5_12_16_snap%03d.h5'%(i))
	print 'bob',i
	gm = hdf['particles/gas']['mass'].as_matrix()
	gtemp = hdf['particles/gas']['temp'].as_matrix()
        my_C = sum(gm[gtemp<1e4])
	hdf.close()
	print 'my name is', i	
except Exception,e:
	print e
	print 'no gas?'
	my_C = 0
my_C = comm.gather(my_C, root=0)
if rank == 0:
	a = np.genfromtxt('/nobackup/afitts/analysis_scripts/mvir_out/Halo848_13_5_12_16_mstar.out')
	np.savetxt('halo%s_coldgas.txt'%(hnum),np.column_stack((a[:,0],my_C)))


