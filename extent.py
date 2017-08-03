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
extent = np.linspace(121,160,40)
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
	pos = pg.readsnap('/nobackupp8/afitts/Gadget-2.0.7/production/mfm%s_13_giz5_12_16_raw_output/snapdir_%03d/snapshot_%03d.0.hdf5'%(hnum,i,i),'pos','gas')/1000
	mass = pg.readsnap('/nobackupp8/afitts/Gadget-2.0.7/production/mfm%s_13_giz5_12_16_raw_output/snapdir_%03d/snapshot_%03d.0.hdf5'%(hnum,i,i),'mass','gas')*1e10/.71
	hdf = pd.HDFStore('/nobackup/afitts/Gadget-2.0.7/production/mfm12596_13_giz5_12_16_raw_output/analysis/dataframes/halo12596_13_giz5_12_16_snap%03d.h5'%(i))
	print 'bob',i
	red = np.float(hdf['props']['redshift'])
	rvir = np.float(hdf['props']['rvir'])*1000/(red+1)
	dm = hdf['particles/dm']['mass'].as_matrix()
	dx = hdf['particles/dm']['x'].as_matrix()
	dy = hdf['particles/dm']['y'].as_matrix()
	dz = hdf['particles/dm']['z'].as_matrix()
	dp = np.column_stack((dx,dy,dz))
	smm = hdf['particles/star']['mass'].as_matrix()
	sx = hdf['particles/star']['x'].as_matrix()
	sy = hdf['particles/star']['y'].as_matrix()
	sz = hdf['particles/star']['z'].as_matrix()
	sp = np.column_stack((sx,sy,sz))
	dsp = np.vstack((dp,sp))
	dsm = np.append(dm,smm)
	#dsc = mikecm(fname = dsp,nofile=True,pmass=dsm)
	dsc = np.zeros(3)
	dsc[0] = np.float(hdf['props']['halox'])
        dsc[1] = np.float(hdf['props']['haloy'])
        dsc[2] = np.float(hdf['props']['haloz'])
	dist = np.sqrt((pos[:,0]-dsc[0])**2+(pos[:,1]-dsc[1])**2+(pos[:,2]-dsc[2])**2)/.71*1000/(red+1)
	if rank == 184:
		print 'bobba',sum(mass[dist<rvir])
	cumumass = np.cumsum(mass[dist.argsort()])
	dist = np.sort(dist)
	my_A = dist[np.where(cumumass<1e8)[0][-1]]
	my_B = dist[np.where(cumumass<1e7)[0][-1]]
        my_C = dist[np.where(cumumass<1e6)[0][-1]]
	hdf.close()
	print 'my name is', i	
except Exception,e:
	print e
	print 'no gas?'
	my_A = 0
	my_B = 0
	my_C = 0
my_A = comm.gather(my_A, root=0)#comm.Allgather( [my_A, MPI.DOUBLE], [A, MPI.DOUBLE] )
my_B = comm.gather(my_B, root=0)
my_C = comm.gather(my_C, root=0)
if rank == 0:
	np.savetxt('halo%s_extended_extent.txt'%(hnum),np.column_stack((my_C,my_B,my_A)))


