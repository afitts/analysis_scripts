import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from annika_shapes import axis
from mikecm import mikecm
import sys


def main_worker(prefix='',f='',DMO=0,rad_max=0.0001,iter_max=0.01):	
	hdf = pd.HDFStore(prefix+f)
	if DMO == 0:
		rhalf = np.float(hdf['props']['rhalf'])
		rvir = np.float(hdf['props']['rvir'])#not z corrected, but h factor
		smass = hdf['particles/star']['mass'].as_matrix()
		d_r = hdf['particles/dm']['r'].as_matrix()
		s_r = hdf['particles/star']['r'].as_matrix()
		xc = np.float(hdf['props']['halox'])
		yc = np.float(hdf['props']['haloy'])
		zc = np.float(hdf['props']['haloz'])
		dx = hdf['particles/dm']['x'].as_matrix()
		dy = hdf['particles/dm']['y'].as_matrix()
		dz = hdf['particles/dm']['z'].as_matrix()
		dp = np.column_stack((dx,dy,dz))
		dm = hdf['particles/dm']['mass'].as_matrix()
		sx = hdf['particles/star']['x'].as_matrix()
		sy = hdf['particles/star']['y'].as_matrix()
		sz = hdf['particles/star']['z'].as_matrix()
		smm = hdf['particles/star']['mass']
		sp = np.column_stack((sx,sy,sz))
		dsp = np.vstack((dp,sp))
		dsm = np.append(dm,smm)
		dsc = mikecm(fname = sp,nofile=True,pmass=smm)
		x = hdf['particles/dm']['x'].as_matrix()-dsc[0]
		y = hdf['particles/dm']['y'].as_matrix()-dsc[1]
		z = hdf['particles/dm']['z'].as_matrix()-dsc[2]
		#x = hdf['particles/dm']['x'].as_matrix()-xc
		#y = hdf['particles/dm']['y'].as_matrix()-yc
		#z = hdf['particles/dm']['z'].as_matrix()-zc
		arr_in = np.vstack((x/.71,y/.71,z/.71)).T
		daxratios = axis(arr_in, rad_max, iter_max,shell=False, axes_out=False, fix_volume=True, quiet=False)
		x = hdf['particles/star']['x'].as_matrix()-dsc[0]
		y = hdf['particles/star']['y'].as_matrix()-dsc[1]
		z = hdf['particles/star']['z'].as_matrix()-dsc[2]
		#x = hdf['particles/star']['x'].as_matrix()-xc
		#y = hdf['particles/star']['y'].as_matrix()-yc
		#z = hdf['particles/star']['z'].as_matrix()-zc
		arr_in = np.vstack((x/.71,y/.71,z/.71)).T
		saxratios = axis(arr_in, rad_max, iter_max,shell=False, axes_out=False, fix_volume=True, quiet=False)
		dvx = hdf['particles/dm']['vx'].as_matrix()
		dvy = hdf['particles/dm']['vy'].as_matrix()
		dvz = hdf['particles/dm']['vz'].as_matrix()
		dvdisp_x = np.sqrt(np.std(dvx[d_r<rad_max])**2)
		dvdisp_y = np.sqrt(np.std(dvy[d_r<rad_max])**2)
		dvdisp_z = np.sqrt(np.std(dvz[d_r<rad_max])**2)
		svx = hdf['particles/star']['vx'].as_matrix()
		svy = hdf['particles/star']['vy'].as_matrix()
		svz = hdf['particles/star']['vz'].as_matrix()
		svdisp_x = np.sqrt(np.std(svx[s_r<rad_max])**2)
		svdisp_y = np.sqrt(np.std(svy[s_r<rad_max])**2)
		svdisp_z = np.sqrt(np.std(svz[s_r<rad_max])**2)
		dm_enc = sum(dm[d_r<rad_max])
		sm_enc = sum(smm[s_r<rad_max])
		if rad_max == 0.001:
			np.savetxt('Halo848_star_pos_within_5kpc.txt',np.column_stack((x[s_r<.005]*1000/.71,y[s_r<.005]*1000/.71,z[s_r<.005]*1000/.71)),header = '(0) X (kpc) (1) Y (2) Z')
	else:
		print "Wrong file!!!"
		mstar = 0
		rad_max = 0
		axratios = [0,0,0]
	hdf.close() 
	return rad_max,daxratios[0],daxratios[1],saxratios[0],saxratios[1],dm_enc,sm_enc,dvdisp_x,dvdisp_y,dvdisp_z,svdisp_x,svdisp_y,svdisp_z

iter_max = 0.1#np.float(sys.argv[1])
hnum = ['10q','10v']#['11707','32503','12596','007','848','796','897','948','1016','20192','20910','32257','10q','10v']
res = ['_11','_11']#['_13','_13','_13','_11','_13','_13','_13','_13','_13','_13','_13','_13','_11','_11','_13']
snum = [600,600]#[184,184,184,184,184,184,184,184,184,184,184,184,600,600]
rad = np.array([0.00025,0.0005,0.001])
for i in np.arange(len(hnum)):
	daxratios = np.zeros((len(rad),2))
	saxratios = np.zeros((len(rad),2))
	rad_max = np.zeros(len(rad))
	dm_enc = np.zeros(len(rad))
	sm_enc = np.zeros(len(rad))
	dvdisp_x = np.zeros(len(rad))
	dvdisp_y = np.zeros(len(rad))
	dvdisp_z = np.zeros(len(rad))
	svdisp_x = np.zeros(len(rad))
	svdisp_y = np.zeros(len(rad))
	svdisp_z = np.zeros(len(rad))
	for j in np.arange(len(rad)):
		rad_max[j],daxratios[j,0],daxratios[j,1],saxratios[j,0],saxratios[j,1],dm_enc[j],sm_enc[j],dvdisp_x[j],dvdisp_y[j],dvdisp_z[j],svdisp_x[j],svdisp_y[j],svdisp_z[j]=main_worker('/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz5_12_16_raw_output/analysis/dataframes/'%(hnum[i],res[i]),'halo%s%s_giz5_12_16_snap%03d.h5'%(hnum[i],res[i],snum[i]),0,rad[j],iter_max)
	np.savetxt('Halo%s_Info.txt'%hnum[i],np.column_stack((rad_max*1e6,daxratios[:,1],daxratios[:,0],saxratios[:,1],saxratios[:,0],dm_enc,sm_enc,dvdisp_x,dvdisp_y,dvdisp_z,svdisp_x,svdisp_y,svdisp_z)),fmt = '%.3e',header = '(0) Radius (pc) (1) DM b/a (2) DM c/a (3) Stellar b/a (4) Stellar c/a (4) DM mass enclosed (M_sun) (5) Stellar mass enclosed (6) DM vel disp x (km/s) (7) DM vel disp y (8) DM vel disp z (9) Stellar vel disp x (10) Stellar vel disp y (11) Stellar vel disp z')
