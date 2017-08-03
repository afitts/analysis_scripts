import numpy as np
import pandas as pd

def subhalo_counter(hnum,res,snap,f):
	global tot8,tot7,vis8,vis7,gas8,gas7
	subs = np.genfromtxt('halo%s_contam95_halo_catalogue.txt'%hnum,skip_header=1)
	hdf = pd.HDFStore('/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz5_12_16_raw_output/analysis/dataframes/halo%s%s_giz5_12_16_snap%s.h5'%(hnum,res,hnum,res,snap))
	rvir = np.float(hdf['props']['rvir'])*1e3*0.71
	print subs[(subs[:,9]<rvir)&(subs[:,9]>.1*rvir)&(subs[:,3]!=0),2]
	tot8 = np.append(tot8,len(subs[(subs[:,9]<rvir)&(subs[:,9]>.1*rvir)&(subs[:,2]>1e8)]))
	tot7 = np.append(tot7,len(subs[(subs[:,9]<rvir)&(subs[:,9]>.1*rvir)&(subs[:,2]>1e7)]))
	vis8 = np.append(vis8,len(subs[(subs[:,9]<rvir)&(subs[:,9]>.1*rvir)&(subs[:,2]>1e8)&(subs[:,3]!=0)]))
	vis7 = np.append(vis7,len(subs[(subs[:,9]<rvir)&(subs[:,9]>.1*rvir)&(subs[:,3]!=0)]))
        gas8 = np.append(gas8,len(subs[(subs[:,9]<rvir)&(subs[:,9]>.1*rvir)&(subs[:,2]>1e8)&(subs[:,4]!=0)]))
        gas7 = np.append(gas7,len(subs[(subs[:,9]<rvir)&(subs[:,9]>.1*rvir)&(subs[:,2]>1e7)&(subs[:,4]!=0)]))
	hdf.close()
	

f = open('subhalo_12_count.txt','w')
f.write('(0) Halo #	(1) # of subhalos (M_halo>1e8)	(2) # of subhalos (M_halo>1e7)	(3) # of galaxies (M_halo>1e8)	(4) # of galaxies (M_halo>1e7)	(3) # of gaseous halos (M_halo>1e8)  (4) # of gaseous halos (M_halo>1e7)\n')
hnum =['11707','32503','12596','007','848','796','20192','20910','897','948','1016','32257','10q','10v','1084']#
res = ['_13','_13','_13','_11','_13','_13','_13','_13','_13','_13','_13','_13','_11','_11','_13']
ver = ['5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16','5_12_16']
dmo = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
snum = [185,185,185,185,185,185,185,185,185,185,185,185,601,601,185]
tot8 = np.array([])
tot7 = np.array([])
vis8 = np.array([])
vis7 = np.array([])
gas8 = np.array([])
gas7 = np.array([])
for i in np.arange(len(hnum)):
	print hnum[i]
	subhalo_counter(hnum[i],res[i],snum[i]-1,f)
print len(hnum),len(tot8),sum(vis7)
tbl = np.column_stack((hnum,tot8,tot7,vis8,vis7,gas8,gas7))
#np.savetxt(f,tbl,fmt = '%s',delimiter='\t')
f.close()
