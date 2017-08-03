import numpy as np
import sys
import glob

def find_halo(pathname,x,y,z):
	numfiles = 64
	globalmin = 5
	i = 184
	for j in np.arange(numfiles): #Bring together all the AHF files for the snapshot
		print j
        	temph = glob.glob(pathname+'analysis/ahf_snap%03d/ahf.snap_%03d.%04d.z*.*.AHF_halos'%(i,i,j))
        	temph = str(temph).strip('[]').replace("'","")
        	h = np.genfromtxt(temph)
		dist = np.sqrt((h[:,5]-x)**2+(h[:,6]-y)**2+(h[:,7]-z)**2)
		mindist = min(dist)
		minarg = np.argmin(dist)
		if mindist< globalmin:
			globalmin = mindist
			minfile= j
	print 'min dist is %.3f, halo is in file %d'%(globalmin,minfile)	
		

if __name__ == "__main__":

	x = np.float(sys.argv[1])
	y = np.float(sys.argv[2])
	z = np.float(sys.argv[3])
	pathname = '/nobackup/afitts/Gadget-2.0.7/production/mfm848_14_giz5_12_16_raw_output/'
	find_halo(pathname,x,y,z)
