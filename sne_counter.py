import numpy as np
import sys

snum = np.float(sys.argv[1])
vers = ['','SNe']#['','1','a','b','c','d','e','f','g','h','i','j']
sne = np.zeros(len(vers))

for i in np.arange(len(vers)):
	a = np.genfromtxt('/nobackup/afitts/Gadget-2.0.7/production/mfm1016%s_13_giz5_12_16_raw_output/SNeIIheating.txt'%(vers[i]))
	atime = np.genfromtxt('output_times.txt')
	sne[i] = sum(a[(a[:,0]>atime[snum-1])&(a[:,0]<atime[snum]),3])
	print vers[i],sne[i]

print sne
