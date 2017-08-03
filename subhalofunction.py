import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as co
from matplotlib import rcParams
import sys
import time
import glob

rcParams['lines.linewidth'] = 4
rcParams['axes.linewidth'] = 2
rcParams['xtick.major.size'] = 10
rcParams['xtick.minor.size'] = 5
rcParams['xtick.major.width'] = 2
rcParams['xtick.minor.width'] = 2
rcParams['ytick.major.size'] = 10
rcParams['ytick.minor.size'] = 5
rcParams['ytick.major.width'] = 2
rcParams['ytick.minor.width'] = 2
rcParams['font.size'] = 20
rcParams['xtick.labelsize']= '16'
rcParams['ytick.labelsize']= '16'
rcParams['savefig.bbox'] = 'tight' 
	

def get_subhalos(pathname,hnum,bincount,dmocor=1,low=1):
	if hnum == '10q' or hnum == '10v':
		i = 600
	else:
		i = 184
	numcores = 16
	switch = 0
	if dmocor != 1:
		typ = 'dmo'
	else:
		typ = ''
	if low != 1:
		res = '_low'
	else:
		res = ''
	### For a subhalo function for each individual host dwarf ###
	#for j in np.arange(numcores): #Bring together all the AHF files for the snapshot
	#	temph = glob.glob(pathname+'ahf_snap%03d/ahf.snap_%03d.%04d.z*.*.AHF_halos'%(i,i,j)) 
	#	temph = str(temph).strip('[]').replace("'","")
	#	h = np.genfromtxt(temph)
	#	if switch == 0:
	#		halo = h
	#		switch = 1
	#	if switch == 1:
	#		try:
	#			halo = np.vstack((halo,h))
	#		except:
	#			print "nothing there"
	#hosthalo = np.genfromtxt(pathname+'halo%smerger_hist.txt'%hnum)
	#hosthalo = hosthalo[0,1]
	#subhalos = halo[halo[:,1]==hosthalo,3]*dmocor

	### For a halo function for the entire run (above a certain contamination threshold) ###

	halo = np.genfromtxt('halo%s_contam95%s_halo_catalogue%s.txt'%(hnum,typ,res),skip_header=1)
	subhalos = halo[halo[:,3]>0,3]
	bins = np.logspace(np.log10(1e3),np.log10(5e7),bincount)
	subhalos ,mass_bins = np.histogram(subhalos,bins=bins)
	mass_bins = 10**(np.log10(mass_bins[1:])-np.log10(bins[1]/bins[0])/2)
	return subhalos, mass_bins
	
class subhalofunction(object):
	def __init__(self):
		self.fig = plt.figure()
		self.sub = self.fig.add_subplot(111)
		my_cmap=plt.get_cmap('plasma')
		self.sm = cm.ScalarMappable(cmap=my_cmap,norm=co.Normalize(vmin=5.004, vmax=7.176))
		self.sm._A = []
		self.sub.set_xlabel(r'Mass $(\mathrm{M}_\odot)$')
		self.sub.set_ylabel('dN/dM (>M$_\star$)')
		self.sub.set_xlim(1e3,1e7)
		self.sub.set_ylim(1e-1,2e0)

	def add_line(self,subhalos,mass_bins,low=0):
		self.sub.loglog(mass_bins, (sum(subhalos)-np.cumsum(subhalos))/14, color='k')
		if low == 1:
			self.sub.loglog(mass_bins,(sum(subhalos)-np.cumsum(subhalos))/14, color='r')

	def add_dmoline(self,subhalos,mass_bins,low=0):
		self.sub.loglog(mass_bins, subhalos, color='k',linestyle='--')
		if low == 1:
			self.sub.loglog(mass_bins, subhalos, color='r',linestyle='--')

	def save(self,hnum,date):
		self.fig.savefig('galaxy%s_function_%s.pdf'%(hnum,date),transparent=True)

def main():
	hnum =['11707','32503','12596','007','848','796','20192','20910','897','948','1016','32257']#,'10q','10v','1084']
	res = ['_13','_13','_13','_11','_13','_13','_13','_13','_13','_13','_13','_13','_11','_11','_13']
	lowres = ['','','','','','','','','','','','','_11','_11','']
	#hnum = sys.argv[1]
	#res = sys.argv[2]
	bincount = 30
	date = time.strftime("%m_%d_%Y")
	subfuncplot = subhalofunction()
	tot_hysub = np.zeros(bincount-1)
	tot_dmosub = np.zeros(bincount-1)
	tot_lowhysub = np.zeros(bincount-1)
	tot_lowdmosub = np.zeros(bincount-1)
	for i in np.arange(len(hnum)):
		#subfuncplot = subhalofunction()
		pathname = '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz5_12_16_raw_output/analysis/'%(hnum[i],res[i])
		hysub,hymass_bin = get_subhalos(pathname,hnum[i],bincount,1,1)
		tot_hysub = tot_hysub+hysub
		pathname = '/nobackup/afitts/Gadget-2.0.7/production/mfm%s%s_giz5_12_16_raw_output/analysis/'%(hnum[i],lowres[i])
		lowhysub,lowhymass_bin = get_subhalos(pathname,hnum[i],bincount,1,0)
		tot_lowhysub = tot_lowhysub+lowhysub
		#pathname = '/nobackup/afitts/Gadget-2.0.7/production/gizdm%s%s_raw_output/analysis/'%(hnum[i],res[i])
		#dmosub,dmomass_bin = get_subhalos(pathname,hnum[i],bincount,0.83120300751,1)
		#tot_dmosub = tot_dmosub+dmosub
		#pathname = '/nobackup/afitts/Gadget-2.0.7/production/gizdm%s%s_raw_output/analysis/'%(hnum[i],lowres[i])
		#lowdmosub,lowdmomass_bin = get_subhalos(pathname,hnum[i],bincount,0.83120300751,0)
		#tot_lowdmosub = tot_lowdmosub+lowdmosub
		print hnum[i]
		#subfuncplot.add_line(hysub,hymass_bin)
		#subfuncplot.add_dmoline(dmosub,dmomass_bin)
		#subfuncplot.add_line(lowhysub,lowhymass_bin,1)
		#subfuncplot.add_dmoline(lowdmosub,lowdmomass_bin,1)
		#subfuncplot.save(hnum[i],date)
	subfuncplot.add_line(tot_hysub,hymass_bin)
	#subfuncplot.add_dmoline(tot_dmosub,dmomass_bin)
	subfuncplot.add_line(tot_lowhysub,lowhymass_bin,1)
	#subfuncplot.add_dmoline(tot_lowdmosub,lowdmomass_bin,1)
	subfuncplot.save('all',date)

main()
