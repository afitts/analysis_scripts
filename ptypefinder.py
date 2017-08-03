import pandas as pd
import sys

snum = int(sys.argv[1])

my_cols = ['id', 'type', 'x', 'y', 'z', 'vx', 'vy', 'vz']

a = pd.read_csv('ahf.snap_184.00%.02d.z0.000.AHF_particles'%snum,names=my_cols, sep='\t')
x = map(int,a['id'][1].split()) #Number of particles in halo
y = a['id'][a['type']==2].index.values #df of lowres dm particles of type 2


print x[0]

print len(y[y<x[0]])

