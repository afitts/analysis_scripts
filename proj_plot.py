import yt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import yt.units as units
from numpy import array
import sys

"""                                                                                                                                                    
    This script plots the 2D projections over x,y or z axis 
"""

usage = "usage:  python {0} <output directory>[zoom boxsize in kpccm]".format(sys.argv[0].split('/')[-1])

if len(sys.argv) < 2:
       print usage
       sys.exit(1337)
###using wmap7 cosmology
h=.71
outdir=sys.argv[1]
if not outdir.endswith('/'):    outdir = outdir + '/'

if len(sys.argv)<3:
       zoombox=25000#183.168*2 #rvir in kpccm/h 
else:
       zoombox=2*float(sys.argv[2])*h  #in kpccm/h 
###snapdir
#fns=['/export/home/vrobless/rouge/sidm_mw_disk_zoom/GVD_Z13_nodisk_833/snaps/snapdir_152/snapshot_152.0.hdf5']
fns=['/nobackup/afitts/Gadget-2.0.7/production/mfm848_14_giz5_12_16_raw_output/snapdir_026/snapshot_026.0.hdf5']

#### MW rockstar 833 z=0
mw1cen_cdm=array([12642.68096457, 12133.94379679, 12446.52161890])#[28.0590859375,24.75334765625,25.03334765625])
#mw1cen_cdm=array([28.052732,25.140259,24.641251]) #cdm l13 rockstar mpc/h l13 833 
cen= [mw1cen_cdm]
label=['z=4.9657']

#axes_pad = 0.00
#p.set_font_size(28)
#p.annotate_text([215,270],label[i], coord_system='plot',text_args={'color':'white','fontsize':'50'})
#plt.savefig(outdir+'sidm_proj_833_l13_300.png',facecolor='k',edgecolor=None,transparent=True,bbox_inches='tight',pad_inches=0.0)

fig = plt.figure()
unit_base = {'UnitLength_in_cm': 3.08568e+21,'UnitMass_in_g':1.989e+43,'UnitVelocity_in_cm_per_s' : 100000}
# See http://matplotlib.org/mpl_toolkits/axes_grid/api/axes_grid_api.html
# These choices of keyword arguments produce a four panel plot with a single
# shared narrow colorbar on the right hand side of the multipanel plot. Axes
# labels are drawn for all plots since we're slicing along different directions
# for each plot.
#grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                #nrows_ncols = (1, 1),
                #axes_pad = 0.00,
                #label_mode = "L",
                #share_all = True,
                #cbar_location="right",
                #cbar_mode="single",
                #cbar_size="3%",
                #cbar_pad="0%")

fields=[('PartType0','Density')]
for i, fn in enumerate(fns):
    # Load the data and create a single plot
    ds = yt.load(fn,unit_base=unit_base,n_ref=4,over_refine_factor=2,index_ptype="PartType0") # load data
    ad = ds.all_data()
    #center=ds.arr(cen[i],'kpccm/h') 
    new_box_size = ds.quan(zoombox,'kpccm/h')
    sp=ds.sphere(mw1cen_cdm,6)
    #sp.add_field(('gas','particle_position_x'),function=sp['gas','x'],units='kpccm/h')
    #sp.add_field(('gas','particle_position_y'),function=sp['gas','y'],units='kpccm/h')
    if (i==0)|(i==1):
        #p = yt.ParticleProjectionPlot(ds, 'z',center=mw1cen_cdm,width=6,data_source=sp,weight_field=('deposit', 'PartType0_smoothed_density'))
        #p.set_cmap(field=('gas','density'), cmap='bone_r') 
        p = yt.ProjectionPlot(ds,'z',('deposit','PartType0_smoothed_density'),data_source=sp,center=mw1cen_cdm,width=6)
        #p.hide_colorbar()
        #p.hide_axes()

    # Ensure the colorbar limits match for all plots
        #p.set_unit('particle_mass', 'Msun')
        #p.set_zlim('particle_mass', 5e4, 2e9) 
        #p.annotate_scale(draw_inset_box=False)
        #p.set_font_size(28)
        #p.annotate_text([215,270],label[i], coord_system='plot',text_args={'color':'white','fontsize':'50'})
        #p.set_axes_unit('kpc')
 
        # This forces the ProjectionPlot to redraw itself on the AxesGrid axes.                                                                               
    #plot = p.plots[('PartType0','Density')]
    #plot = p.plots[fields[i]]
    #plot.figure = fig
    #plot.axes = grid[i].axes
    #plot.cax = grid.cbar_axes[i]
    #plot.axes.patch.set_facecolor('white') 

    # Finally, this actually redraws the plot.
    #p._setup_plots()
print 'hi'
p.save(outdir+'yt2dhist_proj_848_z5_7_17_2017.pdf')
plt.savefig(outdir+'yt2dhist_proj_848_z5.pdf')#,transparent=True,bbox_inches='tight',pad_inches=0.0) 
#plt.savefig(outdir+'sidm_proj_833_l13_300.png',facecolor='k',edgecolor=None,transparent=True,bbox_inches='tight',pad_inches=0.0)
