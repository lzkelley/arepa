import plotting.images.dm_images as dm_images
import matplotlib.pyplot as plt
import simread.readsubfHDF5 as subf
import numpy as np

npixels=720

for run_base in ['dm_runs']:
    for snapnum in range(7):
        cat = subf.subfind_catalog( '../'+run_base+'/', snapnum, keysel=['GroupPos'] )
        center = cat.GroupPos[0,:]
        image = dm_images.dm_image( '../dm_runs/output/', snapnum, center=center, xrange=[-5.0,5.0], yrange=[-5.0,5.0], pixels=npixels , cosmo_wrap=True)
        
        print np.min(image), np.min(image[np.isfinite(image)]), np.max(image), np.median(image) 
        
        fig,ax=plt.subplots( figsize=(1,1) )
        ax.imshow(image, vmin=0.0, vmax=4.5)
        fig.subplots_adjust( left=0.0, bottom=0.0, top=1.0, right=1.0)
        ax.axis('off')
        
        filename=run_base+'_dm_box_image_snap_{:.0f}_{:.0f}.png'.format( snapnum, npixels)
        print filename
        
        fig.savefig( './plots/'+filename, dpi=npixels )



