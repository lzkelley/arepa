import plotting.images.gas_images as gas_images
import matplotlib.pyplot as plt
import simread.readsubfHDF5 as subf
import numpy as np
import glob


npixels=720

for run_base in [#'explicit_feedback_256', \
                 'explicit_feedback_256_soft', \
#                 'explicit_feedback_512', \
                 'explicit_feedback_512_soft_amd']:         #gas_prod1_runs', 'gas_cool_runs', 'gas_adiab_runs', 'dm_runs']:
    print "\n\n\n"
    print '../'+run_base+'/output/snap*'
    snaplist = glob.glob( '../'+run_base+'/output/snap*')
    n_snaps = len(snaplist)
    print snaplist
    for snapnum in [10]:        #range(n_snaps):        #range(n_snaps):        #range(13):
        cat = subf.subfind_catalog( '../'+run_base+'/', snapnum, keysel=['GroupPos', 'GroupLenType'] )
        for groupnr in range(np.min([10, len(cat.GroupLenType[:,0])] ) ):
            if cat.GroupLenType[groupnr,4] > 500:
                center = cat.GroupPos[groupnr,:]

                fig_,(ax1_,ax2_,ax3_)=plt.subplots( 1, 3, figsize=(3,1) )
    
                image,massmap = gas_images.sfr_image( dir = '../'+run_base+'/output/', snapnum = snapnum, center=center, xrange=[-0.01,0.01], yrange=[-0.01,0.01], pixels=npixels , cosmo_wrap=True, massmap=True, projaxis=0)
                print "\n\n\n"
                print image.shape
                print np.min(image), np.min(image[np.isfinite(image)]), np.max(image), np.median(image) 
                print "\n\n\n"
                ax1_.imshow(image, vmin=0.25, vmax=6.5 )



                image,massmap = gas_images.sfr_image( dir = '../'+run_base+'/output/', snapnum = snapnum, center=center, xrange=[-0.01,0.01], yrange=[-0.01,0.01], pixels=npixels , cosmo_wrap=True, massmap=True, projaxis=1)
                print "\n\n\n"
                print image.shape
                print np.min(image), np.min(image[np.isfinite(image)]), np.max(image), np.median(image) 
                print "\n\n\n"
                ax2_.imshow(image, vmin=0.25, vmax=6.5 )


                image,massmap = gas_images.sfr_image( dir = '../'+run_base+'/output/', snapnum = snapnum, center=center, xrange=[-0.01,0.01], yrange=[-0.01,0.01], pixels=npixels , cosmo_wrap=True, massmap=True, projaxis=2)
                print "\n\n\n"
                print image.shape
                print np.min(image), np.min(image[np.isfinite(image)]), np.max(image), np.median(image) 
                print "\n\n\n"
                ax3_.imshow(image, vmin=0.25, vmax=6.5 )

                fig_.subplots_adjust( left=0.0, bottom=0.0, top=1.0, right=1.0, hspace=0, wspace=0)
                ax1_.axis('off')
                ax2_.axis('off')
                ax3_.axis('off')
                

                filename='{:s}_sfr_central_image_snap_{:.0f}_group_{:.0f}_res_{:.0f}.png'.format( run_base, snapnum, groupnr, npixels)
                print filename
                
                fig_.savefig( './plots/'+filename, dpi=npixels )




