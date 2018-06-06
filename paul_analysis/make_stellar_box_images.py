import plotting.images.stellar_images as stellar_images
import matplotlib.pyplot as plt
import simread.readsubfHDF5 as subf
import numpy as np
import glob

npixels=720

for run_base in [#'explicit_feedback_256', \
                 'explicit_feedback_256_soft', \
#                 'explicit_feedback_256_soft_old_dens', \
#                 'explicit_feedback_512', \
                 #'explicit_feedback_512_soft_amd']:         #gas_prod1_runs', 'gas_cool_runs', 'gas_adiab_runs', 'dm_runs']:
                ]:
    print "\n\n\n"
    print '../'+run_base+'/output/snap*'
    snaplist = glob.glob( '../'+run_base+'/output/snap*')
    n_snaps = len(snaplist)
    print snaplist
    for snapnum in [11]:        #range(n_snaps):        #range(13):
        cat = subf.subfind_catalog( '../'+run_base+'/', snapnum, keysel=['GroupPos', 'GroupLenType'] )

        for groupnr in range(np.min([10, len(cat.GroupLenType[:,0])] ) ):
          if cat.GroupLenType[groupnr,4] > 500:
            center = cat.GroupPos[groupnr,:]
            print center



            fig,(ax1,ax2,ax3)=plt.subplots( 1, 3, figsize=(3,1) )

            image  = stellar_images.stellar_image( '../'+run_base+'/output/', snapnum, center=center, \
                                                           xrange=[-0.01,0.01], yrange=[-0.01,0.01], pixels=npixels , \
                                                           cosmo_wrap=True, massmap=True, dust=True, cosmo=True,\
                                                           maxden=1.5e6, dynrange=1e3, projaxis=0, unit_length_mpc=True)

            ax1.imshow(image  )

            image  = stellar_images.stellar_image( '../'+run_base+'/output/', snapnum, center=center, \
                                                           xrange=[-0.01,0.01], yrange=[-0.01,0.01], pixels=npixels , \
                                                           cosmo_wrap=True, massmap=True, dust=True, cosmo=True,\
                                                           maxden=1.5e6, dynrange=1e3, projaxis=1, unit_length_mpc=True)
            ax2.imshow(image  )

            image  = stellar_images.stellar_image( '../'+run_base+'/output/', snapnum, center=center, \
                                                           xrange=[-0.01,0.01], yrange=[-0.01,0.01], pixels=npixels , \
                                                           cosmo_wrap=True, massmap=True, dust=True, cosmo=True,\
                                                           maxden=1.5e6, dynrange=1e3, projaxis=2, unit_length_mpc=True)
            ax3.imshow(image  )

            fig.subplots_adjust( left=0.0, bottom=0.0, top=1.0, right=1.0, wspace=0,hspace=0)
            for ax in [ax1,ax2,ax3]:
                ax.axis('off')

            filename='{:s}_stellar_central_image_WITH_DUST_snap_{:.0f}_group_{:.0f}_res_{:.0f}.png'.format( run_base, snapnum, groupnr, npixels)
            print filename

            fig.savefig( './plots/'+filename, dpi=npixels )






            if True:
                fig,(ax1,ax2,ax3)=plt.subplots( 1, 3, figsize=(3,1) )
    
                image  = stellar_images.stellar_image( '../'+run_base+'/output/', snapnum, center=center, \
                                                               xrange=[-0.01,0.01], yrange=[-0.01,0.01], pixels=npixels , \
                                                               cosmo_wrap=True, massmap=True, dust=False, cosmo=True,\
                                                               maxden=1.5e6, dynrange=1e3, projaxis=0)
    
                ax1.imshow(image  )
    
                image  = stellar_images.stellar_image( '../'+run_base+'/output/', snapnum, center=center, \
                                                               xrange=[-0.01,0.01], yrange=[-0.01,0.01], pixels=npixels , \
                                                               cosmo_wrap=True, massmap=True, dust=False, cosmo=True,\
                                                               maxden=1.5e6, dynrange=1e3, projaxis=1)
                ax2.imshow(image  )
    
                image  = stellar_images.stellar_image( '../'+run_base+'/output/', snapnum, center=center, \
                                                               xrange=[-0.01,0.01], yrange=[-0.01,0.01], pixels=npixels , \
                                                               cosmo_wrap=True, massmap=True, dust=False, cosmo=True,\
                                                               maxden=1.5e6, dynrange=1e3, projaxis=2)
                ax3.imshow(image  )
    
                fig.subplots_adjust( left=0.0, bottom=0.0, top=1.0, right=1.0, wspace=0,hspace=0)
                for ax in [ax1,ax2,ax3]:
                    ax.axis('off')
    
                filename='{:s}_stellar_central_image_snap_{:.0f}_group_{:.0f}_res_{:.0f}.png'.format( run_base, snapnum, groupnr, npixels)
                print filename
    
                fig.savefig( './plots/'+filename, dpi=npixels )












